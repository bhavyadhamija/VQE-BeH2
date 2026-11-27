"""
qubits.py — from molecular integrals to a qubit Hamiltonian, exact diagonalisation and VQE.

  * frozen-core active space
  * Jordan–Wigner mapping to Pauli strings (spin-orbital k -> qubit k, interleaved alpha/beta)
  * sparse qubit matrices, exact (FCI) ground state in the right particle-number sector
  * UCCSD ansatz as a product of exponentials, optimised with exact adjoint gradients

Pauli strings are stored in "XZ form": an integer pair (x, z) meaning  prod_k X_k^x_k Z_k^z_k.
"""
import itertools
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import expm_multiply, eigsh
from scipy.optimize import minimize


def popcount(a):
    a = np.asarray(a, dtype=np.int64)
    c = np.zeros_like(a)
    while np.any(a):
        c += a & 1
        a = a >> 1
    return c


def _pc(n):
    return bin(n).count("1")


# ---------------------------------------------------------------- Pauli algebra
def op_mul(A, B):
    out = {}
    for (x1, z1), c1 in A.items():
        for (x2, z2), c2 in B.items():
            sign = -1 if _pc(z1 & x2) % 2 else 1          # Z^z1 X^x2 = (-1)^{z1.x2} X^x2 Z^z1
            key = (x1 ^ x2, z1 ^ z2)
            out[key] = out.get(key, 0) + sign * c1 * c2
    return out


def op_add(A, B, s=1.0):
    out = dict(A)
    for k, v in B.items():
        out[k] = out.get(k, 0) + s * v
    return out


def creation(j):
    low = (1 << j) - 1
    return {(1 << j, low): 0.5, (1 << j, low | (1 << j)): 0.5}      # Z_<j X_j (1+Z_j)/2


def annihilation(j):
    low = (1 << j) - 1
    return {(1 << j, low): 0.5, (1 << j, low | (1 << j)): -0.5}     # Z_<j X_j (1-Z_j)/2


def clean(A, tol=1e-12):
    return {k: v for k, v in A.items() if abs(v) > tol}


def to_labels(A, n):
    """Convert XZ form to standard Pauli labels (leftmost char = highest qubit) with complex coeffs."""
    out = {}
    for (x, z), c in A.items():
        lab, coef = [], c
        for q in range(n - 1, -1, -1):
            xb, zb = (x >> q) & 1, (z >> q) & 1
            if xb and zb:
                lab.append("Y"); coef *= -1j                  # XZ = -iY
            elif xb:
                lab.append("X")
            elif zb:
                lab.append("Z")
            else:
                lab.append("I")
        out["".join(lab)] = out.get("".join(lab), 0) + coef
    return {k: v for k, v in out.items() if abs(v) > 1e-10}


def to_sparse(A, n):
    """Sparse 2^n x 2^n matrix of an XZ-form operator: <b^x| X^x Z^z |b> = (-1)^{z.b}."""
    dim = 1 << n
    b = np.arange(dim, dtype=np.int64)
    by_x = {}
    for (x, z), c in A.items():
        by_x.setdefault(x, []).append((z, c))
    rows, cols, vals = [], [], []
    for x, terms in by_x.items():
        v = np.zeros(dim, dtype=complex)
        for z, c in terms:
            v += c * (1 - 2 * (popcount(z & b) % 2))
        nz = np.abs(v) > 1e-14
        rows.append((b ^ x)[nz]); cols.append(b[nz]); vals.append(v[nz])
    M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(dim, dim))
    if np.max(np.abs(M.data.imag)) < 1e-12:
        M = sp.csr_matrix((M.data.real, M.indices, M.indptr), shape=M.shape)
    return M


# ---------------------------------------------------------------- active space
def active_space_integrals(h_mo, eri_mo, e_nuc, core, active):
    """Frozen-core energy and effective integrals (chemists' notation) for the active orbitals."""
    e_core = e_nuc
    for c in core:
        e_core += 2 * h_mo[c, c]
    for c in core:
        for d in core:
            e_core += 2 * eri_mo[c, c, d, d] - eri_mo[c, d, d, c]
    h_eff = h_mo.copy()
    for c in core:
        h_eff += 2 * eri_mo[:, :, c, c] - eri_mo[:, c, c, :]
    A = np.ix_(active, active)
    return e_core, h_eff[A], eri_mo[np.ix_(active, active, active, active)]


def fermion_hamiltonian(e_core, h, g):
    """H = e_core + sum h_pq a+_ps a_qs + 1/2 sum (pq|rs) a+_ps a+_rt a_st a_qs  (JW, XZ form)."""
    n = h.shape[0]
    nq = 2 * n
    cr = [creation(k) for k in range(nq)]
    an = [annihilation(k) for k in range(nq)]
    H = {(0, 0): e_core}
    for p, q in itertools.product(range(n), repeat=2):
        if abs(h[p, q]) < 1e-12:
            continue
        for s in range(2):
            H = op_add(H, op_mul(cr[2 * p + s], an[2 * q + s]), h[p, q])
    for p, q, r, s_ in itertools.product(range(n), repeat=4):
        v = g[p, q, r, s_]
        if abs(v) < 1e-12:
            continue
        for s, t in itertools.product(range(2), repeat=2):
            P, Q, R, S = 2 * p + s, 2 * q + s, 2 * r + t, 2 * s_ + t
            if P == R or Q == S:
                continue
            term = op_mul(op_mul(cr[P], cr[R]), op_mul(an[S], an[Q]))
            H = op_add(H, term, 0.5 * v)
    return clean(H), nq


def one_body_operator(t, e0=0.0):
    n = t.shape[0]
    O = {(0, 0): e0}
    for p, q in itertools.product(range(n), repeat=2):
        for s in range(2):
            O = op_add(O, op_mul(creation(2 * p + s), annihilation(2 * q + s)), t[p, q])
    return clean(O)


def sector(nq, n_alpha, n_beta):
    b = np.arange(1 << nq, dtype=np.int64)
    amask = sum(1 << k for k in range(0, nq, 2)); bmask = sum(1 << k for k in range(1, nq, 2))
    return np.where((popcount(b & amask) == n_alpha) & (popcount(b & bmask) == n_beta))[0]


def exact_ground_state(Hs, idx):
    sub = Hs[idx][:, idx].toarray()
    w, v = np.linalg.eigh(sub)
    psi = np.zeros(Hs.shape[0]); psi[idx] = v[:, 0]
    return w[0], psi


def qubitwise_groups(labels):
    groups = []
    for lab in labels:
        for g in groups:
            if all(a == b or a == "I" or b == "I" for other in g for a, b in zip(lab, other)):
                g.append(lab); break
        else:
            groups.append([lab])
    return groups


# ---------------------------------------------------------------- UCCSD
def uccsd_generators(n_spatial, n_occ_spatial):
    """Spin-conserving single and double excitations from the closed-shell reference.
    Returns list of (label, sparse anti-Hermitian generator T - T^dagger)."""
    nq = 2 * n_spatial
    occ = list(range(2 * n_occ_spatial)); vir = list(range(2 * n_occ_spatial, nq))
    spin = lambda k: k % 2
    cr = [creation(k) for k in range(nq)]; an = [annihilation(k) for k in range(nq)]
    gens = []
    for i, j in itertools.combinations(occ, 2):
        for a, b in itertools.combinations(vir, 2):
            if sorted([spin(i), spin(j)]) != sorted([spin(a), spin(b)]):
                continue
            T = op_mul(op_mul(cr[a], cr[b]), op_mul(an[j], an[i]))
            G = clean(op_add(T, _dagger(T), -1.0))
            if G:
                gens.append((f"d {i},{j}->{a},{b}", to_sparse(G, nq).real.tocsr()))
    for i in occ:
        for a in vir:
            if spin(i) != spin(a):
                continue
            T = op_mul(cr[a], an[i])
            G = clean(op_add(T, _dagger(T), -1.0))
            gens.append((f"s {i}->{a}", to_sparse(G, nq).real.tocsr()))
    return gens


def _dagger(A):
    """Hermitian conjugate in XZ form: (X^x Z^z)^dag = Z^z X^x = (-1)^{x.z} X^x Z^z."""
    return {(x, z): np.conj(c) * (-1 if _pc(x & z) % 2 else 1) for (x, z), c in A.items()}


def uccsd_state(theta, gens, ref):
    psi = ref.copy()
    for t, (_, G) in zip(theta, gens):
        psi = expm_multiply(t * G, psi)
    return psi


def vqe_uccsd(Hs, gens, ref, theta0=None, tol=1e-10, maxiter=500):
    """Minimise <psi(theta)|H|psi(theta)>, with the exact gradient from one forward and one backward sweep."""
    K = len(gens)
    history = []

    def energy_grad(theta):
        psi = uccsd_state(theta, gens, ref)
        lam = Hs @ psi
        E = float(psi @ lam)
        g = np.zeros(K)
        phi = psi.copy()
        for k in range(K - 1, -1, -1):
            G = gens[k][1]
            g[k] = 2 * float(lam @ (G @ phi))
            phi = expm_multiply(-theta[k] * G, phi)
            lam = expm_multiply(-theta[k] * G, lam)
        history.append(E)
        return E, g

    theta0 = np.zeros(K) if theta0 is None else theta0
    res = minimize(energy_grad, theta0, jac=True, method="L-BFGS-B", options={"maxiter": maxiter, "ftol": tol, "gtol": 1e-8})
    return res, history
