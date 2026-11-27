"""
qchem.py — a small, self-contained quantum-chemistry engine in NumPy.

  * STO-3G basis for H, Be, O (Hehre, Stewart & Pople 1969)
  * one- and two-electron integrals with the McMurchie–Davidson scheme
  * restricted Hartree–Fock (RHF) with DIIS

Everything is in atomic units (bohr, hartree).
"""
import itertools
import math
import numpy as np
from scipy.special import gammainc, gamma as gamma_fn

ANGSTROM = 1.0 / 0.529177210903          # bohr per angstrom
Z = {"H": 1, "Be": 4, "O": 8}

# ---------------------------------------------------------------------------
# STO-3G: contraction for a Slater orbital with exponent 1, scaled by zeta^2
# ---------------------------------------------------------------------------
_1S = ([2.227660584, 0.405771156, 0.109818036], [0.154328967, 0.535328142, 0.444634542])
_2SP = ([0.994202990, 0.231031333, 0.075138568],
        [-0.099967230, 0.399512826, 0.700115469],      # 2s
        [0.155916275, 0.607683719, 0.391957393])       # 2p
ZETA = {"H": {"1s": 1.24}, "Be": {"1s": 3.68, "2sp": 1.15}, "O": {"1s": 7.66, "2sp": 2.25}}


class BasisFunction:
    """Contracted Cartesian Gaussian  N * x^l y^m z^n * sum_i d_i exp(-a_i r^2)."""

    def __init__(self, origin, shell, exps, coefs, atom, label):
        self.origin = np.asarray(origin, float)
        self.shell = shell
        self.exps = np.asarray(exps, float)
        self.coefs = np.asarray(coefs, float)
        self.atom = atom
        self.label = label
        l, m, n = shell
        L = l + m + n
        df = lambda k: 1 if k <= 0 else k * df(k - 2)
        self.norm = (2 / np.pi) ** 0.75 * 2 ** L * self.exps ** ((2 * L + 3) / 4) / math.sqrt(
            df(2 * l - 1) * df(2 * m - 1) * df(2 * n - 1))
        # normalise the contraction
        s = 0.0
        for i, j in itertools.product(range(len(self.exps)), repeat=2):
            s += (self.norm[i] * self.norm[j] * self.coefs[i] * self.coefs[j]
                  * _overlap_prim(self.exps[i], shell, self.origin, self.exps[j], shell, self.origin))
        self.coefs = self.coefs / math.sqrt(s)


def sto3g(atoms, scale=1.0):
    """atoms = [(symbol, xyz_bohr), ...]; scale multiplies every exponent by scale**2
    (used for the virial-scaling experiment). Returns a list of BasisFunction."""
    basis = []
    for ia, (sym, xyz) in enumerate(atoms):
        z1 = ZETA[sym]["1s"]
        basis.append(BasisFunction(xyz, (0, 0, 0), np.array(_1S[0]) * z1 ** 2 * scale ** 2, _1S[1], ia, f"{sym}{ia} 1s"))
        if "2sp" in ZETA[sym]:
            z2 = ZETA[sym]["2sp"]
            e = np.array(_2SP[0]) * z2 ** 2 * scale ** 2
            basis.append(BasisFunction(xyz, (0, 0, 0), e, _2SP[1], ia, f"{sym}{ia} 2s"))
            for k, ax in enumerate("xyz"):
                sh = [0, 0, 0]; sh[k] = 1
                basis.append(BasisFunction(xyz, tuple(sh), e, _2SP[2], ia, f"{sym}{ia} 2p{ax}"))
    return basis


# ---------------------------------------------------------------------------
# McMurchie–Davidson building blocks
# ---------------------------------------------------------------------------
def boys(n, T):
    if T < 1e-12:
        return 1.0 / (2 * n + 1)
    return gamma_fn(n + 0.5) * gammainc(n + 0.5, T) / (2.0 * T ** (n + 0.5))


def hermite_E(i, j, Qx, a, b):
    """Table E[t] for t = 0..i+j  (Hermite expansion of a 1-D Gaussian product)."""
    p = a + b; q = a * b / p
    XPA = -b / p * Qx; XPB = a / p * Qx           # Qx = A - B
    Etab = np.zeros((i + 1, j + 1, i + j + 2))
    Etab[0, 0, 0] = math.exp(-q * Qx * Qx)
    for ii in range(i + 1):
        for jj in range(j + 1):
            if ii == 0 and jj == 0:
                continue
            for t in range(ii + jj + 1):
                if ii > 0:
                    v = XPA * Etab[ii - 1, jj, t]
                    if t > 0: v += Etab[ii - 1, jj, t - 1] / (2 * p)
                    v += (t + 1) * Etab[ii - 1, jj, t + 1]
                else:
                    v = XPB * Etab[ii, jj - 1, t]
                    if t > 0: v += Etab[ii, jj - 1, t - 1] / (2 * p)
                    v += (t + 1) * Etab[ii, jj - 1, t + 1]
                Etab[ii, jj, t] = v
    return Etab[i, j, : i + j + 1]


def hermite_R(L, p, PC):
    """R[t,u,v] (n = 0) for t+u+v <= L, via the standard downward recursion."""
    X, Y, Zc = PC
    T = p * (X * X + Y * Y + Zc * Zc)
    Rn = np.zeros((L + 1, L + 1, L + 1, L + 1))       # Rn[n, t, u, v]
    for n in range(L + 1):
        Rn[n, 0, 0, 0] = (-2 * p) ** n * boys(n, T)
    for n in range(L - 1, -1, -1):
        for t in range(L - n + 1):
            for u in range(L - n - t + 1):
                for v in range(L - n - t - u + 1):
                    if t + u + v == 0:
                        continue
                    if t > 0:
                        val = X * Rn[n + 1, t - 1, u, v]
                        if t > 1: val += (t - 1) * Rn[n + 1, t - 2, u, v]
                    elif u > 0:
                        val = Y * Rn[n + 1, t, u - 1, v]
                        if u > 1: val += (u - 1) * Rn[n + 1, t, u - 2, v]
                    else:
                        val = Zc * Rn[n + 1, t, u, v - 1]
                        if v > 1: val += (v - 1) * Rn[n + 1, t, u, v - 2]
                    Rn[n, t, u, v] = val
    return Rn[0]


def _overlap_prim(a, s1, A, b, s2, B):
    p = a + b
    out = (np.pi / p) ** 1.5
    for k in range(3):
        out *= hermite_E(s1[k], s2[k], A[k] - B[k], a, b)[0]
    return out


def _kinetic_prim(a, s1, A, b, s2, B):
    l2, m2, n2 = s2
    term0 = b * (2 * (l2 + m2 + n2) + 3) * _overlap_prim(a, s1, A, b, s2, B)
    term1 = 0.0
    for k in range(3):
        up = list(s2); up[k] += 2
        term1 += _overlap_prim(a, s1, A, b, tuple(up), B)
    term1 *= -2 * b * b
    term2 = 0.0
    for k in range(3):
        if s2[k] >= 2:
            dn = list(s2); dn[k] -= 2
            term2 += s2[k] * (s2[k] - 1) * _overlap_prim(a, s1, A, b, tuple(dn), B)
    term2 *= -0.5
    return term0 + term1 + term2


def _pair(a, s1, A, b, s2, B):
    p = a + b
    P = (a * A + b * B) / p
    Ex = hermite_E(s1[0], s2[0], A[0] - B[0], a, b)
    Ey = hermite_E(s1[1], s2[1], A[1] - B[1], a, b)
    Ez = hermite_E(s1[2], s2[2], A[2] - B[2], a, b)
    return p, P, Ex, Ey, Ez


def _nuclear_prim(pair, C):
    p, P, Ex, Ey, Ez = pair
    L = len(Ex) + len(Ey) + len(Ez) - 3
    R = hermite_R(L, p, P - C)
    val = 0.0
    for t, et in enumerate(Ex):
        for u, eu in enumerate(Ey):
            for v, ev in enumerate(Ez):
                val += et * eu * ev * R[t, u, v]
    return 2 * np.pi / p * val


def _eri_prim(pr1, pr2):
    p, P, Ex1, Ey1, Ez1 = pr1
    q, Q, Ex2, Ey2, Ez2 = pr2
    alpha = p * q / (p + q)
    L = len(Ex1) + len(Ey1) + len(Ez1) + len(Ex2) + len(Ey2) + len(Ez2) - 6
    R = hermite_R(L, alpha, P - Q)
    s1 = np.einsum("t,u,v->tuv", Ex1, Ey1, Ez1)
    sgn = lambda n: (-1.0) ** np.arange(n)
    s2 = np.einsum("t,u,v->tuv", Ex2 * sgn(len(Ex2)), Ey2 * sgn(len(Ey2)), Ez2 * sgn(len(Ez2)))
    val = 0.0
    for t, u, v in itertools.product(*map(range, s1.shape)):
        if s1[t, u, v] == 0.0:
            continue
        sub = R[t:t + s2.shape[0], u:u + s2.shape[1], v:v + s2.shape[2]]
        val += s1[t, u, v] * np.sum(s2 * sub)
    return 2 * np.pi ** 2.5 / (p * q * math.sqrt(p + q)) * val


# ---------------------------------------------------------------------------
# Contracted integrals
# ---------------------------------------------------------------------------
def _prims(bf):
    return [(bf.exps[i], bf.norm[i] * bf.coefs[i]) for i in range(len(bf.exps))]


def one_electron(basis, atoms):
    """Return S, T, V_ne (n x n) matrices."""
    n = len(basis)
    S = np.zeros((n, n)); T = np.zeros((n, n)); V = np.zeros((n, n))
    centres = [(Z[s], np.asarray(x, float)) for s, x in atoms]
    for i in range(n):
        for j in range(i, n):
            A, B = basis[i], basis[j]
            s = t = v = 0.0
            for a, ca in _prims(A):
                for b, cb in _prims(B):
                    c = ca * cb
                    s += c * _overlap_prim(a, A.shell, A.origin, b, B.shell, B.origin)
                    t += c * _kinetic_prim(a, A.shell, A.origin, b, B.shell, B.origin)
                    pr = _pair(a, A.shell, A.origin, b, B.shell, B.origin)
                    for Zc, C in centres:
                        v -= c * Zc * _nuclear_prim(pr, C)
            S[i, j] = S[j, i] = s; T[i, j] = T[j, i] = t; V[i, j] = V[j, i] = v
    return S, T, V


def two_electron(basis):
    """(ij|kl) in chemists' notation, using 8-fold symmetry."""
    n = len(basis)
    pairs = {}
    for i in range(n):
        for j in range(n):
            A, B = basis[i], basis[j]
            pairs[i, j] = [(ca * cb, _pair(a, A.shell, A.origin, b, B.shell, B.origin))
                           for a, ca in _prims(A) for b, cb in _prims(B)]
    eri = np.zeros((n, n, n, n))
    for i in range(n):
        for j in range(i + 1):
            ij = i * (i + 1) // 2 + j
            for k in range(n):
                for l in range(k + 1):
                    kl = k * (k + 1) // 2 + l
                    if ij < kl:
                        continue
                    val = 0.0
                    for c1, p1 in pairs[i, j]:
                        for c2, p2 in pairs[k, l]:
                            val += c1 * c2 * _eri_prim(p1, p2)
                    for (a, b, c, d) in [(i, j, k, l), (j, i, k, l), (i, j, l, k), (j, i, l, k),
                                         (k, l, i, j), (l, k, i, j), (k, l, j, i), (l, k, j, i)]:
                        eri[a, b, c, d] = val
    return eri


def nuclear_repulsion(atoms):
    e = 0.0
    for (s1, x1), (s2, x2) in itertools.combinations(atoms, 2):
        e += Z[s1] * Z[s2] / np.linalg.norm(np.subtract(x1, x2))
    return e


# ---------------------------------------------------------------------------
# Restricted Hartree–Fock
# ---------------------------------------------------------------------------
def rhf(S, Hcore, eri, n_elec, e_nuc, max_iter=200, tol=1e-10, verbose=False):
    nocc = n_elec // 2
    w, U = np.linalg.eigh(S)
    X = U @ np.diag(w ** -0.5) @ U.T                    # symmetric orthogonalisation
    eps, C = np.linalg.eigh(X.T @ Hcore @ X)
    C = X @ C
    D = 2 * C[:, :nocc] @ C[:, :nocc].T
    focks, errs = [], []
    E_old = 0.0
    for it in range(max_iter):
        J = np.einsum("pqrs,rs->pq", eri, D)
        K = np.einsum("prqs,rs->pq", eri, D)
        F = Hcore + J - 0.5 * K
        E_el = 0.5 * np.sum(D * (Hcore + F))
        err = F @ D @ S - S @ D @ F
        focks.append(F); errs.append(err)
        focks, errs = focks[-8:], errs[-8:]
        if len(focks) > 1:                               # DIIS extrapolation
            m = len(focks)
            Bm = -np.ones((m + 1, m + 1)); Bm[-1, -1] = 0
            for a in range(m):
                for b in range(m):
                    Bm[a, b] = np.sum(errs[a] * errs[b])
            rhs = np.zeros(m + 1); rhs[-1] = -1
            cfs = np.linalg.lstsq(Bm, rhs, rcond=None)[0][:-1]
            F = sum(c * f for c, f in zip(cfs, focks))
        eps, Cp = np.linalg.eigh(X.T @ F @ X)
        C = X @ Cp
        D_new = 2 * C[:, :nocc] @ C[:, :nocc].T
        if verbose:
            print(it, E_el + e_nuc)
        if abs(E_el - E_old) < tol and np.max(np.abs(D_new - D)) < 1e-8:
            D = D_new
            break
        D, E_old = D_new, E_el
    J = np.einsum("pqrs,rs->pq", eri, D)
    K = np.einsum("prqs,rs->pq", eri, D)
    F = Hcore + J - 0.5 * K
    E_el = 0.5 * np.sum(D * (Hcore + F))
    eps, Cp = np.linalg.eigh(X.T @ F @ X)
    C = X @ Cp
    return E_el + e_nuc, C, eps, D


def molecule_integrals(atoms, scale=1.0):
    basis = sto3g(atoms, scale)
    S, T, V = one_electron(basis, atoms)
    eri = two_electron(basis)
    return dict(basis=basis, S=S, T=T, V=V, eri=eri, e_nuc=nuclear_repulsion(atoms))
