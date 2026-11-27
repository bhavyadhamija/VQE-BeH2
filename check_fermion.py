"""Independent check: build the active-space Hamiltonian from explicit sparse fermion operators (Kronecker products)
and compare with the Pauli-algebra construction used in qubits.py."""
import itertools, numpy as np, scipy.sparse as sp
from beh2_vqe import Point, CORE, ACTIVE, ANGSTROM
from qubits import active_space_integrals, fermion_hamiltonian, to_sparse
p = Point(1.326 * ANGSTROM)
e_core, h, g = active_space_integrals(p.h, p.g, p.e_nuc, CORE, ACTIVE)
n = 2 * h.shape[0]
I2 = sp.identity(2, format="csr"); Zm = sp.csr_matrix(np.diag([1., -1.])); lower = sp.csr_matrix(np.array([[0., 1.], [0., 0.]]))  # |0><1| lowers |1> -> |0>
def kron_list(ops):
    out = sp.identity(1, format="csr")
    for o in ops: out = sp.kron(out, o, format="csr")
    return out
def a(k):   # qubit k is the k-th least significant bit; matrices ordered highest qubit first
    ops = [I2] * n
    for q in range(k): ops[n - 1 - q] = Zm
    ops[n - 1 - k] = lower
    return kron_list(ops)
A = [a(k) for k in range(n)]; Ad = [x.T.tocsr() for x in A]
H = e_core * sp.identity(1 << n, format="csr")
for pq in itertools.product(range(n // 2), repeat=2):
    for s in range(2):
        H = H + h[pq] * (Ad[2 * pq[0] + s] @ A[2 * pq[1] + s])
for p_, q, r, s_ in itertools.product(range(n // 2), repeat=4):
    if abs(g[p_, q, r, s_]) < 1e-12: continue
    for s, t in itertools.product(range(2), repeat=2):
        H = H + 0.5 * g[p_, q, r, s_] * (Ad[2*p_+s] @ Ad[2*r+t] @ A[2*s_+t] @ A[2*q+s])
Hp = to_sparse(fermion_hamiltonian(e_core, h, g)[0], n)
print("max |H_kron - H_pauli| =", abs(H - Hp).max())
