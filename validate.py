import time, numpy as np
from qchem import molecule_integrals, rhf

def run(atoms, nel, name):
    t0 = time.time()
    m = molecule_integrals(atoms)
    E, C, eps, D = rhf(m["S"], m["T"] + m["V"], m["eri"], nel, m["e_nuc"])
    print(f"{name}: E_RHF = {E:.10f} Ha   ({len(m['basis'])} basis functions, {time.time()-t0:.1f} s)")
    return m, E, C

# H2, R = 1.4 bohr (Szabo & Ostlund): E_RHF = -1.1167, E_FCI = -1.1373
m, E, C = run([("H", (0, 0, 0)), ("H", (0, 0, 1.4))], 2, "H2")
# 2x2 FCI for H2 in the MO basis
h = C.T @ (m["T"] + m["V"]) @ C
g = np.einsum("pi,qj,rk,sl,pqrs->ijkl", C, C, C, C, m["eri"])
E0 = 2*h[0,0] + g[0,0,0,0]; E1 = 2*h[1,1] + g[1,1,1,1]; K = g[0,1,0,1]
Hci = np.array([[E0, K], [K, E1]])
print(f"H2: E_FCI = {np.linalg.eigvalsh(Hci)[0] + m['e_nuc']:.10f} Ha")

# H2O, Crawford programming project geometry (bohr): E_SCF = -74.942079928192
geo = [("O", (0.0, -0.143225816552, 0.0)), ("H", (1.638036840407, 1.136548822547, 0.0)),
       ("H", (-1.638036840407, 1.136548822547, 0.0))]
run(geo, 10, "H2O")
