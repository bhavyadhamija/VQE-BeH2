"""
beh2_vqe.py — BeH2 ground-state energy curve with Hartree–Fock, exact diagonalisation (FCI)
and VQE-UCCSD, plus a molecular virial-theorem analysis.

Geometry: linear H–Be–H, Be at the origin, both Be–H distances equal to R (symmetric stretch).
Basis: STO-3G (7 spatial orbitals, 6 electrons). Active space: Be 1s frozen -> 6 orbitals,
4 electrons -> 12 qubits (Jordan–Wigner).

The VQE is a classical statevector simulation of the UCCSD circuit (no quantum hardware).

Usage:
    python beh2_vqe.py                 # full run: curve + virial + scaling test (a few minutes)
    python beh2_vqe.py --quick         # 5 geometries only
Outputs (in results/): results_curve.csv, results_equilibrium.json, beh2_pauli_hamiltonian.json
"""
import argparse
import json
import os
import time

import numpy as np
from scipy.optimize import minimize_scalar

from qchem import molecule_integrals, rhf, ANGSTROM
from qubits import (active_space_integrals, fermion_hamiltonian, one_body_operator, to_sparse, to_labels,
                    sector, exact_ground_state, qubitwise_groups, uccsd_generators, uccsd_state, vqe_uccsd)

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
CORE, ACTIVE = [0], list(range(1, 7))
N_ELEC = 6


def geometry(R_bohr):
    return [("Be", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, R_bohr)), ("H", (0.0, 0.0, -R_bohr))]


class Point:
    """Everything at one geometry R (bohr) and basis scale eta (exponents x eta^2)."""

    def __init__(self, R, eta=1.0):
        self.R, self.eta = R, eta
        m = molecule_integrals(geometry(R), scale=eta)
        self.E_hf, C, self.eps, D = rhf(m["S"], m["T"] + m["V"], m["eri"], N_ELEC, m["e_nuc"])
        self.T_hf = float(np.sum(D * m["T"]))
        h = C.T @ (m["T"] + m["V"]) @ C
        t = C.T @ m["T"] @ C
        g = np.einsum("pi,qj,rk,sl,pqrs->ijkl", C, C, C, C, m["eri"], optimize=True)
        self.h, self.t, self.g, self.e_nuc = h, t, g, m["e_nuc"]

    def full_fci(self):
        H, nq = fermion_hamiltonian(self.e_nuc, self.h, self.g)
        Hs = to_sparse(H, nq)
        E, psi = exact_ground_state(Hs, sector(nq, 3, 3))
        Ts = to_sparse(one_body_operator(self.t), nq)
        self.E_fci, self.T_fci = E, float(psi @ Ts @ psi)
        return E

    def active(self):
        e_core, h_a, g_a = active_space_integrals(self.h, self.g, self.e_nuc, CORE, ACTIVE)
        H, nq = fermion_hamiltonian(e_core, h_a, g_a)
        self.H_pauli, self.nq = H, nq
        self.Hs = to_sparse(H, nq)
        t_core = 2 * sum(self.t[c, c] for c in CORE)
        self.Ts = to_sparse(one_body_operator(self.t[np.ix_(ACTIVE, ACTIVE)], e0=t_core), nq)
        self.E_cas, psi = exact_ground_state(self.Hs, sector(nq, 2, 2))
        self.T_cas = float(psi @ self.Ts @ psi)
        self.ref = np.zeros(1 << nq); self.ref[0b1111] = 1.0
        return self.E_cas

    def vqe(self, gens, theta0=None):
        res, hist = vqe_uccsd(self.Hs, gens, self.ref, theta0=theta0, tol=1e-12, maxiter=1000)
        psi = uccsd_state(res.x, gens, self.ref)
        self.E_vqe, self.T_vqe, self.theta, self.hist = res.fun, float(psi @ self.Ts @ psi), res.x, hist
        return res.fun


def solve(R, eta=1.0, gens=None, theta0=None, vqe=True):
    p = Point(R, eta)
    p.full_fci(); p.active()
    if vqe:
        p.vqe(gens, theta0)
    return p


def virial_terms(E, T, R, dEdR):
    V = E - T
    return dict(T=T, V=V, ratio=-V / T, lhs=2 * T + V, rhs=-R * dEdR, residual=2 * T + V + R * dEdR)


def main():
    os.makedirs(OUT, exist_ok=True)
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    t_start = time.time()
    gens = uccsd_generators(6, 2)
    d = 0.01                                    # bohr, for finite differences in R
    de = 0.002                                  # for the derivative with respect to basis scale eta

    Rs_A = [1.0, 1.3, 1.6, 2.2, 3.0] if args.quick else [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.5, 2.8, 3.0]
    rows, theta = [], None
    print(f"{'R/A':>5} {'E_HF':>13} {'E_FCI':>13} {'E_CAS':>13} {'E_VQE':>13} {'VQE-CAS':>9}")
    for RA in Rs_A:
        R = RA * ANGSTROM
        p = solve(R, gens=gens, theta0=theta)
        theta = p.theta
        pm = solve(R - d, gens=gens, theta0=theta)
        pp = solve(R + d, gens=gens, theta0=theta)
        e_lo = solve(R, eta=1 - de, gens=gens, theta0=theta)
        e_hi = solve(R, eta=1 + de, gens=gens, theta0=theta)
        row = dict(R_A=RA, R_bohr=R, E_HF=p.E_hf, E_FCI=p.E_fci, E_CAS=p.E_cas, E_VQE=p.E_vqe,
                   VQE_minus_CAS=p.E_vqe - p.E_cas, n_iter=len(p.hist))
        for name, E, T, Em, Ep, El, Eh in [
                ("HF", p.E_hf, p.T_hf, pm.E_hf, pp.E_hf, e_lo.E_hf, e_hi.E_hf),
                ("FCI", p.E_fci, p.T_fci, pm.E_fci, pp.E_fci, e_lo.E_fci, e_hi.E_fci),
                ("VQE", p.E_vqe, p.T_vqe, pm.E_vqe, pp.E_vqe, e_lo.E_vqe, e_hi.E_vqe)]:
            vt = virial_terms(E, T, R, (Ep - Em) / (2 * d))
            vt["eta_dEdeta"] = (Eh - El) / (2 * de)          # predicted residual: eta * dE/deta at eta = 1
            for k, v in vt.items():
                row[f"{name}_{k}"] = v
        rows.append(row)
        print(f"{RA:5.2f} {p.E_hf:13.8f} {p.E_fci:13.8f} {p.E_cas:13.8f} {p.E_vqe:13.8f} {p.E_vqe - p.E_cas:9.2e}")

    # equilibrium distances (minimise each method's curve)
    def emin(method):
        f = lambda RA: getattr(solve(RA * ANGSTROM, vqe=(method == "E_vqe"), gens=gens, theta0=theta), method)
        r = minimize_scalar(f, bounds=(1.20, 1.45), method="bounded", options={"xatol": 1e-4})
        return r.x, r.fun
    Re = {m: emin(m) for m in ["E_hf", "E_fci", "E_vqe"]}
    print("equilibrium R (A) and energy:", {k: (round(v[0], 4), round(v[1], 8)) for k, v in Re.items()})

    # detailed record at the FCI equilibrium
    RA_e = Re["E_fci"][0]; Rb_e = RA_e * ANGSTROM
    pe = solve(Rb_e, gens=gens, theta0=theta)
    labels = to_labels(pe.H_pauli, pe.nq)
    groups = qubitwise_groups(list(labels))
    json.dump({k: [float(np.real(v)), float(np.imag(v))] for k, v in sorted(labels.items())},
              open(os.path.join(OUT, "beh2_pauli_hamiltonian.json"), "w"), indent=0)

    # basis-scaling test at the FCI equilibrium: optimise eta, then check the virial theorem
    def best_eta(Rb):
        r = minimize_scalar(lambda e: Point(Rb, e).full_fci(), bounds=(0.9, 1.1), method="bounded",
                            options={"xatol": 1e-10})
        return r.x, r.fun
    eta_s, E_s = best_eta(Rb_e)
    ps = Point(Rb_e, eta_s); ps.full_fci()
    # dG/dR with Richardson extrapolation (steps d and 2d) to remove the O(d^2) finite-difference error
    g1 = (best_eta(Rb_e + d)[1] - best_eta(Rb_e - d)[1]) / (2 * d)
    g2 = (best_eta(Rb_e + 2 * d)[1] - best_eta(Rb_e - 2 * d)[1]) / (4 * d)
    scaled = virial_terms(ps.E_fci, ps.T_fci, Rb_e, (4 * g1 - g2) / 3)
    scaled["eta"] = eta_s
    scaled["residual_without_extrapolation"] = 2 * ps.T_fci + (ps.E_fci - ps.T_fci) + Rb_e * g1

    eq = dict(R_eq_A={k: v[0] for k, v in Re.items()}, E_at_R_eq={k: v[1] for k, v in Re.items()},
              at_R_eq_FCI=dict(R_A=RA_e, E_HF=pe.E_hf, E_FCI_full=pe.E_fci, E_FCI_active=pe.E_cas,
                               E_VQE=pe.E_vqe, frozen_core_error=pe.E_cas - pe.E_fci,
                               VQE_minus_active_FCI=pe.E_vqe - pe.E_cas, correlation_energy=pe.E_fci - pe.E_hf,
                               qubits=pe.nq, pauli_terms=len(labels), qwc_groups=len(groups),
                               uccsd_parameters=len(gens), vqe_energy_history=pe.hist,
                               orbital_energies=pe.eps.tolist()),
              scaled_basis_virial_FCI=scaled, runtime_s=time.time() - t_start)
    json.dump(eq, open(os.path.join(OUT, "results_equilibrium.json"), "w"), indent=1)

    import csv
    with open(os.path.join(OUT, "results_curve.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(json.dumps({k: v for k, v in eq["at_R_eq_FCI"].items() if k != "vqe_energy_history"}, indent=1))
    print("scaled-basis virial (FCI):", json.dumps(scaled, indent=1))
    print(f"done in {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    main()
