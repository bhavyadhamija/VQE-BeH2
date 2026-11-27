"""Plots from results/results_curve.csv and results/results_equilibrium.json (run beh2_vqe.py first)."""
import csv
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAVY, PURPLE, BLUE, GREY, RED = "#1E2A4A", "#5B4BB7", "#2F6FDE", "#6B7280", "#B03A2E"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
rows = list(csv.DictReader(open(os.path.join(OUT, "results_curve.csv"))))
col = lambda k: np.array([float(r[k]) for r in rows])
eq = json.load(open(os.path.join(OUT, "results_equilibrium.json")))
R = col("R_A")

# 1. energy curve + VQE error
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 7.2), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
a1.plot(R, col("E_HF"), "o-", color=GREY, label="Hartree–Fock")
a1.plot(R, col("E_FCI"), "-", color=NAVY, lw=2.2, label="Exact (full CI, 14 qubits)")
a1.plot(R, col("E_VQE"), "o", color=PURPLE, mfc="none", ms=8, mew=1.6, label="VQE-UCCSD (12 qubits, frozen core)")
a1.axvline(eq["R_eq_A"]["E_fci"], color=GREY, ls=":", lw=1)
a1.text(eq["R_eq_A"]["E_fci"] + 0.03, -15.25, f"R$_e$ = {eq['R_eq_A']['E_fci']:.3f} Å", color=GREY)
a1.set_ylabel("Energy (Hartree)")
a1.set_title("BeH$_2$ symmetric stretch, STO-3G")
a1.legend(frameon=False)
a2.semilogy(R, col("VQE_minus_CAS"), "o-", color=PURPLE, label="VQE − exact (same 12-qubit Hamiltonian)")
a2.semilogy(R, col("E_CAS") - col("E_FCI"), "s--", color=BLUE, ms=4, label="frozen-core error")
a2.axhline(1.6e-3, color=RED, ls="--", lw=1)
a2.text(R[0], 1.9e-3, "chemical accuracy (1.6 mHa)", color=RED, fontsize=9)
a2.set_xlabel("Be–H distance R (Å)")
a2.set_ylabel("Error (Hartree)")
a2.legend(frameon=False, fontsize=9, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "beh2_energy_curve.png"), dpi=200)

# 2. virial theorem
fig, (b1, b2) = plt.subplots(1, 2, figsize=(11, 4.4))
for name, c, mk in [("HF", GREY, "o"), ("FCI", NAVY, "s"), ("VQE", PURPLE, "^")]:
    b1.plot(R, col(f"{name}_ratio"), mk + "-", color=c, label=name, ms=5)
    b2.plot(col(f"{name}_eta_dEdeta"), col(f"{name}_residual"), mk, color=c, label=name, ms=6, mfc="none" if name == "VQE" else c)
sv = eq["scaled_basis_virial_FCI"]
b1.plot([eq["at_R_eq_FCI"]["R_A"]], [sv["ratio"]], "*", color=RED, ms=14, label=f"FCI, optimally scaled basis (η = {sv['eta']:.4f})")
b1.axhline(2, color=GREY, ls="--", lw=1)
b1.set_xlabel("Be–H distance R (Å)"); b1.set_ylabel("−⟨V⟩/⟨T⟩")
b1.set_title("Virial ratio (exact theory: 2 at equilibrium)")
b1.legend(frameon=False, fontsize=9)
lim = [min(col("HF_eta_dEdeta").min(), col("FCI_residual").min()) - 0.1, max(col("HF_eta_dEdeta").max(), col("HF_residual").max()) + 0.1]
b2.plot(lim, lim, color=GREY, lw=1, ls="--", label="y = x")
b2.set_xlabel("η ∂E/∂η at η = 1   (Hartree)")
b2.set_ylabel("2⟨T⟩ + ⟨V⟩ + R dE/dR   (Hartree)")
b2.set_title("Virial residual = basis-scaling derivative")
b2.legend(frameon=False, fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "beh2_virial.png"), dpi=200)

# 3. VQE convergence at equilibrium
h = np.array(eq["at_R_eq_FCI"]["vqe_energy_history"])
Ecas = eq["at_R_eq_FCI"]["E_FCI_active"]
fig, c1 = plt.subplots(figsize=(6.5, 3.8))
c1.semilogy(np.arange(1, len(h) + 1), h - Ecas, "o-", color=PURPLE, ms=4)
c1.axhline(1.6e-3, color=RED, ls="--", lw=1); c1.text(1, 1.9e-3, "chemical accuracy", color=RED, fontsize=9)
c1.axhline(eq["at_R_eq_FCI"]["E_HF"] - Ecas, color=GREY, ls=":", lw=1)
c1.text(len(h) * 0.55, (eq["at_R_eq_FCI"]["E_HF"] - Ecas) * 1.15, "Hartree–Fock start", color=GREY, fontsize=9)
c1.set_xlabel("energy + gradient evaluation"); c1.set_ylabel("E(θ) − exact  (Hartree)")
c1.set_title(f"VQE-UCCSD at R = {eq['at_R_eq_FCI']['R_A']:.3f} Å ({eq['at_R_eq_FCI']['uccsd_parameters']} parameters)")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "beh2_vqe_convergence.png"), dpi=200)
print("plots written")
