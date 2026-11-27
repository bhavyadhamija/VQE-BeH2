# BeH₂ with the Variational Quantum Eigensolver and the Molecular Virial Theorem

Ground-state potential-energy curve of BeH₂ computed three ways: Hartree–Fock, exact diagonalisation (full CI) and a VQE with a UCCSD ansatz. The project also checks the molecular virial theorem 2⟨T⟩ + ⟨V⟩ = −R·dE/dR for every method.

Everything is written from scratch in NumPy/SciPy, from Gaussian integrals to the qubit Hamiltonian and the VQE. No PySCF, Qiskit or PennyLane is needed. The VQE is a noiseless statevector simulation of the quantum circuit on an ordinary computer.

## Results

| Quantity | Value |
|---|---|
| Qubits (Jordan–Wigner, Be 1s frozen) | 12 |
| Pauli terms in the Hamiltonian | 327 (101 qubit-wise commuting groups) |
| UCCSD parameters | 92 (76 doubles + 16 singles) |
| Equilibrium Be–H distance: HF / exact / VQE | 1.2905 / 1.3165 / 1.3159 Å |
| Energy at equilibrium: HF / exact / VQE | −15.56080 / −15.59525 / −15.59455 Ha |
| Correlation energy | −0.0344 Ha (−21.6 kcal/mol) |
| VQE error vs exact (same 12 qubits), at equilibrium | 0.36 mHa |
| VQE within chemical accuracy (1.6 mHa) | 0.9 – 1.8 Å |
| Largest VQE error | 8.9 mHa at 2.5 Å |
| Virial ratio −⟨V⟩/⟨T⟩ at equilibrium (standard STO-3G) | 1.986 |
| Virial ratio with optimally scaled basis (η = 0.9931) | 1.99943 (residual 1×10⁻⁶ Ha) |

![Energy curve](results/beh2_energy_curve.png)

**VQE accuracy.** Near equilibrium the VQE recovers about 99% of the active-space correlation energy. When the bonds are stretched past about 1.8 Å the ground state becomes strongly multi-configurational. A single-layer Trotterised UCCSD circuit then can't represent it exactly, so the error grows. This is a known limit of the ansatz, consistent with published UCCSD results for BeH₂, and not an optimiser failure.

![Virial theorem](results/beh2_virial.png)

**Virial theorem.** In a basis set with fixed exponents the theorem is not satisfied exactly. The leftover term equals η·∂E/∂η, the change in energy if every basis function were uniformly compressed or expanded:

```
2⟨T⟩ + ⟨V⟩ + R·dE/dR = η·∂E/∂η |η=1
```

The right panel confirms this identity at all 14 geometries for HF and full CI. The VQE follows it only approximately because UCCSD is not fully variational, so the virial mismatch is a second indicator of where the ansatz falls short. Optimising the single scale factor η (Löwdin scaling) restores the theorem to 10⁻⁶ Ha.

![VQE convergence](results/beh2_vqe_convergence.png)

## Method

1. **Geometry:** linear H–Be–H, symmetric stretch, R = 0.9–3.0 Å (14 points).
2. **Basis:** STO-3G (7 basis functions, 6 electrons).
3. **Integrals:** overlap, kinetic, nuclear attraction and electron repulsion via the McMurchie–Davidson scheme.
4. **Hartree–Fock:** restricted HF with DIIS.
5. **Active space:** Be 1s frozen, leaving 4 electrons in 6 orbitals (12 spin orbitals = 12 qubits).
6. **Qubit Hamiltonian:** Jordan–Wigner mapping to Pauli strings.
7. **Exact reference:** diagonalisation in the correct electron-number sector, plus a full 14-qubit FCI without the frozen core.
8. **VQE:** UCCSD from the Hartree–Fock state, L-BFGS-B with exact adjoint gradients, warm-started along the curve.
9. **Virial analysis:** ⟨T⟩ and ⟨V⟩ from density matrices, dE/dR and ∂E/∂η by central finite differences, Richardson extrapolation for the scaled-basis test.

A full explanation of every step, the derivation of the virial identity and a discussion of the numbers is in [docs/WRITEUP.md](docs/WRITEUP.md).

## Validation

| Test | This code | Reference |
|---|---|---|
| H₂, R = 1.4 bohr, RHF | −1.1167142 Ha | −1.1167 (Szabo & Ostlund) |
| H₂, R = 1.4 bohr, FCI | −1.1372759 Ha | −1.1373 (Szabo & Ostlund) |
| H₂O, RHF (Crawford geometry) | −74.94207955 Ha | −74.94207993 (Crawford) |
| BeH₂, R = 1.3 Å, RHF | −15.56127750 Ha | −15.56127803 (MPQC) |
| BeH₂, R = 1.3 Å, force on H | 0.00415797 Ha/bohr | 0.00415795 (MPQC) |

The Pauli-string Hamiltonian also matches an independent build from explicit fermion operators to 4×10⁻¹⁴ (`check_fermion.py`).

## Repository layout

```
.
├── qchem.py            STO-3G basis, McMurchie–Davidson integrals, restricted Hartree–Fock
├── qubits.py           active space, Jordan–Wigner, exact diagonalisation, UCCSD-VQE
├── beh2_vqe.py         the study: energy curve, equilibria, virial analysis, scaled-basis test
├── make_plots.py       figures from the results files
├── validate.py         H2 and H2O reference checks
├── check_fermion.py    independent rebuild of the qubit Hamiltonian
├── results/
│   ├── results_curve.csv            all energies and virial terms at each R
│   ├── results_equilibrium.json     equilibrium data, VQE history, scaled-basis test
│   ├── beh2_pauli_hamiltonian.json  12-qubit Hamiltonian at equilibrium (label -> [re, im])
│   └── *.png                        figures
└── docs/WRITEUP.md     detailed explanation
```

## How to run

```bash
pip install -r requirements.txt
python validate.py          # H2 and H2O checks, a few seconds
python beh2_vqe.py          # full study, about 7 minutes (--quick: 5 geometries, about 3 minutes)
python make_plots.py        # figures into results/
python check_fermion.py     # independent Hamiltonian check
```

Running `beh2_vqe.py` overwrites the files in `results/`.

The saved Hamiltonian can be loaded into Qiskit with:

```python
import json
from qiskit.quantum_info import SparsePauliOp
terms = json.load(open("results/beh2_pauli_hamiltonian.json"))
H = SparsePauliOp.from_list([(k, complex(*v)) for k, v in terms.items()])
```

## Limitations

- STO-3G is a minimal basis, so absolute energies are far from the real molecule. The comparison between methods is what matters.
- Only the symmetric linear stretch is studied.
- The VQE is noiseless. On real hardware each of the 101 measurement groups needs many shots, and gate noise would raise the energy.

## References

- A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry* (Dover, 1996).
- W. J. Hehre, R. F. Stewart and J. A. Pople, *J. Chem. Phys.* 51, 2657 (1969): STO-3G.
- L. E. McMurchie and E. R. Davidson, *J. Comput. Phys.* 26, 218 (1978): integral scheme.
- P.-O. Löwdin, *J. Mol. Spectrosc.* 3, 46 (1959): scaling and the virial theorem.
- A. Peruzzo et al., *Nat. Commun.* 5, 4213 (2014): the variational quantum eigensolver.
- K. Sugisaki et al., *Phys. Chem. Chem. Phys.* (2022): UCCSD and full-CI energies for BeH₂/STO-3G.
- T. D. Crawford, Programming Projects (Project #3): H₂O reference energy.
- MPQC validation output for BeH₂ SCF/STO-3G.

## Author

Bhavya Dhamija, BITS Pilani
