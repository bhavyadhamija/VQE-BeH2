# BeH₂ with VQE and the Virial Theorem: Detailed Write-up

## Results in brief

BeH₂ (linear, symmetric stretch) in the STO-3G basis. All integrals and every simulation are computed from scratch in Python. The VQE runs on an exact statevector simulator, so there is no shot noise or hardware noise.

| Quantity | Value |
| --- | --- |
| Qubits after Jordan–Wigner (frozen Be 1s core) | 12 |
| Pauli terms in the qubit Hamiltonian | 327 (101 qubit-wise commuting measurement groups) |
| UCCSD parameters | 92 (76 doubles + 16 singles) |
| Equilibrium Be–H distance: Hartree–Fock / exact (FCI) / VQE | 1.2905 / 1.3165 / 1.3159 Å |
| Ground-state energy at R_e: HF / FCI / VQE | −15.56080 / −15.59525 / −15.59455 Ha |
| Correlation energy (FCI − HF) | −0.0344 Ha (−21.6 kcal/mol) |
| VQE error vs exact on the same 12 qubits, at R_e | 0.36 mHa (0.22 kcal/mol) |
| Range where VQE is within chemical accuracy (1.6 mHa) | 0.9 – 1.8 Å |
| Largest VQE error on the curve | 8.9 mHa at 2.5 Å |
| Virial ratio −⟨V⟩/⟨T⟩ near R_e, standard STO-3G | 1.986 |
| Virial ratio with the optimally scaled basis (η = 0.9931) | 1.99943; virial residual 1×10⁻⁶ Ha |

Three things to take away:

1. **VQE works well near equilibrium.** At R_e it recovers 99% of the correlation energy available in the 12-qubit active space (98% of the full-CI correlation energy, counting the frozen core). It gives the bond length to within 0.0006 Å of the exact answer.
2. **It degrades when bonds are stretched.** Beyond about 1.8 Å the ground state becomes strongly multi-configurational, and single-layer Trotterised UCCSD cannot represent it exactly. That is a property of the ansatz, not a bug in the optimiser. The gradient was checked against finite differences, and reordering the operators or adding a second layer did not close the gap.
3. **The virial theorem is a sharp diagnostic.** With a fixed basis set it is *not* satisfied exactly. The leftover term equals exactly how much the energy would change if the basis were uniformly compressed or stretched. Optimise that one scale factor and the theorem holds to 10⁻⁶ Ha.

## The method, step by step

Every step from the molecule to the VQE energy is implemented from scratch in NumPy/SciPy. Units are atomic (bohr, hartree); 1 Å = 1.8897 bohr.

1. **Geometry.** Linear H–Be–H with Be at the origin and both Be–H distances equal to R (the symmetric stretch). R runs from 0.9 to 3.0 Å at 14 points.
2. **Basis set.** STO-3G: Be gets 1s, 2s, 2pₓ, 2pᵧ and 2p_z; each H gets 1s. That makes 7 basis functions, each a fixed combination of 3 Gaussians.
3. **Integrals.** Overlap S, kinetic T, nuclear attraction V and electron repulsion (pq|rs), computed with the McMurchie–Davidson method (Hermite Gaussians plus the Boys function). Kinetic and potential matrices are kept separate because the virial theorem needs ⟨T⟩ and ⟨V⟩ separately.
4. **Hartree–Fock (RHF).** A self-consistent field loop with DIIS acceleration gives 7 molecular orbitals; 3 are doubly occupied by the 6 electrons. Orbital energies at equilibrium (Ha):
   - Be 1s: −4.519
   - occupied σg and σu: −0.460 and −0.424
   - empty degenerate Be 2pπ pair: +0.211, +0.211
   - empty σ*: +0.468 and +0.963
5. **Active space.** The Be 1s core orbital (−4.5 Ha) is frozen, leaving 4 electrons in 6 orbitals. That is 12 spin orbitals and therefore 12 qubits. Freezing the core changes the energy by only 3.4 × 10⁻⁴ Ha.
6. **Qubit Hamiltonian.** The second-quantised Hamiltonian is mapped to qubits with the Jordan–Wigner transformation: spin orbital k becomes qubit k, with α and β interleaved. The result has 327 Pauli strings, which can be measured in 101 qubit-wise commuting groups.
7. **Exact reference (FCI).** The qubit Hamiltonian is diagonalised exactly in the 4-electron, S_z = 0 sector (225 states). A full 14-qubit FCI without the frozen core gives the true STO-3G answer.
8. **VQE ansatz (UCCSD).** The circuit starts from the Hartree–Fock bit string |000000001111⟩. It applies one exponential exp(θₖ(τₖ − τₖ†)) for each spin-conserving excitation: 76 doubles, then 16 singles, giving 92 parameters.
9. **Optimisation.** The energy is minimised with L-BFGS-B. Exact gradients come from one forward and one backward sweep through the circuit (the adjoint method), checked against finite differences to 7 digits. Each geometry is warm-started from the previous geometry's angles.
10. **Virial analysis.** At every R:
    - ⟨T⟩ and ⟨V⟩ are computed for HF, FCI and VQE.
    - dE/dR comes from central differences (±0.01 bohr).
    - ∂E/∂η comes from rescaling every basis exponent by η².
    - At equilibrium, η is optimised and the theorem is re-checked.

The VQE is a noiseless statevector simulation of the circuit on an ordinary computer; no quantum hardware or shot noise is involved.

## Validation

The chemistry engine reproduces every published value it was tested against to within 6 × 10⁻⁷ Ha. The remaining digits come from rounding in the basis-set coefficients.

| Test | This code | Published | Difference |
| --- | --- | --- | --- |
| H₂, R = 1.4 bohr, RHF energy (Ha) | −1.1167142 | −1.1167 (Szabo & Ostlund) | < 10⁻⁷ |
| H₂, R = 1.4 bohr, FCI energy (Ha) | −1.1372759 | −1.1373 (Szabo & Ostlund) | < 10⁻⁷ |
| H₂O, Crawford tutorial geometry, RHF (Ha) | −74.94207955 | −74.94207993 (Crawford) | 3.8 × 10⁻⁷ |
| BeH₂, R = 1.3 Å, RHF energy (Ha) | −15.56127750 | −15.56127803 ([MPQC reference output](https://molecuilder.de/browser/ThirdParty/mpqc_open/src/bin/mpqc/validate/ref/basis1_beh2scfsto3gd2h.out)) | 5.4 × 10⁻⁷ |
| BeH₂, R = 1.3 Å, force on each H (Ha/bohr) | 0.00415797 | 0.00415795 (MPQC) | 2 × 10⁻⁸ |

**Internal checks:**

- The qubit Hamiltonian's energy on the Hartree–Fock bit string equals the RHF energy exactly. This checks the whole Jordan–Wigner mapping.
- An independent construction of H from explicit fermion operators (sparse Kronecker products, `check_fermion.py`) agrees with the Pauli-string version to 4 × 10⁻¹⁴.
- The adjoint VQE gradient agrees with a finite-difference gradient to 7 digits (0.0224151 vs 0.0224151).

**The UCCSD accuracy is consistent with the literature.** A Trotterised UCCSD circuit is not exactly the same as the textbook UCCSD exponential, and its energy depends slightly on the order of the operators. For BeH₂ in STO-3G, UCCSD energies are reported a few × 10⁻⁴ Ha above full CI near the equilibrium region, e.g. −15.77875 vs −15.77917 Ha ([Sugisaki et al., PCCP 2022, supplementary data](https://www.rsc.org/suppdata/d1/cp/d1cp04318h/d1cp04318h1.pdf)). This work finds 3.6 × 10⁻⁴ Ha at equilibrium. Reversing the operator order (singles before doubles) and adding a second circuit layer did not lower it. That confirms this is the ansatz's limit and not an optimiser failure.

## Energy curve

![Energy curve](../results/beh2_energy_curve.png)

The three curves have the same shape:

- **Hartree–Fock** lies above the exact curve everywhere. The gap (the correlation energy) is 0.034 Ha at equilibrium and grows to 0.31 Ha at 3.0 Å, because one determinant cannot describe two breaking bonds.
- **The exact energy** has its minimum at R = 1.317 Å (E = −15.595246 Ha). Hartree–Fock puts the minimum at 1.291 Å and VQE at 1.316 Å. For comparison, VQE studies of BeH₂ commonly use a near-equilibrium distance of 1.3264 Å ([e.g. arXiv:2310.09418](https://arxiv.org/pdf/2310.09418)).

The VQE error (lower panel) is measured against the exact energy of the *same* 12-qubit Hamiltonian, so it isolates the circuit and the optimiser. Near equilibrium the error is 0.35 mHa, a quarter of chemical accuracy. It grows as the bonds stretch, because the circuit starts from Hartree–Fock and UCCSD includes only single and double excitations. Freezing the Be 1s core adds a separate, nearly constant 0.33–0.35 mHa at every distance.

![VQE convergence](../results/beh2_vqe_convergence.png)

## Molecular virial theorem

### What the theorem says

For the exact electronic wavefunction of a molecule with the nuclei clamped at separation R, the kinetic energy ⟨T⟩ and the total potential energy ⟨V⟩ satisfy

**2⟨T⟩ + ⟨V⟩ = −R · dE/dR**

Here ⟨V⟩ covers electron–nuclear attraction, electron–electron repulsion and nuclear–nuclear repulsion.

This is the molecular generalisation of the particle-in-a-box and atomic results. It follows because every Coulomb term scales as 1/length and the kinetic energy as 1/length². At the equilibrium geometry dE/dR = 0, so it reduces to the textbook form 2⟨T⟩ + ⟨V⟩ = 0, i.e. **−⟨V⟩/⟨T⟩ = 2 and E = −⟨T⟩**. Away from equilibrium the right-hand side is not zero, and its sign says whether the bond wants to shrink or stretch.

The code computes each piece directly:

- ⟨T⟩ from the one-particle density matrix contracted with the kinetic-energy integrals.
- ⟨V⟩ = E − ⟨T⟩.
- dE/dR by central finite differences (step 0.01 bohr), each displaced point re-solved from scratch.

This is done for the Hartree–Fock, exact (FCI) and VQE wavefunctions at all 14 bond lengths.

### Why it is not satisfied in a finite basis, and the exact correction

The theorem assumes the wavefunction is free to rescale itself in space. A Gaussian basis set with fixed exponents cannot do that, so the identity picks up an extra term. Scaling every exponent α → η²α is equivalent to stretching all electron coordinates by 1/η. Following that through gives an identity that holds *exactly* for any variational method in a fixed basis:

**2⟨T⟩ + ⟨V⟩ + R · dE/dR = η · ∂E/∂η |η=1**

The right-hand side measures how much the energy would drop if the basis were uniformly compressed or expanded. The code evaluates it independently, by re-running the whole calculation with η = 1 ± 0.002.

![Virial theorem](../results/beh2_virial.png)

The right panel confirms the identity at every geometry: residual and scaling derivative fall on the line y = x. They agree to 1.2×10⁻⁴ Ha for HF and FCI, which is the size of the finite-difference error.

### Selected numbers (exact FCI, standard STO-3G)

| R (Å) | −⟨V⟩/⟨T⟩ | 2⟨T⟩+⟨V⟩ | −R·dE/dR | Residual | η∂E/∂η |
| --- | --- | --- | --- | --- | --- |
| 1.0 | 1.9502 | 0.8110 | 0.8954 | −0.0844 | −0.0843 |
| 1.3 | 1.9860 | 0.2221 | 0.0319 | 0.1902 | 0.1902 |
| 1.6 | 1.9934 | 0.1033 | −0.3815 | 0.4848 | 0.4848 |
| 2.2 | 1.9704 | 0.4694 | −0.4680 | 0.9374 | 0.9374 |
| 3.0 | 1.9394 | 0.9885 | −0.0086 | 0.9971 | 0.9971 |

All energies are in hartree. Residual = (2⟨T⟩+⟨V⟩) − (−R·dE/dR).

How to read the table:

- **Near equilibrium the ratio is 1.986, not 2.** The scaling derivative is positive, meaning the basis is slightly too compact: the energy would still go down if every orbital were expanded by about 0.7%.
- **When stretched the residual grows to about 1 Ha.** At 3.0 Å the molecule is nearly Be + 2H. STO-3G hydrogen uses ζ = 1.24, tuned for H inside molecules; a free H atom wants ζ = 1.0. The basis is far too compact there, so the scaling derivative is large.
- **HF and FCI follow the same identity.** Both are variational in the same space, so both satisfy it.
- **VQE follows it only approximately.** The mismatch is 0.0007 Ha near equilibrium and grows to 0.015 Ha at 2.5 Å. The identity requires the energy to be fully stationary with respect to every allowed change of the wavefunction. UCCSD is stationary only within its own restricted parameter space. So the virial mismatch is a second, independent indicator of where the ansatz is falling short, and it lines up with the energy error curve.

### Restoring the theorem: optimally scaled basis

At the FCI equilibrium distance the code minimises the FCI energy over the single scale factor η (the Löwdin scaling procedure). The optimum is **η = 0.99307**, a 0.7% expansion of every orbital. At that basis:

- −⟨V⟩/⟨T⟩ = **1.99943**.
- 2⟨T⟩+⟨V⟩ = 0.0089296 Ha and −R·dE/dR = 0.0089286 Ha, so the virial residual is **9.9×10⁻⁷ Ha**.

The remaining 0.009 Ha on each side is real, not an error. Rescaling the basis shifts the minimum slightly, so this geometry is no longer exactly at the bottom of the curve and there is a small force. The theorem holds in its full, force-corrected form.

One numerical detail matters here. With a plain central difference the residual came out as 3.4×10⁻⁵ Ha. The error shrank fourfold each time the step was halved, which is the signature of an O(d²) truncation error rather than physics. Richardson extrapolation of two step sizes, (4·g(d) − g(2d))/3, removes the leading error and gives the 10⁻⁶ value above.

**Link to the particle-in-a-box:** this is the same physics in a harder setting. For the box, ⟨V⟩ = 0 and the virial content is all in the wall force: 2⟨T⟩ = L · (−dE/dL). Here the R·dE/dR term plays the role of the wall: it is the work done when the nuclear framework is stretched.

## The code

Six Python files, about 830 lines in total. Everything needs only NumPy, SciPy and Matplotlib, and the full run takes 6–8 minutes on a laptop.

| File | What it does |
| --- | --- |
| `qchem.py` | STO-3G basis; overlap, kinetic, nuclear-attraction and electron-repulsion integrals (McMurchie–Davidson); restricted Hartree–Fock with DIIS |
| `qubits.py` | Frozen-core active space; Jordan–Wigner mapping to Pauli strings; sparse qubit matrices; exact ground state; UCCSD generators; VQE with adjoint gradients |
| `beh2_vqe.py` | The experiment: energy curve, equilibrium distances, virial analysis, scaled-basis test; writes the results files |
| `make_plots.py` | Three figures from the results files |
| `validate.py` | Reproduces the H₂ and H₂O textbook energies |
| `check_fermion.py` | Rebuilds the qubit Hamiltonian independently and compares |

### qchem.py, function by function

- **`sto3g(atoms, scale)`** builds the basis functions. Each atom's exponents are the universal STO-3G values times ζ², where ζ is the atom's Slater exponent (H 1.24; Be 1s 3.68, 2sp 1.15). The `scale` argument multiplies every exponent by η² for the virial test.
- **`hermite_E`** and **`hermite_R`** are the two McMurchie–Davidson building blocks. The first expands a product of two Gaussians in Hermite Gaussians. The second gives the Coulomb integrals of Hermite Gaussians through the Boys function.
- **`one_electron`** returns S, T and V_ne.
- **`two_electron`** returns (pq|rs), computing only the 1 of every 8 symmetry-related integrals that is unique.
- **`rhf`** does symmetric orthogonalisation, then an SCF loop with DIIS (least-squares form, so it never fails once converged). It returns the energy, the orbitals, the orbital energies and the density matrix.

### qubits.py, function by function

- **Pauli algebra.**
  - A Pauli string is a pair of bit masks (x, z) meaning ∏ Xₖ^xₖ Zₖ^zₖ.
  - `op_mul` multiplies two operators, tracking the sign from moving Z past X.
  - `creation(j)` = Z₀…Zⱼ₋₁ · Xⱼ(1+Zⱼ)/2 is the Jordan–Wigner aⱼ†; `annihilation(j)` is aⱼ.
- **`active_space_integrals`** folds the frozen core into a constant energy and an effective one-electron matrix.
- **`fermion_hamiltonian`** builds H = E_core + Σ h_pq a†_p a_q + ½ Σ (pq|rs) a†_p a†_r a_s a_q over spins, as Pauli strings.
- **`to_labels`** gives Qiskit-style labels such as "IIZZXXYY…" (rightmost character = qubit 0).
- **`to_sparse`** gives the 4096 × 4096 matrix; each Pauli string is a signed permutation, built with bit operations.
- **`sector`** and **`exact_ground_state`** diagonalise in the right electron-number sector.
- **`uccsd_generators`** lists every spin-conserving single and double excitation τ and returns the anti-Hermitian matrix τ − τ† for each.
- **`vqe_uccsd`** applies exp(θₖGₖ) one generator at a time. It gets the energy and all 92 gradient components from one forward and one backward sweep, then runs L-BFGS-B.

### beh2_vqe.py, the experiment

- **`Point(R, eta)`** computes integrals, RHF, the MO-basis integrals, and the separate kinetic matrix at one geometry and basis scale.
- **`.full_fci()`** runs the 14-qubit exact calculation.
- **`.active()`** builds the 12-qubit Hamiltonian and its exact energy.
- **`.vqe()`** runs VQE.
- **`main()`** loops over 14 distances. At each one it also solves at R ± 0.01 bohr (for dE/dR) and at η = 1 ± 0.002 (for ∂E/∂η). It then finds each method's equilibrium distance, runs the scaled-basis virial test, and writes the files to `results/`.

## Assumptions and limitations

The numbers are exact for the model chosen. The model itself is minimal, and the VQE is simulated.

- **Minimal basis.** STO-3G has one function per occupied atomic shell (plus Be 2p), so absolute energies are far from the true molecule. For comparison, the Hartree–Fock limit for BeH₂ is about −15.77 Ha, against −15.56 Ha here. The trends and method comparisons are what matter.
- **Symmetric stretch only.** Both Be–H bonds change together and the molecule stays linear; bending and asymmetric stretching are not studied.
- **Frozen core.** Be 1s is frozen for the qubit calculations, at a cost of 0.33–0.35 mHa, which is measured, not assumed.
- **Noiseless simulation.** Energies are exact expectation values. On hardware, each of the 101 measurement groups needs many shots, and gate noise would raise the energy.
- **Ansatz.** A single Trotter step of UCCSD from Hartree–Fock. It is accurate near equilibrium but degrades when bonds break (up to 8.9 mHa at 2.5 Å), where a multireference starting point or an adaptive ansatz would be needed.
- **Finite differences.** dE/dR and ∂E/∂η come from central differences. Their error is about 10⁻⁴ Ha for the fixed-basis checks; the scaled-basis test removes it by Richardson extrapolation.

## Sources

- K. Sugisaki, T. Kato, Y. Minato, K. Okuwaki, Y. Mochizuki, "Variational quantum eigensolver simulations with the multireference unitary coupled cluster ansatz: a case study of the C2v quasi-reaction pathway of beryllium insertion into a H₂ molecule," *Phys. Chem. Chem. Phys.* (2022), [supplementary information](https://www.rsc.org/suppdata/d1/cp/d1cp04318h/d1cp04318h1.pdf). UCCSD vs full-CI energies for BeH₂/STO-3G and Trotter-ordering effects.
- [MPQC validation output, BeH₂ SCF/STO-3G at 1.3 Å](https://molecuilder.de/browser/ThirdParty/mpqc_open/src/bin/mpqc/validate/ref/basis1_beh2scfsto3gd2h.out). Reference energy and gradient.
- W. J. Hehre, R. F. Stewart, J. A. Pople, *J. Chem. Phys.* 51, 2657 (1969). STO-3G contraction parameters.
- [arXiv:2310.09418](https://arxiv.org/pdf/2310.09418). BeH₂ near-equilibrium distance of 1.3264 Å used in VQE studies.
- [NASA NTRS 19670058410](https://ntrs.nasa.gov/citations/19670058410). Hartree–Fock limit estimate for BeH₂ (about −15.773 Ha).
- A. Szabo and N. S. Ostlund, *Modern Quantum Chemistry* (Dover, 1996). H₂ STO-3G reference energies.
- T. D. Crawford, Programming Projects (Project #3, Hartree–Fock). H₂O STO-3G reference energy.
- P.-O. Löwdin, "Scaling problem, virial theorem, and connected relations in quantum mechanics," *J. Mol. Spectrosc.* 3, 46 (1959). The scaling form of the virial theorem.
