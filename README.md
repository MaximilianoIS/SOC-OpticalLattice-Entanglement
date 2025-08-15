# Spin Dynamics Simulation: OAT & Exact Diagonalization

This repository contains C++ implementations of spin dynamics simulations for models described in the main paper and related works.  
So far, this includes:

- One-Axis Twisting (OAT) in the fully symmetric Dicke manifold, implemented following Eq. (9) of the main paper.

- Full Exact Diagonalization (ED) of the Bose–Hubbard model, implemented following the tutorial by Zhang & Dong (Eur. J. Phys. 31, 591 (2010)).

Further equations and models from the main paper will be incorporated over time, along with analysis of entanglement and spin-squeezing properties.

---

## Overview

### **1. OAT Hamiltonian (`H_oat`)**
- Implements the reduced model for large systems in the **S = L/2** Dicke manifold.
- Maps XXZ couplings to OAT parameters via the main paper’s Eq. (5).
- Hamiltonian:
  \[
  H_{\text{OAT}} = \frac{2(J_\parallel + J_\perp)}{L-1} S^2 + \frac{2(J_\parallel - J_\perp)}{L-1} S_x^2 + 2 J_\Omega S_x
  \]
- Uses collective spin operators (`Sx`, `Sy`, `Sz`) in the Dicke basis to evolve the state.

### **2. Exact Diagonalization (ED)**
- Implements full ED of the Bose–Hubbard model.
- Follows Zhang’s approach:
  - Generate Fock basis for \( N \) bosons on \( M \) sites.
  - Construct the Hamiltonian using efficient hashing & sparse matrix storage.
  - Diagonalize using **ARPACK** (`eigs`) for ground state & low-lying excitations.
- Computes observables:
  - Condensate fraction \( f_c \) (Penrose–Onsager criterion)
  - Off-diagonal long-range order (ODLRO)
  - On-site number fluctuations

---

## Features

- **Two modeling approaches**:
  - **OAT** → Efficient for large \( L \) with symmetric spin states.
  - **ED** → Full solution for small–medium \( M \), exact but more computationally expensive.
- **Time evolution**: Uses eigen-decomposition for exact unitary evolution.
- **Observable computation**: Contrast \( C_x, C_z \), condensate fraction, correlations, variances.
- **Parameter mapping**: XXZ model → OAT couplings via analytical formulas.

---

## ⚙️ Requirements

- **C++17 or newer**
- [Armadillo](http://arma.sourceforge.net/) (linear algebra)
- [ARPACK](https://www.caam.rice.edu/software/ARPACK/) (eigenvalue problems)

---

## ▶️ Compilation

### **Build**
```bash
mkdir build && cd build
cmake ..
make
```

### **Run**
```bash
./oat_simulation         # Runs OAT model
./ed_simulation          # Runs full exact diagonalization
```

---

## 📂 Project Structure

```
├── src/
│   ├── H_OAT.cpp              # OAT Hamiltonian & evolution
│   ├── spdm.cpp               # Exact diagonalization of Bose–Hubbard
│   └── ...
├── data/                    # Output CSV files
├── plots/                   # Scripts to visualize results
├── CMakeLists.txt
└── README.md
```

---

##  Example Output

### OAT Contrast (Main Paper Fig. 2d parameters)
```
L = 8, J/U = 50, ϕ = 0.5, Ω/J = 15, tmax = 200
```
Produces CSV columns:
```
t, Cx, Cz
```
which can be plotted to reproduce the contrast dynamics.

### ED Observables
Outputs:
- `fc_vs_UJ.csv` – condensate fraction vs \( U/J \)
- `rho_far_vs_UJ.csv` – farthest-site correlation
- `sigma_vs_UJ.csv` – on-site number variance

---

## 📚 References

- Main paper: *Phys. Rev. Lett. 126, 113401* – Eq. (9) for OAT Hamiltonian.  
- Tutorial: J. M. Zhang, R. X. Dong, *Eur. J. Phys.* **31**, 591 (2010) – Basis generation & ED method.  
- Code structure inspired by Zhang’s hashing & sparse-matrix techniques.
