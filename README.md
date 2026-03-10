# Universal Multifractality in TAI Transition

This repository contains Python simulations of the **Topological Anderson Insulator (TAI)** phase within the Haldane model on a honeycomb lattice. It explores how different types of quantum chaos (diagonal vs. off-diagonal disorder) affect topological phase transitions and the emergence of multifractal wavefunctions.

## 📌 Project Overview
The Topological Anderson Insulator (TAI) is an anomalous quantum phase where adding disorder to a trivial band insulator *creates* a topological state rather than destroying it. 

This project simulates the real-space evolution of the Local Chern Marker $C(r)$ to map the topological phase diagram. Specifically, it tests the theoretical hypothesis that the TAI phase is strictly dependent on the *type* of disorder applied to the lattice. We contrast standard **Anderson (diagonal) disorder** against a novel **Magnetic Phase (off-diagonal) disorder**.

## ⚛️ Theoretical Background
In the clean Haldane model, the topological state is a "tug-of-war" between two parameters:
1. **Staggered Mass ($M$):** The on-site energy difference between Sublattices A and B. This drives the system toward a trivial Band Insulator (BI).
2. **Complex Hopping ($t_2$):** A synthetic, microscopic magnetic flux on next-nearest-neighbor hopping. This breaks time-reversal symmetry and drives the system toward a topological Chern Insulator (CI).

When $M$ dominates, the system is trivial. However, adding disorder causes second-order quantum scattering, which renormalizes the effective mass ($\overline{M}$). 

### The Core Experiment: Diagonal vs. Off-Diagonal Disorder
Using the Self-Consistent Born Approximation, the renormalized mass depends on the Pauli matrix algebra of the specific disorder:

* **Anderson Disorder (Scalar On-Site):** Modeled as random voltages added to the atoms ($\sigma_0$). Because $\sigma_0 \sigma_z \sigma_0 = \sigma_z$, the scattering provides a **negative** correction to the mass. The gap shrinks, the bands invert, and the TAI phase emerges.
* **Magnetic Phase Disorder (Complex Hopping):** Modeled as random Peierls phases ($e^{i\phi}$) applied to the hopping bridges ($\sigma_x, \sigma_y$). Because off-diagonal matrices anti-commute with the mass matrix ($\sigma_x \sigma_z \sigma_x = -\sigma_z$), the scattering provides a **positive** correction to the mass. Destructive quantum interference isolates the sublattices, widening the gap and strictly forbidding the TAI transition.

## 📊 Key Results

### 1. The TAI Phase Diagram (Anderson Disorder)
*(Drag and drop your magenta Anderson Phase Diagram image here)*
> **Description:** Sweeping Staggered Mass ($M$) vs. Disorder Strength ($W$). The emergence of the magenta TAI dome proves that scalar diagonal noise successfully shrinks the bandgap to invert the bands, shifting the system from BI $\rightarrow$ TAI $\rightarrow$ Anderson Insulator (AI).

### 2. The Trivial Phase Diagram (Magnetic Phase Disorder)
*(Drag and drop your cyan Magnetic Phase Diagram image here)*
> **Description:** Applying random magnetic phases to the hopping bridges perfectly destroys the topological engine. The effective mass is driven wider, preventing band inversion. The system transitions directly from BI $\rightarrow$ AI, skipping the topological and multifractal critical states entirely.

### 3. Real-Space Local Chern Marker Evolution
*(Drag and drop your side-by-side real-space lattice plots here)*
> **Description:** Direct spatial visualization of the Local Chern Marker $C(r)$. 
> * **Left Panel ($M=1.00$):** A clean Chern Insulator gradually destroyed by chaos.
> * **Right Panel ($M=1.89$):** A trivial Band Insulator blooming into a percolating TAI network under moderate Anderson disorder, before finally succumbing to Anderson Localization.

## 🚀 How to Run
The codebase relies on standard Python scientific libraries.
1. Clone the repository: `git clone https://github.com/Yukoku115/Universal-Multifractality-in-TAI-transition.git`
2. Install dependencies: `pip install numpy matplotlib scipy`
3. Run the main script to generate the Hamiltonian matrices and plot the real-space or phase diagram sweeps. Adjust the `disorder_type` parameter to `"Anderson"` or `"Magnetic"` to replicate the comparative results.

## 📚 References
1. Haldane, F. D. M. (1988). *Model for a Quantum Hall Effect without Landau Levels*. Physical Review Letters.
2. Song, J., Liu, H., Jiang, H., Sun, Q.-f., & Xie, X. C. (2012). *Dependence of topological Anderson insulator on the type of disorder*. Physical Review B, 85(19), 195125.
