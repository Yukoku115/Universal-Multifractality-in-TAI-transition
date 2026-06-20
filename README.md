# Multifractal Analysis of the Topological Anderson Insulator (TAI) Transition

This repository contains tight-binding simulation pipelines to compute topological phase diagrams and perform multifractal scaling analysis on 2D honeycomb lattices, utilizing the **Haldane model** (spinless) and the **Kane-Mele model** (spinful). 

The simulations use the real-space **Local Chern Marker (LCM)** and the **Local Spin Chern Marker (LSCM)** to map phase boundaries under Anderson (electrostatic) disorder and Rashba spin-orbit coupling (SOC) disorder. At the phase boundaries, a 2D box-counting algorithm characterizes the multifractal dimensions to identify the universality classes of the metal-insulator transition.

---

## Directory Structure

```directory
├── final_plots/           # High-resolution phase diagrams & multifractal curves (Phase 3 results)
├── archive/               # Legacy files consolidated to keep the root directory clean
│   ├── old_plots/         # Deprecated plots from Phase 1 & Phase 2
│   ├── old_scripts/       # Legacy python scripts and helper tools
│   └── intermediate_runs/ # Intermediate plots and trial sweeps
├── LCM_ploter/            # Core physics library
│   ├── __init__.py        # Exported functions
│   ├── core.py            # Hamiltonian builders (Haldane & Kane-Mele) and topological markers
│   ├── multifractal.py    # Box-counting algorithms and parabolic gamma fitting
│   └── plotting.py        # Shared plotting routines for phase diagrams & spectra
├── main_haldane.py        # Sweep script for the spinless Haldane model
├── main_km.py             # Sweep script for the spinful Kane-Mele model
├── check_progress.py      # Helper utility to inspect checkpoint progress
└── nanoribbon.py          # Helper utility for nanoribbon edge-state calculations
```

---

## Physics and Models

### 1. The Models
* **Haldane Model:** A spinless model on a honeycomb lattice that breaks time-reversal symmetry (TRS) locally using next-nearest-neighbor (NNN) hopping phases $\phi = \pm \pi/2$ while maintaining zero net macroscopic magnetic flux.
* **Kane-Mele Model:** A spinful model consisting of two time-reversal partners of the Haldane model. The NNN hopping acts as intrinsic spin-orbit coupling (SOC). It preserves TRS and exhibits the Quantum Spin Hall (QSH) effect.

### 2. Helicity Projection (Prodan's Method)
When Rashba SOC ($t_R$) is introduced, spin $S_z$ is no longer conserved in the lab frame due to spin-flipping terms. To calculate the topological invariant (the Spin Chern Number), the code uses **Emil Prodan's generalized method**:
1. The lab spin operator $S_z$ is projected onto the occupied valence band: $P_{\text{sub}} = P S_z P$.
2. Diagonalizing $P_{\text{sub}}$ yields eigenvalues (the spin spectrum) that are split by a gap around $0$.
3. We separate the occupied space into two orthogonal helicity projection operators based on the sign of the eigenvalues: $P_+$ (positive helicity, $s_n > 0.1$) and $P_-$ (negative helicity, $s_n < -0.1$).
4. The Local Spin Chern Marker is calculated using the gapped projection operator $P_+$:
   $$c_s(\mathbf{r}) = -2\pi i \langle \mathbf{r} | [P_+ x P_+, P_+ y P_+] | \mathbf{r} \rangle$$

### 3. Relativistic Spin Disorder Scaling
Substrate roughness and lattice defects cause spatial fluctuations in the vertical electric field $E_z(\mathbf{r})$, leading to random fluctuations in the Rashba hopping terms (bond-dependent spin disorder $W_R$). The code scales this disorder perturbative relative to the on-site Anderson disorder $W_A$:
$$W_R = \frac{0.1}{\sqrt{3}} W_A$$
* **The 0.1 Factor:** Reflects the relativistic scaling of spin-orbit coupling (typically orders of magnitude smaller than electrostatic hopping energies).
* **The $1/\sqrt{3}$ Factor:** Coordination normalization. Since each site on the honeycomb lattice has $z = 3$ nearest-neighbor bonds, the variance of the independent bond fluctuations adds up. Dividing the bond amplitude by $\sqrt{3}$ ensures the total variance felt by a scattering electron is normalized relative to the Anderson disorder.

### 4. Multifractal Exponent ($\gamma$) and Universality Classes
At the critical boundary of the transition, the wavefunctions exhibit multifractal scaling. We calculate the mass exponent $\tau(q) = D_q(q-1)$ and fit the anomalous dimension $\gamma$ using the parabolic approximation:
$$\tau(q) = d(q-1) + \gamma q(1-q) \quad (\text{for } q \le 4)$$
* **Unitary Class** ($t_R = 0$, spin conserved): Theoretical limit $\gamma \approx 0.26$.
* **Symplectic Class** ($t_R > 0$, spin mixed): Theoretical limit $\gamma \approx 0.17$.

---

## Installation & Setup

### Dependencies
The codebase requires Python 3.8+ and the following scientific packages:
* `numpy`
* `scipy`
* `matplotlib`
* `joblib` (for multi-core parallelization)

Install them via `pip`:
```bash
pip install numpy scipy matplotlib joblib
```

---

## Running the Simulations

### 1. Mapping Phase Diagrams
To sweep the staggered mass $M$ and disorder $W$ to map out the topological phases:
* **Haldane Model:**
  ```bash
  python main_haldane.py
  ```
* **Kane-Mele Model:**
  ```bash
  python main_km.py
  ```
*Both scripts support automatic checkpointing. If interrupted, they will automatically save progress to a `.npz` file and resume from the last completed row upon restart.*

### 2. Analyzing Checkpoint Progress
To check the number of completed rows and parameters inside a checkpoint file:
```bash
python check_progress.py
```

### 3. Parameters Configuration
Key physical and numerical parameters can be configured directly at the top of the `main_km.py` and `main_haldane.py` files:
* `Nx_val, Ny_val`: Dimensions of the honeycomb lattice (e.g., $15 \times 15$ or $20 \times 20$).
* `grid_size`: Resolution of the $(M, W)$ phase diagram (e.g., $30 \times 30$ or $50 \times 50$ points).
* `tR_base`: Baseline uniform Rashba SOC strength.
* `N_JOBS`: Number of parallel CPU cores to use for the joblib workers (set to `-1` for all available cores).
