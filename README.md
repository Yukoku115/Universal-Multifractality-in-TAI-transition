# Universal Multifractality in Topological Anderson Insulator Transition

A comprehensive computational framework for studying **Topological Anderson Insulators (TAI)** and multifractal wave functions in the Haldane and Kane-Mele models on honeycomb lattices. This project investigates how different types of disorder—Anderson (diagonal) vs. Magnetic Phase (off-diagonal)—drive topological phase transitions and the emergence of multifractal eigenstates.

## 🎯 Project Overview

The **Topological Anderson Insulator** is a paradigmatic example of a disorder-induced topological phase: adding disorder to a trivial band insulator can *create* rather than destroy topological states. This repository provides:

- **Phase diagram computations** for multiple quantum models (Haldane, Kane-Mele)
- **Multifractal analysis** of wave function statistics across the TAI transition
- **Real-space visualization** of local topological markers
- **Nanoribbon geometry support** for edge state calculations
- **Checkpoint/resume capability** for long-running calculations

## 🏗️ Project Structure

```
├── main_haldane.py               # Haldane model phase diagrams (Anderson vs. Magnetic disorder)
├── main_km.py                    # Kane-Mele model with Rashba spin-orbit coupling
├── nanoribbon.py                 # Band structure analysis for zigzag/armchair nanoribbons
├── check_progress.py             # Progress monitoring and partial visualization
├──
│
├── LCM_ploter/                   # Core computational library
│   ├── core.py                   # Hamiltonian construction, lattice generation
│   ├── multifractal.py           # Box-counting multifractal analysis (D_q spectrum)
│   └── plotting.py               # Visualization utilities
│
├── LCM_phase/                    # Pre-computed phase diagrams and analysis
│   ├── LCM_ver2_vectorised_noplot
│   └── LCM_ver2_verctorised.py
│
├── plots/                        # Output directory for generated plots
├── LCM_plots/                    # LCM-specific visualizations
└── papers/                       # Reference literature
```

## ⚛️ Physical Models

### 1. **Haldane Model (Spinless)**
The Haldane Hamiltonian on a honeycomb lattice includes:
- **Nearest-neighbor hopping** ($t_1 = 1$): connects A and B sublattices
- **Next-nearest-neighbor hopping with flux** ($t_2 = 1/3$): breaks time-reversal symmetry
- **Staggered mass** ($M$): on-site energy difference between sublattices
- **Disorder** ($W$): varies by type (Anderson or Magnetic Phase)

**Key parameters:** Sweep $(M, W)$ space to map the phase diagram in terms of Chern marker $C(r)$.

### 2. **Kane-Mele Model (Spinful)**
Extends Haldane with spin degrees of freedom:
- Intrinsic spin-orbit coupling with opposite signs for spin-up/down
- Rashba spin-orbit coupling ($t_R$) on nearest-neighbor bonds
- Optional Rashba disorder for further control
- Computes **spin Chern number** $C_s = (C_{\uparrow} - C_{\downarrow})/2$

**Key parameters:** $M$, $W_{\text{Anderson}}$, $W_{\text{Rashba}}$, $t_R$.

### 3. **Nanoribbon Geometries**
Support for:
- **Zigzag nanoribbons:** A and B staggered edges with topological end states
- **Armchair nanoribbons:** 4-atom repeating unit with different band structure

## 🔬 Analysis Capabilities

### Phase Diagram Computation
- **Parallel CPU computation** using `joblib` for grid sweeps
- **Checkpoint/resume** system preserves progress across sessions
- **Bulk masking** to avoid edge artifacts
- Configurable grid resolution (typically $50 \times 50$ or larger)

### Multifractal Analysis
Computes the full multifractal spectrum $D_q$ via 2D spatial box-counting:
- **Probability moments** $P_q(\lambda)$ across box sizes
- **Extraction of scaling exponent** $\tau(q)$ from log-log fits
- **Generalized dimensions** $D_q = \tau(q) / (q - 1)$
- Supports arbitrary moment orders $q \in \mathbb{R}$

### Real-Space Visualization
- Local Chern marker evolution across $(M, W)$ space
- Wavefunction probability density heatmaps
- Band structure along high-symmetry paths (nanoribbons)

## 🚀 Quick Start

### Installation
```bash
# Clone the repository
git clone <repository-url>
cd "Multifractal in TAI transition"

# Install dependencies
pip install numpy scipy matplotlib joblib
```

### Running Phase Diagram Calculations

**Haldane Model (Anderson Disorder):**
```bash
python main_haldane.py
```
Computes a $50 \times 50$ phase diagram for $M \in [1.3, 2.8]$ and $W \in [0.0, 9.5]$.  
Saves checkpoint to `phase_diagram_checkpoint.npz`.

**Kane-Mele Model:**
```bash
python main_km.py
```
Similar grid sweep, saving to `km_phase_diagram_checkpoint.npz`.  
Edit `sweep_type` in the script to toggle between Anderson and Rashba disorder.

### Monitoring Progress
```bash
python check_progress.py
```
Generates `partial_plot.png` showing completed points and current progress.

### Nanoribbon Band Structure
```bash
python nanoribbon.py
```
Plots band structure for zigzag or armchair geometries with interactive sliders for $M$ and $t_2$.

## 🛠️ Configuration

All main scripts expose parameters at the top:

```python
Nx_val, Ny_val = 20, 20              # System size
grid_size      = 50                  # Phase diagram resolution
W_vals         = np.linspace(0.0, 9.5, grid_size)   # Disorder range
M_vals         = np.linspace(1.3, 2.8, grid_size)   # Mass range
N_JOBS         = 4                   # Parallel CPU cores
```

Modify these values, then re-run the script. **Note:** Changing grid parameters invalidates old checkpoints automatically.

## 📊 Output Files

| File | Description |
|------|-------------|
| `phase_diagram_checkpoint.npz` | Computed Haldane phase diagram with progress tracking |
| `km_phase_diagram_checkpoint.npz` | Computed Kane-Mele phase diagram |
| `multifractal_checkpoint.npz` | Cached multifractal D_q spectra |
| `partial_plot.png` | Visual progress report |
| `plots/` | All generated figures (heatmaps, band structures, etc.) |

## 🔍 Core Library: LCM_ploter

### `core.py`
- `build_honeycomb(Nx, Ny)` — Generate finite honeycomb lattice
- `precompute_geometry(x, y, sublattices)` — Cache bond connectivity and topology
- `build_H_from_geometry(...)` — Construct Hamiltonian with disorder
- `compute_row(i, M, W_vals, ...)` — Parallel computation kernel for one parameter slice

### `multifractal.py`
- `compute_Pq(psi, atom_x, atom_y, box_sizes, q_vals, L)` — Extract probability moments
- `extract_tau_Dq(mean_Pq, log_lambda, q_vals)` — Fit D_q spectrum
- Implements equations 3–6 from multifractal literature

### `plotting.py`
- `plot_phase_diagram(phase_map, W_vals, M_vals, ...)` — Generate publication-quality heatmaps
- Customizable colormaps (tai_cmap: blue→cyan→white→magenta→pink)

## 📚 Physical Interpretation

### Anderson Disorder (Diagonal)
- **Mechanism:** On-site random potential ($\sigma_0$ Pauli matrix)
- **Effect on mass:** Negative correction → gap shrinks
- **Result:** Band inversion in trivial phase → TAI phase emerges
- **Signature:** Topological dome in phase diagram

### Magnetic Phase Disorder (Off-Diagonal)
- **Mechanism:** Random Peierls phases on hopping terms ($\sigma_x, \sigma_y$)
- **Effect on mass:** Positive correction → gap widens
- **Result:** Destructive interference prevents topological state
- **Signature:** No TAI dome; direct BI → AI transition

### Multifractal Critical Point
At the TAI transition, eigenstates exhibit **self-similar** structure:
- Wave function is neither extended nor localized
- Fractal dimensions $D_q$ vary with moment order $q$
- Maximum multifractality at criticality

## 📖 References

1. **Haldane, F. D. M.** (1988). *Model for a Quantum Hall Effect without Landau Levels*. Physical Review Letters, **61**(15), 2015.

2. **Song, J., Liu, H., Jiang, H., Sun, Q.-f., & Xie, X. C.** (2012). *Dependence of topological Anderson insulator on the type of disorder*. Physical Review B, **85**(19), 195125.

3. **Prodan, E.** (2011). *Topological Anderson insulators*. Physical Review B, **80**(8), 085119.

4. **Kane, C. L., & Mele, E. J.** (2005). *Quantum Spin Hall Effect in Graphene*. Physical Review Letters, **95**(22), 226801.

5. **Evers, F., & Mirlin, A. D.** (2008). *Anderson transitions*. Reviews of Modern Physics, **80**(3), 1355.

## ⚙️ Requirements

- **Python 3.7+**
- **NumPy** — matrix operations and Fourier transforms
- **SciPy** — eigenvalue solvers, spatial distance computations
- **Matplotlib** — plotting and visualization
- **Joblib** — parallel job scheduling

## 📝 Notes

- **Large lattices** ($Nx, Ny > 25$) increase both memory and computation time quadratically.
- **High-resolution grids** ($> 100 \times 100$) may require days of computation on multi-core systems.
- Checkpoints are **compatible** with resumed runs—modify parameters and re-run to extend calculations.
- For production runs, consider submitting to a **cluster scheduler** (SLURM/PBS) with job arrays.

## 👤 Author

Yuko K.

## 📄 License

[Specify your chosen license, e.g., MIT, GPL-3.0, etc.]
