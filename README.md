# Topological Phase Diagram Simulations

Two scripts for computing phase diagrams of 2D topological insulators on a honeycomb lattice, using the Local Chern Marker (LCM) and multifractal analysis.

---

## Files

- **`main_haldane.py`** — Haldane model phase diagram + multifractal analysis
- **`main_km.py`** — Kane-Mele model phase diagram (Spin Chern Number)
- **`LCM_ploter.py`** — shared library (lattice building, Hamiltonian construction, plotting)

---

## What Each Script Does

### `main_haldane.py`
Sweeps a 2D grid of (staggered mass `M`, Anderson disorder `W`) and computes the Local Chern Marker for each point to map out the phase boundaries between:
- **CI** — Chern Insulator
- **TAI** — Topological Anderson Insulator
- **AI** — Anderson Insulator
- **BI** — Band Insulator

After the phase diagram, it runs a **multifractal analysis** (generalised dimension D2) along the phase boundaries to characterise the localisation of eigenstates near E = 0.

### `main_km.py`
Same sweep structure, but uses the **Kane-Mele Hamiltonian** (2N × 2N, with spin-up/down sectors). Computes the **Spin Chern Number** C_s = (C↑ − C↓) / 2. Supports both Anderson disorder and Rashba spin disorder sweeps, controlled by `sweep_type`.

---

## Dependencies

```
numpy
scipy
joblib
matplotlib
```

Install with:
```bash
pip install numpy scipy joblib matplotlib
```

---

## Usage

```bash
python main_haldane.py
python main_km.py
```

Both scripts support **checkpointing** — if interrupted, they resume from where they left off automatically.

---

## Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `Nx_val, Ny_val` | Lattice size | 20 × 20 |
| `grid_size` | Phase diagram resolution (grid_size²points) | 50 |
| `W_vals` | Disorder strength range | 0.0 – 9.5 |
| `M_vals` | Staggered mass range | 1.3 – 2.8 |
| `N_JOBS` | Parallel CPU cores (-1 = all) | 4 |
| `tR_base` | Rashba coupling (KM only) | 0.4 |
| `sweep_type` | `"anderson"` or `"rashba"` (KM only) | `"anderson"` |

---

## Outputs

| File | Description |
|------|-------------|
| `phase_diagram_checkpoint.npz` | Haldane checkpoint |
| `km_phase_diagram_checkpoint.npz` | Kane-Mele checkpoint |
| `multifractal_checkpoint.npz` | Multifractal checkpoint |
| `multifractal_analysis.png` | D2 map + spectra plots |
| `km_phase_diagram_{sweep_type}.png` | Kane-Mele phase diagram |