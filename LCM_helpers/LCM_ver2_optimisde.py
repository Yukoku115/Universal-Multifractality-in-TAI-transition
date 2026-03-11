import os
import sys
import time
import numpy as np
from collections import deque
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from joblib import Parallel, delayed

# ==============================================================================
# Colormap & Plot Styling
# ==============================================================================
tai_cmap = LinearSegmentedColormap.from_list(
    "tai_cmap",
    [
        (0.00, "#08306b"),
        (0.25, "#00e5ff"),
        (0.50, "#ffffff"),
        (0.75, "#ff00ff"),
        (1.00, "#5e00ff"),
    ]
)
mpl.colormaps.register(tai_cmap, name="tai_cmap", force=True)

plt.rcParams.update({
    'font.size': 12,
    'font.family': 'serif',
    'mathtext.fontset': 'stix',
    'axes.linewidth': 1.2,
    'lines.linewidth': 1.5,
    'figure.dpi': 150
})

# ==============================================================================
# Lattice
# ==============================================================================
_NNN_VECS_POS = np.array([
    [ np.sqrt(3),   0.0],
    [-np.sqrt(3)/2, 1.5],
    [-np.sqrt(3)/2,-1.5],
])
_NNN_VECS_NEG = -_NNN_VECS_POS


def build_honeycomb(Nx, Ny, a=1.0):
    """Generates a finite honeycomb lattice."""
    a1 = np.array([np.sqrt(3) * a, 0])
    a2 = np.array([np.sqrt(3) / 2 * a, 3 / 2 * a])

    i_idx, j_idx = np.meshgrid(np.arange(Nx), np.arange(Ny), indexing='ij')
    R = (i_idx[..., np.newaxis] * a1 + j_idx[..., np.newaxis] * a2).reshape(-1, 2)

    x_coords = np.empty(2 * len(R))
    y_coords = np.empty(2 * len(R))
    x_coords[0::2] = R[:, 0];  y_coords[0::2] = R[:, 1]
    x_coords[1::2] = R[:, 0];  y_coords[1::2] = R[:, 1] + a

    sublattices = np.empty(2 * len(R), dtype='U1')
    sublattices[0::2] = 'A'
    sublattices[1::2] = 'B'

    return x_coords, y_coords, sublattices


# ==============================================================================
# Geometry precomputation  —  called ONCE before the sweep
# ==============================================================================
def precompute_geometry(x, y, sublattices):
    coords      = np.column_stack((x, y))
    dist_matrix = cdist(coords, coords)
    diff_vecs   = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]

    nn_mask  = np.isclose(dist_matrix, 1.0,        atol=1e-3)
    nnn_mask = np.isclose(dist_matrix, np.sqrt(3), atol=1e-3)

    sign_matrix = np.zeros((len(x), len(x)), dtype=float)
    for v in _NNN_VECS_POS:
        sign_matrix[np.all(np.isclose(diff_vecs, v, atol=1e-3), axis=2)] =  1.0
    for v in _NNN_VECS_NEG:
        sign_matrix[np.all(np.isclose(diff_vecs, v, atol=1e-3), axis=2)] = -1.0
    sign_matrix[sublattices == 'B', :] *= -1

    active_nnn = nnn_mask & (sign_matrix != 0)

    return nn_mask, active_nnn, sign_matrix


# ==============================================================================
# Hamiltonian  —  geometry passed in, only physics recomputed per (M,W)
# ==============================================================================
def build_H_from_geometry(x, y, sublattices, nn_mask, active_nnn, sign_matrix,
                           t1=1.0, t2=1/3, M=0.0, W=0.0, disorder_type="Anderson"):
    N = len(x)
    H = np.zeros((N, N), dtype=complex)

    phase_matrix = np.zeros((N, N), dtype=float)
    if disorder_type == "Magnetic" and W > 0:
        raw          = np.random.uniform(-W / 2, W / 2, (N, N))
        phase_matrix = (raw - raw.T) / 2

    H[nn_mask]    = -t1 * np.exp(1j * phase_matrix[nn_mask])
    H[active_nnn] = (-t2
                     * np.exp(1j * sign_matrix[active_nnn] * np.pi / 2)
                     * np.exp(1j * phase_matrix[active_nnn]))

    diag = np.where(sublattices == 'A', M, -M)
    if disorder_type == "Anderson" and W > 0:
        diag = diag + np.random.uniform(-W / 2, W / 2, N)
    H[np.arange(N), np.arange(N)] += diag

    return H


# ==============================================================================
# Local Chern Marker
# ==============================================================================
def calc_local_chern(H, x, y, E_Fermi=0.0):
    evals, evecs = np.linalg.eigh(H)
    occ = evecs[:, evals < E_Fermi]
    P   = occ @ occ.conj().T

    PxP = (P * x[np.newaxis, :]) @ P
    PyP = (P * y[np.newaxis, :]) @ P

    return np.real(np.diag(-2j * np.pi * (PxP @ PyP - PyP @ PxP)))


# ==============================================================================
# Per-row worker  —  used by joblib
# ==============================================================================
def compute_row(i, M_val, W_vals, done_row, x, y, sub,
                nn_mask, active_nnn, sign_matrix, bulk_mask):
    row = np.full(len(W_vals), np.nan)
    for j, W_val in enumerate(W_vals):
        if done_row[j]:
            continue
        H       = build_H_from_geometry(x, y, sub, nn_mask, active_nnn, sign_matrix,
                                        M=M_val, W=W_val, disorder_type="Anderson")
        C_local = calc_local_chern(H, x, y)
        row[j]  = float(np.mean(C_local[bulk_mask]))
    return i, row


# ==============================================================================
# Parameters  —  edit these
# ==============================================================================
Nx_val, Ny_val = 30, 30
grid_size      = 50
W_vals         = np.linspace(0.0, 9.5, grid_size)
M_vals         = np.linspace(1.3, 2.8, grid_size)
CHECKPOINT     = "phase_diagram_checkpoint.npz"
N_JOBS         = 4   # -1 = all available CPU cores; set to e.g. 4 to limit

# ==============================================================================
# Build lattice, bulk mask, and geometry  —  all done once
# ==============================================================================
x_lat, y_lat, sub_lat = build_honeycomb(Nx=Nx_val, Ny=Ny_val)

dist_from_center = np.sqrt((x_lat - x_lat.mean())**2 + (y_lat - y_lat.mean())**2)
bulk_mask        = dist_from_center < np.max(dist_from_center) * 0.40

print("Precomputing lattice geometry...", end=" ", flush=True)
nn_mask, active_nnn, sign_matrix = precompute_geometry(x_lat, y_lat, sub_lat)
print("done.")

# ==============================================================================
# Checkpoint: load existing progress or start fresh
# ==============================================================================
if os.path.exists(CHECKPOINT):
    data = np.load(CHECKPOINT)

    grid_matches = (
        data["W_vals"].shape == W_vals.shape   and
        data["M_vals"].shape == M_vals.shape   and
        np.allclose(W_vals, data["W_vals"])    and
        np.allclose(M_vals, data["M_vals"])
    )

    if not grid_matches:
        print("Grid parameters have changed — discarding old checkpoint and starting fresh.")
        phase_map = np.full((len(M_vals), len(W_vals)), np.nan, dtype=float)
        done_mask = np.zeros((len(M_vals), len(W_vals)), dtype=bool)
    else:
        phase_map = data["phase_map"]
        done_mask = data["done_mask"].astype(bool)
        completed = int(done_mask.sum())
        print(f"Resuming from checkpoint — {completed}/{grid_size**2} points already done.")
else:
    phase_map = np.full((len(M_vals), len(W_vals)), np.nan, dtype=float)
    done_mask = np.zeros((len(M_vals), len(W_vals)), dtype=bool)
    print("No checkpoint found — starting fresh.")

# ==============================================================================
# Main sweep  —  parallelised row-by-row, checkpoint after each batch
# ==============================================================================
pending_rows  = [i for i in range(len(M_vals)) if not done_mask[i].all()]
total_points  = grid_size ** 2

# Batch size matches N_JOBS (or actual core count if N_JOBS=-1)
BATCH_SIZE    = N_JOBS if N_JOBS > 0 else (os.cpu_count() or 4)

# Rolling window for ETA: stores (points_completed, seconds) per batch
# Capped at 10 entries so the rate reflects recent pace, not early slow batches
timing_window = deque(maxlen=10)

print(f"Computing {len(pending_rows)} remaining rows (out of {len(M_vals)}) "
      f"using {N_JOBS if N_JOBS != -1 else 'all'} CPU cores...\n")

for batch_start in range(0, len(pending_rows), BATCH_SIZE):

    batch         = pending_rows[batch_start : batch_start + BATCH_SIZE]
    t_batch_start = time.time()
    n_batches     = (len(pending_rows) + BATCH_SIZE - 1) // BATCH_SIZE
    batch_num     = batch_start // BATCH_SIZE + 1
    sys.stdout.write(f"\r  Batch {batch_num}/{n_batches} — rows {batch[0]}–{batch[-1]} running...{' '*20}")
    sys.stdout.flush()

    results = Parallel(n_jobs=N_JOBS)(
        delayed(compute_row)(
            i, M_vals[i], W_vals, done_mask[i],
            x_lat, y_lat, sub_lat,
            nn_mask, active_nnn, sign_matrix, bulk_mask
        )
        for i in batch
    )

    batch_points = 0
    for i, row in results:
        for j, val in enumerate(row):
            if not np.isnan(val):
                phase_map[i, j] = val
                done_mask[i, j] = True
                batch_points    += 1

    np.savez(CHECKPOINT, W_vals=W_vals, M_vals=M_vals,
             phase_map=phase_map, done_mask=done_mask)

    # ── Rolling-window ETA ────────────────────────────────────────────────
    timing_window.append((batch_points, time.time() - t_batch_start))

    completed  = int(done_mask.sum())
    remaining  = total_points - completed
    total_t    = sum(t for _, t in timing_window)
    total_pts  = sum(p for p, _ in timing_window)

    if total_pts > 0 and total_t > 0:
        rate    = total_pts / total_t               # points/sec, recent average
        eta_sec = remaining / rate
        h, rem  = divmod(int(eta_sec), 3600)
        m, s    = divmod(rem, 60)
        eta_str = f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"
    else:
        eta_str = "--:--"

    progress = completed / total_points
    bar      = '=' * int(30 * progress) + '-' * (30 - int(30 * progress))
    sys.stdout.write(
        f"\r[{bar}] {completed}/{total_points} ({progress*100:.1f}%)  "
        f"ETA {eta_str}   "
    )
    sys.stdout.flush()

print("\nDone. Generating plot...")

# ==============================================================================
# Plot
# ==============================================================================
fig, ax = plt.subplots(figsize=(5, 5))
mesh = ax.imshow(
    phase_map,
    origin="lower",
    extent=[W_vals.min(), W_vals.max(), M_vals.min(), M_vals.max()],
    cmap="tai_cmap",
    vmin=-0.2, vmax=1.2,
    aspect="auto"
)
ax.set_box_aspect(1)
ax.set_title("Topological Phase Diagram")
ax.set_xlabel("Disorder Strength ($W$)")
ax.set_ylabel("Staggered Mass ($M$)")
fig.colorbar(mesh, ax=ax, label='Bulk Chern Number')

fig.text(
    0.5, -0.02,
    f"Honeycomb lattice: $N_x = {Nx_val}$, $N_y = {Ny_val}$ "
    f"({2 * Nx_val * Ny_val} sites) · Grid: {grid_size}×{grid_size}",
    ha='center', va='top', fontsize=9, color='0.4',
    fontstyle='italic'
)

plt.tight_layout()
plt.show()
