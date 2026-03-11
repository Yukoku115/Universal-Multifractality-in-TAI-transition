import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

tai_cmap = LinearSegmentedColormap.from_list(
    "tai_cmap",
    [
        (0.00, "#3831f6"),
        (0.25, "#00e5ff"),
        (0.50, "#ffffff"),
        (0.75, "#9211ee"),
        (1.00, "#fe0ea6"),
    ]
)
mpl.colormaps.register(tai_cmap, name="tai_cmap", force=True)

# ── Load checkpoint ───────────────────────────────────────────────────────────
data      = np.load("phase_diagram_checkpoint.npz", allow_pickle=True)
W_vals    = data["W_vals"]
M_vals    = data["M_vals"]
phase_map = data["phase_map"]
done_mask = data["done_mask"]

# Grid info
grid_W      = len(W_vals)
grid_M      = len(M_vals)
n_completed = int(done_mask.sum())
n_total     = grid_W * grid_M

# System size from checkpoint (saved by LCM_vectorised.py)
Nx      = int(data["Nx"])      if "Nx"      in data else "?"
Ny      = int(data["Ny"])      if "Ny"      in data else "?"
N_sites = int(data["N_sites"]) if "N_sites" in data else "?"
N_bulk  = int(data["N_bulk"])  if "N_bulk"  in data else "?"

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5.5, 5.5))
fig.subplots_adjust(bottom=0.12, top=0.93)


mesh = ax.imshow(
    phase_map,
    origin="lower",
    extent=[W_vals.min(), W_vals.max(), M_vals.min(), M_vals.max()],
    cmap="tai_cmap",
    vmin=-0.2, vmax=1.2,
    aspect="auto"
)
ax.set_box_aspect(1)
ax.set_title("Topological Phase Diagram", pad=10)
ax.set_xlabel("Disorder Strength ($W$)")
ax.set_ylabel("Staggered Mass ($M$)")
fig.colorbar(mesh, ax=ax, label='Bulk Local Chern Marker $\\langle C \\rangle$')

# ── Parameter caption ─────────────────────────────────────────────────────────
caption = (
    f"System: ${Nx} \\times {Ny}$ unit cells "
    f"($N_{{\\mathrm{{sites}}}} = {N_sites}$, $N_{{\\mathrm{{bulk}}}} = {N_bulk}$)  |  "
    f"Grid: ${grid_M} \\times {grid_W}$  "
    f"($W \\in [{W_vals.min():.2f},\\,{W_vals.max():.2f}]$, "
    f"$M \\in [{M_vals.min():.2f},\\,{M_vals.max():.2f}]$)\n"
    f"Completed: {n_completed}/{n_total}"
)

fig.text(
    0.5, 0.02,
    caption,
    ha='center', va='bottom',
    fontsize=8, fontstyle='italic',
    color='#333333',
    transform=fig.transFigure,
)


plt.savefig("phase_diagram.png", dpi=200, bbox_inches='tight')
plt.show()
