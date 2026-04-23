import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

# Same colormap as LCM_ploter.py
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

plt.rcParams.update({
    'font.size': 12,
    'font.family': 'serif',
    'mathtext.fontset': 'stix',
    'axes.linewidth': 1.2,
    'figure.dpi': 150
})

data = np.load('km_phase_diagram_checkpoint.npz')
phase_map = data['phase_map'].copy()
done_mask = data['done_mask']
W_vals = data['W_vals']
M_vals = data['M_vals']

phase_map[~done_mask] = np.nan

completed = int(done_mask.sum())
total = done_mask.size

fig, ax = plt.subplots(figsize=(6.5, 6))
mesh = ax.imshow(
    phase_map,
    origin="lower",
    extent=[W_vals.min(), W_vals.max(), M_vals.min(), M_vals.max()],
    cmap="tai_cmap",
    vmin=0, vmax=1,
    aspect="auto"
)
ax.set_box_aspect(1)
cbar = fig.colorbar(mesh, ax=ax, pad=0.04)
cbar.set_label(r'Spin Chern Marker $C_s$', labelpad=15)
ax.set_xlabel("Disorder Strength ($W$)")
ax.set_ylabel("Staggered Mass ($M$)")
ax.set_title(f'Partial plot — {completed}/{total} points done', pad=12)

fig.tight_layout()
plt.savefig('partial_plot.png', dpi=200, bbox_inches='tight')
print(f'Saved partial_plot.png  ({completed}/{total} points, {100*completed/total:.1f}%)')