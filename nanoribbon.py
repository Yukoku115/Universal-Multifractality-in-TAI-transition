import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# --- 1. The Universal Physics Engine ---
def precompute_matrices(lattice_type, Ny=12):
    """Precomputes the spatial matrices for either ribbon type."""
    if lattice_type == 'zigzag':
        a_trans = np.array([np.sqrt(3), 0])
        a2 = np.array([np.sqrt(3)/2, 1.5])
        
        coords, sublattices = [], []
        for j in range(Ny):
            R = j * a2
            coords.append(R + np.array([0, 0]));   sublattices.append('A')
            coords.append(R + np.array([0, 1.0])); sublattices.append('B')
            
        # Topological vectors for Zigzag orientation, the three directions an e can hop to its NNN
        v_pos_A = [a_trans, a2 - a_trans, -a2]
        
    elif lattice_type == 'armchair':
        a_trans = np.array([3.0, 0])
        a_stack = np.array([0, np.sqrt(3)])
        
        coords, sublattices = [], []
        for j in range(Ny):
            R = j * a_stack
            # Armchair requires a 4-atom repeating unit cell
            coords.append(R + np.array([0, 0]));                 sublattices.append('A')
            coords.append(R + np.array([0.5, np.sqrt(3)/2]));    sublattices.append('B')
            coords.append(R + np.array([1.5, np.sqrt(3)/2]));    sublattices.append('A')
            coords.append(R + np.array([2.0, 0]));               sublattices.append('B')
            
        # Topological vectors for Armchair orientation (rotated 90 degrees)
        v_pos_A = [np.array([1.5, -np.sqrt(3)/2]), np.array([0, np.sqrt(3)]), np.array([-1.5, -np.sqrt(3)/2])]

    v_neg_A = [-v for v in v_pos_A]
    coords = np.array(coords)
    N_atoms = len(coords)
    k_values = np.linspace(-np.pi, np.pi, 80)

    # Base Matrix Skeletons
    H_M_base = np.zeros((N_atoms, N_atoms), dtype=complex)
    for i in range(N_atoms):
        H_M_base[i, i] = 1.0 if sublattices[i] == 'A' else -1.0

    H_t1_list, H_t2_list = [], []
    for k in k_values:
        H_t1_k = np.zeros((N_atoms, N_atoms), dtype=complex)
        H_t2_k = np.zeros((N_atoms, N_atoms), dtype=complex)
        
        for i in range(N_atoms):# i is the target atom
            for j in range(N_atoms):# j is the source atom
                for p in [-1, 0, 1]:
                    dist_vec = coords[i] - (coords[j] + p * a_trans)
                    dist = np.linalg.norm(dist_vec)
                    
                    if np.isclose(dist, 1.0, atol=1e-3):
                        H_t1_k[i, j] += -1.0 * np.exp(1j * k * p)
                        
                    if np.isclose(dist, np.sqrt(3), atol=1e-3):
                        sign = 0
                        for v in v_pos_A:
                            if np.allclose(dist_vec, v, atol=1e-3): sign = 1
                        for v in v_neg_A:
                            if np.allclose(dist_vec, v, atol=1e-3): sign = -1
                        if sublattices[i] == 'B': sign *= -1 
                            
                        if sign != 0:
                            H_t2_k[i, j] += -1.0 * np.exp(1j * sign * np.pi/2) * np.exp(1j * k * p)
                            
        H_t1_list.append(H_t1_k)
        H_t2_list.append(H_t2_k)

    return k_values, N_atoms, H_M_base, H_t1_list, H_t2_list

def fast_calc(M, t2, k_vals, N_at, H_M, H_t1, H_t2):
    """Calculates bands using precomputed matrices."""
    all_energies = np.zeros((len(k_vals), N_at))
    for idx in range(len(k_vals)):
        H_total = H_t1[idx] + (t2 * H_t2[idx]) + (M * H_M)
        all_energies[idx, :] = np.linalg.eigvalsh(H_total)
    return all_energies

# --- 2. User Input ---
print("\n" + "="*50)
print("HALDANE MODEL NANORIBBON SIMULATOR")
print("="*50)
choice = input("Enter ribbon type to plot ('zigzag', 'armchair', or 'both'): ").strip().lower()

if choice not in ['zigzag', 'armchair', 'both']:
    print("Invalid choice. Defaulting to 'both'.")
    choice = 'both'

# --- 3. UI and Plotting Setup ---
init_M = 1.48
init_t2 = 0.33

configs = []
if choice in ['zigzag', 'both']:
    print("Precomputing Zigzag matrices...")
    configs.append(('Zigzag', *precompute_matrices('zigzag')))
if choice in ['armchair', 'both']:
    print("Precomputing Armchair matrices...")
    configs.append(('Armchair', *precompute_matrices('armchair')))

# Create dynamic figure size based on choice
fig, axes = plt.subplots(1, len(configs), figsize=(6 * len(configs), 7))
if len(configs) == 1: axes = [axes] # Ensure axes is iterable
plt.subplots_adjust(bottom=0.3)

plot_data = [] # Store lines for quick updating
for ax, (name, k_vals, N_at, H_M, H_t1, H_t2) in zip(axes, configs):
    bands = fast_calc(init_M, init_t2, k_vals, N_at, H_M, H_t1, H_t2)
    lines = []
    for n in range(bands.shape[1]):
        line, = ax.plot(k_vals, bands[:, n], color='black', linewidth=1)
        lines.append(line)
    
    ax.set_title(f"{name} Ribbon")
    ax.set_xlabel("Momentum (k)")
    ax.set_ylabel("Energy")
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(-4, 4)
    ax.set_xticks([-np.pi, 0, np.pi], [r'$-\pi/a$', '0', r'$\pi/a$'])
    ax.grid(True, linestyle='--', alpha=0.5)
    
    plot_data.append((lines, k_vals, N_at, H_M, H_t1, H_t2))

# Sliders
ax_M = plt.axes([0.2, 0.15, 0.65, 0.03])
ax_t2 = plt.axes([0.2, 0.08, 0.65, 0.03])
slider_M = Slider(ax_M, 'Mass (M)', 0.0, 3.0, valinit=init_M)
slider_t2 = Slider(ax_t2, 'Hopping (t2)', 0.0, 1.0, valinit=init_t2)

def update(val):
    new_M, new_t2 = slider_M.val, slider_t2.val
    for lines, k_vals, N_at, H_M, H_t1, H_t2 in plot_data:
        new_bands = fast_calc(new_M, new_t2, k_vals, N_at, H_M, H_t1, H_t2)
        for n in range(new_bands.shape[1]):
            lines[n].set_ydata(new_bands[:, n])
    fig.canvas.draw_idle()

slider_M.on_changed(update)
slider_t2.on_changed(update)

print("Opening interactive plot...")
plt.show()