import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

def build_honeycomb(Nx, Ny, a=1.0):
    """Generates a finite honeycomb lattice."""
    a1 = np.array([np.sqrt(3) * a, 0])
    a2 = np.array([np.sqrt(3)/2 * a, 3/2 * a])
    delta_A = np.array([0, 0])
    delta_B = np.array([0, a])
    
    x_coords, y_coords, sublattices = [], [], []
    
    for i in range(Nx):
        for j in range(Ny):
            R = i * a1 + j * a2
            
            pos_A = R + delta_A
            x_coords.append(pos_A[0]); y_coords.append(pos_A[1]); sublattices.append('A')
            
            pos_B = R + delta_B
            x_coords.append(pos_B[0]); y_coords.append(pos_B[1]); sublattices.append('B')
            
    return np.array(x_coords), np.array(y_coords), np.array(sublattices)

def build_H_real_space(x, y, sublattices, t1=1.0, t2=1/3, M=0.0, W=0.0):
    """Builds the complete Haldane Hamiltonian with t2 and Disorder in real space."""
    N = len(x)
    H = np.zeros((N, N), dtype=complex)
    coords = np.column_stack((x, y))
    dist_matrix = cdist(coords, coords)#this is a 288 x 288 matrix
    
    # Haldane Vectors for t2 phase
    v_pos_A = [np.array([np.sqrt(3), 0]), np.array([-np.sqrt(3)/2, 1.5]), np.array([-np.sqrt(3)/2, -1.5])]
    v_neg_A = [-v for v in v_pos_A]
    
    # 1. Fill the Hopping Terms
    for i in range(N):
        for j in range(N):
            dist = dist_matrix[i, j]
            
            # Nearest-Neighbor (t1)
            if np.isclose(dist, 1.0, atol=1e-3):
                H[i, j] += -t1
                
            # Next-Nearest-Neighbor (t2)
            elif np.isclose(dist, np.sqrt(3), atol=1e-3):
                dist_vec = coords[i] - coords[j]
                sign = 0
                for v in v_pos_A:
                    if np.allclose(dist_vec, v, atol=1e-3): sign = 1
                for v in v_neg_A:
                    if np.allclose(dist_vec, v, atol=1e-3): sign = -1
                
                if sublattices[i] == 'B': sign *= -1
                    
                if sign != 0:
                    H[i, j] += -t2 * np.exp(1j * sign * np.pi/2)

    # 2. Add Mass (M) and Disorder (W)
    for i in range(N):
        H[i, i] += M if sublattices[i] == 'A' else -M
        if W > 0:
            H[i, i] += np.random.uniform(-W/2, W/2)
            
    return H

def calc_local_chern(H, x, y, E_Fermi=0.0):
    """
    Calculates the Local Chern Marker C(r_i) for every atom in the lattice.
    """
    N = len(x)
    # 1. Solve the Hamiltonian
    evals, evecs = np.linalg.eigh(H)
    
    # 2. Create the Projector Matrix (P) for occupied states (below Fermi Energy)
    occupied = evals < E_Fermi
    # P = sum(|psi><psi|) for all occupied states
    P = evecs[:, occupied] @ evecs[:, occupied].conj().T
    
    # 3. Create Position Matrices
    X = np.diag(x)
    Y = np.diag(y)
    
    # 4. The Local Chern Marker Math
    # C = -2 * pi * i * [PXP, PYP]
    PxP = P @ X @ P
    PyP = P @ Y @ P
    
    Commutator = PxP @ PyP - PyP @ PxP
    C_matrix = -2j * np.pi * Commutator
    
    # The local value for each atom is the diagonal of the resulting matrix
    return np.real(np.diag(C_matrix))

# ==========================================
# Run the Simulation
# ==========================================

# 1. Build a larger flake (12x12 gives a good bulk area)
Nx, Ny = 12, 12
x, y, sub = build_honeycomb(Nx, Ny)

# 2. Build Hamiltonian (Parameters for Topological Phase: M=1.0, W=0)
print(f"Building {len(x)}x{len(x)} Hamiltonian...")
H = build_H_real_space(x, y, sub, t1=1.0, t2=1/3, M=1.0, W=0.0)

# 3. Calculate Local Chern Marker
print("Calculating Local Chern Marker...")
C_local = calc_local_chern(H, x, y)

# 4. Plot the Color Map
plt.figure(figsize=(8, 7))
# Scatter plot where the color 'c' is mapped to the Local Chern value
scatter = plt.scatter(x, y, c=C_local, cmap='coolwarm', s=80, edgecolors='black', linewidth=0.5, vmin=-1, vmax=1)

plt.colorbar(scatter, label='Local Chern Marker $C(r)$')
plt.title(f"Topological Phase (M=1.0, W=0)\nBulk $C(r) \\approx 1$")
plt.xlabel("x position")
plt.ylabel("y position")
plt.axis('equal')
plt.grid(True, linestyle='--', alpha=0.3)
plt.show()