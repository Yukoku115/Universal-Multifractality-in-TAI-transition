import numpy as np
from scipy.spatial.distance import cdist

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
# Geometry (used by both Haldane and Kane-Mele)
# ==============================================================================
def precompute_geometry(x, y, sublattices):
    """
    Precomputes lattice geometry matrices needed for Hamiltonian construction.

    Returns
    -------
    nn_mask      : (N,N) bool   — nearest-neighbour pair mask
    active_nnn   : (N,N) bool   — next-nearest-neighbour pair mask (with correct sign)
    sign_matrix  : (N,N) float  — +/-1 Haldane phase signs for NNN hoppings
    dx_hat       : (N,N) float  — x-component of unit vectors along NN bonds
    dy_hat       : (N,N) float  — y-component of unit vectors along NN bonds
    """
    coords      = np.column_stack((x, y))
    dist_matrix = cdist(coords, coords)
    diff_vecs   = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]

    nn_mask  = np.isclose(dist_matrix, 1.0,        atol=1e-3)
    nnn_mask = np.isclose(dist_matrix, np.sqrt(3), atol=1e-3)

    # Haldane NNN sign matrix
    sign_matrix = np.zeros((len(x), len(x)), dtype=float)
    for v in _NNN_VECS_POS:
        sign_matrix[np.all(np.isclose(diff_vecs, v, atol=1e-3), axis=2)] =  1.0
    for v in _NNN_VECS_NEG:
        sign_matrix[np.all(np.isclose(diff_vecs, v, atol=1e-3), axis=2)] = -1.0
    sign_matrix[sublattices == 'B', :] *= -1

    active_nnn = nnn_mask & (sign_matrix != 0)

    # NN bond unit vectors (needed for Rashba coupling)
    dist_safe = np.where(nn_mask, dist_matrix, 1.0)   # avoid div-by-zero
    dx_hat = diff_vecs[:, :, 0] / dist_safe
    dy_hat = diff_vecs[:, :, 1] / dist_safe

    return nn_mask, active_nnn, sign_matrix, dx_hat, dy_hat

# ==============================================================================
# Haldane (spinless) Hamiltonian  —  unchanged physics, updated signature
# ==============================================================================
def build_H_from_geometry(x, y, sublattices, nn_mask, active_nnn, sign_matrix,
                           t1=1.0, t2=1/3, M=0.0, W=0.0, disorder_type="Anderson",
                           dx_hat=None, dy_hat=None):
    """Builds the spinless N x N Haldane Hamiltonian."""
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
# Kane-Mele (spinful) Hamiltonian  —  2N x 2N
# ==============================================================================
def build_H_KaneMele(x, y, sublattices, nn_mask, active_nnn, sign_matrix,
                     dx_hat, dy_hat,
                     t1=1.0, t2=1/3, M=0.0,
                     tR_base=0.0, W_anderson=0.0, W_rashba=0.0):
    """
    Builds the 2N x 2N Kane-Mele Hamiltonian.

    The basis order is [spin-up sites 0..N-1 | spin-down sites 0..N-1].

    Parameters
    ----------
    tR_base   : float — uniform Rashba coupling strength
    W_anderson: float — Anderson on-site disorder amplitude
    W_rashba  : float — amplitude of random bond-disorder on the Rashba coupling
    """
    N = len(x)

    # ── Diagonal (on-site) terms ──────────────────────────────────────────────
    diag = np.where(sublattices == 'A', M, -M)
    if W_anderson > 0:
        diag = diag + np.random.uniform(-W_anderson / 2, W_anderson / 2, N)

    # ── Spin-conserving blocks (N x N each) ───────────────────────────────────
    H_uu = np.zeros((N, N), dtype=complex)
    H_dd = np.zeros((N, N), dtype=complex)

    # NN hopping (same for both spins)
    H_uu[nn_mask] = -t1
    H_dd[nn_mask] = -t1

    # Intrinsic SOC: matches Haldane convention (original code: -i*t2*nu_ij)
    # H_uu must equal H_Haldane in the NNN sector so that C_up = +1 in the
    # topological phase, and H_dd = H_uu* (time-reversal) gives C_down = -1.
    H_uu[active_nnn] = -1j * t2 * sign_matrix[active_nnn]   # same as Haldane
    H_dd[active_nnn] = +1j * t2 * sign_matrix[active_nnn]   # TR partner

    # On-site energies
    H_uu[np.arange(N), np.arange(N)] = diag
    H_dd[np.arange(N), np.arange(N)] = diag

    # ── Spin-flipping block: Rashba SOC ───────────────────────────────────────
    H_ud = np.zeros((N, N), dtype=complex)   # up <- down
    H_du = np.zeros((N, N), dtype=complex)   # down <- up

    # Bond-disordered Rashba coupling strength (symmetric: t_R^{ij} = t_R^{ji})
    if W_rashba > 0:
        tR_rand = np.random.uniform(-W_rashba / 2, W_rashba / 2, (N, N))
        tR_rand = (tR_rand + tR_rand.T) / 2
        tR_mat  = tR_base + tR_rand
    else:
        tR_mat = np.full((N, N), tR_base)

    # Matrix elements: i*t_R*(dx_hat ± i*dy_hat)
    H_ud[nn_mask] = 1j * tR_mat[nn_mask] * (dx_hat[nn_mask] - 1j * dy_hat[nn_mask])
    H_du[nn_mask] = 1j * tR_mat[nn_mask] * (dx_hat[nn_mask] + 1j * dy_hat[nn_mask])

    # ── Assemble full 2N x 2N matrix ─────────────────────────────────────────
    H_KM = np.block([
        [H_uu, H_ud],
        [H_du, H_dd]
    ])

    return H_KM

# ==============================================================================
# Local Chern Marker — spinless
# ==============================================================================
def calc_local_chern(H, x, y, E_Fermi=0.0):
    """Returns the local Chern marker for each site (spinless N x N system)."""
    evals, evecs = np.linalg.eigh(H)
    occ = evecs[:, evals < E_Fermi]
    P   = occ @ occ.conj().T

    PxP = (P * x[np.newaxis, :]) @ P
    PyP = (P * y[np.newaxis, :]) @ P

    return np.real(np.diag(-2j * np.pi * (PxP @ PyP - PyP @ PxP)))

# ==============================================================================
# Local Spin Chern Marker — spinful (Kane-Mele)
# ==============================================================================
def calc_local_spin_chern(H_2N, x, y, E_Fermi=0.0):
    """
    Returns the local Spin Chern Marker C_s = (C_up - C_down) / 2 per site.

    The total Chern number C_up + C_down is exactly zero in a time-reversal
    symmetric system. The spin Chern number C_s = (C_up - C_down) / 2
    is the correct topological invariant for the Quantum Spin Hall phase.
    """
    x_2N = np.concatenate([x, x])
    y_2N = np.concatenate([y, y])

    evals, evecs = np.linalg.eigh(H_2N)
    occ = evecs[:, evals < E_Fermi]
    P   = occ @ occ.conj().T

    PxP = (P * x_2N[np.newaxis, :]) @ P
    PyP = (P * y_2N[np.newaxis, :]) @ P

    C_2N = np.real(np.diag(-2j * np.pi * (PxP @ PyP - PyP @ PxP)))

    N = len(x)
    # C_up in the first N entries, C_down in the last N entries
    return (C_2N[:N] - C_2N[N:]) / 2.0

# ==============================================================================
# Per-row workers for parallelisation
# ==============================================================================
def compute_row(i, M_val, W_vals, done_row, x, y, sub,
                nn_mask, active_nnn, sign_matrix, bulk_mask,
                dx_hat=None, dy_hat=None):
    """Spinless Haldane row worker."""
    row = np.full(len(W_vals), np.nan)
    for j, W_val in enumerate(W_vals):
        if done_row[j]:
            continue
        H       = build_H_from_geometry(x, y, sub, nn_mask, active_nnn, sign_matrix,
                                        M=M_val, W=W_val, disorder_type="Anderson")
        C_local = calc_local_chern(H, x, y)
        row[j]  = float(np.mean(C_local[bulk_mask]))
    return i, row

def compute_row_km(i, M_val, W_vals, done_row, x, y, sub,
                   nn_mask, active_nnn, sign_matrix, dx_hat, dy_hat,
                   bulk_mask, tR_base=0.3, W_anderson_fixed=0.0, W_rashba_fixed=0.0,
                   sweep_type="anderson"):
    """
    Kane-Mele row worker.

    Sweeps over W_vals for a fixed M_val. The parameter swept is determined
    by `sweep_type`:
      - 'anderson': W_vals represents on-site Anderson disorder.
                    W_rashba is held constant at W_rashba_fixed.
      - 'rashba'  : W_vals represents bond Rashba disorder (spin disorder).
                    W_anderson is held constant at W_anderson_fixed.
    """
    row = np.full(len(W_vals), np.nan)
    for j, W_val in enumerate(W_vals):
        if done_row[j]:
            continue

        if sweep_type == "anderson":
            w_and = W_val
            w_rash = W_rashba_fixed
        elif sweep_type == "rashba":
            w_and = W_anderson_fixed
            w_rash = W_val
        else:
            raise ValueError("sweep_type must be 'anderson' or 'rashba'")

        H_KM    = build_H_KaneMele(
                      x, y, sub, nn_mask, active_nnn, sign_matrix, dx_hat, dy_hat,
                      M=M_val, tR_base=tR_base, W_anderson=w_and, W_rashba=w_rash)
        C_spin  = calc_local_spin_chern(H_KM, x, y)
        row[j]  = float(np.mean(C_spin[bulk_mask]))
    return i, row

