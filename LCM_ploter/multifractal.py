import numpy as np
from scipy.optimize import curve_fit

# ==============================================================================
# Multifractal D_q via TRUE 2D spatial box-counting  (Paper Eqs. 3-6)
# ==============================================================================
def compute_Pq(psi, atom_x, atom_y, box_sizes, q_vals, L):
    """
    Compute the generalised IPR  P_q(λ)  for one eigenstate.

    Returns a dict  { q : array of P_q values, one per box size }
    and the corresponding  log(λ)  array, where λ = l / L.

    This is the *per-realisation* building block; the caller should
    arithmetic-average ⟨P_q⟩ over realisations before extracting τ(q).
    """
    prob = np.abs(psi)**2
    prob /= prob.sum()

    x_min, x_max = atom_x.min(), atom_x.max()
    y_min, y_max = atom_y.min(), atom_y.max()

    log_lambda = []
    Pq_arrays  = {q: [] for q in q_vals}

    for l in box_sizes:
        x_edges = np.arange(x_min, x_max + l, l)
        y_edges = np.arange(y_min, y_max + l, l)

        ix = np.searchsorted(x_edges, atom_x, side='right') - 1
        iy = np.searchsorted(y_edges, atom_y, side='right') - 1

        n_xcells = len(x_edges) - 1
        box_id   = iy * n_xcells + ix

        # Vectorised box summation via np.bincount
        mu = np.bincount(box_id, weights=prob)
        mu = mu[mu > 0]

        log_lambda.append(np.log(l / L))          # λ = l / L  (Eq. 5)

        for q in q_vals:
            if q == 1:
                Pq = -np.sum(mu * np.log(mu))      # Shannon entropy
            else:
                Pq = np.sum(mu**q)
            Pq_arrays[q].append(Pq)

    log_lambda = np.array(log_lambda)
    for q in q_vals:
        Pq_arrays[q] = np.array(Pq_arrays[q])

    return Pq_arrays, log_lambda


def extract_tau_Dq(mean_Pq, log_lambda, q_vals):
    """
    Given ensemble-averaged ⟨P_q(λ)⟩ and log(λ), extract τ(q) and D_q.

    τ(q) = slope of  log⟨P_q⟩  vs  log(λ)        (Eq. 5)
    D_q  = τ(q) / (q − 1)                          (Eq. 6)
    """
    tau, Dq = {}, {}
    for q in q_vals:
        lP    = np.log(np.maximum(mean_Pq[q], 1e-300))
        slope = np.polyfit(log_lambda, lP, 1)[0]
        tau[q] = slope
        Dq[q]  = slope / (q - 1) if q != 1 else slope
    return tau, Dq


# ---------- legacy wrapper (single-shot, kept for quick D2 scans) ----------
def compute_Dq(psi, atom_x, atom_y, box_sizes, q_vals, L=None):
    """
    Single-realisation convenience wrapper.
    If L is None, estimate L from the sample extent.
    """
    if L is None:
        Lx = atom_x.max() - atom_x.min()
        Ly = atom_y.max() - atom_y.min()
        L  = max(Lx, Ly)
    Pq, log_lam = compute_Pq(psi, atom_x, atom_y, box_sizes, q_vals, L)
    return extract_tau_Dq(Pq, log_lam, q_vals)


# ==============================================================================
# Parabolic approximation  τ(q) = d(q−1) + γ q(1−q)   (Eq. 7)
# ==============================================================================
def parabolic_tau(q, d, gamma):
    return d * (q - 1) + gamma * q * (1 - q)

def fit_gamma(q_vals, tau_vals, d=2, q_max=4.0):
    """
    Fit the anomalous exponent γ using only q ≤ q_max,
    where the parabolic approximation is physically valid.
    """
    q = np.array(q_vals)
    t = np.array(tau_vals)
    mask = q <= q_max
    if mask.sum() < 2:
        mask = np.ones_like(q, dtype=bool)  # fallback
    popt, _ = curve_fit(
        lambda qq, g: parabolic_tau(qq, d, g),
        q[mask], t[mask], p0=[0.3], bounds=(0, 3)
    )
    return popt[0]
