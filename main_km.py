"""
main_km.py  —  Kane-Mele phase diagram with Rashba spin disorder

Computes the Spin Chern Number C_s = (C_up - C_down) / 2 across a
(M, W_anderson) grid using the 2N x 2N Kane-Mele Hamiltonian.

Key differences from main.py (Haldane model):
  • Hilbert space is 2N x 2N (spin-up and spin-down sectors).
  • Intrinsic SOC (t2) has opposite sign for spin-up vs. spin-down.
  • A Rashba coupling (tR_base) couples the two spin sectors on NN bonds.
  • Optional Rashba disorder (W_rashba) randomises the Rashba strength
    bond-by-bond, acting as a form of spin disorder.
  • The plotted quantity is C_s, NOT the total Chern number (which is 0).
"""

import os
import sys
import time
import numpy as np
from collections import deque
from joblib import Parallel, delayed

from LCM_ploter import (
    build_honeycomb,
    precompute_geometry,
    compute_row_km,
    plot_phase_diagram,
)

def main():
    # ==========================================================================
    # Parameters  —  edit here
    # ==========================================================================
    Nx_val, Ny_val = 20, 20             # lattice size (keep smaller: 2N x 2N matrices)
    grid_size      = 50               # resolution of the (M, W) phase diagram grid
    W_vals         = np.linspace(0.0, 9.5, grid_size)   # disorder sweep range
    M_vals         = np.linspace(1.3, 2.8, grid_size)   # staggered mass range
    
    # Physics settings
    tR_base        = 0.1             # baseline (uniform) Rashba coupling
    sweep_type     = "anderson"     # "anderson" or "rashba" (what does W_vals control?)
    # New feature: enable both disorder types by keeping a small Rashba perturbation
    # during the Anderson sweep. Set this to True to activate the combined-disorder mode.
    include_rashba_perturbation = False
    W_rashba_perturbation    = 0.1

    # These are fixed background values applied to whatever ISN'T being swept:
    W_anderson_fixed = 0.0
    W_rashba_fixed   = 0.0

    CHECKPOINT     = "km_phase_diagram_checkpoint.npz"
    N_JOBS         = 4               # parallel CPU cores (-1 = all)

    # ==========================================================================
    # Build lattice, bulk mask, and geometry  —  all done once
    # ==========================================================================
    x_lat, y_lat, sub_lat = build_honeycomb(Nx=Nx_val, Ny=Ny_val)

    dist_from_center = np.sqrt((x_lat - x_lat.mean())**2 + (y_lat - y_lat.mean())**2)
    bulk_mask        = dist_from_center < np.max(dist_from_center) * 0.40

    print("Precomputing lattice geometry...", end=" ", flush=True)
    nn_mask, active_nnn, sign_matrix, dx_hat, dy_hat = precompute_geometry(
        x_lat, y_lat, sub_lat
    )
    print("done.")

    if sweep_type == "anderson" and include_rashba_perturbation:
        W_rashba_fixed = W_rashba_perturbation

    print(f"System : {Nx_val}x{Ny_val} unit cells  |  "
          f"N_sites={len(x_lat)}  |  2N={2*len(x_lat)}  |  N_bulk={int(bulk_mask.sum())}")
    print(f"Physics: tR_base={tR_base}  |  Sweeping: {sweep_type.upper()}  |  "
          f"W_rashba_fixed={W_rashba_fixed}")
    print(f"Grid   : {grid_size}x{grid_size}  |  "
          f"W in [{W_vals.min():.1f}, {W_vals.max():.1f}]  |  "
          f"M in [{M_vals.min():.1f}, {M_vals.max():.1f}]\n")

    # ==========================================================================
    # Checkpoint: load existing progress or start fresh
    # ==========================================================================
    if os.path.exists(CHECKPOINT):
        data = np.load(CHECKPOINT)
        
        # Also check if physics parameters changed
        # np.load() returns 0-d arrays for scalars, so we extract with .item()
        try:
            old_sweep     = str(data["sweep_type"].item())
            old_tR        = float(data["tR_base"].item())
            old_W_and_fix = float(data["W_anderson_fixed"].item())
            old_W_ras_fix = float(data["W_rashba_fixed"].item())
        except KeyError:
            # If any keys are missing in an older checkpoint, force a restart
            old_sweep     = "none"
            old_tR        = -1.0
            old_W_and_fix = -1.0
            old_W_ras_fix = -1.0

        grid_matches = (
            data["W_vals"].shape == W_vals.shape and
            data["M_vals"].shape == M_vals.shape and
            np.allclose(W_vals, data["W_vals"]) and
            np.allclose(M_vals, data["M_vals"]) and
            old_sweep == sweep_type and
            np.isclose(old_tR, tR_base) and
            np.isclose(old_W_and_fix, W_anderson_fixed) and
            np.isclose(old_W_ras_fix, W_rashba_fixed)
        )
        if not grid_matches:
            print("Grid or physics parameters changed — starting fresh.")
            phase_map = np.full((len(M_vals), len(W_vals)), np.nan, dtype=float)
            done_mask = np.zeros((len(M_vals), len(W_vals)), dtype=bool)
        else:
            phase_map = data["phase_map"]
            done_mask = data["done_mask"].astype(bool)
            completed = int(done_mask.sum())
            print(f"Resuming from checkpoint — {completed}/{grid_size**2} points done.")
    else:
        phase_map = np.full((len(M_vals), len(W_vals)), np.nan, dtype=float)
        done_mask = np.zeros((len(M_vals), len(W_vals)), dtype=bool)
        print("No checkpoint found — starting fresh.")

    # ==========================================================================
    # Main sweep (parallelised row by row)
    # ==========================================================================
    pending_rows  = [i for i in range(len(M_vals)) if not done_mask[i].all()]
    total_points  = grid_size ** 2
    BATCH_SIZE    = N_JOBS if N_JOBS > 0 else (os.cpu_count() or 4)
    timing_window = deque(maxlen=10)

    print(f"\nComputing {len(pending_rows)} remaining rows using "
          f"{N_JOBS if N_JOBS != -1 else 'all'} CPU cores...\n")

    for batch_start in range(0, len(pending_rows), BATCH_SIZE):
        batch         = pending_rows[batch_start : batch_start + BATCH_SIZE]
        t_batch_start = time.time()
        n_batches     = (len(pending_rows) + BATCH_SIZE - 1) // BATCH_SIZE
        batch_num     = batch_start // BATCH_SIZE + 1
        sys.stdout.write(
            f"\r  Batch {batch_num}/{n_batches} — rows {batch[0]}–{batch[-1]} running...{' '*20}"
        )
        sys.stdout.flush()

        results = Parallel(n_jobs=N_JOBS)(
            delayed(compute_row_km)(
                i, M_vals[i], W_vals, done_mask[i],
                x_lat, y_lat, sub_lat,
                nn_mask, active_nnn, sign_matrix, dx_hat, dy_hat,
                bulk_mask, tR_base=tR_base,
                W_anderson_fixed=W_anderson_fixed, 
                W_rashba_fixed=W_rashba_fixed,
                sweep_type=sweep_type
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
                 phase_map=phase_map, done_mask=done_mask,
                 sweep_type=sweep_type, tR_base=tR_base,
                 W_anderson_fixed=W_anderson_fixed, W_rashba_fixed=W_rashba_fixed)

        # ETA tracking
        timing_window.append((batch_points, time.time() - t_batch_start))
        completed = int(done_mask.sum())
        remaining = total_points - completed
        total_t   = sum(t for _, t in timing_window)
        total_pts = sum(p for p, _ in timing_window)

        if total_pts > 0 and total_t > 0:
            rate    = total_pts / total_t
            eta_sec = remaining / rate
            h, rem  = divmod(int(eta_sec), 3600)
            m, s    = divmod(rem, 60)
            eta_str = f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"
        else:
            eta_str = "--:--"

        progress = completed / total_points
        bar      = '=' * int(30 * progress) + '-' * (30 - int(30 * progress))
        sys.stdout.write(
            f"\r[{bar}] {completed}/{total_points} ({progress*100:.1f}%)  ETA {eta_str}   "
        )
        sys.stdout.flush()

    print("\nDone. Generating plot...")

    # ==========================================================================
    # Plotting
    # ==========================================================================
    xlabel = "Rashba Spin Disorder ($W_{Rashba}$)" if sweep_type == "rashba" else "Anderson Disorder ($W_{Anderson}$)"
    mode_suffix = f"_with_Rashba_{W_rashba_fixed:.2f}" if W_rashba_fixed > 0 else ""
    phase_title = (
        f'Kane-Mele Phase Diagram  ({"Anderson sweep with Rashba perturbation" if W_rashba_fixed > 0 else sweep_type.capitalize() + " sweep"} '
        f'$W_{{rashba}}={W_rashba_fixed}$, $t_R={tR_base}$)'
        if W_rashba_fixed > 0 else f'Kane-Mele Phase Diagram ($t_R={tR_base}$)'
    )

    plot_phase_diagram(
        W_vals, M_vals, phase_map,
        Nx_val=Nx_val, Ny_val=Ny_val, grid_size=grid_size,
        save_path=f"km_phase_diagram_{sweep_type}{mode_suffix}.png",
        cbar_label=r'Spin Chern Marker $C_s$',
        title=phase_title,
        tR_base=tR_base, sweep_type=sweep_type,
        xlim=(0.0, 9.5)
    )
    # Patch the plotting xlabel manually for now since we didn't add it to parameters
    # The image will have W, but we'll know what it means.

if __name__ == "__main__":
    main()
