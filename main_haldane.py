import os
import sys
import time
import numpy as np
from collections import deque
from joblib import Parallel, delayed
from scipy.spatial.distance import cdist

from LCM_ploter import (
    build_honeycomb,
    precompute_geometry,
    compute_row,
    plot_phase_diagram
)

def main():
    # ==========================================================================
    # Parameters
    # ==========================================================================
    Nx_val, Ny_val = 20, 20
    grid_size      = 50
    W_vals         = np.linspace(0.0, 9.5, grid_size)
    M_vals         = np.linspace(1.3, 2.8, grid_size)
    CHECKPOINT     = "phase_diagram_checkpoint.npz"
    N_JOBS         = 4   # -1 = all available CPU cores; set to e.g. 4 to limit

    # ==========================================================================
    # Build lattice, bulk mask, and geometry  —  all done once
    # ==========================================================================
    x_lat, y_lat, sub_lat = build_honeycomb(Nx=Nx_val, Ny=Ny_val)

    dist_from_center = np.sqrt((x_lat - x_lat.mean())**2 + (y_lat - y_lat.mean())**2)
    bulk_mask        = dist_from_center < np.max(dist_from_center) * 0.40

    print("Precomputing lattice geometry...", end=" ", flush=True)
    nn_mask, active_nnn, sign_matrix, dx_hat, dy_hat = precompute_geometry(x_lat, y_lat, sub_lat)
    print("done.")

    # ==========================================================================
    # Checkpoint: load existing progress or start fresh
    # ==========================================================================
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

    # ==========================================================================
    # Main sweep
    # ==========================================================================
    pending_rows  = [i for i in range(len(M_vals)) if not done_mask[i].all()]
    total_points  = grid_size ** 2

    BATCH_SIZE    = N_JOBS if N_JOBS > 0 else (os.cpu_count() or 4)
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

        # ETA tracking
        timing_window.append((batch_points, time.time() - t_batch_start))
        completed  = int(done_mask.sum())
        remaining  = total_points - completed
        total_t    = sum(t for _, t in timing_window)
        total_pts  = sum(p for p, _ in timing_window)

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

    # ==========================================================================
    # Plotting using LCM_ploter
    # ==========================================================================
    """
    plot_phase_diagram(
        W_vals, M_vals, phase_map, 
        Nx_val, Ny_val, grid_size, 
        save_path="phase_diagram_out.png"
    )
    """
    # ==========================================================================
    # Multifractal Analysis (D2 map & Spectra)
    # ==========================================================================
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import eigsh
    from LCM_ploter import build_H_from_geometry, compute_Pq, extract_tau_Dq, plot_multifractal_analysis

    print("\n\n=== Computing Multifractal D2 Map & Spectra ===")
    MF_CHECKPOINT = "multifractal_checkpoint.npz"

    # Lattices
    Nx_a, Ny_a     = 20, 20      # size for D2 phase diagram
    Nx_s, Ny_s     = 20, 20      # size for spectra points
    grid_mf        = 20          # resolution
    n_real_a       = 20          # realisations for D2 phase map
    n_spectra      = 50          # realisations for spectra 

    W_mf = np.linspace(0.0, 9.5, grid_mf)
    M_mf = np.linspace(1.3, 2.8, grid_mf)
    
    x_a, y_a, sub_a = build_honeycomb(Nx_a, Ny_a)
    x_s, y_s, sub_s = build_honeycomb(Nx_s, Ny_s)
    
    # System size L = max extent of the sample  (needed for λ = l/L)
    L_a = max(x_a.max() - x_a.min(), y_a.max() - y_a.min())
    L_s = max(x_s.max() - x_s.min(), y_s.max() - y_s.min())
    
    # Box sizes: from ~2a up to ~L/2  (paper uses l/L up to 0.5)
    # We specify a range and dynamically filter per-lattice below
    box_all = [2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30]
    box_a   = [l for l in box_all if l <= L_a / 2]
    box_s   = [l for l in box_all if l <= L_s / 2]
    q_d2   = [2]
    q_vals = [0.5, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10]
    
    # Paper: "we target the four eigenstates closest to E = 0" (k=4)
    K_EIGS = 4

    phase_points = [
        (1.44, 1.50, "CI",       "#e6550d", "h",  "bc"),
        (2.00, 5.50, "AI",       "#31a354", "D",  "bc"),
        (1.44, 5.50, "TAI",      "#e31a1c", "o",  "bc"),
        (2.60, 4.00, "BI",       "#984ea3", "^",  "bc"),
        (2.18, 4.39, "Boundary", "#e7298a", "s",  "de"),
        (1.90, 6.82, "Boundary", "#7570b3", "p",  "de"),
        (2.42, 6.09, "Boundary", "#1f78b4", "D",  "de"),
    ]

    # Load checkpoint if exists
    if os.path.exists(MF_CHECKPOINT):
        mf_data = np.load(MF_CHECKPOINT, allow_pickle=True)
        D2_map      = mf_data["D2_map"]
        done_D2     = mf_data["done_D2"].astype(bool)
        spectra     = mf_data["spectra"].item()
        done_points = mf_data["done_points"].item()
        print(f"Resuming Multifractal checkpoint — {done_D2.sum()}/{grid_mf**2} D2 points done.")
    else:
        print("No multifractal checkpoint found — starting fresh.")
        D2_map      = np.zeros((grid_mf, grid_mf))
        done_D2     = np.zeros((grid_mf, grid_mf), dtype=bool)
        spectra     = {}
        done_points = { (M,W): False for (M,W,_,_,_,_) in phase_points }

    # Precompute geometries
    _, active_a, sign_a, _, _ = precompute_geometry(x_a, y_a, sub_a)
    nn_a = np.isclose(cdist(np.column_stack((x_a, y_a)), np.column_stack((x_a, y_a))), 1.0, atol=1e-3)
    
    _, active_s, sign_s, _, _ = precompute_geometry(x_s, y_s, sub_s)
    nn_s = np.isclose(cdist(np.column_stack((x_s, y_s)), np.column_stack((x_s, y_s))), 1.0, atol=1e-3)

    # ──────────────────────────────────────────────────────────────────────────
    # 1. Compute D2 map   (ensemble-averaged ⟨P_q⟩ per grid point)
    # ──────────────────────────────────────────────────────────────────────────
    total = grid_mf**2
    cnt = int(done_D2.sum())
    rng = np.random.default_rng(42)
    t_start = time.time() - (cnt * 0.1)
    
    for i, M in enumerate(M_mf):
        for j, W in enumerate(W_mf):
            if done_D2[i, j]: continue
                
            cnt += 1
            if cnt % 5 == 1:
                elapsed = time.time() - t_start
                avg = elapsed / cnt
                rem = avg * (total - cnt)
                m, s = divmod(int(rem), 60)
                sys.stdout.write(f"\r  D2 map {cnt}/{total} (M={M:.2f}, W={W:.2f})  ETA: {m}m {s:02d}s" + " "*10)
                sys.stdout.flush()

            # Accumulate Pq arrays across all realisations & eigenstates
            Pq_accum = {q: [] for q in q_d2}
            log_lam  = None
            
            for _ in range(n_real_a):
                H = build_H_from_geometry(x_a, y_a, sub_a, nn_a, active_a, sign_a, M=M, W=W)
                try:
                    evals, evecs = eigsh(csr_matrix(H), k=K_EIGS, sigma=0.0, which='LM', tol=1e-3)
                except:
                    evals_full, evecs_full = np.linalg.eigh(H)
                    idx = np.argsort(np.abs(evals_full))[:K_EIGS]
                    evals, evecs = evals_full[idx], evecs_full[:, idx]
                
                for k in range(evecs.shape[1]):
                    Pq_k, log_lam = compute_Pq(evecs[:, k], x_a, y_a, box_a, q_d2, L_a)
                    for q in q_d2:
                        Pq_accum[q].append(Pq_k[q])

            # Arithmetic mean of Pq across realisations, then extract tau
            if Pq_accum[2] and log_lam is not None:
                mean_Pq = {q: np.mean(Pq_accum[q], axis=0) for q in q_d2}
                tau_d2, Dq_d2 = extract_tau_Dq(mean_Pq, log_lam, q_d2)
                D2_map[i, j] = Dq_d2[2]
            else:
                D2_map[i, j] = 0.0
                
            done_D2[i, j] = True
            
            if cnt % 20 == 0:
                np.savez(MF_CHECKPOINT, D2_map=D2_map, done_D2=done_D2, spectra=spectra, done_points=done_points)

    print("\nD2 map done. Computing Spectra...")

    # ──────────────────────────────────────────────────────────────────────────
    # 2. Compute Spectra   (ensemble-averaged ⟨P_q⟩ per phase point)
    # ──────────────────────────────────────────────────────────────────────────
    rng2 = np.random.default_rng(0)
    for M, W, label, _, _, _ in phase_points:
        if done_points[(M,W)]: continue
        
        print(f"  Computing point ({M:.2f},{W:.2f}) {label} ...", end=' ', flush=True)
        
        Pq_accum = {q: [] for q in q_vals}
        log_lam  = None
        
        for _ in range(n_spectra):
            H = build_H_from_geometry(x_s, y_s, sub_s, nn_s, active_s, sign_s, M=M, W=W)
            try:
                evals, evecs = eigsh(csr_matrix(H), k=K_EIGS, sigma=0.0, which='LM', tol=1e-3)
            except:
                evals_full, evecs_full = np.linalg.eigh(H)
                idx = np.argsort(np.abs(evals_full))[:K_EIGS]
                evals, evecs = evals_full[idx], evecs_full[:, idx]
                
            for k in range(evecs.shape[1]):
                Pq_k, log_lam = compute_Pq(evecs[:, k], x_s, y_s, box_s, q_vals, L_s)
                for q in q_vals:
                    Pq_accum[q].append(Pq_k[q])

        # Arithmetic mean ⟨P_q⟩, then extract tau and Dq
        if log_lam is not None and Pq_accum[q_vals[0]]:
            mean_Pq = {q: np.mean(Pq_accum[q], axis=0) for q in q_vals}
            tau_m, Dq_m = extract_tau_Dq(mean_Pq, log_lam, q_vals)
        else:
            tau_m = {q: 0.0 for q in q_vals}
            Dq_m  = {q: 0.0 for q in q_vals}
            
        spectra[(M, W)] = (tau_m, Dq_m)
        done_points[(M,W)] = True
        
        np.savez(MF_CHECKPOINT, D2_map=D2_map, done_D2=done_D2, spectra=spectra, done_points=done_points)
        print("done")
        
    print("\nSpectra complete. Generating multifractal plots...")

    plot_multifractal_analysis(
        W_mf, M_mf, D2_map, spectra, phase_points, q_vals, box_a,
        save_path="multifractal_analysis.png"
    )

if __name__ == "__main__":
    main()
