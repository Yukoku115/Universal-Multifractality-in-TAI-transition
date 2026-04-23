import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap

# Setup colormap
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

paper_cmap = LinearSegmentedColormap.from_list('paper_d2', [
    (0.00, '#08306b'),   # deep blue      D2 ~ 0
    (0.15, '#084594'),   # extended deep blue
    (0.30, '#2171b5'),   # mid blue
    (0.45, '#6baed6'),   # light blue
    (0.55, '#deebf7'),   # very pale blue D2 ~ 1
    (0.65, '#ffffb2'),   # pale yellow
    (0.80, '#fecc5c'),   # yellow
    (0.92, '#f03b20'),   # orange-red
    (1.00, '#bd0026'),   # dark red       D2 ~ 2
])
mpl.colormaps.register(paper_cmap, name="paper_d2", force=True)

def setup_plot_style():
    plt.rcParams.update({
        'font.size': 12,
        'font.family': 'serif',
        'mathtext.fontset': 'stix',
        'axes.linewidth': 1.2,
        'lines.linewidth': 1.5,
        'figure.dpi': 150
    })

def plot_phase_diagram(W_vals, M_vals, phase_map,
                       Nx_val=None, Ny_val=None, grid_size=None,
                       save_path=None, show_plot=True,
                       vmin=None, vmax=None,
                       cbar_label='Bulk Chern Number',
                       title='Topological Phase Diagram',
                       tR_base=None, sweep_type=None,
                       xlim=None, ylim=None):
    """
    Plot the phase diagram.

    Parameters
    ----------
    vmin, vmax : float or None
        Colour scale limits. If None (default), they are set symmetrically
        around zero using the 99th-percentile absolute value of the data,
        which works well for both the Haldane (0/1) and Kane-Mele (±C_s) cases.
    cbar_label : str
        Label on the colourbar.
    title : str
        Axes title.
    """
    setup_plot_style()
    """
    # Auto-scale: make a symmetric range about 0 so both signs show clearly
    if vmin is None or vmax is None:
        finite = phase_map[np.isfinite(phase_map)]
        if len(finite):
            abs_max = np.percentile(np.abs(finite), 99)
            abs_max = max(abs_max, 0.05)   # floor so the map is never collapsed
        else:
            abs_max = 1.0
        vmin = vmin if vmin is not None else -abs_max
        vmax = vmax if vmax is not None else  abs_max
    """
    fig, ax = plt.subplots(figsize=(6.5, 6))
    mesh = ax.imshow(
        phase_map,
        origin="lower",
        extent=[W_vals.min(), W_vals.max(), M_vals.min(), M_vals.max()],
        cmap="tai_cmap",
        vmin=vmin, vmax=vmax,
        aspect="auto"
    )
    ax.set_box_aspect(1)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_title(title, pad=12)
    ax.set_xlabel("Disorder Strength ($W$)")
    ax.set_ylabel("Staggered Mass ($M$)")
    
    # Add padding so the colorbar and its label don't collide with the title
    cbar = fig.colorbar(mesh, ax=ax, pad=0.04)
    cbar.set_label(cbar_label, labelpad=15)

    if Nx_val and Ny_val and grid_size:
        caption = (f"Honeycomb lattice: $N_x = {Nx_val}$, $N_y = {Ny_val}$ "
                   f"({2 * Nx_val * Ny_val} sites) · Grid: {grid_size}×{grid_size}")
        if tR_base is not None and sweep_type is not None:
             caption += f" · $t_R = {tR_base}$ · Sweep: {sweep_type.capitalize()}"
             
        caption_text = fig.text(
            0.5, 0.02,
            caption,
            ha='center', va='bottom', fontsize=9, color='0.4',
            fontstyle='italic'
        )
    else:
        caption_text = None
    
    if caption_text:
        fig.set_layout_engine(layout='tight', rect=[0, 0.05, 1, 1])
    else:
        fig.set_layout_engine('tight')
    if save_path:
        if caption_text:
            plt.savefig(save_path, dpi=200, bbox_inches='tight', bbox_extra_artists=[caption_text])
        else:
            plt.savefig(save_path, dpi=200, bbox_inches='tight')
        
    if show_plot:
        plt.show()
    else:
        plt.close(fig)

def plot_multifractal_analysis(W_vals, M_vals, D2_map, spectra, phase_points, q_vals, box_a=None, 
                               save_path="multifractal_analysis.png", show_plot=True):
    """
    Plots the large 2x3 D2 phase diagram and 4 tau/Dq spectra subplots.
    """
    from .multifractal import fit_gamma, parabolic_tau
    
    setup_plot_style()
    
    fig = plt.figure(figsize=(15, 8))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.9, 1, 1],
                          hspace=0.45, wspace=0.32,
                          left=0.06, right=0.97, top=0.91, bottom=0.10)
    ax_a = fig.add_subplot(gs[:,0])
    ax_b = fig.add_subplot(gs[0,1])
    ax_c = fig.add_subplot(gs[1,1])
    ax_d = fig.add_subplot(gs[0,2])
    ax_e = fig.add_subplot(gs[1,2])

    # ── (a) D2 map ────────────────────────────────────────────────────────────
    D2_plot = np.clip(D2_map, 0, 2)
    im = ax_a.pcolormesh(W_vals, M_vals, D2_plot, cmap='paper_d2', shading='auto', vmin=0, vmax=2)
    cbar = fig.colorbar(im, ax=ax_a, pad=0.02)
    cbar.set_label(r'$D_2$', fontsize=12)

    for M, W, label, color, marker, panel in phase_points:
        if M_vals.min() <= M <= M_vals.max() and W_vals.min() <= W <= W_vals.max():
            ax_a.plot(W, M, marker=marker, color=color, ms=8, mec='white', mew=0.8, zorder=5)

    # Note: Text tags were hardcoded to the specific W, M range from the plot, dynamically moving or
    # hiding them if they sit outside the bounds prevents overlapping errors.
    if M_vals.max() >= 2.6: ax_a.text(1.5, 2.55, 'BI',  fontsize=12, fontweight='bold', color='white')
    if M_vals.max() >= 2.2 and W_vals.max() >= 8.0: ax_a.text(7.5, 2.10, 'AI',  fontsize=12, fontweight='bold', color='white')
    if M_vals.max() >= 1.5: ax_a.text(1.5, 1.42, 'CI',  fontsize=12, fontweight='bold', color='white')
    if M_vals.max() >= 1.9 and W_vals.max() >= 5.0: ax_a.text(4.5, 1.80, 'TAI', fontsize=12, fontweight='bold', color='white')
    
    ax_a.set_xlabel(r'$W$', fontsize=13)
    ax_a.set_ylabel(r'$M$', fontsize=13)
    ax_a.set_title(r'(a)  $D_2$', fontsize=12, loc='left')
    ax_a.set_aspect('auto') 

    # ── (b,c,d,e) Spectra ─────────────────────────────────────────────────────
    q_arr = np.array(q_vals)
    q_fit = np.linspace(0.1, 10, 200)

    for ax_tau, ax_Dq, panels in [(ax_b, ax_c, 'bc'), (ax_d, ax_e, 'de')]:
        for M, W, label, color, marker, panel in phase_points:
            if panel != panels: continue
            if (M, W) not in spectra: continue
                
            tau_d, Dq_d = spectra[(M,W)]
            tau_arr = np.array([tau_d[q] for q in q_vals])
            Dq_arr  = np.array([Dq_d[q]  for q in q_vals])
            lbl = f"({M},{W})"
            
            ax_tau.plot(q_arr, tau_arr, marker=marker, color=color, ms=5, label=lbl)
            ax_Dq.plot(q_arr, Dq_arr,  marker=marker, color=color, ms=5, label=lbl)
            
            try:
                qf = [q for q in q_vals if q >= 0.5]
                tf = [tau_d[q] for q in qf]
                gamma = fit_gamma(qf, tf)
                ax_tau.plot(q_fit, parabolic_tau(q_fit, 2, gamma), '--', color=color, lw=1.0, alpha=0.6)
            except Exception:
                pass

        ax_tau.axhline(0, color='k', lw=0.5, ls=':')
        ax_tau.set_ylabel(r'$\tau(q)$', fontsize=12)
        ax_tau.legend(fontsize=7.5, framealpha=0.7)
        ax_tau.set_xlim(0, 10)
        ax_tau.set_xticks([0, 2, 4, 6, 8, 10])

        ax_Dq.axhline(0, color='k', lw=0.5, ls=':')
        ax_Dq.axhline(2, color='k', lw=0.5, ls=':')
        ax_Dq.set_xlabel(r'$q$', fontsize=12)
        ax_Dq.set_ylabel(r'$D_q$', fontsize=12)
        ax_Dq.set_ylim(0, 2.5)
        ax_Dq.set_xlim(0, 10)
        ax_Dq.set_xticks([0, 2, 4, 6, 8, 10])

    ax_b.set_title(r'(b)', fontsize=12, loc='left')
    ax_c.set_title(r'(c)', fontsize=12, loc='left')
    ax_d.set_title(r'(d)', fontsize=12, loc='left')
    ax_e.set_title(r'(e)', fontsize=12, loc='left')

    fig.suptitle(r"Multifractal analysis (Anderson disorder, $t_2=\frac{1}{3}$)", fontsize=13)

    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
        
    if show_plot:
        plt.show()
    else:
        plt.close(fig)
