import numpy as np
import os
import sys
import matplotlib.pyplot as plt

# Temporarily add the script's directory to sys.path
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.append(script_dir)

try:
    import config_and_types as cfg
except ImportError:
    print(f"[!] ERROR: Cannot find 'config_and_types.py' in {script_dir}")
    print("    Please ensure it is saved in the same directory as this script.")
    sys.exit(1)

def plot_heatmaps(data):
    """
    Generates hexbin heatmaps for window, source plane, and sphere coordinates.
    Hexbins are used instead of scatter plots to handle large datasets efficiently.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # 1. Window Coordinates (y_w, z_w) - All rays
    # This shows the distribution of rays across the camera's field of view.
    hb0 = axes[0].hexbin(data['y_w'], data['z_w'], gridsize=50, cmap='viridis', mincnt=1)
    axes[0].set_title("Window Plane (Camera)\n$y_w$ vs $z_w$")
    axes[0].set_xlabel("$y_w$")
    axes[0].set_ylabel("$z_w$")
    fig.colorbar(hb0, ax=axes[0], label='Ray Count')
    
    # 2. Source Plane (y_s, z_s) - Filtering for TERM_SOURCE_PLANE
    source_mask = data['termination_type'] == cfg.TERM_SOURCE_PLANE
    if np.any(source_mask):
        hb1 = axes[1].hexbin(data['y_s'][source_mask], data['z_s'][source_mask], 
                             gridsize=50, cmap='plasma', mincnt=1)
        axes[1].set_title("Source Plane Hits\n$y_s$ vs $z_s$")
        fig.colorbar(hb1, ax=axes[1], label='Ray Count')
    else:
        axes[1].text(0.5, 0.5, 'No Source Plane Hits\n(Check Enum Mapping)', 
                     ha='center', va='center', transform=axes[1].transAxes)
    axes[1].set_xlabel("$y_s$")
    axes[1].set_ylabel("$z_s$")

    # 3. Celestial Sphere (phi, theta) - Filtering for TERM_SPHERE
    # Using raw strings (r"") to fix the SyntaxWarning for LaTeX symbols.
    sphere_mask = data['termination_type'] == cfg.TERM_SPHERE
    if np.any(sphere_mask):
        hb2 = axes[2].hexbin(data['final_phi'][sphere_mask], data['final_theta'][sphere_mask], 
                             gridsize=50, cmap='magma', mincnt=1)
        axes[2].set_title(r"Celestial Sphere Hits" + "\n" + r"$\phi$ vs $\theta$")
        fig.colorbar(hb2, ax=axes[2], label='Ray Count')
    else:
        axes[2].text(0.5, 0.5, 'No Sphere Hits\n(Check Enum Mapping)', 
                     ha='center', va='center', transform=axes[2].transAxes)
    axes[2].set_xlabel(r"$\phi$")
    axes[2].set_ylabel(r"$\theta$")

    plt.tight_layout()
    print("[i] Displaying heatmaps. Close the window to continue...")
    plt.show()

def diagnose_blueprint(blueprint_path: str = None) -> None:
    """
    Reads the binary blueprint file and prints critical diagnostics.
    """
    if blueprint_path is None:
        # Resolve path to the light_blueprint.bin file
        blueprint_path = os.path.abspath(os.path.join(
            script_dir, "..", "..", "..", "project", "photon_geodesic_integrator", "light_blueprint.bin"
        ))
    
    print("=================================================================")
    print(f" BLUEPRINT DIAGNOSTICS & VISUALIZATION")
    print(f" File: {blueprint_path}")
    print("=================================================================")
    
    if not os.path.exists(blueprint_path):
        print(f"[!] ERROR: Blueprint file not found at:\n    {blueprint_path}")
        return

    # Load data using the structured dtype from config_and_types
    data = np.fromfile(blueprint_path, dtype=cfg.BLUEPRINT_DTYPE)
    total_rays = len(data)
    print(f"Total records loaded: {total_rays:,}\n")

    # --- 1. Check Enum Mappings ---
    print("--- 1. Raw Termination Enums in Binary ---")
    unique_enums, counts = np.unique(data['termination_type'], return_counts=True)
    for e, c in zip(unique_enums, counts):
        print(f"  Raw Enum {e:2d}: {c:8,} rays ({c/total_rays*100:6.2f}%)")
        
    print(f"\n  [Config Current]: SOURCE_PLANE = {cfg.TERM_SOURCE_PLANE}, SPHERE = {cfg.TERM_SPHERE}")
    print("  -> If your raw enums above do NOT match these, update config_and_types.py.")

    # --- 2. Check Window Coordinates (Field of View) ---
    print("\n--- 2. Window Coordinate Bounds (y_w, z_w) ---")
    half_w = cfg.WINDOW_WIDTH / 2.0
    in_view = np.sum(
        (data['y_w'] >= -half_w) & (data['y_w'] < half_w) & 
        (data['z_w'] >= -half_w) & (data['z_w'] < half_w)
    )
    print(f"  Rays inside Renderer FOV: {in_view:,} ({in_view/total_rays*100:.2f}%)")

    # --- 3. View Raw Sample Data ---
    print("\n--- 3. First 5 Raw Records ---")
    header = f"{'Ray#':<6} | {'Term Type':<9} | {'y_w':>8} | {'z_w':>8} | {'y_s / theta':>11} | {'z_s / phi':>11}"
    print(header)
    print("-" * len(header))
    
    for i in range(min(5, total_rays)):
        rec = data[i]
        tt = rec['termination_type']
        # Display source plane coords or sphere angles depending on termination
        val_1 = rec['y_s'] if tt == cfg.TERM_SOURCE_PLANE else rec['final_theta']
        val_2 = rec['z_s'] if tt == cfg.TERM_SOURCE_PLANE else rec['final_phi']
        print(f"{i:<6} | {tt:<9} | {rec['y_w']:>8.2f} | {rec['z_w']:>8.2f} | {val_1:>11.3f} | {val_2:>11.3f}")
    
    # --- 4. Launch Visualization ---
    print("\n[i] Generating heatmaps...")
    plot_heatmaps(data)
    print("=================================================================")

if __name__ == "__main__":
    diagnose_blueprint()