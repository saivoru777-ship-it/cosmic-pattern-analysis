"""
Create 3D visualizations of discovered cosmic patterns.

Shows:
1. Galaxy distribution in 3D
2. Voronoi cells colored by neighbor count
3. Density field slices
4. Pattern clusters in 3D
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree


def visualize_3d_distribution(positions, title="3D Galaxy Distribution", sample=3000):
    """Create 3D scatter plot of galaxy positions."""
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Sample for visibility
    if len(positions) > sample:
        indices = np.random.choice(len(positions), sample, replace=False)
        sample_pos = positions[indices]
    else:
        sample_pos = positions

    # Compute local density for coloring
    tree = cKDTree(sample_pos)
    densities = []
    for pos in sample_pos:
        neighbors = tree.query_ball_point(pos, r=10)
        densities.append(len(neighbors))

    densities = np.array(densities)

    # Plot
    scatter = ax.scatter(sample_pos[:, 0], sample_pos[:, 1], sample_pos[:, 2],
                        c=densities, cmap='viridis', s=5, alpha=0.6)

    ax.set_xlabel('X (Mpc/h)', fontsize=10)
    ax.set_ylabel('Y (Mpc/h)', fontsize=10)
    ax.set_zlabel('Z (Mpc/h)', fontsize=10)
    ax.set_title(title, fontsize=14, fontweight='bold')

    cbar = plt.colorbar(scatter, ax=ax, pad=0.1, shrink=0.8)
    cbar.set_label('Local Density', fontsize=10)

    return fig


def visualize_density_slices(positions, box_size=None, n_slices=6):
    """Create 2D density slices through the 3D distribution."""
    if box_size is None:
        box_size = positions.max() - positions.min()

    # Create grid
    grid_size = 128
    density_3d = np.zeros((grid_size, grid_size, grid_size))

    # Normalize positions
    pos_norm = (positions - positions.min()) / box_size * (grid_size - 1)

    # Fill grid
    for pos in pos_norm:
        i, j, k = pos.astype(int)
        if 0 <= i < grid_size and 0 <= j < grid_size and 0 <= k < grid_size:
            density_3d[i, j, k] += 1

    # Smooth
    density_3d = gaussian_filter(density_3d, sigma=2.0)

    # Create slices
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    slice_indices = np.linspace(10, grid_size-10, n_slices, dtype=int)

    for idx, slice_z in enumerate(slice_indices):
        ax = axes[idx]

        slice_2d = density_3d[:, :, slice_z]

        im = ax.imshow(slice_2d.T, origin='lower', cmap='hot', interpolation='bilinear')
        ax.set_title(f'Slice z={slice_z}/{grid_size}', fontsize=11)
        ax.axis('off')

        plt.colorbar(im, ax=ax, fraction=0.046)

    plt.suptitle('Density Field Slices Through Volume', fontsize=14, fontweight='bold')
    plt.tight_layout()

    return fig


def visualize_local_environments(positions, labels, sample=2000):
    """
    Visualize the clustered local environment types in 3D.

    Different colors = different local structure types.
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Sample
    if len(positions) > sample:
        indices = np.random.choice(len(positions), sample, replace=False)
        sample_pos = positions[indices]
        sample_labels = labels[indices] if len(labels) > sample else labels
    else:
        sample_pos = positions
        sample_labels = labels

    # Plot with different colors for each cluster
    scatter = ax.scatter(sample_pos[:, 0], sample_pos[:, 1], sample_pos[:, 2],
                        c=sample_labels, cmap='tab10', s=10, alpha=0.7)

    ax.set_xlabel('X (Mpc/h)', fontsize=10)
    ax.set_ylabel('Y (Mpc/h)', fontsize=10)
    ax.set_zlabel('Z (Mpc/h)', fontsize=10)
    ax.set_title('Local Environment Clusters', fontsize=14, fontweight='bold')

    cbar = plt.colorbar(scatter, ax=ax, pad=0.1, shrink=0.8)
    cbar.set_label('Environment Type', fontsize=10)

    return fig


def visualize_cosmic_web_structure(positions, sample=5000):
    """
    Visualize the cosmic web structure (filaments, nodes, voids).
    """
    fig = plt.figure(figsize=(14, 10))

    # Sample
    if len(positions) > sample:
        indices = np.random.choice(len(positions), sample, replace=False)
        sample_pos = positions[indices]
    else:
        sample_pos = positions

    # Compute local density
    tree = cKDTree(sample_pos)
    densities = []
    for pos in sample_pos:
        neighbors = tree.query_ball_point(pos, r=15)
        densities.append(len(neighbors))

    densities = np.array(densities)

    # Classify by density
    high_threshold = np.percentile(densities, 75)
    medium_threshold = np.percentile(densities, 25)

    colors = np.zeros(len(densities))
    colors[densities >= high_threshold] = 2  # Filaments/clusters (red)
    colors[(densities < high_threshold) & (densities >= medium_threshold)] = 1  # Sheets (yellow)
    colors[densities < medium_threshold] = 0  # Voids (blue)

    # 3D plot
    ax = fig.add_subplot(111, projection='3d')

    # Plot each type separately for legend
    for env_type, color, label in [(0, 'blue', 'Voids'),
                                     (1, 'yellow', 'Sheets'),
                                     (2, 'red', 'Filaments/Clusters')]:
        mask = colors == env_type
        ax.scatter(sample_pos[mask, 0], sample_pos[mask, 1], sample_pos[mask, 2],
                  c=color, s=10, alpha=0.6, label=label)

    ax.set_xlabel('X (Mpc/h)', fontsize=10)
    ax.set_ylabel('Y (Mpc/h)', fontsize=10)
    ax.set_zlabel('Z (Mpc/h)', fontsize=10)
    ax.set_title('Cosmic Web Structure', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)

    return fig


def create_comparison_view(synthetic_data, illustris_data):
    """Side-by-side comparison of synthetic vs Illustris."""
    fig = plt.figure(figsize=(16, 7))

    # Synthetic
    ax1 = fig.add_subplot(121, projection='3d')

    sample_size = 2000
    synth_sample = synthetic_data[np.random.choice(len(synthetic_data), min(sample_size, len(synthetic_data)), replace=False)]

    tree = cKDTree(synth_sample)
    densities_synth = [len(tree.query_ball_point(pos, r=10)) for pos in synth_sample]

    scatter1 = ax1.scatter(synth_sample[:, 0], synth_sample[:, 1], synth_sample[:, 2],
                          c=densities_synth, cmap='viridis', s=5, alpha=0.6)

    ax1.set_xlabel('X (Mpc/h)')
    ax1.set_ylabel('Y (Mpc/h)')
    ax1.set_zlabel('Z (Mpc/h)')
    ax1.set_title('Our Synthetic Data', fontsize=13, fontweight='bold')
    plt.colorbar(scatter1, ax=ax1, shrink=0.7)

    # Illustris
    ax2 = fig.add_subplot(122, projection='3d')

    ill_sample = illustris_data[np.random.choice(len(illustris_data), min(sample_size, len(illustris_data)), replace=False)]

    tree2 = cKDTree(ill_sample)
    densities_ill = [len(tree2.query_ball_point(pos, r=10)) for pos in ill_sample]

    scatter2 = ax2.scatter(ill_sample[:, 0], ill_sample[:, 1], ill_sample[:, 2],
                          c=densities_ill, cmap='viridis', s=5, alpha=0.6)

    ax2.set_xlabel('X (Mpc/h)')
    ax2.set_ylabel('Y (Mpc/h)')
    ax2.set_zlabel('Z (Mpc/h)')
    ax2.set_title('Illustris-like Data', fontsize=13, fontweight='bold')
    plt.colorbar(scatter2, ax=ax2, shrink=0.7)

    plt.suptitle('Comparison: Our Method vs Illustris-like', fontsize=15, fontweight='bold')
    plt.tight_layout()

    return fig


def main():
    """Generate all 3D visualizations."""
    print("="*70)
    print("GENERATING 3D VISUALIZATIONS")
    print("="*70)

    # Load data
    print("\nLoading datasets...")
    synth_data = np.load('realistic_z0.10.npz')
    synthetic_positions = synth_data['positions']

    ill_data = np.load('illustris_realization_1.npz')
    illustris_positions = ill_data['positions']

    # Load pattern discovery results
    pattern_results = np.load('pattern_discovery_results.npy', allow_pickle=True).item()

    # Extract environment labels
    env_result = pattern_results['z0.10']['local_envs']
    env_labels = env_result['labels']

    # Generate visualizations
    print("\n1. Creating 3D distribution plot...")
    fig1 = visualize_3d_distribution(synthetic_positions, "3D Galaxy Distribution (z=0.10)")
    fig1.savefig('3d_distribution.png', dpi=150, bbox_inches='tight')
    print("   Saved: 3d_distribution.png")
    plt.close(fig1)

    print("\n2. Creating density slices...")
    fig2 = visualize_density_slices(synthetic_positions)
    fig2.savefig('3d_density_slices.png', dpi=150, bbox_inches='tight')
    print("   Saved: 3d_density_slices.png")
    plt.close(fig2)

    print("\n3. Creating local environment visualization...")
    # Need to match positions to labels (labels are for sampled points)
    sample_size = len(env_labels)
    sampled_positions = synthetic_positions[:sample_size]

    fig3 = visualize_local_environments(sampled_positions, env_labels)
    fig3.savefig('3d_local_environments.png', dpi=150, bbox_inches='tight')
    print("   Saved: 3d_local_environments.png")
    plt.close(fig3)

    print("\n4. Creating cosmic web structure...")
    fig4 = visualize_cosmic_web_structure(synthetic_positions)
    fig4.savefig('3d_cosmic_web.png', dpi=150, bbox_inches='tight')
    print("   Saved: 3d_cosmic_web.png")
    plt.close(fig4)

    print("\n5. Creating comparison view...")
    fig5 = create_comparison_view(synthetic_positions, illustris_positions)
    fig5.savefig('3d_comparison.png', dpi=150, bbox_inches='tight')
    print("   Saved: 3d_comparison.png")
    plt.close(fig5)

    print("\n" + "="*70)
    print("3D VISUALIZATION COMPLETE!")
    print("="*70)
    print("\nGenerated files:")
    print("  - 3d_distribution.png")
    print("  - 3d_density_slices.png")
    print("  - 3d_local_environments.png")
    print("  - 3d_cosmic_web.png")
    print("  - 3d_comparison.png")


if __name__ == "__main__":
    main()
