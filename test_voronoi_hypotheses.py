"""
Test specific hypotheses about the Voronoi neighbor count anomaly.

We found ~25 neighbors instead of expected ~15.5 for random Poisson.

HYPOTHESES TO TEST:
H1: Filamentary structure causes high neighbor counts
H2: It's a sampling artifact (grid-based generation)
H3: It's a real ΛCDM property (will appear in Illustris too)
H4: There's a geometric constraint beyond gravity
H5: Void regions vs filament regions have different counts

Method: Statistical comparison across:
- Our synthetic data
- Illustris-like data
- Pure random Poisson distribution
- Filament-only vs void-only subsamples
"""

import numpy as np
from scipy.spatial import Voronoi, cKDTree
from scipy.stats import ks_2samp, mannwhitneyu
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def compute_voronoi_neighbors(positions, sample_size=2000):
    """
    Compute Voronoi neighbor count distribution.

    Neighbors = cells sharing a face, counted via ridge_points adjacency.
    Note: vor.regions gives vertex indices (not neighbors) — do not use
    len(region) as a neighbor count, that was the original bug giving ~25.

    Returns both summary statistics and full distribution.
    """
    if len(positions) > sample_size:
        indices = np.random.choice(len(positions), sample_size, replace=False)
        sample_points = positions[indices]
    else:
        sample_points = positions

    print(f"    Computing Voronoi for {len(sample_points)} points...")

    try:
        vor = Voronoi(sample_points)

        # Count face-adjacent neighbors via ridge_points (correct method).
        # Each entry in ridge_points is a pair (p1, p2) sharing a Voronoi face.
        neighbor_counts = np.zeros(len(sample_points), dtype=int)
        for p1, p2 in vor.ridge_points:
            neighbor_counts[p1] += 1
            neighbor_counts[p2] += 1

        return {
            'mean': np.mean(neighbor_counts),
            'std': np.std(neighbor_counts),
            'median': np.median(neighbor_counts),
            'distribution': neighbor_counts,
            'n_cells': len(neighbor_counts)
        }

    except Exception as e:
        print(f"    Voronoi failed: {e}")
        return None


def classify_environment(positions, percentile_high=75, percentile_low=25):
    """
    Classify points as filament (high density) or void (low density).

    Method: Compute local density, split into quantiles.
    """
    print("    Classifying environments...")

    tree = cKDTree(positions)

    # Compute local density (number of neighbors within fixed radius)
    r_probe = 10.0  # Mpc/h
    neighbor_counts = tree.query_ball_point(positions, r=r_probe)
    densities = np.array([len(n) for n in neighbor_counts])

    # Split into high/medium/low density
    high_threshold = np.percentile(densities, percentile_high)
    low_threshold = np.percentile(densities, percentile_low)

    high_density_mask = densities >= high_threshold  # Filaments/clusters
    low_density_mask = densities <= low_threshold   # Voids

    return {
        'high_density_indices': np.where(high_density_mask)[0],
        'low_density_indices': np.where(low_density_mask)[0],
        'densities': densities
    }


def test_hypothesis_1_filaments(positions):
    """
    H1: Filamentary structure causes high neighbor counts.

    Test: Compare Voronoi neighbors in high-density vs low-density regions.
    If H1 is true, filaments should have higher neighbor counts.
    """
    print("\n" + "="*70)
    print("HYPOTHESIS 1: Filamentary Structure Effect")
    print("="*70)

    # Classify environments
    env_class = classify_environment(positions)

    high_indices = env_class['high_density_indices']
    low_indices = env_class['low_density_indices']

    print(f"\n  High-density points (filaments): {len(high_indices)}")
    print(f"  Low-density points (voids): {len(low_indices)}")

    # Compute Voronoi for each subset
    print("\n  Voronoi in HIGH-density regions:")
    high_voronoi = compute_voronoi_neighbors(positions[high_indices], sample_size=1000)

    print("\n  Voronoi in LOW-density regions:")
    low_voronoi = compute_voronoi_neighbors(positions[low_indices], sample_size=1000)

    if high_voronoi is not None and low_voronoi is not None:
        print(f"\n  Results:")
        print(f"    High-density mean neighbors: {high_voronoi['mean']:.2f} ± {high_voronoi['std']:.2f}")
        print(f"    Low-density mean neighbors: {low_voronoi['mean']:.2f} ± {low_voronoi['std']:.2f}")

        # Statistical test
        stat, pval = mannwhitneyu(high_voronoi['distribution'], low_voronoi['distribution'])

        print(f"\n  Mann-Whitney U test:")
        print(f"    U-statistic: {stat:.2f}")
        print(f"    p-value: {pval:.4e}")

        if pval < 0.01:
            print(f"    ✅ SIGNIFICANT DIFFERENCE (p < 0.01)")
            print(f"    → H1 SUPPORTED: Filaments have different Voronoi structure!")
        else:
            print(f"    ❌ No significant difference (p >= 0.01)")
            print(f"    → H1 NOT SUPPORTED")

        return {
            'high': high_voronoi,
            'low': low_voronoi,
            'p_value': pval,
            'supported': pval < 0.01
        }

    return None


def test_hypothesis_2_sampling(synthetic_data, illustris_data):
    """
    H2: High neighbor count is a sampling artifact.

    Test: Compare grid-generated synthetic data vs Illustris-like data.
    If H2 is true, they should differ significantly.
    """
    print("\n" + "="*70)
    print("HYPOTHESIS 2: Sampling Artifact")
    print("="*70)

    print("\n  Computing Voronoi for SYNTHETIC data...")
    synthetic_voronoi = compute_voronoi_neighbors(synthetic_data)

    print("\n  Computing Voronoi for ILLUSTRIS-like data...")
    illustris_voronoi = compute_voronoi_neighbors(illustris_data)

    if synthetic_voronoi is not None and illustris_voronoi is not None:
        print(f"\n  Results:")
        print(f"    Synthetic mean neighbors: {synthetic_voronoi['mean']:.2f} ± {synthetic_voronoi['std']:.2f}")
        print(f"    Illustris mean neighbors: {illustris_voronoi['mean']:.2f} ± {illustris_voronoi['std']:.2f}")

        # Statistical test
        stat, pval = ks_2samp(synthetic_voronoi['distribution'], illustris_voronoi['distribution'])

        print(f"\n  Kolmogorov-Smirnov test:")
        print(f"    KS-statistic: {stat:.4f}")
        print(f"    p-value: {pval:.4e}")

        if pval < 0.01:
            print(f"    ✅ SIGNIFICANT DIFFERENCE (p < 0.01)")
            print(f"    → H2 SUPPORTED: Sampling method matters!")
        else:
            print(f"    ❌ No significant difference (p >= 0.01)")
            print(f"    → H2 NOT SUPPORTED: Both methods give same result")

        return {
            'synthetic': synthetic_voronoi,
            'illustris': illustris_voronoi,
            'p_value': pval,
            'supported': pval < 0.01
        }

    return None


def test_hypothesis_3_lcdm_property():
    """
    H3: ~25 neighbors is a real ΛCDM property.

    Test: Compare against pure random Poisson distribution.
    If H3 is true, ALL structured data should show ~25, random should show ~15.
    """
    print("\n" + "="*70)
    print("HYPOTHESIS 3: Real ΛCDM Property")
    print("="*70)

    # Generate pure random Poisson
    print("\n  Generating pure random Poisson distribution...")
    box_size = 100
    n_points = 10000
    random_positions = np.random.uniform(0, box_size, (n_points, 3))

    print("\n  Computing Voronoi for RANDOM distribution...")
    random_voronoi = compute_voronoi_neighbors(random_positions, sample_size=2000)

    # Load one of our datasets
    data = np.load('realistic_z0.10.npz')
    structured_positions = data['positions']

    print("\n  Computing Voronoi for STRUCTURED distribution...")
    structured_voronoi = compute_voronoi_neighbors(structured_positions, sample_size=2000)

    if random_voronoi is not None and structured_voronoi is not None:
        print(f"\n  Results:")
        print(f"    Random Poisson mean: {random_voronoi['mean']:.2f} ± {random_voronoi['std']:.2f}")
        print(f"    Structured (ΛCDM) mean: {structured_voronoi['mean']:.2f} ± {structured_voronoi['std']:.2f}")

        # Expected for random 3D Poisson: ~15.5
        print(f"\n  Comparison to theory:")
        print(f"    Expected random 3D Poisson: ~15.5")
        print(f"    Observed random: {random_voronoi['mean']:.2f}")
        print(f"    Observed structured: {structured_voronoi['mean']:.2f}")

        # Statistical test
        stat, pval = ks_2samp(random_voronoi['distribution'], structured_voronoi['distribution'])

        print(f"\n  Kolmogorov-Smirnov test:")
        print(f"    KS-statistic: {stat:.4f}")
        print(f"    p-value: {pval:.4e}")

        if pval < 0.01 and structured_voronoi['mean'] > random_voronoi['mean']:
            print(f"    ✅ SIGNIFICANT DIFFERENCE (p < 0.01)")
            print(f"    → H3 SUPPORTED: ΛCDM structure has more neighbors than random!")
        else:
            print(f"    ❌ Not significantly different")

        return {
            'random': random_voronoi,
            'structured': structured_voronoi,
            'p_value': pval,
            'supported': pval < 0.01 and structured_voronoi['mean'] > random_voronoi['mean']
        }

    return None


def test_hypothesis_4_geometric_constraint(all_datasets):
    """
    H4: There's a universal geometric constraint.

    Test: If neighbor count is remarkably consistent across ALL datasets
    (different methods, different parameters), it suggests a fundamental principle.
    """
    print("\n" + "="*70)
    print("HYPOTHESIS 4: Universal Geometric Constraint")
    print("="*70)

    print("\n  Computing Voronoi for ALL datasets...")

    results = []

    for name, positions in all_datasets.items():
        print(f"\n  Dataset: {name}")
        voronoi_stats = compute_voronoi_neighbors(positions, sample_size=1500)

        if voronoi_stats is not None:
            results.append({
                'name': name,
                'mean': voronoi_stats['mean'],
                'std': voronoi_stats['std']
            })

    if len(results) > 0:
        means = [r['mean'] for r in results]
        overall_mean = np.mean(means)
        overall_std = np.std(means)

        print(f"\n  Results across all datasets:")
        for r in results:
            print(f"    {r['name']:30s}: {r['mean']:.2f} ± {r['std']:.2f}")

        print(f"\n  Overall statistics:")
        print(f"    Mean of means: {overall_mean:.2f}")
        print(f"    Std of means: {overall_std:.2f}")
        print(f"    Coefficient of variation: {overall_std/overall_mean:.2%}")

        if overall_std / overall_mean < 0.10:  # Less than 10% variation
            print(f"\n    ✅ REMARKABLY CONSISTENT (<10% variation)")
            print(f"    → H4 SUPPORTED: Suggests universal geometric principle!")
        else:
            print(f"\n    ❌ Significant variation (>10%)")
            print(f"    → H4 NOT SUPPORTED")

        return {
            'results': results,
            'overall_mean': overall_mean,
            'overall_std': overall_std,
            'cv': overall_std / overall_mean,
            'supported': overall_std / overall_mean < 0.10
        }

    return None


def visualize_voronoi_comparison(results_dict):
    """Create comprehensive visualization of Voronoi analysis."""
    print("\n  Generating visualization...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Plot 1: Neighbor distribution comparison
    ax1 = axes[0, 0]

    if 'h3_results' in results_dict and results_dict['h3_results'] is not None:
        h3 = results_dict['h3_results']

        random_dist = h3['random']['distribution']
        structured_dist = h3['structured']['distribution']

        bins = np.arange(0, 50, 2)
        ax1.hist(random_dist, bins=bins, alpha=0.6, label='Random Poisson', density=True, edgecolor='black')
        ax1.hist(structured_dist, bins=bins, alpha=0.6, label='ΛCDM Structure', density=True, edgecolor='black')

        ax1.axvline(15.5, color='red', linestyle='--', label='Expected Random', linewidth=2)
        ax1.set_xlabel('Number of Neighbors', fontsize=11)
        ax1.set_ylabel('Probability Density', fontsize=11)
        ax1.set_title('Voronoi Neighbor Distribution', fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

    # Plot 2: Environment comparison
    ax2 = axes[0, 1]

    if 'h1_results' in results_dict and results_dict['h1_results'] is not None:
        h1 = results_dict['h1_results']

        high_dist = h1['high']['distribution']
        low_dist = h1['low']['distribution']

        bins = np.arange(0, 50, 2)
        ax2.hist(high_dist, bins=bins, alpha=0.6, label='High Density (Filaments)', density=True, edgecolor='black')
        ax2.hist(low_dist, bins=bins, alpha=0.6, label='Low Density (Voids)', density=True, edgecolor='black')

        ax2.set_xlabel('Number of Neighbors', fontsize=11)
        ax2.set_ylabel('Probability Density', fontsize=11)
        ax2.set_title('Voronoi by Environment Type', fontweight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

    # Plot 3: Cross-dataset comparison
    ax3 = axes[1, 0]

    if 'h4_results' in results_dict and results_dict['h4_results'] is not None:
        h4 = results_dict['h4_results']

        names = [r['name'] for r in h4['results']]
        means = [r['mean'] for r in h4['results']]
        stds = [r['std'] for r in h4['results']]

        x = np.arange(len(names))
        ax3.bar(x, means, yerr=stds, alpha=0.7, capsize=5, edgecolor='black')
        ax3.set_xticks(x)
        ax3.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
        ax3.set_ylabel('Mean Neighbors', fontsize=11)
        ax3.set_title('Voronoi Consistency Across Datasets', fontweight='bold')
        ax3.axhline(y=15.5, color='red', linestyle='--', label='Expected Random', alpha=0.7)
        ax3.legend()
        ax3.grid(True, alpha=0.3, axis='y')

    # Plot 4: Summary statistics
    ax4 = axes[1, 1]
    ax4.axis('off')

    # Create summary text
    summary_text = "HYPOTHESIS TESTING SUMMARY\n"
    summary_text += "="*40 + "\n\n"

    if 'h1_results' in results_dict and results_dict['h1_results']:
        h1_result = "✅ SUPPORTED" if results_dict['h1_results']['supported'] else "❌ NOT SUPPORTED"
        summary_text += f"H1: Filament Effect\n  {h1_result}\n\n"

    if 'h2_results' in results_dict and results_dict['h2_results']:
        h2_result = "✅ SUPPORTED" if results_dict['h2_results']['supported'] else "❌ NOT SUPPORTED"
        summary_text += f"H2: Sampling Artifact\n  {h2_result}\n\n"

    if 'h3_results' in results_dict and results_dict['h3_results']:
        h3_result = "✅ SUPPORTED" if results_dict['h3_results']['supported'] else "❌ NOT SUPPORTED"
        summary_text += f"H3: Real ΛCDM Property\n  {h3_result}\n\n"

    if 'h4_results' in results_dict and results_dict['h4_results']:
        h4_result = "✅ SUPPORTED" if results_dict['h4_results']['supported'] else "❌ NOT SUPPORTED"
        summary_text += f"H4: Geometric Constraint\n  {h4_result}\n\n"

    ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes,
             fontsize=11, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plt.savefig('voronoi_hypothesis_tests.png', dpi=150, bbox_inches='tight')
    print("  Saved: voronoi_hypothesis_tests.png")
    plt.close()


def main():
    """Run all hypothesis tests."""
    print("="*70)
    print("VORONOI ANOMALY: HYPOTHESIS TESTING")
    print("="*70)

    results = {}

    # Load datasets
    print("\nLoading datasets...")

    # Our synthetic data
    synth_data = np.load('realistic_z0.10.npz')
    synthetic_positions = synth_data['positions']

    # Illustris-like data
    ill_data = np.load('illustris_realization_1.npz')
    illustris_positions = ill_data['positions']

    # All datasets for H4
    all_datasets = {
        'Synthetic z=0.10': synthetic_positions,
        'Illustris-like #1': illustris_positions,
    }

    # Load more for robustness
    for i in range(2, 6):
        try:
            data = np.load(f'illustris_realization_{i}.npz')
            all_datasets[f'Illustris-like #{i}'] = data['positions']
        except:
            pass

    try:
        for z in ['0.25', '0.45', '0.65']:
            data = np.load(f'realistic_z{z}.npz')
            all_datasets[f'Synthetic z={z}'] = data['positions']
    except:
        pass

    print(f"Loaded {len(all_datasets)} datasets for comparison")

    # Run tests
    results['h1_results'] = test_hypothesis_1_filaments(synthetic_positions)
    results['h2_results'] = test_hypothesis_2_sampling(synthetic_positions, illustris_positions)
    results['h3_results'] = test_hypothesis_3_lcdm_property()
    results['h4_results'] = test_hypothesis_4_geometric_constraint(all_datasets)

    # Visualize
    print("\n" + "="*70)
    print("GENERATING VISUALIZATIONS")
    print("="*70)
    visualize_voronoi_comparison(results)

    # Save results
    np.save('voronoi_hypothesis_results.npy', results, allow_pickle=True)
    print("\nSaved: voronoi_hypothesis_results.npy")

    print("\n" + "="*70)
    print("HYPOTHESIS TESTING COMPLETE!")
    print("="*70)

    return results


if __name__ == "__main__":
    results = main()
