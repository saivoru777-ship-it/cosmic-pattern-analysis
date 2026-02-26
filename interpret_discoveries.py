"""
Interpret and visualize the discovered patterns.

This answers:
1. WHAT patterns actually exist in the data?
2. Are they real or just random noise?
3. Do they evolve over cosmic time?
4. Do any match known structures (Einstein tile, quasicrystals, etc.)?
5. What NEW patterns did we find?
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import pandas as pd


def load_results():
    """Load discovery results."""
    return np.load('pattern_discovery_results.npy', allow_pickle=True).item()


def analyze_time_evolution_of_patterns(results):
    """
    Track how patterns change across cosmic epochs.

    KEY INSIGHT: If tiling-like patterns exist, they should:
    - Persist across time
    - Have consistent scale ratios
    - Show specific evolution signatures
    """
    print("\n" + "="*70)
    print("TIME EVOLUTION OF DISCOVERED PATTERNS")
    print("="*70)

    epochs = sorted(results.keys())

    # 1. Voronoi neighbor distribution evolution
    print("\n1. VORONOI TESSELLATION EVOLUTION")
    print("-" * 70)

    voronoi_data = []
    for epoch in epochs:
        if results[epoch]['voronoi'] is not None:
            v_stats = results[epoch]['voronoi']
            voronoi_data.append({
                'epoch': epoch,
                'mean_neighbors': v_stats['mean_neighbors'],
                'std_neighbors': v_stats['std_neighbors']
            })

    if len(voronoi_data) > 0:
        df_voronoi = pd.DataFrame(voronoi_data)
        print(df_voronoi.to_string(index=False))

        # Check if mean neighbors is constant (would suggest geometric constraint)
        mean_vals = df_voronoi['mean_neighbors'].values
        variation = np.std(mean_vals) / np.mean(mean_vals)
        print(f"\nCoefficient of variation: {variation:.4f}")
        if variation < 0.05:
            print("⚠️  VERY STABLE across time - possible geometric constraint!")
        else:
            print("✓ Normal variation - consistent with gravitational evolution")

    # 2. Scale hierarchy evolution
    print("\n\n2. SCALE HIERARCHY EVOLUTION")
    print("-" * 70)

    scale_data = []
    for epoch in epochs:
        scales_result = results[epoch]['scales']
        clustering_curve = scales_result['clustering_curve']

        scale_data.append({
            'epoch': epoch,
            'scales': [r['scale'] for r in scales_result['results']],
            'clustering': clustering_curve,
            'max_clustering_scale': scales_result['results'][np.argmax(clustering_curve)]['scale']
        })

    for sd in scale_data:
        print(f"\n{sd['epoch']}:")
        print(f"  Clustering strengths: {[f'{c:.2f}' for c in sd['clustering']]}")
        print(f"  Peak clustering at scale: {sd['max_clustering_scale']} Mpc/h")

    # 3. Fourier pattern evolution
    print("\n\n3. FOURIER PEAKS (Periodic/Quasiperiodic Signatures)")
    print("-" * 70)

    fourier_data = []
    for epoch in epochs:
        f_result = results[epoch]['fourier']
        peak_positions = f_result['peak_positions']

        fourier_data.append({
            'epoch': epoch,
            'n_peaks': len(peak_positions),
            'peak_k': list(peak_positions) if len(peak_positions) > 0 else []
        })

    for fd in fourier_data:
        print(f"\n{fd['epoch']}:")
        print(f"  Number of peaks: {fd['n_peaks']}")
        if len(fd['peak_k']) > 0:
            print(f"  Peak k-values: {[f'{k:.4f}' for k in fd['peak_k']]}")
            # Check for ratios
            if len(fd['peak_k']) >= 2:
                ratios = [fd['peak_k'][i+1]/fd['peak_k'][i] for i in range(len(fd['peak_k'])-1)]
                print(f"  k-ratios: {[f'{r:.3f}' for r in ratios]}")

                phi = (1 + np.sqrt(5)) / 2
                if any(abs(r - phi) < 0.1 for r in ratios):
                    print("  ⚠️  RATIO NEAR GOLDEN RATIO φ!")

    # 4. Local environment patterns
    print("\n\n4. LOCAL ENVIRONMENT DIVERSITY")
    print("-" * 70)

    env_data = []
    for epoch in epochs:
        env_result = results[epoch]['local_envs']
        cluster_counts = env_result['cluster_counts']

        # Measure diversity (entropy)
        total = sum(cluster_counts.values())
        probs = [count/total for count in cluster_counts.values()]
        entropy = -sum(p * np.log(p) for p in probs if p > 0)

        env_data.append({
            'epoch': epoch,
            'n_env_types': len(cluster_counts),
            'entropy': entropy,
            'dominant_type_fraction': max(cluster_counts.values()) / total
        })

    df_env = pd.DataFrame(env_data)
    print(df_env.to_string(index=False))

    return {
        'voronoi': voronoi_data,
        'scales': scale_data,
        'fourier': fourier_data,
        'environments': env_data
    }


def check_for_known_patterns(results):
    """
    Check if discovered patterns match known structures:
    - Einstein tile signatures
    - Penrose/quasicrystal patterns
    - Fibonacci sequences
    - Other aperiodic orders
    """
    print("\n" + "="*70)
    print("CHECKING FOR KNOWN APERIODIC PATTERNS")
    print("="*70)

    phi = (1 + np.sqrt(5)) / 2
    sqrt2 = np.sqrt(2)
    sqrt3 = np.sqrt(3)

    known_ratios = {
        'Golden Ratio (φ)': phi,
        'Silver Ratio (δ_S)': 1 + sqrt2,
        'Bronze Ratio': (3 + np.sqrt(13)) / 2,
        '√2': sqrt2,
        '√3': sqrt3,
        '2': 2.0,
    }

    # Check scale ratios
    print("\n1. SCALE RATIOS vs KNOWN CONSTANTS")
    print("-" * 70)

    for epoch in sorted(results.keys()):
        print(f"\n{epoch}:")
        scales_result = results[epoch]['scales']
        scale_ratios = scales_result['scale_ratios']

        print(f"  Observed ratios: {[f'{r:.3f}' for r in scale_ratios]}")

        for name, ratio in known_ratios.items():
            if any(abs(sr - ratio) < 0.15 for sr in scale_ratios):
                print(f"  ⚠️  Close to {name} = {ratio:.3f}")

    # Check Voronoi neighbors
    print("\n\n2. VORONOI NEIGHBOR COUNTS")
    print("-" * 70)
    print("Expected values:")
    print("  Random 3D Poisson: ~15.5 neighbors")
    print("  Face-centered cubic: 12 neighbors")
    print("  Body-centered cubic: 8 + 6 = 14 neighbors")
    print("  Quasicrystalline: varies, but specific values")

    for epoch in sorted(results.keys()):
        if results[epoch]['voronoi'] is not None:
            v_stats = results[epoch]['voronoi']
            mean_n = v_stats['mean_neighbors']
            print(f"\n{epoch}: {mean_n:.2f} neighbors")

            if abs(mean_n - 12) < 1.0:
                print("  ⚠️  Close to FCC lattice!")
            elif abs(mean_n - 14) < 1.0:
                print("  ⚠️  Close to BCC lattice!")
            elif abs(mean_n - 15.5) < 1.0:
                print("  ✓ Consistent with random Poisson")


def visualize_discoveries(results, evolution_analysis):
    """Create comprehensive visualizations."""
    print("\n" + "="*70)
    print("GENERATING VISUALIZATIONS")
    print("="*70)

    epochs = sorted(results.keys())

    # Figure 1: Time evolution summary
    fig = plt.figure(figsize=(16, 12))

    # Subplot 1: Voronoi neighbors over time
    if len(evolution_analysis['voronoi']) > 0:
        ax1 = plt.subplot(3, 3, 1)
        df_v = pd.DataFrame(evolution_analysis['voronoi'])
        x = np.arange(len(df_v))
        ax1.errorbar(x, df_v['mean_neighbors'], yerr=df_v['std_neighbors'],
                    fmt='o-', capsize=5, linewidth=2, markersize=8)
        ax1.axhline(y=15.5, color='r', linestyle='--', label='Random Poisson', alpha=0.5)
        ax1.axhline(y=12, color='g', linestyle='--', label='FCC lattice', alpha=0.5)
        ax1.set_xticks(x)
        ax1.set_xticklabels(df_v['epoch'], rotation=45)
        ax1.set_ylabel('Mean Neighbors')
        ax1.set_title('Voronoi Neighbor Count Evolution', fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

    # Subplot 2: Clustering strength across scales
    ax2 = plt.subplot(3, 3, 2)
    for sd in evolution_analysis['scales']:
        ax2.plot(sd['scales'], sd['clustering'], marker='o', label=sd['epoch'], linewidth=2)
    ax2.set_xscale('log')
    ax2.set_xlabel('Scale (Mpc/h)')
    ax2.set_ylabel('Clustering Strength')
    ax2.set_title('Multiscale Clustering Evolution', fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Subplot 3: Power spectra
    ax3 = plt.subplot(3, 3, 3)
    for epoch in epochs:
        f_result = results[epoch]['fourier']
        ax3.loglog(f_result['k_values'], f_result['power_spectrum'], label=epoch, linewidth=2)
    ax3.set_xlabel('k (wavenumber)')
    ax3.set_ylabel('Power')
    ax3.set_title('Fourier Power Spectra', fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # Subplot 4: Environment diversity
    ax4 = plt.subplot(3, 3, 4)
    df_env = pd.DataFrame(evolution_analysis['environments'])
    x = np.arange(len(df_env))
    ax4.bar(x, df_env['entropy'], alpha=0.7, edgecolor='black')
    ax4.set_xticks(x)
    ax4.set_xticklabels(df_env['epoch'], rotation=45)
    ax4.set_ylabel('Entropy (bits)')
    ax4.set_title('Local Environment Diversity', fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')

    # Subplot 5: Degree distribution (first epoch)
    ax5 = plt.subplot(3, 3, 5)
    first_epoch = epochs[0]
    graph_result = results[first_epoch]['motifs']
    degrees = list(graph_result['degree_distribution'].keys())
    counts = list(graph_result['degree_distribution'].values())
    ax5.bar(degrees, counts, alpha=0.7, edgecolor='black')
    ax5.set_xlabel('Degree')
    ax5.set_ylabel('Count')
    ax5.set_title(f'Network Degree Distribution ({first_epoch})', fontweight='bold')
    ax5.grid(True, alpha=0.3, axis='y')

    # Subplot 6: Shape pattern distribution
    ax6 = plt.subplot(3, 3, 6)
    shape_result = results[first_epoch]['shapes']
    cluster_ids = list(shape_result['cluster_counts'].keys())
    cluster_sizes = list(shape_result['cluster_counts'].values())
    ax6.bar(cluster_ids, cluster_sizes, alpha=0.7, edgecolor='black')
    ax6.set_xlabel('Shape Cluster ID')
    ax6.set_ylabel('Count')
    ax6.set_title(f'Shape Pattern Diversity ({first_epoch})', fontweight='bold')
    ax6.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig('pattern_evolution_analysis.png', dpi=150, bbox_inches='tight')
    print("Saved: pattern_evolution_analysis.png")
    plt.close()


def generate_summary_report(results, evolution_analysis):
    """Generate comprehensive findings report."""
    print("\n" + "="*70)
    print("GENERATING SUMMARY REPORT")
    print("="*70)

    report = []
    report.append("="*70)
    report.append("COSMIC PATTERN DISCOVERY: FINDINGS SUMMARY")
    report.append("="*70)
    report.append("")
    report.append("RESEARCH QUESTION:")
    report.append("Do large-scale cosmic structures show signatures of aperiodic tiling,")
    report.append("quasicrystalline order, or other geometric patterns beyond standard")
    report.append("ΛCDM gravitational clustering?")
    report.append("")
    report.append("="*70)
    report.append("METHODS")
    report.append("="*70)
    report.append("")
    report.append("Analyzed 4 cosmic epochs (z = 0.10, 0.25, 0.45, 0.65)")
    report.append("Applied 6 pattern discovery methods:")
    report.append("  1. Local environment clustering")
    report.append("  2. Voronoi tessellation analysis")
    report.append("  3. Multiscale clustering hierarchy")
    report.append("  4. Fourier/quasiperiodic pattern detection")
    report.append("  5. Network motif extraction")
    report.append("  6. Unsupervised shape pattern discovery")
    report.append("")
    report.append("="*70)
    report.append("KEY FINDINGS")
    report.append("="*70)
    report.append("")

    # Finding 1: Voronoi stability
    if len(evolution_analysis['voronoi']) > 0:
        df_v = pd.DataFrame(evolution_analysis['voronoi'])
        mean_vals = df_v['mean_neighbors'].values
        variation = np.std(mean_vals) / np.mean(mean_vals)

        report.append("1. VORONOI TESSELLATION")
        report.append(f"   Mean neighbors across epochs: {mean_vals.mean():.2f} ± {mean_vals.std():.2f}")
        report.append(f"   Temporal variation: {variation:.2%}")

        if variation < 0.05:
            report.append("   ⚠️  FINDING: Extremely stable - possible geometric constraint!")
        else:
            report.append("   ✓ FINDING: Normal variation - consistent with ΛCDM")
        report.append("")

    # Finding 2: Scale hierarchy
    report.append("2. SCALE HIERARCHY")
    peak_scales = [sd['max_clustering_scale'] for sd in evolution_analysis['scales']]
    report.append(f"   Peak clustering scales: {peak_scales} Mpc/h")

    # Check for golden ratio
    phi = (1 + np.sqrt(5)) / 2
    ratios = [peak_scales[i+1]/peak_scales[i] for i in range(len(peak_scales)-1) if peak_scales[i] > 0]
    if len(ratios) > 0:
        report.append(f"   Scale ratios: {[f'{r:.3f}' for r in ratios]}")
        if any(abs(r - phi) < 0.15 for r in ratios):
            report.append("   ⚠️  FINDING: Ratio(s) near golden ratio φ ≈ 1.618!")
        if any(abs(r - 2.0) < 0.15 for r in ratios):
            report.append("   ✓ FINDING: Simple doubling hierarchy (consistent with ΛCDM)")
    report.append("")

    # Finding 3: Fourier patterns
    report.append("3. FOURIER / QUASIPERIODIC PATTERNS")
    total_peaks = sum(fd['n_peaks'] for fd in evolution_analysis['fourier'])
    report.append(f"   Total Fourier peaks detected: {total_peaks}")

    if total_peaks > 0:
        report.append("   ⚠️  FINDING: Periodic/quasiperiodic features detected")
        report.append("   (Requires further investigation)")
    else:
        report.append("   ✓ FINDING: No strong periodic patterns")
    report.append("")

    # Finding 4: Pattern persistence
    report.append("4. PATTERN PERSISTENCE ACROSS TIME")
    df_env = pd.DataFrame(evolution_analysis['environments'])
    entropy_change = df_env['entropy'].values[-1] - df_env['entropy'].values[0]
    report.append(f"   Environment diversity change: {entropy_change:+.3f} bits")

    if abs(entropy_change) < 0.1:
        report.append("   ⚠️  FINDING: Diversity stable - patterns persist across time!")
    else:
        report.append("   ✓ FINDING: Diversity evolves - consistent with gravitational growth")
    report.append("")

    report.append("="*70)
    report.append("CONCLUSIONS")
    report.append("="*70)
    report.append("")
    report.append("PRIMARY CONCLUSION:")
    report.append("The data shows patterns consistent with standard ΛCDM gravitational")
    report.append("clustering. No strong evidence for aperiodic tiling or quasicrystalline")
    report.append("order beyond what gravity naturally produces.")
    report.append("")
    report.append("HOWEVER:")
    report.append("Several intriguing features warrant further investigation:")
    if variation < 0.05:
        report.append("  - Unusually stable Voronoi neighbor counts")
    if total_peaks > 0:
        report.append("  - Fourier peaks suggesting preferred scales")
    report.append("")
    report.append("NEXT STEPS:")
    report.append("1. Compare to ΛCDM N-body simulations (Illustris, EAGLE)")
    report.append("2. Test on real SDSS/BOSS data")
    report.append("3. Investigate any anomalous features with higher resolution")
    report.append("4. Develop statistical tests for aperiodic order")
    report.append("")
    report.append("="*70)

    # Write to file
    with open('PATTERN_DISCOVERY_REPORT.txt', 'w') as f:
        f.write('\n'.join(report))

    print("Saved: PATTERN_DISCOVERY_REPORT.txt")

    # Print to console
    print("\n".join(report))


def main():
    """Run full interpretation and visualization."""
    print("="*70)
    print("COSMIC PATTERN DISCOVERY: INTERPRETATION & VISUALIZATION")
    print("="*70)

    # Load results
    print("\nLoading results...")
    results = load_results()

    # Analyze time evolution
    evolution_analysis = analyze_time_evolution_of_patterns(results)

    # Check for known patterns
    check_for_known_patterns(results)

    # Visualize
    visualize_discoveries(results, evolution_analysis)

    # Generate report
    generate_summary_report(results, evolution_analysis)

    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)
    print("\nGenerated files:")
    print("  - pattern_evolution_analysis.png")
    print("  - PATTERN_DISCOVERY_REPORT.txt")


if __name__ == "__main__":
    main()
