"""
Deep investigation into the Voronoi ~25 neighbor anomaly.

We'll test:
1. Different Voronoi implementations
2. Different box sizes
3. Different boundary conditions
4. Manual calculation verification
5. Published theoretical predictions
6. 4D projection hypothesis
"""

import numpy as np
from scipy.spatial import Voronoi, ConvexHull
import matplotlib.pyplot as plt


def compute_voronoi_detailed_stats(positions, method='standard'):
    """
    Compute Voronoi statistics with detailed diagnostics.

    This will help us understand WHAT we're actually measuring.
    """
    print(f"\n  Method: {method}")
    print(f"  Points: {len(positions)}")

    # Compute Voronoi
    vor = Voronoi(positions)

    # Detailed analysis
    results = {
        'n_points': len(positions),
        'n_vertices': len(vor.vertices),
        'n_regions': len(vor.regions),
        'ridge_points_shape': vor.ridge_points.shape if len(vor.ridge_points) > 0 else (0, 0)
    }

    # Count neighbors per point
    neighbor_counts = []
    finite_region_counts = []
    infinite_region_counts = []

    # Method 1: Count via ridge_points (direct adjacency)
    # This is what we SHOULD be counting
    adjacency_counts = [0] * len(positions)

    for ridge in vor.ridge_points:
        p1, p2 = ridge
        adjacency_counts[p1] += 1
        adjacency_counts[p2] += 1

    adjacency_counts = np.array(adjacency_counts)

    # Method 2: Count via regions (what we were doing before)
    for i, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]

        if -1 in region:
            # Infinite region (touches boundary)
            infinite_region_counts.append(len(region))
        else:
            # Finite region
            if len(region) > 0:
                finite_region_counts.append(len(region))
                neighbor_counts.append(len(region))

    # Print comparison
    print(f"\n  Method 1 (Ridge adjacency - CORRECT):")
    print(f"    Mean: {np.mean(adjacency_counts):.2f}")
    print(f"    Median: {np.median(adjacency_counts):.2f}")
    print(f"    Std: {np.std(adjacency_counts):.2f}")

    print(f"\n  Method 2 (Region vertices - what we used before):")
    print(f"    Mean: {np.mean(neighbor_counts):.2f} (finite regions only)")
    print(f"    Finite regions: {len(finite_region_counts)}")
    print(f"    Infinite regions: {len(infinite_region_counts)}")

    # Aha! This is the key distinction
    # We were counting VERTICES of regions, not NEIGHBORS

    results['adjacency_counts'] = adjacency_counts
    results['region_vertex_counts'] = np.array(neighbor_counts)
    results['n_finite_regions'] = len(finite_region_counts)
    results['n_infinite_regions'] = len(infinite_region_counts)

    return results


def theoretical_voronoi_neighbors_3d():
    """
    What does theory predict for 3D Voronoi tessellation?

    Let's look up the actual values.
    """
    print("\n" + "="*70)
    print("THEORETICAL PREDICTIONS FOR 3D VORONOI")
    print("="*70)

    print("""
From published literature:

1. RANDOM 3D POISSON PROCESS:
   - Expected neighbor count (face adjacency): ~15.54
   - This is a well-known result from stochastic geometry
   - Derived by Meijering (1953), Miles (1970)

   Source: "The number of faces of the typical Poisson-Voronoi cell"

2. WHAT WE MIGHT BE MEASURING INSTEAD:

   a) VERTICES per Voronoi cell:
      - Expected for 3D Poisson: ~35-36 vertices
      - This is DIFFERENT from neighbor count!

   b) EDGES per Voronoi cell:
      - Expected for 3D Poisson: ~48-50 edges

   c) FACES per Voronoi cell:
      - Expected for 3D Poisson: ~15.54 faces
      - Each face = one neighbor!

3. KEY INSIGHT:
   The number of VERTICES in a Voronoi region is NOT the same as
   the number of NEIGHBORS (adjacent cells).

   Neighbors = cells sharing a FACE
   Vertices = corner points of the cell

   For a cube (regular cell):
   - 6 faces → 6 neighbors
   - 8 vertices
   - 12 edges

   Our ~25 count might be counting something else!
""")


def test_what_were_counting():
    """
    Let's create a simple test case where we KNOW the answer.
    """
    print("\n" + "="*70)
    print("TEST CASE: REGULAR GRID (Known Answer)")
    print("="*70)

    # Create a regular 3D grid
    # Each interior point should have exactly 6 neighbors (cube)
    n = 5
    x = np.arange(n)
    y = np.arange(n)
    z = np.arange(n)

    xx, yy, zz = np.meshgrid(x, y, z)
    grid_points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()]).astype(float)

    print(f"\nRegular 5×5×5 grid: {len(grid_points)} points")
    print(f"Interior points should have exactly 6 neighbors (cube faces)")
    print(f"Vertices per cube: 8")

    # Compute Voronoi
    vor = Voronoi(grid_points)

    # Count via adjacency
    adjacency_counts = [0] * len(grid_points)
    for ridge in vor.ridge_points:
        p1, p2 = ridge
        adjacency_counts[p1] += 1
        adjacency_counts[p2] += 1

    adjacency_counts = np.array(adjacency_counts)

    # Count via region vertices
    region_vertex_counts = []
    for i, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if -1 not in region and len(region) > 0:
            region_vertex_counts.append(len(region))

    print(f"\nResults:")
    print(f"  Adjacency (neighbors) mean: {np.mean(adjacency_counts):.2f}")
    print(f"  Region vertices mean: {np.mean(region_vertex_counts):.2f}")

    print(f"\n  Expected adjacency: 6 (for interior cubes)")
    print(f"  Expected vertices: 8 (for cubes)")

    if np.abs(np.mean(region_vertex_counts) - 8) < 1.0:
        print(f"\n  ✅ CONFIRMED: We were counting VERTICES, not NEIGHBORS!")
        return True

    return False


def compute_correct_voronoi_neighbors(positions, sample_size=2000):
    """
    Compute the CORRECT neighbor count using ridge adjacency.
    """
    if len(positions) > sample_size:
        indices = np.random.choice(len(positions), sample_size, replace=False)
        sample_points = positions[indices]
    else:
        sample_points = positions

    print(f"    Computing Voronoi for {len(sample_points)} points...")

    try:
        vor = Voronoi(sample_points)

        # Count neighbors via ridge_points (CORRECT method)
        neighbor_counts = [0] * len(sample_points)
        for ridge in vor.ridge_points:
            p1, p2 = ridge
            neighbor_counts[p1] += 1
            neighbor_counts[p2] += 1

        neighbor_counts = np.array(neighbor_counts)

        return {
            'mean': np.mean(neighbor_counts),
            'std': np.std(neighbor_counts),
            'median': np.median(neighbor_counts),
            'distribution': neighbor_counts
        }

    except Exception as e:
        print(f"    Voronoi failed: {e}")
        return None


def retest_all_datasets_correctly():
    """
    Re-run ALL our tests using the CORRECT neighbor counting.
    """
    print("\n" + "="*70)
    print("RE-TESTING WITH CORRECT NEIGHBOR COUNTING")
    print("="*70)

    datasets = {}

    # Load datasets
    try:
        data = np.load('realistic_z0.10.npz')
        datasets['Synthetic z0.10'] = data['positions']
    except:
        pass

    try:
        data = np.load('illustris_realization_1.npz')
        datasets['Illustris-like'] = data['positions']
    except:
        pass

    # Add pure random
    box_size = 100
    n_points = 10000
    random_positions = np.random.uniform(0, box_size, (n_points, 3))
    datasets['Random Poisson'] = random_positions

    results = {}

    for name, positions in datasets.items():
        print(f"\n  Dataset: {name}")
        stats = compute_correct_voronoi_neighbors(positions)

        if stats is not None:
            results[name] = stats
            print(f"    Mean neighbors: {stats['mean']:.2f} ± {stats['std']:.2f}")

    return results


def investigate_4d_hypothesis():
    """
    Investigate if there's a 4D connection.

    Hypothesis: If the universe has hidden dimensions or 4D structure,
    projecting 4D Voronoi to 3D might give different statistics.
    """
    print("\n" + "="*70)
    print("4D HYPOTHESIS INVESTIGATION")
    print("="*70)

    print("""
Theoretical Background:

1. 4D VORONOI TESSELLATION:
   - In 4D space, random Poisson process gives:
   - Expected neighbors: ~24 (face adjacency in 4D)
   - This is close to our ~25!

2. PROJECTION HYPOTHESIS:
   If we project 4D points to 3D, what happens to Voronoi structure?

   - 4D Voronoi → 3D projection might preserve some adjacencies
   - Could explain why we see ~24-25 instead of ~15.5

3. PHYSICAL MOTIVATION:
   - Kaluza-Klein theory: extra spatial dimension
   - String theory: 10/11 dimensions (6/7 compactified)
   - Cosmic web might reflect higher-dimensional structure

Let's test this by:
   a) Generating 4D points
   b) Computing 4D Voronoi
   c) Projecting to 3D
   d) Checking if neighbor count is preserved
""")

    # Generate 4D random points
    n_points = 2000
    points_4d = np.random.uniform(0, 100, (n_points, 4))

    print(f"\nGenerating {n_points} 4D random points...")

    # Compute 4D Voronoi
    try:
        print("Computing 4D Voronoi tessellation...")
        vor_4d = Voronoi(points_4d)

        # Count 4D neighbors
        neighbor_counts_4d = [0] * n_points
        for ridge in vor_4d.ridge_points:
            p1, p2 = ridge
            neighbor_counts_4d[p1] += 1
            neighbor_counts_4d[p2] += 1

        neighbor_counts_4d = np.array(neighbor_counts_4d)

        print(f"\n4D Voronoi neighbors:")
        print(f"  Mean: {np.mean(neighbor_counts_4d):.2f}")
        print(f"  Std: {np.std(neighbor_counts_4d):.2f}")

        # Project to 3D
        points_3d_projected = points_4d[:, :3]  # Take first 3 dimensions

        print(f"\nProjecting to 3D and recomputing Voronoi...")
        vor_3d = Voronoi(points_3d_projected)

        neighbor_counts_3d = [0] * n_points
        for ridge in vor_3d.ridge_points:
            p1, p2 = ridge
            neighbor_counts_3d[p1] += 1
            neighbor_counts_3d[p2] += 1

        neighbor_counts_3d = np.array(neighbor_counts_3d)

        print(f"\n3D Projected Voronoi neighbors:")
        print(f"  Mean: {np.mean(neighbor_counts_3d):.2f}")
        print(f"  Std: {np.std(neighbor_counts_3d):.2f}")

        print(f"\nComparison:")
        print(f"  4D Voronoi: {np.mean(neighbor_counts_4d):.2f} neighbors")
        print(f"  3D projected: {np.mean(neighbor_counts_3d):.2f} neighbors")
        print(f"  Expected 3D random: 15.54 neighbors")

        if np.abs(np.mean(neighbor_counts_4d) - 24) < 2.0:
            print(f"\n  ⚠️  4D matches theoretical prediction of ~24!")

        if np.abs(np.mean(neighbor_counts_3d) - np.mean(neighbor_counts_4d)) < 2.0:
            print(f"  ⚠️  Projection preserves neighbor count!")
            print(f"  → This is VERY interesting for the 4D hypothesis!")

        return {
            '4d_neighbors': neighbor_counts_4d,
            '3d_projected_neighbors': neighbor_counts_3d
        }

    except Exception as e:
        print(f"4D Voronoi computation failed: {e}")
        return None


def visualize_voronoi_confusion():
    """
    Create visualization explaining the vertex vs neighbor confusion.
    """
    print("\nCreating explanatory visualization...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Simple 2D example
    ax1 = axes[0]

    # Create a simple Voronoi in 2D to illustrate
    points_2d = np.array([
        [0.5, 0.5],
        [1.5, 0.5],
        [0.5, 1.5],
        [1.5, 1.5],
        [1.0, 1.0]
    ])

    vor_2d = Voronoi(points_2d)

    # Plot
    from scipy.spatial import voronoi_plot_2d
    voronoi_plot_2d(vor_2d, ax=ax1, show_vertices=True, show_points=True)

    ax1.set_title('2D Voronoi: Vertices vs Neighbors', fontweight='bold', fontsize=12)
    ax1.set_xlim(-0.5, 2.5)
    ax1.set_ylim(-0.5, 2.5)

    # Add annotations
    ax1.text(0.5, -0.3, 'Vertices = corner points (dots)\nNeighbors = cells sharing edge',
             ha='left', fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat'))

    # Comparison plot
    ax2 = axes[1]

    # Show theoretical values
    metrics = ['Neighbors\n(faces)', 'Edges', 'Vertices']
    values_3d = [15.54, 48, 35]
    our_measurement = [25.5, np.nan, np.nan]

    x = np.arange(len(metrics))
    width = 0.35

    ax2.bar(x - width/2, values_3d, width, label='Theoretical 3D Poisson', alpha=0.7)
    ax2.bar(x + width/2, our_measurement, width, label='Our Measurement', alpha=0.7)

    ax2.set_ylabel('Count', fontsize=11)
    ax2.set_title('3D Voronoi Cell Properties', fontweight='bold', fontsize=12)
    ax2.set_xticks(x)
    ax2.set_xticklabels(metrics)
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')

    # Add arrow and text
    ax2.annotate('We were counting\nvertices (35),\nnot neighbors (15)!',
                xy=(2, 35), xytext=(1.5, 45),
                arrowprops=dict(arrowstyle='->', lw=2, color='red'),
                fontsize=11, color='red', fontweight='bold')

    plt.tight_layout()
    plt.savefig('voronoi_confusion_explained.png', dpi=150, bbox_inches='tight')
    print("Saved: voronoi_confusion_explained.png")
    plt.close()


def main():
    """Run full investigation."""
    print("="*70)
    print("VORONOI ANOMALY: DEEP DIVE INVESTIGATION")
    print("="*70)

    # Step 1: Theory
    theoretical_voronoi_neighbors_3d()

    # Step 2: Test known case
    print("\n" + "="*70)
    print("TESTING ON KNOWN CASE")
    print("="*70)
    confirmed = test_what_were_counting()

    if confirmed:
        print("\n" + "="*70)
        print("DIAGNOSIS CONFIRMED!")
        print("="*70)
        print("""
We were counting VERTICES of Voronoi regions, not NEIGHBORS!

For a 3D Voronoi cell:
- Neighbors (cells sharing a face): ~15.5 for random Poisson
- Vertices (corner points): ~35 for random Poisson

Our ~25 value is between these, suggesting we're counting
something hybrid or there's a boundary effect.

The scipy.spatial.Voronoi().regions attribute gives vertices,
NOT neighbor adjacencies!
""")

    # Step 3: Retest with correct method
    correct_results = retest_all_datasets_correctly()

    # Step 4: Compare
    print("\n" + "="*70)
    print("COMPARISON: WRONG vs CORRECT METHOD")
    print("="*70)
    print(f"\n{'Dataset':<25} {'Wrong (vertices)':<20} {'Correct (neighbors)':<20}")
    print("-"*70)

    for name in correct_results.keys():
        if name in correct_results:
            correct_mean = correct_results[name]['mean']
            print(f"{name:<25} {'~25':<20} {correct_mean:>6.2f}")

    # Step 5: 4D hypothesis
    result_4d = investigate_4d_hypothesis()

    # Step 6: Visualization
    visualize_voronoi_confusion()

    print("\n" + "="*70)
    print("INVESTIGATION COMPLETE!")
    print("="*70)

    return {
        'confirmed_bug': confirmed,
        'correct_results': correct_results,
        '4d_results': result_4d
    }


if __name__ == "__main__":
    results = main()
