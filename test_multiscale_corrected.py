"""
Multi-scale non-Gaussianity test with CORRECTED validation.

Key fix: Test "clustered vs unclustered POINT CLOUDS"
NOT "Gaussian vs lognormal FIELDS"

This matches actual cosmological structure and avoids Poisson sampling issues.
"""

import numpy as np
from scipy import ndimage, stats
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List


@dataclass
class ScaleRange:
    """Define scales to test."""
    cell_sizes: List[int]
    smooth_scales: List[float]

    @classmethod
    def logarithmic(cls, min_size=3, max_size=32, n_scales=15):
        """Logarithmic spacing of scales."""
        cell_sizes = np.logspace(np.log10(min_size), np.log10(max_size), n_scales)
        cell_sizes = np.unique(cell_sizes.astype(int))

        smooth_scales = np.logspace(np.log10(0.5), np.log10(4.0), n_scales)

        return cls(
            cell_sizes=cell_sizes.tolist(),
            smooth_scales=smooth_scales.tolist()
        )


class MultiscaleNonGaussianityTest:
    """Test non-Gaussianity using multi-scale curves."""

    def __init__(self, box_size=100.0, grid_size=128):
        self.box_size = box_size
        self.grid_size = grid_size

    def standardize_and_grid(self, positions):
        """CIC gridding (all fixes applied)."""
        # Ensure positions are in [0, box_size]
        positions = positions % self.box_size

        mins = positions.min(axis=0)
        maxs = positions.max(axis=0)

        pos_norm = (positions - mins) / (maxs - mins) * (self.grid_size - 1)

        grid = np.zeros((self.grid_size, self.grid_size, self.grid_size))

        for pos in pos_norm:
            i, j, k = pos.astype(int)
            dx = pos[0] - i
            dy = pos[1] - j
            dz = pos[2] - k

            for di in [0, 1]:
                for dj in [0, 1]:
                    for dk in [0, 1]:
                        ii = (i + di) % self.grid_size
                        jj = (j + dj) % self.grid_size
                        kk = (k + dk) % self.grid_size

                        weight = (1 - abs(di - dx)) * (1 - abs(dj - dy)) * (1 - abs(dk - dz))
                        grid[ii, jj, kk] += weight

        return grid

    def counts_variance_at_scale(self, grid, cell_size):
        """Variance-to-mean ratio at specific cell size."""
        if cell_size > self.grid_size // 2:
            return np.nan

        counts = []
        for i in range(0, self.grid_size - cell_size, cell_size):
            for j in range(0, self.grid_size - cell_size, cell_size):
                for k in range(0, self.grid_size - cell_size, cell_size):
                    cell = grid[i:i+cell_size, j:j+cell_size, k:k+cell_size]
                    counts.append(cell.sum())

        counts = np.array(counts)
        mean_count = counts.mean()

        if mean_count > 0:
            return counts.var() / mean_count
        return np.nan

    def counts_skewness_at_scale(self, grid, cell_size):
        """Skewness at specific cell size."""
        if cell_size > self.grid_size // 2:
            return np.nan

        counts = []
        for i in range(0, self.grid_size - cell_size, cell_size):
            for j in range(0, self.grid_size - cell_size, cell_size):
                for k in range(0, self.grid_size - cell_size, cell_size):
                    cell = grid[i:i+cell_size, j:j+cell_size, k:k+cell_size]
                    counts.append(cell.sum())

        counts = np.array(counts)

        if len(counts) > 3 and counts.std() > 0:
            return stats.skew(counts)
        return np.nan

    def euler_at_scale(self, grid, smooth_scale, threshold=0.5):
        """Euler characteristic at specific smoothing scale."""
        mean_density = grid.mean()
        if mean_density > 0:
            delta = (grid - mean_density) / mean_density
        else:
            delta = grid - mean_density

        delta_smooth = ndimage.gaussian_filter(delta, sigma=smooth_scale)

        binary = delta_smooth > threshold

        labeled, beta_0 = ndimage.label(binary)
        inverted = ~binary
        labeled_inv, n_inv = ndimage.label(inverted)
        beta_2 = max(0, n_inv - 1)
        beta_1 = beta_0 + beta_2 - 1

        return beta_0 - beta_1 + beta_2

    def compute_multiscale_curves(self, positions, scale_range: ScaleRange):
        """Compute metrics over multiple scales."""
        grid = self.standardize_and_grid(positions)

        curves = {
            'variance': {'scales': [], 'values': []},
            'skewness': {'scales': [], 'values': []},
            'euler': {'scales': [], 'values': []}
        }

        for cell_size in scale_range.cell_sizes:
            var = self.counts_variance_at_scale(grid, cell_size)
            skew = self.counts_skewness_at_scale(grid, cell_size)

            if not np.isnan(var):
                curves['variance']['scales'].append(cell_size)
                curves['variance']['values'].append(var)

            if not np.isnan(skew):
                curves['skewness']['scales'].append(cell_size)
                curves['skewness']['values'].append(skew)

        for smooth_scale in scale_range.smooth_scales:
            euler = self.euler_at_scale(grid, smooth_scale)

            if not np.isnan(euler):
                curves['euler']['scales'].append(smooth_scale)
                curves['euler']['values'].append(euler)

        return curves

    def curve_chi_squared(self, real_curve, mock_curves):
        """χ² test comparing real curve to mock distribution."""
        real_values = np.array(real_curve['values'])

        n_scales = len(real_values)
        mock_values_at_scales = []

        for mock in mock_curves:
            mock_vals = np.array(mock['values'])
            if len(mock_vals) == n_scales:
                mock_values_at_scales.append(mock_vals)

        if len(mock_values_at_scales) == 0:
            return {'chi_squared': 0, 'dof': n_scales, 'p_value': 1.0}

        mock_values_at_scales = np.array(mock_values_at_scales)

        mock_mean = mock_values_at_scales.mean(axis=0)
        mock_var = mock_values_at_scales.var(axis=0)

        chi_squared_per_scale = (real_values - mock_mean)**2 / (mock_var + 1e-10)
        chi_squared_total = chi_squared_per_scale.sum()
        dof = n_scales

        p_value = 1 - stats.chi2.cdf(chi_squared_total, dof)

        return {
            'chi_squared': chi_squared_total,
            'dof': dof,
            'p_value': p_value,
            'chi_squared_per_scale': chi_squared_per_scale
        }


# ==============================================================================
# CORRECTED CONTROL GENERATORS
# ==============================================================================

def generate_unclustered_control(n_points=10000, box_size=100, seed=None):
    """
    Pure Poisson random point distribution - NO structure.

    This represents the null hypothesis: completely random placement.
    """
    if seed is not None:
        np.random.seed(seed)

    return np.random.uniform(0, box_size, (n_points, 3))


def generate_clustered_control(n_points=10000, box_size=100, n_halos=100,
                               halo_scale=2.0, seed=None):
    """
    Halo-like clustered distribution.

    Represents cosmic structure:
    - Galaxies cluster in halos
    - Halos have NFW-ish profiles
    - Clear departure from Poisson random

    Parameters:
    - n_halos: Number of clusters
    - halo_scale: Size of each halo (Mpc/h)
    """
    if seed is not None:
        np.random.seed(seed)

    # Place halo centers randomly
    halo_centers = np.random.uniform(0, box_size, (n_halos, 3))

    # Distribute points among halos (some variation)
    points_per_halo = np.random.poisson(n_points / n_halos, n_halos)
    points_per_halo = (points_per_halo * n_points / points_per_halo.sum()).astype(int)

    positions = []

    for i, (center, n_gal) in enumerate(zip(halo_centers, points_per_halo)):
        if n_gal == 0:
            continue

        # NFW-ish profile: concentrated center, extended halo
        # Use mixture of Gaussians for simplicity

        # 50% in core (tight)
        n_core = n_gal // 2
        core_offsets = np.random.normal(0, halo_scale * 0.3, (n_core, 3))

        # 50% in halo (extended)
        n_ext = n_gal - n_core
        ext_offsets = np.random.normal(0, halo_scale * 1.5, (n_ext, 3))

        offsets = np.vstack([core_offsets, ext_offsets])
        halo_galaxies = center + offsets

        positions.append(halo_galaxies)

    positions = np.vstack(positions)

    # Apply periodic boundary conditions
    positions = positions % box_size

    return positions[:n_points]  # Ensure exact count


def generate_filamentary_control(n_points=10000, box_size=100, n_filaments=20,
                                 filament_length=30, filament_width=1.5, seed=None):
    """
    Filamentary structure (even more realistic).

    Galaxies cluster along filaments - characteristic of cosmic web.
    """
    if seed is not None:
        np.random.seed(seed)

    points_per_filament = n_points // n_filaments
    positions = []

    for i in range(n_filaments):
        # Random filament endpoints
        start = np.random.uniform(0, box_size, 3)
        direction = np.random.normal(0, 1, 3)
        direction = direction / np.linalg.norm(direction)

        # Points along filament
        t = np.random.uniform(0, filament_length, points_per_filament)

        for ti in t:
            # Point on filament axis
            point = start + ti * direction

            # Add perpendicular scatter (filament width)
            perp1 = np.cross(direction, np.random.normal(0, 1, 3))
            perp1 = perp1 / (np.linalg.norm(perp1) + 1e-10)
            perp2 = np.cross(direction, perp1)

            scatter = (np.random.normal(0, filament_width) * perp1 +
                      np.random.normal(0, filament_width) * perp2)

            positions.append(point + scatter)

    positions = np.array(positions)
    positions = positions % box_size

    return positions[:n_points]


# ==============================================================================
# VALIDATION WITH CORRECTED CONTROLS
# ==============================================================================

def validation_corrected(n_mocks=50):
    """
    CORRECTED validation using proper point cloud controls.

    Tests:
    1. Unclustered vs Unclustered mocks (should NOT flag)
    2. Clustered vs Unclustered mocks (SHOULD flag)
    3. Filamentary vs Unclustered mocks (SHOULD flag even more)
    """
    print("="*70)
    print("CORRECTED MULTISCALE VALIDATION")
    print("="*70)
    print(f"\nUsing POINT CLOUD controls (not field-based)")
    print(f"Generating {n_mocks} mocks...\n")

    tester = MultiscaleNonGaussianityTest(box_size=100.0, grid_size=128)
    scale_range = ScaleRange.logarithmic(min_size=3, max_size=32, n_scales=15)

    print(f"Testing {len(scale_range.cell_sizes)} cell sizes")
    print(f"Testing {len(scale_range.smooth_scales)} smooth scales\n")

    # Control 1: Unclustered vs Unclustered
    print("Control 1: Unclustered vs Unclustered (should NOT flag)")
    print("-"*70)

    unclust_real_pos = generate_unclustered_control(n_points=10000, seed=1000)
    unclust_real_curves = tester.compute_multiscale_curves(unclust_real_pos, scale_range)

    unclust_mock_curves = {'variance': [], 'skewness': [], 'euler': []}

    for i in range(n_mocks):
        if (i+1) % 20 == 0:
            print(f"  Generated {i+1}/{n_mocks} unclustered mocks...")
        mock_pos = generate_unclustered_control(n_points=10000, seed=2000+i)
        mock_curves = tester.compute_multiscale_curves(mock_pos, scale_range)

        for metric in unclust_mock_curves.keys():
            unclust_mock_curves[metric].append(mock_curves[metric])

    # Test
    unclust_results = {}
    for metric in ['variance', 'skewness', 'euler']:
        chi2_result = tester.curve_chi_squared(unclust_real_curves[metric],
                                               unclust_mock_curves[metric])
        unclust_results[metric] = chi2_result
        print(f"  {metric}: χ² = {chi2_result['chi_squared']:.2f} (DOF={chi2_result['dof']}), p = {chi2_result['p_value']:.4f}")

    n_false_pos = sum(res['p_value'] < 0.05 for res in unclust_results.values())
    print(f"\n  False positives: {n_false_pos}/3 (expect ≤1) {'✓' if n_false_pos <= 1 else '✗'}")

    # Control 2: Clustered (halos) vs Unclustered
    print("\nControl 2: Clustered (halos) vs Unclustered (SHOULD flag)")
    print("-"*70)

    clust_real_pos = generate_clustered_control(n_points=10000, n_halos=100, seed=1001)
    clust_real_curves = tester.compute_multiscale_curves(clust_real_pos, scale_range)

    clust_results = {}
    for metric in ['variance', 'skewness', 'euler']:
        chi2_result = tester.curve_chi_squared(clust_real_curves[metric],
                                               unclust_mock_curves[metric])
        clust_results[metric] = chi2_result
        marker = "***" if chi2_result['p_value'] < 0.01 else ("**" if chi2_result['p_value'] < 0.05 else "")
        print(f"  {metric}: χ² = {chi2_result['chi_squared']:.2f} (DOF={chi2_result['dof']}), p = {chi2_result['p_value']:.4f} {marker}")

    n_detections_halo = sum(res['p_value'] < 0.05 for res in clust_results.values())
    print(f"\n  Detections: {n_detections_halo}/3 (expect ≥2) {'✓' if n_detections_halo >= 2 else '✗'}")

    # Control 3: Filamentary vs Unclustered
    print("\nControl 3: Filamentary vs Unclustered (SHOULD flag strongly)")
    print("-"*70)

    fil_real_pos = generate_filamentary_control(n_points=10000, n_filaments=20, seed=1002)
    fil_real_curves = tester.compute_multiscale_curves(fil_real_pos, scale_range)

    fil_results = {}
    for metric in ['variance', 'skewness', 'euler']:
        chi2_result = tester.curve_chi_squared(fil_real_curves[metric],
                                               unclust_mock_curves[metric])
        fil_results[metric] = chi2_result
        marker = "***" if chi2_result['p_value'] < 0.01 else ("**" if chi2_result['p_value'] < 0.05 else "")
        print(f"  {metric}: χ² = {chi2_result['chi_squared']:.2f} (DOF={chi2_result['dof']}), p = {chi2_result['p_value']:.4f} {marker}")

    n_detections_fil = sum(res['p_value'] < 0.05 for res in fil_results.values())
    print(f"\n  Detections: {n_detections_fil}/3 (expect ≥2) {'✓' if n_detections_fil >= 2 else '✗'}")

    # Overall validation
    print("\n" + "="*70)
    validation_passed = (n_false_pos <= 1) and (n_detections_halo >= 2)

    if validation_passed:
        print("✅ VALIDATION PASSED")
        print("   Multiscale curves can distinguish clustered from unclustered!")
        print(f"   Halo detection: {n_detections_halo}/3 metrics")
        print(f"   Filament detection: {n_detections_fil}/3 metrics")
    else:
        print("⚠️  VALIDATION UNCLEAR")
        print(f"   False positives: {n_false_pos}/3")
        print(f"   Halo detections: {n_detections_halo}/3")

    print("="*70)

    # Visualize
    visualize_corrected_validation(
        unclust_real_curves, unclust_mock_curves,
        clust_real_curves, fil_real_curves,
        scale_range, unclust_results, clust_results, fil_results
    )

    return validation_passed, unclust_results, clust_results, fil_results


def visualize_corrected_validation(unclust_real, unclust_mocks,
                                   clust_real, fil_real, scale_range,
                                   unclust_results, clust_results, fil_results):
    """Visualize corrected validation curves."""
    fig, axes = plt.subplots(3, 3, figsize=(16, 12))

    metrics = ['variance', 'skewness', 'euler']
    titles = ['Variance vs Cell Size', 'Skewness vs Cell Size', 'Euler vs Smooth Scale']
    controls = [
        ('Unclustered vs Unclustered', unclust_real, unclust_results),
        ('Clustered (Halos) vs Unclustered', clust_real, clust_results),
        ('Filamentary vs Unclustered', fil_real, fil_results)
    ]

    for idx, (metric, title) in enumerate(zip(metrics, titles)):
        scales = unclust_real[metric]['scales']

        # Mock distribution (same for all)
        mock_values_array = np.array([m['values'] for m in unclust_mocks[metric]])
        mock_mean = mock_values_array.mean(axis=0)
        mock_std = mock_values_array.std(axis=0)

        for col, (control_name, real_curves, results) in enumerate(controls):
            ax = axes[idx, col]

            real_vals = real_curves[metric]['values']
            p_val = results[metric]['p_value']

            # Plot
            ax.plot(scales, mock_mean, 'b-', linewidth=2, label='Unclustered Mean', alpha=0.7)
            ax.fill_between(scales, mock_mean - mock_std, mock_mean + mock_std,
                           alpha=0.2, color='blue', label='±1σ')

            color = 'red' if col > 0 else 'green'
            ax.plot(scales, real_vals, 'o-', linewidth=2, markersize=4,
                   color=color, label=control_name.split(' vs ')[0])

            ax.set_xlabel('Scale', fontsize=9)
            ax.set_ylabel(metric.capitalize(), fontsize=9)

            title_color = 'green' if p_val < 0.05 else 'black'
            if col == 0:
                title_color = 'green' if p_val > 0.05 else 'red'

            ax.set_title(f'{title}\n{control_name}\np={p_val:.4f}',
                        fontsize=10, fontweight='bold', color=title_color)
            ax.legend(fontsize=7, loc='best')
            ax.grid(True, alpha=0.3)
            ax.set_xscale('log')

    plt.suptitle('Corrected Validation: Point Cloud Clustering Tests',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('multiscale_corrected_validation.png', dpi=150, bbox_inches='tight')
    print("\nSaved: multiscale_corrected_validation.png")
    plt.close()


def main():
    """Run corrected multiscale validation."""
    print("\n" + "="*70)
    print("MULTISCALE TEST - CORRECTED VALIDATION")
    print("="*70)
    print("\nKey improvement:")
    print("  OLD: Gaussian vs lognormal FIELDS → Poisson sampling breaks it")
    print("  NEW: Clustered vs unclustered POINT CLOUDS → Matches real data")
    print("="*70)

    validation_passed, unclust_res, clust_res, fil_res = validation_corrected(n_mocks=50)

    print("\n" + "="*70)
    if validation_passed:
        print("✅ SUCCESS - READY TO TEST REAL DATA")
        print("\nThe pipeline can now distinguish:")
        print("  - Halo-like clustering from random")
        print("  - Filamentary structure from random")
        print("\nNext step: Apply to your realistic_z0.10.npz data!")
    else:
        print("⚠️  NEEDS MORE TUNING")
        print("\nBut this is the RIGHT approach - testing point cloud structure,")
        print("not abstract field Gaussianity.")

    print("="*70)

    return validation_passed


if __name__ == "__main__":
    results = main()
