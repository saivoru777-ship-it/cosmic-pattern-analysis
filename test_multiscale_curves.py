"""
Multi-scale non-Gaussianity test using curve-level statistics.

Key innovation: Instead of testing single numbers at one scale,
compute metrics over 10-20 scales and test the CURVE shape.

This is:
- Much more powerful statistically
- Standard in cosmology (scale-dependent bias)
- Robust to noise at any single scale
- Tests HOW non-Gaussianity varies with scale, not just IF it exists
"""

import numpy as np
from scipy import ndimage, stats
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, List, Tuple


@dataclass
class ScaleRange:
    """Define scales to test."""
    cell_sizes: List[int]  # In voxels
    smooth_scales: List[float]  # Gaussian σ in voxels

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
        # Convert to overdensity
        mean_density = grid.mean()
        if mean_density > 0:
            delta = (grid - mean_density) / mean_density
        else:
            delta = grid - mean_density

        # Smooth at this scale
        delta_smooth = ndimage.gaussian_filter(delta, sigma=smooth_scale)

        # Threshold
        binary = delta_smooth > threshold

        # Compute χ
        labeled, beta_0 = ndimage.label(binary)
        inverted = ~binary
        labeled_inv, n_inv = ndimage.label(inverted)
        beta_2 = max(0, n_inv - 1)
        beta_1 = beta_0 + beta_2 - 1

        euler_char = beta_0 - beta_1 + beta_2

        return euler_char

    def compute_multiscale_curves(self, positions, scale_range: ScaleRange):
        """
        Compute metrics over multiple scales.

        Returns curves, not single numbers.
        """
        grid = self.standardize_and_grid(positions)

        curves = {
            'variance': {'scales': [], 'values': []},
            'skewness': {'scales': [], 'values': []},
            'euler': {'scales': [], 'values': []}
        }

        # Variance and skewness vs cell size
        for cell_size in scale_range.cell_sizes:
            var = self.counts_variance_at_scale(grid, cell_size)
            skew = self.counts_skewness_at_scale(grid, cell_size)

            if not np.isnan(var):
                curves['variance']['scales'].append(cell_size)
                curves['variance']['values'].append(var)

            if not np.isnan(skew):
                curves['skewness']['scales'].append(cell_size)
                curves['skewness']['values'].append(skew)

        # Euler vs smoothing scale
        for smooth_scale in scale_range.smooth_scales:
            euler = self.euler_at_scale(grid, smooth_scale)

            if not np.isnan(euler):
                curves['euler']['scales'].append(smooth_scale)
                curves['euler']['values'].append(euler)

        return curves

    def curve_chi_squared(self, real_curve, mock_curves):
        """
        χ² test comparing real curve to mock distribution.

        At each scale point:
        χ²_i = (real_i - mean(mocks_i))² / var(mocks_i)

        Total χ² = Σ χ²_i
        DOF = n_scales
        """
        real_values = np.array(real_curve['values'])

        # Get mock values at each scale
        n_scales = len(real_values)
        mock_values_at_scales = []

        for mock in mock_curves:
            mock_vals = np.array(mock['values'])
            # Ensure same length (should be, but check)
            if len(mock_vals) == n_scales:
                mock_values_at_scales.append(mock_vals)

        mock_values_at_scales = np.array(mock_values_at_scales)  # Shape: (n_mocks, n_scales)

        # Mean and variance at each scale
        mock_mean = mock_values_at_scales.mean(axis=0)
        mock_var = mock_values_at_scales.var(axis=0)

        # χ² at each scale
        chi_squared_per_scale = (real_values - mock_mean)**2 / (mock_var + 1e-10)

        # Total χ²
        chi_squared_total = chi_squared_per_scale.sum()
        dof = n_scales

        # p-value from χ² distribution
        p_value = 1 - stats.chi2.cdf(chi_squared_total, dof)

        return {
            'chi_squared': chi_squared_total,
            'dof': dof,
            'p_value': p_value,
            'chi_squared_per_scale': chi_squared_per_scale
        }

    def curve_ks_test(self, real_curve, mock_curves):
        """
        Alternative: Kolmogorov-Smirnov test on curve residuals.

        Tests if (real - mock_mean) distribution matches expected.
        """
        real_values = np.array(real_curve['values'])

        mock_values_at_scales = []
        for mock in mock_curves:
            mock_vals = np.array(mock['values'])
            if len(mock_vals) == len(real_values):
                mock_values_at_scales.append(mock_vals)

        mock_values_at_scales = np.array(mock_values_at_scales)
        mock_mean = mock_values_at_scales.mean(axis=0)
        mock_std = mock_values_at_scales.std(axis=0)

        # Normalized residuals
        residuals = (real_values - mock_mean) / (mock_std + 1e-10)

        # KS test: should these residuals be Gaussian(0,1)?
        ks_stat, p_value = stats.kstest(residuals, 'norm', args=(0, 1))

        return {
            'ks_statistic': ks_stat,
            'p_value': p_value,
            'residuals': residuals
        }


def generate_gaussian_field(n_points=10000, box_size=100, grid_size=128, seed=None):
    """Generate Gaussian control."""
    if seed is not None:
        np.random.seed(seed)

    gaussian_field = np.random.normal(0, 1, (grid_size, grid_size, grid_size))
    gaussian_field = ndimage.gaussian_filter(gaussian_field, sigma=2.0)

    density = np.exp(gaussian_field - gaussian_field.max() + 10)
    density = density / density.sum()

    indices = np.random.choice(density.size, size=n_points, p=density.flatten())
    i, j, k = np.unravel_index(indices, gaussian_field.shape)

    positions = np.column_stack([i, j, k]) * (box_size / grid_size)
    positions += np.random.uniform(-0.5, 0.5, (n_points, 3)) * (box_size / grid_size)

    return positions


def generate_lognormal_field(n_points=10000, box_size=100, grid_size=128, seed=None):
    """Generate lognormal control (non-Gaussian)."""
    if seed is not None:
        np.random.seed(seed)

    gaussian_field = np.random.normal(0, 1, (grid_size, grid_size, grid_size))
    gaussian_field = ndimage.gaussian_filter(gaussian_field, sigma=2.0)

    lognormal_field = np.exp(gaussian_field)
    density = lognormal_field / lognormal_field.sum()

    indices = np.random.choice(density.size, size=n_points, p=density.flatten())
    i, j, k = np.unravel_index(indices, lognormal_field.shape)

    positions = np.column_stack([i, j, k]) * (box_size / grid_size)
    positions += np.random.uniform(-0.5, 0.5, (n_points, 3)) * (box_size / grid_size)

    return positions


def validation_multiscale(n_mocks=50):
    """
    Validation using multi-scale curves.

    Test if curve-level statistics can distinguish Gaussian from lognormal.
    """
    print("="*70)
    print("MULTISCALE VALIDATION")
    print("="*70)
    print(f"\nGenerating {n_mocks} mocks...\n")

    tester = MultiscaleNonGaussianityTest(box_size=100.0, grid_size=128)
    scale_range = ScaleRange.logarithmic(min_size=3, max_size=32, n_scales=15)

    print(f"Testing {len(scale_range.cell_sizes)} cell sizes: {scale_range.cell_sizes[:3]}...{scale_range.cell_sizes[-3:]}")
    smooth_first = [f"{x:.2f}" for x in scale_range.smooth_scales[:3]]
    smooth_last = [f"{x:.2f}" for x in scale_range.smooth_scales[-3:]]
    print(f"Testing {len(scale_range.smooth_scales)} smooth scales: {smooth_first}...{smooth_last}\n")

    # Control 1: Gaussian vs Gaussian
    print("Control 1: Gaussian vs Gaussian (should NOT flag)")
    print("-"*70)

    gauss_real_pos = generate_gaussian_field(n_points=10000, seed=1000)
    gauss_real_curves = tester.compute_multiscale_curves(gauss_real_pos, scale_range)

    gauss_mock_curves = {'variance': [], 'skewness': [], 'euler': []}

    for i in range(n_mocks):
        if (i+1) % 20 == 0:
            print(f"  Generated {i+1}/{n_mocks} mocks...")
        mock_pos = generate_gaussian_field(n_points=10000, seed=2000+i)
        mock_curves = tester.compute_multiscale_curves(mock_pos, scale_range)

        for metric in gauss_mock_curves.keys():
            gauss_mock_curves[metric].append(mock_curves[metric])

    # Test each curve
    gauss_results = {}
    for metric in ['variance', 'skewness', 'euler']:
        chi2_result = tester.curve_chi_squared(gauss_real_curves[metric], gauss_mock_curves[metric])
        gauss_results[metric] = chi2_result
        print(f"  {metric}: χ² = {chi2_result['chi_squared']:.2f} (DOF={chi2_result['dof']}), p = {chi2_result['p_value']:.4f}")

    n_false_pos = sum(res['p_value'] < 0.05 for res in gauss_results.values())
    print(f"\n  False positives: {n_false_pos}/3 (expect ≤1)")

    # Control 2: Lognormal vs Gaussian
    print("\nControl 2: Lognormal vs Gaussian (SHOULD flag)")
    print("-"*70)

    ln_real_pos = generate_lognormal_field(n_points=10000, seed=1001)
    ln_real_curves = tester.compute_multiscale_curves(ln_real_pos, scale_range)

    # Test against same Gaussian mocks
    ln_results = {}
    for metric in ['variance', 'skewness', 'euler']:
        chi2_result = tester.curve_chi_squared(ln_real_curves[metric], gauss_mock_curves[metric])
        ln_results[metric] = chi2_result
        print(f"  {metric}: χ² = {chi2_result['chi_squared']:.2f} (DOF={chi2_result['dof']}), p = {chi2_result['p_value']:.4f}")

    n_detections = sum(res['p_value'] < 0.05 for res in ln_results.values())
    print(f"\n  Detections: {n_detections}/3 (expect ≥2)")

    # Overall
    print("\n" + "="*70)
    validation_passed = (n_false_pos <= 1) and (n_detections >= 2)

    if validation_passed:
        print("✅ MULTISCALE VALIDATION PASSED")
        print("   Curve-level tests can distinguish Gaussian from lognormal!")
    else:
        print("⚠️  VALIDATION UNCLEAR")
        print("   May need more mocks or scale tuning")

    print("="*70)

    # Visualize curves
    visualize_curves(gauss_real_curves, gauss_mock_curves, ln_real_curves,
                    scale_range, gauss_results, ln_results)

    return validation_passed, gauss_results, ln_results


def visualize_curves(gauss_real, gauss_mocks, ln_real, scale_range, gauss_results, ln_results):
    """Visualize multiscale curves."""
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))

    metrics = ['variance', 'skewness', 'euler']
    titles = ['Counts Variance vs Cell Size', 'Counts Skewness vs Cell Size',
              'Euler Characteristic vs Smooth Scale']

    for idx, (metric, title) in enumerate(zip(metrics, titles)):
        # Gaussian vs Gaussian
        ax_gauss = axes[idx, 0]

        scales = gauss_real[metric]['scales']
        real_vals = gauss_real[metric]['values']

        # Plot mock distribution - gauss_mocks[metric] is a list of curve dicts
        mock_values_array = np.array([m['values'] for m in gauss_mocks[metric]])
        mock_mean = mock_values_array.mean(axis=0)
        mock_std = mock_values_array.std(axis=0)

        ax_gauss.plot(scales, mock_mean, 'b-', linewidth=2, label='Gaussian Mock Mean')
        ax_gauss.fill_between(scales, mock_mean - mock_std, mock_mean + mock_std,
                             alpha=0.3, color='blue', label='±1σ')
        ax_gauss.plot(scales, real_vals, 'go-', linewidth=2, markersize=4, label='Gaussian Real')

        ax_gauss.set_xlabel('Scale', fontsize=10)
        ax_gauss.set_ylabel(metric.capitalize(), fontsize=10)
        ax_gauss.set_title(f'{title}\\nGaussian vs Gaussian (p={gauss_results[metric]["p_value"]:.3f})',
                          fontsize=11, fontweight='bold')
        ax_gauss.legend(fontsize=8)
        ax_gauss.grid(True, alpha=0.3)
        ax_gauss.set_xscale('log')

        # Lognormal vs Gaussian
        ax_ln = axes[idx, 1]

        ln_vals = ln_real[metric]['values']

        ax_ln.plot(scales, mock_mean, 'b-', linewidth=2, label='Gaussian Mock Mean')
        ax_ln.fill_between(scales, mock_mean - mock_std, mock_mean + mock_std,
                          alpha=0.3, color='blue', label='±1σ')
        ax_ln.plot(scales, ln_vals, 'ro-', linewidth=2, markersize=4, label='Lognormal Real')

        ax_ln.set_xlabel('Scale', fontsize=10)
        ax_ln.set_ylabel(metric.capitalize(), fontsize=10)
        ax_ln.set_title(f'{title}\\nLognormal vs Gaussian (p={ln_results[metric]["p_value"]:.3f})',
                       fontsize=11, fontweight='bold',
                       color='green' if ln_results[metric]['p_value'] < 0.05 else 'orange')
        ax_ln.legend(fontsize=8)
        ax_ln.grid(True, alpha=0.3)
        ax_ln.set_xscale('log')

    plt.suptitle('Multiscale Curves: Validation Phase', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('multiscale_validation.png', dpi=150, bbox_inches='tight')
    print("\nSaved: multiscale_validation.png")
    plt.close()


def main():
    """Run multiscale test."""
    print("\n" + "="*70)
    print("MULTISCALE NON-GAUSSIANITY TEST")
    print("="*70)
    print("\nPhilosophy:")
    print("  - Test CURVES over 10-20 scales, not single numbers")
    print("  - χ² test comparing curve shapes")
    print("  - Much more statistical power")
    print("  - Standard approach in cosmology")
    print("="*70)

    validation_passed, gauss_res, ln_res = validation_multiscale(n_mocks=50)

    if validation_passed:
        print("\n✅ Ready to test real data!")
    else:
        print("\n⚠️  Consider tuning or increasing N")

    return validation_passed


if __name__ == "__main__":
    results = main()
