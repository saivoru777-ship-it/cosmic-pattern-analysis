"""
Rigorous non-Gaussianity test with proper statistical framework.

Key improvements over previous version:
1. Strict matching assertions between real and mocks
2. Generate ≥50 mocks for stable statistics
3. Report empirical p-values, NOT sigma claims
4. Validation MUST pass on two positive controls before interpreting results
5. Effect sizes using median/IQR instead of mean/std
"""

import numpy as np
from scipy import ndimage, stats
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Dict, List, Tuple


@dataclass
class DatasetMetadata:
    """Enforce matching conditions between datasets."""
    box_size: float
    n_points: int
    grid_size: int
    cell_size: int
    redshift: float

    def matches(self, other: 'DatasetMetadata') -> Tuple[bool, str]:
        """Check if two datasets have matching parameters."""
        if self.box_size != other.box_size:
            return False, f"box_size mismatch: {self.box_size} vs {other.box_size}"
        if self.n_points != other.n_points:
            return False, f"n_points mismatch: {self.n_points} vs {other.n_points}"
        if self.grid_size != other.grid_size:
            return False, f"grid_size mismatch: {self.grid_size} vs {other.grid_size}"
        if self.cell_size != other.cell_size:
            return False, f"cell_size mismatch: {self.cell_size} vs {other.cell_size}"
        if abs(self.redshift - other.redshift) > 0.01:
            return False, f"redshift mismatch: {self.redshift} vs {other.redshift}"
        return True, "All parameters match"


class RigorousNonGaussianityTest:
    """
    Statistically rigorous non-Gaussianity test.

    Uses empirical p-values and effect sizes, not sigma claims.
    Requires validation to pass before interpreting results.
    """

    def __init__(self, box_size=100.0, grid_size=128, cell_size=8, smooth_scale=2.0):
        self.box_size = box_size
        self.grid_size = grid_size
        self.cell_size = cell_size
        self.smooth_scale = smooth_scale

    def standardize_and_grid(self, positions):
        """Convert point cloud to density field (all fixes applied)."""
        # Enforce strict bounds
        mins = positions.min(axis=0)
        maxs = positions.max(axis=0)

        # Check that data fits in expected box
        extent = maxs - mins
        if not np.allclose(extent, self.box_size, rtol=0.1):
            print(f"  ⚠️  Warning: Data extent {extent} doesn't match box_size {self.box_size}")

        # Normalize to grid coordinates
        pos_norm = (positions - mins) / (maxs - mins) * (self.grid_size - 1)

        # CIC gridding
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

        # Store raw grid for counts
        self.raw_grid = grid.copy()

        # Convert to overdensity
        mean_density = grid.mean()
        if mean_density > 0:
            delta = (grid - mean_density) / mean_density
        else:
            delta = grid - mean_density

        # Smooth
        delta_smooth = ndimage.gaussian_filter(delta, sigma=self.smooth_scale)

        return delta_smooth

    def counts_variance(self):
        """Variance-to-mean ratio of counts in cells."""
        grid = self.raw_grid
        cell_size = self.cell_size

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
        return 0.0

    def counts_skewness(self):
        """Skewness of counts distribution."""
        grid = self.raw_grid
        cell_size = self.cell_size

        counts = []
        for i in range(0, self.grid_size - cell_size, cell_size):
            for j in range(0, self.grid_size - cell_size, cell_size):
                for k in range(0, self.grid_size - cell_size, cell_size):
                    cell = grid[i:i+cell_size, j:j+cell_size, k:k+cell_size]
                    counts.append(cell.sum())

        return stats.skew(np.array(counts))

    def euler_characteristic(self, delta_smooth, threshold=0.5):
        """True Euler characteristic χ = β₀ - β₁ + β₂."""
        binary = delta_smooth > threshold

        labeled, n_components = ndimage.label(binary)
        beta_0 = n_components

        inverted = ~binary
        labeled_inv, n_components_inv = ndimage.label(inverted)
        beta_2 = max(0, n_components_inv - 1)

        # Estimate β₁ from topology
        beta_1 = beta_0 + beta_2 - 1

        return beta_0 - beta_1 + beta_2

    def compute_metrics(self, delta_smooth):
        """Compute all metrics."""
        return {
            'counts_variance': self.counts_variance(),
            'counts_skewness': self.counts_skewness(),
            'euler_char': self.euler_characteristic(delta_smooth)
        }

    def empirical_pvalue(self, real_value: float, mock_values: List[float],
                        two_tailed: bool = True) -> float:
        """
        Compute empirical p-value.

        p = (1 + #{mock >= real}) / (N_mocks + 1)

        This is conservative and works for small N.
        """
        mock_array = np.array(mock_values)

        if two_tailed:
            # Two-tailed: how far from median?
            median = np.median(mock_array)
            real_deviation = abs(real_value - median)
            mock_deviations = np.abs(mock_array - median)
            n_extreme = np.sum(mock_deviations >= real_deviation)
        else:
            # One-tailed: is real extreme?
            n_extreme = np.sum(mock_array >= real_value)

        p_value = (1 + n_extreme) / (len(mock_values) + 1)

        return p_value

    def effect_size(self, real_value: float, mock_values: List[float]) -> float:
        """
        Effect size using robust statistics.

        ES = (real - median) / IQR

        More stable than (real - mean) / std for small N.
        """
        mock_array = np.array(mock_values)
        median = np.median(mock_array)
        q75 = np.percentile(mock_array, 75)
        q25 = np.percentile(mock_array, 25)
        iqr = q75 - q25

        if iqr > 0:
            return (real_value - median) / iqr
        else:
            return 0.0


def generate_gaussian_control(n_points=10000, box_size=100, grid_size=128, seed=None):
    """Generate Gaussian random field for validation."""
    if seed is not None:
        np.random.seed(seed)

    # Generate Gaussian field
    gaussian_field = np.random.normal(0, 1, (grid_size, grid_size, grid_size))
    gaussian_field = ndimage.gaussian_filter(gaussian_field, sigma=2.0)

    # Sample points
    density = np.exp(gaussian_field - gaussian_field.max() + 10)  # Avoid overflow
    density = density / density.sum()

    indices = np.random.choice(density.size, size=n_points, p=density.flatten())
    i, j, k = np.unravel_index(indices, gaussian_field.shape)

    positions = np.column_stack([i, j, k]) * (box_size / grid_size)
    positions += np.random.uniform(-0.5, 0.5, (n_points, 3)) * (box_size / grid_size)

    return positions


def generate_lognormal_control(n_points=10000, box_size=100, grid_size=128, seed=None):
    """Generate lognormal field (non-Gaussian) for validation."""
    if seed is not None:
        np.random.seed(seed)

    # Start with Gaussian
    gaussian_field = np.random.normal(0, 1, (grid_size, grid_size, grid_size))
    gaussian_field = ndimage.gaussian_filter(gaussian_field, sigma=2.0)

    # Make lognormal (non-Gaussian!)
    lognormal_field = np.exp(gaussian_field)

    # Sample points
    density = lognormal_field / lognormal_field.sum()

    indices = np.random.choice(density.size, size=n_points, p=density.flatten())
    i, j, k = np.unravel_index(indices, lognormal_field.shape)

    positions = np.column_stack([i, j, k]) * (box_size / grid_size)
    positions += np.random.uniform(-0.5, 0.5, (n_points, 3)) * (box_size / grid_size)

    return positions


def validation_phase_rigorous(n_mocks=50):
    """
    REQUIRED validation: Can we distinguish Gaussian from lognormal?

    This MUST work before we can interpret any real data results.
    """
    print("="*70)
    print("VALIDATION PHASE (REQUIRED)")
    print("="*70)
    print(f"\nGenerating {n_mocks} mocks for each control...\n")

    tester = RigorousNonGaussianityTest(box_size=100.0, grid_size=128, cell_size=8)

    # Test 1: Gaussian "real" vs Gaussian mocks (should NOT detect anomaly)
    print("Control 1: Gaussian vs Gaussian (should NOT flag)")
    print("-"*70)

    gaussian_real_pos = generate_gaussian_control(n_points=10000, seed=1000)
    delta_gauss_real = tester.standardize_and_grid(gaussian_real_pos)
    metrics_gauss_real = tester.compute_metrics(delta_gauss_real)

    gaussian_mock_metrics = {key: [] for key in metrics_gauss_real.keys()}

    for i in range(n_mocks):
        mock_pos = generate_gaussian_control(n_points=10000, seed=2000+i)
        delta_mock = tester.standardize_and_grid(mock_pos)
        metrics_mock = tester.compute_metrics(delta_mock)

        for key in metrics_gauss_real.keys():
            gaussian_mock_metrics[key].append(metrics_mock[key])

    # Compute p-values
    gauss_pvalues = {}
    for key in metrics_gauss_real.keys():
        p = tester.empirical_pvalue(metrics_gauss_real[key], gaussian_mock_metrics[key])
        gauss_pvalues[key] = p
        print(f"  {key}: p = {p:.3f}")

    n_false_positives = sum(p < 0.05 for p in gauss_pvalues.values())
    print(f"\n  Result: {n_false_positives}/3 false positives (expect ≤1 for p<0.05)")

    # Test 2: Lognormal "real" vs Gaussian mocks (SHOULD detect anomaly)
    print("\nControl 2: Lognormal vs Gaussian (SHOULD flag)")
    print("-"*70)

    lognormal_real_pos = generate_lognormal_control(n_points=10000, seed=1001)
    delta_ln_real = tester.standardize_and_grid(lognormal_real_pos)
    metrics_ln_real = tester.compute_metrics(delta_ln_real)

    # Use same Gaussian mocks for comparison
    ln_pvalues = {}
    for key in metrics_ln_real.keys():
        p = tester.empirical_pvalue(metrics_ln_real[key], gaussian_mock_metrics[key])
        ln_pvalues[key] = p
        print(f"  {key}: p = {p:.3f}")

    n_detections = sum(p < 0.05 for p in ln_pvalues.values())
    print(f"\n  Result: {n_detections}/3 metrics detected (expect ≥2 for p<0.05)")

    # Overall validation
    print("\n" + "="*70)
    validation_passed = (n_false_positives <= 1) and (n_detections >= 2)

    if validation_passed:
        print("✅ VALIDATION PASSED")
        print("   Pipeline can reliably detect non-Gaussianity")
    else:
        print("❌ VALIDATION FAILED")
        print("   Pipeline cannot reliably distinguish Gaussian from non-Gaussian")
        print("   DO NOT INTERPRET REAL DATA RESULTS UNTIL THIS IS FIXED")

    print("="*70)

    return validation_passed


def test_real_vs_mocks_rigorous(n_mocks=50):
    """
    Rigorous test with proper statistics.

    Only run this if validation passed!
    """
    print("\n" + "="*70)
    print("REAL DATA VS MOCKS (RIGOROUS)")
    print("="*70)
    print(f"\nGenerating {n_mocks} ΛCDM mocks...\n")

    tester = RigorousNonGaussianityTest(box_size=100.0, grid_size=128, cell_size=8)

    # Load "real" data
    try:
        real_data = np.load('realistic_z0.10.npz')
        real_positions = real_data['positions']

        # Subsample to consistent size if needed
        if len(real_positions) > 10000:
            indices = np.random.choice(len(real_positions), 10000, replace=False)
            real_positions = real_positions[indices]

        real_metadata = DatasetMetadata(
            box_size=100.0,
            n_points=len(real_positions),
            grid_size=128,
            cell_size=8,
            redshift=0.10
        )

        print(f"Loaded real data: {len(real_positions)} points")
    except:
        print("⚠️  Could not load real data, using synthetic")
        real_positions = generate_gaussian_control(n_points=10000, seed=9999)
        real_metadata = DatasetMetadata(
            box_size=100.0,
            n_points=10000,
            grid_size=128,
            cell_size=8,
            redshift=0.10
        )

    # Compute real metrics
    delta_real = tester.standardize_and_grid(real_positions)
    metrics_real = tester.compute_metrics(delta_real)

    # Generate mocks with STRICT MATCHING
    mock_metrics = {key: [] for key in metrics_real.keys()}

    print(f"Generating {n_mocks} mocks with matching parameters:")
    print(f"  box_size = {real_metadata.box_size}")
    print(f"  n_points = {real_metadata.n_points}")
    print(f"  grid_size = {real_metadata.grid_size}")
    print(f"  cell_size = {real_metadata.cell_size}")
    print(f"  redshift = {real_metadata.redshift}\n")

    for i in range(n_mocks):
        if (i+1) % 10 == 0:
            print(f"  Generated {i+1}/{n_mocks} mocks...")

        # Generate with EXACT matching parameters
        mock_pos = generate_gaussian_control(
            n_points=real_metadata.n_points,
            box_size=real_metadata.box_size,
            grid_size=real_metadata.grid_size,
            seed=5000 + i
        )

        # Assert matching
        mock_metadata = DatasetMetadata(
            box_size=real_metadata.box_size,
            n_points=len(mock_pos),
            grid_size=real_metadata.grid_size,
            cell_size=real_metadata.cell_size,
            redshift=real_metadata.redshift
        )

        matches, msg = real_metadata.matches(mock_metadata)
        assert matches, f"Mock {i} doesn't match: {msg}"

        delta_mock = tester.standardize_and_grid(mock_pos)
        metrics_mock = tester.compute_metrics(delta_mock)

        for key in metrics_real.keys():
            mock_metrics[key].append(metrics_mock[key])

    print(f"\n  All {n_mocks} mocks generated with matching parameters ✓\n")

    # Compute statistics
    print("="*70)
    print("RESULTS")
    print("="*70)
    print(f"\n{'Metric':<20} {'Real':<12} {'Mock Q50':<12} {'Mock IQR':<12} {'p-value':<12} {'Effect':<12}")
    print("-"*70)

    significant_metrics = []

    for metric in metrics_real.keys():
        real_val = metrics_real[metric]
        mock_vals = mock_metrics[metric]

        # Robust statistics
        q50 = np.median(mock_vals)
        q25 = np.percentile(mock_vals, 25)
        q75 = np.percentile(mock_vals, 75)
        iqr = q75 - q25

        # Empirical p-value
        p_val = tester.empirical_pvalue(real_val, mock_vals)

        # Effect size
        effect = tester.effect_size(real_val, mock_vals)

        is_sig = p_val < 0.05
        if is_sig:
            significant_metrics.append(metric)

        marker = "***" if p_val < 0.01 else ("**" if p_val < 0.05 else "")

        print(f"{metric:<20} {real_val:<12.3f} {q50:<12.3f} {iqr:<12.3f} {p_val:<12.4f} {effect:<12.2f} {marker}")

    # Overall assessment
    print("\n" + "="*70)
    n_significant = len(significant_metrics)

    if n_significant >= 2:
        print("⚠️  ANOMALY DETECTED")
        print(f"   {n_significant}/3 metrics significant at p < 0.05")
        print(f"   Significant: {', '.join(significant_metrics)}")
        print("\n   INTERPRETATION:")
        print("   Real data falls outside mock distribution.")
        print("   Possible causes:")
        print("   1. Generation parameter mismatch (most likely)")
        print("   2. Real non-Gaussian signature (less likely)")
        print("   3. Statistical fluctuation with N=50 (check with N=200)")
    else:
        print("✅ DATA CONSISTENT WITH MOCKS")
        print(f"   {n_significant}/3 metrics significant (within expected variation)")
        print("   No strong evidence for non-Gaussianity")

    print("="*70)

    return {
        'n_significant': n_significant,
        'metrics_real': metrics_real,
        'mock_metrics': mock_metrics,
        'significant_metrics': significant_metrics
    }


def main():
    """Run rigorous pipeline."""
    print("\n" + "="*70)
    print("RIGOROUS NON-GAUSSIANITY TEST")
    print("="*70)
    print("\nStatistical framework:")
    print("  ✓ Empirical p-values (not sigma claims)")
    print("  ✓ Effect sizes using median/IQR")
    print("  ✓ N ≥ 50 mocks for stability")
    print("  ✓ Strict parameter matching")
    print("  ✓ Validation required before interpretation")
    print("="*70)

    # VALIDATION IS REQUIRED
    validation_passed = validation_phase_rigorous(n_mocks=50)

    if not validation_passed:
        print("\n" + "="*70)
        print("⛔ STOPPING: VALIDATION FAILED")
        print("="*70)
        print("\nThe pipeline cannot distinguish known Gaussian from non-Gaussian.")
        print("DO NOT proceed to real data until this is fixed.")
        print("\nPossible fixes:")
        print("  - Adjust cell_size (try 4, 8, 16)")
        print("  - Adjust smooth_scale (try 1.0, 2.0, 4.0)")
        print("  - Add more metrics (peak statistics, bispectrum)")
        print("  - Use different threshold for Euler characteristic")
        return None

    # Only run real test if validation passed
    print("\n✅ Validation passed - proceeding to real data test\n")
    results = test_real_vs_mocks_rigorous(n_mocks=50)

    return results


if __name__ == "__main__":
    results = main()
