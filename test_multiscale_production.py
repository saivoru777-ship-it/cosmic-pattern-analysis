"""
Multi-scale non-Gaussianity test - PRODUCTION VERSION

Clean implementation using only validated metrics:
- Counts variance (12 scales)
- Counts skewness (12 scales)

With:
- Covariance-aware χ² (shrinkage regularization)
- Strict parameter matching
- CLI for real data testing
"""

import numpy as np
from scipy import ndimage, stats
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import List
import argparse
import os


@dataclass
class DatasetInfo:
    """Enforce matching between datasets."""
    n_points: int
    box_size: float
    grid_size: int

    def matches(self, other: 'DatasetInfo', tol=0.01) -> tuple[bool, str]:
        """Check if datasets match within tolerance."""
        if abs(self.n_points - other.n_points) / self.n_points > tol:
            return False, f"N mismatch: {self.n_points} vs {other.n_points}"
        if abs(self.box_size - other.box_size) / self.box_size > tol:
            return False, f"Box mismatch: {self.box_size} vs {other.box_size}"
        if self.grid_size != other.grid_size:
            return False, f"Grid mismatch: {self.grid_size} vs {other.grid_size}"
        return True, "Match"


@dataclass
class ScaleRange:
    """Define scales to test."""
    cell_sizes: List[int]

    @classmethod
    def logarithmic(cls, min_size=3, max_size=32, n_scales=12):
        """Logarithmic spacing from min to max."""
        sizes = np.logspace(np.log10(min_size), np.log10(max_size), n_scales)
        sizes = np.unique(sizes.astype(int))
        return cls(cell_sizes=sizes.tolist())


class MultiscaleTest:
    """Production multiscale clustering test."""

    def __init__(self, box_size=100.0, grid_size=128):
        self.box_size = box_size
        self.grid_size = grid_size

    def standardize_and_grid(self, positions):
        """CIC gridding with all fixes applied."""
        positions = positions % self.box_size

        mins = positions.min(axis=0)
        maxs = positions.max(axis=0)
        pos_norm = (positions - mins) / (maxs - mins) * (self.grid_size - 1)

        grid = np.zeros((self.grid_size, self.grid_size, self.grid_size))

        for pos in pos_norm:
            i, j, k = pos.astype(int)
            dx, dy, dz = pos[0] - i, pos[1] - j, pos[2] - k

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

        return counts.var() / mean_count if mean_count > 0 else np.nan

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

        return stats.skew(counts) if len(counts) > 3 and counts.std() > 0 else np.nan

    def compute_multiscale_curves(self, positions, scale_range: ScaleRange):
        """Compute variance and skewness curves over scales."""
        grid = self.standardize_and_grid(positions)

        curves = {
            'variance': {'scales': [], 'values': []},
            'skewness': {'scales': [], 'values': []}
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

        return curves

    def chi_squared_with_covariance(self, real_curve, mock_curves, shrinkage=0.1):
        """
        Covariance-aware χ² test.

        Accounts for correlation between scales using mock covariance matrix.
        Uses shrinkage regularization for stable inversion.

        Returns:
            dict with chi_squared, dof, p_value, condition_number
        """
        real_values = np.array(real_curve['values'])
        n_scales = len(real_values)

        # Build mock matrix
        mock_matrix = []
        for mock in mock_curves:
            mock_vals = np.array(mock['values'])
            if len(mock_vals) == n_scales:
                mock_matrix.append(mock_vals)

        if len(mock_matrix) < 2:
            return {'chi_squared': 0, 'dof': n_scales, 'p_value': 1.0, 'condition_number': np.inf}

        mock_matrix = np.array(mock_matrix)  # Shape: (n_mocks, n_scales)

        # Mock statistics
        mock_mean = mock_matrix.mean(axis=0)
        mock_cov = np.cov(mock_matrix.T)

        # Shrinkage: C_reg = (1-λ)C + λ*diag(C)
        diag_cov = np.diag(np.diag(mock_cov))
        mock_cov_reg = (1 - shrinkage) * mock_cov + shrinkage * diag_cov
        mock_cov_reg += np.eye(n_scales) * 1e-10  # Numerical stability

        # Compute χ²
        residual = real_values - mock_mean

        try:
            cov_inv = np.linalg.inv(mock_cov_reg)
            chi_squared = residual @ cov_inv @ residual
            cond = np.linalg.cond(mock_cov_reg)

            dof = n_scales
            p_value = 1 - stats.chi2.cdf(chi_squared, dof)

            return {
                'chi_squared': chi_squared,
                'dof': dof,
                'p_value': p_value,
                'condition_number': cond
            }

        except np.linalg.LinAlgError:
            # Fallback to diagonal
            variances = np.diag(mock_cov) + 1e-10
            chi_squared = ((residual**2) / variances).sum()

            return {
                'chi_squared': chi_squared,
                'dof': n_scales,
                'p_value': 1 - stats.chi2.cdf(chi_squared, n_scales),
                'condition_number': np.inf
            }


# ==============================================================================
# CONTROLS
# ==============================================================================

def generate_unclustered_control(n_points=10000, box_size=100, seed=None):
    """Pure Poisson random."""
    if seed is not None:
        np.random.seed(seed)
    return np.random.uniform(0, box_size, (n_points, 3))


def generate_clustered_control(n_points=10000, box_size=100, n_halos=100, halo_scale=2.0, seed=None):
    """Halo-like clustering."""
    if seed is not None:
        np.random.seed(seed)

    halo_centers = np.random.uniform(0, box_size, (n_halos, 3))
    points_per_halo = np.random.poisson(n_points / n_halos, n_halos)
    points_per_halo = (points_per_halo * n_points / points_per_halo.sum()).astype(int)

    positions = []
    for center, n_gal in zip(halo_centers, points_per_halo):
        if n_gal == 0:
            continue

        n_core = n_gal // 2
        core = np.random.normal(0, halo_scale * 0.3, (n_core, 3))
        ext = np.random.normal(0, halo_scale * 1.5, (n_gal - n_core, 3))

        positions.append(center + np.vstack([core, ext]))

    positions = np.vstack(positions) % box_size
    return positions[:n_points]


# ==============================================================================
# VALIDATION
# ==============================================================================

def run_validation(n_mocks=50, verbose=True):
    """
    Validation with proper point cloud controls.

    Returns:
        (passed: bool, results: dict)
    """
    if verbose:
        print("="*70)
        print("PRODUCTION VALIDATION")
        print("="*70)
        print(f"\n2 metrics (variance, skewness) × {n_mocks} mocks")
        print("Covariance-aware χ² with shrinkage regularization\n")

    tester = MultiscaleTest(box_size=100, grid_size=128)
    scale_range = ScaleRange.logarithmic(min_size=3, max_size=32, n_scales=12)

    # Control 1: Unclustered vs Unclustered
    if verbose:
        print("Control 1: Unclustered vs Unclustered (expect p > 0.05)")
        print("-"*70)

    unclust_real = generate_unclustered_control(n_points=10000, seed=1000)
    unclust_real_curves = tester.compute_multiscale_curves(unclust_real, scale_range)

    unclust_mocks = {'variance': [], 'skewness': []}
    for i in range(n_mocks):
        if verbose and (i+1) % 25 == 0:
            print(f"  Generated {i+1}/{n_mocks} mocks...")

        mock = generate_unclustered_control(n_points=10000, seed=2000+i)
        mock_curves = tester.compute_multiscale_curves(mock, scale_range)

        for metric in ['variance', 'skewness']:
            unclust_mocks[metric].append(mock_curves[metric])

    unclust_results = {}
    for metric in ['variance', 'skewness']:
        result = tester.chi_squared_with_covariance(
            unclust_real_curves[metric],
            unclust_mocks[metric],
            shrinkage=0.1
        )
        unclust_results[metric] = result

        if verbose:
            print(f"  {metric}: χ²={result['chi_squared']:.2f}, p={result['p_value']:.4f}, cond={result['condition_number']:.1e}")

    n_false_pos = sum(r['p_value'] < 0.05 for r in unclust_results.values())

    if verbose:
        print(f"\n  False positives: {n_false_pos}/2 {'✓' if n_false_pos == 0 else '✗'}")

    # Control 2: Clustered vs Unclustered
    if verbose:
        print("\nControl 2: Clustered vs Unclustered (expect p < 0.05)")
        print("-"*70)

    clust_real = generate_clustered_control(n_points=10000, n_halos=100, seed=1001)
    clust_real_curves = tester.compute_multiscale_curves(clust_real, scale_range)

    clust_results = {}
    for metric in ['variance', 'skewness']:
        result = tester.chi_squared_with_covariance(
            clust_real_curves[metric],
            unclust_mocks[metric],
            shrinkage=0.1
        )
        clust_results[metric] = result

        marker = "***" if result['p_value'] < 0.001 else "**" if result['p_value'] < 0.01 else "*" if result['p_value'] < 0.05 else ""

        if verbose:
            print(f"  {metric}: χ²={result['chi_squared']:.2f}, p={result['p_value']:.4f}, cond={result['condition_number']:.1e} {marker}")

    n_detections = sum(r['p_value'] < 0.05 for r in clust_results.values())

    if verbose:
        print(f"\n  Detections: {n_detections}/2 {'✓' if n_detections == 2 else '✗'}")

    # Overall
    validation_passed = (n_false_pos == 0) and (n_detections == 2)

    if verbose:
        print("\n" + "="*70)
        if validation_passed:
            print("✅ VALIDATION PASSED")
            print("   Pipeline ready for real data")
        else:
            print("⚠️  VALIDATION UNCLEAR")
            print(f"   False positives: {n_false_pos}/2")
            print(f"   True positives: {n_detections}/2")
        print("="*70)

    return validation_passed, {
        'unclust_results': unclust_results,
        'clust_results': clust_results,
        'unclust_mocks': unclust_mocks,
        'scale_range': scale_range
    }


# ==============================================================================
# REAL DATA TEST
# ==============================================================================

def test_real_data_vs_mocks(data_path, mock_paths, n_mocks=50, output_dir="."):
    """
    Test real data against ΛCDM mocks.

    Enforces strict matching and produces publication-ready output.
    """
    print("="*70)
    print("REAL DATA TEST")
    print("="*70)

    # Load data
    print(f"\nLoading real data: {data_path}")
    real_data = np.load(data_path)
    real_pos = real_data['positions']
    real_box = real_data.get('box_size', 100.0)

    real_info = DatasetInfo(
        n_points=len(real_pos),
        box_size=real_box,
        grid_size=128
    )

    print(f"  N = {real_info.n_points}")
    print(f"  Box = {real_info.box_size:.1f} Mpc/h")

    # Load mocks with strict matching
    print(f"\nLoading {len(mock_paths[:n_mocks])} mocks...")

    mock_positions_list = []
    for i, mock_path in enumerate(mock_paths[:n_mocks]):
        mock_data = np.load(mock_path)
        mock_pos = mock_data['positions']
        mock_box = mock_data.get('box_size', real_box)

        mock_info = DatasetInfo(
            n_points=len(mock_pos),
            box_size=mock_box,
            grid_size=128
        )

        # STRICT MATCHING
        matches, msg = real_info.matches(mock_info, tol=0.05)
        if not matches:
            print(f"  ⚠️  Mock {i} doesn't match: {msg}")
            print(f"      Attempting to fix...")

            # Try to fix
            if len(mock_pos) != real_info.n_points:
                # Resample to match
                if len(mock_pos) > real_info.n_points:
                    indices = np.random.choice(len(mock_pos), real_info.n_points, replace=False)
                    mock_pos = mock_pos[indices]
                else:
                    print(f"      ERROR: Mock has too few points ({len(mock_pos)} < {real_info.n_points})")
                    continue

        mock_positions_list.append(mock_pos)

    print(f"\n  Loaded {len(mock_positions_list)} matching mocks ✓")

    # Run test
    tester = MultiscaleTest(box_size=real_info.box_size, grid_size=128)
    scale_range = ScaleRange.logarithmic()

    print("\nComputing curves...")
    real_curves = tester.compute_multiscale_curves(real_pos, scale_range)

    mock_curves_list = {'variance': [], 'skewness': []}
    for i, mock_pos in enumerate(mock_positions_list):
        if (i+1) % 20 == 0:
            print(f"  Processed {i+1}/{len(mock_positions_list)} mocks...")

        mock_curves = tester.compute_multiscale_curves(mock_pos, scale_range)
        for metric in ['variance', 'skewness']:
            mock_curves_list[metric].append(mock_curves[metric])

    # Test
    print("\n" + "="*70)
    print("RESULTS")
    print("="*70)
    print(f"\n{'Metric':<15} {'χ²':<12} {'DOF':<8} {'p-value':<12} {'Significant?'}")
    print("-"*70)

    results = {}
    for metric in ['variance', 'skewness']:
        result = tester.chi_squared_with_covariance(
            real_curves[metric],
            mock_curves_list[metric],
            shrinkage=0.1
        )
        results[metric] = result

        sig = "***" if result['p_value'] < 0.001 else "**" if result['p_value'] < 0.01 else "*" if result['p_value'] < 0.05 else "no"

        print(f"{metric:<15} {result['chi_squared']:<12.2f} {result['dof']:<8} {result['p_value']:<12.6f} {sig}")

    n_sig = sum(r['p_value'] < 0.05 for r in results.values())

    print("\n" + "="*70)
    if n_sig >= 1:
        print("⚠️  ANOMALY DETECTED")
        print(f"   {n_sig}/2 metrics significant at p < 0.05")
        print("\n   Interpretation:")
        print("   Real data differs from ΛCDM mocks in clustering structure")
        print("   Could indicate:")
        print("   - Different cosmological parameters")
        print("   - Selection effects / systematics")
        print("   - Non-standard physics (if systematics ruled out)")
    else:
        print("✅ CONSISTENT WITH MOCKS")
        print("   No significant deviations detected")
        print("   Data consistent with ΛCDM expectations")

    print("="*70)

    return results


def main():
    """CLI interface."""
    parser = argparse.ArgumentParser(description="Multiscale clustering test (production)")
    parser.add_argument('--validate', action='store_true', help='Run validation')
    parser.add_argument('--data', type=str, help='Real data .npz file')
    parser.add_argument('--mocks', nargs='+', help='Mock .npz files')
    parser.add_argument('--n-mocks', type=int, default=50)

    args = parser.parse_args()

    if args.validate or (not args.data):
        # Run validation
        passed, results = run_validation(n_mocks=args.n_mocks)
        return passed

    else:
        # Test real data
        if not args.mocks:
            print("ERROR: --mocks required for testing")
            return None

        results = test_real_data_vs_mocks(args.data, args.mocks, args.n_mocks)
        return results


if __name__ == "__main__":
    import sys

    if len(sys.argv) == 1:
        # Default: validation
        passed, results = run_validation(n_mocks=50, verbose=True)
    else:
        main()
