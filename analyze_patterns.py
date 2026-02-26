"""
Analyze synthetic universe patterns using topological data analysis.

Tests:
1. Persistent homology (detect hierarchical structure)
2. Scale hierarchy detection (look for golden ratio or other preferred ratios)
3. Power spectrum analysis
4. Motif counting
"""

import numpy as np
from ripser import ripser
from scipy.ndimage import gaussian_filter
from scipy.stats import ks_2samp, pearsonr
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN


class PatternAnalyzer:
    """Analyze patterns for tiling-like signatures."""

    def __init__(self, field):
        """
        Parameters:
        - field: 2D numpy array representing density field
        """
        self.field = field
        self.grid_size = field.shape[0]

    def extract_high_density_points(self, threshold_percentile=80, max_points=2000):
        """
        Extract high-density points for topological analysis.

        Parameters:
        - threshold_percentile: keep points above this density percentile
        - max_points: maximum number of points to keep (for computational efficiency)
        """
        threshold = np.percentile(self.field, threshold_percentile)
        y, x = np.where(self.field > threshold)

        # Add density as third dimension for weighted analysis
        densities = self.field[y, x]

        points = np.column_stack([x, y])

        # Subsample if too many points
        if len(points) > max_points:
            indices = np.random.choice(len(points), max_points, replace=False)
            points = points[indices]

        return points

    def compute_persistence(self, points, max_dim=1):
        """
        Compute persistent homology.

        Returns:
        - diagrams: persistence diagrams for H0 (components), H1 (loops)
        """
        print(f"  Computing persistence for {len(points)} points...")

        # Normalize points to unit square for numerical stability
        points_normalized = points / self.grid_size

        # Compute persistence
        result = ripser(points_normalized, maxdim=max_dim)

        return result['dgms']

    def analyze_persistence_statistics(self, diagrams):
        """
        Extract statistical features from persistence diagrams.

        Key metrics:
        - Persistence distribution (birth-death)
        - Number of significant features
        - Scale ratios between features
        """
        stats = {}

        for dim, dgm in enumerate(diagrams):
            # Remove infinite persistence (connected component at H0)
            dgm_finite = dgm[dgm[:, 1] < np.inf]

            if len(dgm_finite) == 0:
                continue

            # Persistence = death - birth
            persistence = dgm_finite[:, 1] - dgm_finite[:, 0]

            # Sort by persistence
            sorted_persistence = np.sort(persistence)[::-1]

            # Store statistics
            stats[f'H{dim}_num_features'] = len(dgm_finite)
            stats[f'H{dim}_total_persistence'] = np.sum(persistence)
            stats[f'H{dim}_max_persistence'] = np.max(persistence) if len(persistence) > 0 else 0
            stats[f'H{dim}_mean_persistence'] = np.mean(persistence) if len(persistence) > 0 else 0
            stats[f'H{dim}_std_persistence'] = np.std(persistence) if len(persistence) > 0 else 0

            # Check for scale hierarchy (ratios between consecutive features)
            if len(sorted_persistence) >= 3:
                ratios = sorted_persistence[:-1] / sorted_persistence[1:]
                stats[f'H{dim}_scale_ratios_mean'] = np.mean(ratios)
                stats[f'H{dim}_scale_ratios_std'] = np.std(ratios)

                # Check if ratios cluster around golden ratio (phi ≈ 1.618)
                phi = (1 + np.sqrt(5)) / 2
                ratios_near_phi = np.sum(np.abs(ratios - phi) < 0.2)
                stats[f'H{dim}_ratios_near_phi'] = ratios_near_phi / len(ratios) if len(ratios) > 0 else 0

        return stats

    def compute_power_spectrum(self):
        """
        Compute 2D power spectrum P(k).
        Look for extra structure beyond smooth power law.
        """
        # FFT
        field_fft = np.fft.fft2(self.field)
        power = np.abs(field_fft)**2

        # Radial averaging
        k = np.fft.fftfreq(self.grid_size)
        kx, ky = np.meshgrid(k, k)
        k_mag = np.sqrt(kx**2 + ky**2)

        # Bin by k
        k_bins = np.linspace(0, 0.5, 50)
        k_centers = (k_bins[:-1] + k_bins[1:]) / 2
        power_binned = []

        for i in range(len(k_bins) - 1):
            mask = (k_mag >= k_bins[i]) & (k_mag < k_bins[i+1])
            if np.any(mask):
                power_binned.append(np.mean(power[mask]))
            else:
                power_binned.append(0)

        power_binned = np.array(power_binned)

        return k_centers, power_binned

    def detect_scale_hierarchy(self, smoothing_scales=[2, 4, 8, 16, 32]):
        """
        Detect if structure exists at discrete hierarchical scales.

        Method: Smooth field at multiple scales, count peaks, check for ratio patterns.
        """
        peak_counts = []
        phi = (1 + np.sqrt(5)) / 2

        for sigma in smoothing_scales:
            smoothed = gaussian_filter(self.field, sigma=sigma)

            # Count peaks (local maxima)
            from scipy.ndimage import maximum_filter
            local_max = (smoothed == maximum_filter(smoothed, size=5))
            n_peaks = np.sum(local_max)

            peak_counts.append(n_peaks)

        peak_counts = np.array(peak_counts)

        # Check for exponential relationship (power-law scaling)
        log_scales = np.log(smoothing_scales)
        log_counts = np.log(peak_counts + 1)  # +1 to avoid log(0)

        # Fit power law
        slope, intercept = np.polyfit(log_scales, log_counts, 1)

        # Check scale ratios
        scale_ratios = np.array(smoothing_scales[1:]) / np.array(smoothing_scales[:-1])

        stats = {
            'peak_counts': peak_counts,
            'smoothing_scales': smoothing_scales,
            'power_law_slope': slope,
            'scale_ratios': scale_ratios,
            'scale_ratios_mean': np.mean(scale_ratios),
            'scale_ratios_std': np.std(scale_ratios),
        }

        return stats

    def morphology_classification(self, smoothing_scale=5):
        """
        Classify each point as node, filament, sheet, or void.
        Based on Hessian eigenvalue method.
        """
        # Smooth field
        smoothed = gaussian_filter(self.field, sigma=smoothing_scale)

        # Compute gradients
        grad_y, grad_x = np.gradient(smoothed)

        # Compute Hessian
        grad_yy, grad_yx = np.gradient(grad_y)
        grad_xy, grad_xx = np.gradient(grad_x)

        # Eigenvalues at each point
        morphology = np.zeros_like(smoothed, dtype=int)

        for i in range(self.grid_size):
            for j in range(self.grid_size):
                hessian = np.array([
                    [grad_xx[i, j], grad_xy[i, j]],
                    [grad_yx[i, j], grad_yy[i, j]]
                ])

                eigenvalues = np.linalg.eigvalsh(hessian)
                eigenvalues = np.sort(eigenvalues)[::-1]  # Sort descending

                # Classification
                if eigenvalues[0] < 0 and eigenvalues[1] < 0:
                    morphology[i, j] = 0  # Peak/Node
                elif eigenvalues[0] < 0 and eigenvalues[1] >= 0:
                    morphology[i, j] = 1  # Filament
                elif eigenvalues[0] >= 0 and eigenvalues[1] >= 0:
                    morphology[i, j] = 2  # Void

        # Count each type
        n_nodes = np.sum(morphology == 0)
        n_filaments = np.sum(morphology == 1)
        n_voids = np.sum(morphology == 2)

        total = n_nodes + n_filaments + n_voids
        if total == 0:
            total = 1  # Avoid division by zero

        stats = {
            'n_nodes': n_nodes,
            'n_filaments': n_filaments,
            'n_voids': n_voids,
            'ratio_nodes': n_nodes / total,
            'ratio_filaments': n_filaments / total,
            'ratio_voids': n_voids / total,
        }

        return stats, morphology


def compare_patterns(fields_dict):
    """
    Run full comparison across all patterns.

    Parameters:
    - fields_dict: dictionary of {name: field} pairs
    """
    results = {}

    print("\n" + "="*60)
    print("ANALYZING PATTERNS FOR TILING SIGNATURES")
    print("="*60)

    for name, field in fields_dict.items():
        print(f"\n\nAnalyzing: {name}")
        print("-" * 50)

        analyzer = PatternAnalyzer(field)

        # 1. Persistent homology
        print("1. Extracting high-density points...")
        points = analyzer.extract_high_density_points()

        print("2. Computing persistent homology...")
        diagrams = analyzer.compute_persistence(points)

        print("3. Analyzing persistence statistics...")
        persistence_stats = analyzer.analyze_persistence_statistics(diagrams)

        # 2. Power spectrum
        print("4. Computing power spectrum...")
        k_centers, power = analyzer.compute_power_spectrum()

        # 3. Scale hierarchy
        print("5. Detecting scale hierarchy...")
        hierarchy_stats = analyzer.detect_scale_hierarchy()

        # 4. Morphology
        print("6. Classifying morphology...")
        morphology_stats, morphology_map = analyzer.morphology_classification()

        # Store results
        results[name] = {
            'persistence_stats': persistence_stats,
            'power_spectrum': (k_centers, power),
            'hierarchy_stats': hierarchy_stats,
            'morphology_stats': morphology_stats,
            'diagrams': diagrams,
            'points': points,
        }

    return results


def main():
    """Run full analysis pipeline."""
    print("Loading synthetic universes...")

    fields = {
        'Random (ΛCDM-like)': np.load('random_field.npy'),
        'Evolved (Gravity)': np.load('evolved_field.npy'),
        'Substitution Tiling': np.load('tiling_field.npy'),
        'Penrose-Inspired': np.load('penrose_field.npy'),
    }

    # Run analysis
    results = compare_patterns(fields)

    # Save results
    print("\n\nSaving results...")
    np.save('analysis_results.npy', results, allow_pickle=True)
    print("Results saved to: analysis_results.npy")

    return results


if __name__ == "__main__":
    results = main()
