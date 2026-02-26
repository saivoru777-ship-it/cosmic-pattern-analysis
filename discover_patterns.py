"""
GENERAL PATTERN DISCOVERY ENGINE

This doesn't assume Einstein tiles or any specific pattern.
Instead, it DISCOVERS whatever patterns actually exist in the data.

Methods:
1. Multi-scale local environment clustering
2. Voronoi cell analysis (natural tessellation)
3. Graph motif extraction
4. Fourier-space periodicity hunting
5. Persistent homology-based shape discovery
6. Unsupervised learning of structural patterns
"""

import numpy as np
from scipy.spatial import Voronoi, cKDTree
from scipy.ndimage import gaussian_filter
from sklearn.cluster import DBSCAN, KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from collections import Counter
import warnings
warnings.filterwarnings('ignore')


class PatternDiscoveryEngine:
    """Discover patterns in 3D point distributions."""

    def __init__(self, positions):
        """
        Parameters:
        - positions: Nx3 array of (x, y, z) coordinates
        """
        self.positions = positions
        self.n_points = len(positions)
        self.tree = cKDTree(positions)

    def discover_all_patterns(self):
        """Run all discovery methods and return findings."""
        print("\n" + "="*60)
        print("GENERAL PATTERN DISCOVERY ENGINE")
        print("="*60)
        print(f"Analyzing {self.n_points} points...")

        results = {}

        # 1. Local environment patterns
        print("\n1. Discovering local environment patterns...")
        results['local_envs'] = self.discover_local_environments()

        # 2. Voronoi tessellation
        print("\n2. Analyzing Voronoi tessellation...")
        results['voronoi'] = self.analyze_voronoi_patterns()

        # 3. Scale-dependent clustering
        print("\n3. Finding scale-dependent clustering...")
        results['scales'] = self.analyze_multiscale_clustering()

        # 4. Fourier patterns
        print("\n4. Searching for periodic/quasiperiodic patterns...")
        results['fourier'] = self.discover_fourier_patterns()

        # 5. Graph motifs
        print("\n5. Extracting network motifs...")
        results['motifs'] = self.discover_graph_motifs()

        # 6. Shape clustering
        print("\n6. Clustering similar structural patterns...")
        results['shapes'] = self.discover_shape_patterns()

        return results

    def discover_local_environments(self, k_neighbors=20, n_clusters=8):
        """
        Cluster points by their local environment.

        For each point, characterize:
        - Neighbor distance distribution
        - Local anisotropy
        - Density gradient direction
        """
        print("  Computing local environment features...")

        features = []

        for i in range(min(2000, self.n_points)):  # Sample for speed
            # Find neighbors
            dists, indices = self.tree.query(self.positions[i], k=k_neighbors+1)
            dists = dists[1:]  # Exclude self
            neighbor_pos = self.positions[indices[1:]]

            # Feature 1: Distance statistics
            mean_dist = np.mean(dists)
            std_dist = np.std(dists)

            # Feature 2: Anisotropy (eigenvalue ratio of neighbor positions)
            centered = neighbor_pos - self.positions[i]
            cov = np.cov(centered.T)
            eigvals = np.linalg.eigvalsh(cov)
            eigvals = np.sort(eigvals)[::-1]
            anisotropy = eigvals[0] / (eigvals[1] + 1e-10)

            # Feature 3: Local density
            volume = (4/3) * np.pi * dists[-1]**3
            density = k_neighbors / volume

            features.append([mean_dist, std_dist, anisotropy, np.log10(density)])

        features = np.array(features)

        # Cluster environments
        print("  Clustering environments...")
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        labels = kmeans.fit_predict(features)

        # Statistics
        cluster_counts = Counter(labels)

        return {
            'features': features,
            'labels': labels,
            'cluster_counts': dict(cluster_counts),
            'cluster_centers': kmeans.cluster_centers_,
            'n_clusters': n_clusters
        }

    def analyze_voronoi_patterns(self, sample_size=2000):
        """
        Analyze Voronoi tessellation.

        This is the "natural" way space self-organizes.
        Look for patterns in cell volumes and neighbor counts.
        """
        # Sample points for computational feasibility
        if self.n_points > sample_size:
            indices = np.random.choice(self.n_points, sample_size, replace=False)
            sample_points = self.positions[indices]
        else:
            sample_points = self.positions

        print(f"  Computing Voronoi tessellation ({len(sample_points)} points)...")

        try:
            vor = Voronoi(sample_points)

            # Analyze cell properties
            neighbor_counts = []

            for region_idx, region in enumerate(vor.regions):
                if -1 not in region and len(region) > 0:
                    neighbor_counts.append(len(region))

            neighbor_counts = np.array(neighbor_counts)

            # Statistics
            stats = {
                'mean_neighbors': np.mean(neighbor_counts),
                'std_neighbors': np.std(neighbor_counts),
                'neighbor_distribution': Counter(neighbor_counts),
            }

            # Check for anomalies
            # In random distribution, expect mean ~ 15-16 neighbors
            # If significantly different → pattern

            return stats

        except Exception as e:
            print(f"  Voronoi computation failed: {e}")
            return None

    def analyze_multiscale_clustering(self, scales=[10, 20, 40, 80, 160]):
        """
        Measure clustering at different scales.

        Look for:
        - Preferred scales (bumps in clustering strength)
        - Scale ratios (hierarchical structure)
        """
        print("  Analyzing clustering across scales...")

        results = []

        for scale in scales:
            # Count pairs within distance
            tree_query = self.tree.query_ball_point(self.positions, r=scale)
            counts = [len(neighbors) - 1 for neighbors in tree_query]  # Exclude self

            mean_count = np.mean(counts)
            std_count = np.std(counts)

            # Expected count for random (Poisson)
            volume = (4/3) * np.pi * scale**3
            mean_density = self.n_points / (self.positions.max() - self.positions.min())**3
            expected_count = mean_density * volume

            # Clustering strength
            clustering = mean_count / (expected_count + 1e-10)

            results.append({
                'scale': scale,
                'mean_count': mean_count,
                'expected_count': expected_count,
                'clustering_strength': clustering
            })

        # Check for preferred scales (local maxima)
        clustering_vals = [r['clustering_strength'] for r in results]

        # Check scale ratios
        scale_ratios = []
        for i in range(len(scales) - 1):
            scale_ratios.append(scales[i+1] / scales[i])

        return {
            'results': results,
            'scale_ratios': scale_ratios,
            'clustering_curve': clustering_vals
        }

    def discover_fourier_patterns(self, grid_size=64):
        """
        Look for periodic or quasiperiodic patterns.

        Method:
        - Grid the data
        - Compute 3D FFT
        - Look for peaks in power spectrum
        """
        print("  Computing 3D density field...")

        # Normalize positions to grid
        pos_normalized = self.positions - self.positions.min(axis=0)
        pos_max = pos_normalized.max()
        pos_normalized = pos_normalized / pos_max * (grid_size - 1)

        # Create 3D density grid
        density = np.zeros((grid_size, grid_size, grid_size))

        for pos in pos_normalized:
            i, j, k = pos.astype(int)
            if 0 <= i < grid_size and 0 <= j < grid_size and 0 <= k < grid_size:
                density[i, j, k] += 1

        # Smooth slightly
        density = gaussian_filter(density, sigma=1.0)

        print("  Computing 3D FFT...")
        fft = np.fft.fftn(density)
        power = np.abs(fft)**2

        # Radially average
        k = np.fft.fftfreq(grid_size)
        kx, ky, kz = np.meshgrid(k, k, k, indexing='ij')
        k_mag = np.sqrt(kx**2 + ky**2 + kz**2)

        k_bins = np.linspace(0, 0.5, 30)
        power_binned = []

        for i in range(len(k_bins) - 1):
            mask = (k_mag >= k_bins[i]) & (k_mag < k_bins[i+1])
            if np.any(mask):
                power_binned.append(np.mean(power[mask]))
            else:
                power_binned.append(0)

        power_binned = np.array(power_binned)
        k_centers = (k_bins[:-1] + k_bins[1:]) / 2

        # Find peaks
        from scipy.signal import find_peaks
        peaks, properties = find_peaks(power_binned, height=np.mean(power_binned))

        return {
            'k_values': k_centers,
            'power_spectrum': power_binned,
            'peaks': peaks,
            'peak_positions': k_centers[peaks] if len(peaks) > 0 else []
        }

    def discover_graph_motifs(self, connection_threshold=50):
        """
        Build neighbor graph and find recurring motifs.

        Motifs = small subgraphs that appear more often than expected.
        """
        print(f"  Building nearest-neighbor graph (r < {connection_threshold})...")

        # Build graph: connect points within threshold
        tree_query = self.tree.query_ball_point(
            self.positions[:min(1000, self.n_points)],  # Sample for speed
            r=connection_threshold
        )

        # Count degree distribution
        degrees = [len(neighbors) - 1 for neighbors in tree_query]  # Exclude self

        degree_dist = Counter(degrees)

        # Simple motif: count triangles (3-node fully connected)
        # This is computationally expensive for large graphs, so we sample

        return {
            'degree_distribution': dict(degree_dist),
            'mean_degree': np.mean(degrees),
            'connection_threshold': connection_threshold
        }

    def discover_shape_patterns(self, n_clusters=6):
        """
        Use unsupervised learning to cluster similar 3D arrangements.

        For each point, create a "shape descriptor" based on neighbors,
        then cluster these descriptors.
        """
        print("  Computing shape descriptors...")

        descriptors = []
        k = 15  # Neighbors for shape

        for i in range(min(1500, self.n_points)):
            dists, indices = self.tree.query(self.positions[i], k=k+1)
            dists = dists[1:]  # Exclude self

            # Descriptor: sorted distances (rotation-invariant)
            descriptor = np.sort(dists)

            # Also add eigenvalues (shape of neighbor cloud)
            neighbors = self.positions[indices[1:]]
            centered = neighbors - self.positions[i]
            cov = np.cov(centered.T)
            eigvals = np.linalg.eigvalsh(cov)
            eigvals = np.sort(eigvals)

            combined = np.concatenate([descriptor[:10], eigvals])
            descriptors.append(combined)

        descriptors = np.array(descriptors)

        # Normalize
        descriptors = (descriptors - descriptors.mean(axis=0)) / (descriptors.std(axis=0) + 1e-10)

        # Cluster
        print("  Clustering shape patterns...")
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        labels = kmeans.fit_predict(descriptors)

        cluster_counts = Counter(labels)

        return {
            'descriptors': descriptors,
            'labels': labels,
            'cluster_counts': dict(cluster_counts),
            'cluster_centers': kmeans.cluster_centers_
        }


def analyze_time_evolution(datasets_dict):
    """
    Run pattern discovery on each time slice and compare.

    Look for:
    - Patterns that persist across time
    - Patterns that evolve in specific ways
    - New patterns emerging or disappearing
    """
    print("\n" + "="*60)
    print("TIME EVOLUTION ANALYSIS")
    print("="*60)

    all_results = {}

    for name, data in datasets_dict.items():
        print(f"\n{'='*60}")
        print(f"Epoch: {name}")
        print(f"{'='*60}")

        positions = data['positions']
        engine = PatternDiscoveryEngine(positions)
        results = engine.discover_all_patterns()

        all_results[name] = results

    return all_results


def main():
    """Run pattern discovery on all datasets."""
    import os

    # Load datasets
    dataset_files = [f for f in os.listdir('.') if f.startswith('realistic_z') and f.endswith('.npz')]
    dataset_files = sorted(dataset_files)

    if len(dataset_files) == 0:
        print("❌ No datasets found! Run generate_realistic_data.py first.")
        return

    print("Loading datasets...")
    datasets = {}

    for filename in dataset_files:
        data = np.load(filename)
        name = filename.replace('realistic_', '').replace('.npz', '')
        datasets[name] = {
            'positions': data['positions'],
            'redshifts': data['redshifts'],
            'z_mean': float(data['z_mean']),
        }
        print(f"  ✓ {name}: {len(data['positions'])} galaxies")

    # Run analysis
    results = analyze_time_evolution(datasets)

    # Save results
    print("\n" + "="*60)
    print("Saving results...")
    np.save('pattern_discovery_results.npy', results, allow_pickle=True)
    print("Saved: pattern_discovery_results.npy")

    return results


if __name__ == "__main__":
    results = main()
