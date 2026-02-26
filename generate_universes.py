"""
Generate synthetic "universe" patterns for testing tiling hypothesis.

We'll create three types of patterns:
1. Random Gaussian field (ΛCDM-like initial conditions)
2. Gravitationally evolved field (approximate structure formation)
3. Aperiodic tiling-influenced pattern (what the hypothesis predicts)
"""

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt


class UniverseGenerator:
    """Generate different synthetic cosmic structure patterns."""

    def __init__(self, grid_size=512, box_size=1000):
        """
        Parameters:
        - grid_size: number of grid points per dimension
        - box_size: physical size in Mpc/h
        """
        self.grid_size = grid_size
        self.box_size = box_size
        self.dx = box_size / grid_size

    def generate_random_field(self, power_law_index=-2.5, seed=None):
        """
        Generate random Gaussian field with power-law spectrum.
        This mimics ΛCDM initial conditions.

        P(k) ~ k^n where n ~ -2 to -3 for CDM
        """
        if seed is not None:
            np.random.seed(seed)

        # Create k-space grid
        k = np.fft.fftfreq(self.grid_size, d=self.dx)
        kx, ky = np.meshgrid(k, k)
        k_mag = np.sqrt(kx**2 + ky**2)
        k_mag[0, 0] = 1  # Avoid division by zero

        # Power spectrum P(k) ~ k^n
        P_k = k_mag**(power_law_index)

        # Generate random phases
        phase = 2 * np.pi * np.random.random((self.grid_size, self.grid_size))

        # Create field in Fourier space
        field_k = np.sqrt(P_k) * np.exp(1j * phase)

        # Transform to real space
        field = np.fft.ifft2(field_k).real

        # Normalize
        field = (field - field.min()) / (field.max() - field.min())

        return field

    def evolve_gravity(self, initial_field, evolution_factor=2.0):
        """
        Approximate gravitational evolution.
        High-density regions get denser (rich get richer).
        """
        # Simple approximation: amplify high-density regions
        field = initial_field.copy()

        # Apply non-linear evolution (overdensities grow faster)
        field = field**evolution_factor

        # Smooth slightly to represent gravitational smoothing
        field = gaussian_filter(field, sigma=2)

        # Normalize
        field = (field - field.min()) / (field.max() - field.min())

        return field

    def generate_fibonacci_spiral_points(self, n_points=2000):
        """
        Generate points using Fibonacci/golden spiral.
        This creates hierarchical structure with golden ratio scaling.
        """
        phi = (1 + np.sqrt(5)) / 2  # Golden ratio

        points = []
        for i in range(n_points):
            # Fibonacci spiral
            theta = 2 * np.pi * i / phi
            r = np.sqrt(i / n_points) * self.grid_size / 2

            x = self.grid_size / 2 + r * np.cos(theta)
            y = self.grid_size / 2 + r * np.sin(theta)

            if 0 <= x < self.grid_size and 0 <= y < self.grid_size:
                points.append([x, y])

        return np.array(points)

    def generate_substitution_tiling_pattern(self, n_iterations=5, seed=None):
        """
        Generate a pattern inspired by substitution tiling rules.
        Uses hierarchical structure with scale ratios.
        """
        if seed is not None:
            np.random.seed(seed)

        field = np.zeros((self.grid_size, self.grid_size))

        phi = (1 + np.sqrt(5)) / 2  # Golden ratio

        # Start with base scale
        base_scale = self.grid_size / 20

        # Hierarchical levels with golden ratio scaling
        for level in range(n_iterations):
            scale = base_scale * (phi ** level)
            n_centers = int(200 / (phi ** level))  # Fewer centers at larger scales

            # Random centers
            centers_x = np.random.randint(0, self.grid_size, n_centers)
            centers_y = np.random.randint(0, self.grid_size, n_centers)

            # Add Gaussian peaks at each center
            for cx, cy in zip(centers_x, centers_y):
                y, x = np.ogrid[:self.grid_size, :self.grid_size]
                distance = np.sqrt((x - cx)**2 + (y - cy)**2)
                peak = np.exp(-(distance**2) / (2 * scale**2))
                field += peak * (1.0 / (level + 1))  # Weight by level

        # Add connectivity (filaments) between peaks
        # Smooth to create web-like structure
        field = gaussian_filter(field, sigma=3)

        # Normalize
        field = (field - field.min()) / (field.max() - field.min())

        return field

    def generate_penrose_inspired_pattern(self, seed=None):
        """
        Generate pattern inspired by Penrose tiling structure.
        Uses golden ratio subdivisions.
        """
        if seed is not None:
            np.random.seed(seed)

        field = np.zeros((self.grid_size, self.grid_size))
        phi = (1 + np.sqrt(5)) / 2

        # Create rhombus-like structures at multiple scales
        scales = [self.grid_size / (phi**i) for i in range(2, 8)]

        for scale in scales:
            n_structures = int(self.grid_size**2 / scale**2 / 2)

            for _ in range(n_structures):
                cx = np.random.randint(0, self.grid_size)
                cy = np.random.randint(0, self.grid_size)
                angle = np.random.random() * 2 * np.pi

                # Create elongated structure (filament-like)
                y, x = np.ogrid[:self.grid_size, :self.grid_size]
                dx = x - cx
                dy = y - cy

                # Rotate
                dx_rot = dx * np.cos(angle) - dy * np.sin(angle)
                dy_rot = dx * np.sin(angle) + dy * np.cos(angle)

                # Anisotropic Gaussian (creates filament)
                structure = np.exp(-(dx_rot**2) / (2 * scale**2) -
                                  (dy_rot**2) / (2 * (scale/phi)**2))

                field += structure

        field = gaussian_filter(field, sigma=2)
        field = (field - field.min()) / (field.max() - field.min())

        return field


def main():
    """Generate all patterns and save."""
    print("Generating synthetic universes...")

    gen = UniverseGenerator(grid_size=512, box_size=1000)

    # Generate different patterns
    print("1. Generating random Gaussian field (ΛCDM-like)...")
    random_field = gen.generate_random_field(seed=42)

    print("2. Generating gravitationally evolved field...")
    evolved_field = gen.evolve_gravity(random_field, evolution_factor=1.5)

    print("3. Generating substitution-tiling pattern...")
    tiling_field = gen.generate_substitution_tiling_pattern(seed=42)

    print("4. Generating Penrose-inspired pattern...")
    penrose_field = gen.generate_penrose_inspired_pattern(seed=42)

    # Save all fields
    np.save('random_field.npy', random_field)
    np.save('evolved_field.npy', evolved_field)
    np.save('tiling_field.npy', tiling_field)
    np.save('penrose_field.npy', penrose_field)

    # Visualize
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))

    axes[0, 0].imshow(random_field, cmap='viridis', origin='lower')
    axes[0, 0].set_title('Random Gaussian Field\n(ΛCDM initial conditions)', fontsize=12)
    axes[0, 0].axis('off')

    axes[0, 1].imshow(evolved_field, cmap='viridis', origin='lower')
    axes[0, 1].set_title('Gravitationally Evolved Field\n(Approximate ΛCDM structure)', fontsize=12)
    axes[0, 1].axis('off')

    axes[1, 0].imshow(tiling_field, cmap='viridis', origin='lower')
    axes[1, 0].set_title('Substitution Tiling Pattern\n(Hypothesis: hierarchical w/ golden ratio)', fontsize=12)
    axes[1, 0].axis('off')

    axes[1, 1].imshow(penrose_field, cmap='viridis', origin='lower')
    axes[1, 1].set_title('Penrose-Inspired Pattern\n(Aperiodic structure)', fontsize=12)
    axes[1, 1].axis('off')

    plt.tight_layout()
    plt.savefig('synthetic_universes.png', dpi=150, bbox_inches='tight')
    print("\nSaved visualization: synthetic_universes.png")
    print("Saved data files: random_field.npy, evolved_field.npy, tiling_field.npy, penrose_field.npy")

    return random_field, evolved_field, tiling_field, penrose_field


if __name__ == "__main__":
    main()
