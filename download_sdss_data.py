"""
Download real SDSS galaxy data for pattern analysis.

We'll use the SDSS SkyServer SQL API to get:
- Galaxy positions (RA, Dec, redshift)
- Multiple redshift slices for time evolution
- Spectroscopic data (more reliable than photometric)
"""

import numpy as np
import requests
import time
import json
from urllib.parse import urlencode


class SDSSDownloader:
    """Download galaxy data from SDSS."""

    def __init__(self):
        self.base_url = "http://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"

    def build_query(self, z_min, z_max, max_galaxies=50000):
        """
        Build SQL query for SDSS galaxies in redshift range.

        Parameters:
        - z_min, z_max: redshift range
        - max_galaxies: maximum number to download
        """
        query = f"""
        SELECT TOP {max_galaxies}
            p.objid,
            p.ra,
            p.dec,
            s.z as redshift,
            s.zErr as redshift_error,
            p.petroMag_r as r_mag,
            p.type
        FROM PhotoObj AS p
        JOIN SpecObj AS s ON s.bestobjid = p.objid
        WHERE
            s.z BETWEEN {z_min} AND {z_max}
            AND s.zErr < 0.01
            AND s.zWarning = 0
            AND p.type = 3
            AND p.petroMag_r < 18.0
            AND p.clean = 1
        ORDER BY p.ra
        """
        return query

    def execute_query(self, query, format='csv'):
        """Execute SQL query on SDSS server."""
        params = {
            'cmd': query,
            'format': format
        }

        print(f"Executing query...")
        print(f"URL: {self.base_url}")

        try:
            response = requests.get(self.base_url, params=params, timeout=300)
            response.raise_for_status()
            return response.text
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")
            return None

    def download_redshift_slice(self, z_min, z_max, max_galaxies=50000):
        """Download galaxies in specific redshift range."""
        print(f"\nDownloading galaxies: z = {z_min:.2f} to {z_max:.2f}")
        print(f"Maximum galaxies: {max_galaxies}")

        query = self.build_query(z_min, z_max, max_galaxies)
        result = self.execute_query(query)

        if result is None:
            print("Download failed!")
            return None

        # Parse CSV
        lines = result.strip().split('\n')

        if len(lines) < 2:
            print("No data returned!")
            return None

        # Get headers
        headers = lines[0].split(',')
        print(f"Headers: {headers}")

        # Parse data
        data = []
        for line in lines[1:]:
            if line.strip():
                values = line.split(',')
                if len(values) == len(headers):
                    data.append(values)

        if len(data) == 0:
            print("No galaxies found in this range!")
            return None

        # Convert to numpy array
        data_array = np.array(data)

        print(f"Downloaded {len(data)} galaxies")

        return {
            'headers': headers,
            'data': data_array,
            'z_min': z_min,
            'z_max': z_max,
            'n_galaxies': len(data)
        }

    def convert_to_comoving_coordinates(self, ra, dec, z):
        """
        Convert RA, Dec, redshift to comoving 3D coordinates.

        Using flat ΛCDM cosmology:
        - H0 = 70 km/s/Mpc
        - Ωm = 0.3
        - ΩΛ = 0.7

        This is approximate but good enough for pattern analysis.
        """
        from scipy.integrate import quad

        # Hubble constant
        H0 = 70.0  # km/s/Mpc
        c = 299792.458  # km/s

        # Cosmological parameters
        Om = 0.3
        OL = 0.7

        # Comoving distance integral
        def E(z):
            return np.sqrt(Om * (1 + z)**3 + OL)

        # For vectorization
        if np.isscalar(z):
            z_array = np.array([z])
        else:
            z_array = np.array(z)

        # Compute comoving distance for each redshift
        dc = np.zeros_like(z_array, dtype=float)
        for i, zi in enumerate(z_array):
            if zi > 0:
                integral, _ = quad(lambda x: 1/E(x), 0, zi)
                dc[i] = (c / H0) * integral
            else:
                dc[i] = 0

        # Convert RA, Dec to radians
        ra_rad = np.deg2rad(ra)
        dec_rad = np.deg2rad(dec)

        # Convert to Cartesian coordinates
        x = dc * np.cos(dec_rad) * np.cos(ra_rad)
        y = dc * np.cos(dec_rad) * np.sin(ra_rad)
        z = dc * np.sin(dec_rad)

        return x, y, z


def download_time_slices():
    """
    Download multiple redshift slices for time evolution analysis.

    Redshift ranges:
    - z ~ 0.05-0.1: Local universe (today)
    - z ~ 0.2-0.3: ~2.5 billion years ago
    - z ~ 0.4-0.5: ~4.5 billion years ago
    - z ~ 0.6-0.7: ~6 billion years ago
    """
    downloader = SDSSDownloader()

    time_slices = [
        {'name': 'z0.05-0.10', 'z_min': 0.05, 'z_max': 0.10, 'max': 50000},
        {'name': 'z0.20-0.25', 'z_min': 0.20, 'z_max': 0.25, 'max': 50000},
        {'name': 'z0.40-0.45', 'z_min': 0.40, 'z_max': 0.45, 'max': 50000},
        {'name': 'z0.60-0.65', 'z_min': 0.60, 'z_max': 0.65, 'max': 30000},
    ]

    results = {}

    for slice_info in time_slices:
        print(f"\n{'='*60}")
        print(f"TIME SLICE: {slice_info['name']}")
        print(f"{'='*60}")

        data = downloader.download_redshift_slice(
            slice_info['z_min'],
            slice_info['z_max'],
            slice_info['max']
        )

        if data is not None:
            results[slice_info['name']] = data

            # Save to file
            filename = f"sdss_{slice_info['name']}.npz"
            np.savez_compressed(
                filename,
                headers=data['headers'],
                data=data['data'],
                z_min=data['z_min'],
                z_max=data['z_max'],
                n_galaxies=data['n_galaxies']
            )
            print(f"Saved: {filename}")

        # Be nice to the server
        time.sleep(2)

    return results


def convert_all_to_3d():
    """Convert all downloaded slices to 3D comoving coordinates."""
    import os

    downloader = SDSSDownloader()

    slice_files = [
        'sdss_z0.05-0.10.npz',
        'sdss_z0.20-0.25.npz',
        'sdss_z0.40-0.45.npz',
        'sdss_z0.60-0.65.npz',
    ]

    print("\n" + "="*60)
    print("CONVERTING TO 3D COMOVING COORDINATES")
    print("="*60)

    for filename in slice_files:
        if not os.path.exists(filename):
            print(f"Skipping {filename} (not found)")
            continue

        print(f"\nProcessing: {filename}")

        # Load data
        loaded = np.load(filename, allow_pickle=True)
        data = loaded['data']
        headers = loaded['headers']

        # Find column indices
        headers_list = list(headers)
        ra_idx = headers_list.index('ra')
        dec_idx = headers_list.index('dec')
        z_idx = headers_list.index('redshift')

        # Extract coordinates
        ra = data[:, ra_idx].astype(float)
        dec = data[:, dec_idx].astype(float)
        z = data[:, z_idx].astype(float)

        print(f"  Galaxies: {len(ra)}")
        print(f"  RA range: {ra.min():.2f} to {ra.max():.2f}")
        print(f"  Dec range: {dec.min():.2f} to {dec.max():.2f}")
        print(f"  z range: {z.min():.4f} to {z.max():.4f}")

        # Convert to 3D
        print(f"  Converting to comoving coordinates...")
        x, y, z_coord = downloader.convert_to_comoving_coordinates(ra, dec, z)

        # Save 3D positions
        output_file = filename.replace('.npz', '_3d.npz')
        np.savez_compressed(
            output_file,
            x=x,
            y=y,
            z=z_coord,
            ra=ra,
            dec=dec,
            redshift=z
        )

        print(f"  Saved: {output_file}")
        print(f"  3D extent:")
        print(f"    x: {x.min():.1f} to {x.max():.1f} Mpc/h")
        print(f"    y: {y.min():.1f} to {y.max():.1f} Mpc/h")
        print(f"    z: {z_coord.min():.1f} to {z_coord.max():.1f} Mpc/h")


def main():
    """Main download workflow."""
    print("="*60)
    print("SDSS GALAXY DATA DOWNLOADER")
    print("="*60)
    print("\nThis will download real galaxy data from SDSS DR17")
    print("for pattern analysis across cosmic time.\n")

    # Download data
    print("Step 1: Downloading redshift slices...")
    results = download_time_slices()

    if len(results) == 0:
        print("\n⚠️  No data downloaded! Check internet connection and SDSS server status.")
        print("Continuing with conversion of any existing data...")

    # Convert to 3D
    print("\nStep 2: Converting to 3D coordinates...")
    convert_all_to_3d()

    print("\n" + "="*60)
    print("DOWNLOAD COMPLETE!")
    print("="*60)
    print("\nDownloaded files:")
    import os
    for f in os.listdir('.'):
        if f.startswith('sdss_') and f.endswith('.npz'):
            print(f"  - {f}")


if __name__ == "__main__":
    main()
