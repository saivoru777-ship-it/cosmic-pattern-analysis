# cosmic-pattern-analysis

Started as an attempt to find Einstein "hat" tile signatures in galaxy distributions.
The idea: if large-scale cosmic structure has aperiodic, non-random geometry, it should
show up in clustering statistics across scales. That hypothesis didn't hold — the patterns
are consistent with standard ΛCDM — but the validation process forced me to build something
more useful: a general multiscale clustering detector for sparse 3D point data.

The core question the tool answers: given a cloud of 3D points, are they randomly placed
or do they have hierarchical clustering structure? It works on galaxies, but there's no
reason it couldn't work on seismic data, neural spike locations, disease cases, or anything
else that's a sparse point process in 3D space.

---

## How it works

The pipeline takes positions of points in a 3D box, grids them using cloud-in-cell
interpolation, then computes variance and skewness of counts across 12 logarithmic scales
(from ~3 to 32 voxels). Those curves get compared against a mock distribution using a
covariance-aware χ² test with Ledoit-Wolf style shrinkage regularization.

Two metrics ended up in the final version: **counts variance** and **counts skewness**.
A third metric (Euler characteristic) got dropped — there was a formula bug that made it
return 1 regardless of input, and fixing it wasn't worth it when the other two already work.

The shrinkage parameter (λ=0.1) keeps the covariance matrix condition number around 10³.
With only 5 mocks and 12 scales the matrix blows up without regularization.

---

## Files

```
test_multiscale_production.py    main pipeline, start here
generate_realistic_data.py       ΛCDM-like synthetic data
generate_universes.py            synthetic universe variations
download_illustris.py            pull Illustris simulation data
download_sdss_data.py            pull SDSS galaxy catalog
visualize_3d_patterns.py         3D structure plots
analyze_patterns.py              general pattern analysis
voronoi_deep_dive.py             Voronoi tessellation analysis
test_voronoi_hypotheses.py       4 specific geometric hypothesis tests
```

The `test_multiscale_*.py` files are earlier iterations — `test_multiscale_production.py`
is the one that actually validates correctly.

---

## Quick start

Validate the pipeline first (~2 min on a laptop):

```bash
python test_multiscale_production.py --validate
```

Expected output:

```
Control 1 (unclustered vs unclustered): 0/2 false positives
Control 2 (clustered vs unclustered):   2/2 detections
VALIDATION PASSED
```

Test your own data against mocks:

```bash
python test_multiscale_production.py \
  --data your_data.npz \
  --mocks mock_1.npz mock_2.npz mock_3.npz \
  --n-mocks 50
```

Input format: `.npz` with a `positions` key (N×3 array in Mpc/h) and optionally `box_size`.
The box sizes need to match between data and mocks — if they don't, the χ² values will be
enormous and meaningless (learned that the hard way).

---

## Where things stand

The pipeline validates cleanly. The last real run compared a z=0.1 galaxy sample against
5 Illustris mocks and got χ² ~350,000 — which looks dramatic but is just a box size
mismatch (data at 100 Mpc/h, mocks at 75 Mpc/h). The "anomaly" was a systematic.

Next step is regenerating mocks at the correct box size and running a clean comparison.

There's also an open question in the Voronoi analysis: every dataset tested — random,
clustered, and real — consistently shows around 25 neighbors per Voronoi cell. The
expected value for a 3D Poisson process is ~15.5. Probably a boundary condition artifact
in the tessellation, but hasn't been tracked down yet.

---

## Dependencies

```bash
pip install numpy scipy matplotlib
```

Optional, for downloading real data:

```bash
pip install requests astropy
```

---

## Background

The Einstein "hat" monotile was published in 2023 — a single shape that tiles the plane
without periodic repetition. The idea of looking for similar aperiodic structure in the
cosmic web was worth testing. It didn't pan out, but the validation-first approach to
building the detection pipeline turned out to be the more interesting result.
