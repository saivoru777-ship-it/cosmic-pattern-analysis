# Archived: cosmic-pattern-analysis V1

**Status:** Archived exploratory project. Do not submit.

## What this project was

An exploratory investigation of whether aperiodic tiling (Einstein hat monotile)
signatures exist in cosmic large-scale structure. The project generated synthetic
ΛCDM-like galaxy distributions and compared them against tiling-inspired
alternatives using power spectra, persistent homology, and multiscale
count-variance statistics.

## Why it was archived

A thorough review (September 2026) identified critical issues:

1. **Voronoi 25.5 anomaly**: Confirmed as a bug. `test_voronoi_hypotheses.py`
   already documents that `len(region)` counts vertices, not neighbors. The
   corrected `ridge_points` approach converges to ~15.5 as expected.

2. **Cosmic evolution result**: Each redshift uses `seed=42+i` — four independent
   random draws, not snapshots of one evolving universe. The `-0.030 bits`
   entropy stability is trivially expected from matched generators.

3. **Code/manuscript mismatches**: The paper claims 3D H0/H1/H2 persistent
   homology; the code runs 2D H0/H1. The paper claims empirical Monte Carlo
   p-values; the code uses `chi2.cdf`. The paper claims Eisenstein-Hu transfer
   function; the code uses `k^n * exp(-(k/k_cut)^2)`.

4. **"Illustris" validation**: `download_illustris.py` generates another custom
   synthetic field, not actual Illustris data. Two custom generators compared
   to each other.

5. **Engineered tiling controls**: The substitution-tiling generator explicitly
   uses `scale = base_scale * phi^level` — golden ratio baked in by construction.

6. **Conceptual flaw**: The Hat tile's vertex diffraction can show periodic
   (sixfold) structure (Kaplan et al. 2024). "Aperiodic → universal Fourier
   signature" is not generally valid.

## What emerged from the failure

The failure exposed a deeper methodological question: what information does
multiscale count variance actually contain, and where does it fail to distinguish
generating processes? This question became a new project:

**spatial-fingerprint-identifiability** — a systematic characterization of the
discriminative limits, degeneracies, and failure modes of multiscale
count-variance fingerprints across diverse point-process families.

## Files

- `paper/main.tex` — the original manuscript (DO NOT SUBMIT)
- `CLAIMS_MAP.md` — the novelty/claims analysis that motivated the archive
- All `.py` scripts and `.npy` data — preserved for provenance
