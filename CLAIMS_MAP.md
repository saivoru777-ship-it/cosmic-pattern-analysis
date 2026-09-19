# Claims Map: cosmic-pattern-analysis vs biospatial-regimes

## Purpose
Map every proposed claim for both papers against: (1) prior literature ownership,
(2) overlap between manuscripts, (3) required experiment, (4) falsifiability.

---

## BIOSPATIAL-REGIMES: What it should uniquely own

### Claim B1: VMR descriptors define a structured regime space for biological spatial mechanisms
- **Prior art:** Detto & Muller-Landau 2013 (scalewise variance for process inference); Greig-Smith 1952 (multi-scale quadrats); Baddeley et al. 2015 (spatial point pattern methodology)
- **What's genuinely new:** The specific 6-descriptor extraction from VMR ratio curves + PCA embedding across 7 biological mechanisms. The descriptor set itself (peak height, threshold scale, width, slopes, integrated excess) is a concrete contribution, even if the underlying variance-vs-scale idea is old.
- **Overlap with cosmic:** If cosmic also builds a descriptor-based classifier, they collide.
- **Required fix:** Cite Detto & Muller-Landau, Greig-Smith, Torquato. Reframe as "we operationalize the scalewise-variance idea as a specific descriptor set for biological point processes" rather than "we introduce VMR fingerprints."
- **Status:** KEEP as biospatial's core contribution. Must add ~15-20 references.

### Claim B2: Clustering-positive mechanisms are separable; sub-Poisson mechanisms are degenerate
- **Prior art:** Well-known that VMR ≈ 1 for both Poisson and regular processes. The sub-Poisson degeneracy is expected.
- **What's genuinely new:** The specific demonstration with these 7 mechanisms + parameter sweeps showing exactly WHERE separation holds and fails.
- **Required fix:** Acknowledge that sub-Poisson degeneracy is expected from VMR's mathematical definition. The finding is "here's the exact boundary," not "we discovered a limitation."
- **Falsifiability:** Would be falsified if adding anisotropy or pair correlation descriptors resolves the null cloud. SHOULD TEST THIS.

### Claim B3: Empirical fingerprints land in biologically interpretable regions
- **Prior art:** None directly — this specific projection is new.
- **What's genuinely new:** Yes, this is novel.
- **Required fix:** Restate 61% out-of-regime more carefully. Current: "suggests biology exceeds single-mechanism complexity." Better: "61% lie outside the simulated reference distribution. Possible explanations include mechanism mixtures, multi-parameter interactions, observation noise, domain shift, or model misspecification." Then test which explanation holds.
- **Falsifiability:** Would be weakened if random noise perturbation of simulated fingerprints produces similar out-of-regime rates.

### Claim B4: Mixture recovery fails (errors up to 0.999)
- **Prior art:** N/A — this is a negative result about their own method.
- **What's genuinely new:** Yes. This is valuable and should be PROMOTED, not buried.
- **Required fix:** Feature this as a key result. It tells users: "don't interpret KDE affinities as mixture fractions."

### Claim B5: Silhouette score 0.172 → "well-separated"
- **Prior art:** Standard interpretation: 0.172 is WEAK separation.
- **Required fix:** Stop calling this "well-separated." Report it honestly. The within-clustering-manifold silhouette (excluding the null cloud) would be more informative. Compute and report both.

### Claim B6: Anisotropy descriptor computed but excluded from embedding
- **Issue:** Layered patterns are exactly where anisotropy matters. Excluding it, then calling the layered/random/exclusion collapse a "fundamental VMR limitation" is misleading — it might be an induced limitation.
- **Required fix:** Run the embedding WITH anisotropy. If layered separates, report that and restructure the paper. If it still doesn't, that's a stronger negative result.

---

## COSMIC-PATTERN-ANALYSIS: What it should uniquely own

### The new identity: "Limits and Identifiability of Multiscale Count-Variance Fingerprints"

This paper becomes the STRESS TEST. Not "VMR works" but "where VMR fails, why, and what you need instead."

### Claim C1: VMR(R), g(r), and S(k) carry substantially the same second-order information
- **Prior art:** Torquato & Stillinger 2003; Torquato 2018 review; mathematical statistics textbooks. Number variance is determined by the pair correlation function / structure factor under stationarity.
- **What's new:** Explicit empirical demonstration across diverse point-process families, with quantified redundancy (e.g., mutual information, prediction error when one predicts the other).
- **Required experiment:** For each process family, compute all three. Show that a classifier using VMR alone vs. g(r) alone vs. S(k) alone achieves similar discrimination. Then show where they DIFFER (finite windows, boundary effects, discretization, anisotropy).
- **Falsifiability:** If VMR consistently adds information beyond g(r)+S(k), the "same information" claim fails. That would also be interesting.

### Claim C2: VMR fingerprints cannot uniquely identify generating mechanisms (equifinality)
- **Prior art:** Ecological literature explicitly warns about this (e.g., Wiegand et al. 2021 in Ecological Processes). Different processes → same pattern; same process → different patterns.
- **What's new:** Systematic demonstration across a carefully chosen process zoo, with explicit construction of equifinal pairs (two different mechanisms producing indistinguishable VMR fingerprints).
- **Required experiment:** For each pair of process families, find parameter settings that produce statistically indistinguishable VMR curves. Report the "confusion matrix" of the fingerprint.
- **Falsifiability:** Would be weakened if NO equifinal pairs exist (unlikely, but would strengthen VMR's claim).

### Claim C3: Specific failure modes under observation effects
- **Prior art:** Survey-window literature in cosmology (Beutler et al. 2021); Hawat et al. for finite-sample S(k) estimation.
- **What's genuinely new:** Systematic characterization of how thinning, coordinate error, finite windows, holes/masks, anisotropic windows, variable density, and dimensionality affect VMR fingerprint shape and classification accuracy.
- **Required experiment:** Take known point processes. Apply each degradation at multiple severity levels. Measure classification accuracy vs. degradation. Find the regime where VMR becomes unreliable.
- **Falsifiability:** If VMR is robust to all tested degradations, that's a positive result (strengthens the method). If it fails under specific conditions, that maps the reliability boundary.

### Claim C4: Different forms of aperiodic order produce dramatically different fingerprints
- **Prior art:** Kaplan et al. 2024 (Hat has periodic diffraction); Spectre addendum (non-periodic diffraction); Oğuz et al. 2017 (quasicrystal hyperuniformity).
- **What's genuinely new:** Using VMR + S(k) + persistent homology to distinguish Hat, Spectre, Penrose, and Ammann-Beenker vertex sets — point processes with known but dramatically different long-range order.
- **Required experiment:** Generate actual aperiodic tilings from their mathematical constructions. Compute VMR, S(k), persistence diagrams. Show which tilings are distinguishable and which are degenerate under each summary statistic.
- **Falsifiability:** If VMR can't distinguish any of them, that's a failure mode. If it distinguishes Hat from Spectre despite both being aperiodic, that's interesting. If persistent homology succeeds where VMR fails (or vice versa), that ranks the methods.

### Claim C5: Higher-order statistics (persistent homology, bispectrum) resolve equifinal pairs
- **Prior art:** Persistent homology in cosmology (Feldbrugge et al. 2019, Biagetti et al. 2021); bispectrum as beyond-P(k) in cosmology is well-established.
- **What's new:** Systematic comparison of which equifinal pairs under VMR become distinguishable under higher-order statistics. This is the "what do you need beyond second order?" question.
- **Required experiment:** For each equifinal pair from C2, test whether persistent homology, nearest-neighbor distributions, or higher-order count statistics break the degeneracy.
- **Falsifiability:** If higher-order statistics don't help, the equifinality is "deep." If they do, we've mapped exactly when to upgrade from VMR to more expensive methods.

### Claim C6: Injection-recovery establishes detection limits
- **Prior art:** Standard in cosmology (e.g., BAO reconstruction). Standard in biology (spike-in controls).
- **What's new:** Applied to VMR fingerprint classification — mixture parameter sweep showing detection sensitivity vs. contamination fraction.
- **Required experiment:** Superposition or thinning-based mixture of process families at varying fractions f. Measure classification accuracy vs. f. Report the minimum detectable fraction for each pair.
- **Mathematical note:** Point patterns can't be mixed as (1-f)X_A + f*X_B the way fields can. Must use thinning-plus-superposition or probabilistic point replacement. Define this explicitly.

---

## CLAIMS THAT SHOULD BE DELETED FROM BOTH PAPERS

| Claim | Reason |
|-------|--------|
| "VMR fingerprints are a new method" | Greig-Smith 1952, Detto & Muller-Landau 2013, Torquato & Stillinger 2003 |
| "Flat fingerprint = emergent; peaked = imposed" | Poisson is flat without emergence; self-organized vegetation produces peaks |
| "Universal organizing mechanism across nm to Mpc" | Second-order statistical similarity ≠ shared mechanism |
| Voronoi 25.5 anomaly | Confirmed bug |
| Cosmic evolution frozen by z=0.65 | Independent random seeds, not evolved snapshots |
| Illustris validation | Two synthetic generators compared to each other |
| "mechanism identification" | Replace with "discrimination among specified generative model families" |

---

## CLAIMS THAT NEED VOCABULARY CHANGES

| Current | Better | Why |
|---------|--------|-----|
| "VMR fingerprints identify mechanisms" | "VMR fingerprints discriminate among specified process families" | Equifinality; identification requires exclusion of alternatives |
| "mediation" (in later papers) | "conditional attenuation" or "statistical mediation under explicit assumptions" | Causal mediation requires stronger assumptions than observational association |
| "well-separated" (silhouette 0.172) | "distinguishable within the clustering manifold; degenerate in the null cloud" | 0.172 is weak overall |
| "61% out-of-regime suggests biology is complex" | "61% lie outside the simulated reference distribution" then test explanations | Multiple explanations, not just one |

---

## THE PAPER COSMIC SHOULD BECOME

### Title: "Limits and Identifiability of Multiscale Count-Variance Fingerprints for Spatial Point Patterns"

### Structure:
1. **Introduction:** Scalewise variance has a long history (cite properly). VMR fingerprints operationalize this idea. But what are their actual discriminative limits? When are they redundant with g(r) and S(k)? When do they fail? No systematic stress test exists.

2. **Mathematical framework:** VMR ↔ g(r) ↔ S(k) equivalence under stationarity. What each statistic uniquely captures under finite windows, discretization, and non-stationarity.

3. **Point process benchmark zoo:** Homogeneous Poisson; inhomogeneous Poisson; Thomas cluster; log-Gaussian Cox; Matérn hard-core; perturbed lattice; determinantal; Penrose vertices; Hat vertices; Spectre vertices; Ammann-Beenker vertices. All at matched intensity and window.

4. **Experiment 1 — Redundancy test:** VMR vs. g(r) vs. S(k) classification accuracy. Ablation: what does each add?

5. **Experiment 2 — Equifinality map:** Systematic search for process pairs with indistinguishable VMR but different higher-order structure.

6. **Experiment 3 — Observation degradation:** Thinning, noise, boundary clipping, anisotropy, density variation. Classification accuracy vs. degradation severity.

7. **Experiment 4 — Aperiodic adversarial test:** Hat vs. Spectre vs. Penrose vs. Ammann-Beenker. Which are VMR-distinguishable? Which require S(k) or persistent homology?

8. **Experiment 5 — Injection-recovery:** Detection sensitivity for thinning-plus-superposition mixtures. Minimum detectable fraction per process pair.

9. **Discussion:** What VMR can reliably tell you. What it can't. When to upgrade to higher-order methods. Implications for cross-domain application.

---

## RELATIONSHIP BETWEEN THE TWO PAPERS

| | biospatial-regimes | cosmic-pattern-analysis |
|-|---|---|
| **Identity** | Application atlas | Foundational stress test |
| **Question** | Where do biological mechanisms sit in VMR descriptor space? | What information does VMR actually contain, and where does it fail? |
| **Process families** | 7 biological mechanisms | 11+ mathematical process families (including aperiodic) |
| **Novel contribution** | Specific descriptor set + biological regime projections | Equifinality map + redundancy quantification + degradation limits |
| **Key result type** | "Here's where your data lands" | "Here's when to trust where it lands" |
| **Dependency** | Can stand alone, but stronger with cosmic as foundation | Makes biospatial's claims defensible |

Cosmic should be published FIRST or simultaneously. It provides the methodological foundation that biospatial's claims rest on.

---

## NON-NEGOTIABLE CITATIONS TO ADD

### For cosmic:
- Greig-Smith 1952 (multi-scale quadrats)
- Baddeley, Rubak & Turner 2015 (spatial point patterns methodology)
- Detto & Muller-Landau 2013 (scalewise variance)
- Torquato & Stillinger 2003 (number variance, hyperuniformity)
- Torquato 2018 (hyperuniformity review)
- Hawat et al. 2022 (structure factor estimation from finite samples)
- Myllymäki et al. 2017 (global envelope tests)
- Smith et al. 2024 Combinatorial Theory (Hat tile — use peer-reviewed, not arXiv)
- Kaplan, O'Keeffe & Treacy 2024 (Hat diffraction is periodic)
- Oğuz et al. 2017 (quasicrystal hyperuniformity)
- Wiegand et al. / Velázquez et al. 2016 (equifinality in spatial patterns)
- Functional point-pattern classification (Mateu & Schoenberg 2023 or similar)
- Feldbrugge et al. 2019 or Biagetti et al. 2021 (persistent homology in cosmology)

### For biospatial-regimes:
- All of the above foundational references
- Diggle 2003 (already cited)
- Ripley 1977 (already cited)
- Illian et al. 2008 (Statistical Analysis and Modelling of Spatial Point Patterns)
- Wiegand & Moloney 2014 (Handbook of Spatial Point-Pattern Analysis in Ecology)
