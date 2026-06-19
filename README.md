# ADHD Functional Brain Network Analysis

Graph-theoretic analysis of resting-state fMRI connectivity to investigate whether ADHD is associated with altered small-world network topology.

---

## Overview

This project applies graph theory to functional connectivity matrices derived from resting-state fMRI to compare large-scale brain network organization between children with ADHD and healthy controls. Five graph metrics are computed across six proportional thresholds, with permutation testing used to assess group differences — both with and without regression of mean functional connectivity.

**Key finding:** ADHD is associated with shorter characteristic path lengths and reduced clustering at higher network densities, suggesting a shift toward more random (rather than regular) network topology — not a breakdown of network architecture.

---

## Repository Structure

```
.
├── project_ADHD_Folasewa.m        # Main analysis script
├── ADHD_connectivity.mat          # Preprocessed ADHD connectivity matrices (200×200×831)
├── Control_connectivity.mat       # Preprocessed control connectivity matrices (200×200×382)
├── report/
│   └── Network_Neuroscience_Mini_Project_Folasewa.pdf
└── README.md
```

> **Note:** `ADHD_connectivity.mat` and `Control_connectivity.mat` are not tracked in this repository due to file size. See [Data](#data) for instructions on how to reproduce them.

---

## Dependencies

- **MATLAB** (R2021a or later recommended)
- **Brain Connectivity Toolbox (BCT)** — must be on your MATLAB path
  - Download: https://sites.google.com/site/bctnet/
  - BCT functions used: `threshold_proportional`, `clustering_coef_wu`, `weight_conversion`, `distance_wei`, `charpath`, `efficiency_wei`, `randmio_und_connected`

---

## Data

Raw data is sourced from the **Healthy Brain Network (HBN)** dataset, released as part of the WiDS 2025 Datathon on Kaggle:

> https://www.kaggle.com/competitions/widsdatathon2025/data

### Preprocessing (commented out in script)

The raw dataset contains flattened pairwise connectivity vectors (Pearson correlation, 19,901 features per participant) alongside clinical labels. Preprocessing steps:

1. Merge connectivity data with ADHD outcome labels on `participant_id`
2. Split by `ADHD_Outcome` into patient (n=831) and control (n=382) groups
3. Reconstruct each flattened vector into a full 200×200 symmetric connectivity matrix
4. Save as `ADHD_connectivity.mat` (`fmri_ADHD`: 200×200×831) and `Control_connectivity.mat` (`fmri_control`: 200×200×382)

These steps are preserved as commented code at the top of `project_ADHD_Folasewa.m`.

---

## Pipeline

```
Load .mat files
      │
      ▼
For each threshold (5%, 10%, 15%, 20%, 25%, 30%):
      │
      ├── Generate 50 random networks (template: first ADHD subject)
      │       └── Compute reference C_rand and L_rand
      │
      ├── For each subject (ADHD + Control):
      │       ├── Threshold matrix (proportional)
      │       ├── Clustering coefficient (clustering_coef_wu)
      │       ├── Characteristic path length (charpath)
      │       ├── Global efficiency (efficiency_wei)
      │       ├── Local efficiency (efficiency_wei, local)
      │       └── Small-worldness σ = (C/C_rand) / (L/L_rand)
      │
      ├── Regress out mean FC from each metric (within-group OLS)
      │
      └── Permutation testing (10,000 permutations, two-tailed)
              ├── With FC regression
              └── Without FC regression
```

---

## Running the Analysis

1. Clone this repository
2. Add BCT to your MATLAB path:
   ```matlab
   addpath('/path/to/BCT')
   ```
3. Place `ADHD_connectivity.mat` and `Control_connectivity.mat` in the working directory (or regenerate from raw data using the commented preprocessing section)
4. Open and run `project_ADHD_Folasewa.m`

Results are printed to the MATLAB console, organized by threshold and condition (with/without FC regression).

---

## Results Summary

| Metric | Finding |
|---|---|
| Characteristic path length | Significantly shorter in ADHD across **all** thresholds (post-FC regression) |
| Clustering coefficient | Significantly lower in ADHD at 25% and 30% thresholds |
| Global efficiency | Significantly higher in ADHD at 5% and 10% thresholds |
| Local efficiency | No significant group differences at any threshold |
| Small-worldness (σ) | Significant only at 30% threshold (p = 0.0241); both groups exhibit σ > 1 |

Full results with p-values are reported in the [project report](report/Network_Neuroscience_Mini_Project_Folasewa.pdf).

---

## Methods Notes

**Proportional thresholding** was applied at six densities (5%–30%) to ensure consistent edge counts across subjects and thresholds.

**Random network generation:** Rather than generating random networks per subject (computationally prohibitive at this sample size), 50 random networks were generated from a single template graph per threshold using `randmio_und_connected`. Reference clustering and path length values were averaged across these 50 networks to compute σ.

**Mean FC regression:** To control for individual differences in overall connectivity strength, mean FC was regressed out of each graph metric within each group using OLS, with group means added back to the residuals before group comparison.

**Permutation testing:** Group differences were assessed with 10,000 permutations (random label shuffling), with two-tailed p-values computed against the resulting null distribution.

---

## Citation

If you use or build on this work, please cite the data source:

> Healthy Brain Network (HBN) / WiDS Datathon 2025. https://www.kaggle.com/competitions/widsdatathon2025/data

---

## Author

**Folasewa** — Network Neuroscience Mini Project
