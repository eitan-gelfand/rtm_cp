# Altered use of stimulus history in face processing in congenital prosopagnosia

## Overview

This repository contains the analysis workspace for the congenital prosopagnosia (CP) face-discrimination study described in the manuscript *Altered use of stimulus history in face processing in congenital prosopagnosia*.

The pipeline examines how face discrimination is shaped by:

- diagnostic group (`TD` vs `CP`)
- face race (`Own-Race` vs `Other-Race`)
- regression-to-the-mean condition (`biasp` vs `biasm`)
- age
- recent and accumulated stimulus history

The repository includes:

- Weibull psychometric fitting
- d-prime / signal-detection analysis
- trial-level GLMM analysis
- bootstrap follow-up for peak-age estimates
- histogram and range-based descriptive plots
- RTM correlation analyses
- updating-model analysis

The manuscript contains the full theoretical and methodological interpretation. The README is intended only as a compact guide to the repository structure and execution workflow.

**Manuscript Contact:** Bat-sheva Hadad - bhadad@edu.haifa.ac.il
**Repository / Analysis Contact:** Eitan Gelfand - eitan.gelfand@gmail.com

---

## Quick Start

### 1. Restore the R environment

This project uses `renv` for R package version control.

```r
renv::restore()
```

### 2. Install Python dependencies

The d-prime workflow is implemented in Python.

```bash
pip install -r requirements.txt
```

### 3. Font setup

Shared plotting helpers are located in `R/`. Most figure scripts use:

```r
source("R/setup_fonts.R")
setup_fonts()
```

The helper selects a Times-like serif font when available and otherwise falls back to a system serif font.

### 4. Execution order

Some analyses are independent, but the main derived-metric workflow is:

1. Run the Weibull fitting script.
2. Run the d-prime calculation script.
3. Run downstream model-based and descriptive analyses that use those derived outputs.

---

## File Organization

### Main directories

- `weibull_analysis/` - subject-level psychometric fitting and Weibull-derived metrics
- `d-prime/` - workflow for d-prime, criterion, and related summary outputs
- `glmm_analysis/` - trial-level GLMM analysis and bootstrap follow-up
- `rtm_corr/` - correlation analyses linking RTM-related measures
- `histograms/` - grouped summary plots for d-prime and PSE
- `range_analysis/` - accuracy-by-range descriptive plots
- `trail_updating_analysis/` - updating-model analysis of recent vs accumulated history
- `R/` - shared helpers for paths, fonts, and plotting theme
- `data/` - input and derived CSV files used by the analysis scripts
- `output/` - generated figures and summary outputs


### Shared R helpers

- `R/paths.R` - centralized project paths
- `R/theme_pub.R` - shared publication plotting theme
- `R/setup_fonts.R` - font setup for figures

---

## Main Scripts

### Derived metrics

- `weibull_analysis/Weibull_CP_experiment.Rmd`
  - Fits Weibull functions by subject and condition.
  - Generates psychometric parameters and prediction outputs used in downstream summaries.

- `d-prime/dprime_calculation.py`
  - Computes d-prime, criterion, and related signal-detection outputs.
  - Writes derived CSV files used by downstream scripts.

### Main model-based analyses

- `glmm_analysis/GLMM_CP.R`
  - Main trial-level GLMM of accuracy.
  - Tests effects of Group, Race, Regression condition, and age.

- `glmm_analysis/GLMM_CP_bootstrap.R`
  - Bootstrap follow-up to the main GLMM.
  - Evaluates stability of peak-age estimates from the modeled trajectories.

- `trail_updating_analysis/last_trail_updating_analysis.R`
  - Compares updating models based on the long-term average (`t-inf`) and previous trial (`t-1`).

### Follow-up and descriptive analyses

- `hisotgrams/dprime_histograms.R`
  - Subject-level d-prime grouped summaries by Group, Race, and Regression condition.

- `hisotgrams/dprime_young_sub_hist.R`
  - Younger-subgroup version of the d-prime grouped summary.

- `hisotgrams/pse_histograms.R`
  - Grouped PSE summaries based on the Weibull outputs.

- `hisotgrams/dprime_cfmt_group.R`
  - d-prime grouped summary split by CFMT subgroup.

- `range_analysis/accXrange_CP_bias.R`
  - Accuracy as a function of morph range, grouped by Race within Bias condition.

- `range_analysis/accXrange_CP_race.R`
  - Accuracy as a function of morph range, grouped by Bias condition within Race.

- `rtm_corr/rtm_correlation_cp.R`
  - Correlates RTM magnitude with Bias- accuracy.

- `rtm_corr/dprime_rtm_corr.R`
  - Correlates RTM magnitude with d-prime.

---

## Workflow

### Stage 1: Generate derived metrics

Run:

- `weibull_analysis/Weibull_CP_experiment.Rmd`
- `d-prime/dprime_calculation.py`

These scripts generate the derived psychometric and signal-detection measures used later in the pipeline.

### Stage 2: Run the main inferential models

Run:

- `glmm_analysis/GLMM_CP.R`
- `glmm_analysis/GLMM_CP_bootstrap.R`
- `trail_updating_analysis/last_trail_updating_analysis.R`

### Stage 3: Run descriptive and follow-up plots

Run as needed:

- scripts in `hisotgrams/`
- scripts in `range_analysis/`
- scripts in `rtm_corr/`

These scripts are mostly parallel and can be run according to the figure or metric of interest.

---

## Data Inputs

The scripts expect CSV inputs in `data/`. The main files used across the repository are:

- `data/full_data_cp.csv` - main trial-level accuracy dataset
- `data/updating_data_cp.csv` - updating-analysis dataset
- `data/dprime_results.csv` - d-prime summaries generated by the Python workflow
- `data/dprime_results_with_range.csv` - d-prime values retaining the Range dimension
- `data/weibull_metric_results.csv` - subject-level Weibull outputs

Additional intermediate files are generated by the scripts themselves.

### Common variables

Most scripts use combinations of:

- `Subject`
- `Age`
- `Group`
- `ExperimentName`
- `Regression`
- `Range`
- `ACC`

---

## Reproducibility

### R

```r
renv::restore()
```

### Python

```bash
pip install -r requirements.txt
```

The repository uses local path helpers in `R/paths.R`, so scripts can be run either from the project root or from their own subdirectories.

---
