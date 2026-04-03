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
- SPSS repeated-measures ANOVA follow-up
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
- `SPSS/` - SPSS repeated-measures ANOVA syntax and short documentation
- `rtm_corr/` - correlation analyses linking RTM-related measures
- `hisotgrams/` - grouped summary plots for d-prime and PSE
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
  - Main outputs:
    - `data/weibull_trial_predictions.csv`
    - `data/weibull_metric_results.csv`
    - `output/plots/Weibull_groups_by_race.png`
    - `output/plots/Weibull_groups_by_bias.png`

- `d-prime/dprime_calculation.py`
  - Computes d-prime, criterion, and related signal-detection outputs.
  - Writes derived CSV files used by downstream scripts.
  - Main outputs:
    - `data/dprime_results_with_range.csv`
    - `data/dprime_results.csv`
    - `output/plots/Dprime_summary.png`
    - `output/plots/CR_summary.png`
    - `output/plots/Bprime_summary.png`

### Main model-based analyses

- `glmm_analysis/GLMM_CP.R`
  - Main trial-level GLMM of accuracy.
  - Tests effects of Group, Race, Regression condition, and age.
  - Main visualization output:
    - `output/plots/glmm_cp_.png`

- `glmm_analysis/GLMM_CP_bootstrap.R`
  - Bootstrap follow-up to the main GLMM.
  - Evaluates stability of peak-age estimates from the modeled trajectories.
  - Main outputs:
    - bootstrap summary CSV/RDS files written to `output/`

- `trail_updating_analysis/last_trail_updating_analysis.R`
  - Compares updating models based on the long-term average (`t-inf`) and previous trial (`t-1`).
  - Main outputs:
    - `output/updating_beta_results.csv`
    - `output/last_trail_updating_analysis.png`

### Follow-up and descriptive analyses

- `SPSS/IBM SPSS Statistics - Code.sps`
  - SPSS repeated-measures ANOVA syntax.
  - Use after generating the d-prime and Weibull/PSE-derived metrics.
  - See `SPSS/SPSS_repeated_measures_ANOVA.md` for a short explanation of the model and factors.

- `hisotgrams/dprime_histograms.R`
  - Subject-level d-prime grouped summaries by Group, Race, and Regression condition.
  - Main visualization output:
    - `output/plots/dprime_histogram.png`

- `hisotgrams/dprime_young_sub_hist.R`
  - Younger-subgroup version of the d-prime grouped summary.
  - Main visualization output:
    - `output/plots/dprime_histogram_young.png`

- `hisotgrams/pse_histograms.R`
  - Grouped PSE summaries based on the Weibull outputs.
  - Main visualization output:
    - `output/plots/pse_histogram.png`

- `hisotgrams/dprime_cfmt_group.R`
  - d-prime grouped summary split by CFMT subgroup.
  - Main visualization output:
    - `output/plots/dprime_cfmt_group.png`

- `range_analysis/accXrange_CP_bias.R`
  - Accuracy as a function of morph range, grouped by Race within Bias condition.
  - Main visualization output:
    - `output/plots/acc_range_CP_bias.png`

- `range_analysis/accXrange_CP_race.R`
  - Accuracy as a function of morph range, grouped by Bias condition within Race.
  - Main visualization output:
    - `output/plots/acc_range_CP_race.png`

- `rtm_corr/rtm_correlation_cp.R`
  - Correlates RTM magnitude with Bias- accuracy.
  - Main visualization output:
    - `output/plots/Correlation_biasM_vs_magnitude_CP.png`

- `rtm_corr/dprime_rtm_corr.R`
  - Correlates RTM magnitude with d-prime.
  - Main visualization output:
    - `output/plots/Correlation_dprime_biasM_vs_magnitude_CP.png`

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

Dependency note:

- `glmm_analysis/GLMM_CP_bootstrap.R` is a follow-up to the main GLMM and should be run after `glmm_analysis/GLMM_CP.R`

### Stage 3: Run descriptive and follow-up plots

Run as needed:

- files in `SPSS/` repeated-measures ANOVA
- scripts in `hisotgrams/`
- scripts in `range_analysis/`
- scripts in `rtm_corr/`

Dependency note:

- The SPSS analysis should be run only after Stage 1, because it uses the derived d-prime and Weibull/PSE metrics.
- `SPSS/IBM SPSS Statistics - Code.sps` contains the syntax, and `SPSS/SPSS_repeated_measures_ANOVA.md` gives the short explanation.

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

### Subject exclusions

- In the Weibull analysis, one CP subject: (829) was excluded because below-chance performance did not allow a stable Weibull fit within a reasonable range; accordingly, two matched control (1704, 1700) subjects were also excluded.


### Main parameters

| Parameter | Typical values | Role in the analysis |
|-----------|----------------|----------------------|
| `Subject` | participant ID | Repeated-measures unit; random intercept / grouping unit across analyses |
| `Age` | continuous years | Main age-related predictor in the GLMM and age-based follow-up analyses |
| `Group` | `TD`, `CP` | Between-subject diagnostic grouping |
| `ExperimentName` | `Caucasian`, `Asian` | Face-race condition; interpreted as Own-Race vs Other-Race |
| `Regression` | `biasp`, `biasm`, sometimes `Null` | RTM-related condition; `biasp` and `biasm` are the main discrimination conditions, while `Null` is used in the d-prime workflow for false-alarm estimation |
| `Range` | morph-distance levels | Stimulus-difficulty / perceptual-distance variable |
| `ACC` | `0`, `1` | Trial-level accuracy outcome |
| `CFMTgroup` | typically `L`, `H` | TD ability subgrouping used in selected follow-up summaries |
| `New_tinf` | numeric | Updating-model term reflecting accumulated / long-term stimulus history |
| `New_t1` | numeric | Updating-model term reflecting recent-history contribution from the previous trial |

### Main visualization outputs

The repository generates publication-oriented figures in `output/plots/` and `output/`. Main examples include:

- GLMM age-trajectory figure: `output/plots/glmm_cp_.png`
- Weibull summary figures: `output/plots/Weibull_groups_by_race.png`, `output/plots/Weibull_groups_by_bias.png`
- Histogram summaries: `output/plots/dprime_histogram.png`, `output/plots/dprime_histogram_young.png`, `output/plots/pse_histogram.png`, `output/plots/dprime_cfmt_group.png`
- Range-based descriptive figures: `output/plots/acc_range_CP_bias.png`, `output/plots/acc_range_CP_race.png`
- Correlation figures: `output/plots/Correlation_biasM_vs_magnitude_CP.png`, `output/plots/Correlation_dprime_biasM_vs_magnitude_CP.png`
- Updating-model figure: `output/last_trail_updating_analysis.png`

---


## Reproducibility

### R

```r
renv::snapshot()  # Capture current state
renv::restore()   # Restore to locked state
```

Repository R version stored in `renv.lock`:

- `R 4.5.2`

### Python

```bash
pip install -r requirements.txt
```

Local Python version in this environment:

- `Python 3.14.3`

Key Python packages are pinned in `requirements.txt`.

The repository uses local path helpers in `R/paths.R`, so scripts can be run either from the project root or from their own subdirectories.

---

## For More Context

The analytical methods, theoretical background, hypotheses, and full discussion of results are presented in the associated manuscript. This repository is intended to support **reproducibility and transparency** of the quantitative analyses.

**To understand the research questions and scientific context, please consult the authors.**

---

## License

This repository is licensed under the **MIT License**.

---

**Last Updated:** April 2026
