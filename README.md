# VineMETIC

Reproducibility material for the manuscript "Analysis of multivariate
event times under informative censoring using vine copula" by Xinyuan
Chen, Yiwei Li, and Qian M. Zhou (revision & resubmission, *Lifetime
Data Analysis*).

The proposed stage-wise vine-copula PMLE methodology itself is
implemented in the R package
[VineSurvIC](https://github.com/michellezhou2009/VineSurvIC), which
this repository depends on. This repo contains only the harness code
needed to reproduce the manuscript's three simulation studies: the
data-generating configurations, the simulation drivers, and the
scripts that gather raw output into the tables and figures reported in
the manuscript and its supplementary material.

## Repository structure

- `settings/` -- data-generating configuration for each simulation
  study (`sim1_setting.R`, `sim2_setting.R`, `sim3_setting.R`, plus
  their pre-built `.rds` equivalents).
  - Sim I: heterogeneous copula families with covariates (Gumbel
    C(1,3), Clayton C(2,3), Gaussian C(1,2|3)).
  - Sim II: same DGP as Sim I, fit under four copula-family
    misspecification scenarios.
  - Sim III: symmetric trivariate Clayton DGP without covariates, at
    two true tau levels, comparing VineSurvIC's vine-PMLE against the
    Li (2020) pseudo-MLE benchmark.
- `drivers/` -- simulation drivers that generate data, fit models, and
  save one `.rds` per replication under `raw_results/`.
  - `run_sim1.R`, `run_sim2.R`, `run_sim3.R`
  - `compute_tau12_sim3.R` -- post-processing step run after
    `run_sim3.R` to add the `tau12` (T1,T2-given-D) estimate to each
    Sim III replication file.
- `gather/` -- aggregation scripts that read `raw_results/` and write
  the tables/figures reported in the manuscript and supplement, as
  `.csv` files (and, for Sim II, `.pdf`/`.png` figures under
  `figures/`).
  - `gather_sim1.R` -> Table 1 (main text) and Web Table S1 (supp)
  - `gather_sim2.R` -> Table 2 (main text) and Figure 2 (main text),
    Web Figures S2-S4 (supp)
  - `gather_sim3.R` -> Web Table (Sim III, supp)
- `utils/` -- shared helper functions used by the settings/drivers
  scripts (censoring-rate calibration, tau12 quadrature, the Li (2020)
  pseudo-MLE fitting routine), trimmed down from the original
  simulation harness to only what the current code calls.

## Reproducing the results

All scripts assume the working directory is the repository root.

```r
# 1. Build the data-generating settings (only needed once; .rds files
#    are already included)
source("settings/sim1_setting.R")
source("settings/sim2_setting.R")
source("settings/sim3_setting.R")

# 2. Run the simulations (writes raw_results/simN/...)
source("drivers/run_sim1.R")
source("drivers/run_sim2.R")
source("drivers/run_sim3.R")
source("drivers/compute_tau12_sim3.R")  # after run_sim3.R

# 3. Gather results into tables/figures (writes gather/ and figures/)
source("gather/gather_sim1.R")
source("gather/gather_sim2.R")
source("gather/gather_sim3.R")
```

## Dependencies

- [VineSurvIC](https://github.com/michellezhou2009/VineSurvIC)
- VineCopula, survival, trust, dplyr, tidyr, ggplot2, purrr

## Acknowledgments

`utils/li2020_fit.R` implements the Li (2020) pseudo-MLE benchmark used
in Simulation Study III. Its `npest.star_12()` function is revised from
the function of the same name in the
[`bvic`](https://github.com/dli-stats/bvic) package, provided in:

Li, D., Hu, X. J., & Wang, R. (2023). Evaluating Association Between
Two Event Times with Observations Subject to Informative Censoring.
*Journal of the American Statistical Association*, 118(542), 1282-1294.
https://doi.org/10.1080/01621459.2021.1990766
