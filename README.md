# Nigeria Diphtheria Outbreak: Kano Line-List Sub-Analysis

R scripts for the Kano case-level sub-analysis supplementing "Epidemiology, Mortality, and
Intervention Strategies for Diphtheria Control: A Comprehensive Analysis of Nigeria's
2022-2026 Outbreak" (Musa BM et al.).

Rebuilds the outbreak analysis from a 22,379-row Kano line list (2022-05 to 2024-05),
replacing several sitrep-based estimates in the manuscript with case-level analysis:
a real (not back-calculated) epidemic curve, direct + Bayesian hierarchical age-specific
case fatality rates, a stable EpiEstim Rt estimate, and a Bayesian SEIR transmission model
(custom Metropolis-Hastings MCMC in base R -- no Stan/rstanarm dependency required to
rerun these scripts).

## Running

Scripts run in order from the project root (working directory containing `data/`,
`figures/`, `scripts/`):

```
Rscript scripts/00_setup.R
Rscript scripts/01_clean_linelist.R
Rscript scripts/02_epidemic_curve.R
Rscript scripts/03_age_cfr.R
Rscript scripts/04_rt_estimation.R
Rscript scripts/05_bayesian_transmission_model.R
Rscript scripts/06_vaccination_scenarios.R
Rscript scripts/07_policy_impact.R
Rscript scripts/08_geographic_distribution.R
```

Each script reads its inputs from `data/` and writes intermediate CSVs back to `data/`
and figures to `figures/`; run them in order the first time.

Input data (the raw Kano line list, and small reference CSVs transcribed from the
manuscript's own tables) are not included in this repository.

## Methodology notes

- `05_bayesian_transmission_model.R`'s header comment explains the SEIR model structure,
  why it's fit to a 27-week window rather than the full 91-week series, and how the
  effective population parameter (N_eff) is a fitted scale, not Kano's real population.
- Bayesian components use custom Metropolis-Hastings sampling rather than Stan/rstanarm,
  for portability across environments without a working C++/Stan toolchain.
