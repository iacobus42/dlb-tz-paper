# R Code for Tz/Dz/Az Analysis in DLB

The `code` folder contains all the code used in the analysis as Quarto file (`.qmd`) and RMarkdown (`.Rmd`) files. Rendered HTML files of the output are included as well (`.html`).

Files are:

- `find_treated` - find people who ever had a dispensing claim for Tz/Dz/Az, tamsulosin, or one of the 5ARIs
- `find_treated_demographics_enrollments` - find the enrollment periods, demographics for these users
- `find_LBD_dates` - find enrollees diagnosed with LDB and the first diagnosis date
- `find_luts_procs` - find LUTS-specific diagnoses and procedures
- `reduce_treated_enrollments` - apply exclusion rules
- `build_model_data_treated` - extract variables for candidate users to build propensity score model
- `fit_psm` - fit the propensity score models
- `match` generate the propensity score matches, uses `match.cpp` for C++ speed in the matching process
- `analysis/analysis_functions.R` - contains functions for use in the analysis shared between the three comparisons
- `analysis/tz-tam`, `analysis/tz-5ari`, `analysis/tam-5ari` - contains the primary analyses reported in the paper
- `analysis/make_figures.R` - makes the figures, numbers for tables in the paper
- `sensitivity-analyses/any_rx_number/analysis` - repeats the primary analysis in an intent-to-treat format
- `sensitivity-analyses/bph-restricted/fit_psm`, `sensitivity-analyses/bph-restricted/match`, `sensitivity-analyses/bph-restricted/analysis` - repeat the analysis but require a diagnosis of BPH to be included
- `sensitivity-analyses/time-varying-effects/analysis` - generalizes the primary analysis to allow time-interacted effects
- `sensitivity-analyses/washout/fit_psm`, `sensitivity-analyses/washout/match`, `sensitivity-analyses/washout/analysis` - repeat the primary analysis but delay the start time of follow-up
