# restraintsfalls

Research compendium investigating the association between physical restraint use and falls in nursing homes, using longitudinal quality monitoring data from the province of Trento (Italy), 2016–2025.

## Overview

Data come from **Indicare Salute Lab**, a quality monitoring system managed by UPIPA (Unione Provinciale Istituzioni Per l'Assistenza). The analysis covers facility-level aggregated indicators collected annually across nursing homes in the province.

The main research question is whether physical restraint use is associated with a reduction in falls — and whether any such association is clinically meaningful.

## Project structure

This project is organized as an R package for consistency and reproducibility:
  
```
R/               Functions (data retrieval, utilities)
data-raw/        Scripts that produce the analysis dataset
inst/stan/       Stan model files
inst/extdata/    Pre-computed results (posteriors, cleaned data)
tests/           Unit tests
vignettes/       Reproducible analysis (accessible after install)
manuscript/      Journal submission version (not in package)
```

## Setup

1. Clone the repository
2. Open in RStudio/Positron
3. Run `renv::restore()` to install dependencies
4. Set `UID` and `PWD` environment variables in `.Renviron` for database access

Alternatively, install directly from GitHub to access functions and the analysis vignette:

```r
renv::install("upipa/restraints-falls-nh")
vignette("analysis", "restraintsfalls")
```

## Methods

Bayesian binomial models estimated with Stan, linking the logit of restraint prevalence to the logit of fall rates, with a temporal trend component. Models are compared against an independence baseline.

## Status

Work in progress. Target journal: *BMC Geriatrics*.
