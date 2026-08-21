# restraintsfalls

Research compendium investigating the association between physical
restraint use and falls in nursing homes, using longitudinal quality
monitoring data from the province of Trento (Italy), 2016–2025.

## Overview

Data come from **Indicare Salute Lab**, a quality monitoring system
managed by UPIPA (Unione Provinciale Istituzioni Per l’Assistenza).
Observations are facility × care-sector × month combinations, pairing
the restraint indicator (2.1, point prevalence) with three nested fall
indicators: total falls (1.1), falls with any outcome (1.3), and falls
with major outcome such as fractures (1.5).

The main research question is whether physical restraint use is
associated with a reduction in falls — and whether any such association
is clinically meaningful.

## Key findings

The restraint–fall association is **strongly heterogeneous** across
facility–sector pairs and **not generalizable**:

- Of 47 facility–sector pairs, ~15 show credible benefit, ~8 credible
  harm, and ~24 are uncertain — the effect changes sign across
  facilities.
- For a new, unobserved facility, the posterior predictive probability
  of any benefit is only ~61%.
- There is no meaningful evidence that restraint reduces **major-outcome
  falls** (the clinically important ones): the average effect is
  compatible with zero.
- On minor-outcome falls, more pairs show credible harm (10) than
  credible benefit (6).

Effectiveness of restraint as a fall-prevention strategy should be
evaluated case by case, based on the characteristics of each facility,
and confirmed with individual-level data.

## Project structure

This project is organized as an R package for consistency and
reproducibility:

    R/               Functions (data retrieval, utilities)
    data-raw/        Model development document + scripts that produce the dataset
    inst/stan/       Stan model files
    inst/extdata/    Pre-computed results (posteriors, cleaned data)
    tests/           Unit tests
    vignettes/       Analysis vignette (accessible after install)
    manuscript/      Journal submission version (not in package)

## Setup

1.  Clone the repository
2.  Open in RStudio/Positron
3.  Run
    [`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html)
    to install dependencies
4.  Set `UID` and `PWD` environment variables in `.Renviron` for
    database access

Alternatively, install directly from GitHub to access functions and the
analysis vignette:

``` r

renv::install("upipa/restraints-falls-nh")
vignette("analysis", "restraintsfalls")
```

## Methods

Bayesian hierarchical models estimated with Stan (cmdstanr). The final
model uses a collapsed likelihood — falls conditional on the observed
restraint count, $`Y \mid X \sim \text{Bin}(n_y, \bar p)`$ with
$`\bar p = (X/n_x)\,q + (1 - X/n_x)\,r`$ — with the nested fall outcomes
factorized exactly as Binomial(fall) × Multinomial(outcome \| fall).
Restraint prevalence follows a per-pair quadratic time trend
(orthogonalized, standardized); fall rates are partially pooled across
facility–sector pairs via Beta hyperpriors, and outcome shares via a
non-centered logistic-normal hierarchy. Model development proceeded
through 12 iterations with Bayes factors (bridge sampling) and posterior
predictive checks; the full workflow is documented in
`data-raw/model-development.qmd`.

## Status

Analysis complete. Manuscript in preparation. Target journal: *BMC
Geriatrics*.
