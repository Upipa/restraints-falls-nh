# Changelog

## restraintsfalls (development version)

- New
  [`fit_model()`](https://upipa.github.io/restraints-falls-nh/reference/fit_model.md)
  compiles and fits the final multi-outcome model (M12) with cmdstanr,
  reproducing the published fit by default (with a `quick` mode for fast
  checks).
- New
  [`load_fit()`](https://upipa.github.io/restraints-falls-nh/reference/load_fit.md)
  loads the posterior draws of the published fit shipped with the
  package (light version, without the per-observation predictive
  replications).
- New
  [`prepare_model_data()`](https://upipa.github.io/restraints-falls-nh/reference/prepare_model_data.md)
  builds the complete Stan data list for the final multi-outcome model
  from the analysis dataset.
