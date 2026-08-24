# restraintsfalls (development version)

* New `fit_model()` compiles and fits the final multi-outcome model (M12) with cmdstanr, reproducing the published fit by default (with a `quick` mode for fast checks).
* New `prepare_model_data()` builds the complete Stan data list for the final multi-outcome model from the analysis dataset.
