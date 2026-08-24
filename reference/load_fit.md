# Load the posterior draws of the final model

Loads the posterior draws of the final multi-outcome model (M12),
shipped with the package as a light `draws_df` object: all parameters
and contrasts (4,000 draws, 4 chains, seed 919214260), without the bulky
per-observation posterior predictive replications.

## Usage

``` r
load_fit()
```

## Value

A `draws_df` object with the posterior draws of the published fit.

## Details

The full `CmdStanMCMC` object is too large to ship (202 MB, mostly the
predictive replications used for PPC). To reproduce the analysis
end-to-end, use
[`prepare_model_data()`](https://upipa.github.io/restraints-falls-nh/reference/prepare_model_data.md)
and
[`fit_model()`](https://upipa.github.io/restraints-falls-nh/reference/fit_model.md).

## Examples

``` r
if (FALSE) { # \dontrun{
draws <- load_fit()
posterior::summarise_draws(
  posterior::subset_draws(draws, variable = c("delta_fall", "rr_fall"))
)
} # }
```
