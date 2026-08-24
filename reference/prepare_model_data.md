# Prepare the data list for the multi-outcome model

Builds the complete Stan data list for the final multi-outcome model
(multinomial fall outcomes + per-pair quadratic restraint trend) from
the analysis dataset. The function reproduces the preparation pipeline
used for the published analysis: pivoting the four indicators, building
the multinomial outcome cells by successive differences (the fall
indicators are nested, 1.5 within 1.3 within 1.1), splitting records
into observation groups, and computing the orthogonalized, standardized
quadratic time term.

## Usage

``` r
prepare_model_data(data = restraints_falls)
```

## Arguments

- data:

  A data frame with the structure of
  [restraints_falls](https://upipa.github.io/restraints-falls-nh/reference/restraints_falls.md):
  columns `anno`, `mese`, `ente`, `id_settore`, `indicatore`, `n`, `d`.
  Defaults to the packaged dataset.

## Value

A list with two elements:

- `stan_data`: the complete data list for the Stan model
  `restraint_multinomial_fac_sec_quad.stan`.

- `meta`: metadata (pair labels, number of pairs `K`, time
  standardization constants, group sizes).

## Details

Records with nesting violations (`n_1.3 > n_1.1` or `n_1.5 > n_1.3` — 59
records in the packaged dataset, probable data-entry errors) fall back
to the binomial groups, which use only the total fall count. Row
ordering may differ from `data-raw/model-development.qmd`; the model is
identical.

## Examples

``` r
prep <- prepare_model_data()
str(prep$meta)
#> List of 6
#>  $ pair_labels: chr [1:47] "10_191" "10_192" "10_193" "10_194" ...
#>  $ K          : int 47
#>  $ t_mean     : num 2021
#>  $ t_sd       : num 2.94
#>  $ t2_sd      : num 0.834
#>  $ group_sizes: Named int [1:5] 2766 2174 350 349 70
#>   ..- attr(*, "names")= chr [1:5] "N_x" "N_multi" "N_bin" "N_mm" ...
```
