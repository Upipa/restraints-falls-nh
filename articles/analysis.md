# Does Physical Restraint Use Prevent Falls in Nursing Homes?

## Background

Physical restraint use in nursing homes is a common but controversial
practice. One of the most frequently cited justifications is fall
prevention. This analysis examines whether facility-level restraint
prevalence is associated with lower fall rates, using longitudinal data
from 2016 to 2025 across nursing homes in the province of Trento, Italy.

## Data

Data come from the **Indicare Salute Lab** quality monitoring system.
Each observation represents a facility × care-sector × month
combination, pairing the restraint indicator (2.1, point prevalence on
the index day) with three nested fall indicators: total falls (1.1),
falls with any outcome (1.3), and falls with major outcome such as
fractures (1.5). The nesting ($`1.5 \subseteq 1.3 \subseteq 1.1`$)
allows the outcome cells to be modeled as a single multinomial vector
over residents.

## Model

The final model is a Bayesian hierarchical model with a collapsed
likelihood. Conditional on the observed restraint count $`X`$ in a
facility–sector pair $`k`$:

``` math
X \sim \text{Binomial}(n_x, p_k), \qquad
Y \mid X \sim \text{Binomial}\!\left(n_y,\; \bar p_k\right), \quad
\bar p_k = \frac{X}{n_x}\, q_k + \left(1 - \frac{X}{n_x}\right) r_k
```

where $`q_k`$ and $`r_k`$ are the fall probabilities for restrained and
non-restrained residents. Fall outcomes are factorized exactly as
Binomial(fall) × Multinomial(outcome \| fall), with outcome cells (no
outcome / minor / major). Restraint prevalence follows a per-pair
quadratic time trend on the logit scale (orthogonalized and
standardized),
$`\text{logit}(p_k) = \beta_k + \alpha_k t + \gamma_k t^2`$; fall rates
are partially pooled across pairs via Beta hyperpriors, and outcome
shares via a non-centered logistic-normal hierarchy.

``` r

# The posterior draws of the final model ship with the package (light
# version: all parameters and contrasts, without the bulky per-observation
# predictive replications):
draws <- load_fit()
```

``` r

# Key contrasts, computed live from the shipped draws:
posterior::summarise_draws(
  posterior::subset_draws(
    draws,
    variable = c("delta_fall", "rr_fall", "delta_major", "rr_major")
  )
)
```

    # A tibble: 4 × 10
      variable        mean   median      sd     mad       q5      q95  rhat ess_bulk
      <chr>          <dbl>    <dbl>   <dbl>   <dbl>    <dbl>    <dbl> <dbl>    <dbl>
    1 delta_fall  -1.96e-2 -1.99e-2 0.00862 0.00853 -0.0333  -0.00493  1.01     444.
    2 rr_fall      7.70e-1  7.62e-1 0.0938  0.0899   0.629    0.937    1.01     446.
    3 delta_major -7.68e-4 -7.94e-4 0.00152 0.00154 -0.00317  0.00175  1.01     430.
    4 rr_major     8.72e-1  8.06e-1 0.391   0.346    0.381    1.59     1.01     430.
    # ℹ 1 more variable: ess_tail <dbl>

``` r

# The full CmdStanMCMC object is too large to ship (202 MB, mostly the
# predictive replications used for PPC). To reproduce the analysis
# end-to-end from the packaged dataset:
prep <- prepare_model_data()
fit <- fit_model(prepared = prep) # full fit (~6 h); use quick = TRUE for a fast check

# The complete development workflow (all 12 model iterations, Bayes factors,
# posterior predictive checks) is in data-raw/model-development.qmd of the
# source repository.
```

## Results

The association between restraint and falls is **strongly
heterogeneous** across facility–sector pairs. Of the 47 pairs, about 15
show credible benefit (restrained residents fall less), about 8 show
credible harm, and the rest are uncertain — the effect changes sign
across facilities. For a new, unobserved facility–sector pair, the
posterior predictive probability of any benefit is only about 61%.

For **major-outcome falls** (the clinically important ones), the average
effect is compatible with zero and predictive evidence of benefit is
weak (P ≈ 0.60). On falls with minor outcomes, more pairs show credible
harm (10) than credible benefit (6).

## Conclusions

The effect of physical restraint on falls is real but not generalizable:
its sign and magnitude vary unpredictably across facilities. There is no
meaningful evidence that restraint reduces the falls that matter most —
those with major outcomes. Effectiveness as a fall-prevention strategy
should be evaluated case by case, based on the characteristics of each
facility, and confirmed with individual-level data.
