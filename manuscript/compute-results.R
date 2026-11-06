# Compute all manuscript numbers from the M12 posterior draws.
# Sourced by the first inline R expression in the YAML abstract (knitr
# evaluates YAML inline code before running any chunk) and, guarded, by the
# `setup` chunk. Working directory is manuscript/ when rendering.

library(restraintsfalls)

draws <- load_fit()
dm <- posterior::as_draws_matrix(draws)

getv <- function(v) as.numeric(posterior::subset_draws(dm, variable = v))
getm <- function(v) as.matrix(posterior::subset_draws(dm, variable = v))

# Key contrasts (generated quantities: means across the observed pairs)
delta_fall <- getv("delta_fall")
rr_fall <- getv("rr_fall")
delta_major <- getv("delta_major")
rr_major <- getv("rr_major")
delta_minor <- getm("q_bar")[, 2] - getm("r_bar")[, 2]


# Population-level contrasts via the Beta hyperparameters (means of the
# hyperpriors; the in-text estimands — realized per-pair values are in the
# supplement table)
a_q <- getv("a_q")
b_q <- getv("b_q")
a_r <- getv("a_r")
b_r <- getv("b_r")
pred_q <- a_q / (a_q + b_q)
pred_r <- a_r / (a_r + b_r)
pop_delta_fall <- pred_q - pred_r
pop_rr_fall <- pred_q / pred_r
pop_base_fall <- pred_r
p_pop_neg <- mean(pop_delta_fall < 0)

# Posterior predictive draws for a new, unobserved pair
set.seed(7301)
ndraw <- length(a_q)
q_new <- rbeta(ndraw, a_q, b_q)
r_new <- rbeta(ndraw, a_r, b_r)
p_new_benefit <- mean(q_new < r_new)
rr_new_med <- median(q_new / r_new)

# New-pair outcome shares from the logistic-normal hierarchy
# (rowwise max via vectorized pmax, no apply)
softmax_rows <- function(m) {
  mx <- do.call(pmax, as.data.frame(m))
  e <- exp(m - mx)
  e / rowSums(e)
}
s_q_new <- softmax_rows(
  getm("mu_q") +
    getm("sigma_q") * matrix(rnorm(3 * ndraw), ndraw, 3)
)
s_r_new <- softmax_rows(
  getm("mu_r") +
    getm("sigma_r") * matrix(rnorm(3 * ndraw), ndraw, 3)
)

# Per-resident predictive contrasts for a new pair, per outcome cell
delta_minor_new <- q_new * s_q_new[, 2] - r_new * s_r_new[, 2]
delta_major_new <- q_new * s_q_new[, 3] - r_new * s_r_new[, 3]
p_minor_new_red <- mean(delta_minor_new < 0)
p_major_new_red <- mean(delta_major_new < 0)

# Bayes factors for a population-level directional effect on the mean:
# H1 = {E[q_new | theta] < E[r_new | theta]} vs H0 = {>=}, i.e. the sign of the
# population mean difference. Under the prior, restrained and non-restrained
# rates are exchangeable, so the prior odds are exactly 1 and the Bayes factor
# equals the posterior odds of H1.
# Total falls: E[q_new | theta] is the Beta hyperprior mean -> p_pop_neg above.
p_mean_fall <- p_pop_neg
bf_mean_fall <- p_mean_fall / (1 - p_mean_fall)

softmax3 <- function(a) {
  mx <- pmax(a[,, 1], a[,, 2], a[,, 3])
  e <- exp(sweep(a, 1:2, mx))
  sweep(e, 1:2, e[,, 1] + e[,, 2] + e[,, 3], `/`)
}
# Outcome cells: E[rate | theta] = pred * E[s[m] | theta]; E[s[m]] under the
# logistic-normal has no closed form -> Monte Carlo over eta draws.
set.seed(9157)
Msim <- 2000
# Returns the population mean rates (restrained and non-restrained) per draw
cell_mean_rates <- function(m) {
  muq <- getm("mu_q")
  siq <- getm("sigma_q")
  mur <- getm("mu_r")
  sir <- getm("sigma_r")
  idx <- rep(1:3, each = Msim)
  eq <- array(
    rnorm(ndraw * Msim * 3, muq[, idx], siq[, idx]),
    c(ndraw, Msim, 3)
  )
  er <- array(
    rnorm(ndraw * Msim * 3, mur[, idx], sir[, idx]),
    c(ndraw, Msim, 3)
  )
  list(
    q = pred_q * rowMeans(softmax3(eq)[,, m]),
    r = pred_r * rowMeans(softmax3(er)[,, m])
  )
}
rates_minor <- cell_mean_rates(2)
rates_major <- cell_mean_rates(3)
ed_minor <- rates_minor$q - rates_minor$r
ed_major <- rates_major$q - rates_major$r

p_mean_minor <- mean(ed_minor < 0)
bf_mean_minor <- p_mean_minor / (1 - p_mean_minor)
p_mean_major <- mean(ed_major < 0)
bf_mean_major <- p_mean_major / (1 - p_mean_major)
rates_none <- cell_mean_rates(1)
ed_none <- rates_none$q - rates_none$r
p_mean_none <- mean(ed_none < 0)
bf_mean_none <- p_mean_none / (1 - p_mean_none)

# Population risk ratios per cell (report medians: the mean does not exist)
pop_rr_none <- rates_none$q / rates_none$r
pop_rr_minor <- rates_minor$q / rates_minor$r
pop_rr_major <- rates_major$q / rates_major$r
pop_base_major <- rates_major$r

# Realized (across observed pairs) delta for no-outcome falls
# (q_bar/r_bar columns are [total, minor, major])
delta_none <- delta_fall - delta_minor

# Share-of-outcome contrasts conditional on a fall having occurred (new pair)
p_share_minor_incr <- mean(s_q_new[, 2] > s_r_new[, 2])
p_share_major_incr <- mean(s_q_new[, 3] > s_r_new[, 3])

# Per-pair deltas: total falls and minor-outcome falls
# (in the flattened s_q/s_r draws, pair k varies fastest within outcome cell m:
#  cell m of pair k is column (m - 1) * K + k)
K <- ncol(getm("q_fall"))
delta_pairs <- getm("q_fall") - getm("r_fall")
delta_minor_pairs <- getm("q_fall") *
  getm("s_q")[, K + 1:K] -
  getm("r_fall") * getm("s_r")[, K + 1:K]

count_pairs <- function(dmat) {
  p_neg <- colMeans(dmat < 0)
  c(
    benefit = sum(p_neg >= 0.95),
    harm = sum(p_neg <= 0.05),
    uncertain = sum(p_neg > 0.05 & p_neg < 0.95)
  )
}
n_fall <- count_pairs(delta_pairs)
n_minor <- count_pairs(delta_minor_pairs)
max_harm_delta <- max(colMeans(delta_pairs)[colMeans(delta_pairs < 0) <= 0.05])

# Between-pair heterogeneity of fall rates (logit-scale SD, posterior mean)
# Rowwise SD vectorized: sqrt(mean squared deviation) with Bessel correction
row_sd <- function(m) {
  lq <- qlogis(m)
  sqrt(rowSums((lq - rowMeans(lq))^2) / (ncol(m) - 1))
}
sd_logit_q <- mean(row_sd(getm("q_fall")))
sd_logit_r <- mean(row_sd(getm("r_fall")))

# Formatting helpers (95% HDI is the bayestestR default)
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_hdi <- function(v, digits = 3) {
  h <- bayestestR::hdi(v)
  stringr::str_glue("[{fmt(h$CI_low, digits)}, {fmt(h$CI_high, digits)}]")
}
