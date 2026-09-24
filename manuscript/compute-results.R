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

# Population-level contrast via the Beta hyperparameters
a_q <- getv("a_q")
b_q <- getv("b_q")
a_r <- getv("a_r")
b_r <- getv("b_r")
p_pop_neg <- mean(a_q / (a_q + b_q) - a_r / (a_r + b_r) < 0)

# Posterior predictive draws for a new, unobserved pair
set.seed(7301)
ndraw <- length(a_q)
q_new <- rbeta(ndraw, a_q, b_q)
r_new <- rbeta(ndraw, a_r, b_r)
p_new_benefit <- mean(q_new < r_new)
rr_new_med <- median(q_new / r_new)

# New-pair outcome shares from the logistic-normal hierarchy
softmax_rows <- function(m) {
  t(apply(m, 1, function(x) exp(x - max(x)) / sum(exp(x - max(x)))))
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
sd_logit_q <- mean(apply(getm("q_fall"), 1, function(v) sd(qlogis(v))))
sd_logit_r <- mean(apply(getm("r_fall"), 1, function(v) sd(qlogis(v))))

# Formatting helpers (95% HDI is the bayestestR default)
fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
fmt_hdi <- function(v, digits = 3) {
  h <- bayestestR::hdi(v)
  paste0("[", fmt(h$CI_low, digits), ", ", fmt(h$CI_high, digits), "]")
}
