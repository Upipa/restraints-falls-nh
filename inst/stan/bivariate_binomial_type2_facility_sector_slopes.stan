// Bivariate Binomial Type II — facility×sector intercepts AND slopes (M9)
// COLLAPSED LIKELIHOOD (no convolutions)
//
// logit(p_k) = beta_k + alpha_k * t_i
//   beta_k:  intercept for facility×sector pair k (k = 1..K)
//   alpha_k: time slope for facility×sector pair k
//   q, r:    global fall probabilities (restrained / not restrained)
//
// The grouping unit is the (ente, id_settore) pair — sector 25 of facility A
// is a different physical unit than sector 25 of facility B.
//
// Likelihood approximation (as M7): falls conditional on the observed
// restraint count follow a single binomial over the full fall denominator:
//   X ~ Bin(n_x, p)                                   (exact)
//   Y | X ~ Bin(n_y, p_bar),  p_bar = (X/n_x)*q + (1 - X/n_x)*r
//
// Prior: beta_k, alpha_k ~ Normal(0, 10) (fully normalized);
//        q, r Uniform(0,1) implicit.
// NOTE: explicit lpmf/lpdf everywhere to retain normalizing constants for BF.

functions {
  real collapsed_bi2_lpmf(array[] int xy, int n_x, int n_y,
                          real p, real q, real r) {
    int x = xy[1];
    int y = xy[2];
    real lp = binomial_lpmf(x | n_x, p);
    real frac = (n_x > 0) ? x * 1.0 / n_x : 0.0;
    real p_bar = frac * q + (1 - frac) * r;
    return lp + binomial_lpmf(y | n_y, p_bar);
  }
}

data {
  // Both observed
  int<lower=0> N_both;
  array[N_both] int<lower=0> n_x_both;
  array[N_both] int<lower=0> n_y_both;
  array[N_both] int<lower=0> x_both;
  array[N_both] int<lower=0> y_both;
  array[N_both] int<lower=1> fac_sec_both;   // facility×sector index
  array[N_both] real t_both;

  // Only restraint observed
  int<lower=0> N_x_only;
  array[N_x_only] int<lower=0> n_x_only;
  array[N_x_only] int<lower=0> x_only;
  array[N_x_only] int<lower=1> fac_sec_x_only;
  array[N_x_only] real t_x_only;

  // Only falls observed
  int<lower=0> N_y_only;
  array[N_y_only] int<lower=0> n_y_only;
  array[N_y_only] int<lower=0> y_only;
  array[N_y_only] int<lower=1> fac_sec_y_only;
  array[N_y_only] real t_y_only;

  int<lower=1> K;  // number of facility×sector pairs
}

parameters {
  vector[K] beta;              // facility×sector intercepts on logit scale
  vector[K] alpha;             // facility×sector time slopes
  real<lower=0, upper=1> q;    // P(falls | restrained)
  real<lower=0, upper=1> r;    // P(falls | not restrained)
}

model {
  // Priors — fully normalized for Bayes factor computation
  target += normal_lpdf(beta  | 0, 10);
  target += normal_lpdf(alpha | 0, 10);

  // Likelihood — both observed
  for (i in 1:N_both) {
    int k = fac_sec_both[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_both[i]);
    array[2] int xy = {x_both[i], y_both[i]};
    target += collapsed_bi2_lpmf(xy | n_x_both[i], n_y_both[i], p_i, q, r);
  }

  // Likelihood — only restraint
  for (i in 1:N_x_only) {
    int k = fac_sec_x_only[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_x_only[i]);
    target += binomial_lpmf(x_only[i] | n_x_only[i], p_i);
  }

  // Likelihood — only falls
  for (i in 1:N_y_only) {
    int k = fac_sec_y_only[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_y_only[i]);
    real prob_fall_i = p_i * q + (1 - p_i) * r;
    target += binomial_lpmf(y_only[i] | n_y_only[i], prob_fall_i);
  }
}

generated quantities {
  real delta_qr = q - r;
  real rr = q / r;

  // Posterior predictive simulation — both observed
  array[N_both] int x_rep;
  array[N_both] int y_rep;
  for (i in 1:N_both) {
    int k = fac_sec_both[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_both[i]);
    int xr = binomial_rng(n_x_both[i], p_i);
    real frac = (n_x_both[i] > 0) ? xr * 1.0 / n_x_both[i] : 0.0;
    real p_bar = frac * q + (1 - frac) * r;
    x_rep[i] = xr;
    y_rep[i] = binomial_rng(n_y_both[i], p_bar);
  }

  // Posterior predictive simulation — only restraint
  array[N_x_only] int x_only_rep;
  for (i in 1:N_x_only) {
    int k = fac_sec_x_only[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_x_only[i]);
    x_only_rep[i] = binomial_rng(n_x_only[i], p_i);
  }

  // Posterior predictive simulation — only falls
  array[N_y_only] int y_only_rep;
  for (i in 1:N_y_only) {
    int k = fac_sec_y_only[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_y_only[i]);
    real prob_fall_i = p_i * q + (1 - p_i) * r;
    y_only_rep[i] = binomial_rng(n_y_only[i], prob_fall_i);
  }
}
