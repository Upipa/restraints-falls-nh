// Bivariate Binomial Type II — facility×sector quadratic trends (M11)
// COLLAPSED LIKELIHOOD + HIERARCHICAL q, r per facility×sector pair
//
// logit(p_k) = beta_k + alpha_k * t_i + gamma_k * t2_i
//   beta_k:  intercept for facility×sector pair k (k = 1..K)
//   alpha_k: linear time slope for pair k
//   gamma_k: quadratic time coefficient for pair k
//   q_k ~ Beta(a_q, b_q), r_k ~ Beta(a_r, b_r)  (hierarchical, per pair)
//   a_q, b_q, a_r, b_r ~ Exponential(1) (hyperpriors, independent)
//
// t2_i is the quadratic term ORTHOGONALIZED against 1 and t (Gram–Schmidt via
// regression of t^2 on t) and STANDARDIZED to unit pooled sd — computed in R
// on the z-score scale. Same construction as M8.
//
// Extends M10 (the current reference) by letting each facility×sector pair's
// restraint trend curve over time, aiming to capture the residual X-side
// variability that the linear trend misses.
//
// Likelihood (collapsed, as M7–M10):
//   X ~ Bin(n_x, p)                                   (exact)
//   Y | X ~ Bin(n_y, p_bar),  p_bar = (X/n_x)*q_k + (1 - X/n_x)*r_k
//
// Prior: beta_k, alpha_k ~ Normal(0, 10); gamma_k ~ Normal(0, 5) — still
//        weakly informative, but on the standardized t2 scale (sd 1) it
//        excludes the absurd curvatures allowed by sd = 10 and improves the
//        NUTS geometry. Hyperpriors Exponential(1) explicit; Beta hierarchy
//        explicit.
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
  array[N_both] real t2_both;                // orthogonalized quadratic term

  // Only restraint observed
  int<lower=0> N_x_only;
  array[N_x_only] int<lower=0> n_x_only;
  array[N_x_only] int<lower=0> x_only;
  array[N_x_only] int<lower=1> fac_sec_x_only;
  array[N_x_only] real t_x_only;
  array[N_x_only] real t2_x_only;

  // Only falls observed
  int<lower=0> N_y_only;
  array[N_y_only] int<lower=0> n_y_only;
  array[N_y_only] int<lower=0> y_only;
  array[N_y_only] int<lower=1> fac_sec_y_only;
  array[N_y_only] real t_y_only;
  array[N_y_only] real t2_y_only;

  int<lower=1> K;  // number of facility×sector pairs
}

parameters {
  vector[K] beta;                 // facility×sector intercepts on logit scale
  vector[K] alpha;                // facility×sector linear time slopes
  vector[K] gamma;                // facility×sector quadratic coefficients
  vector<lower=0, upper=1>[K] q;  // P(falls | restrained), per pair
  vector<lower=0, upper=1>[K] r;  // P(falls | not restrained), per pair
  real<lower=0> a_q;              // Beta hyperparameters for q
  real<lower=0> b_q;
  real<lower=0> a_r;              // Beta hyperparameters for r
  real<lower=0> b_r;
}

model {
  // Priors — fully normalized for Bayes factor computation
  target += normal_lpdf(beta  | 0, 10);
  target += normal_lpdf(alpha | 0, 10);
  target += normal_lpdf(gamma | 0, 5);

  // Hyperpriors and hierarchical priors
  target += exponential_lpdf(a_q | 1);
  target += exponential_lpdf(b_q | 1);
  target += exponential_lpdf(a_r | 1);
  target += exponential_lpdf(b_r | 1);
  target += beta_lpdf(q | a_q, b_q);
  target += beta_lpdf(r | a_r, b_r);

  // Likelihood — both observed
  for (i in 1:N_both) {
    int k = fac_sec_both[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_both[i] + gamma[k] * t2_both[i]);
    array[2] int xy = {x_both[i], y_both[i]};
    target += collapsed_bi2_lpmf(xy | n_x_both[i], n_y_both[i], p_i, q[k], r[k]);
  }

  // Likelihood — only restraint
  for (i in 1:N_x_only) {
    int k = fac_sec_x_only[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_x_only[i] + gamma[k] * t2_x_only[i]);
    target += binomial_lpmf(x_only[i] | n_x_only[i], p_i);
  }

  // Likelihood — only falls
  for (i in 1:N_y_only) {
    int k = fac_sec_y_only[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_y_only[i] + gamma[k] * t2_y_only[i]);
    real prob_fall_i = p_i * q[k] + (1 - p_i) * r[k];
    target += binomial_lpmf(y_only[i] | n_y_only[i], prob_fall_i);
  }
}

generated quantities {
  // Population-level summaries (comparable to the global q, r of M7–M9)
  real q_bar = mean(q);
  real r_bar = mean(r);
  real delta_qr = q_bar - r_bar;
  real rr = q_bar / r_bar;

  // Posterior predictive simulation — both observed
  array[N_both] int x_rep;
  array[N_both] int y_rep;
  for (i in 1:N_both) {
    int k = fac_sec_both[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_both[i] + gamma[k] * t2_both[i]);
    int xr = binomial_rng(n_x_both[i], p_i);
    real frac = (n_x_both[i] > 0) ? xr * 1.0 / n_x_both[i] : 0.0;
    real p_bar = frac * q[k] + (1 - frac) * r[k];
    x_rep[i] = xr;
    y_rep[i] = binomial_rng(n_y_both[i], p_bar);
  }

  // Posterior predictive simulation — only restraint
  array[N_x_only] int x_only_rep;
  for (i in 1:N_x_only) {
    int k = fac_sec_x_only[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_x_only[i] + gamma[k] * t2_x_only[i]);
    x_only_rep[i] = binomial_rng(n_x_only[i], p_i);
  }

  // Posterior predictive simulation — only falls
  array[N_y_only] int y_only_rep;
  for (i in 1:N_y_only) {
    int k = fac_sec_y_only[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_y_only[i] + gamma[k] * t2_y_only[i]);
    real prob_fall_i = p_i * q[k] + (1 - p_i) * r[k];
    y_only_rep[i] = binomial_rng(n_y_only[i], prob_fall_i);
  }
}
