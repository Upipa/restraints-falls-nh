// Bivariate Binomial Type II — facility intercepts + facility-specific time slopes
//
// logit(p_ij) = beta_j + alpha_j * t_i
//   beta_j:  facility-specific intercept (j = 1..J)
//   alpha_j: facility-specific time slope
//   q, r:    global fall probabilities (restrained / not restrained)
//
// Prior: beta_j, alpha_j ~ Normal(0, 10) independent (fully normalized);
//        q, r Uniform(0,1) implicit.
//
// NOTE: All log-probability contributions use explicit lpmf/lpdf
// to retain normalizing constants for Bayes factor computation.

functions {
  real bivariate_binomial_type1_lpmf(array[] int xy, int n, real p, real q, real r) {
    int x = xy[1];
    int y = xy[2];

    int m = max(0, x + y - n);
    int M = min(x, y);

    real log_base = x * log(p) + (n - x) * log1m(p)
                  + x * log1m(q) + y * log(r) + (n - x - y) * log1m(r);

    real log_ratio = log(q) + log1m(r) - log1m(q) - log(r);

    int n_terms = M - m + 1;
    vector[n_terms] log_terms;
    for (i in 1:n_terms) {
      int d = m + i - 1;
      log_terms[i] = lgamma(n + 1)
                   - lgamma(n - x - y + d + 1)
                   - lgamma(x - d + 1)
                   - lgamma(y - d + 1)
                   - lgamma(d + 1)
                   + d * log_ratio;
    }

    return log_base + log_sum_exp(log_terms);
  }

  real bivariate_binomial_type2_lpmf(array[] int xy, int n_x, int n_y,
                                      real p, real q, real r) {
    int x = xy[1];
    int y = xy[2];
    int n_min = min(n_x, n_y);
    int excess_x = n_x - n_min;
    int excess_y = n_y - n_min;
    real prob_fall = p * q + (1 - p) * r;

    if (excess_x == 0 && excess_y == 0) {
      return bivariate_binomial_type1_lpmf(xy | n_min, p, q, r);
    }

    if (excess_y > 0) {
      int y1_lo = max(0, y - excess_y);
      int y1_hi = min(y, n_min);
      int n_terms = y1_hi - y1_lo + 1;
      vector[n_terms] log_terms;
      for (i in 1:n_terms) {
        int y1 = y1_lo + i - 1;
        int y2 = y - y1;
        array[2] int xy1 = {x, y1};
        log_terms[i] = bivariate_binomial_type1_lpmf(xy1 | n_min, p, q, r)
                      + binomial_lpmf(y2 | excess_y, prob_fall);
      }
      return log_sum_exp(log_terms);
    } else {
      int x1_lo = max(0, x - excess_x);
      int x1_hi = min(x, n_min);
      int n_terms = x1_hi - x1_lo + 1;
      vector[n_terms] log_terms;
      for (i in 1:n_terms) {
        int x1 = x1_lo + i - 1;
        int x2 = x - x1;
        array[2] int xy1 = {x1, y};
        log_terms[i] = bivariate_binomial_type1_lpmf(xy1 | n_min, p, q, r)
                      + binomial_lpmf(x2 | excess_x, p);
      }
      return log_sum_exp(log_terms);
    }
  }
}

data {
  // Both observed
  int<lower=0> N_both;
  array[N_both] int<lower=0> n_x_both;
  array[N_both] int<lower=0> n_y_both;
  array[N_both] int<lower=0> x_both;
  array[N_both] int<lower=0> y_both;
  array[N_both] int<lower=1> ente_both;
  array[N_both] real t_both;

  // Only restraint observed
  int<lower=0> N_x_only;
  array[N_x_only] int<lower=0> n_x_only;
  array[N_x_only] int<lower=0> x_only;
  array[N_x_only] int<lower=1> ente_x_only;
  array[N_x_only] real t_x_only;

  // Only falls observed
  int<lower=0> N_y_only;
  array[N_y_only] int<lower=0> n_y_only;
  array[N_y_only] int<lower=0> y_only;
  array[N_y_only] int<lower=1> ente_y_only;
  array[N_y_only] real t_y_only;

  int<lower=1> J;  // number of facilities
}

parameters {
  array[J] real beta;            // facility-specific intercepts on logit scale
  array[J] real alpha;           // facility-specific time slopes
  real<lower=0, upper=1> q;     // P(falls | restrained)
  real<lower=0, upper=1> r;     // P(falls | not restrained)
}

model {
  // Priors — fully normalized for Bayes factor computation
  target += normal_lpdf(beta  | 0, 10);
  target += normal_lpdf(alpha | 0, 10);

  // Likelihood — both observed
  for (i in 1:N_both) {
    int j = ente_both[i];
    real p_i = inv_logit(beta[j] + alpha[j] * t_both[i]);
    array[2] int xy = {x_both[i], y_both[i]};
    target += bivariate_binomial_type2_lpmf(xy | n_x_both[i], n_y_both[i], p_i, q, r);
  }

  // Likelihood — only restraint
  for (i in 1:N_x_only) {
    int j = ente_x_only[i];
    real p_i = inv_logit(beta[j] + alpha[j] * t_x_only[i]);
    target += binomial_lpmf(x_only[i] | n_x_only[i], p_i);
  }

  // Likelihood — only falls
  for (i in 1:N_y_only) {
    int j = ente_y_only[i];
    real p_i = inv_logit(beta[j] + alpha[j] * t_y_only[i]);
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
    int j = ente_both[i];
    real p_i = inv_logit(beta[j] + alpha[j] * t_both[i]);
    real prob_fall_i = p_i * q + (1 - p_i) * r;
    int n_min = min(n_x_both[i], n_y_both[i]);
    int excess_x = n_x_both[i] - n_min;
    int excess_y = n_y_both[i] - n_min;

    int x1 = binomial_rng(n_min, p_i);
    int y1_from_restrained = binomial_rng(x1, q);
    int y1_from_free = binomial_rng(n_min - x1, r);
    int y1 = y1_from_restrained + y1_from_free;

    int x2 = (excess_x > 0) ? binomial_rng(excess_x, p_i) : 0;
    int y2 = (excess_y > 0) ? binomial_rng(excess_y, prob_fall_i) : 0;

    x_rep[i] = x1 + x2;
    y_rep[i] = y1 + y2;
  }

  // Posterior predictive simulation — only restraint
  array[N_x_only] int x_only_rep;
  for (i in 1:N_x_only) {
    int j = ente_x_only[i];
    real p_i = inv_logit(beta[j] + alpha[j] * t_x_only[i]);
    x_only_rep[i] = binomial_rng(n_x_only[i], p_i);
  }

  // Posterior predictive simulation — only falls
  array[N_y_only] int y_only_rep;
  for (i in 1:N_y_only) {
    int j = ente_y_only[i];
    real p_i = inv_logit(beta[j] + alpha[j] * t_y_only[i]);
    real prob_fall_i = p_i * q + (1 - p_i) * r;
    y_only_rep[i] = binomial_rng(n_y_only[i], prob_fall_i);
  }
}
