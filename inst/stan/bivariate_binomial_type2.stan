// Bivariate Binomial Type II — with partial observations
//
// BI2(n_x, n_y, p, q, r) decomposes as:
//   (X1, Y1) ~ BI1(min(n_x, n_y), p, q, r)    [joint part]
//   X2 ~ Bin(n_x - min(n_x, n_y), p)           [excess X, independent]
//   Y2 ~ Bin(n_y - min(n_x, n_y), p*q+(1-p)*r) [excess Y, independent]
//   X = X1 + X2,  Y = Y1 + Y2
//
// The PMF is computed via convolution (log-sum-exp over the split).
// Since |n_x - n_y| is typically 0-4, the convolution has very few terms.
//
// For partial observations (only X or only Y observed):
//   Marginal of X ~ Bin(n_x, p)
//   Marginal of Y ~ Bin(n_y, p*q + (1-p)*r)
//
// NOTE: All log-probability contributions use explicit lpmf (not ~)
// to retain normalizing constants, enabling Bayes factor computation.

functions {
  // Log-PMF of BI1 (joint distribution with common n)
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

  // Log-PMF of BI2 via convolution
  real bivariate_binomial_type2_lpmf(array[] int xy, int n_x, int n_y,
                                      real p, real q, real r) {
    int x = xy[1];
    int y = xy[2];
    int n_min = min(n_x, n_y);
    int excess_x = n_x - n_min;
    int excess_y = n_y - n_min;
    real prob_fall = p * q + (1 - p) * r;

    if (excess_x == 0 && excess_y == 0) {
      // n_x == n_y: pure BI1
      return bivariate_binomial_type1_lpmf(xy | n_min, p, q, r);
    }

    if (excess_y > 0) {
      // n_x <= n_y: X = X1, Y = Y1 + Y2
      // Convolve over y1: P(X=x, Y=y) = sum_{y1} P_BI1(x, y1) * P_Bin(y-y1)
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
      // n_x > n_y: X = X1 + X2, Y = Y1
      // Convolve over x1: P(X=x, Y=y) = sum_{x1} P_BI1(x1, y) * P_Bin(x-x1)
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

  // Only restraint observed
  int<lower=0> N_x_only;
  array[N_x_only] int<lower=0> n_x_only;
  array[N_x_only] int<lower=0> x_only;

  // Only falls observed
  int<lower=0> N_y_only;
  array[N_y_only] int<lower=0> n_y_only;
  array[N_y_only] int<lower=0> y_only;
}

parameters {
  real<lower=0, upper=1> p;   // P(restrained)
  real<lower=0, upper=1> q;   // P(falls | restrained)
  real<lower=0, upper=1> r;   // P(falls | not restrained)
}

model {
  // Likelihood — both observed (joint BI2)
  for (i in 1:N_both) {
    array[2] int xy = {x_both[i], y_both[i]};
    target += bivariate_binomial_type2_lpmf(xy | n_x_both[i], n_y_both[i], p, q, r);
  }

  // Likelihood — only restraint (marginal: X ~ Bin(n_x, p))
  for (i in 1:N_x_only) {
    target += binomial_lpmf(x_only[i] | n_x_only[i], p);
  }

  // Likelihood — only falls (marginal: Y ~ Bin(n_y, p*q + (1-p)*r))
  for (i in 1:N_y_only) {
    target += binomial_lpmf(y_only[i] | n_y_only[i], p * q + (1 - p) * r);
  }
}

generated quantities {
  real prob_fall = p * q + (1 - p) * r;
  real delta_qr = q - r;
  real rr = q / r;

  // Posterior predictive simulation — both observed
  array[N_both] int x_rep;
  array[N_both] int y_rep;
  for (i in 1:N_both) {
    int n_min = min(n_x_both[i], n_y_both[i]);
    int excess_x = n_x_both[i] - n_min;
    int excess_y = n_y_both[i] - n_min;

    // Simulate from BI1(n_min, p, q, r) via decomposition
    int x1 = binomial_rng(n_min, p);
    int y1_from_restrained = binomial_rng(x1, q);
    int y1_from_free = binomial_rng(n_min - x1, r);
    int y1 = y1_from_restrained + y1_from_free;

    // Simulate excess
    int x2 = (excess_x > 0) ? binomial_rng(excess_x, p) : 0;
    int y2 = (excess_y > 0) ? binomial_rng(excess_y, prob_fall) : 0;

    x_rep[i] = x1 + x2;
    y_rep[i] = y1 + y2;
  }

  // Posterior predictive simulation — only restraint
  array[N_x_only] int x_only_rep;
  for (i in 1:N_x_only) {
    x_only_rep[i] = binomial_rng(n_x_only[i], p);
  }

  // Posterior predictive simulation — only falls
  array[N_y_only] int y_only_rep;
  for (i in 1:N_y_only) {
    y_only_rep[i] = binomial_rng(n_y_only[i], prob_fall);
  }
}
