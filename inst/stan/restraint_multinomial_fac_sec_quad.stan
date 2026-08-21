// Restraint + fall outcomes — facility×sector quadratic trends (M12)
// FACTORIZED FORM: Binomial(fall) x Multinomial(outcome | fall)
//
// Joint model over restraint (2.1) and all fall outcomes (1.1, 1.3, 1.5).
// The 4-cell multinomial over residents is factorized exactly:
//   Level 1: Y ~ Bin(n_y, p_fall)         — fall vs not (M11 structure)
//   Level 2: (c1, c2, c3) | Y ~ Multinomial(Y, s_bar) — outcome composition
// with cells (fall without outcome, minor, major) = (n_1.1-n_1.3, n_1.3-n_1.5,
// n_1.5) — the indicators are nested (1.5 ⊆ 1.3 ⊆ 1.1).
//
// The factorization is EXACT (standard multinomial property); the linear
// transition in the restraint fraction frac = X/n_x applies at both levels:
//   p_fall = frac*q_fall[k] + (1-frac)*r_fall[k]
//   s_bar  = [frac*q_fall[k]*s_q[k] + (1-frac)*r_fall[k]*s_r[k]] / p_fall
// (note: the exact level-2 weights involve the fall probabilities too).
//
// Hierarchy:
//   q_fall[k] ~ Beta(a_q, b_q),  r_fall[k] ~ Beta(a_r, b_r)   (as M10/M11)
//   a_q, b_q, a_r, b_r ~ Exponential(1)
// Outcome shares: LOGISTIC-NORMAL hierarchy (no simplex boundaries):
//   s_q[k] = softmax(eta_q[k]),  eta_q[k] ~ Normal(mu_q, sigma_q)
//   s_r[k] = softmax(eta_r[k]),  eta_r[k] ~ Normal(mu_r, sigma_r)
//   mu ~ Normal(0, 2), sigma ~ Exponential(1)
// The Dirichlet-on-simplex version had pathological boundary geometry (100%
// max-treedepth saturation at ~9 s/iter); the softmax parameterization keeps
// every parameter natively unconstrained. M12 needs no Bayes factors, so the
// model uses the ~ notation throughout (constants have zero gradient cost).
//
// logit(p_k) = beta_k + alpha_k * t_i + gamma_k * t2_i   (restraint, as M11)
//
// Observation groups:
//   - restraint block: X ~ Bin(n_x, p_k)  [both records first, then x_only]
//   - multi: restraint + complete outcomes → Bin(Y) + Multinomial(cells | Y)
//   - bin:   restraint + only total falls (or nesting violations) → Bin(Y)
//   - mm:    complete outcomes without restraint → marginal versions (p_k mix)
//   - bm:    only total falls without restraint → marginal Bin

data {
  int<lower=1> K;          // number of facility×sector pairs
  int<lower=0> N_both;     // 2524 — for GQ arrays (compute_fit compatibility)
  int<lower=0> N_x_only;   // 242
  int<lower=0> N_y_only;   // 419

  // Restraint block: both-records first (original both order), then x_only
  int<lower=0> N_x;
  array[N_x] int<lower=0> n_x;
  array[N_x] int<lower=0> x;
  array[N_x] int<lower=1> fs_x;
  array[N_x] real t_x;
  array[N_x] real t2_x;

  // Both with complete outcome vector
  int<lower=0> N_multi;
  array[N_multi] int<lower=0> n_y_multi;
  array[N_multi] int<lower=0> x_multi;    // restraint count (for frac)
  array[N_multi] int<lower=0> nx_multi;   // restraint denominator (for frac)
  array[N_multi, 4] int<lower=0> c_multi; // (no outcome, minor, major, not fallen)
  array[N_multi] int<lower=1> fs_multi;
  array[N_multi] int<lower=1> pos_multi;  // position within both (1..N_both)

  // Both with only total falls
  int<lower=0> N_bin;
  array[N_bin] int<lower=0> n_y_bin;
  array[N_bin] int<lower=0> y_bin;
  array[N_bin] int<lower=0> x_bin;
  array[N_bin] int<lower=0> nx_bin;
  array[N_bin] int<lower=1> fs_bin;
  array[N_bin] int<lower=1> pos_bin;      // position within both

  // Complete outcomes without restraint
  int<lower=0> N_mm;
  array[N_mm] int<lower=0> n_y_mm;
  array[N_mm, 4] int<lower=0> c_mm;
  array[N_mm] int<lower=1> fs_mm;
  array[N_mm] real t_mm;
  array[N_mm] real t2_mm;
  array[N_mm] int<lower=1> pos_mm;        // position within y_only

  // Only total falls without restraint
  int<lower=0> N_bm;
  array[N_bm] int<lower=0> n_y_bm;
  array[N_bm] int<lower=0> y_bm;
  array[N_bm] int<lower=1> fs_bm;
  array[N_bm] real t_bm;
  array[N_bm] real t2_bm;
  array[N_bm] int<lower=1> pos_bm;        // position within y_only
}

transformed data {
  // Total falls per record (sum of the three outcome cells)
  array[N_multi] int<lower=0> y_multi;
  array[N_mm] int<lower=0> y_mm;
  for (i in 1:N_multi) {
    y_multi[i] = c_multi[i, 1] + c_multi[i, 2] + c_multi[i, 3];
  }
  for (i in 1:N_mm) {
    y_mm[i] = c_mm[i, 1] + c_mm[i, 2] + c_mm[i, 3];
  }
}

parameters {
  vector[K] beta;                     // facility×sector intercepts (logit)
  vector[K] alpha;                    // facility×sector linear slopes
  vector[K] gamma;                    // facility×sector quadratic coefficients
  vector<lower=0, upper=1>[K] q_fall; // P(fall | restrained), per pair
  vector<lower=0, upper=1>[K] r_fall; // P(fall | not restrained), per pair
  real<lower=0> a_q;                  // Beta hyperparameters for q_fall
  real<lower=0> b_q;
  real<lower=0> a_r;                  // Beta hyperparameters for r_fall
  real<lower=0> b_r;
  array[K] vector[3] eta_q_raw;       // non-centered log-odds (restrained)
  array[K] vector[3] eta_r_raw;       // non-centered log-odds (not restrained)
  vector[3] mu_q;                     // population mean of eta_q
  vector<lower=0>[3] sigma_q;         // between-pair sd of eta_q
  vector[3] mu_r;                     // population mean of eta_r
  vector<lower=0>[3] sigma_r;         // between-pair sd of eta_r
}

transformed parameters {
  // Non-centered parameterization (exact reparameterization — breaks the
  // hierarchical funnel between mu/sigma and the pair-level log-odds)
  array[K] vector[3] s_q;             // outcome composition if restrained
  array[K] vector[3] s_r;             // outcome composition if not restrained
  for (k in 1:K) {
    s_q[k] = softmax(mu_q + sigma_q .* eta_q_raw[k]);
    s_r[k] = softmax(mu_r + sigma_r .* eta_r_raw[k]);
  }
}

model {
  // Priors (~ notation: constants dropped, no Bayes factor needed for M12)
  beta ~ normal(0, 10);
  alpha ~ normal(0, 10);
  gamma ~ normal(0, 5);
  a_q ~ exponential(1);
  b_q ~ exponential(1);
  a_r ~ exponential(1);
  b_r ~ exponential(1);
  q_fall ~ beta(a_q, b_q);
  r_fall ~ beta(a_r, b_r);
  mu_q ~ normal(0, 2);
  mu_r ~ normal(0, 2);
  sigma_q ~ exponential(1);
  sigma_r ~ exponential(1);
  for (k in 1:K) {
    eta_q_raw[k] ~ std_normal();
    eta_r_raw[k] ~ std_normal();
  }

  // Likelihood — restraint
  for (i in 1:N_x) {
    int k = fs_x[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_x[i] + gamma[k] * t2_x[i]);
    x[i] ~ binomial(n_x[i], p_i);
  }

  // Likelihood — both, complete outcomes (factorized)
  for (i in 1:N_multi) {
    int k = fs_multi[i];
    real frac = (nx_multi[i] > 0) ? x_multi[i] * 1.0 / nx_multi[i] : 0.0;
    real pf = frac * q_fall[k] + (1 - frac) * r_fall[k];
    y_multi[i] ~ binomial(n_y_multi[i], pf);
    if (y_multi[i] > 0) {
      vector[3] s_bar = (frac * q_fall[k] * s_q[k]
                         + (1 - frac) * r_fall[k] * s_r[k]) / pf;
      c_multi[i, 1:3] ~ multinomial(s_bar);
    }
  }

  // Likelihood — both, only total falls
  for (i in 1:N_bin) {
    int k = fs_bin[i];
    real frac = (nx_bin[i] > 0) ? x_bin[i] * 1.0 / nx_bin[i] : 0.0;
    real pf = frac * q_fall[k] + (1 - frac) * r_fall[k];
    y_bin[i] ~ binomial(n_y_bin[i], pf);
  }

  // Likelihood — complete outcomes without restraint (marginal)
  for (i in 1:N_mm) {
    int k = fs_mm[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_mm[i] + gamma[k] * t2_mm[i]);
    real pf = p_i * q_fall[k] + (1 - p_i) * r_fall[k];
    y_mm[i] ~ binomial(n_y_mm[i], pf);
    if (y_mm[i] > 0) {
      vector[3] s_mix = (p_i * q_fall[k] * s_q[k]
                         + (1 - p_i) * r_fall[k] * s_r[k]) / pf;
      c_mm[i, 1:3] ~ multinomial(s_mix);
    }
  }

  // Likelihood — only total falls without restraint (marginal)
  for (i in 1:N_bm) {
    int k = fs_bm[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_bm[i] + gamma[k] * t2_bm[i]);
    real pf = p_i * q_fall[k] + (1 - p_i) * r_fall[k];
    y_bm[i] ~ binomial(n_y_bm[i], pf);
  }
}

generated quantities {
  // Population-level 4-cell probability vectors (means across pairs)
  vector[4] q_bar = rep_vector(0, 4);
  vector[4] r_bar = rep_vector(0, 4);
  for (k in 1:K) {
    q_bar[1:3] += q_fall[k] * s_q[k];
    q_bar[4]   += 1 - q_fall[k];
    r_bar[1:3] += r_fall[k] * s_r[k];
    r_bar[4]   += 1 - r_fall[k];
  }
  q_bar /= K;
  r_bar /= K;

  // Key contrasts: total falls and major-outcome falls
  real fall_q = 1 - q_bar[4];
  real fall_r = 1 - r_bar[4];
  real delta_fall = fall_q - fall_r;
  real rr_fall = fall_q / fall_r;
  real delta_major = q_bar[3] - r_bar[3];
  real rr_major = q_bar[3] / r_bar[3];

  // PPC — arrays in the ORIGINAL stan_data ordering (compute_fit compatible)
  array[N_both] int x_rep;
  array[N_both] int y_rep;
  array[N_x_only] int x_only_rep;
  array[N_y_only] int y_only_rep;

  // Restraint: first N_both entries are both-records, then x_only
  for (i in 1:N_x) {
    int k = fs_x[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_x[i] + gamma[k] * t2_x[i]);
    int xr = binomial_rng(n_x[i], p_i);
    if (i <= N_both) {
      x_rep[i] = xr;
    } else {
      x_only_rep[i - N_both] = xr;
    }
  }

  // Both, complete outcomes: total falls from the binomial level
  for (i in 1:N_multi) {
    int k = fs_multi[i];
    real frac = (nx_multi[i] > 0) ? x_multi[i] * 1.0 / nx_multi[i] : 0.0;
    real pf = frac * q_fall[k] + (1 - frac) * r_fall[k];
    y_rep[pos_multi[i]] = binomial_rng(n_y_multi[i], pf);
  }

  // Both, only total falls
  for (i in 1:N_bin) {
    int k = fs_bin[i];
    real frac = (nx_bin[i] > 0) ? x_bin[i] * 1.0 / nx_bin[i] : 0.0;
    real pf = frac * q_fall[k] + (1 - frac) * r_fall[k];
    y_rep[pos_bin[i]] = binomial_rng(n_y_bin[i], pf);
  }

  // Marginal, complete outcomes
  for (i in 1:N_mm) {
    int k = fs_mm[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_mm[i] + gamma[k] * t2_mm[i]);
    real pf = p_i * q_fall[k] + (1 - p_i) * r_fall[k];
    y_only_rep[pos_mm[i]] = binomial_rng(n_y_mm[i], pf);
  }

  // Marginal, only total falls
  for (i in 1:N_bm) {
    int k = fs_bm[i];
    real p_i = inv_logit(beta[k] + alpha[k] * t_bm[i] + gamma[k] * t2_bm[i]);
    real pf = p_i * q_fall[k] + (1 - p_i) * r_fall[k];
    y_only_rep[pos_bm[i]] = binomial_rng(n_y_bm[i], pf);
  }
}
