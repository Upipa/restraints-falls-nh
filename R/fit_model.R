#' Fit the multi-outcome restraint-falls model
#'
#' Compiles and fits the final Bayesian hierarchical model (M12) with NUTS
#' via cmdstanr: multinomial fall outcomes (total / minor / major) linked to
#' the observed restraint count, with per-pair quadratic restraint trends and
#' hierarchical fall rates and outcome shares.
#'
#' @param data A data frame with the structure of [restraints_falls]. Ignored
#'   if `prepared` is supplied.
#' @param prepared Optional output of [prepare_model_data()]; computed from
#'   `data` if not supplied.
#' @param seed Random seed passed to the sampler. The default reproduces the
#'   published fit.
#' @param chains Number of MCMC chains.
#' @param iter_warmup Warmup iterations per chain.
#' @param iter_sampling Sampling iterations per chain.
#' @param adapt_delta Target acceptance probability. The default (0.95)
#'   reproduces the published fit exactly. The Stan default (0.8) is likely
#'   sufficient and considerably faster — worth trying if you do not need an
#'   exact reproduction.
#' @param quick If `TRUE`, runs a single short chain (200 warmup + 100
#'   sampling iterations) to verify that the model compiles and samples,
#'   before committing to the full fit. Overrides `chains`, `iter_warmup`
#'   and `iter_sampling`.
#' @param ... Further arguments passed to the `sample` method of the
#'   compiled model.
#'
#' @details Fitting the full model takes several hours (about 6 h for the
#'   published fit) and requires a C++ toolchain (Rtools on Windows) to
#'   compile the Stan model. cmdstanr is not on CRAN; install it with
#'   `install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))`.
#'
#' @return A `CmdStanMCMC` object with the posterior draws.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- fit_model(quick = TRUE) # quick check (minutes)
#' fit <- fit_model() # full fit, reproduces the published analysis (~6 h)
#' }
fit_model <- function(
  data = restraints_falls,
  prepared = prepare_model_data(data),
  seed = 919214260,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  quick = FALSE,
  ...
) {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "Il pacchetto 'cmdstanr' non e' installato. Installalo con:\n",
      '  install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))',
      call. = FALSE
    )
  }
  if (quick) {
    chains <- 1
    iter_warmup <- 200
    iter_sampling <- 100
  }
  model_file <- system.file(
    "stan",
    "restraint_multinomial_fac_sec_quad.stan",
    package = "restraintsfalls"
  )
  if (!file.exists(model_file)) {
    stop("File Stan non trovato nel pacchetto: ", model_file, call. = FALSE)
  }
  mod <- cmdstanr::cmdstan_model(model_file)
  mod$sample(
    data = prepared$stan_data,
    seed = seed,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    ...
  )
}
