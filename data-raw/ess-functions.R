# =====================================================================
# ESS bulk e tail — Vehtari, Gelman, Simpson, Carpenter, Bürkner (2021)
# "Rank-normalization, folding, and localization: An improved R-hat"
#
# Implementazione allineata all'algoritmo di posterior:::.ess (verificata
# contro posterior::ess_bulk / ess_tail sui draw di questo progetto).
# Input: vettore numerico (una catena) oppure matrice draw × catene.
# =====================================================================

# Normalizza l'input in una matrice numerica draw × catene,
# eliminando classi particolari (draws_array, tibble, ...)
as_chain_matrix <- function(x) {
  if (is.null(dim(x))) {
    return(matrix(as.numeric(x), ncol = 1))
  }
  d <- dim(x)
  if (length(d) == 3 && d[3] == 1) {
    d <- d[1:2] # draws_array con una sola variabile: iterazioni x catene
  }
  if (length(d) != 2) {
    stop("Fornire un vettore (una catena) o una matrice draw x catene.")
  }
  matrix(as.numeric(x), nrow = d[1], ncol = d[2])
}

# Motore comune: ESS da autocorrelazione con sequenza iniziale positiva e
# monotona di Geyer, cappio a S*log10(S). Rispecchia posterior:::.ess.
ess_engine <- function(x) {
  x <- as_chain_matrix(x)
  n <- nrow(x) # draw per catena
  m <- ncol(x) # catene
  S <- n * m
  if (n < 3) {
    return(NA_real_)
  }

  # Autocovarianza per catena (acf usa il denominatore n, come in Stan)
  acov <- vapply(
    seq_len(m),
    function(j) {
      v <- x[, j]
      a <- as.numeric(acf(v, plot = FALSE, lag.max = n - 1)$acf)
      a * sum((v - mean(v))^2) / n
    },
    numeric(n)
  )
  acov_means <- rowMeans(acov) # media per LAG tra le catene

  mean_var <- acov_means[1] * n / (n - 1) # varianza within (unbiased)
  var_plus <- mean_var * (n - 1) / n
  if (m > 1) {
    var_plus <- var_plus + var(colMeans(x)) # + varianza between (B/n)
  }
  if (var_plus <= 0) {
    return(NA_real_)
  }

  rho_at <- function(lag) 1 - (mean_var - acov_means[lag + 1]) / var_plus

  # Sequenza iniziale positiva di Geyer (coppie: lag 0+1, 2+3, ...)
  rho_t <- numeric(n)
  t <- 0
  rho_even <- 1 # lag 0
  rho_odd <- rho_at(1) # lag 1
  rho_t[1] <- rho_even
  rho_t[2] <- rho_odd
  while (t < n - 5 && (rho_even + rho_odd > 0)) {
    t <- t + 2
    rho_even <- rho_at(t)
    rho_odd <- rho_at(t + 1)
    if ((rho_even + rho_odd) >= 0) {
      rho_t[t + 1] <- rho_even
      rho_t[t + 2] <- rho_odd
    }
  }
  max_t <- t
  if (rho_even > 0) {
    rho_t[max_t + 1] <- rho_even
  }

  # Sequenza monotona: la somma di ciascuna coppia non deve superare la
  # precedente; in caso, la coppia corrente è sostituita dalla media
  t <- 0
  while (t <= max_t - 4) {
    t <- t + 2
    if (rho_t[t + 1] + rho_t[t + 2] > rho_t[t - 1] + rho_t[t]) {
      rho_t[t + 1] <- (rho_t[t - 1] + rho_t[t]) / 2
      rho_t[t + 2] <- rho_t[t + 1]
    }
  }

  tau <- -1 + 2 * sum(rho_t[1:max_t]) + rho_t[max_t + 1]
  tau <- max(tau, 1 / log10(S)) # cappio: ess <= S*log10(S)
  S / tau
}

# ESS bulk: rank-normalizzazione (ranghi sul pool di tutte le catene,
# trasformazione di Blom a scala normale), poi ESS da autocorrelazione.
ess_bulk <- function(x) {
  x <- as_chain_matrix(x)
  r <- rank(as.numeric(x), ties.method = "average")
  z <- qnorm((r - 3 / 8) / (length(r) + 1 / 4)) # Blom: (r - 3/8)/(n + 1/4)
  ess_engine(matrix(z, nrow = nrow(x), ncol = ncol(x)))
}

# ESS tail: ESS delle indicatrici "draw <= q5" e "draw <= q95"
# (già su scala libera, nessuna rank-normalizzazione), minimo delle due.
ess_tail <- function(x, probs = c(0.05, 0.95)) {
  x <- as_chain_matrix(x)
  qs <- quantile(as.numeric(x), probs = probs)
  min(ess_engine(x <= qs[1]), ess_engine(x <= qs[2]))
}

# ESS per quantità derivate calcolate su un draws_df, preservando la
# struttura a catene (come fa fit$summary()). Richiede le colonne
# .chain/.iteration (o chain/iteration dopo janitor::clean_names()).
# Uso:
#   fit$draws(vars, format = "df") |> as_tibble() |> mutate(...) |>
#     ess_derived(c("delta", "ratio"))
ess_derived <- function(df, vars) {
  chain_col <- intersect(c(".chain", "chain"), names(df))[1]
  iter_col <- intersect(c(".iteration", "iteration"), names(df))[1]
  if (is.na(chain_col) || is.na(iter_col)) {
    stop("Servono le colonne .chain e .iteration (usa format = \"df\").")
  }
  purrr::map_dfr(vars, function(v) {
    m <- df |>
      dplyr::select(dplyr::all_of(c(iter_col, chain_col, v))) |>
      tidyr::pivot_wider(
        names_from = dplyr::all_of(chain_col),
        values_from = dplyr::all_of(v)
      ) |>
      dplyr::arrange(.data[[iter_col]]) |>
      dplyr::select(-dplyr::all_of(iter_col)) |>
      as.matrix()
    tibble::tibble(
      variable = v,
      ess_bulk = posterior::ess_bulk(m),
      ess_tail = posterior::ess_tail(m)
    )
  })
}
