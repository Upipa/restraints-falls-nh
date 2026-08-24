#' Prepare the data list for the multi-outcome model
#'
#' Builds the complete Stan data list for the final multi-outcome model
#' (multinomial fall outcomes + per-pair quadratic restraint trend) from the
#' analysis dataset. The function reproduces the preparation pipeline used
#' for the published analysis: pivoting the four indicators, building the
#' multinomial outcome cells by successive differences (the fall indicators
#' are nested, 1.5 within 1.3 within 1.1), splitting records into observation
#' groups, and computing the orthogonalized, standardized quadratic time term.
#'
#' @param data A data frame with the structure of [restraints_falls]: columns
#'   `anno`, `mese`, `ente`, `id_settore`, `indicatore`, `n`, `d`. Defaults
#'   to the packaged dataset.
#'
#' @return A list with two elements:
#'   * `stan_data`: the complete data list for the Stan model
#'     `restraint_multinomial_fac_sec_quad.stan`.
#'   * `meta`: metadata (pair labels, number of pairs `K`, time
#'     standardization constants, group sizes).
#'
#' @details Records with nesting violations (`n_1.3 > n_1.1` or
#'   `n_1.5 > n_1.3` — 59 records in the packaged dataset, probable data-entry
#'   errors) fall back to the binomial groups, which use only the total fall
#'   count. Row ordering may differ from `data-raw/model-development.qmd`;
#'   the model is identical.
#'
#' @export
#'
#' @examples
#' prep <- prepare_model_data()
#' str(prep$meta)
prepare_model_data <- function(data = restraints_falls) {
  required <- c("anno", "mese", "ente", "id_settore", "indicatore", "n", "d")
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols)) {
    stop(
      "Colonne mancanti in `data`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Pivot: una riga per ente x settore x mese, quattro indicatori
  df <- data[, required]
  df$ind <- substr(df$indicatore, 1, 3)
  df$indicatore <- NULL
  wide <- tidyr::pivot_wider(
    df,
    names_from = "ind",
    values_from = c("n", "d"),
    names_sep = "_"
  )

  n_11 <- wide[["n_1.1"]]
  d_11 <- wide[["d_1.1"]]
  n_13 <- wide[["n_1.3"]]
  n_15 <- wide[["n_1.5"]]
  n_21 <- wide[["n_2.1"]]
  d_21 <- wide[["d_2.1"]]

  has_x <- !is.na(n_21) & !is.na(d_21)
  has_11 <- !is.na(n_11) & !is.na(d_11)
  complete <- has_11 &
    !is.na(n_13) &
    !is.na(n_15) &
    (n_13 <= n_11) &
    (n_15 <= n_13)
  keep <- has_x | has_11

  # Tempo: z-score sul pool delle osservazioni utilizzate
  t_raw <- wide$anno + (wide$mese - 0.5) / 12
  t_mean <- mean(t_raw[keep])
  t_sd <- sd(t_raw[keep])
  t <- (t_raw - t_mean) / t_sd

  # Termine quadratico ortogonalizzato (Gram-Schmidt) e standardizzato
  tz <- t[keep]
  orth <- lm(I(tz^2) ~ tz)
  t2_sd <- sd(residuals(orth))
  t2 <- (t^2 - predict(orth, newdata = data.frame(tz = t))) / t2_sd

  # Mappa coppia (ente, id_settore) -> indice intero
  pair_labels <- sort(unique(paste(
    wide$ente[keep],
    wide$id_settore[keep],
    sep = "_"
  )))
  pair_map <- setNames(seq_along(pair_labels), pair_labels)
  k <- unname(pair_map[paste(wide$ente, wide$id_settore, sep = "_")])

  # Gruppi di osservazioni
  idx_both <- which(has_x & has_11)
  idx_x_only <- which(has_x & !has_11)
  idx_y_only <- which(!has_x & has_11)
  complete_b <- complete[idx_both]
  complete_y <- complete[idx_y_only]

  i_multi <- idx_both[complete_b]
  i_bin <- idx_both[!complete_b]
  i_mm <- idx_y_only[complete_y]
  i_bm <- idx_y_only[!complete_y]
  i_restr <- c(idx_both, idx_x_only) # both prima, poi x_only (convenzione GQ)

  # Celle multinomiali: (senza esito, esito minore, esito maggiore, non caduto)
  mk_cells <- function(i) {
    cbind(
      as.integer(n_11[i] - n_13[i]),
      as.integer(n_13[i] - n_15[i]),
      as.integer(n_15[i]),
      as.integer(d_11[i] - n_11[i])
    )
  }

  stan_data <- list(
    K = length(pair_labels),
    N_both = length(idx_both),
    N_x_only = length(idx_x_only),
    N_y_only = length(idx_y_only),
    N_x = length(i_restr),
    n_x = as.integer(d_21[i_restr]),
    x = as.integer(n_21[i_restr]),
    fs_x = k[i_restr],
    t_x = t[i_restr],
    t2_x = t2[i_restr],
    N_multi = length(i_multi),
    n_y_multi = as.integer(d_11[i_multi]),
    x_multi = as.integer(n_21[i_multi]),
    nx_multi = as.integer(d_21[i_multi]),
    c_multi = mk_cells(i_multi),
    fs_multi = k[i_multi],
    pos_multi = which(complete_b),
    N_bin = length(i_bin),
    n_y_bin = as.integer(d_11[i_bin]),
    y_bin = as.integer(n_11[i_bin]),
    x_bin = as.integer(n_21[i_bin]),
    nx_bin = as.integer(d_21[i_bin]),
    fs_bin = k[i_bin],
    pos_bin = which(!complete_b),
    N_mm = length(i_mm),
    n_y_mm = as.integer(d_11[i_mm]),
    c_mm = mk_cells(i_mm),
    fs_mm = k[i_mm],
    t_mm = t[i_mm],
    t2_mm = t2[i_mm],
    pos_mm = which(complete_y),
    N_bm = length(i_bm),
    n_y_bm = as.integer(d_11[i_bm]),
    y_bm = as.integer(n_11[i_bm]),
    fs_bm = k[i_bm],
    t_bm = t[i_bm],
    t2_bm = t2[i_bm],
    pos_bm = which(!complete_y)
  )

  meta <- list(
    pair_labels = pair_labels,
    K = length(pair_labels),
    t_mean = t_mean,
    t_sd = t_sd,
    t2_sd = t2_sd,
    group_sizes = c(
      N_x = length(i_restr),
      N_multi = length(i_multi),
      N_bin = length(i_bin),
      N_mm = length(i_mm),
      N_bm = length(i_bm)
    )
  )

  list(stan_data = stan_data, meta = meta)
}
