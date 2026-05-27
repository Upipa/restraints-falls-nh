#' Retrieve and clean restraint-fall data from the database
#'
#' Connects to the Indicare Salute Lab database, extracts indicators related to
#' physical restraints and falls (2016-2025), removes outliers using the IQR
#' method, and anonymizes facility identifiers.
#'
#' @param con A `DBIConnection` object. If `NULL` (the default), a new
#'   connection is opened via [db_connect()] and closed on exit.
#'
#' @return A tibble with columns: `ente` (anonymized integer ID), `anno`,
#'   `indicatore`, `n`, `d`, `id_settore`, `data`.
#'
#' @export
get_data <- function(con = NULL) {
  close_on_exit <- is.null(con)
  if (close_on_exit) {
    con <- db_connect()
    on.exit(DBI::dbDisconnect(con))
  }

  dplyr::tbl(con, "indicare_salute") |>
    dplyr::filter(anno %in% 2016:2025) |>
    dplyr::collect() |>
    tidyr::drop_na(n, d) |>
    dplyr::filter(
      stringr::str_starts(indicatore, "1.1|1.3|1.5|2.1")
    ) |>
    dplyr::mutate(
      data = lubridate::ym(paste(anno, mese, sep = "-"))
    ) |>
    dplyr::group_by(ente, indicatore, id_settore) |>
    dplyr::mutate(
      outlier_n = !dplyr::between(
        n,
        stats::quantile(n, 0.25) - 1.5 * stats::IQR(n),
        stats::quantile(n, 0.75) + 1.5 * stats::IQR(n)
      ),
      outlier_d = !dplyr::between(
        d,
        stats::quantile(d, 0.25) - 1.5 * stats::IQR(d),
        stats::quantile(d, 0.75) + 1.5 * stats::IQR(d)
      )
    ) |>
    dplyr::filter(!outlier_n, !outlier_d) |>
    dplyr::ungroup() |>
    dplyr::select(-outlier_n, -outlier_d) |>
    dplyr::mutate(ente = as.integer(factor(ente)))
}
