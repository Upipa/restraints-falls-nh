library(tidyverse)
source("R/db_connect.R")

con <- db_connect()

cleaned_data <- tbl(con, "indicare_salute") |>
  filter(anno %in% 2016:2025) |>
  collect() |>
  drop_na(n, d) |>
  filter(
    str_starts(indicatore, "1.1|1.3|1.5|2.1")
  ) |>
  mutate(
    data = ym(str_c(anno, "-", mese))
  ) |>
  group_by(ente, indicatore, id_settore) |>
  mutate(
    outlier_n = !between(
      n,
      quantile(n, 0.25) - 1.5 * IQR(n),
      quantile(n, 0.75) + 1.5 * IQR(n)
    ),
    outlier_d = !between(
      d,
      quantile(d, 0.25) - 1.5 * IQR(d),
      quantile(d, 0.75) + 1.5 * IQR(d)
    )
  ) |>
  filter(!outlier_n, !outlier_d) |>
  mutate(ente = as.integer(factor(ente))) |>
  select(-outlier_n, -outlier_d)

DBI::dbDisconnect(con)

write_rds(cleaned_data, "data/cleaned_data.rds")
