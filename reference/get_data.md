# Retrieve and clean restraint-fall data from the database

Connects to the Indicare Salute Lab database, extracts indicators
related to physical restraints and falls (2016-2025), removes outliers
using the IQR method, and anonymizes facility identifiers.

## Usage

``` r
get_data(con = NULL)
```

## Arguments

- con:

  A `DBIConnection` object. If `NULL` (the default), a new connection is
  opened via
  [`db_connect()`](https://upipa.github.io/restraints-falls-nh/reference/db_connect.md)
  and closed on exit.

## Value

A tibble with columns: `ente` (anonymized integer ID), `anno`,
`indicatore`, `n`, `d`, `id_settore`, `data`.
