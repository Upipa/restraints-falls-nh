# Connect to the Indicare Salute Lab database

Opens a connection to the `IndicareSaluteLab` Azure SQL database via
ODBC. Credentials are read from environment variables by default and can
be set in `.Renviron`.

## Usage

``` r
db_connect(uid = Sys.getenv("UID"), pwd = Sys.getenv("PWD"))
```

## Arguments

- uid:

  Character. Database username. Defaults to the `UID` environment
  variable.

- pwd:

  Character. Database password. Defaults to the `PWD` environment
  variable.

## Value

A `DBIConnection` object.

## Examples

``` r
if (FALSE) { # \dontrun{
con <- db_connect()
DBI::dbListTables(con)
DBI::dbDisconnect(con)
} # }
```
