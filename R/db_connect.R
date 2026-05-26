#' Connect to the Indicare Salute database
#'
#' Opens a connection to the `IndicareSaluteLab` Azure SQL database via ODBC.
#' Credentials are read from environment variables by default and can be set
#' in `.Renviron`.
#'
#' @param uid Character. Database username. Defaults to the `UID` environment variable.
#' @param pwd Character. Database password. Defaults to the `PWD` environment variable.
#'
#' @return A `DBIConnection` object.
#'
#' @examples
#' \dontrun{
#' con <- db_connect()
#' DBI::dbListTables(con)
#' DBI::dbDisconnect(con)
#' }
db_connect <- function(uid = Sys.getenv("UID"), pwd = Sys.getenv("PWD")) {
  DBI::dbConnect(
    odbc::odbc(),
    Driver = "ODBC Driver 18 for SQL Server",
    Server = "upipa-acs.database.windows.net",
    Database = "IndicareSaluteLab",
    UID = uid,
    PWD = pwd,
    Port = 1433
  )
}
