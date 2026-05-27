# Prepare the package dataset
# Run this script with: source("data-raw/prepare_dataset.R")

devtools::load_all()

restraints_falls <- get_data()

usethis::use_data(restraints_falls, overwrite = TRUE)
