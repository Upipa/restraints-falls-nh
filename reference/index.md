# Package index

## Data

Functions for retrieving and preparing the analysis dataset

- [`db_connect()`](https://upipa.github.io/restraints-falls-nh/reference/db_connect.md)
  : Connect to the Indicare Salute Lab database
- [`get_data()`](https://upipa.github.io/restraints-falls-nh/reference/get_data.md)
  : Retrieve and clean restraint-fall data from the database

## Dataset

The cleaned and anonymized analysis dataset

- [`restraints_falls`](https://upipa.github.io/restraints-falls-nh/reference/restraints_falls.md)
  : Restraint and fall indicators in nursing homes

## Model

Prepare the data, fit the final multi-outcome model, load the fit

- [`prepare_model_data()`](https://upipa.github.io/restraints-falls-nh/reference/prepare_model_data.md)
  : Prepare the data list for the multi-outcome model
- [`fit_model()`](https://upipa.github.io/restraints-falls-nh/reference/fit_model.md)
  : Fit the multi-outcome restraint-falls model
- [`load_fit()`](https://upipa.github.io/restraints-falls-nh/reference/load_fit.md)
  : Load the posterior draws of the final model
