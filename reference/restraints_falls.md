# Restraint and fall indicators in nursing homes

Facility-level data on physical restraint prevalence and fall rates from
nursing homes in the province of Trento, Italy (2016-2025). Extracted
from the Indicare Salute Lab quality monitoring system and anonymized.

## Usage

``` r
restraints_falls
```

## Format

A tibble with 11,878 rows and 8 columns:

- anno:

  Year of observation

- mese:

  Month of observation

- ente:

  Anonymized facility identifier (integer)

- indicatore:

  Indicator code. One of:

  1.1

  :   Residenti caduti nel periodo di riferimento (residents who fell
      during the reference period)

  1.3

  :   Cadute con esito (falls with injury)

  1.5

  :   Cadute con esito maggiore (falls with major injury)

  2.1

  :   Persone soggette a contenzione fisica (residents subject to
      physical restraint)

- id_settore:

  Care unit identifier within the facility

- n:

  Numerator. Meaning depends on `indicatore`:

  1.1

  :   Number of residents with at least one registered fall during the
      reference period

  1.3

  :   Number of falls with injury (minor + major) during the reference
      period

  1.5

  :   Number of falls with major injury during the reference period

  2.1

  :   Number of residents with an active physical restraint prescription
      on the index day

- d:

  Denominator. Meaning depends on `indicatore`:

  1.1

  :   Number of residents present during the reference period

  1.3

  :   Total number of registered falls during the reference period

  1.5

  :   Total number of registered falls during the reference period

  2.1

  :   Number of residents present on the index day

- data:

  Date (year-month, as Date object)

## Source

Indicare Salute Lab, UPIPA (Unione Provinciale Istituzioni Per
l'Assistenza), Trento, Italy.

## Details

The index day (giorno indice) for indicator 2.1 is a randomly chosen day
within the reference month, common to all participating facilities
(i.e., all facilities collect data on the same day, communicated
centrally by UPIPA). The actual date is not recorded in the dataset.

Note that `n/d` represents different quantities depending on the
indicator: for 1.1 and 2.1 it is a prevalence (proportion of residents),
while for 1.3 and 1.5 it is a proportion of falls with a given severity.
