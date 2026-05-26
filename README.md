# restraints-falls-nh

Research project investigating the association between physical restraint use and falls in nursing homes, using longitudinal quality monitoring data from the province of Trento (Italy), 2016–2025.

## Overview

Data come from **Indicare Salute**, a quality monitoring system managed by UPIPA (Unione Provinciale Istituzioni Per l'Assistenza). The analysis covers facility-level aggregated indicators collected monthly across nursing homes in the province.

The main research question is whether physical restraint use is associated with a reduction in falls — and whether any such association is clinically meaningful.

## Methods

Bayesian binomial models estimated with Stan, linking the logit of restraint prevalence to the logit of fall rates, with a temporal trend component. Models are compared against an independence baseline.

## Status

Work in progress. Target journal: *BMC Geriatrics*.
