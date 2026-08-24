# Rende disponibile renv e punta alla renv del progetto padre
# (restraints-falls-nh), così il rendering Quarto vede la library del
# progetto (pacchetto restraintsfalls e dipendenze).
# NOTA: source() da solo attiverebbe una renv vuota con root in manuscript/;
# renv::load("..") riporta i .libPaths al progetto corretto.
source("../renv/activate.R")
renv::load(project = "..")
