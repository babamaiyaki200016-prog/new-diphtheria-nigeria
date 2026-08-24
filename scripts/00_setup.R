# ------------------------------------------------------------------------------
# 00_setup.R
# Nigeria Diphtheria Outbreak Modeling -- package setup
#
# Run once. Installs anything missing, then loads every package used across
# the pipeline (01-08). Individual scripts assume this has already been run.
# ------------------------------------------------------------------------------

required_pkgs <- c(
  "tidyverse",   # data wrangling + ggplot2
  "lubridate",   # dates
  "readxl",      # read the Kano line list
  "binom",       # Wilson CIs for CFR
  "EpiEstim",    # Rt estimation
  "deSolve",     # ODE solver for the SEIR model
  "scales",      # plot label formatting
  "knitr"        # simple table printing
)

missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

invisible(lapply(required_pkgs, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# Scripts are always run with the working directory set to the project root
# (see README), so paths below are relative to getwd().
paths <- list(
  data      = "data",
  scripts   = "scripts",
  figures   = "figures",
  reports   = "reports",
  manuscript = "manuscript"
)

for (p in paths) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

cat("Setup complete. Working directory:", getwd(), "\n")
