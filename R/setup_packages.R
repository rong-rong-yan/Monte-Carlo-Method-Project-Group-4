############################################################
# setup_packages.R
# One-time install of all packages required by the pipeline.
# Run interactively:
#   Rscript R/setup_packages.R
# or source() it from an R session.
############################################################

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "yaml",
  "readr",
  "pROC",
  "glmnet",
  "xgboost",
  "fastICA",
  "nnls",
  "BiocManager"   # needed to fetch Biobase
)

is_installed <- function(pkg) {
  pkg %in% rownames(installed.packages())
}

is_loadable <- function(pkg) {
  res <- tryCatch(
    suppressMessages(suppressWarnings(
      requireNamespace(pkg, quietly = TRUE)
    )),
    error = function(e) FALSE
  )
  isTRUE(res)
}

############################################################
# 1. CRAN packages
############################################################

to_install <- cran_pkgs[!sapply(cran_pkgs, is_installed)]

if (length(to_install) > 0) {
  cat("\nInstalling CRAN packages:\n")
  print(to_install)
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

############################################################
# 2. Biobase from Bioconductor
#    NMF depends on Biobase, which is NOT on CRAN.
#    Must be installed via BiocManager.
############################################################

if (!is_installed("Biobase")) {
  cat("\nInstalling Biobase from Bioconductor...\n")
  BiocManager::install("Biobase", update = FALSE, ask = FALSE)
}

############################################################
# 3. NMF from CRAN (after Biobase is in place)
############################################################

if (!is_installed("NMF")) {
  cat("\nInstalling NMF from CRAN...\n")
  install.packages("NMF", repos = "https://cloud.r-project.org")
}

############################################################
# 4. Verify everything is installed AND loadable.
#    "Installed but not loadable" usually means a dep is missing.
############################################################

all_pkgs <- c(cran_pkgs, "Biobase", "NMF")

status_df <- data.frame(
  package   = all_pkgs,
  installed = sapply(all_pkgs, is_installed),
  loadable  = sapply(all_pkgs, is_loadable),
  stringsAsFactors = FALSE,
  row.names = NULL
)

cat("\nInstall check:\n")
print(status_df)

if (any(!status_df$loadable)) {
  cat("\nWARNING: some packages are not loadable. ",
      "Look at rows where loadable = FALSE above.\n", sep = "")
} else {
  cat("\nAll packages installed and loadable.\n")
}