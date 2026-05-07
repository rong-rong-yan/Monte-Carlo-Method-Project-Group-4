############################################################
# setup_packages.R
# One-time install of all packages required by the pipeline.
#
# Usage:
#   Rscript R/setup_packages.R
#
# This script is idempotent: it skips packages that are already
# installed, installs Bioconductor and CRAN packages in the
# correct order, verifies every package can be loaded, and
# fails loudly with a non-zero exit code if anything is broken.
############################################################

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "scales",        # used by Visualization/visualize_abalation.R
  "yaml",
  "readr",
  "pROC",
  "glmnet",
  "xgboost",
  "fastICA",
  "nnls",
  "BiocManager"    # needed to fetch Biobase
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
# 2. Biobase from Bioconductor (required by NMF)
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
# 4. Verify and report
############################################################

all_pkgs <- c(cran_pkgs, "Biobase", "NMF")

status_df <- data.frame(
  package   = all_pkgs,
  installed = sapply(all_pkgs, is_installed),
  loadable  = sapply(all_pkgs, is_loadable),
  version   = sapply(all_pkgs, function(p) {
    if (is_installed(p)) {
      as.character(packageVersion(p))
    } else {
      NA_character_
    }
  }),
  stringsAsFactors = FALSE,
  row.names = NULL
)

cat("\nR version: ", R.version.string, "\n", sep = "")
cat("Platform:  ", R.version$platform,   "\n\n", sep = "")
cat("Install check:\n")
print(status_df)

############################################################
# 5. Fail loudly if anything is broken
############################################################

broken <- status_df$package[!status_df$loadable]

if (length(broken) > 0) {
  cat(
    "\nERROR: the following packages are not loadable:\n  ",
    paste(broken, collapse = ", "),
    "\n\nThe pipeline will not run until these are fixed.\n",
    "Common causes:\n",
    "  - missing system libraries (a C/C++/Fortran compiler on Mac/Linux)\n",
    "  - a failed Bioconductor install (Biobase blocking NMF)\n",
    "  - a network failure during install.packages()\n\n",
    "Try re-running this script. If it still fails, check the\n",
    "install.packages() output above for the underlying error.\n",
    sep = ""
  )
  quit(status = 1)
}

cat("\nAll packages installed and loadable. Environment is ready.\n")
