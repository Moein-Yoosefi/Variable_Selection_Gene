# Installation Script for Variable Selection Gene Project
# Run this script to install all required and optional packages

cat("=== Installing Required Packages for Variable Selection Gene ===\n\n")

# Function to check and install packages
install_if_missing <- function(package_name, repo = "CRAN") {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    cat(sprintf("Installing %s from %s...\n", package_name, repo))

    if (repo == "CRAN") {
      install.packages(package_name, dependencies = TRUE)
    } else if (repo == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(package_name, update = FALSE)
    }

    cat(sprintf("%s installed successfully!\n\n", package_name))
  } else {
    cat(sprintf("%s is already installed.\n", package_name))
  }
}

# Core CRAN packages
cat("Installing core CRAN packages...\n")
cran_packages <- c(
  "glmnet",
  "randomForest",
  "caret",
  "MASS",
  "ggplot2",
  "pheatmap"
)

for (pkg in cran_packages) {
  install_if_missing(pkg, "CRAN")
}

cat("\n")

# Optional CRAN packages
cat("Installing optional CRAN packages...\n")
optional_cran <- c("Boruta")

for (pkg in optional_cran) {
  tryCatch({
    install_if_missing(pkg, "CRAN")
  }, error = function(e) {
    cat(sprintf("Warning: Could not install optional package %s\n", pkg))
  })
}

cat("\n")

# Bioconductor packages
cat("Installing Bioconductor packages (optional)...\n")
bioc_packages <- c("preprocessCore", "impute")

for (pkg in bioc_packages) {
  tryCatch({
    install_if_missing(pkg, "Bioconductor")
  }, error = function(e) {
    cat(sprintf("Warning: Could not install Bioconductor package %s\n", pkg))
    cat("You can skip this if you don't need advanced normalization features.\n")
  })
}

cat("\n=== Installation Complete ===\n")
cat("All required packages have been installed.\n")
cat("You can now run the example analysis with: source('examples/example_analysis.R')\n")

# Print session info
cat("\n=== Session Info ===\n")
sessionInfo()
