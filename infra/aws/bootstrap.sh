#!/usr/bin/env bash
# Environment build for the AWS validation run (Ubuntu 24.04 "noble").
# Installs R 4.5 + all simulation packages as fast apt binaries via r2u,
# then a CPU TensorFlow/keras environment. Idempotent; fails loudly.
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive

echo "=== system packages ==="
sudo apt-get update -y
sudo apt-get install -y --no-install-recommends \
  ca-certificates curl gnupg git rsync tmux htop python3-venv

echo "=== r2u (CRAN apt binaries + bspm) ==="
curl -fsSL https://raw.githubusercontent.com/eddelbuettel/r2u/master/inst/scripts/add_cranapt_noble.sh | sudo bash
sudo apt-get install -y r-base-core

echo "=== R packages (binary via r2u) ==="
sudo Rscript -e '
pkgs <- c("here","tidyverse","MASS","Matrix","psych","rpart","ipred",
          "randomForest","gbm","glmnet","nnet","ranger","xgboost",
          "SuperLearner","WeightIt","dbarts","survey","Hmisc","SimDesign",
          "qs2","reticulate","tensorflow","keras")
install.packages(pkgs)
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) stop("MISSING PACKAGES: ", paste(missing, collapse=", "))
cat("all", length(pkgs), "packages installed\n")
'

echo "=== TensorFlow (CPU) for keras ==="
# The simulation uses the keras v2 API ('keras' CRAN package). install_keras
# creates its own virtualenv that reticulate auto-discovers.
Rscript -e 'keras::install_keras(method = "virtualenv")'
Rscript -e '
library(keras)
ok <- tryCatch(is_keras_available(), error = function(e) FALSE)
if (!isTRUE(ok)) stop("KERAS/TENSORFLOW NOT AVAILABLE — inspect install_keras output")
cat("keras backend OK; TF version:", as.character(tensorflow::tf_version()), "\n")
'

echo "=== smoke: every Analyse dependency loads ==="
Rscript -e '
suppressMessages({
  library(tidyverse); library(MASS); library(psych); library(rpart); library(ipred)
  library(randomForest); library(gbm); library(glmnet); library(nnet)
  library(ranger); library(xgboost); library(SuperLearner); library(WeightIt)
  library(dbarts); library(survey); library(Hmisc); library(SimDesign)
})
cat("all libraries load\n"); print(sessionInfo()$R.version$version.string)
'
echo "=== bootstrap complete ==="
