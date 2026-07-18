#!/usr/bin/env bash
# Environment build for the AWS validation run (Ubuntu 24.04 "noble").
# Installs R 4.5 + all simulation packages as fast apt binaries via r2u,
# then a CPU TensorFlow/keras environment. Idempotent; fails loudly.
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive

echo "=== system packages ==="
sudo apt-get update -y
sudo apt-get install -y --no-install-recommends \
  ca-certificates curl gnupg git rsync tmux htop python3-venv python3-pip python3-dev

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
# The simulation uses the keras v2 API ('keras' CRAN package). Lessons from
# the 2026-07 AWS validation, encoded here:
# - keras::install_keras() pins TF 2.15, which has no wheels for Ubuntu
#   24.04's Python 3.12 — build the venv manually with TF 2.16.
# - Python 3.12 removed distutils; the setuptools .pth shim is NOT processed
#   by R's embedded interpreter, so TF import fails under reticulate even
#   though it works at the CLI. Use Python 3.11 (deadsnakes), where distutils
#   is stdlib.
# - tf-keras (the legacy-API shim TF 2.16 needs) must be version-pinned or it
#   drags in an incompatible protobuf.
# - Runners need TF_USE_LEGACY_KERAS=1 at runtime (run_validation.sh sets it).
sudo add-apt-repository -y ppa:deadsnakes/ppa
sudo apt-get install -y -q python3.11 python3.11-venv python3.11-dev
python3.11 -m venv "$HOME/.virtualenvs/r-tensorflow"
"$HOME/.virtualenvs/r-tensorflow/bin/pip" install -q -U pip
"$HOME/.virtualenvs/r-tensorflow/bin/pip" install -q \
  "tensorflow-cpu==2.16.2" "tf-keras==2.16.0" "protobuf>=3.20.3,<5"
TF_USE_LEGACY_KERAS=1 Rscript -e '
library(keras)
ok <- tryCatch(is_keras_available(), error = function(e) FALSE)
if (!isTRUE(ok)) stop("KERAS/TENSORFLOW NOT AVAILABLE — inspect install_keras output")
cat("keras backend OK; TF version:", as.character(tensorflow::tf_version()), "\n")
'

echo "=== smoke: every Analyse dependency loads ==="
TF_USE_LEGACY_KERAS=1 Rscript -e '
suppressMessages({
  library(tidyverse); library(MASS); library(psych); library(rpart); library(ipred)
  library(randomForest); library(gbm); library(glmnet); library(nnet)
  library(ranger); library(xgboost); library(SuperLearner); library(WeightIt)
  library(dbarts); library(survey); library(Hmisc); library(SimDesign)
})
cat("all libraries load\n"); print(sessionInfo()$R.version$version.string)
'
echo "=== bootstrap complete ==="
