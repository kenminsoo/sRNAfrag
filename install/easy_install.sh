#!/usr/bin/bash

# First start with R packages

Rscript install.R

echo "R packages installed."

# Create conda  environment under name sRNAfrag
eval "$(conda shell.bash hook)"
conda create --name sRNAfrag --file requirements_conda.txt
conda activate sRNAfrag

echo "conda installable packages, installed"

# Pip installable packages
pip install -r requirements.txt

echo "pip installable packages, installed"

pwd > install_wd.txt