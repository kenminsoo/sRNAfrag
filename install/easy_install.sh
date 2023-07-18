#!/usr/bin/bash

# First start with R packages

Rscript install.R

echo "R packages installed."

# Create conda  environment under name sRNAfrag
eval "$(conda shell.bash hook)"
conda create --name sRNAfrag --file requirements_conda.txt
conda activate sRNAfrag

conda install -c bioconda -y bowtie=1.3.1
conda install -c bioconda -y bowtie2=2.5.1
conda install -c bioconda -y umi_tools=1.1.4
conda install -c bioconda -y adapterremoval=2.3.3
conda install -c bioconda -y bedtools=2.30.0
conda install -c bioconda -y fastqc=0.12.1
conda install -c bioconda -y fastx_toolkit=0.0.14
conda install -c bioconda -y cutadapt=4.4
conda install -c bioconda -y hisat2=2.2.1
conda install -c bioconda -y bedtools=2.30.0

echo "conda installable packages, installed"

pwd > install_wd.txt