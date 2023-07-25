#!/usr/bin/bash

# First start with R packages

mkdir sRNA_packages_R

Rscript install.R

cwd=$(pwd)

export R_LIBS=$cwd/sRNA_packages_R

echo "R packages installed."

# Create conda  environment under name sRNAfrag
eval "$(conda shell.bash hook)"
conda config --set channel_priority disabled
conda env create -f requirements_conda.yml
conda activate sRNAfrag

conda install -c bioconda -y bowtie=1.3.1
conda install -c bioconda -y bowtie2=2.5.1
conda install -c bioconda -y umi_tools=1.1.4
conda install -c bioconda -y adapterremoval=2.3.3
conda install -c bioconda -y bedtools=2.30.0
conda install -c bioconda -y fastqc=0.12.1
conda install -c bioconda -y fastx_toolkit=0.0.14
conda install -c bioconda -y hisat2=2.2.1
conda install -c bioconda -y bedtools=2.30.0
conda install -c bioconda -y subread=2.0.6
conda install -c bioconda -y samtools=1.17

pip install -f requirements_pip.txt

echo "conda installable packages, installed"

pwd > install_wd.txt

echo "creating test yaml file"

python create_test_yaml.py

echo "moving yaml into main dir. will replace old template file."

mv sRNA_frag_config.yaml ..

echo "moving out of $(pwd)"

cd ..

python sRNA_fragment.py