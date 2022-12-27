#!/bin/bash

module load lang/Anaconda3/2022.05

echo "please enter location of text file with conda packages to install separated by newlines"
read packages

echo "please enter the conda environment you'd like to install into"
read environment

eval "$(conda shell.bash hook)"
conda activate $environment

while read p; do
	conda install -y $p
done <$packages
