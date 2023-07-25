#!/usr/bin/bash
pwd > install_wd.txt

cwd=$(pwd)

export R_LIBS=$cwd/sRNA_packages_R

mkdir test_results

echo "creating test yaml file"

python create_test_yaml.py

echo "moving yaml into main dir. will replace old template file."

mv sRNA_frag_config.yaml ..

echo "moving out of $(pwd)"

cd ..

python sRNA_fragment.py