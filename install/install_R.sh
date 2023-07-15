#!/usr/bin/bash
Rscript -e 'install.packages("devtools")'

while IFS=" " read -r package version; 
do 
  Rscript -e "devtools::install_version('"$package"', version='"$version"')"; 
done < "requirements_R.txt"

Rscript -e 'devtools::install_github("omarwagih/ggseqlogo")'