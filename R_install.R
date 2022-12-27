#the following script will install packages required for any pipeline
#files are read from a text file of packages delmited by newlines
#will not handle errors and is meant to be a part of pipes that have had dependencies resolved previously

if (!require('magrittr', quietly = TRUE))
install.packages('magrittr', repos = "http://cran.us.r-project.org")
library('magrittr')

local_dir <- '/home/kenmn/R/x86_64-pc-linux-gnu-library/4.1'

if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager", lib=local_dir)
BiocManager::install(version = "3.14", lib=local_dir)

library(BiocManager, lib.loc=local_dir)

#print("enter path of txt file with r packages to install")
#packages <- readline()

con <- file('R_packages2.txt', "r")

already_installed <- list()
installed <- list()
manual <- list()

#helpful if on cluster working with unwritable R library, recreate bioconductor on home storage)

while(TRUE) {
	line = readLines(con, n =1)
	if (length(line) == 0) {
		break
}
	#check if package installed and removes white space
	
	package_name <- gsub(" ", "", line)

	status <- sapply(package_name, require, character.only=TRUE)
	if (status == TRUE){
		print("package has already been installed")
		already_installed <- append(already_installed, line)
		next
}
	install.packages(gsub(" ", "", line), repos = "http://cran.us.r-project.org")
	
	status <- sapply(package_name, require, character.only=TRUE)
		
	if (status == TRUE){
		print(paste(line, "installed"))
		installed <- append(installed, line)
		next
} #now check if false
	else if (status == FALSE){
		print(paste(line, "not installed, retrying with bioconductor"))
		
		BiocManager::install(gsub(" ", "", line))

		status_bio <- sapply(package_name, require, character.only=TRUE)

		if (status_bio == FALSE){
                	print(paste(line, "please reinstall manually"))
			manual <- append(manual, line)
} 		else {
			installed <- append(installed, line)

}
}
}

