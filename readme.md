# sRNAfrag - A Pipeline to Analyze Fragmentation in sRNA-seq data 

```
MIT License

Copyright (c) [2023] [Ken Nakatsu]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## This repository is associated:
Ken Nakatsu, Mayumi Jijiwa, Vedbar Khadka, Masaki Nasu, Youping Deng, sRNAfrag: a pipeline and suite of tools to analyze fragmentation in small RNA sequencing data, Briefings in Bioinformatics, Volume 25, Issue 1, January 2024, bbad515, https://doi.org/10.1093/bib/bbad515

## Overview
Welcome to our sRNAfrag github repo. sRNAfrag is a pipeline that analyzes small RNA fragmentation in small RNA-seq data. On this page, installation and basic usage can be found. Please consult the wiki attached to this repo for more information. The wiki can be found by navigating to the top bar on the github page and clicking on the "Wiki" name which has the label of the book.

## Installation & Environment
sRNAfrag was developed with R version >4.1.3 and Python 3.10.0. It was developed on MacOS Monterey 12.2.1 M1 64GB, and tested on the same computer. 

### **Installing with Conda**

#### *Compatible Systems*
MAC-OS x86_64

UBUNTU x86_64

IF you are working on an ARM infrastructure on MACOS, please download the x86_64 version of Conda. Unfortunately, if you are working on an ARM infrastructure on Linux, then it will not work. If your HPC has conda support, you should be able to install the dependencies of sRNAfrag. 

#### Prerequisite
Install R, please consult their website. [R](https://cran.r-project.org/)

Install conda, please consult their website. [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

The quick install process takes approximately 20 minutes and will test for successful installation. Run the following commands to install.
```
git clone https://github.com/kenminsoo/sRNAfrag
cd sRNAfrag/install
sh easy_install.sh
```
This creates a conda environment that has all the required package. This makes it easy to use as all functions can be accessed from the command line in accordance with the documentation. 

Please check the `install/test_results/out/tables` directory after the install. There should be 7 tables. 

#### Note about R!
Sometimes, your package folder will not be in R lib path. If this occurs, enter this:
```
cd install
cwd=$(pwd)
export R_LIBS=$cwd/sRNA_packages_R
```

This will add it back to your libPaths, and then you can run your pipeline as needed. 

### Installing Dependencies Manually

#### **PYTHON**

Users should use pip python package manager to download all packages. To do so, they can run the following command.

```
cd install
pip install -f requirements_pip.txt
```

Behavior: Should install all packages. Ensure that you have python installed. In general, pip comes with your installation of Python. 

#### **R**

To accomodate for HPC systems that do not allow users to have write access to the packages that R will use, sRNAfrag creates its own directory for R packages. If the user chooses to install packages manually, they can do the following:

Create a folder of sRNAfrag packages with the install directory:

1. First, create a directory with the following name.

```
cd install
mkdir sRNA_packages_R
```

2. Then, run
```
Rscript install.R
```

3. Add the directory to your R Path such that R searches for libraries within your custom directory.

```
cwd=$(pwd)
export R_LIBS=$cwd/sRNA_packages_R
```

Behavior: Will install all R packages within the install directory.

Install R packages to the standard user library directory:

1. Run 
```
Rscript install_write.R
```

Behavior: Will simply install R packages normally. 

#### **BIO TOOLS** (i.e. Bowtie)

Please add all tools to your path.

[bowtie=1.3.1](https://bowtie-bio.sourceforge.net/manual.shtml) - REQUIRED

[bedtools=2.30.0](https://bedtools.readthedocs.io/en/latest/content/installation.html) - REQUIRED

[fastx_toolkit=0.0.14](http://hannonlab.cshl.edu/fastx_toolkit/download.html) - REQUIRED

[subread=2.0.6](https://subread.sourceforge.net/SubreadUsersGuide.pdf) - REQUIRED

[samtools=1.17](http://www.htslib.org/) - REQUIRED

[umi_tools=1.1.4](https://umi-tools.readthedocs.io/en/latest/INSTALL.html) - OPTIONAL (iff need to remove UMIs)

[adapterremoval=2.3.3](https://adapterremoval.readthedocs.io/en/stable/installation.html) OPTIONAL (iff need to remove adapters)

[fastqc=0.12.1](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - OPTIONAL (iff need to run quality check module)

[hisat2=2.2.1](http://daehwankimlab.github.io/hisat2/download/) - OPTIONAL (iff you want to use a certain command in the scripts)

[bowtie2=2.5.1]() - OPTIONAL (iff you want to use a certain command in the scripts)

[ViennaRNA](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/install.html) - OPTIONAL (iff you want to create secondary stuctures with RNAFold)

## Testing The Pipeline
If sRNAfrag is installed using `quick_install.sh` it will run this test automatically. 

If you would like to test this pipeline after manually installing dependencies or running in the singularity container:

```
cd install
sh test_pipeline.sh
```

*THIS WILL OVERRIDE YOUR CURRENT sRNA_frag_config.yaml FILE!*

## Basic Usage
1. Modify the sRNA_frag_config.yaml file to fit your needs. 

2. If you would like to find out of space maps (an important component, in my opinion) please build an index with your reference genome.
```
bowtie-build <ref_genome> <prefix>
```

3. Add the path to the YAML file in the P2 section: `/ref/index/loc/prefix`

4. Ensure that all fastq files are gzipped.
```
cd sample_dir
gzip *.fastq
```

5. Then run
```
cd sRNAfrag  \\ installed dir
python sRNA_fragment.py \\ Runs pipeline
```

## Config Variables

NOTE! All directories must not be ended with a `/` (i.e. `X/Y/Z` and not `X/Y/Z/`)
NOTE!! Make sure there is a space between the colon and your entry.
NOTE!!! While not required, it is best to put strings in quotations.

module_options:

  P1:

    bool: Set as true of fasle depending on if you want to run this module.

    prefix: A prefix for all IDs assigned to fragments.

    umi_removal: Uses regex and the unique method from UMI tools. 
      bool: Set true if want to remove Qiagen UMI.
      regex: ".+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12}).+"
    
    adapter_removal: Uses AdapterRemoval to remove adapters.
      bool: Set true if want to remove adapters.
      sequence: Provide adapter sequence. 

    annotation_options: This takes in the processed annotation files that can be found in the scripts directory. 
      location: "Give the full path to the location of the processed annotation file"
    
    trim5p3p: Trim n5 or n3 nucleotides off the 5p or 3p end during bowtie alignment.
      bool: Set True if want to trim off 5p or 3p end. 
      p5: n5
      p3: n3

    min: Minimum fragment length (integer)
    max: Maximum fragment length (integer)

    fastqc: Run the fastQC quality check modules on first 2MB of data. 

      bool: Set True if want to run quality checks.

      pause: Set True if want to pause pipeline to check if QC looks OK. 
      Note: Will automatically stop if it detects adapter presence if adapter removal was set to true. 

  ### S1:

    bool: Set as True or False depending on if you want to run this module. Note that P2 will not run without this. 

  ### P2: 

    bool: Set as True or False depending on if you want to run this module.

    annotation_file: Give an annotation file for any transcripts you want to look for in OUT-OF-SPACE maps. Used only if find_out_bool is true. 
    NOTE: The annotation file is anti-joined with the filtered/processed annotation file. It should have the same primary key such that the pipeline can filter out matches that are of the correct biotype OR it should be an annotation file that has been filtered for the biotypes you are NOT looking for. 
    
    find_out_bool: Set True if want to find out-of-space maps. (Non-biotype matches)

    built_index_location: Give the locaton of the built index. Must be done prior to running pipeline. 

    plot_every_source: Plot both cluster regions and counts (relative to loci) for every potential source.

    look_for: The attribute in the annotation_file to look for. (i.e. transcript_id or biotype)

    col2_name: The value in the second column to look for. (i.e. exon or transcript)

  SUMMARY: Generate a summary report.
  NOTE: Output is an HTML file that is dependent upon the file structure. If you'd like to save the HTML report, please save it as a PDF from your browser. 
  
  delete_working: Delete the working directory. Would reccomend since it takes up a lot of space. 

### **system_options:**

    num_cores: Maximum number of cores to use or jobs to run.

    mem_free: Keep this amount of gigabytes free (i.e. 12G)

### **dir_locations:**

    working_dir: Path to working directory. Please give full path. Will create if doesn't exist.

    sample_dir: Path to sample directory. Please give full path and ensure samples are gzipped with extension, "fastq.gz".

    out_dir: Path to out directory. Please give full path. Will create if doesn't exist.

## Interested in Collaborating?

Note that development will occur in another repository soon to be uploaded. Development areas are marked by (DEVELOP FLAG) and will slowly be worked on. 

We welcome contributions to sRNAfrag! If you have a feature or improvement you'd like to add, here's how you can initiate a pull request on GitHub:

1. **Fork the Repository**: Start by forking the sRNAfrag repository to your own GitHub account.

2. **Create a New Branch**: In your forked repository, create a new branch for your feature. This helps to keep different changes separate and organized.

3. **Develop Your Feature**:
   - **New Python File**: Develop your feature in a new Python file. Ensure that the code is well-commented and adheres to the project's coding standards.
   - **README for Your Feature**: Create a README file specifically for your feature. This README should include:
     - The purpose of the feature.
     - Any new packages or dependencies required for your feature.
     - An example command or usage guide demonstrating how to use the feature.

4. **Test Your Code**: Before submitting a pull request, thoroughly test your code to ensure it works as expected and does not introduce any new bugs.

5. **Update Documentation**: If your feature adds new functionality, update the main `README.md` or the wiki to reflect this change.

6. **Create a Pull Request**:
   - Commit your changes to your branch and push them to your forked repository.
   - Navigate to the original sRNAfrag repository on GitHub.
   - Click on the 'Pull Requests' tab and then the 'New Pull Request' button.
   - Choose your fork and branch, and then submit the pull request.

7. **Code Review**:
   - Once your pull request is submitted, it will be reviewed by the repository maintainers.
   - Engage in the review process, addressing any feedback or changes requested by the maintainers.

8. **Pull Request Approval and Merge**:
   - If your pull request is approved, the repository maintainers will merge it into the main branch of the sRNAfrag repository.

9. **Celebrate Your Contribution**: Congratulations, you've successfully contributed to the sRNAfrag project!

We look forward to seeing your innovative ideas and contributions to the sRNAfrag project!

If you'd like to help with restructuring this project into a more organized framework (which is something I'd like to do, but do not have the ability to do on my own), please feel free to contact us. Reach out to Ken Nakatsu at knakats@emory.edu for further discussions and collaboration opportunities.
