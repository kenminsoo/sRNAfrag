from basics import *
from fragmentation import *
from gtf_modifiers import *
from conversion_tools import *
from gtf_groundtruth import *
import glob
import yaml
import os
import datetime

# CONFIG START #

with open("sRNA_frag_config.yaml", "r") as file:
    config_vars = yaml.safe_load(file)

# WARNINGS #
##--!!USE ABSOLUTE PATHS!!--##
##--!!FASTQ MUST BE GZ COMPRESSED!!##
# WARNINGS #

# Working directory - All intermediate files will be deleted
working_dir = config_vars["dir_locations"]["working_dir"]

# Sample directory
sample_dir = config_vars["dir_locations"]["sample_dir"]

# Out directory with all final files
out_dir = config_vars["dir_locations"]["out_dir"]

# Annotation file to fragment, should be filtered for one biotype
annotation_file = config_vars["module_options"]["P1"]["annotation_options"]["location"]

# ID Prefix
id_prefix = config_vars["module_options"]["P1"]["prefix"]

# Remove UMIs, True or False
remove_umis_bool = config_vars["module_options"]["P1"]["umi_removal"]["bool"]

# UMI Pattern
umi_pattern = config_vars["module_options"]["P1"]["umi_removal"]["regex"]

# Remove Adapters
remove_adapters_bool = config_vars["module_options"]["P1"]["adapter_removal"]["bool"]

# Adapter Settings
adapter = config_vars["module_options"]["P1"]["adapter_removal"]["sequence"]
min_length = config_vars["module_options"]["P1"]["min"]
max_length = config_vars["module_options"]["P1"]["max"]

# Quality Check?
fastqc_bool = config_vars["module_options"]["P1"]["fastqc"]["bool"]
fastqc_pause = config_vars["module_options"]["P1"]["fastqc"]["pause"]

# Alignment Option
trim_bool = config_vars["module_options"]["P1"]["trim5p3p"]["bool"]
trim_5p = config_vars["module_options"]["P1"]["trim5p3p"]["p5"]
trim_3p = config_vars["module_options"]["P1"]["trim5p3p"]["p3"]

# Keep This Amount of Memory Free (G)
mem_free = config_vars["system_options"]["mem_free"]

# Max number of cores
max_cores = config_vars["system_options"]["num_cores"]

# CONFIG END #

## Directory Setup ##

os.system("mkdir " + working_dir + "; mkdir " + out_dir)
os.system("rm " + working_dir + "/pipeline_summary.csv")
os.system("touch " + working_dir + "/pipeline_summary.csv")

## Move files to working
os.system('cd ' + sample_dir + ';\
            cp *fastq.gz ' + working_dir)

# Make the sample sheet
os.system('cd ' + working_dir + ';\
          rm samples_adapters.txt; \
          touch samples_adapters.txt; \
          for f in *.fastq.gz ; do\
                stripped=${f%%.fastq.gz};\
                echo $stripped >> samples_adapters.txt;\
                done')

# Count number of samples
with open(working_dir + "/samples_adapters.txt", "r") as samples:
    n = 0
    for line in samples:
        n += 1
    
if n == 0:
    raise ValueError("No Samples Detected. Are your files gzipped?")

# Create progress sheet
os.system('cd ' + working_dir + ';\
    echo JOB,STATUS > pipeline_summary.csv')

## FASTQC Module 1 ##
if fastqc_bool == True:
    os.system('cd ' + working_dir + ';\
    cat samples_adapters.txt | parallel "gunzip -c {}.fastq.gz | head -n 1000000 > subset_{}.fastq";\
    cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "fastqc subset_{}.fastq";\
    mkdir ' + out_dir + "/fastqc_pretrim" + ';\
    cat samples_adapters.txt | parallel "mv subset_{}_fastqc.html ' +  out_dir + "/fastqc_pretrim" +'";\
    cat samples_adapters.txt | parallel "mv subset_{}_fastqc.zip ' +  out_dir + "/fastqc_pretrim" +'";\
    echo PRETRIM FASTQC,COMPLETE >> pipeline_summary.csv')

    # Check that the command ran completely
    progress_track = pd.read_csv(working_dir + "/pipeline_summary.csv")

    status = list(progress_track.loc[progress_track['JOB'] == "PRETRIM FASTQC"]["STATUS"])[0]

    if status == "COMPLETE":
        
        # Check that the output dir has been populated
        
        fastqc_pre = os.listdir(out_dir + "/fastqc_pretrim")

        if len(fastqc_pre) > 0:
            print("OK")
        else:
            raise ValueError("FASTQC module did not run properly")

    else:
        raise ValueError("FASTQC Module Did not Run Properly")
    

    if fastqc_pause == True:
        print("Pausing to allow user to check fastQC report.")
        print("Press any key to continue.")
        print("Use ctrl-C to cancel pipeline.")
        input()

## Remove UMI from files ##
## Deduplication BAM files ##
if remove_umis_bool == True:
    os.system('cd ' + working_dir + '; \
              cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "samtools sort {}.bam -o sorted_{}.bam;"\
              cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "mv sorted_{}.bam {}.bam";\
              cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "samtools index {}.bam";\
              cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free)  + 'G --memfree ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "umi_tools dedup --method=unique --stdin={}.bam --log={}.logfile > dedup_{}.bam";\
              cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "mv dedup_{}.bam {}.bam";\
              echo UMI DEDUP,COMPLETE >> pipeline_summary.csv')

    # Ensure that bam files exist
    # There is a case where this module will fail, not producing all files
    files = glob.glob(working_dir + "/*.bam")

    if len(files) != n:
        raise ValueError("UMI Dedup failed.")

if remove_adapters_bool == True:
    os.system('cd ' + working_dir + '; \
            cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "AdapterRemoval --file1 {}.fastq.gz --adapter1 ' + adapter + ' --basename trimmed_{}.fastq.gz --gzip --trimqualities --minlength ' + str(min_length) + '";\
            cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "rm {}.fastq.gz";\
            cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "mv trimmed_{}.fastq.gz.truncated.gz trimmed_{}.fastq.gz";\
            cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "rm trimmed_{}.fastq.gz.discarded.gz";\
            cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "rm trimmed_{}.fastq.gz.settings";\
            cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "mv trimmed_{}.fastq.gz {}.fastq.gz";\
                echo ADAPTER REMOVAL,COMPLETE >> pipeline_summary.csv')

    progress_track = pd.read_csv(working_dir + "/pipeline_summary.csv")

    status = list(progress_track.loc[progress_track['JOB'] == "ADAPTER REMOVAL"]["STATUS"])[0]

    if status == "COMPLETE":
        print("OK")
    else:
        raise ValueError("Adapter Removal failed.")

## FASTQC MODULE ##
if remove_adapters_bool == True:
    if fastqc_bool == True:
        os.system('cd ' + working_dir + ';\
        cat samples_adapters.txt | parallel "gunzip -c {}.fastq.gz | head -n 1000000 > subset_{}.fastq";\
        cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "fastqc subset_{}.fastq";\
        mkdir ' + out_dir + "/fastqc_posttrim" + ';\
        cat samples_adapters.txt | parallel "mv subset_{}_fastqc.html ' +  out_dir + "/fastqc_posttrim" +'";\
        cat samples_adapters.txt | parallel "mv subset_{}_fastqc.zip ' +  out_dir + "/fastqc_posttrim" +'";\
        echo POSTTRIM FASTQC,COMPLETE >> pipeline_summary.csv')

        # Check that the command ran completely
        progress_track = pd.read_csv(working_dir + "/pipeline_summary.csv")

        status = list(progress_track.loc[progress_track['JOB'] == "POSTTRIM FASTQC"]["STATUS"])[0]

        if status == "COMPLETE":
            
            # Check that the output dir has been populated
            
            fastqc_pre = os.listdir(out_dir + "/fastqc_posttrim")

            if len(fastqc_pre) > 0:
                print("OK")
            else:
                raise ValueError("FASTQC module did not run properly")

        else:
            raise ValueError("FASTQC Module Did not Run Properly")

        # Check if adapter removal was successful
        with open(working_dir + "/samples_adapters.txt", "r") as samples:
            for line in samples:
                sample = line
                sample = sample.replace("\n", "")

                extract_zip = out_dir + "/fastqc_posttrim/subset_" + sample + "_fastqc.zip"

                os.system("cd " + out_dir + ";\
                          unzip " + extract_zip)

                adapter_status = pd.read_csv(out_dir + "/subset_" + sample + "_fastqc/summary.txt", sep = "\t")

                status = list(adapter_status.loc[adapter_status["Basic Statistics"] == "Adapter Content"].iloc[:,0])[0]

                if status == "PASS":
                    print("Adapter QC Pass")
                else:
                    raise ValueError("Adapter QC Failed")

    if fastqc_pause == True:
        print("Pausing to allow user to check fastQC report.")
        print("Press any key to continue.")
        print("Use ctrl-C to cancel pipeline.")
        input()

gtf_attribute_to_fasta(annotation_file, working_dir + "/transcriptome.fa", "sequence", "transcript_id", pipeline = True)

generate_index_bowtie(working_dir + "/transcriptome.fa", working_dir + "/transcriptome")

## Alignment Module ##
if trim_bool:
    os.system('cd ' + working_dir + '; \
                cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "bowtie -x ' +  working_dir + "/transcriptome" + ' {}.fastq.gz -k 101 --best --strata -v 0 -S {}.sam --reorder --trim5 ' + str(trim_5p) + ' --trim3 ' + str(trim_5p) + '";\
                cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "samtools view {}.sam -b -o {}.bam";\
                cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "rm {}.sam";\
                echo ALIGNMENT MODULE,COMPLETE >> pipeline_summary.csv')
else:
    os.system('cd ' + working_dir + '; \
                cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "bowtie -x ' + working_dir + "/transcriptome" + ' {}.fastq.gz -k 101 --best --strata -v 0 -S {}.sam --reorder";\
                cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "samtools view {}.sam -b -o {}.bam";\
                cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "rm {}.sam";\
                echo ALIGNMENT MODULE,COMPLETE >> pipeline_summary.csv')

progress_track = pd.read_csv(working_dir + "/pipeline_summary.csv")

status = list(progress_track.loc[progress_track['JOB'] == "ALIGNMENT MODULE"]["STATUS"])[0]

if status == "COMPLETE":
    print("OK")
else:
    raise ValueError("Alignment Module failed.")

## Deduplication BAM files ##
if remove_umis_bool == True:
    os.system('cd ' + working_dir + '; \
              cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "samtools index {}.bam";\
              cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free)  + 'G --memfree ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "umi_tools dedup --method=unique --stdin={}.bam --log={}.logfile > dedup_{}.bam";\
              cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "mv dedup_{}.bam {}.bam";\
              cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "samtools sort {}.bam -o sorted_{}.bam;"\
              cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "mv sorted_{}.bam {}.bam";\
              cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "samtools index {}.bam;"\
              echo UMI DEDUP,COMPLETE >> pipeline_summary.csv')

    # Ensure that bam files exist
    # There is a case where this module will fail, not producing all files
    files = glob.glob(working_dir + "/*.bam")

    if len(files) != n:
        raise ValueError("UMI Dedup failed.")

## Fragmentation Specific Module ##

# Create Lookup Table
create_lookup(annotation_file, "sequence", "transcript_id", min_length, max_length, working_dir + "/sRNA_frag_lookup.csv", id_prefix)
os.system('cd ' + working_dir + '; \
          echo LOOKUP TABLE,COMPLETE >> pipeline_summary.csv')
# Index

# Query for Range
os.system('cd ' + working_dir + "; \
          for f in *.bam; do\
    stripped=${f%.bam}; \
    rm ${stripped}.tsv; \
    echo sequence" + '"\t"' + "${stripped} > ${stripped}.tsv" + '; \
    done')

# Some syntax error exists so just make hidden shell script for this
os.system('cp .change.sh ' + working_dir + ";\
          cd " + working_dir + ";\
          sh .change.sh;\
          echo FRAGMENT COUNTS EXTRACTED,COMPLETE >> pipeline_summary.csv")

# Ensure Finished
progress_track = pd.read_csv(working_dir + "/pipeline_summary.csv")

status = list(progress_track.loc[progress_track['JOB'] == "FRAGMENT COUNTS EXTRACTED"]["STATUS"])[0]

if status == "COMPLETE":
    print("OK")
else:
    raise ValueError("Count extraction failed.")


# Get count info
os.system('cd ' + working_dir + '; \
          touch num_reads_bam.csv ;\
    echo "sample,n" > num_reads_bam.csv;\
    for f in *.bam; do\
    stripped=${f%.bam};\
    n_reads=$(samtools flagstat $f | head -n 1 | cut -f1 -d" ");\
    echo ${stripped},${n_reads} >> num_reads_bam.csv;\
    done;\
    echo RUN STATS,COMPLETE >> pipeline_summary.csv;\
          mv num_reads_bam.csv ' + out_dir)

# Combine files
os.system('cp .combine.py ' + working_dir + ";\
          cd " + working_dir + ";\
            python .combine.py;\
                echo COUNTS COMBINED >> pipeline_summary.csv; \
                    mv sRNA_frag_counts.csv " + out_dir + ";\
                    echo COUNTS MERGED,COMPLETE >> pipeline_summary.csv")