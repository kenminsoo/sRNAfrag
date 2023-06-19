from gtf_scripts.base.basics import *
from fragmentation import *
import os
from gtf_scripts.base.gtf_modifiers import *
from gtf_scripts.base.conversion_tools import *
from gtf_scripts.base.gtf_groundtruth import *
import glob
import time

# CONFIG START #

# WARNINGS #
##--!!USE ABSOLUTE PATHS!!--##
##--!!FASTQ MUST BE GZ COMPRESSED!!##
# WARNINGS #

# Working directory - All intermediate files will be deleted
working_dir = "/Volumes/Extreme_SSD/rerun_cell_lines/working"

# Sample directory
sample_dir = "/Volumes/Extreme_SSD/rerun_cell_lines/samples"

# Reference genome location
build_bool = False
reference_genome = ""

# Alignment index name or location
index_name = "/Users/kenminsoo/Desktop/Projects/JABSOM/Nakatsu/snoRNA_Breast/smallRNA_structure/output/adaptertrim/hg38"

# Out directory with all final files
out_dir = "/Volumes/Extreme_SSD/rerun_cell_lines/out"

# Annotation file to fragment, should be filtered for one biotype
annotation_file = "/Users/kenminsoo/Desktop/Projects/JABSOM/Nakatsu/projects/tools/snoRNA_std_box_v2.gtf"
add_sequences_bool = False

# Remove UMIs, True or False
remove_umis_bool = True

# UMI Pattern
umi_pattern = ".+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12}).+"

# Remove Adapters
remove_adapters_bool = True

# Adapter Settings
adapter = "AGATCGGAAGAG"
min_length = 15
max_length = 45

# Quality Check?
fastqc_bool = True

# Keep This Amount of Memory Free (G)
mem_free = 12

# Max number of cores
max_cores = 3

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

timing = []
st = time.time()

## FASTQC Module 1 ##
if fastqc_bool == True:
    os.system('cd ' + working_dir + ';\
    cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "fastqc {}.fastq.gz";\
    mkdir ' + out_dir + "/fastqc_pretrim" + ';\
    cat samples_adapters.txt | parallel "mv {}_fastqc.html ' +  out_dir + "/fastqc_pretrim" +'";\
    cat samples_adapters.txt | parallel "mv {}_fastqc.zip ' +  out_dir + "/fastqc_pretrim" +'";\
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

et = time.time()

timing.append(["FASTQC_Pre", str(et - st)])

## Remove UMI from files ##

st = time.time()

if remove_umis_bool == True:
    os.system('cd ' + working_dir + ';\
    cat samples_adapters.txt | parallel -I ,,,  --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "umi_tools extract --extract-method=regex --bc-pattern=' + "'"+ umi_pattern  + "'" + ' -I ,,,.fastq.gz -S processed.,,,.fastq.gz";\
        cd ' + working_dir + "; \
            for f in processed* ; do \
                new_file=${f#processed.}; \
                mv $f $new_file;\
                done; \
        echo UMI REMOVAL,COMPLETE >> pipeline_summary.csv")

    progress_track = pd.read_csv(working_dir + "/pipeline_summary.csv")

    status = list(progress_track.loc[progress_track['JOB'] == "UMI REMOVAL"]["STATUS"])[0]

    if status == "COMPLETE":
        print("OK")
    else:
        raise ValueError("UMI extraction failed.")


et = time.time()

timing.append(["UMI_Extraction", str(et - st)])
## Remove Adapters ##

st = time.time()

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

et = time.time()

timing.append(["Adapter Removal", str(et - st)])

## FASTQC MODULE ##

st = time.time()

if remove_adapters_bool == True:
    if fastqc_bool == True:
        os.system('cd ' + working_dir + ';\
        cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "fastqc {}.fastq.gz";\
        mkdir ' + out_dir + "/fastqc_posttrim" + ';\
        cat samples_adapters.txt | parallel "mv {}_fastqc.html ' +  out_dir + "/fastqc_posttrim" +'";\
        cat samples_adapters.txt | parallel "mv {}_fastqc.zip ' +  out_dir + "/fastqc_posttrim" +'";\
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

                extract_zip = out_dir + "/fastqc_posttrim/" + sample + "_fastqc.zip"

                os.system("cd " + out_dir + "/fastqc_posttrim;\
                unzip -o " + extract_zip)

                adapter_status = pd.read_csv(out_dir + "/fastqc_posttrim/" + sample + "_fastqc/summary.txt", sep = "\t")

                status = list(adapter_status.loc[adapter_status["Basic Statistics"] == "Adapter Content"]["PASS"])[0]

                if status == "PASS":
                    print("Adapter QC Pass")
                else:
                    raise ValueError("Adapter QC Failed")

et = time.time()

timing.append(["FASTQC_Post", str(et - st)])

## Alignment Module ##

st = time.time()

os.system('cd ' + working_dir + '; \
            cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "bowtie -x ' + index_name + ' {}.fastq.gz -k 101 --best --strata -v 0 -S {}.sam --reorder";\
            cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "samtools view {}.sam -b -o {}.bam";\
            cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "samtools sort {}.bam -o sorted_{}.bam";\
            cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "mv sorted_{}.bam {}.bam";\
            cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "rm {}.sam";\
            echo ALIGNMENT MODULE,COMPLETE >> pipeline_summary.csv')

progress_track = pd.read_csv(working_dir + "/pipeline_summary.csv")

status = list(progress_track.loc[progress_track['JOB'] == "ALIGNMENT MODULE"]["STATUS"])[0]

if status == "COMPLETE":
    print("OK")
else:
    raise ValueError("Alignment Module failed.")

et = time.time()

timing.append(["Alignment", str(et - st)])

st = time.time()

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

et = time.time()

timing.append(["UMI Deduplication", str(et - st)])

## Fragmentation Specific Module ##

# Add sequence if needed !!! ENSURE SEQUENCE IS LABELED AS SEQUENCE
# !!! ENSURE PRIMARY KEY IS CALLED TRANSCRIPT ID

st = time.time()

gtf_to_bed(annotation_file, working_dir + "/filtering.bed")
os.system('cd ' + working_dir + '; \
          echo FILTERING BED,COMPLETE >> pipeline_summary.csv')

et = time.time()

timing.append(["Bed Creation", str(et - st)])

st = time.time()

# Create Lookup Table
create_lookup(annotation_file, "sequence", "new_tid", min_length, max_length, working_dir + "/snoRNA_frag_lookup.csv")
os.system('cd ' + working_dir + '; \
          echo LOOKUP TABLE,COMPLETE >> pipeline_summary.csv')

et = time.time()

timing.append(["Lookup Creation", str(et - st)])

st = time.time()

# Index
os.system('cd ' + working_dir + '; \
          cat samples_adapters.txt | parallel --memsuspend ' + str(mem_free) + 'G -j ' + str(max_cores) + ' --progress "samtools index {}.bam";\
          echo BAM FILES INDEXED,COMPLETE >> pipeline_summary.csv')

et = time.time()

timing.append(["BAM indexing", str(et - st)])

st = time.time()

# Query for Range
os.system('cd ' + working_dir + "; \
          for f in *.bam; do\
    stripped=${f%.bam}; \
    rm ${stripped}.csv; \
    echo sequence,${stripped} > ${stripped}.csv; \
    samtools view -M -L " + working_dir + "/filtering.bed" + " $f | cut -f 10,2,15 | grep -P 'XM:i:[0-9]{1,2}$' | cut -f 1,2 | awk -v OFS='_' '{print $1, $2}' | sort | uniq -c | awk -v OFS=',' '{print $2, $1}' >> ${stripped}.csv; \
done;\
          echo FRAGMENT COUNTS EXTRACTED,COMPLETE >> pipeline_summary.csv")

# Ensure Finished
progress_track = pd.read_csv(working_dir + "/pipeline_summary.csv")

status = list(progress_track.loc[progress_track['JOB'] == "FRAGMENT COUNTS EXTRACTED"]["STATUS"])[0]

if status == "COMPLETE":
    print("OK")
else:
    raise ValueError("Count extraction failed.")

et = time.time()

timing.append(["Fragment Extraction", str(et - st)])

st = time.time()

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

et = time.time()

timing.append(["Bam Stats", str(et - st)])

st = time.time()

# Combine files
os.system('cp combine.py ' + working_dir + ";\
          cd " + working_dir + ";\
            python combine.py;\
                echo COUNTS COMBINED >> pipeline_summary.csv; \
                    mv snoRNA_frag_counts.csv " + out_dir + ";\
                    echo COUNTS MERGED,COMPLETE >> pipeline_summary.csv")

et = time.time()

timing.append(["Combine Counts", str(et - st)])

timing_df = pd.DataFrame(timing, columns = ["Job", "Time"])

timing_df.to_csv("time.csv", index = False)