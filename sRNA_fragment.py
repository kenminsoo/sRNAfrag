import os
import yaml
import time
import pandas as pd
from basics import *
from datetime import timedelta
from collections import defaultdict
# Import the variables
with open("sRNA_frag_config.yaml", "r") as file:
    config_vars = yaml.safe_load(file)

## -- Config Variables -- ##

P1_bool = config_vars['module_options']['P1']['bool']

S1_bool = config_vars['module_options']['S1']['bool']

P2_bool = config_vars['module_options']['P2']['bool']

out_dir = config_vars['dir_locations']['out_dir']
working_dir = config_vars["dir_locations"]["working_dir"]
sample_dir =config_vars["dir_locations"]["sample_dir"]

delete_working = config_vars['module_options']["delete_working"]

make_summary = config_vars["module_options"]["SUMMARY"]

# Confirm that dir exists
if os.path.isdir(out_dir) == False:
    raise ValueError("Out dir does not exist")
if os.path.isdir(working_dir)  == False:
    raise ValueError("Working dir does not exist")
if os.path.isdir(sample_dir) == False:
    raise ValueError("Sample directory does not exist")

## -- Config Variables -- ##

if P1_bool == True:
    os.system("python sRNA_fragment_P1.py")

if S1_bool == True:
    os.system("cp .S1_figures.R " + out_dir + ";\
            cd " + out_dir + ";\
            Rscript .S1_figures.R")

if P2_bool == True:
    # Copy bash script for edge case where there are over 8000 fragments detected
    os.system("mkdir " + working_dir + "/aligned_sams")
    os.system("cp .add_prefix_outspace.sh " + working_dir + "/aligned_sams")

    os.system("cp sRNA_frag_config.yaml " + out_dir + ";\
            cp sRNA_fragment_P2.py " + out_dir + ";\
            cp gtf_groundtruth.py " + out_dir + ";\
            cp gtf_modifiers.py " + out_dir + ";\
            cp conversion_tools.py " + out_dir + ";\
            cp alias_work.py " + out_dir + ";\
            cp basics.py " + out_dir + ";\
                cd " + out_dir + ";\
                    python sRNA_fragment_P2.py")

# Clean up Module
if P1_bool == True:
    os.system("cd " + out_dir + ";\
            mkdir tables;\
            mv sRNA_frag_counts.csv tables;\
            mv num_reads_bam.csv tables")
    
if S1_bool == True:
    os.system("cd " + out_dir + ";\
            mkdir figures;\
            mv S1_std_loci.jpeg figures;\
            mv S1_motifs_counts.jpeg figures;\
            mv S1_freq_motifs.jpeg figures;\
            mv S1_five-mer_5p.png figures;\
            mv S1_five-mer_3p.png figures;\
            mv S1_correction_density.jpeg figures;\
            mv filtered_corrected_counts.csv tables;\
            mv filtered_counts.csv tables")


if P2_bool == True:
    os.system("cd " + out_dir + ";\
            mkdir int_files;\
            mkdir sequences;\
            mv sRNA_frag_config.yaml int_files;\
            mv *.py int_files;\
            mv *.R int_files;\
            mv merged_counts.csv tables;\
            mv annotated_ref_table.csv tables;\
            mv filtered_sequences.fa sequences;\
            mv ref_table.csv tables;\
            mv P2_cluster_start_end_disagreement.jpeg figures;\
            mv alternate_peaks.csv tables;\
            mv cluster_peak_relationship_table.csv tables;\
            mv P2_Filter_Passing_Hist.jpeg figures")
    

if delete_working == True:
    os.system("rm -rf " + working_dir)

if make_summary == True:
    os.system("cp summary.md " + out_dir + ";\
            cd " + out_dir + ";\
            python -m markdown summary.md > summary.html")