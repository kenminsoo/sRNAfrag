import os
import yaml
import time
import pandas as pd
from basics import *
from datetime import timedelta
from collections import defaultdict
# Import the variables

species = ["hs"]
biotypes = ["snrna"]

base = "sRNA_frag_config"

for organ in species:
    for biotype in biotypes:

        run_name = organ + "_" + biotype

        os.system("mv " + base + "_" + organ + "_" + biotype + ".yaml sRNA_frag_config.yaml") 

        with open("sRNA_frag_config.yaml", "r") as file:
            config_vars = yaml.safe_load(file)

        ## -- Config Variables -- ##

        P0_bool = config_vars['module_options']['P0']['bool']

        P1_bool = config_vars['module_options']['P1']['bool']

        S1_bool = config_vars['module_options']['S1']['bool']

        P2_bool = config_vars['module_options']['P2']['bool']

        out_dir = config_vars['dir_locations']['out_dir']
        working_dir = config_vars["dir_locations"]["working_dir"]

        build_bool = config_vars["module_options"]["P1"]["build_index"]["bool"]
        reference_genome = config_vars["module_options"]["P1"]["build_index"]["reference_location"]

        delete_working = config_vars['module_options']["delete_working"]

        make_summary = config_vars["module_options"]["SUMMARY"]

        ## -- Config Variables -- ##

        # Build an index if needed

        if build_bool == True:
            if reference_genome is None or reference_genome == "":
                raise ValueError("Reference location required if build index is TRUE.")
            elif os.path.isfile(reference_genome) == False:
                raise ValueError("Invalid path for build_index reference_location. (check YAML config)")

            os.system("cd " + working_dir + ";\
                    bowtie-build " + reference_genome + " fragmentation_pipeline")

        if P0_bool == True:
            os.system("python sRNA_fragment_P0_basic.py")

        overall_timing = timing()

        if P1_bool == True:
            os.system("python sRNA_fragment_P1_v2_timed.py")

        overall_timing.take_time("P1")

        if S1_bool == True:
            os.system("cp S1_figures.R " + out_dir + ";\
                    cd " + out_dir + ";\
                    Rscript S1_figures.R")
            
        overall_timing.take_time("S1")

        if P2_bool == True:
            # Copy bash script for edge case where there are over 8000 fragments detected
            os.system("mkdir " + working_dir + "/aligned_sams")
            os.system("cp .add_prefix_outspace.sh " + working_dir + "/aligned_sams")

            os.system("cp sRNA_frag_config.yaml " + out_dir + ";\
                    cp sRNA_fragment_P2_timed.py " + out_dir + ";\
                    cp gtf_groundtruth.py " + out_dir + ";\
                    cp gtf_modifiers.py " + out_dir + ";\
                    cp conversion_tools.py " + out_dir + ";\
                    cp alias_work.py " + out_dir + ";\
                    cp basics.py " + out_dir + ";\
                        cd " + out_dir + ";\
                            python sRNA_fragment_P2_timed.py")
            
        overall_timing.take_time("P2")
        overall_timing.export(run_name + "_overall")

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