import os
import yaml
# Import the variables

with open("snoRNA_frag_config.yaml", "r") as file:
    config_vars = yaml.safe_load(file)

## -- Config Variables -- ##

P0_bool = config_vars['module_options']['P0']['bool']

P1_bool = config_vars['module_options']['P1']['bool']

S1_bool = config_vars['module_options']['S1']['bool']

P2_bool = config_vars['module_options']['P2']['bool']

out_dir = config_vars['dir_locations']['out_dir']

## -- Config Variables -- ##

if P0_bool == True:
    os.system("python snoRNA_fragment_P0_basic.py")

if P1_bool == False:
    os.system("python snoRNA_fragment_P1.py")

if S1_bool == False:
    os.system("cp S1_figures.R " + out_dir + ";\
              cd " + out_dir + ";\
              Rscript S1_figures.R")

if P2_bool == True:
    os.system("cp snoRNA_frag_config.yaml " + out_dir + ";\
              cp snoRNA_fragment_P2.py " + out_dir + ";\
              cp gtf_scripts/base/gtf_groundtruth.py " + out_dir + ";\
              cp gtf_scripts/base/gtf_modifiers.py " + out_dir + ";\
              cp gtf_scripts/base/conversion_tools.py " + out_dir + ";\
              cp gtf_scripts/base/alias_work.py " + out_dir + ";\
              cp gtf_scripts/base/basics.py " + out_dir + ";\
                cd " + out_dir + ";\
                    python snoRNA_fragment_P2.py")
    
# Clean up Module

if P2_bool == True:
    os.system("cd " + out_dir + ";\
            rm snoRNA_frag_config.yaml;\
            rm *.py;\
            rm *.R")