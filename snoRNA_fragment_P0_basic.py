from gtf_scripts.base.alias_work import *
from gtf_scripts.base.conversion_tools import *
from gtf_scripts.base.gtf_descriptors import *
from gtf_scripts.base.gtf_generation import *
from gtf_scripts.base.gtf_modifiers import *
import os

print("Do you need to convert GFF file to GTF?")
print("Y/N")

gff_bool = input()

working_file = []

# Change GFF format to GTF if needed
if gff_bool == "N" or "n":
    print("Moving on.")
else:
    print("Input path to GFF file to convert to GTF.")
    path_gff = input()
    if os.path.isfile(path_gff) == False:
        raise ValueError("Path is incorrect.")
    print("Input output file path & name. i.e. /home/new.gtf")
    output_gtf = input()
    working_file.append(output_gtf)

    gff3_to_gtf(path_gff, output_gtf)

print("Do you need to convert primary keys in your GTF to transcript_id?")
print("(Probably yes)")
print("Y/N")

convert_keys_bool = input()

# Convert primary key to transcript_id
if convert_keys_bool == "Y" or "y":
    conv_dict = my_dictionary()
    print("What is the primary key in your gtf file?")

    prim_key = input()

    conv_dict.add(prim_key, "transcript_id")

    if len(working_file) == 0:
        print("Input path to GTF to standarize.")
        path_gtf = input()

        if os.path.isfile(path_gtf) == False:
            raise ValueError("Path is incorrect.")
    
    else:
        path_gtf = working_file[-1]

    print("Input output file path & name. i.e. /home/new.gtf")
    output_gtf = input()
    working_file.append(output_gtf)

    standardize_attributes(path_gtf, output_gtf, conv_dict)

print("Do you need to add sequences to your GTF File?")
print("Y/N")

add_bool = input()

if add_bool == "Y" or "y":
    if len(working_file) == 0:
        print("Input path to GTF to add sequences.")
        path_gtf = input()

        if os.path.isfile(path_gtf) == False:
            raise ValueError("Path is incorrect.")

    else:
        path_gtf = working_file[-1]

    print("Input output file path & name. i.e. /home/new.gtf")
    output_gtf = input()
    working_file.append(output_gtf)

    print("Input path to a reference genome.")
    ref_genome_i = input()

    add_sequence_gtf(path_gtf, ref_genome_i, output_gtf, "sequence")