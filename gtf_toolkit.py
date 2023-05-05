import numpy as np
import pandas as pd
from pybedtools import BedTool
import os
import time
from basics import *
#from pybedtools import BedTool #we need a new inst

# List of functions in this file and descriptions

# Create dictionary
# This is a class to create a dictionary
# And allows for easy adding

## -- Basic Functions -- ##
""" class my_dictionary(dict):
 
  # __init__ function
  def __init__(self):
    self = dict()
 
  # Function to add key:value
  def add(self, key, value):
    self[key] = value

# Calculate the hamming distance between two sequences
def hamming(seq1, seq2):
    i = 0
    counter = 0
    for letter in seq1:

        if letter != seq2[i]:
            counter += 1
        i += 1

    return counter

# Translate RNA if ever needed
def rna_trans(transcript):
    
    RNA = ""
    
    for letter in transcript:
        if letter == "T":
            RNA = RNA + "U"
        else:
            RNA = RNA + letter
    
    return RNA """

## -- Basic Functions -- ##

# turn gff3 into gtf file
""" def gff3_to_gtf(gff3_file, output_name):
    with open(gff3_file, "r") as gff3, open(output_name, "w") as new:
        for line in gff3:

            temp_line = line.strip()
            temp_line = temp_line.replace("=", ' "')
            temp_line = temp_line.replace(";", '"; ')
            temp_line = temp_line + '"'
            new.write(temp_line + "\n") """

# add a feature that extracts the sequence from the genome
# Note that lower case letters represent masked regions
# It only takes a gtf file 
# BedTools recognizes this
""" def add_sequence_legacy(gtf_file, ref_genome, output_name):
    start_time = time.time()
    
    fasta = ref_genome
    # extract chromosome and location
    with open(gtf_file, "r") as gtf, open(output_name, "w") as new:
        for line in gtf:
            transcript_info = BedTool(line, from_string=True)

            # Ensure that strandedness is on for sequence extract

            transcript_seq = transcript_info.sequence(fi = fasta, s = True)
            sequence = open(transcript_seq.seqfn).read()

            extracted_sequence = sequence.strip().split(sep = "\n")[1]

            temp_line = line.strip()

            new.write(temp_line + '; bed_sequence "' + extracted_sequence + '"\n')

    print("--- %s seconds ---" % (time.time() - start_time)) """


""" # faster version of the above - less safe and assumes that 
# lines will correspond to each other i.e. 
# 1 in gtf = 1 i in the output
# should be true
# NOTE Bed = gtf, it can accept both i think
def add_sequence_fast(bed_file, ref_genome):
    
    fasta = ref_genome

    with open(bed_file, "r") as bed:

        bed_file = bed.read()
        transcript_info = BedTool(bed_file, from_string=True)
        transcript_seq = transcript_info.sequence(fi = fasta, s = True)
        sequence = open(transcript_seq.seqfn).read()

        seq_list = sequence.strip().split(sep = "\n")

    return seq_list

def add_sequence_gtf(gtf_file, ref_genome, output_name, attribute_name = "bed_sequence"):
    
    seq_list = add_sequence_fast(gtf_file, ref_genome)

    with open(gtf_file, "r") as bed, open(output_name, "w") as new:
        
        iter = 1
        
        for line in bed:
            extracted_sequence = seq_list[iter]

            temp_line = line.strip()

            new.write(temp_line + '; ' + attribute_name + ' "' + extracted_sequence + '"\n')

            iter += 2
 """
# we define a new file type, bedseq format
""" def add_sequence_bed(bed_file, ref_genome, output_name):

    seq_list = add_sequence_fast(bed_file, ref_genome)

    with open(bed_file, "r") as gtf, open(output_name, "w") as new:

        iter = 1

        for line in gtf:
            extracted_sequence = seq_list[iter]

            temp_line = line.strip()

            new.write(temp_line + "\t" + extracted_sequence + "\n")

            iter += 2 """

""" # Input gtf file, if sequence exists, this will compare the extracted one (with bedtools) to see if any errors exist in the annotation
def compare_sequence(gtf_file, output_name, sequence_feature1, sequence_feature2):
    # store data in here

    convert_df = []

    with open(gtf_file, "r") as gtf:
        for line in gtf:
            temp_line = line
            # split
            split_line = line.strip().split(sep = "\t")
            # get only sequence info
            attributes = split_line[-1].split(sep = ";")
            # strip
            attributes = list(map(str.strip, attributes))
            # combine
            attributes_split = [item.split(sep = " ") for item in attributes]
            # combine
            attributes_split = sum(attributes_split, [])
            
            # this is done because combining annotation files often results in
            # different index numbers for features

            # find index of first feature and pull sequence
            index1 = attributes_split.index(sequence_feature1)
            seq1 = attributes_split[index1 + 1].replace('"', "").upper()

            # find index of second feature and pull sequence
            index2 = attributes_split.index(sequence_feature2)
            seq2 = attributes_split[index2 + 1].replace('"', "").upper()

            # check if same length
            if len(seq1) != len(seq2):
                dist = "error"
            else:
                dist = hamming(seq1, seq2)
            temp_data = [temp_line.strip(), dist, seq1, seq2]
            convert_df.append(temp_data)
    df_for_output = pd.DataFrame(convert_df)
    df_for_output.to_csv(output_name) """

# Next, we want to turn a tsv file to a gtf
# We need to provide a list that refers to elements of a gtf
# In the tsv


""" extract_list = [12, "snoDB2023", "transcript_ID", 13, 14, ".", 15, ".",]

def tsv_to_gtf(tsv, out_name, extract_list):
    with open(tsv, "r") as tsv_file, open(out_name, "w") as new_file:
        
        iter = 0 
        for line in tsv_file:
            # skip first line
            if iter == 0:
                iter += 1
                continue

            features = line.split(sep = "\t")
            # an ugly way to merge everything
            new_file.write("\t".join([features[extract_list[0]], extract_list[1], 
            extract_list[2], features[extract_list[3]], 
            features[extract_list[4]], extract_list[5], features[extract_list[6]], 
            extract_list[7], "; ".join(['biotype "snoRNA"', ("sequence " + '"' + str(features[22]) + '"'),
            ("transcript_id " + '"' + str(features[0]) + '/' + str(features[16]) + '"')])]) + "\n") """

# Combines reference genomes with
# Different chr aliases
""" def ref_combine(fa1_name, fa2_name, out_name, reference):
    
    reference_nomen = pd.read_csv(reference, sep = "\t")
    reference_nomen["Source"] = ""
    reference_nomen["Exist"] = False

    add_seq = False
    number_columns = len(reference_nomen.columns)
    
    with open(fa1_name, "r") as fa1, open(fa2_name, "r") as fa2, open(out_name, "w") as new:
        for line in fa1:
            # Identify if marks chromosome
            if line[0:1] == ">":
                # See what column it is in

                splt = line[1:].split()[0]

                for column in reference_nomen:
                    iters = 0
                    if splt in list(reference_nomen[column]):
                        # Find where chr is in the dataframe
                        index_num = reference_nomen.index[reference_nomen[column] == splt].to_list()[0]
                        
                        # Mark as existing
                        if reference_nomen.iloc[index_num, -1] == False:
                            # Set it to true if it was false
                            reference_nomen.iloc[index_num, -1] = True
                            reference_nomen.iloc[index_num, -2] = fa1_name
                            add_seq = True
                            new.write(">" + str(reference_nomen.iloc[index_num, 0] + "\n"))
                            break
                        
                        else:
                            # if it is true we do not add sequence
                            add_seq = False
                            break

                    else:
                        print("not in this " + column)
                        iters += 1

                        # If it's never found, ignore
                        # Perhaps we can skip and add to summary report after?
                        if iters == number_columns:
                            raise ValueError("Chromosome Not found")
            else:
                if add_seq == True:
                    new.write(line)
                else:
                    continue

        for line in fa2:
            # Identify if marks chromosome
            if line[0:1] == ">":
                # See what column it is in

                splt = line[1:].split()[0]

                for column in reference_nomen:
                    iters = 0
                    if splt in list(reference_nomen[column]):
                        # Find where chr is in the dataframe
                        index_num = reference_nomen.index[reference_nomen[column] == splt].to_list()[0]
                        
                        # Mark as existing
                        if reference_nomen.iloc[index_num, -1] == False:
                            # Set it to true if it was false
                            reference_nomen.iloc[index_num, -1] = True
                            reference_nomen.iloc[index_num, -2] = fa2_name
                            add_seq = True
                            new.write(">" + str(reference_nomen.iloc[index_num, 0] + "\n"))
                            break
                        
                        else:
                            # if it is true we do not add sequence
                            add_seq = False
                            break

                    else:
                        print("not in this " + column)
                        iters += 1

                        # If it's never found, ignore
                        # Perhaps we can skip and add to summary report after?
                        if iters == number_columns:
                            raise ValueError("Chromosome Not found")

            else:
                if add_seq == True:
                    new.write(line)
                else:
                    continue

    # Extract the informational file
    reference_nomen.to_csv("merge_info.csv") """

""" # Standardize gtf chr names
# When different chr aliases are used
def gtf_naming_stan(gtf_file, reference, out_name):
    reference_nomen = pd.read_csv(reference, sep = "\t")

    num_columns = len(reference_nomen.columns)

    with open(gtf_file, "r") as gtf, open(out_name, "w") as new:
        
        for line in gtf:
            
            modify_line = line.split(sep = "\t")

            chromosome_name = modify_line[0]

            chr_not_found = True

            while chr_not_found == True:
                for column in reference_nomen:
                    if chromosome_name in list(reference_nomen[column]):
                        index_num = reference_nomen.index[reference_nomen[column] == chromosome_name].to_list()[0]
                        chr_not_found = False
                        break
                    else:
                        continue

                if chr_not_found == False:
                    break
                elif chr_not_found == True:
                    raise ValueError(chromosome_name + " not in alaias!")

            modify_line[0] = reference_nomen.iloc[index_num, 0]
            new.write("\t".join(modify_line)) """

# commented out 
# This is a quick script that users can use if their 
# chromosome names have been altered and needs 
# modification before replacement
""" def gtf_naming_stan_temp(gtf_file, reference, out_name):
    reference_nomen = pd.read_csv(reference, sep = "\t")

    with open(gtf_file, "r") as gtf, open(out_name, "w") as new:
        
        for line in gtf:
            
            modify_line = line.split(sep = "\t")

            chromosome_name = modify_line[0]

            if "CHR_" in chromosome_name:
                split = chromosome_name.split(sep = "_")
                joining = split[1:]

                chromosome_name = "_".join(joining)

            chr_not_found = True

            while chr_not_found == True:
                for column in reference_nomen:
                    if chromosome_name in list(reference_nomen[column]):
                        index_num = reference_nomen.index[reference_nomen[column] == chromosome_name].to_list()[0]
                        chr_not_found = False
                        break
                    else:
                        continue

                if chr_not_found == False:
                    break
                elif chr_not_found == True:
                    raise ValueError(chromosome_name + " not in alaias!")

            modify_line[0] = reference_nomen.iloc[index_num, 0]
            new.write("\t".join(modify_line)) """


#gtf_naming_stan_temp("testing_env/final/homo_sapies.GRCh38.gtf", "testing_env/hg38.chromAlias.txt", "RC_Fixed_Chr.gtf")

# Change bed to gtf loci format
""" def zero_to_one(gtf_file, out_name):
    with open(gtf_file, "r") as gtf, open(out_name, "w") as new:
        for line in gtf:
            modify_line = line.split(sep="\t")
            modify_line[3] = int(modify_line[3])
            modify_line[3] = modify_line[3] + 1
            modify_line[3] = str(modify_line[3])

            new.write("\t".join(modify_line))

# Change from gtf to bed loci format
def one_to_zero(gtf_file, out_name):
    with open(gtf_file, "r") as gtf, open(out_name, "w") as new:
        for line in gtf:
            modify_line = line.split(sep="\t")
            modify_line[3] = int(modify_line[3])
            modify_line[3] = modify_line[3] - 1
            modify_line[3] = str(modify_line[3])

            new.write("\t".join(modify_line)) """

""" # Modify the exon value to something else
def gtf_change_middle(gtf_file, out_name, change_value):
    with open(gtf_file, "r") as gtf, open(out_name, "w") as new:
        for line in gtf:
            modify_line = line.split(sep="\t")
            modify_line[2] = str(change_value)

            new.write("\t".join(modify_line)) """

# Filter by attribute
# Will take out whatever attribute you put in
""" def filter_gtf(input_gtf, output_name, filter_by, value, skip = False, num = 1):
    with open(input_gtf, "r") as annotation, open(output_name, "w") as new:
        # skip n number of lines
        iter = 1
        for line in annotation:
            if skip == True: 
                if iter <= num:
                    iter += 1
                    continue
            # split
            split_line = line.strip().split(sep = "\t")
            # get only sequence info
            attributes = split_line[-1].split(sep = ";")
            # strip
            attributes = list(map(str.strip, attributes))
            # combine
            attributes_split = [item.split(sep = " ") for item in attributes]
            # combine
            attributes_split = sum(attributes_split, [])
            # this is done because combining annotation files often results in
            # different index numbers for features
            # find index of first feature and pull sequence
            index1 = attributes_split.index(filter_by)
            type_feature = attributes_split[index1 + 1]

            if type_feature == value:
                continue
            else: 
                new.write(line) """

# Will filter by a column
# 0 index
""" def select_column(input_gtf, output_name, col_number, value):
    with open(input_gtf, "r") as gtf, open(output_name, "w") as new:
        for line in gtf:
            sep = line.split(sep = "\t")

            if sep[col_number] == value:
                new.write(line)

            else:
                continue

def filter_column(input_gtf, output_name, col_number, value):
    with open(input_gtf, "r") as gtf, open(output_name, "w") as new:
        for line in gtf:
            sep = line.split(sep = "\t")

            if sep[col_number] == value:
                continue

            else:
                new.write(line) """

# Count the number of times a certain field appears
""" def countby_field(input_gtf, output_name, field_index):
    
    return_data = my_dictionary()
    
    with open(input_gtf, "r") as gtf, open(output_name, "w") as new:
        for line in gtf:
            split_line = line.strip().split(sep = "\t")
            value = split_line[field_index]

            if value not in return_data:
                return_data.add(value, 0)
            
            return_data[value] += 1
        
        new.write(str(return_data)) """

# Count the number of times a certain attribute appears 
""" def countby_attribute(input_gtf, output_name, countby_value, skip = False, num = 1):
    
    return_data = my_dictionary()

    with open(input_gtf, "r") as gtf, open(output_name, "w") as new:
        
        iter = 1
        for line in gtf:
            if skip == True: 
                if iter <= num:
                    iter += 1
                    continue
            split_line = line.strip().split(sep = "\t")
            # get only sequence info
            attributes = split_line[-1].split(sep = ";")
            # strip
            attributes = list(map(str.strip, attributes))
            # combine
            attributes_split = [item.split(sep = " ") for item in attributes]
            # combine
            attributes_split = sum(attributes_split, [])
            # this is done because combining annotation files often results in
            # different index numbers for features
            # find index of first feature and pull sequence
            index1 = attributes_split.index(countby_value)
            type_feature = attributes_split[index1 + 1]

            if type_feature not in return_data:
                return_data.add(type_feature, 0)

            return_data[type_feature] += 1

        new.write(str(return_data)) """


# Merges overlapping sequences
# Able to be offset by a certain amount to extend to regions
# And clusters
""" def merge_overlaps(gtf_file, output_tsv, output_gtf, offset):
    # Assume gtf format
    column_names = ["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
    gtf_data = pd.read_csv(gtf_file, sep = "\t", names = column_names)
    
    gtf_data = gtf_data.sort_values("start").groupby(["chr", "strand"])

    keys_df = gtf_data.groups.keys()

    # loop through all chromosomes and strands
    stored_data = my_dictionary()
    stored_data.add("overall", 0)

    bigger_data = pd.DataFrame(columns=column_names)

    for key in keys_df:
        stored_data.add("".join(key), 0)

        temp_data = gtf_data.get_group(key)
        temp_data = pd.DataFrame(temp_data)

        n = -1
        n_name = 1
        any_trues = True


        # Calculates if any overlap exist
        # And how deep the overlap goes down
        # If sum overlaps = 1
        # It means the next transcript overlaps with it
        # Can it deal with nested overlaps? 
        # ----   ----  <--- the expression that encapsulates
        #    -----          my face while working with annotation files
        # We want to smash together while maintaining annotation info

        while any_trues == True:

            temp_data["start_next" + str(n_name)] = temp_data["start"].shift(n)

            temp_data["difference_" + str(n_name)] = temp_data["start_next" + str(n_name)] - temp_data["end"] - offset

            temp_data["overlap_" + str(n_name)] = temp_data["difference_" + str(n_name)] < 0

            if True in list(temp_data["overlap_" + str(n_name)]):
                n += -1
                n_name += 1
                continue
            else:
                any_trues = False

        temp_data["Sum_Overlaps"] = temp_data[list(temp_data.filter(regex='overlap'))].sum(axis=1)
        
        temp_data["overlapping?"] = temp_data["Sum_Overlaps"] > 0
        
        smaller_data = temp_data.iloc[:, 0:9] 

        smaller_data["Sum_Overlaps"] = temp_data["Sum_Overlaps"]

        smaller_data["overlapping"] = temp_data["overlapping?"]

        bigger_data = pd.concat([bigger_data, smaller_data])

    bigger_data.to_csv(output_tsv, sep = "\t", header = False, index = False, quotechar= "~")

    # Now merge

    with open(output_tsv, "r") as gtf, open(output_gtf, "w") as new:
        
        skip_n = 0

        end_indicies = []
        start_indicies = []
        combined_attribute = my_dictionary()

        for line in gtf:
            # Separate data
            sep_data = line.split(sep = "\t")
            
            # skip the number of lines obtained from below
            # We use this section to construct a new merged line
            if skip_n != 0:
                new_chr = sep_data[0]
                new_strand = sep_data[6]

                # Initiate finishing procedure if chr or strand dont match
                # I think this is wrong--we need to fix this
                # This needs to add the previous FIRST
                # Then restart and if it needs to start skipping
                # Then we add the root chrs, strands, etc.
                if root_chr != new_chr or root_strand != new_strand:
                    if sep_data[-1] == "False\n":
                        new.write("\t".join(sep_data[0:9]) + "\n")
                    else:
                        # Establish the chromosome and strand
                        root_chr = sep_data[0]
                        root_strand = sep_data[6]

                        start_indicies.append(sep_data[3])
                        end_indicies.append(sep_data[4])

                        # get only attribute info
                        attributes = sep_data[8].split(sep = ";")
                        # strip
                        attributes = list(map(str.strip, attributes))
                        # combine
                        attributes_split = [item.split(sep = " ") for item in attributes]
                        # combine
                        attributes_split = sum(attributes_split, [])

                        # even will be the feature name
                        # odd will be the feature value
                        i = 0
                        for entry in attributes_split:
                            if i % 2 == 0:
                                # it is even and thus a feature name
                                combined_attribute.add(entry.strip(), '"')
                                i += 1
                            else:
                                combined_attribute[attributes_split[i - 1].strip()] = combined_attribute[attributes_split[i - 1]] + entry.replace('"', "")

                                i += 1     

                    chosen_start = min(start_indicies)
                    chosen_end = max(end_indicies)

                    # we have stored our attribute info in a dictionary
                    # Keys :)
                    add_atribute = []

                    for key in combined_attribute:
                        add_atribute.append(key + " ")
                        add_atribute.append(combined_attribute[key] + '"; ')

                    # can't forget the wonderful newline!
                    add_atribute.append("\n")

                    final_entry = [sep_data[0], "merged", "exon", chosen_start, chosen_end, sep_data[5], 
                        sep_data[6], sep_data[7], "".join(add_atribute)]

                    new.write("\t".join(final_entry))

                    # Reset
                    end_indicies = []
                    start_indicies = []
                    combined_attribute = my_dictionary()

                    # this breaks out of loop and starts the skipping process
                    # if needed

                    skip_n = int(float(sep_data[-2]))
                    skip_n_end = int(sep_data[4])

                    continue

                skip_n += -1

                # Add indicies
                start_indicies.append(sep_data[3])
                end_indicies.append(sep_data[4])
                # If the next entry end point is over the skip
                # root, then we add
                current_line_end = int(sep_data[4])
                
                if int(current_line_end) > int(skip_n_end):
                    skip_n = int(float(sep_data[-2]))
                    skip_n_end = sep_data[4]

                # get only attribute info
                attributes = sep_data[8].split(sep = ";")
                # strip
                attributes = list(map(str.strip, attributes))
                # combine
                attributes_split = [item.split(sep = " ") for item in attributes]
                # combine
                attributes_split = sum(attributes_split, [])

                # even will be the feature name
                # odd will be the feature value
                i = 0
                for entry in attributes_split:
                    if i % 2 == 0:
                        # it is even and thus a feature name
                        # make sure to add only if it is not in there already
                        if entry.strip() in combined_attribute:
                            i += 1
                            continue
                        else:
                            combined_attribute.add(entry.strip(), '"')
                            i += 1
                    else:
                        combined_attribute[attributes_split[i - 1].strip()] = combined_attribute[attributes_split[i - 1].strip()] + "=!=" + entry.replace('"', "")
                        i += 1


                # Final action => Add to annotation
                if skip_n == 0:
                    chosen_start = min(start_indicies)
                    chosen_end = max(end_indicies)

                    # we have stored our attribute info in a dictionary
                    # Keys :)
                    add_atribute = []

                    for key in combined_attribute:
                        add_atribute.append(key + " ")
                        add_atribute.append(combined_attribute[key] + '"; ')

                    # can't forget the wonderful newline!
                    add_atribute.append("\n")

                    final_entry = [sep_data[0], "merged", "exon", chosen_start, chosen_end, sep_data[5], 
                        sep_data[6], sep_data[7], "".join(add_atribute)]

                    new.write("\t".join(final_entry))

                    # Reset
                    end_indicies = []
                    start_indicies = []
                    combined_attribute = my_dictionary()
                
                continue

            elif sep_data[-1] == "False\n":
                new.write("\t".join(sep_data[0:9]) + "\n")
            else:
                # Establish the chromosome and strand
                root_chr = sep_data[0]
                root_strand = sep_data[6]

                skip_n = int(float(sep_data[-2]))
                skip_n_end = int(sep_data[4])

                start_indicies.append(sep_data[3])
                end_indicies.append(sep_data[4])

                # get only attribute info
                attributes = sep_data[8].split(sep = ";")
                # strip
                attributes = list(map(str.strip, attributes))
                # combine
                attributes_split = [item.split(sep = " ") for item in attributes]
                # combine
                attributes_split = sum(attributes_split, [])

                # even will be the feature name
                # odd will be the feature value
                i = 0
                for entry in attributes_split:
                    if i % 2 == 0:
                        # it is even and thus a feature name
                        combined_attribute.add(entry.strip(), '"')
                        i += 1
                    else:
                        combined_attribute[attributes_split[i - 1].strip()] = combined_attribute[attributes_split[i - 1]] + entry.replace('"', "")

                        i += 1

    return [stored_data] """

""" # This requires HISAT2 to be installed and a index to be built
# Obtains transcript-copy-id which marks different locations in genome
def align_hisat(gtf_file, output_gtf, index_name):
    with open(gtf_file, "r") as gtf, open("extracted_seq.fa", "w") as new:
        
        added_transcripts = my_dictionary()
        # Create a fasta file with sequences
        # and transcript names
        for line in gtf:
                sep_data = line.split(sep = "\t")
                # get only attribute info
                attributes = sep_data[8].split(sep = ";")
                # strip
                attributes = list(map(str.strip, attributes))
                # combine
                attributes_split = [item.split(sep = " ") for item in attributes]
                # combine
                attributes_split = sum(attributes_split, [])

                transcript_id_index = attributes_split.index("transcript_id")
                transcript_id = attributes_split[transcript_id_index + 1]

                sequence_index = attributes_split.index("sequence")
                sequence = attributes_split[sequence_index + 1]
                
                transcript_id_clean = transcript_id.replace('"', "")

                sequence_clean = sequence.replace('"', "")

                biotype_index = attributes_split.index("biotype")
                biotype = attributes_split[biotype_index + 1]

                biotype_clean = biotype.replace('"', "")

                database_id = sep_data[1]

                if transcript_id_clean not in added_transcripts:
                    added_transcripts.add(transcript_id_clean, sequence_clean)

                    new.write(">" + transcript_id_clean + "--" + biotype_clean + "--" + database_id + "\n" + sequence_clean + "\n")
                else:
                    if added_transcripts[transcript_id_clean].upper() == sequence_clean.upper():
                        continue
                    else:
                        raise ValueError("Value of sequence does not match previously added transcript")

    # Now align with HISAT2 and create new, remapped annotations
    # Expected 

    os.system("hisat2 -x " + index_name + " -f extracted_seq.fa" + 
    " --mp 10000,10000  --no-softclip --rfg 10000,10000 --no-spliced-alignment -a -S seq.sam")

    with open("seq.sam", "r") as sam, open(output_gtf, "w") as new:
        
        num_transcripts = my_dictionary()

        # If flag = 4, then add to unaligned transcripts
        columns = ["transcript_info", "sequence"]
        unaligned_transcripts = pd.DataFrame(columns = columns)

        # I think this could be useful one day
        flag_key = {4:"no_alignments", 16:"-", 256:"not_primary", 276:"-"}
        
        for line in sam:
            # Skip comment lines
            if line[0] == "@":
                continue
            else:
                # Read in the sam file format
                sep_line = line.split(sep = "\t")
                
                # If unassigned, mark and then continue
                if sep_line[1] == '4':
                    unaligned_data = pd.DataFrame([{"transcript_info":sep_line, "sequence":sep_line[9]}])
                    unaligned_transcripts = pd.concat[unaligned_transcripts, unaligned_data]
                    continue

                # Mark if on sense or missense
                elif sep_line[1] == '256' or sep_line[1] == '0':
                    strand = "+"
                elif sep_line[1] == '16' or sep_line[1] == '272':
                    strand = "-"

                # Extract transcript_id, biotype, and database

                feature_one = sep_line[0]
                feature_one_sep = feature_one.split(sep = "--")

                transcript_id = feature_one_sep[0]
                biotype = feature_one_sep[1]
                database = feature_one_sep[2]

                if transcript_id not in num_transcripts:
                    # Add to dictionary
                    num_transcripts.add(transcript_id, 1)

                else: 
                    num_transcripts[transcript_id] += 1

                transcript_copy_id = transcript_id + "_" + str(num_transcripts[transcript_id])

                start_loci = sep_line[3]
                length_of_sequence = len(sep_line[9])
                end_loci = str(int(start_loci) + length_of_sequence - 1)

                sequence = added_transcripts[transcript_id]

                seq_type = "exon"

                chromosome = sep_line[2]

                attributes = ('transcript_id "' + transcript_id + '"; ' + 'transcript_copy_id "' + 
                            transcript_copy_id + '"; ' + 'sequence "' + sequence + '"; ' + 'biotype "' +
                            biotype + '"')
                feature_list = [chromosome, database, seq_type, start_loci, end_loci, ".", strand, ".", attributes]

                new.write("\t".join(feature_list) + "\n") """

# This function will take a gtf file and assign
# biotype of a transcript based upon the origin database.
# It will do this given a key
# It should also reorder
""" def key_biotype_gtf(gtf_file, gtf_output, a_dict):
    with open(gtf_file, "r") as gtf, open(gtf_output, "w") as new:
        for line in gtf:
            sep_line = line.split(sep = "\t")
            # get only attribute info
            attributes = sep_line[8].split(sep = ";")
            # strip
            attributes = list(map(str.strip, attributes))
            # combine
            attributes_split = [item.split(sep = " ") for item in attributes]
            # combine
            attributes_split = sum(attributes_split, [])

            data_origin = sep_line[1]

            biotype = a_dict[data_origin]

            attributes_split.append("biotype")

            attributes_split.append('"' + biotype + '"')

            index_ti = attributes_split.index("transcript_id")
            try:
                index_tcp = attributes_split.index("transcript_copy_id")
                tcp_present = True
            except:
                tcp_present = False
            index_seq = attributes_split.index("sequence")
            index_biotype = attributes_split.index("biotype")

            ti = attributes_split[index_ti + 1] 
            if tcp_present == True:
                tcp = attributes_split[index_tcp + 1] 
            seq = attributes_split[index_seq + 1]
            biotype = attributes_split[index_biotype + 1] 

            if tcp_present == True:
                final_entry = ["transcript_id " + ti + "; ","transcript_copy_id " + tcp + "; ", 
                        "sequence " + seq + "; ", "biotype " + biotype ]
            elif tcp_present == False:
                final_entry = ["transcript_id " + ti + "; ", 
                        "sequence " + seq + "; ", "biotype " + biotype ]

            new.write("\t".join(sep_line[0:8] + ["".join(final_entry)]) + "\n")
 """
""" # Given a gtf file and a dictionary that gives attribute identifier
# relationships
# i.e. {"name":"transcript_id"}
# Add copy ID with HISAT
# The key exists within the old dataset
def standardize_attributes(gtf_file, gtf_output, a_dict):
    with open(gtf_file, "r") as gtf, open(gtf_output, "w") as new:
        for line in gtf:
            sep_line = line.split(sep = "\t")
            # get only attribute info
            attributes = sep_line[8].split(sep = ";")
            # strip
            attributes = list(map(str.strip, attributes))
            # combine
            attributes_split = [item.split(sep = " ") for item in attributes]
            # combine
            attributes_split = sum(attributes_split, [])

            new_attributes = []

            for key in a_dict:
                the_index = attributes_split.index(key)
                new_attributes.append(a_dict[key])
                new_attributes.append(attributes_split[the_index + 1])

            iter = 0

            attribute_combined = ""
            length_attributes = len(new_attributes)

            for item in new_attributes:
                if iter % 2 == 0:
                    attribute_combined = attribute_combined + item + " "

                else:
                    if iter + 1 == length_attributes:
                        attribute_combined = attribute_combined  + item 
                        continue
                    else:
                        attribute_combined = attribute_combined  + item + '; '

                iter += 1
            
            sep_line[8] = attribute_combined
            
            new.write("\t".join(sep_line) + "\n")

    return gtf_output """

#def sequence_length(gtf, field):

""" def gtf_to_bed(gtf, outputname, attribute_to_name = False):
    with open(gtf, "r") as gtf, open(outputname, "w") as new:
        for line in gtf:
            modify_line = line.split(sep="\t")
            modify_line[3] = int(modify_line[3])
            modify_line[3] = modify_line[3] - 1
            modify_line[3] = str(modify_line[3])

            bed_line = []

            bed_line.append(modify_line[0])
            bed_line.append(modify_line[3])
            bed_line.append(modify_line[4])
            # obtain name info
            if attribute_to_name != False:
                attributes = modify_line[8].split(sep = ";")
                # strip
                attributes = list(map(str.strip, attributes))
                # combine
                attributes_split = [item.split(sep = " ") for item in attributes]
                # combine
                attributes_split = sum(attributes_split, [])
                index_id = attributes_split.index(attribute_to_name)
                bed_line.append(attributes_split[index_id + 1])
            else:
                bed_line.append(".")
            # add the length of sequence here later
            bed_line.append(".")
            bed_line.append(modify_line[6])


            new.write("\t".join(bed_line) + "\n") """

# extract info from a gtf file into a tsv for downstream analysis
#def gtf_to_tsv(gtf_file, output_name, dict_extract):



#===Working Playground===#
#biotype_cov = {"pirnadb_v1_7_6":"piRNA", "UCSC":"rRNA", "GtRNAdb_v2_0":"tRNA_full", "miRBase_v22_1_mature_miRNA":"miRNA"}
#key_biotype_gtf("final/human_all.gtf", "final/human_all_biotype.gtf", biotype_cov)

# We then use command line to combine files

# zero_to_one("final/human_all_biotype.gtf", "human_all_biotype.gtf")

# merge_overlaps("synthetic/SPRMT_hg38_031323.gtf", "SPRMT_merge.tsv", "SPRMT_merge_031323.gtf", 0)

# add_sequence("SPRMT_merge_031323.gtf", "hg38_complete.fa", "synthetic/SPRMT_merge_031323_seq.gtf")

# Obtain miRNA in gtf format

# gff3_to_gtf("synthetic/hsa-miRNA.gff3", "synthetic/hsa-miRNA.gtf")

# select_column("synthetic/hsa-miRNA.gtf", "synthetic/hsa-pre-miRNA.gtf", 2, "miRNA_primary_transcript")

# add_sequence("hsa-pre-miRNA.gtf", "hg38_complete.fa", "hsa_pre_mirna.gtf")

#filter_column("SPRMT_merge_031323_seq.gtf", "SPRMT_merge_031323_seq_remove_merged.gtf", 1, "merged")
#select_column("SPRMT_merge_031323_seq.gtf", "SPRMT_merge_031323_seq_merged_only.gtf", 1, "merged")

# pre miRNA processing

#the_gtf = "hsa-pre-miRNA.gtf"

#standardize_attributes(the_gtf, "new.gtf", {"Name":"transcript_id"})
#add_sequence("new.gtf", "hg38_complete.fa", "new_seq.gtf", "sequence")
#key_biotype_gtf("new_seq.gtf", "hsa-pre_miRNA.gtf", {".":"pre_miRNA"})
#align_hisat("hsa-pre_miRNA.gtf", "hsa_pre_miRNA.gtf", "hg38comp")
#add_sequence("hsa_pre_miRNA.gtf", "hg38_complete.fa", "hsa_pre_miRNA_final.gtf")

#standardize_attributes("SPRMT_merge_031323_seq.gtf", "SPRMT_merge_031323_seq_manatee.gtf", {"transcript_copy_id":"gene_id","transcript_id":"gene_name","biotype":"gene_biotype"})

# Generate unique sequences
#select_column("RC_Fixed_Chr.gtf", "RNA_Central_NCexons.gtf", 2, "noncoding_exon")
#gtf_to_bed("RNA_Central_NCexons.gtf", "RNA_Central_NCexons.bed")
#add_sequence_bed("RNA_Central_NCexons.bed", "testing_env/hg38_complete.fa", "RNA_Central.bedseq")


# add_sequence("RNA_Central_NCexons.gtf", "testing_env/hg38_complete.fa", "RC_NCexons_Seqs.gtf", "sequence")

#import 

# tRNA fragments
# gtf_naming_stan("testing_env/human_tRNA_derived.gtf", "testing_env/hg38.chromAlias.txt", "tRNA_derived_chr.gtf")
#zero_to_one("tRNA_derived_chr.gtf", "tRNA_coord_fix.gtf")
#add_sequence_gtf("tRNA_coord_fix.gtf", "testing_env/hg38_complete.fa", "tRNA_derived_seq.gtf", "bed_sequence")
#compare_sequence("tRNA_derived_seq.gtf", "trna_check.csv", "sequence", "bed_sequence")