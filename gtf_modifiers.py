from basics import *

# Description:
# Functions to modify gtf files.
# - Add sequence
# - Change the exon attribute to something elsse
# - filter gtf files by attribute
# - filter or select entries by column values
# - merge overlaps
# - generate new annotation with sequences

# add a feature that extracts the sequence from the genome
# Note that lower case letters represent masked regions
# It only takes a gtf file 
# BedTools recognizes this
def add_sequence_legacy(gtf_file, ref_genome, output_name):
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

    print("--- %s seconds ---" % (time.time() - start_time))

# faster version of the above - less safe and assumes that 
# lines will correspond to each other i.e. 
# 1 in gtf = 1 i in the output
# should be true
# NOTE Bed = gtf, it can accept both
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

# we define a new file type, bedseq format
def add_sequence_bed(bed_file, ref_genome, output_name):

    seq_list = add_sequence_fast(bed_file, ref_genome)

    with open(bed_file, "r") as gtf, open(output_name, "w") as new:

        iter = 1

        for line in gtf:
            extracted_sequence = seq_list[iter]

            temp_line = line.strip()

            new.write(temp_line + "\t" + extracted_sequence + "\n")

            iter += 2

# Modify the exon value to something else
def gtf_change_middle(gtf_file, out_name, change_value):
    with open(gtf_file, "r") as gtf, open(out_name, "w") as new:
        for line in gtf:
            modify_line = line.split(sep="\t")
            modify_line[2] = str(change_value)

            new.write("\t".join(modify_line))

# Filter by attribute
# Will take out whatever attribute you put in
def filter_gtf(input_gtf, output_name, filter_by, value, skip = False, num = 1):
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
            type_feature = type_feature.replace('"', "")
            if type_feature == value:
                continue
            else: 
                new.write(line)

# Will filter by a column
# 0 index
def select_column(input_gtf, output_name, col_number, value):
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
                new.write(line)

# Merge Overlapping Transcripts
def merge_overlaps(gtf_file, output_tsv, output_gtf, offset):
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

                    x = 1

                    for key in combined_attribute:
                        if x == attribute_length:
                            add_atribute.append(key + " ")
                            add_atribute.append(combined_attribute[key] + '"')
                        else:
                            add_atribute.append(key + " ")
                            add_atribute.append(combined_attribute[key] + '"; ')

                        x += 1

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

                    attribute_length = len(combined_attribute)

                    x = 1

                    for key in combined_attribute:
                        if x == attribute_length:
                            add_atribute.append(key + " ")
                            add_atribute.append(combined_attribute[key] + '"')
                        else:
                            add_atribute.append(key + " ")
                            add_atribute.append(combined_attribute[key] + '"; ')

                        x += 1


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

    return [stored_data]

# This function will take a gtf file and assign
# biotype of a transcript based upon the origin database.
# It will do this given a key
# It should also reorder - make it so that the attributes are not hardcoded
def key_biotype_gtf(gtf_file, gtf_output, a_dict):
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

# Given a gtf file and a dictionary that gives attribute identifier
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

    return gtf_output

# Bin Annotation File Into Segment