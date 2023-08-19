from basics import *
from gtf_modifiers import *

# Description:
# Functions that works to convert files from one file type or indexing system to another.
# GFF3 to GTF
# TSV to GTF
# Zero to One Index
# One to Zero Index
# GTF to BED

# turn gff3 into gtf file
def gff3_to_gtf(gff3_file, output_name):
    with open(gff3_file, "r") as gff3, open(output_name, "w") as new:
        for line in gff3:

            if line[0] == "#":
                continue

            temp_line = line.strip()
            temp_line = temp_line.replace("=", ' "')
            temp_line = temp_line.replace(";", '"; ')
            temp_line = temp_line + '"'
            new.write(temp_line + "\n")

# tsv to gtf
# The TSV file have the following:
# Chromsome number
# A database source
# Name of the third column
# Start
# End
# na
# strand
# na
# attributes as a list using another dictionary of key:value

# the default dictionary for a query
extract_dictionary = {"chr":0, 
                      "source":"", 
                      "feature":"", 
                      "start":0, 
                      "end":0, 
                      "score":".", 
                      "strand":0, 
                      "frame":".", 
                      "attributes":my_dictionary()}

def tsv_to_gtf(tsv, out_name, extract_dictionary=extract_dictionary, fill_dict = True, skip_lines = 1, header = True, csv = False):
    with open(tsv, "r") as tsv_file, open(out_name, "w") as new_file:
        
        # fill the extraction dictionary
        if fill_dict == True:
            # keep track of what we can find in the tsv
            binary_dictionary = {"chr":False, 
                      "source":False, 
                      "feature":False, 
                      "start":False, 
                      "end":False, 
                      "score":True, 
                      "strand":False, 
                      "frame":True, 
                      "attributes":False}
            # first look through the header 
            if header == True:
                header_row = tsv_file.readline().strip()

                # Allow for csv files to be used as well
                if csv == False:
                    header_split = header_row.split(sep = "\t")
                elif csv == True:
                    header_split = header_row.split(sep = ",")

                header_split = [x.lower() for x in header_split]

                # automatically fill if we can find a matching index in the 
                # tsv file
                for key in binary_dictionary:
                    if binary_dictionary[key] == True:
                        continue
                    else:
                        # try to find an index where it exists
                        try:
                            header_split.index(key)
                        except:
                            binary_dictionary[key] == False
                        else:
                            extract_dictionary[key] = header_split.index(key)
                            binary_dictionary[key] = True

            # now have the user identify what columns they want
            # to use to fill in the sttributes of the
            # gtf file

            for key in binary_dictionary:
                
                # add attributes differently
                if key == "attributes":
                    attribute_complete = False

                    while attribute_complete == False:
                        print("--- Addition of Attributes (Features in the final column of gtf file) ---")

                        print("Add your gtf attributes:")

                        print("Add the attribute name")
                        attribute_name = input("String, no space: ")


                        print("Pick either a column (index) or a constant string value.")
                        attribute_value = input("please input column integer or constant string (NO NUMERIC):")
                        # convert to numeric if is a column index
                        try:
                            int(attribute_value)
                        except:
                            attribute_value = attribute_value
                        else:
                            attribute_value = int(attribute_value)
                        
                        extract_dictionary[key].add(attribute_name, attribute_value)

                        print("Would you like to add more attributes?")
                        valueok = input("y/n: ")

                        if valueok.strip() in ["y", "yes", "Y", "Yes"]:
                            continue
                        else:
                            print("Are you sure you're done?")
                            trulyok = input("y/n :")
                            if trulyok.strip().lower() in ["y", "ye", "yes"]:
                                attribute_complete = True
                            else:
                                attribute_complete = False
                                
                elif binary_dictionary[key] == True:
                    print(key + ": The column or value to be Used: " + str(extract_dictionary[key]))
                    if type(extract_dictionary[key]) == int:
                        print("name in tsv: " + header_split[extract_dictionary[key]])
                    valueok = input("Is this Ok? y/n: ")

                    if valueok.strip() in ["y", "yes", "Y", "Yes"]:
                        continue
                    # if userwants to change, then allow for change
                    else:
                        print("This is the " + key + " column")
                        print("What would you like to change it to?")
                        print("Pick either a column (index) or a constant string value.")

                        if header == True:
                            i = 0 
                            for item in header_split:
                                print(str(i) + " :" + item)
                                i += 1
                            extract_dictionary[key] = input("please input column integer or constant string (NO NUMERIC): ")
    
                        # now convert to integer if it is one
                        try:
                            int(extract_dictionary[key])
                        except:
                            binary_dictionary[key] == True
                        else:
                            extract_dictionary[key] = int(extract_dictionary[key])
                            binary_dictionary[key] == True

                else:
                    if header == True:
                        i = 0 
                        for item in header_split:
                            print(str(i) + " :" + item)
                            i += 1

                        print("This is the " + key + " column")
                        print("What would you like to have it as?")
                        print("Pick either a column (index) or a constant string value.")
                        extract_dictionary[key] = input("please input column integer or constant string (NO NUMERIC): ")

                        try:
                            int(extract_dictionary[key])
                        except:
                            binary_dictionary[key] == True
                        else:
                            extract_dictionary[key] = int(extract_dictionary[key])
                            binary_dictionary[key] == True
                    
        iter = 1
        for line in tsv_file:
            # skip n line
            if iter <= skip_lines:
                iter += 1
                continue

            if csv == False:
                features = line.split(sep = "\t")
            elif csv == True:
                features = line.split(sep = ",")

            # build the attributes value
            attributes_list = []
            for key in extract_dictionary["attributes"]:
                entry_key = key
                entry_index = extract_dictionary["attributes"][key]

                if type(entry_index) == int:
                    entry_value = '"' + features[entry_index] + '"'
                else:
                    entry_value = '"' + entry_index + '"'

                attributes_list.append(entry_key + " " + entry_value)
            
            # an ugly way to merge everything
            new_file.write("\t".join([features[extract_dictionary["chr"]], 
                                    extract_dictionary["source"], 
                                    extract_dictionary["feature"], 
                                    features[extract_dictionary["start"]], 
                                    features[extract_dictionary["end"]], 
                                    extract_dictionary["score"], 
                                    features[extract_dictionary["strand"]], 
                                    extract_dictionary["frame"], 
                                    "; ".join(attributes_list)]) + "\n")

# Convert ex

# Zero index to One - i.e. bed to gtf
def zero_to_one(gtf_file, out_name):
    with open(gtf_file, "r") as gtf, open(out_name, "w") as new:
        for line in gtf:
            modify_line = line.split(sep="\t")
            modify_line[3] = int(modify_line[3])
            modify_line[3] = modify_line[3] + 1
            modify_line[3] = str(modify_line[3])

            new.write("\t".join(modify_line))

# One index to zero - i.e. gtf to bed
def one_to_zero(gtf_file, out_name):
    with open(gtf_file, "r") as gtf, open(out_name, "w") as new:
        for line in gtf:
            modify_line = line.split(sep="\t")
            modify_line[3] = int(modify_line[3])
            modify_line[3] = modify_line[3] - 1
            modify_line[3] = str(modify_line[3])

            new.write("\t".join(modify_line))

# Gtf to BED format (True bed format)
def gtf_to_bed(gtf, outputname, attribute_to_name = False):
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


            new.write("\t".join(bed_line) + "\n")

# BED to GTF
# single database source
def bed_to_gtf(bed, output, source, biotype):
    with open(bed, "r") as bed, open(output, "w") as new:
        for line in bed:
            sep = line.split(sep = "\t")

            chr = str(sep[0])

            if "chr" not in chr:
                chr = "chr" + chr

            start = str(int(sep[1]) + 1)
            end = str(sep[2])

            id = sep[3]

            strand = sep[5].strip()
            strand = strand.replace("\n", "")

            attributes = ["transcript_id " + '"'+id+'"', "biotype " + '"'+biotype+'"']

            new_line = [chr.strip(), source.strip(), "exon", start.strip(), end.strip(), ".", strand, ".", "; ".join(attributes)]

            new.write("\t".join(new_line)+"\n")

def gtf_to_fasta(gtf, output, ref_genome, primary_key):
    seq_list = add_sequence_fast(gtf, ref_genome)
    
    with open(gtf, "r") as gtf, open(output, "w") as new:
        

        iter = 1

        for line in gtf:
            extracted_sequence = seq_list[iter]

            sep = separate_gtf_line(line)

            columns = sep[0]
            attributes = sep[1]

            primary_index = attributes.index(primary_key)
            primary = attributes[primary_index + 1]

            new.write(">" + primary + "\n" + extracted_sequence + "\n")

            iter += 2

# Turn to fasta and generate an alias file for all transcripts that have alias
def gtf_attribute_to_fasta(gtf, output, attribute, primary_key, pipeline = False):
    names = defaultdict(set)
    if pipeline:
        with open(gtf, "r") as gtf, open(output, "w") as new, open(output + ".gtf", "w") as new_gtf:
            iter = 1

            for line in gtf:

                sep = separate_gtf_line(line)

                columns = sep[0]
                attributes = sep[1]

                seq_index = attributes.index(attribute)
                extracted_sequence = attributes[seq_index + 1]
                extracted_sequence = extracted_sequence.upper()

                primary_index = attributes.index(primary_key)
                primary = attributes[primary_index + 1]

                if extracted_sequence in names:
                    names[extracted_sequence].add(primary)
                    # generate a deduplicated gtf
                else:
                    new.write(">" + primary.replace('"', "") + "\n" + extracted_sequence.replace('"', "") + "\n")
                    new_gtf.write(line)

                    names[extracted_sequence].add(primary)


            iter += 2
    else:
        with open(gtf, "r") as gtf, open(output, "w") as new:
            iter = 1

            for line in gtf:

                sep = separate_gtf_line(line)

                columns = sep[0]
                attributes = sep[1]

                seq_index = attributes.index(attribute)
                extracted_sequence = attributes[seq_index + 1]
                extracted_sequence = extracted_sequence.upper()

                primary_index = attributes.index(primary_key)
                primary = attributes[primary_index + 1]

                if extracted_sequence in names:
                    names[extracted_sequence].add(primary)
                    # generate a deduplicated gtf
                else:
                    new.write(">" + primary.strip('"') + "\n" + extracted_sequence.strip('"') + "\n")

                    names[extracted_sequence].add(primary)

    # generate alias for user reference
    alias_table_name = output + "_alias.csv"
    with open(alias_table_name, "w") as alias:
        alias.write("primary_id,alias\n")
        for entry in names:
            if len(names[entry]) != 1:
                # extract names
                k = 0
                for entry in names[entry]:
                    if k == 0:
                        first_name = entry
                        k += 1
                    else:
                        entry_out = first_name + "," + entry
                        alias.write(entry_out + "\n")




def fasta_to_tsv(fasta, output, header_name):
    with open(fasta, "r") as fasta, open(output, "w") as new:

        new.write(header_name + "\t" + "sequence")

        for line in fasta:
            if line[0] == ">":
                entry = line
                entry_name = line.replace(">", "")
                entry_name = entry_name.replace("\n", "")

                new.write("\n" + entry_name + "\t")

            else:
                entry = line.replace("\n", "")
                new.write(entry)

def tsv_to_fasta(tsv, output, key_col, seq_col, delim = "\t", skip = 1):
    with open(output, "w") as new, open(tsv, "r") as tsv:
        
        i = 1
        for line in tsv:

            if i <= skip:
                i += 1
                continue

            separate = line.split(sep = delim)

            new.write(">" + separate[key_col].replace('"', "") + "\n")

            new.write(separate[seq_col].replace('"', "") + "\n")

if __name__ == '__main__':
  fire.Fire()