from basics import *
from MINTplates import *
import numpy as np
from collections import defaultdict

# takes in a gtf with sequences

def create_lookup(gtf_file, seq_key, source_key,start, end, output, id_pre = "sdr_snoRNAdb2023_", flank = 5):
    lookup_dict = my_dictionary()

    constructor_headers = ["sequence", "sources", "num", "5p_end", "n_5p", "avg_start", "std", "flanking5p", "all_same5p","flanking3p", "all_same3p"]
    df_constructor = []

    with open(gtf_file, "r") as gtf, open(output, "w") as new:
        for line in gtf:
            sepped = separate_gtf_line(line)

            main_cols = sepped[0]
            attributes = sepped[1]

            strand = main_cols[6]

            seq_index = attributes.index(seq_key)

            seq = attributes[seq_index + 1]
            seq = seq.replace('"', "")
            seq = seq.strip(" ")
            seq = seq.upper()

            source_index = attributes.index(source_key)

            source_base = attributes[source_index + 1]
            source_base = source_base.replace('"', "")
            source_base = source_base + "__" + strand

            seq_length = len(seq)
            
            # get each length
            for k in range(start - 1, end):
                
                for i in range(0, seq_length - k):
                    segment = seq[i:i + k + 1]

                    source = str(i + 1) + ";" + str(i + k + 1) + "__" + source_base + "__.__" + str(i/seq_length)

                    # Get 3p and 5p flanking sequences, auto six

                    p5_flank = seq[(i - flank):i]

                    p3_flank = seq[i + k + 1:i + k + 1 + flank]

                    if segment not in lookup_dict:

                        if i == 0:

                            lookup_dict.add(segment, [[source], 1, 1, [i/seq_length], [p5_flank], 1, [p3_flank], 1])
                        
                        else:

                            lookup_dict.add(segment, [[source], 1, 0, [i/seq_length], [p5_flank], 1, [p3_flank], 1])

                    else:

                        lookup_dict[segment][0].append(source)

                        lookup_dict[segment][1] = lookup_dict[segment][1] + 1

                        lookup_dict[segment][3].append(i/seq_length)

                        lookup_dict[segment][4].append(p5_flank)

                        lookup_dict[segment][6].append(p3_flank)

                        if i == 0:

                            lookup_dict[segment][2] = lookup_dict[segment][2] + 1
    
    for key in lookup_dict:

        # Skip if there are more than 300 sources
        # Likely represents sequence of low complexity
        if lookup_dict[key][1] > 300:
            continue

        entry = []

        entry.append(key)

        entry.append(">".join(lookup_dict[key][0]))
        entry.append(str(lookup_dict[key][1]))

        avg = np.mean(lookup_dict[key][3])
        std = np.std(lookup_dict[key][3])

        if lookup_dict[key][2] > 0:
            entry.append(True)
        else:
            entry.append(False)

        entry.append(str(lookup_dict[key][2]))

        entry.append(str(avg))
        entry.append(str(std))

        # check if all 5p flanks are the same
        if all(element == lookup_dict[key][4][0] for element in lookup_dict[key][4]):
            entry.append(max(lookup_dict[key][4]))
            entry.append(str(1))
        else:
            entry.append(max(lookup_dict[key][4]))
            entry.append(str(0))

        # check if all 3p flanks are the same
        if all(element == lookup_dict[key][6][0] for element in lookup_dict[key][6]):
            entry.append(max(lookup_dict[key][6]))
            entry.append(str(1))
        else:
            entry.append(max(lookup_dict[key][6]))
            entry.append(str(0))

        df_constructor.append(entry)

    lookup_table = pd.DataFrame(df_constructor, columns = constructor_headers)

    lookup_table["sequence"] = lookup_table["sequence"].str.upper()

    lookup_table["flanking5p"] = lookup_table["flanking5p"].str.upper()

    lookup_table["flanking3p"] = lookup_table["flanking3p"].str.upper()

    lookup_table = lookup_table.sort_values(by = ["sequence"])

    lookup_table = lookup_table.loc[~lookup_table["sequence"].str.contains("N"),:]

    lookup_table["ID"] = lookup_table["sequence"].apply(lambda x: encode_sequence(x, id_pre))

    lookup_table.to_csv(output, index = False)

def create_lookup_selective(gtf_file, seq_key, source_key,min_length, max_length, output, lookup_builder, id_pre = "sdr_snoRNAdb2023_", flank = 5):
    lookup_dict = my_dictionary()
    selector = defaultdict(set)

    constructor_headers = ["sequence", "sources", "num", "5p_end", "n_5p", "avg_start", "std", "flanking5p", "all_same5p","flanking3p", "all_same3p"]
    df_constructor = []

    with open(lookup_builder, "r") as lookup_build:
        for line in lookup_build:
            splitted_line = line.split(sep = "\t")

            the_number = splitted_line[2]

            the_number = the_number.replace("M", "")

            the_number = int(the_number)
            
            loci = splitted_line[1]
            
            loci = int(loci)

            selector[splitted_line[0].replace('"', "")].add((loci, the_number))

    with open(gtf_file, "r") as gtf, open(output, "w") as new:
        for line in gtf:
            sepped = separate_gtf_line(line)

            main_cols = sepped[0]
            attributes = sepped[1]

            strand = main_cols[6]

            seq_index = attributes.index(seq_key)

            seq = attributes[seq_index + 1]
            seq = seq.replace('"', "")
            seq = seq.strip(" ")

            source_index = attributes.index(source_key)

            source_base = attributes[source_index + 1]
            source_base = source_base.replace('"', "")
            
            selector_key_in_gtf = source_base
            
            source_base = source_base + "__" + strand

            seq_length = len(seq)

            # make it s.t. we extract each info
            query_set_list = selector[selector_key_in_gtf]
            

            for item in query_set_list:
                start = item[0]
                start = start - 1
                length = item[1]
                end = start + length

                segment = seq[start:end]
                segment = segment.strip(" ")
                segment = segment.upper()

                if len(segment) < min_length:
                    continue
                elif len(segment) > max_length:
                    continue

                source = str(start + 1) + ";" + str(end + 1) + "__" + source_base + "__.__" + str(start/seq_length)

                p5_flank = seq[(start - flank):start]
                p3_flank = seq[end:end + flank]

                if segment in lookup_dict:

                    lookup_dict[segment][0].append(source)

                    lookup_dict[segment][1] = lookup_dict[segment][1] + 1

                    lookup_dict[segment][3].append(start/seq_length)

                    lookup_dict[segment][4].append(p5_flank)

                    lookup_dict[segment][6].append(p3_flank)

                    if start == 0:

                        lookup_dict[segment][2] = lookup_dict[segment][2] + 1

                else:

                    if start == 0:
                        lookup_dict.add(segment, [[source], 1, 1, [start/seq_length], [p5_flank], 1, [p3_flank], 1])

                    else:
                        lookup_dict.add(segment, [[source], 1, 0, [start/seq_length], [p5_flank], 1, [p3_flank], 1])

    for key in lookup_dict:

        # Skip if there are more than 300 sources
        # Likely represents sequence of low complexity
        if lookup_dict[key][1] > 300:
            continue

        entry = []

        entry.append(key)

        entry.append(">".join(lookup_dict[key][0]))
        entry.append(str(lookup_dict[key][1]))

        avg = np.mean(lookup_dict[key][3])
        std = np.std(lookup_dict[key][3])

        if lookup_dict[key][2] > 0:
            entry.append(True)
        else:
            entry.append(False)

        entry.append(str(lookup_dict[key][2]))

        entry.append(str(avg))
        entry.append(str(std))

        # check if all 5p flanks are the same
        if all(element == lookup_dict[key][4][0] for element in lookup_dict[key][4]):
            entry.append(max(lookup_dict[key][4]))
            entry.append(str(1))
        else:
            entry.append(max(lookup_dict[key][4]))
            entry.append(str(0))

        # check if all 3p flanks are the same
        if all(element == lookup_dict[key][6][0] for element in lookup_dict[key][6]):
            entry.append(max(lookup_dict[key][6]))
            entry.append(str(1))
        else:
            entry.append(max(lookup_dict[key][6]))
            entry.append(str(0))

        df_constructor.append(entry)

    lookup_table = pd.DataFrame(df_constructor, columns = constructor_headers)

    lookup_table["sequence"] = lookup_table["sequence"].str.upper()

    lookup_table["flanking5p"] = lookup_table["flanking5p"].str.upper()

    lookup_table["flanking3p"] = lookup_table["flanking3p"].str.upper()

    lookup_table = lookup_table.sort_values(by = ["sequence"])

    lookup_table = lookup_table.loc[~lookup_table["sequence"].str.contains("N"),:]

    lookup_table["ID"] = lookup_table["sequence"].apply(lambda x: encode_sequence(x, id_pre))
    
    lookup_table.to_csv(output, index = False)