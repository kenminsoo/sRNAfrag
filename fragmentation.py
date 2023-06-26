from basics import *
from MINTplates import *
import numpy as np

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

            source_index = attributes.index(source_key)

            source_base = attributes[source_index + 1]
            source_base = source_base.replace('"', "")
            source_base = source_base + "_" + strand

            seq_length = len(seq)
            
            # get each length
            for k in range(start - 1, end):
                
                for i in range(0, seq_length - k):
                    segment = seq[i:i + k + 1]

                    source = str(i + 1) + ";" + str(i + k + 1) + "_" + source_base + "_._" + str(i/seq_length)

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

    lookup_table["ID"] = lookup_table["sequence"].apply(lambda x: encode_sequence(x, id_pre))

    lookup_table.to_csv(output, index = False)