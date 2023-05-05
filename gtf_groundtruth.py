from basics import *
from gtf_modifiers import *
from alias_work import *

# Description:
# Create a ground truth database that
# represents how different transcripts
# Will map to each other with varying
# Mismatch levels

reference_genome = "hg38_std.fa"

# information dictionary should be in the format
# note zero index for column
# {info:[gen=0 or attribute=1, column index or name if attribute]}
pd.set_option("display.max.colwidth", 10000)

def generate_groundtruth(gtf, ref, new_dir, primary_key, information_dict):
    seq_list = add_sequence_fast(gtf, ref)
    
    os.system("mkdir " + new_dir)

    with open(gtf, "r") as gtf:
        i = 1
        # Generate the PRIMARY table
        # Minimum Info: Primary Key: Sequence

        primary_table = []
        primary_table_headers = [primary_key, "sequence", "length"]

        for key in information_dict:
            primary_table_headers.append(key)

        for line in gtf:
            table_row = []
            
            # split the line
            line_data = separate_gtf_line(line)

            general_data = line_data[0]
            attribute_data = line_data[1]
            
            # extract the sequence
            sequence = seq_list[i]

            # get the index for the primary key
            primary_key_index = attribute_data.index(primary_key)
            primary_key_value = attribute_data[primary_key_index + 1]
            primary_key_value = primary_key_value.replace('"', '')
            # add to row
            table_row.append(primary_key_value)

            # add the sequence
            table_row.append(sequence)

            # add the length
            table_row.append(len(sequence))

            # now we can add extra data
            for key in information_dict:
                # go through the general data
                if information_dict[key][0] == 0:
                    value = general_data[information_dict[key][1]]
                    table_row.append(value)

                # go through the attribute data 
                elif information_dict[key][0] == 1:
                    value_index = attribute_data.index(information_dict[key][1])
                    value = attribute_data[value_index + 1]
                    value = value.replace('"', "")
                    table_row.append(value)

                else:
                    raise ValueError("There is an issue with your information dictionary")
            
            # now we add this line to the primary table
            primary_table.append(table_row)
        
            i += 2

    primary_df = pd.DataFrame(primary_table, columns=primary_table_headers)

    # filter out too short transcripts
    primary_df = primary_df.loc[primary_df["length"] > 10]

    primary_df.to_csv(new_dir + "/" + "primary.csv", index = False)

    # Begin the fasta creation 

    fasta_entry = ">" + primary_df["transcript_id"] + "SPACER" + primary_df["sequence"]


    with open(new_dir + "/sequences.fasta", "w") as f:
        stringed = fasta_entry.to_string(header = False, index = False)
        
        f.write(stringed)

        os.system("sed 's/SPACER/\n/g' " + new_dir + "/sequences.fasta")

    # Begin the alignment section


generate_groundtruth("RNA_Central_noLNC_std.gtf", reference_genome, "RNA_central", "transcript_id", {"biotype":[1, "biotype"], "source":[1, "databases"]})

## Standardize 
#ref_chr_select("testing_env/hg38_complete.fa", "hg38_std.fa")
#gtf_chr_select("RNA_Central_Exon_std.gtf", "RNA_Central_noLNC_std.gtf")