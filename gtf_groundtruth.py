from basics import *
from gtf_modifiers import *
from alias_work import *
import glob

# Description:
# Create a ground truth database that
# represents how different transcripts
# Will map to each other with varying
# Mismatch levels

reference_genome = "hg38_std.fa"

# information dictionary should be in the format
# note zero index for column
# {info:[gen=0 or attribute=1, column index or name if attribute]}
pd.set_option("display.max.colwidth", 10000000)

def generate_groundtruth(gtf, ref, new_dir, primary_key, information_dict, hisat = False, num_mismatch = 2, input_fasta = "sequences.fasta ", p = 1):
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
    df_tooshort = primary_df.loc[primary_df["length"] <= 10]
    primary_df = primary_df.loc[primary_df["length"] > 10]

    # filter out too long transcripts
    df_toolong =  primary_df.loc[primary_df["length"] >= 1000]
    primary_df = primary_df.loc[primary_df["length"] < 1000]

    # make capital
    primary_df["sequence"] = primary_df["sequence"].str.upper()

    # filter out transcripts with N in sequence
    df_n = primary_df.loc[primary_df["sequence"].str.contains("N")]
    primary_df = primary_df.loc[~(primary_df["sequence"].str.contains("N"))]

    # remove duplicates
    primary_df = eliminate_duplicates(primary_df, information_dict, new_dir)

    primary_df.to_csv(new_dir + "/" + "primary.csv", index = False)
    df_tooshort.to_csv(new_dir + "/tooshort.csv", index = False)
    df_toolong.to_csv(new_dir + "/toolong.csv", index = False)
    df_n.to_csv(new_dir + "/n_exist.csv", index = False)

    # Begin the fasta creation 

    fasta_entry = ">" + primary_df["transcript_id"] + "SPACER" + primary_df["sequence"]


    with open(new_dir + "/sequences.fasta", "w") as f:
        stringed = fasta_entry.to_string(header = False, index = False)
        stringed = stringed.replace(" ", "")
        stringed = stringed.replace("SPACER", "\n")
        f.write(stringed)
    # Begin the alignment section


    # Check if index is built
    built = index_exist(new_dir)

    # 1) Check if HISAT2 is installed

    if hisat == True:

        try:
            os.system("hisat2")

        except:
            raise ModuleNotFoundError("HISAT2 is not installed.")
        
        # 2) Check if indexed

        # Build the index if it doesn't exist
        # Give the user the choice to rebuild if needed
        if built == True:
            print("Index has previously been built")
            print("Do you want to rebuild the index?")

            rebuild = input("y/n")

            if rebuild == "y":
                generate_index_hisat(new_dir)
            else:
                print("Continuing")
        else:
            generate_index_hisat(new_dir)

        # 3) Begin the alignment
        hisat_align(new_dir, num_mismatch, input_fasta)

    # Use bowtie by default => better results
    else:
        try:
            os.system("bowtie2")

        except:
            raise ModuleNotFoundError("bowtie2 is not installed.")
        
        if built == True:
            print("Index has previously been built")
            print("Do you want to rebuild the index?")

            rebuild = input("y/n")

            if rebuild == "y":
                generate_index_bowtie2(new_dir)
            else:
                print("Continuing")
        else:
            generate_index_bowtie2(new_dir)

        # Align

        bowtie_align(new_dir, num_mismatch, input_fasta)

    # Parse the SAM file to create the mapping
    parse_sam(new_dir, num_mismatch)

    # Create the table
    create_table(new_dir, num_mismatch)

# Create the primary table
def create_primary_table(gtf, ref, primary_key, information_dict, new_dir):
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

    return primary_df

# Generate an index
def generate_index_hisat(new_dir):
    os.system("hisat2-build " + new_dir + "/" + "sequences.fasta" + " " + new_dir + "/" + "seqs")

def generate_index_bowtie2(new_dir):
    os.system("bowtie2-build " + new_dir + "/" + "sequences.fasta" + " " + new_dir + "/" + "seqs")

def generate_index_bowtie(fasta, outname):
    os.system("bowtie-build " + fasta + " " + outname)

# Check if index exists
def index_exist(new_dir):
    files = glob.glob(new_dir + "/" + "seqs.*")

    num_index = len(files)

    if num_index >= 6:
        return True
    else:
        return False
    

# Align
def hisat_align(new_dir, num_mismatch = 2, input_fasta = "sequences.fasta ", p = 1):
    os.system("hisat2 --rfg 50000,50000 --rdg 50000,50000 --np 50000 --no-softclip --mp 20,20 -k 9223372036854775807 --max-seeds 9223372036854775807 --rna-strandness F --no-spliced-alignment --secondary -p " + str(p) + " --score-min L," +
               str(num_mismatch*-20) + ",0 -x " + new_dir + "/seqs -f " + input_fasta +
               " -S " + new_dir + "/sequences.sam")
    
index = "/seqs"

def bowtie_align(new_dir, num_mismatch = 2, input_fasta = "sequences.fasta ", p = 1, index = index):

    if index == "/seqs":
        index2 = new_dir + index
    else:
        index2 = index

    os.system("bowtie2 --ma 0 --rfg 50000,50000 --rdg 50000,50000 --np 50000 -L 10 -R 20 --mp 20,20 --norc -p " + str(p) + " --score-min L," +
               str(num_mismatch*-20) + ",0 -a -x " + index2 + " -U " + input_fasta + " -f " +
               "-S " + new_dir + "/sequences.sam")


def bowtie_align_pipeline(index_name, working_dir, fasta):
    os.system('cd ' + working_dir + '; \
    bowtie -f -x ' + index_name + ' ' + fasta +' -k 101 --best --strata -v 0 -S lookup_filtered.sam --reorder --norc')


# Parse SAM File
def parse_sam(new_dir, num_mismatch = 2):
    for i in range(0,num_mismatch + 1):
        os.system("fgrep 'AS:i:" + str(i * -20) + "' " + new_dir + "/sequences.sam > " + new_dir + "/seq_" + str(i) + "_mm.sam")

# Create Table
def create_table(new_dir, num_mismatch = 2):
    for i in range(0, num_mismatch + 1):
        # import dataset
        sam_df = pd.read_csv(new_dir + "/seq_" + str(i) + "_mm.sam", sep = "\t", header = None, usecols=[0,1])

        sam_df.to_csv(new_dir + "/mapping_" + str(i) + "_mm.sam")

# Generate aliases and remove duplicates

def eliminate_duplicates(df, information_dict, new_dir):

    # replace within the first function with the primary df
    dataset = df

    dataset["sequence"] = dataset["sequence"].str.upper()

    # Create the alias table


    ## Store sequence info
    sequence_tid = my_dictionary()

    alias_table = []

    # Honestly such a bad habit to iterate over rows...
    for row in dataset.iterrows():
        new_row = []
        
        # Add if it doesn't exist
        if row[1][1] not in sequence_tid:
            # Add seq to dict
            sequence_tid.add(row[1][1], "")

            # Add TID to dict
            sequence_tid[row[1][1]] = row[1][0]

            # If it doesn't exist, we map to self
            new_row.append(row[1][0])
            new_row.append(row[1][0])

            # Add in information regarding the other transcript

            i = 3
            for key in information_dict:
                new_row.append(row[1][i])
                i += 1

        # If already there, create a new row with existing tid
        else:
            # add the primary key to the main table
            new_row.append(sequence_tid[row[1][1]])
            new_row.append(row[1][0])

            # Add in information about the other transcript
            
            i = 3
            for key in information_dict:
                new_row.append(row[1][i])
                i += 1

        alias_table.append(new_row)

    # Assemble the headers
    header = ["primary_transcript_id", "alias_transcript_id"] + list(information_dict.keys())

    # Export the alias table
    alias_as = pd.DataFrame(alias_table, columns=header)

    alias_as.to_csv(new_dir + "/alias.csv",index = False)

    return df.drop_duplicates(subset = ["sequence"], keep = "first")

# Compare DF to TSV
def eliminate_duplicates_tsv(df, tsv, new_dir, tsv_source, df_source):

    # replace within the first function with the primary df
    dataset = df

    dataset2 = pd.read_csv(tsv, sep = "\t")

    dataset["sequence"] = dataset["sequence"].str.upper()

    # Create the alias table


    ## Store sequence info
    sequence_tid = my_dictionary()

    alias_table = []

    # Honestly such a bad habit to iterate over rows...
    for row in dataset.iterrows():
        new_row = []
        
        # Add if it doesn't exist
        if row[1][1] not in sequence_tid:
            # Add seq to dict
            sequence_tid.add(row[1][1], "")

            # Add TID to dict
            sequence_tid[row[1][1]] = row[1][0]

            # If it doesn't exist, we map to self
            new_row.append(row[1][0])
            new_row.append(row[1][0])

            # Add in information regarding the other transcript
            new_row.append(df_source)

        # If already there, create a new row with existing tid
        else:
            # add the primary key to the main table
            new_row.append(sequence_tid[row[1][1]])
            new_row.append(row[1][0])

            # Add in information about the other transcript
            new_row.append(df_source)

        alias_table.append(new_row)

    for row in dataset2.iterrows():
        new_row = []
        
        # Add if it doesn't exist
        if row[1][1] not in sequence_tid:
            # Add seq to dict
            sequence_tid.add(row[1][1], "")

            # Add TID to dict
            sequence_tid[row[1][1]] = row[1][0]

            # If it doesn't exist, we map to self
            new_row.append(row[1][0])
            new_row.append(row[1][0])

            # Add in information regarding the other transcript

            new_row.append(tsv_source)

        # If already there, create a new row with existing tid
        else:
            # add the primary key to the main table
            new_row.append(sequence_tid[row[1][1]])
            new_row.append(row[1][0])

            # Add in information about the other transcript
            
            new_row.append(tsv_source)

        alias_table.append(new_row)

    # Assemble the headers
    header = ["primary_transcript_id", "alias_transcript_id", "source"]

    # Export the alias table
    alias_as = pd.DataFrame(alias_table, columns=header)

    alias_as.to_csv(new_dir + "/alias.csv",index = False)

    return df.drop_duplicates(subset = ["sequence"], keep = "first")

def generate_fasta(primary_df, new_dir):
    os.system("mkdir " + new_dir)
    
    fasta_entry = ">" + primary_df["transcript_id"] + "SPACER" + primary_df["sequence"]

    with open(new_dir + "/sequences.fasta", "w") as f:
        stringed = fasta_entry.to_string(header = False, index = False)
        stringed = stringed.replace(" ", "")
        stringed = stringed.replace("SPACER", "\n")
        f.write(stringed)

if __name__ == '__main__':
  fire.Fire()