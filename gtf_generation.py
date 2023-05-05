from basics import *

# Description:
# Generate new annotations.
# - Align transcripts with HISAT2 to create new annotations. 

# This requires HISAT2 to be installed and a index to be built
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

                new.write("\t".join(feature_list) + "\n")