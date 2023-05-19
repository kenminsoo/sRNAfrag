from basics import *

# Description:
# We will add the following functions
# 1 vs. 0 coord
# fasta gtf comp.

# take an input fasta
# compare it with extracted sequences

def compare_fasta_gtf(gtf, fasta, key, output, sequence_attribute, new_dir):
    
    os.system("mkdir " + new_dir)
    
    with open(gtf, "r") as gtf, open(fasta, "r") as fasta, open(output, "w") as new:
        # first go through the annotation file to find internal mismatches

        internal_sequences = my_dictionary()

        internal_mismatchs = []

        for line in gtf:
            sep = separate_gtf_line(line)
            columns = sep[0]
            attributes = sep[1]

            seq_index = attributes.index(sequence_attribute)
            seq = attributes[seq_index + 1]

            key_index = attributes.index(key)
            gtf_key = attributes[key_index + 1]

            if gtf_key not in internal_sequences:
                internal_sequences.add(gtf_key,seq)

            elif internal_sequences[gtf_key] == seq:
                continue

            elif internal_sequences[gtf_key] != seq:
                mismatch_info = [gtf_key, seq, internal_sequences[gtf_key]]

                internal_mismatchs.append(mismatch_info)

        if len(internal_mismatchs) > 0:
            internal_mismatch_pd = pd.DataFrame(internal_mismatchs, columns = ["Key", "Key_Seq", "Recorded_Seq"])

            internal_mismatch_pd.to_csv(new_dir + "/internal_mismatch.csv",index = False)
            
            raise ValueError("Primary key sequences do not match.")
            
    # now we go through comparing each primary key from the fasta
    # to that of the annotation

        external_mismatchs = []

        for line in fasta:
            if line[0] == ">":
                new_line = line[1:]

                key = new_line.strip()

            else:
                sequence = line

                gtf_sequence = internal_sequences[key]

                if sequence == gtf_sequence:
                    continue

                else:
                    mismatch_info_ext = [key, sequence, gtf_sequence]
                    
                    external_mismatchs.append(mismatch_info_ext)

        mismatch_info_ext_pd = pd.DataFrame(mismatch_info_ext, columns=["key", "fasta_seq", "gtf_seq"])

        mismatch_info_ext_pd.to_csv(new_dir + "/gtf_mismatch.csv", index = False)
