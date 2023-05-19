from basics import *
from alias_work import *
from gtf_descriptors import *
from gtf_modifiers import *
from gtf_generation import *
from conversion_tools import *
from gtf_groundtruth import *

# complete
#tsv_to_gtf("testing_env/snoDB_All_V2.0.tsv", "testing.tsv", extract_dictionary=extract_dictionary)


# fixed
#merge_overlaps("overlap_testing.gtf", "overlap.tsv", "overlap_output.gtf", 0)

# create tRNA from BED
# do at a later time



# Merge and create a big annotation

# 1) Turn everything into gtf format

# piRNAdb

# needed to take out a ; from the function
## key_biotype_gtf("piRNAs/pirnadb.v1_7_6.hg38.gtf", "piRNAs/pirnadb.gtf", {"piRNAdb":"piRNA"})

# Standardize the entries to biotype & tid
##standardize_attributes("piRNAs/pirnadb.gtf", "piRNAs/pirnadb_v1.gtf", {"piRNA_code": "transcript_id", "biotype" :"biotype"})
 

# pirbase
## bed_to_gtf("piRNAs/hsa.align_pirbase.bed", "piRNAs/pirbase.gtf", "piRbase", "piRNA")

# piRNABank

# Manually we need to do this'
""" cov = {"Plus":"+", "Minus":"-"}
with open("piRNAs/human_pir_piRNABank.txt", "r") as pibank, open("piRNAs/piRNAbank.fasta", "w") as fasta, open("piRNAs/piRNAbank.gtf", "w") as gtf:
    for line in pibank:
        if line[0] == ">":
            #extract information
            real_line = line[1:]

            sep = real_line.split(sep = ":")

            name_of_transcript = sep[0].split(sep = "|")[0]

            biotype = "piRNA"

            start  = str(sep[2])
            end = str(sep[3])
            chr = sep[1]
            chr = "chr" + str(chr)

            strand = sep[4].replace("\n", "")

            strand = cov[strand]

            attribute = ["transcript_id " + '"' + name_of_transcript + '"', "biotype " + '"' + biotype + '"']

            new_line = [chr, "piRNABank", "piRNA", start, end, ".", strand, ".", "; ".join(attribute)]

            gtf.write("\t".join(new_line )+ "\n")

        else:
            sequence = line.replace("U", "T")

            entry = ">" + name_of_transcript + "\n" + sequence

            fasta.write(entry)
 """

# compare pirbase and pirnadb

#gtf_naming_stan("piRNAs/two_piRNAs_v3.gtf", "testing_env/hg38.chromAlias.txt", "piRNAs/two_piRNAs_v4.gtf")
#gtf_chr_select("piRNAs/two_piRNAs.gtf", "piRNAs/two_piRNAs_v2.gtf", chr)

generate_groundtruth("piRNAs/two_piRNAs_v2.gtf", "testing_env/hg38_complete.fa", "piRNAs", "transcript_id", {"source":[0, 1]}, p = 16)
generate_groundtruth("piRNAs/two_piRNAs_v2.gtf", "testing_env/hg38_complete.fa", "piRNAs_piRNAbank", "transcript_id", {"source":[0, 1]}, input_fasta="piRNAbank.fasta", p = 16)


# Add sequences to all
#add_sequence_gtf("piRNAs/piRNAbank.gtf", "testing_env/hg38_complete.fa", "piRNAs/piRNAbank_v2.gtf", "sequence")
#add_sequence_gtf("piRNAs/pirnadb_v1.gtf", "testing_env/hg38_complete.fa", "piRNAs/pirnadb_v2.gtf", "sequence")
#add_sequence_gtf("piRNAs/pirbase.gtf", "testing_env/hg38_complete.fa", "piRNAs/pirbase_v2.gtf", "sequence")


# Convert tRNA
#bed_to_gtf("hg38-tRNAs.bed", "tRNAs.gtf", "tRNAbase", "tRNA")
#bin_gtf("tRNAs.gtf", "binned_tRNAs.gtf", 10, "biotype")

#select_column("binned_tRNAs.gtf", "plus_trna.gtf", 6, "+")
#select_column("binned_tRNAs.gtf", "neg_trna.gtf", 6, "-")