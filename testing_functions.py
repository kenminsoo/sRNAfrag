from basics import *
from alias_work import *
from gtf_descriptors import *
from gtf_modifiers import *
from gtf_generation import *
from conversion_tools import *
from gtf_groundtruth import *
from fragmentation import *

# complete
#tsv_to_gtf("testing_env/snoDB_All_V2.0.tsv", "testing.tsv", extract_dictionary=extract_dictionary)

# Working directory - All intermediate files will be deleted
working_dir = "/Volumes/Extreme_SSD/snoRNA_frag_exosome/working"

# Sample directory
sample_dir = "/Volumes/Extreme_SSD/snoRNA_frag_exosome/raw_data"

# Reference genome location
build_bool = False
reference_genome = ""

# Alignment index name or location
index_name = "/Users/kenminsoo/Desktop/Projects/JABSOM/Nakatsu/snoRNA_Breast/smallRNA_structure/output/adaptertrim/hg38"

# Out directory with all final files
out_dir = "/Volumes/Extreme_SSD/snoRNA_frag_exosome/out"

# Annotation file to fragment, should be filtered for one biotype
annotation_file = "/Users/kenminsoo/Desktop/Projects/JABSOM/Nakatsu/projects/tools/snoRNA_std_box_v2.gtf"
add_sequences_bool = False

# Remove UMIs, True or False
remove_umis_bool = True

# UMI Pattern
umi_pattern = ".+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12}).+"

# Remove Adapters
remove_adapters_bool = True

# Adapter Settings
adapter = "AGATCGGAAGAG"
min_length = 15
max_length = 45

# Quality Check?
fastqc_bool = False

# Max cores
max_cores = 12

# Keep This Amount of Memory Free
mem_free = 16
max_cores = 6
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

#generate_groundtruth("piRNAs/two_piRNAs_v2.gtf", "testing_env/hg38_complete.fa", "piRNAs", "transcript_id", {"source":[0, 1]}, p = 16)
#generate_groundtruth("piRNAs/two_piRNAs_v2.gtf", "testing_env/hg38_complete.fa", "piRNAs_piRNAbank", "transcript_id", {"source":[0, 1]}, input_fasta="piRNAbank.fasta", p = 16)


# Add sequences to all
#add_sequence_gtf("piRNAs/piRNAbank.gtf", "testing_env/hg38_complete.fa", "piRNAs/piRNAbank_v2.gtf", "sequence")
#add_sequence_gtf("piRNAs/pirnadb_v1.gtf", "testing_env/hg38_complete.fa", "piRNAs/pirnadb_v2.gtf", "sequence")
#add_sequence_gtf("piRNAs/pirbase.gtf", "testing_env/hg38_complete.fa", "piRNAs/pirbase_v2.gtf", "sequence")


# Convert tRNA
#bed_to_gtf("hg38-tRNAs.bed", "tRNAs.gtf", "tRNAbase", "tRNA")
#bin_gtf("tRNAs.gtf", "binned_tRNAs.gtf", 10, "biotype")

#select_column("binned_tRNAs.gtf", "plus_trna.gtf", 6, "+")
#select_column("binned_tRNAs.gtf", "neg_trna.gtf", 6, "-")
## Deduplication BAM files ##
# Create Lookup Table
# create_lookup(annotation_file, "sequence", "new_tid", min_length, max_length, working_dir + "/snoRNA_frag_lookup.csv")