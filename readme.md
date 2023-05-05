# Features to Implement & Fixes to Make

## Fixes

## Features to Implement
1) If sequence exists from external source, validate 1-indexing (or 0)
2) Create a way to externally validate sequences pulled from the genome
3) Create a more interactive way to fill the parameters for some functions
4) Eventually add the ability to generate a SQL relational database that has information on secondary structure, mismatch mappings, and a way to merge similar sequences. 
5) AASRA alignment if the differences are represented as 2-3 n.t. 
6) Merge sequences as one if the difference between databases <= 1. 

# Tools For Working Environment and GTF Files

## Contents
1) Tools for managing environments
2) Tools for working with gtf files 

## Example Usage for (2)

We begin with a gff file, example.gff and we want to filter lncRNAs from db3, label biotypes, extract sequences, and cluster ncRNAs that are either overlapping or within 1kb of each other. 

First convert to gtf. 

gff3_to_gtf(example.gff, example.gtf) => example.gtf in the gtf format

Then, standarize the chromosome names to UCSC format which is the default mode of the function (customizable will come later).

gtf_naming_stan(example.gtf, hg38.alias.txt, example_standardized.gtf)

Then we want to add biotypes as an attribute to the gtf file. 

key_biotype_gtf(example_standardized.gtf, example_bio_stan.gtf, {db1:rRNA, db2:tRNA, db3:piRNA, db4:lncRNA})

Then, we want to filter out all lncRNA biotypes from our annotation file, and skip the first line due it containing a comment. 

filter_gtf(example_bio_stan.gtf, example_filtered.gtf, biotype, lncRNA, skip = True, num = 1)

Then we want to take all unique sequences and transcript ids in this file, and align it to a reference genome in  order to ensure that no alternate loci exist in other haplotypes. This helps to reduce bias in downstream analysis when counting features. 

align_hisat(example_filtered.gtf, alt_example_filtered.gtf, humangenome_index)

Then we want to merge all overlapping sequences that are within 1kb of each other in order to create ncRNA gene clusters. 

merge_overlaps(alt_example_filtered.gtf, information.tsv, merged_example.gtf, 1000)

Then we want to extract sequences for our newly clustered genes, along with other things. 

add_sequence(merged_example.gtf, hg38.fa, seq_merged_example.gtf)

We now have a merged annotation file with 1kb reach to annotate reads. 

## Tools for Managing Environments

Tools for managing environments consists of three scripts.
a) Conda install 
    Takes a text file of conda packages and installs them
b) Java install 
    Takes a text file of java packages and installs them
c) R install
    Installs R packages
    
Purpose: The previous scripts will hopefully one day be used in combination with a tool that pulls packages used in publications to install them with ease. It will allow for the easy recreation of environments, complete with versions. This will help save time when setting up environments to recreate analyses and to improve reproducability in the computational community.

## Tools for working with gtf files

I hope to eventually turn this python script for working with gtf files into a package.
It consists of the following functions. 

### Tool Categories
1) alias_work.py - to work with chromosomal aliases
2) gtf_descriptors.py - to describe annotation files, i.e. the number of times certain attributes appear
3) gtf_generation.py - to create new annotation files from sequence
4) gtf_combining_metrics.py - to compare how merging annotation files would work, # overlaps, total number, conflicts
5) basics.py - basic functions
6) gtf_modifiers.py - modify existing gtf files, i.e. changing keys, adding parent-child relationships
7) gtf_groundtruth.py - generate a ground truth dataset of different levels of mismatching of sequences
8) gtf_tracker.py - tracks how any annotation files will be changed so that at the end, one can obtain metrics about what was done
9) conversion_tools.py - convert files to gtf, change index system


a) hamming(seq1, seq2) 
    Takes two sequences and computes the hamming distance. Is mainly used for calculating if two extracted sequences are the same. 
    
b) rna_trans(transcript)
    Takes a DNA sequence and transcribes it to RNA.

c) gff3_to_gtf(gff3_file, output_name)
    Takes a gff3 file and converts it to a gtf file.
    
d) add_sequence_legacy(gtf_file, ref_genome, output_name)
    Takes a gtf file and a reference genome. It will output a gtf file with sequences extracted from the gtf file. This is a slower version of the proceeding function, but is safer in that it works line by line. It is about 200x slower though. 

e) gtf_to_bed(gtf, outputname, attribute_to_name = False)
    Takes a gtf file and turns it into a BED file. The user can specify a certain attribute that they want to keep in the name column of bed files. 

f) compare_sequence(gtf_file, output_name, sequence_feature1, sequence_feature2)
    It takes a gtf_file and expects it to have two features in its attribute column with sequence information. This is typically used to compare one seqeuence feature, possibly from a database of sequences, to a second sequence feature which is extracted by a tool such as the one above. Works to confirm that the sequence in the database is the same as the loci being referred to by the annotation file. 
    
g) tsv_to_gtf(tsv, out_name, extract_list)
    An interactive walk through converting a tsv to a gtf file. 
    
h) ref_combine(fa1_name, fa2_name, out_name, reference)
    This function takes two reference genomes in the fasta format. It also requires a reference text file of aliases for chromosomes. This function can be used to merge UCSC reference genomes with NCBI ones, for example. It is useful mainly for converting alternative haplotypes. Note that is a slow function.
    
i) gtf_naming_stan(gtf_file, reference, out_name)
    This function takes a gtf file and a reference text file of aliases for chromosomes. This then changes all gtf naming of chromosomes to UCSC format. In a future version we may have it so that the user can specify what column of the reference text file to convert to. However, oftentimes there are empty entries which can make this process hard. 

j) zero_to_one(gtf_file, out_name)
    This changes the index from zero based to one based. (add one to start)
    
k) one_to_zero(gtf_file, out_name)
    This changes the index from one based to zero based. (subtract one from start)
    
l) gtf_change_middle(gtf_file, out_name, change_value)
    This modifies the "Exon" value. 

m) filter_gtf(input_gtf, output_name, filter_by, value, skip = False, num = 1)
    This functions takes in an input gtf file and outputs a filtered file. This should will filter OUT a certain attribute in the final column of a gtf file. 

n) countby_field(input_gtf, output_name, field_index)
    This function takes in an input gtf file and outputs how many times a certain field appears. This can be used to count the number annotations that appear on each chromsome. Count by column.

o) countby_attribute(input_gtf, output_name, countby_value, skip = False, num = 1)
    This function takes in a gtf file and counts by a certain attribute. 
    
p) merge_overlaps(gtf_file, output_tsv, output_gtf, offset)
    This function takes in a gtf file and outputs two files. It outputs a tsv that shows the number of transcripts that each annotation overlaps with. The offset parameter extends the overlapping query region and allows one to cluster genes that are located in a certain number of basepairs from a annotation. 

q) align_hisat(gtf_file, output_gtf, index_name)
    This function takes in a gtf file with the following fields:
    1) transcript_id
    2) biotype
    3) sequence
    And requires hisat2 installed on the OS in either the working environment or path. It aligns sequences in the gtf file with the the specific reference genome which should be built prior to running this. index_name takes in the name of this reference genome. It outputs a gtf file with the 1-indexed loci and infers end loci by the following formula: start loci + length - 1. Hisat aligns with the following command: hisat2 -f <gen_fasta> --mp 10000,10000 --no-softclip --rfg 10000,10000 --no-spliced-alignment -a
    Output gtf attributes:
    1) transcript_id
    2) transcript_copy_id
    3) biotype
    4) sequence

r) key_biotype_gtf(gtf_file, gtf_output, a_dict)
    This function takes in a gtf file and a dictionary in the following format:
        a_dict = {database;biotype}
    This only works in specific cases when one pulls data from different databases and wants to merge them to together. 

s) select_column(input_gtf, output_name, col_number, value)
    This function takes in a gtf file, the index of the column (0-based indexing), the output name of the file, and the value that you want to keep. It is a select function, so for example select all that are in. 
    
t) filter_column(input_gtf, output_name, col_number, value)
    This function takes in a gtf file, the index of the column (0-based indexing), the output name of the file, and the value that you want to filter out. It is a filter function. 
    
u) standardize_attributes(gtf_file, gtf_output, a_dict)
    Takes a gtf and a dictionary. The keys in the dictionary should be attribute names from the initial gtf files. The value associated should be what you want to change it to. 
    name => transcript_id should be {"name":"transcript_Id"}
    Key should be in the gtf attributes column. 

v) add_sequence_bed(bed_file, ref_genome, output_name)
    Takes a bed file and a reference genome. It will output a bedseq file which adds an additional sequence column to use in downstream purposes. It will also put the length of the sequence into the score section. 

w) add_sequence_gtf(gtf_file, ref_genome, output_name, attribute_name = "bed_sequence")
    Takes a gtf file and a reference genome. It will output a gtf file with sequences extracted from the gtf file. Works in seconds compared to hours. 

x) gtf_to_tsv(gtf_file, output_name, dict_extract)
    Takes a gtf and a dictionary in the following format:
    {"attributes":["list", "of", "attributes", "to extract"], "columns": [0 index columns to extract]}

y) ref_chr_select(ref_genome, out_name,chr_list) or gtf_chr_select(gtf, out_name,chr_list)
    Both take a list of chromosomes and filters them. Default is the standard human chromosomes without haplotypes. 