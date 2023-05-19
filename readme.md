# Features to Implement & Fixes to Make

## Fixes

## Features to Implement
1) If sequence exists from external source, validate 1-indexing (or 0)
2) Create a way to externally validate sequences pulled from the genome
3) Create a more interactive way to fill the parameters for some functions
4) Eventually add the ability to generate a SQL relational database that has information on secondary structure, mismatch mappings, and a way to merge similar sequences. 

# Tools For Working Environment and GTF Files

## Contents
1) Tools for managing environments
2) Tools for working with gtf files 

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

## 1) alias_work.py - to work with chromosomal aliases

### `ref_combine`

Combine reference genomes.

### Parameters:
- `fa1_name` (str): The file name of the first reference genome in FASTA format.
- `fa2_name` (str): The file name of the second reference genome in FASTA format.
- `out_name` (str): The desired output file name for the combined reference genome.
- `reference` (str): The file name of a reference containing chromosome aliases and information in tabular format.

### Usage:
```python
ref_combine(fa1_name, fa2_name, out_name, reference)
```

---

### `gtf_naming_stan`

Standardize GTF (Gene Transfer Format) chromosome names when different aliases are used.

### Parameters:
- `gtf_file` (str): The file name of the GTF file to be standardized.
- `reference` (str): The file name of a reference containing chromosome aliases and information in tabular format.
- `out_name` (str): The desired output file name for the standardized GTF file.

### Usage:
```python
gtf_naming_stan(gtf_file, reference, out_name)
```

---

### `ref_chr_select`

Filter a reference genome based on a list of desired chromosomes.

### Parameters:
- `ref_genome` (str): The file name of the reference genome to be filtered.
- `out_name` (str): The desired output file name for the filtered reference genome.
- `chr_list` (list, optional): List of chromosome names to be included in the filtered output. Defaults to `['chr1', 'chr10', 'chr11', ...]`.

### Usage:
```python
ref_chr_select(ref_genome, out_name, chr_list=chr)
```

---

### `gtf_chr_select`

Filter a GTF file based on a list of desired chromosomes.

### Parameters:
- `gtf` (str): The file name of the GTF file to be filtered.
- `out_name` (str): The desired output file name for the filtered GTF file.
- `chr_list` (list, optional): List of chromosome names to be included in the filtered output. Defaults to `['chr1', 'chr10', 'chr11', ...]`.

### Usage:
```python
gtf_chr_select(gtf, out_name, chr_list=chr)
```

---

### `bam_chr_extract`

Extract chromosome names from a BAM file.

### Parameters:
- `bam` (str): The file name of the BAM file.
- `out` (str): The desired output file name to store the extracted chromosome names.

### Returns:
- `chroms` (list): A list of extracted chromosome names.

### Usage:
```python
chroms = bam_chr_extract(bam, out)
```

## 2) gtf_descriptors.py - to describe annotation files, i.e. the number of times certain attributes appear

### `compare_sequence`

Compare internal sequence with pulled sequence in a GTF file.

### Parameters:
- `gtf_file` (str): The file name of the GTF file to be analyzed.
- `output_name` (str): The desired output file name to store the comparison results.
- `sequence_feature1` (str): The first sequence feature to compare.
- `sequence_feature2` (str): The second sequence feature to compare.

### Usage:
```python
compare_sequence(gtf_file, output_name, sequence_feature1, sequence_feature2)
```

---

### `countby_field`

Count the number of occurrences of a certain column value in a GTF file.

### Parameters:
- `input_gtf` (str): The file name of the input GTF file.
- `output_name` (str): The desired output file name to store the count results.
- `field_index` (int): The index of the column to be counted.

### Returns:
- `return_data` (my_dictionary): A custom dictionary object containing the count results.

### Usage:
```python
return_data = countby_field(input_gtf, output_name, field_index)
```

---

### `countby_attribute`

Count the number of occurrences of an attribute in a GTF file.

### Parameters:
- `input_gtf` (str): The file name of the input GTF file.
- `output_name` (str): The desired output file name to store the count results.
- `countby_value` (str): The attribute value to be counted.
- `skip` (bool, optional): Flag to skip a certain number of lines in the input GTF file. Defaults to `False`.
- `num` (int, optional): The number of lines to skip if `skip` is set to `True`. Defaults to `1`.

### Returns:
- `return_data` (my_dictionary): A custom dictionary object containing the count results.

### Usage:
```python
return_data = countby_attribute(input_gtf, output_name, countby_value, skip=False, num=1)
```

## 3) gtf_generation.py - to create new annotation files from sequence

### `align_hisat_gtf`

This function aligns transcripts with HISAT2 and generates new annotations. It takes an input GTF file, aligns the transcripts using HISAT2 with a specified index, and creates a new GTF file with the aligned annotations.

### Parameters:
- `gtf_file` (str): The path to the input GTF file.
- `output_gtf` (str): The path to the output GTF file.
- `index_name` (str): The name of the HISAT2 index.

### Usage:
```python
align_hisat_gtf(gtf_file, output_gtf, index_name)
```

---

### `generate_from_fasta`

This function aligns transcripts from a FASTA file using HISAT2 and generates new annotations. It takes an input FASTA file, aligns the transcripts using HISAT2 with a specified index, and creates a new GTF file with the aligned annotations. Intermediate files are stored in the specified directory.

### Parameters:
- `fasta` (str): The path to the input FASTA file.
- `output_gtf` (str): The path to the output GTF file.
- `index_name` (str): The name of the HISAT2 index.
- `new_dir` (str): The name of the directory to store intermediate files.

### Usage:
```python
generate_from_fasta(fasta, output_gtf, index_name, new_dir)
```

---

### `bin_gtf`

This function bins a GTF file into smaller regions. It divides each sequence in the input GTF file into equal-sized bins and creates new entries for each bin. The resulting binned GTF file is stored in the output file. The `n` parameter specifies the number of bins, and the `addition_term` is a term that is appended to the number slice due to it taking on the remainder of bases after the length is divided by `n`.

### Parameters:
- `gtf` (str): The path to the input GTF file.
- `output` (str): The path to the output GTF file.
- `n` (int): The number of bins to divide each sequence into.
- `addition_term` (str): A term that is appended to the number slice.

### Usage:
```python
bin_gtf(gtf, output, n, addition_term)
```

---

### `countby_field`

This function counts the number of occurrences of values in a specific column of a GTF file. It reads the input GTF file, counts the occurrences of values in the specified column (identified by `field_index`), and writes the count data to the output file.

### Parameters:
- `input_gtf` (str): The path to the input GTF file.
- `output_name` (str): The path to the output file for storing the count data.
- `field_index` (int): The index of the field/column to count occurrences.

### Usage:
```python
countby_field(input_gtf, output_name, field_index)
```

### `countby_attribute`

This function counts the number of occurences of a particular attribute in a GTF file. The index is found after specification of the `countby_value` which should be in the 9th column of a gtf file. 

### Parameters:
- `input_gtf` (str): The path to the input GTF file.
- `output_name` (str): The path to the output file for storing the count data.
- `countby_value` (str): The attribute to count occurrences of.
- `skip` (bool): (Optional) Whether to skip a certain number of lines at the beginning of the file. Default is False.
- `num` (int): (Optional) The number of lines to skip at the beginning of the file if `skip` is True

### Usage:
```python
countby_attribute(input_gtf, output_name, countby_value, skip, num)
```

## 4) gtf_modifiers.py - modify existing gtf files, i.e. changing keys, adding parent-child relationships

### `add_sequence_gtf`

This function adds the sequence information from a reference genome to a GTF file. It uses the `add_sequence_fast` function to extract the sequences and appends them as a new attribute with the specified attribute name in the output file.

### Parameters
- `gtf_file` (str): The path to the input GTF file.
- `ref_genome` (str): The path to the reference genome file in FASTA format.
- `output_name` (str): The desired name and path of the output file.
- `attribute_name` (str, optional): The name of the attribute to be added for the extracted sequences. Default is "bed_sequence".

### Usage
```python
add_sequence_gtf(gtf_file, ref_genome, output_name, attribute_name)
```

---

### `add_sequence_bed`

This function adds the sequence information from a reference genome to a BED file. It uses the `add_sequence_fast` function to extract the sequences and appends them as an **additional column** in the output file.

### Parameters
- `bed_file` (str): The path to the input BED file.
- `ref_genome` (str): The path to the reference genome file in FASTA format.
- `output_name` (str): The desired name and path of the output file.

### Usage
```python
add_sequence_bed(bed_file, ref_genome, output_name)
```

---

### `gtf_change_middle`

This function modifies the value of the "exon" attribute in a GTF file to a new specified value. It replaces the original value with the desired value.

### Parameters
 - `gtf_file` (str): The path of the input GTF file.
 - `out_name` (str): The desired output name for the file.
 - `change_value` (str): The value to replace the "exon" attribute column with. 

 ### Usage
 ```python
 gtf_change_middle(gtf_file, out_name, change_value)
 ```

---

### `filter_gtf`

Filters a GTF file based on a specified feature and value.

### Parameters:
- `input_gtf` (str): Path to the input GTF file.
- `output_name` (str): Path to the output file.
- `filter_by` (str): The feature to filter by.
- `value` (str): The value to filter for.
- `skip` (bool, optional): If `True`, skips a specified number of lines. Default is `False`.
- `num` (int, optional): The number of lines to skip. Default is `1`.

### Usage:
```python
filter_gtf(input_gtf, output_name, filter_by, value, skip=False, num=1)
```

---

### `select_column`

Selects rows from a GTF file where a specific column matches a given value.

### Parameters:
- `input_gtf` (str): Path to the input GTF file.
- `output_name` (str): Path to the output file.
- `col_number` (int): The index of the column to select.
- `value` (str): The value to match in the specified column.

### Usage:
```python
select_column(input_gtf, output_name, col_number, value)
```

---

### `filter_column`

Filters out rows from a GTF file where a specific column matches a given value.

#### Parameters:
- `input_gtf` (str): Path to the input GTF file.
- `output_name` (str): Path to the output file.
- `col_number` (int): The index of the column to filter.
- `value` (str): The value to filter out in the specified column.

#### Usage:
```python
filter_column(input_gtf, output_name, col_number, value)
```

---

### `merge_overlaps`

Merges overlapping transcripts in a GTF file.

#### Parameters:
- `gtf_file` (str): Path to the input GTF file.
- `output_tsv` (str): Path to the output TSV file.
- `output_gtf` (str): Path to the output GTF file.
- `offset` (int): Number of nucleotides to extend beyond a transcript to merge. (i.e. merge x and y if within `offset` n.t.)

#### Usage:
```python
merge_overlaps(gtf_file, output_tsv, output_gtf, offset)
```

---

### `key_biotype_gtf`

Assigns biotypes to transcripts in a GTF file based on a dictionary mapping.

### Parameters:
- `gtf_file` (str): Path to the input GTF file.
- `gtf_output` (str): Path to the output GTF file.
- `a_dict` (dict): Dictionary mapping data origins to biotypes.

### Usage:
```python
key_biotype_gtf(gtf_file, gtf_output, a_dict)
```

---

### `standardize_attributes`

Standardizes attributes in a GTF file based on a dictionary mapping.

### Parameters:
- `gtf_file` (str): Path to the input GTF file.
- `gtf_output` (str): Path to the output GTF file.
- `a_dict` (dict): Dictionary mapping attribute identifiers to new attribute names.

### Usage:
```python
standardize_attributes(gtf_file, gtf_output, a_dict)
```

## 5) gtf_groundtruth.py - generate a ground truth dataset of different levels of mismatching of sequences

### `generate_groundtruth`

Generates a ground truth database that represents how different transcripts will map to each other with varying mismatch levels.

### Parameters:

- `gtf` (str): Path to the GTF file containing transcript information.
- `ref` (str): Path to the reference genome file.
- `new_dir` (str): Path to the directory where the output files will be saved.
- `primary_key` (str): The primary key column name for the primary table.
- `information_dict` (dict): Dictionary specifying additional information to include in the primary table. The keys are the column names, and the values are lists in the format `[gen_or_attribute, column_index_or_name]`, where `gen_or_attribute` is either 0 (for general data) or 1 (for attribute data), and `column_index_or_name` is the index or name of the column. i.e. `{"biotype":[1, "biotype"], "source":[1, "databases"]}`
- `hisat` (bool, optional): Whether to use HISAT2 for alignment. Defaults to False.
- `num_mismatch` (int, optional): The maximum number of mismatches allowed during alignment. Defaults to 2.  **Each increase will make this run much slower**
- `input_fasta` (str, optional): Path to the input FASTA file for alignment. Defaults to "sequences.fasta". Change to align with another fasta. Keys are what would be read names in the fasta. 
- `p` (int, optional): Number of cores to use.

### Usage:

```python
generate_groundtruth("transcripts.gtf", "reference_genome.fa", "output_dir", "transcript_id", {"gene_name": [1, "gene_name"], "gene_type":[0, "miRNA"]}, hisat=False, num_mismatch=3, input_fasta="input.fasta", p=2)
```

---

### `generate_fasta`

Generate a FASTA file from a primary_df with the following fields:
transcript_id
sequence
length

```python
generate_fasta(primary_df, new_dir)
```

- `primary_df` (str): DataFrame containing the primary table.
- `new_dir` (str): Path to the directory where output files should be saved.

## 6) conversion_tools.py - convert files to gtf, change index system

### `gff3_to_gtf`

Converts a GFF3 file to a GTF file.

### Parameters:
- `gff3_file` (str): Path to the input GFF3 file.
- `output_name` (str): Path to the output GTF file.

### Usage:
```python
gff3_to_gtf(gff3_file, output_name)
```

---

### `tsv_to_gtf`

Converts a TSV file to a GTF file.

### Parameters:
- `tsv` (str): Path to the input TSV file.
- `out_name` (str): Path to the output GTF file.
- `extract_dictionary` (dict): Dictionary specifying additional information to include in the primary table. The keys are the column names, and the values are lists in the format `[gen_or_attribute, column_index_or_name]`, where `gen_or_attribute` is either 0 (for general data) or 1 (for attribute data), and `column_index_or_name` is the index or name of the column. i.e. `{"biotype":[1, "biotype"], "source":[1, "databases"]}`
- `fill_dict` (bool): Whether to fill the extraction dictionary automatically based on the TSV file header. (default: `True`)
- `skip_lines` (int): Number of lines to skip at the beginning of the TSV file. (default: `1`)
- `header` (bool): Whether the TSV file has a header row. (default: `True`)
- `csv` (bool): Whether the TSV file is in CSV format. (default: `False`)

**Usage:**
```python
tsv_to_gtf(tsv, out_name, extract_dictionary=extract_dictionary, fill_dict=True, skip_lines=1, header=True, csv=False)
```

---

### `zero_to_one`

Converts a file with zero-based indices to one-based indices.

### Parameters:
- `gtf_file` (str): Path to the input file with zero-based indices.
- `out_name` (str): Path to the output file with one-based indices.

### Usage:
```python
zero_to_one(gtf_file, out_name)
```

---

### `one_to_zero`

Converts a file with one-based indices to zero-based indices.

### Parameters:
- `gtf_file` (str): Path to the input file with one-based indices.
- `out_name` (str): Path to the output file with zero-based indices.

### Usage:
```python
one_to_zero(gtf_file, out_name)
```

---

### `gtf_to_bed`

Converts a GTF file to BED format.

### Parameters:
- `gtf` (str): Path to the input GTF file.
- `outputname` (str): Path to the output BED file.
- `attribute_to_name` (str): Attribute name to be used as the BED name column. (default: `False`)

### Usage:
```python
gtf_to_bed(gtf, outputname, attribute_to_name=False)
```

---

### `bed_to_gtf`

Converts a BED file to GTF format.

### Parameters:
- `bed` (str): Path to the input BED file.
- `output` (str): Path to the output GTF file.
- `source` (str): Database source to be used in the GTF file.
- `biotype` (str): Biotype to be used in the GTF file.

### Usage:
```python
bed_to_gtf(bed, output, source, biotype)
```

---

### `gtf_to_fasta`

Converts a GTF file to FASTA format.

### Parameters:
- `gtf` (str): Path to the input GTF file.
- `output` (str): Path to the output FASTA file.
- `ref_genome` (str): Reference genome file.
- `primary_key` (str): A primary key in the attributes column of the GTF file to set as the read name in the fasta file.

### Usage:
```python
gtf_to_fasta(gtf, output, ref_genome, primary_key)
```
---

### `fasta_to_tsv`

Converts a fasta file to a two column tsv file with read name and sequence. 

### Parameters:
- `fasta` (str): Path to input fasta file.
- `output` (str): Path to output generated tsv file.
- `header_name` (str): Name of the non-sequence column.

## 7) gtf_verifiers.py
