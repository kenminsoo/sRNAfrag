from basics import *

# Description:
# Tools for working with chromosome aliases in references and annotations.
# - Combine reference genomes (Missing alias/haplotype)
#    - Not useful, might pivot away from using alt. haplotypes
# - Stadardize names between different aliases 
#   - From different sources

# Combine reference genomes
# documented
def ref_combine(fa1_name, fa2_name, out_name, reference):
    
    reference_nomen = pd.read_csv(reference, sep = "\t")
    reference_nomen["Source"] = ""
    reference_nomen["Exist"] = False

    add_seq = False
    number_columns = len(reference_nomen.columns)
    
    with open(fa1_name, "r") as fa1, open(fa2_name, "r") as fa2, open(out_name, "w") as new:
        for line in fa1:
            # Identify if marks chromosome
            if line[0:1] == ">":
                # See what column it is in

                splt = line[1:].split()[0]

                for column in reference_nomen:
                    iters = 0
                    if splt in list(reference_nomen[column]):
                        # Find where chr is in the dataframe
                        index_num = reference_nomen.index[reference_nomen[column] == splt].to_list()[0]
                        
                        # Mark as existing
                        if reference_nomen.iloc[index_num, -1] == False:
                            # Set it to true if it was false
                            reference_nomen.iloc[index_num, -1] = True
                            reference_nomen.iloc[index_num, -2] = fa1_name
                            add_seq = True
                            new.write(">" + str(reference_nomen.iloc[index_num, 0] + "\n"))
                            break
                        
                        else:
                            # if it is true we do not add sequence
                            add_seq = False
                            break

                    else:
                        print("not in this " + column)
                        iters += 1

                        # If it's never found, ignore
                        # Perhaps we can skip and add to summary report after?
                        if iters == number_columns:
                            raise ValueError("Chromosome Not found")
            else:
                if add_seq == True:
                    new.write(line)
                else:
                    continue

        for line in fa2:
            # Identify if marks chromosome
            if line[0:1] == ">":
                # See what column it is in

                splt = line[1:].split()[0]

                for column in reference_nomen:
                    iters = 0
                    if splt in list(reference_nomen[column]):
                        # Find where chr is in the dataframe
                        index_num = reference_nomen.index[reference_nomen[column] == splt].to_list()[0]
                        
                        # Mark as existing
                        if reference_nomen.iloc[index_num, -1] == False:
                            # Set it to true if it was false
                            reference_nomen.iloc[index_num, -1] = True
                            reference_nomen.iloc[index_num, -2] = fa2_name
                            add_seq = True
                            new.write(">" + str(reference_nomen.iloc[index_num, 0] + "\n"))
                            break
                        
                        else:
                            # if it is true we do not add sequence
                            add_seq = False
                            break

                    else:
                        print("not in this " + column)
                        iters += 1

                        # If it's never found, ignore
                        # Perhaps we can skip and add to summary report after?
                        if iters == number_columns:
                            raise ValueError("Chromosome Not found")

            else:
                if add_seq == True:
                    new.write(line)
                else:
                    continue

    # Extract the informational file
    reference_nomen.to_csv("merge_info.csv")

# Standardize gtf chr names
# When different chr aliases are used
def gtf_naming_stan(gtf_file, reference, out_name):
    reference_nomen = pd.read_csv(reference, sep = "\t")

    num_columns = len(reference_nomen.columns)

    with open(gtf_file, "r") as gtf, open(out_name, "w") as new:
        
        for line in gtf:
            
            modify_line = line.split(sep = "\t")

            chromosome_name = modify_line[0]

            chr_not_found = True

            while chr_not_found == True:
                for column in reference_nomen:
                    if chromosome_name in list(reference_nomen[column]):
                        index_num = reference_nomen.index[reference_nomen[column] == chromosome_name].to_list()[0]
                        chr_not_found = False
                        break
                    else:
                        continue

                if chr_not_found == False:
                    break
                elif chr_not_found == True:
                    raise ValueError(chromosome_name + " not in alaias!")

            modify_line[0] = reference_nomen.iloc[index_num, 0]
            new.write("\t".join(modify_line))

# takes in a reference genome and chromosome list
# and filters them
chr = ['chr1', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrM', 'chrX', 'chrY']

def ref_chr_select(ref_genome, out_name,chr_list = chr):
    switch = False
    with open(ref_genome, "r") as ref, open(out_name, "w") as new:
        for line in ref:
            if line[0] == ">":
                if line[1:].strip() in chr_list:
                    new.write(line)
                    switch = True
                else:
                    switch = False
            elif switch == True:
                new.write(line)
            else:
                continue

# takes in a gtf and chromosome list
# and filters them
def gtf_chr_select(gtf, out_name,chr_list = chr):
    switch = False
    with open(gtf, "r") as gtf, open(out_name, "w") as new:
        for line in gtf:
            split = line.split(sep = "\t")
            chr = split[0]

            if chr in chr_list:
                new.write(line)
            else:
                continue

# Extracting chromosome names 
def bam_chr_extract(bam, out):
    os.system("samtools view -H " + bam + " | grep @SQ | cut -f 2 | sed 's/SN://g' > " + out)

    chroms = []

    with open(out, "r") as chrs:
        for line in chrs:
            chroms.append(line)

    return chroms

# if used with fire, will print out a list
def fasta_chr_extract(fasta):
    chrs = []
    
    with open(fasta, "r") as seqs:
        for line in seqs:
            if line[0] == ">":
                sep_line = line.split(sep = " ")

                chr_name = sep_line[0]

                chrs.append(chr_name.strip(">").strip("\n"))

    print(chrs)

if __name__ == '__main__':
  fire.Fire()