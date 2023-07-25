import pandas as pd
import os
import re
from Bio.Seq import Seq

files = os.listdir()

csvs = []

rev_tab = {"A":"T", "T":"A", "C":"G", "G":"C"}

def rev_comp(seq):
    rev = seq[::-1]

    seq_obj = Seq(rev)

    comp = seq_obj.complement()

    return str(comp)


with open("samples_adapters.txt", "r") as samples:
    for line in samples:
        line_s = line.replace("\n", "")

        csv = line_s + ".tsv"

        csvs.append(csv)

i = 1
for csv in csvs:
    cur = pd.read_csv(csv, sep = "\t")
    
    if i == 1:
        combined = cur
        i += 1

    else:
        combined = combined.merge(cur, how = "outer", on = "sequence")

combined = combined.fillna(0)

combined[['flag','sequence']] = combined["sequence"].str.split('_', expand = True)

combined = combined.loc[combined['flag'] != "16", :]

lookup_table = pd.read_csv("sRNA_frag_lookup.csv")

new_table = combined.merge(lookup_table, how = "left", on = "sequence")

new_table['length'] = new_table["sequence"].str.len()

new_table.to_csv("sRNA_frag_counts.csv", index = False)
