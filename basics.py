import numpy as np
import pandas as pd
from pybedtools import BedTool
import os
import fire
from Bio.Seq import Seq
from collections import defaultdict
import time
from datetime import timedelta

# Descrition:
# Functions that are helpful for other scripts. 
# - Hamming Distance
# - Dictionary
# - DNA -> RNA

## -- Basic Functions -- ##
class my_dictionary(dict):
 
  # __init__ function
  def __init__(self):
    self = dict()
 
  # Function to add key:value
  def add(self, key, value):
    self[key] = value

# Calculate the hamming distance between two sequences
def hamming(seq1, seq2):
    i = 0
    counter = 0
    for letter in seq1:

        if letter != seq2[i]:
            counter += 1
        i += 1

    return counter

# Translate RNA if ever needed
def rna_trans(transcript):
    
    RNA = ""
    
    for letter in transcript:
        if letter == "T":
            RNA = RNA + "U"
        else:
            RNA = RNA + letter
    
    return RNA

def separate_gtf_line(line):
    split_line = line.split(sep = "\t")
    # attributes
    attributes = split_line[-1].split(sep = ";")
    # strip
    attributes = list(map(str.strip, attributes))
    # combine
    attributes_split = [item.split(sep = " ") for item in attributes]
    # combine
    attributes_split = sum(attributes_split, [])

    return_split_line = split_line[0:8]
    return [return_split_line, attributes_split]

def remove_comments(gtf,out):
    with open(gtf, "r") as gtf_f, open(out, "w") as new:
        for line in gtf_f:
            if line[0] == "#":
                continue
            else:
                new.write(line)

def rev_comp(seq):
    rev = seq[::-1]

    seq_obj = Seq(rev)

    comp = seq_obj.complement()

    return str(comp)

# create a timing class
class timing(object):
    # do init
    def __init__(self):
        self.timing = defaultdict(list)
        self.start_time = time.monotonic()

    def take_time(self, module_name):
        self.timing["module"].append(module_name)
        self.timing["timing"].append(timedelta(seconds=time.monotonic() - self.start_time))

    def export(self, run_name):
        pd.DataFrame.from_dict(self.timing).to_csv("timing_" + run_name + ".csv", index = False)

## -- Basic Functions -- ##

if __name__ == '__main__':
  fire.Fire()