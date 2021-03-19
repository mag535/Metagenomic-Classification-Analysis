# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 12:51:51 2020

@author: Melissa Gray
"""

#%% VAIRABLES & IMPORTS
from Bio import Entrez
import re
from time import sleep

Entrez.email = "mag535@drexel.edu"
### This is the accession number start with NZ
# REF_SEQ_ID = "NZ_CP046528.1" 


#%% FUNCTIONS

def gb_id_map(REF_SEQ_ID):
    handle = Entrez.efetch(db="nucleotide", id=REF_SEQ_ID, retmode="xml")
    #print(handle.read())
    Assembly_ID = None
    res = handle.read()
    
    no_line = open("no_line.txt", "w+")
    many_lines = open("many_lines.txt", "w+")
    
    found = re.findall(".*GCF_.*", res)
    
    if len(found) == 0:
        no_line.write(REF_SEQ_ID + "\n")
    elif len(found) == 1:
        Assembly_ID = re.split("[<>]", found[0])[2]
    elif len(found) > 1:
        many_lines.write(REF_SEQ_ID + "\n")
        for item in found:
            many_lines.write("\t" + item + "\n")
    
    no_line.close()
    many_lines.close()
    return REF_SEQ_ID, Assembly_ID

#%% MAIN

if __name__ == '__main__':
    no_match = open("From_Proteus/no_matches.txt", "r", encoding='utf-8')
    test = open("test_id_match.txt", "w+")
    
    for line in no_match:
        parts = line.split("|")
        nz = parts[-1].split(" ")[0]
        ref_seq_id, a_id = gb_id_map(nz)
        test.write(ref_seq_id + " : " + str(a_id) + "\n")
        sleep(0.35)
    
    no_match.close()
    test.close()
    
