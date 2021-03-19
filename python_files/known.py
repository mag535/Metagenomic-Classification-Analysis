# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 16:43:01 2020

@author: Melissa Gray
"""

## Array of Known

#%% VARIABLES & IMPORTS

import numpy as np
import pandas as pd
import pickle
import sys
import os.path
#import cross_validation.py as cv


#%% FUNCTIONS FROM OTHER PY FILES

def count_read_match(classify_res, output_path):
    nodes_dic = pickle.load(open("/mnt/HA/groups/rosenGrp/mag535/nodes_dic.p", 'rb'))
    count_dic = {}
    with open(classify_res, 'r') as cr:
        for line in cr:
            status = line.split('\t')[0].strip()
            # CHECK UNCLASSIFIED
            if status == 'U':
                if "Unclassified" not in count_dic:
                    count_dic["Unclassified"] = 1
                else:
                    count_dic["Unclassified"] += 1
            
            # CHECK CLASSIFIED
            else:
                splt_line = line.split("\t")
                real = int(splt_line[1].strip().split(":")[-1])
                pred = splt_line[2].strip()
                
                # TRACEBACK
                res = trace_back_taxa(pred, nodes_dic)
                no_species = True
                for tid, rank in res:
                    if rank == 'species':
                        no_species = False
                        pred = int(tid)
                        break
                
                # LOG HOW MANY CAN'T BE TRACED BACK TO SPECIES LEVEL
                if no_species:
                    if "No_Species" not in count_dic:
                        count_dic["No_Species"] = 1
                    else:
                        count_dic["No_Species"] += 1
                    continue
                
                # CHECK IF MATCHED
                if pred == real:
                    if "Matched" not in count_dic:
                        count_dic['Matched'] = 1
                    else:
                        count_dic["Matched"] += 1
                else:
                    if "Incorrect" not in count_dic:
                        count_dic["Incorrect"] = 1
                    else:
                        count_dic["Incorrect"] += 1
    return count_dic


#%% FUNCTIONS

def trace_back_taxa(taxid, nodes_dic):
    if taxid not in nodes_dic:
        #### write to a file
        with open("no_traceback.txt", 'w') as file:
            file.write(taxid + '\n')
        return []
    res = []
    cur_id = taxid
    while cur_id != "1":
        parent_id, rank = nodes_dic[cur_id]
        res.append((cur_id, rank))
        cur_id = parent_id
    return res

def get_known_species(library_headers, nodes_dic):
    known_species = set()
    with open(library_headers, 'r') as headers:
        for line in headers:
            taxid = line.split("|")[1]
            # Traceback to species level and add to known_species
            res = trace_back_taxa(taxid, nodes_dic)
            for tid, rank in res:
                if rank == "species":
                    known_species.add(int(tid))
    return known_species

def array_of_known(classify_results, library_headers, nodes_dic):
    array_of_known = []
    
    known_species = get_known_species(library_headers, nodes_dic)
    print("Known Species ->")
    for species in known_species:
        print(species)
    print()
    
    uncl = 0
    no_spec = 0
    others = 0
    
    with open(classify_results, 'r') as cr:
        for line in cr:
            # Set Up and append True Species
            line_array = []
            splt_line = line.split("\t")
            true_species = int(splt_line[1].strip().split(":")[-1])
            line_array.append(true_species)
            
            # Check if Pred_Species is Unclassified
            status = splt_line[0].strip()
            if status == 'U':
                line_array.append(-2)
                uncl += 1
            else:
                # Check if Pred_Species can be traced back to species level
                pred_species = splt_line[2].strip()
                res = trace_back_taxa(pred_species, nodes_dic)
                no_species = True
                for tid, rank in res:
                    if rank == 'species':
                        no_species = False
                        pred_species = int(tid)
                        others += 1
                        break
                # Append -1 for No_Species or the Species level taxid found after traceback
                if no_species:
                    line_array.append(-1)
                    no_spec += 1
                else:
                    line_array.append(pred_species)
            
            # Check if Pred_Species is in Library header file of Kraken Database
            if true_species in known_species:
                line_array.append(1)
            else:
                line_array.append(0)
            
            # Append line array to known_array
            array_of_known.append(line_array)
    return array_of_known, uncl, no_spec, others

def calculate_match_accuracy(array_of_known, classify_results, output_path):
    # import cv.count_read_match(classify_res, output_path)
    # total up matches by accessing count_dic["Matched"]
    ### total_matches = count_read_match(classify_results, output_path)["Matched"]
    # total up known by counting numbers of 1's in column three of array_of_known
    total_matches = 0
    total_known = 0
    total_known_matched = 0
    total_known_classified = 0
    unknown_to_model = 0
    for row in array_of_known:
        if row[1] == row[0]:
            total_matches += 1
        if row[2] == 1:
            total_known += 1
            if row[1] == row[0]:
                total_known_matched += 1
            if row[1] != -2:
                total_known_classified += 1
        elif (row[2] == 0) and (row[1] == -2):
            unknown_to_model += 1
    # total reads by counting number of lines in classify_results
    total_reads = 0
    with open(classify_results, "r") as file:
        for line in file:
            total_reads += 1
    return total_matches, total_known, total_known_matched, total_known_classified, total_reads, unknown_to_model


#%%% MAIN

if __name__ == '__main__':
    cls_res, lib_headers, op, year, fold = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    
    nodes_dic = pickle.load(open("/mnt/HA/groups/rosenGrp/mag535/nodes_dic.p", 'rb'))
    aok, u, ns, o = array_of_known(cls_res, lib_headers, nodes_dic)
    print("Year :", year)
    print("Classification Result :", cls_res)
    print("Unclassified :", u)
    print("No Speices :", ns)
    print("Others :", o)
    
    tm, tk, tkm, tkc, tr, utm = calculate_match_accuracy(aok, cls_res, op)
    print("Total Correct Matches :", tm)
    print("Total Known Matches :", tkm)
    print("Total Known :", tk)
    print("Total Known Classified :", tkc)
    print("Total Unknown To Model :", utm)
    print("Total Reads :", tr)
    print("Known Accuracy (Total Known Matches / Total Known) :", (tkm / tk))
    print("Known Classified Accuracy (Total Correct Matches / Total Known Classified) :", (tm / tkc))
    print("Raw Accuracy (Total Correct Matches / Total Reads) :", (tm / tr))
    print("Real Unknown (Total Unknown To Model / Unclassified) :", (utm / u))
    
    aok_array = np.array(aok)
    aok_df = pd.DataFrame(aok_array, columns=['True Species', 'Predicted Species', 'In Training Data'])
    aok_df.to_csv(os.path.join(op, "{}_fold{}_array_of_known_species.csv".format(year, fold)))
    
    
    