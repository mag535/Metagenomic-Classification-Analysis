# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 10:39:45 2020

@author: Melissa Gray
"""

#%% IMPORTS & VARIABLES

import sys
import os.path
import glob
import pickle
import numpy as np


#%% FUNCTIONS & CLASSES

def create_nz_map(library_headers, output_path):
    # output_path should be /lustre/scratch/mag535/cross_validation_files/
    nz_map = {}
    with open(library_headers, 'r') as headers:
        for line in headers:
            nz_accession = line.split(' ')[0].split('|')[2]
            nz_map[nz_accession] = line
    pickle.dump(nz_map, open(os.path.join(output_path, 'nz_map.p'), 'wb'))
    return

def get_file_list(root_dir):
    # root_dir should look like this: /lustre/scratch/mag535/year_exp/{year}/{fold}
    file_list = []
    cls_list = glob.glob('{}/*'.format(root_dir))
    for cls_dir in cls_list:
        seq_list = glob.glob('{}/*.fna'.format(cls_dir))
        for seq_path in seq_list:
            file_list.append(seq_path)
    return file_list

def read_in_file(file_name, nz_map, fw):
    flag = False
    abnormal = 0
    with open(file_name, 'r') as f:
        for line in f:
            if line[0] == '>':
                nz_accession = line[1:].split(' ')[0]
                #print(nz_accession)
                if nz_accession not in nz_map:
                    #return False
                    abnormal += 1
                    flag = False
                else:
                    flag = True
                    fw.write(nz_map[nz_accession])
            elif flag:
                fw.write(line)
    return abnormal

def create_library_file(year, fold, cored_dir, library_headers, output_path):
    # cored_dir should look like this: '/lustre/scratch/mag535/year_exp/{}/{}'
    try:
        nz_map = pickle.load(open(os.path.join(output_path, 'nz_map.p'), 'rb'))
    except FileNotFoundError:
        create_nz_map(library_headers, output_path)
        nz_map = pickle.load(open(os.path.join(output_path, 'nz_map.p'), 'rb'))
    
    root_dir = cored_dir.format(year, fold)
    file_list = get_file_list(root_dir)
    output_f = os.path.join(output_path, '{}_kraken/training_data_{}_{}.fna'.format(year, year, fold))
    fw = open(output_f, 'w+')
    abnormal_counter = 0
    for file_name in file_list:
        tmp_abnormal = read_in_file(file_name, nz_map, fw)
        abnormal_counter += tmp_abnormal
    print('Abnormal Counter:', abnormal_counter)
    print('Length of File List:', len(file_list))
    return

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

#%% MAIN

if __name__ == '__main__':
    '''
    year, fold, cored_dir, library_headers, output_path = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
    create_library_file(year, fold, cored_dir, library_headers, output_path)
    '''
    
    classify_results, output_path = sys.argv[1], sys.argv[2]
    cd = count_read_match(classify_results, output_path)
    for t in cd:
        print("{} : {}".format(t, cd[t]))
    