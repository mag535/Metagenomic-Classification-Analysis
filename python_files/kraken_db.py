# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 10:00:51 2020

@author: Melissa Gray
"""

#%% IMPORTS & VARIABLES

import sys
import os.path
import random
import pickle
import numpy as np


nodes_dic = {}
names_dic = {}

#%% FUNCTIONS & CLASSES

def create_prelim(library, output_path):
    lib = open(library, "r", encoding='utf-8')
    prelim = open(os.path.join(output_path, "prelim_map.txt"), "w+")
    
    for line in lib:
        if line[0] == ">":
            accession = line.split("|")[-1].split(" ")[0]
            taxid = line.split("|")[1]
            prelim.write('TAXID\tkraken:taxid|' + taxid +"|"+ accession +"|\t" + taxid + "\n")
    
    lib.close()
    prelim.close()
    return

def create_header_file(library, output_path):
    lib = open(library, "r", encoding='utf-8')
    headers = open(os.path.join(output_path, "library_headers.fna"), "w+")
    
    for line in lib:
        if line[0] == ">":
            headers.write(line)
            
    
    lib.close()
    headers.close()
    return

def create_fa(library, output_path):
    lib = open(library, "r", encoding='utf-8')
    content = lib.read()
    lib.close()
    
    max_rand = len(content)
    nucleotides = ["G", "C", "T", "A", "x"]
    seq = open(os.path.join(output_path, "generated_test.fa"), "w+")
    
    i=0
    while i < 1000:
        i += 1      # continues outer loop
        
        #GET THE READS
        reads = ''
        r = random.randint(0, max_rand)
        while content[r] not in nucleotides:
            r = random.randint(0, max_rand)
        start = r
        
        reads += content[r:r+210]
        
        problem = False
        new_lines = 0
        for letter in reads:
            if letter == "\n":
                new_lines += 1
            elif letter not in nucleotides:
                i -= 1                  # redo this run in outer loop
                problem = True              # signal to skip finding taxid
                break
        '''
        if problem:
            continue
        else:
            #reads.replace("\n", "")
            j=0
            while j < new_lines:
                j += 1
                if content[r+210+j] == "\n":
                    j -= 1
                elif content[r+210+j] in nucleotides:
                    reads += content[r+210+j]
                else:
                    i -= 1
                    problem = True
                    break
        '''
        if problem:
            continue
        
        #GET THE TAXID
        found = False
        taxid = ''
        b = -1
        for k in range(start, -1, -1):
            if content[k] == "_":
                b=k
                while content[b] != "|":
                    b -= 1
                taxid += content[b-1]
                found = True
            
            if found:
                h = b-2
                while content[h] != "|":
                    taxid += content[h]
                    h -= 1
                break
        
        # WRITE
        seq.write(">test_" + taxid[::-1] + "\n")
        seq.write(reads + "\n")
    
    seq.close()
    return

def create_manifest(db_path):
    prelim = open(os.path.join(db_path, "prelim_map.txt"), "r", encoding='utf-8')
    old = open("/lustre/scratch/mag535/Kraken/bac_kraken_db/library/bacteria/manifest.txt", "r", encoding='utf-8')
    
    nz_dict = pickle.load(open("/mnt/HA/groups/rosenGrp/mag535/kraken_misc_files/nz_gcf.p", "rb"))
    mani_lines = []
    
    for line in prelim:
        nz = line.split("|")[-1].split(" ")[0]
        gcf = nz_dict[nz]
        for line in old:
            if gcf in line:
                mani_lines.append(line)
                break
    
    prelim.close()
    old.close()
    mani = open(os.path.join(db_path, "manifest.txt"), "w+")
    
    for line in mani_lines:
        mani.write(line)
    
    mani.close()
    return

def count_taxid(fasta):
    taxid_count = {}
    taxids = set()
    
    with open(fasta, "r") as file:
        for line in file:
            parse = line.split("\t")
            if parse[0] == "U":
                if "Unclassified" not in taxids:
                    taxids.add("Unclassified")
                    taxid_count["Unclassified"] = 0
                taxid_count["Unclassified"] += 1
            elif parse[0] == "C":
                if parse[2] not in taxids:
                    taxids.add(parse[2])
                    taxid_count[parse[2]] = 0
                taxid_count[parse[2]] += 1
        
    if "Unclassified" not in taxids:
        taxids.add("Unclassified")
        taxid_count["Unclassified"] = 0
    
    pickle.dump(taxid_count, open("/lustre/scratch/mag535/taxid_count_1999.p", "wb"))
    
    total = np.sum([val for val in taxid_count.values()])
    print("Year: 1999")
    print("SUM:", total)
    print("Expected Sum: 74796007")
    print("These sums are equal:", (total==74796007))
    
    '''
    with open("taxid_count_1999.txt", "w") as output:
        for taxid in taxid_count.keys():
            output.write(str(taxid) + " : " + str(taxid_count[taxid]) + "\n")
    '''
    return

def create_recursive_dictionaries(names_dmp, nodes_dmp):
    names_dic = {}          # saves as {taxid : scientific name} dictionary
    names_log = set()       # saves taxids that have more than one scientific name (or else the previous scientific names would be overwritten)
    nodes_dic = {}          # saves as a {taxid : (parent, rank)} dictionary, where rank is for taxid, not parent
    
    with open(names_dmp, "r") as name_file:
        for line in name_file:
            if "scientific name" in line:
                parse = [info.strip() for info in line.split("|")]
                if parse[0] in names_dic.keys():
                    names_log.add((parse[0], names_dic[parse[0]]))
                    names_log.add((parse[0], parse[1]))
                names_dic[parse[0]] = parse[1]
    pickle.dump(names_dic, open("names_dic.p", "wb"))
    
    with open(nodes_dmp, "r") as node_file:
        for line in node_file:
            if line[0] != "0":
                parse = [info.strip() for info in line.split("|")]
                nodes_dic[str(parse[0])] = (parse[1], parse[2])
    pickle.dump(nodes_dic, open("nodes_dic.p", "wb"))
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

def taxa_frequency(classification_file, level):
    frequency = {}
    #levels = ['strain', 'subspecies', 'species', 'subgenus', 'genus', 'subfamily', 'family', 'suborder', 'order', 'subclass', 'class', 'subphylum', 'subdivision', 'phylum', 'division', 'superkingdom']
    file = open(classification_file, "r", encoding='utf-8')
    
    for line in file:
        parse = [info.strip() for info in line.split("\t")]
        if parse[0] == "U":
            if "Unclassified" not in frequency:
                frequency["Unclassified"] = {}
            if parse[2] not in frequency["Unclassified"]:
                frequency["Unclassified"][parse[2]] = 0
            frequency["Unclassified"][parse[2]] += 1
        else:
            res = trace_back_taxa(parse[2], nodes_dic)
            for p in res:
                if p[1] not in frequency:
                    frequency[p[1]] = {}
                if p[0] not in frequency[p[1]]:
                    frequency[p[1]][p[0]] = 0
                frequency[p[1]][p[0]] += 1
    
    file.close()
    return frequency

def print_abundance(year, level, input_path, output_path):
    tmp = pickle.load(open(os.path.join(input_path, "{}_taxa_frequency.p".format(year)), 'rb'))
    names_dic = pickle.load(open(os.path.join(input_path, "names_dic.p"), "rb"))
    total = 0
    with open(os.path.join(output_path, "{}_{}.txt".format(year, level)), 'w') as file:
        for p in tmp[level]:
            total += tmp[level][p]
            file.write(p + ": " + str(tmp[level][p]) + ", " + names_dic[p] + '\n')
        file.write("Unclassified: " + str(tmp['Unclassified']['0']) + '\n')
        total += tmp['Unclassified']['0']
        file.write("TOTAL: " + str(total) + '\n')
    return

def find_species_ID(headers_file):
    species_taxid = []
    nodes_dic = pickle.load(open("/mnt/HA/groups/rosenGrp/mag535/nodes_dic.p", "rb"))
    with open(headers_file, 'r') as headers:
        for line in headers:
            no_species = True
            cur_taxid = line.split("|")[1]
            res = trace_back_taxa(cur_taxid, nodes_dic)
            
            for taxid, rank in res:
                if rank == 'species':
                    no_species = False
                    species_taxid.append(taxid)
                    break
            
            if no_species:
                species_taxid.append(-1)
    return species_taxid

def create_testing_data(PARENT_DIR, YEAR, FOLD):
    output = '{}/{}_kraken/combined_data_{}_fold{}.fna'.format(PARENT_DIR, YEAR, YEAR, FOLD)
    fw = open(output, 'w+')
    
    for year in range(1999, YEAR+1):
        if not os.path.isdir(os.path.join(PARENT_DIR, "{}_kraken".format(year))):
            continue
        for fold_idx in range(1, 6):
            if fold_idx == FOLD:
                continue
            file = '{}/{}_kraken/training_data_{}_fold{}.fna'.format(PARENT_DIR, year, year, fold_idx)
            with open(file, 'r') as f:
                fw.write(f.read())
    return


#%% MAIN

if __name__ == "__main__":
    PARENT_DIR, YEAR, FOLD = sys.argv[1], sys.argv[2], sys.argv[3]
    
    #create_prelim(library, output_path)
    
    #create_fa(library, output_path)
    
    #count_taxid(fasta)
    
    '''
    # only need to run once; if file does not exist, then run to create
    file_path, starter_taxid = sys.argv[1], sys.argv[2]
    test_recursive_parent_search(file_path, starter_taxid)
    '''
    '''
    classification_file, output_name, level = sys.argv[1], sys.argv[2], sys.argv[3]
    nodes_dic = pickle.load(open("nodes_dic.p", "rb"))
    frequency = taxa_frequency(classification_file, level)
    pickle.dump(frequency, open(output_name, "wb"))
    '''
    #print_abundance(2020, 'phylum', "C:\\Users\\milkg\\OneDrive\\Desktop\\From_Proteus", "C:\\Users\\milkg\\OneDrive\\Desktop\\From_Proteus")
    
    #species_list = find_species_ID(headers_file)
    
    create_testing_data(PARENT_DIR, int(YEAR), int(FOLD))
    
    