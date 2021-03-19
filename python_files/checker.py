# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 10:36:49 2020

@author: Melissa Gray
"""

#%%% VARIABLES & IPMORTS
import sys
import os.path

#%% FUNCTIONS & CLASSES

def check_header_ending(header_file):
    info_dic = {"complete genome" : 0, "complete sequence" : 0, "neither" : 0, "no comma" : 0}
    line_dic = {"complete genome" : list(), "complete sequence" : list(), "neither" : list(), "no comma" : list()}
    
    file = open(header_file, "r")
    
    for line in file:
        if  "," in line:
            if "complete genome" in line.split(",")[-1]:
                info_dic["complete genome"] += 1
                line_dic["complete genome"].append(line)
            elif "complete sequence" in line.split(",")[-1]:
                info_dic["complete sequence"] += 1
                line_dic["complete sequence"].append(line)
            else:
                info_dic["neither"] += 1
                line_dic["neither"].append(line)
        else:
            info_dic["no comma"] += 1
            line_dic["no comma"].append(line)
    
    file.close()
    return info_dic, line_dic

def find_header_key_words(header_file):
    key_nums = {"chromosome" : 0, "plasmid" : 0, "neither" : 0}
    key_lines = {"chromosome" : list(), "plasmid" : list(), "neither" : list()}
    
    file = open(header_file, "r")
    
    for line in file:
        if "chromosome" in line:
            key_nums["chromosome"] += 1
            key_lines["chromosome"].append(line)
        elif "plasmid" in line:
            key_nums["plasmid"] += 1
            key_lines["plasmid"].append(line)
        else:
            key_nums["neither"] += 1
            key_lines["neither"].append(line)
    
    file.close()
    return key_nums, key_lines

def do_endings(header_file):
    checked, sorted_lines = check_header_ending(header_file)
    
    output = open("/lustre/scratch/mag535/output_header_check.txt", "w+")
    for key in checked:
        output.write(key + " : " + str(checked[key]) + "\n")
    output.close()
    
    output2  = open("/lustre/scratch/mag535/genome_headers.txt", "w+")
    for line in sorted_lines['complete genome']:
        output2.write(line)
    output2.close()
    
    output3  = open("/lustre/scratch/mag535/sequence_headers.txt", "w+")
    for line in sorted_lines['complete sequence']:
        output3.write(line)
    output3.close()
    
    output4  = open("/lustre/scratch/mag535/neither_headers.txt", "w+")
    for line in sorted_lines['neither']:
        output4.write(line)
    output4.close()
    
    output5  = open("/lustre/scratch/mag535/no_comma_headers.txt", "w+")
    for line in sorted_lines['no comma']:
        output5.write(line)
    output5.close()
    return

def do_key_words(header_file):
    nums, kl = find_header_key_words(header_file)
    
    output = open("/lustre/scratch/mag535/output_key_header_check.txt", "w+")
    for key in nums:
        output.write(key + " : " + str(nums[key]) + "\n")
    output.close()
    
    output2 = open("/lustre/scratch/mag535/chromosome_headers.txt", "w+")
    for line in kl["chromosome"]:
        output2.write(line)
    output2.close()
    
    output3 = open("/lustre/scratch/mag535/plasmid_headers.txt", "w+")
    for line in kl["plasmid"]:
        output3.write(line)
    output3.close()
    
    output4 = open("/lustre/scratch/mag535/neither_headers.txt", "w+")
    for line in kl["neither"]:
        output4.write(line)
    output4.close()
    return

def check_genomes(year_dir, num):
    file = open("/mnt/HA/groups/rosenGrp/mag535/{}/library.fna".format(year_dir), "r", encoding='utf-8')
    temp = []
    
    counter = 0
    for line in file:
        if line[0] == ">":
            counter += 1
            if counter == num:
                temp.append(line)
        else:
            if counter == num:
                temp.append(line)
            elif counter > num:
                break
    
    file.close()
    origin = open("/lustre/scratch/mag535/Kraken/bac_kraken_db/library/bacteria/library.fna", "r", encoding='utf-8')
    
    flag = False
    checks = 0
    stop = len(temp)
    
    for line in origin:
        if line[0] == ">" and line == temp[0]:
            print("The headers match.")
            checks += 1
            flag = True
        else:
            if flag:
                if line == temp[checks]:
                    checks += 1
                    if checks == stop:
                        print("The genomes match")
                        break
    
    origin.close()
    return


#%% MAIN

if __name__ == "__main__":
    year_dir, num = sys.argv[1], sys.argv[2]
    
    check_genomes(year_dir, int(num))
    
    