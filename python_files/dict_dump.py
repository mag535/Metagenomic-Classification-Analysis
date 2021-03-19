# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 12:39:40 2020

@author: Melissa Gray
"""
#%%
import pickle

ID_dic = {}

#%%

if __name__ == "__main__":
    dic = open("/mnt/HA/groups/rosenGrp/mag535/kraken_misc_files/nz_gcf_completed_dict.txt", "r", encoding="utf-8")
    
    for line in dic:
        vals = line.split(" : ")
        ID_dic[vals[0].strip()] = vals[1].strip()
        
    dic.close()
    pickle.dump(ID_dic, open("/lustre/scratch/mag535/nz_gcf.p", "wb"))