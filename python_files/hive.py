# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 11:51:32 2020

@author: Melissa Gray
"""

#%% VARS & IMPORTS
from summary import Summary
import os.path
import sys


#%% MAIN

if __name__ == '__main__':
    path = "/lustre/scratch/mag535/Kraken/bac_kraken_db/library/bacteria"
    Sum = Summary(os.path.join(path, "assembly_summary.txt"))
    
    #inputs
    '''
    input_years = (sys.argv[1], sys.argv[2])
    output_dir = sys.argv[3]
    '''
    time_range = (sys.argv[1], sys.argv[2])
    
    #dictionary creation and check
    '''
    Sum.dictionary_ID()
    Sum.check_final_dictionary("latest", 'Complete Genome')
    '''
    
    # sub_library.fna creation
    '''
    ids = Sum.find_kraken_IDs(input_years)
    name = os.path.join("/mnt/HA/groups/rosenGrp/mag535/{}".format(output_dir, "library.fna")
    Sum.gather_data(ids, name)
    '''
    
    # still check data inside gather_data()
    year_dir = time_range[0].split("/")[-1]
    Sum.check_data(time_range, year_dir)
    