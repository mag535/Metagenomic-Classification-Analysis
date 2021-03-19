# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 13:21:57 2020

@author: Melissa Gray
"""


#%% VARIABLES & IMPORTS

import sys
import os.path
import glob
import pandas as pd

...

#%% FUNCTIONS

def find_date(time_list=("-1/-1/-1", "-1/-1/-1"), name=''): #end_year=-1, start_year=-1):
    # Gathers all the records from during and/or before the year specified
    ref_df = pd.read_csv('/ifs/groups/rosenGrp/mag535/python_files/reference')
    ref_df['Seq_rel_date'] = pd.to_datetime(pd.read_csv('/ifs/groups/rosenGrp/mag535/python_files/reference')['Seq_rel_date'])
    
    if all(time_list) == "-1/-1/-1":
        # no dates
        time_df = ref_df
    elif time_list[1] == "-1/-1/-1" and time_list[0] != "-1/-1/-1":
        # End date, no start date
        # before and including that time
        time_df = ref_df[ref_df['Seq_rel_date'] <= pd.to_datetime(time_list[0])]
    elif time_list[1] != "-1/-1/-1" and time_list[0] == "-1/-1/-1":
        # Start date, no end date
        # that time and after
        time_df = ref_df[ref_df['Seq_rel_date'] >= pd.to_datetime(time_list[1])]
    elif time_list[1] != "-1/-1/-1" and time_list[0] != "-1/-1/-1":
        # Date range
        # get time after and including start date
        temp_df = ref_df[ref_df['Seq_rel_date'] >= pd.to_datetime(time_list[1])]
        # get time before (not including) end date
        time_df = temp_df[temp_df['Seq_rel_date'] < pd.to_datetime(time_list[0])]
    
    if name != '':
        file_name = os.path.join('/ifs/groups/rosenGrp/mag535/python_files', name + '.csv')
        time_df.sort_values(by='Seq_rel_date').to_csv(file_name, index=False, header=True)
        print("\'{}.csv\' is saved.".format(file_name))
    return time_df

def find_kraken_IDs(times):
    tmp_df = find_date(times)
    gcf_list = tmp_df['Assembly_accession']
    return gcf_list

def find_gcf_files(gcf_list, clark_db):
    genome_paths = []
    not_found = []
    for gcf in gcf_list:
        #print(gcf)
        gcf_files = glob.glob('{}/Bacteria/{}*'.format(clark_db, gcf))
        if len(gcf_files) == 0:
            not_found.append(gcf)
        else:
            for gcf_path in gcf_files:
                genome_paths.append(gcf_path)
    return genome_paths, not_found

def create_clark_sub_database(time_range, clark_db, output_path):
    year = str(time_range[0].split("/")[-1])
    
    gcf_list = find_kraken_IDs(time_range)
    genome_paths, not_found = find_gcf_files(gcf_list, clark_db)
    with open(os.path.join(output_path, "sub_db_scrap/{}_genome_paths.txt".format(year)), "w+") as file:
        for gp in genome_paths:
            file.write(gp + "\n")
    
    if len(not_found) > 0:
        with open(os.path.join(output_path, "sub_db_scrap/{}_not_found.txt".format(year)), "w+") as nf:
            for gcf in not_found:
                nf.write(gcf + "\n")
    return

#%% MAIN

if __name__ == '__main__':
    end, start, clark_db, output_path = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    
    time_range = (end, start)
    create_clark_sub_database(time_range, clark_db, output_path)