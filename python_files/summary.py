# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 10:48:41 2020

@author: Melissa Gray
"""

#%% IMPORTS & VARS
import pandas as pd
import numpy as np
import re
import os.path
import pickle
from Bio import Entrez
from time import sleep


Entrez.email = "mag535@drexel.edu"


#%% CLASSES

class Summary():
    def __init__(self, file):
        self.file = file
        ### gonna make reference file and use that
        return
    
    # SETTERS
    def set_file_name(self, fn):
        self.file = fn
        return
    # GETTERS
    def get_file(self):
        f = open(self.file, encoding='utf-8')
        file_content = f.readlines()
        f.close()
        return file_content
    
    # FUNCTIONS
    def split_and_strip(self, line):
        lst_1 = re.split("\t", line, 15)
        lst_1.pop()
        return lst_1
    
    def get_categories(self, content):
        cat_line = ''
        for line in content:
            if re.search("^#.*", line) and (".txt" not in line):
                cat_line = line[2:]
                break
        categories = [cat.strip().capitalize() for cat in re.split("\t+", cat_line, 15)]
        categories.pop()
        return categories
    
    def parse_content(self, content):
        categories = self.get_categories(content)
        
        # numpy array -> pandas dataframe
        parsed = []
        for line in content:
            if not re.search("^#.*|^@.*", line):
                new_line = self.split_and_strip(line)
                parsed.append(new_line)
        
        ref_df = pd.DataFrame(np.array(parsed), columns=categories)
        ref_df['Seq_rel_date'] = pd.to_datetime(ref_df['Seq_rel_date'])
        return ref_df
    
    def save_df(self, df, df_name='reference'):
        df.to_csv(df_name, index=False, header=True)
        return
    
    def find_date(self, time_list=("-1/-1/-1", "-1/-1/-1"), name=''): #end_year=-1, start_year=-1):
        # Gathers all the records from during and/or before the year specified
        try:
            ref_df = pd.read_csv('reference')
        except FileNotFoundError:
            self.setup()
            ref_df = pd.read_csv('reference')
        ref_df['Seq_rel_date'] = pd.to_datetime(pd.read_csv('reference')['Seq_rel_date'])
        
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
            time_df.sort_values(by='Seq_rel_date').to_csv(name + '.csv', index=False, header=True)
            print("\'{}.csv\' is saved.".format(name))
        return time_df
    
    def find_kraken_IDs(self, times):
        tmp_df = self.find_date(times)
        gcf_list = tmp_df['Assembly_accession']
        return gcf_list
    
    def check_final_dictionary(self, version, level):
        version_df = pd.DataFrame()
        self.set_file_name('/mnt/HA/groups/rosenGrp/mag535/kraken_misc_files/assembly_summary.txt')
        old_df = self.parse_content(self.get_file())
        
        version_df = old_df[old_df["Version_status"] == version]
        final_df = version_df[version_df["Assembly_level"] == level]
        target = final_df['Assembly_accession']
        
        d = open("/lustre/scratch/mag535/nz_gcf.p", "rb")
        ID_dic = pickle.load(d)
        d.close()
        gcf_dic_set = set([val.strip() for val in ID_dic.values()])
        
        missing_gcf = open("/lustre/scratch/mag535/missing_gcf.txt", "w")
        for gcf in target:
            if gcf not in gcf_dic_set:
                missing_gcf.write(gcf + "\n")
        missing_gcf.close()
        return
        
    def gb_id_map(self, REF_SEQ_ID):
        handle = Entrez.efetch(db="nucleotide", id=REF_SEQ_ID, retmode="xml")
        Assembly_ID = None
        res = handle.read()
        #print(res, end="\n")
        no_line = open("/lustre/scratch/mag535/no_line.txt", "w+")
        many_lines = open("/lustre/scratch/mag535/many_lines.txt", "w+")
        
        found = re.findall(".*GCF_.*", res)
        #print(found, end="\n")
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
    
    def dictionary_ID(self):
        ###Entrez.email = 'mag535@drexel.edu'
        
        # cols => ["Assembly_accession", "Taxid", "Organism_name", "Infraspecific_name"]
        # cols => [         (0),           (5),         (7),                (8)]
        try:
             ref = pd.read_csv('reference')
        except FileNotFoundError:
            self.setup()
            ref = pd.read_csv('reference')
        
        ID_dic = {}
        headers = open("/lustre/scratch/mag535/Kraken/bac_kraken_db/library/bacteria/library_headers.fna", "r", encoding='utf-8')
        no_matches = open("/lustre/scratch/mag535/no_matches.txt", "w+")
        output = open('/lustre/scratch/mag535/nz_gcf_dict.txt', 'w+')
        
        for line in headers:
            parts = line.split("|")
            taxid = int(parts[1])
            name = re.split(" ", parts[-1], 1)[-1]
            nz = parts[-1].split(" ")[0]
            
            temp_df = ref[ref['Taxid'] == taxid]
            
            if temp_df.shape[0] == 1:
                ID_dic[nz] = temp_df.iloc[0, 0]
                output.write(nz + " : " + temp_df.iloc[0, 0] + "\n")
            elif temp_df.shape[0] == 0:
                no_matches.write(line)
            elif temp_df.shape[0] > 1:
                flag = False
                for ind in temp_df.index:
                    strain = ''
                    try:
                        strain = temp_df.loc[ind, "Infraspecific_name"].split("strain=")[-1]
                    except AttributeError:
                        strain = ''
                    if (len(strain) > 0) and (strain in temp_df.loc[ind, "Organism_name"]):
                        if temp_df.loc[ind, "Organism_name"] in name:
                            ID_dic[nz] = temp_df.loc[ind, "Assembly_accession"]
                            output.write(nz + " : " + temp_df.loc[ind, "Assembly_accession"] + "\n")
                            flag = True
                            break
                    elif (len(strain) == 0):
                        if (temp_df.loc[ind, "Organism_name"] in name):
                            ID_dic[nz] = temp_df.loc[ind, "Assembly_accession"]
                            output.write(nz + " : " + temp_df.loc[ind, "Assembly_accession"] + "\n")
                            flag = True
                            break
                    else:
                        if (temp_df.loc[ind, "Organism_name"] in name) and (strain in name):
                            ID_dic[nz] = temp_df.loc[ind, "Assembly_accession"]
                            output.write(nz + " : " + temp_df.loc[ind, "Assembly_accession"] + "\n")
                            flag = True
                            break
                if flag == False:
                    no_matches.write(line)
        
        '''
        no_matches.seek(0)
        no_mat = no_matches.readlines()
        for line in no_mat:
            parts = line.split("|")
            nz = parts[-1].split(" ")[0]
            print("~~", nz)
            ref_seq_id, a_id = self.gb_id_map(nz)
            ID_dic[nz] = ref_seq_id
            output.write(ref_seq_id + " : " + str(a_id))
            sleep(3.5)
        '''
        
        headers.close()
        no_matches.close()
        output.close()
        
        '''
        d = open("/lustre/scratch/mag535/nz_gcf.p", "wb" )
        pickle.dump(ID_dic, d)
        d.close()
        '''
        return
    
    def gather_data(self, gcf_list, output_name):
        d = open("/mnt/HA/groups/rosenGrp/mag535/kraken_misc_files/nz_gcf.p", "rb" )
        pre_ma = pickle.load(d)
        d.close()
        target = set()
        for item in gcf_list:
            for key, val in pre_ma.items():
                if item == val:
                    target.add(key)
                    break
        
        library_file = open("/lustre/scratch/mag535/Kraken/bac_kraken_db/library/bacteria/library.fna", "r", encoding='utf-8')
        output = open(output_name, "w+", encoding='utf-8')
        
        flag = False
        for idx, con in enumerate(library_file):
            if con[0] == ">":
                nz = con.split("|")[-1].split(" ")[0]
                if nz in target:
                    output.write(con)
                    flag = True
                else:
                    flag = False
            else:
                if flag:
                    output.write(con)
        
        library_file.close()
        output.close()
        return #self.check_data()
    
    def check_data(self, times, year_dir):
        accession_dic = pickle.load(open("/mnt/HA/groups/rosenGrp/mag535/kraken_misc_files/nz_gcf.p", "rb"))
        gcf_df = self.find_kraken_IDs(times)
        gcf = gcf_df.sample().iloc[0]
        
        nz = ''
        nz_list = []
        for key, val in accession_dic.items():
            if val in gcf_df.values:
                nz_list.append(key)
                if val == gcf:
                    nz = key
        
        # searches
        file = open("/mnt/HA/groups/rosenGrp/mag535/{}/library.fna".format(year_dir), "r", encoding='utf-8')
        temp = []
        missing_nz = []
        
        flag = False
        for line in file:
            if line[0] == ">":
                accession = line.split("|")[-1].split(" ")[0]
                if accession not in nz_list:
                    missing_nz.append(accession)
                if nz in line:
                    flag = True
                    temp.append(line)
                else:
                    flag = False
            else:
                if flag:
                    temp.append(line)
        
        file.close()
        original = open("/lustre/scratch/mag535/Kraken/bac_kraken_db/library/bacteria/library.fna", "r", encoding='utf-8')
        genome = []
        
        flag = False
        for line in original:
            if line[0] == ">":
                if nz in line:
                    flag = True
                    genome.append(line)
                else:
                    if flag:
                        break
                    else:
                        flag = False
            else:
                if flag:
                    genome.append(line)
                            
        
        original.close()
        
        # info
        print(year_dir)
        print(times)
        
        # checks
        if len(temp) == 0:
            print("Couldn't find nz accession in sub_library.fna file")
        elif len(genome) == 0:
            print("Couldn't find NZ accession in main library.fna file")
        elif len(temp) == len(genome):
            if all([temp[i]==genome[i] for i in range(len(temp))]):
                print("The header and genome match")
            else:
                print("The lines in the genomes don't match")
        else:
            print("The length of the genomes don't match")
        
        # missing nz's
        print("MISSING Accessions:")
        if len(missing_nz) == 0:
            print("No missing accessions")
        else:
            for nz in missing_nz:
                print(nz)
        
        return
    
    def setup(self):
        ref_df = self.parse_content(self.get_file()).sort_values(by='Seq_rel_date')
        self.save_df(ref_df)
        return


#%% MAIN

if __name__ == '__main__':
    '''
    path = '' #'Kraken/bac_kraken_db/library/bacteria'
    Sum = Summary(os.path.join(path, "assembly_summary.txt"))
    
    Sum.find_date(("11/22/2003", "10/20/2001"), "birth_year_test")
    '''