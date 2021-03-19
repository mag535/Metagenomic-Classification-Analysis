# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 13:43:15 2020

@author: Melissa Gray
"""

import pandas as pd
import numpy as np
import os
import pickle
import sys
import matplotlib.pyplot as plt
from scipy.spatial import distance
import seaborn as sns


def create_bar_graphs(year_list, level, unclassified='y'):
    main_year_dic = {}
    uncl_dic = {}
    for year in year_list:
        year_dic = pickle.load(open("{}_taxa_frequency.p".format(year), "rb"))
        level_dic = year_dic[level]
        uncl_dic[year] = {}
        uncl_dic[year]["Unclassified"] = year_dic["Unclassified"]['0']
        main_year_dic[year] = level_dic
    year_df = pd.DataFrame.from_dict(main_year_dic, orient='columns').fillna(0)
    uncl_df = pd.DataFrame.from_dict(uncl_dic, orient='columns')
    
    # SORTING
    cls_sum_df = year_df.sum() #sum of classified
    print("cls sum:\n", cls_sum_df, end='\n\n')
    all_df = year_df.append(uncl_df)
    total_sum_df = all_df.sum() # sum of all
    print("total sum:\n", total_sum_df, end='\n\n')
    cls_percent_df = cls_sum_df.divide(total_sum_df).multiply(100) # percentage of classified
    print("cls percent:\n", cls_percent_df, end='\n\n')
    
    if unclassified.lower() == 'y':
        year_df = year_df.append(uncl_df)
    
    for col in year_df.columns:
        year_df[col] = year_df[col].divide(year_df[col].sum()).multiply(100)    # now a percentage
    
    others = []
    if unclassified.lower() == 'y':
        threshold = 5
    else:
        threshold = 10
    for col in year_df.columns:
        year_df.sort_values(by=col, axis=0, ascending=False, inplace=True)
        below_series = year_df[year_df[col] < threshold].loc[:, col]
        others.append(below_series.sum())
        for ind in below_series.index:
            year_df.at[ind, col] = 0
    oth = pd.DataFrame.from_dict({'Others':others}, orient='index', columns=year_df.columns)
    year_df = pd.concat([year_df, oth])
    
    
    # GRAPHING
    colors = ['#fc995b', #coral
              '#fcde5b', #pale yellow
              '#ceff3b', #yellow-green
              '#27b044', #green
              '#8cffee', #teal
              '#91caff', #pale blue
              '#bd9cff', #lavender
              '#fa9cff', #pink
              '#693800', #brown
              '#000000', #black
              '#d45200', #peachy-orange
              '#857532', #pale gold
              '#c9ffd5', #mint
              '#14594f', #turquiose
              '#5500ff', #neon purple
              ]
    fig = plt.subplots(figsize=(12,7))
    ind = np.arange(len(year_list))
    bars = []
    tmp_df = pd.DataFrame()
    o = "#01024d" #deep blue for others
    u = "#6e0000" #deep red for others
    c = 0
    for taxid in year_df.index:
        if (np.sum(year_df.loc[taxid,:]) != 0):
            p = -1
            if tmp_df.empty:
                if taxid == 'Others':
                    p = plt.bar(ind, list(year_df.loc[taxid,:]), color=o)
                elif taxid == 'Unclassified':
                    p = plt.bar(ind, list(year_df.loc[taxid,:]), color=u)
                else:
                    p = plt.bar(ind, list(year_df.loc[taxid,:]), color=colors[c])
                bars.append(p)
            else:
                if taxid == 'Others':
                    p = plt.bar(ind, list(year_df.loc[taxid,:]), bottom=list(tmp_df.sum()), color=o)
                elif taxid == 'Unclassified':
                    p = plt.bar(ind, list(year_df.loc[taxid,:]), bottom=list(tmp_df.sum()), color=u)
                else:
                    p = plt.bar(ind, list(year_df.loc[taxid,:]), bottom=list(tmp_df.sum()), color=colors[c])
                bars.append(p)
            tmp_df = tmp_df.append(year_df.loc[taxid,:])
            c += 1
    
    names_dic = pickle.load(open("names_dic.p", 'rb'))
    names_list = []
    for i in tmp_df.index:
        if (i == 'Others') or (i == 'Unclassified'):
            names_list.append(i)
        else:
            names_list.append(names_dic[i])
    
    #PLOT PERCENTAGE of classified as curved line
    line_dots = plt.plot(ind, list(cls_percent_df), 'd', color='#14b582') #dark teal
    
    legend_markers = [p[0] for p in bars]
    legend_markers.append(line_dots[0])
    names_list.append('Percent Classified')
    
    if unclassified == 'y':
        title = '{} Level'.format(level.capitalize())
        filename = 'Test_{}.png'.format(level)
    else:
        title = '{} Level excluding Unclassified'.format(level.capitalize())
        filename = 'Test_{}_wout.png'.format(level)
    plt.ylabel("Presence (percentage)", fontdict={'fontsize':16, 'fontweight':'semibold', 'style':'italic'})
    plt.yticks(np.arange(0,101,10))
    plt.xlabel("Year", fontdict={'fontsize':16, 'fontweight':'semibold', 'style':'italic'})
    plt.xticks(ind, year_list)
    plt.title(title, fontdict={'fontsize':24, 'fontweight':'bold'})
    plt.legend(legend_markers, names_list, bbox_to_anchor=(1, 1), prop={'size':16})
    
    plt.savefig(filename)
    return

def bray_curt_graph(year_list, level):
    main_year_dic = {}
    year_array = []
    year_array_2 = []
    uncl_dic = {}
    for year in year_list:
        year_dic = pickle.load(open("C:\\Users\\milkg\\OneDrive\\Desktop\\From_Proteus\\{}_taxa_frequency.p".format(year), "rb"))
        level_dic = year_dic[level]
        uncl_dic[year] = {}
        uncl_dic[year]['Unclassified'] = year_dic['Unclassified']['0']
        main_year_dic[year] = level_dic
    year_df = pd.DataFrame.from_dict(main_year_dic, orient='columns').fillna(0)
    uncl_df = pd.DataFrame.from_dict(uncl_dic, orient='columns')
    
    year_df_2 = year_df.append(uncl_df).fillna(0)
    for col in year_df.columns:
        year_array.append([val for val in year_df_2.loc[:, col]])
        year_array_2.append([val for val in year_df.loc[:, col]])
    m_one_array = np.array(year_array)
    
    # Method 2 Normalization
    m_two_array = []
    for i in range(len(year_array_2)):
        row_sum = np.sum(year_array_2[i])
        new_row = np.multiply(np.divide(year_array_2[i], row_sum), 100)
        m_two_array.append(list(new_row))
    m_two_array = np.array(m_two_array)
    
    bray_curt_one = []
    bray_curt_two = []
    labels = []
    for i in range(len(m_one_array)):
        if i == len(m_one_array) - 1:
            break
        else:
            bc = distance.pdist([m_one_array[i], m_one_array[i+1]], 'braycurtis')
            bray_curt_one.append(bc)
            bc2 = distance.pdist([m_two_array[i], m_two_array[i+1]], 'braycurtis')
            bray_curt_two.append(bc2)
            labels.append("-".join([str(year_list[i]), str(year_list[i+1])]))
    
    x = np.arange(0,len(labels), 1)
    line = plt.plot(x, bray_curt_one)
    line2 = plt.plot(x, bray_curt_two)
    
    plt.xticks(x, labels)
    plt.xlabel("Pair of Years")
    plt.ylabel("Bray Curtis Distance")
    plt.title("{} Level".format(level))
    plt.legend([line[0], line2[0]], ['Method 1', 'Method 2 Normalization'], bbox_to_anchor=(1, 0.8))
    plt.savefig("{}_bray_curtis_distance.png".format(level))
    return

def create_triangle_bray_curti(year_list, level):
    #import year dictionaries of level
    level_dic = {}
    for year in year_list:
        year_dic = pickle.load(open("C:\\Users\\milkg\\OneDrive\\Desktop\\From_Proteus\\{}_taxa_frequency.p".format(year), "rb"))
        level_dic[year] = year_dic[level]
    
    level_df = pd.DataFrame.from_dict(level_dic, orient='columns').fillna(0)
    level_array = []
    for col in level_df.columns:
        level_array.append([val for val in level_df.loc[:, col]])
    
    bray_curtis = distance.pdist(level_array, 'braycurtis')
    
    i=0
    prev = 0
    chart_dic = {}
    for year in year_list:
        if year not in chart_dic.keys():
            chart_dic[year] = []
        for j in range(i+1):
            chart_dic[year].append(np.nan)
        ####
        for k in range(prev, prev + (len(year_list)-(i+1))):
            chart_dic[year].append(bray_curtis[k])
        prev += len(year_list)-(i+1)
        i += 1
    
    chart_df = pd.DataFrame.from_dict(chart_dic, orient='index', columns=year_list)
    
    chart_df.fillna('--').to_csv('C:\\Users\\milkg\\EESI\\Summary\\{}-{}_{}_bray-curtis_triangular_distance_chart.csv'.format(year_list[0], year_list[-1], level))
    
    filled_chart_df = chart_df.fillna(0)
    
    plt.title("Bray Curtis Distance Between Sub-Databases\n{} Level".format(level.capitalize()))
    heat_ax = sns.heatmap(filled_chart_df, cmap='YlOrBr')
    plt.savefig('triangular_bray-curtis_distance_{}'.format(level))
    
    return

def unclassified(year_list):
    uncl_dic = {}
    for year in year_list:
        year_dic = pickle.load(open("C:\\Users\\milkg\\OneDrive\\Desktop\\From_Proteus\\{}_taxa_frequency.p".format(year), "rb"))
        uncl_dic[year] = year_dic['Unclassified']
    
    labels = []
    bray_curts = []
    years = [year for year in uncl_dic.keys()]
    for i in range(len(years)):
        if i == len(years)-1 :
            break
        else:
            labels.append(str(years[i]) + '-' + str(years[i+1]))
            bc = distance.pdist([[uncl_dic[years[i]]['0']], [uncl_dic[years[i+1]]['0']]], 'braycurtis')
            bray_curts.append(bc[0])
    
    x = np.arange(0,len(labels))
    plt.plot(x, bray_curts)
    plt.xticks(x, labels)
    plt.yticks(np.arange(0,1.1, .1))
    plt.show()
    return

def genome_relative_abundance(year, level, output='graph'):
    names_dic = pickle.load(open('C:\\Users\\milkg\\OneDrive\\Desktop\\From_Proteus\\names_dic.p', 'rb'))
    level_dic = pickle.load(open('C:\\Users\\milkg\\OneDrive\\Desktop\\From_Proteus\\{}_taxa_frequency.p'.format(year), 'rb'))[level]
    count_sum = np.sum([val for val in level_dic.values()])
    
    # condensing
    condensed_level_dic = {}
    taxids = [key for key in level_dic.keys()]
    percentages = np.divide([val for val in level_dic.values()], count_sum)
    percentages = zip(tuple(taxids), tuple(percentages))
    exclude = []
    others = 0
    if output.lower() == 'graph':
        for taxid, percent in percentages:
            if percent < 0.05:
                exclude.append((taxid, percent))
                others += percent
            else:
                condensed_level_dic[names_dic[taxid]] = percent
    elif output.lower() == 'sheet':
        for taxid, percent in percentages:
            condensed_level_dic[names_dic[taxid]] = percent
    
    keys = [key for key in condensed_level_dic.keys()]
    sorted_keys = sorted(keys)
    sorted_condensed_level_dic = {}
    for skey in sorted_keys:
        sorted_condensed_level_dic[skey] = condensed_level_dic[skey]
    
    if output.lower() == 'graph':
        if len(exclude) > 1:
            sorted_condensed_level_dic['Others'] = others
        elif len(exclude) == 1:
            sorted_condensed_level_dic[names_dic[exclude[0][0]]] = exclude[0][1]
    
    # graphing
    if output.lower() == 'graph':
        x = np.arange(0, len([key for key in sorted_condensed_level_dic.keys()]))
        y = [val for val in sorted_condensed_level_dic.values()]
        plt.bar(x, y)
        
        plt.xlabel('Genomes')
        plt.xticks(x, [key for key in sorted_condensed_level_dic.keys()])
        plt.ylabel('Relative Abundance')
        plt.yticks(np.arange(0, 1.1, 0.1))
        plt.title("Relative Abundance of Genomes\nClassified by the {} Sub-Database for {} Level".format(year, level.capitalize()))
        
        plt.show()
        
    elif output.lower() == 'sheet':
        return pd.DataFrame.from_dict(sorted_condensed_level_dic, orient='index', columns=['Relative Abundance'])
    
    '''
    with open('C:\\Users\\milkg\\EESI\\Summary\\relative_abundance_{}_{}.txt'.format(level, year), 'w') as file:
        for key in sorted_condensed_level_dic:
            file.write(key + " : " + str(sorted_condensed_level_dic[key]) + "\n")
    '''
    return

def calculate_relative_abundance(year_list, level):
    excel_name = 'Relative_Abundance_Classified_{}.xlsx'.format(level)
    with pd.ExcelWriter(excel_name) as writer:
        for year in year_list:
            year_df = genome_relative_abundance(year, level, 'sheet')
            year_df.to_excel(writer, sheet_name=str(year))
    return

#%% MAIN

if __name__ == "__main__":
    #path = sys.argv[1]
    
    #create_prelim(library, output_path)
    
    #create_manifest(path)
    
    '''
    level, uncl, year_list = sys.argv[1], sys.argv[2], sys.argv[3:]
    create_bar_graphs(year_list, level, uncl)
    '''
    
    #bray_curt_graph([1999, 2004, 2009, 2014, 2019, 2020], 'phylum')
    #arr = np.array([[1,2,3], [4,5,6], [4,2,1], [4,5,2]])
    
    '''
    level, year_list = sys.argv[1], sys.argv[2:]
    bray_curt_graph(year_list, level)
    '''
    
    #genome_relative_abundance(2009, 'phylum')
    
    #create_triangle_bray_curti([1999, 2004, 2009, 2014, 2019, 2020], 'species')
    
    calculate_relative_abundance([1999, 2004, 2009, 2014, 2019, 2020], 'species')