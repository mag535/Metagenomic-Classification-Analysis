# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 10:05:54 2021

@author: Melissa Gray

Graphing for CLARK classification files.

What to calculate:
    - relative abundance per taxanomy level (specifically phylum and species) for each year
    - triangle bray curtis distance
    - taxid occurance w/ traceback (as a dict)
    - create traceback dictionary (?)
    
What to graph/print:
    - text file of relative abundance
    - bar graph of relative abundance
        - one that includes unclassified
        - one that includes relative of abundance of unclassified as scatter plot points
"""



#%% VARIABLES & IMPORTS

import pickle, os, sys # all built-in modules: https://docs.python.org/3/py-modindex.html
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
#import matplotlib.colors as mcolors

import pandas as pd
import numpy as np
from scipy.spatial import distance
import seaborn as sns

# color dictionaries for graphing
rank_colors = {
    "Unclassified": "#6e0000", # deep red
    "Others": "#01024d", # deep blue
    "Percent Classified": "aliceblue",
    "phylum": { #base colors
        "Spirochaetes": "yellow",
        "Proteobacteria": "magenta",
        "Chlamydiae": "red",
        "Aquificae": "cyan",
        "Bacteroidetes": "blue",
        "Firmicutes": "green",
        "Verrucomicrobia": "grey"
         },
    "class": { #tab;eau palette
        "Spirochaetia": "darkkhaki",
        "Epsilonproteobacteria": "tab:pink",
        "Chlamydiia": "tab:red",
        "Aquificae": "tab:cyan",
        "Bacteroidia": "tab:blue",
        "Clostridia": "tab:green",
        "Negativicutes": "tab:olive",
        "Verrucomicrobiae" : "tab:gray",
        "Bacilli": "tab:orange",
        "Gammaproteobacteria": "tab:purple",
        "Betaproteobacteria": "tab:brown"
        },
    "order": { #CSS colors
        "Spirochaetales": "gold",
        "Campylobacterales": "darkmagenta",
        "Chlamydiales": "darkred",
        "Aquificales": "darkcyan",
        "Bacteroidales": "darkblue",
        "Verrucomicrobiales": "darkgrey",
        "Clostridiales": "darkgreen",
        "Acidaminococcales": "forestgreen"
        },
    "family": { #CSS colors
        "Spirochaetaceae" : "khaki",
        "Helicobacteraceae": "lightpink",
        "Chlamydiaceae": "lightcoral",
        "Aquificaceae": "lightcyan",
        "Bacteroidaceae": "cornflowerblue",
        "Tannerellaceae": "lightblue",
        "Akkermansiaceae": "lightgrey",
        "Lachnospiraceae": "lightgreen",
        "Odoribacteraceae": "skyblue",
        "Rikenellaceae": "plum",
        "Barnesiellaceae": "lavender",
        "Ruminococcaceae": "springgreen",
        "Acidaminococcaceae": "lime",
        "Oscillospiraceae": "mintcream",
        "Porphyromonadaceae": "aqumarine",
        "Prevotellaceae": "thistle"
        },
    "genus": { #CSS colors
        "Treponema" : "goldenrod",
        "Helicobacter": "mediumvioletred",
        "Chlamydia": "firebrick",
        "Aquifex": "teal",
        "Bacteroides": "slateblue",
        "Phocaeicola": "darkslateblue",
        "Roseburia": "mediumseagreen",
        "Odoribacter": "midnightblue",
        "Alistipes": "royalblue",
        "Faecalibacterium": "lightseagreen",
        "Phascolarctobacterium": "turquoise",
        "Parabacteroides": "dodgerblue",
        "Akkermansia": "gainsboro",
        "Dysosmobacter": "aquamarine",
        "Blautia": "greenyellow",
        "Porphyromonas": "aqua"
        },
    "species" : { #CSS colors
        "Treponema pallidum": "orange",
        "Helicobacter pylori": "crimson",
        "Chlamydia pneumoniae": "maroon",
        "Aquifex aeolicus": "teal",
        "Bacteroides thetaiotaomicron": "steelblue",
        "Bacteroides fragilis": "#4a2500", # deep brown
        "Phocaeicola vulgatus": "#14b582", # dark teal
        "Roseburia hominis": "palegreen",
        "Odoribacter splanchnicus": "powderblue",
        "Bacteroides helcogenes": "#ffd885", # pale orange
        "Bacteroides ovatus": "darkturquoise",
        "Bacteroides sp. A1C1": "indigo",
        "Bacteroides caccae": "darkslategrey",
        "Bacteroides uniformis": "darkslateblue",
        "Phocaeicola dorei": "mediumorchid",
        "Faecalibacterium prausnitzii": "lawngreen",
        "Bacteroides xylanisolvens": "saddlebrown",
        "Alistipes communis": "brown",
        "Bacteroides sp. CACC 737": "plum",
        "Porphyromonas gingivalis": "tomato",
        "Alistipes shahii": "darkorange",
        "Bacteroides sp. M10": "rosybrown",
        "Bacteroides sp. HF-5141": "chocolate"
        }
    }


#%% FUNCTIONS

# The Maths
def trace_back(taxid, nodes_dic, output_path=""):
    if taxid not in nodes_dic:
        #### write to a file
        with open(os.path.join(output_path, "no_traceback.txt"), 'w') as file:
            file.write(taxid + '\n')
        return []
    res = []
    cur_id = taxid
    while cur_id != "1":
        parent_id, rank = nodes_dic[cur_id]
        res.append((cur_id, rank))
        cur_id = parent_id
    return res

def kraken2_taxa_frequency(year, input_path="", output_path="", to_text=0):
    nodes_dic = pickle.load(open(os.path.join(output_path, "nodes_dic.p"), "rb"))
    frequency = {}
    with open(os.path.join(input_path, "kraken2_{}_real_classify_output.txt".format(year)), "r", encoding='utf-8') as file:
        for line in file:
            parse = [info.strip() for info in line.split("\t")]
            if parse[0] == "U":
                # -2 will signify unclassified, like in array-of-known
                if "Unclassified" not in frequency:
                    frequency["Unclassified"] = {}
                if -2 not in frequency["Unclassified"]:
                    frequency["Unclassified"][-2] = 0
                frequency["Unclassified"][-2] += 1
            else:
                res = trace_back(parse[2], nodes_dic, output_path)
                for p in res:
                    if p[1] not in frequency:
                        frequency[p[1]] = {}
                    if p[0] not in frequency[p[1]]:
                        frequency[p[1]][p[0]] = 0
                    frequency[p[1]][p[0]] += 1
    pickle.dump(frequency, open(os.path.join(output_path, "kraken2_{}_taxa_frequency.p".format(year)), "wb"))
    
    if to_text == 1:
        with open(os.path.join(output_path, "kraken2_{}_taxa_frequency.txt".format(year)), "w") as output:
            for rank in frequency.keys():
                output.write(rank + "\n")
                for taxa,freq in frequency[rank].items():
                    output.write("\t" + str(taxa) + " : " + str(freq) + "\n")
    return frequency

def clark_taxa_occurance(year, input_path="", output_path="", to_text=0):
    # creates a dictionary for only one year's taxa frequency
    nodes_dic = pickle.load(open(os.path.join(output_path, "nodes_dic.p"), "rb"))
    print("nodes loaded")
    frequency = {}
    """
    dictionary is this  -->
    
    {Rank : 
     {taxid : count}
     }
    
    To Acess: frequency[Rank][taxid] = count
    """
    
    results = pd.read_csv(os.path.join(input_path, "clark_{}_bac_real_classify_output.csv".format(year)))
    print("results read")
    for ind in results.index:
        assignment = results.loc[ind, " Assignment"]
        if np.isnan(assignment):
            if "Unclassified" not in frequency:
                frequency["Unclassified"] = {}
            # -2 will signify unclassified, like in array-of-known
            if -2 not in frequency["Unclassified"]:
                frequency["Unclassified"][-2] = 0
            frequency["Unclassified"][-2] += 1
        else:
            res = trace_back(str(int(assignment)), nodes_dic, output_path)
            for p in res:
                taxid, rank = p
                if rank not in frequency:
                    frequency[rank] = {}
                if taxid not in frequency[rank]:
                    frequency[rank][taxid] = 0
                frequency[rank][taxid] += 1
    print("taxa counted")
    pickle.dump(frequency, open(os.path.join(output_path, "clark_{}_taxa_frequency.p".format(year)), "wb"))
    print("frequency dic -> pickle file")
    
    if to_text == 1:
        with open(os.path.join(output_path, "clark_{}_taxa_frequency.txt".format(year)), "w") as output:
            for rank in frequency.keys():
                output.write(rank + "\n")
                for taxa,freq in frequency[rank].items():
                    output.write("\t" + str(taxa) + " : " + str(freq) + "\n")
        print("frequency written")
    return frequency

# The Graphs
def _export_excel(classified_df, unclassified_df, names_dic, rank, who, output_path=""):
    # whole number -> fraction -> percent
    classified_relative_abundance = classified_df.divide(classified_df.sum()).multiply(100)
    all_df = (classified_df.append(unclassified_df)).fillna(0)
    all_relative_abundance = all_df.divide(all_df.sum()).multiply(100)
    # change indices from taxid to names
    for ind in classified_relative_abundance.index:
        if ind == -2:
            classified_relative_abundance.rename(index={ind : "Unclassified"}, inplace=True)
        else:
            classified_relative_abundance.rename(index={ind : names_dic[str(ind)]}, inplace=True)
    for ind in all_relative_abundance.index:
        if ind == -2:
            all_relative_abundance.rename(index={ind : "Unclassified"}, inplace=True)
        else:
            all_relative_abundance.rename(index={ind : names_dic[str(ind)]}, inplace=True)
    # export both DataFrames
    with pd.ExcelWriter(os.path.join(output_path, "{}_{}_relative_abundance_percentage.xlsx".format(who, rank)), engine='xlsxwriter') as writer:
        classified_relative_abundance.to_excel(writer, sheet_name='Classified Only')
        all_relative_abundance.to_excel(writer, sheet_name='All')
    print("ra -> excel")
    return all_relative_abundance.loc["Unclassified", :]

def _group_others(relative_abundance_df, threshold):
    others = []
    for col in relative_abundance_df.columns:
        relative_abundance_df.sort_values(by=col, axis=0, ascending=False, inplace=True)
        too_small = relative_abundance_df[relative_abundance_df[col] < threshold].loc[:, col]
        others.append(too_small.sum())
        for ind in too_small.index:
            relative_abundance_df.at[ind, col] = 0
    others_df = pd.DataFrame.from_dict({"Others" : others}, orient="index", columns=relative_abundance_df.columns)
    print("small nums -> others")
    return others_df

def plot_relative_abundance(year_list, rank, classifier, output_path="", input_path="", unclassified=0, excel=0):
    # bar graph
    # unclassified=to include unclassified as a bar : 1=yes, 0=no
    # % classified is included as plot points
    # excel=to export dataframe as an excel sheet : 1=yes, 0=no
    names_dic = pickle.load(open(os.path.join(output_path, "names_dic.p"), "rb"))
    collective_rank_dic = {}
    collective_unclassified_dic = {}
    
    who = ""
    if classifier.lower() == "c":
        who = "clark"
    elif classifier.lower() == "k":
        who = "kraken2"
    
    # GATHER DATA
    for year in year_list:
        # just in case
        try:
            year_dic = pickle.load(open(os.path.join(output_path, "{}_{}_taxa_frequency.p".format(who, year)), "rb"))
            print("freq loaded")
        except FileNotFoundError:
            if classifier.lower() == "c":
                year_dic = clark_taxa_occurance(year, input_path, output_path, 0)
            elif classifier.lower() == "k":
                year_dic = kraken2_taxa_frequency(year, input_path, output_path, 0)
            print("freq created")
        
        collective_rank_dic[year] = year_dic[rank]
        collective_unclassified_dic[year] = year_dic["Unclassified"]
    collective_rank_df = pd.DataFrame.from_dict(collective_rank_dic, orient="columns").fillna(0)
    collective_unclassified_df = pd.DataFrame.from_dict(collective_unclassified_dic, orient="columns")
    
    # DATA -> EXCEL
    if excel == 1:
        unclassified_percentage = _export_excel(collective_rank_df, collective_unclassified_df, names_dic, rank, who, output_path)
    
    # CONDENSING & SORTING
    threshold = 0
    if unclassified == 0:
        if classifier == "k":
            threshold = 5.0 # %
        elif classifier == "c":
            threshold = 2.0 # %
    elif unclassified == 1:
        if classifier == "k":
            threshold = 3.0 # %
        elif classifier == "c":
            threshold = 1.0 # %
        # add unclassified
        collective_rank_df = collective_rank_df.append(collective_unclassified_df) #.fillna(0)
        print("uncl added to collective")
    
    if (excel == 0) and (unclassified == 0):
        all_df = collective_rank_df.append(collective_unclassified_df)
        unclassified_percentage = all_df.divide(all_df.sum()).multiply(100).loc[-2, :] # whole number -> fraction -> percent -> get unclassified
    collective_rank_df = collective_rank_df.divide(collective_rank_df.sum()).multiply(100) # whole number -> fraction -> percent
    print("rel abun calculated")
    
    others_df = _group_others(collective_rank_df, threshold)
    collective_rank_df = pd.concat([collective_rank_df, others_df])
    print("others grouped & added to collective")
    
    # GRAPHING
    
    names_from_dic = list(names_dic.values())
    ids_from_dic = list(names_dic.keys())
    
    colors = rank_colors[rank]
    oth = rank_colors["Others"]
    uncl = rank_colors["Unclassified"]
    c = 0 #counter for what color to use for which bar
    x_range = np.arange(len(year_list))
    bars = [] #to keep track of the stacked bars for the graph
    tmp_df = pd.DataFrame() #to keep track of the "new bottom" for the stacked bars in the graph
    
    # The Stacked Bar Graphs
    fig, ax = plt.subplots(figsize=(17,7)) # the figure to plot on
    for taxa in collective_rank_df.index:
        if (np.sum(collective_rank_df.loc[taxa, :]) != 0):
            p = -1 #the bar graph plot
            if tmp_df.empty: #aka if just starting to create stacked bars
                if taxa == "Others":
                    p = ax.bar(x_range, 
                              collective_rank_df.loc[taxa, :].tolist(),
                              color=oth,
                              edgecolor="#000000",
                              linewidth=1,
                              zorder=-1)
                elif taxa == -2:
                    p = ax.bar(x_range,
                               collective_rank_df.loc[taxa, :].tolist(),
                               color=uncl,
                               edgecolor="#000000",
                               linewidth=1,
                               zorder=-1)
                else:
                    p = ax.bar(x_range,
                               collective_rank_df.loc[taxa, :].tolist(),
                               color=colors[names_from_dic[ids_from_dic.index(taxa)]],
                               edgecolor="#000000",
                               linewidth=1,
                               zorder=-1)
            else: #aka if adding more bars to the stack
                if taxa == "Others":
                    p = ax.bar(x_range,
                               collective_rank_df.loc[taxa, :].tolist(),
                               bottom=tmp_df.sum().tolist(),
                               color=oth,
                               edgecolor="#000000",
                               linewidth=1,
                               zorder=-1)
                elif taxa == -2:
                    p = ax.bar(x_range,
                               collective_rank_df.loc[taxa, :].tolist(),
                               bottom=tmp_df.sum().tolist(),
                               color=uncl,
                               edgecolor="#000000",
                               linewidth=1,
                               zorder=-1)
                else:
                    p = ax.bar(x_range,
                               collective_rank_df.loc[taxa, :].tolist(),
                               bottom=tmp_df.sum().tolist(),
                               color=colors[names_from_dic[ids_from_dic.index(taxa)]],
                               edgecolor="#000000",
                               linewidth=1,
                               zorder=-1)
            bars.append(p)
            tmp_df = tmp_df.append(collective_rank_df.loc[taxa, :])
            c += 1
    print("bars stacked")
    
    # for the legend
    names_list = []
    names_to_id = {}
    for i in tmp_df.index:
        if (i == "Others"):
            names_list.append(i)
        elif (i == -2) or (i == "Unclassified"):
            names_list.append("Unclassified")
        else:
            names_list.append(names_dic[i])
            names_to_id[names_dic[i]] = i
    
    # The Percentage of Classified
    if (excel == 0) and (unclassified == 1):
        unclassified_percentage = np.array(collective_rank_df.loc[-2, :].tolist())
    line_dots = ax.scatter(x_range,
                           100 - unclassified_percentage,
                           c=rank_colors["Percent Classified"],
                           marker="D",
                           linewidths=2.0,
                           edgecolor="#000000",
                           zorder=1)
    
    # graphing the legend
    legend_markers = [p[0] for p in bars]
    legend_markers.append(line_dots)
    names_list.append('Percent Classified')
    
    # graph/file name setup
    if who == "clark":
        who = who.upper()
    elif who == "kraken2":
        who = who.capitalize()
    if unclassified == 1:
        title = "{} Relative Abundance of {} Level Genomes".format(who, rank.capitalize())
        subtitle = "including Unclassified"
        filename = '{}_{}_TEST.png'.format(who.lower(), rank)
    else:
        title = "{} Relative Abundance of {} Level Genomes".format(who, rank.capitalize())
        subtitle = "excluding Unclassified"
        filename = '{}_{}_TEST_wout.png'.format(who.lower(), rank)
    # adding details to the graph
    plt.xlabel('Year (Sub-Database)',
               fontdict={"font" : 'Courier New',
                         "weight" : 'semibold',
                         "size" : 16,
                         "style" : 'italic'})
    plt.xticks(x_range, year_list)
    plt.ylabel('Relative Abundance (%)',
               fontdict={"font" : 'Courier New',
                         "weight" : 'semibold',
                         "size" : 16,
                         "style" : 'italic'})
    plt.yticks(np.arange(0, 110, 10))
    plt.suptitle(title, font='Courier New', weight='bold', size=24) # main title
    plt.title(subtitle,
              fontdict={"font" : 'Courier New',
                        "weight" : 'semibold',
                        "size" : 18, }) # subtitle
    fig.legend(legend_markers,
               names_list,
               #scatterpoints=1,
               bbox_to_anchor=(1.1, 0.9),
               borderaxespad=0.1,
               prop={'size':14}
               )
    
    #plt.show()
    print("ra graphed")
    plt.savefig(os.path.join(output_path, filename), bbox_inches="tight")
    print("ra -> png")
    return

def plot_triangular_bray_curtis(year_list, rank, classifier, output_path="", input_path="", unclassified=0):
    # heatmap
    # uses taxa frequency, separate by rank
    rank_dic = {}
    who = ""
    classifier = classifier.lower()
    if classifier == "c":
        who = "clark"
    elif classifier == "k":
        who = "kraken2"
    else:
        print("Warning: An incorrect value was passed into argument 'classifier'.")
    
    for year in year_list:
        # just in case
        try:
            year_dic = pickle.load(open(os.path.join(output_path, "{}_{}_taxa_frequency.p".format(who, year)), "rb"))
            print("dic loaded")
        except FileNotFoundError:
            if classifier == "c":
                year_dic = clark_taxa_occurance(year, input_path, output_path, 0)
            elif classifier == "k":
                year_dic = kraken2_taxa_frequency(year, input_path, output_path, 0)
            print("dic created")
        
        rank_dic[year] = year_dic[rank]
        if unclassified == 1:
            rank_dic[year][-2] = year_dic["Unclassified"][-2]
    
    rank_df = pd.DataFrame.from_dict(rank_dic, orient='columns').fillna(0)
    rank_array = []
    for col in rank_df.columns:
        rank_array.append(rank_df.loc[:, col].tolist())
    
    bray_curtis = distance.pdist(rank_array, 'braycurtis')
    print("bray curtis calculated")
    
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
    print("triangular sorted")
    
    chart_df = pd.DataFrame.from_dict(chart_dic, orient='index', columns=year_list)
    
    title1 = ""
    title2 = ""
    filename = ""
    if who == "clark":
        title1 = "%s Bray Curtis Distance Between Sub-Databases" %who.upper()
    elif who == "kraken2":
        title1 = "%s Bray Curtis Distance Between Sub-Databases" %who.capitalize()
    else:
        print("Warning: Variable 'who' was not set properly.")
    if unclassified == 0:
        title2 = "%s Level, Classified" %rank.capitalize()
        filename = "{}_{}_classified_triangular_bray-curtis_distance".format(who, rank)
    elif unclassified == 1:
        title2 = "%s Level, Unclassified" %rank.capitalize()
        filename = "{}_{}_unclassified_triangular_bray-curtis_distance".format(who, rank)
    else:
        print("Warning: An incorrect value was passed into argument 'unclassified'.")
    
    chart_df.fillna("--").to_csv(os.path.join(output_path, filename + ".csv"))
    print("chart -> csv")
    
    filled_chart_df = chart_df.fillna(0)
    
    plt.suptitle(title1)
    plt.title(title2)
    heat_ax = sns.heatmap(filled_chart_df, cmap='YlOrBr')
    plt.savefig(os.path.join(output_path, filename + ".png"))
    print("graph -> png")
    return

#%% MAIN

if __name__ == "__main__":
    # input from slurm script
    ip, op, level, who, uncl, excel, plot, yy = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8:]
    
    # creating dictionaries, calculations, and graphs
    
    #res = "C:\\Users\\milkg\\OneDrive\\Desktop\\From_Picotte\\clark_bac_generated_classify_output_1.csv"
    #op = "C:\\Users\\milkg\\OneDrive\\Desktop"
    #ip = "C:\\Users\\milkg\\OneDrive\\Desktop\\From_Picotte"
    #level = "phylum"
    #yy = [1999]
    
    '''
    if who.lower() == "c":
        clark_taxa_occurance(yy[0], ip, op, 1)
    elif who.lower() == "k":
        kraken2_taxa_frequency(yy[0], ip, op, 1)
    '''
    #kraken2_taxa_frequency(1999, ip, op, 1)
    
    if int(plot) == 1:
        plot_relative_abundance(yy, level, who, op, ip, int(uncl), int(excel))
    elif int(plot) == 2:
        plot_triangular_bray_curtis(yy, level, who, op, ip, int(uncl))
    
    