#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 10:16:47 2019

@author: ptruong
"""

import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import seaborn as sns
from triqlerParser import *
from triqlerProcessor import * 
from parseReformattedPSSS3equDecoy import *   

def pivotTriqlerFormattedData(triqlerInput = "../data/PSSS3_triqlerFormatted_nonShared.csv", path = "wideint.pkl"):
    """
    triqlerInput - data in triqler input format.
    path - path to check for data and output data.
    condition - run or condition
    """
    result_file = Path(path)
    if not result_file.exists():
        long = pd.read_table(triqlerInput)
        long["species"] = long.proteins.str[-5:]
        long["runs"] = long.run.str[-3:]
        long["decoy"] = (long.proteins.str[:5] == "decoy")
        wideInt = long.pivot_table(index=['proteins','decoy','species','peptide'], columns=["condition","runs"], values='intensity')
        wideInt.to_pickle(result_file)
    else:
        print(str(result_file) + " already exists")

long.head(30)

######################
# COMPARE S02 V S06 ##
######################

########################################
# READIN PSSS3 top 3 peptide intensity #
########################################

psss3 = pd.read_pickle("wideint_full.pkl")

# Remove decoy
psss3 = psss3[psss3.index.get_level_values("decoy") == False]


psss3.loc[:, ("S02", slice(None))]
psss3.loc[:, ("S02", "R01")]

runs = ["R01", "R02", "R03", "R04", "R05"]
samples = ["S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10"]

# Create pickled pivot table
new_df = pd.DataFrame()
for i in runs:
    for j in samples:
        import time
        start = time.time()
        top3 = psss3.loc[:, (j, i)].groupby(level=[0]).apply(lambda x: x.nlargest(3)).reset_index(level=0, drop=True)
        end = time.time()
        print(end-start) #45s
        
        start = time.time()
        sumOfTop3 = top3.groupby(level=[0]).apply(lambda x: x.sum())
        end = time.time()
        print(start-end) # #9s
        
        if new_df.empty:
            new_df = new_df.append(pd.DataFrame(sumOfTop3))
        else:
            new_df = new_df.join(sumOfTop3)
            
df= new_df.reset_index()
df["protSpec"] = df.proteins
df["proteins"] =  df.protSpec.str.split("_", expand = True)[0]
df["species"] = df.protSpec.str.split("_", expand = True)[1]
df.drop(["protSpec"], axis = 1, inplace = True)
df = df.set_index(["proteins", "species"])

df.to_pickle("top3IntensePeptides_full.pkl")
df = pd.read_pickle("top3IntensePeptides_full.pkl")

# melt to sns format
df = df.div(df.sum(axis=1),axis=0) # Normalize df

df_reshape = pd.DataFrame()
for i in runs:
    for j in samples:
        print("concatinating " + i + " " + j)
        df_tmp = pd.DataFrame(df.loc[:,(j, i)])
        df_tmp["intensity"] = df_tmp.values
        df_tmp["sample"] = [df_tmp.columns[0][0]] * len(df_tmp)
        df_tmp["run"] = [df_tmp.columns[0][1]] * len(df_tmp)
        df_tmp = df_tmp.drop(j, axis=1, level=0)
        df_tmp.columns = df_tmp.columns.droplevel(level = 1)
        df_reshape = pd.concat([df_reshape, df_tmp])

prot = df_reshape.index.get_level_values("proteins").values
spec = df_reshape.index.get_level_values("species").values
intensity = df_reshape.loc[:,("intensity")].values
samples = df_reshape.loc[:,("sample")].values
runs = df_reshape.loc[:,("run")].values

df_reshaped = pd.DataFrame(np.array([prot, spec, intensity, samples, runs]).T, columns = ["proteins", "species", "intensities", "samples", "runs"])
df_reshaped = df_reshaped.convert_objects(convert_numeric = True)

df_reshaped.loc[df_reshaped['species'] == "CAEEL"].intensities.sum()

a_mix = np.array([0.5] * 10)
a_mix[0] += 0.0001
c_mix = np.array([0.5, 0.25, 0.125, 0.0625, 0.031, 0.0155, 0.008, 0.004, 0.002, 0.0])
h_mix = np.array([0.0, 0.25, 0.375, 0.4375, 0.469, 0.4845, 0.492, 0.496, 0.498, 0.5])
a_mix = a_mix/a_mix.sum()
c_mix = c_mix/c_mix.sum()
h_mix = h_mix/h_mix.sum()
mix = {'A thaliana':a_mix, 'C elegans': c_mix, 'H sapiens': h_mix}



fig, ax = plt.subplots(figsize=(13, 8))
#sns.violinplot(data=sampgroup.loc[sampgroup['Species'] == "H sapiens"], x="Sample",y="Expression", hue="Method",split=True)
#sns.boxplot(data=df_reshaped.loc[df_reshaped['species'] == "HUMAN"], x="samples",y="intensities",fliersize=0)
sns.boxplot(data=df_reshaped.loc[df_reshaped['proteins'] == "U4PC44"], x="samples",y="intensities",fliersize=0)
for i in range(len(h_mix)):
    ax.plot([i-.2, i+.2], [h_mix[i], h_mix[i]], '-', color='red', lw=5)
ax.set_ylim(0,0.26)
ax.set_ylabel("Fraction of protein total abundance per sample")
sns.despine(offset=10, trim=True)
plt.title("Homo sapiens")  
############################
# READIN PSSS3 PG.Quantity #
############################
        
df = pd.read_csv("../data/500-PSSS3-raw-reformatted_dropna_dropdup_decoy_nonShared_again.csv", sep ="\t")
#Treshold without removing NaN
FDR_treshold = 0.01
df = df[df["qvalue"].fillna(-1)<FDR_treshold]

decoy = df["decoy"].apply(lambda x: "decoy_" if x == True else "")
df["protein"] = decoy+df["protein"]

df_reform = getTriqlerOutputFormat_from_reformattedRawPSSS3(df, noDecoy = False)
psss3 = df_reform

psss3.reset_index(level = 0, inplace = True)
psss3 = addSpecieProteinColumns(psss3, startRun = "S01:S01_R01", endRun = "S10:S10_R05")

repl_cols = ["S01","S02","S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10"]

psss3 = renameTriqler2Vital_nonDetailed(psss3)

normalized_protein = psss3.div(psss3.sum(axis=1), axis = 0)

##################
# READIN Triqler #
##################

triqler = table2Df("proteins.1vs2.tsv", peptideSeperator = True)

# FDR_tresholding
triqler.sort_values(by = "protein_id_posterior_error_prob", inplace = True)
triqler["FDR"] = triqler["protein_id_posterior_error_prob"].expanding().mean()
triqler = triqler[triqler["FDR"] < 0.01]
triqler = triqler2Numeric(triqler)
triqler = getTriqlerDetailed(triqler)
triqler = renameTriqler2Vital(at_t)


long.head()
