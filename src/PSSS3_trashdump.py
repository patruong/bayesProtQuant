#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:44:33 2019

@author: ptruong
"""

from triqlerParser import *
from triqlerProcessor import * 
from utils import * 

import os
import pandas as pd
import numpy as np
from scipy import stats
import subprocess
%matplotlib inline
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set(font_scale=2)
sns.set_style("whitegrid", {'axes.grid' : False})#, 'grid.color': '.1', 'grid.linestyle': u'-'})


############
## triqler #
############

triqler = table2Df("PSSS3_triqler_output_proteins.1vs2.tsv", peptideSeperator = True)
#triqler2 = table2Df("PSSS3_triqler_output_proteins.1vs3.tsv", peptideSeperator = True)
#triqler3 = table2Df("PSSS3_triqler_output_proteins.2vs3.tsv", peptideSeperator = True)

df_triqler = addSpecieProteinColumns(triqler, startRun = "S01:S01_R01", endRun = "S10:S10_R05")

runs = getRuns(triqler, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
conditionRuns = getConditionRuns(runs)
unique_run = conditionRuns.run.unique()
run_id = getConditionRunId(triqler, unique_run[0])

# All runs
at_t = df_triqler[df_triqler["specie"] == "ARATH"]
ce_t = df_triqler[df_triqler["specie"] == "CAEEL"]
hs_t = df_triqler[df_triqler["specie"] == "HUMAN"]

# Specfic runs
df_triqler[["id", "protein", "specie"] + run_id.tolist()]

############
## psss3 ###
############

def spectronaut_extract(species):
    ss = psss3_raw[psss3_raw['PG.Organisms'] == species]
    ss = ss.groupby(['PG.ProteinAccessions', 'R.Replicate', 'R.Condition'])
    ss = ss.apply(lambda x: x['PG.Quantity'].iloc[0]).unstack().unstack()
    # ss = ss.apply(lambda row: row.where(row>0,np.min(row[row>0])/2.0),axis=1)
    ss.columns = ss.columns.get_level_values(0)
    ss.columns = [c.split("-")[-1] for c in ss.columns] #<--- you lose which run it belongs too lol.
    return ss

psss3_raw = pd.read_table('500-PSSS3-precursor_Report.xls')
psss3 = pd.read_table('500-PSSS3-precursor_Report.proteins.tsv')

def addId2psss3(psss3):
    """
    Input:
        psss file
    """
    firstSplit = psss3["R.FileName"].str.split("_", n = 6, expand = True)
    secondSplit = firstSplit[2].str.split("-", n = 3, expand = True)
    
    run = firstSplit[4]
    samples = secondSplit[2]
    identity = samples + ":" + samples + "_" + run
    psss3["id"] = identity
    return psss3

#psss3 = addId2psss3(psss3)

df_psss3 = addSpecieProteinColumns(psss3, startRun = "S01:S01_R01", endRun = "S10:S10_R05")

runs = getRuns(psss3, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
conditionRuns = getConditionRuns(runs)
unique_run = conditionRuns.run.unique()
run_id = getConditionRunId(psss3, unique_run[0])

# All runs
at_s = df_psss3[df_psss3["specie"] == "ARATH"]
ce_s = df_psss3[df_psss3["specie"] == "CAEEL"]
hs_s = df_psss3[df_psss3["specie"] == "HUMAN"]



at_s = renameTriqler2Vital(at_s)
ce_s = renameTriqler2Vital(ce_s)
hs_s = renameTriqler2Vital(hs_s)

# Specfic runs
df_psss3[["id", "protein", "specie"] + run_id.tolist()]


# All runs
at_s = psss3[psss3["PG.Organisms"] == 'Arabidopsis thaliana']
ce_s = psss3[psss3["PG.Organisms"] == 'Caenorhabditis elegans OX=6239']
hs_s = psss3[psss3["PG.Organisms"] == 'Homo sapiens']

at_s_raw = spectronaut_extract('Arabidopsis thaliana')
ce_s_raw = spectronaut_extract('Caenorhabditis elegans OX=6239')
hs_s_raw = spectronaut_extract('Homo sapiens')



#at_s_raw = psss3_raw[psss3_raw["PG.Organisms"] == 'Arabidopsis thaliana']
#ce_s_raw = psss3_raw[psss3_raw["PG.Organisms"] == 'Caenorhabditis elegans OX=6239']
#hs_s_raw = psss3_raw[psss3_raw["PG.Organisms"] == 'Homo sapiens']



# Specific runs



###################################################
## Count detected proteins per species and method #
###################################################
protein_count = pd.DataFrame(columns = ['species', 'method'])
protein_count['species'] = (['ARATH']*len(at_t) + ['CAEEL']*len(ce_t) + ["HUMAN"]*len(hs_t) + 
             ['ARATH']*len(at_s) + ['CAEEL']*len(ce_s) + ["HUMAN"]*len(hs_s))
protein_count['method'] = (['Triqler']*len(at_t) + ['Triqler']*len(ce_t) + ['Triqler']*len(hs_t) + 
             ['Spectronaut']*len(at_s) + ['Spectronaut']*len(ce_s) + ['Spectronaut']*len(hs_s))

fig, ax = plt.subplots(figsize=(15, 9))
sns.countplot(x='species', hue='method', data=protein_count)
sns.despine(offset=10, trim=True);
plt.show()

################################
## Per-protein-total-abundance #
################################

#SPECTRONAUT:
#per-protein-total-abundance A_thaliana Spectronaut:
pp_at_s = at_s.sum(axis=0).values / len(at_s)
#per-protein-total-abundance C_elegans Spectronaut:
pp_ce_s = ce_s.sum(axis=0).values / len(ce_s)
#per-protein-total-abundance H_sapiens Spectronaut:
pp_hs_s = hs_s.sum(axis=0).values / len(hs_s)

# IS THIS REALLY PER PROTEIN ABUNDANCE?

#DIFFACTO:
#per-protein-total-abundance A_thaliana Diffacto:
pp_at_d = at_d[repl_cols].sum(axis=0).values / len(at_d)
#per-protein-total-abundance C_elegans Diffacto:
pp_ce_d = ce_d[repl_cols].sum(axis=0).values / len(ce_d)
#per-protein-total-abundance H_sapiens Diffacto:
pp_hs_d = hs_d[repl_cols].sum(axis=0).values / len(hs_d)

######################################
# Normalize to per-protein abundance #
######################################

at_s = at_s.drop(["protein", "id", "specie"], axis = 1)
ce_s = ce_s.drop(["protein", "id", "specie"], axis = 1)
hs_s = hs_s.drop(["protein", "id", "specie"], axis = 1)

#normalized_protein_at_s = at_s.div(at_s.sum(axis=1), axis=0)
normalized_protein_at_s = at_s / at_s.sum()
normalized_protein_ce_s = ce_s.div(ce_s.sum(axis=1), axis=0)
normalized_protein_ce_s = ce_s / ce_s.sum()
normalized_protein_hs_s = hs_s.div(hs_s.sum(axis=1), axis=0)
normalized_protein_hs_s = hs_s / hs_s.sum()

sampgroup_at_s,sampgroup_hs_s,sampgroup_ce_s = pd.DataFrame(),pd.DataFrame(),pd.DataFrame()

for s in repl_cols_s:
    sampgroup_at_s[s] = normalized_protein_at_s[s].sum(1)
    sampgroup_ce_s[s] = normalized_protein_ce_s[s].sum(1)
    sampgroup_hs_s[s] = normalized_protein_hs_s[s].sum(1)
    
sampgroup_at_s['Species'] = ['A thaliana'] *len(normalized_protein_at_s)
sampgroup_ce_s['Species'] = ['C elegans']  *len(normalized_protein_ce_s)
sampgroup_hs_s['Species'] = ['H sapiens']  *len(normalized_protein_hs_s)
sampgroup_at_s['Method']  = ['Spectronaut']*len(normalized_protein_at_s)
sampgroup_ce_s['Method']  = ['Spectronaut']*len(normalized_protein_ce_s)
sampgroup_hs_s['Method']  = ['Spectronaut']*len(normalized_protein_hs_s)

sampgroup_spec = pd.concat([sampgroup_at_s,sampgroup_ce_s,sampgroup_hs_s])

######################################

repl_cols_s = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]

normalized_protein_at_s = at_s.div(at_s.sum(axis=1), axis=0)
normalized_protein_ce_s = ce_s.div(ce_s.sum(axis=1), axis=0)
normalized_protein_hs_s = hs_s.div(hs_s.sum(axis=1), axis=0)

sampgroup_at_s,sampgroup_hs_s,sampgroup_ce_s = pd.DataFrame(),pd.DataFrame(),pd.DataFrame()

for s in repl_cols_s:
    sampgroup_at_s[s] = normalized_protein_at_s[s].sum(1)
    sampgroup_ce_s[s] = normalized_protein_ce_s[s].sum(1)
    sampgroup_hs_s[s] = normalized_protein_hs_s[s].sum(1)
    
sampgroup_at_s['Species'] = ['A thaliana'] *len(normalized_protein_at_s)
sampgroup_ce_s['Species'] = ['C elegans']  *len(normalized_protein_ce_s)
sampgroup_hs_s['Species'] = ['H sapiens']  *len(normalized_protein_hs_s)
sampgroup_at_s['Method']  = ['Spectronaut']*len(normalized_protein_at_s)
sampgroup_ce_s['Method']  = ['Spectronaut']*len(normalized_protein_ce_s)
sampgroup_hs_s['Method']  = ['Spectronaut']*len(normalized_protein_hs_s)

sampgroup_spec = pd.concat([sampgroup_at_s,sampgroup_ce_s,sampgroup_hs_s])

sampgroup = sampgroup_spec

sampgroup.where((pd.notnull(sampgroup)), np.NaN,inplace=True)
sampgroup.dropna(inplace=True)

#nn = sampgroup[repl_cols_s][sampgroup[repl_cols_s] > 0.0].min(numeric_only=True,axis=1)/2.0
#nnn = pd.concat([nn]*sampgroup.shape[1],axis=1)
#sampgroup = sampgroup.where(sampgroup!=0.0, nnn)
        

sampgroup['Protein'] = range(len(sampgroup))

unmeltsamp = sampgroup.copy()
sampgroup=pd.melt(sampgroup, id_vars=['Method', 'Species','Protein'], value_vars=repl_cols_s)
sampgroup.rename(columns={'variable':'Sample','value':'Expression'}, inplace=True)

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
sns.boxplot(data=sampgroup.loc[sampgroup['Species'] == "H sapiens"], x="Sample",y="Expression", hue="Method",fliersize=0)
for i in range(len(h_mix)):
    ax.plot([i-.2, i+.2], [h_mix[i], h_mix[i]], '-', color='red', lw=5)
ax.set_ylim(0,0.15)
ax.set_ylabel("Fraction of protein total abundance per sample")
sns.despine(offset=10, trim=True)
plt.title("Homo sapiens")
fig.savefig("sapiens_sample.png")

