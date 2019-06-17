#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 10:33:15 2019

@author: ptruong
"""



from triqlerParser import *
from triqlerProcessor import * 
from utils import * 
from parseReformattedPSSS3equDecoy import * # Needs to tidy and place functions in right place!
from corrPlot import * #correlation plotting


import pandas as pd
import numpy as np
from scipy import stats
import subprocess
%matplotlib inline
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set(font_scale=2)
sns.set_style("whitegrid", {'axes.grid' : False})#, 'grid.color': '.1', 'grid.linestyle': u'-'})

##########
# PSSS3 ##
##########
df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup_decoy.csv", sep ="\t")
df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup_decoy_nonShared_again.csv", sep ="\t")
FDR_treshold = 0.01
df = df[df["qvalue"].fillna(-1)<FDR_treshold]
decoy = df["decoy"].apply(lambda x: "decoy_" if x == True else "")
df["protein"] = decoy+df["protein"]

df_reform = getTriqlerOutputFormat_from_reformattedRawPSSS3(df, noDecoy = False)

df_reform_temp = df_reform
df_reform = df_reform_temp

#df_reform = df_reform.T.fillna(df_reform.mean(axis=1)).T
#df_reform = df_reform.T.fillna(df_reform.median(axis=1)).T
#df_reform = df_reform.T.fillna(df_reform.min(axis=1)).T
#df_reform = df_reform.T.fillna(df_reform.max(axis=1)).T
#df_reform = df_reform.dropna()
#psss3 = df_reform

df_reform = df_reform.fillna(df_reform.mean().mean())
df_reform = df_reform.fillna(df_reform.median().median())
df_reform = df_reform.fillna(df_reform.min().min())
df_reform = df_reform.fillna(df_reform.max().max())

df_reform = df_reform.dropna()
psss3 = df_reform

psss3.reset_index(level=0, inplace = True)
df_psss3 = addSpecieProteinColumns(psss3, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
#df_psss3 = df_psss3.dropna()

at_s = df_psss3[df_psss3["specie"] == "ARATH"]
ce_s = df_psss3[df_psss3["specie"] == "CAEEL"]
hs_s = df_psss3[df_psss3["specie"] == "HUMAN"]

at_s = renameTriqler2Vital_nonDetailed(at_s)
ce_s = renameTriqler2Vital_nonDetailed(ce_s)
hs_s = renameTriqler2Vital_nonDetailed(hs_s)

#########################
##### pure data sets ####
#########################

psss3 = pd.read_table('500-PSSS3-precursor_Report_Lukas.xls')
psss3_2 = pd.read_table('500-PSSS3-precursor_Report_new.xls')
psss3_3 = pd.read_table('500-PSSS3-precursor_Report.xls')

def spectronaut_extract(species):
    ss = psss3[psss3['PG.Organisms'] == species]
    ss = ss.groupby(['PG.ProteinAccessions', 'R.Replicate', 'R.Condition'])
    ss = ss.apply(lambda x: x['PG.Quantity'].iloc[0]).unstack().unstack()
    # ss = ss.apply(lambda row: row.where(row>0,np.min(row[row>0])/2.0),axis=1)
    ss.columns = ss.columns.get_level_values(0)
    ss.columns = [c.split("-")[-1] for c in ss.columns]
    return ss

at_s = spectronaut_extract('Arabidopsis thaliana')
ce_s = spectronaut_extract('Caenorhabditis elegans OX=6239')
hs_s = spectronaut_extract('Homo sapiens')

###########
# triqler #
###########


triqler = table2Df("proteins.1vs2_noShared.tsv", peptideSeperator = True)
#triqler = table2Df("proteins.1vs2.tsv", peptideSeperator = True)
triqler = triqler.dropna()
triqler = triqler2Numeric(triqler)
triqler = getTriqlerDetailed(triqler)


triqler = table2Df("PSSS3_triqler_output_proteins.1vs10.tsv", peptideSeperator = True)
triqler = triqler2Numeric(triqler)
triqler = getTriqlerDetailed(triqler)
#triqler = triqler[triqler["q_value"] < 0.01] #FDR filter

# FDR on posterior_error_prob
triqler.sort_values(by = "posterior_error_prob", inplace = True)
triqler["FDR"] = triqler["posterior_error_prob"].expanding().mean()
triqler = triqler[triqler["FDR"] < 0.01]

# FDR on protein_id_posterior_error_prob
triqler.sort_values(by = "protein_id_posterior_error_prob", inplace = True)
triqler["FDR"] = triqler["protein_id_posterior_error_prob"].expanding().mean()
triqler = triqler[triqler["FDR"] < 0.01]


#############
# Triqler ###
#############

#FDR = getTriqler_FDR(triqler, 0.05)

PEP = findFDR(triqler, target_FDR = .01, lower_PEP = 0.0, upper_PEP =1)
triqler = triqler[triqler["protein_id_posterior_error_prob"] < PEP]
decoyMatches = sum(triqler[triqler["protein_id_posterior_error_prob"] < PEP].decoy)


#df_triqler = addSpecieProteinColumns(triqler, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
df_triqler = triqler
runs = getRuns(triqler, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
conditionRuns = getConditionRuns(runs)
unique_run = conditionRuns.run.unique()
run_id = getConditionRunId(triqler, unique_run[0])

# All runs
at_t = df_triqler[df_triqler["specie"] == "ARATH"]
ce_t = df_triqler[df_triqler["specie"] == "CAEEL"]
hs_t = df_triqler[df_triqler["specie"] == "HUMAN"]

# DROP species and drop peptides <------------------------
at_t = renameTriqler2Vital(at_t, diff_exp = 1.0)
ce_t = renameTriqler2Vital(ce_t, diff_exp = 1.0)
hs_t = renameTriqler2Vital(hs_t, diff_exp = 1.0)


at_t = pd.read_csv('diffacto_A_thaliana.tsv', delimiter='\t')
ce_t = pd.read_csv('diffacto_C_elegans.tsv', delimiter='\t')
hs_t = pd.read_csv('diffacto_H_sapiens.tsv', delimiter='\t')
#at_dq = at_d[at_d["S/N"]>-20]
at_t = at_t[at_t["S/N"]>-20]
hs_t = hs_t[hs_t["S/N"]>-20]
ce_t = ce_t[ce_t["S/N"]>-20]




protein_count = pd.DataFrame(columns=['species', 'method'])
protein_count['species'] = ['A_thaliana']*len(at_s) + ['C_elegans']*len(ce_s) + ['H_sapiens']*len(hs_s) + \
                           ['A_thaliana']*len(at_t) + ['C_elegans']*len(ce_t) + ['H_sapiens']*len(hs_t) 
protein_count['method']  = ['Spectronaut']*(len(at_s)+len(ce_s)+len(hs_s)) + \
                           ['Triqler']   *(len(at_t)+len(ce_t)+len(hs_t))
fig, ax = plt.subplots(figsize=(15, 9))
sns.countplot(x='species', hue='method', data=protein_count)
sns.despine(offset=10, trim=True);
plt.show()

col_order = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
#col_order = ["S01", "S02"]

def countRightlySorted(df,col_order=col_order):
    df["lt"] = [0.]*len(df)
    df["gt"] = [0.]*len(df)
    df["eq"] = [0.]*len(df)
    for left in range(0,len(col_order)-1):
        for right in range(left+1,len(col_order)):
            df["lt"] += (df[col_order[left]]<df[col_order[right]])
            df["gt"] += (df[col_order[left]]>df[col_order[right]])
            df["eq"] += (df[col_order[left]]==df[col_order[right]])
    df["tot"] = df["lt"] + df["gt"] + df["eq"]
    
    
from scipy.stats.mstats import gmean
missingValue=np.array([0.])


sampgroup_at_s,sampgroup_hs_s,sampgroup_ce_s = pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
for s in col_order:
    sampgroup_at_s[s] = at_s[s].sum(1)
    sampgroup_ce_s[s] = ce_s[s].sum(1)
    sampgroup_hs_s[s] = hs_s[s].sum(1)
sampgroup_at_s['Species'] = ['A thaliana']*len(at_s)
sampgroup_ce_s['Species'] = ['C elegans']*len(ce_s)
sampgroup_hs_s['Species'] = ['H sapiens']*len(hs_s)
sampgroup_at_s['Method'] = ['Spectronaut']*len(at_s)
sampgroup_ce_s['Method'] = ['Spectronaut']*len(ce_s)
sampgroup_hs_s['Method'] = ['Spectronaut']*len(hs_s)

sampgroup_spec = pd.concat([sampgroup_at_s,sampgroup_ce_s,sampgroup_hs_s])


repl_cols_d = [c for c in at_t.columns if "-" in c]

#DIFFACTO:

at_t['Species'] = ['A thaliana']*len(at_t)
at_t['Method'] = ['Triqler']*len(at_t)
ce_t['Species'] = ['C elegans']*len(ce_t)
ce_t['Method'] = ['Triqler']*len(ce_t)
hs_t['Species'] = ['H sapiens']*len(hs_t)
hs_t['Method'] = ['Triqler']*len(hs_t)

diffac = pd.concat([at_t,ce_t,hs_t])

sampgroup = diffac.loc[:,['Species','Method']]
for s in col_order:
    sampgroup[s] = gmean(diffac.loc[:,diffac.columns.str.contains(s)],axis=1)

sampgroup=sampgroup.append(sampgroup_spec)
#missingValue=0
sampgroup.where((pd.notnull(sampgroup)), np.NaN,inplace=True)
#sampgroup.dropna(inplace=True)
        
#sampgroup.where((pd.notnull(sampgroup)), missingValue, inplace=True)

sampgroup['Protein'] = range(len(sampgroup))

unmeltsamp = sampgroup.copy()
#sampgroup=pd.melt(sampgroup, id_vars=['Method', 'Species','Protein'], value_vars=col_order)
#sampgroup.rename(columns={'variable':'Sample','value':'Expression'}, inplace=True)

countRightlySorted(unmeltsamp)


def plotSorted(species,compareOperator,figname=None):
    fig, ax = plt.subplots(figsize=(10, 6))
    for meth in unmeltsamp['Method'].unique():
        sns.distplot(unmeltsamp[(unmeltsamp['Method']==meth) & (unmeltsamp['Species']==species)][compareOperator],
#             hist_kws={"cumulative":-1, "linewidth": 3, "alpha": 1,"histtype": "step"},
             hist_kws={"linewidth": 3, "alpha": 0.5},
             bins=range(10,46),kde=False, ax=ax, label=meth)
    ax.legend()
    ax.set_title(species + " (" + compareOperator + ")")
    #fig.legend(handles=handles, labels=labels)
#    ax.set_xlim(25,45)
#    ax.set_ylim(0,150)
    sns.despine(offset=10, trim=True)
    ax.set_ylabel("Number of proteins")
    ax.set_xlabel("Correctly sorted")
    if figname:
        fig.savefig(figname)
        
        
#plotSorted('H sapiens',"lt","Hsap_lt.png")
#plotSorted('C elegans',"gt","Cele_gt.png")
plotSorted('H sapiens',"lt")
plotSorted('C elegans',"gt")
plotSorted('A thaliana',"eq")


plotSorted('H sapiens',"gt")
plotSorted('C elegans',"gt")
plotSorted('A thaliana',"gt")

plotSorted('A thaliana',"lt")
plotSorted('C elegans',"eq")
plotSorted('C elegans',"tot")



