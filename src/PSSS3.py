#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 11:07:34 2019

@author: ptruong
"""


from triqlerParser import *
from triqlerProcessor import * 
from utils import * 
from parseReformattedPSSS3equDecoy import * # Needs to tidy and place functions in right place!
from corrPlot import * #correlation plotting
#from diffExp import * 
from qvalue import *

import os
import pandas as pd
import numpy as np
from scipy import stats
import scipy as sp
import subprocess
%matplotlib inline
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set(font_scale=2)
sns.set_style("whitegrid", {'axes.grid' : False})#, 'grid.color': '.1', 'grid.linestyle': u'-'})


############
## triqler #
############


triqler = table2Df("proteins.1vs2.tsv", peptideSeperator = True)
triqler = table2Df("proteins.1vs2_multiPriorCorrect.tsv", peptideSeperator = True)
triqler = table2Df("multiPrior_groupNorm_version2_beforeChange_imputedDiffsGroups.tsv", peptideSeperator = True)


triqler = table2Df("proteins.1vs2_noShared.tsv", peptideSeperator = True)
triqler = table2Df("proteins.1vs2_filteredOn_PGQvalue_EGQvalue_PGCscore.tsv", peptideSeperator = True)
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



#####


#####
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



def meltTriqler2SNSformat(df, specie, diff_exp_prob = 1.0, normalize = True):
    """
    Takes at_t, hs_t, ce_t and melt them to SNS df.
    """
    if normalize == True:
        df.iloc[:,6:-5] = df.iloc[:, 6:-5].div(df.iloc[:, 6:-5].sum(axis=1), axis=0)
    diff_exp_prob = str(diff_exp_prob)
    df = df.drop(["id", "specie", "q_value", "posterior_error_prob",
                      "num_peptides", "protein_id_posterior_error_prob",
                      "log2_fold_change", "diff_exp_prob_"+diff_exp_prob,
                      "peptides", "specie", "decoy", "id"], axis = 1).melt("protein")
    df = df.rename(index = str, columns={"value":"Expression", "variable":"id", "protein":"Protein"})
    #df["protein"] + "~" +  df["variable"]
    df["Method"] = ["Triqler"] * len(df)
    df['Species'] = [specie] *len(df)
    df["Sample"] = df["id"].str.split(":", expand = True)[0]
    return df

sampgroup_at_t = meltTriqler2SNSformat(at_t, "A thaliana")
sampgroup_hs_t = meltTriqler2SNSformat(hs_t, "H sapiens")
sampgroup_ce_t = meltTriqler2SNSformat(ce_t, "C elegans")
sampgroup_spec_t = pd.concat([sampgroup_at_t,sampgroup_ce_t,sampgroup_hs_t])

# DROP species and drop peptides <------------------------
at_t = renameTriqler2Vital(at_t)
ce_t = renameTriqler2Vital(ce_t)
hs_t = renameTriqler2Vital(hs_t)

at_t = 10**at_t
ce_t = 10**ce_t
hs_t = 10**hs_t

# Drop S01 or S10
at_t = at_t.drop(["S01", "S10"], axis = 1)
ce_t = ce_t.drop(["S01", "S10"], axis = 1)
hs_t = hs_t.drop(["S01", "S10"], axis = 1)


at_s[at_s.index == "P51413"].T.plot(title = "spectronaut")
at_t[at_t.index == "P51413"].T.plot(title = "triqler")

ce_t[ce_t.index == "Q23551"].T.plot()

ce_s[ce_s.index == "P48053"].T.plot(title = "spectronaut")


ce_t["S02"].idxmin()
ce_t["S02"].min()
ce_t[ce_t.index == "Q9N337"].T.plot()


##########################
## psss3  Matthew data ###
##########################
#psss3_raw = pd.read_table("S500-PSSS3-equ decoy_Report.xls")
psss3 = pd.read_table('500-PSSS3-precursor_Report.proteins.tsv')
df_psss3 = addSpecieProteinColumns(psss3, startRun = "S01:S01_R01", endRun = "S10:S10_R05")

at_s = df_psss3[df_psss3["specie"] == "ARATH"]
ce_s = df_psss3[df_psss3["specie"] == "CAEEL"]
hs_s = df_psss3[df_psss3["specie"] == "HUMAN"]

at_s = renameTriqler2Vital_nonDetailed(at_s)
ce_s = renameTriqler2Vital_nonDetailed(ce_s)
hs_s = renameTriqler2Vital_nonDetailed(hs_s)

########################
## psss3  unfiltered ###
########################

df = pd.read_csv("500-PSSS3-equ decoy_Report.xls", sep ="\t")

df = parsePSSS3equDecoy(filename = '500-PSSS3-equ decoy_Report_nonShared.xls')
df = parsePSSS3equDecoy(filename = '500-PSSS3-equ decoy_Report_filteredOn_PGQvalue_EGQvalue_PGCscore.xls', index = True)
#df = parsePSSS3equDecoy(filename = '500-PSSS3-equ decoy_Report.xls') # SLOW
#df = parsePSSS3equDecoy("500-PSSS3-raw-reformatted_dropna_dropdup_decoy.csv") # SLOW
#df2 = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t")
#df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup_droppedCols.csv", sep = "\t")
df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup_decoy.csv", sep ="\t")
df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup_decoy_sharedCAEEL.csv", sep ="\t")
df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup_decoy_nonShared.csv", sep ="\t")
df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup_decoy_nonShared_again.csv", sep ="\t")
df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup_decoy_filteredOn_PQValue_EGQvalue_PGCscore.csv", sep ="\t")

#df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t")
df = pd.read_csv("500-PSSS3-precursor_Report.xls", sep = "\t")
df = pd.read_csv("500-PSSS3-precursor_Report_reformatted.xls", sep = "\t")


#df.fillna(-1, inplace = True)
#Treshold without removing NaN
FDR_treshold = 0.01
df = df[df["qvalue"].fillna(-1)<FDR_treshold]
#df.replace(-1, np.NaN, inplace = True)

#########################
# COUNT NAN TABLE #######
#########################
def rawSpecies(df, specie = "Caenorhabditis elegans OX=6239"):
    """
    Slices raw data (500-PSSS3-equ decoy_Report.xls) into 
    """
    species = df[df["PG.Organisms"] == specie]
    return species

def samplePercNaN(df, specie = "Caenorhabditis elegans OX=6239", sample = "S01"):
    """
    Get sample percentage NaNs.
    """
    spec = rawSpecies(df = df, specie = specie)
    nonDecoy = spec[spec["EG.IsDecoy"] == False]    
    sample = sample
    sampleDat = nonDecoy[nonDecoy["R.Condition"].str[-3:] == sample]    
    perc_NaN = sampleDat["PG.Quantity"].isna().sum() / (sampleDat["PG.Quantity"].count() + sampleDat["PG.Quantity"].isna().sum())
    return perc_NaN

def allSamplesPercNaN(df, specie = "Caenorhabditis elegans OX=6239"):
    repl_cols = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
    percNaN_list = []
    for i in repl_cols:
        percNaN = samplePercNaN(df, specie = specie, sample = i)
        percNaN_list.append(percNaN)
    percentage_df = pd.DataFrame(np.array([repl_cols, np.around(percNaN_list, 4)]).T,  columns = ["sample", "percentage_NaN"])
    return percentage_df

def computeRawPercNaNTable(df):
    """
    Compute perc NaN for all species in the raw data (df = pd.read_csv("500-PSSS3-equ decoy_Report.xls", sep ="\t"))
    """
    species = df["PG.Organisms"].unique()
    percentageList = []
    specieList = []
    for i in species:
        percentages = allSamplesPercNaN(df, specie = i)
        percentageList.append(percentages)
        specieList.append(i)
    
    valuesList = []
    samples = percentageList[0]["sample"]
    for i in range(len(percentageList)):
        valuesList.append(percentageList[i]["percentage_NaN"].values)
    percTable = pd.DataFrame(np.array([samples.values] + valuesList).T, columns = ["sample"] + list(species))
    return percTable


def show_values_on_bars(axs):
    def _show_on_single_plot(ax):        
        for p in ax.patches:
            _x = p.get_x() + p.get_width() / 2
            _y = p.get_y() + p.get_height()
            value = '{:.3f}'.format(p.get_height())
            ax.text(_x, _y, value, ha="center") 

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _show_on_single_plot(ax)
    else:
        _show_on_single_plot(axs)
        
percTable = computeRawPercNaNTable(df)


meltedPercTable = percTable.melt("sample")
meltedPercTable = meltedPercTable.convert_objects(convert_numeric=True)

        
fig, ax = plt.subplots(figsize=(13*2.5, 8*2.5))
#sns.violinplot(data=sampgroup.loc[sampgroup['Species'] == "H sapiens"], x="Sample",y="Expression", hue="Method",split=True)
sns.barplot(data=meltedPercTable, x="sample",y="value", hue="variable")
ax.set_ylim(0,1.5)
ax.set_ylabel("Fraction of NaN")
#ax1= g.axes[0]

ax.axhline(1, ls='--', color = "r")
show_values_on_bars(ax)

sns.despine(offset=10, trim=True)
plt.title("Homo sapiens")  
fig.savefig("sapiens_sample.png")

#nonDecoy["FG.NormalizedMS2PeakArea"].isna().sum()# Peptide precursor quantity
#nonDecoy["PG.Cscore"].isna().sum()
#nonDecoy["EG.Cscore"].isna().sum()



#########################
## INTERSECTION METHOD ##
#########################
"""
Choose only intersecting proteins.
"""

overlappingProteins = list(set(df["protein"]).intersection(set(triqler["id"])))

# complete subset of proteins
df = df[df["protein"].isin(overlappingProteins)]
triqler = triqler[triqler["id"].isin(overlappingProteins)]

# sampled subset of proteins
import random
n_samples = 3500
sampleProteins = random.sample(overlappingProteins, n_samples)
df = df[df["protein"].isin(sampleProteins)]
triqler = triqler[triqler["id"].isin(sampleProteins)]

#####################

decoy = df["decoy"].apply(lambda x: "decoy_" if x == True else "")
df["protein"] = decoy+df["protein"]

df_reform = getTriqlerOutputFormat_from_reformattedRawPSSS3(df, noDecoy = False)

#####################
# Global Imputation #
#####################

df_reform = df_reform.fillna(df_reform.mean().mean())
df_reform = df_reform.fillna(df_reform.median().median())
df_reform = df_reform.fillna(df_reform.min().min())
df_reform = df_reform.fillna(df_reform.max().max())

df_reform = df_reform.dropna()
psss3 = df_reform.copy()


#######################
# row-wise imputation ###
#######################

#df_reform = df_reform.T
#df_reform = df_reform.fillna(df_reform.mean())
#df_reform = df_reform.fillna(df_reform.median())
#df_reform = df_reform.fillna(df_reform.min())
#df_reform = df_reform.fillna(df_reform.max())
#df_reform = df_reform.T

df_reform = df_reform.T.fillna(df_reform.mean(axis=1)).T
df_reform = df_reform.T.fillna(df_reform.median(axis=1)).T
df_reform = df_reform.T.fillna(df_reform.min(axis=1)).T
df_reform = df_reform.T.fillna(df_reform.max(axis=1)).T

df_reform = df_reform.dropna()
psss3 = df_reform

## TEST
#df_t = df_reform[df_reform.index == "V6CM07_CAEEL"]
#prettyPrint(df_reform[df_reform.index == "V6CM07_CAEEL"])
#df_reform[df_reform.index == "V6CM07_CAEEL"]
#df_t.T.fillna(df_t.mean(axis=1))

#df_t2 = df_reform[df_reform.index == "V6CLU8_CAEEL"]
#prettyPrint(df_reform[df_reform.index == "V6CLU8_CAEEL"])
#df_t2.T.fillna(df_t2.mean(axis=1))

#df_t3 = df_reform[df_reform.index == "V6CLV0_CAEEL"]
#prettyPrint(df_reform[df_reform.index == "V6CLV0_CAEEL"])
#df_t3.T.fillna(df_t3.mean(axis=1))
##collisionEnergy=

##############################
# sample based imputation ####
##############################

df_multiIdx = triqlerFormat2multiIndex(df_reform)
import time
start = time.time()
df_multiIdx = multiIdxFillna_min(df_multiIdx) #mean imputation #<--- IMPUTATION STEP
end = time.time()
tot = end - start
print(tot)

df_multiIdx = df_multiIdx.dropna(thresh = len(df_multiIdx.iloc[0]))

psss3 = multiIndex2TriqlerFormat(df_multiIdx)
###########################
### END OF IMPUTATION #####
###########################

psss3.reset_index(level=0, inplace = True)
df_psss3 = addSpecieProteinColumns(psss3, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
#df_psss3 = df_psss3.dropna()

at_s = df_psss3[df_psss3["specie"] == "ARATH"]
ce_s = df_psss3[df_psss3["specie"] == "CAEEL"]
hs_s = df_psss3[df_psss3["specie"] == "HUMAN"]

# ToDo remove noise from at_s, ce_s, hs_s IF not enough samples detected in runs remove protein.

# NORMALIZATION
def meltSpectronaut2SNSformat(df, specie, normalize = True):
    """
    Takes at_s, hs_s, ce_s and melt them to SNS df.
    """
    if normalize == True:
        df.iloc[:, 3:] = df.iloc[:, 3:].div(df.iloc[:, 3:].sum(axis=1), axis=0)
    df = df.drop(["id", "specie"], axis = 1).melt("protein")
    df = df.rename(index = str, columns={"value":"Expression", "variable":"id", "protein":"Protein"})
    #df["protein"] + "~" +  df["variable"]
    df["Method"] = ["Spectronaut"] * len(df)
    df['Species'] = [specie] *len(df)
    df["Sample"] = df["id"].str.split(":", expand = True)[0]
    return df

sampgroup_at_s = meltSpectronaut2SNSformat(at_s, "A thaliana")
sampgroup_hs_s = meltSpectronaut2SNSformat(hs_s, "H sapiens")
sampgroup_ce_s =   meltSpectronaut2SNSformat(ce_s, "C elegans")
sampgroup_spec_s = pd.concat([sampgroup_at_s,sampgroup_ce_s,sampgroup_hs_s])

at_s.iloc[:,3:] = np.log10(at_s.iloc[:,3:])
ce_s.iloc[:,3:] = np.log10(ce_s.iloc[:,3:])
hs_s.iloc[:,3:] = np.log10(hs_s.iloc[:,3:])

repl_cols = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]

# Dropped S01 and S10
repl_cols = ["S02","S03","S04","S05","S06","S07","S08","S09"]

'''
count_table = (at_s[at_s["id"].str.slice(0,5) != "decoy"].isna().sum()/at_s[at_s["id"].str.slice(0,5) != "decoy"].count())

repl_list = []
perc_list = []
for i in repl_cols:
    repl_list.append(i)
    perc_avg = count_table[count_table.index.str.slice(0,3) == i].sum() / 5
    perc_list.append(perc_avg)
count_df = pd.DataFrame(np.array([repl_list, np.around(perc_list, 4)]).T, columns = ["sample", "percentage_NaN"])

count_table = (ce_s[ce_s["id"].str.slice(0,5) != "decoy"].isna().sum()/ce_s[ce_s["id"].str.slice(0,5) != "decoy"].count())

ce_s[ce_s["id"].str.slice(0,5) != "decoy"].count()

repl_list = []
perc_list = []
for i in repl_cols:
    repl_list.append(i)
    perc_avg = count_table[count_table.index.str.slice(0,3) == i].sum() / 5
    perc_list.append(perc_avg)
count_df = pd.DataFrame(np.array([repl_list, np.around(perc_list, 4)]).T,  columns = ["sample", "percentage_NaN"])

count_table = (hs_s[hs_s["id"].str.slice(0,5) != "decoy"].isna().sum()/hs_s[hs_s["id"].str.slice(0,5) != "decoy"].count())

repl_list = []
perc_list = []
for i in repl_cols:
    repl_list.append(i)
    perc_avg = count_table[count_table.index.str.slice(0,3) == i].sum() / 5
    perc_list.append(perc_avg)
count_df = pd.DataFrame(np.array([repl_list, np.around(perc_list, 4)]).T,  columns = ["sample", "percentage_NaN"])
'''

at_s = renameTriqler2Vital_nonDetailed(at_s)
ce_s = renameTriqler2Vital_nonDetailed(ce_s)
hs_s = renameTriqler2Vital_nonDetailed(hs_s)


####################################################
# DATA EXPLORATION ... NOTE DO NOT IMPUTE BEFORE ###
####################################################
all_t = renameTriqler2Vital(df_triqler, diff_exp = 1.0).dropna()
all_s = renameTriqler2Vital_nonDetailed(df_psss3)
all_s

samples = all_s.columns.unique()

sample = "S01"

all_s[all_s["S01"].isna().sum(axis = 1) == 5]




# Get diffExp probabilities TRIQLER

diffExpDir = "triqlerOutput_noShared_largeScaleOptimizations_minSamp15"
triqler = table2Df(diffExpDir+"/" + "proteins.3vs7.tsv", peptideSeperator = True)
triqler = triqler2Numeric(triqler)
triqler = getTriqlerDetailed(triqler)

diff_exp_treshold = str(1.0)
FDR_treshold = 0.05
diff_t = triqler["diff_exp_prob_"+diff_exp_treshold]
diff_t = pd.to_numeric(diff_t)
diff_t = (diff_t > (1-FDR_treshold))
count_t  = diff_t.sum()

# Get diffExp spectronaut
all_s = renameTriqler2Vital_nonDetailed(df_psss3)
#val1 = hs_s["S03"]
#val2 = hs_s["S07"]

val1 = all_s["S03"]
val2 = all_s["S07"]

diff_s = diffExp(val1, val2)
diff_s = pd.DataFrame(diff_s)
diff_s = diff_s.dropna()
count_s = diff_s.sum()

comp_df = neatifyDiffExpComps(diffExpComps)
comp_df.to_csv("differentialExpressionComparison_ARATH_minSamp15.csv") #Can read in this file instead!

############3

all_s = renameTriqler2Vital_nonDetailed(df_psss3)
val1 = all_s[spectronautComps[0][0]]
val2 = all_s[spectronautComps[0][1]]
pVals = stats.ttest_ind(val1.T, val2.T).pvalue
pVals = pd.DataFrame(pVals)
qVals = qvalues(pVals)
#diff_s = (stats.ttest_ind(val1.T, val2.T).pvalue > 0.05) 
diff_s = (qVals < 0.05)
count_s = diff_s.sum()
#ax.set_xticks(45)
# ToDo
# - Split up the analysis between species.

sns.set(style="whitegrid")
tips = sns.load_dataset("tips")
ax = sns.barplot(x="day", y="total_bill", data=tips)
####################


# Drop S01 or S10
at_s = at_s.drop(["S01", "S10"], axis = 1)
ce_s = ce_s.drop(["S01", "S10"], axis = 1)
hs_s = hs_s.drop(["S01", "S10"], axis = 1)


def imputeFullSamples(df, impute_method = "min"):
    """
    Drops all non-full samples and impute the rest.
    df - dataframe with proteins as rows and multiple samples as columns.
    
    example:
        - imputeFullSamples(at_s)
    """
    df = df
    samples = ["S0"+str(i) for i in range(1,10)]
    samples.append("S10")
    for i in samples:
        df = df[df[i].sum(axis = 1) != 0]
    if impute_method == "mean":
        df = df.T.fillna(df.mean(axis=1)).T
    elif impute_method == "median":
        df = df.T.fillna(df.median(axis=1)).T
    elif impute_method == "min":
        df = df.T.fillna(df.min(axis=1)).T
    elif impute_method == "max":
        df = df.T.fillna(df.max(axis=1)).T
    else:
        raise Exception("Invalid impute_method!")
    return df

impute_method = "min"
at_s = imputeFullSamples(at_s, impute_method = impute_method)
ce_s = imputeFullSamples(ce_s, impute_method = impute_method)
hs_s = imputeFullSamples(hs_s, impute_method = impute_method)


at_s


#### TEST
#ce_s[ce_s.index == "V6CM07"]["S01"].sum(axis=1)
#ce_s[ce_s.index == "V6CLV0"]["S01"].sum(axis=1)
#ce_s[ce_s.index == "V6CLU8"]["S01"].sum(axis=1)

#ce_s[ce_s.index == "V6CM07"]["S10"].sum(axis=1)
#ce_s[ce_s.index == "V6CLV0"]["S10"].sum(axis=1)
#ce_s[ce_s.index == "V6CLU8"]["S10"].sum(axis=1)

###################################################
## Count detected proteins per species and method #
###################################################
protein_count = pd.DataFrame(columns = ['species', 'method'])
protein_count['species'] = (['A thaliana']*len(at_t) + ['C elegans']*len(ce_t) + ["H sapiens"]*len(hs_t) + 
             ['A thaliana']*len(at_s) + ['C elegans']*len(ce_s) + ["H sapiens"]*len(hs_s))
protein_count['method'] = (['Triqler']*len(at_t) + ['Triqler']*len(ce_t) + ['Triqler']*len(hs_t) + 
             ['Spectronaut']*len(at_s) + ['Spectronaut']*len(ce_s) + ['Spectronaut']*len(hs_s))

fig, ax = plt.subplots(figsize=(15, 9))
sns.countplot(x='species', hue='method', data=protein_count)
sns.despine(offset=10, trim=True);
plt.show()

##################################################
#                                               ##
#   ########################################    ##
#   ###     COMPARISONS                 ####    ##
#   ########################################    ##
#                                               ##    
##################################################

######################################
# Normalize to per-protein abundance #
######################################

### TEST
ce_s[ce_s.index == "A0A0K3AQU9"]#["S05"]
prettyPrint(ce_s[ce_s.index == "A0A0K3AQU9"])#["S05"]

ce_s[ce_s.index == "A0A0K3ARN9"]
prettyPrint(ce_s[ce_s.index == "A0A0K3ARN9"])

ce_s[ce_s.index == "V6CLV0"].sum(axis=1)

####

repl_cols = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]

# Dopped S01 and S10
repl_cols = ["S02","S03","S04","S05","S06","S07","S08","S09"]

ce_s.T/ce_s.T.sum()

ce_s.T.sum()
ce_s.sum(axis = 1)

normalized_protein_at_s = at_s.div(at_s.sum(axis=1), axis=0)
normalized_protein_ce_s = ce_s.div(ce_s.sum(axis=1), axis=0)
normalized_protein_hs_s = hs_s.div(hs_s.sum(axis=1), axis=0)

#normalized_protein_at_s = at_s.div(at_s.prod(axis=1), axis=0)
#normalized_protein_ce_s = ce_s.div(ce_s.prod(axis=1), axis=0)
#normalized_protein_hs_s = hs_s.div(hs_s.prod(axis=1), axis=0)

##########################################
## Merge normalized spectronaut species ##
##########################################

sampgroup_at_s,sampgroup_hs_s,sampgroup_ce_s = pd.DataFrame(),pd.DataFrame(),pd.DataFrame()

def adjustPandasNaNSum(df, replicate = "S01", maxNaN = 0):
    """
    Pandas sum of NaN = 0, this
    
    maxNaN - how many NaN per sample to allow.
    """
    repl = replicate
    replicates = df[repl]
    replicates = replicates[replicates.isna().sum(axis = 1) <= maxNaN ]
    replicates['total'] = replicates.sum(axis=1)    
    replicates.loc[replicates[repl].isnull().all(1),'total']=np.nan
    return replicates['total']

######################################################
## FIX THIS SSSSUUUUUM STUFFFY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##################################################################

maxNaN = 6
for s in repl_cols: #<----------------------------------------- Sums of stuff, we need to drop before summing up stuff!
    #sampgroup_at_s[s] = normalized_protein_at_s[s].sum(1)
    #sampgroup_ce_s[s] = normalized_protein_ce_s[s].sum(1)
    #sampgroup_hs_s[s] = normalized_protein_hs_s[s].sum(1)
    sampgroup_at_s[s] = adjustPandasNaNSum(normalized_protein_at_s, s, 5)
    sampgroup_ce_s[s] = adjustPandasNaNSum(normalized_protein_ce_s, s, 5)
    sampgroup_hs_s[s] = adjustPandasNaNSum(normalized_protein_hs_s, s, 5)

#for s in repl_cols:
#sampgroup_at_s[s] = normalized_protein_at_s[s].stack().reset_index().drop("level_1", axis = 1).set_index("protein")[0]
#s="S01"
#s_1 = normalized_protein_hs_s["S01"].stack().reset_index().drop("level_1", axis = 1).set_index("protein").rename(index = str, columns={0:"S01"})
#s_2 = normalized_protein_hs_s["S02"].stack().reset_index().drop("level_1", axis = 1).set_index("protein").rename(index = str, columns={0:"S02"})
#s_3 = normalized_protein_hs_s["S03"].stack().reset_index().drop("level_1", axis = 1).set_index("protein").rename(index = str, columns={0:"S03"})

#pd.concat([s_1, s_)


sampgroup_at_s['Species'] = ['A thaliana'] *len(sampgroup_at_s)
sampgroup_ce_s['Species'] = ['C elegans']  *len(sampgroup_ce_s)
sampgroup_hs_s['Species'] = ['H sapiens']  *len(sampgroup_hs_s)
sampgroup_at_s['Method']  = ['Spectronaut']*len(sampgroup_at_s)
sampgroup_ce_s['Method']  = ['Spectronaut']*len(sampgroup_ce_s)
sampgroup_hs_s['Method']  = ['Spectronaut']*len(sampgroup_hs_s)

sampgroup_spec_s = pd.concat([sampgroup_at_s,sampgroup_ce_s,sampgroup_hs_s])

######################################
## Merge normalized Triqler species ##
######################################

normalized_protein_at_t = at_t.div(at_t.sum(axis=1), axis=0)
normalized_protein_ce_t = ce_t.div(ce_t.sum(axis=1), axis=0)
normalized_protein_hs_t = hs_t.div(hs_t.sum(axis=1), axis=0)

#normalized_protein_at_s = at_s.div(at_s.prod(axis=1), axis=0)
#normalized_protein_ce_s = ce_s.div(ce_s.prod(axis=1), axis=0)
#normalized_protein_hs_s = hs_s.div(hs_s.prod(axis=1), axis=0)
sampgroup_at_t,sampgroup_hs_t,sampgroup_ce_t = pd.DataFrame(),pd.DataFrame(),pd.DataFrame()

for t in repl_cols:
    sampgroup_at_t[t] = normalized_protein_at_t[t].sum(1)
    sampgroup_ce_t[t] = normalized_protein_ce_t[t].sum(1)
    sampgroup_hs_t[t] = normalized_protein_hs_t[t].sum(1)
    
sampgroup_at_t['Species'] = ['A thaliana'] *len(normalized_protein_at_t)
sampgroup_ce_t['Species'] = ['C elegans']  *len(normalized_protein_ce_t)
sampgroup_hs_t['Species'] = ['H sapiens']  *len(normalized_protein_hs_t)
sampgroup_at_t['Method']  = ['Triqler']*len(normalized_protein_at_t)
sampgroup_ce_t['Method']  = ['Triqler']*len(normalized_protein_ce_t)
sampgroup_hs_t['Method']  = ['Triqler']*len(normalized_protein_hs_t)

sampgroup_spec_t = pd.concat([sampgroup_at_t,sampgroup_ce_t,sampgroup_hs_t])


###############################
# Merge both analysis methods #
###############################
sampgroup = sampgroup_spec_t
sampgroup = sampgroup_spec_s
sampgroup = sampgroup_spec_t.append(sampgroup_spec_s)


#sampgroup.where((pd.notnull(sampgroup)), np.NaN,inplace=True)
#sampgroup.dropna(inplace=True)


# Check count for each sample after dropna.
test = sampgroup[sampgroup["Method"] == "Spectronaut"]
(test.dropna()["Species"] == "H sapiens").sum()
(test.dropna()["Species"] == "A thaliana").sum()
(test.dropna()["Species"] == "C elegans").sum()

##################
# LATE IMPUTE ####
##################

# IMPLEMENT LATE IMPUTE HERE

#nn = sampgroup[repl_cols_s][sampgroup[repl_cols_s] > 0.0].min(numeric_only=True,axis=1)/2.0
#nnn = pd.concat([nn]*sampgroup.shape[1],axis=1)
#sampgroup = sampgroup.where(sampgroup!=0.0, nnn)
        

#sampgroup['Protein'] = range(len(sampgroup))
sampgroup["Protein"] = sampgroup.index

unmeltsamp = sampgroup.copy()
sampgroup=pd.melt(sampgroup, id_vars=['Method', 'Species','Protein'], value_vars=repl_cols)
sampgroup.rename(columns={'variable':'Sample','value':'Expression'}, inplace=True)

###########
# FILLNA ##
###########

sampgroup = unmeltsamp
sampgroup = pd.melt(sampgroup, id_vars=['Method', 'Species','Protein'], value_vars=repl_cols)
sampgroup.rename(columns={'variable':'Sample','value':'Expression'}, inplace=True)

sampgroup = sampgroup.fillna(sampgroup["Expression"].mean())
sampgroup = sampgroup.fillna(sampgroup["Expression"].median())
sampgroup = sampgroup.fillna(sampgroup["Expression"].min())

a_mix = np.array([0.5] * 10)
a_mix[0] += 0.0001
c_mix = np.array([0.5, 0.25, 0.125, 0.0625, 0.031, 0.0155, 0.008, 0.004, 0.002, 0.0])
h_mix = np.array([0.0, 0.25, 0.375, 0.4375, 0.469, 0.4845, 0.492, 0.496, 0.498, 0.5])
a_mix = a_mix/a_mix.sum()
c_mix = c_mix/c_mix.sum()
h_mix = h_mix/h_mix.sum()
mix = {'A thaliana':a_mix, 'C elegans': c_mix, 'H sapiens': h_mix}

# Dropped S01 and S10
a_mix = np.array([0.5] * 8)
a_mix[0] += 0.0001
c_mix = np.array([0.25, 0.125, 0.0625, 0.031, 0.0155, 0.008, 0.004, 0.002])
h_mix = np.array([0.25, 0.375, 0.4375, 0.469, 0.4845, 0.492, 0.496, 0.498])
a_mix = a_mix/a_mix.sum()
c_mix = c_mix/c_mix.sum()
h_mix = h_mix/h_mix.sum()
mix = {'A thaliana':a_mix, 'C elegans': c_mix, 'H sapiens': h_mix}



fig, ax = plt.subplots(figsize=(13, 8))
#sns.violinplot(data=sampgroup.loc[sampgroup['Species'] == "H sapiens"], x="Sample",y="Expression", hue="Method",split=True)
sns.boxplot(data=sampgroup.loc[sampgroup['Species'] == "H sapiens"], x="Sample",y="Expression", hue="Method",fliersize=0)
for i in range(len(h_mix)):
    ax.plot([i-.2, i+.2], [h_mix[i], h_mix[i]], '-', color='red', lw=5)
ax.set_ylim(0,0.26)
ax.set_ylabel("Fraction of protein total abundance per sample")
sns.despine(offset=10, trim=True)
plt.title("Homo sapiens")  
fig.savefig("sapiens_sample.png")
fig.savefig("TETETETE.png")

fig, ax = plt.subplots(figsize=(13, 8))
#sns.violinplot(data=sampgroup.loc[sampgroup['Species'] == "H sapiens"], x="Sample",y="Expression", hue="Method",split=True)
sns.boxplot(data=sampgroup.loc[sampgroup['Species'] == "C elegans"], x="Sample",y="Expression", hue="Method",fliersize=0)
for i in range(len(c_mix)):
    ax.plot([i-.2, i+.2], [c_mix[i], c_mix[i]], '-', color='red', lw=5)
ax.set_ylim(0,0.9)
ax.set_ylabel("Fraction of protein total abundance per sample")
plt.title('Caenorhabditis elegans')
sns.despine(offset=10, trim=True)

fig, ax = plt.subplots(figsize=(15, 9)) 
#sns.violinplot(data=sampgroup.loc[sampgroup['Species'] == "H sapiens"], x="Sample",y="Expression", hue="Method",split=True)
sns.boxplot(data=sampgroup.loc[sampgroup['Species'] == "A thaliana"], x="Sample",y="Expression", hue="Method",fliersize=0)
for i in range(len(a_mix)):
    ax.plot([i-.2, i+.2], [a_mix[i], a_mix[i]], '-', color='red', lw=5)
ax.set_ylim(0,0.2)
ax.set_ylabel("Fraction of protein total abundance per sample")
plt.title('Arabidopsis thaliana')
sns.despine(offset=10, trim=True)


#############################
## CORRELATION COEFFICIENTS #
#############################
# CORRELATION PLOTTING FUNCTIONS.
from scipy.stats.stats import pearsonr
from numpy import log10

def getMixtures(trunc = False):
    if trunc == False:
        a_mix = np.array([0.5] * 10)
        a_mix[0] += 0.0001
        c_mix = np.array([0.5, 0.25, 0.125, 0.0625, 0.031, 0.0155, 0.008, 0.004, 0.002, 0.0])
        h_mix = np.array([0.0, 0.25, 0.375, 0.4375, 0.469, 0.4845, 0.492, 0.496, 0.498, 0.5])
        a_mix = a_mix/a_mix.sum()
        c_mix = c_mix/c_mix.sum()
        h_mix = h_mix/h_mix.sum()
        mix = {'A thaliana':a_mix, 'C elegans': c_mix, 'H sapiens': h_mix}
        return mix, h_mix, c_mix, a_mix
    else:
        # Dropped S01 and S10
        a_mix = np.array([0.5] * 8)
        a_mix[0] += 0.0001
        c_mix = np.array([0.25, 0.125, 0.0625, 0.031, 0.0155, 0.008, 0.004, 0.002])
        h_mix = np.array([0.25, 0.375, 0.4375, 0.469, 0.4845, 0.492, 0.496, 0.498])
        a_mix = a_mix/a_mix.sum()
        c_mix = c_mix/c_mix.sum()
        h_mix = h_mix/h_mix.sum()
        mix = {'A thaliana':a_mix, 'C elegans': c_mix, 'H sapiens': h_mix}
        return mix, h_mix, c_mix, a_mix

def my_corr(row):
    mix, h_mix, c_mix, a_mix = getMixtures(trunc = True)
    
    mymix = mix[row['Species']]
    row = row.loc[repl_cols]
    ixs = (row != 0.0)
    a = np.array(row.loc[ixs])
    b = mymix[ixs]
    return pearsonr(a,b)[0]

def p_r_plot(species,figname=None):
    fig, ax = plt.subplots(figsize=(10, 6))
    N=5000
    for meth in unmeltsamp['Method'].unique():
        sns.distplot(unmeltsamp[(unmeltsamp['Method']==meth) & (unmeltsamp['Species']==species)]['Correlation, R'],
             hist_kws={"cumulative":-1, "linewidth": 3, "alpha": 1,"histtype": "step"},
             bins=N, kde=False, ax=ax, label=meth, norm_hist = True)
    ax.set_xlim(-0.9,1.0)
    ax.set_ylim(0,1.2)
    ax.set_title(species)
    sns.despine(offset=10, trim=True)
    ax.set_ylabel("Number of proteins")
    ax.set_xlabel("Cumulative samples "r"$\leq$" " correlation coefficient")
    plt.legend()
    if figname:
        fig.savefig(figname)


unmeltsamp['Correlation, R'] = unmeltsamp.apply(my_corr, axis=1)
unmeltsamp.dropna(inplace=True)
#unmeltsamp.isna().sum().sum()
unmeltsamp.isnull().sum().sum()



p_r_plot('H sapiens')
p_r_plot('C elegans')
p_r_plot('A thaliana')

spec = 'H sapiens'
unmeltsamp[unmeltsamp["Species"] == spec ]["Correlation, R"].min()
unmeltsamp[unmeltsamp["Species"] == spec ]["Correlation, R"].max()
