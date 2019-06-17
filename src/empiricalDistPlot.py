#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 10:44:45 2019

@author: ptruong
"""

import numpy as np
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt

import os
import shutil


"""
df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup_decoy.csv", sep = "\t")

df = pd.read_csv("500-PSSS3-equ decoy_Report_filteredOn_PGQvalue_EGQvalue_PGCscore_EGCscore.xls", sep = "\t")
run_ids = df.run_id.unique()
run_ids.sort()

os.mkdir("distributions_allDist")
#f = open("distributions_log/percNaN.csv", "w")

#i = 44
fig, ax = plt.subplots(figsize=(23.4, 16.54))
first_iter = True
kde = None
count = 0
for i in range(len(run_ids)):
    print("Generating plots for " + run_ids[i])
    #df[df["run_id"] == run_ids[i]]
    #data = df.proteinQuantity
    data = df[df["run_id"] == run_ids[i]].proteinQuantity
    data = np.log(data)
    percNaN = data.isna().sum()/(data.isna().sum() + data.count())
    print("Percent NaN: " + str(percNaN))
    #f.write(run_ids[i] + "\t" + str(percNaN))
    x = data.dropna()
    
    #imputeConst = data.mean()
    #x = data.fillna(imputeConst)
    sns.distplot(x, bins = 1000, hist = False, kde_kws = {"color":"red", "alpha":0.2}, ax = ax).set_title(run_ids[i] + " PercNaN: " + str(percNaN))    
    #sns.distplot(x, bins = 1000, hist = False, kde_kws = {"color":"red", "alpha":0.2}, ax = ax).set_title("All sample distributions")
#    if first_iter == True:
#        kde = sns.distplot(x).get_lines()[0].get_data()[1]
#        count += 1
#    else:
#        kde += sns.distplot(x).get_lines()[0].get_data()[1]
#        count += 1
#kde_mu = kde / count * 5
    fig.savefig("distributions_allDist/"+run_ids[i]+".png")
    #fig.clf()
#fig.savefig("distributionsTEST/"+"allDist_log.png")
    #density = gaussian_kde(data)
    
    
# ToDO make average and confidence interval for the distributions...
    
dat = sns.distplot(x).get_lines()[0].get_data()[1]
dat2 = sns.distplot(x).get_lines()[0].get_data()[1] 
len(df)
"""

################################
###### ON PEPTIDE ##############
################################

def calcPercNaN(data):
    percNaN = data.isna().sum()/(data.isna().sum() + data.count())
    return percNaN

def getConditions(df):
    conditions = df["R.Condition"].str[-3:].unique()
    return conditions

def getSampleGroupDistributions(df, conditions, variable = "FG.NormalizedMS2PeakArea", outputDest = "distributions", figureSize = (23.4, 16.54)):
    """
    Function generates empirical distributions and kde for all sample groups in the PSSS3 dataset.
    
    variables = {"FG.NormalizedMS2PeakArea", "EG.Cscore", "EG.PEP", etc...}
    """
    outputDest = outputDest 
    
    if os.path.isdir(outputDest) != True:
        os.mkdir(outputDest)
    
    fig, ax = plt.subplots(figsize=figureSize)
    alpha = 1
    for i in range(len(conditions)):
        fig, ax = plt.subplots(figsize=figureSize)
        print("Generating plots for " + conditions[i])
        sampleGroup = df[df["R.Condition"].str[-3:] == conditions[i]]     
        data = sampleGroup[variable]
        data = data.dropna()
        data = np.log(data)
        percNaN = calcPercNaN(data)
        print("Percent NaN: " + str(percNaN))
        sns.distplot(data, bins = 1000, hist = True, 
                     kde_kws = {"color":"red", "alpha":alpha}, 
                     ax = ax).set_title(conditions[i] + " PercNaN: " + str(percNaN))
        fig.savefig(outputDest  + "/" + conditions[i] + ".png")
        #fig.clf()    

def getSampleGroupDistributionsAllPlot(df, conditions, variable = "FG.NormalizedMS2PeakArea", alpha = 0.2, outputDest = "distributions", hist = False, figureSize = (23.4, 16.54)):
    """
    Function for generating all kde distributions in one plot. Lower alpha gives transparence where kde are less dense.
    """
    outputDest = outputDest 
    if os.path.isdir(outputDest) != True:
        os.mkdir(outputDest)
    
    fig, ax = plt.subplots(figsize=figureSize)
    alpha = alpha
    for i in range(len(conditions)):
        print("Generating data for " + conditions[i])
        sampleGroup = df[df["R.Condition"].str[-3:] == conditions[i]]     
        data = sampleGroup[variable]
        data = data.dropna()
        data = np.log(data)
        sns.distplot(data, bins = 1000, hist = hist, 
                     kde_kws = {"color":"red", "alpha":alpha},
                     hist_kws = {"color":"blue", "alpha":alpha},
                     ax = ax).set_title(conditions[i] + " PercNaN: " + str(percNaN))
    if hist == False:
        figname = outputDest  + "/" + "All samples" + ".png"
    else:
        figname = outputDest  + "/" + "All samples" + "Hist.png"
    fig.savefig(figname)
            

if __name__ == "__main__":
    df = pd.read_csv("500-PSSS3-equ decoy_Report_filteredOn_PGQvalue_EGQvalue_PGCscore_EGCscore.xls", sep = "\t")
    conditions = getConditions(df)
    var = "FG.NormalizedMS2PeakArea"
    getSampleGroupDistributions(df, conditions, variable = var, 
                                outputDest = "distributions")
    getSampleGroupDistributionsAllPlot(df, conditions, variable = var, alpha = 0.2, 
                                       outputDest = "distributions")
    getSampleGroupDistributionsAllPlot(df, conditions, variable = var, alpha = 0.2, 
                                       outputDest = "distributions", hist = True)