#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 16:18:48 2019

@author: ptruong
"""
from triqlerParser import *
from triqlerProcessor import * 
from utils import * 
from parseReformattedPSSS3equDecoy import * # Needs to tidy and place functions in right place!
from corrPlot import * #correlation plotting
from qvalue import *

import os
import pandas as pd
import numpy as np
from scipy import stats
import scipy as sp
import subprocess
#%matplotlib inline
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set(font_scale=2)
sns.set_style("whitegrid", {'axes.grid' : False})#, 'grid.color': '.1', 'grid.linestyle': u'-'})

##############
# Diff Exp ###
##############

def getComparisons(fileDir = "triqlerOutput_noShared_largeScaleOptimizations", filePrefix = "proteins", fileSuffix = ".tsv"):
    """
    Get all comparisons from triqler outputs.
    """
    fileDir = fileDir 
    startFrom = len(filePrefix) + 1
    endAt = len(fileSuffix)
    dirFiles = os.listdir(fileDir)
    comparisons = []
    for i in dirFiles:
        comparisons.append(i[startFrom:-endAt].split("vs"))
    return comparisons

def getComparisonLabels(comparisons):
    """
    Take comparisons list from getComparisons to get comparison labels.
    """
    comps = []
    for i in comparisons:
        comps.append(i[0]+"vs"+i[1])
    return comps

def getTriqlerComparisons(comparisons, prefix = "proteins", suffix = "tsv"):
    """
    Construct triqler comparisons names.
    """
    triqlerComps = []
    for i in comparisons:
        triqlerComps.append(prefix+"."+i[0]+"vs"+i[1]+"."+suffix)
    return triqlerComps

def getTriqlerComparisonsFromFileDir(fileDir = "triqlerOutput_noShared_largeScaleOptimizations"):
    """
    Get triqler comparisons from file directory.
    """
    dirFiles = os.listdir(fileDir)
    return dirFiles

def getSpectronautComparison(comparisons):
    """
    Get Spectronaut comparisons from comparisons list.
    """
    spectronautComps = []
    for i in comparisons:
        if i[0] == "10":
           bufferStr1 = ""
        else:
            bufferStr1 = "0"
        if i[1] == "10":
            bufferStr2 = ""
        else:
            bufferStr2 = "0"
        i1 = "S"+bufferStr1+i[0]
        i2 = "S"+bufferStr2+i[1]
#        i1 = "S0"+i[0]
#        i2 = "S0"+i[1]
        spectronautComps.append([i1, i2])
    return spectronautComps

def processComparisons(triqlerOutputDir, PSSS3, comparisons, triqlerComps, spectronautComps, FDRTreshold, specie, logSpectronaut):  
    """
    Get differential expressions from  triqler and PSSS3.
    
    triqlerOutputDir - directory of triqler outputs.
    PSSS3 = addSpecieProteinColumns(psss3, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
    """
    diffExpDir = triqlerOutputDir
    if logSpectronaut == True:
        PSSS3.iloc[:,3:] = np.log10(PSSS3.iloc[:,3:])
    all_s = PSSS3
    diffExpComps = []
    counter = 0
    print("loading files from " + diffExpDir)
    for i in range(len(comparisons)):
        print(str(counter) + "of" + str(len(comparisons)))
        #Triqler Part
        triqler = table2Df(diffExpDir+"/" + triqlerComps[i], peptideSeperator = True)
        triqler = triqler2Numeric(triqler)
        triqler = getTriqlerDetailed(triqler)
        if specie == "ARATH":
            triqler = triqler[triqler.specie == "ARATH"]
        elif specie == "CAEEL":
            triqler = triqler[triqler.specie == "CAEEL"]
        elif specie == "HUMAN":
            triqler = triqler[triqler.specie == "HUMAN"]
        #diff_exp_treshold = str(diffExpTreshold)
        FDR_treshold = FDRTreshold
        #diff_t = triqler["diff_exp_prob_"+diff_exp_treshold]
        diff_t = triqler["q_value"]
        diff_t = pd.to_numeric(diff_t)
#        diff_t = (diff_t > (1-FDR_treshold))  # diff_t over 1-treshold. We want 95% or more prob of being diff exp.
        #diff_t = (diff_t < 0.05)
        diff_t = (diff_t < FDR_treshold)
        count_t  = diff_t.sum()
    
        #Spectronaut Part
        #all_s = renameTriqler2Vital_nonDetailed(df_psss3)
        if False:
            val1 = all_s[spectronautComps[i][0]]
            val2 = all_s[spectronautComps[i][1]]
            
            diff_s = diffExp(val1, val2)
            diff_s = pd.DataFrame(diff_s)
            diff_s = diff_s.dropna()
            count_s = diff_s.sum()
            
        if specie == "ARATH":
            PSSS3 = PSSS3[PSSS3.specie == "ARATH"]
            all_s = renameTriqler2Vital_nonDetailed(PSSS3)
        elif specie == "CAEEL":
            PSSS3 = PSSS3[PSSS3.specie == "CAEEL"]
            all_s = renameTriqler2Vital_nonDetailed(PSSS3)
        elif specie == "HUMAN":
            PSSS3 = PSSS3[PSSS3.specie == "HUMAN"]
            all_s = renameTriqler2Vital_nonDetailed(PSSS3)
        else:
            all_s = renameTriqler2Vital_nonDetailed(PSSS3)
            
        val1 = all_s[spectronautComps[i][0]]
        val2 = all_s[spectronautComps[i][1]]
        pVals = stats.ttest_ind(val1.T, val2.T).pvalue
        pVals = pd.DataFrame(pVals)
        qVals = qvalues(pVals)
        #diff_s = (stats.ttest_ind(val1.T, val2.T).pvalue > 0.05) 
        diff_s = (qVals < FDR_treshold) # treshold by Q-value
        count_s = diff_s.sum()
        
        diffExpComp = [count_t, count_s]
        diffExpComps.append(diffExpComp)
        counter+=1
    return diffExpComps

def neatifyDiffExpComps(diffExpComps, comparisons):
    """
    neatify diffEcpComps output.
    """
    comps = getComparisonLabels(comparisons)
    valueList1 = []
    valueList2 = []
    sampleList1 = []
    sampleList2 = []
    for i in range(len(diffExpComps)):
        valueList1.append(diffExpComps[i][0])
        valueList2.append(diffExpComps[i][1][0])
        sampleList1.append(int(comps[i].split("vs")[0]))
        sampleList2.append(int(comps[i].split("vs")[1]))
    comp_df = pd.DataFrame([valueList1, valueList2, sampleList1, sampleList2], index = ["triqler", "spectronaut", "sample1", "sample2"])
    comp_df = comp_df.T
    #comp_df = comp_df.sort_index()
    comp_df = comp_df.sort_values(["sample1", "sample2"])
    #comp_df.to_csv("differentialExpressionComparison.csv") #Can read in this file instead!
    return comp_df


########################
# plotting comparisons #
########################

def meltCompDf(comp_df):
    """
    melt the comp_df to seaborn plotting format.
    
    NOTE: should lookup how pd.melt function works!
    """
    allVals = []
    for i in range(len(comp_df)):
        value_t = [str(comp_df.sample1[i])+"vs"+str(comp_df.sample2[i]), comp_df.triqler[i], "triqler", comp_df.sample1[i], comp_df.sample2[i]]
        value_s = [str(comp_df.sample1[i])+"vs"+str(comp_df.sample2[i]), comp_df.spectronaut[i], "spectronaut", comp_df.sample1[i], comp_df.sample2[i]]
        allVals.append(value_t)
        allVals.append(value_s)
    
    all_df = pd.DataFrame(allVals, columns = ["comparison","diffExp", "method", "sample1", "sample2"])
    all_df = all_df.sort_values(["sample1", "sample2"])
    #d.melt(comp_df, id_vars=['triqler', 'spectronaut'])
    return all_df


def barplotComps(meltedDf, title, figname = None):    
    """
    Bar diff expression comparison plot between methods.
    """
    import matplotlib.pyplot as plt
    a4_dims = (11.7*1.5, 8.27*1.5)
    #df = mylib.load_data()
    fig, ax = plt.subplots(figsize=a4_dims)
    plt.xticks(rotation=90)
    sns.barplot(x="comparison", y="diffExp", hue="method", data=meltedDf, ax = ax).set_title(title)
    if figname:
        fig.savefig(figname)

"""
#def getDifferentialComps:
comparisons = getComparisons(fileDir = "triqlerOutput_noShared_largeScaleOptimizations", filePrefix = "proteins", fileSuffix = ".tsv")
triqlerComps  = getTriqlerComparisons(comparisons, prefix = "proteins", suffix = "tsv")
spectronautComps = getSpectronautComparison(comparisons)

diffExpDir = "triqlerOutput_noShared_largeScaleOptimizations_minSamp15"
diffExpTreshold = 1.0 # Integer to be converted to String to match triqler output labelling.
FDRTreshold = 0.05
specie = None
diffExpComps = processComparisons(triqlerOutputDir = diffExpDir,
                   PSSS3 = df_psss3,
                   comparisons = comparisons,
                   triqlerComps = triqlerComps,
                   spectronautComps = spectronautComps,
                   FDRTreshold = FDRTreshold,
                   specie = "ARATH")
all_df = meltCompDf(comp_df)
barplotComps(all_df, "ARATH")
"""