#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 10:41:19 2019

@author: ptruong
"""
import sys 

from proteinSearch import *

from triqlerParser import *
from triqlerProcessor import * 
from utils import * 
from parseReformattedPSSS3equDecoy import * # Needs to tidy and place functions in right place!
from proc import *
#from corrPlot import * #correlation plotting

import os
import argparse
import pandas as pd
import numpy as np
from scipy import stats
import scipy as sp
import subprocess
import matplotlib.pyplot as plt 
import seaborn as sns

def comparisonPlot(sampgroup, specie, hardRange = True, truncated = False, title = "plot"):
    
    homoSapiens = ["hs", "HS", "homo sapiens", "Homo Sapiens", "Homo sapeins", 
                   "homosapiens", "HomoSapiens", "HUMAN", "human",
                   "H sapeins", "h sapiens", "H Sapeins"]
    caenorhabditisElegans = ["ce", "CE", "Caenorhabditis elegans", "Caenorhabditis Elegans",
                             "caenorhabditis elegans", "caenorhabditiselegans", "CaenorhabditisElegans",
                             "CAEEL", "caeel",
                             "C elegans", "c elegans", "C Elegans"]
    arabidopsisThaliana = ["at", "AT", "Arabidopsis thaliana", "Arabidopsis Thaliana",
                           "arabidopsis thaliana", "arabidopsisthaliana", "ArabidopsisThaliana",
                           "ARATH", "arath",
                           "A thaliana", "a thaliana", "A Thaliana"]
    
    species = homoSapiens + caenorhabditisElegans + arabidopsisThaliana
    
    if specie not in species:
        raise Exception("Incorrect specie!")
    elif specie in homoSapiens:
        a_mix, c_mix, h_mix, mix = getMixtures(truncated = truncated)
        specie = "H sapiens"
        species_mixture = h_mix
    elif specie in caenorhabditisElegans:
        a_mix, c_mix, h_mix, mix = getMixtures(truncated = truncated)
        specie = "C elegans"
        species_mixture = c_mix
    elif specie in arabidopsisThaliana:
        a_mix, c_mix, h_mix, mix = getMixtures(truncated = truncated)
        specie = "A thaliana"
        species_mixture = a_mix
    

    fig, ax = plt.subplots(figsize=(13, 8))
    #sns.violinplot(data=sampgroup.loc[sampgroup['Species'] == "H sapiens"], x="Sample",y="Expression", hue="Method",split=True)
    sns.boxplot(data=sampgroup.loc[sampgroup['Species'] == specie], x="Sample",y="Expression", hue="Method",fliersize=0)
    for i in range(len(species_mixture)):
        ax.plot([i-.2, i+.2], [species_mixture[i], species_mixture[i]], '-', color='red', lw=5)
    #ax.set_ylim(0,0.26)
    
    if hardRange == False:
        nonOutlierMax = sampgroup.loc[sampgroup['Species'] == specie]["Expression"].quantile([0.99]).values[0]
        ax.set_ylim(0, nonOutlierMax)
    else:
        if specie == "H sapiens":
            ax.set_ylim(0,0.26)
        elif specie == "C elegans":
            ax.set_ylim(0,0.9)
        elif specie == "A thaliana":
            ax.set_ylim(0,0.2)
    ax.set_ylabel("Fraction of protein total abundance per sample")
    sns.despine(offset=10, trim=True)
    plt.title(specie + " " + title)  
    fig.savefig(specie+ "_" + title+".png")

#############################
## CORRELATION COEFFICIENTS #
#############################
# CORRELATION PLOTTING FUNCTIONS.
from scipy.stats.stats import pearsonr
from numpy import log10

def my_corr(row):
    repl_cols = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
    a_mix, c_mix, h_mix, mix = getMixtures()
    mymix = mix[row['Species']]
    row = row.loc[repl_cols]
    ixs = (row != 0.0)
    a = np.array(row.loc[ixs])
    b = mymix[ixs]
    return pearsonr(a,b)[0]

def p_r_plot(unmeltsamp, species,figname=None):
    unmeltsamp['Correlation, R'] = unmeltsamp.apply(my_corr, axis=1)
    fig, ax = plt.subplots(figsize=(10, 6))
    N=5000
    for meth in unmeltsamp['Method'].unique():
        sns.distplot(unmeltsamp[(unmeltsamp['Method']==meth) & (unmeltsamp['Species']==species)]['Correlation, R'],
             hist_kws={"cumulative":-1, "linewidth": 3, "alpha": 1,"histtype": "step"},
             bins=N, kde=False, ax=ax, label=meth) #, norm_hist = )
    ax.set_xlim(-0.8,1.0)
    #ax.set_ylim(0,1.2)
    ax.set_title(species)
    sns.despine(offset=10, trim=True)
    ax.set_ylabel("Number of proteins")
    ax.set_xlabel("Cumulative samples "r"$\leq$" " correlation coefficient")
    plt.legend()
    #if figname:
    fig.savefig(figname)

#######################################
# GENERATE DIFF EXPRESSION COMPARISON #
#######################################

def generateDiffExpPlot(fileDirectory, df_psss3, FDR_treshold, logSpectronaut, figName = "plot"):
        comparisons = getComparisons(fileDir = fileDirectory, filePrefix = "proteins", fileSuffix = ".tsv")
        comparisons.remove([''])
        triqlerComps  = getTriqlerComparisons(comparisons, prefix = "proteins", suffix = "tsv")
        spectronautComps = getSpectronautComparison(comparisons)
        
        # Input directory of triqler outputs    
        diffExpDir = fileDirectory
        FDRTreshold = FDR_treshold
        specieList = [None, "ARATH", "CAEEL", "HUMAN"]
        for k in specieList:
            specie = k
            if specie == None:
                print("Generating AllSpecies DiffExp...")
            else:
                print("Generating " + specie + " DiffExp...")
            diffExpComps = processComparisons(triqlerOutputDir = diffExpDir,
                               PSSS3 = df_psss3,
                               comparisons = comparisons,
                               triqlerComps = triqlerComps,
                               spectronautComps = spectronautComps,
                               FDRTreshold = FDRTreshold,
                               specie = specie,
                               logSpectronaut = logSpectronaut)
            comp_df = neatifyDiffExpComps(diffExpComps, comparisons)
            all_df = meltCompDf(comp_df)
            if specie == None:
                specieTitle = "AllSpec"
            else:
                specieTitle = specie
            diffExpTitle = specieTitle + figName
            
            barplotComps(all_df, diffExpTitle, 
                         fileDirectory + "/results/" + specieTitle+'_diffExpPlot'+figName+".png")

