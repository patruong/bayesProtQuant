#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 15:44:43 2019

@author: ptruong
"""

import pandas as pd
import seaborn as sns
import numpy as np

def logPeptideIntensityDistPlot(specie, title):
    intensities = df[df["species"] == specie].intensity
    ax = sns.distplot(np.log10(intensities.values), norm_hist=False)
    ax.set(xlabel = "log10(peptide intensity)", ylabel = "Probability")
    ax.set_title(specie + " " + "density")
    plt.savefig(specie + "distPlot.png")
    plt.clf()    

def logPeptideIntensityDistPlotSpecSamp(df, specie = None, sample = None):
    """
    plot distplot for species and sample
    """
    specie = specie 
    sample = sample
    
    if specie == None:
        df_spec = df
    else:
        df_spec = df[df["PG.Organisms"] == specie]
    if sample == None:
        df_spec_samp = df_spec
    else:
        df_spec_samp = df_spec[df_spec["R.Condition"].str[-3:] == sample]
    print("Plotting " + specie + ":" + sample)
    intensities = df_spec_samp["FG.NormalizedMS2PeakArea"]
    ax = sns.distplot(np.log10(intensities.values), norm_hist=False)
    ax.set(xlabel = "log10(peptide intensity)", ylabel = "Probability")
    ax.set_title(specie + " " + sample + " " + "density")
    ax.set_ylim((0, 1))
    ax.set_xlim((0, 10))
    plt.savefig(specie + "_" + sample + "_" + "distPlot.png")
    plt.clf()
    
def logPeptideIntensityDistPlotSpec(df, specie = None):
    """
    plot distplot for species.
    """
    specie = specie 
    if specie == None:
        df_spec = df
    else:
        df_spec = df[df["PG.Organisms"] == specie]
    intensities = df_spec["FG.NormalizedMS2PeakArea"]
    if specie == None:
        print("Plotting all samples")
    else:
        print("Plotting " + specie)
    ax = sns.distplot(np.log10(intensities.values), norm_hist=False)
    ax.set(xlabel = "log10(peptide intensity)", ylabel = "Probability")
    ax.set_ylim((0, 1))
    ax.set_xlim((0, 10))
    if specie == None:
        ax.set_title("all density")
        plt.savefig("all_distPlot.png")
    else:
        ax.set_title(specie + " " + "density")
        plt.savefig(specie + "distPlot.png")
    plt.clf()
def logPeptideIntensityDistPlotSamp(df, sample = None):
    """
    plot distplot for sample.
    """

    sample = sample
    if sample == None:
        df_samp = df
    else:
        df_samp = df[df["R.Condition"].str[-3:] == sample]    
    intensities = df_samp["FG.NormalizedMS2PeakArea"]
    if sample == None:
        print("All samples")
    else:
        print("Plotting " + sample)
    ax = sns.distplot(np.log10(intensities.values), norm_hist=False)
    ax.set(xlabel = "log10(peptide intensity)", ylabel = "Probability")
    ax.set_ylim((0, 1))
    ax.set_xlim((0, 10))
    if sample == None:
        ax.set_title("All density")
        plt.savefig("All_distplot.png")
    else:
        ax.set_title(sample + " " + "density")
        plt.savefig(sample + "distPlot.png")
    plt.clf()
if __name__ == "__main__":
    ### NOTE SHOULD NOT BE LIKE THIS NEED TO DROP JUNK PEPTIDES...
    """
    df = pd.read_csv("PSSS3_triqlerFormatted_nonShared.csv", sep = "\t")
    df["species"] = df.proteins.str[-5:]
    df["decoy"] = df.proteins.str[:5]
    df = df[df["decoy"] != "decoy"]
    specList = df["species"].unique()
    
    df = df.drop(["run", "charge", "searchScore", "peptide"], axis = 1)
    
    df_spec = df[df["species"] == "CAEEL"]
    df_spec = df_spec[df_spec["condition"] == "S10"]    
    intensities = df_spec.intensity
    ax = sns.distplot(np.log10(intensities.values))
    
    
    logPeptideIntensityDistPlot(specList[0], "CElegans")
    logPeptideIntensityDistPlot(specList[1], "HSapiens")
    logPeptideIntensityDistPlot(specList[2], "AThaliana")
    """
    FDR_treshold = 0.01
    df = pd.read_csv("500-PSSS3-equ decoy_Report.xls", sep ="\t")
    df = df[~df["EG.Qvalue"].isna()]
    df_f = df[df["EG.Qvalue"] < FDR_treshold]
    
    species = df["PG.Organisms"].unique()
    repl_cols = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
    
    for specie in species:
        logPeptideIntensityDistPlotSpec(df_f, specie)
    for sample in repl_cols:
        logPeptideIntensityDistPlotSamp(df_f, sample)    
    for sample in repl_cols:
        for specie in species:
            logPeptideIntensityDistPlotSpecSamp(df_f, specie, sample)
    logPeptideIntensityDistPlotSamp(df_f) # All data  
    #logPeptideIntensityDistPlotSpec(df_f, species[0])    
    #logPeptideIntensityDistPlotSamp(df_f, "S01")    
    
