#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 10:36:42 2019

@author: ptruong
"""

import sys 

from proteinSearch import *

from triqlerParser import *
from triqlerProcessor import * 
from utils import * 
from parseReformattedPSSS3equDecoy import * # Needs to tidy and place functions in right place!
#from corrPlot import * #correlation plotting

import os
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from scipy import stats
import scipy as sp
import subprocess
import matplotlib.pyplot as plt 
import seaborn as sns


#############
# Triqler ###
#############

def processTriqler(triqlerFile = "proteins.1vs2.tsv", FDR_treshold = 0.01):
    """
    Parse and process triqler to analysis format.
    """    
    triqler = table2Df(triqlerFile, peptideSeperator = True)
    triqler = triqler.dropna()
    triqler = triqler2Numeric(triqler)
    triqler = getTriqlerDetailed(triqler)
    
    
    # FDR on protein_id_posterior_error_prob
    triqler.sort_values(by = "protein_id_posterior_error_prob", inplace = True)
    triqler["FDR"] = triqler["protein_id_posterior_error_prob"].expanding().mean()
    triqler = triqler[triqler["FDR"] < FDR_treshold]
    
    return triqler

def splitTriqlerBySpecies(triqler, exponential = 2, truncated = False): #transform = np.exp2):
    """
    Split triqler into species dataframes.
    
    example transform = np.exp2 or transform = np.exp
    truncated - can be subset of repl_cols
    """    
    # All runs
    at_t = triqler[triqler["specie"] == "ARATH"]
    ce_t = triqler[triqler["specie"] == "CAEEL"]
    hs_t = triqler[triqler["specie"] == "HUMAN"]
        
    # DROP species and drop peptides <------------------------
    at_t = renameTriqler2Vital(at_t)
    ce_t = renameTriqler2Vital(ce_t)
    hs_t = renameTriqler2Vital(hs_t)
    
    # transform 
    #if transform != False:
    #    at_t = transform(at_t)
    #    ce_t = transform(ce_t)
    #    hs_t = transform(hs_t)
    if truncated != False:
        if truncated == True:
            print("Dropping [S01 S02] in Triqler split...")
            at_t = at_t.drop(["S01", "S10"], axis = 1)
            ce_t = ce_t.drop(["S01", "S10"], axis = 1)
            hs_t = hs_t.drop(["S01", "S10"], axis = 1)
        else:
            print("Dropping [" + " ".join(truncated) +"] in Triqler split...")
            at_t = at_t.drop(truncateDrop(truncated), axis = 1)
            ce_t = ce_t.drop(truncateDrop(truncated), axis = 1)
            hs_t = hs_t.drop(truncateDrop(truncated), axis = 1)            
    if exponential != False:
        print("Before exponential sample...")
        print(hs_t.iloc[2582,:].values)
        print("Taking the " + str(exponential) + " exponential of the Triqler values...")
        at_t = exponential**(at_t)
        ce_t = exponential**(ce_t)
        hs_t = exponential**(hs_t)
        print("After exponential sample...")
        print(hs_t.iloc[2582,:].values)
    return at_t, ce_t, hs_t

############################
# Spectronaut unfiltered ###
############################

def processSpectronaut(spectronautFile = "500-PSSS3-raw-reformatted_dropna_dropdup_decoy.csv",
                       FDR_treshold = 0.01, impute = "mean", global_impute = True):
    """
    Parse and process reformatted spectronaut to analysis format. The data is pre-reformatted 
    because the formatting procedure takes some time.
    """
    df = pd.read_csv(spectronautFile, sep ="\t")
    
    # FDR filtereing
    df = df[df["qvalue"].fillna(-1)<FDR_treshold]
    
    # Add decoy prefix to decoy proteins (to keep the same format in triqler and spectronaut the same)
    decoy = df["decoy"].apply(lambda x: "decoy_" if x == True else "")
    df["protein"] = decoy+df["protein"]
    
    # Reformat spectronaut to triqler format
    df_reform = getTriqlerOutputFormat_from_reformattedRawPSSS3(df, noDecoy = False)
    #df_reform = df_reform.drop_duplicates()
    if global_impute == True:
        if impute in ["mean", "average", "Mean", "Average"]:
            df_reform = df_reform.fillna(df_reform.mean().mean())
        elif impute in ["median", "Median"]:        
            df_reform = df_reform.fillna(df_reform.median().median())
        elif impute in ["min","Min", "minimal", "Minimal"]:
            df_reform = df_reform.fillna(df_reform.min().min())
        elif impute in ["max", "Max", "maximal", "Maximal"]:
            df_reform = df_reform.fillna(df_reform.max().max())
        elif impute == None:
            pass
        else:
            raise Exception("Must fill in imputation method!")
    
    else:
        if impute in ["mean", "average", "Mean", "Average"]:
            df_reform = df_reform.T.fillna(df_reform.mean(axis=1)).T
        elif impute in ["median", "Median"]:        
            df_reform = df_reform.T.fillna(df_reform.median(axis=1)).T
        elif impute in ["min","Min", "minimal", "Minimal"]:
            df_reform = df_reform.T.fillna(df_reform.min(axis=1)).T
        elif impute in ["max", "Max", "maximal", "Maximal"]:
            df_reform = df_reform.T.fillna(df_reform.max(axis=1)).T
        elif impute == None:
            pass
        else:
            raise Exception("Must fill in imputation method!")
    
    # ToDo: Implement columnswise impute?????
    
    #df_reform = df_reform.dropna()
    df_reform.reset_index(level=0, inplace = True)
    psss3 = addSpecieProteinColumns(df_reform, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
    #psss3 = psss3.dropna()
    return psss3


def splitSpectronautBySpecies(psss3, truncated = False):
    """
    Split Spectronaut into species dataframes.
    """     
    at_s = psss3[psss3["specie"] == "ARATH"]
    ce_s = psss3[psss3["specie"] == "CAEEL"]
    hs_s = psss3[psss3["specie"] == "HUMAN"]
    
    at_s = renameTriqler2Vital_nonDetailed(at_s)
    ce_s = renameTriqler2Vital_nonDetailed(ce_s)
    hs_s = renameTriqler2Vital_nonDetailed(hs_s)
    
    if truncated != False:
        if truncated == True:
            print("Dropping [S01 S02] in Spectronaut split...")
            at_s = at_s.drop(["S01", "S10"], axis = 1)
            ce_s = ce_s.drop(["S01", "S10"], axis = 1)
            hs_s = hs_s.drop(["S01", "S10"], axis = 1)
        else:
            print("Dropping [" + " ".join(truncated) +"] in Spectronaut split...")
            at_s = at_s.drop(truncateDrop(truncated), axis = 1)
            ce_s = ce_s.drop(truncateDrop(truncated), axis = 1)
            hs_s = hs_s.drop(truncateDrop(truncated), axis = 1)   
    return at_s, ce_s, hs_s

def normalizePerProteinAbundance(df):
    """
    Takes a species dataframe and normalizes per protein abundance.
    """    
    normalized_protein = df.div(df.sum(axis=1), axis=0)
    return normalized_protein

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

def mergeAllSpecies(at_s, ce_s, hs_s, at_t, ce_t, hs_t, maxNaN = 6, melt = True, truncated=  False):
    """
    Merges and normalizes alll the different species data frames into plotting format.
    maxNaN - maximal amount of NaN in a sample to keep the protein.
    Note: truncated can take subset of repl_cols instead of boolean.
    """
    if truncated != False:
        if truncated == True:
            repl_cols = ["S02","S03","S04","S05","S06","S07","S08","S09"]
            print("Replicates list: " + "[" + " ".join(repl_cols) +"]")
        else:
            repl_cols = truncated
            print("Replicates list: " + "[" + " ".join(repl_cols) +"]")
    else:
        repl_cols = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
        print("Replicates list: " + "[" + " ".join(repl_cols) +"]")
    normalized_protein_at_s = normalizePerProteinAbundance(at_s)
    normalized_protein_ce_s = normalizePerProteinAbundance(ce_s)
    normalized_protein_hs_s = normalizePerProteinAbundance(hs_s)

    sampgroup_at_s,sampgroup_hs_s,sampgroup_ce_s = pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
    
    for s in repl_cols:
        sampgroup_at_s[s] = adjustPandasNaNSum(normalized_protein_at_s, s, maxNaN)
        sampgroup_ce_s[s] = adjustPandasNaNSum(normalized_protein_ce_s, s, maxNaN)
        sampgroup_hs_s[s] = adjustPandasNaNSum(normalized_protein_hs_s, s, maxNaN)
        
    sampgroup_at_s['Species'] = ['A thaliana'] *len(normalized_protein_at_s)
    sampgroup_ce_s['Species'] = ['C elegans']  *len(normalized_protein_ce_s)
    sampgroup_hs_s['Species'] = ['H sapiens']  *len(normalized_protein_hs_s)
    sampgroup_at_s['Method']  = ['Spectronaut']*len(normalized_protein_at_s)
    sampgroup_ce_s['Method']  = ['Spectronaut']*len(normalized_protein_ce_s)
    sampgroup_hs_s['Method']  = ['Spectronaut']*len(normalized_protein_hs_s)
    
    sampgroup_spec_s = pd.concat([sampgroup_at_s,sampgroup_ce_s,sampgroup_hs_s])
    
    
    ######################################
    ## Merge normalized Triqler species ##
    ######################################
    
    normalized_protein_at_t = normalizePerProteinAbundance(at_t)
    normalized_protein_ce_t = normalizePerProteinAbundance(ce_t)
    normalized_protein_hs_t = normalizePerProteinAbundance(hs_t)
    
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
    
    #nn = sampgroup[repl_cols_s][sampgroup[repl_cols_s] > 0.0].min(numeric_only=True,axis=1)/2.0
    #nnn = pd.concat([nn]*sampgroup.shape[1],axis=1)
    #sampgroup = sampgroup.where(sampgroup!=0.0, nnn)
            
    
    sampgroup['Protein'] = range(len(sampgroup))
    #sampgroup['Protein'] = sampgroup.index
    #unmeltsamp = sampgroup.copy()
    if melt == True:
        sampgroup=pd.melt(sampgroup, id_vars=['Method', 'Species','Protein'], value_vars=repl_cols)
        sampgroup.rename(columns={'variable':'Sample','value':'Expression'}, inplace=True)
        return sampgroup
    else:
        return sampgroup

def mergeSpecies(at, ce, hs, method_name = "Spectronaut"):
    """
    Merges and normalizes alll the different species data frames into plotting format. 
    """
    repl_cols = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
    
    normalized_protein_at = normalizePerProteinAbundance(at)
    normalized_protein_ce = normalizePerProteinAbundance(ce)
    normalized_protein_hs = normalizePerProteinAbundance(hs)

    sampgroup_at,sampgroup_hs,sampgroup_ce = pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
    
    for s in repl_cols:
        sampgroup_at[s] = normalized_protein_at[s].sum(1)
        sampgroup_ce[s] = normalized_protein_ce[s].sum(1)
        sampgroup_hs[s] = normalized_protein_hs[s].sum(1)
        
    sampgroup_at['Species'] = ['A thaliana'] *len(normalized_protein_at)
    sampgroup_ce['Species'] = ['C elegans']  *len(normalized_protein_ce)
    sampgroup_hs['Species'] = ['H sapiens']  *len(normalized_protein_hs)
    sampgroup_at['Method']  = [method_name]*len(normalized_protein_at)
    sampgroup_ce['Method']  = [method_name]*len(normalized_protein_ce)
    sampgroup_hs['Method']  = [method_name]*len(normalized_protein_hs)
    
    sampgroup_spec = pd.concat([sampgroup_at,sampgroup_ce,sampgroup_hs])
    return sampgroup_spec

repl_cols = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
find = ["S02", "S03", "S05"] #<---------------------------------- HERE CKEEP WORKING

def findReplIdx(replicates):
    """
    returns replicates indices
    """
    allReplicates = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
    idxs = np.where(np.isin(allReplicates, replicates))
    return idxs[0]

def takeReplicatesMixtureAmount(mix, replicates):
    """
    return mixture corresponding the selected replicates.
    - mix - mixture from getMixtures
    - allReplicates - ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
    - find - subset of allReplicates
    """
    allReplicates = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
    mixture = np.take(mix, findReplIdx(replicates))
    return mixture

def truncateDrop(replicates):
    """
    Function to find which replicates to drop.
    """
    allReplicates = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
    dropRepl = np.take(allReplicates, np.where(~np.isin(allReplicates, replicates))).tolist()[0]
    return dropRepl

def getMixtures(truncated = False):
    """
    Defines the true mixture levels.
"""
    a_mix = np.array([0.5] * 10)
    a_mix[0] += 0.0001
    c_mix = np.array([0.5, 0.25, 0.125, 0.0625, 0.031, 0.0155, 0.008, 0.004, 0.002, 0.0])
    h_mix = np.array([0.0, 0.25, 0.375, 0.4375, 0.469, 0.4845, 0.492, 0.496, 0.498, 0.5])
    if truncated != False:
        if truncated == True:
            print("Getting mixtures for [S02 S03 S04 S05 S06 S07 S08 S09]")
            a_mix = a_mix[1:-1]
            c_mix = c_mix[1:-1]
            h_mix = h_mix[1:-1]
            a_mix = a_mix/a_mix.sum()
            c_mix = c_mix/c_mix.sum()
            h_mix = h_mix/h_mix.sum()
            mix = {'A thaliana':a_mix, 'C elegans': c_mix, 'H sapiens': h_mix}
        else:
            print("Getting mixtures for [" + " ".join(truncated) +"]")
            a_mix = takeReplicatesMixtureAmount(a_mix, truncated)
            c_mix = takeReplicatesMixtureAmount(c_mix, truncated)
            h_mix = takeReplicatesMixtureAmount(h_mix, truncated)
            a_mix = a_mix/a_mix.sum()
            c_mix = c_mix/c_mix.sum()
            h_mix = h_mix/h_mix.sum()
            mix = {'A thaliana':a_mix, 'C elegans': c_mix, 'H sapiens': h_mix}
            
    else:
        print("Getting mixtures for [S01 S02 S03 S04 S05 S06 S07 S08 S09 S10]")
        a_mix = a_mix/a_mix.sum()
        c_mix = c_mix/c_mix.sum()
        h_mix = h_mix/h_mix.sum()
        mix = {'A thaliana':a_mix, 'C elegans': c_mix, 'H sapiens': h_mix}
    
    return a_mix, c_mix, h_mix, mix


#####################
# ADJUST NON SHARED #
#####################
def adjustNonShared(libraryDirectory = "library/LKaell/", 
                 humanLib = "uniprot_sprot_2017-10-25_HUMAN_ISOFORMS.fasta",
                 arathLib = "uniprot_sprot_2018-01-24_ARATH_ISOFORMS.fasta",
                 caeelLib = "uniprot-elegans-filtered-organism__Caenorhabditis+elegans.fasta",
                 filename = "500-PSSS3-equ decoy_Report.xls",
                 outputname = "500-PSSS3-equ decoy_Report_nonShared_20190507.xls"):
    findUniqueProteins(inputName = filename, outputName = "proteins.csv")
    uniprot = getPSSS3Uniprot(libraryDirectory = libraryDirectory, 
                              humanLib = humanLib,
                              arathLib = arathLib,
                              caeelLib = caeelLib)
    singleProteins = findSingleProteins(filename="proteins.csv")
    shared, nonShared = findShared(filename = "proteins.csv") 
    
    filename = filename
    outputname = outputname
    
    f = open(filename, "r")
    new_file = open(outputname, "w")
    #line = f.readline().split("\t")
    #i = 0
    for line in f:
        #i+=1
        row = line.split("\t")
        if row[3] in nonShared:
            row[3] = row[3].split(";")[0]       
        new_line = "\t".join(row)
        new_file.write(new_line)
    
    ################################################################
    # Code for fixing laste 15 lines                               #
    # unknown why last 15 lines wont get appended on previous code #
    ################################################################
    
    def tail(filename, lines):
        import sys
        import os
        
        bufsize = 8192
        
        lines = lines #int(sys.argv[1])
        fname = filename #sys.argv[2]
        fsize = os.stat(fname).st_size
        
        iter = 0
        #with open(sys.argv[2]) as f:
        with open(fname) as f:
            if bufsize > fsize:
                bufsize = fsize-1
            data = []
            while True:
                iter +=1
                f.seek(fsize-bufsize*iter)
                data.extend(f.readlines())
                if len(data) >= lines or f.tell() == 0:
                    output = (''.join(data[-lines:]))
                    return output
                    break
        #return tailLines
    
    tailLines = tail(filename, 16)
    tailLines = tailLines.split("\n")#.split("\t")
    new_file = open(outputname, "a")
    for line in tailLines:
        if line == "":
            continue
        row = line.split("\t")
        if row[3] in nonShared:
            row[3] = row[3].split(";")[0]
        new_line = "\t".join(row)
        new_file.write(new_line)

#########################################
# PICKLE USED IN computeDiffExp script ##
#########################################
def removePickle(pickle_path = "pickled"):
    """
    Need to remove Pickle to use other parameters!
    """
    if not Path(pickle_path).exists():
        raise Exception ("No pickled folder exists!")
    answ = input("Delete pickled folder " +"(" + str(Path(pickle_path))+ ")" + " with contents? [Y/n]")
    while True:
        if answ == "Y":
            shutil.rmtree("pickled")
        elif answ == "n":
            break
        else:
            print("Y (delete) or n (don't delete)...")
            
def read_pickle_Triqler(triqlerFile, protein_id_fdr_treshold, pickle_path = "../data/pickled/"):
    result_triqler = Path(pickle_path + "triqler.pkl")
    if not result_triqler.exists():
        triqler = processTriqler(triqlerFile = triqlerFile, FDR_treshold = protein_id_fdr_treshold)
        if not Path(pickle_path).exists():
            os.mkdir(pickle_path)
        triqler.to_pickle(pickle_path + "triqler.pkl")
        pickle_params = open(pickle_path + "triqlerParams.txt", "w")
        pickle_params.write("triqlerFile: " + triqlerFile + "\n")
        pickle_params.write("protein_id_fdr_treshold: " + str(protein_id_fdr_treshold) + "\n")
        pickle_params.close()
    else:
        f = open(pickle_path + "triqlerParams.txt", "r")
        print("Reading pickled Triqler...")
        print("Triqler parameters")
        for i in f:
            print(i,  end = "")
        triqler = pd.read_pickle(result_triqler)
    return triqler 

def read_pickle_Spectronaut(spectronautFile, protein_id_fdr_treshold, impute, global_impute, pickle_path = "../data/pickled/"):
    result_spectronaut = Path(pickle_path + "spectronaut.pkl")
    if not result_spectronaut.exists():
        spectronaut = processSpectronaut(spectronautFile = spectronautFile, FDR_treshold = protein_id_fdr_treshold, impute = impute, global_impute = global_impute)
        if not Path(pickle_path).exists():
            os.mkdir(pickle_path)
        pickle_params = open(pickle_path + "spectronautParams.txt", "w")
        pickle_params.write("spectronautFile: " + spectronautFile + "\n")
        pickle_params.write("protein_id_fdr_treshold: " + str(protein_id_fdr_treshold) + "\n")
        pickle_params.write("Impute: " + str(impute) + "\n")
        pickle_params.write("Global_impute: " + str(global_impute))
        pickle_params.close()
        spectronaut.to_pickle(pickle_path + "spectronaut.pkl")
    else:
        f = open(pickle_path + "spectronautParams.txt", "r")
        print("Reading pickled Spectronaut...")
        print("Spectronaut parameters")
        for i in f:
            print(i, end = "")
        spectronaut = pd.read_pickle(result_spectronaut)
    return spectronaut