#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 10:02:39 2019

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
import pandas as pd
import numpy as np
from scipy import stats
import scipy as sp
import subprocess
import matplotlib.pyplot as plt 
import seaborn as sns

########
# MAIN #
########

def main():
    parser = argparse.ArgumentParser(prog='PROG')
    subparsers = parser.add_subparsers(help = 'sub-command help', dest = "subparser_name")
    
    # create parser for the preprocess shared command
    parser_preprocess = subparsers.add_parser('preprocessShared', help = "Adjusts Shared peptides. If peptides are shared but unique, they become adjusted so that they are non-shared.")
    parser_preprocess.add_argument("--input", metavar = "INPUT_FILE", help = "PSSS3 file.")
    parser_preprocess.add_argument("--output", metavar = "OUTPUT_NAME", help = "Output name.")
    parser_preprocess.add_argument("--libDir", metavar = "LIB_DIR", help = "Library directory.")
    parser_preprocess.add_argument("--humanLib", metavar = "HUMAN_LIB", help =".fasta file for HUMAN proteins.")
    parser_preprocess.add_argument("--arathLib", metavar = "ARATH_LIB", help = ".fasta file for ARATH proteins.")
    parser_preprocess.add_argument("--caeelLib", metavar = "CAEEL_LIB", help = ".fasta file for CAEEL proteins.")
    
    # parser for convertSpectronaut2Triqler
    parser_convertSpectronaut2Triqler = subparsers.add_parser("convertSpectronaut2Triqler", help = "Converts PSSS3 to triqler input format")
    parser_convertSpectronaut2Triqler.add_argument("--input", metavar = "INPUT_FILE", help = "Input Spectronaut file (PSSS3)") 
    parser_convertSpectronaut2Triqler.add_argument("--output", metavar = "OUTPUT_NAME", help = "Output name.")

    # parser for callTriqler
    parser_triqler = subparsers.add_parser("triqler", help = "Calls triqler")
    parser_triqler.add_argument("--input", metavar = "INPUT_FILE", help = "triqler input file.")
    parser_triqler.add_argument("--output", metavar = "OUTPUT_NAME", help = "Output name.")
    parser_triqler.add_argument("--fold_change_eval", metavar = "F", help = "log2 fold change evaluation treshold. (default = 1.0)")
    parser_triqler.add_argument("--decoy_pattern", metavar = "P", help = "Prefix for decoy proteins.")
    parser_triqler.add_argument("--min_samples", metavar = "N", help = "Minimum number of samples a peptide need to be quantified in. (default: 2)")
    
    # parser for convertPSSS3toReformattedPSSS3
    parser_convertPSSS3toReformattedPSSS3 = subparsers.add_parser("reformattedPSSS3", help = "Reformattes the PSSS3 to comparisonPlot and diffExp format. Procedure takes some time!")
    parser_convertPSSS3toReformattedPSSS3.add_argument("--input", metavar = "INPUT_FILE", help = "PSSS3 file.")
    parser_convertPSSS3toReformattedPSSS3.add_argument("--output", metavar = "OUTPUT_NAME", help = "Output name.")
    
    # parser for main_singleImpute
    parser_mainSingle = subparsers.add_parser("comparisonPlot", help = "Generate comparison plots.")
    parser_mainSingle.add_argument("--PSSS3_input", metavar = "PSSS3_INPUT_FILE", default = None, help = "Input reformatted PSSS3 file.")
    parser_mainSingle.add_argument("--triqler_input", metavar = "TRIQLER_RESULTS", default = None, help = "<x>vs<y> triqler result file.")
    parser_mainSingle.add_argument("--impute_method", metavar = "IMPUTE_STRING", default = None, help = 'Imputation methods to be used. "min" = minimal, "max" = maximal, "mean" = mean, "median" = median, None = no impute.')
    parser_mainSingle.add_argument("--fdr_treshold", metavar = "F", default = 0.01, help = "FDR treshold for comparison.")
    parser_mainSingle.add_argument("--sample_list", metavar = "LIST", nargs="+", default = ["S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10"], help = "List of samples to compare.")
    parser_mainSingle.add_argument("--triqlerExponential", metavar = "N", default = False, help = "Take the exponent of triqler values to unlog. Default == False; keep non-exponential values.")
    parser_mainSingle.add_argument("--global_impute", metavar = "BOOLEAN", default = True, help = "Impute globally.")
    parser_mainSingle.add_argument("--maxNaN", metavar = "N", default = 6, help = "N NaNs to include protein in comparison.")
    parser_mainSingle.add_argument("--outputSuffix", metavar = "OUTPUT_SUFFIX", default = "plot", help = "Output plot suffices.")
    # parser for mainDiffExp
    parser_mainDiffExp = subparsers.add_parser("diffExp", help = "Generate diff expression comparison plot.")
    parser_mainDiffExp.add_argument("--resultDir", metavar = "TRIQLER_RESULT_DIR", help = "Directory of triqler results.")
    parser_mainDiffExp.add_argument("--PSSS3", metavar = "PSSS3_INPUT_FILE", help = "Reformatted PSSS3 file.")
    parser_mainDiffExp.add_argument("--impute_method", metavar = "IMPUTE_STRING", default = None, help = 'Imputation methods to be used. "min" = minimal, "max" = maximal, "mean" = mean, "median" = median, None = no impute.')
    parser_mainDiffExp.add_argument("--fdr_treshold", metavar = "F", default = 0.01, help = "FDR treshold for comparison.")
    parser_mainDiffExp.add_argument("--global_impute", metavar = "BOOLEAN", default = True, help = "Impute globally.")
    parser_mainDiffExp.add_argument("--figName", metavar = "STRING", default = "plot", help = "Output suffix.")
   
    args = parser.parse_args()
    if args.subparser_name == 'preprocessShared':
        adjustNonShared(libraryDirectory = args.libDir,
                        humanLib = args.humanLib,
                        arathLib = args.arathLib,
                        caeelLib = args.caeelLib,
                        filename = args.input,
                        outputname = args.output)
    elif args.subparser_name == 'convertSpectronaut2Triqler':
        convertSpectronaut2Triqler(inputArg = args.input, 
                                   outputName = args.output)        
    elif args.subparser_name == 'triqler':
        callTriqler(inputArg = args.input,
                    outputName = args.output,
                    fold_change_eval = args.fold_change_eval,
                    decoy_pattern = args.decoy_pattern)
    elif args.subparser_name == 'comparisonPlot':
        main_singleImpute(impute_method = args.impute_method,
                          FDR_treshold = float(args.fdr_treshold),
                          truncated = args.sample_list,
                          triqlerExponential = float(args.triqlerExponential),
                          global_impute = args.global_impute,
                          maxNaN = int(float(args.maxNaN)),
                          PSSS3Input = args.PSSS3_input,
                          triqlerInput = args.triqler_input,
                          title_suffix = args.outputSuffix)
    elif args.subparser_name == 'diffExp':
        mainDiffExp(fileDir = args.resultDir,
                    PSSS3Input = args.PSSS3,
                    impute_method = args.impute_method,
                    FDR_treshold = float(args.fdr_treshold),
                    global_impute = args.global_impute,
                    figName = args.figName)
    
def testImputeDf():
    return pd.DataFrame([[1, np.nan, np.nan, np.nan], 
                               [np.nan, 2, 2, np.nan], [3, 3, 3, np.nan],
                               [np.nan, 4, 4, 4], 
                               [np.nan, np.nan, np.nan, np.nan]])
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
        elif imput == None:
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

def getProteinCountPlot(at_t, ce_t, hs_t, at_s, ce_s, hs_s):
    protein_count = pd.DataFrame(columns = ['species', 'method'])
    protein_count['species'] = (['A thaliana']*len(at_t) + ['C elegans']*len(ce_t) + ["H sapiens"]*len(hs_t) + 
                 ['A thaliana']*len(at_s) + ['C elegans']*len(ce_s) + ["H sapiens"]*len(hs_s))
    protein_count['method'] = (['Triqler']*len(at_t) + ['Triqler']*len(ce_t) + ['Triqler']*len(hs_t) + 
                 ['Spectronaut']*len(at_s) + ['Spectronaut']*len(ce_s) + ['Spectronaut']*len(hs_s))

    fig, ax = plt.subplots(figsize=(15, 9))
    sns.countplot(x='species', hue='method', data=protein_count)
    sns.despine(offset=10, trim=True);
    fig.savefig("proteinCount.png")

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


##############
### MAIN #####
##############

def main_singleImpute(impute_method = "min", FDR_treshold = 0.01, truncated = False, triqlerExponential = 2, global_impute = True, maxNaN = 6,
                      PSSS3Input = "/home/ptruong/Git/libvar/data/PSSS3/500-PSSS3-equ decoy_Report_new/500-PSSS3-raw-reformatted_dropna_dropdup_decoy.csv",
                    triqlerInput = "/home/ptruong/Git/libvar/data/PSSS3/500-PSSS3-equ decoy_Report_new/triqlerOutput/20190312_unfiltered/proteins.1vs2.tsv",
                    title_suffix = "plot"):
    print("Starting single Impute procedure...")
    if impute_method == None:
        print("No Imputation!")
    else:
        print("Imputation method: " + impute_method)
    if impute_method == None:
        pass
    else:
        print("Global impute: " + str(global_impute))
    print("Processing Spectronaut file: " + PSSS3Input)
    spectronaut =  processSpectronaut(spectronautFile = PSSS3Input,
                       FDR_treshold = FDR_treshold, impute = impute_method,
                       global_impute = global_impute)
    
    print("Splitting Spectronaut file by species...")
    at_s, ce_s, hs_s = splitSpectronautBySpecies(spectronaut, truncated = truncated)
    
    print("Processing Triqler file:" + triqlerInput)
    triqler = processTriqler(triqlerFile = triqlerInput, FDR_treshold = FDR_treshold)
    print("Splitting Triqler file by species...")
    at_t, ce_t, hs_t = splitTriqlerBySpecies(triqler, exponential = triqlerExponential, truncated = truncated)
    
    print("Merging samples")
    sampgroup = mergeAllSpecies(at_s, ce_s, hs_s, at_t, ce_t, hs_t, maxNaN = maxNaN, melt = True, truncated = truncated)
    unmeltsamp = mergeAllSpecies(at_s, ce_s, hs_s, at_t, ce_t, hs_t, maxNaN = maxNaN, melt = False, truncated = truncated)
    
    print("Generating protein count plot...")
    #getProteinCountPlot(at_t, ce_t, hs_t, at_s, ce_s, hs_s)
    print("Generating comparison plots...")
    comparisonPlot(sampgroup, "H sapeins",  hardRange = False, truncated = truncated, title = title_suffix)
    comparisonPlot(sampgroup, "C elegans",  hardRange = False, truncated = truncated, title = title_suffix)
    comparisonPlot(sampgroup, "A thaliana",  hardRange = False, truncated = truncated, title = title_suffix)
    
    #print("Generating correlation treshold plots...")
    #p_r_plot(unmeltsamp, 'H sapiens', "Hsapiens_CorrTresh.png")
    #p_r_plot(unmeltsamp,'C elegans', "Celegans_CorrTresh.png")
    #p_r_plot(unmeltsamp,'A thaliana', "Athaliana_CorrTresh.png")
    print("Done!")
    
def mainMultiImpute(triqlerExponential = 2,
                    PSSS3Input = "/home/ptruong/Git/libvar/data/PSSS3/500-PSSS3-equ decoy_Report_new/500-PSSS3-raw-reformatted_dropna_dropdup_decoy.csv",
                    triqlerInput = "/home/ptruong/Git/libvar/data/PSSS3/500-PSSS3-equ decoy_Report_new/triqlerOutput/20190312_unfiltered/proteins.1vs2.tsv"):
    def spectronautImpute(impute_method, FDR_treshold = 0.01):
        spectronaut =  processSpectronaut(spectronautFile = PSSS3Input,
                           FDR_treshold = FDR_treshold, impute = impute_method)
        at_s, ce_s, hs_s = splitSpectronautBySpecies(spectronaut)
        return at_s, ce_s, hs_s
    
    triqler = processTriqler(triqlerFile = triqlerInput, FDR_treshold = 0.01)
    
    #imputation_methods =  ["min", "max", "mean", "median"]
    imputation_methods =  ["min", "mean", "median"]
    at_s = []
    ce_s = []
    hs_s = []
    
    at_t = []
    ce_t = []
    hs_t = []
    for i in imputation_methods:
        at_s_i, ce_s_i, hs_s_i = spectronautImpute(impute_method = i, FDR_treshold = 0.01)
        at_s.append(at_s_i)
        ce_s.append(ce_s_i)
        hs_s.append(hs_s_i)
    
        at_t_i, ce_t_i, hs_t_i = splitTriqlerBySpecies(triqler, exponential = triqlerExponential)
        at_t.append(at_t_i)
        ce_t.append(ce_t_i)
        hs_t.append(hs_t_i)
        
    at_s = pd.concat(at_s)
    ce_s = pd.concat(ce_s)
    hs_s = pd.concat(hs_s)
    
    at_t = pd.concat(at_t)
    ce_t = pd.concat(ce_t)
    hs_t = pd.concat(hs_t)
    
    #at_t, ce_t, hs_t = splitTriqlerBySpecies(triqler)
    
    sampgroup = mergeAllSpecies(at_s, ce_s, hs_s, at_t, ce_t, hs_t, melt = True)
    unmeltsamp = mergeAllSpecies(at_s, ce_s, hs_s, at_t, ce_t, hs_t, melt = False)
    
    
    getProteinCountPlot(at_t, ce_t, hs_t, at_s, ce_s, hs_s)
    comparisonPlot(sampgroup, "H sapeins",  hardRange = True)
    comparisonPlot(sampgroup, "C elegans",  hardRange = True)
    comparisonPlot(sampgroup, "A thaliana",  hardRange = True)

    #p_r_plot(unmeltsamp, 'H sapiens', "Hsapiens_CorrTresh.png")
    #p_r_plot(unmeltsamp,'C elegans', "Celegans_CorrTresh.png")
    #p_r_plot(unmeltsamp,'A thaliana', "Athaliana_CorrTresh.png")

"""
def mainMultiImpute(PSSS3Input = "/home/ptruong/Git/libvar/data/PSSS3/500-PSSS3-equ decoy_Report_new/500-PSSS3-raw-reformatted_dropna_dropdup_decoy.csv",
                    triqlerInput = "/home/ptruong/Git/libvar/data/PSSS3/500-PSSS3-equ decoy_Report_new/triqlerOutput/20190312_unfiltered/proteins.1vs2.tsv"):
    def spectronautImpute(impute_method, FDR_treshold = 0.01):
        spectronaut =  processSpectronaut(spectronautFile = PSSS3Input,
                           FDR_treshold = FDR_treshold, impute = impute_method)
        at_s, ce_s, hs_s = splitSpectronautBySpecies(spectronaut)
        return at_s, ce_s, hs_s
    
    triqler = processTriqler(triqlerFile = triqlerInput, FDR_treshold = 0.01)
    
    #imputation_methods =  ["min", "max", "mean", "median"]
    imputation_methods =  ["min", "mean", "median"]
    at_s = []
    ce_s = []
    hs_s = []
    
    for i in imputation_methods:
        at_s_i, ce_s_i, hs_s_i = spectronautImpute(impute_method = i, FDR_treshold = 0.01)
        at_s.append(at_s_i)
        ce_s.append(ce_s_i)
        hs_s.append(hs_s_i)
    
    at_t, ce_t, hs_t = splitTriqlerBySpecies(triqler)
    
    sampgroup_s = []    
    for i in range(len(imputation_methods)):
        sampgroup_i =  mergeSpecies(at_s[i], ce_s[i], hs_s[i], method_name = "Spectronaut_"+imputation_methods[i])
        sampgroup_s.append(sampgroup_i)
        
    sampgroup = mergeAllSpecies(at_s, ce_s, hs_s, at_t, ce_t, hs_t, melt = True)
    unmeltsamp = mergeAllSpecies(at_s, ce_s, hs_s, at_t, ce_t, hs_t, melt = False)

    #getProteinCountPlot(at_t, ce_t, hs_t, at_s, ce_s, hs_s)
    comparisonPlot(sampgroup, "H sapeins",  hardRange = True)
    comparisonPlot(sampgroup, "C elegans",  hardRange = True)
    comparisonPlot(sampgroup, "A thaliana",  hardRange = True)

    #p_r_plot(unmeltsamp, 'H sapiens', "Hsapiens_CorrTresh.png")
    #p_r_plot(unmeltsamp,'C elegans', "Celegans_CorrTresh.png")
    #p_r_plot(unmeltsamp,'A thaliana', "Athaliana_CorrTresh.png")
"""

def mainDiffExp(fileDir, PSSS3Input, impute_method = None, FDR_treshold = 0.01, global_impute = True, figName = "plot"):
    """
    Generates amount of differential expression plot.
    
    fileDir should be map containing all triqler generated outputs.
    PSSS3Input - PSSS3 file (parsed with parsePSSS3equDecoy() or convertPSSS3toReformattedPSSS3()).
    """
    print("Performining analysis diff Exp analysis")
    fileDir = fileDir
    try:
        dirName = fileDir + "/results"
        os.mkdir(dirName)
        print("Directory " + dirName + " created!")
    except:
        print("Directory " + dirName + "already exists!")
        
    print("Analyzing " + fileDir)
    fileDirectory = fileDir
    #print("params")
    #print("unlogTriqler: " + str(unlogTriqler))
    #print("logSpectronaut: " + str(logSpectronaut))
    #print("truncate: " + str(dropS01S10))
    #print("FDR_treshold: " + str(FDR_treshold))

    #basefig = ""
    #if unlogTriqler == True:
    #    basefig += "_unlogTriqler"
    #if logSpectronaut == True:
    #    basefig += "_logSpectronaut"
    #if dropS01S10 == True:
    #    basefig += "_trunc"
    #basefig += "_FDR" + str(FDR_treshold).replace(".", "")
    
    df_psss3 = processSpectronaut(spectronautFile = PSSS3Input,
                                  FDR_treshold = FDR_treshold, impute = impute_method,
                                  global_impute = global_impute)
    ###########################
    # DiffExp COMPARISONS #####
    ###########################
    print("Generating DiffExp plot...")
    generateDiffExpPlot(fileDirectory = fileDirectory,
                        df_psss3 = df_psss3,
                        FDR_treshold = FDR_treshold,
                        logSpectronaut = logSpectronaut,
                        figName = figName)
##############
# Converters #
##############

def convertSpectronaut2Triqler(inputArg = "500-PSSS3-equ\ decoy_Report.xls", outputName = "PSSS3_triqlerFormatted.csv"):
    """
    Calls convert_spectronaut_to_triqler_input.py to convert PSSS3 file to triqler input formatted file.
    """
    inputArg = inputArg
    outputName = outputName
    command = "python2 convert_spectronaut_to_triqler_input.py" + " " + inputArg + " " + outputName
    os.system(command)

def callTriqler(inputArg = "PSSS3_triqlerFormatted.csv", 
                out_file = "triqlerOutput/protein.tsv",
                fold_change_eval = 1.0,
                decoy_pattern = "decoy_",
                min_samples = 2):
    inputArg = "PSSS3_triqlerFormatted.csv"
    command = ("python -m triqler" + " --out_file " + out_file + " " + 
               "--fold_change_eval " + str(fold_change_eval) + " " +
               "--decoy_pattern " + decoy_pattern + " " +
               "--min_samples " + str(min_samples) + " " + inputArg)
    os.system(command)

def convertPSSS3toReformattedPSSS3(inputArg = '500-PSSS3-equ decoy_Report.xls', outputName = "500-PSSS3-raw-reformatted_dropna_dropdup_decoy.csv"):
    inputArg = inputArg
    outputName = outputName
    df = parsePSSS3equDecoy(filename = inputArg) # VERY SLOW
    df = df.drop_duplicates()
    df.to_csv(outputName, sep = "\t", index = False)
    del df    

if __name__ == "__main__":
    main()

    
    
    
