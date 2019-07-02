#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 10:02:39 2019

@author: ptruong
"""
import sys 

from proteinSearch import *

from log2FC import *
from triqlerParser import *
from triqlerProcessor import * 
from utils import * 
from parseReformattedPSSS3equDecoy import * # Needs to tidy and place functions in right place!
from proc import *
from plotting import * 
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
    
    # parser for logfc between samples
    parser_log2FC = subparsers.add_parser("log2fc", help = 'Computes log2FC between A and B')
    parser_log2FC.add_argument("--a", metavar = "F", help = "Value A for log2fc(A,B)")
    parser_log2FC.add_argument("--b", metavar = "F", help = "Value B For log2fc(A,B)")
    
    # parser for get the "True" log2fc between samples.
    parser_log2FC_true = subparsers.add_parser("true_fc", help = 'Get the "True" log2fc between samples.') 
    parser_log2FC_true.add_argument("--s1", metavar = "F", help = "Sample 1-10 for S01,..., S10")
    parser_log2FC_true.add_argument("--s2", metavar = "F", help = "Sample 1-10 for S01,..., S10")
    
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
    elif args.subparser_name == "log2fc":
        print(log2fc(a = float(args.a), b = float(args.b)))
    elif args.subparser_name == "true_fc":
        print(get_fc_between_samples(s1 = int(args.s1), s2 = int(args.s2)))
        
def testImputeDf():
    return pd.DataFrame([[1, np.nan, np.nan, np.nan], 
                               [np.nan, 2, 2, np.nan], [3, 3, 3, np.nan],
                               [np.nan, 4, 4, 4], 
                               [np.nan, np.nan, np.nan, np.nan]])


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

##################
# computeDiffExp #
##################

if __name__ == "__main__":
    main()

    
    
    
