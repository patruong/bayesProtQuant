#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 10:33:30 2019

@author: ptruong
"""

import os
import shutil
from pathlib import Path
import pandas as pd
from scipy import stats
import scipy as sp
from qvalue import *
from proc import *
import argparse


def main():
    parser = argparse.ArgumentParser(prog='PROG')
    subparsers = parser.add_subparsers(help = 'sub-command help', dest = "subparser_name")
    
    # create parser for ttest triqler
    parser_ttest_triqler = subparsers.add_parser("ttest_triqler", help = "Computes differentially expressed protein from triqler results with t-test.")
    parser_ttest_triqler.add_argument("--input", metavar = "INPUT_FILE", help = "triqler results <x>vs<y>.tsv file.")
    parser_ttest_triqler.add_argument("--output", metavar = "OUTPUT_FILE", default = False, help = "output file name.") 
    parser_ttest_triqler.add_argument("--proteinFDR", metavar = "F", default = 0.01, help = "protein id fdr treshold.")
    parser_ttest_triqler.add_argument("--specie", metavar = "S", help = "One of following species, ARATH, HUMAN, CAEEL.")
    parser_ttest_triqler.add_argument("--sample1", metavar = "S", help = "One of following [S01, S02 ..., S10]")
    parser_ttest_triqler.add_argument("--sample2", metavar = "S", help = "One of following [S01, S02 ..., S10]")
    parser_ttest_triqler.add_argument("--diffExpFDR", metavar = "F", default = 0.05, help = "Q-value treshold for amount of differentially expressed proteins.")
    parser_ttest_triqler.add_argument("--exponential", metavar = "I", default = False, help = "Take exponential of triqler absolute protein quantifications") 
    
    # create parser for ttest spectronaut
    parser_ttest_spectronaut = subparsers.add_parser("ttest_spectronaut", help = 'Computes differentially expressed protein from spectronaut results with t-test.')
    parser_ttest_spectronaut.add_argument("--input", metavar = "INPUT_FILE", help = "reformatted psss3 file.")
    parser_ttest_spectronaut.add_argument("--output", metavar = "OUTPUT_FILE", default = False, help = "output file name.")
    parser_ttest_spectronaut.add_argument("--proteinFDR", metavar = "F", default = 0.01, help = "protein id fdr treshold.")
    parser_ttest_spectronaut.add_argument("--specie", metavar = "S", help = "One of following species, ARATH, HUMAN, CAEEL.")
    parser_ttest_spectronaut.add_argument("--sample1", metavar = "S", help = "One of following [S01, S02 ..., S10]")
    parser_ttest_spectronaut.add_argument("--sample2", metavar = "S", help = "One of following [S01, S02 ..., S10]")
    parser_ttest_spectronaut.add_argument("--diffExpFDR", metavar = "F", default = 0.05, help = "Q-value treshold for amount of differentially expressed proteins.")
    parser_ttest_spectronaut.add_argument("--impute", metavar = "S", default = None, help = "Imputation method; mean, median, min, max") 
    parser_ttest_spectronaut.add_argument("--global_impute", metavar = "Boolean", default = True, help = "Impute globally (True) or locally (False).") 
    
    # create parser for log fold change triqler
    parser_lfc_triqler = subparsers.add_parser("lfc_triqler", help = 'Computes differentially expressed protein from triqler results with log fold change q-values.')
    parser_lfc_triqler.add_argument("--input", metavar = "INPUT_FILE", help = "triqler results <x>vs<y>.tsv file.")
    parser_lfc_triqler.add_argument("--output", metavar = "OUTPUT_FILE", default = False, help = "output file name.")
    parser_lfc_triqler.add_argument("--proteinFDR", metavar = "F", default = 0.01, help = "protein id fdr treshold.")
    parser_lfc_triqler.add_argument("--specie", metavar = "S", help = "One of following species, ARATH, HUMAN, CAEEL.")
    parser_lfc_triqler.add_argument("--diffExpFDR", metavar = "F", default = 0.05, help = "Q-value treshold for amount of differentially expressed proteins.")
    
    # create parser for removing pickle folder.
    parser_removePickle = subparsers.add_parser("removePickle", help = 'Removes pickled results. Is required to reset parameters for t-test computations.')

    args = parser.parse_args()
    if args.subparser_name == "ttest_triqler":
        n_diffExp = ttest_triqler(triqlerFile = args.input,
                                  protein_id_fdr_treshold = float(args.proteinFDR),
                                  specie = args.specie,
                                  sample1 = args.sample1,
                                  sample2 = args.sample2,
                                  fdr_treshold = float(args.diffExpFDR),
                                  exponential = int(args.exponential))
        if not args.output == False:
            writeResults(n_diffExp, output_file = args.output)
        else:
            print("Differentially expressed proteins: " + n_diffExp)
    elif args.subparser_name == "ttest_spectronaut":
        n_diffExp = ttest_spectronaut(spectronautFile = args.input,
                                      protein_id_fdr_treshold = float(args.proteinFDR),
                                      specie = args.specie, 
                                      sample1 = args.sample1,
                                      sample2 = args.sample2,
                                      fdr_treshold = float(args.diffExpFDR),
                                      impute = args.impute,
                                      global_impute = args.global_impute)
        if not args.output == False:
            writeResults(n_diffExp, output_file = args.output)
        else:
            print("Differentially expressed proteins: " + n_diffExp)
    elif args.subparser_name == "lfc_triqler":
        n_diffExp = get_qvalue_diff_exp_triqler(triqler_result = args.input,
                                                protein_id_fdr = float(args.proteinFDR),
                                                specie = args.specie,
                                                diff_exp_fdr_treshold = float(args.diffExpFDR))
        if not args.output == False:
            writeResults(n_diffExp, output_file = args.output)
        else:
            print("Differentially expressed proteins: " + n_diffExp)
    elif args.subparser_name == "removePickle":
        removePickle()

def removePickle():
    """
    Need to remove Pickle to use other parameters!
    """
    if not Path("pickled").exists():
        raise Exception ("No pickled folder exists!")
    answ = input("Delete pickled folder with contents? [Y/n]")
    while True:
        if answ == "Y":
            shutil.rmtree("pickled")
        elif answ == "n":
            break
        else:
            print("Y (delete) or n (don't delete)...")
    
def readTriqler(triqlerFile, FDR_treshold):
    result_triqler = Path("pickled/triqler.pkl")
    if not result_triqler.exists():
        triqler = processTriqler(triqlerFile = triqlerFile, FDR_treshold = FDR_treshold)
        if not Path("pickled").exists():
            os.mkdir("pickled")
        triqler.to_pickle("pickled/triqler.pkl")
        pickle_params = open("pickled/triqlerParams.txt", "w")
        pickle_params.write("triqlerFile: " + triqlerFile + "\n")
        pickle_params.write("FDR_treshold: " + str(FDR_treshold) + "\n")
        pickle_params.close()
    else:
        f = open("pickled/triqlerParams.txt", "r")
        print("Reading pickled Triqler...")
        print("Triqler parameters")
        for i in f:
            print(i,  end = "")
        triqler = pd.read_pickle(result_triqler)
    return triqler 

def readSpectronaut(spectronautFile, FDR_treshold, impute, global_impute):
    result_spectronaut = Path("pickled/spectronaut.pkl")
    if not result_spectronaut.exists():
        spectronaut = processSpectronaut(spectronautFile = spectronautFile, FDR_treshold = FDR_treshold, impute = impute, global_impute = global_impute)
        if not Path("pickled").exists():
            os.mkdir("pickled")
        pickle_params = open("pickled/spectronautParams.txt", "w")
        pickle_params.write("spectronautFile: " + spectronautFile + "\n")
        pickle_params.write("FDR_treshold: " + str(FDR_treshold) + "\n")
        pickle_params.write("Impute: " + str(impute) + "\n")
        pickle_params.write("Global_impute: " + str(global_impute))
        pickle_params.close()
        spectronaut.to_pickle("pickled/spectronaut.pkl")
    else:
        f = open("pickled/spectronautParams.txt", "r")
        print("Reading pickled Spectronaut...")
        print("Spectronaut parameters")
        for i in f:
            print(i, end = "")
        spectronaut = pd.read_pickle(result_spectronaut)
    return spectronaut

def computeDiffExp(sample1, sample2):
    """
    t-test differential expression.
    """
    if not (sample1.index.tolist() == sample2.index.tolist()):
        raise Exception("matching error!")
    pVals = []
    for i in range(len(sample1)):
        val1 = sample1.iloc[i].values
        val1 = val1[np.logical_not(np.isnan(val1))]
        val2 = sample2.iloc[i].values
        val2 = val2[np.logical_not(np.isnan(val2))]
        pVal = stats.ttest_ind(val1, val2).pvalue
        pVals.append(pVal)
    pVals = np.array(pVals)
    prots = sample1.index
    pVals = pd.DataFrame(pVals, index = prots)#, columns = ["S02"+"_"+"S06"])
    qVals = qvalues(pVals)
    return qVals

"""
triqlerFile = "../data/triqlerOutput_noShared_largeScaleOptimizations_minSamp15/proteins.1vs2.tsv",
             FDR_treshold = 0.01, # Protein ID treshold.
             exponential = False,
             truncated =  False,
             spectronautFile = "../data/500-PSSS3-raw-reformatted_dropna_dropdup_decoy_nonShared_again.csv",
             impute = None,
             global_impute = False):     
"""

    
###########################
# q-value triqler diffExp #
###########################
"""
triqler_result = "proteins.2vs6.tsv"
FDR_treshold = 0.01 # protein id FDR
specie = ["ARATH", "CAEEL", "HUMAN"]
"""



##################
# main functions #
##################
"""
def readData(triqlerFile = "../data/triqlerOutput_noShared_largeScaleOptimizations_minSamp15/proteins.1vs2.tsv",
             FDR_treshold = 0.01, # Protein ID treshold.
             exponential = False,
             truncated =  False,
             spectronautFile = "../data/500-PSSS3-raw-reformatted_dropna_dropdup_decoy_nonShared_again.csv",
             impute = None,
             global_impute = False):     
    triqler = readTriqler(triqlerFile = triqlerFile, FDR_treshold = FDR_treshold)
    spectronaut = readSpectronaut(spectronautFile = spectronautFile, FDR_treshold = FDR_treshold, impute = impute, global_impute = global_impute)
    at_t, ce_t, hs_t = splitTriqlerBySpecies(triqler, exponential = exponential, truncated = truncated)
    at_s, ce_s, hs_s = splitSpectronautBySpecies(spectronaut, truncated = truncated)
    return triqler, spectronaut, at_t, ce_t, hs_t, at_s, ce_s, hs_s
"""
def ttest_triqler(triqlerFile, protein_id_fdr_treshold, specie, sample1, sample2, fdr_treshold, exponential):
    if specie not in ["HUMAN", "ARATH", "CAEEL"]:
        raise Exception ("species should be HUMAN, ARATH or CAEEL!")
    triqler = readTriqler(triqlerFile, protein_id_fdr_treshold)
    at_t, ce_t, hs_t = splitTriqlerBySpecies(triqler, exponential = exponential, truncated = False)
    if specie == "HUMAN":
        df = hs_t
    elif specie == "ARATH":
        df = at_t
    elif specie == "CAEEL":
        df = ce_t
    s1 = df[sample1]
    s2 = df[sample2]
    qVals = computeDiffExp(s1, s2)
    n_diffExp = (qVals < fdr_treshold).sum().values[0]
    return n_diffExp

def ttest_spectronaut(spectronautFile, protein_id_fdr_treshold, specie, sample1, sample2, fdr_treshold, impute, global_impute):
    if specie not in ["HUMAN", "ARATH", "CAEEL"]:
        raise Exception ("species should be HUMAN, ARATH or CAEEL!")
    spectronaut = readSpectronaut(spectronautFile = spectronautFile, FDR_treshold = protein_id_fdr_treshold, impute = impute, global_impute = global_impute)
    at_s, ce_s, hs_s = splitSpectronautBySpecies(spectronaut, truncated = False)
    if specie == "HUMAN":
        df = hs_s
    elif specie == "ARATH":
        df = at_s
    elif specie == "CAEEL":
        df = ce_s
    s1 = df[sample1]
    s2 = df[sample2]
    qVals = computeDiffExp(s1, s2)
    n_diffExp = (qVals < fdr_treshold).sum().values[0]
    return n_diffExp

def get_qvalue_diff_exp_triqler(triqler_result, protein_id_fdr, specie, diff_exp_fdr_treshold):
    """
    triqler_result = triqler <x>vs<y>.tsv file (e.g. "proteins.2vs6.tsv")
    protein_id_fdr = protein id fdr tresholding.
    specie = ["ARATH", "CAEEL", "HUMAN"]
    """
    triqler = table2Df(triqler_result, peptideSeperator = True)
    triqler = triqler2Numeric(triqler)
    triqler = getTriqlerDetailed(triqler)
    
    triqler.sort_values(by = "protein_id_posterior_error_prob", inplace = True)
    triqler["FDR"] = triqler["protein_id_posterior_error_prob"].expanding().mean()
    triqler = triqler[triqler["FDR"] < protein_id_fdr]
    
    #df_triqler = addSpecieProteinColumns(triqler, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
    df_triqler = triqler
    runs = getRuns(triqler, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
    conditionRuns = getConditionRuns(runs)
    unique_run = conditionRuns.run.unique()
    run_id = getConditionRunId(triqler, unique_run[0])
    
    n_diffExp = (df_triqler[df_triqler["specie"] == specie].q_value < diff_exp_fdr_treshold).sum().values[0]
    return n_diffExp

def writeResults(n_diffExp, output_file = "result.csv"):
    print("Differentially expressed proteins: " + str(n_diffExp))
    f = open(output_file, "w")
    f.write("n_diffExp\t"+str(n_diffExp))
    f.close()

if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    main()







