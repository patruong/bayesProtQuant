#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 15:45:41 2019

@author: ptruong
"""
from plotting import *
from proc import *
from func import *

s1 = "S02"
s2 = "S06"
logFunc = np.log10
side = "two"
fc_treshold = 0.8
fdr_treshold = 0.05
col_diffExp = "red"
col_not_diffExp = "blue"
compute_FDR = compute_logP #NOTE: function computes -log(FDR)

# NOTE CHECK pi0 for q-value computations!


# Read data
def getTriqler_pickle(triqlerFile, protein_id_fdr_treshold, pickle_path = "../data/pickled/"):
    triqler = read_pickle_Triqler(triqlerFile, protein_id_fdr_treshold, pickle_path = pickle_path)
    at_t, ce_t, hs_t = splitTriqlerBySpecies(triqler, exponential = False, truncated = False)
    return triqler, at_t, ce_t, hs_t

def getSpectronaut_pickle(spectronautFile, protein_id_fdr_treshold, impute = None, global_impute = True, pickle_path = "../data/pickled/"):
    spectronaut =  read_pickle_Spectronaut(spectronautFile, protein_id_fdr_treshold, impute = None, global_impute = True, pickle_path = pickle_path)
    at_s, ce_s, hs_s = splitSpectronautBySpecies(spectronaut, truncated = False)
    return spectronaut, at_s, ce_s, hs_s

def compute_nDiff_exp(sample1, sample2, fc_treshold = 1.0, fdr_treshold = 0.05, compute_FDR = compute_logQ, col_diffExp = "red", col_not_diffExp = "blue"):
    """
    eg. 
    df = hs_t
    sample1 = df[s1]
    sample2 = df[s2]
    """
    log2fc_vals = get_point_estimate_log2fc(sample1, sample2)
    logFDR = compute_FDR(sample1, sample2, logFunc)
    logFDR = logFDR.squeeze()
    log2fc_vals = log2fc_vals.reindex(index = logFDR.index)
    df_volcano_format = volcano_df_format(log2fc_vals, logFDR, logFunc = logFunc, side = side, fc_treshold = fc_treshold, p_treshold =  fdr_treshold, col_diffExp = col_diffExp, col_not_diffExp = col_not_diffExp)
    
    n_diffExp = n_diffExp_df(df_volcano_format, col_diffExp)
    return n_diffExp

sample1 = hs_t["S02"]
sample2 = hs_t["S06"]

n = 100 #each p takes ~3s and q takes ~17s
fdrs = np.linspace(0,0.2, n+1)  
delta = fdrs[1]-fdrs[0]
fdrs = np.append(fdrs, np.linspace(0.2 + delta, 1, int(n/10 + 1)))

diffExp_array = np.array([])
count = 0
for fdr in fdrs:
    print(str(count) + " of " + str(len(fdrs)) + " bins.")
    fdr_treshold = fdr
    n_diffExp = compute_nDiff_exp(sample1, sample2, fc_treshold = fc_treshold, fdr_treshold = fdr_treshold, compute_FDR = compute_FDR, col_diffExp = col_diffExp, col_not_diffExp = col_not_diffExp )
    diffExp_array = np.append(diffExp_array, n_diffExp)
    count += 1


 n_diffExp = compute_nDiff_exp(sample1, sample2, fc_treshold = fc_treshold, fdr_treshold = fdr_treshold, compute_FDR = compute_FDR, col_diffExp = col_diffExp, col_not_diffExp = col_not_diffExp )


# ACTUALLY ONLY NEED df_volcano_format to compute, dont need to do all the time consuming stuff....
 


