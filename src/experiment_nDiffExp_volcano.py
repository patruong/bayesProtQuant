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

sample1 = hs_t["S02"]
sample2 = hs_t["S06"]
 
log2fc_vals = get_point_estimate_log2fc(sample1, sample2)
logFDR = compute_FDR(sample1, sample2, logFunc)
logFDR = logFDR.squeeze()
log2fc_vals = log2fc_vals.reindex(index = logFDR.index)

n = 10000 #each p takes ~3s and q takes ~17s
fdrs = np.linspace(0,1, n+1)
count = 0
for fdr in fdrs:  
    print(str(count) + " of " + str(len(fdrs)) + " bins.")
    df_volcano_format = volcano_df_format(log2fc_vals, logFDR, logFunc = logFunc, side = side, fc_treshold = fc_treshold, p_treshold =  fdr, col_diffExp = col_diffExp, col_not_diffExp = col_not_diffExp)
    n_diffExp = n_diffExp_df(df_volcano_format, col_diffExp)
    count += 1

