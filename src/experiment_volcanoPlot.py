#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 17:15:18 2019

@author: ptruong
"""

from plotting import *
from proc import *
from func import *


#############
# Triqler  ##
#############

def get_triqler_volcano_df_pVals(triqlerFile, protein_id_fdr_treshold, specie, s1, s2, logFunc = np.log2, side = "two", fc_treshold = 1.0, fc_dev = 0.2, p_treshold = 0.05, col_diffExp = "red", col_not_diffExp = "blue", pickle_path = "../data/pickled/"):
    """
    s1,s2 = S01...S10
    """
    # Read data
    triqler = read_pickle_Triqler(triqlerFile, protein_id_fdr_treshold, pickle_path = pickle_path)
    at_t, ce_t, hs_t = splitTriqlerBySpecies(triqler, exponential = False, truncated = False)
    
    if specie == "HUMAN":
        df = hs_t
    elif specie == "ARATH":
        df = at_t
    elif specie == "CAEEL":
        df = ce_t
    #define samples
    sample1 = df[s1]
    sample2 = df[s2]
    log2fc_vals = get_point_estimate_log2fc(sample1, sample2)
    #log2P = compute_log2P(sample1, sample2)
    #log10P = compute_log10P(sample1, sample2)
    logP = compute_logP(sample1, sample2, logFunc)
    #df_volcano_format = volcano_df_format(log2fc_vals, log2P[0], side = side, fc_treshold = fc_treshold, p_treshold = p_treshold, col_diffExp = col_diffExp, col_not_diffExp = col_not_diffExp)
    df_volcano_format = volcano_df_format(log2fc_vals, logP[0], logFunc = logFunc, side = side, fc_treshold = fc_treshold, fc_dev = fc_dev, p_treshold = p_treshold, col_diffExp = col_diffExp, col_not_diffExp = col_not_diffExp)
    return df_volcano_format

def get_triqler_volcano_df_qVals(triqlerFile, protein_id_fdr_treshold, specie, s1, s2, logFunc = np.log2, side = "two", fc_treshold = 1.0, fc_dev = 0.2, p_treshold = 0.05, col_diffExp = "red", col_not_diffExp = "blue", pickle_path = "../data/pickled/"):
    """
    s1,s2 = S01...S10
    """
    # Read data
    triqler = read_pickle_Triqler(triqlerFile, protein_id_fdr_treshold, pickle_path = pickle_path)
    at_t, ce_t, hs_t = splitTriqlerBySpecies(triqler, exponential = False, truncated = False)
    
    if specie == "HUMAN":
        df = hs_t
    elif specie == "ARATH":
        df = at_t
    elif specie == "CAEEL":
        df = ce_t
    #define samples
    sample1 = df[s1]
    sample2 = df[s2]
    log2fc_vals = get_point_estimate_log2fc(sample1, sample2)
    logQ = compute_logQ(sample1, sample2, logFunc)
    log2fc_vals = log2fc_vals.reindex(index = logQ.index)
    df_volcano_format = volcano_df_format(log2fc_vals, logQ, logFunc = logFunc, side = side, fc_treshold = fc_treshold, fc_dev = fc_dev, p_treshold = p_treshold, col_diffExp = col_diffExp, col_not_diffExp = col_not_diffExp)
    return df_volcano_format


###############
# Spectronaut #
###############
    
def get_spectronaut_volcano_df_pVals(spectronautFile, protein_id_fdr_treshold, specie, s1, s2, logFunc = np.log2, side = "two", fc_treshold = 1.0, fc_dev = 0.2, p_treshold = 0.05, col_diffExp = "red", col_not_diffExp = "blue", pickle_path = "../data/pickled/"):
    """
    s1,s2 = S01...S10
    """
    spectronaut =  read_pickle_Spectronaut(spectronautFile, protein_id_fdr_treshold, impute = None, global_impute = True, pickle_path = pickle_path)
    at_s, ce_s, hs_s = splitSpectronautBySpecies(spectronaut, truncated = False)

    if specie == "HUMAN":
        df = hs_s
    elif specie == "ARATH":
        df = at_s
    elif specie == "CAEEL":
        df = ce_s        
    #define samples
    sample1 = df[s1]
    sample2 = df[s2]    
    log2fc_vals = get_point_estimate_log2fc(sample1, sample2)
    #log2P = compute_log2P(sample1, sample2)    
    #log10P = compute_log10P(sample1, sample2)
    logP = compute_logP(sample1, sample2, logFunc)    
    #df_volcano_format = volcano_df_format(log2fc_vals, log2P[0], side = side, fc_treshold = fc_treshold, p_treshold = p_treshold, col_diffExp = col_diffExp, col_not_diffExp = col_not_diffExp)
    df_volcano_format = volcano_df_format(log2fc_vals, logP[0], logFunc = logFunc, side = side, fc_treshold = fc_treshold, fc_dev = fc_dev, p_treshold = p_treshold, col_diffExp = col_diffExp, col_not_diffExp = col_not_diffExp)
    return df_volcano_format

def get_spectronaut_volcano_df_qVals(spectronautFile, protein_id_fdr_treshold, specie, s1, s2, logFunc = np.log2, side = "two", fc_treshold = 1.0, fc_dev = 0.2, p_treshold = 0.05, col_diffExp = "red", col_not_diffExp = "blue", pickle_path = "../data/pickled/"):
    """
    s1,s2 = S01...S10
    """
    spectronaut =  read_pickle_Spectronaut(spectronautFile, protein_id_fdr_treshold, impute = None, global_impute = True, pickle_path = pickle_path)
    at_s, ce_s, hs_s = splitSpectronautBySpecies(spectronaut, truncated = False)

    if specie == "HUMAN":
        df = hs_s
    elif specie == "ARATH":
        df = at_s
    elif specie == "CAEEL":
        df = ce_s        
    #define samples
    sample1 = df[s1]
    sample2 = df[s2]    
    log2fc_vals = get_point_estimate_log2fc(sample1, sample2)
    logQ = compute_logQ(sample1, sample2, logFunc)
    log2fc_vals = log2fc_vals.reindex(index = logQ.index)
    df_volcano_format = volcano_df_format(log2fc_vals, logQ, logFunc = logFunc, side = side, fc_treshold = fc_treshold, fc_dev = fc_dev, p_treshold = p_treshold, col_diffExp = col_diffExp, col_not_diffExp = col_not_diffExp)
    return df_volcano_format



if __name__ == "__main__":
    # plot params
    plotWidth = 1000
    plotHeight = 1000
    x_label = "log2(FC)"
    y_label = "-log10(q-value)"
    
    # params    
    triqlerFile = "tmp"
    spectronautFile = "tmp"
    protein_id_fdr_treshold = 0.01
    fc_treshold = 0.95
    fc_dev = fc_treshold/2
    p_treshold = 0.05
    side = "fc"
    specie = "HUMAN"
    s1 = "S02"
    s2 = "S06"
    logFunc = np.log10 #Adjust y-label to match this!
    logP_treshold = -logFunc(p_treshold)

    triqler_diffExp_col = "firebrick"
    triqler_not_diffExp_col = "gold" 
    spectronaut_diffExp_col = "green"
    spectronaut_not_diffExp_col = "steelblue"
#    
#    triqler_diffExp_col = "steelblue"
#    triqler_not_diffExp_col = "steelblue" 
#    spectronaut_diffExp_col = "green"
#    spectronaut_not_diffExp_col = "green"
    
    title_triqler = "Triqler| " + "proteid_id_fdr: " + str(protein_id_fdr_treshold) + "| fc_treshold: " + str(fc_treshold) + "| p_treshold " + str(p_treshold) + "| specie: " + specie + " " + s1 + ":" + s2
    title_spec = "Spectronaut| " + "proteid_id_fdr: " + str(protein_id_fdr_treshold) + "| fc_treshold: " + str(fc_treshold) + "| p_treshold " + str(p_treshold) + "| specie: " + specie + " " + s1 + ":" + s2
    title_overlay = "Overlay| " + "proteid_id_fdr: " + str(protein_id_fdr_treshold) + "| fc_treshold: " + str(fc_treshold) + "| p_treshold " + str(p_treshold) + "| specie: " + specie + " " + s1 + ":" + s2
    
    outputFile_triqler = "volcano_triqler_"+specie+"_"+s1+"_"+s2 + ".html"
    outputFile_spec = "volcano_spectronaut_"+specie+"_"+s1+"_"+s2+ ".html"
    outputFile_overlay = "volcano_overlay_"+specie+"_"+s1+"_"+s2+ ".html"
    

    
    df_triqler = get_triqler_volcano_df_pVals(triqlerFile, protein_id_fdr_treshold, specie, s1, s2, 
                                              logFunc = logFunc, side = side, fc_treshold = fc_treshold, fc_dev = fc_dev,
                                              p_treshold = p_treshold, col_diffExp = triqler_diffExp_col,
                                              col_not_diffExp = triqler_not_diffExp_col)
#    df_triqler = get_triqler_volcano_df_qVals(triqlerFile, protein_id_fdr_treshold, specie, s1, s2, 
#                                              logFunc = logFunc, side = side, fc_treshold = fc_treshold, fc_dev = fc_dev,
#                                              p_treshold = p_treshold, col_diffExp = triqler_diffExp_col,
#                                              col_not_diffExp = triqler_not_diffExp_col)
   
    volcanoPlot(df_triqler,
                 vertical_line_r = fc_treshold, vertical_line_l = -fc_treshold, horizontal_line = logP_treshold,
                 plotWidth = plotWidth, plotHeight = plotHeight, title = title_triqler,
                x_label = x_label, y_label = y_label, outputFile = outputFile_triqler)

    

    
    df_spec = get_spectronaut_volcano_df_pVals(spectronautFile, protein_id_fdr_treshold, specie, s1, s2,
                                               logFunc = logFunc, side = side, fc_treshold = fc_treshold, fc_dev = fc_dev,
                                               p_treshold = p_treshold, col_diffExp = spectronaut_diffExp_col,
                                               col_not_diffExp = spectronaut_not_diffExp_col)
#    df_spec = get_spectronaut_volcano_df_qVals(spectronautFile, protein_id_fdr_treshold, specie, s1, s2,
#                                               logFunc = logFunc, side = side, fc_treshold = fc_treshold, fc_dev = fc_dev
#                                               p_treshold = p_treshold, col_diffExp = spectronaut_diffExp_col,
#                                               col_not_diffExp = spectronaut_not_diffExp_col)
#    
    volcanoPlot(df_spec, 
                vertical_line_r = fc_treshold, vertical_line_l = -fc_treshold, horizontal_line = logP_treshold,
                plotWidth = plotWidth, plotHeight = plotHeight, title = title_spec,
                x_label = x_label, y_label = y_label, outputFile = outputFile_spec)
    
    volcanoPlot_spectronaut_triqler(df_triqler, df_spec, 
                                    vertical_line_r = fc_treshold, vertical_line_l = -fc_treshold, horizontal_line = logP_treshold,
                                    plot_width = plotWidth, plot_height = plotHeight, title = title_overlay,
                                    x_label = x_label, y_label = y_label, outputFile = outputFile_overlay)
    
    
    df_triqler.