#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 16:58:14 2019

@author: ptruong
"""


from plotting import *
from proc import *
from func import *
from log2FC import * 

#############
# Triqler  ##
#############

def readinTriqler(triqlerFile, protein_id_fdr_treshold, pickle_path = "../data/pickled/"):
    triqler = read_pickle_Triqler(triqlerFile, protein_id_fdr_treshold, pickle_path = pickle_path)
    at_t, ce_t, hs_t = splitTriqlerBySpecies(triqler, exponential = False, truncated = False)
    return at_t, ce_t, hs_t

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

def df_triqler_volcano_qVals(df, s1, s2, logFunc = np.log2, side = "two", fc_treshold = 1.0, fc_dev = 0.2, p_treshold = 0.05, col_diffExp = "red", col_not_diffExp = "blue"):
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
    
def readinSpectronaut(spectronautFile, protein_id_fdr_treshold, pickle_path = "../data/pickled/"):    
    spectronaut =  read_pickle_Spectronaut(spectronautFile, protein_id_fdr_treshold, impute = None, global_impute = True, pickle_path = pickle_path)
    at_s, ce_s, hs_s = splitSpectronautBySpecies(spectronaut, truncated = False)
    return at_s, ce_s, hs_s

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

def df_spectronaut_volcano_qVals(df, s1, s2, logFunc = np.log2, side = "two", fc_treshold = 1.0, fc_dev = 0.2, p_treshold = 0.05, col_diffExp = "red", col_not_diffExp = "blue"):
    sample1 = df[s1]
    sample2 = df[s2]    
    log2fc_vals = get_point_estimate_log2fc(sample1, sample2)
    logQ = compute_logQ(sample1, sample2, logFunc)
    log2fc_vals = log2fc_vals.reindex(index = logQ.index)
    df_volcano_format = volcano_df_format(log2fc_vals, logQ, logFunc = logFunc, side = side, fc_treshold = fc_treshold, fc_dev = fc_dev, p_treshold = p_treshold, col_diffExp = col_diffExp, col_not_diffExp = col_not_diffExp)
    return df_volcano_format

###############
# computation #
###############
    
def n_diffExp(df):
    """
    df - df from get_triqler_volcano_df_pVals() or get_spectronaut_volcano_df_pVals or qVals equivalent.
    """
    n = len(df[df.diffExp == True])
    return n


def diffExpRatio(df, includeNaN = False):
    """
    get the ratio of differential expression.
    """
    diffExp = n_diffExp(df)
    if includeNaN == False:
        ratio = diffExp / df.diffExp.count()
    else:
        ratio = diffExp / len(df)
    return ratio

def getComparisons(type = str):
    if type == str:
        samples = ["S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10"]
    elif type == int:
        samples = [i for i in range(1,11)]
    compList = []
    for i in samples:
        for j in samples[1:]:
            if ((i,j) in compList or (j,i) in compList or (i,j) == (i,i)):
                pass
            else:
                compList.append((i,j))
    return compList

def getFCs():
    comparisons = getComparisons(type = int)    
    comparisonKeys = getComparisons(type = str)
    FCs = {}
    for i in range(len(comparisons)):
        s1,s2 = comparisons[i]
        s1Key, s2Key = comparisonKeys[i]
        key = s1Key + ":" + s2Key
        FCs.update({key : get_fc_between_samples(s1,s2)})
    return FCs

if __name__ == "__main__":
    import time
    start = time.time()
    # plot params
    plotWidth = 1000
    plotHeight = 1000
    x_label = "log2(FC)"
    y_label = "-log10(q-value)"
    
    # params    
    triqlerFile = "tmp"
    spectronautFile = "tmp"
    protein_id_fdr_treshold = 0.01
#    
#    specie = "CAEEL"
#    s1 = "S02"
#    s2 = "S06"
    logFunc = np.log10 #Adjust y-label to match this!
    
#    triqler_diffExp_col = "firebrick"
#    triqler_not_diffExp_col = "gold" 
#    spectronaut_diffExp_col = "green"
#    spectronaut_not_diffExp_col = "steelblue"
#    
    triqler_diffExp_col = "steelblue"
    triqler_not_diffExp_col = "steelblue" 
    spectronaut_diffExp_col = "green"
    spectronaut_not_diffExp_col = "green"
    
    at_t, ce_t, hs_t = readinTriqler(triqlerFile, protein_id_fdr_treshold, pickle_path = "../data/pickled/")
    at_s, ce_s, hs_s = readinSpectronaut(spectronautFile, protein_id_fdr_treshold, pickle_path = "../data/pickled/")
    
    #fc_treshold = -4.00
    #fc_dev = fc_treshold/2
    p_treshold = 0.05
    logP_treshold = -logFunc(p_treshold)
    side = "fc"
#    vertical_line_r = fc_treshold + fc_dev
#    vertical_line_l = fc_treshold - fc_dev
    
#    title_triqler = "Triqler| " + "proteid_id_fdr: " + str(protein_id_fdr_treshold) + "| fc_treshold: " + str(fc_treshold) + "| p_treshold " + str(p_treshold) + "| specie: " + specie + " " + s1 + ":" + s2
#    title_spec = "Spectronaut| " + "proteid_id_fdr: " + str(protein_id_fdr_treshold) + "| fc_treshold: " + str(fc_treshold) + "| p_treshold " + str(p_treshold) + "| specie: " + specie + " " + s1 + ":" + s2
#    title_overlay = "Overlay| " + "proteid_id_fdr: " + str(protein_id_fdr_treshold) + "| fc_treshold: " + str(fc_treshold) + "| p_treshold " + str(p_treshold) + "| specie: " + specie + " " + s1 + ":" + s2
#    
#    outputFile_triqler = "volcano_triqler_"+specie+"_"+s1+"_"+s2 + ".html"
#    outputFile_spec = "volcano_spectronaut_"+specie+"_"+s1+"_"+s2+ ".html"
#    outputFile_overlay = "volcano_overlay_"+specie+"_"+s1+"_"+s2+ ".html"
    
    
    diffExp_dict = {}
    FCs = getFCs()
    species = ["ARATH", "CAEEL", "HUMAN"]
    totalCount = len(FCs)*len(species)
    count = 0
    for key in FCs.keys():
        for specie in species:
            print("Progress " + str(count) + "/" + str(totalCount) + ", Time elapsed: " , end = "")
            print(time.time()-start)
            s1, s2 = key.split(":")
            fc_treshold = FCs[key][specie]
            if fc_treshold > 9:
                fc_treshold = 9 #set a limit to log2FC(x) = np.log2(0.5) - np.log2(0.001)
            elif fc_treshold < -9:
                fc_treshold = -9
            fc_dev = fc_treshold/2
            side = "two"
            vertical_line_r = fc_treshold 
            vertical_line_l = -fc_treshold
            if specie == "HUMAN":
                df_t = hs_t
                df_s = hs_s
            elif specie == "ARATH":
                df_t = at_t
                df_s = at_s
                #side = "two" # ARATH trueFC = 0.0
                fc_treshold = 0.4 # Hard-coded treshold for ARATH
            elif specie == "CAEEL":
                df_t = ce_t
                df_s = ce_s
            df_t = df_triqler_volcano_qVals(df_t, s1, s2, logFunc = logFunc, side = side, 
                                          fc_treshold = fc_treshold, fc_dev = fc_dev,
                                          p_treshold = p_treshold, col_diffExp = triqler_diffExp_col,
                                          col_not_diffExp = triqler_not_diffExp_col)
            df_s = df_spectronaut_volcano_qVals(df_s, s1, s2, logFunc = logFunc, side = side,
                                              fc_treshold = fc_treshold, fc_dev = fc_dev,
                                              p_treshold = p_treshold, col_diffExp = spectronaut_diffExp_col,
                                              col_not_diffExp = spectronaut_not_diffExp_col)                
            n_diffExp_proteins_t = n_diffExp(df_t)
            diffExp_ratio_t = diffExpRatio(df_t)
            n_diffExp_proteins_s = n_diffExp(df_s)
            diffExp_ratio_s = diffExpRatio(df_s)
            diffExp_dict.update({"triqler" + "_" + key + "_" + specie : (n_diffExp_proteins_t, diffExp_ratio_t)})
            diffExp_dict.update({"spectronaut" + "_" + key + "_" + specie : (n_diffExp_proteins_s, diffExp_ratio_s)})

            title_triqler = ("Triqler" + 
                             "| " + specie + " " + s1 + ":" + s2 +
                             "| protein_fdr: " + str(protein_id_fdr_treshold) + 
                             "| fc: " + str(fc_treshold) + 
                             "| q_tresh: " + str(p_treshold) + 
                             "| DE: " + str(n_diffExp_proteins_t) + 
                             "| DE-ratio: " + str(diffExp_ratio_t)
                             )
            title_spec = ("Spectronaut" +
                          "| " + specie + " " + s1 + ":" + s2 + 
                          "| protein_fdr: " + str(protein_id_fdr_treshold) + 
                          "| fc: " + str(fc_treshold) + 
                          "| q_treshold: " + str(p_treshold) +
                          "| DE: " + str(n_diffExp_proteins_s) + 
                          "| DE ratio: " + str(diffExp_ratio_s)                          
                          )
            title_overlay = ("Overlay" +
                             "| " + specie + " " + s1 + ":" + s2 +
                             "| protein_fdr: " + str(protein_id_fdr_treshold) + 
                             "| fc: " + str(fc_treshold) + 
                             "| q_treshold: " + str(p_treshold) +
                             "| DE(t): " + str(n_diffExp_proteins_t) + 
                             "| DEr(t): " + str(diffExp_ratio_t) +
                             "| DE(s): " + str(n_diffExp_proteins_s) +
                             "| DEr(s): " + str(diffExp_ratio_s)
                             )                
     
            outputFile_triqler = "volcano_triqler_"+specie+"_"+s1+"_"+s2 + ".html"
            outputFile_spec = "volcano_spectronaut_"+specie+"_"+s1+"_"+s2+ ".html"
            outputFile_overlay = "volcano_overlay_"+specie+"_"+s1+"_"+s2+ ".html"
            
            volcanoPlot(df_t,
                        vertical_line_r = vertical_line_r, vertical_line_l = vertical_line_l, horizontal_line = logP_treshold,
                        plotWidth = plotWidth, plotHeight = plotHeight, title = title_triqler,
                        x_label = x_label, y_label = y_label, outputFile = outputFile_triqler)
            volcanoPlot(df_s, 
                        vertical_line_r =  vertical_line_r, vertical_line_l = vertical_line_l, horizontal_line = logP_treshold,
                        plotWidth = plotWidth, plotHeight = plotHeight, title = title_spec,
                        x_label = x_label, y_label = y_label, outputFile = outputFile_spec)
            volcanoPlot_spectronaut_triqler(df_s, df_t, 
                                    vertical_line_r =  vertical_line_r, vertical_line_l = vertical_line_l, horizontal_line = logP_treshold,
                                    plot_width = plotWidth, plot_height = plotHeight, title = title_overlay,
                                    x_label = x_label, y_label = y_label, outputFile = outputFile_overlay)
            count+=1
    # toCSV with pandas
    results = pd.DataFrame.from_dict(diffExp_dict, orient = "index", columns = ["diffExp", "diffExpRatio"])
    results.to_csv("volcanoResults.csv", sep = "\t")
    CSV = "\n".join([k+":"+ "("+",".join(tuple(str(x) for x in v))+")" for k,v in diffExp_dict.items()])
    # to dict formatted .txt file
    f = open("volcanoResultsDict.txt", "w")
    f.write(CSV)
    f.close()
    end = time.time()
    print("Process finished!")
    print("Elapsed time: ", end ="")
    print(end-start)
    
    # Find all fcs...
    # Find all fc_ratios...
    
    
    

#    df_triqler = get_triqler_volcano_df_qVals(triqlerFile, protein_id_fdr_treshold, specie, s1, s2, 
#                                              logFunc = logFunc, side = side, fc_treshold = fc_treshold, fc_dev = fc_dev,
#                                              p_treshold = p_treshold, col_diffExp = triqler_diffExp_col,
#                                              col_not_diffExp = triqler_not_diffExp_col)
#
#    df_spec = get_spectronaut_volcano_df_qVals(spectronautFile, protein_id_fdr_treshold, specie, s1, s2,
#                                               logFunc = logFunc, side = side, fc_treshold = fc_treshold, fc_dev = fc_dev,
#                                               p_treshold = p_treshold, col_diffExp = spectronaut_diffExp_col,
#                                               col_not_diffExp = spectronaut_not_diffExp_col)
#    
#
#    n_diffExp(df_triqler)
#    n_diffExp(df_spec)