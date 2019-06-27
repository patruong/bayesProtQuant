#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 18:02:28 2019

@author: ptruong
"""

import pandas as pd 
import numpy as np

#filename = "500-PSSS3-raw-reformatted_dropna_dropdup.csv"
#df = pd.read_csv(filename, sep ="\t")
    
def getRunId_from_reformattedRawPSSS3(df, noDecoy = True):
    """
    Function returns run_id from reformatted raw PSSS3 data.
    
    example:
        # if not parsed raw data
        df = parsePSSS3equDecoy(filename = '500-PSSS3-equ decoy_Report.xls') #SLOW
        
        # else
        df.to_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t", index = False)
        df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t")    
        run_id = getRunId_from_reformattedRawPSSS3(df)
    """

    #df = df.drop(['Cscore', 'nrStrippedSequences', 'qvalue'], axis = 1)
    if noDecoy == True:
        df.drop(['qvalue'], axis = 1)
    else:
        df.drop(['qvalue', 'decoy'], axis = 1)
    # Get unique run_id
    run_id = df.sort_values(by = "run_id").run_id.unique()
    return run_id

def getProteins_from_reformattedRawPSSS3(df, noDecoy = True):
    """
    Function returns run_id from reformatted raw PSSS3 data.
    
    example:
        # if not parsed raw data
        df = parsePSSS3equDecoy(filename = '500-PSSS3-equ decoy_Report.xls') #SLOW
        
        # else
        df.to_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t", index = False)
        df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t")    
        proteins = getProteins_from_reformattedRawPSSS3(df)

    """
    
    #df = df.drop(['Cscore', 'nrStrippedSequences', 'qvalue'], axis = 1)
    if noDecoy == True:
        df.drop(['qvalue'], axis = 1)
    else:
        df = df.drop(['qvalue', 'decoy'], axis = 1)
    # Get unique proteins
    proteins = df.sort_values(by = "run_id").protein.unique()
    return proteins


def getTriqlerOutputFormat_from_reformattedRawPSSS3(df, noDecoy = True):
    """
    Takes reformatted raw PSSS3 and gives triqler format output
    
    example:
         # if not parsed raw data
        df = parsePSSS3equDecoy(filename = '500-PSSS3-equ decoy_Report.xls') #SLOW
        
        # else
        df.to_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t", index = False)
        df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t")    
        df = getTriqlerOutputFormat_from_reformattedRawPSSS3(df)        
    """
    run_id = getRunId_from_reformattedRawPSSS3(df, noDecoy = noDecoy)
    proteins = getProteins_from_reformattedRawPSSS3(df, noDecoy = noDecoy)
    #df = df.drop(["Cscore", "nrStrippedSequences", "qvalue"], axis = 1)
    if noDecoy == True:
        df = df.drop(["qvalue"], axis = 1)
    else:
        df = df.drop(["qvalue", "decoy"], axis = 1)
    df = df.set_index("run_id")
    
    # Split up the different run_ids
    dfRunIds = [] #same order as run_id
    for i in range(len(run_id)):
        runId = run_id[i]
        df_runID = df[df.index == runId]
        df_runID = df_runID.set_index("protein")
        df_runID.columns = [runId]
        dfRunIds.append(df_runID)
        
    # Merge to triqler format
    df = pd.concat(dfRunIds, axis = 1)
    df = df.drop_duplicates()
    return df
    '''
    df_main = pd.DataFrame()
    for i in dfRunIds:
        df_main = df_main.join(i, how = 'outer')
    return df_main
    '''

def triqlerFormat2multiIndex(df):
    """
    Takes triqler outputted PSSS3 format and converts to multiIndex df
    """
        
    run_id = df.columns
    proteins = df.index
    
    # Split samples and runs
    samples = []
    runs = []
    for i in run_id:
        vals = i.split(":")[1].split("_")
        samples.append(vals[0])
        runs.append(vals[1])
    samples = np.unique(samples)
    runs = np.unique(runs)
    
    # Make header tuple for multilevel dataframe
    headerTuples = []
    for i in samples:
        for j in runs:
            headerTuples.append((i,j))
    
    df_values = df.values.T
    df_zipped = dict(zip(headerTuples, df_values))
    
    df_multiIdx = pd.DataFrame(df_zipped, index=proteins)
    df_multiIdx.columns.names = ('sample','run')
    return df_multiIdx
    
def multiIndex2TriqlerFormat(df_multiIdx):
    """
    Reformat multiIndex df from triqlerFormat2multiIndex(df) to triqler PSSS3 output format.
    """
    levels = df_multiIdx.columns.levels
    labels = df_multiIdx.columns.labels
    headers = []
    for i in range(len(labels[0])):
        header = levels[0][labels[0][i]] + ":" + levels[0][labels[0][i]] + "_" + levels[1][labels[1][i]]
        headers.append(header)
    proteins = df_multiIdx.index
    df_values = df_multiIdx.values
    df_reformatted = pd.DataFrame(df_values, index = proteins, columns = headers)
    return df_reformatted


def multiIdxFillna_mean(df_multiIdx):
    """
    Sample-wise fillna with mean. Input is multiIndex df 
    
    example:
        df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t")
        df_reform = getTriqlerOutputFormat_from_reformattedRawPSSS3(df)
        df_multiIdx = triqlerFormat2multiIndex(df_reform)
        df_imputed = multiIdxFillna_mean(df_multiIdx)
    """
    df_multiIdx = df_multiIdx.transpose()
    for i in df_multiIdx.index.levels[0]:
        df_multiIdx.loc[i].fillna(df_multiIdx.loc[i].mean(), inplace = True)
    df_multiIdx = df_multiIdx.transpose()
    return df_multiIdx

def multiIdxFillna_median(df_multiIdx):
    """
    Sample-wise fillna with median. Input is multiIndex df 
    
    example:
        df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t")
        df_reform = getTriqlerOutputFormat_from_reformattedRawPSSS3(df)
        df_multiIdx = triqlerFormat2multiIndex(df_reform)
        df_imputed = multiIdxFillna_median(df_multiIdx)
    """
    df_multiIdx = df_multiIdx.transpose()
    for i in df_multiIdx.index.levels[0]:
        df_multiIdx.loc[i].fillna(df_multiIdx.loc[i].median(), inplace = True)
    df_multiIdx = df_multiIdx.transpose()
    return df_multiIdx

def multiIdxFillna_min(df_multiIdx):
    """
    Sample-wise fillna with median. Input is multiIndex df 
    
    example:
        df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t")
        df_reform = getTriqlerOutputFormat_from_reformattedRawPSSS3(df)
        df_multiIdx = triqlerFormat2multiIndex(df_reform)
        df_imputed = multiIdxFillna_min(df_multiIdx)
    """
    df_multiIdx = df_multiIdx.transpose()
    for i in df_multiIdx.index.levels[0]:
        df_multiIdx.loc[i].fillna(df_multiIdx.loc[i].min(), inplace = True)
    df_multiIdx = df_multiIdx.transpose()
    return df_multiIdx

def multiIdxFillna_max(df_multiIdx):
    """
    Sample-wise fillna with median. Input is multiIndex df 
    
    example:
        df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t")
        df_reform = getTriqlerOutputFormat_from_reformattedRawPSSS3(df)
        df_multiIdx = triqlerFormat2multiIndex(df_reform)
        df_imputed = multiIdxFillna_max(df_multiIdx)
    """
    df_multiIdx = df_multiIdx.transpose()
    for i in df_multiIdx.index.levels[0]:
        df_multiIdx.loc[i].fillna(df_multiIdx.loc[i].max(), inplace = True)
    df_multiIdx = df_multiIdx.transpose()
    return df_multiIdx

"""
import time
start = time.time()
df = pd.read_csv("500-PSSS3-raw-reformatted_dropna_dropdup.csv", sep = "\t")
df_reform = getTriqlerOutputFormat_from_reformattedRawPSSS3(df)
df_multiIdx = triqlerFormat2multiIndex(df_reform)
df_multiIdx = multiIdxFillna_mean(df_multiIdx)
df_psss3 = multiIndex2TriqlerFormat(df_multiIdx)
end = time.time()
total = end-start
print(total)
"""