#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 19:27:09 2019

@author: ptruong
"""
import pandas as pd
import numpy as np

'''
def parsePSSS3equDecoy(filename = '500-PSSS3-equ decoy_Report.xls'):
    """
    File for parsing PSSS3 unfiltered with decoys and PG.Quantity column. (Very slow)
    """
    #df = pd.read_csv(filename, sep = "\t", usecols = [1,2,3,4,6,7,9])
    df = pd.read_csv(filename, sep = "\t", usecols = [1,2,3,6,7,9])
    #df = df.dropna() #dropna from raw file to make processing feasible <--
    
    # remove all shared proteins
    proteinId = df["PG.ProteinAccessions"].apply(lambda x: np.nan if ";" in x else x)

    ###################
    #df.loc[["PG.Organisms","PG.ProteinAccessions"]].apply(lambda x: np.nan if ";" in x else x)
    df.loc[:,["PG.Organisms", "PG.ProteinAccessions"]].apply(lambda x: np.nan if ";" in x else x)
    df[df["PG.Organisms"] != "Caenorhabditis elegans OX=6239"]["PG.ProteinAccessions"] = df[df["PG.Organisms"] != "Caenorhabditis elegans OX=6239"]["PG.ProteinAccessions"].apply(lambda x: np.nan if ";" in x else x)
    df[df["PG.Organisms"] != "Caenorhabditis elegans OX=6239"]["PG.ProteinAccessions"]
   
    def sharedProteinAdj(x):
        if "Caenorhabditis elegans OX=6239" in x:
            return x.replace("_Caenorhabditis elegans OX=6239", "")
        elif "Homo sapiens" in x:
            return x.replace("_Homo sapiens", "")
        elif "Arabidopsis thaliana" in x:
            return x.replace("_Arabidopsis thaliana", "")
        
    proteinId = (df["PG.ProteinAccessions"] + "_" + df["PG.Organisms"]).apply()

    ###################
    def sharedProtein_Adjustement(row):
        CE = "Caenorhabditis elegans OX=6239"
        if row["PG.Organisms"] ==  CE:
            return row["PG.ProteinAccessions"]
        elif row["PG.Organisms"] != CE:
            if ";" in row["PG.ProteinAccessions"]:
                return np.nan
            else:
                return row["PG.ProteinAccessions"]
                
    df.loc[:,["PG.Organisms", "PG.ProteinAccessions"]].apply(sharedProtein_Adjustement, axis = 1)
    ##############
    
    df["proteinId"] = proteinId
    df = df[pd.notnull(df["proteinId"])]
    
    
    reformatted_df = pd.DataFrame()
    reformatted_df["qvalue"] = df["PG.Qvalue"]
    df.drop(["PG.Qvalue"], inplace = True, axis = 1)
    reformatted_df["proteinQuantity"] = df["PG.Quantity"]
    df.drop(["PG.Quantity"], inplace = True, axis = 1)
    reformatted_df["decoy"] = df["EG.IsDecoy"]
    df.drop(["EG.IsDecoy"], inplace = True, axis = 1)
    specie = df["PG.Organisms"].apply(lambda x: "CAEEL" if x == "Caenorhabditis elegans OX=6239" 
             else "ARATH" if x == "Arabidopsis thaliana" else  "HUMAN" if x == 'Homo sapiens' else x)
    df.drop("PG.Organisms", inplace = True, axis = 1)
    #samples = df["R.FileName"].str.split("-", n=3, expand = True)[2].str.split("_", n=4, expand=True)[0]
    #runs = df["R.FileName"].str.split("-", n=3, expand = True)[2].str.split("_", n=4, expand=True)[2]
    samples = df["R.FileName"].str[-15:-12]
    runs = df["R.FileName"].str[-6:-3]
    df.drop("R.FileName", inplace = True, axis = 1)
    del df
    run_id = samples + ":" + samples + "_" + runs 
    run_id = run_id.apply(lambda x: "S03:S03_R05" if x == 'MHR:MHR_X01' else x) #special case for S03:S03_R05...
    protein = proteinId + "_" + specie
    
    #cScore = df["PG.Cscore"]
    #nrStripped = df["PG.NrOfStrippedSequencesIdentified"]
    
    #decoy = 
    
    reformatted_df["run_id"] = run_id
    reformatted_df["protein"] = protein
    #reformatted_df["Cscore"] = cScore
    #reformatted_df["nrStrippedSequences"] = nrStripped
    
    
    
    
    #reformatted_df = reformatted_df.dropna() # Bascially drop all shared proteins. (= NaN from previous row)
    
    # Drop duplicates
    reformatted_df = reformatted_df.drop_duplicates()

    # Sort by run_id
    # reformatted_df = df.sort_values(by = "run_id")

    return reformatted_df
'''

# input - df: a Dataframe, chunkSize: the chunk size
# output - a list of DataFrame
# purpose - splits the DataFrame into smaller of max size chunkSize (last is smaller)
def splitDataFrameIntoSmaller(df, chunkSize = 10000): 
    listOfDf = list()
    numberChunks = len(df) // chunkSize + 1
    for i in range(numberChunks):
        listOfDf.append(df[i*chunkSize:(i+1)*chunkSize])
    return listOfDf

def parsePSSS3equDecoy(filename = '500-PSSS3-equ decoy_Report.xls', index = False):
    if index == False:
        df = pd.read_csv(filename, sep = "\t", usecols = [1,2,3,6,7,9])
    else:
        df = pd.read_csv(filename, sep = "\t", usecols = [2,3,4,7,8,10])
    chunks = splitDataFrameIntoSmaller(df, chunkSize = 1000000)
    
    reformatted_chunks = []
    for i in range(len(chunks)):
        chunk = chunks[i]
        
        #proteinId = chunk["PG.ProteinAccessions"].apply(lambda x: np.nan if ";" in x else x) 
        
        def sharedProteinAdj(x):
            """
            Remove Shared Peptide, ignore CAEEL peptides. Basic version. 
            Probabily will need to adjust for HS och AT also, depending on database search.
            """
            if "Caenorhabditis elegans OX=6239" in x:
                return x.replace("_Caenorhabditis elegans OX=6239", "")
            elif "Homo sapiens" in x:
                if ";" in x:
                    return np.nan
                else:
                    return x.replace("_Homo sapiens", "")
            elif "Arabidopsis thaliana" in x:
                if ";" in x:
                    return np.nan
                else:
                    return x.replace("_Arabidopsis thaliana", "")
            
        proteinId = (chunk["PG.ProteinAccessions"] + "_" + chunk["PG.Organisms"]).apply(sharedProteinAdj)
        
        chunk["proteinId"] = proteinId
        chunk = chunk[pd.notnull(chunk["proteinId"])]
        
        
        reformatted_chunk = pd.DataFrame()
        reformatted_chunk["qvalue"] = chunk["PG.Qvalue"]
        chunk.drop(["PG.Qvalue"], inplace = True, axis = 1)
        reformatted_chunk["proteinQuantity"] = chunk["PG.Quantity"]
        chunk.drop(["PG.Quantity"], inplace = True, axis = 1)
        reformatted_chunk["decoy"] = chunk["EG.IsDecoy"]
        chunk.drop(["EG.IsDecoy"], inplace = True, axis = 1)
        specie = chunk["PG.Organisms"].apply(lambda x: "CAEEL" if x == "Caenorhabditis elegans OX=6239" 
                 else "ARATH" if x == "Arabidopsis thaliana" else  "HUMAN" if x == 'Homo sapiens' else x)
        chunk.drop("PG.Organisms", inplace = True, axis = 1)
        samples = chunk["R.FileName"].str[-15:-12]
        runs = chunk["R.FileName"].str[-6:-3]
        chunk.drop("R.FileName", inplace = True, axis = 1)
        del chunk
        run_id = samples + ":" + samples + "_" + runs 
        run_id = run_id.apply(lambda x: "S03:S03_R05" if x == 'MHR:MHR_X01' else x) #special case for S03:S03_R05...
        protein = proteinId + "_" + specie
        
        
        reformatted_chunk["run_id"] = run_id
        reformatted_chunk["protein"] = protein
        
        # Drop duplicates
        reformatted_chunk = reformatted_chunk.drop_duplicates()
        
        reformatted_chunks.append(reformatted_chunk)
        
    df = pd.concat(reformatted_chunks)
    df = df.drop_duplicates()
    return df
