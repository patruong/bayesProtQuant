#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 13:17:05 2019

@author: ptruong

NOTE: Function naming is a bit off... getAbundances and getAllAbundances acutally gets fractions,
they are used in getFractions() function to get a neat table of fractions.

File should be renamed TriqlerProcessor.py

"""

from triqlerParser import * 
import matplotlib.pyplot as plt
import matplotlib

def formatTriqlerInput(df, startRun = "S01:S01_R01", endRun = "S10:S10_R05"):
    """
    Format the triqler <x>vs<y> output to protein, species and abundances in condition:run
    Input:
        df - dataframe of <x>vs<y>.tsv output from triqler. (dataframe parsing from table2Df function.)
    """
    formatted_df = df.loc[:,"S01:S01_R01":"S10:S10_R05"] 
    formatted_df = formatted_df.apply(pd.to_numeric)
    formatted_df.insert(0, 'protein', df["protein"])
    species = formatted_df["protein"].str[-5:]
    formatted_df.insert(1, "species", species)
    return formatted_df


def getRuns(df, startRun = "S01:S01_R01", endRun = "S10:S10_R05"):
    """
    Get a series containing all the runs from pandas dataframe of <x>vs<y>.tsv. 
    
    example:
        triqler = table2Df("PSSS3_triqler_output_proteins.1vs2.tsv", peptideSeperator = True)
        runs = getRuns(triqler)
    """
    formatted_df = df.loc[:,startRun:endRun] # get the data for the conditions and runs
    runs = formatted_df.columns # get the conditions and runs
    runs = pd.DataFrame(runs, columns = ["runs"])
    return runs

def getConditionRunId(df, run = 'R01', startRun = "S01:S01_R01", endRun = "S10:S10_R05"):
    """
    Input:
        df - dataframe of <x>vs<y>.tsv output from triqler. (dataframe parsing from table2Df function.)
        run - run number as string. 
        startRun - starting column with run (as string)
        endRun - ending column with run (as string)
        
    example:
        triqler = table2Df("PSSS3_triqler_output_proteins.1vs3.tsv", peptideSeperator = True)
        runs = getRuns(triqler, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
        conditionRuns = getConditionRuns(runs)
        unique_run = conditionRuns.run.unique()
        run_id = getConditionRunId(triqler, unique_run[0])
    """
    run = run
    startRun = startRun
    endRun = endRun
    
    runs = getRuns(df, startRun = startRun, endRun = endRun)
    conditionRuns = getConditionRuns(runs)
    mask = conditionRuns[conditionRuns.run == run]
    run_id = (mask.condition + ":" + mask.condition + "_" + mask.run).values
    return run_id


def addSpecieProteinColumns(df, startRun = "S01:S01_R01", endRun = "S10:S10_R05"):
    """
    Split the protein string into species and proteins seperately and adds it into dataframe
    Inputs:
        df - dataframe of <x>vs<y>.tsv output from triqler. (dataframe parsing from table2Df function.)
    
    example:
        triqler = table2Df("PSSS3_triqler_output_proteins.1vs2.tsv", peptideSeperator = True)
        
    """
    runs = getRuns(df, startRun = startRun, endRun = endRun)
    added_df = df[["protein"]+runs.values.transpose().tolist()[0]]
    species = added_df["protein"].str[-5:]
    proteins = added_df["protein"].str[:-6]
    
    added_df["id"] = added_df["protein"]
    added_df["specie"] = species
    added_df["protein"] = proteins
    cols = added_df.columns[-2:].tolist() + added_df.columns[:-2].tolist() 
    added_df = added_df[cols]
    return added_df

def getAbundances(df, run):
    """
    df - PSSS3 dataframe input with species column.
    run - the name of the run columns.
    """
    totalProteinAbundance = df[run].sum()
    fraction = df[run]/totalProteinAbundance
    sampleRun = pd.DataFrame(data = [df["species"], fraction]).transpose()
    uniqueSpecies = sampleRun["species"].unique()
    #print(uniqueSpecies)
    fractionSpecies = []
    for i in uniqueSpecies:
        sumOfSpecies = sampleRun.loc[sampleRun["species"] == i][run].sum()
        fractionSpecies.append(sumOfSpecies)
    fractions = dict(zip(uniqueSpecies, fractionSpecies))
    return fractions

def getConditionRuns(runs):
    """
    Takes the output from getRuns and seperated conditions and runs.
    """
    conditionRuns = runs["runs"].str.split(":", n = 1, expand = True)
    conditionRuns = conditionRuns[1].str.split("_", n = 1, expand = True)
    conditionRuns = pd.DataFrame(conditionRuns.values, columns = ["condition", "run"])
    return conditionRuns

def getAllAbundances(df, startRun = "S01:S01_R01", endRun = "S10:S10_R05"):
    """
    Gives all abundances for all condition:runs from triqler <x>vs<y>.tsv file
    Input:
        df - dataframe of <x>vs<y>.tsv output from triqler. (dataframe parsing from table2Df function.)
    """
    runs = getRuns(df, startRun = startRun, endRun = endRun)
    formatted_df = formatTriqlerInput(df)
    runs = np.array(runs.transpose())[0]
    abundances = []
    for i in range(len(runs)):
        abundance_i = getAbundances(formatted_df, runs[i])
        abundances.append(abundance_i)
    abundances = pd.DataFrame([runs, abundances], index = ["condition:run", "fractions"]).transpose()
    return abundances

def getFractions(df):
    """
    Gives a df with condition:run, condition, run and fractions of the species.
    Input:
        df - dataframe of <x>vs<y>.tsv output from triqler. (dataframe parsing from table2Df function.)
    
    example: 
        triqler = table2Df("PSSS3_triqler_output_proteins.1vs2.tsv", peptideSeperator = True)
        fractions = getFractions(triqler)
    """
    abundances = getAllAbundances(df)
    runs = getRuns(df)
    conditionRuns = getConditionRuns(runs)
    run = conditionRuns["run"]
    condition = conditionRuns["condition"]
    abundances.insert(1, "condition", condition)
    abundances.insert(2, "run", run)
    return abundances

def dictMerge(dictionaryList):
    """
    Function for merging dictionaries.
    Input:
        dictionaryList - list of dictionaries
    """
    d = {}
    for k in dictionaryList[0]: #for k in d.iterkeys()
        d[k] = tuple(d[k] for d in dictionaryList)
    return d

def meanVarFromConditions(fractions, condition):
    """
    Gives mean and standard deviation from one condition from fractions table (from getFractions()).
    Input:
        fractions - df from getFractions()
        condition - specific condition name from fractions["condition"]
    Output:
        dictionary of means and std of the different species.
    example:
        fractions = getFractions(triqler)
        samples = fractions["condition"].unique()
        mean, std = meanVarFromConditions(fractions, samples[0])
    """
    conditionFractions = dictMerge(list(fractions[fractions["condition"] == condition]["fractions"]))
    conditionFractions = pd.DataFrame(conditionFractions)
    mean = conditionFractions.mean().to_dict()
    std = conditionFractions.std().to_dict()
    return mean, std

def meanVarFractions(fractions):
    """
    Gives means and variance for all conditions in fractions dataframe (from getFractions()).
    Input:
        fractions - df from getFractions()
    """
    means = []
    stds = []
    samples = fractions["condition"].unique()
    for i in range(len(samples)):
        mean,std = meanVarFromConditions(fractions, samples[i])
        means.append(mean)
        stds.append(std)
        
    means = dictMerge(means)    
    stds = dictMerge(stds)
    return means, stds

def getAllConditionFractions(fractions):
    """
    Get all the fractional intensities from all the conditions.
    Input:
        fractions - fractions - df from getFractions()
    example:
        fractions = getFractions(triqler)
        conditionsData = getAllConditionFractions(fractions)
    """
    conditions = fractions["condition"].unique()
    conditionsData = []
    for i in range(len(samples)):
        conditionFractions = dictMerge(list(fractions[fractions["condition"] == conditions[i]]["fractions"]))
        #conditionData_i = {conditions[i]:conditionFractions}        
        #conditionsData.append(conditionData_i)
    return conditionsData


def conditionsData2Df(conditionsData):
    """
    Takes conditionsData and reformats to dataframe.
    Input:
        conditionsData - output from getAllConditionsFractions().
    Example:
        triqler = table2Df("PSSS3_triqler_output_proteins.1vs2.tsv", peptideSeperator = True)
        fractions = getFractions(triqler)
        conditionsData = getAllConditionFractions(fractions) #data for plots containing all runs.
    """
    key_list = []
    species_list = []
    val_list = []
    for i in range(len(conditionsData)):
        key = list(conditionsData[i].keys())[0]
        species = list(conditionsData[i][key])
        for j in species:
            val = conditionsData[i][key][j]
            for k in val:
                key_list.append(key)
                species_list.append(j)
                val_list.append(float(k))
    df = pd.DataFrame(np.array([np.array(key_list), np.array(species_list), 
                                np.array(val_list)]).T, columns = ["condition","specie","fraction"])
    df["fraction"] = df["fraction"].apply(pd.to_numeric)
    return df


def fractions2Df(fractions):
    """
    OBS! SOMETHING IS BUGGING IN THIS ONE CHEECK THE LIST Appending thingy!
    
    Takes fractions from getFractions() and seperates the dictionary in "fractions" column into seperate columns.
    Input:
        fractions - fractions - df from getFractions()
    Example:
        triqler = table2Df("PSSS3_triqler_output_proteins.1vs2.tsv", peptideSeperator = True)
        fractions = getFractions(triqler)
        df = fractions2Df(fractions)
    """
    condition_list = []
    run_list = []
    species_list = []
    fraction_list = []
    for i in range(len(fractions)):
        condition = fractions.loc[i]["condition"]
        run = fractions.loc[i]["run"]
        fraction_i = fractions.loc[i]["fractions"]
        keys = list(fraction_i.keys())
        for j in keys:
            val = fraction_i[j]
            condition_list.append(condition)
            run_list.append(run)
            species_list.append(j)
            fraction_list.append(float(val))
    df = pd.DataFrame(np.array([np.array(condition_list), np.array(run_list),
                                np.array(species_list), np.array(val_list)]).T,
        columns = ["condition","run","specie", "fraction"])
    df["fraction"] = df["fraction"].apply(pd.to_numeric)
    return df

def getMixtures():
    """
    Gives the mixture fractions
    """
    conditions = np.array(['S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09',
       'S10'])
    a_mix = np.array([0.5] * 10)
    #a_mix[0] += 0.0001
    c_mix = np.array([0.5, 0.25, 0.125, 0.0625, 0.031, 0.0155, 0.008, 0.004, 0.002, 0.0])
    h_mix = np.array([0.0, 0.25, 0.375, 0.4375, 0.469, 0.4845, 0.492, 0.496, 0.498, 0.5])    
    mixtures = pd.DataFrame([a_mix, c_mix, h_mix], columns = conditions, index = ["HUMAN","CAEEL","ARATH"])
    return mixtures

##########################################
# Cyano PSSS3 specific functions (Vital) #
##########################################
def renameTriqler2Vital_nonDetailed(df, keepRuns = False):
    """
    Rename triqler output to become similar to Vitals format for PSSS3 analysis
    
    Example:
        triqler = table2Df("PSSS3_triqler_output_proteins.1vs2.tsv", peptideSeperator = True)
        df_triqler = addSpecieProteinColumns(triqler, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
        at_t = df_triqler[df_triqler["specie"] == "ARATH"]
        at_t = renameTriqler(at_t)
    """
    if keepRuns == False:
        repl_cols_s = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
        cols_allRun = list(np.array([[i,i,i,i,i] for i in repl_cols_s]).flatten())
        old_cols = df.drop(["id", "specie", "protein"], axis = 1).columns
        renameColsMap = dict(zip(old_cols, cols_allRun))
        df = df.drop(["id", "specie"], axis = 1).set_index(["protein"]).rename(columns = renameColsMap)
    else:
        df = df.drop(["id", "specie"], axis = 1).set_index(["protein"])
    df = df.convert_objects(convert_numeric=True)
    return df

def renameTriqler2Vital(df, keepRuns = False):
    """
    Rename triqler output to become similar to Vitals format for PSSS3 analysis
    
    Example:
        triqler = table2Df("PSSS3_triqler_output_proteins.1vs2.tsv", peptideSeperator = True)
        df_triqler = addSpecieProteinColumns(triqler, startRun = "S01:S01_R01", endRun = "S10:S10_R05")
        at_t = df_triqler[df_triqler["specie"] == "ARATH"]
        at_t = renameTriqler(at_t)
    """
#    old_cols = ['S01:S01_R01', 'S01:S01_R02', 'S01:S01_R03',
#           'S01:S01_R04', 'S01:S01_R05', 'S02:S02_R01', 'S02:S02_R02',
#           'S02:S02_R03', 'S02:S02_R04', 'S02:S02_R05', 'S03:S03_R01',
#           'S03:S03_R02', 'S03:S03_R03', 'S03:S03_R04', 'S03:S03_R05',
#           'S04:S04_R01', 'S04:S04_R02', 'S04:S04_R03', 'S04:S04_R04',
#           'S04:S04_R05', 'S05:S05_R01', 'S05:S05_R02', 'S05:S05_R03',
#           'S05:S05_R04', 'S05:S05_R05', 'S06:S06_R01', 'S06:S06_R02',
#           'S06:S06_R03', 'S06:S06_R04', 'S06:S06_R05', 'S07:S07_R01',
#           'S07:S07_R02', 'S07:S07_R03', 'S07:S07_R04', 'S07:S07_R05',
#           'S08:S08_R01', 'S08:S08_R02', 'S08:S08_R03', 'S08:S08_R04',
#           'S08:S08_R05', 'S09:S09_R01', 'S09:S09_R02', 'S09:S09_R03',
#           'S09:S09_R04', 'S09:S09_R05', 'S10:S10_R01', 'S10:S10_R02',
#           'S10:S10_R03', 'S10:S10_R04', 'S10:S10_R05']
    old_cols = df.columns[df.columns.str[:1] == "S"]
    repl_cols_s = old_cols.str[:3].unique()
    old_cols = ["protein"] + old_cols.tolist()
    if keepRuns == False:
#        repl_cols_s = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]        
        cols_allRun = list(np.array([[i,i,i,i,i] for i in repl_cols_s]).flatten())
        #old_cols = df.drop(["id", "specie", "protein", 'q_value', 'posterior_error_prob', 'num_peptides',
        #   'protein_id_posterior_error_prob', 'log2_fold_change',
        #   'diff_exp_prob_'+str(diff_exp), 'peptides', 'decoy'], axis = 1).columns
        renameColsMap = dict(zip(old_cols, cols_allRun))
        df = df[old_cols].set_index(["protein"]).rename(columns = renameColsMap)  
        #df2 = df.drop(["id", "specie", 'q_value', 'posterior_error_prob', 'num_peptides',
        #   'protein_id_posterior_error_prob', 'log2_fold_change',
        #   'diff_exp_prob_'+str(diff_exp), 'peptides', 'decoy'], axis = 1).set_index(["protein"]).rename(columns = renameColsMap)
    else:
        df = df[old_cols].set_index(["protein"]) 
    df = df.convert_objects(convert_numeric=True)
    return df


#############
## PARSERS ##
#############
'''    
def parsePSSS3equDecoy(filename = '500-PSSS3-equ decoy_Report.xls'):
    """
    File for parsing PSSS3 unfiltered with decoys and PG.Quantity column. (Very slow)
    """
    F = open(filename)
    
    header = F.readline().split("\t")[1:8] #not used, only used to skip first row of readin
    lines = []
    
    organisms = ['Arabidopsis thaliana',
                 'Caenorhabditis elegans OX=6239',
                 'Homo sapiens']
    
    #for i in range(3000):
    for line in F:    
        #line = F.readline()
        element = line.split("\t")[1:8]
        
        # Ignore shared proteins
        if ";" in element[2]:
            continue
        
        reformatted_line = []
        
        #lines.append(element)
        # Make identifier for sample and runs
        sample = element[0].split("_")[2].split("-")[2]
        run = element[0].split("_")[-2]
        run_id = sample + ":" + sample + "_" + run
        reformatted_line.append(run_id)
        
        # Species and Protein
        if element[1] == organisms[0]:
            specie = "ARATH"
        elif element[1] == organisms[1]:
            specie = "CAEEL"
        elif element[1] == organisms[2]:
            specie = "HUMAN"
        protein = element[2] + "_" + specie
        reformatted_line.append(protein)
        
        # PG.Cscore
        reformatted_line.append(float(element[3]))
        
        # PG.NrOfStrippedSequences
        reformatted_line.append(float(element[4]))
        
        # Q-value
        reformatted_line.append(float(element[5]))
        
        # Protein Quantity
        reformatted_line.append(float(element[6]))
        lines.append(reformatted_line)
    
    new_header = ["run_id", "protein", "Cscore", "NStrippedSequences", "Qvalue", "ProteinQuantity"]
    df = pd.DataFrame(lines, columns=new_header)
    return df
'''
#df = psss3.head(1000)

'''
def parsePSSS3equDecoy(filename = '500-PSSS3-equ decoy_Report.xls'):
    """
    File for parsing PSSS3 unfiltered with decoys and PG.Quantity column. (Very slow)
    """
    df = pd.read_csv(filename, sep = "\t", usecols = [1,2,3,4,5,6,7,9])
    df = df.dropna() #dropna from raw file to make processing feasible <--
    
    reformatted_df = pd.DataFrame()
    
    samples = df["R.FileName"].str.split("-", n=3, expand = True)[2].str.split("_", n=4, expand=True)[0]
    runs = df["R.FileName"].str.split("-", n=3, expand = True)[2].str.split("_", n=4, expand=True)[2]
    
    run_id = samples + ":" + samples + "_" + runs
    
    specie = df["PG.Organisms"].apply(lambda x: "CAEEL" if x == "Caenorhabditis elegans OX=6239" 
             else "ARATH" if x == "Arabidopsis thaliana" else  "HUMAN" if x == 'Homo sapiens' else x)
    
    # remove all shared proteins
    proteinId = df["PG.ProteinAccessions"].apply(lambda x: np.nan if ";" in x else x)
    
    protein = proteinId + "_" + specie
    
    cScore = df["PG.Cscore"]
    nrStripped = df["PG.NrOfStrippedSequencesIdentified"]
    qval = df["PG.Qvalue"]
    proteinQuantity = df["PG.Quantity"]
    decoy = df["EG.IsDecoy"]
    
    reformatted_df["run_id"] = run_id
    reformatted_df["protein"] = protein
    reformatted_df["Cscore"] = cScore
    reformatted_df["nrStrippedSequences"] = nrStripped
    reformatted_df["qvalue"] = qval
    reformatted_df["proteinQuantity"] = proteinQuantity
    reformatted_df["decoy"] = decoy
    
    reformatted_df = reformatted_df.dropna() # Bascially drop all shared proteins. (= NaN from previous row)
    
    # Drop duplicates
    reformatted_df = reformatted_df.drop_duplicates()

    # Sort by run_id
    # reformatted_df = df.sort_values(by = "run_id")

    return reformatted_df
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
    df["proteinId"] = proteinId
    df = df[pd.notnull(df["proteinId"])]
    proteinId = df["PG.ProteinAccessions"]
    
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
    run_id = run_id.apply(lambda x: "S03:S04_R05" if x == 'MHR:MHR_X01' else x) #special case for S03:S04_R05...
    
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
    #reformatted_df = reformatted_df.drop_duplicates()

    # Sort by run_id
    # reformatted_df = df.sort_values(by = "run_id")

    return reformatted_df

def minimizePSSS3equDecoy(filename = "500-PSSS3-raw-reformatted_dropna_dropdup.csv"):
    """
    Remove Cscore, nrStrippedSequences, qvalue from reformatted_df from parsePSSS3equDecoy
    """
    df = pd.read_csv(filename, sep ="\t")
    df = df.drop(['Cscore', 'nrStrippedSequences', 'qvalue'], axis = 1)
    return df


# Split the protein columns to specie and protein columns
def getTriqlerDetailed(triqler):
    """
    Input triqler to get detailed columns
    
    Example:    
        triqler = table2Df("PSSS3_triqler_output_proteins.1vs2.tsv", peptideSeperator = True)
        triqler = getTriqlerDetailed(triqler)
    """
    triqler = triqler.convert_objects(convert_numeric = True)    
    specie = triqler["protein"].str[-5:]
    triqler["specie"] = specie
    decoy = triqler["protein"].str[:-6].str.split("_", expand=True)
    triqler["decoy"] = (decoy[0] == "decoy")
    protein_list1 = decoy[decoy[0] != "decoy"][0]
    protein_list1.fillna(value=pd.np.nan, inplace = True)
    protein_list2 = decoy[1]
    protein_list2.fillna(value=pd.np.nan, inplace = True)
    triqler["prot"] = protein_list2
    protein_list_df = (pd.DataFrame(data = protein_list1.values,index = protein_list1.index, columns = ["prot"]))
    tritest = pd.concat([triqler, protein_list_df], axis = 1, join = "outer")["prot"]
    def sjoin(x): return ';'.join(x[x.notnull()].astype(str))
    protein_list = tritest.groupby(level=0, axis=1).apply(lambda x: x.apply(sjoin, axis=1))
    triqler["id"] = triqler["protein"]
    triqler.drop("protein", axis = 1, inplace = True)
    triqler["protein"] = protein_list
    triqler.drop("prot",axis = 1, inplace = True)
    #triqler.columns = list(triqler.columns[-2:]) + list(triqler.columns[-4:-3]) + list(triqler.columns[-3:-2]) + list(triqler.columns[-5:-4]) + list(triqler.columns[:-5])
    return triqler

def getTriqler_FDR(triqler, PEP_treshold = 0.05):
    """
    Take Triqler input and get FDR for the set using protein_id_posterior_error_prob
    """
    correctPSMs = triqler[triqler["decoy"] == False]
    FDR = correctPSMs[correctPSMs["protein_id_posterior_error_prob"] < PEP_treshold]["protein_id_posterior_error_prob"].sum()/len(correctPSMs)
    return FDR

#findFDR
def findFDR(triqler, target_FDR = 0.01, lower_PEP = 0, upper_PEP = 1):
    """
    Functions find the optimal PEP value for the target FDR
    """
    triqler = triqler
    target = target_FDR
    def findFDR_func(x):
        res = getTriqler_FDR(triqler, PEP_treshold = x) - target
        return res
    optPEP= sp.optimize.bisect(findFDR_func, a = lower_PEP, b = upper_PEP)
    return optPEP 

#Convert Triqler to numeric
def triqler2Numeric(triqler):
    """
    Input:
        triqler - triqler dataframe
    example:
    triqler = table2Df("proteins.1vs2.tsv", peptideSeperator = True)
    triqler = triqler2Numeric(triqler)
    """
    triqler[list(triqler.columns[:2])+list(triqler.columns[3:-1])] = triqler[list(triqler.columns[:2])+list(triqler.columns[3:-1])].apply(pd.to_numeric)
    return triqler

##################################################
#### MAIN ########################################
##################################################
##################################################
    
if __name__ == "__main__":
    
    
    PSSS3 = pd.read_csv("500-PSSS3-precursor_Report.proteins.tsv", sep = "\t")
    
    
    triqler = table2Df("PSSS3_triqler_output_proteins.1vs2.tsv", peptideSeperator = True)
    triqler2 = table2Df("PSSS3_triqler_output_proteins.1vs3.tsv", peptideSeperator = True)

    mixtures = getMixtures()
    
    fractions = getFractions(triqler)
    means, stds = meanVarFractions(fractions) #data for mean-var plots
    conditionsData = getAllConditionFractions(fractions) #data for plots containing all runs.
    df = conditionsData2Df(conditionsData)
    
    fractions_PSSS3 = getFractions(PSSS3)
    conditionsData_PSSS3 = getAllConditionFractions(fractions_PSSS3) #data for plots containing all runs.
    df_PSSS3 = conditionsData2Df(conditionsData_PSSS3)
    
    ax = df[df["specie"] == "HUMAN"].boxplot(column = "fraction", by = "condition", return_type='axes')
    
    boxprops_PSSS3 = dict(linestyle='-', linewidth=2, color='r')
    medianprops_PSSS3 = dict(linestyle='-', linewidth=2, color='r')
    bp = df_PSSS3[df_PSSS3["specie"] == "HUMAN"].boxplot(column = "fraction", by = "condition", 
            boxprops = boxprops_PSSS3, medianprops = medianprops_PSSS3, return_type = 'dict', ax = ax) #correct plotting something wrong in code!
    
    # boxplot style adjustments
    [[item.set_linewidth(2) for item in bp[key]['boxes']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['fliers']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['medians']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['means']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['whiskers']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['caps']] for key in bp.keys()]
    
    [[item.set_color('b') for item in bp[key]['boxes']] for key in bp.keys()]
    # seems to have no effect
    [[item.set_color('k') for item in bp[key]['fliers']] for key in bp.keys()]
    [[item.set_color('r') for item in bp[key]['medians']] for key in bp.keys()]
    [[item.set_markerfacecolor('g') for item in bp[key]['means']] for key in bp.keys()]
    [[item.set_color('k') for item in bp[key]['whiskers']] for key in bp.keys()]
    [[item.set_color('b') for item in bp[key]['caps']] for key in bp.keys()]
    
    
    mixtures[mixtures.index == "HUMAN"].boxplot()
    
    plt.ylim(0,1)
    plt.title("Homo Sapiens")
    
    df_PSSS3[df_PSSS3["specie"] == "HUMAN"].boxplot(column = "fraction", by = "condition") #correct plotting something wrong in code!
    
    
    df[df["specie"] == "CAEEL"].boxplot(column = "fraction", by = "condition")
    
    bp = df_PSSS3[df_PSSS3["specie"] == "CAEEL"].boxplot(column = "fraction", by = "condition", 
            boxprops = boxprops_PSSS3, medianprops = medianprops_PSSS3, return_type = 'dict', ax = ax) #correct plotting something wrong in code!
    
    ####
    
    ax = df[df["specie"] == "CAEEL"].boxplot(column = "fraction", by = "condition", return_type='axes')
    
    boxprops_PSSS3 = dict(linestyle='-', linewidth=2, color='r')
    medianprops_PSSS3 = dict(linestyle='-', linewidth=2, color='r')
    bp = df_PSSS3[df_PSSS3["specie"] == "CAEEL"].boxplot(column = "fraction", by = "condition", 
            boxprops = boxprops_PSSS3, medianprops = medianprops_PSSS3, return_type = 'dict', ax = ax) #correct plotting something wrong in code!
    
    # boxplot style adjustments
    [[item.set_linewidth(2) for item in bp[key]['boxes']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['fliers']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['medians']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['means']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['whiskers']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['caps']] for key in bp.keys()]
    
    [[item.set_color('b') for item in bp[key]['boxes']] for key in bp.keys()]
    # seems to have no effect
    [[item.set_color('k') for item in bp[key]['fliers']] for key in bp.keys()]
    [[item.set_color('r') for item in bp[key]['medians']] for key in bp.keys()]
    [[item.set_markerfacecolor('g') for item in bp[key]['means']] for key in bp.keys()]
    [[item.set_color('k') for item in bp[key]['whiskers']] for key in bp.keys()]
    [[item.set_color('b') for item in bp[key]['caps']] for key in bp.keys()]
    
    mixtures[mixtures.index == "CAEEL"].boxplot()
    plt.ylim(0,1)
    plt.title("Caenorhabditis elegans")
    
    df[df["specie"] == "ARATH"].boxplot(column = "fraction", by = "condition")
    
    
    ax = df[df["specie"] == "ARATH"].boxplot(column = "fraction", by = "condition", return_type='axes')
    
    boxprops_PSSS3 = dict(linestyle='-', linewidth=2, color='r')
    medianprops_PSSS3 = dict(linestyle='-', linewidth=2, color='r')
    bp = df_PSSS3[df_PSSS3["specie"] == "ARATH"].boxplot(column = "fraction", by = "condition", 
            boxprops = boxprops_PSSS3, medianprops = medianprops_PSSS3, return_type = 'dict', ax = ax) #correct plotting something wrong in code!
    
    # boxplot style adjustments
    [[item.set_linewidth(2) for item in bp[key]['boxes']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['fliers']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['medians']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['means']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['whiskers']] for key in bp.keys()]
    [[item.set_linewidth(2) for item in bp[key]['caps']] for key in bp.keys()]
    
    [[item.set_color('b') for item in bp[key]['boxes']] for key in bp.keys()]
    # seems to have no effect
    [[item.set_color('k') for item in bp[key]['fliers']] for key in bp.keys()]
    [[item.set_color('r') for item in bp[key]['medians']] for key in bp.keys()]
    [[item.set_markerfacecolor('g') for item in bp[key]['means']] for key in bp.keys()]
    [[item.set_color('k') for item in bp[key]['whiskers']] for key in bp.keys()]
    [[item.set_color('b') for item in bp[key]['caps']] for key in bp.keys()]
    mixtures[mixtures.index == "ARATH"].boxplot()
    plt.ylim(0,1)
    plt.title("Arabidopsis thaliana")
    
    
    import seaborn as sns
    fig, ax = plt.subplots(figsize=(13, 8))
    sns.boxplot(data=sampgroup.loc[sampgroup['Species'] == "H sapiens"], x="Sample",y="Expression", hue="Method",fliersize=0)
    for i in range(len(h_mix)):
        ax.plot([i-.2, i+.2], [h_mix[i], h_mix[i]], '-', color='red', lw=5)
    ax.set_ylim(0,0.15)
    ax.set_ylabel("Fraction of protein total abundance per sample")
    sb.despine(offset=10, trim=True)
    plt.title("Homo sapiens")
    fig.savefig("sapiens_sample.png")
    
    #import seaborn as sns
    #ax = sns.boxplot(x="condition", y="fraction", data=df[df["specie"] == "HUMAN"])
    
    #df = fractions2Df(fractions)
    
    df.boxplot(column = "fraction", by = "condition")
    #temp append con
    conditionFractions.keys()
    conditionFractions.values()
    formatted_df.apply(pd.to_numeric)
    
    df = PSSS3_triqler
    PSSS3_triqler["S01:S01_R01"].plot.line()
    
    y  = PSSS3_triqler["S01:S01_R01"]
    x = y.index
    
    plt.plot(x[:100], y[:100], '-')
    
    500-PSSS3-precursor_Report.proteins.tsv
    
    
    
    
    #???? BELOW??
    a_mix = a_mix/a_mix.sum()
    c_mix = c_mix/c_mix.sum()
    h_mix = h_mix/h_mix.sum()
    mix = {'A thaliana':a_mix, 'C elegans': c_mix, 'H sapiens': h_mix}
    
