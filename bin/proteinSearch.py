#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 10:50:32 2019

@author: ptruong
"""

import os
import pandas as pd
import numpy as np

#######
# CMD #
#######

def findUniqueProteins(inputName = "500-PSSS3-equ\ decoy_Report.xls", outputName = "proteins.csv"):
    """
    Generates csv of unique proteins from PSSS3 file.
    """
    inputName = inputName
    outputName = outputName
    cmd = "cat " + inputName + " "+ "| cut -f 4 | sort | uniq -c >" + " " + outputName
    os.system(cmd)

############
# Database #
############

def getUniProt(filename = "uniprot_sprot_2017-10-25_HUMAN_ISOFORMS.fasta"):
    """
    Creates a df with UniProt protein accessions and corresponding protein.
    """
    f = open(filename, "r")
    libraryAccessions = []
    libraryProteins = []
    librarySequences = []
    peptideSequence = ""
    for line in f:
        if line[0] == ">":
            # if not first run, append to df
            if peptideSequence != "":
                librarySequences.append(peptideSequence)
            peptideSequence = "" # reset the peptide sequence
            
            libraryInfo = line.split("\n")[0].split("|")
            libraryAccession = libraryInfo[1]
            libraryProtein = libraryInfo[2]
            libraryAccessions.append(libraryAccession)
            libraryProteins.append(libraryProtein)
        else:
            peptideSequence += line.split("\n")[0]
    librarySequences.append(peptideSequence) #include last sequence...
    
    data = np.array([libraryAccessions, libraryProteins, librarySequences]).T
    uniprot = pd.DataFrame(data, columns = ["proteinAccession", "protein", "sequence"])
    return uniprot

def getPSSS3Uniprot(libraryDirectory = "library/LKaell/", 
                    humanLib = "uniprot_sprot_2017-10-25_HUMAN_ISOFORMS.fasta",
                    arathLib = "uniprot_sprot_2018-01-24_ARATH_ISOFORMS.fasta",
                    caeelLib = "uniprot-elegans-filtered-organism__Caenorhabditis+elegans.fasta"):
    """
    Get concatenated library for HUMAN, ARATH and CAEEL.
    """
    libraryDirectory = libraryDirectory
    species = ["HUMAN",
               "ARATH",
               "CAEEL"]
    libraryFiles = [humanLib,
                    arathLib,
                    caeelLib]
    uniprot = pd.DataFrame()
    for i in species:
        if i == "HUMAN":
            libraryFile = humanLib
        if i == "ARATH":
            libraryFile = arathLib
        if i == "CAEEL":
            libraryFile = caeelLib
        uniprot_spec = getUniProt(libraryDirectory + libraryFile)
        uniprot_spec["specie"] = [i for j in range(len(uniprot_spec))]
        uniprot = pd.concat([uniprot, uniprot_spec])
    return uniprot
    

def findSingleProteins(filename = "proteins.csv"):
    """
    Find all single protein proteinAccession number. 
    
    filename - .csv generated file from findUniqueProteins output.
    """
    
    filename = filename
    f = open(filename, "r")
    # Find all the singleton proteins.
    singleProteins = []
    for line in f:
        count, proteinAccession = line.split("\n")[0].split()
        count = int(count)
        proteins = proteinAccession.split(";") #make a list to append all proteins and one to append unique proteins.
        num_proteinAccessions = len(proteins)
        if num_proteinAccessions == 1:
            singleProteins.append(proteins[0])
    return singleProteins

def findShared(filename = "proteins.csv"):
    """
    Find all shared and non-shared proteins among those with multi protein accession numbers.
    
    filename - .csv generated file from findUniqueProteins output.
    """
    filename = filename
    f = open(filename, "r")
    singleProteins = findSingleProteins(filename = filename)
    
    nonShared = []
    shared = []
    #allNonSingleProt = []
    
    for line in f:
        count, proteinAccession = line.split("\n")[0].split()
        count = int(count)
        proteins = proteinAccession.split(";")
        num_proteinAccessions = len(proteins)
        if num_proteinAccessions > 1:
            isShared = False
            for i in proteins:
                if i in singleProteins:
                    isShared = True
            #allNonSingleProt.append(proteinAccession)
            if isShared == False:
                nonShared.append(proteinAccession)
            if isShared == True:
                shared.append(proteinAccession)
    return shared, nonShared

def areProteinAccessionSimilar(protein, uniprot, matchLength = 10):
    """
    Check if proteinAccession numbers are similar or same. If require exact match in protein sequences, use
    matchLength = "full" 
    
    protein - shared proteinAccession number (e.g. 'Q9Y4B6;Q9Y4B6-2') NOTE: requires ";" in the string.
    """
    proteins = protein.split(";")
    matchLength = matchLength
    #If require perfect match
    if matchLength == "full":
        matchLength = len(uniprot[uniprot["proteinAccession"] == protein.split(";")[0]].sequence.values[0])
    if matchLength > len(uniprot[uniprot["proteinAccession"] == protein.split(";")[0]].sequence.values[0]):
        raise Exception("matchLength is exceeds protein sequence length! Protein sequence length is " +
                       str(len(uniprot[uniprot["proteinAccession"] == protein.split(";")[0]].sequence.values[0])))
    matchingSeq = ""
    similarProtein = True
    for i in proteins:
    #protein = proteins[0]
        foundSequence = uniprot[uniprot["proteinAccession"] == i]
        sequence = foundSequence["sequence"]
        # Let the first protein sequence set the matching standard.
        if matchingSeq == "":
            matchingSeq = sequence.str[:matchLength].values[0]
        proteinSeq = sequence.str[:matchLength].values[0]
        if proteinSeq != matchingSeq:
            similarProtein = False
            return similarProtein
    return similarProtein 

def findSimilarProteins(nonSharedProteinAccessions, uniprot, matchLength = 10):
    """
    Creates df with nonShared proteins and boolean table of if the proteinAccessions are similar or not.
    
    nonShared - list of shared proteinAccessions (e.g. 'Q9Y4B6;Q9Y4B6-2') NOTE: requires ";" in the string.
    """
    nonShared = nonSharedProteinAccessions
    uniprot = uniprot
    matchLength = matchLength
    df = pd.DataFrame(nonShared, columns = ["proteinAccession"])
    df["similarProteins"] = df["proteinAccession"].apply(lambda x: areProteinAccessionSimilar(x, uniprot = uniprot, matchLength = matchLength))
    return df

def mapNonShared(nonShared):
    """
    Take nonShared proteinAccessions as input and map to first proteinAccession
    """
    concat_nonShared = []
    for i in nonShared:
        proteinAccession = i.split(";")[0]
        concat_nonShared.append(proteinAccession)
    return concat_nonShared 

if __name__ == "__main__":
    #df = pd.read_csv("proteins.csv", sep = "\s", names  = ["count", "proteinAccession"])

    uniprot = getPSSS3Uniprot(libraryDirectory = "library/LKaell/", 
                    humanLib = "uniprot_sprot_2017-10-25_HUMAN_ISOFORMS.fasta",
                    arathLib = "uniprot_sprot_2018-01-24_ARATH_ISOFORMS.fasta",
                    caeelLib = "uniprot-elegans-filtered-organism__Caenorhabditis+elegans.fasta")
    singleProteins = findSingleProteins(filename="proteins.csv")
    shared, nonShared = findShared(filename = "proteins.csv") 
    
    filename = "500-PSSS3-equ decoy_Report.xls"
    outputname = "500-PSSS3-equ decoy_Report_nonShared_20190507.xls"
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
        
    #Shared are kept as shared, nonShared needs to be checked for protein sequence similarity
    #matchLength = 10 #treshold for sequence start match
    #df_similarity = findSimilarProteins(nonSharedProteinAccessions = nonShared, uniprot = uniprot, matchLength = matchLength)
    
    #df_similarity[df_similarity["similarProteins"] == True]["proteinAccession"]


