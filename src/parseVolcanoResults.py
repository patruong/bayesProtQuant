#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 11:07:57 2019

@author: ptruong
"""

import pandas as pd
#df= pd.read_csv("volcanoResults.csv", sep = "\t", index_col = 0)

def readVolcanoResultsDict(filename = "volcanoResultsDict.txt"):
    f = open(filename, "r")
    data = {}
    for i in f:
        i = i.strip("\n")
        method, comp, extra = i.split("_")
        specie, vals = extra.split(":")
        vals = vals.strip("(").strip(")").split(",")
        vals = (float(vals[0]), float(vals[1]))
        data.update({method+"_"+comp+"_"+specie:vals})
    return data

def volcanoResultsDict_to_meltedDf(filename = "volcanoResultsDict.txt"):
    f = open(filename, "r")
    data= {}
    data.update({"method": []})
    data.update({"specie": []})
    data.update({"comparison": []})
    data.update({"diffExp": []})
    data.update({"diffExpRatio": []})
    for i in f:
        row = []
        i = i.strip("\n")
        method, comp, extra = i.split("_")
        specie, vals = extra.split(":")
        vals = vals.strip("(").strip(")").split(",")
        diffExp = float(vals[0])
        deRatio = float(vals[1])
        row.extend([method, specie, comp, diffExp, deRatio])
        for j in range(len(data.keys())): 
            data[list(data.keys())[j]].append(row[j])
    df = pd.DataFrame.from_dict(data)
    return df

#data = readVolcanoResultsDict(filename = "volcanoResultsDict.txt")
#df = volcanoResultsDict_to_meltedDf(filename = "volcanoResultsDict.txt")










