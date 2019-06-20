#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 10:28:26 2019

@author: ptruong
"""

import numpy as np
import pandas as pd


def findRowLengths(fileRows):  
    """
    Counts the unique row lengths in a list of rows.
    
    Input:
        fileRows - list of rows from readFile in triqlerParser.py
    """
    unique_row_length = []
    for i in range(len(fileRows)):
        length = len(fileRows[i])
       # if length == 28:
       #     print(fileRows[i-1])
       #     print(fileRows[i])
       #     print(i)
        if length not in unique_row_length:
            unique_row_length.append(length)
    return unique_row_length

def readFile(filename):
    """
    read in filename and return header and rows of the file as lists.
    """       
    f = open(filename, "r")
    f_header = f.readline().replace("\n","").split("\t")
    f_rows = []
    count = 0
    for i in f:
        count += 1
        #print(i)
        #if count == 10:
        #    break
        #f_row = f.readline().replace("\n","").split("\t") #python 3.6
        f_row = i.replace("\n","").split("\t")
        f_rows.append(f_row)
    return f_header, f_rows
 

def readFile2df(header, rows, peptideSeperator = False):
    """
    Takes header list and rows list as input to create a pandas df.
    
    Input:
        header - header list
        rows - list or rows (list of lists)
        peptideSeperator - specify if peptide column should be seperated, If False default peptide are peptides in list format.
    """
    
    lenHeader = len(header)
    data = []
    for i in range(len(rows)):    
        column = rows[i][:lenHeader-1]
        peptides = rows[i][lenHeader-1:]
        if peptideSeperator == True:
            peptides = "|".join(peptides)
        elif peptideSeperator != False:
            peptides = peptideSeperator.join(peptides)
        else:
            pass
        column.append(peptides)
        data.append(column)
    df = pd.DataFrame(data, columns = header)
    return df

def table2Df(filename, peptideSeperator = False):
    header, rows = readFile(filename)
    df = readFile2df(header, rows, peptideSeperator = peptideSeperator)
    return df


