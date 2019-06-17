#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 12:30:12 2019

@author: ptruong
"""

from triqlerParser import *

iRPG = table2Df("iPRG2016.tsv", peptideSeperator = True)
r1 = table2Df("proteins.1vs2.tsv")
r2 = table2Df("proteins.1vs3.tsv")
r3 = table2Df("proteins.2vs3.tsv")


header, rows = readFile("iPRG2016.tsv")
f_length = findRowLengths(f_rows)

r1_header, r1_row = readFile("proteins.1vs2.tsv")
r1_length = findRowLengths(r1_row)

r2_header, r2_row = readFile("proteins.1vs3.tsv")
r2_length = findRowLengths(r2_row)

r3_header, r3_row = readFile("proteins.2vs3.tsv")
r3_length = findRowLengths(r3_row)