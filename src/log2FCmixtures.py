#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 16:11:26 2019

@author: ptruong
"""
"""
log2fc
http://rstudio-pubs-static.s3.amazonaws.com/13988_bb11d85b79b2436280de434988558140.html
https://www.quora.com/How-do-you-interpret-Log-Fold-Change-with-Gene-Expression
"""

from proc import *
import numpy as np

a_mix, h_mix, c_mix, mix = getMixtures()


a = 10
b = 100
def log2fc(a,b):
    """
    fold change a to b
    """
    fc = b/a
    return np.log2(fc)


#(h_mix[0:-1]/h_mix[1:])-1
#(c_mix[0:-1]/c_mix[1:])-1
#(a_mix[0:-1]/a_mix[1:])-1

# log2-foldchange formula log2FC - log2(B) - log2(A)
np.log2(h_mix[0:-1])-np.log2(h_mix[1:])
np.log2(c_mix[0:-1])-np.log2(c_mix[1:])
np.log2(a_mix[0:-1])-np.log2(a_mix[1:])
