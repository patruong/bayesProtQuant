#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 09:04:59 2019

@author: ptruong
"""
import pandas as pd

def prettyPrint(df):
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(df)
