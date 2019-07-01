#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 20:56:46 2019

@author: ptruong
"""

import pandas as pd
from proc import * 

from os.path import expanduser
usr = open("../usr.txt", "r")
usr = usr.read().split("\n")[:-1]

import plotly
plotly.__version__
plotly.tools.set_credentials_file(username=usr[0], api_key=usr[1])
from plotly.offline import plot
import plotly.plotly as py
import plotly.graph_objs as go


triqler = pd.read_pickle("triqler.pkl")
spectronaut = pd.read_pickle("spectronaut.pkl")
at_t, ce_t, hs_t = splitTriqlerBySpecies(triqler, exponential = False, truncated = False, keepRuns = False)
at_s, ce_s, hs_s = splitSpectronautBySpecies(spectronaut, truncated = False, keepRuns = False)

# normalize 
at_t = at_t.div(at_t.sum(axis=1), axis=0)
ce_t = ce_t.div(ce_t.sum(axis=1), axis=0)
hs_t = hs_t.div(hs_t.sum(axis=1), axis=0)
at_s = at_s.div(at_s.sum(axis=1), axis=0)
ce_s = ce_s.div(ce_s.sum(axis=1), axis=0)
hs_s = hs_s.div(hs_s.sum(axis=1), axis=0)

# get true values
a_mix, c_mix, h_mix, mix = getMixtures()
at_true = pd.DataFrame(np.repeat(a_mix, 5), index = at_t.columns).T/5
ce_true = pd.DataFrame(np.repeat(c_mix, 5), index = ce_t.columns).T/5
hs_true = pd.DataFrame(np.repeat(h_mix, 5), index = hs_t.columns).T/5
#at_t = at_t.rename(index=str, columns = dict(zip(at_t.columns, at_t.columns.str[:3]+":" +at_t.columns.str[-3:])))
#at_s = at_s.rename(index=str, columns = dict(zip(at_s.columns, at_s.columns.str[:3]+":" +at_s.columns.str[-3:])))

def concatenatedDF(df_s, df_t, df_true):
    df_s["method"] = "spectronaut"
    df_s["method_id"] = 1
    df_t["method"] = "triqler"
    df_t["method_id"] = 2
    df_true["method"] = "true"
    df_true["method_id"] = 3
    df = pd.concat([df_t, df_s, df_true])
    return df

df = concatenatedDF(hs_s, hs_t, hs_true)
# parallell coordinates
parDim = []
for i in range((df.columns.str[:1] == "S").sum()): 
    parDim.append(dict(range = [0,1], label = df.columns[i], values = df.iloc[:,i]))
data = [
    go.Parcoords(
        line = dict(color = df['method_id'],
                   colorscale = [[0,'orange'],[0.5,'blue'],[1,'red']]),
        dimensions = parDim
    )
]


layout = go.Layout(
    plot_bgcolor = '#E5E5E5',
    paper_bgcolor = '#E5E5E5'
)

plot(data, filename="hs.html")

# Find max
at_s[at_s.iloc[:,1] == at_s.iloc[:,1].max()]

df[df.iloc[:,36] > 0.8]