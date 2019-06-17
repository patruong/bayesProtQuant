#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 12:50:43 2019

@author: ptruong
"""

dt = pd.DataFrame()
dt['a'] = [1,2, np.nan, 3]
dt['b'] = [4,5, np.nan, 6]
dt
dt.sum(axis=1)

dt['total'] = dt.sum(axis=1) 
dt.loc[dt['a'].isnull() & dt['b'].isnull(),'total']=np.nan

dt['total'] = dt.sum(axis=1)    
dt.loc[dt[['a','b']].isnull().all(1),'total']=np.nan