#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 18:50:31 2023

@author: blim
"""

import numpy as np
import pandas as pd

# Load your DataFrame
df = pd.read_excel('All_RBD.xlsx', index_col=0)
#%%
# Convert each row into a string
df['String'] = df.apply(lambda row: ''.join(map(str, row)), axis=1)
# Start a counter for identical sequences
counter = df['String'].value_counts()

print(counter)
