#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 12:38:49 2023

@author: blim
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 17:25:10 2023

@author: blim
"""

import pandas as pd
import glob

# Prompt input directory and output file name
directory = input("Enter the directory containing the Excel files: ")
output_file_name = input("Enter the output Excel file name: ")

# List all Excel files in the directory
file_list = glob.glob(directory + '/*.xlsx')
dataframes = []

# Iterate over each file and read it into a DataFrame
for file in file_list:
    df = pd.read_excel(file)
    dataframes.append(df)

# Merge all DataFrames into a single DataFrame
merged_df = pd.concat(dataframes, ignore_index=True)
merged_df.to_excel(output_file_name, index=False)
print(merged_df)
