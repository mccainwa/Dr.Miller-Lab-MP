#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 14:05:18 2023

@author: blim
"""

import pandas as pd

#Ask for input and output file path
input_excel_file = input("Enter the input Excel file path: ")
output_fasta_file = input("Enter the output FASTA file path: ")
df = pd.read_excel(input_excel_file)

# Open the output FASTA file
with open(output_fasta_file, 'w') as file:
    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
        header = f'>{row.iloc[0]}'  # Assuming the first column contains file names
        sequence = ''.join(str(aa) for aa in row.iloc[1:])  # Convert values to strings and join them
        
        file.write(f'{header}\n{sequence}\n\n')

print("Data successfully processed and saved to FASTA.")