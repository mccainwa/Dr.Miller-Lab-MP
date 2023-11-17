#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 17:37:53 2023

@author: blim
"""

import pandas as pd

fasta_file = "/home/blim/Documents/Summer2023/HPV/HPV_aligned.fasta"

data = []  # List to store sequence data

with open(fasta_file, 'r') as file:
    lines = file.readlines()

sequence = ""
header = None

for line in lines:
    line = line.strip()
    if line.startswith('>'):
        if header is not None:
            data.append((header, sequence))
            sequence = ""
        header = line[1:]
    else:
        sequence += line

# Append the last sequence
data.append((header, sequence))

# Get the maximum sequence length
max_length = max(len(seq) for _, seq in data)

# Create a dictionary to store the sequence data for each sample
sequence_dict = {}
for i, (header, sequence) in enumerate(data):
    # Pad the sequence with '-' to match the maximum length
    padded_sequence = sequence + '-' * (max_length - len(sequence))
    sequence_dict[f"Sample_{i+1}"] = list(padded_sequence)

# Create a DataFrame from the sequence data
df = pd.DataFrame(sequence_dict).T

new_df = df.iloc[:, 27:315]

# Assuming your DataFrame is named df
new_df = new_df.replace('-', pd.NA)  # Replace '-' with missing values

# Calculate the mode of each column
column_modes = new_df.mode()

# Replace missing values with the mode of each column
new_df = new_df.fillna(column_modes.iloc[0])

# If you want to replace missing values with the mode of each column, excluding the first column
# df.iloc[:, 1:] = df.iloc[:, 1:].fillna(column_modes.iloc[0])


# Save the DataFrame to a CSV file
new_df.to_excel("/home/blim/Documents/Summer2023/HPV_Trial.xlsx", index=False)
