#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 16:02:51 2023

@author: blim
"""

import pandas as pd

input_file = "/home/blim/Documents/Summer2023/HPV/HPV_complete.txt"
output_file = "/home/blim/Documents/Summer2023/HPV/HPV_extracted.txt"
output_excel = "/home/blim/Documents/Summer2023/HPV/HPV_extracted.xlsx"
proteins_to_identify = ["L1"]
protein_data = {}

with open(input_file, 'r') as file:
    lines = file.readlines()

protein_id = None
protein_sequence = []
sample_counter = 1

for line in lines:
    line = line.strip()
    if line.startswith('>'):
        if protein_id and protein_sequence:
            protein_data[protein_id] = ''.join(protein_sequence)
        protein_id = line[1:]
        protein_sequence = []
    else:
        protein_sequence.append(line)

# Create a list to store the formatted sequences
formatted_sequences = []

# Format the sequences with the incremented headers and blank lines
for protein_id, protein_sequence in protein_data.items():
    for protein_to_identify in proteins_to_identify:
        if protein_to_identify in protein_id:
            header = f"L1_{sample_counter}"
            formatted_sequence = f">{header}\n{protein_sequence}\n\n"
            if 'X' not in protein_sequence and 'B' not in protein_sequence and 'J' not in protein_sequence and (500 <= len(protein_sequence) <= 530):
                formatted_sequences.append(formatted_sequence)
                sample_counter += 1
            break

# Write the formatted sequences to the output file
with open(output_file, 'w') as file:
    file.writelines(formatted_sequences)

# Create a DataFrame from the extracted protein data
df = pd.DataFrame(protein_data.items(), columns=['Header', 'Sequence'])

# Add the incremented headers to the DataFrame
df['Header'] = [f"L1_{i+1}" for i in range(len(df))]

# Filter out rows with sequence length not between 500 and 530, and exclude rows with 'X', 'B', or 'J'
df = df[(df['Sequence'].str.len().between(500, 530)) & (~df['Sequence'].str.contains('X|B|J', regex=True))]

# Save the DataFrame to an Excel file
df.to_excel(output_excel, index=False)

