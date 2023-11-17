#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 14:36:17 2023

@author: blim
"""

import pandas as pd
from Bio import SeqIO

# Specify the input and output file
fasta_file = input("Enter the input FASTA file path: ")
excel_file = input("Enter the output Excel file path: ")

records = list(SeqIO.parse(fasta_file, 'fasta'))
data = {'File': []}
sequences = []

# Extract sample names and sequences
for record in records:
    sample_name = record.id
    sequence = str(record.seq)
    data['File'].append(sample_name)
    sequences.append(list(sequence))

# Prepare transposed_sequences
max_length = max(len(seq) for seq in sequences)
padded_sequences = [seq + ['-'] * (max_length - len(seq)) for seq in sequences]
transposed_sequences = list(map(list, zip(*padded_sequences)))

# Add transposed sequences to the data dictionary
for i in range(max_length):
    column_name = f'Amino Acid {i+1}'
    data[column_name] = transposed_sequences[i]

# Create a DataFrame from the data dictionary
df = pd.DataFrame(data)

df.to_excel(excel_file, index=False)

print("Data successfully processed and saved to Excel.")
