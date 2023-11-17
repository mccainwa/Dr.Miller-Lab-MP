#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 12:33:34 2023

@author: blim
"""

import glob
import os

# Directory containing the .txt files
directory = '/home/blim/Documents/Summer2023/DENV/DENV4'

# List all .txt files in the directory
file_list = glob.glob(directory + '/*.txt')

# Define a custom sorting key function to extract the numeric portion of the file names
def get_file_number(file_path):
    file_name = os.path.splitext(os.path.basename(file_path))[0]
    return int(file_name[4:])

# Sort the file list based on the numeric portion of the file names
file_list = sorted(file_list, key=get_file_number)

new_lines = []
for file_path in file_list:
    file_name = os.path.splitext(os.path.basename(file_path))[0]  # Get the file name without extension
    with open(file_path, 'r') as file:
        lines = file.readlines()

    new_lines.append(f'>{file_name}\n')  # Add the file name as the identifier

    for line in lines:
        if not line.startswith('>'):
            new_line = line.replace('nan', '')
            new_lines.append(new_line)

    new_lines.append('\n\n')  # Add an extra newline before moving to the next file

# Write the modified lines to a new file
with open('merged_output.fasta', 'w') as file:
    file.writelines(new_lines)

