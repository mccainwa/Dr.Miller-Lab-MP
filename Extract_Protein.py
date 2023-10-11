#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 13:15:04 2023

@author: blim
"""

import os

parent_directory = "/home/blim/Downloads"
output_directory = "/home/blim/Downloads"

# Check if output directory exists, create if not
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
    
    # List of protein identifiers
proteins = ["L1 protein"] # Loop through each file in the directory

for filename in os.listdir(parent_directory):
    if filename.endswith('.txt'):  
        # Check if file ends with .txt
        print(f"Processing file {filename}")
        filepath = os.path.join(parent_directory, filename)
        with open(filepath, 'r') as file:
            content = file.readlines()          
            
        for protein in proteins:
            output_content = []          
            for i, line in enumerate(content):
                if line.startswith('>') and protein in line:
                    output_content.append(line.strip())  # Add the header line
                    if i + 1 < len(content):  # If there is a next line
                        next_line = content[i + 1].strip()
                        while next_line and not next_line.startswith('>'):
                            output_content.append(next_line)
                            i += 1
                            if i + 1 < len(content):
                                next_line = content[i + 1].strip()
                            else:
                                break
                elif line.startswith('>'):
                    break
                else:
                    output_content.append(line.strip())  # Add the sequence line
            
            if output_content:  # If there's any sequence to write
                # Create new file name and path
                new_filename = os.path.splitext(filename)[0] + "_" + protein + "_processed.txt"
                new_filepath = os.path.join(output_directory, new_filename)                
                
                # Write sequences to new file
                with open(new_filepath, 'w') as new_file:
                    new_file.write('\n'.join(output_content))        
                    
            print(f"Finished processing file {filename}")
        print("Finished processing all files.")