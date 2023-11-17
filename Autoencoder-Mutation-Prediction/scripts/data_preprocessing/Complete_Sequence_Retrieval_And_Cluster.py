#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 12:40:44 2023

@author: blim


"""

import os
import pandas as pd


#User Prompts
parent_directory = input("Enter the parent directory: ")
start_phrase = input("Enter start phrase: ")
end_phrase = input("Enter end phrase: ")
output_directory = input("Enter the desired directory for the output file: ")


#Iterate over subfolders in the parent directory
for subfolder in os.listdir(parent_directory):
    print(subfolder)
    subfolder_path = os.path.join(parent_directory, subfolder)
    
    #Skip if the item in the parent directory is not a subfolder
    if not os.path.isdir(subfolder_path):
        continue

    # Initialize a list to store the results
    results = []
    
    # Iterate over files in the directory
    for filename in os.listdir(subfolder_path):
        if filename.endswith('.txt'):  # Process only text files
            filepath = os.path.join(subfolder_path, filename)
            with open(filepath, 'r') as file:
                # Read the content of the file and remove newline characters
                content = file.read().replace('\n', '')
                
                # Find the start and end indices
                start_index = content.find(start_phrase)
                end_index = content.find(end_phrase)
                
                # Check if both start and end phrases are present
                if start_index != -1 and end_index != -1:
                    # Extract the desired string
                    extracted_string = content[start_index:end_index+len(end_phrase)]
                    
                    # Check if the extracted string contains 'X'
                    if 'X' not in extracted_string and 'J' not in extracted_string:
                        # Remove spaces from the extracted string
                        extracted_string = extracted_string.replace(" ", "")
                        
                        # Split the extracted string into individual letters
                        letters = list(extracted_string)
                        
                        # Store the results
                        result = {
                            'File': os.path.splitext(filename)[0],  # Remove the file extension                     
                            **{f'Letter {i+1}': letter for i, letter in enumerate(letters)}
                        }
                        results.append(result)
    
# Create a DataFrame from the results
df = pd.DataFrame(results)
    
# Save the DataFrame to an Excel file
output_file = f'{subfolder}.xlsx'
output_path = os.path.join(output_directory, output_file)
df.to_excel(output_file, index=False)
print(output_file)
print(output_path)
    
# Read the Excel file
df = pd.read_excel(output_file)

# Concatenate columns for each row into a single string
rows_concatenated = []
for row in df.values:
    header = ">{}".format(row[0])
    sequence = ''.join(map(str, row[1:]))
    row_string = "{}\n{}".format(header, sequence)
    rows_concatenated.append(row_string)
# Output the concatenated strings to a text file
with open(f'{subfolder}.txt', 'w') as file:
    file.write('\n'.join(rows_concatenated))
