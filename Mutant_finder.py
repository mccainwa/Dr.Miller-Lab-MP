#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 18:38:50 2023

@author: blim
"""

import pandas as pd
import re
import pandas as pd
import re

# Read the Excel file
input_file = input("Enter the input Excel file path: ")
df = pd.read_excel(input_file)

# Create a new DataFrame to store the different values
different_values_df = pd.DataFrame(columns=['File', 'Position', 'Mode', 'Different Value', 'Merged'])

# Iterate over each column except the first column
for column in df.columns[1:]:
    # Get the most common value in the column
    mode_value = df[column].mode().values[0]

    # Check if the column is the second column
    if column == 'Column':
        # Remove the word 'Letter' from the values in column 2
        values = df[column].str.replace('Letter', '')
    else:
        values = df[column]

    # Extract the numeric part from the column name using regular expression
    position = re.search(r'\d+', column).group()

    # Create a list to store the new entries for this column
    new_entries = []

    # Iterate over each value in the column, starting from the second row
    for index, value in values.items():
        # Check if the value is different from the mode value
        if value != mode_value:
            # Create the merged value by combining position, mode, and different value
            merged_position = int(position)
            merged = f"{merged_position}{mode_value}{value}"

            # Create a new entry with the file name, position, mode, different value, and merged column
            new_entry = {
                'File': df.at[index, 'File'],
                'Position': position,
                'Mode': mode_value,
                'Different Value': value,
                'Merged': merged
            }

            # Append the new entry to the list
            new_entries.append(new_entry)

    # Append the new entries for this column to the different_values_df
    different_values_df = pd.concat([different_values_df, pd.DataFrame(new_entries)], ignore_index=True)

# Ask for the output Excel file name
output_file_name = input("Enter the output Excel file name: ")

# Save the different values to a new Excel file
different_values_df.to_excel(output_file_name, index=False)


