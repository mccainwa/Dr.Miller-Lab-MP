#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 10:22:45 2023

@author: blim
"""

import glob
import os

input_files = glob.glob('./*.txt')
output_suffix = '_'

for input_file in input_files:
    with open(input_file, 'r') as file:
        contents = file.read().strip().split('\n\n')    
        
    # Extract the input file name without the extension
    file_name = os.path.splitext(os.path.basename(input_file))[0]    
    
    for i, content in enumerate(contents):
        output_file = f'{file_name}{output_suffix}{i+1}.txt'
        with open(output_file, 'w') as file:
            file.write(content)
