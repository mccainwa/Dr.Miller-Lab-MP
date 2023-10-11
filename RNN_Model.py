#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 17:05:51 2023

@author: blim
"""
#%%
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Model, Sequential
from tensorflow.keras.layers import Input, LSTM, Dense

# Load your amino acid sequence data into a DataFrame
df = pd.read_excel('PCA_Proof.xlsx')

# Remove the last underscore and anything that follows it in the column values
df['File'] = df['File'].str.replace(r'_\d+$', '', regex=True)

# Get the amino acid sequences as input data
sequences = df.iloc[:, 1:].astype(str).values

# Create a mapping from amino acid to integer
amino_acids = np.unique(sequences)
amino_acid_to_int = {aa: i + 1 for i, aa in enumerate(amino_acids)}

# Convert the amino acid sequences to one-hot encoded sequences
input_dim = sequences.shape[1]
num_amino_acids = len(amino_acids) + 1
encoded_sequences = np.zeros((len(sequences), input_dim, num_amino_acids), dtype=int)
for i, seq in enumerate(sequences):
    for j, aa in enumerate(seq):
        encoded_sequences[i, j, amino_acid_to_int[aa]] = 1
#%%
# Define the dimensions for input and output sequences
input_dim = 194
output_dim = 194

# Build the RNN model
model = Sequential()
model.add(LSTM(128, input_shape=(None, input_dim)))
model.add(Dense(output_dim, activation='softmax'))

# Compile the model
model.compile(optimizer='adam', loss='categorical_crossentropy')

# Train the model
model.fit(X_train, y_train, epochs=10, batch_size=128, shuffle=True)

# Predict the output sequences
predicted_sequences = model.predict(X_test)