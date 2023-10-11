#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 15:00:47 2023

@author: blim
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Flatten, Reshape, LSTM
import matplotlib.pyplot as plt
from scipy.stats import mode

#%% Load Data and One-Hot Encode
# Load your amino acid sequence data into a DataFrame
df = pd.read_excel('All_RBD.xlsx')

# Check if there are underscores in the 'File' column
if any(df['File'].str.contains(r'_')):
    # Remove the last underscore and anything that follows it in the 'File' column values
    df['File'] = df['File'].str.replace(r'_\d+$', '', regex=True)

# Get the amino acid sequences as input data
sequences = df.iloc[:, 1:].astype(str).values

# Remove rows containing 'X' or 'J' in the sequences
mask = (~(sequences == 'X')).all(axis=1) & (~(sequences == 'J')).all(axis=1)
sequences = sequences[mask]

# Check for "-" values in sequences and replace them with the column mode
column_modes = mode(sequences, axis=0).mode
sequences = np.where(sequences == '-', np.tile(column_modes, (sequences.shape[0], 1)), sequences)

''' 
#add if you want to save specific sequence
df1 = pd.DataFrame(sequences)
df1.to_excel(r'/home/blim/Documents/Summer2023/DENV/HPV_L1_mode.xlsx')
'''

# Create a mapping from amino acid to integer
amino_acids = np.unique(sequences)
amino_acid_to_int = {aa: i for i, aa in enumerate(amino_acids)}

# Convert the amino acid sequences to one-hot encoded sequences
input_dim = sequences.shape[1]
num_amino_acids = len(amino_acids)
encoded_sequences = np.zeros((len(sequences), input_dim, num_amino_acids), dtype=int)
for i, seq in enumerate(sequences):
    for j, aa in enumerate(seq):
        encoded_sequences[i, j, amino_acid_to_int[aa]] = 1
       
#%% Split Data

# Split the data into train and test sets
X_train, X_test = train_test_split(encoded_sequences, test_size=0.25, random_state=42)

#%% Build Model
# Define the dimensions for input and latent space
latent_dim = sequences.shape[1] # Set to the same value as input_dim

# Build the autoencoder model
input_layer = Input(shape=(input_dim, num_amino_acids))
flatten_layer = Flatten()(input_layer)
encoded = Dense(latent_dim, activation='relu')(flatten_layer)
decoded = Dense(input_dim * num_amino_acids, activation='relu')(encoded)
reshaped_output = Reshape((input_dim, num_amino_acids))(decoded)
autoencoder = Model(inputs=input_layer, outputs=reshaped_output)

# Change the activation function of the output layer to 'softmax'
output_layer = Dense(input_dim * num_amino_acids, activation='softmax')(encoded)
reshaped_output = Reshape((input_dim, num_amino_acids))(output_layer)
autoencoder = Model(inputs=input_layer, outputs=reshaped_output)

# Compile and train the autoencoder
autoencoder.compile(optimizer='adam', loss='mse')
history = autoencoder.fit(X_train, X_train, epochs=20, batch_size=64, shuffle=True)

#%% Deciphering my Model

# Obtain the latent representations
encoders = Model(inputs=input_layer, outputs=encoded)
latent_vectors = encoders.predict(encoded_sequences)

#%% Predict Future Sequences and normalize
   
# Predict future sequences using the trained autoencoder
predicted_sequences = autoencoder.predict(X_test)
predicted_sequences = predicted_sequences/np.max(predicted_sequences, axis=(0,2), keepdims=True)

# Decode the predicted sequences back to amino acids
decoded_sequences = np.argmax(predicted_sequences, axis=2)

int_to_amino_acid = {i: aa for i, aa in enumerate(amino_acids)}

decoded_sequences_aa = np.empty_like(decoded_sequences, dtype=object)
for i in np.ndindex(decoded_sequences.shape):
    decoded_sequences_aa[i] = int_to_amino_acid[decoded_sequences[i]]

#%% Calculating the Mean per AA per position, and gathering the top 2 most probable AA's    
 
# Initialize an array to store the mean values
mean_predicted_sequences = np.zeros((sequences.shape[1], encoded_sequences.shape[2]))

# Calculate the mean for each column per index at axis 2 for each sample
for i in range(predicted_sequences.shape[0]):  # Iterate over samples
    for j in range(predicted_sequences.shape[1]):  # Iterate over positions
        column_values = predicted_sequences[i, j, :]
        mean_predicted_sequences[j, :] += column_values

# Divide by the number of samples to get the mean
mean_predicted_sequences /= predicted_sequences.shape[0]
mean_predicted_sequences = mean_predicted_sequences.T

# Find the indices of the top two maximum values in axis=1 for each column
top_two_indices = np.argsort(-mean_predicted_sequences, axis=0)[:2, :]

# Calculate the difference between the top two values in each column
diff = mean_predicted_sequences[top_two_indices[0], np.arange(sequences.shape[1])] - mean_predicted_sequences[top_two_indices[1], np.arange(sequences.shape[1])]

# Find the positions with the smallest difference
min_diff_positions = np.argsort(diff)[:10]

# Create an empty DataFrame to store the top mutations
mutation_df = pd.DataFrame(columns=['Position', 'Global Position', 'Mutation', 'Confidence Score'])

mutated_aa = decoded_sequences_aa[0].copy()
for pos in min_diff_positions:
    aa_index_1 = top_two_indices[0, pos]
    aa_index_2 = top_two_indices[1, pos]
    aa_label_1 = amino_acids[aa_index_1]
    aa_label_2 = amino_acids[aa_index_2] 
    Confidence_score = 1 - diff[pos]
    new_entry = {'Position': pos, 'Global Position': pos + 331, 'Mutation': f"{aa_label_1} -> {aa_label_2}", 'Confidence Score': Confidence_score}
    mutation_df = pd.concat([mutation_df, pd.DataFrame(new_entry, index=[0])], ignore_index=True)
    
mutated_sequence_string = "".join(mutated_aa)
decoded_sequence_string = "".join(decoded_sequences_aa[0])

# Print the top mutations
print(mutation_df)

# Save the mutation DataFrame to Excel
mutation_df.to_excel('top_mutations.xlsx', index=False)

#%% Analysis comparison of autoencoder to input data

# Reconstruct the encoded sequences
reconstructed_sequences = X_test

# Decode the predicted sequences back to amino acids
decoded_reconstructed_sequences = np.argmax(reconstructed_sequences, axis=2)

decoded_reconstructed_sequences_aa = np.empty_like(decoded_reconstructed_sequences, dtype=object)
for i in np.ndindex(decoded_reconstructed_sequences.shape):
    decoded_reconstructed_sequences_aa[i] = int_to_amino_acid[decoded_reconstructed_sequences[i]]

# Iterate over each element in the matrices and count differences
num_diff = 0
for i in range(decoded_reconstructed_sequences_aa.shape[0]):
    for j in range(decoded_reconstructed_sequences_aa.shape[1]):
        if decoded_reconstructed_sequences_aa[i, j] != decoded_sequences_aa[i, j]:
            num_diff += 1

# Stats            
total_elements = decoded_sequences_aa.size
percent_diff = 1 - (num_diff / total_elements)

# Count the number of samples and mutants that are identical to input sequences
num_identical_samples = np.sum(np.all(sequences == np.tile(decoded_sequences_aa[0, :], (sequences.shape[0], 1)), axis=1))
num_identical_mutants = np.sum(np.all(mutated_aa == sequences, axis=1))

print("Number of Differences:", num_diff)
print("Percent Accuracy:", percent_diff)
print(f"Number of identical samples: {num_identical_samples}")
print(f"Number of identical mutants: {num_identical_mutants}")

#%% Graphs

import seaborn as sns

# Calculate the probabilities for each amino acid at each position
amino_acid_probabilities = np.mean(predicted_sequences, axis=0)

# Create a DataFrame to store the probabilities
amino_acid_prob_df = pd.DataFrame(amino_acid_probabilities, columns=amino_acids).T

# Plot the heatmap
plt.figure(figsize=(12, 6))
ax = sns.heatmap(amino_acid_prob_df, cmap='viridis', xticklabels=5, yticklabels=1)
ax.tick_params(right=True, top=True, labelright=True, labeltop=True, rotation=0)
plt.xlabel('Position Number')
plt.ylabel('Amino Acid')
plt.title('Autoencoder Amino Acid Confidence Score Heatmap')
plt.show()
