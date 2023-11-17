#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 10:46:49 2023

@author: ebroni
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import make_scorer, accuracy_score
from tensorflow.keras.models import Model, load_model
from keras.wrappers.scikit_learn import KerasClassifier, KerasRegressor
from tensorflow.keras.layers import Input, Dense, Flatten, Reshape
import matplotlib.pyplot as plt

# Load Data and One-Hot Encode
df = pd.read_excel('DENV1.xlsx')
df['File'] = df['File'].str.replace(r'_\d+$', '', regex=True)
sequences = df.iloc[:, 1:].astype(str).values
amino_acids = np.unique(sequences)
amino_acid_to_int = {aa: i + 1 for i, aa in enumerate(amino_acids)}
input_dim = sequences.shape[1]
num_amino_acids = len(amino_acids) + 1
encoded_sequences = np.zeros((len(sequences), input_dim, num_amino_acids), dtype=int)
for i, seq in enumerate(sequences):
    for j, aa in enumerate(seq):
        encoded_sequences[i, j, amino_acid_to_int[aa]] = 1

# Split Data
X_train, X_test = train_test_split(encoded_sequences, test_size=0.05, random_state=42)

# Define the hyperparameters and their values to search over
param_grid = {
    'latent_dim': [3394, 70, 70],
    'epochs': [10, 20, 30],
    'batch_size': [32, 64, 128]
}

def build_autoencoder(latent_dim):
    input_layer = Input(shape=(input_dim, num_amino_acids))
    flatten_layer = Flatten()(input_layer)
    encoded = Dense(latent_dim, activation='relu')(flatten_layer)
    decoded = Dense(input_dim * num_amino_acids, activation='relu')(encoded)
    reshaped_output = Reshape((input_dim, num_amino_acids))(decoded)
    autoencoder = Model(inputs=input_layer, outputs=reshaped_output)        # Change the activation function of the output layer to 'softmax'
    output_layer = Dense(input_dim * num_amino_acids, activation='softmax')(encoded)
    reshaped_output = Reshape((input_dim, num_amino_acids))(output_layer)
    autoencoder = Model(inputs=input_layer, outputs=reshaped_output)
    autoencoder.compile(optimizer='adam', loss='mse')
    return autoencoder

# Create a KerasRegressor based on the build_autoencoder function
model = KerasRegressor(build_fn=build_autoencoder)

# Use GridSearchCV to search over the hyperparameter grid
grid_search = GridSearchCV(estimator=model, param_grid=param_grid,
                           scoring='neg_mean_squared_error', cv=3)

# Fit the GridSearchCV on the training data
grid_result = grid_search.fit(X_train, X_train)

# Print the best hyperparameters and associated reconstruction error
print("Best Hyperparameters: ", grid_result.best_params_)
print("Best Reconstruction Error: ", -grid_result.best_score_)

# Retrieve the best model
best_model = grid_result.best_estimator_.model

# Train the best model with the entire training data
history = best_model.fit(X_train, X_train, epochs=grid_result.best_params_['epochs'],
                         batch_size=grid_result.best_params_['batch_size'], shuffle=True)

# Save the trained model
model_path = 'autoencoder_model.h5'
best_model.save(model_path)
print("Trained model saved successfully.")

# Load the saved model
loaded_model = load_model(model_path)
print("Trained model loaded successfully.")

# Deciphering my Model
encoders = Model(inputs=loaded_model.input, outputs=loaded_model.layers[-4].output)
latent_vectors = encoders.predict(encoded_sequences)

# Predict Future Sequences and Normalize
predicted_sequences = loaded_model.predict(X_test)
predicted_sequences = predicted_sequences / np.max(predicted_sequences, axis=(0, 2), keepdims=True)

# Decode the predicted sequences back to amino acids
decoded_sequences = np.argmax(predicted_sequences, axis=2)

int_to_amino_acid = {i + 1: aa for i, aa in enumerate(amino_acids)}

decoded_sequences_aa = np.empty_like(decoded_sequences, dtype=object)
for i in np.ndindex(decoded_sequences.shape):
    if decoded_sequences[i] != 0:
        decoded_sequences_aa[i] = int_to_amino_acid[decoded_sequences[i]]
    else:
        decoded_sequences_aa[i] = 'NaN'

# Calculating the Mean per AA per position, and gathering the top 2 most probable AAs
mean_predicted_sequences = np.zeros((70,18))

for i in range(predicted_sequences.shape[0]):
    for j in range(70):
        column_values = predicted_sequences[i, j, :]
        mean_predicted_sequences[j, :] += column_values

mean_predicted_sequences /= predicted_sequences.shape[0]
mean_predicted_sequences = mean_predicted_sequences.T

top_two_indices = np.argsort(-mean_predicted_sequences, axis=0)[:2, :]
top_two_values = mean_predicted_sequences[top_two_indices, np.arange(70)]

diff = top_two_values[0, :] - top_two_values[1, :]
mutated_aa = decoded_sequences_aa[0].copy()
min_diff_positions = np.argsort(diff)[:3]

for pos in min_diff_positions:
    aa_index_1 = top_two_indices[0, pos]
    aa_index_2 = top_two_indices[1, pos]
    aa_label_1 = amino_acids[aa_index_1 - 1] if aa_index_1 != 0 else 'NaN'
    aa_label_2 = amino_acids[aa_index_2 - 1] if aa_index_2 != 0 else 'NaN'
    probability = 1 - diff[pos]
    mutated_aa[pos] = aa_label_2
    print(f"Position {pos}: Probability = {probability}, Mutations: {aa_label_1} -> {aa_label_2}")

mutated_sequence_string = "".join(mutated_aa)
decoded_sequence_string = "".join(decoded_sequences_aa[2])

print("Decoded Sequence with Mutations:")
print(mutated_sequence_string)
print("Decoded Sequence without Mutations:")
print(decoded_sequence_string)

# Analysis comparison of autoencoder to input data
reconstructed_sequences = X_test
decoded_reconstructed_sequences = np.argmax(reconstructed_sequences, axis=2)

decoded_reconstructed_sequences_aa = np.empty_like(decoded_reconstructed_sequences, dtype=object)
for i in np.ndindex(decoded_reconstructed_sequences.shape):
    if decoded_reconstructed_sequences[i] != 0:
        decoded_reconstructed_sequences_aa[i] = int_to_amino_acid[decoded_reconstructed_sequences[i]]
    else:
        decoded_reconstructed_sequences_aa[i] = 'NaN'

num_diff = 0
for i in range(decoded_reconstructed_sequences_aa.shape[0]):
    for j in range(decoded_reconstructed_sequences_aa.shape[1]):
        if decoded_reconstructed_sequences_aa[i, j] != decoded_sequences_aa[i, j]:
            num_diff += 1

total_elements = decoded_sequences_aa.size
percent_diff = 1 - (num_diff / total_elements)

num_identical_samples = np.sum(np.all(sequences == np.tile(decoded_sequences_aa[2, :], (sequences.shape[0], 1)), axis=1))
num_identical_mutants = np.sum(np.all(mutated_aa == sequences, axis=1))

print("Number of Differences:", num_diff)
print("Percent Accuracy:", percent_diff)
