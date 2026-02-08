# Dr.Miller-Lab-MP

---
## Overview

**AEGIS (AutoEncoder-driven Genomic Insight System)** is an unsupervised machine learning framework designed to identify mutation-prone amino acid positions in viral proteins, with a primary focus on the SARS-CoV-2 Spike Receptor Binding Domain (RBD).

AEGIS leverages a deep autoencoder trained on aligned viral protein sequences and applies Monte Carlo dropout at inference time to quantify epistemic uncertainty in amino acid predictions. Positions exhibiting high predictive ambiguity—characterized by elevated entropy and small probability margins between competing residues—are interpreted as being under increased evolutionary pressure.

To further ground predictions in biological realism, AEGIS incorporates position-specific Markov transition modeling to capture observed amino acid substitution dynamics. Probabilistic outputs from the autoencoder and Markov model are fused to produce a ranked list of mutation “hotspots” that reflect both statistical uncertainty and evolutionary plausibility.

Rather than predicting a single deterministic mutation, AEGIS highlights regions of the protein landscape most susceptible to change, making it well-suited for viral surveillance, hypothesis generation, and downstream experimental prioritization.

## High-Level Pipeline

1. Sequence ingestion and cleaning  
2. One-hot encoding of aligned amino acid sequences  
3. Autoencoder training  
4. Monte Carlo dropout inference  
5. Entropy and epistemic uncertainty estimation  
6. Markov transition modeling  
7. Probability fusion and hotspot ranking  
8. Visualization and export of results  

---

## Functionality

Analyzing DNA/RNA viral amino acid sequences using Machine Learning algorithms to predict high mutability or resistibility amino acid positions. 

A pipeline built in Python to query NCBI with a user-input search information and statistically process the retrieved protein FASTA files. The pipeline provides a user interface to NCBI to pull FASTA sequences by asking the user for what they want to pull. The resulting FASTA sequences are stored in an input file, which will then be pre-processed to align the data, remove any sequences that contain ambiguous amino acids, and output a FASTA and Excel file that contains user-specified protein or amino acid sequences. The Excel file will be used as the input into our 3 machine learning models, which will be trained on the data where it can detect hidden patterns. The Mutation prediction function will output the top 3 most probable mutations with 10 iterations, for a total of 30 predictions. Multiple trials can be conducted.

When running the Python code, be sure to specify the desired path on your device for output_file.

## Features

- **Data Retrieval**: Easily pull genomic data from the National Center for Biotechnology Information (NCBI) to create a robust dataset for analysis.

- **Preprocessing**: Streamline your data by filtering out ambiguous amino acids, ensuring a clean dataset for model training.

- **Machine Learning Models**: Train and fine-tune your predictive models using the provided VAE, RNN, LSTM, optimizing with a mean squared error (MSE) loss function. ML algorithm is equipped to tune hyperparameters such as epoch size, batch size, and latent dimensions.

- **Mutation Prediction**: Employ the toolkit to predict mutations in viral genomes based on hidden patterns, ranking the most probable mutations.

## Installation (I got this from our friend chatgpt so could be useful to add but if you could edit this part that would be great)

To get started, clone the repository and install the required dependencies:

```bash
git clone https://github.com/yourrepository/viral-mutation-toolkit.git
cd viral-mutation-toolkit
pip install -r requirements.txt
```

## Usage

Follow these steps to harness the full power of the Viral Mutation Prediction Toolkit:

1. **Data Retrieval**: Use the provided scripts to fetch data from NCBI or integrate your dataset.

2. **Data Preprocessing**: Utilize the data preprocessing functions to clean and format your dataset for training.

3. **Training the Model**: Train your predictive model using the built-in variational autoencoder (VAE). You can fine-tune it to optimize performance.

4. **Mutation Prediction**: Predict viral mutations based on the trained model's insights into hidden mutational patterns.

## Example

```python
# Import the toolkit functions
from viral_mutation_toolkit import data_retrieval, preprocessing, train_model, predict_mutations

# Step 1: Data Retrieval
data = data_retrieval.fetch_data_from_ncbi("Sars-Cov-2")

# Step 2: Data Preprocessing
clean_data = preprocessing.clean_and_format_data(data)

# Step 3: Training the Model
trained_model = train_model.train_vae_model(clean_data)

# Step 4: Mutation Prediction
predictions = predict_mutations.predict_mutations(trained_model, clean_data)

# Display the top predicted mutations
print(predictions)
```



## Data Input

- Input data is provided as an aligned Excel file:
  - Column 0: sequence identifier
  - Columns 1…N: amino acids per position
- Rows containing ambiguous amino acids (`X`, `J`) are removed
- Gap characters (`-`) are replaced with the column-wise modal amino acid
- The reference implementation was developed using:
  - `All_RBD.xlsx` (SARS-CoV-2 Spike RBD)

---

## Model Architecture

### Autoencoder

- Input: one-hot encoded tensor of shape  
  `(sequence_length × number_of_amino_acids)`
- Encoder:
  - Flatten → Dense → Dropout (MC-enabled)
- Decoder:
  - Dense → Reshape → Softmax per position
- Loss function: Mean Squared Error (MSE)
- Optimizer: Adam

Dropout layers remain active during inference to enable Monte Carlo sampling.

---

## Monte Carlo Dropout and Uncertainty Estimation

- The trained autoencoder is evaluated across multiple stochastic forward passes (default: 200)
- Each pass yields a probability distribution over amino acids per position
- Aggregated metrics include:
  - Mean probability distribution
  - Epistemic variance
  - Shannon entropy per residue

These metrics identify positions where the model exhibits consistent uncertainty, a signal associated with mutational flexibility.

---

## Markov Transition Modeling

To incorporate evolutionary context:

- A position-specific Markov transition matrix is constructed from observed amino acid transitions
- Laplace smoothing is applied to stabilize estimates
- The stationary distribution of each matrix represents evolutionary preference
- Monte Carlo probabilities are fused with Markov priors:

\[
p_{\text{fused}} \propto p_{\text{MC}}^{\lambda} \cdot p_{\text{Markov}}^{(1-\lambda)}
\]

This fusion discourages biologically implausible mutation predictions.

---

## Mutation Hotspot Scoring

Each amino acid position is ranked using a composite score derived from:

- Low Top-1 vs Top-2 probability margin
- High Shannon entropy
- (Optionally) high epistemic variance

Higher scores indicate statistically ambiguous and evolutionarily permissive sites.

---

## Outputs

Key outputs include:

- `top_mutations_mc.xlsx` — most ambiguous mutation candidates
- `aegis_mc_markov_fused_hotspots.xlsx` — final ranked mutation hotspots
- `figures/` directory containing:
  - Entropy and epistemic variance plots
  - Mutation probability heatmaps
  - Hotspot rankings and confidence visualizations

---

## Installation

```bash
git clone https://github.com/<your-username>/Dr.Miller-Lab-MP.git
cd Dr.Miller-Lab-MP
pip install -r requirements.txt

**Running the Model**

1. Place your aligned amino acid Excel file in the project directory

2. Update file paths inside the script if necessary

3. Run: python final_Sars_MC_AE.py

Outputs will be written to the working directory.

**Contributing**

We welcome contributions to improve this toolkit! If you'd like to enhance its functionality, feel free to create a pull request.

**Disclaimer**

This software is intended for research and exploratory use only.
It does not provide clinical, diagnostic, or therapeutic predictions.

**Acknowledgments**

Developed in the Miller Lab.
Code authored by Brandon Lim, Emmanuel Broni, and Whelton Miller, PhD.

