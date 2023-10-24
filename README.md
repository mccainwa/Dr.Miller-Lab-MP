# Dr.Miller-Lab-MP

---

# Viral Amino Acid Mutation Prediction Toolkit

## Overview

The Viral Amino Acid Mutation Prediction Toolkit is a comprehensive software package designed to facilitate the extraction, preprocessing, training, and mutation prediction from viral genomic data. It harnesses the power of unsupervised machine learning, including variational autoencoders (VAE), Recurrent Neural Networks (RNN), and Long Short-Term Memory Neural Networks (LSTM), to uncover hidden mutational patterns within viral genomes. This toolkit focuses primarily on Sars-Cov-2 but can be adapted for other viruses.

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

## Contributing

We welcome contributions to improve this toolkit! If you'd like to enhance its functionality, feel free to create a pull request.

