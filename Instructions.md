# INSTRUCTIONS.md  
## AEGIS (AutoEncoder-driven Genomic Insight System) — Instructions

This guide is written to replicate the core results of the AEGIS project. By the end, you will be able to:

- Load aligned amino acid sequence data
- Clean and encode sequences for machine learning
- Train an autoencoder on viral protein sequence variation
- Use Monte Carlo dropout to measure uncertainty at each amino acid position
- Rank predicted “mutation hotspot” positions
- Generate tables and figures that summarize results

> **Important:** This project is not predicting the future with certainty. It identifies positions that look “mutation-prone” because the model is uncertain and multiple amino acids remain plausible at that site.

---

## 0) What you need before starting

### Software
- Python **3.9+** recommended
- I ran all of this code on Spyder

### Data requirement
You will need an **aligned amino acid dataset** in Excel format, similar to `All_RBD.xlsx`.

The script assumes:
- Column 0: sequence names/IDs
- Columns 1…N: amino acids at each aligned position
- Each row is one sequence
- All sequences have the same length (alignment)

There is a script in the files that can take NCBI genes and filter out extraneous data and align the amino acid dataset
---

## 1) Clone the repository

Option A: Git (recommended)
```bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>
