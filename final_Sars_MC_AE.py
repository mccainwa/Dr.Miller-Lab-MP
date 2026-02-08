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
from tensorflow.keras.layers import Input, Dense, Flatten, Reshape, Dropout
import tensorflow as tf
import matplotlib.pyplot as plt
from scipy.stats import mode
from sklearn.model_selection import train_test_split

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

#%% Build Model (with dropout for MC sampling)
latent_dim = sequences.shape[1]  # same as input_dim

input_layer = Input(shape=(input_dim, num_amino_acids))
x = Flatten()(input_layer)
x = Dense(latent_dim, activation='relu')(x)
x = Dropout(0.30, name="mc_dropout")(x)   # <-- NEW: stochastic latent for MC

encoded = x # alias for legacy code below that uses encoded

output = Dense(input_dim * num_amino_acids, activation='softmax')(x)
output = Reshape((input_dim, num_amino_acids))(output)

autoencoder = Model(inputs=input_layer, outputs=output)
autoencoder.compile(optimizer='adam', loss='mse')
history = autoencoder.fit(X_train, X_train, epochs=20, batch_size=64, shuffle=True, verbose=1)

# optional code but I kept it for latent vectors later
encoders = Model(inputs=input_layer, outputs=x)
latent_vectors = encoders.predict(encoded_sequences, verbose=0)

#%% Monte Carlo Dropout Inference (AEGIS)
T = 200                  # increase if you have time/compute (e.g., 500–1000)
eps = 1e-12

mc_pass_means = []       # each item: (L, V) mean over test samples for that pass
for t in range(T):
    preds = autoencoder(X_test, training=True).numpy()     # keep dropout ON
    preds = preds / np.clip(preds.sum(axis=2, keepdims=True), 1e-12, None)
    pass_mean = preds.mean(axis=0)                         # (L, V)
    mc_pass_means.append(pass_mean)

mc_pass_means = np.stack(mc_pass_means, axis=0)            # (T, L, V)

# Aggregate across MC passes
mean_prob_mc = mc_pass_means.mean(axis=0)                  # (L, V)
var_prob_mc  = mc_pass_means.var(axis=0)                   # (L, V)
epi_var_pos  = var_prob_mc.mean(axis=1)                    # (L,) epistemic variance per position
entropy_pos  = -(mean_prob_mc * np.log(np.clip(mean_prob_mc, eps, None))).sum(axis=1)  # (L,)

#%% Markov transitions + fusion with MC mean probabilities
ALPHA = 0.5     # Laplace/Dirichlet smoothing for transitions. can smooth from 0.1-1.0
LAMBDA = 0.7    # blend weight: fused ∝ p_MC^λ * p_Markov^(1-λ). 0.6-0.8 keeps MC dominant but lets Markov help
eps = 1e-12

V = len(amino_acids)
L = sequences.shape[1]
aa_to_idx = {aa:i for i, aa in enumerate(amino_acids)}
idx_to_aa = {i:aa for aa, i in aa_to_idx.items()}
RBD_OFFSET = 331

# Per-position transition matrices Tmat[j] shape (V,V): P(next=b | current=a) at position j
Tmat = np.zeros((L, V, V), dtype=float)
for j in range(L):
    prev = sequences[:-1, j]
    nxt  = sequences[1:,  j]
    for a, b in zip(prev, nxt):
        ia, ib = aa_to_idx.get(a), aa_to_idx.get(b)
        if ia is not None and ib is not None:
            Tmat[j, ia, ib] += 1.0
    # smooth and row-normalize
    Tmat[j] = Tmat[j] + ALPHA
    Tmat[j] = Tmat[j] / (Tmat[j].sum(axis=1, keepdims=True) + eps)

def stationary(P, iters=2000, tol=1e-12):
    # left stationary distribution via power iteration
    x = np.ones(P.shape[0], dtype=float) / P.shape[0]
    for _ in range(iters):
        x_new = x @ P
        if np.linalg.norm(x_new - x, 1) < tol:
            break
        x = x_new
    return x / (x.sum() + eps)

# Stationary distribution per position (L,V)
pi_markov = np.stack([stationary(Tmat[j]) for j in range(L)], axis=0)

# Fuse with MC mean distribution computed above
p_mc = mean_prob_mc / (mean_prob_mc.sum(axis=1, keepdims=True) + eps)   # (L,V)
p_mk = pi_markov  / (pi_markov.sum(axis=1, keepdims=True) + eps)        # (L,V)
p_fused = np.power(p_mc + eps, LAMBDA) * np.power(p_mk + eps, 1.0 - LAMBDA)
p_fused = p_fused / (p_fused.sum(axis=1, keepdims=True) + eps)

# Rank fused distribution by small margin + high entropy (and you can add epi_var_pos if desired)
order = np.argsort(-p_fused, axis=1)
top1_f = order[:, 0]; top2_f = order[:, 1]
p1f = p_fused[np.arange(L), top1_f]; p2f = p_fused[np.arange(L), top2_f]
margin_f = p1f - p2f
entropy_f = -(p_fused * np.log(np.clip(p_fused, eps, None))).sum(axis=1)

# Normalize & score
def minmax(x):
    lo, hi = float(np.min(x)), float(np.max(x))
    return np.zeros_like(x) if hi <= lo + 1e-12 else (x - lo) / (hi - lo)

score_fused = (1.0 - minmax(margin_f)) * (entropy_f / np.log(V))  # add * minmax(epi_var_pos) if you want

# Save top-K fused hotspots
TOPK = 20
rank_idx = np.argsort(-score_fused)[:TOPK]

fused_rows = []
for j in rank_idx:
    fused_rows.append({
        "Position (RBD idx)": int(j),
        "Global Position": int(j + RBD_OFFSET),
        "Top1 (fused)": idx_to_aa[int(top1_f[j])],
        "Top2 (fused)": idx_to_aa[int(top2_f[j])],
        "Top1 Prob (fused)": float(p1f[j]),
        "Top2 Prob (fused)": float(p2f[j]),
        "Margin (fused)": float(margin_f[j]),
        "Entropy (fused)": float(entropy_f[j]),
        "Score (fused)": float(score_fused[j]),
    })

aegis_fused_df = pd.DataFrame(fused_rows)
try:
    aegis_fused_df.to_excel("aegis_mc_markov_fused_hotspots.xlsx", index=False)
except Exception:
    aegis_fused_df.to_csv("aegis_mc_markov_fused_hotspots.csv", index=False)

print("\nSaved: aegis_mc_markov_fused_hotspots.(xlsx|csv)")
    
#%% Deciphering my Model

# Obtain the latent representations
encoders = Model(inputs=input_layer, outputs=encoded)
latent_vectors = encoders.predict(encoded_sequences)

#%% Predict future sequences and normalize

# mean_prob_mc is already a normalized (L, V) distribution across amino acids per position.
L, V = mean_prob_mc.shape
rows = np.arange(L)

# Top1 / Top2 per position from MC-mean distribution
order = np.argsort(-mean_prob_mc, axis=1)
top1_idx = order[:, 0]
top2_idx = order[:, 1]

p1 = mean_prob_mc[rows, top1_idx]
p2 = mean_prob_mc[rows, top2_idx]
margin = p1 - p2

# Optional: decode the "consensus" sequence as the Top1 AA per position
int_to_amino_acid = {i: aa for i, aa in enumerate(amino_acids)}
decoded_consensus = np.array([int_to_amino_acid[i] for i in top1_idx], dtype=object)

# If you want the same 'mutation_df' you had before, but from MC means:
RBD_OFFSET = 331
K = 10  # top K ambiguous sites (smallest margin)
min_diff_positions = np.argsort(margin)[:K]

mutation_rows = []
for pos in min_diff_positions:
    aa1 = int_to_amino_acid[int(top1_idx[pos])]
    aa2 = int_to_amino_acid[int(top2_idx[pos])]
    confidence = 1.0 - float(margin[pos])  # your legacy "1 - margin" score
    mutation_rows.append({
        "Position": int(pos),
        "Global Position": int(pos + RBD_OFFSET),
        "Mutation": f"{aa1} -> {aa2}",
        "Top1 Prob": float(p1[pos]),
        "Top2 Prob": float(p2[pos]),
        "Confidence Score": confidence
    })

mutation_df = pd.DataFrame(mutation_rows)
mutation_df.to_excel("top_mutations_mc.xlsx", index=False)
print(mutation_df)

#%% plots and graphs

#setup for all graphs
import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
plt.rcParams.update({"figure.dpi": 160, "savefig.dpi": 300, "font.size": 10})
os.makedirs("figures", exist_ok=True)
L, V = mean_prob_mc.shape
pos_axis = np.arange(L) + RBD_OFFSET

#%% Entropy strip across positions
plt.figure(figsize=(10, 1.75))
plt.plot(pos_axis, entropy_pos, lw=1.5)
plt.xlabel("Spike residue (global numbering)"); plt.ylabel("Entropy")
plt.title("Shannon entropy per position (MC mean distribution)")
plt.tight_layout(); plt.savefig("figures/entropy_strip.png")


#%% Epistemic variance strip across positions
plt.figure(figsize=(10, 1.75))
plt.plot(pos_axis, epi_var_pos, lw=1.5)
plt.xlabel("Spike residue (global numbering)"); plt.ylabel("Epistemic var")
plt.title("Epistemic variance per position (MC dropout)")
plt.tight_layout(); plt.savefig("figures/epistemic_var_strip.png")


#%%Top-K hotspots bar chart (MC)
def minmax(x):
    lo, hi = float(np.min(x)), float(np.max(x))
    return np.zeros_like(x) if hi <= lo + 1e-12 else (x - lo) / (hi - lo)

entropy_norm = entropy_pos / np.log(V)
margin_norm  = 1.0 - minmax(margin)
epi_norm     = minmax(epi_var_pos)

w_margin, w_entropy, w_epi = 0.5, 0.3, 0.2
score_mc = (margin_norm**w_margin) * (entropy_norm**w_entropy) * (epi_norm**w_epi)

TOPK = 15
rank_idx = np.argsort(-score_mc)[:TOPK]
labels = [f"{pos_axis[i]} ({amino_acids[top1_idx[i]]}→{amino_acids[top2_idx[i]]})" for i in rank_idx]

plt.figure(figsize=(10, 3.5))
plt.bar(range(TOPK), score_mc[rank_idx])
plt.xticks(range(TOPK), labels, rotation=60, ha="right")
plt.ylabel("Composite hotspot score"); plt.title("Top predicted hotspots (MC)")
plt.tight_layout(); plt.savefig("figures/top_hotspots_mc.png")


#%% Margin vs Entropy scatter

plt.figure(figsize=(5, 4))
plt.scatter(margin, entropy_pos, s=12, alpha=0.6)
for pos in rank_idx[:5]:  # annotate a few top sites
    plt.annotate(str(pos_axis[pos]), (margin[pos], entropy_pos[pos]), fontsize=8)
plt.xlabel("Top1 − Top2 margin (MC)"); plt.ylabel("Entropy")
plt.title("Ambiguity vs. dispersion")
plt.tight_layout(); plt.savefig("figures/margin_vs_entropy.png");

#%% Probability heatmap

subset = rank_idx[:10]  # top 10 hotspots
P = mean_prob_mc[subset, :]  # shape: 10 × V

plt.figure(figsize=(10, 3.5))
plt.imshow(P.T, aspect="auto", origin="lower", interpolation="nearest")
plt.yticks(np.arange(V), amino_acids)
plt.xticks(range(len(subset)), [str(pos_axis[i]) for i in subset], rotation=0)
plt.colorbar(label="Probability")
plt.xlabel("Spike residue"); plt.ylabel("Amino acid")
plt.title("Per-AA probabilities at top hotspots (MC mean)")
plt.tight_layout(); plt.savefig("figures/heatmap_probs_top_hotspots.png")

#%% Markov transition heatmap (OK this doesn't really work idk why)

# Average transition matrix across positions
Tavg = Tmat.mean(axis=0)  # V×V
plt.figure(figsize=(5, 4))
plt.imshow(Tavg, aspect="auto", origin="lower", interpolation="nearest")
plt.xticks(np.arange(V), amino_acids, rotation=90); plt.yticks(np.arange(V), amino_acids)
plt.colorbar(label="P(next=b | current=a)")
plt.title("AA→AA transition probabilities (Markov avg)")
plt.tight_layout(); plt.savefig("figures/markov_transition_avg.png");

# Transition matrix for a selected hotspot (first in rank_idx)
j0 = int(rank_idx[0])
plt.figure(figsize=(5, 4))
plt.imshow(Tmat[j0], aspect="auto", origin="lower", interpolation="nearest")
plt.xticks(np.arange(V), amino_acids, rotation=90); plt.yticks(np.arange(V), amino_acids)
plt.colorbar(label="P(next=b | current=a)")
plt.title(f"Markov transitions at position {pos_axis[j0]}")
plt.tight_layout(); plt.savefig(f"figures/markov_transition_pos_{pos_axis[j0]}.png")

#%% MC vs Fused (Markov) rank comparison

# ranks (lower = better)
rank_mc = np.argsort(-score_mc)
rank_fused = np.argsort(-score_fused)
inv_rank_mc = np.empty_like(rank_mc); inv_rank_mc[rank_mc] = np.arange(L)
inv_rank_fused = np.empty_like(rank_fused); inv_rank_fused[rank_fused] = np.arange(L)

plt.figure(figsize=(5,4))
plt.scatter(inv_rank_mc, inv_rank_fused, s=8, alpha=0.5)
plt.xlabel("Rank (MC)"); plt.ylabel("Rank (Fused)")
plt.title("Hotspot rank shift with Markov fusion")
# annotate a few movers
delta = inv_rank_mc - inv_rank_fused
movers = np.argsort(-np.abs(delta))[:5]
for m in movers:
    plt.annotate(str(pos_axis[m]), (inv_rank_mc[m], inv_rank_fused[m]), fontsize=8)
plt.tight_layout(); plt.savefig("figures/rank_shift_mc_vs_fused.png");

#%% Modeling top_mutations_mc

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({"figure.dpi": 160, "savefig.dpi": 300, "font.size": 10})
os.makedirs("figures", exist_ok=True)

# 1) Load
df = pd.read_excel("top_mutations_mc.xlsx")
# ensure expected columns exist
expected = ["Position","Global Position","Mutation","Top1 Prob","Top2 Prob","Confidence Score"]
missing = [c for c in expected if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns in Excel: {missing}")

# clean / types
for c in ["Position","Global Position"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")
for c in ["Top1 Prob","Top2 Prob","Confidence Score"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df = df.dropna(subset=["Global Position","Top1 Prob","Top2 Prob","Confidence Score"]).copy()

# derived
df["Margin (Top1-Top2)"] = df["Top1 Prob"] - df["Top2 Prob"]
df["Top2/Top1 Ratio"] = df["Top2 Prob"] / np.clip(df["Top1 Prob"], 1e-12, None)

# sort variants for top-K displays
TOPK = 15
df_sorted_conf = df.sort_values("Confidence Score", ascending=False).head(TOPK).reset_index(drop=True)
df_sorted_margin_small = df.sort_values("Margin (Top1-Top2)", ascending=True).head(TOPK).reset_index(drop=True)

# 2) Bar chart: Top-K by Confidence Score (your legacy score = 1 - margin)
labels_conf = [f"{int(gp)} ({m})" for gp, m in zip(df_sorted_conf["Global Position"], df_sorted_conf["Mutation"])]
plt.figure(figsize=(10, 3.6))
plt.bar(range(len(df_sorted_conf)), df_sorted_conf["Confidence Score"])
plt.xticks(range(len(df_sorted_conf)), labels_conf, rotation=60, ha="right")
plt.ylabel("Confidence (1 − margin)")
plt.title("Top predicted mutations by Confidence Score (MC)")
plt.tight_layout(); plt.savefig("figures/top_mutations_confidence.png")

# 3) Bar chart: Top-K most ambiguous (smallest margin)
labels_amb = [f"{int(gp)} ({m})" for gp, m in zip(df_sorted_margin_small["Global Position"], df_sorted_margin_small["Mutation"])]
plt.figure(figsize=(10, 3.6))
plt.bar(range(len(df_sorted_margin_small)), df_sorted_margin_small["Margin (Top1-Top2)"])
plt.xticks(range(len(df_sorted_margin_small)), labels_amb, rotation=60, ha="right")
plt.ylabel("Margin (Top1 − Top2)")  # small = more ambiguous
plt.title("Most ambiguous positions (smallest margin)")
plt.tight_layout(); plt.savefig("figures/top_mutations_small_margin.png")

# 4) Scatter across the RBD: Confidence vs Position
plt.figure(figsize=(8, 3.4))
plt.scatter(df["Global Position"], df["Confidence Score"], s=16, alpha=0.7)
plt.xlabel("Spike residue (global numbering)")
plt.ylabel("Confidence (1 − margin)")
plt.title("Confidence vs position (all candidates)")
# annotate a few key sites if present
for key_gp in [484, 501, 502]:
    hit = df.loc[df["Global Position"] == key_gp]
    if not hit.empty:
        x = float(hit["Global Position"].iloc[0])
        y = float(hit["Confidence Score"].iloc[0])
        plt.scatter([x], [y], s=40)
        plt.annotate(str(key_gp), (x, y), fontsize=9, xytext=(4, 4), textcoords="offset points")
plt.tight_layout(); plt.savefig("figures/scatter_confidence_vs_position.png")

# 5) Side-by-side bars for Top1 vs Top2 probs at Top-K ambiguous sites
K = min(10, len(df_sorted_margin_small))
subset = df_sorted_margin_small.iloc[:K]
x = np.arange(K)
w = 0.42
plt.figure(figsize=(10, 3.6))
plt.bar(x - w/2, subset["Top1 Prob"], width=w, label="Top1 Prob")
plt.bar(x + w/2, subset["Top2 Prob"], width=w, label="Top2 Prob")
xt = [f"{int(gp)}\n{mut}" for gp, mut in zip(subset["Global Position"], subset["Mutation"])]
plt.xticks(x, xt, rotation=0)
plt.ylabel("Probability")
plt.title("Top1 vs Top2 probabilities at most ambiguous sites")
plt.legend()
plt.tight_layout(); plt.savefig("figures/top1_top2_probs_ambiguous.png")

# 6) Lollipop chart: Top2/Top1 ratio (the closer to 1, the more ambiguous)
subset_ratio = df.sort_values("Top2/Top1 Ratio", ascending=False).head(TOPK).reset_index(drop=True)
plt.figure(figsize=(10, 3.6))
xp = np.arange(len(subset_ratio))
plt.hlines(y=xp, xmin=0.0, xmax=subset_ratio["Top2/Top1 Ratio"], colors="black", linewidth=1)
plt.plot(subset_ratio["Top2/Top1 Ratio"], xp, "o")
yt = [f"{int(gp)} ({m})" for gp, m in zip(subset_ratio["Global Position"], subset_ratio["Mutation"])]
plt.yticks(xp, yt)
plt.xlabel("Top2 / Top1 probability ratio")
plt.title("Ambiguity by Top2/Top1 ratio (higher = more ambiguous)")
plt.tight_layout(); plt.savefig("figures/lollipop_top2_over_top1.png")

# 7) Export figure-ready tables (optional)
df_sorted_conf.to_excel("figures/table_topK_confidence.xlsx", index=False)
df_sorted_margin_small.to_excel("figures/table_topK_small_margin.xlsx", index=False)

print("Saved figures to: figures/")

#%% this is just for poster. can delete after

# If you didn't store entropy_pos earlier, compute it now:
eps = 1e-12
if 'entropy_pos' not in globals():
    entropy_pos = -(mean_prob_mc * np.log(np.clip(mean_prob_mc, eps, None))).sum(axis=1)

# Make a 1 x L heatmap (row = entropy)
H = entropy_pos[np.newaxis, :]  # shape (1, L)

plt.figure(figsize=(10, 1.8))
im = plt.imshow(H, aspect="auto", origin="lower", interpolation="nearest")
plt.yticks([0], ["Entropy"])
# Show global spike numbering on x-axis sparsely
tick_idx = np.linspace(0, L-1, 12, dtype=int)
plt.xticks(tick_idx, (pos_axis[tick_idx]).astype(int), rotation=0)

cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
cbar.set_label("Shannon entropy (nats)")

plt.xlabel("Spike residue (global numbering)")
plt.title("Per-residue uncertainty from MC-dropout (Shannon entropy)")
plt.tight_layout()
plt.savefig("figures/Fig2_entropy_heatmap.png")
plt.close()
