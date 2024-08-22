# Project Overview

This repository contains data and codes for analyzing and simulating behavioral data from the Intentions Game, including both real and model-simulated data. Repo for the paper: "Self-Other Generalisation Shapes Social Interaction and Is Disrupted in Borderline Personality Disorder"

## Directory Structure

### 1. `Data`
**Real Behavioral Data**

- **intent_dat.csv**: Cleaned behavioural data from the Intentions Game.

### 2. `Data_Simulated`
**Model-Simulated Behavioral Data** for each group:
  - **g1**: Borderline Personality Disorder (BPD)
  - **g2**: Control Group (CON)

- **full_sim_g1.csv**: Full model results for the BPD group (Model M4).
- **full_sim_g2.csv**: Full model results for the Control group (Model M1).
- **data_sim_g1.csv**: Model-simulated actions for BPD group (Model M4).
- **data_sim_g2.csv**: Model-simulated actions for Control group (Model M1).
- **full_sim_g1_M3.csv**: Full model results for the BPD group (Model M3).
- **full_sim_g2_M3.csv**: Full model results for the Control group (Model M3).

### 3. `Models`
**Model Files** to fit and simulate data:

- **M1.m**: Model 1 file.
- **M2.m**: Model 2 file.
- **M3.m**: Model 3 file.
- **M4.m**: Model 4 file.
- **Beta.m**: Model Beta file.
- **Utilities**: Directory containing utility scripts for model simulations.

### 4. `Analysis`
**Scripts for Reproducing Results** from the paper:

- **BPD_CON_Analysis.R**: Core analysis script to produce outcomes as presented in the paper.
- **UtilityFunctions.R**: Custom utility functions used within the core analysis script.
    
