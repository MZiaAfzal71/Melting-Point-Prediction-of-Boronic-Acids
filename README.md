# Melting-Point-Prediction-of-Boronic-Acids

## Boronic Acids Descriptor Analysis

## Overview

This repository contains scripts and data for extracting, processing, and analyzing boronic acids data using various feature extraction techniques and machine learning models. The dataset is sourced from [Organoborons](https://organoborons.com/boronic-acids/page-1.html) and analyzed using multiple descriptors, including Morgan fingerprints, MACCS fingerprints, Mordred descriptors, and Coulomb matrices. The final performance of these descriptors is evaluated using five machine learning models.

## Folder Structure

```
Data Files/
|-- Excel Files/
|   |-- Boronic_pages.xlsx
|   |-- Boronic_Acids_IM.xlsx
|   |-- Boronic_Acids_SMILES.xlsx
|   |-- Cleaned_Boronic_Acids.xlsx
|   |-- Boronic_Bonds_Desc_Boron_En.xlsx
|   |-- Boronic_Morgan_fingerprint.xlsx
|   |-- Boronic_MACCS_fingerprint.xlsx
|   |-- Boronic_Mordred_3D.xlsx
|   |-- Boronic_Mordred_3DC.xlsx
|   |-- CoulombMatrix_BoronicAcids_Desc.xlsx
|-- XGBoost/
|   |-- (Generated model files used for graph generation)
```

## Workflow

### Step 1: Fetching Boronic Acids Links

- **Script:** `GetKeys_7pages_boronicacids_website.py`
- **Function:** Scrapes boronic acid links from the first 7 pages of the website.
- **Output:** `Excel Files/Boronic_pages.xlsx`

### Step 2: Extracting Data

- **Script:** `Get_Data_Boronic_Acids.py`
- **Function:** Extracts names, InChIKey, SMILES (initially empty), and melting points.
- **Output:**
  - `Excel Files/Boronic_Acids_IM.xlsx`
  - `Excel Files/Boronic_Acids_SMILES.xlsx` (now with SMILES filled)

### Step 3: Cleaning Data

- **Script:** `Cleaned_Boronic_Acids.py`
- **Function:**
  - Removes rows with empty melting points.
  - Averages melting points when provided as a list.
- **Output:** `Excel Files/Cleaned_Boronic_Acids.xlsx` (605 molecules retained)

### Step 4: Generating Molecular Descriptors

#### (a) Custom Descriptor

- **Script:** `describemoleculeEnhanced.py`
- **Function:** Generates a 20-feature vector as described in our research paper.
- **Output:** `Excel Files/Boronic_Bonds_Desc_Boron_En.xlsx`

#### (b) Morgan & MACCS Fingerprints

- **Script:** `MorganMACCSFingerprints.py`
- **Function:** Generates molecular fingerprints.
- **Outputs:**
  - `Excel Files/Boronic_Morgan_fingerprint.xlsx`
  - `Excel Files/Boronic_MACCS_fingerprint.xlsx`

#### (c) Mordred Descriptors

- **Script:** `Compute_Mordred_Alkanes.py`
- **Function:** Computes Mordred descriptors.
- **Output:** `Excel Files/Boronic_Mordred_3D.xlsx`
- **Cleaning Step:** `Clean_Mordred_3D_Desc.py`
- **Final Output:** `Excel Files/Boronic_Mordred_3DC.xlsx`

#### (d) Coulomb Matrix Descriptors

- **Script:** `Compute_CoulombMatrix_Desc.py`
- **Function:** Computes the Coulomb matrix descriptor.
- **Output:** `Excel Files/CoulombMatrix_BoronicAcids_Desc.xlsx`

### Step 5: Machine Learning Analysis

- **Feature vectors analyzed:**
  1. Custom descriptor (Boronic\_Bonds\_Desc\_Boron\_En.xlsx)
  2. Morgan fingerprint (Boronic\_Morgan\_fingerprint.xlsx)
  3. MACCS fingerprint (Boronic\_MACCS\_fingerprint.xlsx)
  4. Mordred descriptors (Boronic\_Mordred\_3DC.xlsx)
  5. Coulomb matrix (CoulombMatrix\_BoronicAcids\_Desc.xlsx)

#### (a) Decision Tree

- **Script:** `DecisionTree.py`
- **Function:** Computes Mean Absolute Error (MAE) and R² scores for the selected descriptor.
- **Usage:** Uncomment the file path in the script to analyze different feature vectors.

#### (b) Other ML Models

- **Scripts:**
  - `LightGBM.py`
  - `RandomForest.py`
  - `SupportVectorM.py`
  - `XGBoostModel.py`
- **Function:** Each script applies a different ML model for performance evaluation.

### Step 6: Graph Generation (XGBoost Only)

- **Script:** `XGBoostModel.py`
- **Function:**
  - Generates results saved in `Data Files/XGBoost/`
  - Used for final graph visualization
- **Script:** `GenerateGraphs.py`
- **Function:** Produces publication-ready graphs from XGBoost results.

## Usage

1. Run scripts sequentially as per workflow.
2. For ML evaluation, comment/uncomment relevant descriptor file paths in `DecisionTree.py` or other ML scripts.
3. Use `GenerateGraphs.py` to visualize final results.
