import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors
from rdkit.Chem.rdMolDescriptors import GetMACCSKeysFingerprint
import numpy as np


smiles_file = "Excel Files/Cleaned_Boronic_Acids.xlsx"  # Update with your actual file path
df = pd.read_excel(smiles_file)


# Output files
output_file1 = "Excel Files/Boronic_Morgan_fingerprint.xlsx"
output_file2 = "Excel Files/Boronic_MACCS_fingerprint.xlsx"

# Store results
descMorgan_list = []
descMACCS_list = []

for index, smiles in enumerate(df['SMILES']):
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))  # Convert SMILES to RDKit Mol
    if mol is None:
        print(f"Invalid SMILES at row {index}: {smiles}")
        continue
    try:
        # Generate 3D Conformer
        mol_3d = Chem.AddHs(mol)  # Add hydrogen atoms
        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())  # Generate 3D structure
        AllChem.UFFOptimizeMolecule(mol_3d)  # Optimize 3D geometry

        # Generate a Morgan fingerprint (ECFP4, radius=2)
        morgan_fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=512)

        # Compute MACCS fingerprint (166-bit)
        maccs_fp = GetMACCSKeysFingerprint(mol)

        # Store results in a dictionary
        descMorgan_list.append(np.array(morgan_fp))
        descMACCS_list.append(list(maccs_fp))
    except:
        print(f'{index} is not being processed!')
    print(f'{index} number of molecules processsed!')

# Convert to DataFrame and save results
results_df1 = pd.DataFrame(descMorgan_list)

results_df1.columns = [f"desc_{i}" for i in range(results_df1.shape[1])]

base_cols = ['Name', 'SMILES', 'Melting Point']
temp_df1 = df[base_cols].copy()

# Concatenate the base columns with the descriptor columns
final_df1 = pd.concat([temp_df1, results_df1], axis=1)

final_df1.to_excel(output_file1, index=False)

# Convert to DataFrame and save results
results_df2 = pd.DataFrame(descMACCS_list)

results_df2.columns = [f"desc_{i}" for i in range(results_df2.shape[1])]

# Concatenate the base columns with the descriptor columns
final_df2 = pd.concat([temp_df1, results_df2], axis=1)

final_df2.to_excel(output_file2, index=False)


print(f"Descriptor calculation completed. Results saved to 'Boronic_*_fingerprint.xlsx'.")
