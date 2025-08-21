import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors

smiles_file = "Excel Files/Cleaned_Boronic_Acids.xlsx"  # Update with your actual file path
df = pd.read_excel(smiles_file)

# Initialize Mordred Calculator for 2D and 3D descriptors
calc_3d = Calculator(descriptors)  # Compute both 2D & 3D descriptors

# Store results
descriptors3d_list = []

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

        # Compute 3D Descriptors
        descriptors_3d = calc_3d(mol_3d)

        # Store results in a dictionary
        descriptors3d_list.append({"SMILES": smiles, **dict(descriptors_3d)})
    except:
        print(f'{index} is not being processed!')
    print(f'{index} number of molecules processsed!')

# Convert to DataFrame and save results
results_df = pd.DataFrame(descriptors3d_list)
results_df.to_excel("Excel Files/Boronic_Mordred_3D.xlsx", index=False)

print(f"Descriptor calculation completed. Results saved to 'Boronic_Mordred_3D.xlsx'.")
