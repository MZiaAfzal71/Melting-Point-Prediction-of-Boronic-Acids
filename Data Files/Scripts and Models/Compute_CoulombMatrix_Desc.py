import pandas as pd
import rdkit
from rdkit.Chem import AllChem
from rdkit import Chem
import numpy as np
import ase
from dscribe.descriptors import CoulombMatrix


def smiles_str_to_rdkit_mol(smiles_str: str) -> rdkit.Chem.Mol:
    """
    Convert a SMILES string to an RDKit mol object and generate its 3D structure.

    Args:
        smiles_str (str): A SMILES string representing a molecule.

    Returns:
        rdkit.Chem.Mol: An RDKit mol object representing the molecule with an optimized 3D structure.
    """
    try:
        # Convert SMILES string to RDKit mol object
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            raise ValueError("Invalid SMILES string")

        # Add hydrogens to the molecule
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates using ETKDG
        if AllChem.EmbedMolecule(mol, AllChem.ETKDG()) != 0:
            raise ValueError("Failed to embed molecule")

        # Optimize 3D geometry using UFF force field
        AllChem.UFFOptimizeMolecule(mol)

        return mol
    except Exception as e:
        print(f"Error processing SMILES '{smiles_str}': {e}")
        return None


def ase_atoms_to_coulomb_matrix(ase_atoms: ase.Atoms) -> np.ndarray:
    """
    Convert an ASE Atoms object to a Coulomb matrix.

    Args:
        ase_atoms (ase.Atoms): ASE Atoms object.

    Returns:
        np.ndarray: Flattened Coulomb matrix representation of the molecule.
    """
    try:
        # Create a Coulomb matrix descriptor
        coulomb_matrix = CoulombMatrix(n_atoms_max=ase_atoms.get_global_number_of_atoms())

        # Compute the Coulomb matrix
        col_matrix = coulomb_matrix.create(ase_atoms)

        # Reshape the matrix for easier processing
        return col_matrix.flatten()
    except Exception as e:
        print(f"Error generating Coulomb matrix: {e}")
        return None


# Define file paths (keeping them unchanged)
input_file = 'Excel Files/Cleaned_Boronic_Acids.xlsx'
output_file = 'Excel Files/CoulombMatrix_BoronicAcids_Desc.xlsx'

# Load input Excel file
chem_file = pd.read_excel(input_file)

# Initialize list for storing descriptors
descriptors = []

# Process each molecule in the dataset
for i, sm in chem_file['SMILES'].items():
    try:
        mol = smiles_str_to_rdkit_mol(sm)
        if mol:
            # Convert RDKit mol to ASE Atoms object
            ase_atoms = ase.Atoms(
                numbers=[atom.GetAtomicNum() for atom in mol.GetAtoms()],
                positions=mol.GetConformer().GetPositions()
            )

            # Compute the Coulomb matrix
            col_mat = ase_atoms_to_coulomb_matrix(ase_atoms)

            # Append the descriptor to the list
            if col_mat is not None:
                descriptors.append(col_mat)
            else:
                descriptors.append([None] * 100)  # Placeholder for missing data

        else:
            descriptors.append([None] * 100)  # Placeholder for missing data

    except Exception as e:
        print(f"Error processing molecule {i}: {e}")
        descriptors.append([None] * 100)

    print(f'{i + 1} molecules processed!')

# Convert descriptors list into a DataFrame
desc_df = pd.DataFrame(descriptors)

# Assign meaningful column names
desc_df.columns = [f"desc_{i}" for i in range(desc_df.shape[1])]

# Extract essential base columns
base_cols = ['Name', 'SMILES', 'Melting Point']
base_df = chem_file[base_cols].copy()

# Concatenate base data with descriptors
final_df = pd.concat([base_df, desc_df], axis=1)

# Save the final dataset to Excel
final_df.to_excel(output_file, index=False)

print(f"Processing complete! Data saved to '{output_file}'")

