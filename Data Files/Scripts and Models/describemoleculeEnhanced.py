import pandas as pd
from rdkit import Chem
import networkx as nx
import numpy as np


def get_mol_graph(mol):
    """Construct a graph from an RDKit molecule (with explicit hydrogens).
    Nodes: atom indices.
    Edges: bonds with attributes: bond_type and aromatic flag."""
    G = nx.Graph()
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        bt = bond.GetBondType()
        aromatic = bond.GetIsAromatic()
        G.add_edge(i, j, bond_type=bt, aromatic=aromatic)
    return G


def compute_path_sum(mol, path, G):
    """
    Compute the path sum for a given path.
    For each atom along the path, its atomic number (Z) is added.
    For the bond connecting the previous atom, apply a multiplier:
      - Aromatic bond: 1.5
      - Double bond (non-aromatic): 2.0
      - Triple bond (non-aromatic): 3.0
      - Otherwise: 1.0
    The reference atom (first in the path) is added normally.
    """
    total = 0.0
    for idx in range(len(path)):
        atom = mol.GetAtomWithIdx(path[idx])
        Z = atom.GetAtomicNum()
        if idx == 0:
            total += Z
        else:
            data = G.get_edge_data(path[idx - 1], path[idx])
            multiplier = 1.0
            if data.get("aromatic"):
                multiplier = 1.5
            else:
                bond_type = data.get("bond_type")
                if bond_type == Chem.rdchem.BondType.DOUBLE:
                    multiplier = 2.0
                elif bond_type == Chem.rdchem.BondType.TRIPLE:
                    multiplier = 3.0
            total += multiplier * Z
    return total


def compute_descriptor_vector(mol):
    """
    Compute a descriptor vector from the molecule.
    The molecule is first converted to include explicit hydrogens.
    The longest bond path length from the Boron reference atom is used as max_depth.
    For each atom at a bond distance d (d>=2) from the Boron reference,
    all shortest paths are obtained and the path sum is computed.
    The descriptor vector component for a given d is the cumulative sum over all such paths.
    """
    # Convert to include explicit hydrogens
    mol = Chem.AddHs(mol)

    # Identify Boron as the reference atom (choose the first Boron found)
    boron_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'B']
    if not boron_indices:
        return None
    ref = boron_indices[0]

    # Build graph
    G = get_mol_graph(mol)

    # Compute the shortest path lengths from reference
    lengths = nx.single_source_shortest_path_length(G, ref)
    if not lengths:
        return None
    max_depth = max(lengths.values())

    descriptor = {}  # keys: bond distance; values: cumulative path sum
    for atom_idx, dist in lengths.items():
        # Only consider atoms with distance >= 2
        if dist < 2:
            continue
        # Get all shortest paths between ref and this atom
        paths = list(nx.all_shortest_paths(G, source=ref, target=atom_idx))
        path_sum_total = sum(compute_path_sum(mol, path, G) for path in paths)
        descriptor.setdefault(dist, 0.0)
        descriptor[dist] += path_sum_total

    # Build the descriptor vector for distances 2 to max_depth
    vector = [descriptor.get(d, 0.0) for d in range(2, max_depth + 1)]
    return np.array(vector)


# Example usage:
# 1. Read an Excel file with a column 'smiles' containing SMILES codes for boronic acids.
input_file = 'Excel Files/Cleaned_Boronic_Acids.xlsx'
output_file = 'Excel Files/Boronic_Bonds_Desc_Boron_En.xlsx'
df = pd.read_excel(input_file)

# 2. For each SMILES code, convert it to canonical form and compute the descriptor vector.
descriptor_list = []
for smi in df['SMILES']:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        descriptor_list.append(None)
        continue
    canonical_smi = Chem.MolToSmiles(mol, canonical=True)
    mol = Chem.MolFromSmiles(canonical_smi)

    vec = compute_descriptor_vector(mol)
    descriptor_list.append(vec)

# Convert the list of descriptors (each a list of numbers) into a DataFrame.
# This will create separate columns for each element in the descriptor.
desc_df = pd.DataFrame(descriptor_list)
# Optionally, rename the descriptor columns (e.g., desc_0, desc_1, ...)
desc_df.columns = [f"desc_{i}" for i in range(desc_df.shape[1])]

# Create a new dataframe with the desired columns from the original data:
# 'Name', 'SMILES', and 'boiling_point'
base_cols = ['Name', 'SMILES', 'Melting Point']
new_df = df[base_cols].copy()

# Concatenate the base columns with the descriptor columns
final_df = pd.concat([new_df, desc_df], axis=1)

# Save the new dataframe to an Excel file
final_df.to_excel(output_file, index=False)

print(f'The processing is being finished and the results have been saved in {output_file}')
