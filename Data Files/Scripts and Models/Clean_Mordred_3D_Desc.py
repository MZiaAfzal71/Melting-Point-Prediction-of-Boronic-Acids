import pandas as pd

# Define file paths
input_file = 'Excel Files/Boronic_Mordred_3D.xlsx'       # File containing Mordred descriptors
input_file1 = 'Excel Files/Cleaned_Boronic_Acids.xlsx'   # File containing cleaned boronic acids data
output_file = 'Excel Files/Boronic_Mordred_3DC.xlsx'    # Output file path

# Load Excel files into DataFrames
df1 = pd.read_excel(input_file)   # Data with Mordred descriptors
df2 = pd.read_excel(input_file1)  # Cleaned Boronic Acids data

# Select only numeric columns from df1 (Mordred descriptor data)
df3 = df1.select_dtypes(include=['number'])

# Remove columns where all values are 0 (i.e., keep only meaningful descriptors)
df3 = df3.loc[:, (df3 != 0).any(axis=0)]

# Define base columns from the cleaned Boronic Acids data
base_cols = ['Name', 'SMILES', 'Melting Point']

# Extract relevant base columns from df2
new_df = df2[base_cols].copy()

# Concatenate the base columns with the descriptor columns
final_df = pd.concat([new_df, df3], axis=1)

# Save the final merged DataFrame to an Excel file
final_df.to_excel(output_file, index=False)

print(f"Processing complete! Data saved to '{output_file}'")
