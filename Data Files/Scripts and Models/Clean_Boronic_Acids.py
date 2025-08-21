import pandas as pd

# Load the Excel file
file_path = "Excel Files/Boronic_Acids_SMILES.xlsx"  # Replace with your actual file path
df = pd.read_excel(file_path)

# Function to convert melting point values
def process_melting_point(mp):
    if isinstance(mp, str) and '-' in mp:
        values = list(map(float, mp.split('-')))  # Convert range to list of floats
        return sum(values) / len(values)  # Compute average
    try:
        return float(mp)  # Convert single values directly to float
    except ValueError:
        return None  # Handle non-numeric cases

# Process the dataset
df_cleaned = df[['Name', 'SMILES', 'Melting Point']].copy()
df_cleaned['Melting Point'] = df_cleaned['Melting Point'].apply(process_melting_point)

# Remove rows with empty (NaN) values
df_cleaned.dropna(inplace=True)

# Save to a new Excel file
output_file = "Excel Files/Cleaned_Boronic_Acids.xlsx"
df_cleaned.to_excel(output_file, index=False)

print(f"Cleaned data saved to {output_file}")
