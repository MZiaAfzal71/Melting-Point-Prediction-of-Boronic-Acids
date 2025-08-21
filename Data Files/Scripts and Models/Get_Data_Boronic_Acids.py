import pandas as pd
import requests
from bs4 import BeautifulSoup

def get_pubchem_info(inchikey):
    """
    Fetches the Canonical SMILES representation from PubChem using an InChIKey.

    Parameters:
    - inchikey (str): The InChIKey of the compound.

    Returns:
    - str: The Canonical SMILES string if found, otherwise an empty string.
    """
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/IUPACName,CanonicalSMILES/JSON'

    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for bad responses (4xx, 5xx)
        jstruct_res = response.json()
        return jstruct_res['PropertyTable']['Properties'][0].get('CanonicalSMILES', ' ')
    except (requests.exceptions.RequestException, KeyError, IndexError):
        return ' '  # Return empty string if any error occurs

# Define file paths
input_file = 'Excel Files/Boronic_pages.xlsx'
output_file1 = 'Excel Files/Boronic_Acids_IM.xlsx'
output_file2 = 'Excel Files/Boronic_Acids_SMILES.xlsx'

# Load Boronic pages from Excel
bor_file = pd.read_excel(input_file)

# Dictionary to store extracted Name, InChIKey, and Melting Point data
inch_melting = {'Name': [], 'InChIKey': [], 'Melting Point': []}

# Iterate over each page URL from the input file
for counter, pg in enumerate(bor_file['pages'], start=1):
    searchaddress = f'https://organoborons.com/{pg}'
    response = requests.get(searchaddress)

    if response.status_code == 200:
        soup = BeautifulSoup(response.text, 'html.parser')
        tables = soup.find_all('table')

        try:
            # Extract compound name from the first table
            inch_melting['Name'].append(tables[0].find('h1').text)

            # Extract data from the 4th table (index 3)
            rows = tables[3].find_all('tr')
            for row in rows:
                if 'Melting' in row.text:
                    melting_list = row.text.split()
                    inch_melting['Melting Point'].append(melting_list[2])  # Extract melting point value
                if 'InChIKey' in row.text:
                    inch_list = row.text.split()
                    inch_melting['InChIKey'].append(inch_list[1])  # Extract InChIKey

            print(f'{counter} entries processed successfully!')
        except (IndexError, AttributeError):
            print(f"Warning: Skipping {searchaddress} due to missing data.")
    else:
        print(f"Failed to retrieve {searchaddress} (Status Code: {response.status_code})")

# Ensure all extracted lists have the same length to prevent DataFrame misalignment
min_length = min(len(inch_melting['Name']), len(inch_melting['InChIKey']), len(inch_melting['Melting Point']))
inch_melting = {key: values[:min_length] for key, values in inch_melting.items()}

# Create DataFrame with an empty SMILES column
data = pd.DataFrame({
    'Name': inch_melting['Name'],
    'InChIKey': inch_melting['InChIKey'],
    'SMILES': [''] * len(inch_melting['Name']),  # Empty SMILES initially
    'Melting Point': inch_melting['Melting Point']
})

# Save the initial data (without SMILES) to an Excel file
data.to_excel(output_file1, index=False)

# Fetch SMILES for each compound using the PubChem API
for idx, inchi in enumerate(data['InChIKey']):
    data.loc[idx, 'SMILES'] = get_pubchem_info(inchi)
    print(f'{idx + 1}/{len(data)} compounds processed!')

# Save the final data with SMILES included
data.to_excel(output_file2, index=False)

print(f"Processing complete! Data saved to '{output_file2}'")
