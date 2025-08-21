import pandas as pd
import requests
from bs4 import BeautifulSoup


def extract_boronic_pages(html_content):
    """
    Extracts all hyperlinks containing the keyword 'data' from the given HTML content.

    Parameters:
    - html_content (str): The HTML content of a webpage.

    Returns:
    - list: A list of extracted URLs that contain the word 'data'.
    """
    soup = BeautifulSoup(html_content, 'html.parser')
    hrefs = [a['href'] for a in soup.find_all('a', href=True)]  # Extract all href links
    filtered_hrefs = list(filter(lambda x: 'data' in x, hrefs))  # Keep only links containing 'data'

    return filtered_hrefs


# Dictionary to store extracted page links
bor_pages = {'pages': []}

# Loop through pages 1 to 7 and extract data links
for i in range(1, 8):
    url = f'https://organoborons.com/boronic-acids/page-{i}.html'
    response = requests.get(url)

    if response.status_code == 200:  # Check if the request was successful
        bor_pages['pages'] += extract_boronic_pages(response.text)
    else:
        print(f"Failed to fetch {url} (Status Code: {response.status_code})")

# Convert extracted links into a DataFrame and save to an Excel file
data = pd.DataFrame(bor_pages)
data.to_excel('Excel Files/Boronic_pages.xlsx', index=False)

print("Extraction complete. Data saved to 'Excel Files/Boronic_pages.xlsx'")