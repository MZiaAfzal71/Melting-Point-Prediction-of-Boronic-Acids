from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score
import pandas as pd

# Load data
# input_file = 'Excel Files/CoulombMatrix_BoronicAcids_Desc.xlsx'
# input_file = 'Excel Files/Boronic_Mordred_3DC.xlsx'
# input_file = 'Excel Files/Boronic_Morgan_fingerprint.xlsx'
# input_file = 'Excel Files/Boronic_MACCS_fingerprint.xlsx'
input_file = 'Excel Files/Boronic_Bonds_Desc_Boron_En.xlsx'


chem_file = pd.read_excel(input_file)
chem_file.fillna(0, inplace=True)

X = chem_file.iloc[:, 3:]
y = chem_file['Melting Point']

X_train, X_valid, y_train, y_valid = train_test_split(X, y, train_size=0.8, random_state=42)

# Train Random Forest model
model = RandomForestRegressor( random_state=42)
model.fit(X_train, y_train)

# Predictions
y_pred = model.predict(X)

# Print results
print(f"The mean absolute error is: {mean_absolute_error(y, y_pred)}")
print(f"The R2 score is: {r2_score(y, y_pred)}")