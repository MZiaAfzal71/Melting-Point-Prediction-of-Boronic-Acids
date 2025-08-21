from sklearn.model_selection import train_test_split # To split data
from xgboost import XGBRegressor # The model to learn Fusion points of Alkanes
from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error  # To find errors
import pandas as pd

# input_file = 'Excel Files/CoulombMatrix_BoronicAcids_Desc.xlsx'
# input_file = 'Excel Files/Boronic_Mordred_3DC.xlsx'
# input_file = 'Excel Files/Boronic_Morgan_fingerprint.xlsx'
# input_file = 'Excel Files/Boronic_MACCS_fingerprint.xlsx'
input_file = 'Excel Files/Boronic_Bonds_Desc_Boron_En.xlsx'

# output_file = 'XGBoost Results/CoulombMatrix.xlsx'
# output_file = 'XGBoost Results/Mordred.xlsx'
# output_file = 'XGBoost Results/Morgan.xlsx'
# output_file = 'XGBoost Results/MACCS.xlsx'
output_file = 'XGBoost Results/Descriptor.xlsx'

chem_file = pd.read_excel(input_file)
chem_file.fillna(0, inplace=True)

no_alk = len(chem_file)

X = chem_file.iloc[:, 3:]
y = chem_file['Melting Point']

X_train, X_valid, y_train, y_valid = train_test_split(X, y, train_size=0.8, random_state=42)


# XGBRegressor with respect to SingularValues to measure Melting Points  of hydrocarbons
my_model = XGBRegressor(random_state = 42)


# Fit the model
my_model.fit(X_train, y_train)

y_pred = my_model.predict(X)

# Print results
print(f'The mean absolute error is: {mean_absolute_error(y, y_pred)}')
print(f'The R2 score is {r2_score(y_pred, y)}')

# To save the results in a file
results = pd.DataFrame({'Name' : chem_file['Name'], 'Observed' : chem_file['Melting Point'],
                        'Predicted' : y_pred, 'Difference' : abs(y - y_pred)})
results.to_excel(output_file, index=False)

