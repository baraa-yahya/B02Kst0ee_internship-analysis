import numpy as np
import pandas as pd
from scipy.integrate import quad
import pickle

# Load the pickle file
pkl_file1 = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeSignal_all_years_MC.pkl"
try:
    with open(pkl_file1, 'rb') as f:
        data_file1 = pickle.load(f)
        print("Pickle file loaded successfully.")
except pickle.UnpicklingError as e:
    print(f"Error loading pickle file: {e}")
    exit(1)

# Filter DataFrame based on condition
df1 = data_file1[data_file1['XGBOutput'] > 0.9]

# Specify the column names
q2_true_column = "q2_true"
q2_NB_column = "q2_NB"

# Calculate Integral of q2_NB from 14.0 to its maximum
q2_NB = df1[q2_NB_column]
q2_true = df1[q2_true_column]

q2_max = 20.0
# Function to integrate q2_NB
def integrate_q2_NB(min_value, max_value):
    return quad(lambda x: np.interp(x, q2_NB, q2_NB.values), min_value, max_value, limit=200)[0]

# Calculate q2_NB_integral
q2_NB_integral = integrate_q2_NB(q2_NB.min(), q2_max)
print(f"Q2_NB_integral calculated: {q2_NB_integral}")


# Adjusted function to integrate q2_true with increased limit on subdivisions
def integrate_q2_true(min_value, max_value):
    return quad(lambda x: np.interp(x, q2_true, q2_true.values), min_value, max_value, limit=200)[0]

print(f"q2_true min: {q2_true.min()}, q2_true max: {q2_true.max()} and q2_NB min: {q2_NB.min()}, q2_NB max: {q2_NB.max()}")

# Find the minimum q2_true_min where integrals are equal
def find_q2_true_min_for_equal_integral(df, q2_NB_integral):
    search_range = np.linspace(df[q2_true_column].min(), df[q2_true_column].max(), num=10000)
    for q2_true_min in search_range:
        try:
            integral_q2_true = integrate_q2_true(q2_true_min, q2_max)
            print(f"Trying q2_true_min: {q2_true_min}, Integral_q2_true: {integral_q2_true}, Q2_NB_integral: {q2_NB_integral}")
            if np.isclose(integral_q2_true, q2_NB_integral, atol=1e-3):
                return q2_true_min
        except ValueError:
            continue
    return None

# Example usage
q2_true_min = find_q2_true_min_for_equal_integral(df1, q2_NB_integral)
if q2_true_min is not None:
    print(f"Minimum {q2_true_column} value (q2_true_min) where integrals are equal: {q2_true_min}")
else:
    print(f"No value found where integrals are equal for {q2_true_column}")