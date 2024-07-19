import numpy as np
import pandas as pd

file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/Bu2Kpipiee_all_years_MC.pkl"
data = pd.read_pickle(file_path)
print("Data loaded successfully.")

data = data[data['XGBOutput'] > 0.99]
print(f"Filtered data shape: {data.shape}")

# Check for the column existence
assert 'B_DTF_PV_B_M_0' in data.columns, "Column 'B_DTF_PV_B_M_0' does not exist in data."

# Assuming the column you want to fit is named 'B_DTF_PV_B_M_0'
mass_data_new = data['B_DTF_PV_B_M_0']
weights = data['w_Track_Trig_Kin'] * data['kpipi_weight']

# Compute weighted histogram
nbins = 50
obs_min = mass_data_new.min()
obs_max = mass_data_new.max()

counts, bins = np.histogram(mass_data_new, bins=nbins, range=(obs_min, obs_max), weights=weights)

# Convert counts and bins to DataFrame
histogram_data = pd.DataFrame({
    'bin_left_edge': bins[:-1],
    'bin_right_edge': bins[1:],
    'counts': counts
})

print("Histogram data as DataFrame:")
print(histogram_data.head())
