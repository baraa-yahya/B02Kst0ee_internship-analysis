import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
from math import pi

# Load the pickle file into a DataFrame
pkl_file1 = "/sps/lhcb/yahya/Internship/test_Bu/pkl/Bu2Kpipiee_all_years_MC_new_root_files_no_q2_cut.pkl"
data_file1 = pd.read_pickle(pkl_file1)
data_file1 = data_file1[data_file1['q2_NB']>14.0]
file_name1 = os.path.splitext(os.path.basename(pkl_file1))[0]  # Remove path and file extension
decay, year = file_name1.split('_')[:2]  # Split by underscore and take first two parts
# Specify the column for which you want to create the histogram
print(data_file1.columns)
column_name = "B_DTF_PV_B_M_0"
#column_name = "B_DTF_PV_Kstar_M_0"
#column_name ="ctk_lhcb"

print(data_file1[column_name].count())

# Plot the histogram
plt.figure()
n_bins=50
plot_range = (data_file1[column_name].min(), data_file1[column_name].max())  # Change this to your desired range
plt.hist(data_file1[column_name], n_bins, weights=data_file1['w_Track_Trig_Kin'], range=plot_range, histtype='step', density=True, color='skyblue', edgecolor='blue', label='Without kpipi_weight')
plt.hist(data_file1[column_name], n_bins, weights=data_file1['w_Track_Trig_Kin'] * data_file1['kpipi_weight_right'], range=plot_range, density=True, histtype='step', color='skyblue', edgecolor='red', label='With kpipi_weight')
# Add labels and title
plt.xlabel("$m(K^{+}pi^{-})$ [MeV/$c^{2}$]")
#plt.xlabel("$m(K^{*0}e^{+}e^{-})$ [MeV/$c^{2}$]")
#plt.xlabel(r"$\phi$")  # Use LaTeX formatting for Greek letters

plt.ylabel('Normalized Events')
# Set x-axis limits
plt.xlim(plot_range)

# Add a legend
plt.legend()

# Specify the output file name using decay, year, and column name
output_file = f"{decay}_{year}_{column_name}_histogram.pdf"

# If the file already exists, delete it
if os.path.exists(output_file):
    os.remove(output_file)

# Save the plot as PDF
plt.savefig(output_file)
plt.close()  # Close the plot to avoid duplicate figures being saved
