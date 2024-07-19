import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
from math import pi

# Load the pickle file into a DataFrame
#pkl_file1 = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/Bu2Kpipiee_all_years_MC.pkl"
pkl_file1 = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeSignal_all_years_MC.pkl"
# pkl_file1="/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeFlat_all_years_MC.pkl"
data_file1 = pd.read_pickle(pkl_file1)
df1 = data_file1[data_file1['XGBOutput']>0.9]# Extract decay and year information from the filename
file_name1 = os.path.splitext(os.path.basename(pkl_file1))[0]  # Remove path and file extension
decay, year = file_name1.split('_')[:2]  # Split by underscore and take first two parts
# Specify the column for which you want to create the histogram
print(df1.columns)
column_name = "q2_true"
# column_name = "ctl_lhcb"
print(df1[column_name].count())

# Plot the histogram
plt.figure()
n_bins=60
plot_range = (df1[column_name].min(), df1[column_name].max())  # Change this to your desired range
# plt.hist(df1[column_name], n_bins, weights=df1['w_Track_Trig_Kin'], range=plot_range, density=True, histtype='step', color='skyblue', edgecolor='red', label='With kpipi_weight')
plt.hist(df1[column_name], n_bins, range=plot_range, histtype='step', color='skyblue', edgecolor='blue')
plt.hist(df1['q2_NB'], n_bins, range=plot_range, histtype='step', color='skyblue', edgecolor='red')


plt.xlabel(r"$m(K^{\ast 0}e^{+}e^{-})$ [MeV/$c^{2}$]")
# plt.xlabel(r"$\phi$")  # Use LaTeX formatting for Greek letters
#plt.xlabel(r"$\cos (\theta_{k})$")
#plt.xlabel(r"$\cos (\theta_{k})$")
plt.ylabel('Normalized Events')
# Set x-axis limits
plt.xlim(plot_range)

# Add a legend
# plt.legend()

# Specify the output file name using decay, year, and column name
output_file = f"{decay}_{year}_{column_name}_histogram.pdf"

# If the file already exists, delete it
if os.path.exists(output_file):
    os.remove(output_file)

# Save the plot as PDF
plt.savefig(output_file)
plt.close()  # Close the plot to avoid duplicate figures being saved
