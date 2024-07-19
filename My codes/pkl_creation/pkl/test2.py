import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Load the pickle file into a DataFrame
pkl_file = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/Bu2Kpipiee_all_years_MC.pkl"
df = pd.read_pickle(pkl_file)
df =df[df['XGBOutput']>0.9]
# Extract decay and year information from the filename
file_name = os.path.splitext(pkl_file)[0]  # Remove file extension
decay, year = file_name.split('_')[:2]  # Split by underscore and take first two parts

# Specify the column for which you want to create the histogram
column_name = "B_DTF_PV_Kstar_M_0"

print(df.columns)
# Plot the histogram without plotting it
nbins=50
bin_heights, bin_edges, _ = plt.hist(df[column_name],weights=df['w_Track_Trig_Kin'] * df['kpipi_weight'],density=False, bins=nbins, histtype='step', color='skyblue', edgecolor='black', alpha=0)
bin_heights2, bin_edges2, _ = plt.hist(df[column_name],weights=df['w_Track_Trig_Kin'], density=False,bins=nbins, histtype='step', color='skyblue', edgecolor='black', alpha=0)

# Calculate the bin centers
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2

# Calculate the Poisson error for each bin
bin_errors = np.sqrt(bin_heights)
bin_errors2 = np.sqrt(bin_heights2)
# Plot the histogram with error bars
plt.errorbar(bin_centers, bin_heights, yerr=bin_errors,  fmt='o', color='red', ecolor='black',markersize=3, label=f'{column_name} with kpipi',elinewidth=1, capsize=4)
plt.errorbar(bin_centers2, bin_heights2, yerr=bin_errors2,  fmt='o', color='blue', ecolor='black',markersize=3, label=f'{column_name}',elinewidth=1, capsize=4)
plt.xlabel(column_name)
plt.ylabel('Frequency')
plt.title(f'Histogram of {column_name}')

# Set x-axis limits
plt.xlim(df[column_name].min(), df[column_name].max())

plt.legend()

# Specify the output file name using decay, year, and column name
output_file = f"{decay}_{year}_{column_name}_histogram.pdf"

# If the file already exists, delete it
if os.path.exists(output_file):
    os.remove(output_file)

# Save the plot as PDF
plt.savefig(output_file)

plt.show()
