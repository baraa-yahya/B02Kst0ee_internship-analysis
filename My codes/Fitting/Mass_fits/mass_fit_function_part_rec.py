import sys
import os
import numpy as np
import pandas as pd
import zfit
import matplotlib.pyplot as plt
from math import pi
import b2kstll
from b2kstll.models import pdfs
from b2kstll.models.pdfs import *
from b2kstll.utils import *

# Add the site-packages directory to sys.path
sys.path.append('/sps/lhcb/yahya/setup/miniconda_install/envs/bd2ksteeEnv/lib/python3.11/site-packages')

# Load the data from the provided file
file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/Bu2Kpipiee_all_years_MC.pkl"
data = pd.read_pickle(file_path)
print("Data loaded successfully.")

data = data[data['XGBOutput'] > 0.99]
print(f"Filtered data shape: {data.shape}")

# Check for the column existence
assert 'B_DTF_PV_B_M_0' in data.columns, "Column 'B_DTF_PV_B_M_0' does not exist in data."

# Assuming the column you want to fit is named 'B_DTF_PV_B_M_0'
mass_data_new = data['B_DTF_PV_B_M_0']
print(f"Mass data loaded with {len(mass_data_new)} entries.")
# Define the observable
obs = zfit.Space('mass', limits=(mass_data_new.min(), mass_data_new.max()))
print(f"Observable defined with limits: {obs.limits}")

# Create the EXPSTEP PDF
params = pdfs.paramExpStep("my_pdf_prefix")
expstep_pdf = pdfs.ExpStep(obs, "my_pdf_suffix", params)
print("EXPSTEP PDF created.")

# Create the data
data_ = zfit.data.Data.from_numpy(obs=obs, array=mass_data_new.values, weights = data['w_Track_Trig_Kin'] * data['kpipi_weight']
)
print("Zfit data object created.")

# Create the loss function
nll = zfit.loss.UnbinnedNLL(model=expstep_pdf, data=data_)
print("Loss function created.")

# Create a minimizer
minimizer = zfit.minimize.Minuit()
print("Minimizer created.")
# Perform the minimization
result = minimizer.minimize(nll)
b2kstll.plot.plot_distributions(result, suffix="part_rec_bkg_test", labels=['data'])