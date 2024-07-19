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
import mplhep
# Add the site-packages directory to sys.path
sys.path.append('/sps/lhcb/yahya/setup/miniconda_install/envs/bd2ksteeEnv/lib/python3.11/site-packages')


# Load the data from the provided file
file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeSignal_all_years_MC.pkl"
data = pd.read_pickle(file_path)
data = data[data['XGBOutput'] > 0.99]

# Assuming the column you want to fit is named 'B_DTF_PV_B_M_0'
mass_data = data['B_DTF_PV_B_M_0']

# Define the observable
obs = zfit.Space('mass', limits=(mass_data.min(), mass_data.max()))

# Create the DCB_DG PDF
params = paramDCBDGauss("my_pdf_prefix")
dcb_dg_pdf = pdfs.doubleCBDGauss(obs, "my_pdf_suffix", params)

# Create the data
data_ = zfit.data.Data.from_numpy(obs=obs, array=mass_data.values)

# Create the loss function
nll = zfit.loss.UnbinnedNLL(model=dcb_dg_pdf, data=data_)

# Create a minimizer
minimizer = zfit.minimize.Minuit()

# Perform the minimization
result = minimizer.minimize(nll)

b2kstll.plot.plot_distributions(result, suffix="MC_signal_" , weights = data['w_Track_Trig_Kin'], labels =['gauss a','gauss b','crystal ball'])