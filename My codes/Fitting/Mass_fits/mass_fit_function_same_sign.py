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
import mplhep
from scipy.integrate import simps

# Add the site-packages directory to sys.path
sys.path.append('/sps/lhcb/yahya/setup/miniconda_install/envs/penv/lib/python3.11/site-packages')


# Load the data from the provided file
file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0Jpsi2eeSS_all_years_Data.pkl"
data = pd.read_pickle(file_path)
data = data[data['XGBOutput']>0.8]
# Assuming the column you want to fit is named 'B_DTF_PV_B_M_0'
mass_data = data['B_DTF_PV_B_M_0']
obs = zfit.Space('mass', limits=(mass_data.min(), 6500))
    # Create the EXPSTEP PDF
params = pdfs.paramExpStep("my_pdf_prefix")
expstep_pdf = pdfs.ExpStep(obs, "my_pdf_suffix", params)

    # Create the data
data_ = zfit.data.Data.from_numpy(obs=obs, array=mass_data.values)

    # Create the loss function
nll = zfit.loss.UnbinnedNLL(model=expstep_pdf, data=data_)

    # Create a minimizer
minimizer = zfit.minimize.Minuit()

    # Perform the minimization
result = minimizer.minimize(nll)
print(result)
suffix = f"comb_bkg_fit_XGB_{0}_"
b2kstll.plot.plot_distributions(result, suffix=suffix, labels=['data'])