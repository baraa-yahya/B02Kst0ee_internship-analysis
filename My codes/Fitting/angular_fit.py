import numpy as np
import pandas as pd
import zfit
import b2kstll
from b2kstll.models import pdfs
from b2kstll.models.pdfs import *
from b2kstll.utils import *
import tensorflow as tf
from tensorflow.python.framework.errors_impl import InvalidArgumentError
import sys
from b2kstll.models.pdfs import Angular3D
import yaml
from math import pi
import argparse
import re
#from analysis.efficiency import load_efficiency_model
import b2kstll
import subprocess

# Add the site-packages directory to sys.path
sys.path.append('/sps/lhcb/yahya/setup/miniconda_install/envs/bd2ksteeEnv/lib/python3.11/site-packages')

# Load the data from the pickle file
file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/Bu2Kpipiee_all_years_MC.pkl"
#file_path = "/sps/lhcb/yahya/Internship/test_Bu/pkl/Bu2Kpipiee_all_years_MC_new_root_files_no_q2_cut.pkl"
#file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeSignal_all_years_MC.pkl"
#file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0Jpsi2eeSS_all_years_Data.pkl"
data = pd.read_pickle(file_path)
#data = data[data['XGBOutput']>0.995]
# data = data[data["q2_NB"] > 14.0]
# Define the observable space for each angle variable
ctl = zfit.Space('ctl', limits=(-1.0, 1.0))
ctk = zfit.Space('ctk', limits=(-1.0, 1.0))
phi = zfit.Space('phi', limits=(-np.pi, np.pi))

# Concatenate the angle variables to create the overall observable space
obs = ctl * ctk * phi

# Convert the angle data into Zfit data format
data_ = data[['ctl_lhcb', 'ctk_lhcb', 'phi_lhcb']].to_numpy()
data_zfit = zfit.data.Data.from_numpy(obs=obs, array=data_, weights=data['kpipi_weight']*data['w_Track_Trig_Kin'])

# Build the Angular3D PDF
angular_pdf, _ = Angular3D(angles=obs, pdfID="not")

# Create the loss function
nll = zfit.loss.UnbinnedNLL(model=angular_pdf, data=data_zfit)

# Create a minimizer
minimizer = zfit.minimize.Minuit()

# Perform the minimization
result = minimizer.minimize(nll)
print(result)
# Plot the result
#b2kstll.plot.plot_distributions_using_partial_integrals(angular_pdf, data_zfit, suffix="angular_fit_")
b2kstll.plot.plot_distributions(result,suffix="angular_fit_part_rec_sps_")
