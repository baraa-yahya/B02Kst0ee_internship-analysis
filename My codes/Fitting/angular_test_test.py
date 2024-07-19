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
import b2kstll
import subprocess

# Add the site-packages directory to sys.path
sys.path.append('/sps/lhcb/yahya/setup/miniconda_install/envs/bd2ksteeEnv/lib/python3.11/site-packages')

# Load the data from the pickle file
file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/Bu2Kpipiee_all_years_MC.pkl"
data = pd.read_pickle(file_path)

# Define the observable space for each angle variable
ctl = zfit.Space('ctl', limits=(-1.0, 1.0))
ctk = zfit.Space('ctk', limits=(-1.0, 1.0))
phi = zfit.Space('phi', limits=(-np.pi, np.pi))

# Concatenate the angle variables to create the overall observable space
obs = ctl * ctk * phi

# Loop over different XGBOutput thresholds
xgb_thresholds = [0.1, 0.2, 0.9, 0.99, 0.999]

for threshold in xgb_thresholds:
    print(f"Processing for XGBOutput threshold: {threshold}")
    
    # Filter data based on the current threshold
    filtered_data = data[data['XGBOutput'] > threshold]

    # Convert the angle data into Zfit data format
    data_ = filtered_data[['ctl_lhcb', 'ctk_lhcb', 'phi_lhcb']].to_numpy()
    data_zfit = zfit.data.Data.from_numpy(obs=obs, array=data_)

    # Build the Angular3D PDF
    angular_pdf, _ = Angular3D(angles=obs, pdfID="Flat")

    # Create the loss function
    nll = zfit.loss.UnbinnedNLL(model=angular_pdf, data=data_zfit)

    # Create a minimizer
    minimizer = zfit.minimize.Minuit()

    # Perform the minimization
    result = minimizer.minimize(nll)

    print(f"""{result}""")
    # Plot the result
    suffix = f"angular_fit_XGB_{threshold}_"
    b2kstll.plot.plot_distributions(result, suffix=suffix)
