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
import subprocess
from b2kstll.models.Kpi import KpiPWavePDF, KpiSWavePDF, ReKpiPandSWavePDF, ImKpiPandSWavePDF
from b2kstll.models import angular
from b2kstll.models import spherical_harmonics

from b2kstll.models.acceptance import AccPDF
# Add the site-packages directory to sys.path
sys.path.append('/sps/lhcb/yahya/setup/miniconda_install/envs/bd2ksteeEnv/lib/python3.11/site-packages')

# Load the data from the pickle file
file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeSignal_all_years_MC.pkl"
data = pd.read_pickle(file_path)

# Define the observable space for each angle variable
ctl = zfit.Space('ctl', limits=(-1.0, 1.0))
ctk = zfit.Space('ctk', limits=(-1.0, 1.0))
phi = zfit.Space('phi', limits=(-np.pi, np.pi))

# Concatenate the angle variables to create the overall observable space
obs = ctl * ctk * phi

# Convert the angle data into Zfit data format
data_ = data[['ctl_lhcb', 'ctk_lhcb', 'phi_lhcb']].to_numpy()
data_zfit = zfit.Data.from_numpy(obs=obs, array=data_)

# Define the function for P-wave
def PWave(FL, S3, S4, S5, AFB, S7, S8, S9, costheta_l, costheta_k, phi, backend=tf):
    sintheta_k = backend.sqrt(1.0 - costheta_k * costheta_k)
    sintheta_l = backend.sqrt(1.0 - costheta_l * costheta_l)
    sintheta_2k = (1.0 - costheta_k * costheta_k)
    sintheta_2l = (1.0 - costheta_l * costheta_l)
    sin2theta_k = (2.0 * sintheta_k * costheta_k)
    cos2theta_l = (2.0 * costheta_l * costheta_l - 1.0)
    sin2theta_l = (2.0 * sintheta_l * costheta_l)

    pdf = ((3.0 / 4.0) * (1.0 - FL) * sintheta_2k +
           FL * costheta_k * costheta_k +
           (1.0 / 4.0) * (1.0 - FL) * sintheta_2k * cos2theta_l +
           -1.0 * FL * costheta_k * costheta_k * cos2theta_l +
           S3 * sintheta_2k * sintheta_2l * backend.cos(2.0 * phi) +
           S4 * sin2theta_k * sin2theta_l * backend.cos(phi) +
           S5 * sin2theta_k * sintheta_l * backend.cos(phi) +
           (4.0 / 3.0) * AFB * sintheta_2k * costheta_l +
           S7 * sin2theta_k * sintheta_l * backend.sin(phi) +
           S8 * sin2theta_k * sin2theta_l * backend.sin(phi) +
           S9 * sintheta_2k * sintheta_2l * backend.sin(2.0 * phi))

    return (9.0 / (32.0 * pi)) * pdf

# Define the AngularPDF for P-wave
params = {
    'FL': zfit.Parameter('FL', 0.6, 0, 1),
    'S3': zfit.Parameter('S3', 0.0, -2, 2),
    'S4': zfit.Parameter('S4', 0.0, -2, 2),
    'S5': zfit.Parameter('S5', 0.0, -2, 2),
    'AFB': zfit.Parameter('AFB', 0.0, -2, 2),
    'S7': zfit.Parameter('S7', 0.0, -2, 2),
    'S8': zfit.Parameter('S8', 0.0, -2, 2),
    'S9': zfit.Parameter('S9', 0.0, -2, 2)
}

angular_pdf = angular.AngularPDF(obs=obs, params=params, func=PWave, name='PWave')

# Create the loss function
nll = zfit.loss.UnbinnedNLL(model=angular_pdf, data=data_zfit)

# Create a minimizer
minimizer = zfit.minimize.Minuit()

# Perform the minimization
result = minimizer.minimize(nll)

# Plot the result
b2kstll.plot.plot_distributions(result, suffix="angular_fit_same_sign_")
