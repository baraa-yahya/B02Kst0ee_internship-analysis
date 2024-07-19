import numpy as np
import zfit
import pandas as pd
import b2kstll
from b2kstll.models import legendre
from b2kstll.models.legendre import *
import matplotlib.pyplot as plt
# Load your data from pickle file (adjust path as necessary)
file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeFlat_all_years_MC.pkl"
# file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeFlat_all_years_without_q2_MC.pkl"
data = pd.read_pickle(file_path)
data = data[data['XGBOutput'] > 0.9]
#for ctl

c1L = zfit.Parameter('c1L', 0.1, -1., 1.)
c2L = zfit.Parameter('c2L', 0.1, -1., 1.)
c3L = zfit.Parameter('c3L', 0.1, -1., 1.)
c4L = zfit.Parameter('c4L', 0.1, -1., 1.)
c5L = zfit.Parameter('c5L', 0.1, -1., 1.)
c6L = zfit.Parameter('c6L', 0.1, -1., 1.)
c7L = zfit.Parameter('c7L', 0.1, -1., 1.)
c8L = zfit.Parameter('c8L', 0.1, -1., 1.)
print("Initial parameter values:")
print(f"c1L: {c1L.numpy()}, c2L: {c2L.numpy()}, c3L: {c3L.numpy()}, c4L: {c4L.numpy()}, c5L: {c5L.numpy()}, c6L: {c6L.numpy()}, c7L: {c7L.numpy()}, c8L: {c8L.numpy()}")


ctl_lhcb_data = data['ctl_lhcb']
obs_ctl = zfit.Space('ctl', limits=(-1, 1))
pdf_ctl = legendre.LegendrePolynomialPDF_8(obs=obs_ctl, c1L=c1L, c2L=c2L, c3L=c3L, c4L=c4L, c5L=c5L, c6L=c6L, c7L=c7L, c8L=c8L)

data_ctl_lhcb = zfit.data.Data.from_numpy(obs=obs_ctl, array=ctl_lhcb_data.values)
print(ctl_lhcb_data.count())# Create negative log-likelihood
nll_ctl = zfit.loss.UnbinnedNLL(model=pdf_ctl, data=data_ctl_lhcb)

# Create minimizer
minimizer_ctl = zfit.minimize.Minuit()

# Perform minimization
result_ctl = minimizer_ctl.minimize(nll_ctl)

# b2kstll.plot.plot_distributions(result_ctl, suffix="ctl")



ctk_data = data['ctk_lhcb'].values

c1k = zfit.Parameter('c1k', 0.1, -1., 1.)
c2k = zfit.Parameter('c2k', 0.1, -1., 1.)
c3k = zfit.Parameter('c3k', 0.1, -1., 1.)
c4k = zfit.Parameter('c4k', 0.1, -1., 1.)
c5k = zfit.Parameter('c5k', 0.1, -1., 1.)
c6k = zfit.Parameter('c6k', 0.1, -1., 1.)
c7k = zfit.Parameter('c7k', 0.1, -1., 1.)
c8k = zfit.Parameter('c8k', 0.1, -1., 1.)
ctk = zfit.Space('ctk', limits=(-1, 1))

pdf_ctk = legendre.LegendrePolynomialPDF_8(obs=ctk, c1L=c1k, c2L=c2k, c3L=c3k, c4L=c4k, c5L=c5k, c6L=c6k, c7L=c7k, c8L=c8k)

data_k = zfit.data.Data.from_numpy(obs=ctk, array=ctk_data)

# Create negative log-likelihood
nll_ctk = zfit.loss.UnbinnedNLL(model=pdf_ctk, data=data_k)

# Create minimizer
minimizer_ctk = zfit.minimize.Minuit()

# Perform minimization
result_ctk = minimizer_ctk.minimize(nll_ctk)

# Note: If you want to plot using b2kstll, ensure you have imported it correctly and adjust as necessary
# b2kstll.plot.plot_distributions(result_ctk, suffix="ctk_lhcb_legendre") #labels=['Legendre Polynomial Order 6'])

phi_data = data['phi_lhcb']


# Define parameters
c1phi = zfit.Parameter('c1phi', 0.0, -1.0, 1.0)
c2phi = zfit.Parameter('c2phi', 0.1, -1.0, 1.0)
c3phi = zfit.Parameter('c3phi', 0.0, -1.0, 1.0)
c4phi = zfit.Parameter('c4phi', 0.1, -1.0, 1.0)
c5phi = zfit.Parameter('c5phi', 0.0, -1.0, 1.0)
c6phi = zfit.Parameter('c6phi', 0.1, -1.0, 1.0)
c7phi = zfit.Parameter('c7phi', 0.0, -1.0, 1.0)
c8phi = zfit.Parameter('c8phi', 0.1, -1.0, 1.0)


# Define observable space
phi = zfit.Space('phi', limits=(-np.pi, np.pi))

# Define PDF
pdf_phi = legendre.LegendrePolynomialPDF_8(obs=phi, c1L=c1phi, c2L=c2phi, c3L=c3phi, c4L=c4phi, c5L=c5phi, c6L=c6phi, c7L=c7phi, c8L=c8phi)

# Convert data to zfit format
data_phi = zfit.data.Data.from_numpy(obs=phi, array=phi_data.values)

# Evaluate the PDF with initial parameters
pdf_values = pdf_phi.pdf(data_phi)

# Create negative log-likelihood
nll_phi = zfit.loss.UnbinnedNLL(model=pdf_phi, data=data_phi)

# Evaluate the NLL with initial parameters
initial_loss = nll_phi.value()

# Create minimizer
minimizer_phi = zfit.minimize.Minuit()

# Perform minimization (without callback)
result_phi = minimizer_phi.minimize(nll_phi)

print(result_ctl)
print(result_ctk)
print(result_phi)

# Plotting results
# b2kstll.plot.plot_distributions(result_phi, suffix="phi_lhcb_legendre")


# Get data points from your DataFrame
ctl_lhcb_data = data['ctl_lhcb'].values
ctk_lhcb_data = data['ctk_lhcb'].values
phi_lhcb_data = data['phi_lhcb'].values

# Evaluate PDFs at data points (these are unnormalized probabilities)
pdf_ctl_values = pdf_ctl.pdf(data_ctl_lhcb)
pdf_ctk_values = pdf_ctk.pdf(data_k)
pdf_phi_values = pdf_phi.pdf(data_phi)

# Print or use the evaluated PDF values (these are unnormalized)
print("Unnormalized PDF values for ctl_lhcb_data:")
print(pdf_ctl_values)

print("Unnormalized PDF values for ctk_lhcb_data:")
print(pdf_ctk_values)

print("Unnormalized PDF values for phi_lhcb_data:")
print(pdf_phi_values)

# Calculate acceptance (likelihood without normalization)
acceptance = 1 / (pdf_ctl_values * pdf_ctk_values * pdf_phi_values)

# Print or use the acceptance value
print(f"Acceptance value: {acceptance}")
