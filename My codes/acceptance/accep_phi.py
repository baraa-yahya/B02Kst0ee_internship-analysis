import numpy as np
import zfit
import pandas as pd
import b2kstll
from b2kstll.models.legendre import LegendrePolynomialPDF_8

# Load your data from pickle file (adjust path as necessary)
file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeFlat_all_years_MC.pkl"
try:
    data = pd.read_pickle(file_path)
except FileNotFoundError:
    print(f"File not found: {file_path}")
    raise

data = data[data['XGBOutput'] > 0.9]
phi_data = data['phi_lhcb']

if phi_data.isna().any():
    raise ValueError("phi_data contains NaN values.")

# Define parameters
c1phi = zfit.Parameter('c1phi', 0.0, -1.0, 1.0)
c2phi = zfit.Parameter('c2phi', 0.1, -1.0, 1.0)
c3phi = zfit.Parameter('c3phi', 0.0, -1.0, 1.0)
c4phi = zfit.Parameter('c4phi', 0.1, -1.0, 1.0)
c5phi = zfit.Parameter('c5phi', 0.0, -1.0, 1.0)
c6phi = zfit.Parameter('c6phi', 0.1, -1.0, 1.0)
c7phi = zfit.Parameter('c7phi', 0.0, -1.0, 1.0)
c8phi = zfit.Parameter('c8phi', 0.1, -1.0, 1.0)


# Print initial parameter values
print("Initial parameter values:")
print(f"c1phi: {c1phi.numpy()}, c2phi: {c2phi.numpy()}, c3phi: {c3phi.numpy()}, c4phi: {c4phi.numpy()}, c5phi: {c5phi.numpy()}, c6phi: {c6phi.numpy()}, c7phi: {c7phi.numpy()}, c8phi: {c8phi.numpy()}")

# Define observable space
phi = zfit.Space('phi', limits=(-np.pi, np.pi))

# Define PDF
pdf_phi = LegendrePolynomialPDF_8(obs=phi, c1=c1phi, c2=c2phi, c3=c3phi, c4=c4phi, c5=c5phi, c6=c6phi, c7=c7phi, c8=c8phi)

# Convert data to zfit format
data_phi = zfit.data.Data.from_numpy(obs=phi, array=phi_data.values)

# Evaluate the PDF with initial parameters
pdf_values = pdf_phi.pdf(data_phi)
print("Initial PDF values:", pdf_values.numpy())

# Create negative log-likelihood
nll_phi = zfit.loss.UnbinnedNLL(model=pdf_phi, data=data_phi)

# Evaluate the NLL with initial parameters
initial_loss = nll_phi.value()
print("Initial NLL value:", initial_loss.numpy())

# If initial NLL value is NaN, print details and exit
if np.isnan(initial_loss):
    print("Initial NLL value is NaN. Exiting...")
    exit()

# Create minimizer
minimizer_phi = zfit.minimize.Minuit()

# Perform minimization (without callback)
result_phi = minimizer_phi.minimize(nll_phi)
print(result_phi)

# Plotting results
b2kstll.plot.plot_distributions(result_phi, suffix="phi_lhcb_legendre")
