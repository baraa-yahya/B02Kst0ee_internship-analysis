import numpy as np
import zfit
import pandas as pd
import b2kstll
import matplotlib.pyplot as plt

# Disable LaTeX rendering in Matplotlib
plt.rcParams['text.usetex'] = False

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

# Check data range
phi_min, phi_max = phi_data.min(), phi_data.max()
print(f"Data range: min={phi_min}, max={phi_max}")

# Define parameters
c1phi = zfit.Parameter('c1phi', 0.1, -2.0, 2.0)
c2phi = zfit.Parameter('c2phi', 0.1, -2.0, 2.0)
c3phi = zfit.Parameter('c3phi', 0.0, -2.0, 2.0)

# Print initial parameter values
print("Initial parameter values:")
print(f"c1phi: {c1phi.numpy()}, c2phi: {c2phi.numpy()}, c3phi: {c3phi.numpy()}")

# Define observable space
phi = zfit.Space('phi', limits=(-np.pi, np.pi))

# Define a simpler PDF for testing
class SimpleLegendrePDF(zfit.pdf.ZPDF):
    _N_OBS = 1
    _PARAMS = ['c1phi', 'c2phi', 'c3phi']

    @zfit.supports()
    def _unnormalized_pdf(self, x, params):
        x0 = x[0]
        c1phi = params['c1phi']
        c2phi = params['c2phi']
        c3phi = params['c3phi']
        return 1. + c1phi * x0 + c2phi * 0.5 * (3. * x0**2 - 1) + \
               c3phi * (5. * x0**3 - 3 * x0) / 2.

pdf_phi = SimpleLegendrePDF(obs=phi, c1phi=c1phi, c2phi=c2phi, c3phi=c3phi)

# Convert data to zfit format
data_phi = zfit.data.Data.from_numpy(obs=phi, array=phi_data.values)

# Evaluate the PDF with initial parameters
pdf_values = pdf_phi.pdf(data_phi)
print("Initial PDF values sample:", pdf_values.numpy()[:10])  # Print a sample of PDF values

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

# Perform minimization
result_phi = minimizer_phi.minimize(nll_phi)
print(result_phi)

# Plotting results (assuming `plot_distributions` handles LaTeX correctly)
b2kstll.plot.plot_distributions(result_phi, suffix="phi_lhcb_legendre")
