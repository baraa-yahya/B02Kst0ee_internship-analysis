import numpy as np
import pickle
import zfit
import matplotlib.pyplot as plt
import sys

# Check if all required arguments are provided
if len(sys.argv) < 4:
    print("Usage: python fit_bkg.py <decay_name> <year> <mode>")
    sys.exit(1)

# Extract arguments
decay_name = sys.argv[1]
year = sys.argv[2]
mode = sys.argv[3]

# Load the PKL file based on decay, year, and mode
pkl_filename = f"/sps/lhcb/yahya/Internship/pkl_creation/pkl/{decay_name}_{year}_{mode}.pkl"
with open(pkl_filename, 'rb') as f:
    data = pickle.load(f)

# Convert the column "B_DTF_PV_B_M_0" to a numpy array
b_mass_array = data['B_DTF_PV_B_M_0'].to_numpy()

# Create a histogram in the specified range
B_mass, bins = np.histogram(b_mass_array, bins=100, range=(3600, 7000), density=True)

# Calculate the bin centers
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Convert the histogram data to zfit data
if mode == 'MC':
    obs = zfit.Space('B_mass', limits=(4500, 5750))
else:
    obs = zfit.Space('B_mass', limits=(3600, 7000))
data_zfit = zfit.Data.from_numpy(obs=obs, array=b_mass_array)

# Define the model parameters
if mode == 'MC':
    mu = zfit.Parameter('mu', 5150, 4800, 5250)
    sigma = zfit.Parameter('sigma', 10, 0, 50)
    alphaL = zfit.Parameter('alphaL', 1, 0.5, 5)
    alphaR = zfit.Parameter('alphaR', 6, 2, 10)
    nL = zfit.Parameter('nL', 1, 0, 10)
    nR = zfit.Parameter('nR', 2, 0, 10)
    crystalball_L = zfit.pdf.CrystalBall(obs=obs, mu=mu, sigma=sigma, alpha=alphaL, n=nL)
    crystalball_R = zfit.pdf.CrystalBall(obs=obs, mu=mu, sigma=sigma, alpha=alphaR, n=nR)
    fraction = zfit.Parameter('fraction', 0.5, 0, 1)
    model = zfit.pdf.SumPDF([crystalball_L, crystalball_R], [fraction])
else:
    lambda_param = zfit.Parameter('lambda_param', -0.1, -1.0, 0.0)
    exp_model = zfit.pdf.Exponential(lambda_param, obs=obs)
    model = exp_model

# Create the negative log likelihood
nll = zfit.loss.UnbinnedNLL(model=model, data=data_zfit)

# Create a minimizer
minimizer = zfit.minimize.Minuit()

# Perform the minimization
result = minimizer.minimize(nll)

# Perform HESSE to estimate errors
result.hesse()

# Get the fitted parameters
fitted_parameters = result.params

# Generate filename based on decay, year, and mode
filename = f'{decay_name}_{year}_{mode}_fit.pdf'

# Print the fitted parameters with their errors
print("Fitted parameters:")
for param in fitted_parameters:
    value = fitted_parameters[param]['value']
    try:
        error = fitted_parameters[param]['hesse']['error']
    except KeyError:
        error = "Error not available"
    print(f"{param.name}: {value:.5f} ± {error}")

# Generate points to plot the fitted model
x = np.linspace(3600, 7000, 1000)
y = model.pdf(x, norm_range=obs).numpy()

# Plot the data and the fitted curve
plt.hist(b_mass_array, bins=100, range=(3600, 7000), density=True, alpha=0.6, color='g', label='Data')
plt.plot(x, y, 'r-', label='Fit')
plt.xlabel('B Mass')
plt.ylabel('Probability')
plt.legend()

# Save the plot with the dynamically generated filename
plt.savefig(filename)
# to run the code python name_of_file decay_name year mode==Data or MC
#this code for same sign background and signal decay both mass dist. fit