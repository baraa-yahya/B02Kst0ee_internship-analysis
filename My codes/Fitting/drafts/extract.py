import numpy as np
import pickle
import zfit
import matplotlib.pyplot as plt

# Load the PKL file
with open('/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeSignal_2012_MC.pkl', 'rb') as f:
    data = pickle.load(f)

# Convert the column "B_DTF_PV_B_M_0" to a numpy array
b_mass_array = data['B_DTF_PV_B_M_0'].to_numpy()

# Create a histogram in the specified range
B_mass, bins = np.histogram(b_mass_array, bins=100, range=(4500, 5750), density=True)

# Calculate the bin centers
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Convert the histogram data to zfit data
obs = zfit.Space('B_mass', limits=(4500, 5750))
data_zfit = zfit.Data.from_numpy(obs=obs, array=b_mass_array)

# Define the model parameters with mu shifted to the left
mu = zfit.Parameter('mu', 5150, 4800, 5250)  # Adjusted initial value and limits to shift the peak leftward
sigma = zfit.Parameter('sigma', 10, 0, 600)
alphaL = zfit.Parameter('alphaL', 1.0, 0, 20)
alphaR = zfit.Parameter('alphaR', 1.0, 0, 20)
nL = zfit.Parameter('nL', 5, 0, 200)
nR = zfit.Parameter('nR', 5, 0, 200)

# Define the CrystalBall PDFs with adjusted mu
crystalball_L = zfit.pdf.CrystalBall(obs=obs, mu=mu, sigma=sigma, alpha=alphaL, n=nL)
crystalball_R = zfit.pdf.CrystalBall(obs=obs, mu=mu, sigma=sigma, alpha=alphaR, n=nR)

# Define the fraction parameter and the combined model
fraction = zfit.Parameter('fraction', 0.5, 0, 1)
model = zfit.pdf.SumPDF([crystalball_L, crystalball_R], [fraction])


# Create the negative log likelihood
nll = zfit.loss.UnbinnedNLL(model=model, data=data_zfit)

# Create a minimizer
minimizer = zfit.minimize.Minuit()

# Perform the minimization
result = minimizer.minimize(nll)

# Get the fitted parameters
fitted_parameters = result.params

# Print the fitted parameters with their errors
print("Fitted parameters:")
for param in result.params:
    value = result.params[param]['value']
    try:
        error = result.params[param]['minuit_hesse']['error']
    except KeyError:
        error = "Error not available"
    print(f"{param.name}: {value} ± {error}")

# Generate points to plot the fitted model
x = np.linspace(4500, 5750, 1000)
y = model.pdf(x, norm_range=obs)

# Plot the data and the fitted curve
plt.hist(b_mass_array, bins=100, range=(4500, 5750), density=True, alpha=0.6, color='g', label='Data')
plt.plot(x, y, 'r-', label='Fit')
plt.xlabel('B Mass')
plt.ylabel('Probability')
plt.legend()

# Save the plot
filename = 'fit_zfit4.pdf'
plt.savefig(filename)
plt.show()
