import numpy as np
import pickle
import zfit
import matplotlib.pyplot as plt
import os

# Load the PKL file
with open('/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeSignal_2012_MC.pkl', 'rb') as f:
    data = pickle.load(f)

# Access the column "B_DTF_PV_B_M_0" and create a histogram in the specified range
B_mass, bins = np.histogram(data['B_DTF_PV_B_M_0'], bins=100, range=(4500, 5750), density=True)

# Calculate the bin centers
bin_centers = 0.5 * (bins[:-1] + bins[1:])

# Convert the histogram data to zfit data
obs = zfit.Space('B_mass', limits=(4500, 5750))
data_zfit = zfit.Data.from_numpy(obs=obs, array=bin_centers)

# Define the model
mu = zfit.Parameter('mu', 5250, 4800, 5250)
sigma = zfit.Parameter('sigma', 15, 0, 50)
alphaL = zfit.Parameter('alphaL', 1, 0.5, 2)
alphaR = zfit.Parameter('alphaR', 1, 0.5, 2)
nL = zfit.Parameter('nL', 1, 1, 10)
nR = zfit.Parameter('nR', 1, 1, 10)

crystalball_L = zfit.pdf.CrystalBall(obs=obs, mu=mu, sigma=sigma, alpha=alphaL, n=nL)
crystalball_R = zfit.pdf.CrystalBall(obs=obs, mu=mu, sigma=sigma, alpha=alphaR, n=nR)

fraction = zfit.Parameter('fraction', 0.2, 0, 1)
model = zfit.pdf.SumPDF([crystalball_L, crystalball_R], [fraction])

# Create the negative log likelihood
nll = zfit.loss.UnbinnedNLL(model=model, data=data_zfit)

# Create a minimizer
minimizer = zfit.minimize.Minuit()

# Perform the minimization
result = minimizer.minimize(nll)

# Get the fitted parameters
fitted_parameters = result.params

# Plot the data and the fitted curve
x = np.linspace(4500, 5750, 1000)
y = model.pdf(x, norm_range=obs)

plt.hist(data['B_DTF_PV_B_M_0'], bins=100, range=(4500, 5750), density=True, alpha=0.6, color='g', label='Data')
plt.plot(x, y, 'r-', label='Fit')
plt.xlabel('B Mass')
plt.ylabel('Probability')
plt.legend()

# Save the plot
filename = 'fit_zfi2.pdf'
plt.savefig(filename)
plt.show()
