import uproot
import zfit
import numpy as np
import matplotlib.pyplot as plt


# Load the ROOT file
file = uproot.open("/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2011/Flagged/B02Kst0Jpsi2eeSS-Data-2011-MagDown-Translated-Flagged-XGBOutput-q2BDTOutput.root")

# Access the DecayTree TTree
tree = file["DecayTree"]

# Extract the B_M data
data_df = tree.arrays(tree.keys(),library="pd")
B_M_obs = zfit.Space('B_M', limits=(4500, 7000))
print(data_df.head())
# Convert the B_M_data to zfit Data
data = zfit.Data.from_pandas(data_df, obs=B_M_obs )

# Define the Gaussian model
mean = zfit.Parameter("mean", 5300, 5200, 5400)
sigma = zfit.Parameter("sigma", 20, 10, 50)
gauss = zfit.pdf.Gauss(obs=B_M_obs, mu=mean, sigma=sigma)

# Create the likelihood function
likelihood = zfit.loss.UnbinnedNLL(model=gauss, data=data)

# Minimize the likelihood
minimizer = zfit.minimize.Minuit()
result = minimizer.minimize(likelihood)

# Print the result
print(result)

# Plot the fit result
plt.hist(data_df["B_M"], bins=100, range=(4500, 7000), density=True, alpha=0.5, label='Data')
x = np.linspace(5000, 6000, 1000)
plt.plot(x, gauss.pdf(x).numpy(), label='Fit')
plt.xlabel('B_M')
plt.ylabel('Normalized counts')
plt.legend()
plt.show()
plt.savefig("mass_fit3.pdf")
