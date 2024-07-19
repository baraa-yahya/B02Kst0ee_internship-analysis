import uproot
import zfit
import numpy as np
import matplotlib.pyplot as plt
from zfit import z

obs = zfit.Space('x', limits=(4500, 70000))
data_root = zfit.Data.from_root(path="/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2011/Flagged/B02Kst0Jpsi2eeSS-Data-2011-MagDown-Translated-Flagged-XGBOutput-q2BDTOutput.root", treepath="DecayTree", branches="B_M", obs=obs)
mu = zfit.Parameter("mu", 0, -10, 10)
sigma = zfit.Parameter("sigma", 1, 0, 10)
gauss = zfit.pdf.Gauss(mu=mu, sigma=sigma, obs=data_root)
nll = zfit.loss.UnbinnedNLL(model=gauss, data=data_root)
minimizer = zfit.minimize.Minuit()

result = minimizer.minimize(nll)

plt.hist(data[:, 0], bins=50, density=True, alpha=0.5, label='Data')

x = np.linspace(4500, 70000, 1000)
y = zfit.run(gauss.pdf(x))
plt.plot(x, y, color='red', lw=2, label='Fit')

plt.xlabel('Variable')
plt.ylabel('Normalized Count')
plt.legend()

# Save the plot as a PDF
plt.savefig("mass_fit.pdf")
