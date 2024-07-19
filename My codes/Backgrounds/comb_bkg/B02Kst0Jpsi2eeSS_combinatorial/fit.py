import ROOT
import matplotlib.pyplot as plt
import mplhep
import numpy as np
import zfit
import zfit.z.numpy as znp
from zfit import z  

obs = zfit.Space('x', limits=(4500,7000 ))


## Create numpy data set
mu_true = 0
sigma_true = 5
data_np= "/sps/lhcb/pietrzyk/analysis/samples/Data/B02Kst0Jpsi2eeSS/2011/Flagged/B02Kst0Jpsi2eeSS-Data-2011-MagDown-Translated-Flagged-XGBOutput-q2BDTOutput.root"
#data_np = np.random.normal(
 #   mu_true, 
  #  sigma_true, 
   # size=10000
#)


## Put this data set into a zfit data set
data = zfit.data.Data.from_numpy(
    obs=obs, 
    array=data_np
)

## Create PDF
mu = zfit.Parameter(
    "mu", 
    2.4, 
    -1., 
    5.,
    10, 
) 
sigma = zfit.Parameter(
    "sigma",
    1.3, 
    0, 
    5.,
    10,
) 
gauss = zfit.pdf.Gauss(
    obs=obs, 
    mu=mu, 
    sigma=sigma
)

## Define the fitting technique
nll = zfit.loss.UnbinnedNLL(model=gauss, data=data)  # loss
minimizer = zfit.minimize.Minuit()

## Do the fit and get fit results
minimum = minimizer.minimize(loss=nll)

## Print params and fit results
print(minimum)

## Plot
import matplotlib.pyplot as plt

n_bins = 50
range_ = (4500,7000)
_ = plt.hist(data_np, bins=n_bins, range=range_)


## Compute y = pdf(x), for a given binned x axis
## By default this is normalized to one.
x = np.linspace(*range_, num=1000)
pdf = zfit.run(gauss.pdf(x))

## Plot the pdf
## data_np.shape[0] equivalent to len(data_np)
plt.plot(x, data_np.shape[0] / n_bins * (range_[1] - range_[0]) * pdf)
plt.savefig("easy_tutorial.pdf")