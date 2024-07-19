import numpy as np
import zfit
import pandas as pd
import b2kstll
from b2kstll.models import legendre
from b2kstll.models.legendre import *


# Load your data from pickle file (adjust path as necessary)
file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeFlat_all_years_MC.pkl"
data = pd.read_pickle(file_path)
data = data[data['XGBOutput'] > 0.9]
#for ctk

ctk_data = data['ctk_lhcb'].values

c1k = zfit.Parameter('c1k', 0.1, -1., 1.)
c2k = zfit.Parameter('c2k', 0.1, -1., 1.)
c3k = zfit.Parameter('c3k', 0.1, -1., 1.)
c4k = zfit.Parameter('c4k', 0.1, -1., 1.)
ctk = zfit.Space('ctk', limits=(-1, 1))

pdf_ctk = legendre.LegendrePolynomialPDF_4(obs=ctk, c1L=c1k, c2L=c2k, c3L=c3k, c4L=c4k)

data_k = zfit.data.Data.from_numpy(obs=ctk, array=ctk_data)

# Create negative log-likelihood
nll_ctk = zfit.loss.UnbinnedNLL(model=pdf_ctk, data=data_k)

# Create minimizer
minimizer_ctk = zfit.minimize.Minuit()

# Perform minimization
result_ctk = minimizer_ctk.minimize(nll_ctk)
print(result_ctk)

# Note: If you want to plot using b2kstll, ensure you have imported it correctly and adjust as necessary
b2kstll.plot.plot_distributions(result_ctk, suffix="ctk_lhcb_legendre") #labels=['Legendre Polynomial Order 6'])

