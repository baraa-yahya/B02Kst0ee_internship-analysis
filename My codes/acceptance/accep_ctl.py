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
print(result_ctl)

b2kstll.plot.plot_distributions(result_ctl, suffix="ctl")
