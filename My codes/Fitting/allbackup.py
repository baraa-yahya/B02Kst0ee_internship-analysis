import numpy as np
import pandas as pd
import zfit
import matplotlib.pyplot as plt
import b2kstll
from b2kstll.models import pdfs
from b2kstll.models.pdfs import *
from b2kstll.utils import *
import mplhep
from scipy.integrate import simps

# Load the actual data
file_path_actual = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0Jpsi2ee_all_years_Data.pkl"
data_file = pd.read_pickle(file_path_actual)
XGBOutput_value = 0.9
data_actual = data_file[data_file['XGBOutput']>XGBOutput_value]
mass_data_actual = data_actual['B_DTF_PV_B_M_0']
# obs_actual = zfit.Space('mass', limits=(mass_data_actual.min(), mass_data_actual.max()))
obs_actual = zfit.Space('mass', limits=(mass_data_actual.min(),6500))
data_actual_zfit = zfit.data.Data.from_numpy(obs=obs_actual, array=mass_data_actual.values)
# Perform the fits for the three components (signal, part_rec, comb_bkg)

# Fit for part_rec background
file_path_part_rec = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/Bu2Kpipiee_all_years_MC.pkl"
data_part_rec = pd.read_pickle(file_path_part_rec)
data_part_rec = data_part_rec[data_part_rec['XGBOutput']>XGBOutput_value]
mass_data_part_rec = data_part_rec['B_DTF_PV_B_M_0'] 
obs_part_rec = zfit.Space('mass', limits=(mass_data_part_rec.min(), mass_data_part_rec.max()))
data_part_rec_zfit = zfit.data.Data.from_numpy(obs=obs_actual, array=mass_data_part_rec.values, weights=data_part_rec['kpipi_weight']*data_part_rec['w_Track_Trig_Kin'])

params_part_rec = pdfs.paramExpStep("part_rec_pdf_prefix")
expstep_pdf_part_rec = pdfs.ExpStep(obs_actual, "part_rec_pdf_suffix", params_part_rec)
nll_part_rec = zfit.loss.UnbinnedNLL(model=expstep_pdf_part_rec, data=data_part_rec_zfit)
minimizer_part_rec = zfit.minimize.Minuit()
result_part_rec = minimizer_part_rec.minimize(nll_part_rec)

# Fit for combinatorial background (same_sign)
file_path_comb_bkg = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0Jpsi2eeSS_all_years_Data.pkl"
data_comb_bkg = pd.read_pickle(file_path_comb_bkg)
data_comb_bkg = data_comb_bkg[data_comb_bkg['XGBOutput']>0.6]
mass_data_comb_bkg = data_comb_bkg['B_DTF_PV_B_M_0']
obs_comb_bkg = zfit.Space('mass', limits=(mass_data_comb_bkg.min(), 6500))
data_comb_bkg_zfit = zfit.data.Data.from_numpy(obs=obs_actual, array=mass_data_comb_bkg.values)

params_comb_bkg = pdfs.paramExpStep("comb_bkg_pdf_prefix")
expstep_pdf_comb_bkg = pdfs.ExpStep(obs_actual, "comb_bkg_pdf_suffix", params_comb_bkg)
nll_comb_bkg = zfit.loss.UnbinnedNLL(model=expstep_pdf_comb_bkg, data=data_comb_bkg_zfit)
minimizer_comb_bkg = zfit.minimize.Minuit()
result_comb_bkg = minimizer_comb_bkg.minimize(nll_comb_bkg)

# Fit for signal
file_path_signal = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeSignal_all_years_MC.pkl"
data_signal = pd.read_pickle(file_path_signal)
data_signal = data_signal[data_signal['XGBOutput']>XGBOutput_value]
mass_data_signal = data_signal['B_DTF_PV_B_M_0']
obs_signal = zfit.Space('mass', limits=(mass_data_signal.min(), mass_data_signal.max()))
data_signal_zfit = zfit.data.Data.from_numpy(obs=obs_actual, array=mass_data_signal.values, weights=data_signal['w_Track_Trig_Kin'])

params_signal = pdfs.paramDCBDGauss("signal_pdf_prefix")
dcb_dg_pdf_signal = pdfs.doubleCBDGauss(obs_actual, "signal_pdf_suffix", params_signal)
nll_signal = zfit.loss.UnbinnedNLL(model=dcb_dg_pdf_signal, data=data_signal_zfit)
minimizer_signal = zfit.minimize.Minuit()
result_signal = minimizer_signal.minimize(nll_signal)

# Combine the fitted functions and fix the parameters

# Fix all parameters of each component
for component_params in [params_part_rec, params_comb_bkg, params_signal]:
    for param in component_params.values():
        param.floating = False

yield_rec_bkg = zfit.Parameter("yield_rec_bkg", 2000, 0, 10000, step_size=1)
yield_comb_bkg = zfit.Parameter("yield_comb_bkg", 6000, 0, 10000, step_size=1)
yield_signal = zfit.Parameter("yield_signal", 2000, 0, 10000, step_size=1)


# Extract the fitted PDFs
part_rec_pdf = expstep_pdf_part_rec.create_extended(yield_rec_bkg)
comb_bkg_pdf = expstep_pdf_comb_bkg.create_extended(yield_comb_bkg)
signal_pdf = dcb_dg_pdf_signal.create_extended(yield_signal)


# Define the total PDF as a sum of the individual PDFs with floating fractions
total_pdf = zfit.pdf.SumPDF([part_rec_pdf, comb_bkg_pdf, signal_pdf ])

# Create the loss function for the total PDF
nll_total = zfit.loss.ExtendedUnbinnedNLL(model=total_pdf, data=data_actual_zfit)

# Create a minimizer and perform the minimization
minimizer_total = zfit.minimize.Minuit()
result_total = minimizer_total.minimize(nll_total)

b2kstll.plot.plot_distributions(result_total, suffix="total_fit_", with_blind=True, labels=['Part_rec bkg', 'Comb bkg', 'Signal'] )

print(f"Valid: {result_total.valid} Converged: {result_total.converged}")

# Compute the errors using the hesse method
param_errors = result_total.hesse()

# Extracting the yields and their errors
yield_rec_bkg_value = yield_rec_bkg.numpy()
yield_rec_bkg_error = param_errors[yield_rec_bkg]['error']

yield_comb_bkg_value = yield_comb_bkg.numpy()
yield_comb_bkg_error = param_errors[yield_comb_bkg]['error']

yield_signal_value = yield_signal.numpy()
yield_signal_error = param_errors[yield_signal]['error']

print(f"Part_rec yield: {yield_rec_bkg_value} ± {yield_rec_bkg_error}")
print(f"Comb_bkg yield: {yield_comb_bkg_value} ± {yield_comb_bkg_error}")
