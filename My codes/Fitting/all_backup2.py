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
data_actual = data_file[data_file['XGBOutput']>0.9]
mass_data_actual = data_actual['B_DTF_PV_B_M_0']
# obs_actual = zfit.Space('mass', limits=(mass_data_actual.min(), mass_data_actual.max()))
obs_actual = zfit.Space('mass', limits=(mass_data_actual.min(),mass_data_actual.max()))
data_actual_zfit = zfit.data.Data.from_numpy(obs=obs_actual, array=mass_data_actual.values)
# Perform the fits for the three components (signal, part_rec, comb_bkg)

# Fit for part_rec background
file_path_part_rec = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/Bu2Kpipiee_all_years_MC.pkl"
data_part_rec = pd.read_pickle(file_path_part_rec)
mass_data_part_rec = data_part_rec['B_DTF_PV_B_M_0']
obs_part_rec = zfit.Space('mass', limits=(mass_data_part_rec.min(), mass_data_part_rec.max()))
data_part_rec_zfit = zfit.data.Data.from_numpy(obs=obs_actual, array=mass_data_part_rec.values)

params_part_rec = pdfs.paramExpStep("part_rec_pdf_prefix")
expstep_pdf_part_rec = pdfs.ExpStep(obs_actual, "part_rec_pdf_suffix", params_part_rec)
nll_part_rec = zfit.loss.UnbinnedNLL(model=expstep_pdf_part_rec, data=data_part_rec_zfit)
minimizer_part_rec = zfit.minimize.Minuit()
result_part_rec = minimizer_part_rec.minimize(nll_part_rec)

# Fit for combinatorial background (same_sign)
file_path_comb_bkg = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0Jpsi2eeSS_all_years_Data.pkl"
data_comb_bkg = pd.read_pickle(file_path_comb_bkg)
mass_data_comb_bkg = data_comb_bkg['B_DTF_PV_B_M_0']
obs_comb_bkg = zfit.Space('mass', limits=(mass_data_comb_bkg.min(), mass_data_comb_bkg.max()))
data_comb_bkg_zfit = zfit.data.Data.from_numpy(obs=obs_actual, array=mass_data_comb_bkg.values)

params_comb_bkg = pdfs.paramExpStep("comb_bkg_pdf_prefix")
expstep_pdf_comb_bkg = pdfs.ExpStep(obs_actual, "comb_bkg_pdf_suffix", params_comb_bkg)
nll_comb_bkg = zfit.loss.UnbinnedNLL(model=expstep_pdf_comb_bkg, data=data_comb_bkg_zfit)
minimizer_comb_bkg = zfit.minimize.Minuit()
result_comb_bkg = minimizer_comb_bkg.minimize(nll_comb_bkg)

# Fit for signal
file_path_signal = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeSignal_all_years_MC.pkl"
data_signal = pd.read_pickle(file_path_signal)
mass_data_signal = data_signal['B_DTF_PV_B_M_0']
obs_signal = zfit.Space('mass', limits=(mass_data_signal.min(), mass_data_signal.max()))
data_signal_zfit = zfit.data.Data.from_numpy(obs=obs_actual, array=mass_data_signal.values)

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

b2kstll.plot.plot_distributions_from_fitresult(result_total, suffix="total_fit_", with_pull=False, with_blind=True)

def plot_pulls(data_points, total_pdf_values_for_bins, bin_centers, bin_widths, suffix=""):
    plt.style.use(mplhep.style.LHCb2)
    plt.close("all")

    fig = plt.figure(1)
    ax2 = fig.add_axes([0.15, 0.07, 0.78, 0.12])

    ax2.axhline(0., color='r', linewidth=2)
    ax2.axhline(3., color='r', linewidth=2, linestyle=":")
    ax2.axhline(-3., color='r', linewidth=2, linestyle=":")

    finite = data_points > 0
    res = (data_points - total_pdf_values_for_bins) / np.sqrt(data_points)

    ax2.errorbar(bin_centers, res, xerr=bin_widths, yerr=np.ones(len(res)),
                 color="k", fmt="o", markersize=5, linewidth=2)

    ax2.set_ylim(-5., 5)
    ax2.set_xlim(bin_centers[0] - bin_widths[0], bin_centers[-1] + bin_widths[-1])

    plt.savefig(suffix + "pull_plot.pdf")
    plt.show()

    plt.clf()
    del ax2
    plt.close(fig)

# Define the bins
bins = np.linspace(mass_data_actual.min(), mass_data_actual.max(), 50)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_widths = (bins[1:] - bins[:-1]) / 2

# Evaluate the total PDF at the middle of each bin
total_pdf_values = total_pdf.pdf(bin_centers).numpy()

# Bin the actual data to get the counts in each bin
data_counts, bin_edges = np.histogram(mass_data_actual, bins=bins)

# Calculate the total PDF values for each bin
total_pdf_values_for_bins = total_pdf_values * (len(mass_data_actual) * (bins[1] - bins[0]))
print(total_pdf_values)
# Apply mask to filter out zero-count bins
mask = data_counts > 0
data_points = data_counts[mask]
total_pdf_values_for_bins = total_pdf_values_for_bins[mask]
bin_centers_filtered = bin_centers[mask]
bin_widths_filtered = bin_widths[mask]

# Plot the pulls
plot_pulls(data_points=data_points, 
           total_pdf_values_for_bins=total_pdf_values_for_bins, 
           bin_centers=bin_centers_filtered,
           bin_widths=bin_widths_filtered, 
           suffix="")


import fitz

result = fitz.open()
for pdf in ['total_fit_SumPDF-mass.pdf', 'pull_plot.pdf']:
    with fitz.open(pdf) as mfile:
        result.insert_pdf(mfile)
result.save("result.pdf")

print(f""" valid: {result_total.valid} converged: {result_total.converged}""")
print(f""" part_rec yield: {yield_rec_bkg} comb_bkg yield: {yield_comb_bkg}""")