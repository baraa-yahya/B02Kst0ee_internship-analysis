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
data_actual = data_file[data_file['XGBOutput'] > 0.9]
mass_data_actual = data_actual['B_DTF_PV_B_M_0']
obs_actual = zfit.Space('mass', limits=(mass_data_actual.min(), mass_data_actual.max()))
data_actual_zfit = zfit.data.Data.from_numpy(obs=obs_actual, array=mass_data_actual.values)

# Perform the fits for the three components (signal, part_rec, comb_bkg)
# Part_rec background
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

# Combinatorial background (same_sign)
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

# Signal
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
total_pdf = zfit.pdf.SumPDF([part_rec_pdf, comb_bkg_pdf, signal_pdf])

# Create the loss function for the total PDF
nll_total = zfit.loss.ExtendedUnbinnedNLL(model=total_pdf, data=data_actual_zfit)

# Create a minimizer and perform the minimization
minimizer_total = zfit.minimize.Minuit()
result_total = minimizer_total.minimize(nll_total)

# Define the bins
bins = np.linspace(mass_data_actual.min(), mass_data_actual.max(), 100)
bin_centers = (bins[:-1] + bins[1:]) / 2
bin_widths = (bins[1:] - bins[:-1]) / 2

# Evaluate the total PDF at the middle of each bin
total_pdf_values = total_pdf.pdf(bin_centers).numpy()

# Bin the actual data to get the counts in each bin
data_counts, bin_edges = np.histogram(mass_data_actual, bins=bins)

# Calculate the total PDF values for each bin
total_pdf_values_for_bins = total_pdf_values * (len(mass_data_actual) * (bins[1] - bins[0]))

# Apply mask to filter out zero-count bins
mask = data_counts > 0
data_points = data_counts[mask]
total_pdf_values_for_bins = total_pdf_values_for_bins[mask]
bin_centers_filtered = bin_centers[mask]
bin_widths_filtered = bin_widths[mask]

def plot_with_pulls(result, data_points, total_pdf_values_for_bins, bin_centers_filtered, bin_widths_filtered, suffix=""):
    plt.style.use(mplhep.style.LHCb2)
    plt.close("all")

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]})

    # Plot the fit result on ax1
    b2kstll.plot.plot_distributions_from_fitresult(result, ax=ax1, suffix=suffix, with_pull=False, with_blind=True)

    # Plot the pulls on ax2
    ax2.axhline(0., color='r', linewidth=2)
    ax2.axhline(3., color='r', linewidth=2, linestyle=":")
    ax2.axhline(-3., color='r', linewidth=2, linestyle=":")

    res = (data_points - total_pdf_values_for_bins) / np.sqrt(total_pdf_values_for_bins)

    ax2.errorbar(bin_centers_filtered, res, xerr=bin_widths_filtered, yerr=np.ones(len(res)),
                 color="k", fmt="o", markersize=5, linewidth=2)

    ax2.set_ylim(-10., 10)
    ax2.set_xlim(bin_centers_filtered[0] - bin_widths_filtered[0], bin_centers_filtered[-1] + bin_widths_filtered[-1])

    plt.tight_layout()
    plt.savefig(suffix + "combined_fit_and_pull_plot.pdf")
    plt.show()

def plot_style(result, data_points, total_pdf_values_for_bins, bin_centers_filtered, bin_widths_filtered, suffix="", with_pulls=False):
    plt.style.use(mplhep.style.LHCb2)
    plt.close("all")

    if with_pulls:
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), gridspec_kw={'height_ratios': [3, 1]})

        # Plot the fit result on ax1
        b2kstll.plot.plot_distributions_from_fitresult(result, ax=ax1, suffix=suffix, with_pull=False, with_blind=True)

        # Plot the pulls on ax2
        ax2.axhline(0., color='r', linewidth=2)
        ax2.axhline(3., color='r', linewidth=2, linestyle=":")
        ax2.axhline(-3., color='r', linewidth=2, linestyle=":")

        res = (data_points - total_pdf_values_for_bins) / np.sqrt(total_pdf_values_for_bins)

        ax2.errorbar(bin_centers_filtered, res, xerr=bin_widths_filtered, yerr=np.ones(len(res)),
                     color="k", fmt="o", markersize=5, linewidth=2)

        ax2.set_ylim(-10., 10)
        ax2.set_xlim(bin_centers_filtered[0] - bin_widths_filtered[0], bin_centers_filtered[-1] + bin_widths_filtered[-1])

        plt.tight_layout()
        plt.savefig(suffix + "combined_fit_and_pull_plot.pdf")
        plt.show()

    else:
        b2kstll.plot.plot_distributions_from_fitresult(result, suffix=suffix, with_pull=False, with_blind=True)
        plt.tight_layout()
        plt.savefig(suffix + "fit_plot.pdf")
        plt.show()
plot_style(result_total, data_points, total_pdf_values_for_bins, bin_centers_filtered, bin_widths_filtered, suffix="total_fit_", with_pulls=True)