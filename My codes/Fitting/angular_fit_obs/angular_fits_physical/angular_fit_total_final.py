import numpy as np
import zfit
import tensorflow as tf
import pandas as pd
import b2kstll
from b2kstll.models import legendre
from b2kstll.models.legendre import LegendrePolynomialPDF_8
from math import pi

# Load your data from pickle file (adjust path as necessary)
file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeFlat_all_years_MC.pkl"
data = pd.read_pickle(file_path)
data = data[data['XGBOutput'] > 0.9]


c1L = zfit.Parameter('c1L', 0.1, -1., 1.)
c2L = zfit.Parameter('c2L', 0.1, -1., 1.)
c3L = zfit.Parameter('c3L', 0.1, -1., 1.)
c4L = zfit.Parameter('c4L', 0.1, -1., 1.)
c5L = zfit.Parameter('c5L', 0.1, -1., 1.)
c6L = zfit.Parameter('c6L', 0.1, -1., 1.)
c7L = zfit.Parameter('c7L', 0.1, -1., 1.)
c8L = zfit.Parameter('c8L', 0.1, -1., 1.)

ctl_lhcb_data = data['ctl_lhcb']
obs_ctl = zfit.Space('ctl', limits=(-1, 1))
pdf_ctl = legendre.LegendrePolynomialPDF_8(obs=obs_ctl, c1L=c1L, c2L=c2L, c3L=c3L, c4L=c4L, c5L=c5L, c6L=c6L, c7L=c7L, c8L=c8L)

data_ctl_lhcb = zfit.data.Data.from_numpy(obs=obs_ctl, array=ctl_lhcb_data.values)
nll_ctl = zfit.loss.UnbinnedNLL(model=pdf_ctl, data=data_ctl_lhcb)

# Create minimizer
minimizer_ctl = zfit.minimize.Minuit()

# Perform minimization
result_ctl = minimizer_ctl.minimize(nll_ctl)

# b2kstll.plot.plot_distributions(result_ctl, suffix="ctl")



ctk_data = data['ctk_lhcb'].values

c1k = zfit.Parameter('c1k', 0.1, -1., 1.)
c2k = zfit.Parameter('c2k', 0.1, -1., 1.)
c3k = zfit.Parameter('c3k', 0.1, -1., 1.)
c4k = zfit.Parameter('c4k', 0.1, -1., 1.)
c5k = zfit.Parameter('c5k', 0.1, -1., 1.)
c6k = zfit.Parameter('c6k', 0.1, -1., 1.)
c7k = zfit.Parameter('c7k', 0.1, -1., 1.)
c8k = zfit.Parameter('c8k', 0.1, -1., 1.)
ctk = zfit.Space('ctk', limits=(-1, 1))

pdf_ctk = legendre.LegendrePolynomialPDF_8(obs=ctk, c1L=c1k, c2L=c2k, c3L=c3k, c4L=c4k, c5L=c5k, c6L=c6k, c7L=c7k, c8L=c8k)

data_k = zfit.data.Data.from_numpy(obs=ctk, array=ctk_data)

# Create negative log-likelihood
nll_ctk = zfit.loss.UnbinnedNLL(model=pdf_ctk, data=data_k)

# Create minimizer
minimizer_ctk = zfit.minimize.Minuit()

# Perform minimization
result_ctk = minimizer_ctk.minimize(nll_ctk)

# Note: If you want to plot using b2kstll, ensure you have imported it correctly and adjust as necessary
# b2kstll.plot.plot_distributions(result_ctk, suffix="ctk_lhcb_legendre") #labels=['Legendre Polynomial Order 6'])

phi_data = data['phi_lhcb']


# Define parameters
c1phi = zfit.Parameter('c1phi', 0.0, -1.0, 1.0)
c2phi = zfit.Parameter('c2phi', 0.1, -1.0, 1.0)
c3phi = zfit.Parameter('c3phi', 0.0, -1.0, 1.0)
c4phi = zfit.Parameter('c4phi', 0.1, -1.0, 1.0)
c5phi = zfit.Parameter('c5phi', 0.0, -1.0, 1.0)
c6phi = zfit.Parameter('c6phi', 0.1, -1.0, 1.0)
c7phi = zfit.Parameter('c7phi', 0.0, -1.0, 1.0)
c8phi = zfit.Parameter('c8phi', 0.1, -1.0, 1.0)


# Define observable space
obs_phi = zfit.Space('phi', limits=(-np.pi, np.pi))

# Define PDF
pdf_phi = legendre.LegendrePolynomialPDF_8(obs=obs_phi, c1L=c1phi, c2L=c2phi, c3L=c3phi, c4L=c4phi, c5L=c5phi, c6L=c6phi, c7L=c7phi, c8L=c8phi)

# Convert data to zfit format
data_phi = zfit.data.Data.from_numpy(obs=obs_phi, array=phi_data.values)

# Evaluate the PDF with initial parameters
pdf_values = pdf_phi.pdf(data_phi)

# Create negative log-likelihood
nll_phi = zfit.loss.UnbinnedNLL(model=pdf_phi, data=data_phi)

# Evaluate the NLL with initial parameters
initial_loss = nll_phi.value()

# Create minimizer
minimizer_phi = zfit.minimize.Minuit()

# Perform minimization (without callback)
result_phi = minimizer_phi.minimize(nll_phi)








# Define the P-wave PDF
class PWavePDF(zfit.pdf.ZPDF):
    _N_OBS = 3
    _PARAMS = ['FL', 'S3', 'S4', 'S5', 'AFB', 'S7', 'S8', 'S9']

    @zfit.supports()
    def _unnormalized_pdf(self, x,params):


        costheta_l = x[0]   
        costheta_k = x[1]
        phi = x[2]  

        sintheta_k = tf.sqrt(1.0 - tf.square(costheta_k))
        sintheta_l = tf.sqrt(1.0 - tf.square(costheta_l))

        # Add debug checks for intermediate calculations
        tf.debugging.check_numerics(sintheta_k, "Check sintheta_k for NaNs")
        tf.debugging.check_numerics(sintheta_l, "Check sintheta_l for NaNs")

        sintheta_2k = 1.0 - tf.square(costheta_k)
        sintheta_2l = 1.0 - tf.square(costheta_l)

        sin2theta_k = 2.0 * sintheta_k * costheta_k
        cos2theta_l = 2.0 * tf.square(costheta_l) - 1.0
        sin2theta_l = 2.0 * sintheta_l * costheta_l

        # Additional debug checks
        tf.debugging.check_numerics(sintheta_2k, "Check sintheta_2k for NaNs")
        tf.debugging.check_numerics(sintheta_2l, "Check sintheta_2l for NaNs")
        tf.debugging.check_numerics(sin2theta_k, "Check sin2theta_k for NaNs")
        tf.debugging.check_numerics(cos2theta_l, "Check cos2theta_l for NaNs")
        tf.debugging.check_numerics(sin2theta_l, "Check sin2theta_l for NaNs")

        FL = params['FL']
        S3 = params['S3']
        S4 = params['S4']
        S5 = params['S5']
        AFB = params['AFB']
        S7 = params['S7']
        S8 = params['S8']
        S9 = params['S9']

        pdf = ((3.0 / 4.0) * (1.0 - FL) * sintheta_2k +
            FL * tf.square(costheta_k) +
            (1.0 / 4.0) * (1.0 - FL) * sintheta_2k * cos2theta_l -
            FL * tf.square(costheta_k) * cos2theta_l +
            S3 * sintheta_2k * sintheta_2l * tf.cos(2.0 * phi) +
            S4 * sin2theta_k * sin2theta_l * tf.cos(phi) +
            S5 * sin2theta_k * sintheta_l * tf.cos(phi) +
            (4.0 / 3.0) * AFB * sintheta_2k * costheta_l +
            S7 * sin2theta_k * sintheta_l * tf.sin(phi) +
            S8 * sin2theta_k * sin2theta_l * tf.sin(phi) +
            S9 * sintheta_2k * sintheta_2l * tf.sin(2.0 * phi))

        # Final debug check for pdf
        pdf = tf.debugging.check_numerics(pdf, "Check pdf for NaNs")
        return (9.0 / (32.0 * pi)) * pdf

# Load the data from the pickle file
file_path = "/sps/lhcb/yahya/Internship/pkl_creation/pkl/B02Kst0eeSignal_all_years_MC.pkl"
data = pd.read_pickle(file_path)
df = data[data['XGBOutput'] > 0.9]

# Prepare the data for fitting
data_array_ctl = df['ctl_lhcb'].to_numpy()
data_array_ctk = df['ctk_lhcb'].to_numpy()
data_array_phi = df['phi_lhcb'].to_numpy()

assert not np.isnan(data_array_ctl).any(), "NaN values found in data_array_ctl"
assert not np.isnan(data_array_ctk).any(), "NaN values found in data_array_ctk"
assert not np.isnan(data_array_phi).any(), "NaN values found in data_array_phi"

# Correct definition of observables with multidimensional limits
obs_ctl2 = zfit.Space('ctl_lhcb', limits=(-1, 1))
obs_ctk2 = zfit.Space('ctk_lhcb', limits=(-1, 1))
obs_phi2 = zfit.Space('phi_lhcb', limits=(-pi, pi))

# Combine the spaces into a single multidimensional space
obs2 = obs_ctl2 * obs_ctk2 * obs_phi2

# Create parameters
FL = zfit.Parameter("FL", 0.5, 0, 1)
S3 = zfit.Parameter("S3", 0.0, -1, 1)
S4 = zfit.Parameter("S4", 0.0, -1, 1)
S5 = zfit.Parameter("S5", 0.0, -1, 1)
AFB = zfit.Parameter("AFB", 0.0, -1, 1)
S7 = zfit.Parameter("S7", 0.0, -1, 1)
S8 = zfit.Parameter("S8", 0.0, -1, 1)
S9 = zfit.Parameter("S9", 0.0, -1, 1)

# Create the PDF
pdf_tot = PWavePDF(obs=obs2,FL=FL, S3=S3, S4=S4, S5=S5, AFB=AFB, S7=S7, S8=S8, S9=S9)

# # Acceptance function calculation
ctl_pdf_values = pdf_ctl.pdf(data_array_ctl)
ctk_pdf_values = pdf_ctk.pdf(data_array_ctk)
phi_pdf_values = pdf_phi.pdf(data_array_phi)

# Normalize the acceptance PDFs
norm_ctl = tf.reduce_sum(ctl_pdf_values) / ctl_pdf_values.shape[0]
norm_ctk = tf.reduce_sum(ctk_pdf_values) / ctk_pdf_values.shape[0]
norm_phi = tf.reduce_sum(phi_pdf_values) / phi_pdf_values.shape[0]

# Calculate the acceptance weights
acceptance_tot = 1 / ((ctl_pdf_values / norm_ctl) * (ctk_pdf_values / norm_ctk) * (phi_pdf_values / norm_phi))

print(acceptance_tot)
# acceptance_tot = 1 / (ctl_pdf_values * ctk_pdf_values * phi_pdf_values)

# Create a dataset
zdata = zfit.Data.from_numpy(obs=obs2, array=np.array([data_array_ctl, data_array_ctk, data_array_phi]).T, weights=acceptance_tot*df['w_Track_Trig_Kin'])

# Define a loss function
loss2 = zfit.loss.UnbinnedNLL(model=pdf_tot, data=zdata)

# Choose a minimizer


minimizer2 = zfit.minimize.Minuit()

# Minimize the loss
result_tot = minimizer2.minimize(loss2)
print(result_tot)

param_errors = result_tot.hesse()

# Print the results with errors
params = [FL, S3, S4, S5, AFB, S7, S8, S9]
for param in params:
    value = param.numpy()
    error = param_errors[param]['error']
    print(f"{param.name}: {value} ± {error}")


# result_tot.loss.draw_projection(pdf_tot, zdata)

# # Plot the results
import b2kstll
b2kstll.plot.plot_distributions(model=pdf_tot, data=zdata, suffix="total_fit",with_pulls=False)  