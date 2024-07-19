import zfit
from zfit import z
import tensorflow as tf
import numpy as np
import pandas as pd
from math import pi


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
df = data[data['XGBOutput'] > 0.995]

# Prepare the data for fitting
data_array_ctl = df['ctl_lhcb'].to_numpy()
data_array_ctk = df['ctk_lhcb'].to_numpy()
data_array_phi = df['phi_lhcb'].to_numpy()

assert not np.isnan(data_array_ctl).any(), "NaN values found in data_array_ctl"
assert not np.isnan(data_array_ctk).any(), "NaN values found in data_array_ctk"
assert not np.isnan(data_array_phi).any(), "NaN values found in data_array_phi"

# Correct definition of observables with multidimensional limits
obs_ctl = zfit.Space('ctl_lhcb', limits=(-1, 1))
obs_ctk = zfit.Space('ctk_lhcb', limits=(-1, 1))
obs_phi = zfit.Space('phi_lhcb', limits=(-pi, pi))

# Combine the spaces into a single multidimensional space
obs = obs_ctl * obs_ctk * obs_phi

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
pdf = PWavePDF(obs=obs,FL=FL, S3=S3, S4=S4, S5=S5, AFB=AFB, S7=S7, S8=S8, S9=S9)

# Create a dataset
zdata = zfit.Data.from_numpy(obs=obs, array=np.array([data_array_ctl, data_array_ctk, data_array_phi]).T)

# Define a loss function
loss = zfit.loss.UnbinnedNLL(model=pdf, data=zdata)

# Choose a minimizer


minimizer = zfit.minimize.Minuit()

# Minimize the loss
result = minimizer.minimize(loss)
print(result)
param_errors = result.hesse()

# Print the results with errors
params = [FL, S3, S4, S5, AFB, S7, S8, S9]
for param in params:
    value = param.numpy()
    error = param_errors[param]['error']
    print(f"{param.name}: {value} ± {error}")




# Plot the results
import b2kstll
b2kstll.plot.plot_distributions_using_partial_integrals(model=pdf, data = zdata, suffix ="test_lhcb", with_pulls=False)
# b2kstll.plot.plot_distributions(result, suffix="pwave")