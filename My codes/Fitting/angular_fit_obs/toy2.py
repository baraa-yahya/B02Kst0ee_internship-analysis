import zfit
from zfit import z
import tensorflow as tf
import numpy as np
import progressbar
import flavio
from math import pi
import time

zfit.settings.set_verbosity(-1)

# Initialize Wilson coefficients (Standard Model defaults)
wc = flavio.WilsonCoefficients()

# Specify the q^2 range (for example, q2min = 1 GeV^2, q2max = 6 GeV^2)
q2min = 14.33204  # GeV^2
q2max = 20  # GeV^2

# Compute the observables with the q^2 cut
FL_SM = flavio.np_prediction('<FL>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
AFB_SM = flavio.np_prediction('<AFB>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S3_SM = flavio.np_prediction('<S3>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S4_SM = flavio.np_prediction('<S4>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S5_SM = flavio.np_prediction('<S5>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S7_SM = flavio.np_prediction('<S7>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S8_SM = flavio.np_prediction('<S8>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S9_SM = flavio.np_prediction('<S9>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)

# Print the results
print("Angular observables from SM:")
print(f"FL: {FL_SM}")
print(f"S3: {S3_SM}")
print(f"S4: {S4_SM}")
print(f"S5: {S5_SM}")
print(f"AFB: {AFB_SM}")
print(f"S7: {S7_SM}")
print(f"S8: {S8_SM}")
print(f"S9: {S9_SM}")

class PWavePDF(zfit.pdf.ZPDF):
    _N_OBS = 3
    _PARAMS = ['FL', 'S3', 'S4', 'S5', 'AFB', 'S7', 'S8', 'S9']

    @zfit.supports()
    def _unnormalized_pdf(self, x, params):
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

# Define the parameter space
FL = zfit.Parameter("FL", FL_SM , 0, 1)
S3 = zfit.Parameter("S3", S3_SM, -1, 1)
S4 = zfit.Parameter("S4", S4_SM, -1, 1)
S5 = zfit.Parameter("S5", S5_SM, -1, 1)
AFB = zfit.Parameter("AFB", AFB_SM, -1, 1)
S7 = zfit.Parameter("S7", S7_SM, -1, 1)
S8 = zfit.Parameter("S8", S8_SM, -1, 1)
S9 = zfit.Parameter("S9", S9_SM, -1, 1)

# Define the observable space
obs = zfit.Space("costheta_l", limits=(-1, 1)) * zfit.Space("costheta_k", limits=(-1, 1)) * zfit.Space("phi", limits=(-np.pi, np.pi))

# Create the model instance
model = PWavePDF(obs=obs, FL=FL, S3=S3, S4=S4, S5=S5, AFB=AFB, S7=S7, S8=S8, S9=S9)

# Generate toy data
n_events = 120
sampler = model.create_sampler(n=n_events)

# Define the loss function
loss = zfit.loss.UnbinnedNLL(model, sampler)

from zfit.minimize import DefaultToyStrategy

minimizer = zfit.minimize.Minuit(strategy=DefaultToyStrategy(), verbosity=0, tol=1e-5, use_minuit_grad=True)

fit_results = []
ntoys = 1

start_time = time.time()

params = loss.get_params()

with progressbar.ProgressBar(max_value=ntoys) as bar:
    while len(fit_results) < ntoys:
        # Generate toys
        sampler.resample()  # this is where the sampling happens

        # Randomize initial values
        for param in params:
            param.randomize()  # or smarter, use `set_value` for your own method

        # Minimize the NLL
        result = minimizer.minimize(loss)

        if result.converged:
            # Calculate uncertainties
            result.hesse()
            fit_results.append(result)
            bar.update(len(fit_results))

end_time = time.time()

print(f"Time taken: {end_time - start_time} seconds")
print(fit_results[:10])
