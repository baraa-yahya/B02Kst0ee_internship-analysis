import flavio

# Initialize Wilson coefficients (Standard Model defaults)
wc = flavio.WilsonCoefficients()

# Specify the q^2 range (for example, q2min = 1 GeV^2, q2max = 6 GeV^2)
q2min = 14.33204  # GeV^2
q2max = 20  # GeV^2

# Compute the observables with the q^2 cut
FL = flavio.np_prediction('<FL>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
AFB = flavio.np_prediction('<AFB>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S3 = flavio.np_prediction('<S3>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S4 = flavio.np_prediction('<S4>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S5 = flavio.np_prediction('<S5>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S7 = flavio.np_prediction('<S7>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S8 = flavio.np_prediction('<S8>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
S9 = flavio.np_prediction('<S9>(B0->K*ee)', wc, q2min=q2min, q2max=q2max)
# Print the results
print(f"FL: {FL}")
print(f"S3: {S3}")
print(f"S4: {S4}")
print(f"S5: {S5}")
print(f"AFB: {AFB}")
print(f"S7: {S7}")
print(f"S8: {S8}")
print(f"S9: {S9}")
