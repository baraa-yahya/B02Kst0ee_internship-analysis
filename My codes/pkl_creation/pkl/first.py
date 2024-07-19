import pickle
import matplotlib.pyplot as plt

# 1. Load the pickle file
with open('B02Kst0eeSignal_201_MC.pkl', 'rb') as file:
    data_dict = pickle.load(file)

# Plot histogram of the 'B_DTF_PV_B_M_0' column with connected lines
plt.hist(data_dict['B_DTF_PV_B_M_0'], bins=200, histtype='step', color='skyblue', edgecolor='black')
plt.xlabel('B_DTF_PV_B_M_0')
plt.ylabel('Frequency')
plt.title('Step Histogram of B_DTF_PV_B_M_0')
plt.show()
plt.savefig("first_plot_7.pdf")

# THE PICKLE FILE WORK NOW MAKE IT GENERAL AND START FITTING PLEASE!!!!!!!!!!!!!!