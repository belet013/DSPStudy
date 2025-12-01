import numpy as np
import matplotlib.pyplot as plt
from sigmf import SigMFFile, sigmffile

# # Read in file. We need to tell it what format it is
# samples = np.fromfile('qpsk_in_noise.iq', np.complex64)
# print(samples)

# # Plot the constellation to make sure it looks right
# plt.plot(np.real(samples), np.imag(samples), '.')
# plt.grid(True)
# plt.show()

# Load a dataset
filename = 'qpsk_in_noise'
signal = sigmffile.fromfile(filename)
samples = signal.read_samples().view(np.complex64).flatten()
print(samples[0:10]) # lets look at the first 10 samples

# Get some metadata and all annotations
sample_rate = signal.get_global_field(SigMFFile.SAMPLE_RATE_KEY)
sample_count = signal.sample_count
signal_duration = sample_count / sample_rate