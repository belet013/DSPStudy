import numpy as np
from scipy.signal import firwin2, convolve, fftconvolve, lfilter
import matplotlib.pyplot as plt

# Create a test signal, we'll use a Gaussian noise
sample_rate = 1e6 # Hz
N = 10000 # Samples to simulate
x = np.random.randn(N) + 1j * np.random.randn(N)

# Save PSD of the input signal
PSD_input = 10*np.log10(np.fft.fftshift(np.abs(np.fft.fft(x))**2)/len(x))


# Create an FIR filter, same on as 2nd example above
freqs = [0,100e3, 110e3, 190e3, 200e3, 300e3, 310e3, 500e3]
gains = [1,1,0,0,0.5,0.5,0,0]
h2 = firwin2(101, freqs, gains, fs=sample_rate)

# Apply filter using 4 methods
#x_numpy = np.convolve(h2,x)
#x_scipy = convolve(h2,x)
x_fft_convolve = fftconvolve(x,h2, 'same')
#x_lfilter = lfilter(h2,1,x)

# Look at PSD of the output signal
PSD_output  = 10*np.log10(np.fft.fftshift(np.abs(np.fft.fft(x_fft_convolve))**2)/len(x_fft_convolve))
f = np.linspace(-sample_rate/2/1e6, sample_rate/2/1e6,len(PSD_output))
plt.plot(f,PSD_input, alpha=0.8)
plt.plot(f,PSD_output, alpha=0.8)
plt.xlabel('Frequency [Mhz]')
plt.ylabel('PSD [dB]')
plt.axis([sample_rate/-2/1e6,sample_rate/2/1e6, -40, 20])
plt.legend(['Input', 'Output'], loc = 1)
plt.grid()
plt.show()
# print(x_numpy[0:2])
# print(x_scipy[0:2])
# print(x_fft_convolve[0:2])
# print(x_lfilter[0:2])

# import numpy as np
# from scipy import signal
# import matplotlib.pyplot as plt

# num_taps = 51       # Use an odd nmber of taps
# cut_off = 3000      # Hz
# sample_rate = 32000 # Hz

# # Create a low-pass filter
# h = signal.firwin(num_taps, cut_off, fs=sample_rate)

# # Shift the filter in frequency by multiplying by exp(j*2*pi*f0*t)
# f0 = 10e3               # Shift amount
# Ts = 1.0/sample_rate    # Sample period
# t = np.arange(0.0, Ts*len(h), Ts)   # Time vector. args are (start, stop, step)
# exponential = np.exp(2j*np.pi*f0*t) # Essentially a complex sine wave

# h_band_pass = h * exponential # Applying shift

# # Plot the impulse response
# plt.figure('impulse')
# plt.plot(np.real(h_band_pass), '.-')
# plt.plot(np.imag(h_band_pass), '.-')
# plt.legend(['real', 'imag'], loc=1)

# # Plot the frequency response
# H = np.abs(np.fft.fft(h_band_pass,1024)) # 1024-point fft & magnitude
# H = np.fft.fftshift(H) # Center to 0 Hz
# w = np.linspace(-sample_rate/2,sample_rate/2,len(H)) # X-axis
# plt.figure('freq')
# plt.plot(w,H,'.-')
# plt.xlabel('Frequency [Hz]')
# plt.show()




# # # Plot the impulse response
# # plt.plot(h, '.-')
# # plt.show()