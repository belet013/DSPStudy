import numpy as np
import matplotlib.pyplot as plt

Fs = 300 # Sample rate
Ts = 1/Fs # Sample period
N = 2048 # Number of samples to simulate

t = Ts*np.arange(N)
x = np.exp(1j*2*np.pi*50*t) # 50 Hz sine

# Create complex noise w/ unity power
n = (np.random.randn(N) + 1j*np.random.randn(N)) / np.sqrt(2) 
noise_power = 2
r = x + n * np.sqrt(noise_power)

windowed_x = r * np.hamming(len(x)) # Hamming window applied

PSD = np.abs(np.fft.fft(windowed_x))**2 / (N*Fs)
PSD_log = 10.0*np.log10(PSD)
PSD_shifted = np.fft.fftshift(PSD_log)

f = np.arange(Fs/-2.0, Fs/2.0, Fs/N) # Start, Stop, Step

plt.plot(f, PSD_shifted)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
plt.grid(True)
plt.show()