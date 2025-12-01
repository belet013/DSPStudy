# #FFT Example
# import numpy as np
# import matplotlib.pyplot as plt

# Fs = 1
# N = 100

# t = np.arange(N)
# s = np.sin(0.15 * 2 * np.pi * t)
# s = s * np.hamming(100)
# S = np.fft.fftshift(np.fft.fft(s))
# S_mag = np.abs(S)
# S_phase = np.angle(S)
# f = np.arange(-Fs/2, Fs/2, Fs/N)
# plt.figure(0)
# plt.plot(f, S_mag, '.-')
# plt.figure(1)
# plt.plot(f, S_phase, '.-')
# plt.show()

##Spectogram Example
# import numpy as np
# import matplotlib.pyplot as plt
# sample_rate = 1e6

# # Generate tone plus noise
# t = np.arange(1024*1000)/sample_rate # time vector
# f = 50e3 # freq of tone
# x = np.sin(2*np.pi*f*t) + 0.2*np.random.randn(len(t))

# # simulate the signal above, or use your own signal

# fft_size = 1024
# num_rows = len(x) // fft_size # // is an integer division which rounds down
# spectrogram = np.zeros((num_rows, fft_size))
# for i in range(num_rows):
#     spectrogram[i,:] = 10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(x[i*fft_size:(i+1)*fft_size])))**2)

# plt.imshow(spectrogram, aspect='auto', extent = [sample_rate/-2/1e6, sample_rate/2/1e6, len(x)/sample_rate, 0])
# plt.xlabel("Frequency [MHz]")
# plt.ylabel("Time [s]")
# plt.show()


# FFT w/o np.fft() Example
import numpy as np
import matplotlib.pyplot as plt

def fft(x):
    N = len(x)
    if N == 1:
        return x
    twiddle_factors = np.exp(-2j * np.pi * np.arange(N//2) / N)
    x_even = fft(x[::2])
    x_odd = fft(x[1::2])
    return np.concatenate([x_even + twiddle_factors * x_odd, x_even - twiddle_factors * x_odd])

# Simulating tone and noise
Fs = 1e6 # 1 MHz sample rate
Fo = 0.2e6 # 200 kHz offset from carrier
N = 1024 # number of samples
t = np.arange(N)/Fs # time vector
s = np.exp(2j * np.pi * Fo * t) # complex tone
n = (np.random.randn(N) + 1j * np.random.randn(N)) / np.sqrt(2) # Unity complex white noise 
r = s + n # received signal 0 dB SNR
r = r * np.hamming(N) # window our signal

# Perform the fft, fftshift, and convert to dB
X = fft(r) # equivalent to np.fft.fft(r)
X_shifted = np.roll(X, N//2) # equivalent to np.fft.fftshift()
X_mag = 10*np.log10(np.abs(X_shifted)**2) # convert to dB

# Plot our results
f = np.linspace(Fs/-2, Fs/2, N) / Fs
plt.plot(f, X_mag)
plt.plot(f[np.argmax(X_mag)],np.max(X_mag), 'rx')
plt.grid()
plt.xlabel('Freqency [MHz]')
plt.ylabel('Magnitude [dB]')
plt.show()