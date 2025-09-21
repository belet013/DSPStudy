import numpy as np, matplotlib.pyplot as plt

fs, N = 1_000_000, 4096
rng = np.random.default_rng(0)
x = rng.normal(0, 1, N)

X = np.fft.fft(x)
f = np.fft.fftfreq(N, 1/fs)

# center around 0 Hz and normalize
Xc = np.fft.fftshift(X) / N
fc = np.fft.fftshift(f)

plt.plot(fc, 20*np.log10(np.abs(Xc)+1e-12))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude (dB)')
plt.title('Two-sided FFT (centered)')
plt.grid(True)
plt.tight_layout()
plt.show()
