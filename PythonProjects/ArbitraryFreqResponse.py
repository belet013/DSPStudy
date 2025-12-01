import numpy as np
import matplotlib.pyplot as plt

H = np.hstack((np.zeros(20), np.arange(10)/10, np.zeros(20)))
w = np.linspace(-0.5, 0.5, 50)
h = np.fft.ifftshift(np.fft.ifft(np.fft.ifftshift(H)))

window = np.hamming(len(h))
h = h * window

H_fft = np.fft.fftshift(np.abs(np.fft.fft(h,1024)))


plt.plot(H_fft)
plt.show()

plt.plot(np.real(h))
plt.plot(np.imag(h))
plt.legend(['real','imag'], loc=1)
plt.show()

plt.plot(w,H,'.-')
plt.show()