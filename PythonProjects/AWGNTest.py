import numpy as np
import matplotlib.pyplot as plt

N = 1024
x =np.random.randn(N)
plt.plot(x, '.-')
plt.show()

X = np.fft.fftshift(np.fft.fft(x))
X = X[N//2:]
plt.plot(np.real(X),'.-')
plt.show()