import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import sigmf
from sigmf import SigMFFile

num_symbols = 10000

# x_symbols array will contain complex numbers representing a QPSK
# symbols. Each symbol will be a complex number with a mangitude of 
# 1 and a phase angle corresponding to one of the four QPSK constellation
# points (45,135,225,315)

x_int = np.random.randint(0,4,num_symbols) # random int 0 to 3
x_degrees = x_int*360/4.0 + 45             # 45,135,225,315 degrees
x_radians = x_degrees*np.pi/180            # sin and cos takes in radians
x_symbols = np.cos(x_radians) + 1j*np.sin(x_radians) # produces our QPSK complex nunbers
n = (np.random.randn(num_symbols) + 1j*np.random.randn(num_symbols))/np.sqrt(2) # AWGN w/ unity power
r = x_symbols + n * np.sqrt(0.01) # noise power of 0.01

print(r)
plt.plot(np.real(r),np.imag(r), '.')
plt.grid(True)
plt.show()

# Saving to an IQ file
print(type(r[0])) # Check data type, should show Complex128
r = r.astype(np.complex64) # Convert to Complex64
print(type(r[0])) # Verify it's Complex64
#r.tofile('qpsk_in_noise.iq') # Save to file

# replaced r.tofile('qpsk_in_noise.iq')
# r.tofile('qpsk_in_noise.iq')
r.tofile('qpsk_in_noise.sigmf-data') # replace line above with this one

# create the metadata
meta = SigMFFile(
    data_file='qpsk_in_noise.sigmf-data', # extension is optional
    global_info = {
        SigMFFile.DATATYPE_KEY: 'cf32_le',
        SigMFFile.SAMPLE_RATE_KEY: 8000000,
        SigMFFile.AUTHOR_KEY: 'Your name and/or email',
        SigMFFile.DESCRIPTION_KEY: 'Simulation of qpsk with noise',
        SigMFFile.VERSION_KEY: sigmf.__version__,
    }
)

# create a capture key at time index 0
meta.add_capture(0, metadata={
    SigMFFile.FREQUENCY_KEY: 915000000,
    SigMFFile.DATETIME_KEY: dt.datetime.now(dt.timezone.utc).isoformat(),
}) 


# check for mistakes and write to disk
meta.validate()
meta.tofile('qpsk_in_noise.sigmf-meta') # extension is optional