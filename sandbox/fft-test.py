import numpy as np
import matplotlib.pyplot as plt

steps = np.arange(0,5120,10)

# adc = np.ones(256) * 10
# adc = np.append(adc, np.zeros(256))

adc = np.ones(512) * 10

asd, xFreq = plt.psd(adc, Fs=1e8, NFFT=512, pad_to=512, visible=True)

plt.show()

