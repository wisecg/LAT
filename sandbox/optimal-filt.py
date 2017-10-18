#!/usr/bin/env python
import numpy as np
from scipy import fft, ifft
import matplotlib.pyplot as plt
import pywt
import waveLibs as wl

def main():
    """ Write an optimual filter like CUORE: https://arxiv.org/pdf/1012.3266.pdf
    Or like Sec. 4 of Marco Carretoni's thesis:
         https://cuore.lngs.infn.it/sites/default/files/tesiPhD_MCarrettoni.pdf
    """

    # Load data and template
    npzfile = np.load("./data/optimumInputs.npz")
    rl, tl = npzfile['arr_0'], npzfile['arr_1']
    data, dataTS, dataE, dataST = rl[0], rl[1], rl[2], rl[3]
    temp, tempTS, tempE, tempST = tl[0], tl[1], tl[2], tl[3]

    # Scale and time-align template
    temp = temp * (dataE / tempE)
    tempTS -= tempST - dataST

    # Denoise WF and find maximum
    wp = pywt.WaveletPacket(data=data, wavelet='haar', mode='symmetric',maxlevel=3)
    new_wp = pywt.WaveletPacket(data=None, wavelet='haar', mode='symmetric')
    new_wp['aaa'] = wp['aaa'].data
    newWave = new_wp.reconstruct(update=False)
    idxRise = np.where((dataTS >= 8000) & (dataTS <= 11000))
    maxIdx = np.argmax(newWave[idxRise])
    maxTS = dataTS[idxRise][maxIdx]
    print "maxTS:",maxTS

    # Find transfer function H(w)

    xfftSignal, powSignal = fftWF(temp)

    expo = np.exp(1j * xfftSignal * maxTS*2)  # why does maxTS*2 seem to align it?

    xfftNoise, powNoise = fftWF(data)

    # try a padded version only w/ the noise part of the baseline
    # xfftNoise, powNoise = fftWF(data, 0, maxIdx-10)
    # nPad = (len(powSignal) - len(powNoise))/2
    # powNoise = np.lib.pad(powNoise, (nPad), 'mean')

    # H = expo * np.conjugate(powSignal) / np.power(powNoise,2)
    H = np.conjugate(powSignal) / np.power(powNoise,2)
    # H = np.conjugate(powSignal) / powNoise

    # Inverse transforms
    invTrH = ifftWF(H,len(data))
    invTrSignal = ifftWF(powSignal,len(data))
    invTrExpo = ifftWF(expo,len(data))
    invTrNoise = ifftWF(powNoise,len(data))

    # plots
    fig = plt.figure(figsize=(8,8), facecolor='w')
    a1 = plt.subplot(311)
    a1.set_title("Energy %.1f keV  Start Time %.0f ns  Max Time %.0f ns" % (dataE, dataST, maxTS))
    a1.set_ylabel("ADC [A.U.]",y=0.95, ha='right')
    a1.set_xlabel("Time (ns)",x=0.95, ha='right')
    a1.plot(dataTS,data,color='blue',alpha=0.8,label='Data WF')
    # a1.plot(dataTS,newWave,color='red',label='Lvl3 Haar Denoised')
    a1.plot(tempTS,temp,color='orange',label='Expected (Template) WF')
    a1.axvline(x=maxTS,color='cyan')
    a1.axvline(x=dataST,color='green')
    a1.legend(loc=4)

    a2 = plt.subplot(312)
    # a2.loglog(xfftNoise,powNoise,label='Data FFT')
    # a2.loglog(xfftSignal,powSignal,alpha=0.8,label='Expected FFT')
    a2.semilogy(xfftSignal,expo,color='blue',label='Expo')
    a2.loglog(xfftSignal,powSignal,color='orange',label='Signal')
    a2.semilogy(xfftSignal,powNoise,color='green',label='Noise/Data')

    a2.legend(loc=1)

    a3 = plt.subplot(313)
    # a3.set_ylim(-0.03,0.1)
    a3.plot(dataTS,invTrSignal,color='orange',label='Signal')
    a3.plot(dataTS,invTrExpo,color='blue',label='Expo')
    a3.plot(dataTS,invTrNoise,color='green',label='Noise')
    a3.plot(dataTS,invTrH,color='red',alpha=0.7,label='H')
    a3.legend(loc=1)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.35)
    # plt.show(block=False)
    # plt.show()
    plt.savefig("./plots/optimal-filt.pdf")


def fftWF(data, tLo=0, tHi=0):
    if tHi==0: tHi = len(data)-1
    yf = fft(data[tLo:tHi])
    T = 1e-8 # sampling frequency
    N = len(yf)
    xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
    yret = 2.0/N * np.abs(yf[:N/2]) # FT spectrum w/ absval
    # yret = 2.0/N * yf[:N/2]
    # ypow = np.multiply(yret, yret)
    return xf, yret

def ifftWF(powSpec, n):
    # invTrData = ifft(powSpec, 2500)*2*np.pi*2*np.pi*len(powSpec)
    invTrData = ifft(powSpec, n)
    return invTrData

def generateBaseline(waveFTPow):
	baselineWavePow = np.zeros(len(waveFTPow), dtype=complex)
	for idx in range(5, len(waveFTPow)):
		wp = waveFTPow[idx]
		baselineWavePow[idx] = np.complex(abs(wp*np.random.random_sample()), abs(wp*np.random.random_sample()) )/np.sqrt(2.)

	baselineWave = ifft(baselineWavePow, 2500)*2*np.pi*2*np.pi*len(baselineWavePow)
	# baselineWaveT = baselineWave[50:2016+50]
	# baselineWaveTS = np.linspace(1, 2016, 2016)
	baselineWaveTS = np.linspace(1, len(baselineWave), len(baselineWave))

	return baselineWaveTS, baselineWave



if __name__ == "__main__":
    main()
