#!/usr/local/bin/python
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
import h5py

""" 'Find an Inspiral' Tutorial: https://losc.ligo.org/tutorial06/ """

# Read in data and template
fs = 4096
dataFile = h5py.File('./data/ligo_data_w_signal.hdf5', 'r')
data = dataFile['strain/Strain'][...]
dataFile.close()
time = np.arange(0, 16, 1./fs)

templateFile = h5py.File('./data/ligo_template.hdf5', 'r')
template = templateFile['strain/Strain'][...]
temp_time = np.arange(0, template.size / (1.0*fs), 1./fs)
templateFile.close()

def main():

    fig = plt.figure(figsize=(15,10), facecolor='w')
    p1 = plt.subplot(231)
    p2 = plt.subplot(232)
    p3 = plt.subplot(233)
    p4 = plt.subplot(234)
    p5 = plt.subplot(235)
    p6 = plt.subplot(236)

    # Plot data and template
    p1.plot(time,data)
    p1.set_xlabel('Time (s)')
    p1.set_ylabel('Strain')
    p2.plot(temp_time,template)
    p2.set_xlabel('Time (s)')
    p2.set_ylabel('Strain')
    p2.set_title('Template')

    # Plot ASD of data
    power_data, freq_psd = plt.psd(data[12*fs:], Fs=fs, NFFT=fs, visible=False)
    p3.loglog(freq_psd, np.sqrt(power_data), color='blue')
    p3.set_xlim([20, 2048])
    p3.set_xlabel('Frequency (Hz)')
    p3.set_ylabel('ASD')

    # Plot ASD of template
    power, freq = plt.psd(template, Fs=fs, NFFT=fs, visible=False)
    p3.loglog(freq, np.sqrt(power), color='red')
    p3.grid('on')

    # Apply a bandpass filter to the data and plot it
    (B,A) = sig.butter(4, [80/(fs/2.0), 250/(fs/2.0)], btype='pass')
    data_pass= sig.lfilter(B, A, data)
    p4.plot(time, data_pass)
    p4.set_title('Band passed data')
    p4.set_xlabel('Time (s)')

    # Plot time-domain cross-correlation
    correlated_raw = np.correlate(data, template, 'valid')
    correlated_passed = np.correlate(data_pass, template, 'valid')
    p5.plot(np.arange(0, (correlated_raw.size*1.)/fs, 1.0/fs),correlated_raw)
    p5.set_title('Time domain cross-correlation')
    p5.set_xlabel('Offest between data and template (s)')

    p5.plot(np.arange(0, (correlated_passed.size*1.)/fs, 1.0/fs), correlated_passed,color='red')
    p5.set_xlabel('Offset between data and template (s)')
    # p5.set_title('Band passed time domain cross-correlation')

    # Optimal Matched Filter (freq. domain)

    # Take the FFT of the data
    data_fft = np.fft.fft(data)

    # Pad template and take FFT
    zero_pad = np.zeros(data.size - template.size)
    template_padded = np.append(template, zero_pad)
    template_fft = np.fft.fft(template_padded)

    # Match FFT frequency bins to PSD frequency bins
    power_data, freq_psd = plt.psd(data[12*fs:], Fs=fs, NFFT=fs, visible=False)
    datafreq = np.fft.fftfreq(data.size)*fs
    power_vec = np.interp(datafreq, freq_psd, power_data)

    # Apply the optimal matched filter
    optimal = data_fft * template_fft.conjugate() / power_vec
    optimal_time = 2*np.fft.ifft(optimal)

    # Normalize the matched filter output
    df = np.abs(datafreq[1] - datafreq[0])
    sigmasq = 2*(template_fft * template_fft.conjugate() / power_vec).sum() * df
    sigma = np.sqrt(np.abs(sigmasq))
    SNR = abs(optimal_time) / (sigma)

    # Plot the result
    p6.cla()
    p6.plot(time, SNR)
    p6.set_title('Optimal Matched Filter')
    p6.set_xlabel('Offset time (s)')
    p6.set_ylabel('SNR')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
