#!/usr/local/bin/python
import sys
from ROOT import TFile,gDirectory
import matplotlib.pyplot as plt
from matplotlib import gridspec
import waveLibs as wl
import numpy as np
import scipy.signal as sig
from scipy.ndimage.filters import gaussian_filter

""" This is an adaptation of the LIGO "Find an Inspiral" tutorial,
    applied to MJD data. https://losc.ligo.org/tutorial06/
    Course Material: http://www2.physics.umd.edu/~pshawhan/courses/CGWA/
    Especially Lecture 2: http://www2.physics.umd.edu/~pshawhan/courses/CGWA/ShawhanLecture2.pdf
"""
def main(argv):
    scanSpeed = 0.0001
    opt1, opt2 = "", ""
    intMode, printWF = False, False
    if (len(argv) >= 1): opt1 = argv[0]
    if (len(argv) >= 2): opt2 = argv[1]
    if "-i" in (opt1, opt2):
        intMode = True
        print "Interactive mode selected."
    if "-s" in (opt1, opt2):
        printWF = True
        print "Saving WF plots to current directory."

    # Set input file and cut
    inputFile = TFile("~/project/match-skim/waveletSkimDS5_run23920.root")
    # inputFile = TFile("~/project/v2-processwfs/waveletSkimDS5_90.root")
    waveTree = inputFile.Get("skimTree")
    theCut = inputFile.Get("cutUsedHere").GetTitle()
    theCut += " && waveS5/trapENFCal < 1200 && trapENFCal < 5"
    # theCut += " && trapENFCal > 1.5 && trapENFCal < 2"
    waveTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    waveTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Using cut:\n",theCut,"\n"
    print "Found",nList,"entries passing cuts."

    # Make a figure
    fig = plt.figure(figsize=(15,10), facecolor='w')
    # p1 = plt.subplot(231)
    # p2 = plt.subplot(232)
    # p3 = plt.subplot(233)
    # p4 = plt.subplot(234)
    # p5 = plt.subplot(235)
    # p6 = plt.subplot(236)
    p1 = plt.subplot2grid((2,3), (0,0), colspan=2)
    p3 = plt.subplot2grid((2,3), (0,2))
    p4 = plt.subplot2grid((2,3), (1,0))
    p5 = plt.subplot2grid((2,3), (1,1))
    p6 = plt.subplot2grid((2,3), (1,2))


    plt.show(block=False)

    # Make templates
    samp, r, z, tempE, tempStart, smooth = 5000, 0, 15, 10, 2500, 100  # huge template
    tempHuge, tempHugeTS = wl.MakeSiggenWaveform(samp,r,z,tempE,tempStart,smooth)

    samp, r, z, tempE, tempStart, smooth = 2016, 0, 15, 10, 1000, 100 # regular size template
    # samp, r, z, tempE, tempStart, smooth = 2016, 30, 30, 10, 1000, 100 # regular size temp, slower rise
    # samp, r, z, tempE, tempStart, smooth = 500, 0, 15, 10, 100, 100  # small template
    tempOrig, tempOrigTS = wl.MakeSiggenWaveform(samp,r,z,tempE,tempStart,smooth)


    # Load stuff from DS1 forced acq. runs
    npzfile = np.load("./data/fft_forcedAcqDS1.npz")
    avgPsd,xPsd,avgPwrSpec,xPwrSpec,data_forceAcq,data_fft = npzfile['arr_0'],npzfile['arr_1'],npzfile['arr_2'],npzfile['arr_3'],npzfile['arr_4'],npzfile['arr_5']

    # Sampling frequency - 100MHz
    fs = 1e8

    # Loop over events
    iList = -1
    while(True):
        iList += 1
        if intMode==True and iList!=0:
            value = raw_input()
            if value=='q': break
            if value=='p': iList -= 2  # previous
            if (value.isdigit()): iList = int(value) # go to entry
        if iList >= elist.GetN(): break

        entry = waveTree.GetEntryNumber(iList);
        waveTree.LoadTree(entry)
        waveTree.GetEntry(entry)
        nChans = waveTree.channel.size()
        nWaves = waveTree.MGTWaveforms.size()
        numPass = waveTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = waveTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))

        # Loop over hits passing cuts
        hitList = (iH for iH in xrange(nChans) if waveTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:

            # ------------------------------------------------------------------------

            # Load data
            run = waveTree.run
            chan = waveTree.channel.at(iH)
            dataE = waveTree.trapENFCal.at(iH)
            dataENM = waveTree.trapENM.at(iH)
            dataMaxTime = waveTree.trapENMSample.at(iH)*10. - 4000
            signal = wl.processWaveform(waveTree.MGTWaveforms.at(iH))
            data = signal.GetWaveBLSub()
            time = signal.GetTS()
            print "%d / %d  Run %d  nCh %d  chan %d  trapENF %.1f  len %d" % (iList,nList,run,nChans,chan,dataE,len(data))

            # Load template
            template, temp_time = tempOrig, tempOrigTS
            template = template * (dataENM / tempE)
            # template = wl.wfDerivative(template)

            # Plot data and template
            p1.cla()
            p1.plot(time,data,color='blue')
            p1.set_xlabel('Time (s)')
            p1.set_ylabel('Voltage (arb)')
            p1.set_title("Run %i  Ch %i  E %.2f  ENM %.2f" % (run,chan,dataE,dataENM))

            # p2.cla()
            p1.plot(temp_time,template,color='red')
            # p2.set_xlabel('Time (s)')
            # p2.set_ylabel('Voltage (arb)')
            # p2.set_title('Template')

            # Plot ASD of data and template
            # (amplitude spectral density (ASD), which is the square root of the PSD)
            p3.cla()
            power_data, freq_psd = plt.psd(data, Fs=fs, NFFT=2048, pad_to=2048, visible=False)
            p3.loglog(freq_psd, np.sqrt(power_data), color='blue')

            power_temp, freq = plt.psd(template, Fs=fs, NFFT=2048, pad_to=2048, visible=False) # check this.
            p3.loglog(freq, np.sqrt(power_temp), color='red')

            p3.loglog(xPsd,np.sqrt(avgPsd),color='green') # DS1 forced acq. results

            p3.set_xlabel('Frequency (Hz)')
            p3.set_ylabel('ASD')
            p3.grid('on')

            # Apply a bandpass filter to the data and plot it
            # NOTE: This is tough.  There might be a feature
            #       due to rising edges, between 3e6 and 5e6 Hz.
            #       But that's only from data. Doing it w/ siggen templates yielded squat.
            # B,A = sig.butter(2, [3e6/(1e8/2),5e6/(1e8/2)], btype='bandpass')

            B,A = sig.butter(2,1e6/(1e8/2),btype='lowpass') # nice lowpass filter
            data_LowPass = sig.lfilter(B, A, data)

            B1,A1 = sig.butter(2, [1e5/(1e8/2),1e6/(1e8/2)], btype='bandpass') # attempt at bandpass
            data_BanPass = sig.lfilter(B1, A1, data)

            p4.cla()
            p4.plot(time, data_LowPass, label='lowpass')
            p4.plot(time, data_BanPass, label='bandpass')
            p4.set_title('Band passed data')
            p4.set_xlabel('Time (s)')
            p4.legend(loc=4)

            # Plot time-domain cross-correlation
            # correlated_raw = np.correlate(data, template, 'same') # ligo used 'valid'.  'full' or 'same'
            # correlated_passed = np.correlate(data_pass, template, 'same')
            # p5.plot(np.arange(0, (correlated_raw.size*1.)/fs, 1.0/fs),correlated_raw)
            # p5.set_title('Time domain cross-correlation')
            # p5.set_xlabel('Offest between data and template (s)')
            # p5.plot(np.arange(0, (correlated_passed.size*1.)/fs, 1.0/fs), correlated_passed,color='red')
            # p5.set_xlabel('Offset between data and template (s)')

            # correlated = np.correlate(data,temp,'full')
            # correlated = np.correlate(data_pass,temp,'full')
            # x_corr = np.arange(0, (correlated.size*1.)/fs, 1.0/fs)

            # -- clint -- Do a simple convolution and smooth it
            temp, tempTS = flipAndAlign(tempHuge,tempHugeTS,time,dataMaxTime,dataENM,tempE)
            # smf = gaussian_filter(temp * data,sigma=float( 10 ))
            smf = gaussian_filter(temp * data_LowPass,sigma=float(5))

            p5.cla()
            p5.plot(time,data,color='blue',alpha=0.4)
            p5.plot(time,data_LowPass,color='green')
            p5.plot(tempTS,temp,color='red')
            p5.plot(tempTS,smf,color='red')
            # p5.plot(tempTS,correlated) # use the 'same' option
            # p5.plot(x_corr,correlated)


            # Matched (Optimal) Filter, freq. domain
            # Take the FFT of the data
            # data_fft = np.fft.fft(data)

            # Take the FFT of the low pass data
            data_fft = np.fft.fft(data_LowPass)

            # Pad template and take FFT
            zero_pad = np.zeros(data.size - template.size)
            template_padded = np.append(template, zero_pad)
            template_fft = np.fft.fft(template_padded)

            # Match FFT frequency bins to PSD frequency bins
            # power_noise, freq_psd_noise = plt.psd(data_forceAcq, Fs=fs, NFFT=2048, pad_to=2048, visible=False)
            # power_noise, freq_psd_noise = plt.psd(data[:256], Fs=fs, visible=False) # 512 or 256?
            power_noise, freq_psd_noise = avgPsd, xPsd # from file

            # just for fun, calculate the integral of the 3 different noise spectra
            # forceAcqSingleBL = np.sum(power_noise)
            # data_bl_psd,_ = plt.psd(data[:512], Fs=fs, NFFT=512, visible=False)
            # dataBL = np.sum(data_bl_psd)
            # avgForceAcqSpec = np.sum(avgPsd)

            datafreq = np.fft.fftfreq(data.size) * fs
            power_vec = np.interp(datafreq, freq_psd_noise, power_noise)

            # Apply the optimal matched filter
            optimal = data_fft * template_fft.conjugate() / power_vec
            optimal_time = 2*np.fft.ifft(optimal)

            # Normalize the matched filter output
            df = np.abs(datafreq[1] - datafreq[0]) # freq. bin size
            sigmasq = 2 * (template_fft * template_fft.conjugate() / power_vec).sum() * df
            sigma = np.sqrt(np.abs(sigmasq))
            SNR = abs(optimal_time) / (sigma)

            opE = np.amax(SNR)
            opEE = opE/dataE

            # Plot the result
            p6.cla()
            p6.plot(time, SNR)
            p6.set_title('Optimal Matched Filter, opE %.1f  opE/E %.2f' % (opE,opEE))
            p6.set_xlabel('Offset time (s)')
            p6.set_ylabel('SNR')

            plt.tight_layout()
            plt.pause(scanSpeed)

            # ------------------------------------------------------------------------


def flipAndAlign(tempOrig,tempOrigTS,dataTS,dataMaxTime,dataENM,tempE):
    temp, tempTS = tempOrig, tempOrigTS
    temp = temp * (dataENM / tempE)
    temp = np.flip(temp,0)
    tempMaxTime = np.argmax(temp)*10 # find max after flipping
    tempTS = tempTS - (tempMaxTime - dataMaxTime)
    idx = np.where((tempTS >= dataTS[0]-5) & (tempTS <= dataTS[-1]+5))
    temp, tempTS = temp[idx], tempTS[idx]
    return temp, tempTS


if __name__ == "__main__":
    main(sys.argv[1:])
