#!/usr/bin/env python
"""
    Works with March 14, 2018 version of siggen
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from siggen import PPC
import seaborn as sns
sns.set(style='darkgrid', context='talk')

def main():
    outDir = os.environ['LATDIR']
    #Set up an instance of a PPC detector
    det = PPC("{}/data/conf/p1_new.config".format(outDir), wf_padding=10000)

    # Don't generate WFs beyond this
    print("Max Length (Z): {}, Max Radius (R): {}".format(det.detector_length, det.detector_radius))

    #Plot waveforms from a few different r/z locations
    phi = 0.
    rs = np.linspace(0,34,2)
    zs = np.linspace(5,50,2)
    # amps = [2]
    # Amplitudes from 2 to 20 ADC (1 - 10 keV)
    # amps = [2*i for i in range(1,11)]
    amps = [1,2,3,4,5,6,7,8,9,10]

    # Secret Parameters without wiggles
    det.lp_num = [1.]
    det.lp_den = [1.,-1.95933813 ,0.95992564]
    det.hp_num = [1.0, -1.999640634643256, 0.99964063464325614]
    det.hp_den = [1, -1.9996247480008278, 0.99962475299714171]

    dataList = []
    # fig1, ax1 = plt.subplots(figsize=(10,7))

    for amp in amps:
        for i, r in enumerate(rs):
            for j, z in enumerate(zs):
                if r == 0 and z == 50: continue
                # Current version of SigGen can't generate more than 1000 samples
                wf_proc = det.MakeSimWaveform(r, phi, z, amp, 150, 0.1, 1000)
                # ax1.plot(wf_proc, label = "R = {:0.2f},  Z = {:.2f} ".format(r, z))
                dataMap = {}
                dataMap['amp'] = amp
                dataMap['z'] = z
                dataMap['r'] = r
                # for idx, sample in enumerate(wf_proc):
                    # dataMap['w{}'.format(idx)] = sample
                dataMap['waveform'] = wf_proc.tolist()
                dataList.append(dataMap)

    # ax1.set_xlabel("Time [{} ns steps]".format(det.time_step_size))
    # ax1.set_ylabel("Amplitude [ADC]")
    # ax1.legend()
    # plt.tight_layout()
    # plt.show()

    # Save waveforms for later -- remember, only saving 1000 samples, ~500 before and ~500 after
    df = pd.DataFrame.from_dict(dataList)
    print(df.head())
    df.to_hdf('{}/data/SigGen_WFs_low.h5'.format(outDir), 'skimTree')


if __name__ == "__main__":
    main()
