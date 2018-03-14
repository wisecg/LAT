#!/usr/bin/env python
"""
    Works with March 14, 2018 version of siggen
    Only works with python3!
"""

import numpy as np
import matplotlib.pyplot as plt
from siggen import PPC
import seaborn as sns
sns.set(style='darkgrid', context='talk')
def main():
    #Set up an instance of a PPC detector
    det = PPC("../data/conf/p1_new.config", wf_padding=10000)

    # Don't generate WFs beyond this
    print("Max Length (Z): {}, Max Radius (R): {}".format(det.detector_length, det.detector_radius))

    #Plot waveforms from a few different r/z locations
    phi = 0.
    rs = np.linspace(0,30,3)
    zs = np.linspace(5,50,3)

    # Secret Parameters without wiggles
    det.lp_num = [1.]
    det.lp_den = [1.,-1.95933813 ,0.95992564]
    det.hp_num = [1.0, -1.999640634643256, 0.99964063464325614]
    det.hp_den = [1, -1.9996247480008278, 0.99962475299714171]

    fig1, ax1 = plt.subplots(figsize=(10,7))
    for i, r in enumerate(rs):
        for j, z in enumerate(zs):
            wf_proc = det.MakeSimWaveform(r, phi, z, 100., 500, 0.5, 1000)
            ax1.plot(wf_proc, label = "R = {:0.2f},  Z = {:.2f} ".format(r, z))

    ax1.set_xlabel("Time [{} ns steps]".format(det.time_step_size))
    ax1.set_ylabel("Amplitude [ADC]")
    ax1.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
