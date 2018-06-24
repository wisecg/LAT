#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')

import waveLibs as wl

def main():

    # Redondo note: I calculated the flux using gae = 0.511*10^-10
    # for other values of gae use:
    # FLUX = Table*[gae/(0.511*10^-10)]^2
    gaeR = 5.11e-11
    gaeP = 4.35e-12 # PandaX PRL, Nov 2017 (gae2.py)
    gae = 1
    rat = gaeR/gaeR
    redondoScale = 1e19 * rat**2 # convert table to [cts / (keV cm^2 d)]

    tmpE, tmpF = [], []
    with open("../data/redondoFlux.txt") as f1:
        lines = f1.readlines()[11:]
        for i, line in enumerate(lines[:-1]):
            data = line.split()
            tmpE.append(float(data[0]))
            tmpF.append(float(data[1]))

    axE, axF, axD = [], [], []
    for i in range(len(tmpE)-10):
        axE.append(tmpE[i])
        axF.append(tmpF[i])
        der = (tmpF[i+10] - tmpF[i]) / (tmpE[i+10] - tmpE[i])
        axD.append(der)

    # # profit
    # msList = []
    # msThresh = 15
    # maxtab,_ = wl.peakdet(axD, msThresh)
    # for iMax in range(len(maxtab)):
    #     idx = int(maxtab[iMax][0])
    #     val = maxtab[iMax][1]
    #     # msList.append(dataTS[idx])
    #     print("%d  idx %d  ene %.2f  val %.2f  thresh %.2f" % (iMax, idx, axE[idx], val, msThresh))
    #     plt.axvline(axE[idx], lw=1, c='r', label="%.2f" % axE[idx])

    # plt.plot(axE, axD, '.', ms=1.)
    plt.plot(axE, axF, '.', ms=1.)

    # plt.legend(ncol=3, fontsize=14)
    plt.xlabel("Energy (keV)", ha='right', x=1)
    plt.tight_layout()

    plt.show()



if __name__=="__main__":
    main()