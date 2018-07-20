#!/usr/bin/env python3
import sys
import numpy as np
import waveLibs as wl
import dsi
import matplotlib.pyplot as plt
plt.style.use('../pltReports.mplstyle')

def main():

    plotRes()
    # getSigma_v2()


def getSigma(E, ds=None):
    """ Get the MJ energy resolution.
    This is just used to set an intial guess on the resolution
    because we let sigma float for the peaks.
    """
    if ds==0:   p = [0.147, 0.0173, 0.0003]
    elif ds==1: p = [0.136, 0.0174, 0.00028]
    elif ds==2: p = [0.143, 0.0172, 0.000284]
    elif ds==3: p = [0.162, 0.0172, 0.000297]
    elif ds==4: p = [0.218, 0.015, 0.00035]
    elif ds=='5A': p = [0.2121, 0.01838, 0.00031137]
    elif ds=='5B': p = [0.18148, 0.01690, 0.00031873]
    else:
        # use DS5B numbers for any other DS
        p = [0.18148, 0.01690, 0.00031873]

    return np.sqrt(p[0]**2 + p[1]**2 * E + p[2]**2 * E**2)


def getSigma_v2(E, ds=None, opt=""):
    """ HG resolutions, from the energy unidoc. They are given in terms of FWHM, but I want a single Gaussian width (sigma).
    """
    eRes = {
        "0" :  {"nat": [1.260e-1, 1.790e-2, 2.370e-4], "enr": [1.500e-1, 1.750e-2, 2.820e-4], "both": [1.470e-1, 1.730e-2, 3.000e-4]},
        "1" :  {"nat": [1.470e-1, 1.770e-2, 2.010e-4], "enr": [1.340e-1, 1.750e-2, 2.820e-4], "both": [1.360e-1, 1.740e-2, 2.800e-4]},
        "2" :  {"nat": [1.410e-1, 1.800e-2, 1.680e-4], "enr": [1.420e-1, 1.720e-2, 2.860e-4], "both": [1.430e-1, 1.720e-2, 2.840e-4]},
        "3" :  {"nat": [1.800e-1, 1.820e-2, 2.090e-4], "enr": [1.580e-1, 1.710e-2, 3.090e-4], "both": [1.620e-1, 1.720e-2, 2.970e-4]},
        "4" :  {"nat": [2.140e-1, 1.540e-2, 3.970e-4], "enr": [2.170e-1, 1.490e-2, 3.190e-4], "both": [2.180e-1, 1.500e-2, 3.500e-4]},
        "5A" : {"nat": [2.248e-1, 1.894e-2, 2.794e-4], "enr": [2.660e-1, 2.215e-2, 2.868e-4], "both": [2.592e-1, 2.057e-2, 3.086e-4]},
        "5B" : {"nat": [1.650e-1, 1.760e-2, 2.828e-4], "enr": [1.815e-1, 1.705e-2, 3.153e-4], "both": [1.815e-1, 1.690e-2, 3.187e-4]},
        "5C" : {"nat": [1.565e-1, 1.810e-2, 2.201e-4], "enr": [1.361e-1, 1.740e-2, 2.829e-4], "both": [1.519e-1, 1.718e-2, 2.762e-4]}
    }
    p = eRes[str(ds)][opt]
    return np.sqrt(p[0]**2 + p[1]**2 * E + p[2]**2 * E**2)/2.355


def getExpo(ds):

    f = np.load("%s/data/expo-totals-e95.npz"  % dsi.latSWDir)
    dsExpo, detExpo = f['arr_0'].item(), f['arr_1'].item()
    enrExp, natExp = dsExpo[ds][0], dsExpo[ds][1]
    return enrExp, natExp


def plotRes():

    fig, (p1, p2) = plt.subplots(1, 2, figsize=(9,5))

    dsList = [0,1,2,3,4,"5A","5B","5C"]

    x = np.arange(0, 250, 0.01)
    resWEnr = np.zeros(len(x))
    resWNat = np.zeros(len(x))

    enrTot, natTot = 0, 0
    cmap = plt.cm.get_cmap('jet', len(dsList)+2)

    for i, ds in enumerate(dsList):
        enrE, natE = getExpo(ds)
        enrTot += enrE
        natTot += natE
        resE = getSigma_v2(x, ds, "enr")
        resN = getSigma_v2(x, ds, "nat")
        resWEnr += resE * enrE
        resWNat += resN * natE

        p1.plot(x, resE, c=cmap(i), label="DS%s : %.0f kg-d, Enr" % (str(ds), enrE))
        p2.plot(x, resN, c=cmap(i), label="DS%s : %.0f kg-d, Nat" % (str(ds), natE))

    # you can analytically weight several functions, not hard at all. f = Sum e_i*f_i / Sum e_i
    resWEnr = resWEnr / enrTot
    resWNat = resWNat / natTot

    p1.plot(x, resWEnr, c='r', lw=6, label="Weighted, Enr")
    p2.plot(x, resWNat, c='r', lw=6, label="Weighted, Nat")

    p1.legend(loc=2, fontsize=10, ncol=2)
    p2.legend(loc=2, fontsize=10, ncol=2)
    p1.set_ylim(0.05, 0.2)
    p2.set_ylim(0.05, 0.2)
    p1.set_xlabel("Energy (keV)", ha='right', x=1)
    p2.set_xlabel("Energy (keV)", ha='right', x=1)
    p1.set_ylabel(r"$\mathregular{\sigma(E)}$", ha='right', y=1)
    plt.tight_layout()
    # plt.show()
    plt.savefig("%s/plots/lat-eres.pdf" % dsi.latSWDir)


if __name__=="__main__":
    main()