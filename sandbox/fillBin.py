import random
import numpy as np
import matplotlib.pyplot as plt

def main():
    # manually increment a numpy array histogram
    # similar to TH1::Fill(x)

    mu, sigma = -1, 1
    vals = np.random.normal(mu, sigma, 10000000)

    fLo, fHi, fpb = -10, 10, 0.01
    nbf = int((fHi-fLo)/fpb)+1
    hist = np.zeros(nbf)

    for v in vals:
        if v < fLo or v > fHi:
            continue

        binIdx = int( nbf * (v-fLo)/(fHi-fLo) ) # same formula as TAxis::FindBin
        hist[binIdx] += 1

    x = np.arange(fLo-fpb,fHi,fpb)

    mu = np.mean(vals)
    max = x[np.argmax(hist)]

    fig = plt.figure()
    plt.plot(x, hist, ls='steps', c='b')
    plt.axvline(mu, c='r', label='vals mean:%.4f' % mu)
    plt.axvline(mu, c='g', label='hist mean:%.4f' % max)
    plt.legend(loc=1)

    plt.xlabel("fitSlo")
    plt.show()


if __name__=="__main__":
    main()