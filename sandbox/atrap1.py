#!/usr/local/bin/python
import sys, time
import numpy as np
import matplotlib.pyplot as plt

def main():

    file5 = np.load("./data/ds5exampleWaveform5.npz")
    dataTS, data, data5_denoised, dataENM, dataNoise = file5['arr_0'], file5['arr_1'], file5['arr_2'], file5['arr_3'], file5['arr_5']

    xAmp, xMu, xSig, xTau = dataENM, 10000., 600., -72000.  # xMu should be from GAT data

if __name__ == "__main__":
    main()