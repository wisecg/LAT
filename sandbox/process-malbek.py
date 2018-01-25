#!/usr/bin/env python3
import sys, os, imp, glob
# sys.argv.append("-b")
# ds = imp.load_source('DataSetInfo',os.environ['LATDIR']+'/DataSetInfo.py')
# import matplotlib.pyplot as plt
# import numpy as np

def main():
    explore()


def explore():
    """
    # test1 has really fast pulser waveforms (faster than physics)
    test1_t4.root
        -> fTree (124340 entries)
    test1_t4_wf.root
        -> tree  (124340 entries)

    # test3 is slow(er)
    test3_t4.root
        -> fTree (147962 entries)
    test3_t4_wf.root
        -> tree (147962 entries)

    ThesisAll_newCal.root
        -> there's a buncha stuff here.  I guess it's mainly bkg data.
    t4_wf_all.root

    LinearCrateSpecial.root
    pb_t4_wf.root

    """
    from ROOT import TFile, TTree


def get_threshold_list():
    """ From Graham (gkg@princeton.edu) """
     return [ 706.2,
              575.1,
              447.7,
              313.6,
              229.1,
              207.2,
              224.8,
              180.2]


def iswt(coefficients, wavelet):
    """ From Graham (gkg@princeton.edu) """
    """
      Input parameters:

        coefficients
          approx and detail coefficients, arranged in level value
          exactly as output from swt:
          e.g. [(cA1, cD1), (cA2, cD2), ..., (cAn, cDn)]

        wavelet
          Either the name of a wavelet or a Wavelet object

    """
    output = coefficients[0][0].copy() # Avoid modification of input data

    #num_levels, equivalent to the decomposition level, n
    num_levels = len(coefficients)
    for j in range(num_levels,0,-1):
        step_size = int(math.pow(2, j-1))
        last_index = step_size
        _, cD = coefficients[num_levels - j]
        for first in range(last_index): # 0 to last_index - 1

            # Getting the indices that we will transform
            indices = numpy.arange(first, len(cD), step_size)

            # select the even indices
            even_indices = indices[0::2]
            # select the odd indices
            odd_indices = indices[1::2]

            # perform the inverse dwt on the selected indices,
            # making sure to use periodic boundary conditions
            x1 = pywt.idwt(output[even_indices], cD[even_indices], wavelet, 'per')
            x2 = pywt.idwt(output[odd_indices], cD[odd_indices], wavelet, 'per')

            # perform a circular shift right
            x2 = numpy.roll(x2, 1)

            # average and insert into the correct indices
            output[indices] = (x1 + x2)/2.

    return output


def apply_threshold(output, scaler = 1., input=None):
    """ From Graham (gkg@princeton.edu) """
   """
       output is a list of vectors (cA and cD, approximation
       and detail coefficients) exactly as you would expect
       from swt decomposition.
          e.g. [(cA1, cD1), (cA2, cD2), ..., (cAn, cDn)]

       If input is none, this function will calculate the
       tresholds automatically for each waveform.
       Otherwise it will use the tresholds passed in, assuming
       that the length of the input is the same as the length
       of the output list.
       input looks like:
          [threshold1, threshold2, ..., thresholdn]

       scaler is a tuning parameter that will be multiplied on
       all thresholds.  Default = 1 (0.8?)
   """

   for j in range(len(output)):
      cA, cD = output[j]
      if input is None:
        dev = numpy.median(numpy.abs(cD - numpy.median(cD)))/0.6745
        thresh = math.sqrt(2*math.log(len(cD)))*dev*scaler
      else: thresh = scaler*input[j]
      cD = pywt.thresholding.hard(cD, thresh)
      output[j] = (cA, cD)


def grahamsWFProcessFunction():
    """ From Graham (gkg@princeton.edu) """

    # -----------------------------------------------
    # do something like this for each wf in the tree
    # ------------------------------------------------

    #set wavelet parameters
    wl = pywt.Wavelet('haar')
    levels = 8
    threshold_list = get_threshold_list()

    # get vector data
    waveform_vector = waveform.GetVectorData()

    # perform 8 level SWT
    swt_output = pywt.swt(waveform_vector, wl, level=levels)

    # threshold the SWT coefficients
    apply_threshold(swt_output, 1., threshold_list)

    # inverse transform
    cA_thresh = iswt(swt_output, wl)

    # make a new denoised waveform
    denoised_waveform = ROOT.MGTWaveform()
    wf_sampling = current_waveform.GetSamplingFrequency()
    wf_length = current_waveform.GetLength()
    denoised_waveform.SetSamplingFrequency(wf_sampling)
    denoised_waveform.SetData(cA_thresh[0:wf_length], wf_length)

    # -----------------------------------------------
    # if you want to calculate wpar, you need to find the maximum of
    # the ccd0 power spectrum
    # ------------------------------------------------

    # get the extreme of cd0 power spectrum
    _, cD_zero = swt_output[0]
    wf_length = current_waveform.GetLength()
    wf_sampling = current_waveform.GetSamplingFrequency()
    wf_period = current_waveform.GetSamplingPeriod()
    offset = 5000

    cd_power_waveform = ROOT.MGTWaveform()
    cd_power_waveform.SetSamplingFrequency(wf_sampling)
    cd_power_waveform.SetData(cD_zero[0:wf_length], wf_length)
    cd_power_waveform.SetLength(wf_length)
    cd_power_waveform *= cd_power_waveform

    max_finder.SetFindMaximum(True)
    max_finder.SetLocalMinimumTime(offset)
    max_finder.SetLocalMaximumTime(wf_length*wf_period - offset)
    max_finder.TransformInPlace(cd_power_waveform)
    cd_zero_max = max_finder.GetTheExtremumValue()

    wpar = cd_zero_max]/energy**2

if __name__=="__main__":
    main()