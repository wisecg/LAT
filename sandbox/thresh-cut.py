#!/usr/bin/env python
import sys, itertools
sys.argv.append("-b")
import DataSetInfo as ds
import tinydb as db
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import scipy.special as sp

def main():

    threshCut()


def threshCut():
    """ ./lat-jobs.py -threshCut
        Applies a threshold cut.
    """
    import numpy as np
    import pandas as pd
    import os

    threshCut = 0.9

    for dsNum in range(0,6):

        if dsNum==2: continue

        print("dataset",dsNum)
        df = pd.read_hdf("./data/ThreshDS%d_Processed.h5" % dsNum, 'threshTree')

        goodRuns, badRuns, goodRunErfs = {}, {}, {}
        for chan in df:
            col = int(chan)
            goodRuns[col], badRuns[col], goodRunErfs[col] = [], [], []

            for idx, vals in enumerate(df.loc[:,chan]):

                if np.isnan(vals).any(): continue # skip NaN run ranges where data wasn't collected for the channel

                thr, sigma, hi, lo = vals[0], vals[1], int(vals[2]), int(vals[3])

                if vals[0] <= threshCut:
                    goodRuns[col].append([hi,lo])
                    goodRunErfs[col].append([hi,lo,thr,sigma])
                else:
                    badRuns[col].append([hi,lo])


        # Now make a list of channels to cut in each bad run.

        runList = []
        for ch in badRuns:
            for pair in badRuns[ch]:
                for run in range(pair[0],pair[1]+1):
                    runList.append(run)
        runSet = set(runList)
        runList = sorted(list(runSet))

        runChanPairs = {}
        for run in runList:

            if run not in runChanPairs:
                runChanPairs[run] = []

            for ch in badRuns:
                for pair in badRuns[ch]:
                    if run >= pair[0] and run <= pair[1]:
                        runChanPairs[run].append(ch)
                        # print(run,ch)

        text_file = open("threshCut_v1.txt", "a")
        for run in sorted(runChanPairs):
            sarr = [str(a) for a in runChanPairs[run]]
            text_file.write(str(run) + ' ' + (' '.join(sarr)) + "\n")
        text_file.close()

if __name__=="__main__":
    main()