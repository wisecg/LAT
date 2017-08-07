#!/usr/local/bin/python
"""
===================== LAT2.py =====================

Calibrate energy parameters in LAT data.

Two modes:
    -cal : Scans a calibration range and updates
           parameters in a pandas file.
    -upd : Updates an input file with data from the
           calibration database file.

v1: 07 Aug 2017

========= C. Wiseman (USC), B. Zhu (LANL) =========
"""
import sys
import numpy as np

def main(argv):

    print "======================================="
    print "LAT2 started:",time.strftime('%X %x %Z')
    startT = time.clock()

    # let's get some margs
    cal, upd = 0, 0
    dsNum, subNum, runNum = -1, -1, -1
    pathToInput, pathToOutput, manualInput, manualOutput = ".", ".", "", ""

    if len(argv)==0: return
    for i,opt in enumerate(argv):
        if opt == "-cal":
            cal = True
        if opt == "-r":
            rangeMode, dsNum, subNum = True, int(argv[i+1]), int(argv[i+2])
            print "Scanning DS-%d sub-range %d" % (dsNum, subNum)


if __name__ == "__main__":
    main(sys.argv[1:])