#!/usr/local/bin/python
"""
===================== LAT3.py =====================

PLACEHOLDER CODE.  THIS SHOULD BE A LAT SKIMMER,
WHICH READS DATASETINFO.PY AND PRODUCES CHANNEL-SPECIFIC
FILES WITH RUN CUTS APPLIED.

v1: 07 Aug 2017

========= C. Wiseman (USC), B. Zhu (LANL) =========
"""
import sys, time
import numpy as np

def main(argv):

    print "========================================"
    print "LAT3 started:",time.strftime('%X %x %Z')
    startT = time.clock()
    stopT = time.clock()
    print "Stopped:",time.strftime('%X %x %Z'),"\nProcess time (min):",(stopT - startT)/60

if __name__ == "__main__":
    main(sys.argv[1:])