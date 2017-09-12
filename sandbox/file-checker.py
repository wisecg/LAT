#!/usr/common/usg/software/python/2.7.9/bin/python
from ROOT import TFile, TTree

def main():
    """ Two ideas here:
    1) Scan over all lat-related files, both BG and Cal:
        - skims
        - wave-skims
        - split-skims
        - split-lats
        - merged lats
    and make sure they all exist, all have tree entries, etc.

    2) Scan over calibration run lists from calDB, and see how many
    entries we have per channel.  Could help to find the sweet spot
    of how many cal runs we need to process per range.

    also, should it check if we have events down to 0.7 kev or whatever?
    """
    print "hi"

if __name__ == "__main__":
    main()