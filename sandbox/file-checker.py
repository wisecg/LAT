#!/usr/common/usg/software/python/2.7.9/bin/python
from ROOT import TFile, TTree
import os.path

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

    also, this should be moved to job-panda once it's working well.
    """

    dsMap = {0:75,1:51,2:7,3:24,4:18,5:112}

    for ds in dsMap:
        for sub in range(dsMap[ds]+1):

            # make sure bg skims exist
            fileName = "/global/homes/w/wisecg/project/bg-skim/skimDS%d_%d_low.root" % (ds,sub)
            if not os.path.isfile(fileName):
                print "file not found! name:", fileName
                continue
            f = TFile(fileName)
            t = f.Get("skimTree")
            n = t.GetEntriesFast()
            print "DS %d  sub %d  entries %d" % (ds, sub, n)
            if n==0:
                print "no entries found! file:", fileName
                continue






if __name__ == "__main__":
    main()