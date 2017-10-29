#!/usr/bin/env python
import imp
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')
wl = imp.load_source('waveLibs','../waveLibs.py')
import tinydb as db
"""
tinyDB items are nested dicts:
{"key":key, "vals":vals}

Cal Tables (gives run coverages):
    key: ds[DS]_calIdx.
    vals: {[idx]:[cal lo, cal hi, cov lo, cov hi]}
    Print one with:
    wl.getDBCalTable(dsNum, verbose=True)

Cal Records (calib consts for each channel in each calIdx)
    key: ds[DS]_idx[n]
    vals: {[chan]:[trapENF, fitAmp, latAF, latAFC]}
    - channel list comes from DataSetInfo.py

Cut records:
    key: [Name]_ds[i]_idx[j]_module[k]_[descriptor].
        Names are: "riseNoise", "fitSlo", "bcMax", "pol2", and "pol3"
        idx is the calIdx
        Descriptor is two numbers (energy range), "continuum", or "peak"
            Continuum range is 5-50 kev
            Peak is 236-240
        descriptors: Peak, Continuum, 50_90, 90_130, 130_170, 170_210
    vals: {[chan]:[1%, 5%, 90%, 95%, 99%]}
    - channel list comes from DataSetInfo.py
"""
def main():

    # testWLFunctions()
    # return

    # iterate over the DB
    # nRec = 0
    calDB = db.TinyDB('../calDB.json')
    keys = []
    types = []
    for item in calDB:

        d = dict(item)
        key = d["key"]
        vals = d["vals"]

        tmp = key.split("_")
        tmp = [str(t) for t in tmp]

        if tmp[0]=="fitSlo":
            # print tmp

            print vals

            return


        # nRec += 1
        # if nRec > 10: break

    types = set(types)
    print types


def testWLFunctions():

    # wl.setDBCalTable()
    # wl.getDBKeys()
    # wl.getDBCalRecord("ds1_idx0")
    # wl.getDBCalRecord("ds1_calIdx")
    # wl.getDBCalRecord("fitSlo_ds1_idx55_m1_170_210")
    # wl.delDBRecord("ds1_idx0")
    # wl.getDBCalTable(5)
    # wl.getDBRunCoverage(1,9999)

    cal = ds.CalInfo()

    # get a cal index for a run
    # key, run = "ds1_m1",10770
    # print cal.GetCalIdx(key,run)

    # generate a list of cal runs for a given index
    # key, idx = "ds3_m1",4
    # print cal.GetCalList(key,idx,runLimit=10)

    # generate cal runs for a given dataset
    # key = "ds2_m1"
    # for idx in range(cal.GetIdxs(key)):
        # print cal.GetCalList(key,idx,runLimit=10)

    # generate all possible cal runs
    # for key in cal.GetKeys():
        # print key
        # for idx in range(cal.GetIdxs(key)):
            # print cal.GetCalList(key,idx,runLimit=10)


if __name__=="__main__":
    main()
