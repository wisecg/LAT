# DCR says that you can compare the onBoardE with the trapENFCal.
# For noise events, they tend to fluctuate (onBoardE is like, random or something?)
# For physics events, they tend to have "similar" values.
# So the DIFFERENCE between the two quantities should indicate if an event is physics or not.

# for the ext pulser data, plot 'trapENFCal - onBoardE' and try to see if the low energy events
# from the external pulser have a consistent value.

# Also somehow this should show the detector threshold.

# It also occurs to me that a calibration measurement could have super-precise
# energy thresholds determined.


def getEff():
    """ Efficiency vs. energy, each detector in Test 3"""
    from ROOT import TChain, GATDataSet

    f = plt.figure()
    plt.cla()

    extPulserInfo = calInfo.GetSpecialList()["extPulserInfo"]
    syncChan = wl.getChan(0,10,0) # 672

    dsNum, modNum, calIdx = 0, 1, 33
    calDB = db.TinyDB('../calDB.json')
    pars = db.Query()
    fsD = ds.getDBRecord("fitSlo_ds%d_idx%d_m%d_Peak" % (dsNum, calIdx, modNum), False, calDB, pars)

    bkgIdx = 75 # runs 6887-6963
    thD = ds.getDBRecord("thresh_ds%d_bkgidx%d" % (dsNum,bkgIdx), False, calDB, pars)

    for pIdx in [19,20,21]:
    # for pIdx in [19]:

        runList = calInfo.GetSpecialRuns("extPulser",pIdx)
        attList = extPulserInfo[pIdx][0]
        extChan = extPulserInfo[pIdx][-1]
        fsCut = fsD[extChan][2] # 90% value (used in LAT3)

        effVals, threshVals, trigVals = [], [], []
        eneVals, sloVals, rtVals = [], [], []
        for i, run in enumerate(runList):
            if run in [7225, 7233]:
                continue

            # elogs: "20 Hz, 150 second runs"
            gds = GATDataSet(run)
            runTime = gds.GetRunTime() # sec
            pulseRate = 20 # Hz

            fileList = ds.getLATRunList([run],"%s/lat" % (ds.specialDir))
            latChain = TChain("skimTree")
            for f in fileList:
                latChain.Add("%s/lat/%s" % (ds.specialDir,f))

            tNames = ["Entry$","mH","channel","trapENFCal","fitSlo","den90","den10","threshKeV","threshSigma"]
            theCut = "(channel==%d || channel==%d) && mH==2" % (syncChan, extChan) # enforce correct sync
            tVals = wl.GetVX(latChain,tNames,theCut)
            nPass = len(tVals["Entry$"])

            enfArr = [tVals["trapENFCal"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
            sloArr = [tVals["fitSlo"][i] for i in range(nPass) if tVals["channel"][i]==extChan]
            rtArr  = [tVals["den90"][i] - tVals["den10"][i] for i in range(nPass) if tVals["channel"][i]==extChan]

            if len(enfArr)==0:
                print("Run %d, No hits in channel %d found.  Continuing ..." % (run,extChan))
                continue

            eneVals.extend(enfArr)
            sloVals.extend(sloArr)
            rtVals.extend(rtArr)
