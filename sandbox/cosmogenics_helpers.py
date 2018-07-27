#!/usr/bin/env python
import sys, os, time
import numpy as np
import ROOT
import dsi
det = dsi.DetInfo()
bkg = dsi.BkgInfo()
cal = dsi.CalInfo()

"""

    Some helper functions to grab common detector lists across datasets
    Also outputs a .txt file with run,startTime,runTime

"""

# process_MJD_as_built_mage_results -o /mnt/mjdDisk1/Majorana/users/psz/MAGE/Processed/surf/Pb210/1_0.99_0.99 -c /mnt/mjdDisk1/Majorana/users/psz/MAGE/config/det_config_DS5.json -t 1.0 0.99 0.99 /mnt/mjdDisk1/Majorana/users/psz/MAGE/Simulations/surf/Pb210/MJDem_DUPTFE_Pb210_surf_p051.root


# Module 1 DS0-6: ['C1P1D2', 'C1P1D3', 'C1P1D4', 'C1P2D2', 'C1P2D3', 'C1P3D4', 'C1P5D3', 'C1P6D3', 'C1P7D2', 'C1P7D3']
# Module 1 DS1-6: ['C1P1D2', 'C1P1D3', 'C1P1D4', 'C1P2D2', 'C1P2D3', 'C1P3D2', 'C1P3D3', 'C1P3D4', 'C1P5D3', 'C1P6D1', 'C1P6D3', 'C1P6D4', 'C1P7D2', 'C1P7D3', 'C1P7D4']
# Module 2 DS4-6: ['C2P1D4', 'C2P3D1', 'C2P3D2', 'C2P6D2', 'C2P7D3'] # Rejected C2P5D3 because of noise

#"C==1&&(P==1&&D==2)||(P==1&&D==3)||(P==1&&D==4)||(P==2&&D==2)||(P==2&&D==3)||(P==3&&D==4)||(P==5&&D==3)||(P==6&&D==3)||(P==7&&D==2)||(P==7&&D==3)"

#C==1&&(P==1&&D==2)||(P==1&&D==3)||(P==1&&D==4)||(P==2&&D==2)||(P==2&&D==3)||(P==3&&D==2)||(P==3&&D==3)||(P==3&&D==4)||(P==5&&D==3)||(P==6&&D==1)||(P==6&&D==3)||(P==6&&D==4)||(P==7&&D==2)||(P==7&&D==3)||(P==7&&D==4)

# M1List = [112, 113, 114, 122, 123, 132, 133, 134, 153, 161, 163, 164, 172, 173, 174]
# M2List = [214, 231, 232, 262, 273]

# Exposure Information
# DS-0
# Enriched (kg-d): 426.3625   No cuts 454.604  - PSA95  28.2419  - Burst 0.0000    (Tot: -28.2419)
# Natural (kg-d) : 158.9550   No cuts 164.528  - PSA95 5.5729   - Burst 0.0000    (Tot: -5.5729)

# DS-1
# Enriched (kg-d): 631.5396   No cuts 656.151  - PSA95  18.1773  - Burst 6.4342    (Tot: -24.6116)
# Natural (kg-d) : 29.5235    No cuts 62.863   - PSA95 33.3395  - Burst 0.0000    (Tot: -33.3395)

# DS-2
# Enriched (kg-d): 104.3238   No cuts 105.821  - PSA95  0.2202   - Burst 1.2773    (Tot: -1.4975)
# Natural (kg-d) : 5.1778     No cuts 10.522   - PSA95 5.3444   - Burst 0.0000    (Tot: -5.3444)

# DS-3
# Enriched (kg-d): 334.9536   No cuts 368.563  - PSA95  29.1560  - Burst 4.4531    (Tot: -33.6091)
# Natural (kg-d) : 48.6482    No cuts 81.738   - PSA95 31.2533  - Burst 1.8369    (Tot: -33.0903)

# DS-4
# Enriched (kg-d): 50.7790    No cuts 102.849  - PSA95  52.0701  - Burst 0.0000    (Tot: -52.0701)
# Natural (kg-d) : 23.1491    No cuts 73.844   - PSA95 50.6944  - Burst 0.0000    (Tot: -50.6944)

# DS-5A
# Enriched (kg-d): 798.8447   No cuts 1244.702 - PSA95  435.5196 - Burst 10.3379   (Tot: -445.8575)
# Natural (kg-d) : 337.4632   No cuts 616.633  - PSA95 276.9703 - Burst 2.1998    (Tot: -279.1701)

# DS-5B
# Enriched (kg-d): 624.1460   No cuts 672.765  - PSA95  44.3761  - Burst 4.2432    (Tot: -48.6192)
# Natural (kg-d) : 226.3708   No cuts 335.903  - PSA95 108.9848 - Burst 0.5472    (Tot: -109.5319)

# DS-5C
# Enriched (kg-d): 173.0972   No cuts 174.474  - PSA95  0.9417   - Burst 0.4350    (Tot: -1.3767)
# Natural (kg-d) : 63.3547    No cuts 86.029   - PSA95 20.3367  - Burst 2.3379    (Tot: -22.6747)

# DS-6
# Enriched (kg-d): 1144.5646   No cuts 1359.229 - PSA95  125.5439 - Burst 89.1208   (Tot: -214.6648)
# Natural (kg-d) : 389.5847   No cuts 509.397  - PSA95 117.5180 - Burst 2.2943    (Tot: -119.8123)
#
# Totals for DS: [0, 1, 2, 3, 4, '5A', '5B', '5C', 6]
# Enriched (kg-y): 11.7416
# Natural (kg-y) : 3.5105

def main():
    # grabDetList()
    hDictEnr, hDictNat = {}, {}
    dTypeList = ['final', 'frb', 'fr', 'th', 'lat']
    for ds in range(0,7):
        # fOut = ROOT.TFile('fHistograms_DS{}.root'.format(ds), 'RECREATE')
        # fOut.cd()
        # for dType in dTypeList:
            # hDictEnr[ds], hDictNat[ds] = grabRawSpectra(ds, dType)
            # hDictEnr[ds].Write()
            # hDictNat[ds].Write()
        grabRunList(ds, True)
        # fOut.Close()



def grabRunList(ds=6, bSave=False):
    """
        Scans through all background runs in a dataset to print out
    """
    from ROOT import GATDataSet, GATTimeInfo, TFile, TTree

    runList = []
    for run in bkg.getRunList(ds):
            runList.append(run)
    runList = sorted(runList) # if the cal run isn't before the bkg runs, it's not the first run in this list

    print("Scanning Dataset {}".format(ds))

    # use GDS once just to pull out the path.
    gds = GATDataSet()
    runPath = gds.GetPathToRun(runList[0],GATDataSet.kGatified)
    filePath = '/'.join(runPath.split('/')[:-1])

    # track time per run so we can identify slowdowns on PDSF
    start = time.time()

    startList = []
    timeList = []

    # begin loop over runs
    for idx, run in enumerate(runList):
        # print progress
        f = np.fabs(100*idx/len(runList) % 10)
        if f < 0.5:
            print("%d/%d (run %d) %.1f%% done." % (idx, len(runList), run, 100*idx/len(runList)))

        # make sure file exists and it's not blind before trying to load it
        fname = filePath + "/mjd_run%d.root" % run
        if not os.path.isfile(fname):
            print("Couldn't find run",run,":",fname)
            continue
        if not os.access(fname, os.R_OK):
            print("File is blind:",fname)
            continue
        tf = TFile(fname)
        mjdTree = tf.Get("mjdTree")
        mjdTree.GetEntry(0)

        tInfo = mjdTree.timeinfo
        startTime = tInfo.startTime
        stopTime = tInfo.stopTime
        startList.append(int(startTime))
        # timeList.append(stopTime - startTime) # Start Time
        timeList.append(int(stopTime)) # Stop Time

        # For debugging
        # print(startTime, stopTime, stopTime-startTime)

        # Close files after finished
        tf.Close()

    timeElapsed = time.time()-start

    print(len(runList), len(startList), len(timeList))

    if bSave:
        with open('./data/DS{}_RunTimeList.txt'.format(ds), 'w') as f:
            for run,start,live in zip(runList, startList, timeList):
                f.write('{},{},{}\n'.format(run,start,live))


    print("Elapsed: %.4f sec, %.4f sec/run" % (timeElapsed, timeElapsed/len(runList)))


def grabDetList():
    """
        Scans through good detector lists of a dataset to print out common enriched detectors
    """
    dsList = [0, 1, 2, 3, 4, 5]

    l0 = det.getGoodChanList(0, mod=1, detType='Enr')
    cpd0 = ['C{}P{}D{}'.format(*det.getChanCPD(0, chan)) for chan in l0]
    l1 = det.getGoodChanList(1, mod=1, detType='Enr')
    cpd1 = ['C{}P{}D{}'.format(*det.getChanCPD(1, chan)) for chan in l1]
    l2 = det.getGoodChanList(2, mod=1, detType='Enr')
    cpd2 = ['C{}P{}D{}'.format(*det.getChanCPD(2, chan)) for chan in l2]
    l3 = det.getGoodChanList(3, mod=1, detType='Enr')
    cpd3 = ['C{}P{}D{}'.format(*det.getChanCPD(3, chan)) for chan in l3]
    l5 = det.getGoodChanList(5, mod=1, detType='Enr')
    cpd5 = ['C{}P{}D{}'.format(*det.getChanCPD(5, chan)) for chan in l5]

    a0 = np.intersect1d(cpd0, cpd1, assume_unique=True)
    a01 = np.intersect1d(cpd0, cpd3, assume_unique=True)
    a02 = np.intersect1d(cpd0, cpd5, assume_unique=True)
    a1 = np.intersect1d(cpd1, cpd3, assume_unique=True)
    a2 = np.intersect1d(cpd1, cpd5, assume_unique=True)
    a3 = np.intersect1d(cpd3, cpd5, assume_unique=True)
    a4 = np.intersect1d(a2, a3, assume_unique=True)

    a03 = np.intersect1d(a0, a01, assume_unique=True)
    a04 = np.intersect1d(a01, a02, assume_unique=True)
    a05 = np.intersect1d(a03, a04, assume_unique=True)

    print(cpd0, cpd1, cpd2, cpd3, cpd5)
    print(a4.tolist())
    print(a05.tolist())

    l4 = det.getGoodChanList(4, mod=2, detType='Enr')
    cpd4 = ['C{}P{}D{}'.format(*det.getChanCPD(4, chan)) for chan in l4]

    l52 = det.getGoodChanList(5, mod=2, detType='Enr')
    cpd52 = ['C{}P{}D{}'.format(*det.getChanCPD(5, chan)) for chan in l52]

    a52 = np.intersect1d(cpd4, cpd52, assume_unique=True)
    print(a52.tolist())


def grabRawSpectra(dsNum=6, dType='lat'):
    print('Scanning DS{}, type={}'.format(dsNum, dType))

    skimTree = ROOT.TChain("skimTree")
    if dType == 'lat':
        skimTree.Add('{}/latSkimDS{}_*.root'.format(dsi.latDir, dsNum))
    elif dType == 'th':
        skimTree.Add('{}/th_ds{}_*.root'.format(dsi.cutDir+'/th', dsNum))
    elif dType == 'fr':
        skimTree.Add('{}/fr_ds{}_*.root'.format(dsi.cutDir+'/fr95', dsNum))
    elif dType == 'frb':
        skimTree.Add('{}/frb95_ds{}_*.root'.format(dsi.cutDir+'/frb95', dsNum))
    elif dType == 'final':
        skimTree.Add('{}/final95_DS{}*.root'.format(dsi.cutDir+'/final95', dsNum))

    hEnr = ROOT.TH1D("hEnr_DS{}_{}".format(dsNum, dType), "DS{} Enriched ({})".format(dsNum, dType), 2500, 0, 250)
    hNat = ROOT.TH1D("hNat_DS{}_{}".format(dsNum, dType), "DS{} Natural ({})".format(dsNum, dType), 2500, 0, 250)

    skimTree.Project("hEnr_DS{}_{}".format(dsNum, dType), "trapENFCal", "isEnr && trapENFCal > 0 && trapENFCal < 250")
    skimTree.Project("hNat_DS{}_{}".format(dsNum, dType), "trapENFCal", "isNat && trapENFCal > 0 && trapENFCal < 250")

    return hEnr, hNat


if __name__ == '__main__':
    main()
