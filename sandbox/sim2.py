#!/usr/bin/env python3
import sys, os, shlex, glob, imp
import subprocess as sp

# load LAT libraries
import dsi
wl = imp.load_source('waveLibs',os.environ['LATDIR']+'/waveLibs.py')
cal = dsi.CalInfo()

def main():

    # processSims()
    # testData()
    processTritSim()


def processSims():
    """ Requires Micah's branches MGDO:sim and GAT:forClint to be set.
        I've built them in ~/mgsw/test/
        To enable them, export MGDODIR and GATDIR in ~/lat/envSetup.sh
    """
    # this is not the complete run list, but it could easily be extended if we need more stats.
    runList = [
        1087013,1087021,1087029,1087037,1087045,1087053,1087061,1087069,1087077,1087085,
        1087093,1087101,1087109,1087245,1087381,1087517,1087653,1087789,1087925,1088061,
        1088197,1088333,1088469,1088605,1088741,1088877,1089013,1089149,1089285,1089421,
        1089557,1089693,1089829
        ]
    outDir = "~/project/sims/xyz"
    gatDir = os.environ['GATDIR']
    app = "%s/Apps/process_MJD_as_built_mage_results" % gatDir
    config = "/global/projecta/projectdirs/majorana/users/mbuuck/sim/mageScripts/det_config_for_GAT_pp_by_det_v2.json"

    # best fit to DS3 long cal
    dlThickness = 3.0
    transPoint = 0.95
    transLevel = 0.4

    for run in runList:
        inFile = "/global/projecta/projectdirs/majorana/sim/MJDG41003Sims/MJDemonstrator/linesource/M1CalSource/A224_Z88/MJDemonstrator_linesource_A224_Z88_from_A224_Z88_to_A208_Z81_in_M1CalSource_500000_-%d.root" % run

        cmd = "%s -o %s -c %s -t %.1f %.2f %.1f %s" % (app,outDir,config,dlThickness,transPoint,transLevel,inFile)

        print(run)
        # print(cmd)
        sp.call(shlex.split(cmd))


def testData():
    from ROOT import TFile, TTree, TVector3

    geomFile = "/global/projecta/projectdirs/majorana/sim/MJDG41003Sims/MJD.json"
    fileList = sorted(glob.glob("%s/sims/xyz/*.root" % ds.dataDir))

    for f in fileList[:1]:
        print(f)
        tf = TFile(f)
        simTree = tf.Get("simTree")
        # simTree.Print()

        # for iEnt in range(simTree.GetEntries()):
        for iEnt in range(100):
            simTree.GetEntry(iEnt)
            ae = simTree.fAnalysisEvent
            totE = ae.GetTotalEnergy()*1000
            mH = ae.GetNElements()

            if mH < 2: continue

            sumE, sumEA = 0, 0
            for iH in range(mH):
                ele = ae.GetElement(iH)
                hitE = ele.GetEnergy()*1000
                if hitE < 0: continue
                hitID = ele.GetWaveformID()

                sumE += hitE
                act = ele.GetActiveness()
                hitEA = hitE/act
                sumEA += hitEA

                phys = ele.GetPhysicsProcesses()
                dpos = ele.GetLocalPositions()
                np, nd = phys.size(), dpos.size()

                fE = 0
                print("hit",iH)
                for i in range(np):
                    procName, pctEdep = phys[i][0], phys[i][1]
                    fE += pctEdep
                    edep = hitE * pctEdep
                    x, y, z = dpos[i].X(), dpos[i].Y(), dpos[i].Z()
                    print("  %-6s  f %-5.4f  edep %-8.4f  x %-6.3f  y %-6.3f  z %.3f" % (procName, pctEdep, edep, x, y, z))

                print("hit %d:  det %s  hitE %.2f  act %.2f  hitEA %.2f  fE %.4f" % (iH, hitID, hitE, act, hitEA, fE))

            print("Event %d  mH %d  totE %.2f  sumE %.2f  sumEA %.2f\n" % (iEnt, mH, totE, sumE, sumEA))

        tf.Close()


def processTritSim():
    """ Special job for the tritium spectrum """

    outDir = "~/project/sims/trit"
    gatDir = os.environ['GATDIR']
    app = "%s/Apps/process_MJD_as_built_mage_results" % gatDir
    config = "/global/projecta/projectdirs/majorana/users/mbuuck/sim/mageScripts/det_config_for_GAT_pp_by_det_v2.json"

    # best fit to DS3 long cal
    dlThickness = 3.0
    transPoint = 0.95
    transLevel = 0.4

    inFile = "/global/projecta/projectdirs/majorana/users/bxyzhu/MaGe/Raw/Tritium_p001.root"
    cmd = "%s -o %s -c %s -t %.1f %.2f %.1f %s" % (app,outDir,config,dlThickness,transPoint,transLevel,inFile)
    # print(cmd)
    sp.call(shlex.split(cmd))


if __name__=="__main__":
    main()