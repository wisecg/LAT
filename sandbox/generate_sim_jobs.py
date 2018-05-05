#!/usr/bin/env python
"""  Short script for generating and submitting a buncha MaGe macros for calibration """
import sys, shlex, glob, os, json
import subprocess as sp

jobsubStr = "sbatch slurm-job_pdsf.sh"
# This file is copied from the MaGe output directories
vol_data = json.load(open('{}/data/uniformGssVolumeList.json'.format(os.getenv('LATDIR'))))
volList = vol_data['volumes']

def main():
    # PDSF
    macroDir = '/projecta/projectdirs/majorana/users/bxyzhu/MaGe/macros/surf'
    outDir = '/projecta/projectdirs/majorana/users/bxyzhu/MaGe/Raw/surf'

    # LANL
    # macroDir = '/mnt/mjdDisk1/Majorana/users/psz/MAGE/macros/batch'
    #outDir = '/mnt/mjdDisk1/Majorana/users/psz/MAGE/Simulations/5M'

    # Calibration Jobs
    # genCalMacros(nStart=601, nStop=999, mod=1, macroDir=macroDir, outDir=outDir, nEvents = 500000)

    # Surface Alpha Jobs
    # Surface Sources: DUPTFE, DUCopper
    genSurfaceMacros(sourceName='DUCopper', nStart=1, nStop=250, macroDir=macroDir, outDir=outDir)
    #genAlphaMacros(sourceName='DUPTFE', nStart=1, nStop=1, macroDir=macroDir, outDir=outDir)

    #runAllJobs(sourceName='DUPTFE', mod=1, nStart=1, nStop = 250, macroType='gss', macroDir=macroDir)

def genCalMacros(mod=1, nStart=1, nStop=100, oldAF = False, macroDir='.', outDir='.', nEvents = 500000):
    for nFile in range(nStart, nStop+1):
        f = open('{}/M{}CalSource_A224_Z88_p{:03d}.mac'.format(macroDir, mod, nFile), 'w+')

        f.write("/MG/manager/mglog routine\n")
        f.write("/MG/eventaction/reportingfrequency 10000\n")
        f.write("/MG/manager/seedWithUUID\n")
        f.write("/MG/processes/realm BBdecay\n")
        f.write("/MG/processes/lowenergy true\n")
        f.write("/MG/processes/useAllHP true\n")
        f.write("/MG/geometry/detector MJDemonstrator\n")
        f.write("/MG/geometry/WorldMaterial Vacuum\n")
        f.write("/MG/eventaction/rootschema MCRun\n")
        f.write("/MG/eventaction/rootfilename {}/M{}CalSource_A224_Z88_p{:03d}.root\n".format(outDir, mod, nFile))
        f.write("/MG/io/MCRun/SetSensitiveIDLabelScheme askGeom\n")
        f.write("/MG/io/MCRun/setRunID {:03d}\n".format(nFile))
        f.write("/MG/io/MCRun/useTimeWindow true\n")
        f.write("/MG/io/MCRun/setTimeWindow 86400 second\n")
        f.write("/MG/demonstrator/muonVetoOn false\n")
        if oldAF:
            f.write("/MG/demonstrator/innerCopperOn false\n")
        else:
            f.write("/MG/demonstrator/innerCopperOn true\n")
        if mod==1:
            f.write("/MG/mjdemocalassemblyW/sourceOn 1\n")
        f.write("/run/initialize\n")
        f.write("/MG/generator/select MJDCalibration\n")
        if mod==2:
            f.write("/MG/generator/MJDCalibration/setSourcePos E\n")
        f.write("/MG/generator/MJDCalibration/setA 224\n")
        f.write("/MG/generator/MJDCalibration/setZ 88\n")
        f.write("/grdm/nucleusLimits 208 224 81 88\n")
        f.write("/run/beamOn {}".format(nEvents))
        f.close()


def genSurfaceMacros(sourceName='DUPTFE', nStart=1, nStop=100, macroDir='.', outDir='.', nEvents = 5000000):
    for nFile in range(nStart, nStop+1):
        f = open('{}/MJDem_{}_A210_Z84_gss_p{:03d}.mac'.format(macroDir, sourceName, nFile), 'w+')

        f.write("/MG/manager/mglog routine\n")
        f.write("/MG/eventaction/reportingfrequency 100000\n")
        f.write("/MG/manager/seedWithUUID\n")
        f.write("/MG/processes/realm BBdecay\n")
        f.write("/MG/processes/lowenergy true\n")
        f.write("/MG/processes/useAllHP true\n")
        f.write("/MG/geometry/detector MJDemonstrator\n")
        f.write("/MG/geometry/WorldMaterial Vacuum\n")
        f.write("/MG/demonstrator/muonVetoOn false\n")
        f.write("/MG/demonstrator/innerCopperOn true\n")
        f.write("/MG/processes/useNoHadPhysics true\n")
        f.write("/run/initialize\n")
        f.write("/MG/generator/select GSS\n")
        f.write("/MG/eventaction/rootschema GSS\n")
        f.write("/MG/eventaction/rootfilename {}/MJDem_{}_A210_Z84_gss_p{:03d}.root\n".format(outDir, sourceName, nFile))

        # Add Volumes here
        for vol in volList:
            stripped = vol.rstrip('_0123456789')
            if sourceName == 'DUCopper':
                if stripped.endswith('HVRing77'): f.write("/MG/io/gss/addVolume {}\n".format(vol))
                if stripped.endswith('HollowHexRod'): f.write("/MG/io/gss/addVolume {}\n".format(vol))
                if stripped.endswith('FlexInsulator'): f.write("/MG/io/gss/addVolume {}\n".format(vol))
                if stripped.endswith('CrystalMountingPlate'): f.write("/MG/io/gss/addVolume {}\n".format(vol))
                if stripped.endswith('ContactPin'): f.write("/MG/io/gss/addVolume {}\n".format(vol))
                if stripped.endswith('SpringFEMount'): f.write("/MG/io/gss/addVolume {}\n".format(vol))
                if stripped.endswith('LMFECoverPlate'): f.write("/MG/io/gss/addVolume {}\n".format(vol))
                if stripped.endswith('SpringNut'): f.write("/MG/io/gss/addVolume {}\n".format(vol))
            elif sourceName == 'DUPTFE':
                if stripped.endswith('HVNut'): f.write("/MG/io/gss/addVolume {}\n".format(vol))
                if stripped.endswith('CrystalInsulator'): f.write("/MG/io/gss/addVolume {}\n".format(vol))
                if stripped.endswith('CenterBushing'): f.write("/MG/io/gss/addVolume {}\n".format(vol))

        f.write("/MG/generator/gss/boundvol RadShieldAssembly_001_RadShieldCuInner_001\n")
        if sourceName == 'DUCopper':
            f.write("/MG/io/gss/setMaxIntersections 26\n")
        elif sourceName == 'DUPTFE':
            f.write("/MG/io/gss/setMaxIntersections 28\n")
        f.write("/run/beamOn {}".format(nEvents))
        f.close()


def genAlphaMacros(sourceName='DUPTFE', nStart=1, nStop=100, macroDir='.', outDir='.'):
    # This assumes that the output directory of the gss files is the same as the output directory of the surf files
    # Requires pyROOT!
    import ROOT
    for nFile in range(nStart, nStop+1):
        nEvents = 0
        f = open('{}/MJDem_{}_A210_Z84_surf_p{:03d}.mac'.format(macroDir, sourceName, nFile), 'w+')
        f.write("/MG/manager/mglog routine\n")
        f.write("/MG/eventaction/reportingfrequency 10000\n")
        f.write("/MG/manager/seedWithUUID\n")
        f.write("/MG/processes/realm BBdecay\n")
        f.write("/MG/processes/lowenergy true\n")
        f.write("/MG/processes/useAllHP true\n")
        f.write("/MG/geometry/detector MJDemonstrator\n")
        f.write("/MG/geometry/WorldMaterial Vacuum\n")
        f.write("/MG/eventaction/rootschema MCRun\n")
        f.write("/MG/eventaction/rootfilename {}/MJDem_{}_A210_Z84_surf_p{:03d}.root\n".format(outDir, sourceName, nFile))
        f.write("/MG/io/MCRun/SetSensitiveIDLabelScheme askGeom\n")
        f.write("/MG/io/MCRun/setRunID {:03d}\n".format(nFile))
        f.write("/MG/io/MCRun/useTimeWindow true\n")
        f.write("/MG/io/MCRun/setTimeWindow 86400 second\n")
        f.write("/MG/demonstrator/muonVetoOn false\n")
        f.write("/MG/demonstrator/innerCopperOn true\n")
        f.write("/run/initialize\n")
        f.write("/MG/generator/select RDMiso\n")
        f.write("/gun/energy 0 eV\n")
        f.write("/grdm/ion 210 84 0\n")
        f.write("/grdm/nucleusLimits 206 210 80 84\n")

        # Here we have to find the gss ROOT file and get the number of entries to set for the number of events
        gssFileName = '{}/MJDem_{}_A210_Z84_gss_p{:03d}.root'.format(outDir, sourceName, nFile)
        gssFile = ROOT.TFile(gssFileName)
        tree = gssFile.Get('GSSTree')
        nEvents = tree.GetEntries()
        gssFile.Close()

        f.write("/MG/generator/gsspositionsfile {}\n".format(gssFileName))
        f.write("/run/beamOn {}".format(nEvents))
        f.close()

def runAllJobs(sourceName='Cal', macroType='gss', mod=1, nStart=0, nStop=None, macroDir='.'):
    # Macro type is either gss or surf
    fileList = []
    if sourceName == 'Cal':
        fileList = sorted(glob.glob('{}/M{}CalSource_*.mac'.format(macroDir, mod)))
    else:
        fileList = sorted(glob.glob('{}/MJDem_{}_A210_Z84_{}_*.mac'.format(macroDir, sourceName, macroType)))

    if nStop > len(fileList):
        print("Stop value too large")
        return
    if nStop == None:
        nStop = len(fileList)
    for macroFile in fileList[nStart-1:nStop]:
        print("""{} 'MaGe {}'""".format(jobsubStr, macroFile))
        sh("""{} 'MaGe {}'""".format(jobsubStr, macroFile))


def sh(cmd):
    sp.call(shlex.split(cmd))
    return


if __name__ == "__main__":
    main()
