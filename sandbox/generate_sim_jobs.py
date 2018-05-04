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
    # runAllCalJobs(mod=1, nStart=601, nStop = 999, macroDir=macroDir)

    # Surface Alpha Jobs
    # Surface Sources: DUPTFE, DUCopper
    genSurfaceMacros(nStart=1, nStop=2, mod=1, macroDir=macroDir, outDir=outDir)
    # genAlphaMacros(nStart=1, nStop=999, mod=1, macroDir=macroDir, outDir=outDir, nEvents = 500000)

def genCalMacros(mod=1, nStart=1, nStop=100, oldAF = False, macroDir='.', outDir='/mnt/mjdDisk1/Majorana/users/psz/MAGE/Simulations/5M', nEvents = 500000):
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


def genSurfaceMacros(sourceName='DUPTFE', mod=1, nStart=1, nStop=100, macroDir='.', outDir='/mnt/mjdDisk1/Majorana/users/psz/MAGE/Simulations/5M', nEvents = 5000000):
    for nFile in range(nStart, nStop+1):
        f = open('{}/M{}Alpha_A210_Z84_gss_p{:03d}.mac'.format(macroDir, mod, nFile), 'w+')

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
        f.write("/MG/eventaction/rootfilename {}/M{}Alpha_A210_Z84_gss_p{:03d}.root\n".format(outDir, mod, nFile))

        # Add Volumes here
        for vol in volList:
            stripped = vol.rstrip('_0123456789')
            if sourceName == 'DUCopper':
                if stripped.endswith('HVRing77'): f.write("/MG/io/gss/addVolume {}\n".format(stripped))
                if stripped.endswith('HollowHexRod'): f.write("/MG/io/gss/addVolume {}\n".format(stripped))
                if stripped.endswith('FlexInsulator'): f.write("/MG/io/gss/addVolume {}\n".format(stripped))
                if stripped.endswith('CrystalMountingPlate'): f.write("/MG/io/gss/addVolume {}\n".format(stripped))
                if stripped.endswith('ContactPin'): f.write("/MG/io/gss/addVolume {}\n".format(stripped))
                if stripped.endswith('SpringFEMount'): f.write("/MG/io/gss/addVolume {}\n".format(stripped))
                if stripped.endswith('LMFECoverPlate'): f.write("/MG/io/gss/addVolume {}\n".format(stripped))
                if stripped.endswith('SpringNut'): f.write("/MG/io/gss/addVolume {}\n".format(stripped))
            elif sourceName == 'DUPTFE':
                if stripped.endswith('HVNut'): f.write("/MG/io/gss/addVolume {}\n".format(stripped))
                if stripped.endswith('CrystalInsulator'): f.write("/MG/io/gss/addVolume {}\n".format(stripped))
                if stripped.endswith('CenterBushing'): f.write("/MG/io/gss/addVolume {}\n".format(stripped))

        f.write("/MG/generator/gss/boundvol RadShieldAssembly_001_RadShieldCuInner_001\n")
        if sourceName == 'DUCopper':
            f.write("/MG/io/gss/setMaxIntersections 26\n")
        elif sourceName == 'DUPTFE':
            f.write("/MG/io/gss/setMaxIntersections 20\n")
        f.write("/run/beamOn {}".format(nEvents))
        f.close()


def runAllCalJobs(mod=1, nStart=0, nStop=None, macroDir='.'):
    fileList = sorted(glob.glob('{}/M{}CalSource_*'.format(macroDir, mod)))
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
