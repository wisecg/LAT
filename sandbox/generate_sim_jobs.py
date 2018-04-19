#!/usr/bin/env python
"""  Short script for generating and submitting a buncha MaGe macros for calibration """
import sys, shlex, glob, os
import subprocess as sp

jobsubStr = "sbatch slurm-job_pdsf.sh"

def main():
    # PDSF
    # macroDir = '/projecta/projectdirs/majorana/users/bxyzhu/MaGe/macros/batch'
    # outDir = '/projecta/projectdirs/majorana/users/bxyzhu/MaGe/Raw/batch'

    # LANL
    macroDir = '/mnt/mjdDisk1/Majorana/users/psz/MAGE/macros/batch'
    outDir = '/mnt/mjdDisk1/Majorana/users/psz/MAGE/Simulations/5M'

    # genCalMacros(150, mod=1, macroDir=macroDir, outDir=outDir)
    runAllCalJobs(mod=2, macroDir=macroDir)

def genCalMacros(nFiles, mod=1, oldAF = True, macroDir='.', outDir='/mnt/mjdDisk1/Majorana/users/psz/MAGE/Simulations/5M', nEvents = 500000):
    for nFile in range(1, nFiles+1):
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
