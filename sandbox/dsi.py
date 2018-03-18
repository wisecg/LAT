#!/usr/bin/env python3
import os

latSWDir    = os.environ['LATDIR']
dataDir     = "/global/projecta/projectdirs/majorana/users/wisecg"
bkgDir      = dataDir+"/bkg"     # subfolders: skim waves split lat
calDir      = dataDir+"/cal"     # subfolders: skim waves split lat
specialDir  = dataDir+"/special" # subfolders: skim waves split lat
skimDir     = bkgDir+"/skim"
waveDir     = bkgDir+"/waves"
splitDir    = bkgDir+"/split"
latDir      = bkgDir+"/lat"
calSkimDir  = calDir+"/skim"
calWaveDir  = calDir+"/waves"
calSplitDir = calDir+"/split"
calLatDir   = calDir+"/lat"
pandaDir    = dataDir+"/pandas"

# class bkgInfo:

def parseDSI():
    print("heyy")


if __name__=="__main__":
    parseDSI()