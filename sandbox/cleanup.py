#!/usr/bin/env python3
import sys, os, glob

skimFails = [
    "skim-ds5-116.txt","skim-ds5-121.txt","skim-ds5-115.txt","skim-ds5-120.txt","skim-ds5-113.txt","skim-ds5-114.txt","skim-ds5-117.txt",
    "skim-ds5-118.txt","skim-ds5-119.txt",
    ]

# delete bad files?

# make a new job sub list

# ./skim_mjd_data 1 11 -l -g -d -t 0.7 /global/projecta/projectdirs/majorana/users/wisecg/bkg/skim >& ./logs/skim-ds1-11.txt &

if __name__=="__main__":

    f = open("../jobLists/test.ls","w")

    for logFile in sorted(skimFails):
        tmp = logFile.split("-")
        dsNum = int(tmp[1][-1])
        subNum = int(tmp[2].split(".")[0])
        dub = "-d" if dsNum < 6 else ""

        # yep, no & at the end.
        # cmd = """./skim_mjd_data %d %d -l -g %s -t 0.7 /global/projecta/projectdirs/majorana/users/wisecg/bkg/skim >& ./logs/skim-ds%d-%d.txt\n""" % (dsNum, subNum, dub, dsNum, subNum)
        # f.write(cmd)

        if os.path.isfile("../logs/"+logFile):
            os.remove("../logs/"+logFile)
            print("deleted",logFile)

    f.close()
