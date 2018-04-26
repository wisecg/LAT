#!/usr/bin/env python3
import glob

def main():

    fileList = sorted(glob.glob("../logs/lat*.txt"))
    for f in fileList:
        print(f)

        with open(f) as ftmp:
            fl = ftmp.readlines()

        iLine = -1
        for i, l in enumerate(fl):
            ltmp = l.rstrip()

            if "Detector Thresholds" in ltmp:
                iLine = i
                break

        for l in fl[iLine:]:
            ltmp = l.rstrip()
            print(ltmp)

        # return

if __name__=="__main__":
    main()