#!/usr/bin/env python
import sys, imp, time, glob
sys.argv.append("-b") # kill all interactive crap
wl = imp.load_source('waveLibs','../waveLibs.py')

def main():
    print "hi!"

if __name__=="__main__":
    main()