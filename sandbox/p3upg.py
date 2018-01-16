#!/usr/bin/env python3
import sys, imp
ds = imp.load_source('DataSetInfo','../DataSetInfo.py')

def main():
    print(sys.version)
    print(ds.GetGoodChanList(5))
    # raw_input()
    input()



if __name__=="__main__":
    main()