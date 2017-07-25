#!/usr/local/bin/python
from ROOT import TFile, TTree, TEntryList, gDirectory, TNamed

def main():

    gFile = TFile("~/project/mjddatadir/gatified/mjd_run13071.root")
    bFile = TFile("~/project/mjddatadir/built/OR_run13071.root")
    gatTree = gFile.Get("mjdTree")
    bltTree = bFile.Get("MGTree")
    gatTree.AddFriend(bltTree)

    theCut = "mH > 1 && !EventDC1Bits && Entry$ < 320"

    gatTree.Draw(">>elist", theCut, "entrylist")
    elist = gDirectory.Get("elist")
    gatTree.SetEntryList(elist)
    nList = elist.GetN()
    print "Using cut:\n",theCut
    print "Found",gatTree.GetEntries(),"input events."
    print "Found",nList,"events passing cuts."

    print "Starting event loop ..."
    for iList in range(nList):
        entry = gatTree.GetEntryNumber(iList);
        gatTree.LoadTree(entry)
        gatTree.GetEntry(entry)
        nChans = gatTree.channel.size()

        # Loop over hits passing cuts
        numPass = gatTree.Draw("channel",theCut,"GOFF",1,iList)
        chans = gatTree.GetV1()
        chanList = list(set(int(chans[n]) for n in xrange(numPass)))
        hitList = (iH for iH in xrange(nChans) if gatTree.channel.at(iH) in chanList)  # a 'generator expression'
        for iH in hitList:

            # Load data
            run = gatTree.run
            chan = gatTree.channel.at(iH)
            dataENF = gatTree.trapENFCal.at(iH)
            print "run %d  iList %d  chan %d  enf %.2f" % (run, iList, chan, dataENF)

    print "Copying tree with same cut ..."
    outFile = TFile("~/project/lat/copyTree_test.root", "RECREATE")
    outTree = gatTree.CopyTree("")
    outTree.Write()
    print "Wrote",outTree.GetEntries(),"events."
    cutUsed = TNamed("theCut",theCut)
    cutUsed.Write()



if __name__ == "__main__":
    main()
