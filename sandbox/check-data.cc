// check-data.cc
// Clint Wiseman, USC/MJD
// May 4, 2017
//
// Intended to assist database-based run selection
// and file consistency checks.  Takes a run list from
// DataSetInfo, the runDB, or (future) some text file,
// and checks its suitability for inclusion in an official DS.
//
// DB access requires MkCookie be run recently.

#include <iostream>
#include <bitset>
#include "TFile.h"
#include "TTreeReader.h" // so I can use Form()
#include "TEntryList.h"
#include "TTree.h"
#include "MJTRun.hh"
#include "MJDatabase.hh"
#include "MJDocument.hh"
#include "GATDataSet.hh"
#include "DataSetInfo.hh"

using namespace std;
using namespace MJDB;

void checkData(int dsNum, bool loadTrees, vector<int> list);
vector<int> getDBRunList(int &dsNum, string partNum, string runRank, string dbKey, bool docs);

int main(int argc, char** argv)
{
  // Get user args
  if (argc < 2) {
    cout << "Usage: ./check-data [options]\n"
         << "  Options: \n"
         << "    -b [dsNum] : use a whole DS from DataSetInfo.hh\n"
         << "    -d [partNum] [runRank] : get list from the run database\n"
         << "    -w [whole key] : use a custom key for the runDB\n"
         << "    -i : set to include_docs=true\n"
         << "    -l : load chains and do checks (slow)\n"
         << "    -x : access DB only \n"
         << "  RunDB access examples (-db option):\n"
         << "    partNum = P3LQK, P3KJR, P3LQG, etc.\n"
         << "    runRank = gold, silver, bronze, cal, etc.\n";
    return 1;
  }
  int dsNum=-1;
  bool loadTrees=0, useDB=0, docs=0, onlyDB=0;
  string partNum="",runRank="",dbKey="";
  vector<string> opt(argv+1, argv+argc);
  for (size_t i = 0; i < opt.size(); i++) {
    if (opt[i] == "-b") { dsNum = stoi(opt[i+1]); }
    if (opt[i] == "-d") { useDB=1; partNum=opt[i+1]; runRank=opt[i+2]; }
    if (opt[i] == "-w") { dbKey = opt[i+1]; }
    if (opt[i] == "-i") { docs=1; }
    if (opt[i] == "-l") { loadTrees=1; }
    if (opt[i] == "-x") { onlyDB=1; }
    // I could also add an option to read in a run list from a file ...
    // Also, a run range limit for the database access ...
  }

  // Run the checker routines
  vector<int> runList;
  if (useDB) runList = getDBRunList(dsNum,partNum,runRank,dbKey,docs);
  if (!onlyDB) checkData(dsNum,loadTrees,runList);
}


void checkData(int dsNum, bool loadTrees, vector<int> list)
{
  // Start with a run list
  vector<int> runList;
  if (list.size()==0) {
    cout << "Scanning DS-" << dsNum << endl;
    GATDataSet ds;
    for (int rs = 0; rs < GetDataSetSequences(dsNum); rs++) LoadDataSet(ds, dsNum, rs);
    for (int i=0; i<(int)ds.GetNRuns(); i++) runList.push_back(ds.GetRunNumber(i));
  }
  else {
    cout << "Using run list from DB ...\n";
    runList = list;
  }

  // Loop over runs, accessing individual files
  for (auto run : runList)
  {
    GATDataSet *dsRun = new GATDataSet(run);

    // Check that built data exists
    string bltPath = dsRun->GetPathToRun(run,GATDataSet::kBuilt);
    TFile *bltFile = new TFile(bltPath.c_str());
    if (bltFile->IsZombie()) {
      cout << Form("Run %i has no built file. Skipping ...\n",run);
      continue;
    }

    // Check BBDecay run bit
    MJTRun *runInfo = (MJTRun*)bltFile->Get("run");
    bool bbDecay=0;
    bitset<32> event_type = runInfo->GetRunBits();
    if (event_type.test(1)) bbDecay=1;

    // Check that gat data exists and it's not a blinded (inaccessible) run
    bool deadGat=0;
    string gatPath = dsRun->GetPathToRun(run,GATDataSet::kGatified);
    TFile *gatFile = new TFile(gatPath.c_str());
    if (gatFile->IsZombie()) deadGat=1;
    if (bbDecay && deadGat) cout << Form("Run %i is a blind run, should not be included!\n",run);
    if (deadGat) {
      delete bltFile; delete gatFile;
      continue;
    }

    // Check for startClockTime branch (indicates 4/2017 reprocessing was successful)
    TTree *gat = (TTree*)gatFile->Get("mjdTree");
    TBranch *b1 = gat->FindBranch("startClockTime");
    if (b1==NULL) {
      cout << Form("Run %i doesn't have startClockTime\n",run);
      continue;
    }
    TBranch *b2 = gat->FindBranch("EventDC1Bits");
    if (b2==NULL) {
      cout << Form("Run %i doesn't have EventDC1Bits\n",run);
      continue;
    }

    if (loadTrees)
    {
      // Check if veto data exists (ignore DS4)
      string vetPath;
      TFile *vetFile = NULL;
      bool noVeto=0;
      if (dsNum != 4) {
        vetPath = dsRun->GetPathToRun(run,GATDataSet::kVeto);
        vetFile = new TFile(vetPath.c_str());
      }
      else noVeto=1;
      if (!noVeto && vetFile->IsZombie()) {
        cout << Form("Run %i has no veto data.\n",run);
        noVeto=1;
      }

      // Check we have nonzero entries in all 3 chains
      long nGat = gat->GetEntries();
      TTree *blt = (TTree*)bltFile->Get("MGTree");
      long nBlt = blt->GetEntries();
      long nVet = -1;
      TTree *vet = NULL;
      if (!noVeto) {
        vet = (TTree*)vetFile->Get("vetoTree");
        vet->GetEntries();
      }
      if (nGat==0 || nBlt==0 || nVet==0)
        cout << Form("Run %i, zero entries in chain:  gat %li  built %li  veto %li\n",run,nGat,nBlt,nVet);

      // Check for orphan LG events above 500 keV.
      // Did some tricks to make the loop fast, but it's still pretty slow.
      vector<int> *channel=0;
      vector<double> *trapENFCal=0;
      gat->SetBranchStatus("*",0);
      gat->SetBranchStatus("channel",1);
      gat->SetBranchStatus("trapENFCal",1);
      gat->SetBranchStatus("P",1);
      gat->SetBranchStatus("D",1);
      gat->SetBranchStatus("EventDC1Bits",1);
      gat->SetBranchAddress("channel",&channel);
      gat->SetBranchAddress("trapENFCal",&trapENFCal);
      gat->Draw(">>elist","trapENFCal > 500 && trapENFCal < 15000 && P!=0 && D!=0 && !EventDC1Bits","entrylist GOFF");
      TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
      gat->SetEntryList(elist);
      for (int iList = 0; iList < elist->GetN(); iList++) {
        int iEvent = elist->GetEntry(iList);
        gat->GetEntry(iEvent);
        // cout << Form("Run %i  iL %i  iE %i  trapENF %.0f\n",run,iList,iEvent,trapENFCal->at(0));
        int nHits = channel->size();
        vector<int> hits;
        map<int,int> gainMap;
        for (int iHit = 0; iHit < nHits; iHit++) gainMap[channel->at(iHit)] = (int)iHit;
        for (int iHit = 0; iHit < nHits; iHit++) {
          int chan = channel->at(iHit);
          if (chan%2==1 && !(gainMap.find(chan-1)!=gainMap.end())) {
            cout << Form("Run %i  Orphan LG Event.  Gat Entry %i  Index %i  Chan %i  trapENFCal %.1f\n", run,iEvent,iHit,chan,trapENFCal->at(iHit));
          }
        }
      }
      cout << Form("Run %i  elist %lli\n",run,elist->GetN()); // Just lets you know the app isn't dead
      delete bltFile;
      if (!noVeto) delete vetFile;
    }
    delete gatFile;
    delete dsRun;
  }
}


vector<int> getDBRunList(int &dsNum, string partNum, string runRank, string dbKey, bool docs)
{
  string view = Form("run_rank?key=[\"%s\",\"%s\"]",partNum.c_str(),runRank.c_str());
  if (dbKey!="") view = dbKey;

  // Access DB.
  const string dbString = "mjd_run_database";
  const string dbServer = "mjdb.phy.ornl.gov";
  MJDatabase runDB(&dbString, &dbServer);
  runDB.SetServerScheme("https");
  MJDocument runDoc;
  runDoc.Get_View(runDB,"dbApp",view,docs);  // Last arg is include_docs.  If true, downloads ALL runDB fields.
  string errorMessage;
  if (runDB.GetDBStatus(errorMessage)!=0){
    cout << "Failed to get document.  cURL error: " << runDB.GetDBStatus(errorMessage)
         << " Message: " << errorMessage << endl;
    return vector<int>();
  }
  int nDocs = runDoc["rows"].Length();
  cout << "Found " << nDocs << " run records.\n";
  cout << runDB.GetURL() << endl; // you can check this in a browser
  // runDoc.Dump();

  // Loop over the document.
  vector<int> runList;
  for (int i = 0; i < nDocs; i++) {

    int run = atoi(runDoc["rows"][i]["value"].Value().AsString().c_str());
    runList.push_back(run);

    // Some examples of other things you can do:
    // If you set include_docs=true above, you can pretty much slice and dice the
    // records from the runDB as much as you want -- by runBits (disruptive work, expert mode ...),
    // gat status, if it's calibration or background, etc.

    // runDoc["rows"][i].Dump(); // dump just one document
    // int runInDoc = atoi(runDoc["rows"][i]["doc"]["RunNumber"].Value().AsString().c_str());
    // bool isOpen = ( runDoc["rows"][i]["doc"]["access"].Value().AsString() == "open" );
    // int runBits = atoi( runDoc["rows"][i]["doc"]["orca_run_bits"].Value().AsString().c_str() );
    // int runQuality = atoi( runDoc["rows"][i]["doc"]["RunQualityVal"].Value().AsString().c_str() );
    // string runRank = runDoc["rows"][i]["doc"]["RunRank"].Value().AsString();  // gold, silver, etc
    string timestamp = runDoc["rows"][i]["doc"]["timestamp"].Value().AsString();
    // double runDuration = stod( runDoc["rows"][i]["doc"]["ElapsedTime"].Value().AsString() );

    cout << Form("Run %i  timestamp %s\n",run,timestamp.c_str());
  }

  // Figure out the dataset
  dsNum = FindDataSet(runList[0]);

  return runList;
}
