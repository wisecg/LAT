// This code was used to create channel maps for DataSetInfo.py ... I think.
// C. Wiseman

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "GATDataSet.hh"
#include "DataSetInfo.hh"
#include "MJTChannelMap.hh"
#include "MJTChannelSettings.hh"

using namespace std;

void buildMap(int dsNum);
int GetCPD(string pos);
int GetDetID(string det);

int main(int argc, char** argv)
{
  int dsNum = stoi(argv[1]);
  // cout << "dsNum:" << dsNum << endl;
  buildMap(dsNum);
  // string det = "P42574B";
  // string pos = "C1P7";
  // cout << "det:" << GetDetID(det) << endl;
  // cout << "pos:" << GetCP(pos) << endl;
}


void buildMap(int dsNum)
{
  // NOTE:  You could add this to $GATDIR/Apps/check-data.cc ... It's the same structure.

  cout << "Building channel map ...\n";

  GATDataSet ds;
  vector<int> runList;
  // int dsNum = 1;
  for (int rs = 0; rs <= GetDataSetSequences(dsNum); rs++) LoadDataSet(ds, dsNum, rs);
  for (int i=0; i < (int)ds.GetNRuns(); i++) runList.push_back(ds.GetRunNumber(i));

  map<int,int> chanDetID;
  map<int,int> chanCPD;
  vector<int> pMons;

  for (int i = 0; i < (int)runList.size(); i++)
  // for (int i = 0; i < 2; i++)
  {
    int run = runList[i];
    // cout << run << endl;

    GATDataSet ds;
    string gatPath = ds.GetPathToRun(run,GATDataSet::kGatified);

    ifstream chkFile(gatPath.c_str());
    if (!chkFile.good()) {
      cout << "Couldn't find file: " << gatPath << endl;
      continue;
    }
    TFile *gatFile = new TFile(gatPath.c_str());
    MJTChannelMap *chMap = (MJTChannelMap*)gatFile->Get("ChannelMap");
    MJTChannelSettings *chSet = (MJTChannelSettings*)gatFile->Get("ChannelSettings");

    // build pulser monitor list
    vector<uint32_t> enabledPMs = chMap->GetPulserChanList();
    for (auto ch : enabledPMs)
      if (find(pMons.begin(),pMons.end(),ch)==pMons.end())
        pMons.push_back(ch);


    // build channel maps
    vector<uint32_t> enabledIDs = chSet->GetEnabledIDList();
    for (auto ch : enabledIDs)
    {
      // 1. fill detID map
      string det = chMap->GetDetectorName(ch);
      int thisDetID = GetDetID(det);

      // key doesn't exist, add it
      if (!chanDetID.count(ch))
        chanDetID[ch] = thisDetID;
      else {
        // BAD: key exists, value doesn't match previous value (also not a PM)
        if (thisDetID!=chanDetID[ch] && find(pMons.begin(),pMons.end(),ch)==pMons.end()) {
          cout << Form("Error!  Run %i  channel %i  thisDetID %i  prevDetID %i\n  Exiting ...\n",run,ch,thisDetID,chanDetID[ch]);
          return;
        }
      }

      // 2. fill detCPD map
      string pos = chMap->GetDetectorPos(ch);
      int thisCPD = GetCPD(pos);

      // key doesn't exist, add it
      if (!chanCPD.count(ch))
        chanCPD[ch] = thisCPD;
      else {
        // BAD: key exists, value doesn't match previous value (also not a PM)
        if (thisCPD!=chanCPD[ch] && find(pMons.begin(),pMons.end(),ch)==pMons.end()) {
          cout << Form("Error!  Run %i  channel %i  thisCPD %i  prevDetCPD %i\n  Exiting ...\n",run,ch,thisCPD,chanCPD[ch]);
          return;
        }
      }
    }

    gatFile->Close();
  }

  // Now make sure no values are repeated in either map
  // x.first - key, x.second - value
  for (auto const& x : chanDetID) {
    if (find(pMons.begin(), pMons.end(), x.first)!=pMons.end())
      continue;
    for (auto const& y : chanDetID) {
      if (find(pMons.begin(), pMons.end(), x.first)!=pMons.end())
        continue;
      if (y.second==x.second && y.first!=x.first && abs(y.first-x.first)!=1)
          cout << Form("Error!  key %i and key %i both have the same value: %i\n",x.first,y.first,x.second);
    }
  }
  for (auto const& x : chanCPD) {
    if (find(pMons.begin(), pMons.end(), x.first)!=pMons.end())
      continue;
    for (auto const& y : chanCPD) {
      if (find(pMons.begin(), pMons.end(), x.first)!=pMons.end())
        continue;
      if (y.second==x.second && y.first!=x.first && abs(y.first-x.first)!=1)
          cout << Form("Error!  key %i and key %i both have the same value: %i\n",x.first,y.first,x.second);
    }
  }


  // Finally, print the maps in a nice python dict format.

  cout << "\nds" << dsNum << "DetID = {";
  for (auto const &x : chanDetID) cout << x.first << ":" << x.second << ", ";
  cout << "}\n\n";

  cout << "ds" << dsNum << "CPD = {";
  for (auto const &x : chanCPD) cout << x.first << ":" << x.second << ", ";
  cout << "}\n\n";

  cout << "ds" << dsNum << "PM = [";
  for (auto ch : pMons) cout << ch << ", ";
  cout << "]\n\n";
}


int GetCPD(string pos){
  if (pos == "") return 0; // this will break stoi
  pos.erase(remove_if(pos.begin(), pos.end(), [](char c) { return isalpha(c); } ), pos.end());
  return stoi(pos);
}


int GetDetID(string det)
{
  if (det == "") return 0; // this will break stoi
  string tmp = det;
  tmp.erase(remove_if(tmp.begin(), tmp.end(), [](char c) { return isalpha(c); } ), tmp.end());
  int detID = 0;
  if (det.find("P")!=string::npos) {
    detID = stoi(tmp)*10;
    detID += 1000000;
    if (det.find("A")!=string::npos) detID += 0;
    if (det.find("B")!=string::npos) detID += 1;
    if (det.find("C")!=string::npos) detID += 2;
  }
  else if (det.find("B")!=string::npos){
    detID = stoi(tmp);
    detID += 20000;
  }
  // cout << "det " << det << "  tmp " << tmp << "  detID " << detID << endl;
  return detID;
}
