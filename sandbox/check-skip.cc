// check-skip.cc
// C. Wiseman, USC

#include <cstdio>
#include <fstream>
#include "TFile.h"
#include "TTreeReader.h"
#include "GATDataSet.hh"

using namespace std;

int GetDataSetSequences(int dsNum);

int main()
{
  // load every skim file in this directory
  vector<string> skims;
  for (int i=0; i < 6; i++) {
    for (int j=0; j <= GetDataSetSequences(i); j++) {
      string fileName = Form("skimDS%i_%i.root",i,j);
      ifstream f(fileName.c_str());
      if (!f.good()) cout << "File not found! " << fileName << endl;
      skims.push_back(fileName);
    }
  }

  // for every tree, check the channel vector
  for (auto fileName: skims)
  {
    cout << "Now scanning " << fileName << " ..";
    TFile f(fileName.c_str());
    TTreeReader reader("skimTree", &f);
    TTreeReaderValue< vector<int> > channelIn(reader, "channel");

    while (reader.Next())
    {
      vector<int> chan = (*channelIn);
      for (auto ch : chan) {
        if (ch%2 == 1) continue;
        if ( find(chan.begin(), chan.end(), ch+1) != chan.end() ) {
          cout << Form("Shit.  Event %lli has channel %i and channel %i.\n", reader.GetCurrentEntry(), ch, ch+1);
        }
      }
    }
    cout << ". done.\n";
  }
}

// When LoadDataSet is updated, this must be updated as well.
int GetDataSetSequences(int dsNum)
{
  map<int,int> dsMap = {{0,75},{1,51},{2,7},{3,24},{4,18},{5,112},{6,5}};
  return dsMap[dsNum];
}
