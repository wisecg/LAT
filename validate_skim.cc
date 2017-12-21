// validate_skim.cc
// I. Guinn, UW

#include <cstdio>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <utility>

#include "GATDataSet.hh"
#include "DataSetInfo.hh"

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TTree.h"
#include "TChain.h"
#include "TLeaf.h"
#include "TTreeFormula.h"

using namespace std;

int main(int argc, char** argv) {
  if(argc!=2) {
    cout << "Expected usage: validate_skim path/to/skimfile.root" << endl;
    return 1;
  }

  string skimname(argv[1]);
  TFile skimfile(skimname.c_str());
  TTree* skimtree = dynamic_cast<TTree*>(skimfile.Get("skimTree"));
  if(skimtree==NULL) {
    cout << "Could not find skimTree" << endl;
    return 1;
  }

  TTreeReader skim(skimtree);
  char isrun[10], issmall[10];
  int dsnum=0, subdsnum=0;

  //load the GATDataSet based on the file name
  size_t strpos = skimname.rfind('/');
  if(strpos==string::npos) strpos=0;
  else strpos++;
  // sscanf(skimname.c_str()+strpos, "skimDS%d%[_ a-z A-Z]%d%[_ a-z A-Z].root", &dsnum, isrun, &subdsnum, issmall); // original
  sscanf(skimname.c_str()+strpos, "skimDS%d%[_ a-z A-Z]%d%[_ a-z A-Z]_low.root", &dsnum, isrun, &subdsnum, issmall); // low-e skim
  // sscanf(skimname.c_str()+strpos, "skimDS%d%[_ a-z A-Z]%d%[_ a-z A-Z]_low.root", &dsnum, isrun, &subdsnum, issmall); // waveSkim

  GATDataSet ds;
  if(strcmp(isrun, "_run")==0) ds.AddRunNumber(subdsnum);
  else LoadDataSet(ds, dsnum, subdsnum);

  TChain* gatch = ds.GetChains();
  TTreeReader gat(gatch);

  //set up branches
  TTreeReaderValue<int> runskim(skim, "run");
  TTreeReaderValue<int> ievskim(skim, "iEvent");
  TTreeReaderArray<int> ihit(skim, "iHit");
  TTreeReaderValue<double> rungat(gat, "run");

  // first string is name of branch in skim tree. vector is branches in gatified and/or built trees to compare to
  vector< pair<string, vector<string> > > leafs2compare;
  leafs2compare.emplace_back("run", vector<string>{"fRunNumber"});
  leafs2compare.emplace_back("iEvent", vector<string>{"LocalEntry$"});
  leafs2compare.emplace_back("channel", vector<string>{"channel", "fDigitizerData.GetID()", "fWaveforms.fID"});
  //leafs2compare.emplace_back("index", vector<string>{"index", "fIndex"});
  leafs2compare.emplace_back("clockTime_s", vector<string>{"clockTime/1e9", "fTime/1e9"});
  leafs2compare.emplace_back("tOffset", vector<string>{"tOffset", "fTOffset"});
  //leafs2compare.emplace_back("TMath::Nint(tOffset/10 + clockTime_s*1e8)", vector<string>{"fTimeStamp"});


  vector< pair<TTreeFormula*, vector<TTreeFormula*> > > forms;
  for(auto& formname : leafs2compare) {
    vector<TTreeFormula*> formsgat;
    for(auto& formname2 : formname.second)
    formsgat.emplace_back(new TTreeFormula((formname2+"gat").c_str(), formname2.c_str(), gatch));
    forms.push_back(make_pair(new TTreeFormula((formname.first+"skim").c_str(), formname.first.c_str(), skimtree), formsgat));
  }

  Long64_t treestart = 0;
  gat.SetEntry(0);

  for(auto entry : skim) {
    //find the gat tree entry that matches the current skim tree entry
    while(*runskim != *rungat) {
      treestart+=gat.GetTree()->GetTree()->GetEntries();
      TTreeReader::EEntryStatus st = gat.SetEntry(treestart);
      if(st != 0) {
        cout << "Warning: could not open tree for run " << *runskim << " read status is " << st << " treestart " << treestart << "/" << gat.GetEntries(1) << endl;
        return 0;
      }

      for(auto& formlist : forms)
      for(auto& form : formlist.second)
      form->UpdateFormulaLeaves();
    }

    if(gat.SetEntry(treestart + *ievskim) != 0) {
      cout << "Warning: could not find entry " << *ievskim << " in run " << *rungat << endl;
      return 0;
    }

    for(auto& formlist : forms) {
      Int_t niterskim = formlist.first->GetNdata();

      for(auto& form : formlist.second) {
        Int_t nitergat = form->GetNdata();
        if(nitergat==1) {
          if(form->EvalInstance(0) != formlist.first->EvalInstance(0)) {
            cout << "Warning: skim leaf " << formlist.first->GetTitle()
            << " in entry " << entry << " = "
            << formlist.first->EvalInstance(0) << " does not equal gat leaf "
            << form->GetTitle() << " in entry " << *ievskim
            << " = " << form->EvalInstance(0) << endl;
            //return 0;
          }
        }
        else {
          if(UInt_t(niterskim) != ihit.GetSize()) {
            cout << "Warning: Run " << *runskim <<  " skim leaf " << formlist.first->GetLeaf(0)->GetName() << " does not have one entry for each hit!" << endl;
            return 0;
          }
          else {
            for(size_t i=0; i<ihit.GetSize(); i++) {
              if(formlist.first->EvalInstance(i) != form->EvalInstance(ihit[i])) {
                cout << "Warning: Run " << *runskim << " skim leaf " << formlist.first->GetTitle()
                << " in entry " << entry << " = "
                << formlist.first->EvalInstance(i) << " does not equal gat leaf "
                << form->GetTitle() << " in entry " << *ievskim
                << " = " << form->EvalInstance(ihit[i]) << endl;
                //return 0;
              }
            }
          }
        }
      }
    }
  }
}

