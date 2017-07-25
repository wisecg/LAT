// Identify events in skim data passing
// a given TCut, and append the corresponding waveforms
// from built data into a new file.
// C. Wiseman, 1/18/2017
//         v.2 3/06/2017
//         v.3 5/18/2017 (thanks to Ian w/ waveform stack)
//         v.4 7/25/2017 - added multisampling support

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stack>
#include <array>
#include "TFile.h"
#include "TROOT.h"
#include "TChain.h"
#include "TEntryList.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "GATDataSet.hh"
#include "MGTWaveform.hh"
#include "MJTMSWaveform.hh"
#include "MJTRun.hh"

using namespace std;

void SkimWaveforms(string theCut, string inFile, string outFile);
void TCutSkimmer(string theCut, int dsNum);
void diagnostic();

int main(int argc, char** argv)
{
  // Get some (m)args
  if (argc < 2){
    cout << "Usage: ./wave-skim [options]\n"
         << "       [-s : run TCutSkimmer]\n"
         << "       [-r [dsNum] [subNum] : specify DS and sub-DS]\n"
         << "       [-f [dsNum] [runNum] : specify DS and run num]\n"
         << "       [-p [inPath] [outPath]: file locations]\n"
         << "       [-c : use calibration TCut]\n";
    return 1;
  }
  string inPath="", outPath="";
  int dsNum=-1, subNum=0, run=0;
  bool sw=0, tcs=0, fil=0, cal=0;
  vector<string> opt(argv, argv+argc);
  for (size_t i = 0; i < opt.size(); i++) {
    if (opt[i] == "-s") { sw=0; tcs=1; }
    if (opt[i] == "-f") { sw=1; fil=1; dsNum = stoi(opt[i+1]); run = stoi(opt[i+2]); }
    if (opt[i] == "-r") { sw=1; dsNum = stoi(opt[i+1]); subNum = stoi(opt[i+2]); }
    if (opt[i] == "-p") { inPath = opt[i+1]; outPath = opt[i+2]; }
    if (opt[i] == "-c") { cal=1; }
  }

  // Set cut
  string theCut = "trapENFCal>0.8 && gain==0 && mH==1 && isGood && !muVeto && !isLNFill1 && !isLNFill2 && P!=0 && D!=0 && trapETailMin<0.5"; // DS0-5 standard cut

  // test cut
  // theCut = "trapENFCal>20 && gain==0 && mH==1 && isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0 && P!=0 && D!=0";
  theCut = "trapENFCal > 20 && trapENFCal < 100 && mH==1";

  if (cal) theCut = "trapENFCal>0.5 && trapENFCal<250 && gain==0 && isGood && !muVeto && !isLNFill1 && !isLNFill2 && P!=0 && D!=0 && trapETailMin<0.5"; // calibration file cut

  if (dsNum == 5) theCut += " && channel!=692 && channel!=1232";

  // Set file I/O
  string inFile = Form("%s/skimDS%i_%i_low.root",inPath.c_str(),dsNum,subNum);
  string outFile = Form("%s/waveSkimDS%i_%i.root",outPath.c_str(),dsNum,subNum);
  if (fil) {
    inFile = Form("%s/skimDS%i_run%i_low.root",inPath.c_str(),dsNum,run);
    outFile = Form("%s/waveSkimDS%i_run%i.root",outPath.c_str(),dsNum,run);
  }

  // Go running
  // diagnostic();
  cout << "Scanning DS-" << dsNum << endl;
  if (tcs) TCutSkimmer(theCut, dsNum);
  if (!tcs && sw) SkimWaveforms(theCut, inFile, outFile);
}


void SkimWaveforms(string theCut, string inFile, string outFile)
{
  // Take an input skim file, copy it with a waveform branch appended.
  // NOTE: The copied vectors are NOT resized to contain only entries passing cuts.
  // This is to preserve 1-1 matching with the other vectors in the copied file.

  TChain *skimTree = new TChain("skimTree");
  skimTree->Add(inFile.c_str());
  cout << "Found " << skimTree->GetEntries() << " input skim entries.\n";
  size_t n = skimTree->Draw(">>elist",theCut.c_str(), "entrylist");
  cout << "Draw successful.  Found " << n << " events passing cuts.\n";
  if (n == 0) {
    cout << "No events found passing cuts.  Exiting ...\n";
    return;
  }
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  skimTree->SetEntryList(elist);
  TFile *output = new TFile(outFile.c_str(),"RECREATE");
  cout << "Attempting tree copy ... \n";
  TTree *cutTree = skimTree->CopyTree("");
  cout << "Tree copy successful.\n";
  cutTree->Write("",TObject::kOverwrite);
  TNamed thisCut("theCut",theCut);	// save the cut used into the file.
  thisCut.Write();
  cout << Form("Using this cut:  \n%s  \nWrote %lli entries to the cut tree.\n",theCut.c_str(),cutTree->GetEntries());

  // Add waveforms to the cut tree, keeping only one run in memory at a time
  int run=0, iEvent=0;
  vector<double> *channel=0;
  cutTree->SetBranchAddress("run",&run);
  cutTree->SetBranchAddress("iEvent",&iEvent);
  cutTree->SetBranchAddress("channel",&channel);

  vector<MGTWaveform*> *waveVector=0;
  stack<MGTWaveform*> usedPointers;
  TBranch *waveBranch = cutTree->Branch("MGTWaveforms","vector<MGTWaveform*>",&waveVector,32000,0);

  GATDataSet *ds = new GATDataSet();
  string runPath = "";
  TChain *built = new TChain("MGTree");
  TChain *gat = new TChain("mjdTree");
  TTreeReader bReader(built);
  TTreeReader gReader(gat);
  TTreeReaderValue<TClonesArray> wfBranch(bReader,"fWaveforms");
  TTreeReaderValue<TClonesArray> wfAuxBranch(bReader,"fAuxWaveforms");
  TTreeReaderArray<double> wfChan(gReader,"channel");
  int prevRun=0;
  bool isMS=0, printMS=0;
  cout << "Adding waveforms to the cut tree ...\n";
  for (int i = 0; i < cutTree->GetEntries(); i++)
  {
    cutTree->GetEntry(i);
    if (run != prevRun) {
      built->Reset();
      gat->Reset();
      runPath = ds->GetPathToRun(run,GATDataSet::kBuilt);
      built->Add(runPath.c_str());
      bReader.SetTree(built);

      // detect multisampling
      TDirectory* tdir = gROOT->CurrentDirectory();
      TFile bFile(runPath.c_str());
      MJTRun* runInfo = (MJTRun*)bFile.Get("run");
      isMS = runInfo->GetUseMultisampling();
      if (printMS==0 && isMS==1) {
        cout << "Multisampling detected.\n";
        printMS=1;
      }
      gROOT->cd(tdir->GetPath());

      runPath = ds->GetPathToRun(run,GATDataSet::kGatified);
      gat->Add(runPath.c_str());
      gReader.SetTree(gat);
    }
    bReader.SetEntry(iEvent);
    gReader.SetEntry(iEvent);
    waveVector->resize(0);
    int nWF = (*wfBranch).GetEntries();

    // Figure out which hits made it into the skim file (some are cut by data cleaning)
    // This preserves the 1-1 matching between the skim vectors and the new MGTWaveform vector
    int numPass = cutTree->Draw("channel","","GOFF",1,i);
    double *channelList = cutTree->GetV1();
    vector<double> chanVec;
    for (int j = 0; j < numPass; j++) chanVec.push_back(channelList[j]);

    // Fill the waveform branch
    for (int iWF = 0; iWF < nWF; iWF++)
    {
      if ( find(chanVec.begin(), chanVec.end(), wfChan[iWF]) != chanVec.end() )
      {
        // handle multisampling
        // MGTWaveform *wave = dynamic_cast<MGTWaveform*>((*wfBranch).At(iWF));
        MGTWaveform* wave = NULL;
        if (!isMS) {
          MGTWaveform* reg = dynamic_cast<MGTWaveform*>((*wfBranch).At(iWF));
          wave = reg;
        }
        else {
          MGTWaveform* reg = dynamic_cast<MGTWaveform*>((*wfBranch).At(iWF));
          MGTWaveform* aux = dynamic_cast<MGTWaveform*>((*wfAuxBranch).At(iWF));
          MJTMSWaveform ms(reg,aux);
          wave = dynamic_cast<MGTWaveform*>(&ms);
        }

        // use a stack, don't clone wf's (huge memory leak)
        if(usedPointers.empty()) waveVector->push_back(new MGTWaveform);
      	else {
      	  waveVector->push_back(usedPointers.top());
      	  usedPointers.pop();
      	}
        *(waveVector->back()) = *wave;
        // cout << Form("   wave -> iWF %i  channel %.0f\n", iWF,wfChan[iWF]);
      }
    }
    waveBranch->Fill();

    // send any new pointers into usedPointers so they can be reused
    for (MGTWaveform* waveptr : *waveVector) {
      waveptr->Clear();
      usedPointers.push(waveptr);
    }
    waveVector->clear();

    // update progress and save
    if (i%5000==0 && i!=0) {
      cout << i << " saved, " << 100*i/(double)cutTree->GetEntries() << "% done.\n";
      cutTree->Write("", TObject::kOverwrite);
    }
    prevRun = run;
  }

  // Save and quit
  cutTree->Write("",TObject::kOverwrite);
  output->Close();

  // delete ds;
  cout << "Wrote file: " << outFile << endl;
}

void TCutSkimmer(string theCut, int dsNum)
{
  // There is a simpler version of this routine in data-cleaning.cc
  // if you don't want all the file I/O at the top here.
  // Cut down a skim file by saving only entries that pass the cut
  // into a new file.  (Very useful for low-threshold skim files).

  cout << "Skimming data set " << dsNum << " using this cut: " << theCut << "\n\n";

  // Find all files matching the expression
  string command = Form("ls ~/project/cal-skim/skimDS%i*",dsNum);
  array<char, 128> buffer;
  vector<string> files;
  string str;
  shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);
  if (!pipe) throw runtime_error("popen() failed!");
  while (!feof(pipe.get())) {
    if (fgets(buffer.data(), 128, pipe.get()) != NULL) {
      str = buffer.data();
      str.erase(remove(str.begin(), str.end(), '\n'), str.end());  // strip newline
      files.push_back(str);
    }
  }

  // Apply theCut to each file and make a new one
  for (auto file : files)
  {
    cout << "input: " << file << endl;

    // Get the run range and create an output file
    string range = file.substr(file.find(Form("DS%i",dsNum))+2);
    string ext = "_cgw.root";
    string::size_type i = range.find(ext);
    if (i != string::npos) range.erase(i, ext.length());
    string outFile = "~/project/cal-skim-basic/skim-basicDS" + range + ".root";
    cout << "output: " << outFile << endl;

    // Skim the input file with theCut
    TChain *skim = new TChain("skimTree");
    skim->Add(file.c_str());
    skim->Draw(">>entList",theCut.c_str(),"entrylist GOFF");
    TEntryList *entList = (TEntryList*)gDirectory->Get("entList");
    skim->SetEntryList(entList);
    TFile *f2 = new TFile(outFile.c_str(),"recreate");
    TTree *small = skim->CopyTree("");
    small->Write();
    cout << "Wrote " << small->GetEntries() << " entries.\n";
    TNamed thisCut("cutUsedHere",theCut);	// save the cut used into the file.
    thisCut.Write();
    f2->Close();
  }
}

void diagnostic()
{
  TFile *f = new TFile("./data/waveSkimDS4_test.root");
  TTree *t = (TTree*)f->Get("skimTree");

  vector<MGTWaveform*> *waves=0;
  int iEvent=0, run=0;
  vector<double> *channel=0;
  t->SetBranchAddress("MGTWaveforms",&waves);
  t->SetBranchAddress("channel",&channel);
  t->SetBranchAddress("run",&run);
  t->SetBranchAddress("iEvent",&iEvent);

  for (int i = 0; i < t->GetEntries(); i++)
  {
    t->GetEntry(i);
    cout << Form("i %i  iEvent %i  run %i  size channel %lu  size MGT %lu\n", i,iEvent,run,channel->size(),waves->size());

    if (channel->size() != waves->size()) {
      cout << "     Warning!! Them sizes ain't equal!  Something went wrong filling branches.\n";
      return;
    }
  }
}
