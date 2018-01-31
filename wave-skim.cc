// Identify events in skim data passing
// a given TCut, and append the corresponding waveforms
// from built data into a new file.
// C. Wiseman, 1/18/2017
//         v.2 3/06/2017
//         v.3 5/18/2017 (thanks to Ian w/ waveform stack)
//         v.4 7/25/2017 - added multisampling support
//         v.5 1/18/2018 - added 'no TCut' mode

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
#include "MGWFNonLinearityCorrector.hh"
#include "MJTGretina4DigitizerData.hh"
#include "MJTypes.hh"

using namespace std;

void SkimWaveforms(string theCut, string inFile, string outFile, bool nlc);
void TCutSkimmer(string theCut, int dsNum);
void diagnostic();
void LoadNLCParameters(int ddID, int run, const MGVDigitizerData* dd, bool useTwoPass=true);

// stuff for NLC.  made 'em global because fk it.
map<int, MGWFNonLinearityCorrectionMap*> NLCMaps;
map<int, MGWFNonLinearityCorrectionMap*> NLCMaps2;
string NLCMapDir = "/global/project/projectdirs/majorana/data/production/NLCDB";

int main(int argc, char** argv)
{
  // hey let's get some (m)args
  if (argc < 2){
   cout << "Usage: ./wave-skim [options]\n"
        << "       [-s : run TCutSkimmer]\n"
        << "       [-r [dsNum] [subNum] : specify DS and sub-DS]\n"
        << "       [-f [dsNum] [runNum] : specify DS and run num]\n"
        << "       [-p [inPath] [outPath]: file locations]\n"
        << "       [-c : use calibration TCut]\n"
        << "       [-x : don't apply a data cleaning TCut]\n"
        << "       [-n : do the Radford 2-pass NLC correction]\n";
   return 1;
  }
  string inPath=".", outPath=".";
  int dsNum=-1, subNum=0, run=0;
  bool sw=0, tcs=0, fil=0, cal=0, longCal=0, nlc=0, noCut=0;
  vector<string> opt(argv, argv+argc);
  for (size_t i = 0; i < opt.size(); i++) {
    if (opt[i] == "-s") { sw=0; tcs=1; }
    if (opt[i] == "-f") { sw=1; fil=1; dsNum = stoi(opt[i+1]); run = stoi(opt[i+2]); }
    if (opt[i] == "-r") { sw=1; dsNum = stoi(opt[i+1]); subNum = stoi(opt[i+2]); }
    if (opt[i] == "-p") { inPath = opt[i+1]; outPath = opt[i+2]; }
    if (opt[i] == "-c") { cal=1; }
    if (opt[i] == "-l") { longCal=1; }
    if (opt[i] == "-x") { noCut=1; }
    if (opt[i] == "-n") {
     cout << "Performing Pass-2 Nonlinearity Correction ...\n";
     nlc=1;
    }
  }

  // 0nBB DS0-5 PRL standard bkg data cut (for reference)
  // string theCut = "!(channel==592 && run>=11348 && run<=11488) && !((channel & 0x2B0) == 0x2B0 && run >= 4239 && run <= 4436) && isGood && !wfDCBits && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && !muVeto && mHL==1 && avse>-1 && dcr99<0 && isEnr"

  // Low-E DS0-5 standard bkg data cut:
  string theCut = "!(channel==592 && run>=11348 && run<=11488) && !((channel & 0x2B0) == 0x2B0 && run >= 4239 && run <= 4436) && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0 && P!=0 && D!=0 && isGood && !muVeto && mH==1 && gain==0 && trapENFCal > 0.7";

  // DS0-5 calibration data cut
  if (cal) theCut = "trapENFCal > 0.7 && trapENFCal < 250 && (mH==1 || mH==2) && gain==0 && isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0 && P!=0 && D!=0";

  // Long calibration run cut
  if (longCal) theCut = "trapENFCal > 0.7 && isGood && gain==0 && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0 && P!=0 && D!=0";

  // No cut (for special runs)
  if (noCut) theCut = "";

  // Debug cut
  // theCut = "trapENFCal > 2 && mH==1 && isGood && !muVeto && !(C==1&&isLNFill1) && !(C==2&&isLNFill2) && C!=0&&P!=0&&D!=0";

  cout << "The cut is: " << theCut << endl;

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
  if (!tcs && sw) SkimWaveforms(theCut, inFile, outFile, nlc);
}


void SkimWaveforms(string theCut, string inFile, string outFile, bool nlc)
{
 // Take an input skim file, copy it with a waveform branch appended.
 // NOTE: The copied vectors are NOT resized to contain only entries passing cuts.
 // This is to preserve 1-1 matching with the other vectors in the copied file.
 TChain *skimTree = new TChain("skimTree");
 skimTree->Add(inFile.c_str());
 cout << "Found " << skimTree->GetEntries() << " input skim entries.\n";

 if (theCut != "") {
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
   output->Close();
 }
 else {
   TFile *output = new TFile(outFile.c_str(),"RECREATE");
   TTree *cutTree = skimTree->CopyTree("");
   cutTree->Write("",TObject::kOverwrite);
   cout << Form("No cuts applied.  Wrote %lli entries to the cut tree.\n",cutTree->GetEntries());
   output->Close();
 }
 TFile *output = new TFile(outFile.c_str(),"UPDATE");
 TTree *cutTree = (TTree*)output->Get("skimTree");

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
 TTreeReaderValue<TClonesArray> dd(bReader,"fDigitizerData");
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
       // don't handle multisampling
       // MGTWaveform *wave = dynamic_cast<MGTWaveform*>((*wfBranch).At(iWF));

       // handle multisampling
       MGTWaveform* wave = NULL;
       if (!isMS) {
         MGTWaveform* reg = dynamic_cast<MGTWaveform*>((*wfBranch).At(iWF));
         reg->SetWFEncScheme(MGTWaveform::kDiffVarInt); // it should copy this over, but just in case ...
         wave = reg;
       }
       else if (!nlc) {
         MGTWaveform* reg = dynamic_cast<MGTWaveform*>((*wfBranch).At(iWF)); // downsampled wf
         MGTWaveform* aux = dynamic_cast<MGTWaveform*>((*wfAuxBranch).At(iWF)); // fully sampled wf
         MJTMSWaveform ms(reg,aux);
         ms.SetWFEncScheme(MGTWaveform::kDiffVarInt);
         wave = dynamic_cast<MGTWaveform*>(&ms);
       }

       // do the 2-pass NLC correction from GAT-v01-06.  See below for a reference.
       if (nlc)
       {
         MGVDigitizerData* d = dynamic_cast<MGVDigitizerData*>((*dd).At(iWF));
         int ddID = d->GetID();
         LoadNLCParameters(ddID, run, d); // Adds NLC maps for this detector/run if they don't exist already

         MGWFNonLinearityCorrector* nlc = new MGWFNonLinearityCorrector();
         nlc->SetNLCCourseFineMaps(NLCMaps[ddID], NLCMaps2[ddID]);
         nlc->SetTimeConstant_samples(190); // 1.9 us time constant for Radford time-lagged method

         // In multisampled mode, we only can do the NLC on the fully sampled part of the WF.
         if (isMS) {
             MGTWaveform* reg = dynamic_cast<MGTWaveform*>((*wfBranch).At(iWF));  // downsampled wf
             MGTWaveform* aux = dynamic_cast<MGTWaveform*>((*wfAuxBranch).At(iWF));  // fully sampled wf
             nlc->TransformInPlace(*aux);
             MJTMSWaveform ms(reg,aux);
             ms.SetWFEncScheme(MGTWaveform::kDiffVarInt);
             wave = dynamic_cast<MGTWaveform*>(&ms);
         }
         else {
           nlc->TransformInPlace(*wave); // preferred method
         }
         delete nlc;
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
   if (i%10000==0 && i!=0) {
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
 // here's an example of grabbing an MGVDigitizerData object.
 GATDataSet ds;
 string path = ds.GetPathToRun(12345, GATDataSet::kBuilt);
 TChain *b = new TChain("MGTree");
 b->Add(path.c_str());
 TTreeReader reader(b);
 TTreeReaderValue<TClonesArray> dd(reader,"fDigitizerData");
 int count = 0;
 while (reader.Next())
 {
   int nEnt = (*dd).GetEntries();
   for (int i = 0; i < nEnt; i++)
   {
     MGVDigitizerData* d = dynamic_cast<MGVDigitizerData*>((*dd).At(i));
     d->SmartDump();
   }
   count++;
   if (count > 5) return;
 }
 return;
}


void LoadNLCParameters(int ddID, int run, const MGVDigitizerData* dd, bool useTwoPass)
{
 // Adapted by Clint on 6 Aug. 2017 from GAT-v01-06.
 // (There was no outside-GAT function I could call.)
 // 'useTwoPass' is set to true by default.
 // Direct path to the original file & GAT revision:
 // https://github.com/mppmu/GAT/tree/de49b732d95cf78265f94033b2fc609238742e8e

 MGWFNonLinearityCorrectionMap* map1 = NLCMaps[ddID];

 // This is a hack. Eventually we need to have these maps in a DB and pull them out properly.
 if(map1 == NULL)
 {
   map1 = new MGWFNonLinearityCorrectionMap;
   if((run >= 6000000 && run < 60000000) || run > 60002419) {
     static map< int, bool > seen;
     if(!seen[run]) {
       cout << "GATNonLinearityCorrector: No NLC files for run " << run << endl;
       seen[run] = true;
     }
     NLCMaps[ddID] = map1;
     return;
   }

   int crate = MJUtil::GetCrate(ddID);
   int card = MJUtil::GetCard(ddID);
   int channel = MJUtil::GetChannel(ddID);

   if(useTwoPass)
   {
     const MJTGretina4DigitizerData* g4dd = dynamic_cast<const MJTGretina4DigitizerData*>(dd);
     uint32_t boardSN = 0;
     if(g4dd == NULL) {
       cout << "GATNonLinearityCorrector::LoadParameters("
            << ddID << ", " << run << "): "
      << "Error: couldn't cast dd to MJTGretina4DigitizerData" << endl;
     }
     else boardSN = g4dd->GetBoardSerialNumber();

     // The board labeled SN-021h return 0x221 when probed by ORCA
     // Our NLC folders use the board labels so change this one.
     if(boardSN == 0x221) boardSN = 0x21;

     if(run < 6000000) {
       if(boardSN == 0 && run < 8184) {
         // No boardSN in part of DS0.
         if(card == 4) boardSN = 0x1b;
         else if(card == 5) boardSN = 0x19;
         else if(card == 6) boardSN = 0x22;
         else if(card == 7) boardSN = 0x12;
         else if(card == 8) boardSN = 0x25;
         else if(card == 9) boardSN = 0x26;
         else if(card == 10) boardSN = 0x18;
         else if(card == 11) boardSN = 0x21;
         else {
                 cout << "GATNonLinearityCorrector::LoadParameters("
                      << ddID << ", " << run << "): "
          << "Got unknown DS0 card number " << card << endl;
         }
       }
       if(run < 11339) {
         // card 019h is in slot 5 but Radford didn't make this NLC files,
         // so use the ones from crate 2 slot 15
         if(boardSN == 0x19 && card == 5) {
           crate = 2;
           card = 15;
         }
       }
       if(run >= 11339 && run <= 11396) {
         // in this run range, slot 5 used card 10h, whose NLCs were
         // measured in crate 2 slot 9
         if(boardSN == 0x10 && card == 5) {
           crate = 2;
           card = 9;
         }
       }
       if(run >= 18643) {
         // slot 11 changes from 014h to 029h, but 029h was only
         // calibrated in slot 8
         if(boardSN == 0x29 && card == 11) card = 8;
       }
       if(run >= 18990) {
         // 014h is usually used in slot 13 in this run range,
         // but 014h was only calibrated in slot 11
         if(boardSN == 0x14 && card == 13) card = 11;
       }
     }
     else { // DS4 runs
       if(run == 60002395 || run == 60002396) {
         // 01ah moved from slot 7 to slot 6 for 2 runs.
         if(boardSN == 0x1a && card == 6) card = 7;
       }
     }

     char bsnString[8];
     if(boardSN > 0xf) sprintf(bsnString, "%03xh", boardSN);
     else sprintf(bsnString, "%03d", boardSN);
     char fileName1a[500];
     sprintf(fileName1a,"%s/Boards/%s/c%dslot%d/Crate%d_GRET%d_Ch%d_part1a.dat",
             NLCMapDir.c_str(), bsnString, crate, card, crate, card, channel);
     map1->LoadFromCombinedFile(fileName1a);
     NLCMaps[ddID] = map1;

     MGWFNonLinearityCorrectionMap* map2 = new MGWFNonLinearityCorrectionMap;
     char fileName2a[500];
     sprintf(fileName2a,"%s/Boards/%s/c%dslot%d/Crate%d_GRET%d_Ch%d_part2a.dat",
             NLCMapDir.c_str(), bsnString, crate, card, crate, card, channel);
     map2->LoadFromCombinedFile(fileName2a);
     NLCMaps2[ddID] = map2;
   }

   else {
     string path = NLCMapDir;
     if(run < 11339) path += "/Run0";
     else if(run < 11397) path += "/Run11339";
     else if(run < 18622) path += "/Run11397";
     else if(run < 18643) path += "/Run18623";
     else if(run < 18990) path += "/Run18643";
     else if(run < 19018) path += "/Run18990";
     else if(run < 19502) path += "/Run19018";
     else if(run < 6000000) path += "/Run19502";
     else if(run < 60000000) {
       static map< int, map< int, bool> > seen;
       if(!seen[run][ddID]) {
         cout << "GATNonLinearityCorrector: No NLC files for run " << run << " ddID " << ddID << endl;
         seen[run][ddID] = true;
       }
       return;
     }
     else if(run < 60002395) path += "/Run60000000";
     else if(run < 60002397) path += "/Run60002395";
     else if(run < 60002419) path += "/Run60002397";
     else {
       static map< int, map< int, bool> > seen;
       if(!seen[run][ddID]) {
         cout << "GATNonLinearityCorrector: No NLC files for run " << run << " ddID " << ddID << endl;
         seen[run][ddID] = true;
       }
       return;
     }
     if(crate == 1) {
       char upFileName[500];
       sprintf(upFileName,"%s/Crate%d_GRET%d_Ch%d_up1.dat", path.c_str(), crate, card, channel);
       char downFileName[500];
       sprintf(downFileName,"%s/Crate%d_GRET%d_Ch%d_up0.dat", path.c_str(), crate, card, channel);
       map1->LoadFromUpDownFiles(upFileName, downFileName);
     }
     else if(crate == 2) {
       char combFileName[500];
       sprintf(combFileName,"%s/Crate%d_GRET%d_Ch%d_comb.dat", path.c_str(), crate, card, channel);
       map1->LoadFromCombinedFile(combFileName);
     }
     else {
       cout << "GATNonLinearityCorrector::LoadParameters("
            << ddID << ", " << run << "): Error: got invalid crate number "
            << crate << endl;
     }
     NLCMaps[ddID] = map1;
   }
 }
}
