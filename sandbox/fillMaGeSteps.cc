// Takes Pinghan's code to output all steps from raw MaGe file into a csv
#include "MGTMCEventSteps.hh"
#include "MGTMCStepData.hh"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TProof.h"
#include "TROOT.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TObject.h"
#include "TString.h"
#include <map>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <getopt.h>

using namespace std;

int main(int argc, char** argv)
{
  if(argc != 3) {
    cout << "Usage: " << argv[0] << " [File] [OutName (no .csv)]" << endl;
    return 1;
  }

  string fName = argv[1];
  string fOutName = argv[2];
  cout.precision(15);

  TChain* fTree = new TChain("fTree");
  fTree->Add(Form("%s", fName.c_str()));
  cout << "Added files" << endl;
  cout << "-------" << endl;
  Long64_t nentries = (Long64_t)fTree->GetEntries();
  cout << fTree->GetNtrees() << " trees with " << nentries << " entries " <<endl;
  cout << "-------" << endl;

  //for tree
  MGTMCEventSteps *eventSteps = 0;
  fTree->SetBranchAddress("eventSteps", &eventSteps);

  const MGTMCStepData *step;

  TString StepVolume;
  TObjArray* VolumePtr;
  Double_t timestamp;
  Double_t Edep;
  Int_t ID = 0;
  Int_t ParticleID;
  Int_t ParentID;
  //Int_t EventId;
  Int_t A, Z, L;
  Double_t Weight;
  Double_t fX,fY,fZ;
  Double_t fPx,fPy,fPz,fP, kinE;
  TString Process;
  Int_t TrackID;
  cout << "Total Entries:" << nentries << endl;

  ofstream fevent(Form("%s.csv",fOutName.c_str()));
  fevent << "iEvent" << "," << "iStep" <<"," << "ParticleID" <<"," << "A" <<","<< "Z"<<"," <<"L"<<"," <<"TrackID" <<"," <<  "ParentID" <<","<< "timestamp" << "," << "KineticE" << "," << "Edep" << "," << "fX" << "," << "fY" << "," << "fZ"<<"," <<"fPx"<< "," << "fPy" << "," << "fPz" <<"," << "Weight" <<","<< "StepVolume" << endl;

  for (Int_t i=0;i<nentries;i++)
  {
    fTree->GetEntry(i);

    for (Int_t k=0; k<eventSteps->GetNSteps(); k++)
    {
      step=eventSteps->GetStep(k);
      ParticleID=step->GetParticleID();
      Weight = step->GetTrackWeight();
      A = ((ParticleID-1000000000)/10)%1000;
      Z = ((ParticleID-1000000000)-A*10)/10000;
      //  L = step->GetNuclearL();
      L = (ParticleID - 1000000000)/10000000;
      fX = step->GetX();
      fY = step->GetY();
      fZ = step->GetZ();
      fPx = step->GetPx();
      fPy = step->GetPy();
      fPz = step->GetPz();
      kinE = step->GetKineticE();
      ParentID = step->GetParentTrackID();
      TrackID = step->GetTrackID();
      Edep = step->GetEdep();
      timestamp = step->GetT();
      StepVolume="";
      StepVolume+=step->GetPhysVolName();
      fevent << i << "," << k <<"," << ParticleID <<"," << A <<","<< Z<<"," <<L<<"," <<TrackID <<"," <<  ParentID <<","<< timestamp*1e-9 << ","<< kinE << "," << Edep << ","<< fX << "," << fY << "," << fZ<<"," << fPx<< "," << fPy << "," << fPz <<","<<Weight <<","<<StepVolume << endl;
    }
  }
}
