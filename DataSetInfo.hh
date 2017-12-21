#ifndef __DataSetInfo_hh__
#define __DataSetInfo_hh__

#include <iostream>
#include <map>
#include <cmath>
#include "glob.h"

#include "GATDataSet.hh"

using namespace std;

// ======================================================================
// Contents:
//
// FindDataSet - Returns the DS number of a given run.
// GetDataSetSequences - Returns a map of all the sub-ranges in a DS.
// GetDSRunAndStartTimes - Start time and run time of each DS
// LoadDataSet - Contains run ranges for all DS's, updates a GATDataSet.
// LoadDetectorList - Gives a list of the detector names for each module.
// GetTotalActiveMass - Total active mass for each dataset.
//                      This could be calculated instead of hardcoded in the future.
// LoadActiveMasses - Returns a map of all active masses.
// LoadActiveMassUncertainties - Returns a map of all active mass uncertainties.
// LoadBadDetectorMap - Returns a map of bad (i.e. not biased, unusuable) detectors.
// LoadVetoDetectorMap - Returns a map of veto-only detectors.
// GetChannelSelectionPath - Returns a string with the path to the highest
//                          version of the channel selection files.
// LoadEnrNatMap - quick way to tell if a given detID is enriched (1) or natural (0).
// CheckModule - Given a detector ID, look up which module it lives in.
// GetVetoActiveMass - Modifies total mass to not include veto-only detectors.
// GetENFC - Parameters for corrected trapENFCal
// GetENMC - Parameters for corrected trapENMCal
// GetAvsE - AvsE parameters
// GetDCR* - DCR parameters
// LoadDS4MuonList - Static muon list for DS-4, calculated manually
//                   by $GATDIR/mjd-veto/skim-veto.cc
// GetLNRunCoverage - Given a run and DS number, verify that this run is covered
//                    by the most recent LN Fill Tag.
// LoadLNFillTimes1 - Returns a vector of M1 LN fills.
// LoadLNFillTimes2 - Returns a vector of M2 LN fills.
//
// ======================================================================

// When LoadDataSet is updated, this must be updated as well.
int FindDataSet(int run);
// When LoadDataSet is updated, this must be updated as well.
int GetDataSetSequences(int dsNum);
void GetDSRunAndStartTimes(int dsNum, double &runTime_s, double &startTime0);
map<int, vector<int> > LoadDataSet(GATDataSet& ds, int dsNum, int subNum=-1);
vector<int> LoadDetectorList(int module);
void GetTotalActiveMass(int dsNum, double& m1Total, double& m1Enr, double& m1Nat,
  double& m2Total, double& m2Enr, double& m2Nat);
map<int,double> LoadActiveMasses(int dsNum);
map<int,double> LoadActiveMassUncertainties(int dsNum);
map<int, bool> LoadBadDetectorMap(int dsNum);
map<int,bool> LoadVetoDetectorMap(int dsNum);
std::string GetChannelSelectionPath(int dsNum, int officialVersion = -1);
map<int, bool> LoadEnrNatMap();
int CheckModule(int detID);
void GetVetoActiveMass(map<int,double> actM4Det_g, map<int,bool> detIDIsVetoOnly,
  double &mVeto_M1Total_kg, double &mVeto_M2Total_kg);
double GetENFC(int chan, int dsNum, double trapENF, int run);
double GetENMC(int chan, int dsNum, double trapENM, int run);
double GetAvsE(int chan, double TSCurrent50nsMax, double TSCurrent100nsMax, double TSCurrent200nsMax,
  double trapENF, double trapENFCal, int dsNum, int run);
double GetDCR90(int chan, double nlcblrwfSlope, double trapMax, int dsNum, int run);
double GetDCRCTC90(int chan, double nlcblrwfSlope, double trapE, double trapMax, int dsNum);
double GetDCR85(int chan, double nlcblrwfSlope, double trapMax, int dsNum, int run);
double GetDCR95(int chan, double nlcblrwfSlope, double trapMax, int dsNum, int run);
double GetDCR98(int chan, double nlcblrwfSlope, double trapMax, int dsNum, int run);
double GetDCR99(int chan, double nlcblrwfSlope, double trapMax, int dsNum, int run);
double GetDCR995(int chan, double nlcblrwfSlope, double trapMax, int dsNum, int run);
double GetDCR999(int chan, double nlcblrwfSlope, double trapMax, int dsNum, int run);
bool GetLNRunCoverage(int dsNum, int run);
void LoadDS4MuonList(vector<int> &muRuns, vector<double> &muRunTStarts,
  vector<double> &muTimes,vector<int> &muTypes, vector<double> &muUncert);
void LoadLNFillTimes1(vector<double>& lnFillTimes1, int dsNum);
void LoadLNFillTimes2(vector<double>& lnFillTimes2, int dsNum);

#endif
