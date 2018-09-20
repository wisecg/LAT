#ifndef _GPX_FITTER_HH_
#define _GPX_FITTER_HH_

/*

        When a Wenqin gains enough statistics experience to reach level 36, he evolves into the Great Professor X (GPX)
        GPX is the final form of Wenqin. It can mega evolve into Mega GPX X using a Wenqinite X or Mega GPX Y using a Wenqinite Y

*/

#include <string>
#include <vector>
#include <map>

class RooDataSet;
class RooFitResult;
class RooPlot;
class RooRealVar;
class RooAbsReal;
class RooWorkspace;
class RooAbsPdf;
class RooMinimizer;
class RooMCStudy;
class RooEffProd;
class RooFormulaVar;
class TChain;
class TH1D;
class RooHistPdf;
class GPXFitter {

	public:
		GPXFitter();

		GPXFitter(std::string ds, double fitMin, double fitMax, std::string mode) {fDS = ds; fFitMin = fitMin; fFitMax = fitMax; fMode = mode;}

		virtual ~GPXFitter();

        // Constructs model PDF -- MUST be called after LoadData()!
        virtual void ConstructPDF(bool bNoEff = false, bool bWFMode = false, std::string fCPD="");

        // Do the fit
        // Set minimizer type here also... there's like Minuit, Minuit2, probably more
        // Honestly Minuit2 and Minuit are the same shit, one's just potentially slightly faster
        virtual void DoFit(std::string Minimizer = "Minuit2");
				// virtual void DoFitEff(std::string Minimizer = "Minuit2");

        // Draws and saves a plot of the fit as well as correlation matrix -- default binning is 0.2 keV
        // Binning is simply for visualization!
        void DrawBasic(double binSize = 0.2, bool drawLabels = true, bool drawResid = true, bool drawMatrix = true);

				void DrawModels(double binSize = 0.2);

        // Draws and saves contour plot -- arguments must have same name as in ConstructPDF()!
        // Parameters that become limited will take forever (as in never finish)
        void DrawContour(std::string argN1 = "Tritium", std::string argN2 = "Ge68");

        // This function calculates the energy resolution as a function of energy
        // According to the BDM PRL paper -- https://arxiv.org/abs/1612.00886
        double GetSigma(double energy);

        // Wrapper for returning fit values for variables and their errors
        std::vector<double> GetVar(std::string argN);

        // This function uses RooMCStudy to generate toy MC and then fit the results
        void GenerateMCStudy(std::vector<std::string> argS = {"Tritium"}, int nMC = 5000);

        // This function generates Toy MC data according to the best fit model and saves to a file
        void GenerateToyMC(std::string fileName);

        // Get Fit Result
        RooFitResult *GetFitResult() {return fFitResult;}

				// Get Workspace
        RooWorkspace *GetWorkspace() {return fFitWorkspace;}

        // Load data from skim TChain with a TCut
        // This assumes standard skimTree format
        void LoadChainData(TChain *skimTree, std::string theCut);

				// Load histogram from file
        void LoadHistData(TChain *skimTree, std::string theCut);

        // Creates, draws, and saves Profile Likelihood -- argument must have same name as in ConstructPDF()!
        // This is the ProfileNLL built into RooFit
        std::map<std::string, std::vector<double>> ProfileNLL(std::vector<std::string> argS = {"Tritium"}, double CL = 0.683);

        // Saves fit results into file
        void SaveShit(std::string outfileName = "TestOutput.root");

				// Manually set efficiency histogram
        void SetEfficiency(TH1D* effSpec){fEffSpec = effSpec;}

				// Set an exposure map different from default
        void SetExposureMap(std::map<std::string, std::vector<double>> fExposure){fExposureMap = fExposure;}

        // Sets range for fit
        void SetFitRange(double fitMin, double fitMax);

				// Sets initial value of parameter
        void SetParameter(std::string arg, double val, bool fix = false);

        // Sets prefix of output files
        void SetSavePrefix(std::string savePrefix) {fSavePrefix = savePrefix;}

        // String to append to beginning of saved files
        std::string fSavePrefix;

        double fChiSquare;

	private:
				// Dataset label
				std::string fDS;
				std::string fMode;

				// Fit range -- in keV
				double fFitMin;
				double fFitMax;

				// String used for cuts
				std::string fCutString;

				// Exposure Map
				std::map<std::string, std::vector<double>> fExposureMap;

				// Efficiency Scaling Map
				std::map<std::string, double> fEffMap;

				// Energy
        RooRealVar *fEnergy;
				RooFormulaVar *fEnergyShift;

        // Real dataset (as in real data)
        RooDataSet *fRealData;

        // Toy MC data -- This is for studying systematics...
        RooDataSet *fMCData;
        RooMCStudy *fMCStudy;

        // Total PDF -- should change to RooSimultaneous for simultaneous fits
        RooAbsPdf *fModelPDF;

				// Total PDF -- with efficiencies
				// RooAbsPdf *fModelPDFFinal;
				// RooEffProd *fModelPDFFinal;
    	  TH1D *fEffSpec;
				bool bCustomEff;
				bool fWFMode;
				bool fDetMode;

        // Minimizer
        RooMinimizer *fMinimizer;

        // NLL and ProfileNLL
        RooAbsReal *fNLL;
        RooAbsReal *fProfileNLL;

        // Saved fit result
        RooFitResult *fFitResult;

        // Fit workspace
        RooWorkspace *fFitWorkspace;

};

#endif
