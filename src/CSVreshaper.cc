#include "../interface/CSVreshaper.h"

using namespace std;

// Constructor
CSVreshaper::CSVreshaper(CSVmultiplet<double>& iScaleFactors, TH1F* iHistogram){ 

	// Original (Old) working points
	CSVmultiplet<double> CSVthresholds(0.244, 0.679, 0.898);

	// Feed new and old CSV working point thresholds to the interpolator
	vector<double> newWP; newWP.clear();
	vector<double> oldWP; oldWP.clear();
	newWP.push_back(0.0);																					oldWP.push_back(0.0);
	CSVthresholds.GetLoose();
	iScaleFactors.GetLoose();
	newWP.push_back(GetNewThreshold(*CSVthresholds.GetLoose(),	*iScaleFactors.GetLoose(),	iHistogram));	oldWP.push_back(*CSVthresholds.GetLoose());
	newWP.push_back(GetNewThreshold(*CSVthresholds.GetMedium(),	*iScaleFactors.GetMedium(),	iHistogram));	oldWP.push_back(*CSVthresholds.GetMedium());
	newWP.push_back(GetNewThreshold(*CSVthresholds.GetTight(),	*iScaleFactors.GetTight(), 	iHistogram));	oldWP.push_back(*CSVthresholds.GetTight());
	newWP.push_back(1.001);																					oldWP.push_back(1.001);

	interpolator =  new ROOT::Math::Interpolator(newWP, oldWP, ROOT::Math::Interpolation::kLINEAR);
}

// Destructor
CSVreshaper::~CSVreshaper(){}

// Clean pointers
void CSVreshaper::Clean(){
	if(interpolator != NULL){ delete interpolator; interpolator = NULL; }
}


// If something goes really wrong, inform and quit
void CSVreshaper::ThrowFatalError(string const iMessage){ cerr << "[ERROR]\t" << iMessage << " Cannot continue. Terminating..." << endl; exit(1); }

// Return reshapped CSV value
double CSVreshaper::GetReshapedCSVvalue(double iOriginalCSVvalue){
	return interpolator->Eval(iOriginalCSVvalue);
}

// Return the new CSV WP threshold
double CSVreshaper::GetNewThreshold(double iOriginalThreshold, double iWPscaleFactor, TH1F* iHistogram){

	// Set up what we have and what we want
	double originalThresholdBin	= iHistogram->FindBin(iOriginalThreshold);
	unsigned int lastBin		= iHistogram->GetNbinsX()+1;
	double originalIntegral		= iHistogram->Integral(originalThresholdBin, lastBin);
	double targetIntegral		= originalIntegral * iWPscaleFactor;

	// Find and return the bin low edge whose integral to the last bin provides at least the target integral amount
	for(int bin = lastBin; bin != -1; bin--){
		double newIntegral = iHistogram->Integral(bin, lastBin);
		if(newIntegral >= targetIntegral){ return (iHistogram->GetBinLowEdge(bin)); }
	}

	// Should never get here
	ThrowFatalError("Could not meet target integral.");
	return -99;
}
