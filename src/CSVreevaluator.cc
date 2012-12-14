#include "../interface/CSVreevaluator.h"

using namespace std;
using namespace boost::assign;

// Constructor
CSVreevaluator::CSVreevaluator(string iSampleName, string iEra, double iScaleBC, double iCharmFactor, double iScaleL){

	bottomFlavorReshapers	= NULL;
	charmFlavorReshapers	= NULL;
	lightFlavorReshapers	= NULL;

	// File and histo paths
	string pathToEfficiencyFile = (string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/");
	string histoName			= iSampleName;

    if (iEra == "2011") {
      pathToEfficiencyFile += "mc_btag_efficiency_7TeV.root";
      histoName += "_7TeV";
    }
    else if (iEra == "2012_52x" || iEra == "2012_53x") {
      pathToEfficiencyFile += "mc_btag_efficiency_8TeV.root";
      histoName += "_8TeV";
    }
    else {
      ThrowFatalError("Invalid era");
	}

	efficiencyFile = new TFile(pathToEfficiencyFile.c_str());
	TH1F* bottomFlavorBtagEfficiencyHistogram	= (TH1F*)efficiencyFile->Get(string(histoName + "_hb").c_str());
	TH1F* charmFlavorBtagEfficiencyHistogram	= (TH1F*)efficiencyFile->Get(string(histoName + "_hc").c_str());
	TH1F* lightFlavorBtagEfficiencyHistogram	= (TH1F*)efficiencyFile->Get(string(histoName + "_hl").c_str());


	// Scale factor container (where some "magic" numbers come from)
	CSVscaleFactorContainer scaleFactorContainer(iEra);


	// === Set up CSV reshapers === //
	// === Here, each flavor (b, c, l) gets its own binning in eta-pt space. For each eta-pt bin, a reshaper is created
	// Bottom flavor
	vector<double> bottomFlavorPtBins;	bottomFlavorPtBins	+= 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670;
	vector<double> bottomFlavorEtaBins;	bottomFlavorEtaBins	+= -2.5, 2.5;
	scaleFactorContainer.SetBottomFlavorBins(bottomFlavorEtaBins, bottomFlavorPtBins);

	bottomFlavorReshapers			= new XYmap<CSVreshaper>(bottomFlavorEtaBins, bottomFlavorPtBins);
	for(unsigned int bin = 0; bin < bottomFlavorReshapers->size(); bin++){
		double eta		= bottomFlavorReshapers->GetCenterX(bin);
		double pt		= bottomFlavorReshapers->GetCenterY(bin);
		double looseSF  = scaleFactorContainer.GetBottomFlavorScaleFactor(eta, pt, 'L', iScaleBC);
		double mediumSF = scaleFactorContainer.GetBottomFlavorScaleFactor(eta, pt, 'M', iScaleBC);
		double tightSF  = scaleFactorContainer.GetBottomFlavorScaleFactor(eta, pt, 'T', iScaleBC);
		CSVmultiplet<double> scaleFactors(looseSF, mediumSF, tightSF);
		bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomFlavorBtagEfficiencyHistogram));	
	}

	// Charm flavor
	// Take the same bins as for bottom flavor
	charmFlavorReshapers			= new XYmap<CSVreshaper>(bottomFlavorEtaBins, bottomFlavorPtBins);
	for(unsigned int bin = 0; bin < charmFlavorReshapers->size(); bin++){
		double eta		= charmFlavorReshapers->GetCenterX(bin);
		double pt		= charmFlavorReshapers->GetCenterY(bin);
		double looseSF  = scaleFactorContainer.GetCharmFlavorScaleFactor(eta, pt, 'L', iScaleBC, iCharmFactor);
		double mediumSF = scaleFactorContainer.GetCharmFlavorScaleFactor(eta, pt, 'M', iScaleBC, iCharmFactor);
		double tightSF  = scaleFactorContainer.GetCharmFlavorScaleFactor(eta, pt, 'T', iScaleBC, iCharmFactor);
		CSVmultiplet<double> scaleFactors(looseSF, mediumSF, tightSF);
		charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmFlavorBtagEfficiencyHistogram));	
	}

	// Light flavor
	vector<double> lightFlavorPtBins;	lightFlavorPtBins	+= 20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 670;
	vector<double> lightFlavorEtaBins;	lightFlavorEtaBins	+= 0, 0.5, 1.0, 1.5, 2.5;
	lightFlavorReshapers			= new XYmap<CSVreshaper>(lightFlavorEtaBins, lightFlavorPtBins);
	for(unsigned int bin = 0; bin < lightFlavorReshapers->size(); bin++){
		double eta		= charmFlavorReshapers->GetCenterX(bin);
		double pt		= charmFlavorReshapers->GetCenterY(bin);
		double looseSF  = scaleFactorContainer.GetLightFlavorScaleFactor(eta, pt, 'L', iScaleL);
		double mediumSF = scaleFactorContainer.GetLightFlavorScaleFactor(eta, pt, 'M', iScaleL);
		double tightSF  = scaleFactorContainer.GetLightFlavorScaleFactor(eta, pt, 'T', iScaleL);
		CSVmultiplet<double> scaleFactors(looseSF, mediumSF, tightSF);
		lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightFlavorBtagEfficiencyHistogram));	
	}


	// Clean up
	bottomFlavorBtagEfficiencyHistogram	= NULL;
	charmFlavorBtagEfficiencyHistogram	= NULL;
	lightFlavorBtagEfficiencyHistogram	= NULL;

}


// Destructor
CSVreevaluator::~CSVreevaluator(){
	delete bottomFlavorReshapers;	bottomFlavorReshapers	= NULL;
	delete charmFlavorReshapers;	charmFlavorReshapers	= NULL;
	delete lightFlavorReshapers;	lightFlavorReshapers	= NULL;
	delete efficiencyFile;			efficiencyFile = NULL;
}


// If something goes really wrong, inform and quit
void CSVreevaluator::ThrowFatalError(string const iMessage){ cerr << "[ERROR]\t" << iMessage << " Cannot continue. Terminating..." << endl; exit(1); }


// Return the corrected CSV value
double CSVreevaluator::GetReshapedCSVvalue(double iEta, double iPt, double iOriginalCSVvalue, int iFlavor){

	// CSV should not be >1
	if( iOriginalCSVvalue > 1){ ThrowFatalError("CSV value > 1."); }

	// Negative CSV indicates from problem with the algorithm
	if( iOriginalCSVvalue < 0 ){ return iOriginalCSVvalue; }

	// If CSV value is 1, don't do any reshaping
	if( fabs(iOriginalCSVvalue - 1) < 0.0001 ){ return iOriginalCSVvalue; }

	// === Reshape based on flavor === //
	double result = 0.0;
	switch(abs(iFlavor)){
		case 0:		result = iOriginalCSVvalue;																		break;
		case 4:		result = charmFlavorReshapers->GetObject(iEta, iPt)->GetReshapedCSVvalue(iOriginalCSVvalue);	break; // charm flavor
		case 5:		result = bottomFlavorReshapers->GetObject(iEta, iPt)->GetReshapedCSVvalue(iOriginalCSVvalue);	break; // bottom flavor
		default:	result = lightFlavorReshapers->GetObject(iEta, iPt)->GetReshapedCSVvalue(iOriginalCSVvalue);	break; // light flavor
	}

	return result;
}

