#include "../interface/CSVreevaluator.h"

using namespace std;
using namespace boost::assign;

// Constructor
CSVreevaluator::CSVreevaluator(string iSampleName, string iEra, double iScaleBC, double iCharmFactor, double iScaleL){

	bottomFlavorReshapers	= NULL;
	charmFlavorReshapers	= NULL;
	lightFlavorReshapers	= NULL;

	// File and histo paths
	string pathToEfficiencyFile1 = (string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/");
	string pathToEfficiencyFile2 = (string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/");
    string histoName1 = "";
    string histoName2 = "";

    if (iEra == "2011") {
      pathToEfficiencyFile1 += "mc_btag_efficiency_7TeV.root";
      histoName1 = iSampleName+"_7TeV";
    }
    else {
      pathToEfficiencyFile1 += "mc_btag_efficiency_8TeV.root";
      histoName1 = iSampleName+"_8TeV";
    }
    pathToEfficiencyFile2 += "mc_btag_efficiency_8TeV_53x.root";
    histoName2 = "ttbar_jj";

	efficiencyFile1 = new TFile(pathToEfficiencyFile1.c_str());
    efficiencyFile2 = new TFile(pathToEfficiencyFile2.c_str());

	TH1F* bottomFlavorBtagEfficiencyHistogram	= (TH1F*)efficiencyFile1->Get(string(histoName1 + "_hb").c_str());
	TH1F* charmFlavorBtagEfficiencyHistogram	= (TH1F*)efficiencyFile1->Get(string(histoName1 + "_hc").c_str());
	TH1F* lightFlavorBtagEfficiencyHistogram	= (TH1F*)efficiencyFile1->Get(string(histoName1 + "_hl").c_str());

	TH1F* bottomBtagEffHist_1A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt20to25_eta0p0to1p2_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_1A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt20to25_eta0p0to1p2_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_1A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt20to25_eta0p0to1p2_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_1B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt20to25_eta1p2to2p1_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_1B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt20to25_eta1p2to2p1_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_1B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt20to25_eta1p2to2p1_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_1C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt20to25_eta2p1to2p4_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_1C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt20to25_eta2p1to2p4_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_1C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt20to25_eta2p1to2p4_8TeV_hl").c_str());
    
	TH1F* bottomBtagEffHist_2A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt25to30_eta0p0to1p2_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_2A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt25to30_eta0p0to1p2_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_2A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt25to30_eta0p0to1p2_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_2B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt25to30_eta1p2to2p1_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_2B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt25to30_eta1p2to2p1_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_2B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt25to30_eta1p2to2p1_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_2C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt25to30_eta2p1to2p4_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_2C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt25to30_eta2p1to2p4_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_2C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt25to30_eta2p1to2p4_8TeV_hl").c_str());

	TH1F* bottomBtagEffHist_3A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt30to45_eta0p0to1p2_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_3A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt30to45_eta0p0to1p2_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_3A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt30to45_eta0p0to1p2_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_3B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt30to45_eta1p2to2p1_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_3B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt30to45_eta1p2to2p1_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_3B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt30to45_eta1p2to2p1_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_3C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt30to45_eta2p1to2p4_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_3C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt30to45_eta2p1to2p4_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_3C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt30to45_eta2p1to2p4_8TeV_hl").c_str());

	TH1F* bottomBtagEffHist_4A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt45to75_eta0p0to1p2_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_4A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt45to75_eta0p0to1p2_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_4A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt45to75_eta0p0to1p2_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_4B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt45to75_eta1p2to2p1_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_4B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt45to75_eta1p2to2p1_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_4B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt45to75_eta1p2to2p1_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_4C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt45to75_eta2p1to2p4_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_4C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt45to75_eta2p1to2p4_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_4C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt45to75_eta2p1to2p4_8TeV_hl").c_str());

	TH1F* bottomBtagEffHist_5A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt75to150_eta0p0to1p2_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_5A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt75to150_eta0p0to1p2_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_5A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt75to150_eta0p0to1p2_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_5B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt75to150_eta1p2to2p1_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_5B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt75to150_eta1p2to2p1_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_5B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt75to150_eta1p2to2p1_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_5C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt75to150_eta2p1to2p4_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_5C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt75to150_eta2p1to2p4_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_5C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt75to150_eta2p1to2p4_8TeV_hl").c_str());

	TH1F* bottomBtagEffHist_6A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt150toInf_eta0p0to1p2_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_6A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt150toInf_eta0p0to1p2_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_6A	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt150toInf_eta0p0to1p2_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_6B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt150toInf_eta1p2to2p1_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_6B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt150toInf_eta1p2to2p1_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_6B	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt150toInf_eta1p2to2p1_8TeV_hl").c_str());
	TH1F* bottomBtagEffHist_6C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt150toInf_eta2p1to2p4_8TeV_hb").c_str());
	TH1F* charmBtagEffHist_6C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt150toInf_eta2p1to2p4_8TeV_hc").c_str());
	TH1F* lightBtagEffHist_6C	= (TH1F*)efficiencyFile2->Get(string(histoName2+"_jet_pt150toInf_eta2p1to2p4_8TeV_hl").c_str());
    
	// Scale factor container (where some "magic" numbers come from)
	CSVscaleFactorContainer scaleFactorContainer(iEra);


	// === Set up CSV reshapers === //
	// === Here, each flavor (b, c, l) gets its own binning in eta-pt space. For each eta-pt bin, a reshaper is created
	// Bottom flavor

    vector<double> bottomFlavorPtBins;
    for(int i=0; i<47; ++i) {
      if (i < 20)      bottomFlavorPtBins.push_back(i*5.0); //Bins 5 GeV wide to 100 GeV
      else if (i < 30) bottomFlavorPtBins.push_back(100+10.0*(i-20)); //Bins 10 GeV wide to 200 GeV
      else             bottomFlavorPtBins.push_back(200+50.0*(i-30)); //Bins 50 GeV wide to 1000 GeV
    }
    vector<double> bottomFlavorEtaBins;
    bottomFlavorEtaBins.push_back(0.0);
    bottomFlavorEtaBins.push_back(1.2);
    bottomFlavorEtaBins.push_back(2.1);
    bottomFlavorEtaBins.push_back(2.5);

	scaleFactorContainer.SetBottomFlavorBins(bottomFlavorEtaBins, bottomFlavorPtBins);

	bottomFlavorReshapers			= new XYmap<CSVreshaper>(bottomFlavorEtaBins, bottomFlavorPtBins);
	for(unsigned int bin = 0; bin < bottomFlavorReshapers->size(); bin++){
		double eta		= bottomFlavorReshapers->GetCenterX(bin);
		double pt		= bottomFlavorReshapers->GetCenterY(bin);
		double looseSF  = scaleFactorContainer.GetBottomFlavorScaleFactor(0, pt, 'L', iScaleBC);
		double mediumSF = scaleFactorContainer.GetBottomFlavorScaleFactor(0, pt, 'M', iScaleBC);
		double tightSF  = scaleFactorContainer.GetBottomFlavorScaleFactor(0, pt, 'T', iScaleBC);
		CSVmultiplet<double> scaleFactors(looseSF, mediumSF, tightSF);
        if (iEra == "2011" || iEra == "2012_52x") {
		bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomFlavorBtagEfficiencyHistogram));
        }
        else if (iEra == "2012_53x") {
        if (pt < 20) { // duplicate of pt 20 - 25
          if (fabs(eta) < 1.2)      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_1A));
          else if (fabs(eta) < 2.1) bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_1B));
          else                      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_1C));
        }
        else if (pt < 25) {
          if (fabs(eta) < 1.2)      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_1A));
          else if (fabs(eta) < 2.1) bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_1B));
          else                      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_1C));
        }
        else if (pt < 30) {
          if (fabs(eta) < 1.2)      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_2A));
          else if (fabs(eta) < 2.1) bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_2B));
          else                      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_2C));
        }
        else if (pt < 45) {
          if (fabs(eta) < 1.2)      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_3A));
          else if (fabs(eta) < 2.1) bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_3B));
          else                      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_3C));
        }
        else if (pt < 75) {
          if (fabs(eta) < 1.2)      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_4A));
          else if (fabs(eta) < 2.1) bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_4B));
          else                      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_4C));
        }
        else if (pt < 150) {
          if (fabs(eta) < 1.2)      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_5A));
          else if (fabs(eta) < 2.1) bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_5B));
          else                      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_5C));
        }
        else {
          if (fabs(eta) < 1.2)      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_6A));
          else if (fabs(eta) < 2.1) bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_6B));
          else                      bottomFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, bottomBtagEffHist_6C));
        }
        }
        else {
          ThrowFatalError("Invalid era "+iEra);
        }
	}

	// Charm flavor
	// Take the same bins as for bottom flavor
	charmFlavorReshapers			= new XYmap<CSVreshaper>(bottomFlavorEtaBins, bottomFlavorPtBins);
	for(unsigned int bin = 0; bin < charmFlavorReshapers->size(); bin++){
		double eta		= charmFlavorReshapers->GetCenterX(bin);
		double pt		= charmFlavorReshapers->GetCenterY(bin);
		double looseSF  = scaleFactorContainer.GetCharmFlavorScaleFactor(0, pt, 'L', iScaleBC, iCharmFactor);
		double mediumSF = scaleFactorContainer.GetCharmFlavorScaleFactor(0, pt, 'M', iScaleBC, iCharmFactor);
		double tightSF  = scaleFactorContainer.GetCharmFlavorScaleFactor(0, pt, 'T', iScaleBC, iCharmFactor);
		CSVmultiplet<double> scaleFactors(looseSF, mediumSF, tightSF);
        if (iEra == "2011" || iEra == "2012_52x") {
		charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmFlavorBtagEfficiencyHistogram));	
        }
        else if (iEra == "2012_53x") {
        if (pt < 20) { // duplicate of pt = 20 - 25
          if (fabs(eta) < 1.2)      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_1A));
          else if (fabs(eta) < 2.1) charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_1B));
          else                      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_1C));
        }
        else if (pt < 25) {
          if (fabs(eta) < 1.2)      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_1A));
          else if (fabs(eta) < 2.1) charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_1B));
          else                      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_1C));
        }
        else if (pt < 30) {
          if (fabs(eta) < 1.2)      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_2A));
          else if (fabs(eta) < 2.1) charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_2B));
          else                      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_2C));
        }
        else if (pt < 45) {
          if (fabs(eta) < 1.2)      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_3A));
          else if (fabs(eta) < 2.1) charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_3B));
          else                      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_3C));
        }
        else if (pt < 75) {
          if (fabs(eta) < 1.2)      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_4A));
          else if (fabs(eta) < 2.1) charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_4B));
          else                      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_4C));
        }
        else if (pt < 150) {
          if (fabs(eta) < 1.2)      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_5A));
          else if (fabs(eta) < 2.1) charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_5B));
          else                      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_5C));
        }
        else {
          if (fabs(eta) < 1.2)      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_6A));
          else if (fabs(eta) < 2.1) charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_6B));
          else                      charmFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, charmBtagEffHist_6C));
        }
        }
        else {
          ThrowFatalError("Invalid era "+iEra);
        }
	}

	// Light flavor

    vector<double> lightFlavorPtBins;
    for(int i=0; i<47; ++i) {
      if (i < 20)      lightFlavorPtBins.push_back(i*5.0); //Bins 5 GeV wide to 100 GeV
      else if (i < 30) lightFlavorPtBins.push_back(100+10.0*(i-20)); //Bins 10 GeV wide to 200 GeV
      else             lightFlavorPtBins.push_back(200+50.0*(i-30)); //Bins 50 GeV wide to 1000 GeV
    }
    vector<double> lightFlavorEtaBins;
    lightFlavorEtaBins.push_back(0.0);
    lightFlavorEtaBins.push_back(0.5);
    lightFlavorEtaBins.push_back(1.2);
    lightFlavorEtaBins.push_back(1.5);
    lightFlavorEtaBins.push_back(2.1);
    lightFlavorEtaBins.push_back(2.5);

	lightFlavorReshapers			= new XYmap<CSVreshaper>(lightFlavorEtaBins, lightFlavorPtBins);
	for(unsigned int bin = 0; bin < lightFlavorReshapers->size(); bin++){
		double eta		= lightFlavorReshapers->GetCenterX(bin);
		double pt		= lightFlavorReshapers->GetCenterY(bin);
		double looseSF  = scaleFactorContainer.GetLightFlavorScaleFactor(eta, pt, 'L', iScaleL);
		double mediumSF = scaleFactorContainer.GetLightFlavorScaleFactor(eta, pt, 'M', iScaleL);
		double tightSF  = scaleFactorContainer.GetLightFlavorScaleFactor(eta, pt, 'T', iScaleL);
		CSVmultiplet<double> scaleFactors(looseSF, mediumSF, tightSF);
        if (iEra == "2011" || iEra == "2012_52x") {
          lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightFlavorBtagEfficiencyHistogram));	
        }
        else if (iEra == "2012_53x") {
        if (pt < 20) { // duplicate of pt = 20 - 25
          if (fabs(eta) < 1.2)      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_1A));
          else if (fabs(eta) < 2.1) lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_1B));
          else                      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_1C));
        }
        else if (pt < 25) {
          if (fabs(eta) < 1.2)      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_1A));
          else if (fabs(eta) < 2.1) lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_1B));
          else                      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_1C));
        }
        else if (pt < 30) {
          if (fabs(eta) < 1.2)      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_2A));
          else if (fabs(eta) < 2.1) lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_2B));
          else                      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_2C));
        }
        else if (pt < 45) {
          if (fabs(eta) < 1.2)      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_3A));
          else if (fabs(eta) < 2.1) lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_3B));
          else                      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_3C));
        }
        else if (pt < 75) {
          if (fabs(eta) < 1.2)      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_4A));
          else if (fabs(eta) < 2.1) lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_4B));
          else                      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_4C));
        }
        else if (pt < 150) {
          if (fabs(eta) < 1.2)      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_5A));
          else if (fabs(eta) < 2.1) lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_5B));
          else                      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_5C));
        }
        else {
          if (fabs(eta) < 1.2)      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_6A));
          else if (fabs(eta) < 2.1) lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_6B));
          else                      lightFlavorReshapers->Set(bin, CSVreshaper(scaleFactors, lightBtagEffHist_6C));
        }
        }
        else {
          ThrowFatalError("Invalid era "+iEra);
        }
	}


	// Clean up
	bottomFlavorBtagEfficiencyHistogram	= NULL;
	charmFlavorBtagEfficiencyHistogram	= NULL;
	lightFlavorBtagEfficiencyHistogram	= NULL;

	bottomBtagEffHist_1A = NULL;
	charmBtagEffHist_1A	 = NULL; 
	lightBtagEffHist_1A	 = NULL; 
	bottomBtagEffHist_1B = NULL;	
	charmBtagEffHist_1B	 = NULL; 
	lightBtagEffHist_1B	 = NULL;
	bottomBtagEffHist_1C = NULL;
	charmBtagEffHist_1C	 = NULL; 
	lightBtagEffHist_1C	 = NULL; 
    
	bottomBtagEffHist_2A = NULL;
	charmBtagEffHist_2A	 = NULL; 
	lightBtagEffHist_2A	 = NULL; 
	bottomBtagEffHist_2B = NULL;	
	charmBtagEffHist_2B	 = NULL; 
	lightBtagEffHist_2B	 = NULL;
	bottomBtagEffHist_2C = NULL;
	charmBtagEffHist_2C	 = NULL; 
	lightBtagEffHist_2C	 = NULL; 

	bottomBtagEffHist_3A = NULL;
	charmBtagEffHist_3A	 = NULL; 
	lightBtagEffHist_3A	 = NULL; 
	bottomBtagEffHist_3B = NULL;	
	charmBtagEffHist_3B	 = NULL; 
	lightBtagEffHist_3B	 = NULL;
	bottomBtagEffHist_3C = NULL;
	charmBtagEffHist_3C	 = NULL; 
	lightBtagEffHist_3C	 = NULL; 

	bottomBtagEffHist_4A = NULL;
	charmBtagEffHist_4A	 = NULL; 
	lightBtagEffHist_4A	 = NULL; 
	bottomBtagEffHist_4B = NULL;	
	charmBtagEffHist_4B	 = NULL; 
	lightBtagEffHist_4B	 = NULL;
	bottomBtagEffHist_4C = NULL;
	charmBtagEffHist_4C	 = NULL; 
	lightBtagEffHist_4C	 = NULL; 

	bottomBtagEffHist_5A = NULL;
	charmBtagEffHist_5A	 = NULL; 
	lightBtagEffHist_5A	 = NULL; 
	bottomBtagEffHist_5B = NULL;	
	charmBtagEffHist_5B	 = NULL; 
	lightBtagEffHist_5B	 = NULL;
	bottomBtagEffHist_5C = NULL;
	charmBtagEffHist_5C	 = NULL; 
	lightBtagEffHist_5C	 = NULL; 

	bottomBtagEffHist_6A = NULL;
	charmBtagEffHist_6A	 = NULL; 
	lightBtagEffHist_6A	 = NULL; 
	bottomBtagEffHist_6B = NULL;	
	charmBtagEffHist_6B	 = NULL; 
	lightBtagEffHist_6B	 = NULL;
	bottomBtagEffHist_6C = NULL;
	charmBtagEffHist_6C	 = NULL; 
	lightBtagEffHist_6C	 = NULL; 
}


// Destructor
CSVreevaluator::~CSVreevaluator(){
	delete bottomFlavorReshapers;	bottomFlavorReshapers	= NULL;
	delete charmFlavorReshapers;	charmFlavorReshapers	= NULL;
	delete lightFlavorReshapers;	lightFlavorReshapers	= NULL;
	delete efficiencyFile1;			efficiencyFile1 = NULL;
	delete efficiencyFile2;			efficiencyFile2 = NULL;
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

