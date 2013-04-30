#include "../interface/BEANhelper.h"

using namespace std;

BEANhelper::BEANhelper(){

	isSetUp = false;

	CSVLwp = 0.244;
	CSVMwp = 0.679;
	CSVTwp = 0.898;

	leptonSFfile	= NULL;

    // new files
    doubleMuonTriggerFile = NULL;
    doubleEleTriggerFile = NULL;
    muonEleTriggerFile = NULL;
    
	jetSFfile		= NULL;
	h_b_eff_		= NULL;
	h_c_eff_		= NULL;
	h_l_eff_		= NULL;
	h_o_eff_		= NULL;
	h_ele_SF_       = NULL;
	h_mu_SF_        = NULL;

	// PU reweighing
	puFile			= NULL;
	h_PU_ratio		= NULL;
	h_PUup_ratio	= NULL;
	h_PUdown_ratio	= NULL;

	// CSV reweighting
	for( int iSys=0; iSys<9; iSys++ ){
	  for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = NULL;
	  for( int iPt=0; iPt<3; iPt++ ){
	    for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
	  }
	}
    f_CSVwgt_HF = NULL;
    f_CSVwgt_LF = NULL;

	// CSV reshaping
	sh_				= NULL;
	sh_hfSFUp_		= NULL;
	sh_hfSFDown_	= NULL;
	sh_lfSFUp_		= NULL;
	sh_lfSFDown_	= NULL;

    samplename = "blank";

}


BEANhelper::~BEANhelper(){
	if(h_b_eff_ != NULL){ delete h_b_eff_; h_b_eff_ = NULL; }
	if(h_c_eff_ != NULL){ delete h_c_eff_; h_c_eff_ = NULL; }
	if(h_l_eff_ != NULL){ delete h_l_eff_; h_l_eff_ = NULL; }
	if(h_o_eff_ != NULL){ delete h_o_eff_; h_o_eff_ = NULL; }
	if(h_PU_ratio != NULL){ delete h_PU_ratio; h_PU_ratio = NULL; }
	if(h_PUup_ratio != NULL){ delete h_PUup_ratio; h_PUup_ratio = NULL; }
	if(h_PUdown_ratio != NULL){ delete h_PUdown_ratio; h_PUdown_ratio = NULL; }
	if(h_ele_SF_ != NULL){ delete h_ele_SF_; h_ele_SF_ = NULL; }
	if(h_mu_SF_ != NULL){ delete h_mu_SF_; h_mu_SF_ = NULL; }
	// CSV reweighting

	// CSV reweighting
	for( int iSys=0; iSys<9; iSys++ ){
	  for( int iPt=0; iPt<5; iPt++ ){ if(h_csv_wgt_hf[iSys][iPt] != NULL){ delete h_csv_wgt_hf[iSys][iPt]; h_csv_wgt_hf[iSys][iPt] = NULL;} }
	  for( int iPt=0; iPt<3; iPt++ ){
	    for( int iEta=0; iEta<3; iEta++ ){
	      if(h_csv_wgt_lf[iSys][iPt][iEta] != NULL){ delete h_csv_wgt_lf[iSys][iPt][iEta]; h_csv_wgt_lf[iSys][iPt][iEta] = NULL;}
	    }
	  }
	}

	// CSV reshaping
	if(sh_ != NULL){ delete sh_; sh_ = NULL; }
	if(sh_hfSFUp_ != NULL){ delete sh_hfSFUp_; sh_hfSFUp_ = NULL; }
	if(sh_hfSFDown_ != NULL){ delete sh_hfSFDown_; sh_hfSFDown_ = NULL; }
	if(sh_lfSFUp_ != NULL){ delete sh_lfSFUp_; sh_lfSFUp_ = NULL; }
	if(sh_lfSFDown_ != NULL){ delete sh_lfSFDown_; sh_lfSFDown_ = NULL; }
	// Close files
	if(leptonSFfile != NULL){ leptonSFfile->Close(); leptonSFfile = NULL; }
	if(jetSFfile != NULL){ jetSFfile->Close(); jetSFfile = NULL; }
	if(puFile != NULL){ puFile->Close(); puFile = NULL; }
	if(f_CSVwgt_HF != NULL){ f_CSVwgt_HF->Close(); f_CSVwgt_HF = NULL; }
	if(f_CSVwgt_LF != NULL){ f_CSVwgt_LF->Close(); f_CSVwgt_LF = NULL; }

    if (doubleMuonTriggerFile != NULL){ doubleMuonTriggerFile->Close(); doubleMuonTriggerFile = NULL;}
    if (doubleEleTriggerFile != NULL){ doubleEleTriggerFile->Close(); doubleEleTriggerFile = NULL;}
    if (muonEleTriggerFile != NULL) { muonEleTriggerFile->Close(); muonEleTriggerFile= NULL;}
}

// Set up parameters one by one
void BEANhelper::SetUp(string iEra, int iSampleNumber, bool iIsLJ, bool iIsData, string iDataset, bool iReshapeCSV, bool iPfLeptons = true, string iCollisionDS){
	// Make sure we don't set up more than once
	if(isSetUp){ ThrowFatalError("Trying to set up 'BEANhelper' for the second time. Check your code."); }

	// Bring in the external values
	era				= iEra;
	sampleNumber	= iSampleNumber;
	isLJ			= iIsLJ;
	isData			= iIsData;
	dataset			= iDataset;
	reshapeCSV		= iReshapeCSV;
	usePfLeptons    = iPfLeptons;

	// Error checking here
	if((era != "2011") && (era != "2012_52x") && (era != "2012_53x")){ ThrowFatalError("era set to '" + era + "' but it has to be either 2011, 2012_52x, or 2012_53x"); }
	if(sampleNumber==0){ ThrowFatalError("'sampleNumber' cannot be '0'."); }
	if(dataset.length()==0){ ThrowFatalError("'dataset' is blank."); }

	// Setup PU reweighing
	SetUpPUreweighing(iCollisionDS);

	// Setup CSV reweighting
	SetUpCSVreweighting();

	// Setup CSV reshaping
	SetUpCSVreshaping();

	// Setup jet efficiency scale factors
	SetUpJetSF();

	// Setup lepton efficiency scale factors
	SetUpLeptonSF();

	// Awknowledge setup
	isSetUp = true;

}

// Set up PU reweighing distributions
void BEANhelper::SetUpPUreweighing(string const iCollisionsDS){

	// Do nothing if we're running on collision data
	if(isData){ return; }

	// Get sample name from sample number
	string samplename = GetSampleName();

	// Don't get any of the histograms if the sample is data
	if( !(sampleNumber>0 && samplename!="data" && !isData) && !(sampleNumber<0 && samplename=="data" && isData) ) {
		cout << "sample number  '"	<< sampleNumber	<< "'" << endl;
		cout << "sample name    '"	<< samplename		<< "'" << endl;
		cout << "isData         "	<< isData			<< endl;
		assert (samplename == "sampleNumber, samplename, and isData inconsistent"); 
	}
	if( samplename == "blank") assert (samplename == "Why is samplename still blank?");
	if( sampleNumber<0 || samplename=="data" || isData) return;

	// Get PU file path
	puFile = new TFile ((string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/pu_distributions.root").c_str());

	// Set up mc histo pointer
	TH1D* h_pu_mc = NULL;

	// Figure out appropriate histo names to obtain from file
	string samplename_pu_input = samplename;

	if (samplename == "ttbar_bb" || samplename == "ttbar_cc")	samplename_pu_input = "ttbar";
	if (samplename == "zjets_lowmass")							samplename_pu_input = "zjets";

	// Histo names
	string mc_histo					= "";
	string collisions_histo			= "";
	string collisions_histo_up		= "";
	string collisions_histo_down	= "";

	if( era=="2012_53x" ){ // === 2012 === //

		// Simulation
		mc_histo				= "puDist_MC_8TeV_53x_S10";

		// Collision data
		vector<string> datasetNames;
		string datasets = iCollisionsDS;
		if(datasets == "All" || datasets == "all" || datasets == ""){ datasets = "2012A_13July,2012A_06Aug,2012B_13July,2012C_PR,2012C_24Aug,2012D_PR"; }

		datasets.erase( remove( datasets.begin(), datasets.end(), ' ' ), datasets.end() );
		stringstream datasets_ss(datasets);
		string item;
		while(getline(datasets_ss, item, ',')){
			if(item.length()==0){ continue; }
			datasetNames.push_back(item);
		} 
		
		for(unsigned int h=0; h<datasetNames.size(); h++){
				string histoName		= era + "/" + "puDist_" + datasetNames.at(h) + "_69300xSec";
				string histoUpName		= era + "/" + "puDist_" + datasetNames.at(h) + "_71795xSec";
				string histoDownName	= era + "/" + "puDist_" + datasetNames.at(h) + "_66805xSec";

				TH1D* h_PU_ratio_temp		= (TH1D*)puFile->Get(histoName.c_str());
				TH1D* h_PUup_ratio_temp		= (TH1D*)puFile->Get(histoUpName.c_str());
				TH1D* h_PUdown_ratio_temp	= (TH1D*)puFile->Get(histoDownName.c_str());
				if(h_PU_ratio_temp		== NULL){ cerr << "ERROR: '" + histoName		+ "' cannot be found in the PU distributions file." << endl; exit(1); }
				if(h_PUup_ratio_temp	== NULL){ cerr << "ERROR: '" + histoUpName		+ "' cannot be found in the PU distributions file." << endl; exit(1); }
				if(h_PUdown_ratio_temp	== NULL){ cerr << "ERROR: '" + histoDownName	+ "' cannot be found in the PU distributions file." << endl; exit(1); }

			if(h==0){
				h_PU_ratio		= (TH1D*)h_PU_ratio_temp->Clone();
				h_PUup_ratio	= (TH1D*)h_PUup_ratio_temp->Clone();  
				h_PUdown_ratio	= (TH1D*)h_PUdown_ratio_temp->Clone();
			}else{
				h_PU_ratio		->Add((TH1D*)puFile->Get(histoName.c_str()));
				h_PUup_ratio	->Add((TH1D*)puFile->Get(histoUpName.c_str()));
				h_PUdown_ratio	->Add((TH1D*)puFile->Get(histoDownName.c_str()));
			}
		}

		
	}else if( era=="2012_52x" ){ // === 2012 === //

		if (samplename == "ttbar" || samplename == "wjets" || samplename == "zjets" ||
				samplename == "ttH120_FullSim"|| samplename == "ttH120_FastSim")		samplename_pu_input = samplename + "_";
		else if (samplename == "ttbar_bb" || samplename == "ttbar_cc")					samplename_pu_input = "ttbar_";
		else if (samplename == "ttbarW" || samplename == "ttbarZ")						samplename_pu_input = "ttZorW_";
		else if (samplename == "zjets_lowmass")											samplename_pu_input = "zjets_";
		else																			samplename_pu_input = "";

		// Simulation
		mc_histo				= samplename_pu_input + "Summer2012_pileup_8TeV";

		// Collision data
		collisions_histo		= "pileup_8TeV_69300xSec";
		collisions_histo_up		= "pileup_8TeV_71795xSec";
		collisions_histo_down	= "pileup_8TeV_66805xSec";
		h_PU_ratio		= (TH1D*)puFile->Get(string(era + "/" + collisions_histo).c_str())->Clone();
		h_PUup_ratio	= (TH1D*)puFile->Get(string(era + "/" + collisions_histo_up).c_str())->Clone();
		h_PUdown_ratio	= (TH1D*)puFile->Get(string(era + "/" + collisions_histo_down).c_str())->Clone();

	}else if ( era=="2011") { // === 2011 === //

		string dataset_infix = dataset;
		if( !(dataset.find("SingleMu")!=string::npos || dataset.find("ElectronHad")!=string::npos) ) dataset_infix = "SingleMu";

		if( samplename == "ttH120" || samplename == "ttbarW" || samplename == "ttbarZ" ){
			// Simulation
			if( samplename == "ttH120" ){	mc_histo = "ttH_7TeV_numGenPV";	}
			else{							mc_histo = "ttV_7TeV_numGenPV";	}
		}else{
			// Simulation
			mc_histo = "F2011exp_7TeV";
		}

		// Collision data
		collisions_histo		= "pileup_7TeV_" + dataset_infix + "_68000_true";
		collisions_histo_up		= "pileup_7TeV_" + dataset_infix + "_73440_true";
		collisions_histo_down	= "pileup_7TeV_" + dataset_infix + "_62560_true";
		h_PU_ratio		= (TH1D*)puFile->Get(string(era + "/" + collisions_histo).c_str())->Clone();
		h_PUup_ratio	= (TH1D*)puFile->Get(string(era + "/" + collisions_histo_up).c_str())->Clone();
		h_PUdown_ratio	= (TH1D*)puFile->Get(string(era + "/" + collisions_histo_down).c_str())->Clone();

	}else{ // === Not 2011 or 2012 === //
		cout << "Era set to '" << era << "'" << endl;
		assert (era == "either 2012_52x, 2012_53x, or 2011");
	}

	// Get histos
	h_pu_mc			= (TH1D*)puFile->Get(string(era + "/" + mc_histo).c_str())->Clone();

	// Normalize every histogram to its area
	h_PU_ratio		->Scale( 1./h_PU_ratio		->Integral() );
	h_PUup_ratio	->Scale( 1./h_PUup_ratio	->Integral() );
	h_PUdown_ratio	->Scale( 1./h_PUdown_ratio	->Integral() );
	h_pu_mc			->Scale( 1./h_pu_mc			->Integral() );

	// Divide copy of the numerator by denominator to get ratio
	h_PU_ratio		->Divide( h_pu_mc );
	h_PUup_ratio	->Divide( h_pu_mc );
	h_PUdown_ratio	->Divide( h_pu_mc );

	// Delete pointer to mc histo
	if(h_pu_mc != NULL){ delete h_pu_mc; h_pu_mc = NULL; }

}


// Set up CSV reshaping
void BEANhelper::SetUpCSVreweighting(){

  // Do not set it up if we're running on collision data
  if(isData){ return; }
  f_CSVwgt_HF = new TFile ((string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/csv_rwt_hf_IT.root").c_str());
  f_CSVwgt_LF = new TFile ((string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/csv_rwt_lf_IT.root").c_str());


  // CSV reweighting
  for( int iSys=0; iSys<9; iSys++ ){
    TString syst_csv_suffix_hf = "final";
    TString syst_csv_suffix_lf = "final";
    
    switch( iSys ){
    case 0:
      // this is the nominal case
      break;
    case 1:
      // JESUp
      syst_csv_suffix_hf = "final_JESUp"; syst_csv_suffix_lf = "final_JESUp";
      break;
    case 2:
      // JESDown
      syst_csv_suffix_hf = "final_JESDown"; syst_csv_suffix_lf = "final_JESDown";
      break;
    case 3:
      // purity up
      syst_csv_suffix_hf = "final_LFUp"; syst_csv_suffix_lf = "final_HFUp";
      break;
    case 4:
      // purity down
      syst_csv_suffix_hf = "final_LFDown"; syst_csv_suffix_lf = "final_HFDown";
      break;
    case 5:
      // stats1 up
      syst_csv_suffix_hf = "final_Stats1Up"; syst_csv_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      // stats1 down
      syst_csv_suffix_hf = "final_Stats1Down"; syst_csv_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      // stats2 up
      syst_csv_suffix_hf = "final_Stats2Up"; syst_csv_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      // stats2 down
      syst_csv_suffix_hf = "final_Stats2Down"; syst_csv_suffix_lf = "final_Stats2Down";
      break;
    }

    for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)f_CSVwgt_HF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );

    for( int iPt=0; iPt<3; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)f_CSVwgt_LF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }

}


// Set up CSV reshaping
void BEANhelper::SetUpCSVreshaping(){

	// Do not set it up if we're running on collision data
	if(isData){ return; }

	// Tweak sample name if needed
	string samplename = GetSampleName();
    if (samplename == "ttbar_jj" || samplename == "ttbar_lj" || samplename == "ttbar_ll" ||
        samplename == "ttbar_bb_jj" || samplename == "ttbar_bb_lj" || samplename == "ttbar_bb_ll" ||
        samplename == "ttbar_cc_jj" || samplename == "ttbar_cc_lj" || samplename == "ttbar_cc_ll" ||
        samplename == "ttbar_bb" || samplename == "ttbar_cc")                  samplename = "ttbar";
    if (samplename == "wjets_1p" || samplename == "wjets_2p" || samplename == "wjets_3p" ||
        samplename == "wjets_4p")                                              samplename = "wjets";
	if (samplename == "zjets_1p" || samplename == "zjets_2p" || samplename == "zjets_3p" ||
        samplename == "zjets_4p" || samplename == "zjets_lowmass" ||
        samplename == "zjets_lowmass_1p" || samplename == "zjets_lowmass_2p")  samplename = "zjets";
    if (samplename == "t_schannel_ll")                                         samplename = "t_schannel";
    if (samplename == "t_tchannel_ll")                                         samplename = "t_tchannel";
    if (samplename == "t_tWchannel_jl" || samplename == "t_tWchannel_lj" ||
        samplename == "t_tWchannel_ll")                                        samplename = "t_tWchannel";
    if (samplename == "t_schannel_ll")                                         samplename = "t_schannel";
    if (samplename == "tbar_tchannel_ll")                                      samplename = "tbar_tchannel";
    if (samplename == "tbar_tWchannel_jl" || samplename == "tbar_tWchannel_lj" ||
        samplename == "tbar_tWchannel_ll")                                     samplename = "tbar_tWchannel";

    
	// Set charm scale factor
	double charmFactor = 2.0 - 1.0;

	// Set up CSVreevaluators for each systype
	sh_				= new CSVreevaluator(samplename, era,    0, charmFactor,  0);
	sh_hfSFUp_		= new CSVreevaluator(samplename, era,  1.5, charmFactor,  0);
	sh_hfSFDown_	= new CSVreevaluator(samplename, era, -1.5, charmFactor,  0);
	sh_lfSFUp_		= new CSVreevaluator(samplename, era,    0, charmFactor,  1);
	sh_lfSFDown_	= new CSVreevaluator(samplename, era,    0, charmFactor, -1);

}

// Set up jet efficiency scale factors
void BEANhelper::SetUpJetSF(){

	// Do not set it up if we're running on collision data
	if(isData){ return; }

	string com_suffix = "";
	string filePath = "";
	if( era=="2012_53x" ){
		filePath = string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/mc_btag_efficiency_8TeV_53x.root";
		com_suffix = "_8TeV";
	}else if( era=="2012_52x" ){
        filePath = string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/mc_btag_efficiency_8TeV.root";
		com_suffix = "_8TeV";
	}else if( era=="2011" ){
		filePath = string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/mc_btag_efficiency_7TeV.root";
		com_suffix = "_7TeV";
	}else{ // === Not 2011 or 2012 === //
		cout << "Era set to '" << era << "'" << endl;
		assert (era == "either 2012_52x, 2012_53x, or 2011");
	}

	jetSFfile = new TFile (filePath.c_str());

	string samplename = GetSampleName();
    if (samplename == "ttbar_jj" || samplename == "ttbar_lj" || samplename == "ttbar_ll" ||
        samplename == "ttbar_bb_jj" || samplename == "ttbar_bb_lj" || samplename == "ttbar_bb_ll" ||
        samplename == "ttbar_cc_jj" || samplename == "ttbar_cc_lj" || samplename == "ttbar_cc_ll" ||
        samplename == "ttbar_bb" || samplename == "ttbar_cc")                  samplename = "ttbar";
    if (samplename == "wjets_1p" || samplename == "wjets_2p" || samplename == "wjets_3p" ||
        samplename == "wjets_4p")                                              samplename = "wjets";
	if (samplename == "zjets_1p" || samplename == "zjets_2p" || samplename == "zjets_3p" ||
        samplename == "zjets_4p" || samplename == "zjets_lowmass" ||
        samplename == "zjets_lowmass_1p" || samplename == "zjets_lowmass_2p")  samplename = "zjets";
    if (samplename == "t_schannel_ll")                                         samplename = "t_schannel";
    if (samplename == "t_tchannel_ll")                                         samplename = "t_tchannel";
    if (samplename == "t_tWchannel_jl" || samplename == "t_tWchannel_lj" ||
        samplename == "t_tWchannel_ll")                                        samplename = "t_tWchannel";
    if (samplename == "t_schannel_ll")                                         samplename = "t_schannel";
    if (samplename == "tbar_tchannel_ll")                                      samplename = "tbar_tchannel";
    if (samplename == "tbar_tWchannel_jl" || samplename == "tbar_tWchannel_lj" ||
        samplename == "tbar_tWchannel_ll")                                     samplename = "tbar_tWchannel";

    
	h_b_eff_ = (TH2D*)jetSFfile->Get(string( samplename + com_suffix + "_jet_pt_eta_b_eff" ).c_str())->Clone();
	h_c_eff_ = (TH2D*)jetSFfile->Get(string( samplename + com_suffix + "_jet_pt_eta_c_eff" ).c_str())->Clone();
	h_l_eff_ = (TH2D*)jetSFfile->Get(string( samplename + com_suffix + "_jet_pt_eta_l_eff" ).c_str())->Clone();
	h_o_eff_ = (TH2D*)jetSFfile->Get(string( samplename + com_suffix + "_jet_pt_eta_o_eff" ).c_str())->Clone();

}

// Set up lepton efficiency scale factors
void BEANhelper::SetUpLeptonSF(){

	// Do not set it up if we're running on collision data
	if(isData){ return; }

	string filePath = "";
    
    string oldFilePath = "";
    
    if ( era=="2012_53x" ){
      
      filePath = string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/lepton_SF_8TeV_53x_baseline.root";
      
      oldFilePath = string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/lepton_SF_8TeV.root";

      oldLeptonScaleFactFile = new TFile (oldFilePath.c_str() );

      
	} else if( era=="2012_52x"  ){
      filePath = string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/lepton_SF_8TeV.root";
	}else if( era=="2011" ){
		filePath = string(getenv("CMSSW_BASE")) + "/src/NtupleMaker/BEANmaker/data/lepton_SF_8TeV.root";
	}else{ // === Not 2011 or 2012 === //
		cout << "Era set to '" << era << "'" << endl;
		assert (era == "either 2012_52x, 2012_53x, or 2011");
	}

	leptonSFfile = new TFile (filePath.c_str());


    // Careful! This might not be compatible with 2012_52x anymore.
	if( isLJ ){
	  //h_ele_SF_ = (TH2D*)leptonSFfile->Get(string( "h_ele_pt_eta_full_id_iso_hlt_pass" ).c_str())->Clone();
	  //h_ele_SF_ = (TH2D*)leptonSFfile->Get(string( "ele_pt_eta_full_id_iso_hlt_8TeV" ).c_str())->Clone();
      h_ele_SF_ = (TH2D*)leptonSFfile->Get(string( "TightEleIdIsoSF" ).c_str())->Clone();      
	  h_mu_SF_  = (TH2D*)leptonSFfile->Get(string( "mu_pt_eta_full_id_iso_hlt_8TeV" ).c_str())->Clone();
        
	}else {

	  h_ele_SF_ = (TH2D*)leptonSFfile->Get(string( "TightEleIdIsoSF" ).c_str())->Clone();
	  h_mu_SF_  = (TH2D*)leptonSFfile->Get(string( "mu_pt_eta_full_id_iso_8TeV" ).c_str())->Clone();
	}

	if ( era=="2012_53x") {
          h_doubleMuTrigSF  = (TH2D*) leptonSFfile->Get("TwoMuonTriggerSF")->Clone("DoubleMuSF");
          h_doubleEleTrigSF = (TH2D*) leptonSFfile->Get("TwoEleTriggerSF")->Clone("DoubleEleSF");
          h_muonEleTrigSF   = (TH2D*) leptonSFfile->Get("MuonEleTriggerSF")->Clone("MuonEleSF");


          h_SingleEle_trig_SF_ = (TH2D*) leptonSFfile->Get(string("TightEleTriggerSF").c_str())->Clone();
          
          h_looseMuonSF = (TH2D*)  leptonSFfile->Get(string( "LooseMuonIdIsoSF" ).c_str())->Clone();
          h_looseEleSF =  (TH2D*)  leptonSFfile->Get(string( "LooseEleIdIsoSF" ).c_str())->Clone();

          cout << "Getting lepton sf histos" << endl;
          h_testSingleMuNew = (TH2D*) leptonSFfile->Get("mu_pt_eta_full_id_iso_hlt_8TeV")->Clone("MuNewTest");          
          h_testSingleEleNew = (TH2D*) leptonSFfile->Get("TightEleIdIsoSF")->Clone("EleNewTest");
          h_testTrigSingleEleNew = (TH2D*) leptonSFfile->Get("TightEleTriggerSF")->Clone("EleTrigNewTest");

          h_testSingleMuOld = (TH2D*) oldLeptonScaleFactFile->Get("mu_pt_eta_full_id_iso_hlt_8TeV")->Clone("MuOldTest");
          h_testSingleEleOld = (TH2D*) oldLeptonScaleFactFile->Get("ele_pt_eta_full_id_iso_hlt_8TeV")->Clone("EleOldTest");          
          
          if (h_testSingleMuNew == NULL
              || h_testSingleEleNew == NULL
              || h_testSingleMuOld == NULL
              || h_testSingleEleOld == NULL
              || h_testTrigSingleEleNew == NULL ) {
            ThrowFatalError ("Can't read test histos");
          }
            
          cout << "Got scale factor Histos successfully" << endl;
	}
    

	

	if ( h_ele_SF_ == NULL
	     || h_mu_SF_ == NULL
	     || h_doubleMuTrigSF == NULL
	     || h_doubleEleTrigSF == NULL
	     || h_muonEleTrigSF == NULL
	     ) {

	  std::cout << "Sorry! could not find lepton SF histograms in file " << filePath
		    << "end. Without SF histograms I cannot continue. Crashing..."
		    << endl;

	  ThrowFatalError("NO LEP SF HISTOS");

	}

}


// Return sample name
string BEANhelper::GetSampleName(){

	string samplename = "";

	// 2012 eras
	if (era == "2012_53x" || era == "2012_52x") {
			 if ( sampleNumber<0 ){		samplename = "data";				}
		else if( sampleNumber==2500 ){	samplename = "ttbar";				}
		else if( sampleNumber==2544 ){	samplename = "ttbar_cc";			}
		else if( sampleNumber==2555 ){	samplename = "ttbar_bb";			}
		else if( sampleNumber==2511 ){	samplename = "ttbar_scaleup";		}
		else if( sampleNumber==2510 ){	samplename = "ttbar_scaledown";		}
		else if( sampleNumber==2513 ){	samplename = "ttbar_matchingup";	}
		else if( sampleNumber==2512 ){	samplename = "ttbar_matchingdown";	}
		else if( sampleNumber==2566 ){	samplename = "ttbar_jj";			}
		else if( sampleNumber==2563 ){	samplename = "ttbar_lj";			}
		else if( sampleNumber==2533 ){	samplename = "ttbar_ll";			}
		else if( sampleNumber==2576 ){	samplename = "ttbar_cc_jj";			}
		else if( sampleNumber==2573 ){	samplename = "ttbar_cc_lj";			}
		else if( sampleNumber==2543 ){	samplename = "ttbar_cc_ll";			}
		else if( sampleNumber==2586 ){	samplename = "ttbar_bb_jj";			}
		else if( sampleNumber==2583 ){	samplename = "ttbar_bb_lj";			}
		else if( sampleNumber==2553 ){	samplename = "ttbar_bb_ll";			}
		else if( sampleNumber==2400 ){	samplename = "wjets";				}
		else if( sampleNumber==2401 ){	samplename = "wjets_1p";			}
		else if( sampleNumber==2402 ){	samplename = "wjets_2p";			}
		else if( sampleNumber==2403 ){	samplename = "wjets_3p";			}
		else if( sampleNumber==2404 ){	samplename = "wjets_4p";			}
		else if( sampleNumber==2850 ){	samplename = "zjets_lowmass";		}
		else if( sampleNumber==2851 ){	samplename = "zjets_lowmass_1p";	}
		else if( sampleNumber==2852 ){	samplename = "zjets_lowmass_2p";	}
		else if( sampleNumber==2800 ){	samplename = "zjets";				}
		else if( sampleNumber==2801 ){	samplename = "zjets_1p";			}
		else if( sampleNumber==2802 ){	samplename = "zjets_2p";			}
		else if( sampleNumber==2803 ){	samplename = "zjets_3p";			}
		else if( sampleNumber==2804 ){	samplename = "zjets_4p";			}
		else if( sampleNumber==2600 ){	samplename = "t_schannel";			}
		else if( sampleNumber==2630 ){	samplename = "t_schannel_ll";		}
		else if( sampleNumber==2601 ){	samplename = "tbar_schannel";		}
		else if( sampleNumber==2631 ){	samplename = "tbar_schannel_ll";	}
		else if( sampleNumber==2602 ){	samplename = "t_tchannel";			}
		else if( sampleNumber==2632 ){	samplename = "t_tchannel_ll";		}
		else if( sampleNumber==2603 ){	samplename = "tbar_tchannel";		}
		else if( sampleNumber==2633 ){	samplename = "tbar_tchannel_ll";	}
		else if( sampleNumber==2604 ){	samplename = "t_tWchannel";			}
		else if( sampleNumber==2654 ){	samplename = "t_tWchannel_lj";		}
		else if( sampleNumber==2664 ){	samplename = "t_tWchannel_jl";		}
		else if( sampleNumber==2634 ){	samplename = "t_tWchannel_ll";		}
		else if( sampleNumber==2605 ){	samplename = "tbar_tWchannel";		}
		else if( sampleNumber==2655 ){	samplename = "tbar_tWchannel_lj";	}
		else if( sampleNumber==2665 ){	samplename = "tbar_tWchannel_jl";	}
		else if( sampleNumber==2635 ){	samplename = "tbar_tWchannel_ll";	}
		else if( sampleNumber==2700 ){	samplename = "ww";					}
		else if( sampleNumber==2710 ){	samplename = "www";					}
		else if( sampleNumber==2720 ){	samplename = "wwz";					}
		else if( sampleNumber==2701 ){	samplename = "wz";					}
		else if( sampleNumber==2721 ){	samplename = "wzz";					}
		else if( sampleNumber==2702 ){	samplename = "zz";					}
		else if( sampleNumber==2722 ){	samplename = "zzz";					}
		else if( sampleNumber==2524 ){	samplename = "ttbarW";				}
		else if( sampleNumber==2534 ){	samplename = "ttbarWW";				}
		else if( sampleNumber==2523 ){	samplename = "ttbarZ";				}
		else if( sampleNumber==2525 ){	samplename = "ttbarttbar";			}
		else if( sampleNumber>=9100 && sampleNumber<=9300 ){
				 if (era == "2012_53x"){	samplename = "ttH120_FullSim";	}
			else if (era == "2012_52x"){	samplename = "ttH120_FastSim";	}
			else{ cout << "Era set to '" << era << "'" << endl; assert (era == "either 2012_52x or 2012_53x"); }
		}else if( sampleNumber>=8100 && sampleNumber<=8300 ){
				 if (era == "2012_53x"){	samplename = "ttH120_bb";		}
			else if (era == "2012_52x"){	samplename = "ttH120_FullSim";	}
			else{	cout << "Era set to '" << era << "'" << endl; assert (era == "either 2012_52x or 2012_53x");	}
		}else if ( sampleNumber>=7100 && sampleNumber<=7300 && era == "2012_53x" ){
			samplename = "ttH120_FullSim";
		}else{ cout << "Sample number set to '" << sampleNumber << "'" << endl; assert (samplename == "No good sampleNumber found for 2012_52x or 2012_53x"); }

	}else if (era == "2011") {

			 if( sampleNumber<0 ){							samplename = "data";			}
		else if( sampleNumber==2300 ){						samplename = "zjets";			}
		else if( sampleNumber==2310 ){						samplename = "zjets_lowmass";	}
		else if( sampleNumber==2400 ){						samplename = "wjets";			}
		else if( sampleNumber==2500 ){						samplename = "ttbar";			}
		else if( sampleNumber==2544 ){						samplename = "ttbar_cc";		}
		else if( sampleNumber==2555 ){						samplename = "ttbar_bb";		}
		else if( sampleNumber==2510 ){						samplename = "ttbar_scaleup";	}
		else if( sampleNumber==2511 ){						samplename = "ttbar_scaledown";	}
		else if( sampleNumber==2523 ){						samplename = "ttbarZ";			}
		else if( sampleNumber==2524 ){						samplename = "ttbarW";			}
		else if( sampleNumber==2600 ){						samplename = "singlet";			}
		else if( sampleNumber==2700 ){						samplename = "ww";				}
		else if( sampleNumber==2701 ){						samplename = "wz";				}
		else if( sampleNumber==2702 ){						samplename = "zz";				}
		else if( sampleNumber>=100 && sampleNumber<=300 ){	samplename = "ttH120";			}
		else{  assert (samplename == "No good sampleNumber found for 2011");				}

	}else{
		cout << "Era set to '" << era << "'" << endl;
		assert (era == "either 2012_52x, 2012_53x, or 2011");
	}

	return samplename;
}

// Check that we are set up, otherwise inform and quit
void BEANhelper::CheckSetUp(){ if(!isSetUp){ ThrowFatalError("BEANhelper not yet set up."); } }

// If something goes really wrong, inform and quit
void BEANhelper::ThrowFatalError(string const iMessage){ cerr << "[ERROR]\t" << iMessage << " Cannot continue. Terminating..." << endl; exit(1); }

// Return corrected MET based on *UNCORRECTED* input jet collection
BNmet BEANhelper::GetCorrectedMET(const BNmet& iMET, const BNjetCollection& iJets, const sysType::sysType iSysType){
	CheckSetUp();
	// IMPORTANT! iJets is the *UNCORRECTED* jet collection upon which to base the MET correction

	TLorentzVector newMET(iMET.px, iMET.py, 0, 0);
    float minPt = -14000.0; //14 TeV
    
	// Correct the MET based on jet pT corrections
	for( BNjetCollection::const_iterator Jet = iJets.begin(); Jet != iJets.end(); ++Jet ){

      
		// Get the corrected jet
		BNjet correctedJet = GetCorrectedJet(*Jet, iSysType);

        // Protect against NAN jets
        if (correctedJet.px > minPt && Jet->px > minPt && correctedJet.py > minPt && Jet->py > minPt) {
          // Compute the correction (difference)
          TLorentzVector jetCorrection((correctedJet.px-Jet->px), (correctedJet.py-Jet->py));

          // Subtract difference from MET
          newMET -= jetCorrection;
        }
	}

	// Make a copy of the input MET for output and update xy values
	BNmet result = iMET;
	result.px = newMET.Px();
	result.py = newMET.Py();
	result.pt = newMET.Pt();

	return result;

}

// Return corrected jet
BNjet BEANhelper::GetCorrectedJet(const BNjet& iJet, const sysType::sysType iSysType){
	CheckSetUp();

	// Return original jet if this sample contains collision data
	if(isData){ return iJet; }

	// Factor by which we will scale all the 4-vector components of the jet
	double jetFactor = 1.;

	// Make appopriate correction to the xy components of input jet
	switch(iSysType){
		case sysType::JERup:	jetFactor *= getJERfactor(1, fabs(iJet.eta), iJet.genJetPT, iJet.pt);	 break;
		case sysType::JERdown:	jetFactor *= getJERfactor(-1, fabs(iJet.eta), iJet.genJetPT, iJet.pt);	 break;
		case sysType::JESup:	jetFactor *= (1. + iJet.JESunc);
                                jetFactor *= getJERfactor(0, fabs(iJet.eta), iJet.genJetPT, jetFactor * iJet.pt); break;
		case sysType::JESdown:	jetFactor *= (1. - iJet.JESunc);
                                jetFactor *= getJERfactor(0, fabs(iJet.eta), iJet.genJetPT, jetFactor * iJet.pt); break;
		default:				jetFactor *= getJERfactor(0, fabs(iJet.eta), iJet.genJetPT, iJet.pt);	 break;
	}

	// Make a copy of the input jet for output and update 4-vector values
	BNjet result = iJet;
    
	result.px	    *= jetFactor;
	result.py	    *= jetFactor;
	result.pz	    *= jetFactor;
	result.pt	    *= jetFactor;
    result.et       *= jetFactor;
	result.energy	*= jetFactor;
    result.mass     *= jetFactor;

    // Update CSV value (i.e. reshape it if applicable)
// 	result.btagCombinedSecVertex = GetCSVvalue(result,iSysType);

	return result;
}

// Return collection with corrected jets
BNjetCollection BEANhelper::GetCorrectedJets(const BNjetCollection& iJets, const sysType::sysType iSysType){
	CheckSetUp();
	BNjetCollection result;		
	for( BNjetCollection::const_iterator Jet = iJets.begin(); Jet != iJets.end(); ++Jet ){
		result.push_back(GetCorrectedJet(*Jet, iSysType));
	}
	return result;
}

BNjetCollection BEANhelper::GetUncorrectedJets(const BNjetCollection& iCorrectedJets, const BNjetCollection& iOriginalJets){
	CheckSetUp();

	if(iOriginalJets.size() < iCorrectedJets.size()){ cerr << "ERROR in BEANhelper::GetUncorrectedJets: size of 'iCorrectedJets' is greater than 'iOriginalJets'. Terminating..." << endl; exit(1); }
	BNjetCollection result;

	for( BNjetCollection::const_iterator cJet = iCorrectedJets.begin(); cJet != iCorrectedJets.end(); ++cJet ){
		bool found = false;
		for( BNjetCollection::const_iterator oJet = iOriginalJets.begin(); oJet != iOriginalJets.end(); ++oJet ){
			if(deltaR(cJet->eta, cJet->phi, oJet->eta, oJet->phi) < 0.001){ result.push_back(*oJet); found = true; break; }
		}
		if(!found){ cerr << "ERROR in BEANhelper::GetUncorrectedJets: could not match one of the jets. Check input jet collections. Terminating..." << endl; exit(1); }
	}

	return result;
}

// Return the btag weight for the input jet (product of efficiency and scale factor)
float BEANhelper::GetBtagWeight(const BNjet& iJet, const sysType::sysType iSysType){
	CheckSetUp();

	// Do nothing if this sample contains collision data
	if(isData){ return 1.0; }

	std::vector<double> myEffSF;
	switch(iSysType){
		case sysType::hfSFup:	myEffSF = getEffSF(1, iJet.pt, iJet.eta, iJet.flavour);		break;
		case sysType::hfSFdown:	myEffSF = getEffSF(-1, iJet.pt, iJet.eta, iJet.flavour);	break;
		case sysType::lfSFup:	myEffSF = getEffSF(2, iJet.pt, iJet.eta, iJet.flavour);		break;
		case sysType::lfSFdown:	myEffSF = getEffSF(-2, iJet.pt, iJet.eta, iJet.flavour);	break;
		default:				myEffSF = getEffSF(0, iJet.pt, iJet.eta, iJet.flavour);		break;

	}

	return ((myEffSF[0])*(myEffSF[1]));
}

bool BEANhelper::PassesCSV(const BNjet& iJet, const char iCSVworkingPoint){
	CheckSetUp();

	float csvValue = iJet.btagCombinedSecVertex;

	// CSV b-tagging requirement
	switch(iCSVworkingPoint){
		case 'L':	if(csvValue > CSVLwp){ return true; }	break;
		case 'M':	if(csvValue > CSVMwp){ return true; }	break;
		case 'T':	if(csvValue > CSVTwp){ return true; }	break;
		case '-':	return true;							break;
	}
	return false;
}

double BEANhelper::GetCSVvalue(const BNjet& iJet, const sysType::sysType iSysType){
	CheckSetUp();
	double result = iJet.btagCombinedSecVertex;

	if(isData){	return result; }

	if(reshapeCSV){
		switch(iSysType){
            case sysType::hfSFup:	result = sh_hfSFUp_->GetReshapedCSVvalue(iJet.eta, iJet.pt, iJet.btagCombinedSecVertex, iJet.flavour);      break;
			case sysType::hfSFdown:	result = sh_hfSFDown_->GetReshapedCSVvalue(iJet.eta, iJet.pt, iJet.btagCombinedSecVertex, iJet.flavour);	break;
			case sysType::lfSFup:	result = sh_lfSFUp_->GetReshapedCSVvalue(iJet.eta, iJet.pt, iJet.btagCombinedSecVertex, iJet.flavour);		break;
			case sysType::lfSFdown:	result = sh_lfSFDown_->GetReshapedCSVvalue(iJet.eta, iJet.pt, iJet.btagCombinedSecVertex, iJet.flavour);	break;
            default:				result = sh_->GetReshapedCSVvalue(iJet.eta, iJet.pt, iJet.btagCombinedSecVertex, iJet.flavour);				break;
		}	
	}
    
	return result;
}

// Return whether or not jet passes cuts
bool BEANhelper::IsGoodJet(const BNjet& iJet, const float iMinPt, const float iMaxAbsEta, const jetID::jetID iJetID, const char iCSVworkingPoint){
	CheckSetUp();
	// Transverse momentum requirement
	if(iJet.pt < iMinPt){ return false; }

	// Absolute eta requirement
	if(fabs(iJet.eta) > iMaxAbsEta){ return false; }

	// Jet ID
	switch(iJetID){
		case jetID::jetMinimal:		if(!iJet.jetIDMinimal){ return false; }	break;
		case jetID::jetLooseAOD:	if(!iJet.jetIDLooseAOD){ return false; }	break;
		case jetID::jetLoose:		if(!iJet.jetIDLoose){ return false; }	break;
		case jetID::jetTight:		if(!iJet.jetIDTight){ return false; }	break;
		default:					return false;	break;

	}

	if(!PassesCSV(iJet, iCSVworkingPoint)){ return false; }

	return true;
}

// Return collection with objects passing cuts
BNjetCollection BEANhelper::GetSelectedJets(const BNjetCollection& iJets, const float iMinPt, const float iMaxAbsEta, const jetID::jetID iJetID, const char iCSVwp){
	BNjetCollection result;		
	for( BNjetCollection::const_iterator Jet = iJets.begin(); Jet != iJets.end(); ++Jet ){
		if(IsGoodJet((*Jet), iMinPt, iMaxAbsEta, iJetID, iCSVwp)){ result.push_back(*Jet); }
	}
	return result;
}

float BEANhelper::GetHT(const BNjetCollection& iJets){
	CheckSetUp();
	float result = 0;
	for( BNjetCollection::const_iterator Jet = iJets.begin(); Jet != iJets.end(); ++Jet ){ result += Jet->pt; }
	return result;
}


// Muon relative isolation
float BEANhelper::GetMuonRelIso(const BNmuon& iMuon){ 
	CheckSetUp();
	float result = 9999;
    if (era=="2011") {
      if( usePfLeptons )
        result = ((iMuon.chargedHadronIso + iMuon.neutralHadronIso + iMuon.photonIso)/iMuon.pt);
      else
        result = ((iMuon.trackIso + iMuon.ecalIso + iMuon.hcalIso)/iMuon.pt);
    }
	else if (era=="2012_52x" || era=="2012_53x") {
      if( usePfLeptons )
        result = (((iMuon.pfIsoR04SumChargedHadronPt + max(0.0, iMuon.pfIsoR04SumNeutralHadronEt + iMuon.pfIsoR04SumPhotonEt - 0.5*iMuon.pfIsoR04SumPUPt)))/iMuon.pt);
      else
        result = (((iMuon.pfIsoR04SumChargedHadronPt + max(0.0, iMuon.pfIsoR04SumNeutralHadronEt + iMuon.pfIsoR04SumPhotonEt - 0.5*iMuon.pfIsoR04SumPUPt)))/iMuon.pt);
    }
    else {
      assert (era == "either 2012_52x, 2012_53x, or 2011");
    }

	return result;
}

float BEANhelper::GetMuonSF( const BNmuon& iMuon, const muonID::muonID inputID ){
  CheckSetUp();
  double SF = 1.0;

  if (isData) return SF;

  double usePT = -10;
  double useEta = -10;
  
  switch (inputID) {
    // Tight selection
  case muonID::muonTight:
    usePT = std::min( iMuon.pt, 499. );
    useEta = std::min( fabs(iMuon.eta), 2.09);  
    SF = h_mu_SF_->GetBinContent( h_mu_SF_->FindBin(usePT, useEta) );
    break;
    
    // Treat all other selections as loose
    // Loose SF is closest for these in any case.
  case muonID::muonLoose:
  case muonID::muonSide:
  case muonID::muonPtOnly:
  case muonID::muonPtEtaOnly:
  case muonID::muonPtEtaIsoOnly:
  case muonID::muonPtEtaIsoTrackerOnly:
    
    usePT = std::min( iMuon.pt, 200. );
    useEta = std::min( fabs(iMuon.eta), 2.09);
    // for this histo x = eta y = pt
    SF = h_looseMuonSF->GetBinContent( h_looseMuonSF->FindBin(useEta, usePT) );
    break;
    
  }

  if (SF < 0.005) {
    TString message = Form("Error! Muon SF too small (%f).\nCalled GetMuonSF with pt = %f eta = %f and tight = %d ",
                           SF, iMuon.pt, iMuon.eta, inputID==muonID::muonTight);
    ThrowFatalError(message.Data());

  }

  return SF;
}
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Double Lepton Trigger SF ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// double muon
float BEANhelper::GetDoubleMuonTriggerSF ( const BNmuon & firstMuon, const BNmuon & secondMuon) {

  CheckSetUp();
  if (era != "2012_53x") ThrowFatalError("We can only do trigger SF in 2012_53x. Please use that era.");

  float SF = 1.0;

  if (isData) return SF;

  float useLep1Eta = std::min (fabs(firstMuon.eta), 2.09);
  float useLep2Eta = std::min (fabs(secondMuon.eta), 2.49);

  SF = h_doubleMuTrigSF->GetBinContent(h_doubleMuTrigSF->FindBin(useLep1Eta, useLep2Eta));
    
  
  return SF;

}

// double electron
float BEANhelper::GetDoubleElectronTriggerSF ( const BNelectron & firstEle, const BNelectron & secondEle) {

  CheckSetUp();
  if (era != "2012_53x") ThrowFatalError("We can only do trigger SF in 2012_53x. Please use that era.");
  
  float SF = 1.0;

  if (isData) return SF;

  float useLep1Eta = std::min (fabs(firstEle.eta), 2.09);
  float useLep2Eta = std::min (fabs(secondEle.eta), 2.49);

  SF = h_doubleEleTrigSF->GetBinContent(h_doubleEleTrigSF->FindBin(useLep1Eta, useLep2Eta));
    
  
  return SF;

}


// muon electron
float BEANhelper::GetMuonEleTriggerSF ( const BNmuon & firstMuon, const BNelectron & secondEle) {

  CheckSetUp();
  if (era != "2012_53x") ThrowFatalError("We can only do trigger SF in 2012_53x. Please use that era.");
  float SF = 1.0;

  if (isData) return SF;

  float useLep1Eta = std::min (fabs(firstMuon.eta), 2.49);
  float useLep2Eta = std::min (fabs(secondEle.eta), 2.49);

  SF = h_muonEleTrigSF->GetBinContent(h_muonEleTrigSF->FindBin(useLep1Eta, useLep2Eta));
    
  
  return SF;

}

// test only
float BEANhelper::TestSingleMuonTriggerNew ( const BNmuon & iMuon) {
  //cout << "Inside Test single muon trigger new" << endl;
  CheckSetUp();
  if (era != "2012_53x") ThrowFatalError("We can only do trigger SF in 2012_53x. Please use that era.");
  float SF = 1.0;

  if (isData) return SF;

  
  double usePt = std::min(iMuon.pt, 499.0);
  usePt = std::max (27.0, usePt);

  float useEta = std::min(fabs(iMuon.eta), 2.05) ;

  SF =  h_testSingleMuNew->GetBinContent(h_testSingleMuNew->FindBin(usePt, useEta));
  //cout << "returning..." << endl;
  return SF;

}

float BEANhelper::TestSingleEleTriggerNew ( const BNelectron & iEle) {
  //cout << "Inside Test single ele trigger new" << endl;
  CheckSetUp();
  if (era != "2012_53x") ThrowFatalError("We can only do trigger SF in 2012_53x. Please use that era.");
  float SF = 1.0;

  if (isData) return SF;

  double usePt = std::min(iEle.pt, 499.0);
  usePt = std::max (30.0, usePt);

  // eventually, but not now
  //float useEta = fabs (std::min(iEle.eta, 2.48)) ;
  //float useEta = std::min(iEle.eta, 2.48) ;
  float useEta = iEle.eta;

  SF = h_testSingleEleNew->GetBinContent(h_testSingleEleNew->FindBin(usePt, useEta));
  //cout << "returning..." << endl;
  return SF;

}

float BEANhelper::TestSingleMuonTriggerOld ( const BNmuon & iMuon) {
  //cout << "Inside Test single muon trigger old" << endl;
  CheckSetUp();
  if (era != "2012_53x") ThrowFatalError("We can only do trigger SF in 2012_53x. Please use that era.");
  float SF = 1.0;

  if (isData) return SF;

  
  double usePt = std::min(iMuon.pt, 499.0);
  usePt = std::max (27.0, usePt);

  float useEta = iMuon.eta;

  SF =  h_testSingleMuOld->GetBinContent(h_testSingleMuOld->FindBin(usePt, useEta));
  //cout << "returning..." << endl;
  return SF;

}

float BEANhelper::TestSingleEleTriggerOld ( const BNelectron & iEle) {
  //cout << "Inside Test single ele trigger old" << endl;
  CheckSetUp();
  if (era != "2012_53x") ThrowFatalError("We can only do trigger SF in 2012_53x. Please use that era.");
  float SF = 1.0;

  if (isData) return SF;

  double usePt = std::min(iEle.pt, 499.0);
  usePt = std::max (30.0, usePt);

  float useEta = std::min(iEle.eta, 2.48);

  SF = h_testSingleEleOld->GetBinContent(h_testSingleEleOld->FindBin(usePt, useEta));
  
  //cout << "returning..." << endl;
  return SF;

}


/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
    

// Return whether or not tau passes cuts
bool BEANhelper::IsVLooseTau(const BNtau& iTau){ return IsGoodTau(iTau, tauID::tauVLoose); }
bool BEANhelper::IsLooseTau(const BNtau& iTau){ return IsGoodTau(iTau, tauID::tauLoose); }
bool BEANhelper::IsMediumTau(const BNtau& iTau){ return IsGoodTau(iTau, tauID::tauLoose); }
bool BEANhelper::IsTightTau(const BNtau& iTau){ return IsGoodTau(iTau, tauID::tauTight); }

bool BEANhelper::IsGoodTau(const BNtau& iTau, const tauID::tauID iTauID){
	CheckSetUp();

	// Be skeptical about this muon making it through
	bool passesKinematics	= false;
	bool passesIso			= false;
	bool passesID			= false;

	// Check if this muon is good enough
	switch(iTauID){
		case tauID::tauNonIso:
			passesKinematics		= (iTau.pt >= 20) && (fabs(iTau.eta) <= 2.1);
			passesID				= (iTau.leadingTrackPt >= 5) && (iTau.HPSdecayModeFinding > 0) && (iTau.HPSagainstElectronLoose > 0) && (iTau.HPSagainstMuonLoose > 0);
			passesIso				= true;
			break;
		case tauID::tauVLoose:
			passesKinematics		= (iTau.pt >= 20) && (fabs(iTau.eta) <= 2.1);
			passesID				= (iTau.leadingTrackPt >= 5) && (iTau.HPSdecayModeFinding > 0) && (iTau.HPSagainstElectronLoose > 0) && (iTau.HPSagainstMuonLoose > 0);
			passesIso				= (iTau.HPSbyVLooseCombinedIsolationDeltaBetaCorr > 0);
			break;
		case tauID::tauLoose:
			passesKinematics		= (iTau.pt >= 20) && (fabs(iTau.eta) <= 2.1);
			passesID				= (iTau.leadingTrackPt >= 5) && (iTau.HPSdecayModeFinding > 0) && (iTau.HPSagainstElectronLoose > 0) && (iTau.HPSagainstMuonLoose > 0);
			passesIso				= (iTau.HPSbyLooseCombinedIsolationDeltaBetaCorr > 0);
			break;
		case tauID::tauMedium:
			passesKinematics		= (iTau.pt >= 20) && (fabs(iTau.eta) <= 2.1);
			passesID				= (iTau.leadingTrackPt >= 5) && (iTau.HPSdecayModeFinding > 0) && (iTau.HPSagainstElectronTight > 0) && (iTau.HPSagainstMuonTight > 0);
			passesIso				= (iTau.HPSbyMediumCombinedIsolationDeltaBetaCorr > 0);
			break;
		case tauID::tauTight:
			passesKinematics		= (iTau.pt >= 20) && (fabs(iTau.eta) <= 2.1);
			passesID				= (iTau.leadingTrackPt >= 5) && (iTau.HPSdecayModeFinding > 0) && (iTau.HPSagainstElectronTight > 0) && (iTau.HPSagainstMuonTight > 0);
			passesIso				= (iTau.HPSbyTightCombinedIsolationDeltaBetaCorr > 0);
			break;
	}

	return (passesKinematics && passesID && passesIso);

}

// Return collection with objects passing cuts
BNtauCollection BEANhelper::GetSelectedTaus(const BNtauCollection& iTaus, const tauID::tauID iTauID){
	CheckSetUp();
	BNtauCollection result;
	for( BNtauCollection::const_iterator Tau = iTaus.begin(); Tau != iTaus.end(); ++Tau ){
		if(IsGoodTau((*Tau), iTauID)){ result.push_back(*Tau); }
	}
	return result;
}

BNtau BEANhelper::GetCorrectedTau(const BNtau& iTau, const sysType::sysType iSysType){
	CheckSetUp();

	// Return original tau if this sample contains collision data
	if(isData){ return iTau; }

	// Factor by which we will scale all the 4-vector components of the tau
	double tauFactor = 1.;

	// Make appopriate correction to the components of input tau
	switch(iSysType){
		case sysType::TESup:	tauFactor *= (1 + 0.03); break;
		case sysType::TESdown:	tauFactor *= (1 - 0.03); break;
		default:				tauFactor *= 1;	 break;
	}

	// Make a copy of the input tau for output and update 4-vector values
	BNtau result = iTau;
    
	result.px                *= tauFactor;
	result.py                *= tauFactor;
	result.pz                *= tauFactor;
	result.pt                *= tauFactor;
	result.et                *= tauFactor;
	result.energy            *= tauFactor;
	result.leadingTrackPt    *= tauFactor;
	//result.mass            *= tauFactor;

	return result;
}

// Return collection with corrected taus
BNtauCollection BEANhelper::GetCorrectedTaus(const BNtauCollection& iTaus, const sysType::sysType iSysType){
	CheckSetUp();
	BNtauCollection result;		
	for( BNtauCollection::const_iterator Tau = iTaus.begin(); Tau != iTaus.end(); ++Tau ){
		result.push_back(GetCorrectedTau(*Tau, iSysType));
	}
	return result;
}

// Return whether or not muon passes cuts
bool BEANhelper::IsSideMuon(const BNmuon& iMuon){ return IsGoodMuon(iMuon, muonID::muonSide); }

bool BEANhelper::IsLooseMuon(const BNmuon& iMuon){ return IsGoodMuon(iMuon, muonID::muonLoose); }

bool BEANhelper::IsTightMuon(const BNmuon& iMuon){ return IsGoodMuon(iMuon, muonID::muonTight); }

bool BEANhelper::IsGoodMuon(const BNmuon& iMuon, const muonID::muonID iMuonID){
	CheckSetUp();

	// Set default kinematic thresholds
	float minLooseMuonPt		= 9999;
	float minTightMuonPt		= 9999;
	float maxLooseMuonAbsEta	= 0;
	float maxTightMuonAbsEta	= 0;

    if (era=="2011") {
      minLooseMuonPt		= 10;
      minTightMuonPt		= ( isLJ ) ? 30. : 20.;
      maxLooseMuonAbsEta	= 2.4;
      maxTightMuonAbsEta	= 2.1;
    }
    else if (era=="2012_52x" || era=="2012_53x") {
      minLooseMuonPt		= 10;
      minTightMuonPt		= ( isLJ ) ? 30. : 20.;
      maxLooseMuonAbsEta	= 2.5;
      maxTightMuonAbsEta	= 2.1;
    }
    else {
      assert (era == "either 2012_52x, 2012_53x, or 2011");
    }

	// Be skeptical about this muon making it through
	bool passesKinematics	= false;
	bool passesIso			= false;
	bool passesID			= false;
	bool isPFMuon			= false;
	bool passesTrackerID    = false;

	// Check if this muon is good enough
    if (era=="2011") { 
      switch(iMuonID){
      case muonID::muonSide:
        passesKinematics		= ((iMuon.pt >= minLooseMuonPt) && (fabs(iMuon.eta) <= maxLooseMuonAbsEta));
        passesIso				= (GetMuonRelIso(iMuon) < 0.800);
        passesID				= (iMuon.isGlobalMuon==1);
        break;
      case muonID::muonLoose:
        passesKinematics		= ((iMuon.pt >= minLooseMuonPt) && (fabs(iMuon.eta) <= maxLooseMuonAbsEta));
        passesIso				= (GetMuonRelIso(iMuon) < 0.200);
        passesID				= (iMuon.isGlobalMuon==1);
        break;
      case muonID::muonTight:
        passesKinematics		= ((iMuon.pt >= minTightMuonPt) && (fabs(iMuon.eta) <= maxTightMuonAbsEta));
        passesIso				= (GetMuonRelIso(iMuon) < 0.125);
        passesTrackerID	        = ((iMuon.isTrackerMuon) && (iMuon.numberOfValidTrackerHitsInnerTrack > 10) && (iMuon.pixelLayersWithMeasurement > 0)
                                   && (iMuon.numberOfMatchedStations > 1) && (fabs(iMuon.correctedD0) < 0.02) && (fabs(iMuon.correctedDZ) < 1.));
        passesID				= ((iMuon.isGlobalMuon==1) && (iMuon.isGlobalMuonPromptTight==1) && passesTrackerID);
        break;
      case muonID::muonPtOnly:
      case muonID::muonPtEtaOnly:
      case muonID::muonPtEtaIsoOnly:
      case muonID::muonPtEtaIsoTrackerOnly:
        cout << "Oops, you asked for loose muon selections in 2011. These aren't defined yet" << endl
             << "Sorry... we have to crash now" << endl;
        ThrowFatalError("Asked for a loose muon selection in 2011, which is not possible.");
        break;
        
      }
    }// End of 2011 era
    else if (era=="2012_52x" || era=="2012_53x") {
      switch(iMuonID){
      case muonID::muonSide:
        passesKinematics		= ((iMuon.pt >= minLooseMuonPt) && (fabs(iMuon.eta) <= maxLooseMuonAbsEta));
        passesIso				= (GetMuonRelIso(iMuon) < 0.800);
        isPFMuon				= true;
        passesID				= (((iMuon.isGlobalMuon==1) || (iMuon.isTrackerMuon==1)) && isPFMuon);
        break;
      case muonID::muonLoose:
        passesKinematics		= ((iMuon.pt >= minLooseMuonPt) && (fabs(iMuon.eta) <= maxLooseMuonAbsEta));
        passesIso				= (GetMuonRelIso(iMuon) < 0.200);
        isPFMuon				= true;
        passesID				= (((iMuon.isGlobalMuon==1) || (iMuon.isTrackerMuon==1)) && isPFMuon);
        break;
      case muonID::muonTight:
        passesKinematics		= ((iMuon.pt >= minTightMuonPt) && (fabs(iMuon.eta) <= maxTightMuonAbsEta));
        passesIso				= (GetMuonRelIso(iMuon) < 0.120);
        isPFMuon				= true;
        passesTrackerID	        = ((iMuon.isGlobalMuon==1)
                                   && (iMuon.normalizedChi2 < 10) && (fabs(iMuon.correctedD0Vertex) < 0.2) && (fabs(iMuon.dVzPVz) < 0.5) 
                                   && (iMuon.numberOfLayersWithMeasurement > 5 ) && (iMuon.numberOfValidMuonHits > 0)
                                   && (iMuon.numberOfValidPixelHits > 0) && (iMuon.numberOfMatchedStations > 1));
        
        passesID				= (((iMuon.isGlobalMuon==1) || (iMuon.isTrackerMuon==1)) && isPFMuon && passesTrackerID);
        break;
      case muonID::muonPtOnly:
        passesKinematics		= ((iMuon.pt >= minTightMuonPt));
        passesIso				= true;
        isPFMuon				= true;
        passesTrackerID	        = true;
        
        passesID				= true;
        break;
      case muonID::muonPtEtaOnly:
        passesKinematics		= ((iMuon.pt >= minTightMuonPt) && (fabs(iMuon.eta) <= maxTightMuonAbsEta));
        passesIso				= true;
        isPFMuon				= true;
        passesTrackerID	        = true;
        
        passesID				= true;
        break;
      case muonID::muonPtEtaIsoOnly:
        passesKinematics		= ((iMuon.pt >= minTightMuonPt) && (fabs(iMuon.eta) <= maxTightMuonAbsEta));
        passesIso				= (GetMuonRelIso(iMuon) < 0.120);
        isPFMuon				= true;
        passesTrackerID	        = true;
        
        passesID				= true;
        break;
      case muonID::muonPtEtaIsoTrackerOnly:
        passesKinematics		= ((iMuon.pt >= minTightMuonPt) && (fabs(iMuon.eta) <= maxTightMuonAbsEta));
        passesIso				= (GetMuonRelIso(iMuon) < 0.120);
        isPFMuon				= true;
        passesTrackerID	        = ((iMuon.isGlobalMuon==1)
                                   && (iMuon.normalizedChi2 < 10) && (fabs(iMuon.correctedD0Vertex) < 0.2) && (fabs(iMuon.dVzPVz) < 0.5) 
                                   && (iMuon.numberOfLayersWithMeasurement > 5 ) && (iMuon.numberOfValidMuonHits > 0)
                                   && (iMuon.numberOfValidPixelHits > 0) && (iMuon.numberOfMatchedStations > 1));
        
        passesID				= passesTrackerID;
        break;
      }
	}// End of 2012 era
    else {
      assert (era == "either 2012_52x, 2012_53x, or 2011");
    }

	return (passesKinematics && passesIso && passesID);

}


bool BEANhelper::IsAnyTriggerBitFired ( const vector<string> & targetTriggers, const BNtriggerCollection &  triggerBits) {

    // check to see if you passed the trigger by looping over the bits
  // looking for your bit
  // and see if it is 1
  bool passTrigger = false;

  
  for ( BNtriggerCollection::const_iterator iBit = triggerBits.begin();
        iBit != triggerBits.end();
        iBit ++ ){

    for (vector<string>::const_iterator iTarget = targetTriggers.begin();
         iTarget != targetTriggers.end();
         iTarget++) {

      // if this is the right name and the bit is set to one
      if ( iBit->name.find((*iTarget)) != std::string::npos && iBit->pass==1) {
        passTrigger = true;
      }
      
    }// end for each target

  }// end for each bit

  return passTrigger;
  
}

bool BEANhelper::SingleObjectMatchesAnyTrigger (double recoEta, double recoPhi, const vector<string> & targetTriggers,  const BNtrigobjCollection & triggerObjects) {

  bool matchTrigger = false;
  
  for ( BNtrigobjCollection::const_iterator iObj = triggerObjects.begin();
        iObj != triggerObjects.end();
        iObj ++ ) {
    
    // look for a matching name again
    for (vector<string>::const_iterator iTarget = targetTriggers.begin();
         iTarget != targetTriggers.end();
         iTarget++) {
      
      // if this is the right name 
      if ( iObj->filter.find((*iTarget)) != std::string::npos) {
        double deltaR = sqrt ( ( recoEta - iObj->eta ) * ( recoEta - iObj->eta )
                               + (recoPhi - iObj->phi) * ( recoPhi - iObj->phi)
                               );

        if (deltaR < 0.3)
          matchTrigger = true;
          
      }
      
    }// end for each target

    
  } // end for each object


  return matchTrigger;
  
  
}


bool BEANhelper::DoubleObjectMatchesAnyTrigger (double firstEta, double firstPhi, double secondEta, double secondPhi,
                                                const vector<string> & targetTriggers,  const BNtrigobjCollection & triggerObjects,
                                                bool sumFilterResults) {

  if (triggerDebug) cout << "DoubleObjectMatchesAnyTrigger:Looking for match to\nDoubleObjectMatchesAnyTrigger: input 1: ("
                         << firstEta << ", " << firstPhi
                         << ")\nDoubleObjectMatchesAnyTrigger: input 2: ("
                         << secondEta << ", " << secondPhi << ")"
                         << endl;

  // This is trickier than the other case
  // Need to make sure that you match both objects to the same name
  // of trigger

  // Define some typedefs so that you
  // can write shorter code
  // 
  // typedef pair< bool, bool> pairOfBools;
  typedef pair<string, int> trigNameAndMatchedLegs;
  typedef vector<trigNameAndMatchedLegs> trigNameAndMatchedLegsCollection;

  trigNameAndMatchedLegsCollection myMatches;
  trigNameAndMatchedLegsCollection myCandidates;

  for ( unsigned iName = 0; iName < targetTriggers.size(); iName++){
    myMatches.push_back( trigNameAndMatchedLegs ( targetTriggers[iName], 0 ));
    myCandidates.push_back( trigNameAndMatchedLegs ( targetTriggers[iName], 0 ));    
  }
  

  /////////////////////////////////////////////////////////////////
  
  
  for ( BNtrigobjCollection::const_iterator iObj = triggerObjects.begin();
        iObj != triggerObjects.end();
        iObj ++ ) {

    // this is very verbose
    //if (triggerDebug) {
    //  cout << "DoubleObjectMatchesAnyTrigger: Found an object for filter " << iObj->filter << endl;
    //}
    
    // look for a matching name again
    unsigned numFilters = 0;
    for (trigNameAndMatchedLegsCollection::iterator iTarget = myMatches.begin();
         iTarget != myMatches.end();
         iTarget++) {
      
      // if this is the right name 
      if ( iObj->filter.find((iTarget->first)) != std::string::npos) {
        myCandidates[numFilters].second++;
        if (triggerDebug) {
          cout << "DoubleObjectMatchesAnyTrigger: Found an object for filter "
               << iTarget->first << endl
               << "eta = " << iObj->eta << " phi = " << iObj->phi
               << endl;          
        }
        
        double firstDR = sqrt ( ( firstEta - iObj->eta ) * ( firstEta - iObj->eta )
                               + (firstPhi - iObj->phi) * ( firstPhi - iObj->phi)
                               );

        if (firstDR < 0.3){
          if (triggerDebug)
            cout << "DoubleObjectMatchesAnyTrigger: Found a match to first object with DR = " << firstDR << endl;
          
          iTarget->second++;
        }

        double secondDR = sqrt ( ( secondEta - iObj->eta ) * ( secondEta - iObj->eta )
                               + (secondPhi - iObj->phi) * ( secondPhi - iObj->phi)
                               );

        if (secondDR < 0.3){
          if (triggerDebug)
            cout << "DoubleObjectMatchesAnyTrigger: Found a match to second object with DR = " << secondDR << endl;
          
          iTarget->second++;
        }
          
      }
      numFilters++;
    }// end for each target

    
  } // end for each object


  bool matchedBothLegsToOneTrigger = false;
  int totalNumFiltersFired =0;
  // look at all the triggers whose matches you counted
  // and see if any had >=2 matches.
  // It is possible to have more than two matches because the trigger selects on >=2 matches
  // In the following cases, there are could be more than 2
  //   1. There are more than two trigger muons in the event
  //   2. The two muons from the trigger are within a cone of 0.3, in which
  //      case both muons match both legs
  //
  if (triggerDebug)
    cout << "Done looping over objects, printing results" << endl;

  int numFilters =0;
  for (trigNameAndMatchedLegsCollection::const_iterator iTarget = myMatches.begin();
       iTarget != myMatches.end();
       iTarget++) {

    if (triggerDebug) {
      cout << "Filter  " << iTarget->first << " matches = "
           << iTarget->second << endl
           << "Candidates for filter " << myCandidates[numFilters].first
           << " numCands = " << myCandidates[numFilters].second
           << endl;

    }

    totalNumFiltersFired += iTarget->second;
    
    if (iTarget->second >=2 )
      matchedBothLegsToOneTrigger = true;

    numFilters++;
  }

  // this is false by default
  // only necessary if you are doing MuEG triggers
  if (sumFilterResults) {
    matchedBothLegsToOneTrigger = (totalNumFiltersFired >=2);
  }

  return matchedBothLegsToOneTrigger;
  
  
}




bool BEANhelper::MuonMatchesSingleMuTrigger(const BNmuon& iMuon, const BNtriggerCollection & triggerBits, const BNtrigobjCollection & triggerObjects ){

  // Define your target triggers
  vector<string> singleMuNames;
  singleMuNames.push_back("HLT_IsoMu24_eta2p1_v");
  singleMuNames.push_back("HLT_IsoMu24_v");


  // Now check to see if there is an HLT object
  // that corresponds to your trigger

  bool bitFired = IsAnyTriggerBitFired ( singleMuNames, triggerBits);

  if (!bitFired) return false;

  vector<string> singleMuFilters;
  singleMuFilters.push_back("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15::HLT");
  singleMuFilters.push_back("hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15::HLT");
  
  bool matchTrigger = SingleObjectMatchesAnyTrigger(iMuon.eta, iMuon.phi, singleMuFilters, triggerObjects);
  
  return matchTrigger;
  
}

bool BEANhelper::MuonsMatchDoubleMuTrigger(const BNmuon& iMuon, const BNmuon & jMuon, const BNtriggerCollection & triggerBits, const BNtrigobjCollection & triggerObjects){


  // Define your target triggers
  vector<string> doubleMuNames;
  doubleMuNames.push_back("HLT_Mu17_Mu8_v");
  doubleMuNames.push_back("HLT_Mu17_TkMu8_v");

  bool bitFired = IsAnyTriggerBitFired ( doubleMuNames, triggerBits);

  if (triggerDebug) cout << "MuonsMatchDoubleMuTrigger: Checking to see if trigger fired ... " << bitFired << endl;
  
  if (!bitFired) return false;

  vector<string> doubleMuFilters;
  doubleMuFilters.push_back("hltDiMuonGlb17Glb8DzFiltered0p2::HLT");
  doubleMuFilters.push_back("hltDiMuonGlb17Trk8DzFiltered0p2::HLT");
  

  if (triggerDebug) cout << "MuonsMatchDoubleMuTrigger: Now checking for matches" << endl;
  bool matchTrigger  = DoubleObjectMatchesAnyTrigger ( iMuon.eta, iMuon.phi, jMuon.eta, jMuon.phi, doubleMuFilters, triggerObjects);
  

  return matchTrigger;
  
}


bool BEANhelper::ElectronMatchesSingleEleTrigger(const BNelectron& iEle, const BNtriggerCollection & triggerBits, const BNtrigobjCollection & triggerObjects ){

  // Define your target triggers
  vector<string> singleEleNames;
  singleEleNames.push_back("HLT_Ele27_WP80_v");
  //singleEleNames.push_back("HLT_IsoMu24_v");


  // Now check to see if there is an HLT object
  // that corresponds to your trigger

  bool bitFired = IsAnyTriggerBitFired ( singleEleNames, triggerBits);

  if (!bitFired) return false;

  vector<string> singleEleFilters;
  singleEleFilters.push_back("hltEle27WP80TrackIsoFilter::HLT");
  //singleEleFilters.push_back("hltL3crIsoL1sMu16L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15::HLT");
  
  bool matchTrigger = SingleObjectMatchesAnyTrigger(iEle.eta, iEle.phi, singleEleFilters, triggerObjects);
  
  return matchTrigger;
  
}

bool BEANhelper::ElectronsMatchDoubleEleTrigger(const BNelectron& iEle, const BNelectron & jEle, const BNtriggerCollection & triggerBits, const BNtrigobjCollection & triggerObjects){


  // Define your target triggers
  vector<string> doubleEleNames;
  doubleEleNames.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  //doubleEleNames.push_back("");

  bool bitFired = IsAnyTriggerBitFired ( doubleEleNames, triggerBits);

  if(triggerDebug) cout << ": Checking to see if trigger fired ... " << bitFired << endl;
  
  if (!bitFired) return false;

  vector<string> doubleEleFilters;
  doubleEleFilters.push_back("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ::HLT");
  //doubleEleFilters.push_back("hltDiMuonGlb17Trk8DzFiltered0p2::HLT");
  

  if (triggerDebug) cout << ": Now checking for matches" << endl;
  bool matchTrigger  = DoubleObjectMatchesAnyTrigger ( iEle.eta, iEle.phi, jEle.eta, jEle.phi, doubleEleFilters, triggerObjects);
  

  return matchTrigger;
  
}


bool BEANhelper::MuEGMatchMuEGTrigger(const BNmuon& iMuon, const BNelectron & jEle, const BNtriggerCollection & triggerBits, const BNtrigobjCollection & triggerObjects){


  // Define your target triggers
  vector<string> muEGNames;
  muEGNames.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  muEGNames.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

  bool bitFired = IsAnyTriggerBitFired ( muEGNames, triggerBits);

  if (triggerDebug) cout << ": Checking to see if trigger fired ... " << bitFired << endl;
  
  if (!bitFired) return false;

  vector<string> muEGFilters;

  // Start with the Mu8_Ele17 filters
  muEGFilters.push_back("hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter::HLT");
  muEGFilters.push_back("hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8::HLT");

  if (triggerDebug) cout << ": Now checking for matches" << endl;
  bool match_Mu8_Ele17  = DoubleObjectMatchesAnyTrigger ( iMuon.eta, iMuon.phi, jEle.eta, jEle.phi, muEGFilters, triggerObjects, true);

  muEGFilters.clear();
  muEGFilters.push_back("hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter::HLT");  
  muEGFilters.push_back("hltL1Mu12EG7L3MuFiltered17::HLT");

  bool match_Mu17_Ele8  = DoubleObjectMatchesAnyTrigger ( iMuon.eta, iMuon.phi, jEle.eta, jEle.phi, muEGFilters, triggerObjects, true);

  if (triggerDebug) cout << ": Results Mu8_Ele17 = " << match_Mu8_Ele17 << ", results Mu17_Ele8 = " << match_Mu17_Ele8 << endl;
  
  return (match_Mu8_Ele17 || match_Mu17_Ele8);
  
}




// Return collection with objects passing cuts
BNmuonCollection BEANhelper::GetSelectedMuons(const BNmuonCollection& iMuons, const muonID::muonID iMuonID){
	CheckSetUp();
	BNmuonCollection result;
	for( BNmuonCollection::const_iterator Muon = iMuons.begin(); Muon != iMuons.end(); ++Muon ){
		if(IsGoodMuon((*Muon), iMuonID)){ result.push_back(*Muon); }
	}
	return result;
}

// Electron relative isolation
// Return whether or not muon passes cuts
bool BEANhelper::IsSideElectron(const BNelectron& iElectron){ return IsGoodElectron(iElectron, electronID::electronSide); }

bool BEANhelper::IsLooseElectron(const BNelectron& iElectron){ return IsGoodElectron(iElectron, electronID::electronLoose); }

bool BEANhelper::IsTightElectron(const BNelectron& iElectron){ return IsGoodElectron(iElectron, electronID::electronTight); }

float BEANhelper::GetElectronRelIso(const BNelectron& iElectron){ 
	CheckSetUp();
	float result = 9999;
    if (era == "2011") {
      if( usePfLeptons )
        result = ((iElectron.chargedHadronIso + iElectron.neutralHadronIso + iElectron.photonIso)/iElectron.pt);
      else
        result = ((iElectron.trackIsoDR03 + iElectron.ecalIsoDR03 + iElectron.hcalIsoDR03)/iElectron.pt);
    }
    else if (era == "2012_52x" || era == "2012_53x") {
      if( usePfLeptons ) 
        result = (((iElectron.chargedHadronIso + max(0.0, iElectron.neutralHadronIso + iElectron.photonIso - iElectron.AEffDr03*iElectron.rhoPrime)))/iElectron.pt);
      else
        result = (((iElectron.chargedHadronIso + max(0.0, iElectron.neutralHadronIso + iElectron.photonIso - iElectron.AEffDr03*iElectron.rhoPrime)))/iElectron.pt);
    }
    else {
      assert (era == "either 2012_52x, 2012_53x, or 2011");
    }
	return result;
}

float BEANhelper::GetElectronSF(const BNelectron& iElectron, const electronID::electronID inputID){
  CheckSetUp();
  double SF = 1.0;
  if (isData) return SF;

  double usePT = -10;
  double useEta = -10;
  
  switch (inputID) {

  case electronID::electronTight:
  case electronID::electronTightMinusTrigPresel:
    usePT = std::min (iElectron.pt, 199.0);
    useEta = std::min (fabs(iElectron.eta), 2.39);
    // ID*ISO and trigger are split 
    SF = h_ele_SF_->GetBinContent(h_ele_SF_->FindBin(useEta,usePT));
    if (isLJ)
      SF = SF * h_SingleEle_trig_SF_->GetBinContent(h_SingleEle_trig_SF_->FindBin(useEta, usePT));
    break;

  case electronID::electronLoose:
  case electronID::electronSide:
  case electronID::electronLooseMinusTrigPresel:
    usePT = std::min (iElectron.pt, 199.0);
    usePT = std::max (usePT, 15.1);
    useEta = std::min (fabs(iElectron.eta), 2.39);                                                  
    SF = h_looseEleSF->GetBinContent(h_looseEleSF->FindBin(useEta, usePT));    
    break;

  }

//   if ( !isLJ ) {

//     // current binning is only 20 to 150
//     // so fix  it!
//     double usePT = std::min( iElectron.pt, 149. );
//     usePT = std::max(21.0, usePT);
//     //double useEta = ( iElectron.eta>0. ) ? std::min( 2.09, iElectron.eta ) : std::max( -2.09, iElectron.eta );
//     double useEta = std::min(fabs(iElectron.eta), 2.49);
//     SF = h_ele_SF_->GetBinContent( h_ele_SF_->FindBin(usePT, useEta) );

//   } else {
//     // current binning is only 20 to 150
//     // so fix  it!
//     double usePT = std::min( iElectron.pt, 499. );
//     usePT = std::max(11.0, usePT);
//     //double useEta = ( iElectron.eta>0. ) ? std::min( 2.09, iElectron.eta ) : std::max( -2.09, iElectron.eta );
//     double useEta = iElectron.eta;
    
//     SF = h_ele_SF_->GetBinContent( h_ele_SF_->FindBin(usePT, useEta) );
//   }


  // Do a quick sanity check for the SF.
  // If it is below half a percent, then it is
  // probably way too small, and there was likely a bug
  // reading it from the histo.
  //---------->>>>>> Throw an exception
  if (SF < 0.005) {
    TString message = Form("Error! Electron SF too small (%f).\nCalled GetElectronSF with pt = %f eta = %f and tight = %d ",
                           SF, iElectron.pt, iElectron.eta, inputID==electronID::electronTight);
    ThrowFatalError(message.Data());

  }
  
  return SF;
}


// Electron ID function to keep things tidy
bool BEANhelper::GetElectronIDresult(const BNelectron& iElectron, const electronID::electronID iElectronID){
	CheckSetUp();

	// Electron ID stuff (DADT...)
	int eidHyperTight1MC		= iElectron.eidHyperTight1MC;
	bool eidHyperTight1MC_dec	= ( (eidHyperTight1MC & 1)==1 );
	int eidTight				= iElectron.eidTight;
	bool eidTight_dec			= ( (eidTight & 1)==1 );
	bool eid					= isLJ ? eidHyperTight1MC_dec : eidTight_dec;
	bool d0						= ( fabs(iElectron.correctedD0) < 0.02 );
	bool dZ						= ( fabs(iElectron.correctedDZ) < 1. );
	bool dist					= ( fabs(iElectron.dist)<0.02 );
	bool dcot					= ( fabs(iElectron.dcot)<0.02 );
	bool nlost					= ( iElectron.numberOfLostHits<1 );
    bool no_exp_inner_trkr_hits = ( iElectron.numberOfExpectedInnerHits <= 0 );
    
	// 2012 era-specific
	double mvaID				= iElectron.mvaTrigV0;
	bool passMVAId				= ( mvaID>0.0 ); // For 2012_52x
    bool passMVAId53x          = ( mvaID>0.5 );  // For 2012_53x, tighter selection
	bool d02					= ( fabs(iElectron.correctedD0Vertex) < 0.02 );
	bool d04					= ( fabs(iElectron.correctedD0Vertex) < 0.04 );


    bool myTrigPresel = false;
    if(fabs(iElectron.scEta) < 1.479) {
      // this is the 
      if(iElectron.scSigmaIEtaIEta < 0.014 &&
         iElectron.hadOverEm < 0.15 &&
         iElectron.trackIsoDR03 / iElectron.pt < 0.2 &&
         iElectron.ecalIsoDR03/iElectron.pt < 0.2 &&
         iElectron.hcalIsoDR03/iElectron.pt < 0.2 &&
         iElectron.numberOfLostHits == 0)
        myTrigPresel = true;
    }
    else {
      if(iElectron.scSigmaIEtaIEta < 0.035 &&
         iElectron.hadOverEm < 0.10 &&
         iElectron.trackIsoDR03/iElectron.pt < 0.2 &&
         iElectron.ecalIsoDR03/iElectron.pt < 0.2 &&
         iElectron.hcalIsoDR03/iElectron.pt < 0.2 &&
         iElectron.numberOfLostHits == 0)
        myTrigPresel = true;
    }

	if(era=="2011"){
		bool notConv				= ( !(dist && dcot) && nlost );
		bool id						= ( eid && d0 && dZ && notConv );
		if(iElectronID==electronID::electronSide){			return true; }
		else if(iElectronID==electronID::electronLoose){    return true; }
		else if(iElectronID==electronID::electronTight){	return id; }
	}
    else if(era=="2012_52x"){
		bool notConv				= ( iElectron.passConvVeto );
		bool id						= ( passMVAId && d02 && dZ && notConv );

		if(iElectronID==electronID::electronSide){			return (passMVAId && d04 && notConv); }
		else if(iElectronID==electronID::electronLoose){  	return (passMVAId && d04 && notConv && myTrigPresel); }
		else if(iElectronID==electronID::electronTight){	return (id && myTrigPresel);}
        else if ( iElectronID == electronID::electronTightMinusTrigPresel ){
          return id;
        } else if ( iElectronID == electronID::electronLooseMinusTrigPresel ){
          return (passMVAId && d04 && notConv);
        }
	}
    else if(era=="2012_53x"){
		bool notConv				= ( iElectron.passConvVeto );
		bool id						= ( passMVAId53x && d02 && dZ && notConv );

		if(iElectronID==electronID::electronSide){			return (passMVAId53x && d04 && notConv); }
		else if(iElectronID==electronID::electronLoose){  	return (passMVAId53x && no_exp_inner_trkr_hits && d04 && notConv && myTrigPresel); }
		else if(iElectronID==electronID::electronTight){	return (id && no_exp_inner_trkr_hits && myTrigPresel);}
        else if ( iElectronID == electronID::electronTightMinusTrigPresel ){
          return (id && no_exp_inner_trkr_hits);
        } else if ( iElectronID == electronID::electronLooseMinusTrigPresel ){
          return (passMVAId53x && no_exp_inner_trkr_hits && d04 && notConv);
        }
	}    
    else {
      assert (era == "either 2012_52x, 2012_53x, or 2011");
    }


	return false;
}

// Return whether or not electron passes cuts
bool BEANhelper::IsGoodElectron(const BNelectron& iElectron, const electronID::electronID iElectronID){
	CheckSetUp();

	// Set default kinematic thresholds
	float minLooseElectronPt		= 9999;
	float minTightElectronPt		= 9999;
	float maxLooseElectronAbsEta	= 0;
	float maxTightElectronAbsEta	= 0;

    if (era=="2011") {
      minLooseElectronPt		= 10;
      minTightElectronPt		= ( isLJ ) ? 30. : 20.;
      maxLooseElectronAbsEta	= 2.5;
      maxTightElectronAbsEta	= 2.5;
    }
    else if (era=="2012_52x" || era=="2012_53x") {
      minLooseElectronPt		= 10;
      minTightElectronPt		= ( isLJ ) ? 30. : 20.;
      maxLooseElectronAbsEta	= 2.5;
      maxTightElectronAbsEta	= 2.5;
    }
    else {
	  cout << "Era set to '" << era << "'" << endl;
      assert (era == "either 2012_52x, 2012_53x, or 2011");
    }

	// Be skeptical about this electron making it through
	bool passesKinematics	= false;
	bool passesIso			= false;
	bool passesID			= false;

	// Check if this electron is good enough
	bool inCrack			= ((fabs(iElectron.scEta) > 1.4442) && (fabs(iElectron.scEta) < 1.5660));
    if (era=="2011") {
      switch(iElectronID){
      case electronID::electronSide:
        passesKinematics	= ((iElectron.pt >= minLooseElectronPt) && (fabs(iElectron.eta) <= maxLooseElectronAbsEta) && (!inCrack));
        passesIso			= (GetElectronRelIso(iElectron) < 0.800);
        passesID			= GetElectronIDresult(iElectron, iElectronID);
        break;
        
      case electronID::electronLoose:
      case electronID::electronLooseMinusTrigPresel:
        passesKinematics	= ((iElectron.pt >= minLooseElectronPt) && (fabs(iElectron.eta) <= maxLooseElectronAbsEta) && (!inCrack));
        passesIso			= (GetElectronRelIso(iElectron) < 0.200);
        passesID			= GetElectronIDresult(iElectron, iElectronID);
        break;

      case electronID::electronTight:
      case electronID::electronTightMinusTrigPresel:
        passesKinematics	= ((iElectron.pt >= minTightElectronPt) && (fabs(iElectron.eta) <= maxTightElectronAbsEta) && (!inCrack));
        passesIso			= (GetElectronRelIso(iElectron) < 0.100);
        passesID			= GetElectronIDresult(iElectron, iElectronID);
        break;
      }
    } // End of 2011 era
    else if (era=="2012_52x" || era=="2012_53x") {
      switch(iElectronID){
      case electronID::electronSide:
        passesKinematics	= ((iElectron.pt >= minLooseElectronPt) && (fabs(iElectron.eta) <= maxLooseElectronAbsEta) && (!inCrack));
        passesIso			= (GetElectronRelIso(iElectron) < 0.800);
        passesID			= GetElectronIDresult(iElectron, iElectronID);
        break;
        
      case electronID::electronLoose:
      case electronID::electronLooseMinusTrigPresel:
        passesKinematics	= ((iElectron.pt >= minLooseElectronPt) && (fabs(iElectron.eta) <= maxLooseElectronAbsEta) && (!inCrack));
        passesIso			= (GetElectronRelIso(iElectron) < 0.200);
        passesID			= GetElectronIDresult(iElectron, iElectronID);
        break;
        
      case electronID::electronTight:
      case electronID::electronTightMinusTrigPresel:
        passesKinematics	= ((iElectron.pt >= minTightElectronPt) && (fabs(iElectron.eta) <= maxTightElectronAbsEta) && (!inCrack));
        passesIso			= (GetElectronRelIso(iElectron) < 0.100);
        passesID			= GetElectronIDresult(iElectron, iElectronID);
        break;
      }
	} // End of 2012 era
    else {
	  cout << "Era set to '" << era << "'" << endl;
      assert (era == "either 2012_52x, 2012_53x, or 2011");
    }

	return (passesKinematics && passesIso && passesID);
}

// Return collection with objects passing cuts
BNelectronCollection BEANhelper::GetSelectedElectrons(const BNelectronCollection& iElectrons, const electronID::electronID iElectronID){
	CheckSetUp();
	BNelectronCollection result;
	for( BNelectronCollection::const_iterator Electron = iElectrons.begin(); Electron != iElectrons.end(); ++Electron ){
		if(IsGoodElectron((*Electron), iElectronID)){ result.push_back(*Electron); }
	}
	return result;
}

BNjetCollection BEANhelper::GetCleanJets(const BNjetCollection& iJets, const vector<TLorentzVector>& iToUnmatch, const float iDeltaR){
	CheckSetUp();
	BNjetCollection result;
	for( BNjetCollection::const_iterator Jet = iJets.begin(); Jet != iJets.end(); ++Jet ){
		bool matched = false;
		for(unsigned int m = 0; m < iToUnmatch.size(); m++){
			if(deltaR(Jet->eta, Jet->phi, iToUnmatch.at(m).Eta(), iToUnmatch.at(m).Phi()) < iDeltaR){ matched = true; break; }
		}

		if(!matched){ result.push_back(*Jet); }
	}
	return result;
}

unsigned int BEANhelper::GetNumCSVbtags(const BNjetCollection& iJets, const char iCSVwp){
	CheckSetUp();
	unsigned int result = 0;
	for( BNjetCollection::const_iterator Jet = iJets.begin(); Jet != iJets.end(); ++Jet ){ if(PassesCSV(*Jet, iCSVwp)){ result++; } }
	return result;
}

unsigned int BEANhelper::GetNumNonCSVbtags(const BNjetCollection& iJets, const char iCSVwp){
	CheckSetUp();
	unsigned int result = 0;
	for( BNjetCollection::const_iterator Jet = iJets.begin(); Jet != iJets.end(); ++Jet ){ if(!PassesCSV(*Jet, iCSVwp)){ result++; } }
	return result;
}



// ******************** MCparticle ****************** //
BNmcparticleCollection BEANhelper::GetSelectedMCparticlesByPDGid(const BNmcparticleCollection& iMCparticles, const vector<int> iPDGid){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( find(iPDGid.begin(), iPDGid.end(), MCparticle->id) != iPDGid.end() ){
			result.push_back(*MCparticle);	
		}
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetSelectedMCparticlesByChildPDGid(const BNmcparticleCollection& iMCparticles, const vector<int> iPDGid){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( find(iPDGid.begin(), iPDGid.end(), MCparticle->daughter0Id) != iPDGid.end() || 
			find(iPDGid.begin(), iPDGid.end(), MCparticle->daughter1Id) != iPDGid.end() ){
			result.push_back(*MCparticle);	
		}
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetSelectedMCparticlesByParentPDGid(const BNmcparticleCollection& iMCparticles, const vector<int> iPDGid){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( find(iPDGid.begin(), iPDGid.end(), MCparticle->mother0Id) != iPDGid.end() || 
			find(iPDGid.begin(), iPDGid.end(), MCparticle->mother1Id) != iPDGid.end() ){
			result.push_back(*MCparticle);	
		}
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetSelectedMCparticlesByGrandParentPDGid(const BNmcparticleCollection& iMCparticles, const vector<int> iPDGid){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( find(iPDGid.begin(), iPDGid.end(), MCparticle->grandMother00Id) != iPDGid.end() || 
			find(iPDGid.begin(), iPDGid.end(), MCparticle->grandMother01Id) != iPDGid.end() ||
			find(iPDGid.begin(), iPDGid.end(), MCparticle->grandMother10Id) != iPDGid.end() ||
			find(iPDGid.begin(), iPDGid.end(), MCparticle->grandMother11Id) != iPDGid.end() ){
			result.push_back(*MCparticle);	
		}
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetUnrejectedMCparticlesByPDGid(const BNmcparticleCollection& iMCparticles, const vector<int> iPDGid){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( find(iPDGid.begin(), iPDGid.end(), MCparticle->id) == iPDGid.end() ){
			result.push_back(*MCparticle);	
		}
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetUnrejectedMCparticlesByChildPDGid(const BNmcparticleCollection& iMCparticles, const vector<int> iPDGid){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( find(iPDGid.begin(), iPDGid.end(), MCparticle->daughter0Id) == iPDGid.end() &&
			find(iPDGid.begin(), iPDGid.end(), MCparticle->daughter1Id) == iPDGid.end()
		){
			result.push_back(*MCparticle);	
		}
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetUnrejectedMCparticlesByParentPDGid(const BNmcparticleCollection& iMCparticles, const vector<int> iPDGid){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( find(iPDGid.begin(), iPDGid.end(), MCparticle->mother0Id) == iPDGid.end() &&
			find(iPDGid.begin(), iPDGid.end(), MCparticle->mother1Id) == iPDGid.end()
		){
			result.push_back(*MCparticle);	
		}
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetUnrejectedMCparticlesByGrandParentPDGid(const BNmcparticleCollection& iMCparticles, const vector<int> iPDGid){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( find(iPDGid.begin(), iPDGid.end(), MCparticle->grandMother00Id) == iPDGid.end() &&
			find(iPDGid.begin(), iPDGid.end(), MCparticle->grandMother01Id) == iPDGid.end() &&
			find(iPDGid.begin(), iPDGid.end(), MCparticle->grandMother10Id) == iPDGid.end() &&
			find(iPDGid.begin(), iPDGid.end(), MCparticle->grandMother11Id) == iPDGid.end()
		){
			result.push_back(*MCparticle);	
		}
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetSelectedMCparticlesByChildStatus(const BNmcparticleCollection& iMCparticles, const bool iKeepStatus1, const bool iKeepStatus2, const bool iKeepStatus3){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( (((MCparticle->daughter0Status == 1) && iKeepStatus1) || ((MCparticle->daughter0Status == 2) && iKeepStatus2) || ((MCparticle->daughter0Status == 3) && iKeepStatus3)) ||
			(((MCparticle->daughter1Status == 1) && iKeepStatus1) || ((MCparticle->daughter1Status == 2) && iKeepStatus2) || ((MCparticle->daughter1Status == 3) && iKeepStatus3))
		){ result.push_back(*MCparticle); }
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetSelectedMCparticlesByStatus(const BNmcparticleCollection& iMCparticles, const bool iKeepStatus1, const bool iKeepStatus2, const bool iKeepStatus3){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( ((MCparticle->status == 1) && iKeepStatus1) ||
			((MCparticle->status == 2) && iKeepStatus2) ||
			((MCparticle->status == 3) && iKeepStatus3)
		){ result.push_back(*MCparticle); }
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetSelectedMCparticlesByParentStatus(const BNmcparticleCollection& iMCparticles, const bool iKeepStatus1, const bool iKeepStatus2, const bool iKeepStatus3){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( (((MCparticle->mother0Status == 1) && iKeepStatus1) || ((MCparticle->mother0Status == 2) && iKeepStatus2) || ((MCparticle->mother0Status == 3) && iKeepStatus3)) ||
			(((MCparticle->mother1Status == 1) && iKeepStatus1) || ((MCparticle->mother1Status == 2) && iKeepStatus2) || ((MCparticle->mother1Status == 3) && iKeepStatus3))
		){ result.push_back(*MCparticle); }
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetSelectedMCparticlesByGrandParentStatus(const BNmcparticleCollection& iMCparticles, const bool iKeepStatus1, const bool iKeepStatus2, const bool iKeepStatus3){
	BNmcparticleCollection result;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		if( (((MCparticle->grandMother00Status == 1) && iKeepStatus1) || ((MCparticle->grandMother00Status == 2) && iKeepStatus2) || ((MCparticle->grandMother00Status == 3) && iKeepStatus3)) ||
			(((MCparticle->grandMother01Status == 1) && iKeepStatus1) || ((MCparticle->grandMother01Status == 2) && iKeepStatus2) || ((MCparticle->grandMother01Status == 3) && iKeepStatus3)) ||
			(((MCparticle->grandMother10Status == 1) && iKeepStatus1) || ((MCparticle->grandMother10Status == 2) && iKeepStatus2) || ((MCparticle->grandMother10Status == 3) && iKeepStatus3)) ||
			(((MCparticle->grandMother11Status == 1) && iKeepStatus1) || ((MCparticle->grandMother11Status == 2) && iKeepStatus2) || ((MCparticle->grandMother11Status == 3) && iKeepStatus3))
		){ result.push_back(*MCparticle); }
	}
	return result;
}

BNmcparticleCollection BEANhelper::GetParents(const BNmcparticle& iMCparticle, const BNmcparticleCollection& iMCparticles){
	BNmcparticleCollection result;

	for( BNmcparticleCollection::const_iterator potentialParent = iMCparticles.begin(); potentialParent != iMCparticles.end(); ++potentialParent){

		// Check that at least one of the children is iMCparticle
		if(potentialParent->daughter0Id != iMCparticle.id && potentialParent->daughter1Id != iMCparticle.id){ continue; }

		// Check that at least on of the parents of the iMCparticle is the potential parent
		if(iMCparticle.mother0Id != potentialParent->id && iMCparticle.mother1Id != potentialParent->id){ continue; }

		// Check that 
		if(
			iMCparticle.grandMother00Id != potentialParent->mother0Id &&
			iMCparticle.grandMother01Id != potentialParent->mother0Id &&
			iMCparticle.grandMother10Id != potentialParent->mother0Id &&
			iMCparticle.grandMother11Id != potentialParent->mother0Id &&
			iMCparticle.grandMother00Id != potentialParent->mother1Id &&
			iMCparticle.grandMother01Id != potentialParent->mother1Id &&
			iMCparticle.grandMother10Id != potentialParent->mother1Id &&
			iMCparticle.grandMother11Id != potentialParent->mother1Id
			){ continue; }

			result.push_back(*potentialParent);

	}

	return result;
}

unsigned int BEANhelper::GetNumExtraPartons(const BNmcparticleCollection& iMCparticles){

  int nList = 0;
  unsigned int nPartons = 0;
  for( unsigned i=0; i< iMCparticles.size(); i++ ){
    nList ++;
    int id = iMCparticles.at(i).id;
    int motherID = iMCparticles.at(i).motherId;
    int status = iMCparticles.at(i).status;
//     int mother0Status = iMCparticles.at(i).mother0Status;
//     int mother1Status = iMCparticles.at(i).mother1Status;
    int aid = abs(id);

    if(status == 3){//only status 3 particles (which are listed first)
      if(nList>6){//dont look at first 6 (incomming event)
        if( (aid>0 && aid<6) || aid ==21 || aid ==9){//udscb gluon
          if(abs(motherID) !=23 && abs(motherID) !=24 && abs(motherID) != 25 && abs(motherID)!=6 ){ //not from WZHt
            nPartons++;
          }
        }
      }
    }else{
      break; //break when you run out of stat==3
    }
  }

  return nPartons;

}


hdecayType::hdecayType BEANhelper::GetHdecayType(const BNmcparticleCollection& iMCparticles){

  int numH=0;
  int daughter0id = -99, daughter1id = -99;
  for( unsigned i=0; i< iMCparticles.size(); i++ ){
    int id = iMCparticles.at(i).id;
    int status = iMCparticles.at(i).status;
    if( id!=25 || status!=3 ) continue;

    numH++;

    daughter0id = iMCparticles.at(i).daughter0Id;
    daughter1id = iMCparticles.at(i).daughter1Id;
  }

  if( numH!=1 ) std::cout << "  Expected to find 1 Higgs and found " << numH << " Higg(es) " << std::endl;

  daughter0id = abs(daughter0id);
  daughter1id = abs(daughter1id);

  hdecayType::hdecayType decayType = hdecayType::hvv;
  if( daughter0id==24 && daughter1id==24 )      decayType = hdecayType::hww;
  else if( daughter0id==23 && daughter1id==23 ) decayType = hdecayType::hzz;
  else if( (daughter0id==5 && daughter1id==5) ||
	   (daughter0id==4 && daughter1id==4))  decayType = hdecayType::hbb;
  else if( daughter0id==3  && daughter1id==3  ) decayType = hdecayType::hbb;
  else if( daughter0id==2  && daughter1id==2  ) decayType = hdecayType::hbb;
  else if( daughter0id==1  && daughter1id==1  ) decayType = hdecayType::hbb;
  else if( daughter0id==11 && daughter1id==11 ) decayType = hdecayType::htt;
  else if( daughter0id==13 && daughter1id==13 ) decayType = hdecayType::htt;
  else if( daughter0id==15 && daughter1id==15 ) decayType = hdecayType::htt;
  else if( (daughter0id==23 && daughter1id==22) || 
	   (daughter1id==23 && daughter0id==22) ||
	   (daughter1id==22 && daughter0id==22) ||
	   (daughter1id==21 && daughter0id==21)) decayType = hdecayType::hvv;
  else{
    std::cout << "\t UNCLASSIFIED EVENT !! daughter0id = " << daughter0id << ",\t daughter1id = " << daughter1id << std::endl;
  }

  return decayType;
}

bool BEANhelper::keepHdecayType(const BNmcparticleCollection& iMCparticles, const hdecayType::hdecayType targetType){

  hdecayType::hdecayType eventType = BEANhelper::GetHdecayType( iMCparticles );
  bool keepEvent = ( targetType==eventType );
  return keepEvent;
}


// (Sort of) draws a feynman diagram with the particle id's of the particle itself, daughters, mothers and grandmothers
void BEANhelper::DrawFeynman(const BNmcparticle& iMCparticle){

	// Prepare particle IDs
	stringstream d0;	d0.str("");
	stringstream d1;	d1.str("");
	stringstream p;		p.str("");
	stringstream m0;	m0.str("");
	stringstream m1;	m1.str("");
	stringstream g00;	g00.str("");
	stringstream g01;	g01.str("");
	stringstream g10;	g10.str("");
	stringstream g11;	g11.str("");

											p	<< iMCparticle.id;
	if(     iMCparticle.daughter0Id != -99){ d0	<< iMCparticle.daughter0Id;		}
	if(     iMCparticle.daughter1Id != -99){ d1	<< iMCparticle.daughter1Id;		};
	if(       iMCparticle.mother0Id != -99){ m0	<< iMCparticle.mother0Id;	 	};
	if(       iMCparticle.mother1Id != -99){ m1	<< iMCparticle.mother1Id;	  	};
	if( iMCparticle.grandMother00Id != -99){ g00	<< iMCparticle.grandMother00Id;	};
	if( iMCparticle.grandMother01Id != -99){ g01	<< iMCparticle.grandMother01Id;	};
	if( iMCparticle.grandMother10Id != -99){ g10	<< iMCparticle.grandMother10Id;	};
	if( iMCparticle.grandMother11Id != -99){ g11	<< iMCparticle.grandMother11Id;	};

	// Prepare each line
	stringstream ssout; ssout.str("");
	ssout << string( 2,' ') << setw(5) << setfill(' ') << g00.str()	<< endl;

	ssout << string( 8,' ') << setw(5) << setfill(' ') << m0.str()	<< string( 7,' ') << d0.str() << endl;

	ssout << string( 2,' ') << setw(5) << setfill(' ') << g01.str()	<< string( 6,' '); 
	ssout << (  (iMCparticle.mother0Id != -99) ? "\\":" ") << string(5,' ');
	ssout << ((iMCparticle.daughter0Id != -99) ?  "/":" ") << endl;

	ssout << string(14,' ') << setw(3) << setfill(' ') << p.str()	<< endl;

	ssout << string( 2,' ') << setw(5) << setfill(' ') << g10.str()	<< string( 6,' '); 
	ssout << (  (iMCparticle.mother1Id != -99) ?  "/":" ") << string(5,' ');
	ssout << ((iMCparticle.daughter1Id != -99) ? "\\":" ") << endl;

	ssout << string( 8,' ') << setw(5) << setfill(' ') << m1.str()	<< string( 7,' ') << d1.str() << endl;

	ssout << string( 2,' ') << setw(5) << setfill(' ') << g11.str()	<< endl;

	// Print end result
	cout << "\n" << string(40,'-') << endl << ssout.str() << string(40,'-') << endl;
}

BNmcparticleCollection BEANhelper::GetGenTaus(const BNmcparticleCollection& iMCparticles){
	vector<int> tauIDs; tauIDs.push_back(15); tauIDs.push_back(-15);
	return GetSelectedMCparticlesByPDGid(iMCparticles, tauIDs);
}

BNmcparticleCollection BEANhelper::GetHadronicGenTaus(const BNmcparticleCollection& iMCparticles){
	BNmcparticleCollection result;
	BNmcparticleCollection genTaus = GetGenTaus(iMCparticles);
	for( BNmcparticleCollection::const_iterator genTau = genTaus.begin(); genTau != genTaus.end(); ++genTau){
		if((abs(genTau->daughter0Id) == 16) || (abs(genTau->daughter1Id) == 16)){ result.push_back(*genTau); }
	}

	return result;
}

// Return the visible gentau (neutrino momentum removed)
BNmcparticle BEANhelper::GetVisGenTau(const BNmcparticle& iTau, const BNmcparticleCollection& iMCparticles){
	BNmcparticle result = iTau;
	TLorentzVector tauP4(iTau.px, iTau.py, iTau.pz, iTau.energy);

	vector<int> tauId(1, iTau.id);
	vector<int> tauParent0Id(1, iTau.mother0Id);
	vector<int> tauParent1Id(1, iTau.mother1Id);

	// Get all particles whose parent is a tau
	BNmcparticleCollection particlesWithTauParents						= GetSelectedMCparticlesByParentPDGid(iMCparticles, tauId);
	// Out of those, keep only those whose grandparent0 is a tau parent
	BNmcparticleCollection particlesWithTauParentsMatchedToTauParent0	= GetSelectedMCparticlesByGrandParentPDGid(particlesWithTauParents, tauParent0Id);
	// Out of those, keep only those whose grandparent1 is a tau parent
	BNmcparticleCollection particlesWithTauParentsMatchedToTauParents	= GetSelectedMCparticlesByGrandParentPDGid(particlesWithTauParentsMatchedToTauParent0, tauParent1Id);

	// Subtrack those neutrinos that are tau children
	for(BNmcparticleCollection::const_iterator tauChild = particlesWithTauParentsMatchedToTauParents.begin(); tauChild != particlesWithTauParentsMatchedToTauParents.end(); ++tauChild){
		if( abs(tauChild->id) == 12 || abs(tauChild->id) == 14 || abs(tauChild->id) == 16 ){
			TLorentzVector neutrinoP4(tauChild->px, tauChild->py, tauChild->pz, tauChild->energy);
			tauP4 = tauP4-neutrinoP4;
		}
	}

	result.mass		= tauP4.M();
	result.energy	= tauP4.Energy();
	result.et		= tauP4.Et();
	result.pt		= tauP4.Pt();
	result.px		= tauP4.Px();
	result.py		= tauP4.Py();
	result.pz		= tauP4.Pz();
	result.phi		= tauP4.Phi();
	result.eta		= tauP4.Eta();
	result.theta	= tauP4.Theta();

	return result;
}

// Is this event kept by tautaulepton analysis?
bool BEANhelper::IsTauTauLeptonEvent(const BNtauCollection& iTaus, const BNjetCollection& iJets, const BNelectronCollection& iElectrons, const BNmuonCollection& iMuons, const sysType::sysType iSysType){

	// B-jets: veto out events with 3 or more CSVM-tags
	BNjetCollection correctedJets				= GetCorrectedJets(iJets, iSysType);
	BNjetCollection selectedCorrectedCSVMJets	= GetSelectedJets(correctedJets, 30, 2.4, jetID::jetLoose, 'M');
	if(selectedCorrectedCSVMJets.size() >= 3){ return false; }

	// Leptons: Select events with exactly one 'tight' lepton and zero 'exclusively loose' leptons
	BNelectronCollection looseElectrons		= GetSelectedElectrons(iElectrons, electronID::electronLoose);
	BNelectronCollection tightElectrons		= GetSelectedElectrons(iElectrons, electronID::electronTight);
	//BNelectronCollection exLooseElectrons	= GetDifference(looseElectrons, tightElectrons);
	BNmuonCollection looseMuons				= GetSelectedMuons(iMuons, muonID::muonLoose);
	BNmuonCollection tightMuons				= GetSelectedMuons(iMuons, muonID::muonTight);
	//BNmuonCollection exLooseMuons			= GetDifference(looseMuons, tightMuons);
	if(tightMuons.size() + tightElectrons.size() != 1){ return false; }
	if(looseMuons.size() + looseElectrons.size() != 1){ return false; } // Note 'loose' includes 'tight'

	// Taus: Clean taus which are also reconstructed as leptons (loose, i.e. exLoose + tight)
	BNtauCollection correctedTaus	= GetCorrectedTaus(iTaus, iSysType);
	BNtauCollection muonlessTaus	= GetDifference(correctedTaus, looseMuons, 0.25);
	BNtauCollection leptonlessTaus	= GetDifference(muonlessTaus, looseElectrons, 0.25);
	BNtauCollection nonIsoTaus		= GetSelectedTaus(leptonlessTaus, tauID::tauNonIso);
	BNtauCollection vlooseTaus		= GetSelectedTaus(nonIsoTaus, tauID::tauVLoose);
	if(nonIsoTaus.size() < 2){ return false; }
	if(vlooseTaus.size() < 1){ return false; }
		
	return true;
}

// === PU reweighing === //
double BEANhelper::GetPUweight(const double iNumBX0){
  if (isData) return 1.0;
  return h_PU_ratio->GetBinContent( h_PU_ratio->FindBin( iNumBX0 ) );
 }

double BEANhelper::GetPUweightUp(const double iNumBX0){
  if (isData) return 1.0;
  return h_PUup_ratio->GetBinContent( h_PUup_ratio->FindBin( iNumBX0 ) );
}

double BEANhelper::GetPUweightDown(const double iNumBX0){
  if (isData) return 1.0;
  return h_PUdown_ratio->GetBinContent( h_PUdown_ratio->FindBin( iNumBX0 ) );
}



// === Jet CSV reweighting === //
double BEANhelper::GetCSVweight(const BNjetCollection& iJets, const sysType::sysType iSysType){
  if (isData) return 1.0;

  CheckSetUp();
  // IMPORTANT! iJets is the *SELECTED* jet collection which you use for your analysis

  int iSysHF = 0;
  switch(iSysType){
  case sysType::JESup:	         iSysHF=1; break;
  case sysType::JESdown:         iSysHF=2; break;
  case sysType::CSVLFup:         iSysHF=3; break;
  case sysType::CSVLFdown:       iSysHF=4; break;
  case sysType::CSVHFStats1up:	 iSysHF=5; break;
  case sysType::CSVHFStats1down: iSysHF=6; break;
  case sysType::CSVHFStats2up:	 iSysHF=7; break;
  case sysType::CSVHFStats2down: iSysHF=8; break;
  default : iSysHF = 0; break;
  }

  int iSysLF = 0;
  switch(iSysType){
  case sysType::JESup:	         iSysLF=1; break;
  case sysType::JESdown:         iSysLF=2; break;
  case sysType::CSVHFup:         iSysLF=3; break;
  case sysType::CSVHFdown:       iSysLF=4; break;
  case sysType::CSVLFStats1up:	 iSysLF=5; break;
  case sysType::CSVLFStats1down: iSysLF=6; break;
  case sysType::CSVLFStats2up:	 iSysLF=7; break;
  case sysType::CSVLFStats2down: iSysLF=8; break;
  default : iSysLF = 0; break;
  }

  double csvWgthf = 1.;
  double csvWgtlf = 1.;

  for( BNjetCollection::const_iterator iJet = iJets.begin(); iJet != iJets.end(); ++iJet ){ 

    double csv = iJet->btagCombinedSecVertex;
    double jetPt = iJet->pt;
    double jetAbsEta = fabs( iJet->eta );
    int flavor = abs( iJet->flavour );

    int iPt = -1; int iEta = -1;
    if (jetPt >=29.99 && jetPt<40) iPt = 0;
    else if (jetPt >=40 && jetPt<60) iPt = 1;
    else if (jetPt >=60 && jetPt<100) iPt = 2;
    else if (jetPt >=100 && jetPt<160) iPt = 3;
    else if (jetPt >=160 && jetPt<10000) iPt = 4;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6) iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.41) iEta = 2;

    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

    if (abs(flavor) == 5 || abs(flavor) == 4){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;

      // if( iSysHF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      // 			     flavor, jetPt, iJet->eta, csv, iCSVWgtHF );
    }
    else {
      if (iPt >=2) iPt=2;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;

      // if( iSysLF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      // 			     flavor, jetPt, iJet->eta, csv, iCSVWgtLF );
    }
  }

  double csvWgtTotal = csvWgthf* csvWgtlf;

  return csvWgtTotal;
}

// === Jet CSV reweighting: return lfSF and hfSF separately === //
vdouble BEANhelper::GetCSVweights(const BNjetCollection& iJets, const sysType::sysType iSysType){
  vdouble result;
  
  if (isData) {
    result.clear();
    result.push_back(1.0);
    result.push_back(1.0);
    return result;
  }

  CheckSetUp();
  // IMPORTANT! iJets is the *SELECTED* jet collection which you use for your analysis

  int iSysHF = 0;
  switch(iSysType){
  case sysType::JESup:	         iSysHF=1; break;
  case sysType::JESdown:         iSysHF=2; break;
  case sysType::CSVLFup:         iSysHF=3; break;
  case sysType::CSVLFdown:       iSysHF=4; break;
  case sysType::CSVHFStats1up:	 iSysHF=5; break;
  case sysType::CSVHFStats1down: iSysHF=6; break;
  case sysType::CSVHFStats2up:	 iSysHF=7; break;
  case sysType::CSVHFStats2down: iSysHF=8; break;
  default : iSysHF = 0; break;
  }

  int iSysLF = 0;
  switch(iSysType){
  case sysType::JESup:	         iSysLF=1; break;
  case sysType::JESdown:         iSysLF=2; break;
  case sysType::CSVHFup:         iSysLF=3; break;
  case sysType::CSVHFdown:       iSysLF=4; break;
  case sysType::CSVLFStats1up:	 iSysLF=5; break;
  case sysType::CSVLFStats1down: iSysLF=6; break;
  case sysType::CSVLFStats2up:	 iSysLF=7; break;
  case sysType::CSVLFStats2down: iSysLF=8; break;
  default : iSysLF = 0; break;
  }

  double csvWgthf = 1.;
  double csvWgtlf = 1.;

  for( BNjetCollection::const_iterator iJet = iJets.begin(); iJet != iJets.end(); ++iJet ){ 

    double csv = iJet->btagCombinedSecVertex;
    double jetPt = iJet->pt;
    double jetAbsEta = fabs( iJet->eta );
    int flavor = abs( iJet->flavour );

    int iPt = -1; int iEta = -1;
    if (jetPt >=29.99 && jetPt<40) iPt = 0;
    else if (jetPt >=40 && jetPt<60) iPt = 1;
    else if (jetPt >=60 && jetPt<100) iPt = 2;
    else if (jetPt >=100 && jetPt<160) iPt = 3;
    else if (jetPt >=160 && jetPt<10000) iPt = 4;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6) iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.41) iEta = 2;

    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

    if (abs(flavor) == 5 || abs(flavor) == 4){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;

      // if( iSysHF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      // 			     flavor, jetPt, iJet->eta, csv, iCSVWgtHF );
    }
    else {
      if (iPt >=2) iPt=2;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;

      // if( iSysLF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      // 			     flavor, jetPt, iJet->eta, csv, iCSVWgtLF );
    }
  }

  result.clear();
  result.push_back(csvWgthf);
  result.push_back(csvWgtlf);
  
  return result;
}


// === Top Pt reweighting === //
double BEANhelper::TopPtWeight(double topPt){
  if( topPt<0 ) return 1;

  double p0 = 1.18246e+00;
  double p1 = 4.63312e+02;
  double p2 = 2.10061e-06;

  if( topPt>p1 ) topPt = p1;

  double result = p0 + p2 * topPt * ( topPt - 2 * p1 );
  return result;
}


double BEANhelper::GetTopPtweight(const BNmcparticleCollection& iMCparticles){
  string samplename = GetSampleName();
  bool isTTJets = false;
  if (samplename == "ttbar_jj" || samplename == "ttbar_lj" || samplename == "ttbar_ll" ||
      samplename == "ttbar_bb_jj" || samplename == "ttbar_bb_lj" || samplename == "ttbar_bb_ll" ||
      samplename == "ttbar_cc_jj" || samplename == "ttbar_cc_lj" || samplename == "ttbar_cc_ll" ||
      samplename == "ttbar_bb" || samplename == "ttbar_cc")  isTTJets = true;

  if( !isTTJets ) return 1.;

  double topPt = -1;
  for( unsigned i=0; i< iMCparticles.size(); i++ ){
    if( iMCparticles.at(i).id==6 ){ topPt = iMCparticles.at(i).pt; break; }
  }
  double topPtSF = BEANhelper::TopPtWeight( topPt );
 
  return topPtSF;
}

double BEANhelper::GetTopPtweightUp(const BNmcparticleCollection& iMCparticles){
  string samplename = GetSampleName();
  bool isTTJets = false;
  if (samplename == "ttbar_jj" || samplename == "ttbar_lj" || samplename == "ttbar_ll" ||
      samplename == "ttbar_bb_jj" || samplename == "ttbar_bb_lj" || samplename == "ttbar_bb_ll" ||
      samplename == "ttbar_cc_jj" || samplename == "ttbar_cc_lj" || samplename == "ttbar_cc_ll" ||
      samplename == "ttbar_bb" || samplename == "ttbar_cc")  isTTJets = true;

  if( !isTTJets ) return 1.;

  double topPt = -1;
  for( unsigned i=0; i< iMCparticles.size(); i++ ){
    if( iMCparticles.at(i).id==6 ){ topPt = iMCparticles.at(i).pt; break; }
  }
  double topPtSF = BEANhelper::TopPtWeight( topPt );

  return 2*(topPtSF-1) + 1;
}

double BEANhelper::GetTopPtweightDown(const BNmcparticleCollection& iMCparticles){
  string samplename = GetSampleName();
  bool isTTJets = false;
  if (samplename == "ttbar_jj" || samplename == "ttbar_lj" || samplename == "ttbar_ll" ||
      samplename == "ttbar_bb_jj" || samplename == "ttbar_bb_lj" || samplename == "ttbar_bb_ll" ||
      samplename == "ttbar_cc_jj" || samplename == "ttbar_cc_lj" || samplename == "ttbar_cc_ll" ||
      samplename == "ttbar_bb" || samplename == "ttbar_cc")  isTTJets = true;
  if( !isTTJets ) return 1.;

  return 1.;
}


/*
Function to get the Q^2 scale weight

For the scale systematics samples, those should be normalized to the same expected yield as the nominal ttbar sample (before any cuts), 
e.g. cross section x integrated luminosity / number generated  

The functions Q2ScaleUpWgt and Q2ScaleDownWgt have an average weight that is different than 1.0.  That is, they would change the expected yield of these scale samples.  These numbers, 1.402 and 0.683, correct the systematics samples so that they have the same MC expectation as the nominal ttbar sample without any cuts.

The numbers above were derived running on the ttbar samples and adding up these weights (up and down), and taking the inverse.  This was done without placing any cuts whatsoever on the ttbar samples.  Running on the semi-leptonic, fully-hadronic, and fully-leptonic ttbar samples produce consistent numbers, so the scale factors above were derived running over each of these samples.

Another thing I want to make sure is clear: the Q^2 scale shape uncertainty (using the weights in the BEANs) is only applicable for ttbar samples.  Another procedure or a dedicated sample would be needed to get the systematics for W and Z.  For now, we apply a rate systematic for that uncertainty (it's already included in the datacards), because these are such small uncertainties for most of our analyses, but that may change in the future.
*/

double BEANhelper::GetQ2ScaleUp(const BNevent& event){
  string samplename = GetSampleName();
  bool isTTJets = false;
  if (samplename == "ttbar_jj" || samplename == "ttbar_lj" || samplename == "ttbar_ll" ||
      samplename == "ttbar_bb_jj" || samplename == "ttbar_bb_lj" || samplename == "ttbar_bb_ll" ||
      samplename == "ttbar_cc_jj" || samplename == "ttbar_cc_lj" || samplename == "ttbar_cc_ll" ||
      samplename == "ttbar_bb" || samplename == "ttbar_cc")  isTTJets = true;

  if( !isTTJets ) return 1.;

  return 1.402 * event.Q2ScaleUpWgt;
}

double BEANhelper::GetQ2ScaleDown(const BNevent& event){
  string samplename = GetSampleName();
  bool isTTJets = false;
  if (samplename == "ttbar_jj" || samplename == "ttbar_lj" || samplename == "ttbar_ll" ||
      samplename == "ttbar_bb_jj" || samplename == "ttbar_bb_lj" || samplename == "ttbar_bb_ll" ||
      samplename == "ttbar_cc_jj" || samplename == "ttbar_cc_lj" || samplename == "ttbar_cc_ll" ||
      samplename == "ttbar_bb" || samplename == "ttbar_cc")  isTTJets = true;

  if( !isTTJets ) return 1.;

  return 0.683 * event.Q2ScaleDownWgt;
}



// ******************** From BEANsUtilities.h ****************** //
double BEANhelper::getJERfactor( int returnType, double jetAbsETA, double genjetPT, double recojetPT){

    CheckSetUp();
    string samplename = GetSampleName();
	double factor = 1.;
    
	double scale_JER = 1., scale_JERup = 1., scale_JERdown = 1.;
    double diff_FullSim_FastSim = 0.;
	if( jetAbsETA<0.5 ){ 
		scale_JER = 1.052; scale_JERup = 1.052 + sqrt( 0.012*0.012 + 0.062*0.062 ); scale_JERdown = 1.052 - sqrt( 0.012*0.012 + 0.061*0.061 );
        diff_FullSim_FastSim = 1.35;
	}
	else if( jetAbsETA<1.1 ){ 
		scale_JER = 1.057; scale_JERup = 1.057 + sqrt( 0.012*0.012 + 0.056*0.056 ); scale_JERdown = 1.057 - sqrt( 0.012*0.012 + 0.055*0.055 );
        diff_FullSim_FastSim = 1.54;
	}
	else if( jetAbsETA<1.7 ){ 
		scale_JER = 1.096; scale_JERup = 1.096 + sqrt( 0.017*0.017 + 0.063*0.063 ); scale_JERdown = 1.096 - sqrt( 0.017*0.017 + 0.062*0.062 );
        diff_FullSim_FastSim = 1.97;
	}
	else if( jetAbsETA<2.3 ){ 
		scale_JER = 1.134; scale_JERup = 1.134 + sqrt( 0.035*0.035 + 0.087*0.087 ); scale_JERdown = 1.134 - sqrt( 0.035*0.035 + 0.085*0.085 );
        diff_FullSim_FastSim = 3.12;
	}
	else if( jetAbsETA<5.0 ){ 
		scale_JER = 1.288; scale_JERup = 1.288 + sqrt( 0.127*0.127 + 0.155*0.155 ); scale_JERdown = 1.288 - sqrt( 0.127*0.127 + 0.153*0.153 );
        diff_FullSim_FastSim = 3.12;
	}

	double jetPt_JER = recojetPT;
	double jetPt_JERup = recojetPT;
	double jetPt_JERdown = recojetPT;

    double diff_recojet_genjet = recojetPT - genjetPT;
    if ((era == "2011" && samplename == "ttH120") || (era == "2012_52x" && samplename == "ttH120_FastSim")) {
      diff_recojet_genjet += diff_FullSim_FastSim; }

	if( genjetPT>10. ){
		jetPt_JER = std::max( 0., genjetPT + scale_JER * ( diff_recojet_genjet ) );
		jetPt_JERup = std::max( 0., genjetPT + scale_JERup * ( diff_recojet_genjet ) );
		jetPt_JERdown = std::max( 0., genjetPT + scale_JERdown * ( diff_recojet_genjet ) );
	}

	if( returnType==1 )       factor = jetPt_JERup/recojetPT;
	else if( returnType==-1 ) factor = jetPt_JERdown/recojetPT;
	else                      factor = jetPt_JER/recojetPT;

	if( !(genjetPT>10.) ) factor = 1.;

	return factor;
}



void BEANhelper::getSp(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, float &aplanarity, float &sphericity) {
	//
	// Aplanarity and sphericity
	//

	int nJets = int(jets.size());

	float mxx = lepton.Px()*lepton.Px() + met.Px()*met.Px();
	float myy = lepton.Py()*lepton.Py() + met.Py()*met.Py();
	float mzz = lepton.Pz()*lepton.Pz() + met.Pz()*met.Pz();
	float mxy = lepton.Px()*lepton.Py() + met.Px()*met.Py();
	float mxz = lepton.Px()*lepton.Pz() + met.Px()*met.Pz();
	float myz = lepton.Py()*lepton.Pz() + met.Px()*met.Pz();

	for (int i=0; i<nJets; i++) {
		mxx += jets[i].Px()*jets[i].Px();
		myy += jets[i].Py()*jets[i].Py();
		mzz += jets[i].Pz()*jets[i].Pz();
		mxy += jets[i].Px()*jets[i].Py();
		mxz += jets[i].Px()*jets[i].Pz();
		myz += jets[i].Py()*jets[i].Pz();		
	}
	float sum = mxx + myy + mzz;
	mxx /= sum;
	myy /= sum;
	mzz /= sum;
	mxy /= sum;
	mxz /= sum;
	myz /= sum;

	TMatrix tensor(3,3);
	tensor(0,0) = mxx;
	tensor(1,1) = myy;
	tensor(2,2) = mzz;
	tensor(0,1) = mxy;
	tensor(1,0) = mxy;
	tensor(0,2) = mxz;
	tensor(2,0) = mxz;
	tensor(1,2) = myz;
	tensor(2,1) = myz;
	TVector eigenval(3);
	tensor.EigenVectors(eigenval);

	sphericity = 3.0*(eigenval(1)+eigenval(2))/2.0;
	aplanarity = 3.0*eigenval(2)/2.0;

	return;
}


void BEANhelper::getFox(vecTLorentzVector jets, float &h0, float &h1, float &h2, float &h3, float &h4) {
	//
	// Aplanarity and sphericity
	//

	int visObjects = int(jets.size());

	float eVis = 0.0;
	for (int i=0; i<visObjects; i++) {
		eVis += jets[i].E();
	}

	h0 = 0.0;
	h1 = 0.0;
	h2 = 0.0;
	h3 = 0.0;
	h4 = 0.0;
	for (int i=0; i<visObjects-1; i++) {
		for (int j=i+1; j<visObjects; j++) {
			float costh = cos(jets[i].Angle(jets[j].Vect()));
			float p0 = 1.0;
			float p1 = costh;
			float p2 = 0.5*(3.0*costh*costh - 1.0);
			float p3 = 0.5*(5.0*costh*costh - 3.0*costh);
			float p4 = 0.125*(35.0*costh*costh*costh*costh - 30.0*costh*costh + 3.0);
			float pipj = jets[i].P()*jets[j].P();
			h0 += (pipj/(eVis*eVis))*p0;
			h1 += (pipj/(eVis*eVis))*p1;
			h2 += (pipj/(eVis*eVis))*p2;
			h3 += (pipj/(eVis*eVis))*p3;
			h4 += (pipj/(eVis*eVis))*p4;
		}
	}

	return;
}


vdouble BEANhelper::getEffSF( int returnType, double jetPt, double jetEta, double jetId ){

	vdouble result;

	if (isData) {
		result.clear();
		result.push_back(1.0);
		result.push_back(1.0);
		return result;
	}

	double m_type = 0.;
	if( returnType==-1 )      m_type = -1.;
	else if( returnType==1 )  m_type = 1.;
	else                      m_type = 0.;

	float ptmin[] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
	float ptmax[] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};

	double SFb_error[] = {
		0.0295675,
		0.0295095,
		0.0210867,
		0.0219349,
		0.0227033,
		0.0204062,
		0.0185857,
		0.0256242,
		0.0383341,
		0.0409675,
		0.0420284,
		0.0541299,
		0.0578761,
		0.0655432 };

	double pt  = jetPt;
	double eta = jetEta;
	int flavor = jetId;

	double threshold = 670;
	pt = ( pt>threshold ) ? threshold-0.0000001 : pt;

	double absEta = fabs(eta);

	int use_bin=-1;
	for( int p=0; p<14; p++ ){
		if( pt>ptmin[p] && pt<ptmax[p] ){
			use_bin = p; break;
		}
	}
	if( use_bin<0 ) std::cout << "   ERROR!! use_bin < 0 " << std::endl;


	double SFb = 0.6981*((1.+(0.414063*pt))/(1.+(0.300155*pt)));

	double SFc = SFb;

	SFb = SFb + m_type * SFb_error[use_bin];
	SFc = SFc + m_type * 2* SFb_error[use_bin];

	double SFl = 1.;
	if( returnType==-2 ){ // min
		if( absEta < 0.8 )                      SFl = ((0.972455+(7.51396e-06*pt))+(4.91857e-07*(pt*pt)))+(-1.47661e-09*(pt*(pt*pt)));
		else if( absEta < 1.6 && absEta > 0.8 ) SFl = ((1.02055+(-0.000378856*pt))+(1.49029e-06*(pt*pt)))+(-1.74966e-09*(pt*(pt*pt)));
		else if( absEta < 2.4 && absEta > 1.6 ) SFl = ((0.983476+(-0.000607242*pt))+(3.17997e-06*(pt*pt)))+(-4.01242e-09*(pt*(pt*pt)));
	}
	else if( returnType==2 ){ // max
		if( absEta < 0.8 )                      SFl = ((1.15116+(0.00122657*pt))+(-3.63826e-06*(pt*pt)))+(2.08242e-09*(pt*(pt*pt)));
		else if( absEta < 1.6 && absEta > 0.8 ) SFl = ((1.20146+(0.000359543*pt))+(-1.12866e-06*(pt*pt)))+(6.59918e-10*(pt*(pt*pt)));
		else if( absEta < 2.4 && absEta > 1.6 ) SFl = ((1.18654+(-0.000795808*pt))+(3.69226e-06*(pt*pt)))+(-4.22347e-09*(pt*(pt*pt)));
	}
	else { // mean
		if( absEta < 0.8 )                      SFl = ((1.06182+(0.000617034*pt))+(-1.5732e-06*(pt*pt)))+(3.02909e-10*(pt*(pt*pt)));
		else if( absEta < 1.6 && absEta > 0.8 ) SFl = ((1.111+(-9.64191e-06*pt))+(1.80811e-07*(pt*pt)))+(-5.44868e-10*(pt*(pt*pt)));
		else if( absEta < 2.4 && absEta > 1.6 ) SFl = ((1.08498+(-0.000701422*pt))+(3.43612e-06*(pt*pt)))+(-4.11794e-09*(pt*(pt*pt)));
	}

	double SF=1;

	double tagEff=0;
	if( abs(flavor)==5 ){
		tagEff = h_b_eff_->GetBinContent( h_b_eff_->FindBin( pt, absEta ) );
		SF = SFb;
	}
	else if( abs(flavor)==4 ){
		tagEff = h_c_eff_->GetBinContent( h_c_eff_->FindBin( pt, absEta ) );
		SF = SFc;
	}
	else if( abs(flavor)==1 || abs(flavor)==2 || abs(flavor)==3 || abs(flavor)==21 ){
		tagEff = h_l_eff_->GetBinContent( h_l_eff_->FindBin( pt, absEta ) );
		SF = SFl;
	}
	else {
		tagEff = h_o_eff_->GetBinContent( h_o_eff_->FindBin( pt, absEta ) );
		SF = SFl;
	}


	result.clear();
	result.push_back(tagEff);
	result.push_back(SF);

	return result;
}


bool BEANhelper::ttPlusHeavyKeepEvent( const BNmcparticleCollection& iMCparticles,
                                       const BNjetCollection& iJets ) {

  CheckSetUp();
  string samplename = GetSampleName();

  // validate input
  bool validInput = false;
  if (samplename == "ttbar" || samplename == "ttbar_jj" || samplename == "ttbar_lj" || samplename == "ttbar_ll") validInput = true;
  if (samplename == "ttbar_cc" || samplename == "ttbar_cc_jj" || samplename == "ttbar_cc_lj" || samplename == "ttbar_cc_ll") validInput = true;
  if (samplename == "ttbar_bb" || samplename == "ttbar_bb_jj" || samplename == "ttbar_bb_lj" || samplename == "ttbar_bb_ll") validInput = true;

  if (!validInput ){
    cout << "ttPlusHeavyKeepEvent: could not recognize samplename " << samplename <<"... failing" <<endl;
  }

  bool validEra = false;
  if (era == "2011" ) validEra = true;
  if (era == "2012_52x" || era == "2012_53x" ) validEra = true;

  if (!validEra ){
    cout << "ttPlusHeavyKeepEvent: could not recognize era " << era <<"... failing" <<endl;
  }

  bool keepEvent = false;
  bool debug_ = false;

  bool isWtoCS = false;

  if (debug_) cout << "Num MC particles = " << iMCparticles.size() << std::endl
                   << "Num pfjets " <<  int(iJets.size()) << std::endl;

  // check to see if the event has a c with a parent W

  for( unsigned i=0; i< iMCparticles.size(); i++ ){
    int id = iMCparticles.at(i).id;
    int motherID = iMCparticles.at(i).motherId;
    int grandMotherID = iMCparticles.at(i).grandMotherId;

    if (debug_) std::cout << "Particle " << i << " is " << id << ", has mother " << motherID << " and grandmother " << grandMotherID << std::endl;

    if (debug_) cout <<" Particle " << i << " has id " << id << endl;

    if (era == "2011") {
      if( abs(id)==4  && abs(motherID)==24 ){
        isWtoCS = true;
        break;
      }
    }
    else if (era == "2012_52x" || era == "2012_53x") {
      int daughter0ID = iMCparticles.at(i).daughter0Id;
      int daughter1ID = iMCparticles.at(i).daughter1Id;

      if( abs(id)==24 && ((abs(daughter0ID)==3 && abs(daughter1ID)==4) || (abs(daughter0ID)==4 && abs(daughter1ID)==3)) ){
        isWtoCS = true;
        break;
      }
    }
    else {
      assert (era == "either 2012_52x, 2012_53x, or 2011");
    }    
  }

  bool isBBbarEvent = false;
  bool isCCbarEvent = false;


  bool gotB = false;
  bool gotC = false;

  int numBmomB=0;
  int numBmomT=0;
  int numBbarmomBbar=0;
  int numBbarmomTbar=0;
  int numCmomC=0;
  int numCbarmomCbar=0;

  if (debug_) cout << "Starting loop over pf jet parton ids to see if you have a b" <<endl;

  for( int i=0; i<int(iJets.size()); i++ ){



    int id = iJets.at(i).genPartonId;
    if( id==-99 ) continue;
    int motherID = iJets.at(i).genPartonMotherId;
    int mother0ID =  iJets.at(i).genPartonMother0Id;
    int mother1ID =  iJets.at(i).genPartonMother1Id;



    // check to see if pf jets is from a  b/c and mother is a gluon
    // or, if mother is some light quark

    if (era == "2011") {
      if( abs(id)==5 && ( motherID==21 || abs(motherID)<abs(id) ) ) gotB=true;
      if( abs(id)==4 && ( motherID==21 || abs(motherID)<abs(id) ) ) gotC=true;
    } else if (era == "2012_52x" || era == "2012_53x") {
      if( abs(id)==5 && abs(mother0ID) != 6  && abs(mother1ID) != 6 ) gotB=true;
      // basically, as long as you didn't come from a W, then you are a tt+cc
      if( abs(id)==4 && (abs(mother0ID)==21 || abs(mother1ID) == 21 || abs(mother0ID) < abs(id) || abs(mother1ID) < abs(id))
          && abs(mother0ID)!=24 && abs(mother1ID)!=24 ) gotC=true;
    }
    else {
      assert (era == "either 2012_52x, 2012_53x, or 2011");
    }
    

    if (debug_) std::cout << "Jet index " << i << " is generator id " << id
                          << ", has mother " << motherID << ", mother0ID = " << mother0ID << ", mother1ID = " << mother1ID  << std::endl;

    if (debug_) std::cout << "----------------> Got B = " << gotB << endl;

    // if things are their own mother,
    // where are they from? Does this mean stable?
    if( id==5  && motherID==id ) numBmomB++;
    if( id==-5 && motherID==id ) numBbarmomBbar++;

    if( id==4  && motherID==id ) numCmomC++;
    if( id==-4 && motherID==id ) numCbarmomCbar++;

    if( id==5  && motherID==6  ) numBmomT++;
    if( id==-5 && motherID==-6 ) numBbarmomTbar++;

  }

  // if at least one b from b & one b from t, or if CC, and your jet was not b
  if( ((numBmomB>=1 && numBmomT>=1) || (numBbarmomBbar>=1 && numBbarmomTbar>=1)) && !gotB ){
    if (debug_) std::cout << "No sign of a b jet, but looping over jets again to check"
                          <<std::endl;
    // for each jet that is  b from b
    for( int i=0; i<int(iJets.size()); i++ ){
      if (debug_) cout << "LOOP: i = " << i << endl;
      int id0 = iJets.at(i).genPartonId;
      int motherID0 = iJets.at(i).genPartonMotherId;
      if( !(abs(id0)==5 && motherID0==id0) ) continue;

      if (debug_) std::cout << "Jet index " << i  << " is a bjet, let us see that it is not from top" <<std::endl;
      // for each jet that is b from t
      for( int j=0; j<int(iJets.size()); j++ ){
        if (debug_) cout << "LOOP: j = " << j << endl;
        int id1 = iJets.at(j).genPartonId;
        int motherID1 = iJets.at(j).genPartonMotherId;
        if (debug_) std::cout << "LOOP: id0 = " << id0 << ", motherID0 = " << motherID0
                              << ", id1 = " << id1 << ", motherID1 = " << motherID1 << endl
                              << "continue? = " << !(id1==id0 && abs(motherID1)==6) << endl;
        if( !(id1==id0 && abs(motherID1)==6) ) continue;
        if (debug_) std::cout << "You didn't skip this event!" << endl;
        // if delta r between b from b and b from t is big enough, then b in final state is OK
        double dR = reco::deltaR(iJets.at(i).genPartonEta,
                                 iJets.at(i).genPartonPhi,
                                 iJets.at(j).genPartonEta,
                                 iJets.at(j).genPartonPhi);
        if (debug_) std::cout << "dR = " << dR << endl;
        if (debug_) std::cout << "gotB = " << gotB << endl;
        if( dR>0.3 ){
          gotB = true;
          if (debug_) std::cout << "Found something with dR > 0.3... now gotB = " << gotB << endl;
          break;
        }
        if (debug_) std::cout << "SECOND: b not from top with dR = " << dR << " gotB = " << gotB << endl;
      }
      if( gotB ) break;
    }
  }

  if( (numCmomC>=1 || numCbarmomCbar>=1) && !isWtoCS ){
    gotC = true;
  }

  if( gotB ) isBBbarEvent = true;
  else if( gotC ) isCCbarEvent = true;



  if( (samplename == "ttbar" || samplename == "ttbar_jj" || samplename == "ttbar_lj" || samplename == "ttbar_ll") && !isBBbarEvent && !isCCbarEvent ) keepEvent = true;
  else if( (samplename == "ttbar_bb" || samplename == "ttbar_bb_jj" || samplename == "ttbar_bb_lj" || samplename == "ttbar_bb_ll") &&  isBBbarEvent && !isCCbarEvent ) keepEvent = true;
  else if( (samplename == "ttbar_cc" || samplename == "ttbar_cc_jj" || samplename == "ttbar_cc_lj" || samplename == "ttbar_cc_ll")  && !isBBbarEvent && isCCbarEvent  ) keepEvent = true;

  if (debug_) cout << "Filter result = " << keepEvent << endl
                   << "isBBbarEvent = " << isBBbarEvent << endl
                   << "isCCbarEvent = " << isCCbarEvent << endl
                   << "... will we skip this? " << (!keepEvent) << endl;

  return keepEvent;

}
