#ifndef INC_BEANSUTILITIES
#define INC_BEANSUTILITIES
///////////////////////////////////////////////////////////////////////////////
// File: BEANsUtilities.h
// 
// Purpose:  Commonly used BEAN functions
//
// Created:   15-JUL-2012  Darren Puigh
// History:   new!
// Modified:  
//
///////////////////////////////////////////////////////////////////////////////
// Dependencies (#includes)

#include <iostream>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <algorithm>
#include "TVector.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "Math/Interpolator.h"

#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

//Headers for the data items
#include "ProductArea/BNcollections/interface/BNelectron.h"
#include "ProductArea/BNcollections/interface/BNevent.h"
#include "ProductArea/BNcollections/interface/BNjet.h"
#include "ProductArea/BNcollections/interface/BNmcparticle.h"
#include "ProductArea/BNcollections/interface/BNmet.h"
#include "ProductArea/BNcollections/interface/BNmuon.h"
#include "ProductArea/BNcollections/interface/BNphoton.h"
#include "ProductArea/BNcollections/interface/BNsupercluster.h"
#include "ProductArea/BNcollections/interface/BNtrack.h"
#include "ProductArea/BNcollections/interface/BNtrigger.h"
#include "ProductArea/BNcollections/interface/BNskimbits.h"
#include "ProductArea/BNcollections/interface/BNtrigobj.h"
#include "ProductArea/BNcollections/interface/BNprimaryvertex.h"

#include "NtupleMaker/BEANmaker/interface/BtagWeight.h"
#include "NtupleMaker/BEANmaker/interface/BTagReshaping.h"

#include "DataFormats/Math/interface/deltaR.h"

#endif


typedef std::vector< TLorentzVector > vecTLorentzVector;
typedef std::vector<std::vector<double> > vvdouble;
typedef std::vector<std::vector<std::string> > vvstring;
typedef std::vector<std::vector<int> > vvint;
typedef std::vector<std::string> vstring;
typedef std::vector<double> vdouble;
typedef std::vector<int> vint;

typedef BNeventCollection::const_iterator         EventIter;
typedef BNelectronCollection::const_iterator      EleIter;
typedef BNjetCollection::const_iterator           JetIter;
typedef BNmcparticleCollection::const_iterator    MCparIter;
typedef BNmetCollection::const_iterator           MetIter;
typedef BNmuonCollection::const_iterator          MuonIter;
typedef BNphotonCollection::const_iterator        PhotonIter;
typedef BNprimaryvertexCollection::const_iterator PVIter;
typedef BNskimbitsCollection::const_iterator      SkimBitIter;
typedef BNsuperclusterCollection::const_iterator  SCIter;
typedef BNtrackCollection::const_iterator         TrackIter;
typedef BNtriggerCollection::const_iterator       TrigIter;
typedef BNtrigobjCollection::const_iterator       TrigObjIter;

char * my_pPath = getenv ("CMSSW_BASE");
std::string my_base_dir(my_pPath);
std::string str_eff_file_7TeV = my_base_dir + "/src/NtupleMaker/BEANmaker/data/mc_btag_efficiency_7TeV.root";
std::string str_eff_file_8TeV = my_base_dir + "/src/NtupleMaker/BEANmaker/data/mc_btag_efficiency_8TeV.root";
std::string str_pu_file_7TeV  = my_base_dir + "/src/NtupleMaker/BEANmaker/data/pu_distributions_7TeV.root";
std::string str_pu_file_8TeV  = my_base_dir + "/src/NtupleMaker/BEANmaker/data/pu_distributions_8TeV.root";
std::string str_lep_file_7TeV  = my_base_dir + "/src/NtupleMaker/BEANmaker/data/lepton_SF_8TeV.root";
std::string str_lep_file_8TeV  = my_base_dir + "/src/NtupleMaker/BEANmaker/data/lepton_SF_8TeV.root";

std::string str_csv_file_8TeV = str_eff_file_8TeV;


BTagShapeInterface *sh_;
BTagShapeInterface *sh_hfSFUp_;
BTagShapeInterface *sh_hfSFDown_;
BTagShapeInterface *sh_lfSFUp_;
BTagShapeInterface *sh_lfSFDown_;


TH2D* h_b_eff_;
TH2D* h_c_eff_;
TH2D* h_l_eff_;
TH2D* h_o_eff_;

TH1D* h_PU_ratio_;
TH1D* h_PUup_ratio_;
TH1D* h_PUdown_ratio_;

TH2D* h_ele_SF_;
TH2D* h_mu_SF_;


bool isFastSim_ = false;
bool isLJ_ = false;
std::string era_ = "2011";


using namespace std;

namespace BEANs{
  const double PI=2.0*acos(0.);
  const double TWOPI=2.0*PI;
  const float ETA_LIMIT=15.0;
  const float EPSILON=1.E-10;

  void setMCsample( int insample=2500, bool is8TeV=true, bool isLJ=true, std::string dset="" );

  void electronSelector( const BNelectronCollection &electrons, bool isLJ, std::string era, vint &tightElectrons, vint &looseElectrons, vdouble &tightElectronSF, vdouble &looseElectronSF );
  void muonSelector( const BNmuonCollection &muons, bool isLJ, std::string era, vint &tightMuons, vint &looseMuons, vdouble &tightMuonSF, vdouble &looseMuonSF );
  void jetSelector( const BNjetCollection &pfjets, std::string sysType, std::string era, vint &tightJets, vint &tagJets, vint &untagJets, 
		    std::vector<BTagWeight::JetInfo> &myjetinfo, double csvCut = 0.679 );
  void jetSelectorV2( const BNjetCollection &pfjets, std::string sysType, std::string era, vint &tightJets, vint &tagJets, vint &untagJets, double csvCut = 0.679 );

  void getPUwgt( double input_numPU, double &PU_scale, double &PUup_scale, double &PUdown_scale ); 

  vdouble getEffSF( int returnType, double jetPts, double jetEtas, double jetIds, std::string era);
  double getJERfactor( int returnType, double jetAbsETA, double genjetPT, double recojetPT );
  void getSp(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, float &aplanarity, float &sphericity);
  void getFox(vecTLorentzVector jets,float &h0, float &h1, float &h2, float &h3, float &h4);	
  void getFox_mod2(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, double HT,
		   float &h0_mod2, float &h1_mod2, float &h2_mod2, float &h3_mod2, float &h4_mod2,  float &h5_mod2,  
		   float &h6_mod2, float &h7_mod2, float &h8_mod2, float &h9_mod2, float &h10_mod2 );

  bool ttPlusHeavyKeepEvent( BNmcparticleCollection const &mcparticles,
                        BNjetCollection const &pfjets,
                        TString ttbarType,
                        TString era);

  double reshape_csv( double eta, double pt, double csv, int flavor, std::string sysType );

}



void BEANs::setMCsample( int insample, bool is8TeV, bool isLJ, std::string dset ){

  bool debug = false;
  
  std::string input_eff_file = str_eff_file_7TeV;
  std::string input_lep_file = str_lep_file_7TeV;
  std::string input_pu_file  = str_pu_file_7TeV;
  std::string com_suffix = "_7TeV";
  if( is8TeV ){
    input_eff_file = str_eff_file_8TeV;
    input_lep_file = str_lep_file_8TeV;
    input_pu_file  = str_pu_file_8TeV;
    com_suffix = "_8TeV";
    era_ = "2012";
  }

  if( isLJ ) isLJ_ = true;

  if (debug)
    cout << "setMCsample: Opening eff file " << input_eff_file
         << ", pu file = " << input_pu_file << endl;

  TFile *f_tag_eff_ = new TFile(input_eff_file.c_str());

  if (f_tag_eff_->IsZombie()){
    cout << "Oops! Tried to open file " << input_eff_file
         << ", but it was a zombie. Crashing" << endl;
    assert (f_tag_eff_->IsZombie() == false);
  }

  std::string samplename = "ttbar";
  if( insample==2300 || insample==2310 ) samplename = "zjets";
  else if( insample==2400 ) samplename = "wjets";
  else if( insample==2500 ) samplename = "ttbar";
  else if( insample==2510 ) samplename = "ttbar_scaleup";
  else if( insample==2511 ) samplename = "ttbar_scaledown";
  else if( insample==2523 ) samplename = "ttbarZ";
  else if( insample==2524 ) samplename = "ttbarW";
  else if( insample==2600 ) samplename = "singlet";
  else if( insample==2700 ) samplename = "ww";
  else if( insample==2701 ) samplename = "wz";
  else if( insample==2702 ) samplename = "zz";
  else if( insample>=100 && insample<=140 ) samplename = "ttH120";

  if( is8TeV ){
    // 12Aug2012 - Include all samples - JGWood
    if( insample==2500      ) samplename = "ttbar";
    else if( insample==2523 ) samplename = "ttbarZ";
    else if( insample==2524 ) samplename = "ttbarW";
    else if( insample==2600 ) samplename = "t_schannel";
    else if( insample==2602 ) samplename = "t_tchannel";
    else if( insample==2504 ) samplename = "t_tWchannel";
    else if( insample==2501 ) samplename = "tbar_schannel";
    else if( insample==2503 ) samplename = "tbar_tchannel";
    else if( insample==2505 ) samplename = "tbar_tWchannel";
    else if( insample==2800 || insample==2850 ) samplename = "zjets";
    else if( insample==2400 ) samplename = "wjets";
    else if( insample==2700 ) samplename = "ww";
    else if( insample==2701 ) samplename = "wz";
    else if( insample==2702 ) samplename = "zz";
    else if( insample>=8000 && insample<9000  ) samplename = "ttH120_FullSim";
    else if( insample>=9000 && insample<10000 ) samplename = "ttH120_FastSim";
  }


  if( insample>=9000 && insample<10000 ) isFastSim_ = true;

  if (debug)
    cout << "setMCSample: Looking for histrograms with names like: "
         << std::string( samplename + com_suffix + "_jet_pt_eta_b_eff")
         << endl;

  
  h_b_eff_ = (TH2D*)f_tag_eff_->Get(std::string( samplename + com_suffix + "_jet_pt_eta_b_eff" ).c_str());
  h_c_eff_ = (TH2D*)f_tag_eff_->Get(std::string( samplename + com_suffix + "_jet_pt_eta_c_eff" ).c_str());
  h_l_eff_ = (TH2D*)f_tag_eff_->Get(std::string( samplename + com_suffix + "_jet_pt_eta_l_eff" ).c_str());
  h_o_eff_ = (TH2D*)f_tag_eff_->Get(std::string( samplename + com_suffix + "_jet_pt_eta_o_eff" ).c_str());
  
  bool bHistoOK =  (h_b_eff_ != 0);
  bool cHistoOK =  (h_c_eff_ != 0);
  bool lHistoOK =  (h_l_eff_ != 0);
  bool oHistoOK =  (h_o_eff_ != 0);

  if (debug)
    cout << "setMCSample: bHistoOK = " << bHistoOK << ", cHistoOK = " << cHistoOK << ", lHistoOK = "
         << lHistoOK << ", oHistoOK = " << oHistoOK << endl;
  
  if (!bHistoOK || !cHistoOK || !lHistoOK || !oHistoOK){
    cout << "setMCSample: Error. We are missing one of the required btag eff histograms. "     
         << endl
         << "Wanted histos with names like: "
         << std::string( samplename + com_suffix + "_jet_pt_eta_b_eff")
         << endl;
    assert ( bHistoOK && cHistoOK && lHistoOK && oHistoOK);
  }



  TFile *f_pu_ = new TFile(input_pu_file.c_str());

  TH1D* h_pu_data;
  TH1D* h_pu_data_up;
  TH1D* h_pu_data_down;
  TH1D* h_pu_mc;

  if( is8TeV ){
    h_pu_data      = (TH1D*)f_pu_->Get((std::string("pileup_8TeV_69300xSec")).c_str());
    h_pu_data_up   = (TH1D*)f_pu_->Get((std::string("pileup_8TeV_71795xSec")).c_str());
    h_pu_data_down = (TH1D*)f_pu_->Get((std::string("pileup_8TeV_66805xSec")).c_str());

    std::string mc_pu_input = "Summer2012";
    if( insample==2300 || insample==2400 || insample==2500 || 
	(insample>=8000 && insample<9000) || 
	(insample>=9000 && insample<10000) ) mc_pu_input = std::string(samplename + "_Summer2012");
    else if( insample==2523 || insample==2524 ) mc_pu_input = "ttZorW_Summer2012";

    h_pu_mc = (TH1D*)f_pu_->Get((std::string(mc_pu_input + "_pileup_8TeV")).c_str());
  }
  else{
    if( !(dset.find("SingleMu")!=std::string::npos || dset.find("ElectronHad")!=std::string::npos) ) dset = "SingleMu";

    if( (insample>=100 && insample<=140) || (insample==2523) || (insample==2524) ){
      h_pu_data      = (TH1D*)f_pu_->Get((std::string("pileup_7TeV_" + dset + "_68000_observed")).c_str());
      h_pu_data_up   = (TH1D*)f_pu_->Get((std::string("pileup_7TeV_" + dset + "_73440_observed")).c_str());
      h_pu_data_down = (TH1D*)f_pu_->Get((std::string("pileup_7TeV_" + dset + "_62560_observed")).c_str());

      h_pu_mc = (TH1D*)f_pu_->Get("ttH_7TeV_numGenPV");
    }
    else{
      h_pu_data      = (TH1D*)f_pu_->Get((std::string("pileup_7TeV_" + dset + "_68000_true")).c_str());
      h_pu_data_up   = (TH1D*)f_pu_->Get((std::string("pileup_7TeV_" + dset + "_73440_true")).c_str());
      h_pu_data_down = (TH1D*)f_pu_->Get((std::string("pileup_7TeV_" + dset + "_62560_true")).c_str());

      h_pu_mc = (TH1D*)f_pu_->Get("F2011exp_7TeV");
    }
  }

  h_pu_data->Scale( 1./h_pu_data->Integral() );
  h_pu_data_up->Scale( 1./h_pu_data_up->Integral() );
  h_pu_data_down->Scale( 1./h_pu_data_down->Integral() );

  h_pu_mc->Scale( 1./h_pu_mc->Integral() );

  h_PU_ratio_     = (TH1D*)h_pu_data->Clone();
  h_PUup_ratio_   = (TH1D*)h_pu_data_up->Clone();
  h_PUdown_ratio_ = (TH1D*)h_pu_data_down->Clone();

  h_PU_ratio_->Divide( h_pu_mc );
  h_PUup_ratio_->Divide( h_pu_mc );
  h_PUdown_ratio_->Divide( h_pu_mc );


  TFile *f_lep_ = new TFile(input_lep_file.c_str());
  if( isLJ_ ){
    h_ele_SF_ = (TH2D*)f_lep_->Get(std::string( "ele_pt_eta_full_id_iso_hlt_8TeV" ).c_str());
    h_mu_SF_  = (TH2D*)f_lep_->Get(std::string( "mu_pt_eta_full_id_iso_hlt_8TeV" ).c_str());
  }
  else {
    h_ele_SF_ = (TH2D*)f_lep_->Get(std::string( "ele_pt_eta_full_id_iso_8TeV" ).c_str());
    h_mu_SF_  = (TH2D*)f_lep_->Get(std::string( "mu_pt_eta_full_id_iso_8TeV" ).c_str());
  }



  sh_ = new BTagShapeInterface(std::string(samplename + com_suffix),str_csv_file_8TeV.c_str(),0,0);
  sh_hfSFUp_ = new BTagShapeInterface(std::string(samplename + com_suffix),str_csv_file_8TeV.c_str(),1.5,0);
  sh_hfSFDown_ = new BTagShapeInterface(std::string(samplename + com_suffix),str_csv_file_8TeV.c_str(),-1.5,0);
  sh_lfSFUp_ = new BTagShapeInterface(std::string(samplename + com_suffix),str_csv_file_8TeV.c_str(),0,1);
  sh_lfSFDown_ = new BTagShapeInterface(std::string(samplename + com_suffix),str_csv_file_8TeV.c_str(),0,-1);


}


void BEANs::getPUwgt( double input_numPU, double &PU_scale, double &PUup_scale, double &PUdown_scale ){

  PU_scale     = h_PU_ratio_->GetBinContent( h_PU_ratio_->FindBin( input_numPU ) );
  PUup_scale   = h_PUup_ratio_->GetBinContent( h_PUup_ratio_->FindBin( input_numPU ) );
  PUdown_scale = h_PUdown_ratio_->GetBinContent( h_PUdown_ratio_->FindBin( input_numPU ) );

}


/////////
///
/// Electrons
///
////////
void BEANs::electronSelector( const BNelectronCollection &electrons, bool isLJ, std::string era, vint &tightElectrons, vint &looseElectrons, vdouble &tightElectronSF, vdouble &looseElectronSF ){

  tightElectrons.clear();
  looseElectrons.clear();
  tightElectronSF.clear();
  looseElectronSF.clear();

  era_ = era;
  isLJ_ = isLJ;

  bool is2011 = ( era_.find("2011")!=std::string::npos );
  double tightPt = ( isLJ_ ) ? 30. : 20.;
  double loosePt = 10.;

  if( is2011 ){
    for( int i=0; i<int(electrons.size()); i++ ){
      double eleSCEta = electrons.at(i).scEta;
      double absSCeta = fabs(eleSCEta);
      double eleEta = electrons.at(i).eta;
      double elePt = electrons.at(i).pt;

      bool isCrack = ( (absSCeta>1.4442) && (absSCeta<1.5660) );

      bool kin = ( (elePt>loosePt) && !isCrack && fabs(eleEta)<2.5 );

      if( !kin ) continue;

      double chargedHadronIso = electrons.at(i).chargedHadronIso;
      double neutralHadronIso = electrons.at(i).neutralHadronIso;
      double photonIso = electrons.at(i).photonIso;

      double relIso = ( chargedHadronIso + neutralHadronIso + photonIso ) * 1./elePt;

      bool looseIso = ( relIso < 0.2 );
      bool tightIso = ( relIso < 0.1 );

      int eidHyperTight1MC = electrons.at(i).eidHyperTight1MC;
      bool eidHyperTight1MC_dec = ( (eidHyperTight1MC & 1)==1 );

      int eidTight = electrons.at(i).eidTight;
      bool eidTight_dec = ( (eidTight & 1)==1 );

      bool eid = isLJ_ ? eidHyperTight1MC_dec : eidTight_dec;

      bool d0 = ( fabs(electrons.at(i).correctedD0) < 0.02 );
      bool dZ = ( fabs(electrons.at(i).correctedDZ) < 1. );

      bool dist  = ( fabs(electrons.at(i).dist)<0.02 );
      bool dcot  = ( fabs(electrons.at(i).dcot)<0.02 );
      bool nlost = ( electrons.at(i).numberOfLostHits<1 );
      bool notConv = ( !(dist && dcot) && nlost );

      bool id = ( eid && d0 && dZ && notConv );

      if( kin && looseIso ){
        if( ((elePt>tightPt) && id && tightIso) ){
	  tightElectrons.push_back(i);
	  tightElectronSF.push_back(1.);
	}
        else{
	  looseElectrons.push_back(i);
	  looseElectronSF.push_back(1.);
	}
      }
    }// end electron loop

  } // end if 2011
  else{ // default is 2012 selection
    for( int i=0; i<int(electrons.size()); i++ ){

      double eleSCEta = electrons.at(i).scEta;
      double absSCeta = fabs(eleSCEta);
      double eleEta = electrons.at(i).eta;
      double elePt = electrons.at(i).pt;

      bool isCrack = ( (absSCeta>1.4442) && (absSCeta<1.5660) );

      bool kin = ( (elePt>loosePt) && !isCrack && fabs(eleEta)<2.5 );

      if( !kin ) continue;

      double chargedHadronIso = electrons.at(i).chargedHadronIso;
      double neutralHadronIso = electrons.at(i).neutralHadronIso;
      double photonIso = electrons.at(i).photonIso;

      double AEffDr03 = electrons.at(i).AEffDr03;
      double rhoPrime = electrons.at(i).rhoPrime;

      double relIso_rho   = ( chargedHadronIso + max(0.0, neutralHadronIso + photonIso - AEffDr03*rhoPrime) ) * 1./elePt;

      bool looseIso = ( relIso_rho < 0.2 );
      bool tightIso = ( relIso_rho < 0.1 );

      double mvaID = electrons.at(i).mvaTrigV0;
      bool passMVAId = ( mvaID>0.0 );

      bool d02 = ( fabs(electrons.at(i).correctedD0Vertex) < 0.02 );
      bool d04 = ( fabs(electrons.at(i).correctedD0Vertex) < 0.04 );
      bool dZ = ( fabs(electrons.at(i).correctedDZ) < 1. );

      bool notConv = ( electrons.at(i).passConvVeto );

      bool id = ( passMVAId && d02 && dZ && notConv );

      if( kin && looseIso && passMVAId && d04 && notConv ){
        if( ((elePt>tightPt) && id && tightIso) ){
	  tightElectrons.push_back(i);
	  double usePT = std::min( elePt, 499. );
	  double useEta = ( eleEta>0. ) ? std::min( 2.39, eleEta ) : std::max( -2.39, eleEta );
	  double SF = h_ele_SF_->GetBinContent( h_ele_SF_->FindBin(usePT, useEta) );
	  tightElectronSF.push_back(SF);
	}
        else{
	  looseElectrons.push_back(i);
	  looseElectronSF.push_back(1.);
	}
      }
    }// end electron loop
  }
} //end electronSelector


/////////
///
/// Muons
///
////////
void BEANs::muonSelector( const BNmuonCollection &muons, bool isLJ, std::string era, vint &tightMuons, vint &looseMuons, vdouble &tightMuonSF, vdouble &looseMuonSF ){

  tightMuons.clear();
  looseMuons.clear();
  tightMuonSF.clear();
  looseMuonSF.clear();

  era_ = era;
  isLJ_ = isLJ;

  bool is2011 = ( era_.find("2011")!=std::string::npos );
  double tightPt = ( isLJ_ ) ? 30. : 20.;
  double loosePt = 10.;

  if( is2011 ){
    for( int i=0; i<int(muons.size()); i++ ){
      double muPt  = muons.at(i).pt;
      double muEta = muons.at(i).eta;
      double muAbsEta = fabs(muEta);

      bool kin = ( (muPt>loosePt) && (muAbsEta<2.4) );

      if( !kin ) continue;

      double chargedHadronIso = muons.at(i).chargedHadronIso;
      double neutralHadronIso = muons.at(i).neutralHadronIso;
      double photonIso = muons.at(i).photonIso;

      double relIso = ( chargedHadronIso + neutralHadronIso + photonIso ) * 1./muPt;

      bool looseIso = ( relIso<0.2 );
      bool tightIso = ( relIso<0.125 );

      bool isGlobalMuon = ( muons.at(i).isGlobalMuon==1 );
      bool isTrackerMuon = ( muons.at(i).isTrackerMuon==1 );
      bool isGlobalMuonPromptTight = ( muons.at(i).isGlobalMuonPromptTight==1 );

      bool numTrackHits = ( muons.at(i).numberOfValidTrackerHitsInnerTrack > 10 );
      bool numPixelHits = ( muons.at(i).pixelLayersWithMeasurement>0 );
      bool numberOfMatches = ( muons.at(i).numberOfMatchedStations>1 );

      bool passd0 = ( fabs(muons.at(i).correctedD0)<0.02 );
      bool passdz = ( fabs(muons.at(i).correctedDZ)<1. );

      bool id = ( isTrackerMuon && isGlobalMuonPromptTight && numTrackHits && numPixelHits && numberOfMatches && passd0 && passdz );

      if( kin && isGlobalMuon && looseIso ){
        if( ((muPt>tightPt) && (muAbsEta<2.1) && id && tightIso) ){
	  tightMuons.push_back(i);
	  tightMuonSF.push_back(1.);
	}
        else{
	  looseMuons.push_back(i);
	  looseMuonSF.push_back(1.);
	}
      }
    }// end muon loop
  } // end if 2011
  else{ // default is 2012 selection
    for( int i=0; i<int(muons.size()); i++ ){
      double muPt  = muons.at(i).pt;
      double muEta = muons.at(i).eta;
      double muAbsEta = fabs(muEta);

      bool kin = ( (muPt>loosePt) && (muAbsEta<2.5) );

      if( !kin ) continue;

      double pfIsoR04SumChargedHadronPt = muons.at(i).pfIsoR04SumChargedHadronPt;
      double pfIsoR04SumNeutralHadronEt = muons.at(i).pfIsoR04SumNeutralHadronEt;
      double pfIsoR04SumPhotonEt = muons.at(i).pfIsoR04SumPhotonEt;
      double pfIsoR04SumPUPt = muons.at(i).pfIsoR04SumPUPt;

      double relIso_dBeta = (pfIsoR04SumChargedHadronPt + max(0.0, pfIsoR04SumNeutralHadronEt + pfIsoR04SumPhotonEt - 0.5*pfIsoR04SumPUPt))/muPt;

      bool looseIso = ( relIso_dBeta<0.20 );
      bool tightIso = ( relIso_dBeta<0.12 );

      bool isPFmuon = ( muons.at(i).isPFMuon==1 );
//      bool isPFmuon = true; //Temporary hack for early 52x BEANs that lack this variable... (KPL)
      bool isGlobalMuon = ( muons.at(i).isGlobalMuon==1 );
      bool isTrackerMuon = ( muons.at(i).isTrackerMuon==1 );

      bool numberOfLayersWithMeasurement = ( muons.at(i).numberOfLayersWithMeasurement > 5 );
      bool numberOfValidMuonHits = (muons.at(i).numberOfValidMuonHits > 0);
      bool numberOfValidPixelHits = (muons.at(i).numberOfValidPixelHits > 0);
      bool numberOfMatchedStations = (muons.at(i).numberOfMatchedStations > 1);

      bool passd0 = ( fabs(muons.at(i).correctedD0Vertex) < 0.2 );
      bool dVzPVz = ( fabs(muons.at(i).dVzPVz) < 0.5 );

      bool normChi2 = ( muons.at(i).normalizedChi2 < 10) ;

      bool id = ( isGlobalMuon && normChi2 && passd0 && dVzPVz && 
		  numberOfLayersWithMeasurement && numberOfValidMuonHits && numberOfValidPixelHits && numberOfMatchedStations );

      if( kin && (isGlobalMuon || isTrackerMuon) && looseIso && isPFmuon ){
        if( ((muPt>tightPt) && (muAbsEta<2.1) && id && tightIso) ){
	  tightMuons.push_back(i);
	  double usePT = std::min( muPt, 499. );
	  double useEta = ( muEta>0. ) ? std::min( 2.09, muEta ) : std::max( -2.09, muEta );
	  double SF = h_mu_SF_->GetBinContent( h_mu_SF_->FindBin(usePT, useEta) );
	  tightMuonSF.push_back(SF);
	}
        else{
	  looseMuons.push_back(i);
	  looseMuonSF.push_back(1.);
	}
      }
    }// end muon loop
  }
} //end muonSelector


/////////
///
/// PFJets
///
////////
void BEANs::jetSelector( const BNjetCollection &pfjets, std::string sysType, std::string era, vint &tightJets, vint &tagJets, vint &untagJets, 
			 std::vector<BTagWeight::JetInfo> &myjetinfo, double csvCut ){

  tightJets.clear();
  tagJets.clear();
  untagJets.clear();
  myjetinfo.clear();

  for( int i=0; i<int(pfjets.size()); i++ ){
    double jetPt = pfjets.at(i).pt;
    double jetEta = pfjets.at(i).eta;
    double jetAbsEta = fabs(jetEta);

    bool eta = ( jetAbsEta<2.4 );
    bool jetId  = ( pfjets.at(i).jetIDLoose==1 );

    double factor=1;
    if( sysType.compare("data")!=0 ){
      double genJetPT = pfjets.at(i).genJetPT;
      if( sysType.compare("JERUp")==0 )        factor = getJERfactor(1,jetAbsEta,genJetPT,jetPt);
      else if( sysType.compare("JERDown")==0 ) factor = getJERfactor(-1,jetAbsEta,genJetPT,jetPt);
      else                                     factor = getJERfactor(0,jetAbsEta,genJetPT,jetPt);
    }
    jetPt *= factor;

    double unc = pfjets.at(i).JESunc;
    if( sysType.compare("JESUp")==0 )        jetPt *= (1. + unc);
    else if( sysType.compare("JESDown")==0 ) jetPt *= (1. - unc);

    double csv = pfjets.at(i).btagCombinedSecVertex;
    bool csvM = ( csv>csvCut );
    if( jetPt>30. && eta && jetId ){
      tightJets.push_back(i);
      if( csvM ) tagJets.push_back(i);
      else       untagJets.push_back(i);

      if( sysType.compare("data")!=0 ){
	int flavour = pfjets.at(i).flavour;
	std::vector<double> myEffSF;
	if( sysType.compare("hfSFUp")==0 )        myEffSF = BEANs::getEffSF( 1,  jetPt, jetEta, flavour, era_ );
	else if( sysType.compare("hfSFDown")==0 ) myEffSF = BEANs::getEffSF( -1, jetPt, jetEta, flavour, era_ );
	else if( sysType.compare("lfSFUp")==0 )   myEffSF = BEANs::getEffSF( 2,  jetPt, jetEta, flavour, era_ );
	else if( sysType.compare("lfSFDown")==0 ) myEffSF = BEANs::getEffSF( -2, jetPt, jetEta, flavour, era_ );
	else                                      myEffSF = BEANs::getEffSF( 0,  jetPt, jetEta, flavour, era_ );

	BTagWeight::JetInfo myjet( myEffSF[0], myEffSF[1] );
	myjetinfo.push_back(myjet);
      }
    }
  } // end loop over jets
}




/////////
///
/// PFJets
///
////////
void BEANs::jetSelectorV2( const BNjetCollection &pfjets, std::string sysType, std::string era, vint &tightJets, vint &tagJets, vint &untagJets, double csvCut ){

  tightJets.clear();
  tagJets.clear();
  untagJets.clear();

  for( int i=0; i<int(pfjets.size()); i++ ){
    double jetPt = pfjets.at(i).pt;
    double jetEta = pfjets.at(i).eta;
    double jetAbsEta = fabs(jetEta);

    bool eta = ( jetAbsEta<2.4 );
    bool jetId  = ( pfjets.at(i).jetIDLoose==1 );

    int flavour = pfjets.at(i).flavour;

    double factor=1;
    if( sysType.compare("data")!=0 ){
      double genJetPT = pfjets.at(i).genJetPT;
      if( sysType.compare("JERUp")==0 )        factor = getJERfactor(1,jetAbsEta,genJetPT,jetPt);
      else if( sysType.compare("JERDown")==0 ) factor = getJERfactor(-1,jetAbsEta,genJetPT,jetPt);
      else                                     factor = getJERfactor(0,jetAbsEta,genJetPT,jetPt);
    }
    jetPt *= factor;

    double unc = pfjets.at(i).JESunc;
    if( sysType.compare("JESUp")==0 )        jetPt *= (1. + unc);
    else if( sysType.compare("JESDown")==0 ) jetPt *= (1. - unc);

    double csv_old = pfjets.at(i).btagCombinedSecVertex;
    double csv = BEANs::reshape_csv( jetEta, jetPt, csv_old, flavour, sysType );

    bool csvM = ( csv>csvCut );
    if( jetPt>30. && eta && jetId ){
      tightJets.push_back(i);
      if( csvM ) tagJets.push_back(i);
      else       untagJets.push_back(i);
    }
  } // end loop over jets
}



double BEANs::getJERfactor( int returnType, double jetAbsETA, double genjetPT, double recojetPT){

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
  if( isFastSim_ ) diff_recojet_genjet += diff_FullSim_FastSim;

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



void BEANs::getSp(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, float &aplanarity, float &sphericity) {
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


void BEANs::getFox(vecTLorentzVector jets, float &h0, float &h1, float &h2, float &h3, float &h4) {
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


void BEANs::getFox_mod2(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, double HT,
			float &h0_mod2, float &h1_mod2, float &h2_mod2, float &h3_mod2, float &h4_mod2,  float &h5_mod2,  
			float &h6_mod2, float &h7_mod2, float &h8_mod2, float &h9_mod2, float &h10_mod2 ){

  int nJets = int(jets.size());

  int visObjects = nJets + 2; // change if # leps in event selection changes
    
  float eVis = HT;    
    
  TLorentzVector all_obj_vect[visObjects];
        
  for (int k=0; k<nJets; k++)
    {
      all_obj_vect[k] = jets[k];
    }
    
  all_obj_vect[nJets] = lepton;
  all_obj_vect[nJets + 1] = met;    

  h0_mod2 = 0.0;
  h1_mod2 = 0.0;
  h2_mod2 = 0.0;
  h3_mod2 = 0.0;
  h4_mod2 = 0.0;
  h5_mod2 = 0.0;
  h6_mod2 = 0.0;
  h7_mod2 = 0.0;
  h8_mod2 = 0.0;
  h9_mod2 = 0.0;
  h10_mod2 = 0.0;
    
    
    
  // according to original paper, the object pairs are double-counted, and also objects are "autocorrelated."
  // see Fox + Wolfram original paper in Nuc. Phys. B149 (1979) 413  ---- section 8 is relevant section
    
  for (int i=0; i<visObjects; i++)
    {
      for (int j=0; j<visObjects; j++)
        {
	  double angle1 = atan(all_obj_vect[i].Py() / all_obj_vect[i].Px());
            
	  if (all_obj_vect[i].Px() == -fabs(all_obj_vect[i].Px()))
            {
	      angle1 += 2*asin(1.0);
            }
            
	  double angle2 = atan(all_obj_vect[j].Py() / all_obj_vect[j].Px());
            
	  if (all_obj_vect[j].Px() == -fabs(all_obj_vect[j].Px()))
            {
	      angle2 += 2*asin(1.0);
            }
            
                        
	  double c0 = 1.;
	  double c1 = cos(angle1 - angle2);
	  double c2 = cos(2*(angle1 - angle2));
	  double c3 = cos(3*(angle1 - angle2));
	  double c4 = cos(4*(angle1 - angle2));
	  double c5 = cos(5*(angle1 - angle2));
	  double c6 = cos(6*(angle1 - angle2));
	  double c7 = cos(7*(angle1 - angle2));
	  double c8 = cos(8*(angle1 - angle2));
	  double c9 = cos(9*(angle1 - angle2));
	  double c10 = cos(10*(angle1 - angle2));
                        
	  float pipj = all_obj_vect[i].Perp()*all_obj_vect[j].Perp();
            
	  h0_mod2 += (pipj/(eVis*eVis))*c0;
	  h1_mod2 += (pipj/(eVis*eVis))*c1;
	  h2_mod2 += (pipj/(eVis*eVis))*c2;
	  h3_mod2 += (pipj/(eVis*eVis))*c3;
	  h4_mod2 += (pipj/(eVis*eVis))*c4;
	  h5_mod2 += (pipj/(eVis*eVis))*c5;
	  h6_mod2 += (pipj/(eVis*eVis))*c6;
	  h7_mod2 += (pipj/(eVis*eVis))*c7;
	  h8_mod2 += (pipj/(eVis*eVis))*c8;
	  h9_mod2 += (pipj/(eVis*eVis))*c9;
	  h10_mod2 += (pipj/(eVis*eVis))*c10;
                                
        }
    }

  return;
}



vdouble BEANs::getEffSF( int returnType, double jetPt, double jetEta, double jetId, std::string era){

  bool is2012 = ( era.find("2012")!=std::string::npos );


  bool debug = false;

  if (debug)
    cout << "getEffSF: Called getEffSF with the following parameters" << endl
         << "  returnType = " << returnType << ", jetPt = " << jetPt
         << ", jetEta = " << jetEta << ", jetId = " << jetId << ", era = " << era
         << endl;

  // check the histogram status before continuing


  bool bHistoOK =  (h_b_eff_ != 0);
  bool cHistoOK =  (h_c_eff_ != 0);
  bool lHistoOK =  (h_l_eff_ != 0);
  bool oHistoOK =  (h_o_eff_ != 0);

  if (debug)
    cout << "getEffSF: bHistoOK = " << bHistoOK << ", cHistoOK = " << cHistoOK << ", lHistoOK = "
         << lHistoOK << ", oHistoOK = " << oHistoOK << endl;
  
  if (!bHistoOK || !cHistoOK || !lHistoOK || !oHistoOK){
    cout << "getEffSF: Oops! We are missing one of the required histograms. Refusing to continue... "
         << endl;
    assert ( bHistoOK && cHistoOK && lHistoOK && oHistoOK);
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


  double SFb_central = 0.6981*((1.+(0.414063*pt))/(1.+(0.300155*pt)));

  double SFc_central = SFb_central;

  double SFb = SFb_central + m_type * SFb_error[use_bin];
  double SFc = SFc_central + m_type * 2* SFb_error[use_bin];

  // 2011/12 SFl correction function, JGWood 12Aug2012
  if(is2012){
    SFb = SFb_central + m_type * 1.5 * SFb_error[use_bin];
    SFc = SFc_central + m_type * 1.5 * 2* SFb_error[use_bin];
  }

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

  // 2011/12 SFl correction function, JGWood 12Aug2012
  //   CSVM    1.10422 + -0.000523856*x + 1.14251e-06*x*x   x=jetPt
  if(is2012) SFl *= 1.10422 + -0.000523856*pt + 1.14251e-06*pt*pt;

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


  vdouble result;
  result.clear();
  result.push_back(tagEff);
  result.push_back(SF);

  return result;
}


///////////////////////////////////////////////////////////
//
//   Input: your pf jets before selection and mc particles
//          the kind of tt+X you want (tt+lf, tt+cc, tt+bb)
//          the year (2011 or 2012)
//   Output: a bool saying whether or not to keep the event
//
///////////////////////////////////////////////////////////////


bool BEANs::ttPlusHeavyKeepEvent( BNmcparticleCollection const &mcparticles,
                           BNjetCollection const &pfjets,
                           TString ttbarType,
                           TString era ){

  // validate input
  bool validInput = false;
  if (ttbarType == "ttbar") validInput = true;
  if (ttbarType == "ttbar_bb") validInput = true;
  if (ttbarType == "ttbar_cc") validInput =true;

  if (!validInput ){
    cout << "ttPlusHeavyKeepEvent: could not recognize ttbarType " << ttbarType <<"... failing" <<endl;  
    assert (ttbarType == "ttbar or ttbar_bb or ttbar_cc");
  }


  bool validEra = false;
  if (era == "2011" ) validEra = true;
  if (era == "2012" ) validEra = true;

  if (!validEra ){
    cout << "ttPlusHeavyKeepEvent: could not recognize era " << era <<"... failing" <<endl;  
    assert (era == "2011 or 2012");
  }





  bool keepEvent = false;           
  bool debug_ = false;

  bool isWtoCS = false;

  if (debug_) cout << "Num MC particles = " << mcparticles.size() << std::endl
                   << "Num pfjets " <<  int(pfjets.size()) << std::endl;

   
  // check to see if the event has a c with a 
  // parent W
  for( unsigned i=0; i< mcparticles.size(); i++ ){
    int id = mcparticles.at(i).id;
    int motherID = mcparticles.at(i).motherId;
    int grandMotherID = mcparticles.at(i).grandMotherId;

    if (debug_) std::cout << "Particle " << i << " is " << id << ", has mother " << motherID << " and grandmother " << grandMotherID << std::endl;

    if (debug_) cout <<" Particle " << i << " has id " << id << endl;

    if (era == "2011") {
      if( abs(id)==4  && abs(motherID)==24 ){
	isWtoCS = true;
	break;
      }
    }
    else if (era == "2012") {
      int daughter0ID = mcparticles.at(i).daughter0Id;
      int daughter1ID = mcparticles.at(i).daughter1Id;

      if( abs(id)==24 && ((abs(daughter0ID)==3 && abs(daughter1ID)==4) || (abs(daughter0ID)==4 && abs(daughter1ID)==3)) ){
	isWtoCS = true;
	break;
      }
    }
  }


  bool isBBbarEvent = false;
  bool isCCbarEvent = false;


  bool gotB = false;
  bool gotC = false;

  int numBmomB=0;
  int numBmomT=0;
  int numBmomHiggs=0;
  int numBbarmomBbar=0;
  int numBbarmomTbar=0;
  int numBbarmomHiggs=0;
  int numCmomC=0;
  int numCbarmomCbar=0;
  int numCmomHiggs=0;
  int numCbarmomHiggs=0;

  if (debug_) cout << "Starting loop over pf jet parton ids to see if you have a b" <<endl;

  for( int i=0; i<int(pfjets.size()); i++ ){


          
    int id = pfjets.at(i).genPartonId;
    if( id==-99 ) continue;
    int motherID = pfjets.at(i).genPartonMotherId;
    int mother0ID =  pfjets.at(i).genPartonMother0Id;
    int mother1ID =  pfjets.at(i).genPartonMother1Id;
          


    // check to see if pf jets is from a  b/c and mother is a gluon
    // or, if mother is some light quark

    if (era == "2011") {
      if( abs(id)==5 && ( motherID==21 || abs(motherID)<abs(id) ) ) gotB=true;
      if( abs(id)==4 && ( motherID==21 || abs(motherID)<abs(id) ) ) gotC=true;
    } else if (era == "2012") {
      if( abs(id)==5 && abs(mother0ID) != 6  && abs(mother1ID) != 6 ) gotB=true;
      // basically, as long as you didn't come from a W, then you are a tt+cc
      if( abs(id)==4 && (abs(mother0ID)==21 || abs(mother1ID) == 21 || abs(mother0ID) < abs(id) || abs(mother1ID) < abs(id))
	  && abs(mother0ID)!=24 && abs(mother1ID)!=24 ) gotC=true;
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

    if( id==5 && motherID==25 ) numBmomHiggs++;
    if( id==-5 && motherID==25 ) numBbarmomHiggs++;

    if( id==4 && motherID==25 ) numCmomHiggs++;
    if( id==-4 && motherID==25 ) numCbarmomHiggs++;
  }

  //std::cout << "b->b: " << numBmomB << " bbar->bbar: " << numBbarmomBbar << " t->b: " << numBmomT << " tbar->bbar: " << numBbarmomTbar << " H->b: " << numBmomHiggs << " H->bbar: " << numBbarmomHiggs << std::endl;
  //std::cout << "c->c: " << numCmomC << " cbar->cbar: " << numCbarmomCbar << " H->c: " << numCmomHiggs << " H->cbar: " << numCbarmomHiggs << std::endl;
        

  // if at least one b from b & one b from t, or if CC, and your jet was not b
  if( ((numBmomB>=1 && numBmomT>=1) || (numBbarmomBbar>=1 && numBbarmomTbar>=1)) && !gotB ){
    if (debug_) std::cout << "No sign of a b jet, but looping over jets again to check"
                          <<std::endl;
    // for each jet that is  b from b
    for( int i=0; i<int(pfjets.size()); i++ ){
      if (debug_) cout << "LOOP: i = " << i << endl;
      int id0 = pfjets.at(i).genPartonId;
      int motherID0 = pfjets.at(i).genPartonMotherId;
      if( !(abs(id0)==5 && motherID0==id0) ) continue;

      if (debug_) std::cout << "Jet index " << i  << " is a bjet, let us see that it is not from top" <<std::endl;
      // for each jet that is b from t
      for( int j=0; j<int(pfjets.size()); j++ ){
        if (debug_) cout << "LOOP: j = " << j << endl;
        int id1 = pfjets.at(j).genPartonId;
        int motherID1 = pfjets.at(j).genPartonMotherId;
        if (debug_) std::cout << "LOOP: id0 = " << id0 << ", motherID0 = " << motherID0
                              << ", id1 = " << id1 << ", motherID1 = " << motherID1 << endl
                              << "continue? = " << !(id1==id0 && abs(motherID1)==6) << endl;
                
        if( !(id1==id0 && abs(motherID1)==6) ) continue;
        if (debug_) std::cout << "You didn't skip this event!" << endl;
        // if delta r between b from b and b from t is big enough, then b in final state is OK
        double dR = reco::deltaR(pfjets.at(i).genPartonEta,
                                   pfjets.at(i).genPartonPhi,
                                   pfjets.at(j).genPartonEta,
                                   pfjets.at(j).genPartonPhi);
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


  
  if( (ttbarType == "ttbar") && !isBBbarEvent && !isCCbarEvent ) keepEvent = true;
  else if( (ttbarType == "ttbar_bb") &&  isBBbarEvent && !isCCbarEvent ) keepEvent = true;
  else if( (ttbarType == "ttbar_cc")  && !isBBbarEvent && isCCbarEvent  ) keepEvent = true;
      

  if (debug_) cout << "Filter result = " << keepEvent << endl
                   << "isBBbarEvent = " << isBBbarEvent << endl
                   << "isCCbarEvent = " << isCCbarEvent << endl
                   << "... will we skip this? " << (!keepEvent) << endl;
        
  return keepEvent;

}

double BEANs::reshape_csv( double eta, double pt, double csv, int flavor, std::string sysType ){

  double new_csv = csv;

  if( sysType.compare("data")==0 ) return csv;

  if( sysType.compare("hfSFUp")==0 )        new_csv = sh_hfSFUp_->reshape(eta, pt, csv, flavor);
  else if( sysType.compare("hfSFDown")==0 ) new_csv = sh_hfSFDown_->reshape(eta, pt, csv, flavor);
  else if( sysType.compare("lfSFUp")==0 )   new_csv = sh_lfSFUp_->reshape(eta, pt, csv, flavor);
  else if( sysType.compare("lfSFDown")==0 ) new_csv = sh_lfSFDown_->reshape(eta, pt, csv, flavor);
  else                                      new_csv = sh_->reshape(eta, pt, csv, flavor);

  return new_csv;
}



#endif // INC_BEANSUTILITIES
