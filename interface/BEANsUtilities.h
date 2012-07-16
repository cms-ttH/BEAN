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
std::string eff_base_dir(my_pPath);
std::string str_eff_file = eff_base_dir + "/src/NtupleMaker/BEANmaker/interface/mc_btag_efficiency.root";


TFile *f_tag_eff_ = new TFile(str_eff_file.c_str());
TH2D* h_b_eff_ = (TH2D*)f_tag_eff_->Get("ttH120_h_jet_pt_eta_b_eff");
TH2D* h_c_eff_ = (TH2D*)f_tag_eff_->Get("ttH120_h_jet_pt_eta_c_eff");
TH2D* h_l_eff_ = (TH2D*)f_tag_eff_->Get("ttH120_h_jet_pt_eta_l_eff");
TH2D* h_o_eff_ = (TH2D*)f_tag_eff_->Get("ttH120_h_jet_pt_eta_o_eff");


using namespace std;

namespace BEANs{
  const double PI=2.0*acos(0.);
  const double TWOPI=2.0*PI;
  const float ETA_LIMIT=15.0;
  const float EPSILON=1.E-10;

  void setMCsample( int insample=2500 );

  void electronSelector( BNelectronCollection electrons, std::string era, vint &tightElectrons, vint &looseElectrons );
  void muonSelector( BNmuonCollection muons, std::string era, vint &tightMuons, vint &looseMuons );
  void jetSelector( BNjetCollection pfjets, std::string sysType, vint &tightJets, vint &tagJets, vint &untagJets, 
		    std::vector<BTagWeight::JetInfo> &myjetinfo, double csvCut = 0.679 );

  vdouble getEffSF( int returnType, double jetPts, double jetEtas, double jetIds );
  double getJERfactor( int returnType, double jetAbsETA, double genjetPT, double recojetPT );
  void getSp(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, float &aplanarity, float &sphericity);
  void getFox(vecTLorentzVector jets,float &h0, float &h1, float &h2, float &h3, float &h4);	
  void getFox_mod2(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, double HT,
		   float &h0_mod2, float &h1_mod2, float &h2_mod2, float &h3_mod2, float &h4_mod2,  float &h5_mod2,  
		   float &h6_mod2, float &h7_mod2, float &h8_mod2, float &h9_mod2, float &h10_mod2 );
}



void setMCsample( int insample=2500 ){

  std::string samplename = "ttH120";
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

  h_b_eff_ = (TH2D*)f_tag_eff_->Get(std::string( samplename + "_h_jet_pt_eta_b_eff" ).c_str());
  h_c_eff_ = (TH2D*)f_tag_eff_->Get(std::string( samplename + "_h_jet_pt_eta_c_eff" ).c_str());
  h_l_eff_ = (TH2D*)f_tag_eff_->Get(std::string( samplename + "_h_jet_pt_eta_l_eff" ).c_str());
  h_o_eff_ = (TH2D*)f_tag_eff_->Get(std::string( samplename + "_h_jet_pt_eta_o_eff" ).c_str());

}


/////////
///
/// Electrons
///
////////
void BEANs::electronSelector( BNelectronCollection electrons, std::string era, vint &tightElectrons, vint &looseElectrons ){

  tightElectrons.clear();
  looseElectrons.clear();

  if( era.compare("2011")==0 ){
    for( int i=0; i<int(electrons.size()); i++ ){
      double eleSCEta = electrons.at(i).scEta;
      double absSCeta = fabs(eleSCEta);
      double eleEta = electrons.at(i).eta;
      double eleEt = electrons.at(i).et;

      bool isCrack = ( (absSCeta>1.4442) && (absSCeta<1.5660) );

      bool kin = ( (eleEt>15.) && !isCrack && fabs(eleEta)<2.5 );

      if( !kin ) continue;

      double chargedHadronIso = electrons.at(i).chargedHadronIso;
      double neutralHadronIso = electrons.at(i).neutralHadronIso;
      double photonIso = electrons.at(i).photonIso;

      double relIso = ( chargedHadronIso + neutralHadronIso + photonIso ) * 1./eleEt;

      bool looseIso = ( relIso < 0.2 );
      bool tightIso = ( relIso < 0.1 );

      int eidHyperTight1MC = electrons.at(i).eidHyperTight1MC;
      bool eidHyperTight1MC_dec = ( (eidHyperTight1MC & 1)==1 );

      bool d0 = ( fabs(electrons.at(i).correctedD0) < 0.02 );
      bool dZ = ( fabs(electrons.at(i).correctedDZ) < 1. );

      bool dist  = ( fabs(electrons.at(i).dist)<0.02 );
      bool dcot  = ( fabs(electrons.at(i).dcot)<0.02 );
      bool nlost = ( electrons.at(i).numberOfLostHits<1 );
      bool notConv = ( !(dist && dcot) && nlost );

      bool id = ( eidHyperTight1MC_dec && d0 && dZ && notConv );

      if( kin && looseIso ){
	if( ((eleEt>30.) && id && tightIso) ) tightElectrons.push_back(i);
	else looseElectrons.push_back(i);
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

      bool kin = ( (elePt>20.) && !isCrack && fabs(eleEta)<2.5 );

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
	if( ((elePt>30.) && id && tightIso) ) tightElectrons.push_back(i);
	else looseElectrons.push_back(i);
      }
    }// end electron loop
  }
} //end electronSelector


/////////
///
/// Muons
///
////////
void BEANs::muonSelector( BNmuonCollection muons, std::string era, vint &tightMuons, vint &looseMuons ){

  tightMuons.clear();
  looseMuons.clear();

  if( era.compare("2011")==0 ){
    for( int i=0; i<int(muons.size()); i++ ){
      double muPt  = muons.at(i).pt;
      double muEta = muons.at(i).eta;
      double muAbsEta = fabs(muEta);

      bool kin = ( (muPt>10.) && (muAbsEta<2.4) );

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
	if( ((muPt>30.) && (muAbsEta<2.1) && id && tightIso) ) tightMuons.push_back(i);
	else looseMuons.push_back(i);
      }
    }// end muon loop
  } // end if 2011
  else{ // default is 2012 selection
    for( int i=0; i<int(muons.size()); i++ ){
      double muPt  = muons.at(i).pt;
      double muEta = muons.at(i).eta;
      double muAbsEta = fabs(muEta);

      bool kin = ( (muPt>10.) && (muAbsEta<2.5) );

      if( !kin ) continue;

      double pfIsoR04SumChargedHadronPt = muons.at(i).pfIsoR04SumChargedHadronPt;
      double pfIsoR04SumNeutralHadronEt = muons.at(i).pfIsoR04SumNeutralHadronEt;
      double pfIsoR04SumPhotonEt = muons.at(i).pfIsoR04SumPhotonEt;
      double pfIsoR04SumPUPt = muons.at(i).pfIsoR04SumPUPt;

      double relIso_dBeta = (pfIsoR04SumChargedHadronPt + max(0.0, pfIsoR04SumNeutralHadronEt + pfIsoR04SumPhotonEt - 0.5*pfIsoR04SumPUPt))/muPt;

      bool looseIso = ( relIso_dBeta<0.20 );
      bool tightIso = ( relIso_dBeta<0.12 );

      bool isPFmuon = ( muons.at(i).isPFMuon==1 );
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
	if( ((muPt>30.) && (muAbsEta<2.1) && id && tightIso) ) tightMuons.push_back(i);
	else looseMuons.push_back(i);
      }
    }// end muon loop
  }
} //end muonSelector


/////////
///
/// PFJets
///
////////
void BEANs::jetSelector( BNjetCollection pfjets, std::string sysType, vint &tightJets, vint &tagJets, vint &untagJets, 
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
	if( sysType.compare("hfSFUp")==0 )        myEffSF = BEANs::getEffSF( 1,  jetPt, jetEta, flavour );
	else if( sysType.compare("hfSFDown")==0 ) myEffSF = BEANs::getEffSF( -1, jetPt, jetEta, flavour );
	else if( sysType.compare("lfSFUp")==0 )   myEffSF = BEANs::getEffSF( 2,  jetPt, jetEta, flavour );
	else if( sysType.compare("lfSFDown")==0 ) myEffSF = BEANs::getEffSF( -2, jetPt, jetEta, flavour );
	else                                      myEffSF = BEANs::getEffSF( 0,  jetPt, jetEta, flavour );

	BTagWeight::JetInfo myjet( myEffSF[0], myEffSF[1] );
	myjetinfo.push_back(myjet);
      }
    }
  } // end loop over jets
}



double BEANs::getJERfactor( int returnType, double jetAbsETA, double genjetPT, double recojetPT){

  double factor = 1.;

  double scale_JER = 1., scale_JERup = 1., scale_JERdown = 1.;
  if( jetAbsETA<0.5 ){ 
    scale_JER = 1.052; scale_JERup = 1.052 + sqrt( 0.012*0.012 + 0.062*0.062 ); scale_JERdown = 1.052 - sqrt( 0.012*0.012 + 0.061*0.061 );
  }
  else if( jetAbsETA<1.1 ){ 
    scale_JER = 1.057; scale_JERup = 1.057 + sqrt( 0.012*0.012 + 0.056*0.056 ); scale_JERdown = 1.057 - sqrt( 0.012*0.012 + 0.055*0.055 );
  }
  else if( jetAbsETA<1.7 ){ 
    scale_JER = 1.096; scale_JERup = 1.096 + sqrt( 0.017*0.017 + 0.063*0.063 ); scale_JERdown = 1.096 - sqrt( 0.017*0.017 + 0.062*0.062 );
  }
  else if( jetAbsETA<2.3 ){ 
    scale_JER = 1.134; scale_JERup = 1.134 + sqrt( 0.035*0.035 + 0.087*0.087 ); scale_JERdown = 1.134 - sqrt( 0.035*0.035 + 0.085*0.085 );
  }
  else if( jetAbsETA<5.0 ){ 
    scale_JER = 1.288; scale_JERup = 1.288 + sqrt( 0.127*0.127 + 0.155*0.155 ); scale_JERdown = 1.288 - sqrt( 0.127*0.127 + 0.153*0.153 );
  }

  double jetPt_JER = recojetPT;
  double jetPt_JERup = recojetPT;
  double jetPt_JERdown = recojetPT;

  if( genjetPT>10. ){
    jetPt_JER = std::max( 0., genjetPT + scale_JER * ( recojetPT - genjetPT ) );
    jetPt_JERup = std::max( 0., genjetPT + scale_JERup * ( recojetPT - genjetPT ) );
    jetPt_JERdown = std::max( 0., genjetPT + scale_JERdown * ( recojetPT - genjetPT ) );
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


vdouble BEANs::getEffSF( int returnType, double jetPt, double jetEta, double jetId ){

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


  vdouble result;
  result.clear();
  result.push_back(tagEff);
  result.push_back(SF);

  return result;
}



#endif // INC_BEANSUTILITIES
