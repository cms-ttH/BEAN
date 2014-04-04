#ifndef _BEANhelper_h
#define _BEANhelper_h

#include <iostream>
#include <vector>
#include <map>
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
#include "TH1D.h"
#include "TH2D.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif

#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "CommonTools/Utils/interface/normalizedPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

//Headers for the data items
#include "BEAN/Collections/interface/BNleptonCollection.h"
#include "BEAN/Collections/interface/BNelectron.h"
#include "BEAN/Collections/interface/BNevent.h"
#include "BEAN/Collections/interface/BNgenjet.h"
#include "BEAN/Collections/interface/BNjet.h"
#include "BEAN/Collections/interface/BNmcparticle.h"
#include "BEAN/Collections/interface/BNmet.h"
#include "BEAN/Collections/interface/BNmuon.h"
#include "BEAN/Collections/interface/BNphoton.h"
#include "BEAN/Collections/interface/BNprimaryvertex.h"
#include "BEAN/Collections/interface/BNskimbits.h"
#include "BEAN/Collections/interface/BNsupercluster.h"
#include "BEAN/Collections/interface/BNtau.h"
#include "BEAN/Collections/interface/BNtrack.h"
#include "BEAN/Collections/interface/BNtrigger.h"
#include "BEAN/Collections/interface/BNtrigobj.h"

#include "BEAN/BEANmaker/interface/BtagWeight.h"
#include "BEAN/BEANmaker/interface/CSVreevaluator.h"


#endif

typedef std::map<std::string, std::string> mparams;
typedef std::vector< TLorentzVector > vecTLorentzVector;
typedef std::vector<std::vector<double> > vvdouble;
typedef std::vector<std::vector<std::string> > vvstring;
typedef std::vector<std::vector<int> > vvint;
typedef std::vector<std::string> vstring;
typedef std::vector<double> vdouble;
typedef std::vector<int> vint;

typedef BNeventCollection::const_iterator		  EventIter;
typedef BNelectronCollection::const_iterator	  EleIter;
typedef BNjetCollection::const_iterator			  JetIter;
typedef BNmcparticleCollection::const_iterator	  MCparIter;
typedef BNmetCollection::const_iterator			  MetIter;
typedef BNmuonCollection::const_iterator		  MuonIter;
typedef BNtauCollection::const_iterator			  TauIter;
typedef BNphotonCollection::const_iterator		  PhotonIter;
typedef BNprimaryvertexCollection::const_iterator PVIter;
typedef BNskimbitsCollection::const_iterator	  SkimBitIter;
typedef BNsuperclusterCollection::const_iterator  SCIter;
typedef BNtrackCollection::const_iterator		  TrackIter;
typedef BNtriggerCollection::const_iterator		  TrigIter;
typedef BNtrigobjCollection::const_iterator		  TrigObjIter;

namespace analysisType{ enum analysisType{ LJ, DIL, TauLJ, TauDIL }; }
namespace sysType{enum sysType{NA, JERup, JERdown, JESup, JESdown, hfSFup, hfSFdown, lfSFdown, lfSFup, TESup, TESdown, CSVLFup, CSVLFdown, CSVHFup, CSVHFdown, CSVHFStats1up, CSVHFStats1down, CSVLFStats1up, CSVLFStats1down, CSVHFStats2up, CSVHFStats2down, CSVLFStats2up, CSVLFStats2down, CSVCErr1up, CSVCErr1down, CSVCErr2up, CSVCErr2down }; }
namespace jetID{		enum jetID{			none, jetMinimal, jetLooseAOD, jetLoose, jetTight }; }
namespace tauID{		enum tauID{			tauNonIso, tauVLoose, tauLoose, tauMedium, tauTight }; }
namespace muonID{		enum muonID{		muonSide, muonSideLooseMVA, muonSideTightMVA, muonLoose, muonTight, muonPtOnly, muonPtEtaOnly, muonPtEtaIsoOnly, muonPtEtaIsoTrackerOnly, muonNoCuts }; }
namespace electronID{	enum electronID{	electronSide, electronSideLooseMVA, electronSideTightMVA, electronLoose, electronTight, electronTightMinusTrigPresel, electronLooseMinusTrigPresel, electronNoCuts }; }
namespace hdecayType{	enum hdecayType{ hbb, hcc, hww, hzz, htt, hgg, hjj, hzg }; }



using namespace std;

class BEANhelper{

	// === Functions === //
	public: 
		// Constructor(s) and destructor
		BEANhelper();
		virtual ~BEANhelper();
		
		// Set up BEANhelper
		void SetUp(string, int, const analysisType::analysisType, bool, string, bool, bool, string iCollisionDS="All");

		template <typename BNobject> void PrintInfo(const BNobject&);
		template <typename BNobject> BNmcparticle GetMatchedMCparticle(const BNmcparticleCollection&, const BNobject&, const double) const;
		template <typename BNobject> BNgenjet GetMatchedGenjet(const BNgenjetCollection&, const BNobject&, const double);
		template <typename BNobject> BNmcparticle GetMatchedGentau(const BNmcparticleCollection&, const BNobject&, const double);
		template <typename BNobject> BNmcparticle GetMatchedVisGentau(const BNmcparticleCollection&, const BNobject&, const double);
		template <typename BNobject> BNtau GetMatchedTau(const BNtauCollection&, const BNobject&, const double);
		template <typename BNobject> BNjet GetClosestJet(const BNjetCollection&, const BNobject&, const double);
		template <typename BNobject> BNtrack GetClosestTrack(const BNtrackCollection&, const BNobject&, const double);

		// Union, intersection, difference
		template <typename BNcollection> BNcollection GetSortedByPt(const BNcollection&);
		template <typename BNcollection> BNcollection GetSortedByCSV(const BNcollection&);
		template <typename BNcollection> BNcollection GetUnion(const BNcollection&, const BNcollection&);
		template <typename BNcollection> BNcollection GetIntersection(const BNcollection&, const BNcollection&);
		template <typename BNcollection> BNcollection GetDifference(const BNcollection&, const BNcollection&);
		template <typename BNcollection1, typename BNcollection2> BNcollection1 GetDifference(const BNcollection1&, const BNcollection2&, const double);
		template <typename BNcollection> BNcollection GetSymmetricDifference(const BNcollection&, const BNcollection&);
		template <typename BNcollection> BNcollection GetUnionUnsorted(const BNcollection&, const BNcollection&);


		// Jets and MET
		float GetBtagWeight(const BNjet&, const sysType::sysType iSysType=sysType::NA);
		bool PassesCSV(const BNjet&, const char);
		bool IsGoodJet(const BNjet&, const float, const float, const jetID::jetID, const char);
		BNmet GetCorrectedMET(const BNmet&, const BNjetCollection&, const sysType::sysType iSysType=sysType::NA);
		BNjet GetCorrectedJet(const BNjet&, const sysType::sysType iSysType=sysType::NA);
		BNjetCollection GetSelectedJets(const BNjetCollection&, const float, const float, const jetID::jetID, const char);
		BNjetCollection GetUncorrectedJets(const BNjetCollection&, const BNjetCollection&);
		BNjetCollection GetCorrectedJets(const BNjetCollection&, const sysType::sysType iSysType=sysType::NA);
		BNjetCollection GetCleanJets(const BNjetCollection&, const vector<TLorentzVector>&, const float, std::vector<unsigned int>* jet_indices=0);
		BNjetCollection GetCleanJets(const BNjetCollection&, const BNleptonCollection&, const float);
		BNjetCollection GetCleanJets_cProj(const BNjetCollection&, const BNleptonCollection&, const float);
		unsigned int GetNumCSVbtags(const BNjetCollection&, const char, std::vector<unsigned int>* jet_indices=0);
		unsigned int GetNumNonCSVbtags(const BNjetCollection&, const char, std::vector<unsigned int>* jet_indices=0);
		float GetHT(const BNjetCollection&);

		// Taus
		bool IsVLooseTau(const BNtau&);
		bool IsLooseTau(const BNtau&);
		bool IsMediumTau(const BNtau&);
		bool IsTightTau(const BNtau&);
		bool IsGoodTau(const BNtau&, const tauID::tauID);
		float GetTauSF(const BNtau&);
		BNtauCollection GetSelectedTaus(const BNtauCollection&, const tauID::tauID);
		BNtauCollection GetCorrectedTaus(const BNtauCollection&, const sysType::sysType iSysType=sysType::NA);
		BNtau GetCorrectedTau(const BNtau&, const sysType::sysType iSysType=sysType::NA);
		bool IsTauEvent(const BNtauCollection&, const BNjetCollection&, const BNelectronCollection&, const BNmuonCollection&, const sysType::sysType iSysType=sysType::NA);
		bool IsTauLeptonLeptonEvent(const BNtauCollection&, const BNjetCollection&, const BNelectronCollection&, const BNmuonCollection&, const sysType::sysType iSysType=sysType::NA);
		bool IsTauTauLeptonEvent(const BNtauCollection&, const BNjetCollection&, const BNelectronCollection&, const BNmuonCollection&, const sysType::sysType iSysType=sysType::NA);
  
		// Muons
		bool IsSideMuon(const BNmuon&);
		bool IsSideMuonLooseMVA(const BNmuon&, const BNjetCollection* = 0);
		bool IsSideMuonTightMVA(const BNmuon&, const BNjetCollection* = 0);
		bool IsLooseMuon(const BNmuon&);
		bool IsTightMuon(const BNmuon&);
		bool IsGoodMuon(const BNmuon&, const muonID::muonID, const BNjetCollection* = 0);
		float GetMuonRelIso(const BNmuon&) const;
		float GetMuonSF(const BNmuon&, const muonID::muonID m = muonID::muonTight) const;
		float GetDoubleMuonTriggerSF ( const BNmuon&, const BNmuon& );
		float GetMuonEleTriggerSF ( const BNmuon&, const BNelectron & );
		float GetMuonLepMVA( const BNmuon&, const BNjetCollection* = 0);

  // Test only
  float TestSingleMuonTriggerNew ( const BNmuon& );
  float TestSingleEleTriggerNew ( const BNelectron & );
  float TestSingleMuonTriggerOld ( const BNmuon& );
  float TestSingleEleTriggerOld ( const BNelectron & );
  
  
		BNmuonCollection GetSelectedMuons(const BNmuonCollection&, const muonID::muonID, const BNjetCollection* = 0);

		// Electrons
		bool IsSideElectron(const BNelectron&);
		bool IsSideElectronLooseMVA(const BNelectron&);
		bool IsSideElectronTightMVA(const BNelectron&);
		bool IsLooseElectron(const BNelectron&);
		bool IsTightElectron(const BNelectron&);
		float GetElectronRelIso(const BNelectron&) const;
		float GetElectronSF(const BNelectron&, const electronID::electronID e = electronID::electronTight ) const;
		float GetDoubleElectronTriggerSF (const BNelectron&, const BNelectron&);
		bool GetElectronIDresult(const BNelectron& iElectron, const electronID::electronID) const;
		bool IsGoodElectron(const BNelectron&, const electronID::electronID);
		BNelectronCollection GetSelectedElectrons(const BNelectronCollection&, const electronID::electronID);
		float GetElectronLepMVA(const BNelectron&);

		// General lepton functions
		float GetDBCorrectedRelIsoDR04(const BNlepton& iLepton, const float& dBeta_factor) const;
  
		// MCparticles
		BNmcparticleCollection	GetSelectedMCparticlesByGrandParentPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetSelectedMCparticlesByParentPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetSelectedMCparticlesByPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetSelectedMCparticlesByChildPDGid(const BNmcparticleCollection&, const vector<int>);

		BNmcparticleCollection	GetUnrejectedMCparticlesByGrandParentPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetUnrejectedMCparticlesByParentPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetUnrejectedMCparticlesByPDGid(const BNmcparticleCollection&, const vector<int>) const;
		BNmcparticleCollection	GetUnrejectedMCparticlesByChildPDGid(const BNmcparticleCollection&, const vector<int>);

		BNmcparticleCollection	GetSelectedMCparticlesByGrandParentStatus(const BNmcparticleCollection&, const bool, const bool, const bool);
		BNmcparticleCollection	GetSelectedMCparticlesByParentStatus(const BNmcparticleCollection&, const bool, const bool, const bool);
		BNmcparticleCollection	GetSelectedMCparticlesByStatus(const BNmcparticleCollection&, const bool, const bool, const bool) const;
		BNmcparticleCollection	GetSelectedMCparticlesByChildStatus(const BNmcparticleCollection&, const bool, const bool, const bool);

		BNmcparticleCollection	GetParents(const BNmcparticle&, const BNmcparticleCollection&);
		void					DrawFeynman(const BNmcparticle&);
		BNmcparticleCollection	GetGenTaus(const BNmcparticleCollection&);
		BNmcparticleCollection	GetHadronicGenTaus(const BNmcparticleCollection&);
		BNmcparticle			GetVisGenTau(const BNmcparticle&, const BNmcparticleCollection&);
		bool					ttPlusHeavyKeepEvent(const BNmcparticleCollection&, const BNjetCollection&);
		bool					ttPlusHFKeepEvent(const BNmcparticleCollection&, const BNjetCollection&);
		int						ttPlusBBClassifyEvent(const BNmcparticleCollection&, const BNjetCollection&);
		int						ttPlusCCClassifyEvent(const BNmcparticleCollection&, const BNjetCollection&);
		unsigned int			GetNumExtraPartons(const BNmcparticleCollection&);
		hdecayType::hdecayType	GetHdecayType(const BNmcparticleCollection&);
		bool					keepHdecayType(const BNmcparticleCollection&, const hdecayType::hdecayType);

		// PU reweighing
		double GetPUweight(const double);
		double GetPUweightUp(const double);
		double GetPUweightDown(const double);

		// CSV reweighting
		double GetCSVweight(const BNjetCollection&, const sysType::sysType iSysType=sysType::NA);
		vdouble GetCSVweights(const BNjetCollection&, const sysType::sysType iSysType=sysType::NA);

		// Top PT reweighting
		double TopPtWeight(double);
		double GetTopPtweight(const BNmcparticleCollection&);
		double GetTopPtweightUp(const BNmcparticleCollection&);
		double GetTopPtweightDown(const BNmcparticleCollection&);

		// Q^2 systematic
		double GetQ2ScaleUp(const BNevent&);
		double GetQ2ScaleDown(const BNevent&);

  // Trigger utility functions 
  bool IsAnyTriggerBitFired ( const vector<string> & targetTriggers, const BNtriggerCollection & triggerBits);
  
  bool SingleObjectMatchesAnyTrigger (double recoEta, double recoPhi, const vector<string> & targetTriggers,
									  const BNtrigobjCollection & triggerObjects);
  
  bool DoubleObjectMatchesAnyTrigger (double firstEta, double firstPhi, double secondEta, double secondPhi,
									  const vector<string> & targetTriggers,  const BNtrigobjCollection & triggerObjects,
									  bool sumFilterResults = false);

  
  bool MuonMatchesSingleMuTrigger(const BNmuon& iMuon, const BNtriggerCollection & triggerBits,
											  const BNtrigobjCollection & triggerObjects );

  
  bool MuonsMatchDoubleMuTrigger(const BNmuon& iMuon, const BNmuon & jMuon,
											 const BNtriggerCollection & triggerBits,
											 const BNtrigobjCollection & triggerObjects);

  bool ElectronMatchesSingleEleTrigger(const BNelectron& iEle, const BNtriggerCollection & triggerBits,
									   const BNtrigobjCollection & triggerObjects );

  
  bool ElectronsMatchDoubleEleTrigger(const BNelectron& iEle, const BNelectron & jEle,
									  const BNtriggerCollection & triggerBits,
									  const BNtrigobjCollection & triggerObjects);

  bool MuEGMatchMuEGTrigger(const BNmuon & iEle, const BNelectron & jEle,
							const BNtriggerCollection & triggerBits,
							const BNtrigobjCollection & triggerObjects);
  
  // Trigger utility variable.
  const static bool triggerDebug = false;


  
  

	private:
		inline void ThrowFatalError(const std::string& m) const { cerr << "[ERROR]\t" << m << " Cannot continue. Terminating..." << endl; exit(1); };
		double GetCSVvalue(const BNjet&, const sysType::sysType iSysType=sysType::NA);
		inline void CheckSetUp() const { if(!isSetUp){ ThrowFatalError("BEANhelper not yet set up."); } };
		void SetUpPUreweighing(string const);
		void SetUpCSVreweighting();
		void SetUpCSVreshaping();
		void SetUpJetSF();
		void SetUpLeptonSF();
		string GetSampleName();

		// Parameter management
	private:
		bool isSetUp;
		string era;
		int sampleNumber;
		analysisType::analysisType analysis;
		bool isData;
		string dataset;
		bool reshapeCSV;
		bool usePfLeptons;

		TFile*			jetSFfile;
		TFile*			leptonSFfile;
  TFile*		  doubleMuonTriggerFile;
  TFile*		  doubleEleTriggerFile;
  TFile*		  muonEleTriggerFile;
  TFile*		  oldLeptonScaleFactFile;

		// PU reweighing
		TFile*			puFile;
		TH1D*			h_PU_ratio;
		TH1D*			h_PUup_ratio;
		TH1D*			h_PUdown_ratio;

		// CSV reweighting
		TFile* f_CSVwgt_HF;
		TFile* f_CSVwgt_LF;


		// CSV reshaping
		CSVreevaluator*	  sh_;
		CSVreevaluator*	  sh_hfSFUp_;
		CSVreevaluator*	  sh_hfSFDown_;
		CSVreevaluator*	  sh_lfSFUp_;
		CSVreevaluator*	  sh_lfSFDown_;

		// lepMVA TMVA readers
		TMVA::Reader*	  mu_reader_high_b;
		TMVA::Reader*	  mu_reader_high_e;
		TMVA::Reader*	  mu_reader_low_b;
		TMVA::Reader*	  mu_reader_low_e;
		TMVA::Reader*	  ele_reader_high_cb;
		TMVA::Reader*	  ele_reader_high_fb;
		TMVA::Reader*	  ele_reader_high_ec;
		TMVA::Reader*	  ele_reader_low_cb;
		TMVA::Reader*	  ele_reader_low_fb;
		TMVA::Reader*	  ele_reader_low_ec;

		Float_t varneuRelIso;
		Float_t varchRelIso;
		Float_t varjetDR_in;
		Float_t varjetPtRatio_in;
		Float_t varjetBTagCSV_in;
		Float_t varsip3d;
		Float_t varmvaId;
		Float_t varinnerHits;
		Float_t vardxy;
		Float_t vardz;
   
		// Old functions
	public:
		vdouble getEffSF( int returnType, double jetPts, double jetEtas, double jetIds );
		double getJERfactor( int returnType, double jetAbsETA, double genjetPT, double recojetPT );
		void getSp(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, float &aplanarity, float &sphericity);
		void getFox(vecTLorentzVector jets,float &h0, float &h1, float &h2, float &h3, float &h4);	

	protected:

	private:


	// === Variables === //
	public:

	protected:

	private:
		mparams			params;
		float			CSVLwp, CSVMwp, CSVTwp;

		// CSV reweighting
		TH1D* h_csv_wgt_hf[9][5];
		TH1D* c_csv_wgt_hf[9][5];
		TH1D* h_csv_wgt_lf[9][3][3];


		// Old parameters
		TH2D*			h_b_eff_;
		TH2D*			h_c_eff_;
		TH2D*			h_l_eff_;
		TH2D*			h_o_eff_;
		TH2D*			h_ele_SF_;
		TH2D*			h_SingleEle_trig_SF_;
		TH2D*			h_mu_SF_;
  TH2D*			  h_doubleMuTrigSF;
  TH2D*			  h_doubleEleTrigSF;
  TH2D*			  h_muonEleTrigSF;
  TH2D*			  h_looseMuonSF;
  TH2D*			  h_looseEleSF;
  // testing -- temporary only
  TH2D*			  h_testSingleMuNew;
  TH2D*			  h_testSingleMuOld;
  TH2D*			  h_testSingleEleNew;
  TH2D*			  h_testTrigSingleEleNew;
  TH2D*			  h_testSingleEleOld;
  
		string			samplename;

}; // End of class prototype


// === Return matched BNmcparticle based on minimum deltaR and a deltaR threshold === //
template <typename BNobject> BNmcparticle BEANhelper::GetMatchedMCparticle(const BNmcparticleCollection& iMCparticles, const BNobject& iObject, const double iMaxDeltaR) const
{
	BNmcparticle result;
	double minDeltaR = 999;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
		double thisDeltaR = deltaR(MCparticle->eta, MCparticle->phi, iObject.eta, iObject.phi);
		if((thisDeltaR <= iMaxDeltaR) && (thisDeltaR < minDeltaR)){ result = (*MCparticle); }
		minDeltaR = std::min(minDeltaR,thisDeltaR);
	}
	return result;
}

template <typename BNobject> BNgenjet BEANhelper::GetMatchedGenjet(const BNgenjetCollection& iGenjets, const BNobject& iObject, const double iMaxDeltaR){
	BNgenjet result;
	double minDeltaR = 999;
	for( BNgenjetCollection::const_iterator Genjet = iGenjets.begin(); Genjet != iGenjets.end(); ++Genjet){
		double thisDeltaR = deltaR(Genjet->eta, Genjet->phi, iObject.eta, iObject.phi);	
		if((thisDeltaR <= iMaxDeltaR) && (thisDeltaR < minDeltaR)){ result = (*Genjet); }
		minDeltaR = std::min(minDeltaR,thisDeltaR);
	}
	return result;
}

template <typename BNobject> BNmcparticle BEANhelper::GetMatchedGentau(const BNmcparticleCollection& iMCparticles, const BNobject& iObject, const double iMaxDeltaR){
	BNmcparticle result;
	double minDeltaR = 999;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle){
		if(abs(MCparticle->id) != 15){ continue; }
		double thisDeltaR = deltaR(MCparticle->eta, MCparticle->phi, iObject.eta, iObject.phi);	
		if((thisDeltaR <= iMaxDeltaR) && (thisDeltaR < minDeltaR)){ result = (*MCparticle); }
		minDeltaR = std::min(minDeltaR,thisDeltaR);
	}
	return result;
}

template <typename BNobject> BNmcparticle BEANhelper::GetMatchedVisGentau(const BNmcparticleCollection& iMCparticles, const BNobject& iObject, const double iMaxDeltaR){
	BNmcparticle result;
	double minDeltaR = 999;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle){
		if(abs(MCparticle->id) != 15){ continue; }
		BNmcparticle visGenTau = GetVisGenTau(*MCparticle, iMCparticles);
		double thisDeltaR = deltaR(visGenTau.eta, visGenTau.phi, iObject.eta, iObject.phi);	
		if((thisDeltaR <= iMaxDeltaR) && (thisDeltaR < minDeltaR)){ result = visGenTau; }
		minDeltaR = std::min(minDeltaR,thisDeltaR);
	}
	return result;
}

template <typename BNobject> BNtau BEANhelper::GetMatchedTau(const BNtauCollection& iTaus, const BNobject& iObject, const double iMaxDeltaR){
	BNtau result;
	double minDeltaR = 999;
	for( BNtauCollection::const_iterator Tau = iTaus.begin(); Tau != iTaus.end(); ++Tau){
		double thisDeltaR = deltaR(Tau->eta, Tau->phi, iObject.eta, iObject.phi);	
		if((thisDeltaR <= iMaxDeltaR) && (thisDeltaR < minDeltaR)){ result = (*Tau); }
		minDeltaR = std::min(minDeltaR,thisDeltaR);
	}
	return result;
}

template <typename BNobject> BNjet BEANhelper::GetClosestJet(const BNjetCollection& iJets, const BNobject& iObject, const double iMaxDeltaR){
	BNjet result;
	double minDeltaR = 999;
	for( BNjetCollection::const_iterator Jet = iJets.begin(); Jet != iJets.end(); ++Jet){
		double thisDeltaR = deltaR(Jet->eta, Jet->phi, iObject.eta, iObject.phi);	
		if((thisDeltaR <= iMaxDeltaR) && (thisDeltaR < minDeltaR)){ result = (*Jet); }
		minDeltaR = std::min(minDeltaR,thisDeltaR);
	}
	return result;
}

template <typename BNobject> BNtrack BEANhelper::GetClosestTrack(const BNtrackCollection& iTracks, const BNobject& iObject, const double iMaxDeltaR){
	BNtrack result;
	double minDeltaR = 999;
	for( BNtrackCollection::const_iterator Track = iTracks.begin(); Track != iTracks.end(); ++Track){
		double thisDeltaR = deltaR(Track->eta, Track->phi, iObject.eta, iObject.phi);	
		if((thisDeltaR <= iMaxDeltaR) && (thisDeltaR < minDeltaR)){ result = (*Track); }
		minDeltaR = std::min(minDeltaR,thisDeltaR);
	}
	return result;
}


// === Returned sorted input collection, by descending pT === //
template <typename BNcollection> BNcollection BEANhelper::GetSortedByPt(const BNcollection& iBNcollection){
	BNcollection result;
	BNcollection tempCollection = iBNcollection;

	if (tempCollection.size() == 0) { return tempCollection; }

	while(tempCollection.size() > 1){
		typename BNcollection::iterator largestPtElement = tempCollection.begin();
		for(typename BNcollection::iterator Object = (tempCollection.begin()+1); Object != tempCollection.end(); ++Object ){
			if(Object->pt > largestPtElement->pt){ largestPtElement = Object; }
		}

		result.push_back(*largestPtElement);
		tempCollection.erase(largestPtElement);
	}
	result.push_back(*(tempCollection.begin()));

	return result;
}

// === Returned sorted input collection, by descending CSV === //
template <typename BNcollection> BNcollection BEANhelper::GetSortedByCSV(const BNcollection& iBNcollection){
	BNcollection result;
	BNcollection tempCollection = iBNcollection;

	if (tempCollection.size() == 0) { return tempCollection; }

	while(tempCollection.size() > 1){
		typename BNcollection::iterator largestCSVElement = tempCollection.begin();
		for(typename BNcollection::iterator Object = (tempCollection.begin()+1); Object != tempCollection.end(); ++Object ){
			if(Object->btagCombinedSecVertex > largestCSVElement->btagCombinedSecVertex){ largestCSVElement = Object; }
		}

		result.push_back(*largestCSVElement);
		tempCollection.erase(largestCSVElement);
	}
	result.push_back(*(tempCollection.begin()));

	return result;
}

// === Return the union of the two input collections, removing the overlap === //
template <typename BNcollection> BNcollection BEANhelper::GetUnion(const BNcollection& iBNcollection1, const BNcollection& iBNcollection2){
	// Start off by adding all the objects in the first collection to the result
	BNcollection result = iBNcollection1;

	// Check to see what objects in the second collection need to be added, and do so
	for(typename BNcollection::const_iterator Object2 = iBNcollection2.begin(); Object2 != iBNcollection2.end(); ++Object2 ){

		// Look for Object2 in the results collection
		bool presentInFirstCollection = false;
		for(typename BNcollection::const_iterator Object1 = result.begin(); Object1 != result.end(); ++Object1 ){

			// If two objects match in deltaR, check that they have virtually the same momentum. Throw fatal error
			// if we get two matching objects with different momenta, as this probably mean that we're mixing corrected and uncorrected
			// input collections
			if(deltaR(Object1->eta, Object1->phi, Object2->eta, Object2->phi) < 0.00001){
				presentInFirstCollection = true;
				bool sameMomentum = (fabs(Object1->px - Object2->px) < 0.00001) &&
									(fabs(Object1->py - Object2->py) < 0.00001) &&
									(fabs(Object1->pz - Object2->pz) < 0.00001);
				if(!sameMomentum){ cerr << "ERROR: found two objects with same eta and phi, but different momenta. This may be caused by mixing corrected and uncorrected collections." << endl; 
				cout << setprecision(7) << "Eta1: " << Object1->eta << "\tPhi1: " << Object1->phi << "\tpT1: " << Object1->pt << endl;
				cout << setprecision(7) << "Eta2: " << Object2->eta << "\tPhi2: " << Object2->phi << "\tpT2: " << Object2->pt << endl;
				throw std::logic_error("Inside GetUnion"); }
			}
				
		}

		// If after looping over the results collection, Object2 wasn't found, add it!
		if(!presentInFirstCollection){ result.push_back(*Object2); }
	}

	// Sort by descending pT
	return GetSortedByPt(result);
}

// === Return the union of the two input collections, removing the overlap, with collection1 objects first, then collection2 === //
template <typename BNcollection> BNcollection BEANhelper::GetUnionUnsorted(const BNcollection& iBNcollection1, const BNcollection& iBNcollection2){
	// Start off by adding all the objects in the first collection to the result
	BNcollection result = iBNcollection1;

	// Check to see what objects in the second collection need to be added, and do so
	for(typename BNcollection::const_iterator Object2 = iBNcollection2.begin(); Object2 != iBNcollection2.end(); ++Object2 ){

		// Look for Object2 in the results collection
		bool presentInFirstCollection = false;
		for(typename BNcollection::const_iterator Object1 = result.begin(); Object1 != result.end(); ++Object1 ){

			// If two objects match in deltaR, check that they have virtually the same momentum. Throw fatal error
			// if we get two matching objects with different momenta, as this probably mean that we're mixing corrected and uncorrected
			// input collections
			if(deltaR(Object1->eta, Object1->phi, Object2->eta, Object2->phi) < 0.00001){
				presentInFirstCollection = true;
				bool sameMomentum = (fabs(Object1->px - Object2->px) < 0.00001) &&
									(fabs(Object1->py - Object2->py) < 0.00001) &&
									(fabs(Object1->pz - Object2->pz) < 0.00001);
				if(!sameMomentum){ cerr << "ERROR: found two objects with same eta and phi, but different momenta. This may be caused by mixing corrected and uncorrected collections." << endl; 
				cout << setprecision(7) << "Eta1: " << Object1->eta << "\tPhi1: " << Object1->phi << "\tpT1: " << Object1->pt << endl;
				cout << setprecision(7) << "Eta2: " << Object2->eta << "\tPhi2: " << Object2->phi << "\tpT2: " << Object2->pt << endl;
				throw std::logic_error("Inside GetUnionUnsorted"); }
			}
				
		}

		// If after looping over the results collection, Object2 wasn't found, add it!
		if(!presentInFirstCollection){ result.push_back(*Object2); }
	}

	return result;
}

// === Return the intersection of the two input collections, sorted by descending pT === //
template <typename BNcollection> BNcollection BEANhelper::GetIntersection(const BNcollection& iBNcollection1, const BNcollection& iBNcollection2){
	// Start off with an empty collection
	BNcollection result;

	// Check to see what objects in the second collection need to be added, and do so
	for(typename BNcollection::const_iterator Object2 = iBNcollection2.begin(); Object2 != iBNcollection2.end(); ++Object2 ){

		// Look for Object2 in the results collection
		bool presentInFirstCollection = false;
		for(typename BNcollection::const_iterator Object1 = iBNcollection1.begin(); Object1 != iBNcollection1.end(); ++Object1 ){

			// If two objects match in deltaR, check that they have virtually the same momentum. Throw fatal error
			// if we get two matching objects with different momenta, as this probably mean that we're mixing corrected and uncorrected
			// input collections
			if(deltaR(Object1->eta, Object1->phi, Object2->eta, Object2->phi) < 0.00001){
				presentInFirstCollection = true;
				bool sameMomentum = (fabs(Object1->px - Object2->px) < 0.00001) &&
									(fabs(Object1->py - Object2->py) < 0.00001) &&
									(fabs(Object1->pz - Object2->pz) < 0.00001);
				if(!sameMomentum){ cerr << "ERROR: found two objects with same eta and phi, but different momenta. This may be caused by mixing corrected and uncorrected collections." << endl; 
				cout << setprecision(7) << "Eta1: " << Object1->eta << "\tPhi1: " << Object1->phi << "\tpT1: " << Object1->pt << endl;
				cout << setprecision(7) << "Eta2: " << Object2->eta << "\tPhi2: " << Object2->phi << "\tpT2: " << Object2->pt << endl;
				throw std::logic_error("Inside GetIntersection"); }

				// If found a match, break loop
				break;
			}
				
		}

		// If found match, add it!
		if(presentInFirstCollection){ result.push_back(*Object2); }
	}

	// Sort by descending pT
	return GetSortedByPt(result);
}

// === Return the difference of the two input collections, sorted by descending pT === //
template <typename BNcollection> BNcollection BEANhelper::GetDifference(const BNcollection& iBNcollection1, const BNcollection& iBNcollection2){
	// Start off with an empty collection
	BNcollection result;

	// Check to see what objects in the first collection need to be added, and do so
	for(typename BNcollection::const_iterator Object1 = iBNcollection1.begin(); Object1 != iBNcollection1.end(); ++Object1 ){

		// Look for Object1 in the second collection
		bool presentInSecondCollection = false;
		for(typename BNcollection::const_iterator Object2 = iBNcollection2.begin(); Object2 != iBNcollection2.end(); ++Object2 ){

			// If two objects match in deltaR, check that they have virtually the same momentum. Throw fatal error
			// if we get two matching objects with different momenta, as this probably mean that we're mixing corrected and uncorrected
			// input collections
			if(deltaR(Object1->eta, Object1->phi, Object2->eta, Object2->phi) < 0.00001){
				presentInSecondCollection = true;
				bool sameMomentum = (fabs(Object1->px - Object2->px) < 0.00001) &&
									(fabs(Object1->py - Object2->py) < 0.00001) &&
									(fabs(Object1->pz - Object2->pz) < 0.00001);
				if(!sameMomentum){ cerr << "ERROR: found two objects with same eta and phi, but different momenta. This may be caused by mixing corrected and uncorrected collections." << endl; 
				cout << setprecision(7) << "Eta1: " << Object1->eta << "\tPhi1: " << Object1->phi << "\tpT1: " << Object1->pt << endl;
				cout << setprecision(7) << "Eta2: " << Object2->eta << "\tPhi2: " << Object2->phi << "\tpT2: " << Object2->pt << endl;
				throw std::logic_error("Inside GetDifference"); }

				break;
			}
		}

		// If no match found, add it!
		if(!presentInSecondCollection){ result.push_back(*Object1); }

	}

	// Sort by descending pT
	return GetSortedByPt(result);
}

// === Return the symmetric difference of the two input collections, sorted by descending pT === //
template <typename BNcollection> BNcollection BEANhelper::GetSymmetricDifference(const BNcollection& iBNcollection1, const BNcollection& iBNcollection2){
	return GetUnion(GetDifference(iBNcollection1, iBNcollection2), GetDifference(iBNcollection2, iBNcollection1));
}

// === Return the difference of the two input collections of different type, sorted by descending pT === //
template <typename BNcollection1, typename BNcollection2> BNcollection1 BEANhelper::GetDifference(const BNcollection1& iBNcollection1, const BNcollection2& iBNcollection2, const double iDeltaR){
	// Start off with an empty collection
	BNcollection1 result;

	// Check to see what objects in the first collection need to be added, and do so
	for(typename BNcollection1::const_iterator Object1 = iBNcollection1.begin(); Object1 != iBNcollection1.end(); ++Object1 ){

		// Look for Object1 in the second collection
		bool presentInSecondCollection = false;
		for(typename BNcollection2::const_iterator Object2 = iBNcollection2.begin(); Object2 != iBNcollection2.end(); ++Object2 ){

			// If two objects match in deltaR, check that they have virtually the same momentum. Throw fatal error
			// if we get two matching objects with different momenta, as this probably mean that we're mixing corrected and uncorrected
			// input collections
			if(deltaR(Object1->eta, Object1->phi, Object2->eta, Object2->phi) < iDeltaR){
				presentInSecondCollection = true;
				break;
			}
		}

		// If no match found, add it!
		if(!presentInSecondCollection){ result.push_back(*Object1); }

	}

	// Sort by descending pT
	return GetSortedByPt(result);
}


// === Print basic object info === //
template <typename BNobject> void BEANhelper::PrintInfo(const BNobject& iBNobject){
	cout << setprecision(4) << ">>> pT: " << setfill(' ') << setw(7) << iBNobject.pt << "\teta: " << iBNobject.eta << "\tphi: " << iBNobject.phi << endl;
}

#endif // _BEANhelper_h
