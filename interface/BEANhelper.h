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
#include "ProductArea/BNcollections/interface/BNelectron.h"
#include "ProductArea/BNcollections/interface/BNevent.h"
#include "ProductArea/BNcollections/interface/BNjet.h"
#include "ProductArea/BNcollections/interface/BNtau.h"
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
#include "NtupleMaker/BEANmaker/interface/CSVreevaluator.h"


#endif

typedef std::map<std::string, std::string> mparams;
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
typedef BNtauCollection::const_iterator           TauIter;
typedef BNphotonCollection::const_iterator        PhotonIter;
typedef BNprimaryvertexCollection::const_iterator PVIter;
typedef BNskimbitsCollection::const_iterator      SkimBitIter;
typedef BNsuperclusterCollection::const_iterator  SCIter;
typedef BNtrackCollection::const_iterator         TrackIter;
typedef BNtriggerCollection::const_iterator       TrigIter;
typedef BNtrigobjCollection::const_iterator       TrigObjIter;

namespace sysType{		enum sysType{		NA, JERup, JERdown, JESup, JESdown, hfSFup, hfSFdown, lfSFdown, lfSFup, TESup, TESdown }; }
namespace jetID{		enum jetID{			jetMinimal, jetLooseAOD, jetLoose, jetTight }; }
namespace tauID{		enum tauID{			tauVLoose, tauLoose, tauMedium, tauTight }; }
namespace muonID{		enum muonID{		muonSide, muonLoose, muonTight, muonPtOnly, muonPtEtaOnly, muonPtEtaIsoOnly, muonPtEtaIsoTrackerOnly }; }
namespace electronID{	enum electronID{	electronSide, electronLoose, electronTight }; }

using namespace std;

class BEANhelper{
	public:


	// === Functions === //
	public: 
		// Constructor(s) and destructor
		BEANhelper();
		virtual ~BEANhelper();
		
		// Set up BEANhelper
		void SetUp(string, int, bool, bool, string, bool, bool, string iCollisionDS="All");

		template <typename BNobject> void PrintInfo(const BNobject&);
		template <typename BNobject> BNmcparticle GetMatchedMCparticle(const BNmcparticleCollection&, const BNobject&, const double);

		// Union, intersection, difference
		template <typename BNcollection> BNcollection GetSortedByPt(const BNcollection&);
		template <typename BNcollection> BNcollection GetUnion(const BNcollection&, const BNcollection&);
		template <typename BNcollection> BNcollection GetIntersection(const BNcollection&, const BNcollection&);
		template <typename BNcollection> BNcollection GetDifference(const BNcollection&, const BNcollection&);
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
		BNjetCollection GetCleanJets(const BNjetCollection&, const vector<TLorentzVector>&, const float);
		unsigned int GetNumCSVbtags(const BNjetCollection&, const char);
		unsigned int GetNumNonCSVbtags(const BNjetCollection&, const char);
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

		// Muons
		bool IsSideMuon(const BNmuon&);
		bool IsLooseMuon(const BNmuon&);
		bool IsTightMuon(const BNmuon&);
		bool IsGoodMuon(const BNmuon&, const muonID::muonID);
		float GetMuonRelIso(const BNmuon&);
		float GetMuonSF(const BNmuon&);
		BNmuonCollection GetSelectedMuons(const BNmuonCollection&, const muonID::muonID);

		// Electrons
		bool IsSideElectron(const BNelectron&);
		bool IsLooseElectron(const BNelectron&);
		bool IsTightElectron(const BNelectron&);
		float GetElectronRelIso(const BNelectron&);
		float GetElectronSF(const BNelectron&);
		bool GetElectronIDresult(const BNelectron& iElectron, const electronID::electronID);
		bool IsGoodElectron(const BNelectron&, const electronID::electronID);
		BNelectronCollection GetSelectedElectrons(const BNelectronCollection&, const electronID::electronID);

		// MCparticles
		BNmcparticleCollection	GetSelectedMCparticlesByGrandParentPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetSelectedMCparticlesByParentPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetSelectedMCparticlesByPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetSelectedMCparticlesByChildPDGid(const BNmcparticleCollection&, const vector<int>);

		BNmcparticleCollection	GetUnrejectedMCparticlesByGrandParentPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetUnrejectedMCparticlesByParentPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetUnrejectedMCparticlesByPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetUnrejectedMCparticlesByChildPDGid(const BNmcparticleCollection&, const vector<int>);

		BNmcparticleCollection	GetSelectedMCparticlesByGrandParentStatus(const BNmcparticleCollection&, const bool, const bool, const bool);
		BNmcparticleCollection	GetSelectedMCparticlesByParentStatus(const BNmcparticleCollection&, const bool, const bool, const bool);
		BNmcparticleCollection	GetSelectedMCparticlesByStatus(const BNmcparticleCollection&, const bool, const bool, const bool);
		BNmcparticleCollection	GetSelectedMCparticlesByChildStatus(const BNmcparticleCollection&, const bool, const bool, const bool);

		BNmcparticleCollection	GetParents(const BNmcparticle&, const BNmcparticleCollection&);
		void 					DrawFeynman(const BNmcparticle&);
		BNmcparticleCollection	GetGenTaus(const BNmcparticleCollection&);
		BNmcparticleCollection	GetHadronicGenTaus(const BNmcparticleCollection&);
		BNmcparticle			GetVisGenTau(const BNmcparticle&, const BNmcparticleCollection&);
        bool                    ttPlusHeavyKeepEvent(const BNmcparticleCollection&, const BNjetCollection&);
        unsigned int            GetNumExtraPartons(const BNmcparticleCollection&);
  
		// PU reweighing
		double GetPUweight(const double);
		double GetPUweightUp(const double);
		double GetPUweightDown(const double);

		// Top PT reweighting
		double TopPtWeight(double);
		double GetTopPtweight(const BNmcparticleCollection&);
		double GetTopPtweightUp(const BNmcparticleCollection&);
		double GetTopPtweightDown(const BNmcparticleCollection&);


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
		void ThrowFatalError(const string);
		double GetCSVvalue(const BNjet&, const sysType::sysType iSysType=sysType::NA);
		void CheckSetUp();
		void SetUpPUreweighing(string const);
		void SetUpCSVreshaping();
		void SetUpJetSF();
		void SetUpLeptonSF();
		string GetSampleName();

		// Parameter management
	private:
		bool isSetUp;
		string era;
		int sampleNumber;
		bool isLJ;
		bool isData;
		string dataset;
		bool reshapeCSV;
        bool usePfLeptons;

		TFile*			jetSFfile;
		TFile*			leptonSFfile;

		// PU reweighing
		TFile*			puFile;
		TH1D*			h_PU_ratio;
		TH1D*			h_PUup_ratio;
		TH1D*			h_PUdown_ratio;

		// CSV reshaping
		CSVreevaluator*	  sh_;
		CSVreevaluator*	  sh_hfSFUp_;
		CSVreevaluator*	  sh_hfSFDown_;
		CSVreevaluator*	  sh_lfSFUp_;
		CSVreevaluator*	  sh_lfSFDown_;

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

		// Old parameters
		TH2D*			h_b_eff_;
		TH2D*			h_c_eff_;
		TH2D*			h_l_eff_;
		TH2D*			h_o_eff_;
        TH2D*           h_ele_SF_;
        TH2D*           h_mu_SF_;
        string          samplename;

}; // End of class prototype


// === Return matched BNmcparticle based on minimum deltaR and a deltaR threshold === //
template <typename BNobject> BNmcparticle BEANhelper::GetMatchedMCparticle(const BNmcparticleCollection& iMCparticles, const BNobject& iObject, const double iMaxDeltaR){
	BNmcparticle result;
	double minDeltaR = 999;
	for( BNmcparticleCollection::const_iterator MCparticle = iMCparticles.begin(); MCparticle != iMCparticles.end(); ++MCparticle ){
        double thisDeltaR = deltaR(MCparticle->eta, MCparticle->phi, iObject.eta, iObject.phi);	
		if((thisDeltaR <= iMaxDeltaR) && (thisDeltaR < minDeltaR)){ result = (*MCparticle); }
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
			if(deltaR(Object1->eta, Object1->phi, Object2->eta, Object2->phi) < 0.001){
				presentInFirstCollection = true;
				bool sameMomentum = (fabs(Object1->px - Object2->px) < 0.001) &&
									(fabs(Object1->py - Object2->py) < 0.001) &&
									(fabs(Object1->pz - Object2->pz) < 0.001);
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

// === Print basic object info === //
template <typename BNobject> void BEANhelper::PrintInfo(const BNobject& iBNobject){
	cout << setprecision(4) << ">>> pT: " << setfill(' ') << setw(7) << iBNobject.pt << "\teta: " << iBNobject.eta << "\tphi: " << iBNobject.phi << endl;
}

#endif // _BEANhelper_h
