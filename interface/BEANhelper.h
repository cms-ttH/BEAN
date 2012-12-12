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
//#include "NtupleMaker/BEANmaker/interface/BTagReshaping.h"
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
typedef BNphotonCollection::const_iterator        PhotonIter;
typedef BNprimaryvertexCollection::const_iterator PVIter;
typedef BNskimbitsCollection::const_iterator      SkimBitIter;
typedef BNsuperclusterCollection::const_iterator  SCIter;
typedef BNtrackCollection::const_iterator         TrackIter;
typedef BNtriggerCollection::const_iterator       TrigIter;
typedef BNtrigobjCollection::const_iterator       TrigObjIter;

namespace sysType{		enum sysType{		NA, JERup, JERdown, JESup, JESdown, hfSFup, hfSFdown, lfSFdown, lfSFup }; }
namespace jetID{		enum jetID{			jetMinimal, jetLooseAOD, jetLoose, jetTight }; }
namespace muonID{		enum muonID{		muonSide, muonLoose, muonTight }; }
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
  void SetUp(string, int, bool, bool, string, bool, bool);

		template <typename BNobject> void PrintInfo(const BNobject&);

		// Union, intersection, difference
		template <typename BNcollection> BNcollection GetSortedByPt(const BNcollection&);
		template <typename BNcollection> BNcollection GetUnion(const BNcollection&, const BNcollection&);
		template <typename BNcollection> BNcollection GetIntersection(const BNcollection&, const BNcollection&);
		template <typename BNcollection> BNcollection GetDifference(const BNcollection&, const BNcollection&);
		template <typename BNcollection> BNcollection GetSymmetricDifference(const BNcollection&, const BNcollection&);
		template <typename BNcollection> BNcollection GetUnionUnsorted(const BNcollection&, const BNcollection&);

		// Jets and MET
		float GetBtagWeight(const BNjet&, const sysType::sysType iSysType=sysType::NA);
		bool PassesCSV(const BNjet&, const char, const sysType::sysType iSysType=sysType::NA);
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
		BNmcparticleCollection	GetSelectedMCparticlesByPDGid(const BNmcparticleCollection&, const vector<int>);
		BNmcparticleCollection	GetGenTaus(const BNmcparticleCollection&);
		BNmcparticleCollection	GetHadronicGenTaus(const BNmcparticleCollection&);
		BNmcparticle			GetMatchedMCparticle(const BNmcparticleCollection&, const BNtau&, const double iMaxDeltaR=0.25);
        bool                    ttPlusHeavyKeepEvent(const BNmcparticleCollection&, const BNjetCollection&);  
  
		// PU reweighing
		double GetPUweight(const double);
		double GetPUweightUp(const double);
		double GetPUweightDown(const double);

	private:
		double GetCSVvalue(const BNjet&, const sysType::sysType iSysType=sysType::NA);
		void CheckSetUp();
		void ThrowFatalError(const string);

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

		// CSV reshaping
		/*
		BTagShapeInterface*	  sh_;
		BTagShapeInterface*	  sh_hfSFUp_;
		BTagShapeInterface*	  sh_hfSFDown_;
		BTagShapeInterface*	  sh_lfSFUp_;
		BTagShapeInterface*	  sh_lfSFDown_;
		//*/

		CSVreevaluator*	  sh_;
		CSVreevaluator*	  sh_hfSFUp_;
		CSVreevaluator*	  sh_hfSFDown_;
		CSVreevaluator*	  sh_lfSFUp_;
		CSVreevaluator*	  sh_lfSFDown_;

		// Old functions
	public:
        void setMCsample( int insample=2500, std::string era="", bool isLJ=true, std::string dset="" );
		void electronSelector( const BNelectronCollection &electrons, bool isLJ, std::string era, vint &tightElectrons, vint &looseElectrons );
		void muonSelector( const BNmuonCollection &muons, bool isLJ, std::string era, vint &tightMuons, vint &looseMuons );
		void jetSelector( const BNjetCollection &pfjets, std::string sysType, vint &tightJets, vint &tagJets, vint &untagJets, 
				std::vector<BTagWeight::JetInfo> &myjetinfo, double csvCut = 0.679 );
		void getPUwgt( double input_numPU, double &PU_scale, double &PUup_scale, double &PUdown_scale ); 
		vdouble getEffSF( int returnType, double jetPts, double jetEtas, double jetIds );
		double getJERfactor( int returnType, double jetAbsETA, double genjetPT, double recojetPT );
		void getSp(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, float &aplanarity, float &sphericity);
		void getFox(vecTLorentzVector jets,float &h0, float &h1, float &h2, float &h3, float &h4);	
		//void getFox_mod2(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, double HT,
		//float &h0_mod2, float &h1_mod2, float &h2_mod2, float &h3_mod2, float &h4_mod2,  float &h5_mod2,  
		//float &h6_mod2, float &h7_mod2, float &h8_mod2, float &h9_mod2, float &h10_mod2 );




	protected:

	private:


	// === Variables === //
	public:

	protected:

	private:
		mparams			params;
		float			CSVLwp, CSVMwp, CSVTwp;

		// Old parameters
		char *			my_pPath;
		string			my_base_dir;
		string			str_eff_file_7TeV, str_eff_file_8TeV, str_pu_file_7TeV, str_pu_file_8TeV;
		TH2D*			h_b_eff_;
		TH2D*			h_c_eff_;
		TH2D*			h_l_eff_;
		TH2D*			h_o_eff_;
		TH1D*			h_PU_ratio_;
		TH1D*			h_PUup_ratio_;
		TH1D*			h_PUdown_ratio_;
        TH2D*           h_ele_SF_;
        TH2D*           h_mu_SF_;
		double			PI;
		double			TWOPI;
		float			ETA_LIMIT;
		float			EPSILON;
        string          samplename;

}; // End of class prototype



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
				if(!sameMomentum){ cerr << "ERROR: found two objects with same eta and phi, but different momenta. This may be caused by mixing corrected and uncorrected collections." << endl; exit(1); }
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
			if(deltaR(Object1->eta, Object1->phi, Object2->eta, Object2->phi) < 0.001){
				presentInFirstCollection = true;
				bool sameMomentum = (fabs(Object1->px - Object2->px) < 0.001) &&
									(fabs(Object1->py - Object2->py) < 0.001) &&
									(fabs(Object1->pz - Object2->pz) < 0.001);
				if(!sameMomentum){ cerr << "ERROR: found two objects with same eta and phi, but different momenta. This may be caused by mixing corrected and uncorrected collections." << endl; exit(1); }
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
			if(deltaR(Object1->eta, Object1->phi, Object2->eta, Object2->phi) < 0.001){
				presentInFirstCollection = true;
				bool sameMomentum = (fabs(Object1->px - Object2->px) < 0.001) &&
									(fabs(Object1->py - Object2->py) < 0.001) &&
									(fabs(Object1->pz - Object2->pz) < 0.001);
				if(!sameMomentum){ cerr << "ERROR: found two objects with same eta and phi, but different momenta. This may be caused by mixing corrected and uncorrected collections." << endl; exit(1); }

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
			if(deltaR(Object1->eta, Object1->phi, Object2->eta, Object2->phi) < 0.001){
				presentInSecondCollection = true;
				bool sameMomentum = (fabs(Object1->px - Object2->px) < 0.001) &&
									(fabs(Object1->py - Object2->py) < 0.001) &&
									(fabs(Object1->pz - Object2->pz) < 0.001);
				if(!sameMomentum){ cerr << "ERROR: found two objects with same eta and phi, but different momenta. This may be caused by mixing corrected and uncorrected collections." << endl; exit(1); }

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
