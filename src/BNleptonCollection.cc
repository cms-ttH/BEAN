#include "../interface/BNleptonCollection.h"

using namespace std;

BNleptonCollection::~BNleptonCollection(){
	// Clean up
	for(auto it = begin(); it != end(); ++it){ delete *it; }
}

void BNleptonCollection::sort(){
	for(BNleptonCollection::iterator it1 = begin(); it1 != end(); ++it1){

		BNleptonCollection::iterator largestPtIt = it1;
		for(BNleptonCollection::iterator it2 = (largestPtIt); it2 != end(); ++it2){
			if((*it2)->pt > (*largestPtIt)->pt){ largestPtIt = it2; }
		}

		BNlepton* tempObject = *it1;
		*it1 = *largestPtIt;
		*largestPtIt = tempObject;

	}
}

void BNleptonCollection::push_back(BNelectron* iElectron){
	vector<BNlepton*>::push_back(new BNelectron(*iElectron));
	back()->isElectron = 1;
	back()->isMuon = 0;
}

void BNleptonCollection::push_back(BNmuon* iMuon){
	vector<BNlepton*>::push_back(new BNmuon(*iMuon));
	back()->isElectron = 0;
	back()->isMuon = 1;
}

void BNleptonCollection::push_back(BNelectron const & iElectron){
	vector<BNlepton*>::push_back(new BNelectron(iElectron));
	back()->isElectron = 1;
	back()->isMuon = 0;
}

void BNleptonCollection::push_back(BNmuon const & iMuon){
	vector<BNlepton*>::push_back(new BNmuon(iMuon));
	back()->isElectron = 0;
	back()->isMuon = 1;
}

void BNleptonCollection::push_back(BNelectronCollection const & iElectrons){
	for(BNelectronCollection::const_iterator it = iElectrons.begin(); it != iElectrons.end(); ++it){ push_back(*it); }
}

void BNleptonCollection::push_back(BNmuonCollection const & iMuons){
	for(BNmuonCollection::const_iterator it = iMuons.begin(); it != iMuons.end(); ++it){ push_back(*it); }
}
