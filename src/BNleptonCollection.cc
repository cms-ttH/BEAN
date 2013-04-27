#include "../interface/BNleptonCollection.h"

using namespace std;

BNleptonCollection::~BNleptonCollection(){
	// Clean up
	for(auto it = begin(); it != end(); ++it){ delete *it; }
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
