#include "../interface/BNleptonCollection.h"

using namespace std;

BNleptonCollection::~BNleptonCollection(){
	// Clean up
	for(auto it = begin(); it != end(); ++it){ delete *it; }
}

BNleptonCollection::BNleptonCollection() : std::vector<BNlepton*>() {}

BNleptonCollection::BNleptonCollection(const BNleptonCollection& collection) : std::vector<BNlepton*>() {
  for (auto& particle: collection) {
    this->push_back(particle);
  }
}

BNleptonCollection& BNleptonCollection::operator=(const BNleptonCollection& collection) {

  for (auto it = begin(); it != end(); ++it) { delete *it; }
  clear();

  for (auto& particle: collection) {
    this->push_back(particle);
  }

  return *this;
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

void BNleptonCollection::push_back(BNlepton* iLepton){
	BNlepton* newLepton = NULL;
	if(iLepton->isElectron){
		newLepton = new BNelectron(*((BNelectron*)iLepton));
	}else if(iLepton->isMuon){
		newLepton = new BNmuon(*((BNmuon*)iLepton));
	}else{
		cerr << "ERROR: Unknown type of lepton." << endl; exit(1);	
	}

	vector<BNlepton*>::push_back(newLepton);
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

void BNleptonCollection::push_back(BNleptonCollection const & iLeptons){
	for(BNleptonCollection::const_iterator it = iLeptons.begin(); it != iLeptons.end(); ++it){ push_back(*it); }
}
