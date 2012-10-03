#ifndef ProductArea_BNtau_h
#define ProductArea_BNtau_h

#include <vector>


// a simple class
struct BNtau {

  explicit BNtau(double v, int c): et(v), pt(v), eta(v), phi(v), emFraction(v), leadingTrackPt(v), leadingTrackIpVtdxy(v), leadingTrackIpVtdz(v), leadingTrackIpVtdxyError(v), leadingTrackIpVtdzError(v), leadingTrackVx(v), leadingTrackVy(v), leadingTrackVz(v), leadingTrackValidHits(v), leadingTrackNormChiSqrd(v), numProngs(c), numSignalGammas(c), numSignalNeutrals(c), numSignalPiZeros(c), decayMode(c), charge(c), inTheCracks(c), HPSagainstElectronLoose(c), HPSagainstElectronMVA(c), HPSagainstElectronMedium(c), HPSagainstElectronTight(c), HPSagainstMuonLoose(c), HPSagainstMuonMedium(c), HPSagainstMuonTight(c), HPSbyLooseCombinedIsolationDeltaBetaCorr(c), HPSbyMediumCombinedIsolationDeltaBetaCorr(c), HPSbyTightCombinedIsolationDeltaBetaCorr(c), HPSbyVLooseCombinedIsolationDeltaBetaCorr(c), HPSdecayModeFinding(c), leadingTrackValid(c) {}

  BNtau(): et(-99), pt(-99), eta(-99), phi(-99), emFraction(-99), leadingTrackPt(-99), leadingTrackIpVtdxy(-99), leadingTrackIpVtdz(-99), leadingTrackIpVtdxyError(-99), leadingTrackIpVtdzError(-99), leadingTrackVx(-99), leadingTrackVy(-99), leadingTrackVz(-99), leadingTrackValidHits(-99), leadingTrackNormChiSqrd(-99), numProngs(-99), numSignalGammas(-99), numSignalNeutrals(-99), numSignalPiZeros(-99), decayMode(-99), charge(-99), inTheCracks(-99), HPSagainstElectronLoose(-99), HPSagainstElectronMVA(-99), HPSagainstElectronMedium(-99), HPSagainstElectronTight(-99), HPSagainstMuonLoose(-99), HPSagainstMuonMedium(-99), HPSagainstMuonTight(-99), HPSbyLooseCombinedIsolationDeltaBetaCorr(-99), HPSbyMediumCombinedIsolationDeltaBetaCorr(-99), HPSbyTightCombinedIsolationDeltaBetaCorr(-99), HPSbyVLooseCombinedIsolationDeltaBetaCorr(-99), HPSdecayModeFinding(-99), leadingTrackValid(-99){}

	double px, py, pz, energy, et, pt, eta, phi, emFraction, leadingTrackPt, leadingTrackIpVtdxy, leadingTrackIpVtdz, leadingTrackIpVtdxyError, leadingTrackIpVtdzError, leadingTrackVx, leadingTrackVy, leadingTrackVz, leadingTrackValidHits, leadingTrackNormChiSqrd;
	int numProngs, numSignalGammas, numSignalNeutrals, numSignalPiZeros, decayMode, charge;
	int inTheCracks, HPSagainstElectronLoose, HPSagainstElectronMVA, HPSagainstElectronMedium, HPSagainstElectronTight, HPSagainstMuonLoose, HPSagainstMuonMedium, HPSagainstMuonTight, HPSbyLooseCombinedIsolationDeltaBetaCorr, HPSbyMediumCombinedIsolationDeltaBetaCorr, HPSbyTightCombinedIsolationDeltaBetaCorr, HPSbyVLooseCombinedIsolationDeltaBetaCorr, HPSdecayModeFinding, leadingTrackValid;

};

// this is our new product, it is simply a 
// collection of BNtau held in an std::vector
typedef std::vector<BNtau> BNtauCollection;

#endif
