#ifndef ProductArea_BNtau_h
#define ProductArea_BNtau_h

#include "DataFormats/Math/interface/Point3D.h"

#include <vector>

namespace pat {
  class Tau;
}

struct BNtau {
  explicit BNtau(const pat::Tau&, math::XYZPoint&);

  BNtau() :
    px(-99), py(-99), pz(-99), energy(-99), et(-99), pt(-99), eta(-99),
    phi(-99), emFraction(-99), leadingTrackValid(-99), leadingTrackPt(-99),
    leadingTrackIpVtdxy(-99), leadingTrackIpVtdz(-99),
    leadingTrackIpVtdxyError(-99), leadingTrackIpVtdzError(-99),
    leadingTrackVx(-99), leadingTrackVy(-99), leadingTrackVz(-99),
    leadingTrackValidHits(-99), leadingTrackNormChiSqrd(-99),
    numProngs(-99), numSignalGammas(-99), numSignalNeutrals(-99),
    numSignalPiZeros(-99), decayMode(-99), charge(-99), inTheCracks(-99),
    HPSagainstElectronDeadECAL(-99),
    HPSagainstElectronLoose(-99),
    HPSagainstElectronMedium(-99),
    HPSagainstElectronTight(-99),
    HPSagainstElectronVLooseMVA5(-99),
    HPSagainstElectronLooseMVA5(-99),
    HPSagainstElectronMediumMVA5(-99),
    HPSagainstElectronTightMVA5(-99),
    HPSagainstElectronVTightMVA5(-99),
    HPSagainstElectronMVA5category(-99),
    HPSagainstElectronMVA5raw(-99),
    HPSagainstMuonLoose(-99),
    HPSagainstMuonMedium(-99),
    HPSagainstMuonTight(-99),
    HPSagainstMuonLoose2(-99),
    HPSagainstMuonMedium2(-99),
    HPSagainstMuonTight2(-99),
    HPSagainstMuonLoose3(-99),
    HPSagainstMuonMedium3(-99),
    HPSagainstMuonTight3(-99),
    HPSagainstMuonLooseMVA(-99),
    HPSagainstMuonMediumMVA(-99),
    HPSagainstMuonTightMVA(-99),
    HPSagainstMuonMVAraw(-99),
    HPSbyVLooseIsolationMVA3newDMwLT(-99),
    HPSbyLooseIsolationMVA3newDMwLT(-99),
    HPSbyMediumIsolationMVA3newDMwLT(-99),
    HPSbyTightIsolationMVA3newDMwLT(-99),
    HPSbyVTightIsolationMVA3newDMwLT(-99),
    HPSbyVVTightIsolationMVA3newDMwLT(-99),
    HPSbyIsolationMVA3newDMwLTraw(-99),
    HPSbyVLooseIsolationMVA3newDMwoLT(-99),
    HPSbyLooseIsolationMVA3newDMwoLT(-99),
    HPSbyMediumIsolationMVA3newDMwoLT(-99),
    HPSbyTightIsolationMVA3newDMwoLT(-99),
    HPSbyVTightIsolationMVA3newDMwoLT(-99),
    HPSbyVVTightIsolationMVA3newDMwoLT(-99),
    HPSbyIsolationMVA3newDMwoLTraw(-99),
    HPSbyVLooseIsolationMVA3oldDMwLT(-99),
    HPSbyLooseIsolationMVA3oldDMwLT(-99),
    HPSbyMediumIsolationMVA3oldDMwLT(-99),
    HPSbyTightIsolationMVA3oldDMwLT(-99),
    HPSbyVTightIsolationMVA3oldDMwLT(-99),
    HPSbyVVTightIsolationMVA3oldDMwLT(-99),
    HPSbyIsolationMVA3oldDMwLTraw(-99),
    HPSbyVLooseIsolationMVA3oldDMwoLT(-99),
    HPSbyLooseIsolationMVA3oldDMwoLT(-99),
    HPSbyMediumIsolationMVA3oldDMwoLT(-99),
    HPSbyTightIsolationMVA3oldDMwoLT(-99),
    HPSbyVTightIsolationMVA3oldDMwoLT(-99),
    HPSbyVVTightIsolationMVA3oldDMwoLT(-99),
    HPSbyIsolationMVA3oldDMwoLTraw(-99),
    HPSbyVLooseCombinedIsolationDeltaBetaCorr(-99),
    HPSbyLooseCombinedIsolationDeltaBetaCorr(-99),
    HPSbyMediumCombinedIsolationDeltaBetaCorr(-99),
    HPSbyTightCombinedIsolationDeltaBetaCorr(-99),
    HPSbyCombinedIsolationDeltaBetaCorrRaw(-99),
    HPSbyLooseCombinedIsolationDeltaBetaCorr3Hits(-99),
    HPSbyMediumCombinedIsolationDeltaBetaCorr3Hits(-99),
    HPSbyTightCombinedIsolationDeltaBetaCorr3Hits(-99),
    HPSbyCombinedIsolationDeltaBetaCorrRaw3Hits(-99),
    HPSdecayModeFinding(-99),
    HPSdecayModeFindingNewDMs(-99),
    HPSdecayModeFindingOldDMs(-99) {};

  double px, py, pz, energy, et, pt, eta, phi;
  double emFraction;
  int leadingTrackValid;
  double leadingTrackPt, leadingTrackIpVtdxy, leadingTrackIpVtdz, leadingTrackIpVtdxyError, leadingTrackIpVtdzError, leadingTrackVx, leadingTrackVy, leadingTrackVz, leadingTrackValidHits, leadingTrackNormChiSqrd;
  int numProngs, numSignalGammas, numSignalNeutrals, numSignalPiZeros;
  int decayMode, charge;
  int inTheCracks;

  int HPSagainstElectronDeadECAL;

  int HPSagainstElectronLoose;
  int HPSagainstElectronMedium;
  int HPSagainstElectronTight;

  int HPSagainstElectronVLooseMVA5;
  int HPSagainstElectronLooseMVA5;
  int HPSagainstElectronMediumMVA5;
  int HPSagainstElectronTightMVA5;
  int HPSagainstElectronVTightMVA5;
  int HPSagainstElectronMVA5category;
  double HPSagainstElectronMVA5raw;

  int HPSagainstMuonLoose;
  int HPSagainstMuonMedium;
  int HPSagainstMuonTight;

  int HPSagainstMuonLoose2;
  int HPSagainstMuonMedium2;
  int HPSagainstMuonTight2;

  int HPSagainstMuonLoose3;
  int HPSagainstMuonMedium3;
  int HPSagainstMuonTight3;

  int HPSagainstMuonLooseMVA;
  int HPSagainstMuonMediumMVA;
  int HPSagainstMuonTightMVA;
  double HPSagainstMuonMVAraw;

  int HPSbyVLooseIsolationMVA3newDMwLT;
  int HPSbyLooseIsolationMVA3newDMwLT;
  int HPSbyMediumIsolationMVA3newDMwLT;
  int HPSbyTightIsolationMVA3newDMwLT;
  int HPSbyVTightIsolationMVA3newDMwLT;
  int HPSbyVVTightIsolationMVA3newDMwLT;
  double HPSbyIsolationMVA3newDMwLTraw;

  int HPSbyVLooseIsolationMVA3newDMwoLT;
  int HPSbyLooseIsolationMVA3newDMwoLT;
  int HPSbyMediumIsolationMVA3newDMwoLT;
  int HPSbyTightIsolationMVA3newDMwoLT;
  int HPSbyVTightIsolationMVA3newDMwoLT;
  int HPSbyVVTightIsolationMVA3newDMwoLT;
  double HPSbyIsolationMVA3newDMwoLTraw;

  int HPSbyVLooseIsolationMVA3oldDMwLT;
  int HPSbyLooseIsolationMVA3oldDMwLT;
  int HPSbyMediumIsolationMVA3oldDMwLT;
  int HPSbyTightIsolationMVA3oldDMwLT;
  int HPSbyVTightIsolationMVA3oldDMwLT;
  int HPSbyVVTightIsolationMVA3oldDMwLT;
  double HPSbyIsolationMVA3oldDMwLTraw;

  int HPSbyVLooseIsolationMVA3oldDMwoLT;
  int HPSbyLooseIsolationMVA3oldDMwoLT;
  int HPSbyMediumIsolationMVA3oldDMwoLT;
  int HPSbyTightIsolationMVA3oldDMwoLT;
  int HPSbyVTightIsolationMVA3oldDMwoLT;
  int HPSbyVVTightIsolationMVA3oldDMwoLT;
  double HPSbyIsolationMVA3oldDMwoLTraw;

  int HPSbyVLooseCombinedIsolationDeltaBetaCorr;
  int HPSbyLooseCombinedIsolationDeltaBetaCorr;
  int HPSbyMediumCombinedIsolationDeltaBetaCorr;
  int HPSbyTightCombinedIsolationDeltaBetaCorr;
  double HPSbyCombinedIsolationDeltaBetaCorrRaw;

  int HPSbyLooseCombinedIsolationDeltaBetaCorr3Hits;
  int HPSbyMediumCombinedIsolationDeltaBetaCorr3Hits;
  int HPSbyTightCombinedIsolationDeltaBetaCorr3Hits;
  double HPSbyCombinedIsolationDeltaBetaCorrRaw3Hits;

  int HPSdecayModeFinding;
  int HPSdecayModeFindingNewDMs;
  int HPSdecayModeFindingOldDMs;
};

typedef std::vector<BNtau> BNtauCollection;

#endif
