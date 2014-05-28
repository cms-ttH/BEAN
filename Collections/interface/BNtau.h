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
    HPSagainstMuonLoose(-99),
    HPSagainstMuonMedium(-99),
    HPSagainstMuonTight(-99),
    HPSagainstMuonLoose2(-99),
    HPSagainstMuonMedium2(-99),
    HPSagainstMuonTight2(-99),
    HPSagainstMuonLoose3(-99),
    HPSagainstMuonMedium3(-99),
    HPSagainstMuonTight3(-99),
    HPSagainstElectronTightMVA3(-99),
    HPSagainstElectronTightMVA2(-99), HPSbyLooseIsolationMVA2(-99),
    HPSagainstElectronMedium(-99),
    HPSbyTightIsolationMVA(-99),
    HPSbyLooseCombinedIsolationDeltaBetaCorr(-99),
    HPSagainstElectronLooseMVA3(-99), HPSagainstElectronLooseMVA2(-99),
    HPSagainstElectronTight(-99),
    HPSbyVLooseCombinedIsolationDeltaBetaCorr(-99),
    HPSagainstElectronVTightMVA3(-99), HPSagainstElectronMediumMVA3(-99),
    HPSagainstElectronMediumMVA2(-99), HPSagainstElectronMVA(-99),
    HPSbyMediumCombinedIsolationDeltaBetaCorr(-99),
    HPSagainstElectronVLooseMVA2(-99), HPSagainstElectronLoose(-99),
    HPSbyMediumIsolationMVA(-99),
    HPSbyMediumCombinedIsolationDeltaBetaCorr3Hits(-99),
    HPSbyMediumIsolationMVA2(-99), HPSbyTightIsolationMVA2(-99),
    HPSagainstElectronMVA2category(-99), HPSagainstElectronDeadECAL(-99),
    HPSagainstElectronMVA3category(-99),
    HPSbyTightCombinedIsolationDeltaBetaCorr(-99),
    HPSbyLooseCombinedIsolationDeltaBetaCorr3Hits(-99),
    HPSbyLooseIsolationMVA(-99),
    HPSbyTightCombinedIsolationDeltaBetaCorr3Hits(-99),
    HPSdecayModeFinding(-99), HPSbyIsolationMVAraw(-99),
    HPSbyCombinedIsolationDeltaBetaCorrRaw3Hits(-99),
    HPSbyIsolationMVA2raw(-99), HPSagainstElectronMVA2raw(-99),
    HPSagainstElectronMVA3raw(-99),
    HPSbyCombinedIsolationDeltaBetaCorrRaw(-99) {};

  double px, py, pz, energy, et, pt, eta, phi;
  double emFraction;
  int leadingTrackValid;
  double leadingTrackPt, leadingTrackIpVtdxy, leadingTrackIpVtdz, leadingTrackIpVtdxyError, leadingTrackIpVtdzError, leadingTrackVx, leadingTrackVy, leadingTrackVz, leadingTrackValidHits, leadingTrackNormChiSqrd;
  int numProngs, numSignalGammas, numSignalNeutrals, numSignalPiZeros;
  int decayMode, charge;
  int inTheCracks;

  int HPSagainstMuonLoose;
  int HPSagainstMuonMedium;
  int HPSagainstMuonTight;
  int HPSagainstMuonLoose2;
  int HPSagainstMuonMedium2;
  int HPSagainstMuonTight2;
  int HPSagainstMuonLoose3;
  int HPSagainstMuonMedium3;
  int HPSagainstMuonTight3;

  int HPSagainstElectronTightMVA3, HPSagainstElectronTightMVA2;
  int HPSbyLooseIsolationMVA2;
  int HPSagainstElectronMedium;
  int HPSbyTightIsolationMVA;
  int HPSbyLooseCombinedIsolationDeltaBetaCorr;
  int HPSagainstElectronLooseMVA3;
  int HPSagainstElectronLooseMVA2;
  int HPSagainstElectronTight;
  int HPSbyVLooseCombinedIsolationDeltaBetaCorr;
  int HPSagainstElectronVTightMVA3;
  int HPSagainstElectronMediumMVA3;
  int HPSagainstElectronMediumMVA2;
  int HPSagainstElectronMVA;
  int HPSbyMediumCombinedIsolationDeltaBetaCorr;
  int HPSagainstElectronVLooseMVA2;
  int HPSagainstElectronLoose;
  int HPSbyMediumIsolationMVA;
  int HPSbyMediumCombinedIsolationDeltaBetaCorr3Hits;
  int HPSbyMediumIsolationMVA2;
  int HPSbyTightIsolationMVA2;
  int HPSagainstElectronMVA2category;
  int HPSagainstElectronDeadECAL;
  int HPSagainstElectronMVA3category;
  int HPSbyTightCombinedIsolationDeltaBetaCorr;
  int HPSbyLooseCombinedIsolationDeltaBetaCorr3Hits;
  int HPSbyLooseIsolationMVA;
  int HPSbyTightCombinedIsolationDeltaBetaCorr3Hits;
  int HPSdecayModeFinding;
  double HPSbyIsolationMVAraw;
  double HPSbyCombinedIsolationDeltaBetaCorrRaw3Hits;
  double HPSbyIsolationMVA2raw;
  double HPSagainstElectronMVA2raw;
  double HPSagainstElectronMVA3raw;
  double HPSbyCombinedIsolationDeltaBetaCorrRaw;
};

typedef std::vector<BNtau> BNtauCollection;

#endif
