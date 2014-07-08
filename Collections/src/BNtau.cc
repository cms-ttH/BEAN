#include "DataFormats/PatCandidates/interface/Tau.h"

#include "../interface/BNtau.h"

inline bool
in_cracks(float etaValue)
{
  return (fabs(etaValue) < 0.018 ||
      (fabs(etaValue)>0.423 && fabs(etaValue)<0.461) ||
      (fabs(etaValue)>0.770 && fabs(etaValue)<0.806) ||
      (fabs(etaValue)>1.127 && fabs(etaValue)<1.163) ||
      (fabs(etaValue)>1.460 && fabs(etaValue)<1.558));
}


BNtau::BNtau(const pat::Tau& tau, math::XYZPoint& pv) {
  px= tau.px();
  py= tau.py();
  pz= tau.pz();
  pt= tau.pt();
  energy= tau.energy();
  et= tau.et();
  eta= tau.eta();
  phi= tau.phi();
  numProngs= tau.signalPFChargedHadrCands().size();
  numSignalGammas= tau.signalPFGammaCands().size();
  numSignalNeutrals= tau.signalPFNeutrHadrCands().size();
  numSignalPiZeros= tau.signalPiZeroCandidates().size();
  decayMode= tau.decayMode();
  emFraction= tau.emFraction();
  inTheCracks= in_cracks(tau.eta());

  HPSagainstElectronDeadECAL = tau.tauID("againstElectronDeadECAL");

  HPSagainstElectronLoose = tau.tauID("againstElectronLoose");
  HPSagainstElectronMedium = tau.tauID("againstElectronMedium");
  HPSagainstElectronTight = tau.tauID("againstElectronTight");

  HPSagainstElectronVLooseMVA5 = tau.tauID("againstElectronVLooseMVA5");
  HPSagainstElectronLooseMVA5 = tau.tauID("againstElectronLooseMVA5");
  HPSagainstElectronMediumMVA5 = tau.tauID("againstElectronMediumMVA5");
  HPSagainstElectronTightMVA5 = tau.tauID("againstElectronTightMVA5");
  HPSagainstElectronVTightMVA5 = tau.tauID("againstElectronVTightMVA5");
  HPSagainstElectronMVA5category = tau.tauID("againstElectronMVA5category");
  HPSagainstElectronMVA5raw = tau.tauID("againstElectronMVA5raw");

  HPSagainstMuonLoose = tau.tauID("againstMuonLoose");
  HPSagainstMuonMedium = tau.tauID("againstMuonMedium");
  HPSagainstMuonTight = tau.tauID("againstMuonTight");

  HPSagainstMuonLoose2 = tau.tauID("againstMuonLoose2");
  HPSagainstMuonMedium2 = tau.tauID("againstMuonMedium2");
  HPSagainstMuonTight2 = tau.tauID("againstMuonTight2");

  HPSagainstMuonLoose3 = tau.tauID("againstMuonLoose3");
  // HPSagainstMuonMedium3 = tau.tauID("againstMuonMedium3");
  HPSagainstMuonMedium3 = -99;
  HPSagainstMuonTight3 = tau.tauID("againstMuonTight3");

  HPSagainstMuonLooseMVA = tau.tauID("againstMuonLooseMVA");
  HPSagainstMuonMediumMVA = tau.tauID("againstMuonMediumMVA");
  HPSagainstMuonTightMVA = tau.tauID("againstMuonTightMVA");
  HPSagainstMuonMVAraw = tau.tauID("againstMuonMVAraw");

  HPSbyVLooseIsolationMVA3newDMwLT = tau.tauID("byVLooseIsolationMVA3newDMwLT");
  HPSbyLooseIsolationMVA3newDMwLT = tau.tauID("byLooseIsolationMVA3newDMwLT");
  HPSbyMediumIsolationMVA3newDMwLT = tau.tauID("byMediumIsolationMVA3newDMwLT");
  HPSbyTightIsolationMVA3newDMwLT = tau.tauID("byTightIsolationMVA3newDMwLT");
  HPSbyVTightIsolationMVA3newDMwLT = tau.tauID("byVTightIsolationMVA3newDMwLT");
  HPSbyVVTightIsolationMVA3newDMwLT = tau.tauID("byVVTightIsolationMVA3newDMwLT");
  HPSbyIsolationMVA3newDMwLTraw = tau.tauID("byIsolationMVA3newDMwLTraw");

  HPSbyVLooseIsolationMVA3newDMwoLT = tau.tauID("byVLooseIsolationMVA3newDMwoLT");
  HPSbyLooseIsolationMVA3newDMwoLT = tau.tauID("byLooseIsolationMVA3newDMwoLT");
  HPSbyMediumIsolationMVA3newDMwoLT = tau.tauID("byMediumIsolationMVA3newDMwoLT");
  HPSbyTightIsolationMVA3newDMwoLT = tau.tauID("byTightIsolationMVA3newDMwoLT");
  HPSbyVTightIsolationMVA3newDMwoLT = tau.tauID("byVTightIsolationMVA3newDMwoLT");
  HPSbyVVTightIsolationMVA3newDMwoLT = tau.tauID("byVVTightIsolationMVA3newDMwoLT");
  HPSbyIsolationMVA3newDMwoLTraw = tau.tauID("byIsolationMVA3newDMwoLTraw");

  HPSbyVLooseIsolationMVA3oldDMwLT = tau.tauID("byVLooseIsolationMVA3oldDMwLT");
  HPSbyLooseIsolationMVA3oldDMwLT = tau.tauID("byLooseIsolationMVA3oldDMwLT");
  HPSbyMediumIsolationMVA3oldDMwLT = tau.tauID("byMediumIsolationMVA3oldDMwLT");
  HPSbyTightIsolationMVA3oldDMwLT = tau.tauID("byTightIsolationMVA3oldDMwLT");
  HPSbyVTightIsolationMVA3oldDMwLT = tau.tauID("byVTightIsolationMVA3oldDMwLT");
  HPSbyVVTightIsolationMVA3oldDMwLT = tau.tauID("byVVTightIsolationMVA3oldDMwLT");
  HPSbyIsolationMVA3oldDMwLTraw = tau.tauID("byIsolationMVA3oldDMwLTraw");

  HPSbyVLooseIsolationMVA3oldDMwoLT = tau.tauID("byVLooseIsolationMVA3oldDMwoLT");
  HPSbyLooseIsolationMVA3oldDMwoLT = tau.tauID("byLooseIsolationMVA3oldDMwoLT");
  HPSbyMediumIsolationMVA3oldDMwoLT = tau.tauID("byMediumIsolationMVA3oldDMwoLT");
  HPSbyTightIsolationMVA3oldDMwoLT = tau.tauID("byTightIsolationMVA3oldDMwoLT");
  HPSbyVTightIsolationMVA3oldDMwoLT = tau.tauID("byVTightIsolationMVA3oldDMwoLT");
  HPSbyVVTightIsolationMVA3oldDMwoLT = tau.tauID("byVVTightIsolationMVA3oldDMwoLT");
  HPSbyIsolationMVA3oldDMwoLTraw = tau.tauID("byIsolationMVA3oldDMwoLTraw");

  HPSbyLooseCombinedIsolationDeltaBetaCorr3Hits = tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
  HPSbyMediumCombinedIsolationDeltaBetaCorr3Hits = tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
  HPSbyTightCombinedIsolationDeltaBetaCorr3Hits = tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
  HPSbyCombinedIsolationDeltaBetaCorrRaw3Hits = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");

  HPSdecayModeFinding = tau.tauID("decayModeFinding");
  HPSdecayModeFindingNewDMs = tau.tauID("decayModeFindingNewDMs");
  HPSdecayModeFindingOldDMs = tau.tauID("decayModeFindingOldDMs");

  if(tau.leadPFChargedHadrCand().isNonnull()){
    leadingTrackPt= tau.leadPFChargedHadrCand()->pt();
    charge= tau.leadPFChargedHadrCand()->charge();

    if(tau.leadPFChargedHadrCand()->trackRef().isNonnull()){
      leadingTrackValid= 1;
      leadingTrackIpVtdxy= tau.leadPFChargedHadrCand()->trackRef()->dxy(pv);
      leadingTrackIpVtdz= tau.leadPFChargedHadrCand()->trackRef()->dz(pv);
      leadingTrackIpVtdxyError= tau.leadPFChargedHadrCand()->trackRef()->dxyError();
      leadingTrackIpVtdzError= tau.leadPFChargedHadrCand()->trackRef()->dzError();
      leadingTrackVx= tau.leadPFChargedHadrCand()->trackRef()->vx();
      leadingTrackVy= tau.leadPFChargedHadrCand()->trackRef()->vy();
      leadingTrackVz= tau.leadPFChargedHadrCand()->trackRef()->vz();
      leadingTrackValidHits= tau.leadPFChargedHadrCand()->trackRef()->numberOfValidHits();
      leadingTrackNormChiSqrd= tau.leadPFChargedHadrCand()->trackRef()->normalizedChi2();
    }else{
      leadingTrackValid= 0;
      leadingTrackIpVtdxy= -99;
      leadingTrackIpVtdz= -99;
      leadingTrackIpVtdxyError= -99;
      leadingTrackIpVtdzError= -99;
      leadingTrackVx= -99;
      leadingTrackVy= -99;
      leadingTrackVz= -99;
      leadingTrackValidHits= -99;
      leadingTrackNormChiSqrd= -99;
    }
  }else{
    leadingTrackValid= 0;
    leadingTrackPt= -99;
    charge= -99;
    leadingTrackIpVtdxy= -99;
    leadingTrackIpVtdz= -99;
    leadingTrackIpVtdxyError= -99;
    leadingTrackIpVtdzError= -99;
    leadingTrackVx= -99;
    leadingTrackVy= -99;
    leadingTrackVz= -99;
    leadingTrackValidHits= -99;
    leadingTrackNormChiSqrd= -99;
  }
}
