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
  HPSagainstElectronLooseMVA2 = tau.tauID("againstElectronLooseMVA2");
  HPSagainstElectronLooseMVA3 = tau.tauID("againstElectronLooseMVA3");
  HPSagainstElectronMVA = tau.tauID("againstElectronMVA");
  HPSagainstElectronMVA2category = tau.tauID("againstElectronMVA2category");
  HPSagainstElectronMVA2raw = tau.tauID("againstElectronMVA2raw");
  HPSagainstElectronMVA3category = tau.tauID("againstElectronMVA3category");
  HPSagainstElectronMVA3raw = tau.tauID("againstElectronMVA3raw");
  HPSagainstElectronMedium = tau.tauID("againstElectronMedium");
  HPSagainstElectronMediumMVA2 = tau.tauID("againstElectronMediumMVA2");
  HPSagainstElectronMediumMVA3 = tau.tauID("againstElectronMediumMVA3");
  HPSagainstElectronTight = tau.tauID("againstElectronTight");
  HPSagainstElectronTightMVA2 = tau.tauID("againstElectronTightMVA2");
  HPSagainstElectronTightMVA3 = tau.tauID("againstElectronTightMVA3");
  HPSagainstElectronVLooseMVA2 = tau.tauID("againstElectronVLooseMVA2");
  HPSagainstElectronVTightMVA3 = tau.tauID("againstElectronVTightMVA3");

  HPSagainstMuonLoose = tau.tauID("againstMuonLoose");
  HPSagainstMuonMedium = tau.tauID("againstMuonMedium");
  HPSagainstMuonTight = tau.tauID("againstMuonTight");
  HPSagainstMuonLoose2 = tau.tauID("againstMuonLoose2");
  HPSagainstMuonMedium2 = tau.tauID("againstMuonMedium2");
  HPSagainstMuonTight2 = tau.tauID("againstMuonTight2");
  HPSagainstMuonLoose3 = tau.tauID("againstMuonLoose3");
  HPSagainstMuonMedium3 = tau.tauID("againstMuonMedium3");
  HPSagainstMuonTight3 = tau.tauID("againstMuonTight3");

  HPSbyCombinedIsolationDeltaBetaCorrRaw = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw");
  HPSbyCombinedIsolationDeltaBetaCorrRaw3Hits = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
  HPSbyIsolationMVA2raw = tau.tauID("byIsolationMVA2raw");
  HPSbyLooseCombinedIsolationDeltaBetaCorr = tau.tauID("byLooseCombinedIsolationDeltaBetaCorr");
  HPSbyLooseCombinedIsolationDeltaBetaCorr3Hits = tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
  HPSbyLooseIsolationMVA = tau.tauID("byLooseIsolationMVA");
  HPSbyLooseIsolationMVA2 = tau.tauID("byLooseIsolationMVA2");
  HPSbyMediumCombinedIsolationDeltaBetaCorr = tau.tauID("byMediumCombinedIsolationDeltaBetaCorr");
  HPSbyMediumCombinedIsolationDeltaBetaCorr3Hits = tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
  HPSbyMediumIsolationMVA = tau.tauID("byMediumIsolationMVA");
  HPSbyMediumIsolationMVA2 = tau.tauID("byMediumIsolationMVA2");
  HPSbyTightCombinedIsolationDeltaBetaCorr = tau.tauID("byTightCombinedIsolationDeltaBetaCorr");
  HPSbyTightCombinedIsolationDeltaBetaCorr3Hits = tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
  HPSbyTightIsolationMVA = tau.tauID("byTightIsolationMVA");
  HPSbyTightIsolationMVA2 = tau.tauID("byTightIsolationMVA2");
  HPSbyVLooseCombinedIsolationDeltaBetaCorr = tau.tauID("byVLooseCombinedIsolationDeltaBetaCorr");
  HPSdecayModeFinding = tau.tauID("decayModeFinding");
  HPSbyIsolationMVAraw = tau.tauID("byIsolationMVAraw");

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
