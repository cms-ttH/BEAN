#ifndef ProductArea_BNmuon_h
#define ProductArea_BNmuon_h

#include <vector>
#include "BNlepton.h"


// a simple class
struct BNmuon : BNlepton{

	explicit BNmuon(double v, int c): BNlepton(v, c, 0, 1), trackVetoIsoDR03(v), ecalVetoIsoDR03(v), hcalVetoIsoDR03(v), caloVetoIsoDR03(v), trackIsoDR05(v), ecalIsoDR05(v), hcalIsoDR05(v), caloIsoDR05(v), trackVetoIsoDR05(v), ecalVetoIsoDR05(v), hcalVetoIsoDR05(v), caloVetoIsoDR05(v), hcalE(v), ecalE(v), samNormChi2(v), samPT(v), samEta(v), samPhi(v), samDZ(v), samDZerr(v), samD0(v), samD0bs(v), samD0err(v), comNormChi2(v), comPT(v), comEta(v), comPhi(v), comDZ(v), comDZerr(v), comD0(v), comD0bs(v), comD0err(v), normalizedChi2(v), dVzPVz(v), dB(v), ptErr(v), innerTrackNormChi2(v), pfIsoR03SumChargedHadronPt(v), pfIsoR03SumNeutralHadronEt(v), pfIsoR03SumPhotonEt(v), pfIsoR03SumPUPt(v), pfIsoR04SumChargedHadronPt(v), pfIsoR04SumNeutralHadronEt(v), pfIsoR04SumPhotonEt(v), pfIsoR04SumPUPt(v), timeAtIpInOut(v), timeAtIpInOutErr(v), timeAtIpOutIn(v), timeAtIpOutInErr(v), ecal_time(v), hcal_time(v), ecal_timeError(v), hcal_timeError(v), energy_ecal(v), energy_hcal(v), e3x3_ecal(v), e3x3_hcal(v), energyMax_ecal(v), energyMax_hcal(v), IDGMPTight(c), tkNumValidHits(c), tkCharge(c), samNumValidHits(c), samCharge(c), comNumValidHits(c), comCharge(c), isPFMuon(c), isGoodMuon_1StationTight(c), isGlobalMuon(c), isTrackerMuon(c), isStandAloneMuon(c), isGlobalMuonPromptTight(c), numberOfValidMuonHits(c), numberOfValidTrackerHits(c), numberOfLayersWithMeasurement(c), pixelLayersWithMeasurement(c), numberOfMatches(c), numberOfValidTrackerHitsInnerTrack(c), numberOfValidPixelHits(c), numberOfMatchedStations(c), time_ndof(c), innerTrackPt(c), innerTrackPtError(c) { }

	BNmuon(): BNlepton(0, 1), trackVetoIsoDR03(-99), ecalVetoIsoDR03(-99), hcalVetoIsoDR03(-99), caloVetoIsoDR03(-99), trackIsoDR05(-99), ecalIsoDR05(-99), hcalIsoDR05(-99), caloIsoDR05(-99), trackVetoIsoDR05(-99), ecalVetoIsoDR05(-99), hcalVetoIsoDR05(-99), caloVetoIsoDR05(-99), hcalE(-99), ecalE(-99), samNormChi2(-99), samPT(-99), samEta(-99), samPhi(-99), samDZ(-99), samDZerr(-99), samD0(-99), samD0bs(-99), samD0err(-99), comNormChi2(-99), comPT(-99), comEta(-99), comPhi(-99), comDZ(-99), comDZerr(-99), comD0(-99), comD0bs(-99), comD0err(-99), normalizedChi2(-99), dVzPVz(-99), dB(-99), ptErr(-99), innerTrackNormChi2(-99), pfIsoR03SumChargedHadronPt(-99), pfIsoR03SumNeutralHadronEt(-99), pfIsoR03SumPhotonEt(-99), pfIsoR03SumPUPt(-99), pfIsoR04SumChargedHadronPt(-99), pfIsoR04SumNeutralHadronEt(-99), pfIsoR04SumPhotonEt(-99), pfIsoR04SumPUPt(-99), timeAtIpInOut(-99), timeAtIpInOutErr(-99), timeAtIpOutIn(-99), timeAtIpOutInErr(-99), ecal_time(-99), hcal_time(-99), ecal_timeError(-99), hcal_timeError(-99), energy_ecal(-99), energy_hcal(-99), e3x3_ecal(-99), e3x3_hcal(-99), energyMax_ecal(-99), energyMax_hcal(-99), IDGMPTight(-99), tkNumValidHits(-99), tkCharge(-99), samNumValidHits(-99), samCharge(-99), comNumValidHits(-99), comCharge(-99), isPFMuon(-99), isGoodMuon_1StationTight(-99), isGlobalMuon(-99), isTrackerMuon(-99), isStandAloneMuon(-99), isGlobalMuonPromptTight(-99), numberOfValidMuonHits(-99), numberOfValidTrackerHits(-99), numberOfLayersWithMeasurement(-99), pixelLayersWithMeasurement(-99), numberOfMatches(-99), numberOfValidTrackerHitsInnerTrack(-99), numberOfValidPixelHits(-99), numberOfMatchedStations(-99), time_ndof(-99), innerTrackPt(-99), innerTrackPtError(-99) { }

	double trackVetoIsoDR03, ecalVetoIsoDR03, hcalVetoIsoDR03, caloVetoIsoDR03, trackIsoDR05, ecalIsoDR05, hcalIsoDR05, caloIsoDR05, trackVetoIsoDR05, ecalVetoIsoDR05, hcalVetoIsoDR05, caloVetoIsoDR05, hcalE, ecalE, samNormChi2, samPT, samEta, samPhi, samDZ, samDZerr, samD0, samD0bs, samD0err, comNormChi2, comPT, comEta, comPhi, comDZ, comDZerr, comD0, comD0bs, comD0err, isolationR03emVetoEt, isolationR03hadVetoEt, normalizedChi2, dVzPVz, dB, ptErr, innerTrackNormChi2, pfIsoR03SumChargedHadronPt, pfIsoR03SumNeutralHadronEt, pfIsoR03SumPhotonEt, pfIsoR03SumPUPt, pfIsoR04SumChargedHadronPt, pfIsoR04SumNeutralHadronEt, pfIsoR04SumPhotonEt, pfIsoR04SumPUPt, timeAtIpInOut, timeAtIpInOutErr, timeAtIpOutIn, timeAtIpOutInErr, ecal_time, hcal_time, ecal_timeError, hcal_timeError, energy_ecal, energy_hcal, e3x3_ecal, e3x3_hcal, energyMax_ecal, energyMax_hcal, innerTrackPt, innerTrackPtError;

	int IDGMPTight,tkNumValidHits,tkCharge,samNumValidHits,samCharge,comNumValidHits,comCharge,isPFMuon,isGoodMuon_1StationTight,isGlobalMuon,isTrackerMuon,isStandAloneMuon,isGlobalMuonPromptTight,numberOfValidMuonHits,numberOfValidTrackerHits,numberOfLayersWithMeasurement,pixelLayersWithMeasurement,numberOfMatches,numberOfValidTrackerHitsInnerTrack,numberOfValidPixelHits,numberOfMatchedStations,time_ndof;

};

// this is our new product, it is simply a 
// collection of BNmuon held in an std::vector
typedef std::vector<BNmuon> BNmuonCollection;

#endif
