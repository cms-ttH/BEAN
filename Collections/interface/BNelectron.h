#ifndef ProductArea_BNelectron_h
#define ProductArea_BNelectron_h

#include <vector>
#include "BNlepton.h"


// a simple class
struct BNelectron : BNlepton {

	explicit BNelectron(double v, int c): BNlepton(v, c, 1, 0), gsfEt(v), pIn(v), pOut(v), EscOverPin(v), EseedOverPout(v), hadOverEm(v), hcalIsoDR03depth1(v), hcalIsoDR03depth2(v), trackIsoDR04(v), ecalIsoDR04(v), hcalIsoDR04(v), hcalIsoDR04depth1(v), hcalIsoDR04depth2(v), caloIsoDR04(v), fbrem(v), absInvEMinusInvPin(v), delPhiIn(v), delEtaIn(v), scEnergy(v), scRawEnergy(v), scSigmaEtaEta(v), scSigmaIEtaIEta(v), scE1x5(v), scE2x5Max(v), scE5x5(v), scEt(v), scEta(v), scPhi(v), scZ(v), mva(v), mvaTrigV0(v), mvaNonTrigV0(v), dist(v), dcot(v), convradius(v), convPointX(v), convPointY(v), convPointZ(v), eMax(v), eLeft(v), eRight(v), eTop(v), eBottom(v), e3x3(v), swissCross(v), seedEnergy(v), seedTime(v), swissCrossNoI85(v), swissCrossI85(v), E2overE9NoI85(v), E2overE9I85(v), classification(c), numClusters(c), tkNumValidHits(c), tkCharge(c), gsfCharge(c), isEB(c), isEE(c), isGap(c), isEBEEGap(c), isEBGap(c), isEEGap(c), isEcalDriven(c), isTrackerDriven(c), numberOfLostHits(c), numberOfExpectedInnerHits(c), numberOfValidPixelHits(c), numberOfValidPixelBarrelHits(c), numberOfValidPixelEndcapHits(c), isHEEP(c), isHEEPnoEt(c), seedRecoFlag(c), eidRobustHighEnergy(c), eidRobustLoose(c), eidRobustTight(c), eidLoose(c), eidTight(c), eidVeryLooseMC(c), eidLooseMC(c), eidMediumMC(c), eidTightMC(c), eidSuperTightMC(c), eidHyperTight1MC(c), eidHyperTight2MC(c), eidHyperTight3MC(c), eidHyperTight4MC(c), passConvVeto(c), isGsfCtfScPixChargeConsistent(c){ }

	BNelectron(): BNlepton(1, 0), gsfEt(-99), pIn(-99), pOut(-99), EscOverPin(-99), EseedOverPout(-99), hadOverEm(-99), hcalIsoDR03depth1(-99), hcalIsoDR03depth2(-99), trackIsoDR04(-99), ecalIsoDR04(-99), hcalIsoDR04(-99), hcalIsoDR04depth1(-99), hcalIsoDR04depth2(-99), caloIsoDR04(-99), fbrem(-99), absInvEMinusInvPin(-99), delPhiIn(-99), delEtaIn(-99), scEnergy(-99), scRawEnergy(-99), scSigmaEtaEta(-99), scSigmaIEtaIEta(-99), scE1x5(-99), scE2x5Max(-99), scE5x5(-99), scEt(-99), scEta(-99), scPhi(-99), scZ(-99), mva(-99), mvaTrigV0(-99), mvaNonTrigV0(-99), dist(-99), dcot(-99), convradius(-99), convPointX(-99), convPointY(-99), convPointZ(-99), eMax(-99), eLeft(-99), eRight(-99), eTop(-99), eBottom(-99), e3x3(-99), swissCross(-99), seedEnergy(-99), seedTime(-99), swissCrossNoI85(-99), swissCrossI85(-99), E2overE9NoI85(-99), E2overE9I85(-99), classification(-99), numClusters(-99), tkNumValidHits(-99), tkCharge(-99), gsfCharge(-99), isEB(-99), isEE(-99), isGap(-99), isEBEEGap(-99), isEBGap(-99), isEEGap(-99), isEcalDriven(-99), isTrackerDriven(-99), numberOfLostHits(-99), numberOfExpectedInnerHits(-99), numberOfValidPixelHits(-99), numberOfValidPixelBarrelHits(-99), numberOfValidPixelEndcapHits(-99), isHEEP(-99), isHEEPnoEt(-99), seedRecoFlag(-99), eidRobustHighEnergy(-99), eidRobustLoose(-99), eidRobustTight(-99), eidLoose(-99), eidTight(-99), eidVeryLooseMC(-99), eidLooseMC(-99), eidMediumMC(-99), eidTightMC(-99), eidSuperTightMC(-99), eidHyperTight1MC(-99), eidHyperTight2MC(-99), eidHyperTight3MC(-99), eidHyperTight4MC(-99), passConvVeto(-99), isGsfCtfScPixChargeConsistent(-99){ }

	double gsfEt, pIn, pOut, EscOverPin, EseedOverPout, hadOverEm, hcalIsoDR03depth1, hcalIsoDR03depth2, trackIsoDR04, ecalIsoDR04, hcalIsoDR04, hcalIsoDR04depth1, hcalIsoDR04depth2, caloIsoDR04, fbrem, absInvEMinusInvPin, delPhiIn, delEtaIn, scEnergy, scRawEnergy, scSigmaEtaEta, scSigmaIEtaIEta, scE1x5, scE2x5Max, scE5x5, scEt, scEta, scPhi, scZ, mva, mvaTrigV0, mvaNonTrigV0, dist, dcot, convradius, convPointX, convPointY, convPointZ, eMax, eLeft, eRight, eTop, eBottom, e3x3, swissCross, seedEnergy, seedTime, swissCrossNoI85, swissCrossI85, E2overE9NoI85, E2overE9I85;

	int classification, genGrandMother11Id, numClusters, tkNumValidHits, tkCharge, gsfCharge, isEB, isEE, isGap, isEBEEGap, isEBGap, isEEGap, isEcalDriven, isTrackerDriven, numberOfLostHits, numberOfExpectedInnerHits, numberOfValidPixelHits, numberOfValidPixelBarrelHits, numberOfValidPixelEndcapHits, isHEEP, isHEEPnoEt, seedRecoFlag, eidRobustHighEnergy, eidRobustLoose, eidRobustTight, eidLoose, eidTight, eidVeryLooseMC, eidLooseMC, eidMediumMC, eidTightMC, eidSuperTightMC, eidHyperTight1MC, eidHyperTight2MC, eidHyperTight3MC, eidHyperTight4MC, passConvVeto, isGsfCtfScPixChargeConsistent;

};

// this is our new product, it is simply a 
// collection of BNelectron held in an std::vector
typedef std::vector<BNelectron> BNelectronCollection;

#endif
