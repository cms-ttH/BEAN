#ifndef ProductArea_BNlepton_h
#define ProductArea_BNlepton_h

#include <vector>


// a simple class
struct BNlepton{

	explicit BNlepton(double v, int c, int e=0, int m=0): energy(v), et(v), pt(v), px(v), py(v), pz(v), phi(v), eta(v), theta(v), trackIso(v), ecalIso(v), hcalIso(v), caloIso(v), trackIsoDR03(v), ecalIsoDR03(v), hcalIsoDR03(v), caloIsoDR03(v), genET(v), genPT(v), genPhi(v), genEta(v), genMotherET(v), genMotherPT(v), genMotherPhi(v), genMotherEta(v), vx(v), vy(v), vz(v), tkNormChi2(v), tkPT(v), tkEta(v), tkPhi(v), tkDZ(v), tkDZerr(v), tkD0(v), tkD0bs(v), tkD0err(v), correctedD0(v), correctedD0Vertex(v), correctedDZ(v), particleIso(v), chargedHadronIso(v), neutralHadronIso(v), photonIso(v), puChargedHadronIso(v), chargedHadronIsoDR03(v), neutralHadronIsoDR03(v), photonIsoDR03(v), puChargedHadronIsoDR03(v), chargedHadronIsoDR04(v), neutralHadronIsoDR04(v), photonIsoDR04(v), puChargedHadronIsoDR04(v), rhoPrime(v), AEffDr03(v), AEffDr04(v), IP(v), IPError(v), isElectron(e), isMuon(m), charge(c), genId(c), genCharge(c), genNumberOfMothers(c), genMotherId(c), genMotherCharge(c), genMother0Id(c), genMother1Id(c), genGrandMother00Id(c), genGrandMother01Id(c), genGrandMother10Id(c), genGrandMother11Id(c){ }

	BNlepton(int e=0, int m=0): energy(-99), et(-99), pt(-99), px(-99), py(-99), pz(-99), phi(-99), eta(-99), theta(-99), trackIso(-99), ecalIso(-99), hcalIso(-99), caloIso(-99), trackIsoDR03(-99), ecalIsoDR03(-99), hcalIsoDR03(-99), caloIsoDR03(-99), genET(-99), genPT(-99), genPhi(-99), genEta(-99), genMotherET(-99), genMotherPT(-99), genMotherPhi(-99), genMotherEta(-99), vx(-99), vy(-99), vz(-99), tkNormChi2(-99), tkPT(-99), tkEta(-99), tkPhi(-99), tkDZ(-99), tkDZerr(-99), tkD0(-99), tkD0bs(-99), tkD0err(-99), correctedD0(-99), correctedD0Vertex(-99), correctedDZ(-99), particleIso(-99), chargedHadronIso(-99), neutralHadronIso(-99), photonIso(-99), puChargedHadronIso(-99), chargedHadronIsoDR03(-99), neutralHadronIsoDR03(-99), photonIsoDR03(-99), puChargedHadronIsoDR03(-99), chargedHadronIsoDR04(-99), neutralHadronIsoDR04(-99), photonIsoDR04(-99), puChargedHadronIsoDR04(-99), rhoPrime(-99), AEffDr03(-99), AEffDr04(-99), IP(-99), IPError(-99), isElectron(e), isMuon(m), charge(-99), genId(-99), genCharge(-99), genNumberOfMothers(-99), genMotherId(-99), genMotherCharge(-99), genMother0Id(-99), genMother1Id(-99), genGrandMother00Id(-99), genGrandMother01Id(-99), genGrandMother10Id(-99), genGrandMother11Id(-99){ }

	double energy, et, pt, px, py, pz, phi, eta, theta, trackIso, ecalIso, hcalIso, caloIso, trackIsoDR03, ecalIsoDR03, hcalIsoDR03, caloIsoDR03, genET, genPT, genPhi, genEta, genMotherET, genMotherPT, genMotherPhi, genMotherEta, vx, vy, vz, tkNormChi2, tkPT, tkEta, tkPhi, tkDZ, tkDZerr, tkD0, tkD0bs, tkD0err, correctedD0, correctedD0Vertex, correctedDZ, particleIso, chargedHadronIso, neutralHadronIso, photonIso, puChargedHadronIso, chargedHadronIsoDR03, neutralHadronIsoDR03, photonIsoDR03, puChargedHadronIsoDR03, chargedHadronIsoDR04, neutralHadronIsoDR04, photonIsoDR04, puChargedHadronIsoDR04, rhoPrime, AEffDr03, AEffDr04, IP, IPError;

	int isElectron, isMuon, charge, genId, genCharge, genNumberOfMothers, genMotherId, genMotherCharge, genMother0Id, genMother1Id, genGrandMother00Id, genGrandMother01Id, genGrandMother10Id, genGrandMother11Id;

};

#endif
