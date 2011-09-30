#ifndef ProductArea_BNjet_h
#define ProductArea_BNjet_h

#include <vector>


// a simple class
struct BNjet
{
  explicit BNjet(double v, int c):energy(v),et(v),pt(v),px(v),py(v),pz(v),phi(v),eta(v),theta(v),Upt(v),Uenergy(v),L2pt(v),L2L3pt(v),L2L3respt(v),respt(v),EMfrac(v),Hadfrac(v),charge(v),mass(v),area(v),fHPD(v),approximatefHPD(v),genPartonET(v),genPartonPT(v),genPartonEta(v),genPartonPhi(v),genJetET(v),genJetPT(v),genJetEta(v),genJetPhi(v),btagTChighPur(v),btagTChighEff(v),btagJetProb(v),btagJetBProb(v),btagSoftEle(v),btagSoftMuon(v),btagSoftMuonNoIP(v),btagSecVertex(v),btagSecVertexHighEff(v),btagSecVertexHighPur(v),n90Hits(v),hitsInN90(v),chargedHadronEnergyFraction(v),neutralHadronEnergyFraction(v),chargedEmEnergyFraction(v),neutralEmEnergyFraction(v),fLong(v),fShort(v),etaetaMoment(v),phiphiMoment(v),JESunc(v),flavour(c),Nconst(c),jetIDMinimal(c),jetIDLooseAOD(c),jetIDLoose(c),jetIDTight(c),genPartonId(c),genPartonMotherId(c),genPartonGrandMotherId(c),chargedMultiplicity(c),neutralMultiplicity(c),nconstituents(c),nHit(c) { }
BNjet():energy(-99),et(-99),pt(-99),px(-99),py(-99),pz(-99),phi(-99),eta(-99),theta(-99),Upt(-99),Uenergy(-99),L2pt(-99),L2L3pt(-99),L2L3respt(-99),respt(-99),EMfrac(-99),Hadfrac(-99),charge(-99),mass(-99),area(-99),fHPD(-99),approximatefHPD(-99),genPartonET(-99),genPartonPT(-99),genPartonEta(-99),genPartonPhi(-99),genJetET(-99),genJetPT(-99),genJetEta(-99),genJetPhi(-99),btagTChighPur(-99),btagTChighEff(-99),btagJetProb(-99),btagJetBProb(-99),btagSoftEle(-99),btagSoftMuon(-99),btagSoftMuonNoIP(-99),btagSecVertex(-99),btagSecVertexHighEff(-99),btagSecVertexHighPur(-99),n90Hits(-99),hitsInN90(-99),chargedHadronEnergyFraction(-99),neutralHadronEnergyFraction(-99),chargedEmEnergyFraction(-99),neutralEmEnergyFraction(-99),fLong(-99),fShort(-99),etaetaMoment(-99),phiphiMoment(-99),JESunc(-99),flavour(-99),Nconst(-99),jetIDMinimal(-99),jetIDLooseAOD(-99),jetIDLoose(-99),jetIDTight(-99),genPartonId(-99),genPartonMotherId(-99),genPartonGrandMotherId(-99),chargedMultiplicity(-99),neutralMultiplicity(-99),nconstituents(-99),nHit(-99) { }
  double energy,et,pt,px,py,pz,phi,eta,theta,Upt,Uenergy,L2pt,L2L3pt,L2L3respt,respt,EMfrac,Hadfrac,charge,mass,area,fHPD,approximatefHPD,genPartonET,genPartonPT,genPartonEta,genPartonPhi,genJetET,genJetPT,genJetEta,genJetPhi,btagTChighPur,btagTChighEff,btagJetProb,btagJetBProb,btagSoftEle,btagSoftMuon,btagSoftMuonNoIP,btagSecVertex,btagSecVertexHighEff,btagSecVertexHighPur,n90Hits,hitsInN90,chargedHadronEnergyFraction,neutralHadronEnergyFraction,chargedEmEnergyFraction,neutralEmEnergyFraction,fLong,fShort,etaetaMoment,phiphiMoment,JESunc;
  int    flavour,Nconst,jetIDMinimal,jetIDLooseAOD,jetIDLoose,jetIDTight,genPartonId,genPartonMotherId,genPartonGrandMotherId,chargedMultiplicity,neutralMultiplicity,nconstituents,nHit;
};

// this is our new product, it is simply a 
// collection of BNjet held in an std::vector
typedef std::vector<BNjet> BNjetCollection;

#endif
