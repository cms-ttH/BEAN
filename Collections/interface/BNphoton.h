#ifndef ProductArea_BNphoton_h
#define ProductArea_BNphoton_h

#include <vector>


// a simple class
struct BNphoton
{
  explicit BNphoton(double v, int c):energy(v),et(v),pt(v),px(v),py(v),pz(v),phi(v),eta(v),theta(v),trackIso(v),ecalIso(v),hcalIso(v),caloIso(v),trackIsoHollowConeDR03(v),trackIsoSolidConeDR03(v),ecalIsoDR03(v),hcalIsoDR03(v),caloIsoDR03(v),trackIsoHollowConeDR04(v),trackIsoSolidConeDR04(v),ecalIsoDR04(v),hcalIsoDR04(v),caloIsoDR04(v),hadOverEm(v),sigmaEtaEta(v),sigmaIetaIeta(v),r9(v),scEnergy(v),scRawEnergy(v),scSeedEnergy(v),scEta(v),scPhi(v),scZ(v),genET(v),genPT(v),genPhi(v),genEta(v),genMotherET(v),genMotherPT(v),genMotherPhi(v),genMotherEta(v),eMax(v),eLeft(v),eRight(v),eTop(v),eBottom(v),e3x3(v),swissCross(v),seedEnergy(v),seedTime(v),swissCrossNoI85(v),swissCrossI85(v),E2overE9NoI85(v),E2overE9I85(v),IDTight(c),IDLoose(c),IDLooseEM(c),genId(c),genCharge(c),genMotherId(c),genMotherCharge(c),isEB(c),isEE(c),isGap(c),isEBEEGap(c),isEBGap(c),isEEGap(c),hasPixelSeed(c),seedRecoFlag(c) { }
  BNphoton():energy(-99),et(-99),pt(-99),px(-99),py(-99),pz(-99),phi(-99),eta(-99),theta(-99),trackIso(-99),ecalIso(-99),hcalIso(-99),caloIso(-99),trackIsoHollowConeDR03(-99),trackIsoSolidConeDR03(-99),ecalIsoDR03(-99),hcalIsoDR03(-99),caloIsoDR03(-99),trackIsoHollowConeDR04(-99),trackIsoSolidConeDR04(-99),ecalIsoDR04(-99),hcalIsoDR04(-99),caloIsoDR04(-99),hadOverEm(-99),sigmaEtaEta(-99),sigmaIetaIeta(-99),r9(-99),scEnergy(-99),scRawEnergy(-99),scSeedEnergy(-99),scEta(-99),scPhi(-99),scZ(-99),genET(-99),genPT(-99),genPhi(-99),genEta(-99),genMotherET(-99),genMotherPT(-99),genMotherPhi(-99),genMotherEta(-99),eMax(-99),eLeft(-99),eRight(-99),eTop(-99),eBottom(-99),e3x3(-99),swissCross(-99),seedEnergy(-99),seedTime(-99),swissCrossNoI85(-99),swissCrossI85(-99),E2overE9NoI85(-99),E2overE9I85(-99),IDTight(-99),IDLoose(-99),IDLooseEM(-99),genId(-99),genCharge(-99),genMotherId(-99),genMotherCharge(-99),isEB(-99),isEE(-99),isGap(-99),isEBEEGap(-99),isEBGap(-99),isEEGap(-99),hasPixelSeed(-99),seedRecoFlag(-99) { }
  double energy,et,pt,px,py,pz,phi,eta,theta,trackIso,ecalIso,hcalIso,caloIso,trackIsoHollowConeDR03,trackIsoSolidConeDR03,ecalIsoDR03,hcalIsoDR03,caloIsoDR03,trackIsoHollowConeDR04,trackIsoSolidConeDR04,ecalIsoDR04,hcalIsoDR04,caloIsoDR04,hadOverEm,sigmaEtaEta,sigmaIetaIeta,r9,scEnergy,scRawEnergy,scSeedEnergy,scEta,scPhi,scZ,genET,genPT,genPhi,genEta,genMotherET,genMotherPT,genMotherPhi,genMotherEta,eMax,eLeft,eRight,eTop,eBottom,e3x3,swissCross,seedEnergy,seedTime,swissCrossNoI85,swissCrossI85,E2overE9NoI85,E2overE9I85;
  int    IDTight,IDLoose,IDLooseEM,genId,genCharge,genMotherId,genMotherCharge,isEB,isEE,isGap,isEBEEGap,isEBGap,isEEGap,hasPixelSeed,seedRecoFlag;
};

// this is our new product, it is simply a 
// collection of BNphoton held in an std::vector
typedef std::vector<BNphoton> BNphotonCollection;

#endif
