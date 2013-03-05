#ifndef ProductArea_BNtrack_h
#define ProductArea_BNtrack_h

#include <vector>


// a simple class
struct BNtrack
{
  explicit BNtrack(double v,int c):
  pt(v),
    ptError(v),
    px(v),
    py(v),
    pz(v),
    phi(v),
    eta(v),
    theta(v),
    normChi2(v),
    dZ(v),
    d0(v),
    d0err(v),
    vx(v),
    vy(v),
    vz(v),
    charge(c),
    numValidHits(c),
    isHighPurity(c),
    isGoodPtResolution(c),
    isMatchedPF(c),
    caloEMDeltaRp3(v),
    caloHadDeltaRp3(v),
    caloEMDeltaRp4(v),
    caloHadDeltaRp4(v),
    caloEMDeltaRp5(v),
    caloHadDeltaRp5(v),
    lastHitOuterRadius(v),
    lastHitOuterEta(v),
    lastHitOuterTheta(v),
    lastHitOuterPhi(v),
    depTrkRp3(v),
    depEcalRp3(v),
    depHcalRp3(v),
    depHoRp3(v),
    nTracksRp3(v),
    trackerVetoPtRp3(v),
    emVetoEtRp3(v),
    hadVetoEtRp3(v),     
    hoVetoEtRp3(v),
    depTrkRp5(v),
    depEcalRp5(v),
    depHcalRp5(v),
    depHoRp5(v),
    nTracksRp5(v),
    trackerVetoPtRp5(v),
    emVetoEtRp5(v),
    hadVetoEtRp5(v),     
    hoVetoEtRp5(v),
    nHitsMissingOuter(v),
    nHitsMissingInner(v),
    nHitsMissingMiddle(v) { }

BNtrack():
  pt(-99),
    ptError(-99),
    px(-99),
    py(-99),
    pz(-99),
    phi(-99),
    eta(-99),
    theta(-99),
    normChi2(-99),
    dZ(-99),
    d0(-99),
    d0err(-99),
    vx(-99),
    vy(-99),
    vz(-99),
    charge(-99),
    numValidHits(-99),
    isHighPurity(-99),
    isGoodPtResolution(-99),
    isMatchedPF(-99),
    caloEMDeltaRp3(-99),
    caloHadDeltaRp3(-99),
    caloEMDeltaRp4(-99),
    caloHadDeltaRp4(-99),
    caloEMDeltaRp5(-99),
    caloHadDeltaRp5(-99),
    lastHitOuterRadius(-99),
    lastHitOuterEta(-99),
    lastHitOuterTheta(-99),
    lastHitOuterPhi(-99),
    depTrkRp3(-99),
    depEcalRp3(-99),
    depHcalRp3(-99),
    depHoRp3(-99),
    nTracksRp3(-99),
    trackerVetoPtRp3(-99),
    emVetoEtRp3(-99),
    hadVetoEtRp3(-99),     
    hoVetoEtRp3(-99),
    depTrkRp5(-99),
    depEcalRp5(-99),
    depHcalRp5(-99),
    depHoRp5(-99),
    nTracksRp5(-99),
    trackerVetoPtRp5(-99),
    emVetoEtRp5(-99),
    hadVetoEtRp5(-99),     
    hoVetoEtRp5(-99),
    nHitsMissingOuter(-99),
    nHitsMissingInner(-99),
    nHitsMissingMiddle(-99) { }
  double pt,    
    ptError,
    px,
    py,
    pz,
    phi,
    eta,
    theta,
    normChi2,
    dZ,
    d0,
    d0err,
    vx,
    vy,
    vz;
  int charge,
    numValidHits,
    isHighPurity,
    isGoodPtResolution,
    isMatchedPF;

  //by hand calo calc
  double caloEMDeltaRp3;  
  double caloHadDeltaRp3;  
  double caloEMDeltaRp4;  
  double caloHadDeltaRp4;  
  double caloEMDeltaRp5;  
  double caloHadDeltaRp5;  

  double lastHitOuterRadius;
  double lastHitOuterEta;
  double lastHitOuterTheta;  
  double lastHitOuterPhi;  

  //calo calc from muon pog
  double depTrkRp3;
  double depEcalRp3;
  double depHcalRp3;
  double depHoRp3;
  double nTracksRp3;
  double trackerVetoPtRp3;
  double emVetoEtRp3;
  double hadVetoEtRp3;     
  double hoVetoEtRp3;

  double depTrkRp5;
  double depEcalRp5;
  double depHcalRp5;
  double depHoRp5;
  double nTracksRp5;
  double trackerVetoPtRp5;
  double emVetoEtRp5;
  double hadVetoEtRp5;     
  double hoVetoEtRp5;

  int nHitsMissingOuter;
  int nHitsMissingInner;
  int nHitsMissingMiddle;
  
  std::vector<int> subDetIdHits;
  std::vector<int> rawDetIdHits;
  std::vector<bool> isValidHits;
  std::vector<double> modulePerpHits;
  std::vector<double> moduleZHits;
  std::vector<double> modulePhiHits;

};

// this is our new product, it is simply a 
// collection of BNtrack held in an std::vector
typedef std::vector<BNtrack> BNtrackCollection;

#endif
