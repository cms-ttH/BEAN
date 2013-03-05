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
    isGoodPtResolution;

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
