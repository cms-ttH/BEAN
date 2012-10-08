#ifndef ProductArea_BNtrack_h
#define ProductArea_BNtrack_h

#include <vector>


// a simple class
struct BNtrack
{
  explicit BNtrack(double v, int c):pt(v),px(v),py(v),pz(v),phi(v),eta(v),theta(v),normChi2(v),dZ(v),d0(v),d0err(v),vx(v),vy(v),vz(v),charge(c),numValidHits(c) { }
  BNtrack():pt(-99),px(-99),py(-99),pz(-99),phi(-99),eta(-99),theta(-99),normChi2(-99),dZ(-99),d0(-99),d0err(-99),vx(-99),vy(-99),vz(-99),charge(-99),numValidHits(-99) { }
  double pt,px,py,pz,phi,eta,theta,normChi2,dZ,d0,d0err,vx,vy,vz;
  int    charge,numValidHits;

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
