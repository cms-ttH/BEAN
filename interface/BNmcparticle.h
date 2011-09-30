#ifndef ProductArea_BNmcparticle_h
#define ProductArea_BNmcparticle_h

#include <vector>


// a simple class
struct BNmcparticle
{
  explicit BNmcparticle(double v, int c):energy(v),et(v),pt(v),px(v),py(v),pz(v),phi(v),eta(v),theta(v),mass(v),vx(v),vy(v),vz(v),motherET(v),motherPT(v),motherPhi(v),motherEta(v),grandMotherET(v),grandMotherPT(v),grandMotherPhi(v),grandMotherEta(v),charge(c),id(c),status(c),motherId(c),grandMotherId(c),grandMotherCharge(c) { }
  BNmcparticle():energy(-99),et(-99),pt(-99),px(-99),py(-99),pz(-99),phi(-99),eta(-99),theta(-99),mass(-99),vx(-99),vy(-99),vz(-99),motherET(-99),motherPT(-99),motherPhi(-99),motherEta(-99),grandMotherET(-99),grandMotherPT(-99),grandMotherPhi(-99),grandMotherEta(-99),charge(-99),id(-99),status(-99),motherId(-99),motherCharge(-99),grandMotherId(-99),grandMotherCharge(-99) { }
  double energy,et,pt,px,py,pz,phi,eta,theta,mass,vx,vy,vz,motherET,motherPT,motherPhi,motherEta,grandMotherET,grandMotherPT,grandMotherPhi,grandMotherEta,ggrandMotherET,ggrandMotherPT,ggrandMotherPhi,ggrandMotherEta;
  int    charge,id,status,motherId,motherCharge,grandMotherId,grandMotherCharge;
};

// this is our new product, it is simply a 
// collection of BNmcparticle held in an std::vector
typedef std::vector<BNmcparticle> BNmcparticleCollection;

#endif
