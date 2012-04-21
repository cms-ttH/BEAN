#ifndef ProductArea_BNgenjet_h
#define ProductArea_BNgenjet_h

#include <vector>
#include <string>

// a simple class
struct BNgenjet
{
  explicit BNgenjet(double v, int c):pt(v),eta(v),phi(v),px(v),py(v),pz(v),et(v),energy(v),mass(v),emEnergy(v),hadEnergy(v),invisibleEnergy(v),auxiliaryEnergy(v),charge(c) { }
  BNgenjet():pt(-99),eta(-99),phi(-99),px(-99),py(-99),pz(-99),et(-99),energy(-99),mass(-99),emEnergy(-99),hadEnergy(-99),invisibleEnergy(-99),auxiliaryEnergy(-99),charge(-99) { }
  double pt,eta,phi,px,py,pz,et,energy,mass,emEnergy,hadEnergy,invisibleEnergy,auxiliaryEnergy;
  int    charge;

};

// this is our new product, it is simply a 
// collection of BNgenjet held in an std::vector
typedef std::vector<BNgenjet> BNgenjetCollection;

#endif
