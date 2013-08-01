#ifndef ProductArea_BNtrigobj_h
#define ProductArea_BNtrigobj_h

#include <vector>
#include <string>

// a simple class
struct BNtrigobj
{
  explicit BNtrigobj(double v, int c, std::string n):pt(v),eta(v),phi(v),px(v),py(v),pz(v),et(v),energy(v),etTotal(v),id(c),charge(c),isIsolated(c),isMip(c),isForward(c),isRPC(c),bx(c),filter(n) { }
  BNtrigobj():pt(-99),eta(-99),phi(-99),px(-99),py(-99),pz(-99),et(-99),energy(-99),etTotal(-99),id(-99),charge(-99),isIsolated(-99),isMip(-99),isForward(-99),isRPC(-99),bx(-99),filter("e") { }
  double pt,eta,phi,px,py,pz,et,energy,etTotal;
  int    id,charge,isIsolated,isMip,isForward,isRPC,bx;
  std::string filter;

};

// this is our new product, it is simply a 
// collection of BNtrigobj held in an std::vector
typedef std::vector<BNtrigobj> BNtrigobjCollection;

#endif
