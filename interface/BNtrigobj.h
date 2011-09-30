#ifndef ProductArea_BNtrigobj_h
#define ProductArea_BNtrigobj_h

#include <vector>
#include <string>

// a simple class
struct BNtrigobj
{
  explicit BNtrigobj(double v, std::string n):pt(v),eta(v),phi(v),filter(n) { }
  BNtrigobj():pt(-99),eta(-99),phi(-99),filter("empty") { }
  double pt,eta,phi;
  std::string filter;

};

// this is our new product, it is simply a 
// collection of BNtrigobj held in an std::vector
typedef std::vector<BNtrigobj> BNtrigobjCollection;

#endif
