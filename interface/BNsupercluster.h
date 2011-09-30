#ifndef ProductArea_BNsupercluster_h
#define ProductArea_BNsupercluster_h

#include <vector>


// a simple class
struct BNsupercluster
{
  explicit BNsupercluster(double v):energy(v),et(v),ex(v),ey(v),ez(v),phi(v),eta(v),theta(v) { }
  BNsupercluster():energy(-99),et(-99),ex(-99),ey(-99),ez(-99),phi(-99),eta(-99),theta(-99) { }
  double energy,et,ex,ey,ez,phi,eta,theta;

};

// this is our new product, it is simply a 
// collection of BNsupercluster held in an std::vector
typedef std::vector<BNsupercluster> BNsuperclusterCollection;

#endif
