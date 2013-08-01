#ifndef ProductArea_BNbxlumi_h
#define ProductArea_BNbxlumi_h

#include <vector>

// a simple class
struct BNbxlumi
{
  explicit BNbxlumi(double v):bx_B1_now(v),bx_B2_now(v),bx_LUMI_now(v) { }
  BNbxlumi():bx_B1_now(-99),bx_B2_now(-99),bx_LUMI_now(-99) { }
  double bx_B1_now,bx_B2_now,bx_LUMI_now;

};

// this is our new product, it is simply a 
// collection of BNbxlumi held in an std::vector
typedef std::vector<BNbxlumi> BNbxlumiCollection;

#endif
