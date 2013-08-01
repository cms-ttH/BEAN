#ifndef ProductArea_BNprimaryvertex_h
#define ProductArea_BNprimaryvertex_h

#include <vector>


// a simple class
struct BNprimaryvertex
{
  explicit BNprimaryvertex(double v,int c):x(v),xError(v),y(v),yError(v),z(v),zError(v),rho(v),normalizedChi2(v),ndof(v),isFake(c),isValid(c),tracksSize(c),isGood(c) { }
BNprimaryvertex():x(-99),xError(-99),y(-99),yError(-99),z(-99),zError(-99),rho(-99),normalizedChi2(-99),ndof(-99),isFake(-99),isValid(-99),tracksSize(-99),isGood(-99) { }
  double x,xError,y,yError,z,zError,rho,normalizedChi2,ndof;
  int isFake,isValid,tracksSize,isGood;

};

// this is our new product, it is simply a 
// collection of BNprimaryvertex held in an std::vector
typedef std::vector<BNprimaryvertex> BNprimaryvertexCollection;

#endif
