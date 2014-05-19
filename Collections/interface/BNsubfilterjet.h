#ifndef ProductArea_BNsubfilterjet_h
#define ProductArea_BNsubfilterjet_h

#include <vector>

#include "BEAN/Collections/interface/BNjet.h"

struct BNsubfilterjet
{
	explicit BNsubfilterjet(double v,  int c) :
		fatjet(BNjet(v,c)),
		subjettiness1(v),
		subjettiness2(v),
		subjettiness3(v),
		subjettiness4(v) {};
		
	BNsubfilterjet() :
		fatjet(BNjet()),
		subjettiness1(-99),
		subjettiness2(-99),
		subjettiness3(-99),
		subjettiness4(-99) {};
		
  BNjet fatjet;
	std::vector<BNjet> subjets;
	std::vector<BNjet> filterjets;

	double subjettiness1;
	double subjettiness2;
	double subjettiness3;
	double subjettiness4;
};

typedef std::vector<BNsubfilterjet> BNsubfilterjetCollection;

#endif
