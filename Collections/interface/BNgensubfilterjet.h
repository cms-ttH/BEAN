#ifndef ProductArea_BNgensubfilterjet_h
#define ProductArea_BNgensubfilterjet_h

#include <vector>

#include "BEAN/Collections/interface/BNgenjet.h"

struct BNgensubfilterjet
{
	explicit BNgensubfilterjet(double v,  int c) :
		fatjet(BNgenjet(v,c)),
		subjettiness1(v),
		subjettiness2(v),
		subjettiness3(v),
		subjettiness4(v) {};

	BNgensubfilterjet() :
		fatjet(BNgenjet()),
		subjettiness1(-99),
		subjettiness2(-99),
		subjettiness3(-99),
		subjettiness4(-99) {};
		
	BNgenjet fatjet;
	std::vector<BNgenjet> subjets;
	std::vector<BNgenjet> filterjets;

	double subjettiness1;
	double subjettiness2;
	double subjettiness3;
	double subjettiness4;
};

typedef std::vector<BNgensubfilterjet> BNgensubfilterjetCollection;

#endif
