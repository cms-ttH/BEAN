#ifndef ProductArea_BNgentoptagjet_h
#define ProductArea_BNgentoptagjet_h

#include <vector>

#include "BEAN/Collections/interface/BNgenjet.h"

struct BNgentoptagjet
{
	explicit BNgentoptagjet(double v,  int c) :
		fatjet(BNgenjet(v,c)),
		topjet(BNgenjet(v,c)),
		nonW(BNgenjet(v,c)),
		W(BNgenjet(v,c)),
		W1(BNgenjet(v,c)),
		W2(BNgenjet(v,c)),
		toptag(false),
		subjettiness1(v),
		subjettiness2(v),
		subjettiness3(v),
		subjettiness4(v) {};
		
	BNgentoptagjet() :
		fatjet(BNgenjet()),
		topjet(BNgenjet()),
		nonW(BNgenjet()),
		W(BNgenjet()),
		W1(BNgenjet()),
		W2(BNgenjet()),
		toptag(false),
		subjettiness1(-99),
		subjettiness2(-99),
		subjettiness3(-99),
		subjettiness4(-99) {};
		
  BNgenjet fatjet;
	BNgenjet topjet;
	BNgenjet nonW;
	BNgenjet W;
	BNgenjet W1;
	BNgenjet W2;
  
	bool toptag;
	
	double subjettiness1;
	double subjettiness2;
	double subjettiness3;
	double subjettiness4;
};

typedef std::vector<BNgentoptagjet> BNgentoptagjetCollection;

#endif
