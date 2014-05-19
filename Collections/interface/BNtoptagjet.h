#ifndef ProductArea_BNtoptagjet_h
#define ProductArea_BNtoptagjet_h

#include <vector>

#include "BEAN/Collections/interface/BNjet.h"

struct BNtoptagjet
{
	explicit BNtoptagjet(double v,  int c) :
		fatjet(BNjet(v,c)),
		topjet(BNjet(v,c)),
		nonW(BNjet(v,c)),
		W(BNjet(v,c)),
		W1(BNjet(v,c)),
		W2(BNjet(v,c)),
		toptag(false),
		subjettiness1(v),
		subjettiness2(v),
		subjettiness3(v),
		subjettiness4(v) {};
		
	BNtoptagjet() :
		fatjet(BNjet()),
		topjet(BNjet()),
		nonW(BNjet()),
		W(BNjet()),
		W1(BNjet()),
		W2(BNjet()),
		toptag(false),
		subjettiness1(-99),
		subjettiness2(-99),
		subjettiness3(-99),
		subjettiness4(-99) {};
		
  BNjet fatjet;
	BNjet topjet;
	BNjet nonW;
	BNjet W;
	BNjet W1;
	BNjet W2;
  
	bool toptag;
	
	double subjettiness1;
	double subjettiness2;
	double subjettiness3;
	double subjettiness4;
};

typedef std::vector<BNtoptagjet> BNtoptagjetCollection;

#endif
