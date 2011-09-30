#ifndef ProductArea_BNevent_h
#define ProductArea_BNevent_h

#include <vector>


// a simple class
struct BNevent
{
  explicit BNevent(double v, int c):weight(v),pthat(v),qScale(v),alphaQCD(v),alphaQED(v),scalePDF(v),x1(v),x2(v),xPDF1(v),xPDF2(v),BSx(v),BSy(v),BSz(v),bField(v),run(c),evt(c),lumi(c),sample(c),numPV(c),W0decay(c),W1decay(c),Z0decay(c),Z1decay(c),hcalnoiseLoose(c),hcalnoiseTight(c) { }
  BNevent():weight(-99),pthat(-99),qScale(-99),alphaQCD(-99),alphaQED(-99),scalePDF(-99),x1(-99),x2(-99),xPDF1(-99),xPDF2(-99),BSx(-99),BSy(-99),BSz(-99),bField(-99),run(-99),evt(-99),lumi(-99),sample(-99),numPV(-99),W0decay(-99),W1decay(-99),Z0decay(-99),Z1decay(-99),hcalnoiseLoose(-99),hcalnoiseTight(-99) { }
  double weight,pthat,qScale,alphaQCD,alphaQED,scalePDF,x1,x2,xPDF1,xPDF2,BSx,BSy,BSz,bField;
  int    run,evt,lumi,sample,numPV,W0decay,W1decay,Z0decay,Z1decay,hcalnoiseLoose,hcalnoiseTight;
};

// this is our new product, it is simply a 
// collection of BNevent held in an std::vector
typedef std::vector<BNevent> BNeventCollection;

#endif
