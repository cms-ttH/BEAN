#ifndef ProductArea_BNevent_h
#define ProductArea_BNevent_h

#include <vector>


// a simple class
struct BNevent
{
  explicit BNevent(double v, int c, long l):weight(v),pthat(v),qScale(v),alphaQCD(v),alphaQED(v),scalePDF(v),x1(v),x2(v),xPDF1(v),xPDF2(v),BSx(v),BSy(v),BSz(v),bField(v),instLumi(v),bxLumi(v),FilterOutScrapingFraction(v),sumNVtx(v),run(c),lumi(c),sample(c),numPV(c),W0decay(c),W1decay(c),Z0decay(c),Z1decay(c),H0decay(c),H1decay(c),hcalnoiseLoose(c),hcalnoiseTight(c),GoodVertex(c),FilterOutScraping(c),HBHENoiseFilter(c),CSCLooseHaloId(c),CSCTightHaloId(c),EcalLooseHaloId(c),EcalTightHaloId(c),HcalLooseHaloId(c),HcalTightHaloId(c),GlobalLooseHaloId(c),GlobalTightHaloId(c),LooseId(c),TightId(c),numGenPV(c),nm1(c),n0(c),np1(c),evt(l) { }
  BNevent():weight(-99),pthat(-99),qScale(-99),alphaQCD(-99),alphaQED(-99),scalePDF(-99),x1(-99),x2(-99),xPDF1(-99),xPDF2(-99),BSx(-99),BSy(-99),BSz(-99),bField(-99),instLumi(-99),bxLumi(-99),FilterOutScrapingFraction(-99),sumNVtx(-99),run(-99),lumi(-99),sample(-99),numPV(-99),W0decay(-99),W1decay(-99),Z0decay(-99),Z1decay(-99),H0decay(-99),H1decay(-99),hcalnoiseLoose(-99),hcalnoiseTight(-99),GoodVertex(-99),FilterOutScraping(-99),HBHENoiseFilter(-99),CSCLooseHaloId(-99),CSCTightHaloId(-99),EcalLooseHaloId(-99),EcalTightHaloId(-99),HcalLooseHaloId(-99),HcalTightHaloId(-99),GlobalLooseHaloId(-99),GlobalTightHaloId(-99),LooseId(-99),TightId(-99),numGenPV(-99),nm1(-99),n0(-99),np1(-99),evt(-99) { }
  double weight,pthat,qScale,alphaQCD,alphaQED,scalePDF,x1,x2,xPDF1,xPDF2,BSx,BSy,BSz,bField,instLumi,bxLumi,FilterOutScrapingFraction,sumNVtx;
  int    run,lumi,sample,numPV,W0decay,W1decay,Z0decay,Z1decay,H0decay,H1decay,hcalnoiseLoose,hcalnoiseTight,GoodVertex,FilterOutScraping,HBHENoiseFilter,CSCLooseHaloId,CSCTightHaloId,EcalLooseHaloId,EcalTightHaloId,HcalLooseHaloId,HcalTightHaloId,GlobalLooseHaloId,GlobalTightHaloId,LooseId,TightId,numGenPV,nm1,n0,np1;
  long   evt;
};

// this is our new product, it is simply a 
// collection of BNevent held in an std::vector
typedef std::vector<BNevent> BNeventCollection;

#endif
