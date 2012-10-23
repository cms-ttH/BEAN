#ifndef ProductArea_BNmet_h
#define ProductArea_BNmet_h

#include <vector>


// a simple class
struct BNmet
{
  explicit BNmet(double v):et(v),pt(v),px(v),py(v),phi(v),Upt(v),Uphi(v),NeutralEMFraction(v),NeutralHadEtFraction(v),ChargedEMEtFraction(v),ChargedHadEtFraction(v),MuonEtFraction(v),Type6EtFraction(v),Type7EtFraction(v),genPT(v),genPhi(v),muonCorEx(v),muonCorEy(v),jet20CorEx(v),jet20CorEy(v),jet1CorEx(v),jet1CorEy(v),sumET(v),corSumET(v),mEtSig(v),metSignificance(v),significance(v),sigmaX2(v),sigmaY2(v),sigmaXY(v),sigmaYX(v),maxEtInEmTowers(v),emEtFraction(v),emEtInEB(v),emEtInEE(v),emEtInHF(v),maxEtInHadTowers(v),hadEtFraction(v),hadEtInHB(v),hadEtInHE(v),hadEtInHF(v),hadEtInHO(v),UDeltaPx(v),UDeltaPy(v),UDeltaP(v),Uscale(v),type2corPx(v),type2corPy(v),T2pt(v),T2px(v),T2py(v),T2phi(v),T2sumET(v),pfT1jet1pt(v),pfT1jet1phi(v),pfT1jet6pt(v),pfT1jet6phi(v),pfT1jet10pt(v),pfT1jet10phi(v) { }
  BNmet():et(-99),pt(-99),px(-99),py(-99),phi(-99),Upt(-99),Uphi(-99),NeutralEMFraction(-99),NeutralHadEtFraction(-99),ChargedEMEtFraction(-99),ChargedHadEtFraction(-99),MuonEtFraction(-99),Type6EtFraction(-99),Type7EtFraction(-99),genPT(-99),genPhi(-99),muonCorEx(-99),muonCorEy(-99),jet20CorEx(-99),jet20CorEy(-99),jet1CorEx(-99),jet1CorEy(-99),sumET(-99),corSumET(-99),mEtSig(-99),metSignificance(-99),significance(-99),sigmaX2(-99),sigmaY2(-99),sigmaXY(-99),sigmaYX(-99),maxEtInEmTowers(-99),emEtFraction(-99),emEtInEB(-99),emEtInEE(-99),emEtInHF(-99),maxEtInHadTowers(-99),hadEtFraction(-99),hadEtInHB(-99),hadEtInHE(-99),hadEtInHF(-99),hadEtInHO(-99),UDeltaPx(-99),UDeltaPy(-99),UDeltaP(-99),Uscale(-99),type2corPx(-99),type2corPy(-99),T2pt(-99),T2px(-99),T2py(-99),T2phi(-99),T2sumET(-99),pfT1jet1pt(-99),pfT1jet1phi(-99),pfT1jet6pt(-99),pfT1jet6phi(-99),pfT1jet10pt(-99),pfT1jet10phi(-99) { }
  double et,pt,px,py,phi,Upt,Uphi,NeutralEMFraction,NeutralHadEtFraction,ChargedEMEtFraction,ChargedHadEtFraction,MuonEtFraction,Type6EtFraction,Type7EtFraction,genPT,genPhi,muonCorEx,muonCorEy,jet20CorEx,jet20CorEy,jet1CorEx,jet1CorEy,sumET,corSumET,mEtSig,metSignificance,significance,sigmaX2,sigmaY2,sigmaXY,sigmaYX,maxEtInEmTowers,emEtFraction,emEtInEB,emEtInEE,emEtInHF,maxEtInHadTowers,hadEtFraction,hadEtInHB,hadEtInHE,hadEtInHF,hadEtInHO,UDeltaPx,UDeltaPy,UDeltaP,Uscale,type2corPx,type2corPy,T2pt,T2px,T2py,T2phi,T2sumET,pfT1jet1pt,pfT1jet1phi,pfT1jet6pt,pfT1jet6phi,pfT1jet10pt,pfT1jet10phi;

};

// this is our new product, it is simply a 
// collection of BNmet held in an std::vector
typedef std::vector<BNmet> BNmetCollection;

#endif
