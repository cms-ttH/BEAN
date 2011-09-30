#ifndef ProductArea_BNskimbit_h
#define ProductArea_BNskimbit_h

#include <vector>
#include <string>

// a simple class
struct BNskimbit
{
  explicit BNskimbit(int c):EG(c),JetMETTau(c),HLT_Ele15_LW_L1R(c),HLT_Jet15U(c),HLT_Jet30U(c),HLT_Jet50U(c),HLT_Ele15_SW_L1R(c),HLT_L1Jet15(c),HLT_Jet30(c),HLT_Jet50(c),Ncalojet(c),Npfjet(c),Ntcjet(c),Nele(c),MET30(c),MET100(c),MET150(c),L1T_TechBit_032(c),L1T_TechBit_033(c),L1T_TechBit_040(c),L1T_TechBit_041(c),L1T_TechBit_032_to_043(c),GoodVertex(c),FilterOutScraping(c),FilterOutScrapingFraction(c) { }
BNskimbit():EG(-99),JetMETTau(-99),HLT_Ele15_LW_L1R(-99),HLT_Jet15U(-99),HLT_Jet30U(-99),HLT_Jet50U(-99),HLT_Ele15_SW_L1R(-99),HLT_L1Jet15(-99),HLT_Jet30(-99),HLT_Jet50(-99),Ncalojet(-99),Npfjet(-99),Ntcjet(-99),Nele(-99),MET30(-99),MET100(-99),MET150(-99),L1T_TechBit_032(-99),L1T_TechBit_033(-99),L1T_TechBit_040(-99),L1T_TechBit_041(-99),L1T_TechBit_032_to_043(-99),GoodVertex(-99),FilterOutScraping(-99),FilterOutScrapingFraction(-99) { }
  int    EG,JetMETTau,HLT_Ele15_LW_L1R,HLT_Jet15U,HLT_Jet30U,HLT_Jet50U,HLT_Ele15_SW_L1R,HLT_L1Jet15,HLT_Jet30,HLT_Jet50,Ncalojet,Npfjet,Ntcjet,Nele,MET30,MET100,MET150,L1T_TechBit_032,L1T_TechBit_033,L1T_TechBit_040,L1T_TechBit_041,L1T_TechBit_032_to_043,GoodVertex,FilterOutScraping,FilterOutScrapingFraction;
};

// this is our new product, it is simply a 
// collection of BNskimbits held in an std::vector
typedef std::vector<BNskimbit> BNskimbitsCollection;

#endif
