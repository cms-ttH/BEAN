#ifndef BTAGRESHAPHING_H
#define BTAGRESHAPHING_H
#include <utility>
#include <math.h>

#include <TH1F.h>
#include <TFile.h>
#include <TGraph.h>
#include <vector>
#include <iostream>
#include <TSpline.h>
#include "Math/Interpolator.h"
#define MAXPOINTS 200

/*#if PROJECT_NAME == CMSSW
#include "VHbbAnalysis/VHbbDataFormats/interface/btag_payload_b.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/btag_payload_light.h"
#else*/
//#include "VHbbDataFormats/interface/btag_payload_b.h"
//#include "VHbbDataFormats/interface/btag_payload_light.h"
//#include "btag_payload_b_test.h"
#include "btag_payload_b.h"
#include "btag_payload_light.h"
//#endif

class BTagShape{
	public: 
		BTagShape();
		BTagShape(TFile*, const char*, const std::vector<std::pair<float, float> >&);
		float eval(float);

	private:
		ROOT::Math::Interpolator * m_i;
};

class EtaPtBin{
	public:
		EtaPtBin();
		EtaPtBin(float, float, float, float);
		bool contains(float, float);
		float centerEta();
		float centerPt();

		float etaMin;
		float etaMax;
		float ptMin;
		float ptMax;
};

class BinnedBTagShape{
	public:
		BinnedBTagShape();
		BinnedBTagShape(std::vector<EtaPtBin>&, std::vector< std::vector<std::pair<float, float> > >&, TFile*, const char *);

		float  eval(float, float, float);
		std::vector<BTagShape> m_shapes;
		std::vector<EtaPtBin> m_bins; 

};

class BTagShapeInterface{
	public:
		BTagShapeInterface();
		BTagShapeInterface(std::string, const char *, float, float);

		float reshape(float, float, float, int);

		beff myBeff;
		TFile * m_file; 
		BinnedBTagShape * m_b;
		BinnedBTagShape * m_c;
		BinnedBTagShape * m_l;

};


#endif
