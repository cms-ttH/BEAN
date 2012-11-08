#pragma once

#ifndef BTAG_PAYLOAD_B_H
#define BTAG_PAYLOAD_B_H

#include <cstring>
#include <vector>

using namespace std;

class beff{

	public:
		beff();
		~beff();

		vector<float> ptmin;
		vector<float> ptmax;

		size_t bins;

		//Tagger: CSVL within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
		inline float CSVL_SFb(float x) { if(x > 670) return 0.963552 ; return 1.02658*((1.+(0.0195388*x))/(1.+(0.0209145*x))); }
		vector<float> CSVL_SFb_error;

		//Tagger: CSVM within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
		inline float CSVM_SFb(float x) {  if(x > 670) return 0.961716 ; return 0.6981*((1.+(0.414063*x))/(1.+(0.300155*x)));}
		vector<float> CSVM_SFb_error;

		//Tagger: CSVT within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
		inline float CSVT_SFb(float x) {  if(x > 670) return 0.910543 ; return 0.901615*((1.+(0.552628*x))/(1.+(0.547195*x))); }
		vector<float> CSVT_SFb_error;

};

#endif
