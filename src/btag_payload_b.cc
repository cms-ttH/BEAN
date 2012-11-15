#include "../interface/btag_payload_b.h"

beff::beff(){

	ptmin.push_back(30);		ptmax.push_back(40);
	ptmin.push_back(40);        ptmax.push_back(50);
	ptmin.push_back(50);        ptmax.push_back(60);
	ptmin.push_back(60);        ptmax.push_back(70);
	ptmin.push_back(70);        ptmax.push_back(80);
	ptmin.push_back(80);        ptmax.push_back(100);
	ptmin.push_back(100);       ptmax.push_back(120);
	ptmin.push_back(120);       ptmax.push_back(160);
	ptmin.push_back(160);       ptmax.push_back(210);
	ptmin.push_back(210);       ptmax.push_back(260);
	ptmin.push_back(260);       ptmax.push_back(320);
	ptmin.push_back(320);       ptmax.push_back(400);
	ptmin.push_back(400);       ptmax.push_back(500);
	ptmin.push_back(500);       ptmax.push_back(670);

	bins=14;

	//Tagger: CSVL within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
	CSVL_SFb_error.push_back(0.0188743);
	CSVL_SFb_error.push_back(0.0161816);
	CSVL_SFb_error.push_back(0.0139824);
	CSVL_SFb_error.push_back(0.0152644);
	CSVL_SFb_error.push_back(0.0161226);
	CSVL_SFb_error.push_back(0.0157396);
	CSVL_SFb_error.push_back(0.0161619);
	CSVL_SFb_error.push_back(0.0168747);
	CSVL_SFb_error.push_back(0.0257175);
	CSVL_SFb_error.push_back(0.026424);
	CSVL_SFb_error.push_back(0.0264928);
	CSVL_SFb_error.push_back(0.0315127);
	CSVL_SFb_error.push_back(0.030734);
	CSVL_SFb_error.push_back(0.0438259);


	//Tagger: CSVM within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
	CSVM_SFb_error.push_back(0.0295675);
	CSVM_SFb_error.push_back(0.0295095);
	CSVM_SFb_error.push_back(0.0210867);
	CSVM_SFb_error.push_back(0.0219349);
	CSVM_SFb_error.push_back(0.0227033);
	CSVM_SFb_error.push_back(0.0204062);
	CSVM_SFb_error.push_back(0.0185857);
	CSVM_SFb_error.push_back(0.0256242);
	CSVM_SFb_error.push_back(0.0383341);
	CSVM_SFb_error.push_back(0.0409675);
	CSVM_SFb_error.push_back(0.0420284);
	CSVM_SFb_error.push_back(0.0541299);
	CSVM_SFb_error.push_back(0.0578761);
	CSVM_SFb_error.push_back(0.0655432);

	//Tagger: CSVT within 30 < pt < 670 GeV, abs(eta) < 2.4, x = pt
	CSVT_SFb_error.push_back(0.0364717);
	CSVT_SFb_error.push_back(0.0362281);
	CSVT_SFb_error.push_back(0.0232876);
	CSVT_SFb_error.push_back(0.0249618);
	CSVT_SFb_error.push_back(0.0261482);
	CSVT_SFb_error.push_back(0.0290466);
	CSVT_SFb_error.push_back(0.0300033);
	CSVT_SFb_error.push_back(0.0453252);
	CSVT_SFb_error.push_back(0.0685143);
	CSVT_SFb_error.push_back(0.0653621);
	CSVT_SFb_error.push_back(0.0712586);
	CSVT_SFb_error.push_back(0.094589);
	CSVT_SFb_error.push_back(0.0777011);
	CSVT_SFb_error.push_back(0.0866563);

}

beff::~beff(){}
