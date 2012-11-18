#include "../interface/CSVscaleFactorContainer.h"

using namespace std;


// === Class to store BTagScaleFactor objects === //
// Constructor
CSVscaleFactorContainer::CSVscaleFactorContainer(){
	bottom_scaleFactor_errors = NULL;
}


// Destructor
CSVscaleFactorContainer::~CSVscaleFactorContainer(){
	delete bottom_scaleFactor_errors; bottom_scaleFactor_errors = NULL;
}


// Just a "table" containing a bunch of "magic numbers"
void CSVscaleFactorContainer::SetBottomFlavorBins(vector<double>& iEtaBins, vector<double>& iPtBins){

	bottom_scaleFactor_errors = new XYmap<CSVmultiplet<double> >(iEtaBins, iPtBins);

	// For b and c (From btag_payload_b)
	bottom_scaleFactor_errors->Set(0, 30,  CSVmultiplet<double>(0.0188743, 0.0295675, 0.0364717));
	bottom_scaleFactor_errors->Set(0, 40,  CSVmultiplet<double>(0.0161816, 0.0295095, 0.0362281));
	bottom_scaleFactor_errors->Set(0, 50,  CSVmultiplet<double>(0.0139824, 0.0210867, 0.0232876));
	bottom_scaleFactor_errors->Set(0, 60,  CSVmultiplet<double>(0.0152644, 0.0219349, 0.0249618));
	bottom_scaleFactor_errors->Set(0, 70,  CSVmultiplet<double>(0.0161226, 0.0227033, 0.0261482));
	bottom_scaleFactor_errors->Set(0, 80,  CSVmultiplet<double>(0.0157396, 0.0204062, 0.0290466));
	bottom_scaleFactor_errors->Set(0, 100, CSVmultiplet<double>(0.0161619, 0.0185857, 0.0300033));
	bottom_scaleFactor_errors->Set(0, 120, CSVmultiplet<double>(0.0168747, 0.0256242, 0.0453252));
	bottom_scaleFactor_errors->Set(0, 160, CSVmultiplet<double>(0.0257175, 0.0383341, 0.0685143));
	bottom_scaleFactor_errors->Set(0, 210, CSVmultiplet<double>(0.0264240, 0.0409675, 0.0653621));
	bottom_scaleFactor_errors->Set(0, 260, CSVmultiplet<double>(0.0264928, 0.0420284, 0.0712586));
	bottom_scaleFactor_errors->Set(0, 320, CSVmultiplet<double>(0.0315127, 0.0541299, 0.0945890));
	bottom_scaleFactor_errors->Set(0, 400, CSVmultiplet<double>(0.0307340, 0.0578761, 0.0777011));
	bottom_scaleFactor_errors->Set(0, 500, CSVmultiplet<double>(0.0438259, 0.0655432, 0.0866563));
}


// === Bottom flavor scale factor stuff
double CSVscaleFactorContainer::GetBottomFlavorScaleFactor(double iEta, double iPt, char iWP, double iScaleBC){
	if(bottom_scaleFactor_errors == NULL){ cerr << "[ERROR]\t'bottom_scaleFactor_errors' has not been initialized." << endl; exit(1); }
	double result = 0;

	switch(iWP){
		case 'L':
			if(iPt > 670){	result = 0.963552; }
			else{ 			result = 1.02658*((1.+(0.0195388*iPt))/(1.+(0.0209145*iPt))); }
			result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetLoose())));
			break;
		case 'M':
			if(iPt > 670){	result = 0.961716; }
			else{			result = 0.6981*((1.+(0.414063*iPt))/(1.+(0.300155*iPt))); }
			result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetMedium())));
			break;
		case 'T':
			if(iPt > 670){	result = 0.910543; }
			else{			result = 0.901615*((1.+(0.552628*iPt))/(1.+(0.547195*iPt))); }
			result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetTight())));
			break;
		default:
			cerr << "[ERROR]\tWP '" << iWP << "' invalid, please check. " << __LINE__ << " " << __FILE__ << endl; exit(1);
			break;
	}

	return result;

}


// === Charm flavor scale factor stuff
double CSVscaleFactorContainer::GetCharmFlavorScaleFactor(double iEta, double iPt, char iWP, double iScaleBC, double iCharmFactor){
	if(bottom_scaleFactor_errors == NULL){ cerr << "[ERROR]\t'bottom_scaleFactor_errors' has not been initialized." << endl; exit(1); }
	double result = 0;

	switch(iWP){
		case 'L':
			if(iPt > 670){	result = 0.963552; }
			else{ 			result = 1.02658*((1.+(0.0195388*iPt))/(1.+(0.0209145*iPt))); }
			result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetLoose()))) * (1 + iCharmFactor);
			break;
		case 'M':
			if(iPt > 670){	result = 0.961716; }
			else{			result = 0.6981*((1.+(0.414063*iPt))/(1.+(0.300155*iPt))); }
			result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetMedium()))) * (1 + iCharmFactor);
			break;
		case 'T':
			if(iPt > 670){	result = 0.910543; }
			else{			result = 0.901615*((1.+(0.552628*iPt))/(1.+(0.547195*iPt))); }
			result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetTight()))) * (1 + iCharmFactor);
			break;
		default:
			cerr << "[ERROR]\tWP '" << iWP << "' invalid, please check. " << __LINE__ << " " << __FILE__ << endl; exit(1);
			break;
	}

	return result;

}


// === Light flavor scale factor stuff === //
double CSVscaleFactorContainer::GetLightFlavorScaleFactor(double iEta, double iPt, char iWP, double iScaleL){
	switch(iWP){
		case 'L': return GetMistagCSVloose(iEta, iPt, iScaleL);	break;
		case 'M': return GetMistagCSVmedium(iEta, iPt, iScaleL);	break;
		case 'T': return GetMistagCSVtight(iEta, iPt, iScaleL);	break;
	}
	cerr << "[ERROR]\tWP '" << iWP << "' invalid, please check. " << __LINE__ << " " << __FILE__ << endl; exit(1);
	return -999;
}

double CSVscaleFactorContainer::GetMistagCSVloose(double iEta, double iPt, double iScaleL){
	float mean	= 1.0;
	float min	= 1.0;
	float max	= 1.0;

	if( fabs(iEta) < 0.5 ){
		min		= 0.994425 + (-8.66392e-05 * iPt) + (-3.03813e-08 * pow(iPt,2)) + (-3.52151e-10 * pow(iPt,3));
		max		= 1.156280 + (0.0004376680 * iPt) + (-1.69625e-06 * pow(iPt,2)) + (1.00718e-09  * pow(iPt,3));
		mean	= 1.075360 + (0.0001755060 * iPt) + (-8.63317e-07 * pow(iPt,2)) + (3.27516e-10  * pow(iPt,3));
	}else if( fabs(iEta) < 1.0 && fabs(iEta) > 0.5 ){
		min		= 0.998088 + (6.94916e-05  * iPt) + (-4.82731e-07 * pow(iPt,2)) + (1.63506e-10  * pow(iPt,3));
		max		= 1.158820 + (0.000579711  * iPt) + (-2.12243e-06 * pow(iPt,2)) + (1.53771e-09  * pow(iPt,3));
		mean	= 1.078460 + (0.000324580  * iPt) + (-1.30258e-06 * pow(iPt,2)) + (8.50608e-10  * pow(iPt,3));
	}else if( fabs(iEta) < 1.5 && fabs(iEta) > 1.0 ){
		min		= 1.002940 + (0.000289844  * iPt) + (-7.9845e-07  * pow(iPt,2)) + (5.38525e-10  * pow(iPt,3));
		max		= 1.162920 + (0.000659848  * iPt) + (-2.07868e-06 * pow(iPt,2)) + (1.72763e-09  * pow(iPt,3));
		mean	= 1.082940 + (0.000474818  * iPt) + (-1.43857e-06 * pow(iPt,2)) + (1.13308e-09  * pow(iPt,3));
	}else if( fabs(iEta) < 2.4 && fabs(iEta) > 1.5 ){
		min		= 0.979816 + (0.000138797  * iPt) + (-3.14503e-07 * pow(iPt,2)) + (2.38124e-10  * pow(iPt,3));
		max		= 1.143570 + (0.000208540  * iPt) + (-7.43519e-07 * pow(iPt,2)) + (8.73742e-10  * pow(iPt,3));
		mean	= 1.061700 + (0.000173654  * iPt) + (-5.29009e-07 * pow(iPt,2)) + (5.55931e-10  * pow(iPt,3));
	}

	if( iPt > 670 ){
		min		= 0.925136;
		max		= 1.09972;
		mean	= 1.0124;
	}

	float SFl = mean;
	if( iScaleL > 0 ){		SFl = mean + iScaleL * (max-mean); }
	else if( iScaleL < 0 ){	SFl = mean + iScaleL * (mean-min); }

	return (SFl * (0.979396 + (0.000205898 * iPt) + (2.49868e-07 * pow(iPt,2))));
}

double CSVscaleFactorContainer::GetMistagCSVmedium(double iEta, double iPt, double iScaleL){
	float mean	= 1.0;
	float min	= 1.0;
	float max	= 1.0;

	if( fabs(iEta) < 0.8 ){
		min		= 0.972455 + (7.51396e-06  * iPt) + (4.91857e-07  * pow(iPt,2)) + (-1.47661e-09 * pow(iPt,3));
		max		= 1.151160 + (0.001226570  * iPt) + (-3.63826e-06 * pow(iPt,2)) + (2.08242e-09  * pow(iPt,3));
		mean	= 1.061820 + (0.000617034  * iPt) + (-1.5732e-06  * pow(iPt,2)) + (3.02909e-10  * pow(iPt,3));
	}else if( fabs(iEta) < 1.6 && fabs(iEta) > 0.8 ){
		min		= 1.02055  + (-0.000378856 * iPt) + (1.49029e-06  * pow(iPt,2)) + (-1.74966e-09 * pow(iPt,3));
		max		= 1.20146  + (0.0003595430 * iPt) + (-1.12866e-06 * pow(iPt,2)) + (6.59918e-10  * pow(iPt,3));
		mean	= 1.11100  + (-9.64191e-06 * iPt) + (1.80811e-07  * pow(iPt,2)) + (-5.44868e-10 * pow(iPt,3));
	}else if( fabs(iEta) < 2.4 && fabs(iEta) > 1.6 ){
		min		= 0.983476 + (-0.000607242 * iPt) + (3.17997e-06  * pow(iPt,2)) + (-4.01242e-09 * pow(iPt,3));
		max		= 1.186540 + (-0.000795808 * iPt) + (3.69226e-06  * pow(iPt,2)) + (-4.22347e-09 * pow(iPt,3));
		mean	= 1.084980 + (-0.000701422 * iPt) + (3.43612e-06  * pow(iPt,2)) + (-4.11794e-09 * pow(iPt,3));
	}

	if( iPt > 670 ){
		min		= 0.844346;
		max		= 1.05012;
		mean	= 0.947232;
	}

	float SFl = mean;
	if( iScaleL > 0 ){		SFl = mean + iScaleL * (max-mean); }
	else if( iScaleL < 0 ){	SFl = mean + iScaleL * (mean-min); }

	return (SFl * (1.10422 + (-0.000523856 * iPt) + (1.14251e-06 * pow(iPt,2))));
}

double CSVscaleFactorContainer::GetMistagCSVtight(double iEta, double iPt, double iScaleL){
	float mean	= 1.0;
	float min	= 1.0;
	float max	= 1.0;

	if( fabs(iEta) < 2.4 ){
		min		= 0.899715 + (0.00102278 * iPt) + (-2.46335e-06 * pow(iPt,2)) + (9.71143e-10 * pow(iPt,3));
		max		= 0.997077 + (0.00473953 * iPt) + (-1.34985e-05 * pow(iPt,2)) + (1.0032e-08  * pow(iPt,3));
		mean	= 0.948463 + (0.00288102 * iPt) + (-7.98091e-06 * pow(iPt,2)) + (5.50157e-09 * pow(iPt,3));
	}

	if( iPt > 670 ){
		min		= 0.771264;
		max		= 1.13034;
		mean	= 0.950785;
	}

	float SFl = mean;
	if( iScaleL > 0 ){		SFl = mean + 1.5*iScaleL*(max-mean); }
	else if( iScaleL < 0 ){	SFl = mean + 1.5*iScaleL*(mean-min); }

	return (SFl * (1.19275 + (-0.00191042 * iPt) + (2.92205e-06 * pow(iPt,2))));
}

