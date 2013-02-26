#include "../interface/CSVscaleFactorContainer.h"

using namespace std;


// === Class to store BTagScaleFactor objects === //
// Constructor
CSVscaleFactorContainer::CSVscaleFactorContainer(string iEra){
	bottom_scaleFactor_errors = NULL;
	era				= iEra;

	// Error checking here
	if((era != "2011") && (era != "2012_52x") && (era != "2012_53x")){ cerr << "[ERROR]\t'era' has to be either '2011', '2012_52x' or '2012_53x''." << endl; exit(1); }
}


// Destructor
CSVscaleFactorContainer::~CSVscaleFactorContainer(){
	delete bottom_scaleFactor_errors; bottom_scaleFactor_errors = NULL;
}


// Just a "table" containing a bunch of "magic numbers"
void CSVscaleFactorContainer::SetBottomFlavorBins(vector<double>& iEtaBins, vector<double>& iPtBins){

	bottom_scaleFactor_errors = new XYmap<CSVmultiplet<double> >(iEtaBins, iPtBins);
	// For b and c (From btag_payload_b)
    for(unsigned int bin = 0; bin < bottom_scaleFactor_errors->size(); bin++) {
      //double eta      = bottom_scaleFactor_errors->GetCenterX(bin);
      double pt       = bottom_scaleFactor_errors->GetCenterY(bin);
      
	if(era == "2011" || era == "2012_52x"){
	  if (pt > 500)      bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0438259, 0.0655432, 0.0866563));
	  else if (pt > 400) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0307340, 0.0578761, 0.0777011));
	  else if (pt > 320) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0315127, 0.0541299, 0.0945890));
	  else if (pt > 260) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0264928, 0.0420284, 0.0712586));
	  else if (pt > 210) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0264240, 0.0409675, 0.0653621));
	  else if (pt > 160) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0257175, 0.0383341, 0.0685143));
	  else if (pt > 120) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0168747, 0.0256242, 0.0453252));
	  else if (pt > 100) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0161619, 0.0185857, 0.0300033));
	  else if (pt > 80)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0157396, 0.0204062, 0.0290466));
	  else if (pt > 70)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0161226, 0.0227033, 0.0261482));
	  else if (pt > 60)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0152644, 0.0219349, 0.0249618));
	  else if (pt > 50)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0139824, 0.0210867, 0.0232876));
	  else if (pt > 40)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0161816, 0.0295095, 0.0362281));
	  else if (pt > 30)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0188743, 0.0295675, 0.0364717));
	  else bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0188743, 0.0295675, 0.0364717)); //duplicate of pt > 30
	} 
	else if(era == "2012_53x"){
	  if (pt > 600)      bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0474291, 0.0717567, 0.0898199));
	  else if (pt > 500) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0481008, 0.0718173, 0.0769779));
	  else if (pt > 400) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0392831, 0.0474666, 0.0575776));
	  else if (pt > 320) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0386531, 0.0465748, 0.0580424));
	  else if (pt > 260) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0198805, 0.0248119, 0.0503182));
	  else if (pt > 210) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.019182,  0.0216242, 0.0474663));
	  else if (pt > 160) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0136017, 0.0184615, 0.0295603));
	  else if (pt > 120) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0126209, 0.0229375, 0.031102));
	  else if (pt > 100) bottom_scaleFactor_errors->Set(bin, CSVmultiplet<double>(0.0160836, 0.0240102, 0.0317642));
	  else if (pt > 80)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0168479, 0.0264232, 0.0333786));
	  else if (pt > 70)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0131145, 0.0200453, 0.024608));
	  else if (pt > 60)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0145441, 0.0208719, 0.0303327));
	  else if (pt > 50)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0141137, 0.0230073, 0.0342831));
	  else if (pt > 40)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0120027, 0.0207019, 0.0263491));
	  else if (pt > 30)  bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0126178, 0.0209663, 0.0266907));
	  else bottom_scaleFactor_errors->Set(bin,  CSVmultiplet<double>(0.0484285, 0.0554504, 0.0567059)); //calculated for 20 < pt < 30
	}
	else { assert (era == "is "+era+", should be either 2012_52x, 2012_53x, or 2011"); }
    }// End for(unsigned int bin = 0 
} // End void CSVscaleFactorContainer::SetBottomFlavorBins 

// === Bottom flavor scale factor stuff
double CSVscaleFactorContainer::GetBottomFlavorScaleFactor(double iEta, double iPt, char iWP, double iScaleBC){
	if(bottom_scaleFactor_errors == NULL){ cerr << "[ERROR]\t'bottom_scaleFactor_errors' has not been initialized." << endl; exit(1); }
	double result = 0;
    bool newBottomSF_v0 = false;
    bool newBottomSF_v1 = true;

	if(era == "2011" || era == "2012_52x"){
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
	}
	else if(era == "2012_53x"){
	  switch(iWP){
		case 'L':
          if (newBottomSF_v0) {
            if (iPt > 300)     { result = 1.01129-4.50363e-10*pow(fabs(300+47.5028),3.52651);
              result += (2 * iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetLoose()))); }
            else if (iPt < 30) { result = 1.01129-4.50363e-10*pow(fabs(30+47.5028),3.52651);
              result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetLoose()))); }
            else               { result = 1.01129-4.50363e-10*pow(fabs(iPt+47.5028),3.52651);
              result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetLoose()))); }
          }
          else if (newBottomSF_v1) {
            if (iPt > 240)     { result = 0.994401 + (240<117.72)*(-0.00820274*pow(0.0186967*(240-117.72),3)) + (240>117.72)*(-0.0339154*pow(0.0186967*(240-117.72),2));
              result += (2 * iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetLoose()))); }
            else if (iPt < 30) { result = 0.994401 + (30<117.72)*(-0.00820274*pow(0.0186967*(30-117.72),3)) + (30>117.72)*(-0.0339154*pow(0.0186967*(30-117.72),2));
              result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetLoose()))); }
            else               { result = 0.994401 + (iPt<117.72)*(-0.00820274*pow(0.0186967*(iPt-117.72),3)) + (iPt>117.72)*(-0.0339154*pow(0.0186967*(iPt-117.72),2));
              result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetLoose()))); }
          }
          else {
            if (iPt > 800) { result = 0.981149*((1.+(-0.000713295*800))/(1.+(-0.000703264*800)));
			  result += (2 * iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetLoose()))); }
			else 	       { result = 0.981149*((1.+(-0.000713295*iPt))/(1.+(-0.000703264*iPt)));
			  result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetLoose()))); }
          }
          break;
		case 'M':
          if (newBottomSF_v0) {
            if (iPt > 300)     { result = 0.971175-1.48871e-06*pow(fabs(300-36.4402),2.20671);
              result += (2 * iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetMedium()))); }
            else if (iPt < 30) { result = 0.971175-1.48871e-06*pow(fabs(30-36.4402),2.20671);
              result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetMedium()))); }
            else               { result = 0.971175-1.48871e-06*pow(fabs(iPt-36.4402),2.20671);
              result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetMedium()))); }
          }
          else if (newBottomSF_v1) {
            if (iPt > 240)     { result = 0.970772 + (240<117.72)*(-0.0184160*pow(0.0103070*(240-117.72),3)) + (240>117.72)*(-0.0864376*pow(0.0103070*(240-117.72),2));
              result += (2 * iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetMedium()))); }
            else if (iPt < 30) { result = 0.970772 + (30<117.72)*(-0.0184160*pow(0.0103070*(30-117.72),3)) + (30>117.72)*(-0.0864376*pow(0.0103070*(30-117.72),2));
              result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetMedium()))); }
            else               { result = 0.970772 + (iPt<117.72)*(-0.0184160*pow(0.0103070*(iPt-117.72),3)) + (iPt>117.72)*(-0.0864376*pow(0.0103070*(iPt-117.72),2));
              result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetMedium()))); } 
          }
          else {
            if (iPt > 800) { result = 0.726981*((1.+(0.253238*800))/(1.+(0.188389*800)));
              result += (2 * iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetMedium()))); }
            else           { result = 0.726981*((1.+(0.253238*iPt))/(1.+(0.188389*iPt)));
              result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetMedium()))); }
          }
          break;
		case 'T':
          if (newBottomSF_v0 || newBottomSF_v1) {
            if (iPt > 300)     { result = 0.968551-2.90523e-06*pow(fabs(300-1.73201),2.0);
              result += (2 * iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetTight()))); }
            else if (iPt < 30) { result = 0.968551-2.90523e-06*pow(fabs(30-1.73201),2.0); 
              result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetTight()))); }
            else               { result = 0.968551-2.90523e-06*pow(fabs(iPt-1.73201),2.0); 
              result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetTight()))); }
          }
          else {
            if (iPt > 800) { result = 0.869965*((1.+(0.0335062*800))/(1.+(0.0304598*800)));
			  result += (2 * iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetTight()))); }
			else           { result = 0.869965*((1.+(0.0335062*iPt))/(1.+(0.0304598*iPt)));
			  result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetTight()))); }
          }
          break;
		default:
			cerr << "[ERROR]\tWP '" << iWP << "' invalid, please check. " << __LINE__ << " " << __FILE__ << endl; exit(1);
			break;
	  }
	}
	else { assert (era == "is "+era+", should be either 2012_52x, 2012_53x, or 2011"); }

	return result;

}


// === Charm flavor scale factor stuff
double CSVscaleFactorContainer::GetCharmFlavorScaleFactor(double iEta, double iPt, char iWP, double iScaleBC, double iCharmFactor){
	if(bottom_scaleFactor_errors == NULL){ cerr << "[ERROR]\t'bottom_scaleFactor_errors' has not been initialized." << endl; exit(1); }
	double result = 0;

	if(era == "2011" || era == "2012_52x"){
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
	}
	else if(era == "2012_53x"){
	  double x = 0;
	  //// for pt > 800 GeV: use the SFb value at 800 GeV with twice the quoted uncertainty 
	  if(iPt > 800){ x = 800; iScaleBC *=2; }
	  else x = iPt;

	  switch(iWP){
		case 'L':
		        result = 0.981149*((1.+(-0.000713295*x))/(1.+(-0.000703264*x)));
			result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetLoose()))) * (1 + iCharmFactor); ///iPt->800?
			break;
		case 'M':
			result = 0.726981*((1.+(0.253238*x))/(1.+(0.188389*x)));
			result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetMedium()))) * (1 + iCharmFactor);
			break;
		case 'T':
			result = 0.869965*((1.+(0.0335062*x))/(1.+(0.0304598*x)));
			result += (iScaleBC * (*(bottom_scaleFactor_errors->GetSafeObject(iEta, iPt)->GetTight()))) * (1 + iCharmFactor);
			break;
		default:
			cerr << "[ERROR]\tWP '" << iWP << "' invalid, please check. " << __LINE__ << " " << __FILE__ << endl; exit(1);
			break;
	  }
	}
	else { assert (era == "is "+era+", should be either 2012_52x, 2012_53x, or 2011"); }

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
	double mean	= 1.0;
	double min	= 1.0;
	double max	= 1.0;
	double result = 0;

	if(era == "2011" || era == "2012_52x"){
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

	  double mistagScaleFactorLightFlavor = mean;
	  if( iScaleL > 0 ){		mistagScaleFactorLightFlavor = mean + iScaleL * (max-mean); }
	  else if( iScaleL < 0 ){	mistagScaleFactorLightFlavor = mean + iScaleL * (mean-min); }


	  if(era == "2012_52x"){	result = (mistagScaleFactorLightFlavor * (0.979396 + (0.000205898 * iPt) + (2.49868e-07 * pow(iPt,2)))); }
	  else if (era == "2011"){	result = mistagScaleFactorLightFlavor; }

	}
	else if(era == "2012_53x"){
	  ///// for pt > 800 GeV and 0<abs(eta)<1.5: use the SFlight value at 800 GeV with twice the quoted uncertainty 
	  ///// for pt > 700 GeV and 1.5<abs(eta): use the SFlight value at 700 GeV with twice the quoted uncertainty 
	  double x = 0.0;
	  if ( iPt > 800 && fabs(iEta) < 1.5 ) { x = 800; iScaleL *=2;}
	  else if ( iPt > 700 && fabs(iEta) > 1.5 ) { x = 700; iScaleL *=2;}
	  else x = iPt;
	  
	  if( fabs(iEta) < 0.5 ){
	    mean = ((1.04901+(0.00152181*x))+(-3.43568e-06*(x*x)))+(2.17219e-09*(x*(x*x)));
	    min  = ((0.973773+(0.00103049*x))+(-2.2277e-06*(x*x)))+(1.37208e-09*(x*(x*x)));
	    max  = ((1.12424+(0.00201136*x))+(-4.64021e-06*(x*x)))+(2.97219e-09*(x*(x*x)));
	  }else if( fabs(iEta) < 1.0 && fabs(iEta) > 0.5 ){
	    mean = ((0.991915+(0.00172552*x))+(-3.92652e-06*(x*x)))+(2.56816e-09*(x*(x*x)));
	    min  = ((0.921518+(0.00129098*x))+(-2.86488e-06*(x*x)))+(1.86022e-09*(x*(x*x)));
	    max  = ((1.06231+(0.00215815*x))+(-4.9844e-06*(x*x)))+(3.27623e-09*(x*(x*x)))  ;
	  }else if( fabs(iEta) < 1.5 && fabs(iEta) > 1.0 ){
	    mean = ((0.962127+(0.00192796*x))+(-4.53385e-06*(x*x)))+(3.0605e-09*(x*(x*x))) ;
	    min  = ((0.895419+(0.00153387*x))+(-3.48409e-06*(x*x)))+(2.30899e-09*(x*(x*x)));
	    max  = ((1.02883+(0.00231985*x))+(-5.57924e-06*(x*x)))+(3.81235e-09*(x*(x*x))) ;
	  }else if( fabs(iEta) < 2.4 && fabs(iEta) > 1.5 ){
	    mean = ((1.06121+(0.000332747*x))+(-8.81201e-07*(x*x)))+(7.43896e-10*(x*(x*x))) ;
	    min  = ((0.983607+(0.000196747*x))+(-3.98327e-07*(x*x)))+(2.95764e-10*(x*(x*x)));
	    max  = ((1.1388+(0.000468418*x))+(-1.36341e-06*(x*x)))+(1.19256e-09*(x*(x*x)))  ;
	  }

	  result = mean;
	  if( iScaleL > 0 ){ result = mean + iScaleL * (max-mean); }
	  else if( iScaleL < 0 ){	result = mean + iScaleL * (mean-min); }

	}
	else { assert (era == "is "+era+", should be either 2012_52x, 2012_53x, or 2011"); }

	return result;
}

double CSVscaleFactorContainer::GetMistagCSVmedium(double iEta, double iPt, double iScaleL){
	double mean	= 1.0;
	double min	= 1.0;
	double max	= 1.0;
	double result = 0;

	if(era == "2011" || era == "2012_52x"){
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

	  double mistagScaleFactorLightFlavor = mean;
	  if( iScaleL > 0 ){		mistagScaleFactorLightFlavor = mean + iScaleL * (max-mean); }
	  else if( iScaleL < 0 ){	mistagScaleFactorLightFlavor = mean + iScaleL * (mean-min); }
	  
	  if(era == "2012_52x"){	result = (mistagScaleFactorLightFlavor * (1.10422 + (-0.000523856 * iPt) + (1.14251e-06 * pow(iPt,2)))); }
	  else if (era == "2011"){	result = mistagScaleFactorLightFlavor; }

	}
	else if(era == "2012_53x"){
	  double x = iPt;
	  /// for iPt > 800???

	  if( fabs(iEta) < 0.8 ){
	    mean = ((1.06238+(0.00198635*x))+(-4.89082e-06*(x*x)))+(3.29312e-09*(x*(x*x))) ;
	    min  = ((0.972746+(0.00104424*x))+(-2.36081e-06*(x*x)))+(1.53438e-09*(x*(x*x)));
	    max  = ((1.15201+(0.00292575*x))+(-7.41497e-06*(x*x)))+(5.0512e-09*(x*(x*x)))  ;
	  }else if( fabs(iEta) < 1.6 && fabs(iEta) > 0.8 ){
	    mean = ((1.08048+(0.00110831*x))+(-2.96189e-06*(x*x)))+(2.16266e-09*(x*(x*x)));
	    min  = ((0.9836+(0.000649761*x))+(-1.59773e-06*(x*x)))+(1.14324e-09*(x*(x*x)));
	    max  = ((1.17735+(0.00156533*x))+(-4.32257e-06*(x*x)))+(3.18197e-09*(x*(x*x)));
	  }else if( fabs(iEta) < 2.4 && fabs(iEta) > 1.6 ){
	    mean = ((1.09145+(0.000687171*x))+(-2.45054e-06*(x*x)))+(1.7844e-09*(x*(x*x))) ;
	    min  = ((1.00616+(0.000358884*x))+(-1.23768e-06*(x*x)))+(6.86678e-10*(x*(x*x)));
	    max  = ((1.17671+(0.0010147*x))+(-3.66269e-06*(x*x)))+(2.88425e-09*(x*(x*x)))  ;
	  }
	  
	  result = mean;
	  if( iScaleL > 0 ){	result = mean + iScaleL * (max-mean); }
	  else if( iScaleL < 0 ){	result = mean + iScaleL * (mean-min); }
	  
	}
	else { assert (era == "is "+era+", should be either 2012_52x, 2012_53x, or 2011"); }
	
	return result;
}

double CSVscaleFactorContainer::GetMistagCSVtight(double iEta, double iPt, double iScaleL){
	double mean	= 1.0;
	double min	= 1.0;
	double max	= 1.0;
	double result = 0;

	if(era == "2011" || era == "2012_52x"){
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

	  if(era == "2012_52x") iScaleL *= 1.5;  /// ???
	  double mistagScaleFactorLightFlavor = mean;
	  if( iScaleL > 0 ){		mistagScaleFactorLightFlavor = mean + iScaleL*(max-mean); }
	  else if( iScaleL < 0 ){	mistagScaleFactorLightFlavor = mean + iScaleL*(mean-min); }

	  if(era == "2012_52x"){	result = (mistagScaleFactorLightFlavor * (1.19275 + (-0.00191042 * iPt) + (2.92205e-06 * pow(iPt,2)))); }
	  else if (era == "2011"){	result = mistagScaleFactorLightFlavor; }

	}
	else if(era == "2012_53x"){
	  double x = iPt;
	  /// for iPt > 800???

	  if( fabs(iEta) < 2.4 ){
	    mean = ((1.01739+(0.00283619*x))+(-7.93013e-06*(x*x)))+(5.97491e-09*(x*(x*x))) ;
	    min  = ((0.953587+(0.00124872*x))+(-3.97277e-06*(x*x)))+(3.23466e-09*(x*(x*x)));
	    max  = ((1.08119+(0.00441909*x))+(-1.18764e-05*(x*x)))+(8.71372e-09*(x*(x*x))) ;
	  }

	  // iScaleL *= 1.5;  /// ???
	  result = mean;	  
	  if( iScaleL > 0 ){	result = mean + iScaleL*(max-mean); }
	  else if( iScaleL < 0 ){	result = mean + iScaleL*(mean-min); }
	}
	else { assert (era == "is "+era+", should be either 2012_52x, 2012_53x, or 2011"); }

	return result;
}

