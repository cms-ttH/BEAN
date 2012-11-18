#ifndef BTAGSHAPECONTAINER_H
#define BTAGSHAPECONTAINER_H

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "CSVmultiplet.h"
#include "XYmap.h"

using namespace std;

// === Class to store BTagScaleFactor objects === //
class CSVscaleFactorContainer{
	public: 
		CSVscaleFactorContainer();
		~CSVscaleFactorContainer();
		void SetBottomFlavorBins(vector<double>&, vector<double>&);
		double GetBottomFlavorScaleFactor(double, double, char, double);
		double GetCharmFlavorScaleFactor(double, double, char, double, double);
		double GetLightFlavorScaleFactor(double, double, char, double);

	private:
		XYmap<CSVmultiplet<double> >* bottom_scaleFactor_errors;

		double GetMistagCSVloose(double, double, double);
		double GetMistagCSVmedium(double, double, double);
		double GetMistagCSVtight(double, double, double);

};


#endif
