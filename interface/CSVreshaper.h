#ifndef _CSVreshaper_h
#define _CSVreshaper_h

#include <iostream>
#include <vector>
#include <map>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <algorithm>

#include "Math/Interpolator.h"
#include "CSVmultiplet.h"
#include "TH1F.h"


#endif

using namespace std;

class CSVreshaper{

	// === Functions === //
	public: 
		// Constructor(s) and destructor
		CSVreshaper(CSVmultiplet<double>&, TH1F*);
		virtual ~CSVreshaper();
		void Clean();
		double GetReshapedCSVvalue(double);

	private:
		void ThrowFatalError(const string);
		double GetNewThreshold(double, double, TH1F*);

	protected:



	// === Variables === //
	public:

	protected:
	bool isSetUp;

	private:
	ROOT::Math::Interpolator* interpolator;

}; // End of class prototype



