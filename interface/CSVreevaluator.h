#ifndef _CSVreevaluator_h
#define _CSVreevaluator_h

#include <iostream>
#include <vector>
#include <map>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <algorithm>

#include "TFile.h"

#include "XYmap.h"
#include "CSVscaleFactorContainer.h"
#include "CSVreshaper.h"
#include "boost/assign/std/vector.hpp"

#endif

using namespace std;

class CSVreevaluator{

	// === Functions === //
	public: 
		// Constructor(s) and destructor
		CSVreevaluator(string, unsigned int, double, double, double);
		virtual ~CSVreevaluator();

		double GetReshapedCSVvalue(float, float, float, int);
		
	private:
		void ThrowFatalError(const string);

	protected:



	// === Variables === //
	public:

	protected:

	private:
		TFile* efficiencyFile;
		XYmap<CSVreshaper>* bottomFlavorReshapers;
		XYmap<CSVreshaper>* charmFlavorReshapers;
		XYmap<CSVreshaper>* lightFlavorReshapers;

}; // End of class prototype



