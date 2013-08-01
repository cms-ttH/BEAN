#ifndef XYmap_H
#define XYmap_H

#include <iostream>
#include <stdlib.h>

#include "TH2F.h"
#include <map>

#define MAXARRAYSIZE 1000
#define INFINITESIMUS 0.001

using namespace std;

// === Class to store a collection of XYbins === //
template <typename TObject=int> class XYmap{

	private:
		map<unsigned int, TObject*> storageMap;
		TH2F* xySpace;
		unsigned int nBinsX, nBinsY;

	public:
		XYmap(){ nBinsX = 0; nBinsY = 0; xySpace = NULL; storageMap.clear(); }
		// Constructor
		XYmap(vector<double>& iXvalues, vector<double>& iYvalues){

			// Check that the requested number of bins isn't too large
			if(iXvalues.size() >= MAXARRAYSIZE){ cerr << "[ERROR]\tNumber of x bins exceeds limit (" << MAXARRAYSIZE << ")" << endl; exit(1); }
			if(iYvalues.size() >= MAXARRAYSIZE){ cerr << "[ERROR]\tNumber of y bins eyceeds limit (" << MAXARRAYSIZE << ")" << endl; exit(1); }

			// Convert bin vectors to arrays to initialize TH2F
			double yBinArray[MAXARRAYSIZE];
			double xBinArray[MAXARRAYSIZE];
			// X bins
			// First take care of the bins of interest
			for(unsigned int x = 0; x < iXvalues.size(); x++){				xBinArray[x] = iXvalues.at(x); }
			// Make sure the remaining ones contain monotonically increasing values (otherwise TH2F poops out)
			for(unsigned int x = iXvalues.size(); x < MAXARRAYSIZE; x++){	xBinArray[x] = iXvalues.back() + x; }
			// Y bins
			// First take care of the bins of interest
			for(unsigned int y = 0; y < iYvalues.size(); y++){ yBinArray[y] = iYvalues.at(y); }
			// Make sure the remaining ones contain monotonically increasing values (otherwise TH2F poops out)
			for(unsigned int y = iYvalues.size(); y < MAXARRAYSIZE; y++){	yBinArray[y] = iYvalues.back() + y; }

			// Make TH2F histo to keep track of things
			xySpace = new TH2F("", "", (iXvalues.size()-1), xBinArray, (iYvalues.size()-1), yBinArray);

			// Figure out the map dimensions
			nBinsX = xySpace->GetNbinsX() + 2; 
			nBinsY = xySpace->GetNbinsY() + 2; 

			// Clear storage map
			storageMap.clear();

			// NULL map
			for(unsigned int xIt = 0; xIt < nBinsX; xIt++){
				for(unsigned int yIt = 0; yIt < nBinsY; yIt++){
					storageMap[(xIt+yIt)] = NULL;
				}
			}
		}

		// Destructor
		~XYmap(){
			delete xySpace;
			for(typename map<unsigned int, TObject*>::iterator it = storageMap.begin(); it != storageMap.end(); ++it){ delete it->second; it->second = NULL; }
		}

		// Clear object map
		void Clear(){ for(typename map<unsigned int, TObject*>::iterator it = storageMap.begin(); it != storageMap.end(); ++it){ delete it->second; it->second = NULL; } }


		// Get bin number given x and y
		unsigned int GetBin(double iX, double iY){ return (xySpace->FindBin(iX,iY)); }


		// Return bin properties
		unsigned int xsize(){ return nBinsX; }
		unsigned int ysize(){ return nBinsY; }
		unsigned int size(){ return (xsize()*ysize()); }
		// Return the bin in the safe area (non under/overflow)
		unsigned int GetSafeBin(unsigned int iBin){
			unsigned int bin = iBin;
			// Check if this bin is in left edge
			if( (iBin % xsize()) == 0 ){ bin += 1; }
			// Check if this bin is in right edge
			else if( ((iBin + 1) % xsize()) == 0 ){ bin -= 1; }
			// Check if this bin is in lower edge
			if(iBin < xsize()){ bin += xsize(); }
			// Check if this bin is in upper edge
			else if(iBin >= (size() - xsize())){ bin -= xsize(); }

			return bin;
		}
		unsigned int GetSafeBin(double iX, double iY){ return GetSafeBin(GetBin(iX, iY)); }

		// By bin number
		double GetMinX(unsigned int iBin){
			// Check if this bin is in left edge
			unsigned int bin = GetSafeBin(iBin);
			double infinitesimus = 0.0;
			if(bin != iBin){ infinitesimus = INFINITESIMUS; }

			Int_t binx, biny, binz;
			xySpace->GetBinXYZ(bin, binx, biny, binz);
			return (xySpace->GetXaxis()->GetBinLowEdge(binx))-infinitesimus;
		}
		double GetMaxX(unsigned int iBin){
			// Check if this bin is in right edge
			unsigned int bin = GetSafeBin(iBin);
			double infinitesimus = 0.0;
			if(bin != iBin){ infinitesimus = INFINITESIMUS; }

			Int_t binx, biny, binz;
			xySpace->GetBinXYZ(bin, binx, biny, binz);
			return (xySpace->GetXaxis()->GetBinUpEdge(binx))+infinitesimus;
		}
		double GetMinY(unsigned int iBin){
			// Check if this bin is in lower edge
			unsigned int bin = GetSafeBin(iBin);
			double infinitesimus = 0.0;
			if(bin != iBin){ infinitesimus = INFINITESIMUS; }

			Int_t binx, biny, binz;
			xySpace->GetBinXYZ(bin, binx, biny, binz);
			return (xySpace->GetYaxis()->GetBinLowEdge(biny))-infinitesimus;
		}
		double GetMaxY(unsigned int iBin){
			// Check if this bin is in upper edge
			unsigned int bin = GetSafeBin(iBin);
			double infinitesimus = 0.0;
			if(bin != iBin){ infinitesimus = INFINITESIMUS; }

			Int_t binx, biny, binz;
			xySpace->GetBinXYZ(bin, binx, biny, binz);
			return (xySpace->GetYaxis()->GetBinUpEdge(biny))+infinitesimus;
		}
		double GetCenterX(unsigned int iBin){ return (GetMaxX(iBin) + GetMinX(iBin))/2.0; }
		double GetCenterY(unsigned int iBin){ return (GetMaxY(iBin) + GetMinY(iBin))/2.0; }
		// By coordinates
		double GetCenterX(double iX, double iY){ return GetCenterX(GetBin(iX, iY)); }
		double GetMinX(double iX, double iY){ return GetMinX(GetBin(iX, iY)); }
		double GetMaxX(double iX, double iY){ return GetMaxX(GetBin(iX, iY)); }
		double GetCenterY(double iX, double iY){ return GetCenterY(GetBin(iX, iY)); }
		double GetMinY(double iX, double iY){ return GetMinY(GetBin(iX, iY)); }
		double GetMaxY(double iX, double iY){ return GetMaxY(GetBin(iX, iY)); }



		// Find and return a pointer to the object in the iX,iY bin
		TObject* GetObject(unsigned int iBin){ 
			TObject* result = storageMap[iBin];
			if(result == NULL){ cerr << "[ERROR]\tCalled GetObject for bin '" << iBin << "(" << GetCenterX(iBin) << ", " << GetCenterY(iBin) << ") but object is NULL." << endl; exit(1); }
			//if(result == NULL){ cerr << "[WARNING]\tCalled GetObject for bin '" << iBin << "(" << GetCenterX(iBin) << ", " << GetCenterY(iBin) << ") but object is NULL." << endl; }
			return result;
		}
		TObject* GetObject(double iX, double iY){ return GetObject(GetBin(iX, iY));			}
		TObject* GetSafeObject(double iX, double iY){ return GetObject(GetSafeBin(iX, iY));	}
		TObject* GetSafeObject(unsigned int iBin){ return GetObject(GetSafeBin(iBin));		}


		// Methods to add objects 
		void Set(unsigned int iBin, TObject* iObject){
			if(storageMap[iBin] != NULL){ cerr << "[ERROR]\t Bin '" << iBin << "' (" << GetCenterX(iBin) << ", " << GetCenterY(iBin) << ") already contains an object." << endl; exit(1); }
			storageMap[iBin]	= new TObject(*iObject);
		}
		void Set(unsigned int iBin, TObject iObject){ 
			Set(iBin, &iObject);
		}
		void Set(double iX, double iY, TObject* iObject){ Set(GetBin(iX, iY), iObject); }
		void Set(double iX, double iY, TObject iObject){ Set(GetBin(iX, iY), &iObject); }


};

#endif
