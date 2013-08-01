#ifndef CSVmultiplet_H
#define CSVmultiplet_H

#include <vector>
#include <iostream>
#include <stdlib.h>

using namespace std;


// === Class to store some arbitrary type of object for each CSV working point === //
template <typename TObject> class CSVmultiplet{
	public:
		CSVmultiplet(TObject iLoose, TObject iMedium, TObject iTight) : loose(iLoose), medium(iMedium), tight(iTight){}
		~CSVmultiplet(){}

		TObject* GetLoose(){	return &loose; }
		TObject* GetMedium(){	return &medium; }
		TObject* GetTight(){	return &tight; }

	private:
		TObject loose;
		TObject medium;
		TObject tight;
};


#endif
