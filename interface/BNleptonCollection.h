#ifndef ProductArea_BNleptonCollection_h
#define ProductArea_BNleptonCollection_h

#include <vector>
#include "BNlepton.h"
#include "BNelectron.h"
#include "BNmuon.h"
#include <iostream>


class BNleptonCollection : public std::vector<BNlepton*>{

	public:

		~BNleptonCollection();

		void push_back(BNelectron*);
		void push_back(BNmuon*);
		void push_back(BNelectron const &);
		void push_back(BNmuon const &);
		void push_back(BNelectronCollection const &);
		void push_back(BNmuonCollection const &);

};

#endif
