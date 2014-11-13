#ifndef ProductArea_BNleptonCollection_h
#define ProductArea_BNleptonCollection_h

#include <iostream>
#include <vector>
#include <algorithm>
#include "BNlepton.h"
#include "BNelectron.h"
#include "BNmuon.h"


class BNleptonCollection : public std::vector<BNlepton*>{

	public:

		~BNleptonCollection();
        BNleptonCollection();
        BNleptonCollection(const BNleptonCollection& collection);

        BNleptonCollection& operator=(const BNleptonCollection& collection);

		void clear();
		void sort();
		void push_back(BNelectron*);
		void push_back(BNmuon*);
		void push_back(BNlepton*);
		void push_back(BNelectron const &);
		void push_back(BNmuon const &);
		void push_back(BNelectronCollection const &);
		void push_back(BNmuonCollection const &);
		void push_back(BNleptonCollection const &);

	private:
		bool ptComparison(BNlepton const &, BNlepton const &);


};

#endif
