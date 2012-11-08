#ifndef INC_BTAGWEIGHT
#define INC_BTAGWEIGHT

#include <math.h>
#include <iostream>
#include <vector>
using namespace std; 
class BTagWeight {
 public:
  struct JetInfo {
  JetInfo(float mceff,float datasf) : eff(mceff), sf(datasf) {}
    float eff;
    float sf;
  };

 BTagWeight(int jmin, int jmax) : 
  maxTags(jmax), minTags(jmin) {}

  bool filter(int t);
  double weight(vector<JetInfo> jets, int useMinTags, int useMaxTags);
 private:
  int maxTags;
  int minTags;
 

};



#endif // INC_BTAGWEIGHT
