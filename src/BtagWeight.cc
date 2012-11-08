
#include "../interface/BtagWeight.h"

bool BTagWeight::filter(int t)
{
  return (t >= minTags && t <= maxTags);
}



double BTagWeight::weight(vector<JetInfo> jets, int useMinTags=-1, int useMaxTags=-1 )
{

  if( useMinTags>-1 ) minTags = useMinTags;
  if( useMaxTags>-1 ) maxTags = useMaxTags;

  int njets=jets.size();
  int comb= 1 << njets;
  float pMC=0;
  float pData=0;
  for(int i=0;i < comb; i++){
    float mc=1.;
    float data=1.;
    int ntagged=0;
    for(int j=0;j<njets;j++){
      bool tagged = ((i >> j) & 0x1) == 1;
      if(tagged){
	ntagged++;
	mc*=jets[j].eff;
	data*=jets[j].eff*jets[j].sf;
      }
      else{
	mc*=(1.-jets[j].eff);
	data*=(1.-jets[j].eff*jets[j].sf);
      }
    }       
   
    if(filter(ntagged)){
      //std::cout << njets << " " << mc << " " << data << endl;
      pMC+=mc;
      pData+=data;
    }
  }

  double ratio = ( pMC>0. ) ? pData/pMC : 0;

  return ratio;
}
