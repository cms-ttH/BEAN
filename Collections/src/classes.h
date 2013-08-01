#include "DataFormats/Common/interface/Wrapper.h"

//Add includes for your classes here
#include "BEAN/Collections/interface/SampleProd.h"
#include "BEAN/Collections/interface/BNbxlumi.h"
#include "BEAN/Collections/interface/BNelectron.h"
#include "BEAN/Collections/interface/BNtau.h"
#include "BEAN/Collections/interface/BNevent.h"
#include "BEAN/Collections/interface/BNjet.h"
#include "BEAN/Collections/interface/BNmcparticle.h"
#include "BEAN/Collections/interface/BNmet.h"
#include "BEAN/Collections/interface/BNmuon.h"
#include "BEAN/Collections/interface/BNlepton.h"
#include "BEAN/Collections/interface/BNphoton.h"
#include "BEAN/Collections/interface/BNsupercluster.h"
#include "BEAN/Collections/interface/BNtrack.h"
#include "BEAN/Collections/interface/BNtrigger.h"
#include "BEAN/Collections/interface/BNskimbits.h"
#include "BEAN/Collections/interface/BNtrigobj.h"
#include "BEAN/Collections/interface/BNprimaryvertex.h"
#include "BEAN/Collections/interface/BNgenjet.h"
#include <vector>

namespace {
   struct BEAN_Collections {
      //add 'dummy' Wrapper variable for each class type you put into the Event
     SampleProd dummy0;
     edm::Wrapper<SampleProd> dummy1;
     std::vector<SampleProd> dummy2;
     edm::Wrapper<std::vector<SampleProd> > dummy3;

     BNbxlumi bxlumidummy0;
     edm::Wrapper<BNbxlumi> bxlumidummy1;
     std::vector<BNbxlumi> bxlumidummy2;
     edm::Wrapper<std::vector<BNbxlumi> > bxlumidummy3;

     BNelectron electrondummy0;
     edm::Wrapper<BNelectron> electrondummy1;
     std::vector<BNelectron> electrondummy2;
     edm::Wrapper<std::vector<BNelectron> > electrondummy3;

     BNtau taudummy0;
     edm::Wrapper<BNtau> taudummy1;
     std::vector<BNtau> taudummy2;
     edm::Wrapper<std::vector<BNtau> > taudummy3;

     BNjet jetdummy0;
     edm::Wrapper<BNjet> jetdummy1;
     std::vector<BNjet> jetdummy2;
     edm::Wrapper<std::vector<BNjet> > jetdummy3;

     BNevent eventdummy0;
     edm::Wrapper<BNevent> eventdummy1;
     std::vector<BNevent> eventdummy2;
     edm::Wrapper<std::vector<BNevent> > eventdummy3;

     BNmcparticle mcparticledummy0;
     edm::Wrapper<BNmcparticle> mcparticledummy1;
     std::vector<BNmcparticle> mcparticledummy2;
     edm::Wrapper<std::vector<BNmcparticle> > mcparticledummy3;

     BNmet metdummy0;
     edm::Wrapper<BNmet> metdummy1;
     std::vector<BNmet> metdummy2;
     edm::Wrapper<std::vector<BNmet> > metdummy3;

     BNmuon muondummy0;
     edm::Wrapper<BNmuon> muondummy1;
     std::vector<BNmuon> muondummy2;
     edm::Wrapper<std::vector<BNmuon> > muondummy3;

     BNlepton leptondummy0;
     edm::Wrapper<BNlepton> leptondummy1;
     std::vector<BNlepton> leptondummy2;
     edm::Wrapper<std::vector<BNlepton> > leptondummy3;

     BNphoton photondummy0;
     edm::Wrapper<BNphoton> photondummy1;
     std::vector<BNphoton> photondummy2;
     edm::Wrapper<std::vector<BNphoton> > photondummy3;

     BNsupercluster superclusterdummy0;
     edm::Wrapper<BNsupercluster> superclusterdummy1;
     std::vector<BNsupercluster> superclusterdummy2;
     edm::Wrapper<std::vector<BNsupercluster> > superclusterdummy3;

     BNtrack trackdummy0;
     edm::Wrapper<BNtrack> trackdummy1;
     std::vector<BNtrack> trackdummy2;
     edm::Wrapper<std::vector<BNtrack> > trackdummy3;

     BNtrigger triggerdummy0;
     edm::Wrapper<BNtrigger> triggerdummy1;
     std::vector<BNtrigger> triggerdummy2;
     edm::Wrapper<std::vector<BNtrigger> > triggerdummy3;

     BNskimbit skimbitdummy0;
     edm::Wrapper<BNskimbit> skimbitdummy1;
     std::vector<BNskimbit> skimbitdummy2;
     edm::Wrapper<std::vector<BNskimbit> > skimbitdummy3;

     BNtrigobj trigobjdummy0;
     edm::Wrapper<BNtrigobj> trigobjdummy1;
     std::vector<BNtrigobj> trigobjdummy2;
     edm::Wrapper<std::vector<BNtrigobj> > trigobjdummy3;

     BNprimaryvertex primaryvertexdummy0;
     edm::Wrapper<BNprimaryvertex> primaryvertexdummy1;
     std::vector<BNprimaryvertex> primaryvertexdummy2;
     edm::Wrapper<std::vector<BNprimaryvertex> > primaryvertexdummy3;

     BNgenjet genjetdummy0;
     edm::Wrapper<BNgenjet> genjetdummy1;
     std::vector<BNgenjet> genjetdummy2;
     edm::Wrapper<std::vector<BNgenjet> > genjetdummy3;

/*
    These classes are commented out because they are used more rarely. If you need them, move them
    outside the comments and make the corresponding change in classes_def.xml
      
uncomment_h_here

      edm::Ref<std::vector<SampleProd> > dummy4;
      edm::RefVector<std::vector<SampleProd> > dummy5;
      edm::RefProd<std::vector<SampleProd> > dummy6;
*/

   };
}
