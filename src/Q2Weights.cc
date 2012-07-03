// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "TMath.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "LHAPDF/LHAPDF.h"

//
// class declaration
//

class Q2Weights : public edm::EDProducer{
   public:
      explicit Q2Weights(const edm::ParameterSet&);
      ~Q2Weights();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:

    struct GenJetHigherPt {
      bool operator() (const reco::GenJet& j1, const reco::GenJet& j2) const {
	return j1.pt() > j2.pt();
      };
    };
    struct GenLeptonHigherPt {
      bool operator() (const math::XYZTLorentzVector l1, const math::XYZTLorentzVector l2) const {
	return l1.pt() > l2.pt();
      };
    };
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);


 };

//
// constructors and destructor
//
Q2Weights::Q2Weights(const edm::ParameterSet& iConfig)
{
  LHAPDF::initPDFSet("cteq6ll.LHpdf");
  LHAPDF::usePDFMember(1,0);

  produces<std::vector<double> >();
}


Q2Weights::~Q2Weights()
{
 
 
}

//
// member functions
//
 double deltaR(math::XYZTLorentzVector v1,  math::XYZTLorentzVector v2) {
   Double_t DeltaR = 0.0;
   Double_t DeltaPhi = TMath::Abs(v1.Phi()-v2.Phi());
   if (DeltaPhi>TMath::Pi())
     DeltaPhi = 2*TMath::Pi() - DeltaPhi;
   Double_t DeltaEta = v1.Eta() - v2.Eta();
   DeltaR = sqrt(TMath::Power(DeltaPhi, 2) + TMath::Power(DeltaEta, 2));
   return DeltaR;
 }


// ------------ method called for each event  ------------
void
Q2Weights::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  
  
  float Q_scale;
  float factorization_scale;

  float Nparton;
  
  float x1,x2;
  int id1,id2;

  edm::Handle<LHEEventProduct> lheevent;
  bool hasLHEevent = iEvent.getByLabel("source", lheevent);
  
  double weight_up = 1., weight_down = 1.;

  if( hasLHEevent ){
    edm::Handle<GenParticleCollection> genPart;
    iEvent.getByLabel("genParticles", genPart);
  
    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByLabel("generator", genInfo); 
    const GenEventInfoProduct& genEventInfo = *(genInfo.product());
  

    factorization_scale = lheevent->hepeup().SCALUP;
    Q_scale = genEventInfo.qScale();
  
    const gen::PdfInfo* pdf = genEventInfo.pdf();
    id1=pdf->id.first;
    id2=pdf->id.second;
    x1=pdf->x.first;
    x2=pdf->x.second;
    //std::cout << factorization_scale << "  " << Q_scale << "  " << pdf->scalePDF << std::endl;
  
    auto_ptr<vector<math::XYZTLorentzVector> > genLepton ( new vector<math::XYZTLorentzVector> );
  
    Nparton=0;

    for(size_t i = 0; i < genPart->size(); ++ i) {
      //const GenParticle & p = (*genParticles)[i];
      const Candidate & p = (*genPart)[i];
      math::XYZTLorentzVector lvec(p.p4());
      int id = p.pdgId();
      //cout << i << " " <<id << " " << p.status() << " " << p.px() << " " <<p.py() << " " <<p.pz() << " " << p.energy() << " " << p.mass() << std::endl;
      bool is_finalstatus3 = false;
      if ( p.status() == 3) { 
	if(abs(id) == 11 || abs(id)==13 || abs(id)==15) { 
	  genLepton->push_back(lvec); 
	}
      
	is_finalstatus3 = true;
      
	for(size_t j=0; j < p.numberOfDaughters(); j++){
	  //cout << "  " <<  p.daughter(j)->pdgId() << " " << p.daughter(j)->status() << endl;
	  if(p.daughter(j)->status()==3){
	    is_finalstatus3=false;
	  }
	}
	if(is_finalstatus3) {
	  Nparton++;
	}   
      }
    }  

    double alpha = LHAPDF::alphasPDF(Q_scale);
    double alpha_down = LHAPDF::alphasPDF(Q_scale*0.5);
    double alpha_up = LHAPDF::alphasPDF(Q_scale*2.);

    //alhpa_S
    int num=4;
    weight_up = pow(alpha_up,Nparton-num)/pow(alpha,Nparton-num);
    weight_down = pow(alpha_down,Nparton-num)/pow(alpha,Nparton-num);
  
    //factorization scale pdf

    double pdf1 = LHAPDF::xfx(1,  x1, factorization_scale,  id1);
    double pdf2 = LHAPDF::xfx(1,  x2, factorization_scale,  id2);
    double pdf1_up = LHAPDF::xfx(1,  x1, factorization_scale*2.,  id1);
    double pdf2_up = LHAPDF::xfx(1,  x2, factorization_scale*2.,  id2);
    double pdf1_down = LHAPDF::xfx(1,  x1, factorization_scale*.5,  id1);
    double pdf2_down = LHAPDF::xfx(1,  x2, factorization_scale*.5,  id2);
  
    weight_up*=pdf1_up/pdf1*pdf2_up/pdf2;
    weight_down*=pdf1_down/pdf1*pdf2_down/pdf2; 
 
  } // end if( hasLHEevent )
 
  auto_ptr<std::vector<double> > Q2Weights (new std::vector<double>);

  Q2Weights->push_back(weight_up);

  Q2Weights->push_back(weight_down);

  iEvent.put(Q2Weights);
}

// ------------ method called once each job just before starting event loop  ------------
void 
Q2Weights::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Q2Weights::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
Q2Weights::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
Q2Weights::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Q2Weights::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Q2Weights::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Q2Weights::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Q2Weights);

