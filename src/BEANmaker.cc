// -*- C++ -*-
//
// Package:    BEANmaker
// Class:      BEANmaker
// 
/**\class BEANmaker BEANmaker.cc NtupleMaker/BEANmaker/src/BEANmaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Darren Michael Puigh
//         Created:  Wed Oct 28 18:09:28 CET 2009
// $Id: BEANmaker.cc,v 1.2 2011/10/05 01:28:37 puigh Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "ProductArea/BNcollections/interface/BNelectron.h"
#include "ProductArea/BNcollections/interface/BNevent.h"
#include "ProductArea/BNcollections/interface/BNjet.h"
#include "ProductArea/BNcollections/interface/BNmcparticle.h"
#include "ProductArea/BNcollections/interface/BNmet.h"
#include "ProductArea/BNcollections/interface/BNmuon.h"
#include "ProductArea/BNcollections/interface/BNphoton.h"
#include "ProductArea/BNcollections/interface/BNsupercluster.h"
#include "ProductArea/BNcollections/interface/BNtrack.h"
#include "ProductArea/BNcollections/interface/BNtrigger.h"
#include "ProductArea/BNcollections/interface/BNbxlumi.h"
#include "ProductArea/BNcollections/interface/BNtrigobj.h"
#include "ProductArea/BNcollections/interface/BNprimaryvertex.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerEvmReadoutRecord.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTReadoutCollection.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/METReco/interface/BeamHaloSummary.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "SusyAnalysis/EventSelector/interface/uncorrectionTypeMET.h"

#include "DataFormats/Common/interface/View.h"

#include <vector>
#include <map>
#include <string>



//
// class declaration
//

class BEANmaker : public edm::EDProducer {
   public:
      explicit BEANmaker(const edm::ParameterSet&);
      ~BEANmaker();

   private:
      virtual void beginJob() ;
      virtual void beginRun( edm::Run& , const edm::EventSetup& ) ;
      virtual void endRun( edm::Run& , edm::EventSetup& ) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;


      // ----------member data ---------------------------
  int count_hits( const std::vector<CaloTowerPtr> & towers );
  bool isActive(int word, int bit);
  // Active boards DAQ record bit number:
  // 0 FDL 
  // 1 PSB_0 9 Techn.Triggers for FDL
  // 2 PSB_1 13 Calo data for GTL
  // 3 PSB_2 14 Calo data for GTL
  // 4 PSB_3 15 Calo data for GTL
  // 5 PSB_4 19 M/Q bits for GMT
  // 6 PSB_5 20 M/Q bits for GMT
  // 7 PSB_6 21 M/Q bits for GMT
  // 8 GMT
  enum activeDAQ { FDL=0, PSB9, PSB13, PSB14, PSB15, PSB19, PSB20, PSB21, GMT };

  edm::InputTag eleTag_;
  edm::InputTag calojetTag_;
  edm::InputTag pfjetTag_;
  edm::InputTag calometTag_;
  edm::InputTag pfmetTag_;
  edm::InputTag tcmetTag_;
  edm::InputTag muonTag_;
  edm::InputTag photonTag_;
  edm::InputTag EBsuperclusterTag_;
  edm::InputTag EEsuperclusterTag_;
  edm::InputTag trackTag_;
  edm::InputTag genParticleTag_;
  edm::InputTag triggerSummaryTag_;
  edm::InputTag triggerResultsTag_;
  edm::InputTag gtSource_;
  edm::InputTag pvTag_;
  edm::InputTag dcsTag_;
  edm::TriggerNames trigNames;

  edm::InputTag reducedBarrelRecHitCollection_;
  edm::InputTag reducedEndcapRecHitCollection_;

  std::string hltProcessName_ ;

  double eventWeight_;
  double minTrackPt_;
  double minJetPt_;
  double minSCEt_;
  double minPhotonEt_;
  int sample_;
  double minNDOF_;
  double maxAbsZ_;
  double maxd0_;
  int numtrack_;
  double thresh_;
  bool verbose_;

  std::vector<std::string> hlt_sd_eg_name;
  std::vector<std::string> hlt_sd_jetmettau_name;

  int numEvents, numGoodVertexEvents, numFilterOutScrapingEvents;
  int numHEEPele, numHEEPeleEB, numHEEPeleEE;

  HLTConfigProvider hltConfig_;

};

//
// constants, enums and typedefs
//
typedef std::vector<BNelectron>     BNelectronCollection;
typedef std::vector<BNjet>          BNjetCollection;
typedef std::vector<BNevent>        BNeventCollection;
typedef std::vector<BNmcparticle>   BNmcparticleCollection;
typedef std::vector<BNmet>          BNmetCollection;
typedef std::vector<BNmuon>         BNmuonCollection;
typedef std::vector<BNphoton>       BNphotonCollection;
typedef std::vector<BNsupercluster> BNsuperclusterCollection;
typedef std::vector<BNtrack>        BNtrackCollection;
typedef std::vector<BNtrigger>      BNtriggerCollection;
typedef std::vector<BNbxlumi>      BNbxlumiCollection;
typedef std::vector<BNtrigobj>      BNtrigobjCollection;
typedef std::vector<BNprimaryvertex> BNprimaryvertexCollection;

//
// static data member definitions
//
static const char* kSC        = "corHybridSCandMulti5x5WithPreshower";
static const char* kTrigger   = "HLT";
static const char* kTriggerL1Talgo = "L1Talgo";
static const char* kTriggerL1Ttech = "L1Ttech";
static const char* kBXlumi    = "BXlumi";
static const char* kMCpar     = "MCstatus3";
static const char* kMCele     = "MCeleStatus1";
static const char* kMCmu      = "MCmuStatus1";


//
// constructors and destructor
//
BEANmaker::BEANmaker(const edm::ParameterSet& iConfig):
  eventWeight_(iConfig.getParameter<double>("eventWeight")),
  minTrackPt_(iConfig.getParameter<double>("minTrackPt")),
  minJetPt_(iConfig.getParameter<double>("minJetPt")),
  minSCEt_(iConfig.getParameter<double>("minSCEt")),
  minPhotonEt_(iConfig.getParameter<double>("minPhotonEt")),
  sample_(iConfig.getParameter<int>("sample")),
  minNDOF_(iConfig.getUntrackedParameter<double>("minNDOF",4)),
  maxAbsZ_(iConfig.getUntrackedParameter<double>("maxAbsZ",15)),
  maxd0_(iConfig.getUntrackedParameter<double>("maxd0",2)),
  numtrack_(iConfig.getUntrackedParameter<int>("numtrack",10)),
  thresh_(iConfig.getUntrackedParameter<double>("thresh",0.25)),
  verbose_(iConfig.getParameter<bool>("verbose"))
{

  // Define InputTags 
  eleTag_ = iConfig.getParameter<edm::InputTag>("eleTag");
  calojetTag_ = iConfig.getParameter<edm::InputTag>("calojetTag");
  pfjetTag_ = iConfig.getParameter<edm::InputTag>("pfjetTag");
  calometTag_ = iConfig.getParameter<edm::InputTag>("calometTag");
  pfmetTag_ = iConfig.getParameter<edm::InputTag>("pfmetTag");
  tcmetTag_ = iConfig.getParameter<edm::InputTag>("tcmetTag");
  muonTag_ = iConfig.getParameter<edm::InputTag>("muonTag");
  photonTag_ = iConfig.getParameter<edm::InputTag>("photonTag");
  EBsuperclusterTag_ = iConfig.getParameter<edm::InputTag>("EBsuperclusterTag");
  EEsuperclusterTag_ = iConfig.getParameter<edm::InputTag>("EEsuperclusterTag");
  trackTag_ = iConfig.getParameter<edm::InputTag>("trackTag");
  genParticleTag_ = iConfig.getParameter<edm::InputTag>("genParticleTag");
  triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("triggerResultsTag");
  gtSource_ = iConfig.getParameter<edm::InputTag>("gtSource");
  pvTag_ = iConfig.getParameter<edm::InputTag>("pvTag");
  triggerSummaryTag_ = iConfig.getParameter<edm::InputTag>("triggerSummaryTag");
  dcsTag_ = iConfig.getParameter<edm::InputTag>("dcsTag");

  reducedBarrelRecHitCollection_ = iConfig.getParameter<edm::InputTag>("reducedBarrelRecHitCollection");
  reducedEndcapRecHitCollection_ = iConfig.getParameter<edm::InputTag>("reducedEndcapRecHitCollection");

  hltProcessName_ = iConfig.getParameter<std::string>("hltProcessName");


  // Register products
  produces<BNeventCollection>().setBranchAlias("event");
  produces<BNelectronCollection>(eleTag_.label()).setBranchAlias("electrons");
  produces<BNjetCollection>(calojetTag_.label()).setBranchAlias("calojets");
  produces<BNjetCollection>(pfjetTag_.label()).setBranchAlias("pfjets");
  produces<BNmetCollection>(calometTag_.label()).setBranchAlias("calomet");
  produces<BNmetCollection>(pfmetTag_.label()).setBranchAlias("pfmet");
  produces<BNmetCollection>(tcmetTag_.label()).setBranchAlias("tcmet");
  produces<BNmuonCollection>(muonTag_.label()).setBranchAlias("muons");
  produces<BNphotonCollection>(photonTag_.label()).setBranchAlias("photons");
  produces<BNsuperclusterCollection>(kSC).setBranchAlias("superclusters");
  produces<BNtrackCollection>(trackTag_.label()).setBranchAlias("tracks");
  produces<BNtriggerCollection>(kTrigger).setBranchAlias("trigger");
  produces<BNtriggerCollection>(kTriggerL1Talgo).setBranchAlias("L1Talgo");
  produces<BNtriggerCollection>(kTriggerL1Ttech).setBranchAlias("L1Ttech");
  produces<BNbxlumiCollection>(kBXlumi).setBranchAlias("bxlumi");
  produces<BNmcparticleCollection>(kMCpar).setBranchAlias("mcparticles");
  produces<BNmcparticleCollection>(kMCele).setBranchAlias("mcelectrons");
  produces<BNmcparticleCollection>(kMCmu).setBranchAlias("mcmuons");
  produces<BNtrigobjCollection>().setBranchAlias("trigobj");
  produces<BNprimaryvertexCollection>(pvTag_.label()).setBranchAlias("pvs");


}


BEANmaker::~BEANmaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//



// ------------ method called to produce the data  ------------
void
BEANmaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;


   // Get Object Handles
   edm::Handle<edm::View<pat::Electron> > electronHandle;
   iEvent.getByLabel(eleTag_,electronHandle);
   edm::View<pat::Electron> electrons = *electronHandle;

   edm::Handle<edm::View<pat::Jet> > calojetHandle;
   iEvent.getByLabel(calojetTag_,calojetHandle);
   edm::View<pat::Jet> calojets = *calojetHandle;

   edm::Handle<edm::View<pat::Jet> > pfjetHandle;
   iEvent.getByLabel(pfjetTag_,pfjetHandle);
   edm::View<pat::Jet> pfjets = *pfjetHandle;

   edm::Handle<edm::View<pat::MET> > calometHandle;
   iEvent.getByLabel(calometTag_,calometHandle);

   edm::Handle<edm::View<pat::MET> > pfmetHandle;
   iEvent.getByLabel(pfmetTag_,pfmetHandle);

   edm::Handle<edm::View<pat::MET> > tcmetHandle;
   iEvent.getByLabel(tcmetTag_,tcmetHandle);

   edm::Handle<edm::View<pat::Muon> > muonHandle;
   iEvent.getByLabel(muonTag_,muonHandle);
   edm::View<pat::Muon> muons = *muonHandle;

   edm::Handle<edm::View<pat::Photon> > photonHandle;
   iEvent.getByLabel(photonTag_,photonHandle);
   edm::View<pat::Photon> photons = *photonHandle;

   edm::Handle<reco::TrackCollection > trackHandle;
   iEvent.getByLabel(trackTag_,trackHandle);
   reco::TrackCollection tracks = *trackHandle;

   edm::Handle<reco::SuperClusterCollection > EBsuperclusterHandle;
   iEvent.getByLabel(EBsuperclusterTag_,EBsuperclusterHandle);
   reco::SuperClusterCollection EBsuperclusters = *EBsuperclusterHandle;

   edm::Handle<reco::SuperClusterCollection > EEsuperclusterHandle;
   iEvent.getByLabel(EEsuperclusterTag_,EEsuperclusterHandle);
   reco::SuperClusterCollection EEsuperclusters = *EEsuperclusterHandle;

   edm::Handle<reco::GenParticleCollection > genParticles;
   iEvent.getByLabel(genParticleTag_,genParticles);


   // // remove all cluster tools for now
   EcalClusterLazyTools lazyTools( iEvent, iSetup, reducedBarrelRecHitCollection_, reducedEndcapRecHitCollection_ );
   Handle<EcalRecHitCollection> Brechit;//barrel
   Handle<EcalRecHitCollection> Erechit;//endcap
   iEvent.getByLabel(reducedBarrelRecHitCollection_,Brechit);
   iEvent.getByLabel(reducedEndcapRecHitCollection_,Erechit);
   const EcalRecHitCollection* barrelRecHits= Brechit.product();
   const EcalRecHitCollection* endcapRecHits= Erechit.product();


   // Get Trigger and Event Handles
   edm::ESHandle<L1GtTriggerMenu> menuRcd;
   iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
   const L1GtTriggerMenu* menu = menuRcd.product();

   // Get DCS information
   edm::Handle<DcsStatusCollection> dcsHandle;
   bool gotDCSinfo = iEvent.getByLabel(dcsTag_, dcsHandle);

   double evt_bField;


   ///  Luminosity information
   double instLumi = -1;
   double bxLumi = -1;

   std::auto_ptr<BNbxlumiCollection> bnbxlumis(new BNbxlumiCollection);

   edm::Handle<LumiDetails> d;
   iEvent.getLuminosityBlock().getByLabel("lumiProducer",d); 
   if( (d.isValid()) ){

     instLumi = 0;

     bxLumi = d->lumiValue(LumiDetails::kOCC1,iEvent.bunchCrossing())*6.37;
     //calibrated here, but not corrected in Hz/ub

     for( int i=0;i<3564;++i ){// Loop on bunch crossings

       BNbxlumi MyBXlumi;

       //calibrated here but not corrected, in Hz/ub
       double bx_LUMI_now = d->lumiValue(LumiDetails::kOCC1,i)*6.37;

       MyBXlumi.bx_B1_now = d->lumiBeam1Intensity(i);
       MyBXlumi.bx_B2_now = d->lumiBeam2Intensity(i);
       MyBXlumi.bx_LUMI_now= bx_LUMI_now;

       instLumi += bx_LUMI_now;

       bnbxlumis->push_back(MyBXlumi);
     }
   }



   Handle<std::vector< PileupSummaryInfo > >  PupInfo;
   iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

   std::vector<PileupSummaryInfo>::const_iterator PVI;

   int npv = -1;
   float sum_nvtx = 0;
   int nm1 = -1; int n0 = -1; int np1 = -1;

   if( (PupInfo.isValid()) ){

     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

       int BX = PVI->getBunchCrossing();

       npv = PVI->getPU_NumInteractions();

       sum_nvtx += float(npv);

       if(BX == -1) { 
	 nm1 = PVI->getPU_NumInteractions();
       }
       if(BX == 0) { 
	 n0 = PVI->getPU_NumInteractions();
       }
       if(BX == 1) { 
	 np1 = PVI->getPU_NumInteractions();
       }

     }

     std::cout << "\t nm1 = " << nm1 << std::endl;
     std::cout << "\t n0  = " << n0  << std::endl;
     std::cout << "\t np1 = " << np1 << std::endl;

     std::cout << "\t sum_nvtx = " << sum_nvtx << std::endl;
   }



   if( gotDCSinfo && sample_<0 && dcsHandle->size()>0 ){
     // scale factor = 3.801/18166.0 which are
     // average values taken over a stable two
     // week period
     float currentToBFieldScaleFactor = 2.09237036221512717e-04;
     float current = (*dcsHandle)[0].magnetCurrent();
     evt_bField = current*currentToBFieldScaleFactor;
   }
   else {
     ESHandle<MagneticField> magneticField;
     iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

     evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
   }
   
   math::XYZPoint beamSpotPosition;
   beamSpotPosition.SetCoordinates(0,0,0);

   edm::Handle<reco::BeamSpot> bsHandle;
   iEvent.getByLabel("offlineBeamSpot",bsHandle);

   double BSx=0,BSy=0,BSz=0;
   if( (bsHandle.isValid()) ){
     reco::BeamSpot bs = *bsHandle;
     BSx = bs.x0();
     BSy = bs.y0();
     BSz = bs.z0();
     beamSpotPosition = bsHandle->position();
   }


   edm::Handle<GenEventInfoProduct> genEvtInfo;
   bool hasGenEvtInfo = iEvent.getByLabel( "generator", genEvtInfo );

   double qScale=-1, alphaQCD=-1, alphaQED=-1, pthat=-1, scalePDF=-1, x1=-1, x2=-1, xPDF1=-1, xPDF2=-1;
   if( (hasGenEvtInfo) ){
     qScale = genEvtInfo->qScale();
     alphaQCD = genEvtInfo->alphaQCD();
     alphaQED = genEvtInfo->alphaQED();
     pthat = ( genEvtInfo->hasBinningValues() ? (genEvtInfo->binningValues())[0] : 0.0);
     scalePDF = genEvtInfo->pdf()->scalePDF;
     x1 = genEvtInfo->pdf()->x.first;
     x2 = genEvtInfo->pdf()->x.second;
     xPDF1 = genEvtInfo->pdf()->xPDF.first;
     xPDF2 = genEvtInfo->pdf()->xPDF.second;
   }


   Handle< bool > hcalNoiseFilterHandle;
   iEvent.getByLabel("HBHENoiseFilterResultProducer","HBHENoiseFilterResult", hcalNoiseFilterHandle);

   Bool_t hcalNoiseFilter = false;
   if( hcalNoiseFilterHandle.isValid() ) hcalNoiseFilter = (Bool_t)(*hcalNoiseFilterHandle);
   else std::cout << " NoiseFilter =====> hcalNoiseFilterHandle.isValid()==false " << std::endl;


   Handle<HcalNoiseSummary> HcalNoiseSummaryHandle;
   iEvent.getByLabel("hcalnoise", HcalNoiseSummaryHandle);

   bool passLooseNoiseFilter=false, passTightNoiseFilter=false;
   if( HcalNoiseSummaryHandle.isValid() ){
     passLooseNoiseFilter = HcalNoiseSummaryHandle->passLooseNoiseFilter();
     passTightNoiseFilter = HcalNoiseSummaryHandle->passTightNoiseFilter();
   }
   else std::cout << " NoiseFilter =====> HcalNoiseSummaryHandle.isValid()==false " << std::endl;



   // Get BeamHaloSummary 
   edm::Handle<BeamHaloSummary> TheBeamHaloSummary ;
   iEvent.getByLabel("BeamHaloSummary", TheBeamHaloSummary) ;

   bool passCSCLooseHaloId=false, passCSCTightHaloId=false, passEcalLooseHaloId=false, passEcalTightHaloId=false, passHcalLooseHaloId=false;
   bool passHcalTightHaloId=false, passGlobalLooseHaloId=false, passGlobalTightHaloId=false, passLooseId=false, passTightId=false;
   if( TheBeamHaloSummary.isValid() ) {
     const BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product() );

     if( TheSummary.CSCLooseHaloId() ) passCSCLooseHaloId = true;
     if( TheSummary.CSCTightHaloId() ) passCSCTightHaloId = true;
     if( TheSummary.EcalLooseHaloId() ) passEcalLooseHaloId = true;
     if( TheSummary.EcalTightHaloId() ) passEcalTightHaloId = true;
     if( TheSummary.HcalLooseHaloId() ) passHcalLooseHaloId = true;
     if( TheSummary.HcalTightHaloId() ) passHcalTightHaloId = true;
     if( TheSummary.GlobalLooseHaloId() ) passGlobalLooseHaloId = true;
     if( TheSummary.GlobalTightHaloId() ) passGlobalTightHaloId = true;
     if( TheSummary.LooseId() ) passLooseId = true;
     if( TheSummary.TightId() ) passTightId = true;

   }
   else std::cout << " NoiseFilter =====> TheBeamHaloSummary.isValid()==false " << std::endl;


   // std::cout << "\t hcalNoiseFilter = " << hcalNoiseFilter << std::endl;
   // std::cout << "" << std::endl;

   // std::cout << "\t passLooseNoiseFilter = " << passLooseNoiseFilter << std::endl;
   // std::cout << "\t passTightNoiseFilter = " << passTightNoiseFilter << std::endl;
   // std::cout << "" << std::endl;

   // std::cout << "\t passCSCLooseHaloId = " << passCSCLooseHaloId << std::endl;
   // std::cout << "\t passCSCTightHaloId = " << passCSCTightHaloId << std::endl;
   // std::cout << "\t passEcalLooseHaloId = " << passEcalLooseHaloId << std::endl;
   // std::cout << "\t passEcalTightHaloId = " << passEcalTightHaloId << std::endl;
   // std::cout << "\t passHcalLooseHaloId = " << passHcalLooseHaloId << std::endl;
   // std::cout << "\t passHcalTightHaloId = " << passHcalTightHaloId << std::endl;
   // std::cout << "\t passGlobalLooseHaloId = " << passGlobalLooseHaloId << std::endl;
   // std::cout << "\t passGlobalTightHaloId = " << passGlobalTightHaloId << std::endl;
   // std::cout << "\t passLooseId = " << passLooseId << std::endl;
   // std::cout << "\t passTightId = " << passTightId << std::endl;


   edm::Handle<reco::VertexCollection> vtxHandle;
   iEvent.getByLabel(pvTag_,vtxHandle);
   reco::VertexCollection vtxs = *vtxHandle;


   bool GoodVertex = false;

   int numPVs = 0;
   double PVz = 0;
   math::XYZPoint vertexPosition;
   vertexPosition.SetCoordinates(0,0,0);
   std::auto_ptr<BNprimaryvertexCollection> bnpvs(new BNprimaryvertexCollection);
   if( (vtxHandle.isValid()) ){
     numPVs = vtxHandle->size();

     bool firstPV = true;
     for( reco::VertexCollection::const_iterator vtx = vtxs.begin(); vtx!=vtxs.end(); ++vtx ){

       BNprimaryvertex MyPV;

       MyPV.x = vtx->x();
       MyPV.xError = vtx->xError();
       MyPV.y = vtx->y();
       MyPV.yError = vtx->yError();
       MyPV.z = vtx->z();
       MyPV.zError = vtx->zError();

       MyPV.normalizedChi2 = vtx->normalizedChi2();
       MyPV.isFake = vtx->isFake();
       MyPV.isValid = vtx->isValid();
       MyPV.tracksSize = vtx->tracksSize();

       MyPV.ndof = vtx->ndof();
       MyPV.rho  = vtx->position().rho();

       bool isGood = false;
       if( vtx->ndof() > minNDOF_ &&
	   ( (maxAbsZ_ <= 0) || fabs(vtx->z()) <= maxAbsZ_ ) &&
	   ( (maxd0_ <= 0) || fabs(vtx->position().rho()) <= maxd0_ )
	   ){
	 GoodVertex = true;
	 isGood = true;
       }

       MyPV.isGood = ( isGood ) ? 1 : 0;

       if( firstPV && isGood ){
	 vertexPosition = vtx->position();
	 PVz = vtx->z();
	 firstPV = false;
       }

       bnpvs->push_back(MyPV);
     }
   }

   double weight = -1;
   if(eventWeight_ < 0.0) {
     weight = 1.0;
   } 
   else {
     weight = eventWeight_;
   }


   numEvents++;

   //  Fill the electron collection
   int numele=0;
   int ele_index=0;
   std::auto_ptr<BNelectronCollection> bnelectrons(new BNelectronCollection);
   for( edm::View<pat::Electron>::const_iterator ele = electrons.begin(); ele!=electrons.end(); ++ele ){

     ele_index++;

     double elePin = ele->trackMomentumAtVtx().R();
     double elePout = ele->trackMomentumOut().R();

     BNelectron MyElectron;

     // general kinematic variables
     MyElectron.energy = ele->energy();
     MyElectron.gsfEt = ele->et();
     MyElectron.pt = ele->pt();
     MyElectron.px = ele->px();
     MyElectron.py = ele->py();
     MyElectron.pz = ele->pz();
     MyElectron.phi = ele->phi();
     MyElectron.eta = ele->eta();
     MyElectron.theta = ele->theta();

     MyElectron.charge = ele->charge();
     MyElectron.classification = ele->classification();
     MyElectron.vx = ele->vx();
     MyElectron.vy = ele->vy();
     MyElectron.vz = ele->vz();

     MyElectron.EscOverPin = ele->eSuperClusterOverP();
     MyElectron.EseedOverPout = ele->eSeedClusterOverPout();
     MyElectron.hadOverEm = ele->hadronicOverEm();
     MyElectron.fbrem = fabs(elePin-elePout)/elePin;
     MyElectron.delPhiIn = ele->deltaPhiSuperClusterTrackAtVtx();
     MyElectron.delEtaIn = ele->deltaEtaSuperClusterTrackAtVtx();
     MyElectron.eidRobustHighEnergy = ele->electronID("eidRobustHighEnergy");
     MyElectron.eidRobustLoose = ele->electronID("eidRobustLoose");
     MyElectron.eidRobustTight = ele->electronID("eidRobustTight");
     MyElectron.eidLoose = ele->electronID("eidLoose");
     MyElectron.eidTight = ele->electronID("eidTight");

     MyElectron.trackIso = ele->trackIso();
     MyElectron.ecalIso = ele->ecalIso();
     MyElectron.hcalIso = ele->hcalIso();
     MyElectron.caloIso = ele->caloIso();

     MyElectron.trackIsoDR03 = ele->dr03TkSumPt();
     MyElectron.ecalIsoDR03 = ele->dr03EcalRecHitSumEt();
     MyElectron.hcalIsoDR03 = ele->dr03HcalTowerSumEt();
     MyElectron.hcalIsoDR03depth1 = ele->dr03HcalDepth1TowerSumEt();
     MyElectron.hcalIsoDR03depth2 = ele->dr03HcalDepth2TowerSumEt();
     MyElectron.caloIsoDR03 = ele->dr03EcalRecHitSumEt()+ele->dr03HcalTowerSumEt();

     MyElectron.trackIsoDR04 = ele->dr04TkSumPt();
     MyElectron.ecalIsoDR04 = ele->dr04EcalRecHitSumEt();
     MyElectron.hcalIsoDR04 = ele->dr04HcalTowerSumEt();
     MyElectron.hcalIsoDR04depth1 = ele->dr04HcalDepth1TowerSumEt();
     MyElectron.hcalIsoDR04depth2 = ele->dr04HcalDepth2TowerSumEt();
     MyElectron.caloIsoDR04 = ele->dr04EcalRecHitSumEt()+ele->dr04HcalTowerSumEt();

     MyElectron.mva = ele->mva();
     MyElectron.numberOfLostHits = ele->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
     MyElectron.numberOfExpectedInnerHits = ele->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
     MyElectron.numberOfValidPixelHits = ele->gsfTrack()->trackerExpectedHitsInner().numberOfValidPixelHits();
     MyElectron.numberOfValidPixelBarrelHits = ele->gsfTrack()->trackerExpectedHitsInner().numberOfValidPixelBarrelHits();
     MyElectron.numberOfValidPixelEndcapHits = ele->gsfTrack()->trackerExpectedHitsInner().numberOfValidPixelEndcapHits();

     MyElectron.scEnergy = ele->caloEnergy();
     MyElectron.scSigmaEtaEta = ele->scSigmaEtaEta();
     MyElectron.scSigmaIEtaIEta = ele->scSigmaIEtaIEta();
     MyElectron.scE1x5 = ele->scE1x5();
     MyElectron.scE2x5Max = ele->scE2x5Max();
     MyElectron.scE5x5 = ele->scE5x5();
     MyElectron.numClusters = ele->numberOfBrems();

     double caloenergy = -1;

     if( (ele->superCluster().isAvailable()) ){
       double eMax    = lazyTools.eMax(    *(ele->superCluster()) );
       double eLeft   = lazyTools.eLeft(   *(ele->superCluster()) );
       double eRight  = lazyTools.eRight(  *(ele->superCluster()) );
       double eTop    = lazyTools.eTop(    *(ele->superCluster()) );
       double eBottom = lazyTools.eBottom( *(ele->superCluster()) );
       double e3x3    = lazyTools.e3x3(    *(ele->superCluster()) );
       double swissCross = -99;

       if( eMax>0.000001 ) swissCross = 1 - (eLeft+eRight+eTop+eBottom)/eMax;

       MyElectron.eMax = eMax;
       MyElectron.eLeft = eLeft;
       MyElectron.eRight = eRight;
       MyElectron.eTop = eTop;
       MyElectron.eBottom = eBottom;
       MyElectron.e3x3 = e3x3;
       MyElectron.swissCross = swissCross;

       caloenergy = ele->caloEnergy();

       MyElectron.scEt = ele->caloEnergy() * sin( ele->superCluster()->position().theta() );
       MyElectron.scRawEnergy = ele->superCluster()->rawEnergy();
       MyElectron.scEta = ele->superCluster()->position().eta();
       MyElectron.scPhi = ele->superCluster()->position().phi();
       MyElectron.scZ = ele->superCluster()->position().Z();

       double seedE = -999, seedTime = -999;
       int seedRecoFlag = -999;
       
       if( (ele->isEB()) ){
       	 DetId idEB = EcalClusterTools::getMaximum( ele->superCluster()->hitsAndFractions(), &(*barrelRecHits) ).first;
       	 EcalRecHitCollection::const_iterator thisHitEB = barrelRecHits->find(idEB);

       	 seedE = thisHitEB->energy();
       	 seedTime = thisHitEB->time();
       	 seedRecoFlag = thisHitEB->recoFlag();
       }
       else if( (ele->isEE()) ){
       	 DetId idEE = EcalClusterTools::getMaximum( ele->superCluster()->hitsAndFractions(), &(*endcapRecHits) ).first;
       	 EcalRecHitCollection::const_iterator thisHitEE = endcapRecHits->find(idEE);

       	 seedE = thisHitEE->energy();
       	 seedTime = thisHitEE->time();
       	 seedRecoFlag = thisHitEE->recoFlag();
       }

       MyElectron.seedEnergy   = seedE;
       MyElectron.seedTime     = seedTime;
       MyElectron.seedRecoFlag = seedRecoFlag;

       // MyElectron.swissCrossNoI85 = swissCrossNoI85;
       // MyElectron.swissCrossI85   = swissCrossI85;
       // MyElectron.E2overE9NoI85 = e2overe9NoI85;
       // MyElectron.E2overE9I85 = e2overe9I85;
     }

     double heepEt = ( caloenergy<0 || ele->energy()==0. ) ? 0 : ele->et()/ele->energy()*caloenergy;
     MyElectron.et = heepEt;

     if( (ele->closestCtfTrackRef().isAvailable()) ) MyElectron.tkCharge = ele->closestCtfTrackRef()->charge();
     if( (ele->gsfTrack().isAvailable()) ){
       double tkvx = ele->gsfTrack()->vx();
       double tkvy = ele->gsfTrack()->vy();
       double tkpx = ele->gsfTrack()->px();
       double tkpy = ele->gsfTrack()->py();
       double tkpt = ele->gsfTrack()->pt();

       double ndof = ele->gsfTrack()->ndof();
       if( (ndof!=0) ) MyElectron.tkNormChi2 = ele->gsfTrack()->chi2()/ndof;
       MyElectron.tkPT = ele->gsfTrack()->pt();
       MyElectron.tkEta = ele->gsfTrack()->eta();
       MyElectron.tkPhi = ele->gsfTrack()->phi();
       MyElectron.tkDZ = ele->gsfTrack()->dz();
       MyElectron.tkD0 = ele->gsfTrack()->d0();
       MyElectron.tkD0bs = (tkvx-BSx)*tkpy/tkpt-(tkvy-BSy)*tkpx/tkpt;
       MyElectron.tkD0err = ele->gsfTrack()->d0Error();
       MyElectron.tkNumValidHits = ele->gsfTrack()->numberOfValidHits();
       MyElectron.gsfCharge = ele->gsfTrack()->charge();

       MyElectron.correctedD0 = ele->gsfTrack()->dxy(beamSpotPosition);
       MyElectron.correctedDZ = ele->gsfTrack()->dz(vertexPosition);
     }

     MyElectron.isEB = ele->isEB();
     MyElectron.isEE = ele->isEE();
     MyElectron.isGap = ele->isGap();
     MyElectron.isEBEEGap = ele->isEBEEGap();
     MyElectron.isEBGap = ele->isEBGap();
     MyElectron.isEEGap = ele->isEEGap();
     MyElectron.isEcalDriven = ele->ecalDrivenSeed();
     MyElectron.isTrackerDriven = ele->trackerDrivenSeed();

     if( (ele->genLepton()) ){
       MyElectron.genId = ele->genLepton()->pdgId();
       MyElectron.genET = ele->genLepton()->et();
       MyElectron.genPT = ele->genLepton()->pt();
       MyElectron.genPhi = ele->genLepton()->phi();
       MyElectron.genEta = ele->genLepton()->eta();
       MyElectron.genCharge = ele->genLepton()->charge();

       if( (ele->genLepton()->numberOfMothers()>0) ){
	 const reco::Candidate *ElectronMother = ele->genLepton()->mother();
	 bool staytrapped = true;
	 while( (abs(ElectronMother->pdgId())==11 && staytrapped) ){
	   if( ElectronMother->numberOfMothers()==1 ) ElectronMother = ElectronMother->mother();
	   else staytrapped = false;
	 }
       
	 MyElectron.genMotherId = ElectronMother->pdgId();
	 MyElectron.genMotherET = ElectronMother->et();
	 MyElectron.genMotherPT = ElectronMother->pt();
	 MyElectron.genMotherPhi = ElectronMother->phi();
	 MyElectron.genMotherEta = ElectronMother->eta();
	 MyElectron.genMotherCharge = ElectronMother->charge();
       }
     }

     ConversionFinder convFinder;
     ConversionInfo convInfo = convFinder.getConversionInfo(*ele, trackHandle, evt_bField);
   
     double dist = convInfo.dist();
     double dcot = convInfo.dcot();
     double convradius = convInfo.radiusOfConversion();
     math::XYZPoint convPoint = convInfo.pointOfConversion();

     MyElectron.dist = dist;
     MyElectron.dcot = dcot;
     MyElectron.convradius = convradius;
     MyElectron.convPointX = convPoint.x();
     MyElectron.convPointY = convPoint.y();
     MyElectron.convPointZ = convPoint.z();


     bool cutDelEta=false, cutDelPhi=false, cutSigIeta=false, cutE2x5=false, cutEMhad1=false, cutHad2=false, cutTrackIso=false;

     double absEta = fabs(ele->superCluster()->position().eta());
     double EMhad1 = ele->dr03EcalRecHitSumEt() + ele->dr03HcalDepth1TowerSumEt();
     bool EB = ( absEta < 1.442 );
     bool EE = ( absEta > 1.560 && absEta < 2.5 );

     bool isECAL = ele->ecalDrivenSeed();

     bool cutET = ( heepEt > 25. );
     bool cutEta = ( EB || EE );
     bool cutHoverE = ( ele->hadronicOverEm() < 0.05 );

     if( EB ){
       cutDelEta = ( fabs(ele->deltaEtaSuperClusterTrackAtVtx())<0.005 );
       cutDelPhi = ( fabs(ele->deltaPhiSuperClusterTrackAtVtx())<0.09 );
       cutSigIeta = true;
       cutE2x5 = ( (ele->scE2x5Max()/ele->scE5x5())>0.94 || (ele->scE1x5()/ele->scE5x5())>0.83 );
       cutEMhad1 = ( EMhad1 < (2 + 0.03*heepEt) );
       cutHad2 = true;
       cutTrackIso = ( ele->dr03TkSumPt()<7.5 );
     }
     else if( EE ){
       cutDelEta = ( fabs(ele->deltaEtaSuperClusterTrackAtVtx())<0.007 );
       cutDelPhi = ( fabs(ele->deltaPhiSuperClusterTrackAtVtx())<0.09 );
       cutSigIeta = ( ele->scSigmaIEtaIEta()<0.03 );
       cutE2x5 = true;
       if( heepEt<50 ) cutEMhad1 = ( EMhad1 < 2.5 );
       else cutEMhad1 = ( EMhad1 < (2.5 + 0.03*(heepEt-50)) );
       cutHad2 = ( ele->dr03HcalDepth2TowerSumEt()<0.5 );
       cutTrackIso = ( ele->dr03TkSumPt()<15. );
     }

     bool isHEEPnoEt = ( isECAL && cutEta && cutDelEta && cutDelPhi && cutHoverE &&  
			 cutSigIeta && cutE2x5 && cutEMhad1 && cutHad2 && cutTrackIso );
     bool isHEEP = ( isHEEPnoEt && cutET );


     MyElectron.isHEEP = ( isHEEP ) ? 1 : 0;
     MyElectron.isHEEPnoEt = ( isHEEPnoEt ) ? 1 : 0;

     bnelectrons->push_back(MyElectron);

     if( isHEEP ){
       numele++;
       if( EB ) numHEEPeleEB++;
       else if( EE ) numHEEPeleEE++;
     }       
   }


   numHEEPele += numele;


   double metJES1corrPx = 0;
   double metJES1corrPy = 0;
   double metJES20corrPx = 0;
   double metJES20corrPy = 0;
   double UDeltaPx = 0;
   double UDeltaPy = 0;
   double USumET = 0;



   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl_Calo;
   iSetup.get<JetCorrectionsRecord>().get("AK5Calo",JetCorParColl_Calo); 
   JetCorrectorParameters const & JetCorPar_Calo = (*JetCorParColl_Calo)["Uncertainty"];
   JetCorrectionUncertainty *jecUnc_Calo = new JetCorrectionUncertainty(JetCorPar_Calo);


   //  Fill the calojet collection
   std::auto_ptr<BNjetCollection> bncalojets(new BNjetCollection);
   for( edm::View<pat::Jet>::const_iterator calojet = calojets.begin(); calojet != calojets.end(); ++ calojet ) {

     if( !(calojet->pt()>minJetPt_) ) continue;

     /*
     // Get uncorrected/raw jets
     pat::Jet rawJet = calojet->correctedJet("raw");
     double RawjetE = rawJet.energy();
     double RawjetEt = rawJet.et();
     double RawjetEta = rawJet.eta();
     double Rawjetpt = rawJet.pt();
     double Rawjetpx = rawJet.px();
     double Rawjetpy = rawJet.py();
     double RawjetEMF = rawJet.emEnergyFraction();


     // Perform L2(relative) and L3(absolute) and residual corrections
     JEC_L2_calo->setJetEta(RawjetEta);      JEC_L2_calo->setJetPt(Rawjetpt);
     JEC_L2L3_calo->setJetEta(RawjetEta);    JEC_L2L3_calo->setJetPt(Rawjetpt);
     JEC_L2L3res_calo->setJetEta(RawjetEta); JEC_L2L3res_calo->setJetPt(Rawjetpt);
     JEC_res_calo->setJetEta(calojet->eta()); JEC_res_calo->setJetPt(calojet->pt());

     double L2jetpt      = Rawjetpt * JEC_L2_calo->getCorrection();
     double L2L3jetpt    = Rawjetpt * JEC_L2L3_calo->getCorrection();
     double L2L3resjetpt = Rawjetpt * JEC_L2L3res_calo->getCorrection();
     double respt = calojet->pt() * JEC_res_calo->getCorrection();

     double L2L3jetpx = calojet->px();
     double L2L3jetpy = calojet->py();
 
     // Criteria on jets to be included in met correction
     if( Rawjetpt>1 && RawjetEMF<0.9 && (fabs(calojet->eta())>2.6 || RawjetEMF>0.01) ){
       metJES1corrPx += -(L2L3jetpx - Rawjetpx);
       metJES1corrPy += -(L2L3jetpy - Rawjetpy);
       if( Rawjetpt>20 ){
     	 metJES20corrPx += -(L2L3jetpx - Rawjetpx);
     	 metJES20corrPy += -(L2L3jetpy - Rawjetpy);
	 UDeltaPx += Rawjetpx;
	 UDeltaPy += Rawjetpy;
	 USumET -= RawjetEt;
       }
     }
     */

     jecUnc_Calo->setJetEta(calojet->eta());
     jecUnc_Calo->setJetPt(calojet->pt());// the uncertainty is a function of the corrected pt
     double unc = jecUnc_Calo->getUncertainty(true);


     BNjet MyCalojet;

     // general kinematic variables
     MyCalojet.energy = calojet->energy();
     MyCalojet.et = calojet->et();
     MyCalojet.pt = calojet->pt();
     MyCalojet.px = calojet->px();
     MyCalojet.py = calojet->py();
     MyCalojet.pz = calojet->pz();
     MyCalojet.phi = calojet->phi();
     MyCalojet.eta = calojet->eta();
     MyCalojet.theta = calojet->theta();

     MyCalojet.EMfrac = calojet->emEnergyFraction();
     MyCalojet.Hadfrac = calojet->energyFractionHadronic();
     MyCalojet.charge = calojet->jetCharge();
     MyCalojet.mass = calojet->mass();
     MyCalojet.area = calojet->towersArea();
     MyCalojet.fHPD = calojet->jetID().fHPD;
     MyCalojet.flavour = calojet->partonFlavour();
     MyCalojet.Nconst = calojet->nConstituents();
     MyCalojet.n90Hits = calojet->jetID().n90Hits;
     MyCalojet.approximatefHPD = calojet->jetID().approximatefHPD;
     MyCalojet.hitsInN90 = calojet->jetID().hitsInN90;
     MyCalojet.fLong = calojet->jetID().fLong;
     MyCalojet.fShort = calojet->jetID().fShort;
     MyCalojet.nHit = count_hits( calojet->getCaloConstituents() );
     MyCalojet.etaetaMoment = calojet->etaetaMoment();
     MyCalojet.phiphiMoment = calojet->phiphiMoment();


     // btag variables
     MyCalojet.btagTChighPur = calojet->bDiscriminator("trackCountingHighPurBJetTags");
     MyCalojet.btagTChighEff = calojet->bDiscriminator("trackCountingHighEffBJetTags");
     MyCalojet.btagJetProb = calojet->bDiscriminator("jetProbabilityBJetTags");
     MyCalojet.btagJetBProb = calojet->bDiscriminator("jetBProbabilityBJetTags");
     MyCalojet.btagSoftEle = calojet->bDiscriminator("softElectronBJetTags");
     MyCalojet.btagSoftMuon = calojet->bDiscriminator("softMuonBJetTags");
     MyCalojet.btagSoftMuonNoIP = calojet->bDiscriminator("softMuonNoIPBJetTags");
     MyCalojet.btagSecVertex = calojet->bDiscriminator("simpleSecondaryVertexBJetTags");
     MyCalojet.btagSecVertexHighEff = calojet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
     MyCalojet.btagSecVertexHighPur = calojet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");

     /*
     MyCalojet.Upt = Rawjetpt;
     MyCalojet.Uenergy = RawjetE;
     MyCalojet.L2pt = L2jetpt;
     MyCalojet.L2L3pt = L2L3jetpt;
     MyCalojet.L2L3respt = L2L3resjetpt;
     MyCalojet.respt = respt;
     */
     MyCalojet.JESunc = unc;

     if( (calojet->genJet()) ){ // if there is a matched genjet, fill variables
       MyCalojet.genJetET = calojet->genJet()->et();
       MyCalojet.genJetPT = calojet->genJet()->pt();
       MyCalojet.genJetEta = calojet->genJet()->eta();
       MyCalojet.genJetPhi = calojet->genJet()->phi();
     }
     if( (calojet->genParton()) ){ // if there is a matched parton, fill variables
       MyCalojet.genPartonET = calojet->genParton()->et();
       MyCalojet.genPartonPT = calojet->genParton()->pt();
       MyCalojet.genPartonEta = calojet->genParton()->eta();
       MyCalojet.genPartonPhi = calojet->genParton()->phi();
       MyCalojet.genPartonId = calojet->genParton()->pdgId();
       MyCalojet.genPartonMotherId = calojet->genParton()->mother()->pdgId();
       MyCalojet.genPartonGrandMotherId = calojet->genParton()->mother()->pdgId();
     }


     // std::cout << " jet
     /*
     // // Get jet ID results
     //const pat::Jet& ijet = *calojet;
     const pat::Jet& ijet = rawJet;

     ret.set(false);
     jetIDMinimal(ijet, ret);
     bool minimal = ( ret[str_MINIMAL_EMF] );
     ret.set(false);
     jetIDLooseAOD(ijet, ret);
     bool loose_aod = ( ret[str_LOOSE_AOD_fHPD] && ret[str_LOOSE_AOD_N90Hits] && ret[str_LOOSE_AOD_EMF] );
     ret.set(false);
     jetIDLoose(ijet, ret);
     bool loose = ( ret[str_LOOSE_fHPD] && ret[str_LOOSE_N90Hits] && ret[str_LOOSE_EMF] );
     ret.set(false);
     jetIDTight(ijet, ret);
     bool tight = ( ret[str_TIGHT_fHPD] && ret[str_TIGHT_EMF] );


     MyCalojet.jetIDMinimal = (minimal) ? 1 : 0;
     MyCalojet.jetIDLooseAOD = (loose_aod) ? 1 : 0;
     MyCalojet.jetIDLoose = (loose) ? 1 : 0;
     MyCalojet.jetIDTight = (tight) ? 1 : 0;
     */


     bncalojets->push_back(MyCalojet);
   }




   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl_PF;
   iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl_PF); 
   JetCorrectorParameters const & JetCorPar_PF = (*JetCorParColl_PF)["Uncertainty"];
   JetCorrectionUncertainty *jecUnc_PF = new JetCorrectionUncertainty(JetCorPar_PF);

   //  Fill the pfjet collection
   std::auto_ptr<BNjetCollection> bnpfjets(new BNjetCollection);
   for( edm::View<pat::Jet>::const_iterator pfjet = pfjets.begin(); pfjet != pfjets.end(); ++ pfjet ) {


     /*
     // Get uncorrected/raw jets
     pat::Jet rawJet = pfjet->correctedJet("raw");
     double Rawjetpt = rawJet.pt();
     double Rawjetpx = rawJet.px();
     double Rawjetpy = rawJet.py();

     if( Rawjetpt>1. ){
       metJet1corrPx_pf += -( pfjet->px() - Rawjetpx );
       metJet1corrPy_pf += -( pfjet->py() - Rawjetpy );
       
       if( Rawjetpt>6. ){
	 metJet6corrPx_pf += -( pfjet->px() - Rawjetpx );
	 metJet6corrPy_pf += -( pfjet->py() - Rawjetpy );

	 if( Rawjetpt>10. ){
	   metJet10corrPx_pf += -( pfjet->px() - Rawjetpx );
	   metJet10corrPy_pf += -( pfjet->py() - Rawjetpy );
	 }
       }
     }
     */

     if( !(pfjet->pt()>minJetPt_) ) continue;


     jecUnc_PF->setJetEta(pfjet->eta());
     jecUnc_PF->setJetPt(pfjet->pt()); // here you must use the CORRECTED jet pt
     double unc = jecUnc_PF->getUncertainty(true);


     BNjet MyPfjet;

     // general kinematic variables
     MyPfjet.energy = pfjet->energy();
     MyPfjet.et = pfjet->et();
     MyPfjet.pt = pfjet->pt();
     MyPfjet.px = pfjet->px();
     MyPfjet.py = pfjet->py();
     MyPfjet.pz = pfjet->pz();
     MyPfjet.phi = pfjet->phi();
     MyPfjet.eta = pfjet->eta();
     MyPfjet.theta = pfjet->theta();

     //MyPfjet.EMfrac = pfjet->emEnergyFraction();
     //MyPfjet.Hadfrac = pfjet->energyFractionHadronic();
     MyPfjet.charge = pfjet->jetCharge();
     MyPfjet.mass = pfjet->mass();
     //MyPfjet.area = pfjet->towersArea();
     MyPfjet.fHPD = pfjet->jetID().fHPD;
     MyPfjet.flavour = pfjet->partonFlavour();
     MyPfjet.Nconst = pfjet->nConstituents();
     MyPfjet.n90Hits = pfjet->jetID().n90Hits;
     MyPfjet.approximatefHPD = pfjet->jetID().approximatefHPD;
     MyPfjet.hitsInN90 = pfjet->jetID().hitsInN90;


     // btag variables
     MyPfjet.btagTChighPur = pfjet->bDiscriminator("trackCountingHighPurBJetTags");
     MyPfjet.btagTChighEff = pfjet->bDiscriminator("trackCountingHighEffBJetTags");
     MyPfjet.btagJetProb = pfjet->bDiscriminator("jetProbabilityBJetTags");
     MyPfjet.btagJetBProb = pfjet->bDiscriminator("jetBProbabilityBJetTags");
     MyPfjet.btagSoftEle = pfjet->bDiscriminator("softElectronBJetTags");
     MyPfjet.btagSoftMuon = pfjet->bDiscriminator("softMuonBJetTags");
     MyPfjet.btagSoftMuonNoIP = pfjet->bDiscriminator("softMuonNoIPBJetTags");
     MyPfjet.btagSecVertex = pfjet->bDiscriminator("simpleSecondaryVertexBJetTags");
     MyPfjet.btagSecVertexHighEff = pfjet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
     MyPfjet.btagSecVertexHighPur = pfjet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags");

     /*
     // MyPfjet.Upt = Rawjetpt;
     // MyPfjet.L2pt = L2jetpt;
     // MyPfjet.L2L3pt = L2L3jetpt;
     MyPfjet.respt = respt;
     */
     MyPfjet.JESunc = unc;

     if( (pfjet->genJet()) ){ // if there is a matched genjet, fill variables
       MyPfjet.genJetET = pfjet->genJet()->et();
       MyPfjet.genJetPT = pfjet->genJet()->pt();
       MyPfjet.genJetEta = pfjet->genJet()->eta();
       MyPfjet.genJetPhi = pfjet->genJet()->phi();
     }
     if( (pfjet->genParton()) ){ // if there is a matched parton, fill variables
       MyPfjet.genPartonET = pfjet->genParton()->et();
       MyPfjet.genPartonPT = pfjet->genParton()->pt();
       MyPfjet.genPartonEta = pfjet->genParton()->eta();
       MyPfjet.genPartonPhi = pfjet->genParton()->phi();
       MyPfjet.genPartonId = pfjet->genParton()->pdgId();
       MyPfjet.genPartonMotherId = pfjet->genParton()->mother()->pdgId();
       MyPfjet.genPartonGrandMotherId = pfjet->genParton()->mother()->pdgId();
     }


     // DataFormats/JetReco/interface/PFJet.h
     // http://indico.cern.ch/getFile.py/access?contribId=8&resId=0&materialId=slides&confId=92249
     MyPfjet.chargedHadronEnergyFraction = pfjet->chargedHadronEnergyFraction();
     MyPfjet.neutralHadronEnergyFraction = pfjet->neutralHadronEnergyFraction();
     MyPfjet.chargedEmEnergyFraction = pfjet->chargedEmEnergyFraction();
     MyPfjet.neutralEmEnergyFraction = pfjet->neutralEmEnergyFraction();
     MyPfjet.chargedMultiplicity = pfjet->chargedMultiplicity();
     MyPfjet.neutralMultiplicity = pfjet->neutralMultiplicity();
     MyPfjet.nconstituents = pfjet->numberOfDaughters();


     bool loose = (
		   pfjet->neutralHadronEnergyFraction() < 0.99 &&
		   pfjet->chargedEmEnergyFraction() < 0.99 &&
		   pfjet->neutralEmEnergyFraction() < 0.99 &&
		   pfjet->numberOfDaughters() > 1
		   );

     if( fabs(pfjet->eta())<2.4 ){
       loose = ( loose &&
		 pfjet->chargedHadronEnergyFraction() > 0.0 &&
		 pfjet->chargedMultiplicity() > 0
		 );
     }

     MyPfjet.jetIDLoose = (loose) ? 1 : 0;

     bnpfjets->push_back(MyPfjet);
   }




   //  Fill the muon collection
   std::auto_ptr<BNmuonCollection> bnmuons(new BNmuonCollection);
   for( edm::View<pat::Muon>::const_iterator muon = muons.begin(); muon!=muons.end(); ++muon ){

     BNmuon MyMuon;

     // general kinematic variables
     MyMuon.energy = muon->energy();
     MyMuon.et = muon->et();
     MyMuon.pt = muon->pt();
     MyMuon.px = muon->px();
     MyMuon.py = muon->py();
     MyMuon.pz = muon->pz();
     MyMuon.phi = muon->phi();
     MyMuon.eta = muon->eta();
     MyMuon.theta = muon->theta();


     MyMuon.charge = muon->charge();
     MyMuon.vx = muon->vx();
     MyMuon.vy = muon->vy();
     MyMuon.vz = muon->vz();

     MyMuon.trackIso = muon->trackIso();
     MyMuon.ecalIso = muon->ecalIso();
     MyMuon.hcalIso = muon->hcalIso();
     MyMuon.caloIso = muon->caloIso();

     MyMuon.trackIsoDR03 = muon->isolationR03().sumPt;
     MyMuon.ecalIsoDR03 = muon->isolationR03().emEt;
     MyMuon.hcalIsoDR03 = muon->isolationR03().hadEt;
     MyMuon.caloIsoDR03 = muon->isolationR03().emEt + muon->isolationR03().hadEt;

     MyMuon.trackVetoIsoDR03 = muon->isolationR03().trackerVetoPt;
     MyMuon.ecalVetoIsoDR03 = muon->isolationR03().emVetoEt;
     MyMuon.hcalVetoIsoDR03 = muon->isolationR03().hadVetoEt;
     MyMuon.caloVetoIsoDR03 = muon->isolationR03().emVetoEt + muon->isolationR03().hadVetoEt;

     MyMuon.trackIsoDR05 = muon->isolationR05().sumPt;
     MyMuon.ecalIsoDR05 = muon->isolationR05().emEt;
     MyMuon.hcalIsoDR05 = muon->isolationR05().hadEt;
     MyMuon.caloIsoDR05 = muon->isolationR05().emEt + muon->isolationR05().hadEt;

     MyMuon.trackVetoIsoDR05 = muon->isolationR05().trackerVetoPt;
     MyMuon.ecalVetoIsoDR05 = muon->isolationR05().emVetoEt;
     MyMuon.hcalVetoIsoDR05 = muon->isolationR05().hadVetoEt;
     MyMuon.caloVetoIsoDR05 = muon->isolationR05().emVetoEt + muon->isolationR05().hadVetoEt;

     ////   Removed due to not understood seg faults. Seems ok in bare ROOT.
     ////   Investigate and uncomment in future.
     //MyMuon.hcalE = muon->hcalIsoDeposit()->candEnergy();
     //MyMuon.ecalE = muon->ecalIsoDeposit()->candEnergy();


     MyMuon.IDGMPTight = ( muon->isGood("GlobalMuonPromptTight") ) ? 1 : 0;

     MyMuon.isGlobalMuon = ( muon->isGlobalMuon() ) ? 1 : 0;
     MyMuon.isTrackerMuon = ( muon->isTrackerMuon() ) ? 1 : 0;
     MyMuon.isStandAloneMuon = ( muon->isStandAloneMuon() ) ? 1 : 0;
     MyMuon.isGlobalMuonPromptTight = ( muon->isGood("GlobalMuonPromptTight") ) ? 1 : 0;

     MyMuon.numberOfMatches = muon->numberOfMatches();
     MyMuon.numberOfMatchedStations = muon->numberOfMatchedStations();

     MyMuon.dVzPVz = muon->vz() - PVz;
     MyMuon.dB = muon->dB();

     if( (muon->globalTrack().isAvailable()) ){
       MyMuon.normalizedChi2 = muon->globalTrack()->normalizedChi2();
       MyMuon.numberOfValidMuonHits = muon->globalTrack()->hitPattern().numberOfValidMuonHits();
       MyMuon.numberOfValidTrackerHits = muon->globalTrack()->hitPattern().numberOfValidTrackerHits();
       MyMuon.ptErr = muon->globalTrack()->ptError();
     }
     if( (muon->innerTrack().isAvailable()) ){
       MyMuon.numberOfValidTrackerHitsInnerTrack = muon->innerTrack()->numberOfValidHits();
       MyMuon.pixelLayersWithMeasurement = muon->innerTrack()->hitPattern().pixelLayersWithMeasurement();

       MyMuon.correctedD0 = muon->innerTrack()->dxy(beamSpotPosition);
       MyMuon.correctedDZ = muon->innerTrack()->dz(vertexPosition);
     }

     // Get track muon info
     if( (muon->track().isAvailable()) ){
       double tkvx = muon->track()->vx();
       double tkvy = muon->track()->vy();
       double tkpx = muon->track()->px();
       double tkpy = muon->track()->py();
       double tkpt = muon->track()->pt();

       double ndof = muon->track()->ndof();
       if( (ndof!=0) ) MyMuon.tkNormChi2 = muon->track()->chi2()/ndof;
       MyMuon.tkPT = muon->track()->pt();
       MyMuon.tkEta = muon->track()->eta();
       MyMuon.tkPhi = muon->track()->phi();
       MyMuon.tkDZ = muon->track()->dz();
       MyMuon.tkD0 = muon->track()->d0();
       MyMuon.tkD0bs = (tkvx-BSx)*tkpy/tkpt-(tkvy-BSy)*tkpx/tkpt;
       MyMuon.tkD0err = muon->track()->d0Error();
       MyMuon.tkNumValidHits = muon->track()->numberOfValidHits();
       MyMuon.tkCharge = muon->track()->charge();
     }
     // Get standalone muon info
     if( (muon->standAloneMuon().isAvailable()) ){
       double samvx = muon->standAloneMuon()->vx();
       double samvy = muon->standAloneMuon()->vy();
       double sampx = muon->standAloneMuon()->px();
       double sampy = muon->standAloneMuon()->py();
       double sampt = muon->standAloneMuon()->pt();

       double ndof = muon->standAloneMuon()->ndof();
       if( (ndof!=0) ) MyMuon.samNormChi2 = muon->standAloneMuon()->chi2()/ndof;
       MyMuon.samPT = muon->standAloneMuon()->pt();
       MyMuon.samEta = muon->standAloneMuon()->eta();
       MyMuon.samPhi = muon->standAloneMuon()->phi();
       MyMuon.samDZ = muon->standAloneMuon()->dz();
       MyMuon.samD0 = muon->standAloneMuon()->d0();
       MyMuon.samD0bs = (samvx-BSx)*sampy/sampt-(samvy-BSy)*sampx/sampt;
       MyMuon.samD0err = muon->standAloneMuon()->d0Error();
       MyMuon.samNumValidHits = muon->standAloneMuon()->numberOfValidHits();
       MyMuon.samCharge = muon->standAloneMuon()->charge();
     }
     // Get global muon info
     if( (muon->combinedMuon().isAvailable()) ){
       double comvx = muon->combinedMuon()->vx();
       double comvy = muon->combinedMuon()->vy();
       double compx = muon->combinedMuon()->px();
       double compy = muon->combinedMuon()->py();
       double compt = muon->combinedMuon()->pt();

       double ndof = muon->combinedMuon()->ndof();
       if( (ndof!=0) ) MyMuon.comNormChi2 = muon->combinedMuon()->chi2()/ndof;
       MyMuon.comPT = muon->combinedMuon()->pt();
       MyMuon.comEta = muon->combinedMuon()->eta();
       MyMuon.comPhi = muon->combinedMuon()->phi();
       MyMuon.comDZ = muon->combinedMuon()->dz();
       MyMuon.comD0 = muon->combinedMuon()->d0();
       MyMuon.comD0bs = (comvx-BSx)*compy/compt-(comvy-BSy)*compx/compt;
       MyMuon.comD0err = muon->combinedMuon()->d0Error();
       MyMuon.comNumValidHits = muon->combinedMuon()->numberOfValidHits();
       MyMuon.comCharge = muon->combinedMuon()->charge();
     }

     if( (muon->genLepton()) ){
       MyMuon.genId = muon->genLepton()->pdgId();
       MyMuon.genET = muon->genLepton()->et();
       MyMuon.genPT = muon->genLepton()->pt();
       MyMuon.genPhi = muon->genLepton()->phi();
       MyMuon.genEta = muon->genLepton()->eta();
       MyMuon.genCharge = muon->genLepton()->charge();

       if( (muon->genLepton()->numberOfMothers()>0) ){
	 const reco::Candidate *MuonMother = muon->genLepton()->mother();
	 bool staytrapped = true;
	 while( (abs(MuonMother->pdgId())==13 && staytrapped) ){
	   if( MuonMother->numberOfMothers()==1 ) MuonMother = MuonMother->mother();
	   else staytrapped = false;
	 }
       
	 MyMuon.genMotherId = MuonMother->pdgId();
	 MyMuon.genMotherET = MuonMother->et();
	 MyMuon.genMotherPT = MuonMother->pt();
	 MyMuon.genMotherPhi = MuonMother->phi();
	 MyMuon.genMotherEta = MuonMother->eta();
	 MyMuon.genMotherCharge = MuonMother->charge();
       }
     }

     bnmuons->push_back(MyMuon);
   }



   //  Fill the photon collection
   std::auto_ptr<BNphotonCollection> bnphotons(new BNphotonCollection);
   for( edm::View<pat::Photon>::const_iterator photon = photons.begin(); photon!=photons.end(); ++photon ){

     if( !(photon->et()>minPhotonEt_) ) continue;

     BNphoton MyPhoton;

     // general kinematic variables
     MyPhoton.energy = photon->energy();
     MyPhoton.et = photon->et();
     MyPhoton.pt = photon->pt();
     MyPhoton.px = photon->px();
     MyPhoton.py = photon->py();
     MyPhoton.pz = photon->pz();
     MyPhoton.phi = photon->phi();
     MyPhoton.eta = photon->eta();
     MyPhoton.theta = photon->theta();

     MyPhoton.trackIso = photon->trackIso();
     MyPhoton.ecalIso = photon->ecalIso();
     MyPhoton.hcalIso = photon->hcalIso();
     MyPhoton.caloIso = photon->caloIso();

     MyPhoton.trackIsoHollowConeDR03 = photon->trkSumPtHollowConeDR03();
     MyPhoton.trackIsoSolidConeDR03 = photon->trkSumPtSolidConeDR03();
     MyPhoton.ecalIsoDR03 = photon->ecalRecHitSumEtConeDR03();
     MyPhoton.hcalIsoDR03 = photon->hcalTowerSumEtConeDR03();
     MyPhoton.caloIsoDR03 = photon->ecalRecHitSumEtConeDR03() + photon->hcalTowerSumEtConeDR03();

     MyPhoton.trackIsoHollowConeDR04 = photon->trkSumPtHollowConeDR04();
     MyPhoton.trackIsoSolidConeDR04 = photon->trkSumPtSolidConeDR04();
     MyPhoton.ecalIsoDR04 = photon->ecalRecHitSumEtConeDR04();
     MyPhoton.hcalIsoDR04 = photon->hcalTowerSumEtConeDR04();
     MyPhoton.caloIsoDR04 = photon->ecalRecHitSumEtConeDR04() + photon->hcalTowerSumEtConeDR04();

     MyPhoton.hadOverEm = photon->hadronicOverEm();
     MyPhoton.sigmaEtaEta = photon->sigmaEtaEta();
     MyPhoton.sigmaIetaIeta = photon->sigmaIetaIeta();
     MyPhoton.r9 = photon->r9();

     MyPhoton.IDLooseEM = -1;
     MyPhoton.IDLoose = -1;
     MyPhoton.IDTight = -1;

     MyPhoton.hasPixelSeed = ( photon->hasPixelSeed() ) ? 1 : 0;

     // Get photon supercluster information
     if( (photon->superCluster().isAvailable()) ){
       double eMax    = lazyTools.eMax(    *(photon->superCluster()) );
       double eLeft   = lazyTools.eLeft(   *(photon->superCluster()) );
       double eRight  = lazyTools.eRight(  *(photon->superCluster()) );
       double eTop    = lazyTools.eTop(    *(photon->superCluster()) );
       double eBottom = lazyTools.eBottom( *(photon->superCluster()) );
       double e3x3    = lazyTools.e3x3(    *(photon->superCluster()) );
       double swissCross = -99;

       if( eMax>0.000001 ) swissCross = 1 - (eLeft+eRight+eTop+eBottom)/eMax;

       MyPhoton.eMax = eMax;
       MyPhoton.eLeft = eLeft;
       MyPhoton.eRight = eRight;
       MyPhoton.eTop = eTop;
       MyPhoton.eBottom = eBottom;
       MyPhoton.e3x3 = e3x3;
       MyPhoton.swissCross = swissCross;

       MyPhoton.scEnergy = photon->superCluster()->energy();
       MyPhoton.scRawEnergy = photon->superCluster()->rawEnergy();
       MyPhoton.scEta = photon->superCluster()->position().eta();
       MyPhoton.scPhi = photon->superCluster()->position().phi();
       MyPhoton.scZ = photon->superCluster()->position().Z();

       double seedE = -999, seedTime = -999;
       int seedRecoFlag = -999;
       
       if( (photon->isEB()) ){
       	 DetId idEB = EcalClusterTools::getMaximum( photon->superCluster()->hitsAndFractions(), &(*barrelRecHits) ).first;
       	 EcalRecHitCollection::const_iterator thisHitEB = barrelRecHits->find(idEB);

       	 seedE = thisHitEB->energy();
       	 seedTime = thisHitEB->time();
       	 seedRecoFlag = thisHitEB->recoFlag();
       }
       else if( (photon->isEE()) ){
       	 DetId idEE = EcalClusterTools::getMaximum( photon->superCluster()->hitsAndFractions(), &(*endcapRecHits) ).first;
       	 EcalRecHitCollection::const_iterator thisHitEE = endcapRecHits->find(idEE);

       	 seedE = thisHitEE->energy();
       	 seedTime = thisHitEE->time();
       	 seedRecoFlag = thisHitEE->recoFlag();
       }

       MyPhoton.seedEnergy   = seedE;
       MyPhoton.seedTime     = seedTime;
       MyPhoton.seedRecoFlag = seedRecoFlag;

       // MyPhoton.swissCrossNoI85 = swissCrossNoI85;
       // MyPhoton.swissCrossI85   = swissCrossI85;
       // MyPhoton.E2overE9NoI85 = e2overe9NoI85;
       // MyPhoton.E2overE9I85 = e2overe9I85;
     }

     MyPhoton.isEB = photon->isEB();
     MyPhoton.isEE = photon->isEE();
     MyPhoton.isEBEEGap = photon->isEBEEGap();
     MyPhoton.isEBGap = photon->isEBGap();
     MyPhoton.isEEGap = photon->isEEGap();
     MyPhoton.isGap = ( (photon->isEBEEGap()) || (photon->isEBGap()) || (photon->isEEGap()) );

     if( (photon->genPhoton()) ){
       MyPhoton.genId = photon->genPhoton()->pdgId();
       MyPhoton.genET = photon->genPhoton()->et();
       MyPhoton.genPT = photon->genPhoton()->pt();
       MyPhoton.genPhi = photon->genPhoton()->phi();
       MyPhoton.genEta = photon->genPhoton()->eta();
       MyPhoton.genCharge = photon->genPhoton()->charge();
       MyPhoton.genMotherId = photon->genPhoton()->mother()->pdgId();
       MyPhoton.genMotherET = photon->genPhoton()->mother()->et();
       MyPhoton.genMotherPT = photon->genPhoton()->mother()->pt();
       MyPhoton.genMotherPhi = photon->genPhoton()->mother()->phi();
       MyPhoton.genMotherEta = photon->genPhoton()->mother()->eta();
       MyPhoton.genMotherCharge = photon->genPhoton()->mother()->charge();
     }

     bnphotons->push_back(MyPhoton);
   }


   int numhighpurity=0;
   reco::TrackBase::TrackQuality _trackQuality = reco::TrackBase::qualityByName("highPurity");

   //  Fill the track collection
   std::auto_ptr<BNtrackCollection> bntracks(new BNtrackCollection);
   for(reco::TrackCollection::const_iterator track = tracks.begin(); track!=tracks.end(); ++track ){

     if( track->quality(_trackQuality) ) numhighpurity++;

     if( !(track->pt()>minTrackPt_) ) continue;

     BNtrack MyTrack;

     MyTrack.pt = track->pt();
     MyTrack.px = track->px();
     MyTrack.py = track->py();
     MyTrack.pz = track->pz();
     MyTrack.phi = track->phi();
     MyTrack.eta = track->eta();
     MyTrack.theta = track->theta();

     MyTrack.vx = track->vx();
     MyTrack.vy = track->vy();
     MyTrack.vz = track->vz();
     MyTrack.charge = track->charge();
     MyTrack.numValidHits = track->numberOfValidHits();

     MyTrack.dZ = track->dz();
     MyTrack.d0 = track->d0();
     MyTrack.d0err = track->d0Error();
     double ndof = track->ndof();
     if( (ndof!=0) ) MyTrack.normChi2 = track->chi2()/ndof;


     bntracks->push_back(MyTrack);
   }


   bool FilterOutScraping = false;
   int tracksSize = int(tracks.size());
   double FilterOutScrapingFraction = ( tracksSize!=0 ) ? (double)numhighpurity/(double)tracksSize : 0;

   if( tracksSize>numtrack_ ){
     if( FilterOutScrapingFraction>thresh_ ) FilterOutScraping=true;
   }
   else FilterOutScraping = true;


   //  Fill the supercluster collection
   std::auto_ptr<BNsuperclusterCollection> bnsuperclusters(new BNsuperclusterCollection);
   for(reco::SuperClusterCollection::const_iterator supercluster = EBsuperclusters.begin(); supercluster!=EBsuperclusters.end(); ++supercluster ){

     double sctheta = 2.0*atan(exp(-supercluster->eta()));
     double scet = supercluster->energy()*sin(sctheta);

     if( !(scet>minSCEt_) ) continue;

     BNsupercluster MySC;

     MySC.energy = supercluster->energy();
     MySC.phi = supercluster->phi();
     MySC.eta = supercluster->eta();

     MySC.theta = sctheta;

     MySC.et = supercluster->energy()*sin(sctheta);
     MySC.ex = supercluster->energy()*sin(sctheta)*cos(supercluster->phi());
     MySC.ey = supercluster->energy()*sin(sctheta)*sin(supercluster->phi());
     MySC.ez = supercluster->energy()*cos(sctheta);


     bnsuperclusters->push_back(MySC);
   }
   for(reco::SuperClusterCollection::const_iterator supercluster = EEsuperclusters.begin(); supercluster!=EEsuperclusters.end(); ++supercluster ){

     double sctheta = 2.0*atan(exp(-supercluster->eta()));
     double scet = supercluster->energy()*sin(sctheta);

     if( !(scet>minSCEt_) ) continue;

     BNsupercluster MySC;

     MySC.energy = supercluster->energy();
     MySC.phi = supercluster->phi();
     MySC.eta = supercluster->eta();

     MySC.theta = sctheta;

     MySC.et = supercluster->energy()*sin(sctheta);
     MySC.ex = supercluster->energy()*sin(sctheta)*cos(supercluster->phi());
     MySC.ey = supercluster->energy()*sin(sctheta)*sin(supercluster->phi());
     MySC.ez = supercluster->energy()*cos(sctheta);


     bnsuperclusters->push_back(MySC);
   }



   //  Fill the mc particles collections
   std::auto_ptr<BNmcparticleCollection> bnmcparticles(new BNmcparticleCollection);
   std::auto_ptr<BNmcparticleCollection> bnmcelectrons(new BNmcparticleCollection);
   std::auto_ptr<BNmcparticleCollection> bnmcmuons(new BNmcparticleCollection);
   std::vector<int> vecZdecay;
   std::vector<int> vecWdecay;
   vecZdecay.clear();
   vecWdecay.clear();

   if( (genParticles.isValid() && sample_>=0) ){
     for( size_t k = 0; k < genParticles->size(); k++ ){
       const reco::Candidate & mcParticle = (*genParticles)[k];

       int status = mcParticle.status();
       int pdgId  = mcParticle.pdgId();
       int numdgt = mcParticle.numberOfDaughters();

       BNmcparticle MyMCelectron;
       BNmcparticle MyMCmuon;

       if( (status==1) ){
	 if( fabs(pdgId)==11 ){
	   MyMCelectron.energy = mcParticle.energy();
	   MyMCelectron.et = mcParticle.et();
	   MyMCelectron.pt = mcParticle.pt();
	   MyMCelectron.px = mcParticle.px();
	   MyMCelectron.py = mcParticle.py();
	   MyMCelectron.pz = mcParticle.pz();
	   MyMCelectron.phi = mcParticle.phi();
	   MyMCelectron.eta = mcParticle.eta();
	   MyMCelectron.theta = mcParticle.theta();
	   MyMCelectron.mass = mcParticle.mass();
	   MyMCelectron.vx = mcParticle.vx();
	   MyMCelectron.vy = mcParticle.vy();
	   MyMCelectron.vz = mcParticle.vz();
	   MyMCelectron.charge = mcParticle.charge();
	   MyMCelectron.id = mcParticle.pdgId();
	   MyMCelectron.status = mcParticle.status();

	   if( (mcParticle.numberOfMothers()==1) ){
	     MyMCelectron.motherId = mcParticle.mother()->pdgId();
	     MyMCelectron.motherCharge = mcParticle.mother()->charge();
	     MyMCelectron.motherET = mcParticle.mother()->et();
	     MyMCelectron.motherPT = mcParticle.mother()->pt();
	     MyMCelectron.motherPhi = mcParticle.mother()->phi();
	     MyMCelectron.motherEta = mcParticle.mother()->eta();

	     if( (mcParticle.mother()->numberOfMothers()==1) ){
	       MyMCelectron.grandMotherId = mcParticle.mother()->mother()->pdgId();
	       MyMCelectron.grandMotherCharge = mcParticle.mother()->mother()->charge();
	       MyMCelectron.grandMotherET = mcParticle.mother()->mother()->et();
	       MyMCelectron.grandMotherPT = mcParticle.mother()->mother()->pt();
	       MyMCelectron.grandMotherPhi = mcParticle.mother()->mother()->phi();
	       MyMCelectron.grandMotherEta = mcParticle.mother()->mother()->eta();
	     }
	   }

	   bnmcelectrons->push_back(MyMCelectron);
	 }
	 else if( fabs(pdgId)==13 ){
	   MyMCmuon.energy = mcParticle.energy();
	   MyMCmuon.et = mcParticle.et();
	   MyMCmuon.pt = mcParticle.pt();
	   MyMCmuon.px = mcParticle.px();
	   MyMCmuon.py = mcParticle.py();
	   MyMCmuon.pz = mcParticle.pz();
	   MyMCmuon.phi = mcParticle.phi();
	   MyMCmuon.eta = mcParticle.eta();
	   MyMCmuon.theta = mcParticle.theta();
	   MyMCmuon.mass = mcParticle.mass();
	   MyMCmuon.vx = mcParticle.vx();
	   MyMCmuon.vy = mcParticle.vy();
	   MyMCmuon.vz = mcParticle.vz();
	   MyMCmuon.charge = mcParticle.charge();
	   MyMCmuon.id = mcParticle.pdgId();
	   MyMCmuon.status = mcParticle.status();

	   if( (mcParticle.numberOfMothers()==1) ){
	     MyMCmuon.motherId = mcParticle.mother()->pdgId();
	     MyMCmuon.motherCharge = mcParticle.mother()->charge();
	     MyMCmuon.motherET = mcParticle.mother()->et();
	     MyMCmuon.motherPT = mcParticle.mother()->pt();
	     MyMCmuon.motherPhi = mcParticle.mother()->phi();
	     MyMCmuon.motherEta = mcParticle.mother()->eta();

	     if( (mcParticle.mother()->numberOfMothers()==1) ){
	       MyMCmuon.grandMotherId = mcParticle.mother()->mother()->pdgId();
	       MyMCmuon.grandMotherCharge = mcParticle.mother()->mother()->charge();
	       MyMCmuon.grandMotherET = mcParticle.mother()->mother()->et();
	       MyMCmuon.grandMotherPT = mcParticle.mother()->mother()->pt();
	       MyMCmuon.grandMotherPhi = mcParticle.mother()->mother()->phi();
	       MyMCmuon.grandMotherEta = mcParticle.mother()->mother()->eta();
	     }
	   }

	   bnmcmuons->push_back(MyMCmuon);
	 }
       }

       if( (status!=3) ) continue;

       BNmcparticle MyMCparticle;

       MyMCparticle.energy = mcParticle.energy();
       MyMCparticle.et = mcParticle.et();
       MyMCparticle.pt = mcParticle.pt();
       MyMCparticle.px = mcParticle.px();
       MyMCparticle.py = mcParticle.py();
       MyMCparticle.pz = mcParticle.pz();
       MyMCparticle.phi = mcParticle.phi();
       MyMCparticle.eta = mcParticle.eta();
       MyMCparticle.theta = mcParticle.theta();
       MyMCparticle.mass = mcParticle.mass();
       MyMCparticle.vx = mcParticle.vx();
       MyMCparticle.vy = mcParticle.vy();
       MyMCparticle.vz = mcParticle.vz();
       MyMCparticle.charge = mcParticle.charge();
       MyMCparticle.id = mcParticle.pdgId();
       MyMCparticle.status = mcParticle.status();

       if( (mcParticle.numberOfMothers()==1) ){
	 MyMCparticle.motherId = mcParticle.mother()->pdgId();
	 MyMCparticle.motherCharge = mcParticle.mother()->charge();
	 MyMCparticle.motherET = mcParticle.mother()->et();
	 MyMCparticle.motherPT = mcParticle.mother()->pt();
	 MyMCparticle.motherPhi = mcParticle.mother()->phi();
	 MyMCparticle.motherEta = mcParticle.mother()->eta();

	 if( (mcParticle.mother()->numberOfMothers()==1) ){
	   MyMCparticle.grandMotherId = mcParticle.mother()->mother()->pdgId();
	   MyMCparticle.grandMotherCharge = mcParticle.mother()->mother()->charge();
	   MyMCparticle.grandMotherET = mcParticle.mother()->mother()->et();
	   MyMCparticle.grandMotherPT = mcParticle.mother()->mother()->pt();
	   MyMCparticle.grandMotherPhi = mcParticle.mother()->mother()->phi();
	   MyMCparticle.grandMotherEta = mcParticle.mother()->mother()->eta();
	 }
       }


       bnmcparticles->push_back(MyMCparticle);


       if( (abs(pdgId)==23 && numdgt>=2 && vecZdecay.size()<=2) ){
	 int zdecay = abs(mcParticle.daughter(0)->pdgId());

	 vecZdecay.push_back(zdecay);
       }
       else if( (abs(pdgId)==24 && numdgt>=2 && vecWdecay.size()<=2) ){
	 int dId0 = abs(mcParticle.daughter(0)->pdgId());
	 int dId1 = abs(mcParticle.daughter(1)->pdgId());

	 int wdecay0 = 100*dId0 + dId1;
	 int wdecay1 = 100*dId1 + dId0;

	 int wdecay = ( wdecay0 < wdecay1 ) ? wdecay0 : wdecay1;

	 vecWdecay.push_back(wdecay);
       }

     }
   }


   std::auto_ptr<BNeventCollection> bnevent(new BNeventCollection);
   BNevent MyEvent;


   MyEvent.weight = weight;
   MyEvent.pthat = pthat;
   MyEvent.qScale = qScale;
   MyEvent.alphaQCD = alphaQCD;
   MyEvent.alphaQED = alphaQED;
   MyEvent.scalePDF = scalePDF;
   MyEvent.x1 = x1;
   MyEvent.x2 = x2;
   MyEvent.xPDF1 = xPDF1;
   MyEvent.xPDF2 = xPDF2;
   MyEvent.BSx = BSx;
   MyEvent.BSy = BSy;
   MyEvent.BSz = BSz;
   MyEvent.run = iEvent.id().run();
   MyEvent.evt = iEvent.id().event();
   MyEvent.lumi = iEvent.luminosityBlock();
   MyEvent.instLumi = instLumi;
   MyEvent.bxLumi   = bxLumi;

   MyEvent.GoodVertex = ( GoodVertex ) ? 1 : 0;
   MyEvent.FilterOutScraping = ( FilterOutScraping ) ? 1 : 0;
   MyEvent.FilterOutScrapingFraction = FilterOutScrapingFraction;

   MyEvent.sample = sample_;
   MyEvent.numPV = numPVs;
   MyEvent.hcalnoiseLoose = ( passLooseNoiseFilter ) ? 1 : 0;
   MyEvent.hcalnoiseTight = ( passTightNoiseFilter ) ? 1 : 0;

   MyEvent.HBHENoiseFilter = ( hcalNoiseFilter ) ? 1 : 0;

   MyEvent.CSCLooseHaloId = ( passCSCLooseHaloId ) ? 1 : 0;
   MyEvent.CSCTightHaloId = ( passCSCTightHaloId ) ? 1 : 0;
   MyEvent.EcalLooseHaloId = ( passEcalLooseHaloId ) ? 1 : 0;
   MyEvent.EcalTightHaloId = ( passEcalTightHaloId ) ? 1 : 0;
   MyEvent.HcalLooseHaloId = ( passHcalLooseHaloId ) ? 1 : 0;
   MyEvent.HcalTightHaloId = ( passHcalTightHaloId ) ? 1 : 0;
   MyEvent.GlobalLooseHaloId = ( passGlobalLooseHaloId ) ? 1 : 0;
   MyEvent.GlobalTightHaloId = ( passGlobalTightHaloId ) ? 1 : 0;
   MyEvent.LooseId = ( passLooseId ) ? 1 : 0;
   MyEvent.TightId = ( passTightId ) ? 1 : 0;


   MyEvent.bField = evt_bField;

   if( (vecWdecay.size()>0) ){
     MyEvent.W0decay = vecWdecay[0];
     if( (vecWdecay.size()>1) ){
       MyEvent.W1decay = vecWdecay[1];
     }
   }
   if( (vecZdecay.size()>0) ){
     MyEvent.Z0decay = vecZdecay[0];
     if( (vecZdecay.size()>1) ){
       MyEvent.Z1decay = vecZdecay[1];
     }
   }


   bnevent->push_back(MyEvent);


   // Calo MET
   std::auto_ptr<BNmetCollection> bncalomet(new BNmetCollection);
   BNmet MyCalomet;

   //double calometPT = calometHandle->front().pt();
   MyCalomet.pt = calometHandle->front().pt();
   MyCalomet.px = calometHandle->front().px();
   MyCalomet.py = calometHandle->front().py();
   MyCalomet.phi = calometHandle->front().phi();
   MyCalomet.sumET = calometHandle->front().sumEt();
   MyCalomet.corSumET = calometHandle->front().corSumEt();
   MyCalomet.Upt = calometHandle->front().uncorrectedPt();
   MyCalomet.Uphi = calometHandle->front().uncorrectedPhi();

   MyCalomet.et = calometHandle->front().et();
   MyCalomet.mEtSig = calometHandle->front().mEtSig();
   MyCalomet.metSignificance = calometHandle->front().metSignificance();

   MyCalomet.maxEtInEmTowers = calometHandle->front().maxEtInEmTowers();
   MyCalomet.emEtFraction = calometHandle->front().emEtFraction();
   MyCalomet.emEtInEB = calometHandle->front().emEtInEB();
   MyCalomet.emEtInEE = calometHandle->front().emEtInEE();
   MyCalomet.emEtInHF = calometHandle->front().emEtInHF();
   MyCalomet.maxEtInHadTowers = calometHandle->front().maxEtInHadTowers();
   MyCalomet.hadEtFraction = calometHandle->front().etFractionHadronic();
   MyCalomet.hadEtInHB = calometHandle->front().hadEtInHB();
   MyCalomet.hadEtInHE = calometHandle->front().hadEtInHE();
   MyCalomet.hadEtInHF = calometHandle->front().hadEtInHF();
   MyCalomet.hadEtInHO = calometHandle->front().hadEtInHO();

   if( (calometHandle->front().genMET()) ){
     MyCalomet.genPT = calometHandle->front().genMET()->pt();
     MyCalomet.genPhi = calometHandle->front().genMET()->phi();
   }


   if( calometHandle->front().nCorrections()==0 ){
     std::cout << "  nCorrections()==0, NO CORRECTIONS HAVE BEEN APPLIED!!" << std::endl;
   }

   // For now, simply trusting the muon corrections to MET
   double calometMuoncorrPx = 0;
   double calometMuoncorrPy = 0;
   if( calometHandle->front().nCorrections()==1 ){
     calometMuoncorrPx = calometHandle->front().corEx(pat::MET::uncorrJES);
     calometMuoncorrPy = calometHandle->front().corEy(pat::MET::uncorrJES);
   }
   else if( calometHandle->front().nCorrections()>1 ){
     calometMuoncorrPx = calometHandle->front().corEx(pat::MET::uncorrMUON);
     calometMuoncorrPy = calometHandle->front().corEy(pat::MET::uncorrMUON);
   }

   MyCalomet.muonCorEx  = calometMuoncorrPx;
   MyCalomet.muonCorEy  = calometMuoncorrPy;
   MyCalomet.jet20CorEx = metJES20corrPx;
   MyCalomet.jet20CorEy = metJES20corrPy;
   MyCalomet.jet1CorEx  = metJES1corrPx;
   MyCalomet.jet1CorEy  = metJES1corrPy;


   UDeltaPx += calometHandle->front().uncorrectedPt() * cos(calometHandle->front().uncorrectedPhi()) ;
   UDeltaPy += calometHandle->front().uncorrectedPt() * sin(calometHandle->front().uncorrectedPhi()) ;

   double UDeltaP = sqrt( UDeltaPx*UDeltaPx + UDeltaPy*UDeltaPy );

   double UscaleA = 1.5;
   double UscaleB = 1.8;
   double UscaleC = -0.06;

   double Uscale = UscaleA+UscaleB*exp(UscaleC*UDeltaP);

   double type2corPx = (Uscale-1.)*UDeltaPx;
   double type2corPy = (Uscale-1.)*UDeltaPy;
   double type2corSumET = (Uscale-1.)*USumET;

   double T2px = calometHandle->front().px() + type2corPx;
   double T2py = calometHandle->front().py() + type2corPy;
   double T2pt = sqrt( T2px*T2px + T2py*T2py );
   double T2phi = atan2(T2py,T2px);
   double T2sumET = calometHandle->front().sumEt() + type2corSumET;

   MyCalomet.UDeltaPx = UDeltaPx;
   MyCalomet.UDeltaPy = UDeltaPy;
   MyCalomet.UDeltaP = UDeltaP;
   MyCalomet.Uscale = Uscale;
   MyCalomet.type2corPx = type2corPx;
   MyCalomet.type2corPy = type2corPy;
   MyCalomet.T2pt = T2pt;
   MyCalomet.T2px = T2px;
   MyCalomet.T2py = T2py;
   MyCalomet.T2phi = T2phi;
   MyCalomet.T2sumET = T2sumET;

   bncalomet->push_back(MyCalomet);


   // Pf MET
   std::auto_ptr<BNmetCollection> bnpfmet(new BNmetCollection);
   BNmet MyPfmet;

   MyPfmet.pt = pfmetHandle->front().pt();
   MyPfmet.px = pfmetHandle->front().px();
   MyPfmet.py = pfmetHandle->front().py();
   MyPfmet.phi = pfmetHandle->front().phi();
   MyPfmet.sumET = pfmetHandle->front().sumEt();
   MyPfmet.corSumET = pfmetHandle->front().corSumEt();
   MyPfmet.Upt = pfmetHandle->front().uncorrectedPt();
   MyPfmet.Uphi = pfmetHandle->front().uncorrectedPhi();

   if( (pfmetHandle->front().genMET()) ){
     MyPfmet.genPT = pfmetHandle->front().genMET()->pt();
     MyPfmet.genPhi = pfmetHandle->front().genMET()->phi();
   }

   /*
   double pfMETjet1corrPx = pfmetHandle->front().px() + metJet1corrPx_pf;
   double pfMETjet1corrPy = pfmetHandle->front().py() + metJet1corrPy_pf;
   double pfMETjet6corrPx = pfmetHandle->front().px() + metJet6corrPx_pf;
   double pfMETjet6corrPy = pfmetHandle->front().py() + metJet6corrPy_pf;
   double pfMETjet10corrPx = pfmetHandle->front().px() + metJet10corrPx_pf;
   double pfMETjet10corrPy = pfmetHandle->front().py() + metJet10corrPy_pf;

   double pfMETjet1corrPt  = sqrt( pfMETjet1corrPx*pfMETjet1corrPx   + pfMETjet1corrPy*pfMETjet1corrPy );
   double pfMETjet6corrPt  = sqrt( pfMETjet6corrPx*pfMETjet6corrPx   + pfMETjet6corrPy*pfMETjet6corrPy );
   double pfMETjet10corrPt = sqrt( pfMETjet10corrPx*pfMETjet10corrPx + pfMETjet10corrPy*pfMETjet10corrPy );

   double pfMETjet1corrPhi = atan2( pfMETjet1corrPy,pfMETjet1corrPx);
   double pfMETjet6corrPhi = atan2( pfMETjet6corrPy,pfMETjet6corrPx);
   double pfMETjet10corrPhi = atan2( pfMETjet10corrPy,pfMETjet10corrPx);

   MyPfmet.pfT1jet1pt = pfMETjet1corrPt;
   MyPfmet.pfT1jet6pt = pfMETjet6corrPt;
   MyPfmet.pfT1jet10pt = pfMETjet10corrPt;
   MyPfmet.pfT1jet1phi = pfMETjet1corrPhi;
   MyPfmet.pfT1jet6phi = pfMETjet6corrPhi;
   MyPfmet.pfT1jet10phi = pfMETjet10corrPhi;
   */

   bnpfmet->push_back(MyPfmet);


   // Tc MET
   std::auto_ptr<BNmetCollection> bntcmet(new BNmetCollection);
   BNmet MyTcmet;

   MyTcmet.pt = tcmetHandle->front().pt();
   MyTcmet.px = tcmetHandle->front().px();
   MyTcmet.py = tcmetHandle->front().py();
   MyTcmet.phi = tcmetHandle->front().phi();
   MyTcmet.sumET = tcmetHandle->front().sumEt();
   MyTcmet.corSumET = tcmetHandle->front().corSumEt();
   MyTcmet.Upt = tcmetHandle->front().uncorrectedPt();
   MyTcmet.Uphi = tcmetHandle->front().uncorrectedPhi();

   if( (tcmetHandle->front().genMET()) ){
     MyTcmet.genPT = tcmetHandle->front().genMET()->pt();
     MyTcmet.genPhi = tcmetHandle->front().genMET()->phi();
   }


   bntcmet->push_back(MyTcmet);



   /// HLTrigger/Configuration/python/HLT_8E29_cff.py
   //  hltL1NonIsoHLTNonIsoSingleElectronLWEt15PixelMatchFilter
   edm::Handle< trigger::TriggerEvent > hltHandle;
   iEvent.getByLabel(triggerSummaryTag_, hltHandle);

   std::auto_ptr<BNtrigobjCollection> bntrigobjs(new BNtrigobjCollection);

   std::vector<trigger::size_type> filtKeys;
   std::vector<std::string> filtKeys_string;
   for(unsigned int i=0; i<hltHandle->sizeFilters(); i++) {
     const edm::InputTag filterTag = hltHandle->filterTag(i);
     const std::string filt = filterTag.encode();

     // if( filt.find("Electron")!=std::string::npos ){
     //   std::cout << "  filt string = " << filt << std::endl;
     // }

     //if( filt == electronTriggerFilter_.encode() ) {
     if( ( (filt.find("Electron")!=std::string::npos) ||
	   (filt.find("hltL1IsoL1sMu")!=std::string::npos) ||
	   (filt.find("hltL2IsoL1sMu")!=std::string::npos) ||
	   (filt.find("hltL3IsoL1sMu")!=std::string::npos) ||
	   (filt.find("hltL1fL1sMu")!=std::string::npos) ||
	   (filt.find("hltL2fL1sMu")!=std::string::npos) ||
	   (filt.find("hltL3fL1sMu")!=std::string::npos) ||
	   (filt.find("hltEle32WP70PFMT50PFMTFilter")!=std::string::npos) ||
	   (filt.find("HLTEle65CaloIdVTTrkIdTSequence")!=std::string::npos) )
	 ){
       const trigger::size_type filtIndex = hltHandle->filterIndex(filterTag);
       const std::vector<trigger::size_type>& theseKeys = hltHandle->filterKeys(filtIndex);
       int numKeys = theseKeys.size();
       for(int i=0; i<numKeys; i++) {
	 filtKeys.push_back(theseKeys.at(i));
	 filtKeys_string.push_back(filt);
       }
     }
   }

   const trigger::TriggerObjectCollection& hltObjects = hltHandle->getObjects();
   for(size_t i=0; i<filtKeys.size(); i++) {
     size_t key = filtKeys.at(i);
     const trigger::TriggerObject obj = hltObjects[key];

     BNtrigobj MyTrigobj;

     MyTrigobj.pt  = obj.pt();
     MyTrigobj.eta = obj.eta();
     MyTrigobj.phi = obj.phi();
     MyTrigobj.filter = filtKeys_string[i];

     bntrigobjs->push_back(MyTrigobj);
   }



   //  Fill the skim bits collection
   edm::Handle<edm::TriggerResults> hltresults;
   iEvent.getByLabel(triggerResultsTag_,hltresults);

   std::auto_ptr<BNtriggerCollection> bntrigger(new BNtriggerCollection);

   if(hltresults.isValid()){

     const edm::TriggerNames & trigNames = iEvent.triggerNames(*hltresults);
     unsigned int numTriggers = trigNames.size();

     for( unsigned int hltIndex=0; hltIndex<numTriggers-1; ++hltIndex ){

       std::string currentTrigName = trigNames.triggerName(hltIndex);
       int accept = hltresults->accept(hltIndex);

       int prescale = hltConfig_.prescaleValue(iEvent, iSetup, currentTrigName);

       BNtrigger MyTrigger;

       MyTrigger.pass = accept;
       MyTrigger.prescale = prescale;
       MyTrigger.name = currentTrigName;

       //std::cout << " =====>  name = " << currentTrigName << "\t prescale = " << prescale << "\t pass = " << accept << std::endl; 

       bntrigger->push_back(MyTrigger);

     }
   }


   // open main GT (DAQ) readout record - exit if failed
   Handle<L1GlobalTriggerReadoutRecord> gtReadoutRecord;
   iEvent.getByLabel(gtSource_, gtReadoutRecord);

   std::auto_ptr<BNtriggerCollection> bntriggerl1talgo(new BNtriggerCollection);
   std::auto_ptr<BNtriggerCollection> bntriggerl1ttech(new BNtriggerCollection);

   if (!gtReadoutRecord.isValid()) {
     std::cout << "  error ====> can't find L1GlobalTriggerReadoutRecord with label " << gtSource_.label();
     return;
   }

   // initialize bx's to invalid value
   int gtfeBx = -1;

   // get info from GTFE DAQ record
   const L1GtfeWord& gtfeWord = gtReadoutRecord->gtfeWord();
   gtfeBx = gtfeWord.bxNr();
   int gtfeActiveBoards = gtfeWord.activeBoards();

   // get info from FDL if active (including decision word)
   if( isActive(gtfeActiveBoards,FDL) ) {
     /// get Global Trigger algo and technical triger bit statistics
     const DecisionWord& gtDecisionWord = gtReadoutRecord->decisionWord();
     const TechnicalTriggerWord& gtTTWord = gtReadoutRecord->technicalTriggerWord();

     // L1 algos
     for( CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
       BNtrigger MyTriggerL1TAlgo;
       int algoBitNumber = (algo->second).algoBitNumber();

       MyTriggerL1TAlgo.pass = gtDecisionWord.at(algoBitNumber);
       MyTriggerL1TAlgo.name = (algo->second).algoName();

       bntriggerl1talgo->push_back(MyTriggerL1TAlgo);
     }

     // L1 technical
     char trigname[100];
     int tbitNumber = 0;
     TechnicalTriggerWord::const_iterator GTtbitItr;
     for(GTtbitItr = gtTTWord.begin(); GTtbitItr != gtTTWord.end(); GTtbitItr++) {
       BNtrigger MyTriggerL1TTech;
       int pass_l1t_tech = 0;
       if (*GTtbitItr) pass_l1t_tech = 1;

       sprintf(trigname, "L1_TechBit_%03d", tbitNumber);

       MyTriggerL1TTech.pass = pass_l1t_tech;
       MyTriggerL1TTech.name = trigname;

       bntriggerl1ttech->push_back(MyTriggerL1TTech);

       tbitNumber++; 
     }
   }


   if( GoodVertex ) numGoodVertexEvents++;
   if( FilterOutScraping ) numFilterOutScrapingEvents++;

   // Put the collections into the event
   iEvent.put(bnevent);
   iEvent.put(bnelectrons,eleTag_.label());
   iEvent.put(bncalojets,calojetTag_.label());
   iEvent.put(bnpfjets,pfjetTag_.label());
   iEvent.put(bncalomet,calometTag_.label());
   iEvent.put(bnpfmet,pfmetTag_.label());
   iEvent.put(bntcmet,tcmetTag_.label());
   iEvent.put(bnmuons,muonTag_.label());
   iEvent.put(bnphotons,photonTag_.label());
   iEvent.put(bnsuperclusters,kSC);
   iEvent.put(bntracks,trackTag_.label());
   iEvent.put(bntrigger,kTrigger);
   iEvent.put(bnmcparticles,kMCpar);
   iEvent.put(bnmcelectrons,kMCele);
   iEvent.put(bnmcmuons,kMCmu);
   iEvent.put(bntriggerl1talgo,kTriggerL1Talgo);
   iEvent.put(bntriggerl1ttech,kTriggerL1Ttech);
   iEvent.put(bntrigobjs);
   iEvent.put(bnpvs,pvTag_.label());

   iEvent.put(bnbxlumis,kBXlumi);

   /*
   vParam_L2_calo.clear(); vParam_L2L3_calo.clear(); vParam_L2L3res_calo.clear();
   delete L2JetCorPar_calo;
   delete L3JetCorPar_calo;
   delete ResJetCorPar_calo;
   delete JEC_L2_calo;
   delete JEC_L2L3_calo;
   delete JEC_L2L3res_calo;

   vParam_res_calo.clear();
   delete JEC_res_calo;
   delete jecUnc_calo;

   vParam_res_pf.clear();
   delete ResJetCorPar_pf;
   delete JEC_res_pf;
   delete jecUnc_pf;
   */

   delete jecUnc_Calo;
   delete jecUnc_PF;

}

// ------------ method called once each job just before starting event loop  ------------
void 
BEANmaker::beginJob()
{

  numEvents = 0;
  numGoodVertexEvents = 0;
  numFilterOutScrapingEvents = 0;
  numHEEPele = 0;
  numHEEPeleEB = 0;
  numHEEPeleEE = 0;

}


// ------------ method called once each run just before starting event loop  ------------
void 
BEANmaker::beginRun( edm::Run& run, const edm::EventSetup& c )
{

  bool changed(true);
  if (hltConfig_.init(run,c,hltProcessName_,changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
      // The HLT config has actually changed wrt the previous Run, hence rebook your
      // histograms or do anything else dependent on the revised HLT config
      std::cout << "Initalizing HLTConfigProvider"  << std::endl;
    }
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    std::cout << " HLT config extraction failure with process name " << hltProcessName_ << std::endl;
    // In this case, all access methods will return empty values!
  }

}


// ------------ method called once each run just before starting event loop  ------------
void 
BEANmaker::endRun( edm::Run& run, edm::EventSetup& c )
{

}


// ------------ method called once each job just after ending the event loop  ------------
void 
BEANmaker::endJob() {

  if( verbose_ ){
    std::cout << " ****************************************** " << std::endl;
    std::cout << "    Number of events = " << numEvents << std::endl;
    std::cout << "    Number of GoodVertex events = " << numGoodVertexEvents << " (" << int( (double)numGoodVertexEvents/(double)numEvents * 100 ) << "%)" << std::endl;
    std::cout << "    Number of FilterOutScraping events = " << numFilterOutScrapingEvents << " (" << int( (double)numFilterOutScrapingEvents/(double)numEvents * 100 ) << "%)" << std::endl;
    std::cout << "    Number of HEEP electrons = " << numHEEPele << std::endl;
    std::cout << "    Number of barrel HEEP electrons = " << numHEEPeleEB << std::endl;
    std::cout << "    Number of endcap HEEP electrons = " << numHEEPeleEE << std::endl;
    std::cout << " ****************************************** " << std::endl;
  }
}

bool BEANmaker::isActive(int word, int bit) {
  if( word & (1<<bit) ) return true;
  return false;
}

int BEANmaker::count_hits( const std::vector<CaloTowerPtr> & towers )
{
  int nHit = 0;
  for ( unsigned int iTower = 0; iTower < towers.size() ; ++iTower ) {
    const std::vector<DetId>& cellIDs = towers[iTower]->constituents();  // cell == recHit
    nHit += cellIDs.size();
  }
  return nHit;
}

//define this as a plug-in
DEFINE_FWK_MODULE(BEANmaker);
