//// -*- C++ -*-
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
// $Id: BEANmaker.cc,v 1.1 2012/10/29 21:01:52 nvallsve Exp $
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
#include "ProductArea/BNcollections/interface/BNgenjet.h"
#include "ProductArea/BNcollections/interface/BNmcparticle.h"
#include "ProductArea/BNcollections/interface/BNmet.h"
#include "ProductArea/BNcollections/interface/BNmuon.h"
#include "ProductArea/BNcollections/interface/BNphoton.h"
#include "ProductArea/BNcollections/interface/BNsupercluster.h"
#include "ProductArea/BNcollections/interface/BNtau.h"
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

#include "CondFormats/L1TObjects/interface/L1GtPrescaleFactors.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsAlgoTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsTechTrigRcd.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
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
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "SusyAnalysis/EventSelector/interface/uncorrectionTypeMET.h"

#include "DataFormats/Common/interface/View.h"

#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

#include "CMGTools/External/interface/PileupJetIdentifier.h"

#include "NtupleMaker/BEANmaker/interface/AnglesUtil.h"

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
  bool tauIsInTheCracks(float);
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
  edm::InputTag pfeleTag_;
  edm::InputTag pfeleLooseTag_;
  edm::InputTag calojetTag_;
  edm::InputTag pfjetTag_;
  edm::InputTag genjetTag_;
  //edm::InputTag calometTag_;
  edm::InputTag pfmetTag_;
  edm::InputTag pfmetTag_type1correctedRECO_;
  edm::InputTag pfmetTag_uncorrectedPF2PAT_;
  edm::InputTag pfmetTag_uncorrectedRECO_;
  edm::InputTag tcmetTag_;
  edm::InputTag muonTag_;
  edm::InputTag pfmuonTag_;
  edm::InputTag pfmuonLooseTag_;
  edm::InputTag cocktailmuonTag_;
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
  edm::InputTag tauTag_;
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
  bool fillTrackHitInfo_;

  std::vector<std::string> hlt_sd_eg_name;
  std::vector<std::string> hlt_sd_jetmettau_name;

  int numEvents, numGoodVertexEvents, numFilterOutScrapingEvents;

  HLTConfigProvider hltConfig_;

  const L1GtPrescaleFactors* m_l1GtPfAlgo;
  unsigned long long m_l1GtPfAlgoCacheID;

  const std::vector<std::vector<int> >* m_prescaleFactorsAlgoTrig;

  const TrackerGeometry* m_tracker;
};

//
// constants, enums and typedefs
//
typedef std::vector<BNelectron>     BNelectronCollection;
typedef std::vector<BNjet>          BNjetCollection;
typedef std::vector<BNgenjet>       BNgenjetCollection;
typedef std::vector<BNevent>        BNeventCollection;
typedef std::vector<BNmcparticle>   BNmcparticleCollection;
typedef std::vector<BNmet>          BNmetCollection;
typedef std::vector<BNmuon>         BNmuonCollection;
typedef std::vector<BNphoton>       BNphotonCollection;
typedef std::vector<BNsupercluster> BNsuperclusterCollection;
typedef std::vector<BNtrack>        BNtrackCollection;
typedef std::vector<BNtrigger>      BNtriggerCollection;
typedef std::vector<BNbxlumi>       BNbxlumiCollection;
typedef std::vector<BNtrigobj>      BNtrigobjCollection;
typedef std::vector<BNprimaryvertex> BNprimaryvertexCollection;
typedef std::vector<BNtau>          BNtauCollection;

//
// static data member definitions
//
static const char* kSC        = "corHybridSCandMulti5x5WithPreshower";
static const char* kHLT       = "HLT";
static const char* kL1Talgo   = "L1Talgo";
static const char* kL1Ttech   = "L1Ttech";
static const char* kBXlumi    = "BXlumi";
static const char* kMCpar     = "MCstatus3";
static const char* kMCele     = "MCeleStatus1";
static const char* kMCmu      = "MCmuStatus1";
static const char* kHLTobj    = "HLT";
static const char* kL1EmParticlesIso      = "L1EmParticlesIsolated";
static const char* kL1EmParticlesNonIso   = "L1EmParticlesNonIsolated";
static const char* kL1EtMissParticlesMET  = "L1EtMissParticlesMET";
static const char* kL1EtMissParticlesMHT  = "L1EtMissParticlesMHT";
static const char* kL1JetParticlesCentral = "L1JetParticlesCentral";
static const char* kL1JetParticlesForward = "L1JetParticlesForward";
static const char* kL1JetParticlesTau     = "L1JetParticlesTau";
static const char* kL1MuonParticles       = "L1MuonParticles";


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
  verbose_(iConfig.getParameter<bool>("verbose")),
  fillTrackHitInfo_(iConfig.getParameter<bool>("fillTrackHitInfo"))
{

  // Define InputTags 
  eleTag_ = iConfig.getParameter<edm::InputTag>("eleTag");
  pfeleTag_ = iConfig.getParameter<edm::InputTag>("pfeleTag");
  pfeleLooseTag_ = iConfig.getParameter<edm::InputTag>("pfeleLooseTag");
  calojetTag_ = iConfig.getParameter<edm::InputTag>("calojetTag");
  pfjetTag_ = iConfig.getParameter<edm::InputTag>("pfjetTag");
  genjetTag_ = iConfig.getParameter<edm::InputTag>("genjetTag");
  //calometTag_ = iConfig.getParameter<edm::InputTag>("calometTag");
  pfmetTag_ = iConfig.getParameter<edm::InputTag>("pfmetTag");
  pfmetTag_type1correctedRECO_ = iConfig.getParameter<edm::InputTag>("pfmetTag_type1correctedRECO");
  pfmetTag_uncorrectedPF2PAT_  = iConfig.getParameter<edm::InputTag>("pfmetTag_uncorrectedPF2PAT");
  pfmetTag_uncorrectedRECO_    = iConfig.getParameter<edm::InputTag>("pfmetTag_uncorrectedRECO");
  tcmetTag_ = iConfig.getParameter<edm::InputTag>("tcmetTag");
  muonTag_ = iConfig.getParameter<edm::InputTag>("muonTag");
  pfmuonTag_ = iConfig.getParameter<edm::InputTag>("pfmuonTag");
  pfmuonLooseTag_ = iConfig.getParameter<edm::InputTag>("pfmuonLooseTag");
  cocktailmuonTag_ = iConfig.getParameter<edm::InputTag>("cocktailmuonTag");
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
  tauTag_ = iConfig.getParameter<edm::InputTag>("tauTag");

  reducedBarrelRecHitCollection_ = iConfig.getParameter<edm::InputTag>("reducedBarrelRecHitCollection");
  reducedEndcapRecHitCollection_ = iConfig.getParameter<edm::InputTag>("reducedEndcapRecHitCollection");

  hltProcessName_ = iConfig.getParameter<std::string>("hltProcessName");


  // Register products
  produces<BNeventCollection>().setBranchAlias("event");
  produces<BNelectronCollection>(eleTag_.label()).setBranchAlias("electrons");
  produces<BNjetCollection>(calojetTag_.label()).setBranchAlias("calojets");
  produces<BNjetCollection>(pfjetTag_.label()).setBranchAlias("pfjets");
  produces<BNgenjetCollection>(genjetTag_.label()).setBranchAlias("genjets");
  //produces<BNmetCollection>(calometTag_.label()).setBranchAlias("calomet");
  produces<BNmetCollection>(pfmetTag_.label()).setBranchAlias("pfmet");
  produces<BNmetCollection>(std::string(pfmetTag_type1correctedRECO_.label() + "BN")).setBranchAlias("pfmet_type1correctedRECO");
  produces<BNmetCollection>(std::string(pfmetTag_uncorrectedPF2PAT_.label() + "BN")).setBranchAlias("pfmet_uncorrectedPF2PAT");
  produces<BNmetCollection>(std::string(pfmetTag_uncorrectedRECO_.label() + "BN")).setBranchAlias("pfmet_uncorrectedRECO");
  produces<BNmetCollection>(tcmetTag_.label()).setBranchAlias("tcmet");
  produces<BNmuonCollection>(muonTag_.label()).setBranchAlias("muons");
  produces<BNelectronCollection>(pfeleTag_.label()).setBranchAlias("pfelectrons");
  produces<BNelectronCollection>(pfeleLooseTag_.label()).setBranchAlias("pfelectrons_loose");
  produces<BNmuonCollection>(pfmuonTag_.label()).setBranchAlias("pfmuons");
  produces<BNmuonCollection>(pfmuonLooseTag_.label()).setBranchAlias("pfmuons_loose");
  produces<BNmuonCollection>(cocktailmuonTag_.label()).setBranchAlias("cocktailmuons");
  produces<BNphotonCollection>(photonTag_.label()).setBranchAlias("photons");
  produces<BNsuperclusterCollection>(kSC).setBranchAlias("superclusters");
  produces<BNtrackCollection>(trackTag_.label()).setBranchAlias("tracks");
  produces<BNtriggerCollection>(kHLT).setBranchAlias("trigger");
  produces<BNtriggerCollection>(kL1Talgo).setBranchAlias("L1Talgo");
  produces<BNtriggerCollection>(kL1Ttech).setBranchAlias("L1Ttech");
  produces<BNbxlumiCollection>(kBXlumi).setBranchAlias("bxlumi");
  produces<BNmcparticleCollection>(kMCpar).setBranchAlias("mcparticles");
  produces<BNmcparticleCollection>(kMCele).setBranchAlias("mcelectrons");
  produces<BNmcparticleCollection>(kMCmu).setBranchAlias("mcmuons");
  produces<BNtrigobjCollection>(kHLTobj).setBranchAlias("hltobjs");
  produces<BNtrigobjCollection>(kL1EmParticlesIso).setBranchAlias("l1isoemobjs");
  produces<BNtrigobjCollection>(kL1EmParticlesNonIso).setBranchAlias("l1nonisoemobjs");
  produces<BNtrigobjCollection>(kL1EtMissParticlesMET).setBranchAlias("l1metobjs");
  produces<BNtrigobjCollection>(kL1EtMissParticlesMHT).setBranchAlias("l1mhtobjs");
  produces<BNtrigobjCollection>(kL1JetParticlesCentral).setBranchAlias("l1cenjetobjs");
  produces<BNtrigobjCollection>(kL1JetParticlesForward).setBranchAlias("l1forjetobjs");
  produces<BNtrigobjCollection>(kL1JetParticlesTau).setBranchAlias("l1taujetobjs");
  produces<BNtrigobjCollection>(kL1MuonParticles).setBranchAlias("l1muonobjs");
  produces<BNprimaryvertexCollection>(pvTag_.label()).setBranchAlias("pvs");
  produces<BNtauCollection>(tauTag_.label()).setBranchAlias("taus");

  m_l1GtPfAlgoCacheID = 0ULL;

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

   edm::Handle<edm::View<pat::Electron> > pfelectronHandle;
   iEvent.getByLabel(pfeleTag_,pfelectronHandle);

   edm::Handle<edm::View<pat::Electron> > pfelectronLooseHandle;
   iEvent.getByLabel(pfeleLooseTag_,pfelectronLooseHandle);

   edm::Handle<edm::View<pat::Jet> > calojetHandle;
   iEvent.getByLabel(calojetTag_,calojetHandle);

   edm::Handle<edm::View<pat::Jet> > pfjetHandle;
   iEvent.getByLabel(pfjetTag_,pfjetHandle);

   edm::Handle<reco::GenJetCollection > genjetHandle;
   iEvent.getByLabel(genjetTag_,genjetHandle);

   //edm::Handle<edm::View<pat::MET> > calometHandle;
   //iEvent.getByLabel(calometTag_,calometHandle);

   edm::Handle<edm::View<pat::MET> > pfmetHandle;
   iEvent.getByLabel(pfmetTag_,pfmetHandle);

   edm::Handle<vector<reco::PFMET> > pfmetHandle_type1correctedRECO;
   iEvent.getByLabel(pfmetTag_type1correctedRECO_,pfmetHandle_type1correctedRECO);

   edm::Handle<edm::View<pat::MET> > pfmetHandle_uncorrectedPF2PAT;
   iEvent.getByLabel(pfmetTag_uncorrectedPF2PAT_,pfmetHandle_uncorrectedPF2PAT);

   edm::Handle<vector<reco::PFMET> > pfmetHandle_uncorrectedRECO;
   iEvent.getByLabel(pfmetTag_uncorrectedRECO_,pfmetHandle_uncorrectedRECO);

   edm::Handle<edm::View<pat::MET> > tcmetHandle;
   iEvent.getByLabel(tcmetTag_,tcmetHandle);

   edm::Handle<edm::View<pat::Muon> > muonHandle;
   iEvent.getByLabel(muonTag_,muonHandle);

   edm::Handle<edm::View<pat::Muon> > pfmuonHandle;
   iEvent.getByLabel(pfmuonTag_,pfmuonHandle);

   edm::Handle<edm::View<pat::Muon> > pfmuonLooseHandle;
   iEvent.getByLabel(pfmuonLooseTag_,pfmuonLooseHandle);

   edm::Handle<edm::View<reco::Muon> > cocktailmuonHandle;
   iEvent.getByLabel(cocktailmuonTag_,cocktailmuonHandle);

   edm::Handle<edm::View<pat::Photon> > photonHandle;
   iEvent.getByLabel(photonTag_,photonHandle);

   edm::Handle<reco::TrackCollection > trackHandle;
   iEvent.getByLabel(trackTag_,trackHandle);

   edm::Handle<edm::View<pat::Tau> > tauHandle;
   iEvent.getByLabel(tauTag_,tauHandle);

   edm::Handle<reco::SuperClusterCollection > EBsuperclusterHandle;
   iEvent.getByLabel(EBsuperclusterTag_,EBsuperclusterHandle);

   edm::Handle<reco::SuperClusterCollection > EEsuperclusterHandle;
   iEvent.getByLabel(EEsuperclusterTag_,EEsuperclusterHandle);

   edm::Handle<reco::GenParticleCollection > genParticles;
   iEvent.getByLabel(genParticleTag_,genParticles);

   edm::Handle< double > rhoHandle;
   iEvent.getByLabel(edm::InputTag("kt6PFJets","rho"), rhoHandle);
   double rho_event = *rhoHandle;   

   edm::Handle< double > rhoHandle_CentralChargedPileUp;
   iEvent.getByLabel(edm::InputTag("kt6PFJetsCentralChargedPileUp","rho"), rhoHandle_CentralChargedPileUp);
   double rho_event_CentralChargedPileUp = ( (rhoHandle_CentralChargedPileUp.isValid()) ) ? *rhoHandle_CentralChargedPileUp : -99;   

   edm::Handle< double > rhoHandle_CentralNeutral;
   iEvent.getByLabel(edm::InputTag("kt6PFJetsCentralNeutral","rho"), rhoHandle_CentralNeutral);
   double rho_event_CentralNeutral = ( (rhoHandle_CentralNeutral.isValid()) ) ? *rhoHandle_CentralNeutral : -99;

   edm::Handle< double > rhoHandle_CentralNeutralTight;
   iEvent.getByLabel(edm::InputTag("kt6PFJetsCentralNeutralTight","rho"), rhoHandle_CentralNeutralTight);
   double rho_event_CentralNeutralTight = ( (rhoHandle_CentralNeutralTight.isValid()) ) ? *rhoHandle_CentralNeutralTight : -99;

   edm::Handle<reco::ConversionCollection> conversionsHandle;
   iEvent.getByLabel("allConversions", conversionsHandle);


   bool produceElectron = ( (eleTag_.label() == "none") ) ? false : true;
   bool producePFelectron = ( (pfeleTag_.label() == "none") ) ? false : true;
   bool producePFelectronLoose = ( (pfeleLooseTag_.label() == "none") ) ? false : true;
   bool produceCaloJet = ( (calojetTag_.label() == "none") ) ? false : true;
   bool producePFJet = ( (pfjetTag_.label() == "none") ) ? false : true;
   bool produceGenJet = ( (genjetTag_.label() == "none") ) ? false : true;
   //bool produceCaloMET = ( (calometTag_.label() == "none") ) ? false : true;
   bool producePFMET = ( (pfmetTag_.label() == "none") ) ? false : true;
   bool producePFMET_type1correctedRECO = ( (pfmetTag_type1correctedRECO_.label() == "none") ) ? false : true;
   bool producePFMET_uncorrectedPF2PAT  = ( (pfmetTag_uncorrectedPF2PAT_.label() == "none") ) ? false : true;
   bool producePFMET_uncorrectedRECO    = ( (pfmetTag_uncorrectedRECO_.label() == "none") ) ? false : true;
   bool produceTCMET = ( (tcmetTag_.label() == "none") ) ? false : true;
   bool produceMuon = ( (muonTag_.label() == "none") ) ? false : true;
   bool producePFMuon = ( (pfmuonTag_.label() == "none") ) ? false : true;
   bool producePFMuonLoose = ( (pfmuonLooseTag_.label() == "none") ) ? false : true;
   bool produceCocktailMuon = ( (cocktailmuonTag_.label() == "none") ) ? false : true;
   bool produceTau = ( (tauTag_.label() == "none") ) ? false : true;
   bool producePhoton = ( (photonTag_.label() == "none") ) ? false : true;
   bool produceSCsEB = ( (EBsuperclusterTag_.label() == "none") ) ? false : true;
   bool produceSCsEE = ( (EEsuperclusterTag_.label() == "none") ) ? false : true;
   bool produceTrack = ( (trackTag_.label() == "none") ) ? false : true;
   bool produceGenParticle = ( (genParticleTag_.label() == "none") ) ? false : true;



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

   float sum_nvtx_gen = 0, sum_nvtx_true = 0;
   int npv_gen = -1;
   float npv_true = -1;
   int nm1 = -1, n0 = -1, np1 = -1;
   float nm1_true = -1, n0_true = -1, np1_true = -1;

   if( (PupInfo.isValid()) ){
     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

       int BX = PVI->getBunchCrossing();

       sum_nvtx_gen  += float(PVI->getPU_NumInteractions());
       sum_nvtx_true += float(PVI->getTrueNumInteractions());

       if( BX==0 ){
	 npv_gen = PVI->getPU_NumInteractions();
	 npv_true = PVI->getTrueNumInteractions();
       }

       if(BX == -1) { 
	 nm1 = PVI->getPU_NumInteractions();
	 nm1_true = PVI->getTrueNumInteractions();
       }
       else if(BX == 0) { 
	 n0 = PVI->getPU_NumInteractions();
	 n0_true = PVI->getTrueNumInteractions();
       }
       else if(BX == 1) { 
	 np1 = PVI->getPU_NumInteractions();
	 np1_true = PVI->getTrueNumInteractions();
       }
     }
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


   edm::Handle<std::vector<double> > q2weights;
   bool hasQ2weights = iEvent.getByLabel( "q2weights", q2weights );

   double Q2ScaleUpWgt=1, Q2ScaleDownWgt=1;
   if( (hasQ2weights) ){
     int q2size = q2weights->size();

     if( q2size!=2 ){
       std::cout << " ERROR!! q2weights->size() = " << q2size << std::endl;
     }

     if( q2size>0 ){
       Q2ScaleUpWgt = q2weights->at(0);
       if( q2size>1 ){
	 Q2ScaleDownWgt = q2weights->at(1);
       }
     }

     //std::cout << "  ==> Q2ScaleUpWgt = " << Q2ScaleUpWgt << ",\t Q2ScaleDownWgt = " << Q2ScaleDownWgt << std::endl;
   }


   edm::Handle<GenEventInfoProduct> genEvtInfo;
   bool hasGenEvtInfo = iEvent.getByLabel( "generator", genEvtInfo );

   double qScale=-1, alphaQCD=-1, alphaQED=-1, pthat=-1, scalePDF=-1, x1=-1, x2=-1, xPDF1=-1, xPDF2=-1;
   int id1=-99, id2=-99;
   if( (hasGenEvtInfo) ){
     qScale = genEvtInfo->qScale();
     alphaQCD = genEvtInfo->alphaQCD();
     alphaQED = genEvtInfo->alphaQED();
     pthat = ( genEvtInfo->hasBinningValues() ? (genEvtInfo->binningValues())[0] : 0.0);
     scalePDF = genEvtInfo->pdf()->scalePDF;
     id1 = genEvtInfo->pdf()->id.first;
     id2 = genEvtInfo->pdf()->id.second;
     x1 = genEvtInfo->pdf()->x.first;
     x2 = genEvtInfo->pdf()->x.second;
     xPDF1 = genEvtInfo->pdf()->xPDF.first;
     xPDF2 = genEvtInfo->pdf()->xPDF.second;
   }


   Handle< bool > hcalNoiseFilterHandle;
   iEvent.getByLabel("HBHENoiseFilterResultProducer","HBHENoiseFilterResult", hcalNoiseFilterHandle);

   Bool_t hcalNoiseFilter = false;
   if( hcalNoiseFilterHandle.isValid() ) hcalNoiseFilter = (Bool_t)(*hcalNoiseFilterHandle);
   //else if( sample_<0 ) std::cout << " NoiseFilter =====> hcalNoiseFilterHandle.isValid()==false " << std::endl;


   Handle<HcalNoiseSummary> HcalNoiseSummaryHandle;
   iEvent.getByLabel("hcalnoise", HcalNoiseSummaryHandle);

   bool passLooseNoiseFilter=false, passTightNoiseFilter=false;
   if( HcalNoiseSummaryHandle.isValid() ){
     passLooseNoiseFilter = HcalNoiseSummaryHandle->passLooseNoiseFilter();
     passTightNoiseFilter = HcalNoiseSummaryHandle->passTightNoiseFilter();
   }
   //else if( sample_<0 ) std::cout << " NoiseFilter =====> HcalNoiseSummaryHandle.isValid()==false " << std::endl;



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
   //else std::cout << " NoiseFilter =====> TheBeamHaloSummary.isValid()==false " << std::endl;



   edm::Handle<reco::VertexCollection> vtxHandle;
   iEvent.getByLabel(pvTag_,vtxHandle);
   reco::VertexCollection vtxs = *vtxHandle;


   bool GoodVertex = false;

   int numPVs = 0;
   double PVx=-99,PVy=-99,PVz =-99;
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
	 PVx = vtx->x();
	 PVy = vtx->y();
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



   /////////////////////////////////////////////
   ///////
   ///////   Fill the electron collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNelectronCollection> bnelectrons(new BNelectronCollection);
   if( produceElectron ){
     edm::View<pat::Electron> electrons = *electronHandle;

     for( edm::View<pat::Electron>::const_iterator ele = electrons.begin(); ele!=electrons.end(); ++ele ){

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

       MyElectron.pIn = elePin;
       MyElectron.pOut = elePout;
       MyElectron.EscOverPin = ele->eSuperClusterOverP();
       MyElectron.EseedOverPout = ele->eSeedClusterOverPout();
       MyElectron.hadOverEm = ele->hadronicOverEm();
       MyElectron.fbrem = fabs(elePin-elePout)/elePin;
       MyElectron.absInvEMinusInvPin = fabs( 1/ele->ecalEnergy()-(ele->eSuperClusterOverP()/ele->ecalEnergy()) );
       MyElectron.delPhiIn = ele->deltaPhiSuperClusterTrackAtVtx();
       MyElectron.delEtaIn = ele->deltaEtaSuperClusterTrackAtVtx();

       /*
       MyElectron.eidRobustHighEnergy = ele->electronID("eidRobustHighEnergy");
       MyElectron.eidRobustLoose = ele->electronID("eidRobustLoose");
       MyElectron.eidRobustTight = ele->electronID("eidRobustTight");
       MyElectron.eidLoose = ele->electronID("eidLoose");
       MyElectron.eidTight = ele->electronID("eidTight");
       MyElectron.eidVeryLooseMC = ele->electronID("eidVeryLooseMC");
       MyElectron.eidLooseMC = ele->electronID("eidLooseMC");
       MyElectron.eidMediumMC = ele->electronID("eidMediumMC");
       MyElectron.eidTightMC = ele->electronID("eidTightMC");
       MyElectron.eidSuperTightMC = ele->electronID("eidSuperTightMC");
       MyElectron.eidHyperTight1MC = ele->electronID("eidHyperTight1MC");
       MyElectron.eidHyperTight2MC = ele->electronID("eidHyperTight2MC");
       MyElectron.eidHyperTight3MC = ele->electronID("eidHyperTight3MC");
       MyElectron.eidHyperTight4MC = ele->electronID("eidHyperTight4MC");
       */
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
       MyElectron.mvaTrigV0 = ele->electronID("mvaTrigV0");
       MyElectron.mvaNonTrigV0 = ele->electronID("mvaNonTrigV0");
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

       MyElectron.IP = ele->dB(pat::Electron::PV3D);
       MyElectron.IPError = ele->edB(pat::Electron::PV3D);

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
         MyElectron.correctedD0Vertex = ele->gsfTrack()->dxy(vertexPosition);
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
	 int genId = ele->genLepton()->pdgId();

	 MyElectron.genId = ele->genLepton()->pdgId();
	 MyElectron.genET = ele->genLepton()->et();
	 MyElectron.genPT = ele->genLepton()->pt();
	 MyElectron.genPhi = ele->genLepton()->phi();
	 MyElectron.genEta = ele->genLepton()->eta();
	 MyElectron.genCharge = ele->genLepton()->charge();

	 MyElectron.genNumberOfMothers = ele->genLepton()->numberOfMothers();

	 if( (ele->genLepton()->numberOfMothers()>0) ){
	   const reco::Candidate *ElectronMother = ele->genLepton()->mother();
	   bool staytrapped = true;
	   while( (ElectronMother->pdgId()==genId && staytrapped) ){
	     if( ElectronMother->numberOfMothers()>=1 ) ElectronMother = ElectronMother->mother();
	     else staytrapped = false;
	   }
       
	   MyElectron.genMotherId = ElectronMother->pdgId();
	   MyElectron.genMotherET = ElectronMother->et();
	   MyElectron.genMotherPT = ElectronMother->pt();
	   MyElectron.genMotherPhi = ElectronMother->phi();
	   MyElectron.genMotherEta = ElectronMother->eta();
	   MyElectron.genMotherCharge = ElectronMother->charge();
	 }

	 if( (ele->genLepton()->numberOfMothers()>=1) ){
	   const reco::Candidate *ElectronMother0 = ele->genLepton()->mother(0);
	   const reco::Candidate *ElectronMother1 = 0;

	   if( (ele->genLepton()->numberOfMothers()>=2) ) ElectronMother1 = ele->genLepton()->mother(1);

	   int mother0id = ElectronMother0->pdgId();
	   int mother1id = ( ElectronMother1!=0 ) ? ElectronMother1->pdgId() : -99;

	   bool staytrapped = true;
	   while( (mother0id==genId || mother1id==genId) && staytrapped ){
	     if( mother0id==genId && (ElectronMother0!=0) ){
	       if( ElectronMother0->numberOfMothers()>=1 ){
		 ElectronMother0 = ElectronMother0->mother(0);
		 mother0id = ElectronMother0->pdgId();
		 mother1id = -99;
		 if( ElectronMother0->numberOfMothers()>=2 ){
		   ElectronMother1 = ElectronMother0->mother(1);
		   mother1id = ElectronMother1->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else if( mother1id==genId && (ElectronMother1!=0) ){
	       if( ElectronMother1->numberOfMothers()>=1 ){
		 ElectronMother1 = ElectronMother1->mother(0);
		 mother1id = ElectronMother1->pdgId();
		 mother0id = -99;
		 if( ElectronMother1->numberOfMothers()>=2 ){
		   ElectronMother0 = ElectronMother1->mother(1);
		   mother0id = ElectronMother0->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else staytrapped = false;
	   }

	   if( mother0id!=-99 ){
	     MyElectron.genMother0Id = ElectronMother0->pdgId();

	     if( (ElectronMother0->numberOfMothers()>=1) ){
	       const reco::Candidate *ElectronGrandMother0 = ElectronMother0->mother(0);
	       const reco::Candidate *ElectronGrandMother1 = 0;

	       if( (ElectronMother0->numberOfMothers()>=2) ) ElectronGrandMother1 = ElectronMother0->mother(1);

	       int gmother0id = ElectronGrandMother0->pdgId();
	       int gmother1id = ( ElectronGrandMother1!=0 ) ? ElectronGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother0id || gmother1id==mother0id) && staytrapped ){
		 if( gmother0id==mother0id && (ElectronGrandMother0!=0) ){
		   if( ElectronGrandMother0->numberOfMothers()>=1 ){
		     ElectronGrandMother0 = ElectronGrandMother0->mother(0);
		     gmother0id = ElectronGrandMother0->pdgId();
		     gmother1id = -99;
		     if( ElectronGrandMother0->numberOfMothers()>=2 ){
		       ElectronGrandMother1 = ElectronGrandMother0->mother(1);
		       gmother1id = ElectronGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother0id && (ElectronGrandMother1!=0) ){
		   if( ElectronGrandMother1->numberOfMothers()>=1 ){
		     ElectronGrandMother1 = ElectronGrandMother1->mother(0);
		     gmother1id = ElectronGrandMother1->pdgId();
		     gmother0id = -99;
		     if( ElectronGrandMother1->numberOfMothers()>=2 ){
		       ElectronGrandMother0 = ElectronGrandMother1->mother(1);
		       gmother0id = ElectronGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyElectron.genGrandMother00Id = gmother0id;
	       MyElectron.genGrandMother01Id = gmother1id;
	     }
	   }
	   if( mother1id!=-99 ){
	     MyElectron.genMother1Id = ElectronMother1->pdgId();

	     if( (ElectronMother1->numberOfMothers()>=1) ){
	       const reco::Candidate *ElectronGrandMother0 = ElectronMother1->mother(0);
	       const reco::Candidate *ElectronGrandMother1 = 0;

	       if( (ElectronMother0->numberOfMothers()>=2) ) ElectronGrandMother1 = ElectronMother1->mother(1);

	       int gmother0id = ElectronGrandMother0->pdgId();
	       int gmother1id = ( ElectronGrandMother1!=0 ) ? ElectronGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother1id || gmother1id==mother1id) && staytrapped ){
		 if( gmother0id==mother1id && (ElectronGrandMother0!=0) ){
		   if( ElectronGrandMother0->numberOfMothers()>=1 ){
		     ElectronGrandMother0 = ElectronGrandMother0->mother(0);
		     gmother0id = ElectronGrandMother0->pdgId();
		     gmother1id = -99;
		     if( ElectronGrandMother0->numberOfMothers()>=2 ){
		       ElectronGrandMother1 = ElectronGrandMother0->mother(1);
		       gmother1id = ElectronGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother1id && (ElectronGrandMother1!=0) ){
		   if( ElectronGrandMother1->numberOfMothers()>=1 ){
		     ElectronGrandMother1 = ElectronGrandMother1->mother(0);
		     gmother1id = ElectronGrandMother1->pdgId();
		     gmother0id = -99;
		     if( ElectronGrandMother1->numberOfMothers()>=2 ){
		       ElectronGrandMother0 = ElectronGrandMother1->mother(1);
		       gmother0id = ElectronGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyElectron.genGrandMother10Id = gmother0id;
	       MyElectron.genGrandMother11Id = gmother1id;
	     }
	   }
	 }

       }// end check if genLepton

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
       MyElectron.passConvVeto = ele->passConversionVeto();
      

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

     }
   }



   /////////////////////////////////////////////
   ///////
   ///////   Fill the pfelectron collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNelectronCollection> bnpfelectrons(new BNelectronCollection);
   if( producePFelectron ){
     edm::View<pat::Electron> pfelectrons = *pfelectronHandle;

     for( edm::View<pat::Electron>::const_iterator pfele = pfelectrons.begin(); pfele!=pfelectrons.end(); ++pfele ){

       double pfelePin = pfele->trackMomentumAtVtx().R();
       double pfelePout = pfele->trackMomentumOut().R();

       BNelectron MyPfelectron;

       // general kinematic variables
       MyPfelectron.energy = pfele->energy();
       MyPfelectron.gsfEt = pfele->et();
       MyPfelectron.pt = pfele->pt();
       MyPfelectron.px = pfele->px();
       MyPfelectron.py = pfele->py();
       MyPfelectron.pz = pfele->pz();
       MyPfelectron.phi = pfele->phi();
       MyPfelectron.eta = pfele->eta();
       MyPfelectron.theta = pfele->theta();

       MyPfelectron.charge = pfele->charge();
       MyPfelectron.classification = pfele->classification();
       MyPfelectron.vx = pfele->vx();
       MyPfelectron.vy = pfele->vy();
       MyPfelectron.vz = pfele->vz();

       MyPfelectron.pIn = pfelePin;
       MyPfelectron.pOut = pfelePout;
       MyPfelectron.EscOverPin = pfele->eSuperClusterOverP();
       MyPfelectron.EseedOverPout = pfele->eSeedClusterOverPout();
       MyPfelectron.hadOverEm = pfele->hadronicOverEm();
       MyPfelectron.fbrem = fabs(pfelePin-pfelePout)/pfelePin;
       MyPfelectron.absInvEMinusInvPin = fabs( 1/pfele->ecalEnergy()-(pfele->eSuperClusterOverP()/pfele->ecalEnergy()) );
       MyPfelectron.delPhiIn = pfele->deltaPhiSuperClusterTrackAtVtx();
       MyPfelectron.delEtaIn = pfele->deltaEtaSuperClusterTrackAtVtx();
       /*
       MyPfelectron.eidRobustHighEnergy = pfele->electronID("eidRobustHighEnergy");
       MyPfelectron.eidRobustLoose = pfele->electronID("eidRobustLoose");
       MyPfelectron.eidRobustTight = pfele->electronID("eidRobustTight");
       MyPfelectron.eidLoose = pfele->electronID("eidLoose");
       MyPfelectron.eidTight = pfele->electronID("eidTight");
       MyPfelectron.eidVeryLooseMC = pfele->electronID("eidVeryLooseMC");
       MyPfelectron.eidLooseMC = pfele->electronID("eidLooseMC");
       MyPfelectron.eidMediumMC = pfele->electronID("eidMediumMC");
       MyPfelectron.eidTightMC = pfele->electronID("eidTightMC");
       MyPfelectron.eidSuperTightMC = pfele->electronID("eidSuperTightMC");
       MyPfelectron.eidHyperTight1MC = pfele->electronID("eidHyperTight1MC");
       MyPfelectron.eidHyperTight2MC = pfele->electronID("eidHyperTight2MC");
       MyPfelectron.eidHyperTight3MC = pfele->electronID("eidHyperTight3MC");
       MyPfelectron.eidHyperTight4MC = pfele->electronID("eidHyperTight4MC");
       */
       MyPfelectron.particleIso = pfele->particleIso();     
       MyPfelectron.chargedHadronIso = pfele->chargedHadronIso();     
       MyPfelectron.neutralHadronIso = pfele->neutralHadronIso();     
       MyPfelectron.photonIso = pfele->photonIso();     
       MyPfelectron.puChargedHadronIso = pfele->puChargedHadronIso();     

       MyPfelectron.chargedHadronIsoDR03 = pfele->userIso(0);     
       MyPfelectron.neutralHadronIsoDR03 = pfele->userIso(1);
       MyPfelectron.photonIsoDR03 = pfele->userIso(2);
       MyPfelectron.puChargedHadronIsoDR03 = pfele->userIso(3);

       MyPfelectron.chargedHadronIsoDR04 = pfele->userIso(4);     
       MyPfelectron.neutralHadronIsoDR04 = pfele->userIso(5);
       MyPfelectron.photonIsoDR04 = pfele->userIso(6);
       MyPfelectron.puChargedHadronIsoDR04 = pfele->userIso(7);

       MyPfelectron.rhoPrime = std::max(rho_event, 0.0);

       MyPfelectron.trackIso = pfele->trackIso();
       MyPfelectron.ecalIso = pfele->ecalIso();
       MyPfelectron.hcalIso = pfele->hcalIso();
       MyPfelectron.caloIso = pfele->caloIso();

       MyPfelectron.trackIsoDR03 = pfele->dr03TkSumPt();
       MyPfelectron.ecalIsoDR03 = pfele->dr03EcalRecHitSumEt();
       MyPfelectron.hcalIsoDR03 = pfele->dr03HcalTowerSumEt();
       MyPfelectron.hcalIsoDR03depth1 = pfele->dr03HcalDepth1TowerSumEt();
       MyPfelectron.hcalIsoDR03depth2 = pfele->dr03HcalDepth2TowerSumEt();
       MyPfelectron.caloIsoDR03 = pfele->dr03EcalRecHitSumEt()+pfele->dr03HcalTowerSumEt();

       MyPfelectron.trackIsoDR04 = pfele->dr04TkSumPt();
       MyPfelectron.ecalIsoDR04 = pfele->dr04EcalRecHitSumEt();
       MyPfelectron.hcalIsoDR04 = pfele->dr04HcalTowerSumEt();
       MyPfelectron.hcalIsoDR04depth1 = pfele->dr04HcalDepth1TowerSumEt();
       MyPfelectron.hcalIsoDR04depth2 = pfele->dr04HcalDepth2TowerSumEt();
       MyPfelectron.caloIsoDR04 = pfele->dr04EcalRecHitSumEt()+pfele->dr04HcalTowerSumEt();

       MyPfelectron.mva = pfele->mva();
       MyPfelectron.mvaTrigV0 = pfele->electronID("mvaTrigV0");
       MyPfelectron.mvaNonTrigV0 = pfele->electronID("mvaNonTrigV0");
       MyPfelectron.numberOfLostHits = pfele->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
       MyPfelectron.numberOfExpectedInnerHits = pfele->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
       MyPfelectron.numberOfValidPixelHits = pfele->gsfTrack()->trackerExpectedHitsInner().numberOfValidPixelHits();
       MyPfelectron.numberOfValidPixelBarrelHits = pfele->gsfTrack()->trackerExpectedHitsInner().numberOfValidPixelBarrelHits();
       MyPfelectron.numberOfValidPixelEndcapHits = pfele->gsfTrack()->trackerExpectedHitsInner().numberOfValidPixelEndcapHits();

       MyPfelectron.scEnergy = pfele->caloEnergy();
       MyPfelectron.scSigmaEtaEta = pfele->scSigmaEtaEta();
       MyPfelectron.scSigmaIEtaIEta = pfele->scSigmaIEtaIEta();
       MyPfelectron.scE1x5 = pfele->scE1x5();
       MyPfelectron.scE2x5Max = pfele->scE2x5Max();
       MyPfelectron.scE5x5 = pfele->scE5x5();
       MyPfelectron.numClusters = pfele->numberOfBrems();

       MyPfelectron.IP = pfele->dB(pat::Electron::PV3D);
       MyPfelectron.IPError = pfele->edB(pat::Electron::PV3D);

       double caloenergy = -1;
       double pfeleEta = pfele->eta();

       if( (pfele->superCluster().isAvailable()) ){
	 double eMax    = lazyTools.eMax(    *(pfele->superCluster()) );
	 double eLeft   = lazyTools.eLeft(   *(pfele->superCluster()) );
	 double eRight  = lazyTools.eRight(  *(pfele->superCluster()) );
	 double eTop    = lazyTools.eTop(    *(pfele->superCluster()) );
	 double eBottom = lazyTools.eBottom( *(pfele->superCluster()) );
	 double e3x3    = lazyTools.e3x3(    *(pfele->superCluster()) );
	 double swissCross = -99;

	 if( eMax>0.000001 ) swissCross = 1 - (eLeft+eRight+eTop+eBottom)/eMax;

	 MyPfelectron.eMax = eMax;
	 MyPfelectron.eLeft = eLeft;
	 MyPfelectron.eRight = eRight;
	 MyPfelectron.eTop = eTop;
	 MyPfelectron.eBottom = eBottom;
	 MyPfelectron.e3x3 = e3x3;
	 MyPfelectron.swissCross = swissCross;

	 caloenergy = pfele->caloEnergy();

	 MyPfelectron.scEt = pfele->caloEnergy() * sin( pfele->superCluster()->position().theta() );
	 MyPfelectron.scRawEnergy = pfele->superCluster()->rawEnergy();
	 MyPfelectron.scEta = pfele->superCluster()->position().eta();
	 MyPfelectron.scPhi = pfele->superCluster()->position().phi();
	 MyPfelectron.scZ = pfele->superCluster()->position().Z();

	 pfeleEta = pfele->superCluster()->position().eta();

	 double seedE = -999, seedTime = -999;
	 int seedRecoFlag = -999;
       
	 if( (pfele->isEB()) ){
	   DetId idEB = EcalClusterTools::getMaximum( pfele->superCluster()->hitsAndFractions(), &(*barrelRecHits) ).first;
	   EcalRecHitCollection::const_iterator thisHitEB = barrelRecHits->find(idEB);

	   seedE = thisHitEB->energy();
	   seedTime = thisHitEB->time();
	   seedRecoFlag = thisHitEB->recoFlag();
	 }
	 else if( (pfele->isEE()) ){
	   DetId idEE = EcalClusterTools::getMaximum( pfele->superCluster()->hitsAndFractions(), &(*endcapRecHits) ).first;
	   EcalRecHitCollection::const_iterator thisHitEE = endcapRecHits->find(idEE);

	   seedE = thisHitEE->energy();
	   seedTime = thisHitEE->time();
	   seedRecoFlag = thisHitEE->recoFlag();
	 }

	 MyPfelectron.seedEnergy   = seedE;
	 MyPfelectron.seedTime     = seedTime;
	 MyPfelectron.seedRecoFlag = seedRecoFlag;

	 // MyPfelectron.swissCrossNoI85 = swissCrossNoI85;
	 // MyPfelectron.swissCrossI85   = swissCrossI85;
	 // MyPfelectron.E2overE9NoI85 = e2overe9NoI85;
	 // MyPfelectron.E2overE9I85 = e2overe9I85;
       }

       double heepEt = ( caloenergy<0 || pfele->energy()==0. ) ? 0 : pfele->et()/pfele->energy()*caloenergy;
       MyPfelectron.et = heepEt;

       if( (pfele->closestCtfTrackRef().isAvailable()) ) MyPfelectron.tkCharge = pfele->closestCtfTrackRef()->charge();
       if( (pfele->gsfTrack().isAvailable()) ){
	 double tkvx = pfele->gsfTrack()->vx();
	 double tkvy = pfele->gsfTrack()->vy();
	 double tkpx = pfele->gsfTrack()->px();
	 double tkpy = pfele->gsfTrack()->py();
	 double tkpt = pfele->gsfTrack()->pt();

	 double ndof = pfele->gsfTrack()->ndof();
	 if( (ndof!=0) ) MyPfelectron.tkNormChi2 = pfele->gsfTrack()->chi2()/ndof;
	 MyPfelectron.tkPT = pfele->gsfTrack()->pt();
	 MyPfelectron.tkEta = pfele->gsfTrack()->eta();
	 MyPfelectron.tkPhi = pfele->gsfTrack()->phi();
	 MyPfelectron.tkDZ = pfele->gsfTrack()->dz();
	 MyPfelectron.tkD0 = pfele->gsfTrack()->d0();
	 MyPfelectron.tkD0bs = (tkvx-BSx)*tkpy/tkpt-(tkvy-BSy)*tkpx/tkpt;
	 MyPfelectron.tkD0err = pfele->gsfTrack()->d0Error();
	 MyPfelectron.tkNumValidHits = pfele->gsfTrack()->numberOfValidHits();
	 MyPfelectron.gsfCharge = pfele->gsfTrack()->charge();

	 MyPfelectron.correctedD0 = pfele->gsfTrack()->dxy(beamSpotPosition);
         MyPfelectron.correctedD0Vertex = pfele->gsfTrack()->dxy(vertexPosition);
	 MyPfelectron.correctedDZ = pfele->gsfTrack()->dz(vertexPosition);
       }

       MyPfelectron.isEB = pfele->isEB();
       MyPfelectron.isEE = pfele->isEE();
       MyPfelectron.isGap = pfele->isGap();
       MyPfelectron.isEBEEGap = pfele->isEBEEGap();
       MyPfelectron.isEBGap = pfele->isEBGap();
       MyPfelectron.isEEGap = pfele->isEEGap();
       MyPfelectron.isEcalDriven = pfele->ecalDrivenSeed();
       MyPfelectron.isTrackerDriven = pfele->trackerDrivenSeed();

       if( (pfele->genLepton()) ){
	 int genId = pfele->genLepton()->pdgId();

	 MyPfelectron.genId = pfele->genLepton()->pdgId();
	 MyPfelectron.genET = pfele->genLepton()->et();
	 MyPfelectron.genPT = pfele->genLepton()->pt();
	 MyPfelectron.genPhi = pfele->genLepton()->phi();
	 MyPfelectron.genEta = pfele->genLepton()->eta();
	 MyPfelectron.genCharge = pfele->genLepton()->charge();

	 MyPfelectron.genNumberOfMothers = pfele->genLepton()->numberOfMothers();

	 if( (pfele->genLepton()->numberOfMothers()>0) ){
	   const reco::Candidate *ElectronMother = pfele->genLepton()->mother();
	   bool staytrapped = true;
	   while( (ElectronMother->pdgId()==genId && staytrapped) ){
	     if( ElectronMother->numberOfMothers()>=1 ) ElectronMother = ElectronMother->mother();
	     else staytrapped = false;
	   }
       
	   MyPfelectron.genMotherId = ElectronMother->pdgId();
	   MyPfelectron.genMotherET = ElectronMother->et();
	   MyPfelectron.genMotherPT = ElectronMother->pt();
	   MyPfelectron.genMotherPhi = ElectronMother->phi();
	   MyPfelectron.genMotherEta = ElectronMother->eta();
	   MyPfelectron.genMotherCharge = ElectronMother->charge();
	 }

	 if( (pfele->genLepton()->numberOfMothers()>=1) ){
	   const reco::Candidate *ElectronMother0 = pfele->genLepton()->mother(0);
	   const reco::Candidate *ElectronMother1 = 0;

	   if( (pfele->genLepton()->numberOfMothers()>=2) ) ElectronMother1 = pfele->genLepton()->mother(1);

	   int mother0id = ElectronMother0->pdgId();
	   int mother1id = ( ElectronMother1!=0 ) ? ElectronMother1->pdgId() : -99;

	   bool staytrapped = true;
	   while( (mother0id==genId || mother1id==genId) && staytrapped ){
	     if( mother0id==genId && (ElectronMother0!=0) ){
	       if( ElectronMother0->numberOfMothers()>=1 ){
		 ElectronMother0 = ElectronMother0->mother(0);
		 mother0id = ElectronMother0->pdgId();
		 mother1id = -99;
		 if( ElectronMother0->numberOfMothers()>=2 ){
		   ElectronMother1 = ElectronMother0->mother(1);
		   mother1id = ElectronMother1->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else if( mother1id==genId && (ElectronMother1!=0) ){
	       if( ElectronMother1->numberOfMothers()>=1 ){
		 ElectronMother1 = ElectronMother1->mother(0);
		 mother1id = ElectronMother1->pdgId();
		 mother0id = -99;
		 if( ElectronMother1->numberOfMothers()>=2 ){
		   ElectronMother0 = ElectronMother1->mother(1);
		   mother0id = ElectronMother0->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else staytrapped = false;
	   }

	   if( mother0id!=-99 ){
	     MyPfelectron.genMother0Id = ElectronMother0->pdgId();

	     if( (ElectronMother0->numberOfMothers()>=1) ){
	       const reco::Candidate *ElectronGrandMother0 = ElectronMother0->mother(0);
	       const reco::Candidate *ElectronGrandMother1 = 0;

	       if( (ElectronMother0->numberOfMothers()>=2) ) ElectronGrandMother1 = ElectronMother0->mother(1);

	       int gmother0id = ElectronGrandMother0->pdgId();
	       int gmother1id = ( ElectronGrandMother1!=0 ) ? ElectronGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother0id || gmother1id==mother0id) && staytrapped ){
		 if( gmother0id==mother0id && (ElectronGrandMother0!=0) ){
		   if( ElectronGrandMother0->numberOfMothers()>=1 ){
		     ElectronGrandMother0 = ElectronGrandMother0->mother(0);
		     gmother0id = ElectronGrandMother0->pdgId();
		     gmother1id = -99;
		     if( ElectronGrandMother0->numberOfMothers()>=2 ){
		       ElectronGrandMother1 = ElectronGrandMother0->mother(1);
		       gmother1id = ElectronGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother0id && (ElectronGrandMother1!=0) ){
		   if( ElectronGrandMother1->numberOfMothers()>=1 ){
		     ElectronGrandMother1 = ElectronGrandMother1->mother(0);
		     gmother1id = ElectronGrandMother1->pdgId();
		     gmother0id = -99;
		     if( ElectronGrandMother1->numberOfMothers()>=2 ){
		       ElectronGrandMother0 = ElectronGrandMother1->mother(1);
		       gmother0id = ElectronGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyPfelectron.genGrandMother00Id = gmother0id;
	       MyPfelectron.genGrandMother01Id = gmother1id;
	     }
	   }
	   if( mother1id!=-99 ){
	     MyPfelectron.genMother1Id = ElectronMother1->pdgId();

	     if( (ElectronMother1->numberOfMothers()>=1) ){
	       const reco::Candidate *ElectronGrandMother0 = ElectronMother1->mother(0);
	       const reco::Candidate *ElectronGrandMother1 = 0;

	       if( (ElectronMother0->numberOfMothers()>=2) ) ElectronGrandMother1 = ElectronMother1->mother(1);

	       int gmother0id = ElectronGrandMother0->pdgId();
	       int gmother1id = ( ElectronGrandMother1!=0 ) ? ElectronGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother1id || gmother1id==mother1id) && staytrapped ){
		 if( gmother0id==mother1id && (ElectronGrandMother0!=0) ){
		   if( ElectronGrandMother0->numberOfMothers()>=1 ){
		     ElectronGrandMother0 = ElectronGrandMother0->mother(0);
		     gmother0id = ElectronGrandMother0->pdgId();
		     gmother1id = -99;
		     if( ElectronGrandMother0->numberOfMothers()>=2 ){
		       ElectronGrandMother1 = ElectronGrandMother0->mother(1);
		       gmother1id = ElectronGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother1id && (ElectronGrandMother1!=0) ){
		   if( ElectronGrandMother1->numberOfMothers()>=1 ){
		     ElectronGrandMother1 = ElectronGrandMother1->mother(0);
		     gmother1id = ElectronGrandMother1->pdgId();
		     gmother0id = -99;
		     if( ElectronGrandMother1->numberOfMothers()>=2 ){
		       ElectronGrandMother0 = ElectronGrandMother1->mother(1);
		       gmother0id = ElectronGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyPfelectron.genGrandMother10Id = gmother0id;
	       MyPfelectron.genGrandMother11Id = gmother1id;
	     }
	   }
	 }

       }// end check if genLepton


       ConversionFinder convFinder;
       ConversionInfo convInfo = convFinder.getConversionInfo(*pfele, trackHandle, evt_bField);
   
       double dist = convInfo.dist();
       double dcot = convInfo.dcot();
       double convradius = convInfo.radiusOfConversion();
       math::XYZPoint convPoint = convInfo.pointOfConversion();

       MyPfelectron.dist = dist;
       MyPfelectron.dcot = dcot;
       MyPfelectron.convradius = convradius;
       MyPfelectron.convPointX = convPoint.x();
       MyPfelectron.convPointY = convPoint.y();
       MyPfelectron.convPointZ = convPoint.z();
       MyPfelectron.passConvVeto = pfele->passConversionVeto();
       

       bool cutDelEta=false, cutDelPhi=false, cutSigIeta=false, cutE2x5=false, cutEMhad1=false, cutHad2=false, cutTrackIso=false;

       double absEta = fabs(pfele->superCluster()->position().eta());
       double EMhad1 = pfele->dr03EcalRecHitSumEt() + pfele->dr03HcalDepth1TowerSumEt();
       bool EB = ( absEta < 1.442 );
       bool EE = ( absEta > 1.560 && absEta < 2.5 );

       bool isECAL = pfele->ecalDrivenSeed();

       bool cutET = ( heepEt > 25. );
       bool cutEta = ( EB || EE );
       bool cutHoverE = ( pfele->hadronicOverEm() < 0.05 );

       if( EB ){
	 cutDelEta = ( fabs(pfele->deltaEtaSuperClusterTrackAtVtx())<0.005 );
	 cutDelPhi = ( fabs(pfele->deltaPhiSuperClusterTrackAtVtx())<0.09 );
	 cutSigIeta = true;
	 cutE2x5 = ( (pfele->scE2x5Max()/pfele->scE5x5())>0.94 || (pfele->scE1x5()/pfele->scE5x5())>0.83 );
	 cutEMhad1 = ( EMhad1 < (2 + 0.03*heepEt) );
	 cutHad2 = true;
	 cutTrackIso = ( pfele->dr03TkSumPt()<7.5 );
       }
       else if( EE ){
	 cutDelEta = ( fabs(pfele->deltaEtaSuperClusterTrackAtVtx())<0.007 );
	 cutDelPhi = ( fabs(pfele->deltaPhiSuperClusterTrackAtVtx())<0.09 );
	 cutSigIeta = ( pfele->scSigmaIEtaIEta()<0.03 );
	 cutE2x5 = true;
	 if( heepEt<50 ) cutEMhad1 = ( EMhad1 < 2.5 );
	 else cutEMhad1 = ( EMhad1 < (2.5 + 0.03*(heepEt-50)) );
	 cutHad2 = ( pfele->dr03HcalDepth2TowerSumEt()<0.5 );
	 cutTrackIso = ( pfele->dr03TkSumPt()<15. );
       }

       bool isHEEPnoEt = ( isECAL && cutEta && cutDelEta && cutDelPhi && cutHoverE &&  
			   cutSigIeta && cutE2x5 && cutEMhad1 && cutHad2 && cutTrackIso );
       bool isHEEP = ( isHEEPnoEt && cutET );


       MyPfelectron.isHEEP = ( isHEEP ) ? 1 : 0;
       MyPfelectron.isHEEPnoEt = ( isHEEPnoEt ) ? 1 : 0;

       if( sample_<0 ){
	 MyPfelectron.AEffDr03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, pfeleEta, ElectronEffectiveArea::kEleEAData2012);
	 MyPfelectron.AEffDr04 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, pfeleEta, ElectronEffectiveArea::kEleEAData2012);  
       }
       else {
	 MyPfelectron.AEffDr03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, pfeleEta, ElectronEffectiveArea::kEleEAFall11MC);
	 MyPfelectron.AEffDr04 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, pfeleEta, ElectronEffectiveArea::kEleEAFall11MC); 
       }

       bnpfelectrons->push_back(MyPfelectron);

     }
   }





   /////////////////////////////////////////////
   ///////
   ///////   Fill the pfelectron collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNelectronCollection> bnpfelectronsLoose(new BNelectronCollection);
   if( producePFelectronLoose ){
     edm::View<pat::Electron> pfelectronsLoose = *pfelectronLooseHandle;

     for( edm::View<pat::Electron>::const_iterator pfele = pfelectronsLoose.begin(); pfele!=pfelectronsLoose.end(); ++pfele ){

       double pfelePin = pfele->trackMomentumAtVtx().R();
       double pfelePout = pfele->trackMomentumOut().R();

       BNelectron MyPfelectronLoose;

       // general kinematic variables
       MyPfelectronLoose.energy = pfele->energy();
       MyPfelectronLoose.gsfEt = pfele->et();
       MyPfelectronLoose.pt = pfele->pt();
       MyPfelectronLoose.px = pfele->px();
       MyPfelectronLoose.py = pfele->py();
       MyPfelectronLoose.pz = pfele->pz();
       MyPfelectronLoose.phi = pfele->phi();
       MyPfelectronLoose.eta = pfele->eta();
       MyPfelectronLoose.theta = pfele->theta();

       MyPfelectronLoose.charge = pfele->charge();
       MyPfelectronLoose.classification = pfele->classification();
       MyPfelectronLoose.vx = pfele->vx();
       MyPfelectronLoose.vy = pfele->vy();
       MyPfelectronLoose.vz = pfele->vz();

       MyPfelectronLoose.pIn = pfelePin;
       MyPfelectronLoose.pOut = pfelePout;
       MyPfelectronLoose.EscOverPin = pfele->eSuperClusterOverP();
       MyPfelectronLoose.EseedOverPout = pfele->eSeedClusterOverPout();
       MyPfelectronLoose.hadOverEm = pfele->hadronicOverEm();
       MyPfelectronLoose.fbrem = fabs(pfelePin-pfelePout)/pfelePin;
       MyPfelectronLoose.absInvEMinusInvPin = fabs( 1/pfele->ecalEnergy()-(pfele->eSuperClusterOverP()/pfele->ecalEnergy()) );
       MyPfelectronLoose.delPhiIn = pfele->deltaPhiSuperClusterTrackAtVtx();
       MyPfelectronLoose.delEtaIn = pfele->deltaEtaSuperClusterTrackAtVtx();
       /*
       MyPfelectronLoose.eidRobustHighEnergy = pfele->electronID("eidRobustHighEnergy");
       MyPfelectronLoose.eidRobustLoose = pfele->electronID("eidRobustLoose");
       MyPfelectronLoose.eidRobustTight = pfele->electronID("eidRobustTight");
       MyPfelectronLoose.eidLoose = pfele->electronID("eidLoose");
       MyPfelectronLoose.eidTight = pfele->electronID("eidTight");
       MyPfelectronLoose.eidVeryLooseMC = pfele->electronID("eidVeryLooseMC");
       MyPfelectronLoose.eidLooseMC = pfele->electronID("eidLooseMC");
       MyPfelectronLoose.eidMediumMC = pfele->electronID("eidMediumMC");
       MyPfelectronLoose.eidTightMC = pfele->electronID("eidTightMC");
       MyPfelectronLoose.eidSuperTightMC = pfele->electronID("eidSuperTightMC");
       MyPfelectronLoose.eidHyperTight1MC = pfele->electronID("eidHyperTight1MC");
       MyPfelectronLoose.eidHyperTight2MC = pfele->electronID("eidHyperTight2MC");
       MyPfelectronLoose.eidHyperTight3MC = pfele->electronID("eidHyperTight3MC");
       MyPfelectronLoose.eidHyperTight4MC = pfele->electronID("eidHyperTight4MC");
       */
       MyPfelectronLoose.particleIso = pfele->particleIso();     
       MyPfelectronLoose.chargedHadronIso = pfele->chargedHadronIso();     
       MyPfelectronLoose.neutralHadronIso = pfele->neutralHadronIso();     
       MyPfelectronLoose.photonIso = pfele->photonIso();     
       MyPfelectronLoose.puChargedHadronIso = pfele->puChargedHadronIso();     

       MyPfelectronLoose.chargedHadronIsoDR03 = pfele->userIso(0);     
       MyPfelectronLoose.neutralHadronIsoDR03 = pfele->userIso(1);
       MyPfelectronLoose.photonIsoDR03 = pfele->userIso(2);
       MyPfelectronLoose.puChargedHadronIsoDR03 = pfele->userIso(3);

       MyPfelectronLoose.chargedHadronIsoDR04 = pfele->userIso(4);     
       MyPfelectronLoose.neutralHadronIsoDR04 = pfele->userIso(5);
       MyPfelectronLoose.photonIsoDR04 = pfele->userIso(6);
       MyPfelectronLoose.puChargedHadronIsoDR04 = pfele->userIso(7);

       MyPfelectronLoose.rhoPrime = std::max(rho_event, 0.0);

       MyPfelectronLoose.trackIso = pfele->trackIso();
       MyPfelectronLoose.ecalIso = pfele->ecalIso();
       MyPfelectronLoose.hcalIso = pfele->hcalIso();
       MyPfelectronLoose.caloIso = pfele->caloIso();

       MyPfelectronLoose.trackIsoDR03 = pfele->dr03TkSumPt();
       MyPfelectronLoose.ecalIsoDR03 = pfele->dr03EcalRecHitSumEt();
       MyPfelectronLoose.hcalIsoDR03 = pfele->dr03HcalTowerSumEt();
       MyPfelectronLoose.hcalIsoDR03depth1 = pfele->dr03HcalDepth1TowerSumEt();
       MyPfelectronLoose.hcalIsoDR03depth2 = pfele->dr03HcalDepth2TowerSumEt();
       MyPfelectronLoose.caloIsoDR03 = pfele->dr03EcalRecHitSumEt()+pfele->dr03HcalTowerSumEt();

       MyPfelectronLoose.trackIsoDR04 = pfele->dr04TkSumPt();
       MyPfelectronLoose.ecalIsoDR04 = pfele->dr04EcalRecHitSumEt();
       MyPfelectronLoose.hcalIsoDR04 = pfele->dr04HcalTowerSumEt();
       MyPfelectronLoose.hcalIsoDR04depth1 = pfele->dr04HcalDepth1TowerSumEt();
       MyPfelectronLoose.hcalIsoDR04depth2 = pfele->dr04HcalDepth2TowerSumEt();
       MyPfelectronLoose.caloIsoDR04 = pfele->dr04EcalRecHitSumEt()+pfele->dr04HcalTowerSumEt();

       MyPfelectronLoose.mva = pfele->mva();
       MyPfelectronLoose.mvaTrigV0 = pfele->electronID("mvaTrigV0");
       MyPfelectronLoose.mvaNonTrigV0 = pfele->electronID("mvaNonTrigV0");
       MyPfelectronLoose.numberOfLostHits = pfele->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
       MyPfelectronLoose.numberOfExpectedInnerHits = pfele->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
       MyPfelectronLoose.numberOfValidPixelHits = pfele->gsfTrack()->trackerExpectedHitsInner().numberOfValidPixelHits();
       MyPfelectronLoose.numberOfValidPixelBarrelHits = pfele->gsfTrack()->trackerExpectedHitsInner().numberOfValidPixelBarrelHits();
       MyPfelectronLoose.numberOfValidPixelEndcapHits = pfele->gsfTrack()->trackerExpectedHitsInner().numberOfValidPixelEndcapHits();

       MyPfelectronLoose.scEnergy = pfele->caloEnergy();
       MyPfelectronLoose.scSigmaEtaEta = pfele->scSigmaEtaEta();
       MyPfelectronLoose.scSigmaIEtaIEta = pfele->scSigmaIEtaIEta();
       MyPfelectronLoose.scE1x5 = pfele->scE1x5();
       MyPfelectronLoose.scE2x5Max = pfele->scE2x5Max();
       MyPfelectronLoose.scE5x5 = pfele->scE5x5();
       MyPfelectronLoose.numClusters = pfele->numberOfBrems();

       MyPfelectronLoose.IP = pfele->dB(pat::Electron::PV3D);
       MyPfelectronLoose.IPError = pfele->edB(pat::Electron::PV3D);

       double caloenergy = -1;
       double pfeleEta = pfele->eta();

       if( (pfele->superCluster().isAvailable()) ){
	 double eMax    = lazyTools.eMax(    *(pfele->superCluster()) );
	 double eLeft   = lazyTools.eLeft(   *(pfele->superCluster()) );
	 double eRight  = lazyTools.eRight(  *(pfele->superCluster()) );
	 double eTop    = lazyTools.eTop(    *(pfele->superCluster()) );
	 double eBottom = lazyTools.eBottom( *(pfele->superCluster()) );
	 double e3x3    = lazyTools.e3x3(    *(pfele->superCluster()) );
	 double swissCross = -99;

	 if( eMax>0.000001 ) swissCross = 1 - (eLeft+eRight+eTop+eBottom)/eMax;

	 MyPfelectronLoose.eMax = eMax;
	 MyPfelectronLoose.eLeft = eLeft;
	 MyPfelectronLoose.eRight = eRight;
	 MyPfelectronLoose.eTop = eTop;
	 MyPfelectronLoose.eBottom = eBottom;
	 MyPfelectronLoose.e3x3 = e3x3;
	 MyPfelectronLoose.swissCross = swissCross;

	 caloenergy = pfele->caloEnergy();

	 MyPfelectronLoose.scEt = pfele->caloEnergy() * sin( pfele->superCluster()->position().theta() );
	 MyPfelectronLoose.scRawEnergy = pfele->superCluster()->rawEnergy();
	 MyPfelectronLoose.scEta = pfele->superCluster()->position().eta();
	 MyPfelectronLoose.scPhi = pfele->superCluster()->position().phi();
	 MyPfelectronLoose.scZ = pfele->superCluster()->position().Z();

	 pfeleEta = pfele->superCluster()->position().eta();

	 double seedE = -999, seedTime = -999;
	 int seedRecoFlag = -999;
       
	 if( (pfele->isEB()) ){
	   DetId idEB = EcalClusterTools::getMaximum( pfele->superCluster()->hitsAndFractions(), &(*barrelRecHits) ).first;
	   EcalRecHitCollection::const_iterator thisHitEB = barrelRecHits->find(idEB);

	   seedE = thisHitEB->energy();
	   seedTime = thisHitEB->time();
	   seedRecoFlag = thisHitEB->recoFlag();
	 }
	 else if( (pfele->isEE()) ){
	   DetId idEE = EcalClusterTools::getMaximum( pfele->superCluster()->hitsAndFractions(), &(*endcapRecHits) ).first;
	   EcalRecHitCollection::const_iterator thisHitEE = endcapRecHits->find(idEE);

	   seedE = thisHitEE->energy();
	   seedTime = thisHitEE->time();
	   seedRecoFlag = thisHitEE->recoFlag();
	 }

	 MyPfelectronLoose.seedEnergy   = seedE;
	 MyPfelectronLoose.seedTime     = seedTime;
	 MyPfelectronLoose.seedRecoFlag = seedRecoFlag;

	 // MyPfelectronLoose.swissCrossNoI85 = swissCrossNoI85;
	 // MyPfelectronLoose.swissCrossI85   = swissCrossI85;
	 // MyPfelectronLoose.E2overE9NoI85 = e2overe9NoI85;
	 // MyPfelectronLoose.E2overE9I85 = e2overe9I85;
       }

       double heepEt = ( caloenergy<0 || pfele->energy()==0. ) ? 0 : pfele->et()/pfele->energy()*caloenergy;
       MyPfelectronLoose.et = heepEt;

       if( (pfele->closestCtfTrackRef().isAvailable()) ) MyPfelectronLoose.tkCharge = pfele->closestCtfTrackRef()->charge();
       if( (pfele->gsfTrack().isAvailable()) ){
	 double tkvx = pfele->gsfTrack()->vx();
	 double tkvy = pfele->gsfTrack()->vy();
	 double tkpx = pfele->gsfTrack()->px();
	 double tkpy = pfele->gsfTrack()->py();
	 double tkpt = pfele->gsfTrack()->pt();

	 double ndof = pfele->gsfTrack()->ndof();
	 if( (ndof!=0) ) MyPfelectronLoose.tkNormChi2 = pfele->gsfTrack()->chi2()/ndof;
	 MyPfelectronLoose.tkPT = pfele->gsfTrack()->pt();
	 MyPfelectronLoose.tkEta = pfele->gsfTrack()->eta();
	 MyPfelectronLoose.tkPhi = pfele->gsfTrack()->phi();
	 MyPfelectronLoose.tkDZ = pfele->gsfTrack()->dz();
	 MyPfelectronLoose.tkD0 = pfele->gsfTrack()->d0();
	 MyPfelectronLoose.tkD0bs = (tkvx-BSx)*tkpy/tkpt-(tkvy-BSy)*tkpx/tkpt;
	 MyPfelectronLoose.tkD0err = pfele->gsfTrack()->d0Error();
	 MyPfelectronLoose.tkNumValidHits = pfele->gsfTrack()->numberOfValidHits();
	 MyPfelectronLoose.gsfCharge = pfele->gsfTrack()->charge();

	 MyPfelectronLoose.correctedD0 = pfele->gsfTrack()->dxy(beamSpotPosition);
         MyPfelectronLoose.correctedD0Vertex = pfele->gsfTrack()->dxy(vertexPosition);
	 MyPfelectronLoose.correctedDZ = pfele->gsfTrack()->dz(vertexPosition);
       }

       MyPfelectronLoose.isEB = pfele->isEB();
       MyPfelectronLoose.isEE = pfele->isEE();
       MyPfelectronLoose.isGap = pfele->isGap();
       MyPfelectronLoose.isEBEEGap = pfele->isEBEEGap();
       MyPfelectronLoose.isEBGap = pfele->isEBGap();
       MyPfelectronLoose.isEEGap = pfele->isEEGap();
       MyPfelectronLoose.isEcalDriven = pfele->ecalDrivenSeed();
       MyPfelectronLoose.isTrackerDriven = pfele->trackerDrivenSeed();

       if( (pfele->genLepton()) ){
	 int genId = pfele->genLepton()->pdgId();

	 MyPfelectronLoose.genId = pfele->genLepton()->pdgId();
	 MyPfelectronLoose.genET = pfele->genLepton()->et();
	 MyPfelectronLoose.genPT = pfele->genLepton()->pt();
	 MyPfelectronLoose.genPhi = pfele->genLepton()->phi();
	 MyPfelectronLoose.genEta = pfele->genLepton()->eta();
	 MyPfelectronLoose.genCharge = pfele->genLepton()->charge();

	 MyPfelectronLoose.genNumberOfMothers = pfele->genLepton()->numberOfMothers();

	 if( (pfele->genLepton()->numberOfMothers()>0) ){
	   const reco::Candidate *ElectronMother = pfele->genLepton()->mother();
	   bool staytrapped = true;
	   while( (ElectronMother->pdgId()==genId && staytrapped) ){
	     if( ElectronMother->numberOfMothers()>=1 ) ElectronMother = ElectronMother->mother();
	     else staytrapped = false;
	   }
       
	   MyPfelectronLoose.genMotherId = ElectronMother->pdgId();
	   MyPfelectronLoose.genMotherET = ElectronMother->et();
	   MyPfelectronLoose.genMotherPT = ElectronMother->pt();
	   MyPfelectronLoose.genMotherPhi = ElectronMother->phi();
	   MyPfelectronLoose.genMotherEta = ElectronMother->eta();
	   MyPfelectronLoose.genMotherCharge = ElectronMother->charge();
	 }

	 if( (pfele->genLepton()->numberOfMothers()>=1) ){
	   const reco::Candidate *ElectronMother0 = pfele->genLepton()->mother(0);
	   const reco::Candidate *ElectronMother1 = 0;

	   if( (pfele->genLepton()->numberOfMothers()>=2) ) ElectronMother1 = pfele->genLepton()->mother(1);

	   int mother0id = ElectronMother0->pdgId();
	   int mother1id = ( ElectronMother1!=0 ) ? ElectronMother1->pdgId() : -99;

	   bool staytrapped = true;
	   while( (mother0id==genId || mother1id==genId) && staytrapped ){
	     if( mother0id==genId && (ElectronMother0!=0) ){
	       if( ElectronMother0->numberOfMothers()>=1 ){
		 ElectronMother0 = ElectronMother0->mother(0);
		 mother0id = ElectronMother0->pdgId();
		 mother1id = -99;
		 if( ElectronMother0->numberOfMothers()>=2 ){
		   ElectronMother1 = ElectronMother0->mother(1);
		   mother1id = ElectronMother1->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else if( mother1id==genId && (ElectronMother1!=0) ){
	       if( ElectronMother1->numberOfMothers()>=1 ){
		 ElectronMother1 = ElectronMother1->mother(0);
		 mother1id = ElectronMother1->pdgId();
		 mother0id = -99;
		 if( ElectronMother1->numberOfMothers()>=2 ){
		   ElectronMother0 = ElectronMother1->mother(1);
		   mother0id = ElectronMother0->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else staytrapped = false;
	   }

	   if( mother0id!=-99 ){
	     MyPfelectronLoose.genMother0Id = ElectronMother0->pdgId();

	     if( (ElectronMother0->numberOfMothers()>=1) ){
	       const reco::Candidate *ElectronGrandMother0 = ElectronMother0->mother(0);
	       const reco::Candidate *ElectronGrandMother1 = 0;

	       if( (ElectronMother0->numberOfMothers()>=2) ) ElectronGrandMother1 = ElectronMother0->mother(1);

	       int gmother0id = ElectronGrandMother0->pdgId();
	       int gmother1id = ( ElectronGrandMother1!=0 ) ? ElectronGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother0id || gmother1id==mother0id) && staytrapped ){
		 if( gmother0id==mother0id && (ElectronGrandMother0!=0) ){
		   if( ElectronGrandMother0->numberOfMothers()>=1 ){
		     ElectronGrandMother0 = ElectronGrandMother0->mother(0);
		     gmother0id = ElectronGrandMother0->pdgId();
		     gmother1id = -99;
		     if( ElectronGrandMother0->numberOfMothers()>=2 ){
		       ElectronGrandMother1 = ElectronGrandMother0->mother(1);
		       gmother1id = ElectronGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother0id && (ElectronGrandMother1!=0) ){
		   if( ElectronGrandMother1->numberOfMothers()>=1 ){
		     ElectronGrandMother1 = ElectronGrandMother1->mother(0);
		     gmother1id = ElectronGrandMother1->pdgId();
		     gmother0id = -99;
		     if( ElectronGrandMother1->numberOfMothers()>=2 ){
		       ElectronGrandMother0 = ElectronGrandMother1->mother(1);
		       gmother0id = ElectronGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyPfelectronLoose.genGrandMother00Id = gmother0id;
	       MyPfelectronLoose.genGrandMother01Id = gmother1id;
	     }
	   }
	   if( mother1id!=-99 ){
	     MyPfelectronLoose.genMother1Id = ElectronMother1->pdgId();

	     if( (ElectronMother1->numberOfMothers()>=1) ){
	       const reco::Candidate *ElectronGrandMother0 = ElectronMother1->mother(0);
	       const reco::Candidate *ElectronGrandMother1 = 0;

	       if( (ElectronMother0->numberOfMothers()>=2) ) ElectronGrandMother1 = ElectronMother1->mother(1);

	       int gmother0id = ElectronGrandMother0->pdgId();
	       int gmother1id = ( ElectronGrandMother1!=0 ) ? ElectronGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother1id || gmother1id==mother1id) && staytrapped ){
		 if( gmother0id==mother1id && (ElectronGrandMother0!=0) ){
		   if( ElectronGrandMother0->numberOfMothers()>=1 ){
		     ElectronGrandMother0 = ElectronGrandMother0->mother(0);
		     gmother0id = ElectronGrandMother0->pdgId();
		     gmother1id = -99;
		     if( ElectronGrandMother0->numberOfMothers()>=2 ){
		       ElectronGrandMother1 = ElectronGrandMother0->mother(1);
		       gmother1id = ElectronGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother1id && (ElectronGrandMother1!=0) ){
		   if( ElectronGrandMother1->numberOfMothers()>=1 ){
		     ElectronGrandMother1 = ElectronGrandMother1->mother(0);
		     gmother1id = ElectronGrandMother1->pdgId();
		     gmother0id = -99;
		     if( ElectronGrandMother1->numberOfMothers()>=2 ){
		       ElectronGrandMother0 = ElectronGrandMother1->mother(1);
		       gmother0id = ElectronGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyPfelectronLoose.genGrandMother10Id = gmother0id;
	       MyPfelectronLoose.genGrandMother11Id = gmother1id;
	     }
	   }
	 }

       }// end check if genLepton


       ConversionFinder convFinder;
       ConversionInfo convInfo = convFinder.getConversionInfo(*pfele, trackHandle, evt_bField);
   
       double dist = convInfo.dist();
       double dcot = convInfo.dcot();
       double convradius = convInfo.radiusOfConversion();
       math::XYZPoint convPoint = convInfo.pointOfConversion();

       MyPfelectronLoose.dist = dist;
       MyPfelectronLoose.dcot = dcot;
       MyPfelectronLoose.convradius = convradius;
       MyPfelectronLoose.convPointX = convPoint.x();
       MyPfelectronLoose.convPointY = convPoint.y();
       MyPfelectronLoose.convPointZ = convPoint.z();
       MyPfelectronLoose.passConvVeto = pfele->passConversionVeto();
       

       bool cutDelEta=false, cutDelPhi=false, cutSigIeta=false, cutE2x5=false, cutEMhad1=false, cutHad2=false, cutTrackIso=false;

       double absEta = fabs(pfele->superCluster()->position().eta());
       double EMhad1 = pfele->dr03EcalRecHitSumEt() + pfele->dr03HcalDepth1TowerSumEt();
       bool EB = ( absEta < 1.442 );
       bool EE = ( absEta > 1.560 && absEta < 2.5 );

       bool isECAL = pfele->ecalDrivenSeed();

       bool cutET = ( heepEt > 25. );
       bool cutEta = ( EB || EE );
       bool cutHoverE = ( pfele->hadronicOverEm() < 0.05 );

       if( EB ){
	 cutDelEta = ( fabs(pfele->deltaEtaSuperClusterTrackAtVtx())<0.005 );
	 cutDelPhi = ( fabs(pfele->deltaPhiSuperClusterTrackAtVtx())<0.09 );
	 cutSigIeta = true;
	 cutE2x5 = ( (pfele->scE2x5Max()/pfele->scE5x5())>0.94 || (pfele->scE1x5()/pfele->scE5x5())>0.83 );
	 cutEMhad1 = ( EMhad1 < (2 + 0.03*heepEt) );
	 cutHad2 = true;
	 cutTrackIso = ( pfele->dr03TkSumPt()<7.5 );
       }
       else if( EE ){
	 cutDelEta = ( fabs(pfele->deltaEtaSuperClusterTrackAtVtx())<0.007 );
	 cutDelPhi = ( fabs(pfele->deltaPhiSuperClusterTrackAtVtx())<0.09 );
	 cutSigIeta = ( pfele->scSigmaIEtaIEta()<0.03 );
	 cutE2x5 = true;
	 if( heepEt<50 ) cutEMhad1 = ( EMhad1 < 2.5 );
	 else cutEMhad1 = ( EMhad1 < (2.5 + 0.03*(heepEt-50)) );
	 cutHad2 = ( pfele->dr03HcalDepth2TowerSumEt()<0.5 );
	 cutTrackIso = ( pfele->dr03TkSumPt()<15. );
       }

       bool isHEEPnoEt = ( isECAL && cutEta && cutDelEta && cutDelPhi && cutHoverE &&  
			   cutSigIeta && cutE2x5 && cutEMhad1 && cutHad2 && cutTrackIso );
       bool isHEEP = ( isHEEPnoEt && cutET );


       MyPfelectronLoose.isHEEP = ( isHEEP ) ? 1 : 0;
       MyPfelectronLoose.isHEEPnoEt = ( isHEEPnoEt ) ? 1 : 0;

       if( sample_<0 ){
	 MyPfelectronLoose.AEffDr03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, pfeleEta, ElectronEffectiveArea::kEleEAData2012);
	 MyPfelectronLoose.AEffDr04 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, pfeleEta, ElectronEffectiveArea::kEleEAData2012);  
       }
       else {
	 MyPfelectronLoose.AEffDr03 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, pfeleEta, ElectronEffectiveArea::kEleEAFall11MC);
	 MyPfelectronLoose.AEffDr04 = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04, pfeleEta, ElectronEffectiveArea::kEleEAFall11MC); 
       }

       bnpfelectronsLoose->push_back(MyPfelectronLoose);

     }
   }



   // double metJES1corrPx = 0;
   // double metJES1corrPy = 0;
   // double metJES20corrPx = 0;
   // double metJES20corrPy = 0;
   // double UDeltaPx = 0;
   // double UDeltaPy = 0;
   // double USumET = 0;



   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl_Calo;
   iSetup.get<JetCorrectionsRecord>().get("AK5Calo",JetCorParColl_Calo); 
   JetCorrectorParameters const & JetCorPar_Calo = (*JetCorParColl_Calo)["Uncertainty"];
   JetCorrectionUncertainty *jecUnc_Calo = new JetCorrectionUncertainty(JetCorPar_Calo);



   /////////////////////////////////////////////
   ///////
   ///////   Fill the calojet collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNjetCollection> bncalojets(new BNjetCollection);
   if( produceCaloJet ){
     edm::View<pat::Jet> calojets = *calojetHandle;

     for( edm::View<pat::Jet>::const_iterator calojet = calojets.begin(); calojet != calojets.end(); ++ calojet ) {

       if( !(calojet->pt()>minJetPt_) ) continue;

       double unc = 1., JECuncUp = 1., JECuncDown = 1.; // JEC uncertainties only defined for jets with |eta| < 5.5 and pt > 9 GeV (2011 data)
       if( calojet->pt()>9. && fabs(calojet->eta())<5.0 ){
	 jecUnc_Calo->setJetEta(calojet->eta());
	 jecUnc_Calo->setJetPt(calojet->pt());// the uncertainty is a function of the corrected pt
	 JECuncUp = jecUnc_Calo->getUncertainty(true); //up variation
	 unc = JECuncUp;
	 jecUnc_Calo->setJetEta(calojet->eta());
	 jecUnc_Calo->setJetPt(calojet->pt());// the uncertainty is a function of the corrected pt
	 JECuncDown = jecUnc_Calo->getUncertainty(false); //up variation
       }

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
       MyCalojet.btagCombinedSecVertex = calojet->bDiscriminator("combinedSecondaryVertexBJetTags");
       MyCalojet.btagCombinedSecVertexMVA = calojet->bDiscriminator("combinedSecondaryVertexMVABJetTags");
       MyCalojet.btagSoftMuonByPt = calojet->bDiscriminator("softMuonByPtBJetTags");
       MyCalojet.btagSoftMuonByIP3 = calojet->bDiscriminator("softMuonByIP3dBJetTags");
       MyCalojet.btagSoftElectronByPt = calojet->bDiscriminator("softElectronByPtBJetTags");
       MyCalojet.btagSoftElectronByIP3 = calojet->bDiscriminator("softElectronByIP3dBJetTags");


       MyCalojet.JESunc = unc;
       MyCalojet.JECuncUp = JECuncUp;
       MyCalojet.JECuncDown = JECuncDown;

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

	 int numberOfMothers = calojet->genParton()->numberOfMothers();
	 if( numberOfMothers==1 ){
	   MyCalojet.genPartonMotherId  = calojet->genParton()->mother()->pdgId();
	   MyCalojet.genPartonMother0Id = calojet->genParton()->mother()->pdgId();
	   int numberOfGrandMothers = calojet->genParton()->mother()->numberOfMothers();
	   if( numberOfGrandMothers==1 ) MyCalojet.genPartonGrandMotherId = calojet->genParton()->mother()->mother()->pdgId();
	   else if( numberOfGrandMothers>=2 ){
	     MyCalojet.genPartonGrandMother00Id = calojet->genParton()->mother()->mother(0)->pdgId();
	     MyCalojet.genPartonGrandMother01Id = calojet->genParton()->mother()->mother(1)->pdgId();
	   }

	   int pdgId = calojet->genParton()->pdgId();
	   int motherId = calojet->genParton()->mother()->pdgId();

	   int last_motherID = -99;
	   int last_grandMotherID = -99;
	   if( pdgId==motherId ){
	     const reco::Candidate* new_mother = calojet->genParton()->mother();

	     bool keepGoing = true;
	     while( keepGoing ){
	       int new_motherID = new_mother->pdgId();
	       last_motherID = new_motherID;
	       if( new_mother->numberOfMothers()>0 ){
		 new_mother = new_mother->mother();
		 last_grandMotherID = new_mother->pdgId();
		 if( new_motherID!=pdgId ) keepGoing = false;
	       }
	       else keepGoing = false;
	     }
	   }

	   if( last_motherID!=-99 ){
	     MyCalojet.genPartonMother0Id = last_motherID;
	     if( last_grandMotherID!=-99 ) MyCalojet.genPartonGrandMotherId = last_grandMotherID;
	   }
	 }
	 else if( numberOfMothers>=2 ){
	   MyCalojet.genPartonMother0Id = calojet->genParton()->mother(0)->pdgId();
	   MyCalojet.genPartonMother1Id = calojet->genParton()->mother(1)->pdgId();

	   if( calojet->genParton()->mother(0)->numberOfMothers()==1 ) MyCalojet.genPartonGrandMother00Id = calojet->genParton()->mother(0)->mother()->pdgId();
	   else if( calojet->genParton()->mother(0)->numberOfMothers()>=2 ){
	     MyCalojet.genPartonGrandMother00Id = calojet->genParton()->mother(0)->mother(0)->pdgId();
	     MyCalojet.genPartonGrandMother01Id = calojet->genParton()->mother(0)->mother(1)->pdgId();
	   }

	   if( calojet->genParton()->mother(1)->numberOfMothers()==1 ) MyCalojet.genPartonGrandMother00Id = calojet->genParton()->mother(1)->mother()->pdgId();
	   else if( calojet->genParton()->mother(1)->numberOfMothers()>=2 ){
	     MyCalojet.genPartonGrandMother10Id = calojet->genParton()->mother(1)->mother(0)->pdgId();
	     MyCalojet.genPartonGrandMother11Id = calojet->genParton()->mother(1)->mother(1)->pdgId();
	   }
	 }
       }

       bncalojets->push_back(MyCalojet);
     }
   }



   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl_PF;
   iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl_PF); 
   JetCorrectorParameters const & JetCorPar_PF = (*JetCorParColl_PF)["Uncertainty"];
   JetCorrectionUncertainty *jecUnc_PF = new JetCorrectionUncertainty(JetCorPar_PF);



   /////////////////////////////////////////////
   ///////
   ///////   Fill the pfjet collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNjetCollection> bnpfjets(new BNjetCollection);
   if( producePFJet ){
     edm::View<pat::Jet> pfjets = *pfjetHandle;
     std::vector<PFCandidatePtr> PFJetPart;

     //full
     Handle<ValueMap<float> > puJetIdMVA_full;
     iEvent.getByLabel("puJetMvaChs","fullDiscriminant",puJetIdMVA_full);
     Handle<ValueMap<int> > puJetIdFlag_full;
     iEvent.getByLabel("puJetMvaChs","fullId",puJetIdFlag_full);

     //simple
     Handle<ValueMap<float> > puJetIdMVA_simple;
     iEvent.getByLabel("puJetMvaChs","simpleDiscriminant",puJetIdMVA_simple);
     Handle<ValueMap<int> > puJetIdFlag_simple;
     iEvent.getByLabel("puJetMvaChs","simpleId",puJetIdFlag_simple);

     //cutbased
     Handle<ValueMap<float> > puJetIdMVA_cutbased;
     iEvent.getByLabel("puJetMvaChs","cutbasedDiscriminant",puJetIdMVA_cutbased);
     Handle<ValueMap<int> > puJetIdFlag_cutbased;
     iEvent.getByLabel("puJetMvaChs","cutbasedId",puJetIdFlag_cutbased);

     //puJet ID  MVA inputs
     Handle<ValueMap<StoredPileupJetIdentifier> > vmap;
     iEvent.getByLabel("puJetIdChs", vmap);

     PileupJetIdentifier puIdentifier;

     for( edm::View<pat::Jet>::const_iterator pfjet = pfjets.begin(); pfjet != pfjets.end(); ++ pfjet ) {

       if( !(pfjet->pt()>minJetPt_) ) continue;

       BNjet MyPfjet;

       unsigned int idx = pfjet-pfjets.begin();

       float mva_full  = (*puJetIdMVA_full)[pfjets.refAt(idx)];
       int idflag_full = (*puJetIdFlag_full)[pfjets.refAt(idx)];

       float mva_simple  = (*puJetIdMVA_simple)[pfjets.refAt(idx)];
       int idflag_simple = (*puJetIdFlag_simple)[pfjets.refAt(idx)];

       float mva_cutbased  = (*puJetIdMVA_cutbased)[pfjets.refAt(idx)];
       int idflag_cutbased = (*puJetIdFlag_cutbased)[pfjets.refAt(idx)];

       puIdentifier = (*vmap)[pfjets.refAt(idx)];

       bool passTight_full=false, passMedium_full=false, passLoose_full=false;
       if( ( ( idflag_full & (1 << 0) ) != 0 ) ) passTight_full  = true;
       if( ( ( idflag_full & (1 << 1) ) != 0 ) ) passMedium_full = true;
       if( ( ( idflag_full & (1 << 2) ) != 0 ) ) passLoose_full = true;

       bool passTight_simple=false, passMedium_simple=false, passLoose_simple=false;
       if( ( ( idflag_simple & (1 << 0) ) != 0 ) ) passTight_simple  = true;
       if( ( ( idflag_simple & (1 << 1) ) != 0 ) ) passMedium_simple = true;
       if( ( ( idflag_simple & (1 << 2) ) != 0 ) ) passLoose_simple = true;

       bool passTight_cutbased=false, passMedium_cutbased=false, passLoose_cutbased=false;
       if( ( ( idflag_cutbased & (1 << 0) ) != 0 ) ) passTight_cutbased  = true;
       if( ( ( idflag_cutbased & (1 << 1) ) != 0 ) ) passMedium_cutbased = true;
       if( ( ( idflag_cutbased & (1 << 2) ) != 0 ) ) passLoose_cutbased = true;

       // printf("  ===> jet %d,\t pt = %4.1f,\t eta = %4.2f \n", idx, pfjet->pt(), pfjet->eta() );
       // printf("\t\t\t mva_full = %4.3f,\t idflag_full = %d,\t passTight_full = %d,\t passMedium_full = %d,\t passLoose_full = %d \n",
       // 	      mva_full, idflag_full, (passTight_full)?1:0, (passMedium_full)?1:0, (passLoose_full)?1:0 );
       // printf("\t\t\t mva_simp = %4.3f,\t idflag_simp = %d,\t passTight_simp = %d,\t passMedium_simp = %d,\t passLoose_simp = %d \n",
       // 	      mva_simple, idflag_simple, (passTight_simple)?1:0, (passMedium_simple)?1:0, (passLoose_simple)?1:0 );
       // printf("\t\t\t mva_cutb = %4.3f,\t idflag_cutb = %d,\t passTight_cutb = %d,\t passMedium_cutb = %d,\t passLoose_cutb = %d \n",
       // 	      mva_cutbased, idflag_cutbased, (passTight_cutbased)?1:0, (passMedium_cutbased)?1:0, (passLoose_cutbased)?1:0 );

       MyPfjet.puJetMVA_full     = mva_full;
       MyPfjet.puJetMVA_simple   = mva_simple;
       MyPfjet.puJetMVA_cutbased = mva_cutbased;

       MyPfjet.puJetId_full     = idflag_full;
       MyPfjet.puJetId_simple   = idflag_simple;
       MyPfjet.puJetId_cutbased = idflag_cutbased;

       MyPfjet.puJetId_tight_full     = (passTight_full)?1:0;
       MyPfjet.puJetId_tight_simple   = (passTight_simple)?1:0;
       MyPfjet.puJetId_tight_cutbased = (passTight_cutbased)?1:0;

       MyPfjet.puJetId_medium_full     = (passMedium_full)?1:0;
       MyPfjet.puJetId_medium_simple   = (passMedium_simple)?1:0;
       MyPfjet.puJetId_medium_cutbased = (passMedium_cutbased)?1:0;

       MyPfjet.puJetId_loose_full     = (passLoose_full)?1:0;
       MyPfjet.puJetId_loose_simple   = (passLoose_simple)?1:0;
       MyPfjet.puJetId_loose_cutbased = (passLoose_cutbased)?1:0;



       double rawpt = pfjet->correctedJet(0).pt();

       PFJetPart = pfjet->getPFConstituents();

       double maxCandPt=0;
       double leadCandVx=-99,leadCandVy=-99,leadCandVz=-99;
       for(UInt_t j=0;j<PFJetPart.size();j++){
	 double pTcand = PFJetPart[j]->pt();
	 if( pTcand>maxCandPt ){
	   maxCandPt = pTcand;
	   leadCandVx = PFJetPart[j]->vx();
	   leadCandVy = PFJetPart[j]->vy();
	   leadCandVz = PFJetPart[j]->vz();
	 }
       }

       double leadCandDistFromPV = sqrt( (leadCandVx-PVx)*(leadCandVx-PVx) + (leadCandVy-PVy)*(leadCandVy-PVy) + (leadCandVz-PVz)*(leadCandVz-PVz) );

       MyPfjet.leadCandPt = maxCandPt;

       MyPfjet.leadCandVx = leadCandVx;
       MyPfjet.leadCandVy = leadCandVy;
       MyPfjet.leadCandVz = leadCandVz;
       MyPfjet.leadCandDistFromPV = leadCandDistFromPV;


       MyPfjet.dZ = puIdentifier.dZ();
       MyPfjet.dR2Mean = puIdentifier.dR2Mean();
       MyPfjet.dRMean = puIdentifier.dRMean();
       MyPfjet.frac01 = puIdentifier.frac01();
       MyPfjet.frac02 = puIdentifier.frac02();
       MyPfjet.frac03 = puIdentifier.frac03();
       MyPfjet.frac04 = puIdentifier.frac04();
       MyPfjet.frac05 = puIdentifier.frac05();
       MyPfjet.frac06 = puIdentifier.frac06();
       MyPfjet.frac07 = puIdentifier.frac07(); //Always 0
       MyPfjet.beta = puIdentifier.beta();
       MyPfjet.betaStar = puIdentifier.betaStar();
       MyPfjet.betaClassic = puIdentifier.betaClassic();
       MyPfjet.betaStarClassic = puIdentifier.betaStarClassic();
       MyPfjet.ptD = puIdentifier.ptD();
       MyPfjet.nvtx = puIdentifier.nvtx();
       MyPfjet.d0 = puIdentifier.d0(); //Not declared in class StoredPileupJetIdentifier (PileupJetIdentifier.h) 

       MyPfjet.Upt = rawpt;

       double unc = 1., JECuncUp = 1., JECuncDown = 1.; // JEC uncertainties only defined for jets with |eta| < 5.5 and pt > 9 GeV (2011 data)
       if( pfjet->pt()>9. && fabs(pfjet->eta())<5.0 ){
	 jecUnc_PF->setJetEta(pfjet->eta());
	 jecUnc_PF->setJetPt(pfjet->pt());// the uncertainty is a function of the corrected pt
	 JECuncUp = jecUnc_PF->getUncertainty(true); //up variation
	 unc = JECuncUp;
	 jecUnc_PF->setJetEta(pfjet->eta());
	 jecUnc_PF->setJetPt(pfjet->pt());// the uncertainty is a function of the corrected pt
	 JECuncDown = jecUnc_PF->getUncertainty(false); //up variation
       }

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
       MyPfjet.btagCombinedSecVertex = pfjet->bDiscriminator("combinedSecondaryVertexBJetTags");
       MyPfjet.btagCombinedSecVertexMVA = pfjet->bDiscriminator("combinedSecondaryVertexMVABJetTags");
       MyPfjet.btagSoftMuonByPt = pfjet->bDiscriminator("softMuonByPtBJetTags");
       MyPfjet.btagSoftMuonByIP3 = pfjet->bDiscriminator("softMuonByIP3dBJetTags");
       MyPfjet.btagSoftElectronByPt = pfjet->bDiscriminator("softElectronByPtBJetTags");
       MyPfjet.btagSoftElectronByIP3 = pfjet->bDiscriminator("softElectronByIP3dBJetTags");

       MyPfjet.JESunc = unc;
       MyPfjet.JECuncUp = JECuncUp;
       MyPfjet.JECuncDown = JECuncDown;

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

	 int numberOfMothers = pfjet->genParton()->numberOfMothers();
	 if( numberOfMothers==1 ){
	   MyPfjet.genPartonMotherId  = pfjet->genParton()->mother()->pdgId();
	   MyPfjet.genPartonMother0Id = pfjet->genParton()->mother()->pdgId();
	   int numberOfGrandMothers = pfjet->genParton()->mother()->numberOfMothers();
	   if( numberOfGrandMothers==1 ) MyPfjet.genPartonGrandMotherId = pfjet->genParton()->mother()->mother()->pdgId();
	   else if( numberOfGrandMothers>=2 ){
	     MyPfjet.genPartonGrandMother00Id = pfjet->genParton()->mother()->mother(0)->pdgId();
	     MyPfjet.genPartonGrandMother01Id = pfjet->genParton()->mother()->mother(1)->pdgId();
	   }

	   int pdgId = pfjet->genParton()->pdgId();
	   int motherId = pfjet->genParton()->mother()->pdgId();

	   int last_motherID = -99;
	   int last_grandMotherID = -99;
	   if( pdgId==motherId ){
	     const reco::Candidate* new_mother = pfjet->genParton()->mother();

	     bool keepGoing = true;
	     while( keepGoing ){
	       int new_motherID = new_mother->pdgId();
	       last_motherID = new_motherID;
	       if( new_mother->numberOfMothers()>0 ){
		 new_mother = new_mother->mother();
		 last_grandMotherID = new_mother->pdgId();
		 if( new_motherID!=pdgId ) keepGoing = false;
	       }
	       else keepGoing = false;
	     }
	   }

	   if( last_motherID!=-99 ){
	     MyPfjet.genPartonMother0Id = last_motherID;
	     if( last_grandMotherID!=-99 ) MyPfjet.genPartonGrandMotherId = last_grandMotherID;
	   }
	 }
	 else if( numberOfMothers>=2 ){
	   MyPfjet.genPartonMother0Id = pfjet->genParton()->mother(0)->pdgId();
	   MyPfjet.genPartonMother1Id = pfjet->genParton()->mother(1)->pdgId();

	   if( pfjet->genParton()->mother(0)->numberOfMothers()==1 ) MyPfjet.genPartonGrandMother00Id = pfjet->genParton()->mother(0)->mother()->pdgId();
	   else if( pfjet->genParton()->mother(0)->numberOfMothers()>=2 ){
	     MyPfjet.genPartonGrandMother00Id = pfjet->genParton()->mother(0)->mother(0)->pdgId();
	     MyPfjet.genPartonGrandMother01Id = pfjet->genParton()->mother(0)->mother(1)->pdgId();
	   }

	   if( pfjet->genParton()->mother(1)->numberOfMothers()==1 ) MyPfjet.genPartonGrandMother00Id = pfjet->genParton()->mother(1)->mother()->pdgId();
	   else if( pfjet->genParton()->mother(1)->numberOfMothers()>=2 ){
	     MyPfjet.genPartonGrandMother10Id = pfjet->genParton()->mother(1)->mother(0)->pdgId();
	     MyPfjet.genPartonGrandMother11Id = pfjet->genParton()->mother(1)->mother(1)->pdgId();
	   }
	 }

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
   }




   /////////////////////////////////////////////
   ///////
   ///////   Fill the muon collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNmuonCollection> bnmuons(new BNmuonCollection);
   if( produceMuon ){
     edm::View<pat::Muon> muons = *muonHandle;

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
     
       MyMuon.timeAtIpInOut = muon->time().timeAtIpInOut;
       MyMuon.timeAtIpInOutErr = muon->time().timeAtIpInOutErr;
       MyMuon.timeAtIpOutIn = muon->time().timeAtIpOutIn;
       MyMuon.timeAtIpOutInErr = muon->time().timeAtIpOutInErr;
       MyMuon.time_ndof = muon->time().nDof;

       if( muon->isEnergyValid() ){
	 MyMuon.ecal_time = muon->calEnergy().ecal_time;
	 MyMuon.hcal_time = muon->calEnergy().hcal_time;
	 MyMuon.ecal_timeError = muon->calEnergy().ecal_timeError;
	 MyMuon.hcal_timeError = muon->calEnergy().hcal_timeError;
	 MyMuon.energy_ecal = muon->calEnergy().em;
	 MyMuon.energy_hcal = muon->calEnergy().had;
	 MyMuon.e3x3_ecal = muon->calEnergy().emS9;
	 MyMuon.e3x3_hcal = muon->calEnergy().hadS9;
	 MyMuon.energyMax_ecal = muon->calEnergy().emMax;
	 MyMuon.energyMax_hcal = muon->calEnergy().hadMax;
       }



       MyMuon.IDGMPTight = ( muon->isGood("GlobalMuonPromptTight") ) ? 1 : 0;

       MyMuon.isGlobalMuon = ( muon->isGlobalMuon() ) ? 1 : 0;
       MyMuon.isTrackerMuon = ( muon->isTrackerMuon() ) ? 1 : 0;
       MyMuon.isStandAloneMuon = ( muon->isStandAloneMuon() ) ? 1 : 0;
       MyMuon.isGlobalMuonPromptTight = ( muon->isGood("GlobalMuonPromptTight") ) ? 1 : 0;

       MyMuon.numberOfMatches = muon->numberOfMatches();
       MyMuon.numberOfMatchedStations = muon->numberOfMatchedStations();

       MyMuon.dVzPVz = muon->vz() - PVz;
       MyMuon.dB = muon->dB();

       MyMuon.IP = muon->dB(pat::Muon::PV3D);
       MyMuon.IPError = muon->edB(pat::Muon::PV3D);


       if( (muon->globalTrack().isAvailable()) ){
	 MyMuon.normalizedChi2 = muon->globalTrack()->normalizedChi2();
	 MyMuon.numberOfValidMuonHits = muon->globalTrack()->hitPattern().numberOfValidMuonHits();
	 MyMuon.numberOfValidTrackerHits = muon->globalTrack()->hitPattern().numberOfValidTrackerHits();
         MyMuon.numberOfLayersWithMeasurement = muon->track()->hitPattern().trackerLayersWithMeasurement();
	 MyMuon.ptErr = muon->globalTrack()->ptError();
       }
       if( (muon->innerTrack().isAvailable()) ){
	 MyMuon.numberOfValidTrackerHitsInnerTrack = muon->innerTrack()->numberOfValidHits();
	 MyMuon.pixelLayersWithMeasurement = muon->innerTrack()->hitPattern().pixelLayersWithMeasurement();
         MyMuon.numberOfValidPixelHits = muon->innerTrack()->hitPattern().numberOfValidPixelHits();

	 MyMuon.correctedD0 = muon->innerTrack()->dxy(beamSpotPosition);
         MyMuon.correctedD0Vertex = muon->innerTrack()->dxy(vertexPosition);
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
	 int genId = muon->genLepton()->pdgId();

	 MyMuon.genId = muon->genLepton()->pdgId();
	 MyMuon.genET = muon->genLepton()->et();
	 MyMuon.genPT = muon->genLepton()->pt();
	 MyMuon.genPhi = muon->genLepton()->phi();
	 MyMuon.genEta = muon->genLepton()->eta();
	 MyMuon.genCharge = muon->genLepton()->charge();

	 MyMuon.genNumberOfMothers = muon->genLepton()->numberOfMothers();

	 if( (muon->genLepton()->numberOfMothers()>0) ){
	   const reco::Candidate *MuonMother = muon->genLepton()->mother();
	   bool staytrapped = true;
	   while( (MuonMother->pdgId()==genId && staytrapped) ){
	     if( MuonMother->numberOfMothers()>=1 ) MuonMother = MuonMother->mother();
	     else staytrapped = false;
	   }
       
	   MyMuon.genMotherId = MuonMother->pdgId();
	   MyMuon.genMotherET = MuonMother->et();
	   MyMuon.genMotherPT = MuonMother->pt();
	   MyMuon.genMotherPhi = MuonMother->phi();
	   MyMuon.genMotherEta = MuonMother->eta();
	   MyMuon.genMotherCharge = MuonMother->charge();
	 }

	 if( (muon->genLepton()->numberOfMothers()>=1) ){
	   const reco::Candidate *MuonMother0 = muon->genLepton()->mother(0);
	   const reco::Candidate *MuonMother1 = 0;

	   if( (muon->genLepton()->numberOfMothers()>=2) ) MuonMother1 = muon->genLepton()->mother(1);

	   int mother0id = MuonMother0->pdgId();
	   int mother1id = ( MuonMother1!=0 ) ? MuonMother1->pdgId() : -99;

	   bool staytrapped = true;
	   while( (mother0id==genId || mother1id==genId) && staytrapped ){
	     if( mother0id==genId && (MuonMother0!=0) ){
	       if( MuonMother0->numberOfMothers()>=1 ){
		 MuonMother0 = MuonMother0->mother(0);
		 mother0id = MuonMother0->pdgId();
		 mother1id = -99;
		 if( MuonMother0->numberOfMothers()>=2 ){
		   MuonMother1 = MuonMother0->mother(1);
		   mother1id = MuonMother1->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else if( mother1id==genId && (MuonMother1!=0) ){
	       if( MuonMother1->numberOfMothers()>=1 ){
		 MuonMother1 = MuonMother1->mother(0);
		 mother1id = MuonMother1->pdgId();
		 mother0id = -99;
		 if( MuonMother1->numberOfMothers()>=2 ){
		   MuonMother0 = MuonMother1->mother(1);
		   mother0id = MuonMother0->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else staytrapped = false;
	   }

	   if( mother0id!=-99 ){
	     MyMuon.genMother0Id = MuonMother0->pdgId();

	     if( (MuonMother0->numberOfMothers()>=1) ){
	       const reco::Candidate *MuonGrandMother0 = MuonMother0->mother(0);
	       const reco::Candidate *MuonGrandMother1 = 0;

	       if( (MuonMother0->numberOfMothers()>=2) ) MuonGrandMother1 = MuonMother0->mother(1);

	       int gmother0id = MuonGrandMother0->pdgId();
	       int gmother1id = ( MuonGrandMother1!=0 ) ? MuonGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother0id || gmother1id==mother0id) && staytrapped ){
		 if( gmother0id==mother0id && (MuonGrandMother0!=0) ){
		   if( MuonGrandMother0->numberOfMothers()>=1 ){
		     MuonGrandMother0 = MuonGrandMother0->mother(0);
		     gmother0id = MuonGrandMother0->pdgId();
		     gmother1id = -99;
		     if( MuonGrandMother0->numberOfMothers()>=2 ){
		       MuonGrandMother1 = MuonGrandMother0->mother(1);
		       gmother1id = MuonGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother0id && (MuonGrandMother1!=0) ){
		   if( MuonGrandMother1->numberOfMothers()>=1 ){
		     MuonGrandMother1 = MuonGrandMother1->mother(0);
		     gmother1id = MuonGrandMother1->pdgId();
		     gmother0id = -99;
		     if( MuonGrandMother1->numberOfMothers()>=2 ){
		       MuonGrandMother0 = MuonGrandMother1->mother(1);
		       gmother0id = MuonGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyMuon.genGrandMother00Id = gmother0id;
	       MyMuon.genGrandMother01Id = gmother1id;
	     }
	   }
	   if( mother1id!=-99 ){
	     MyMuon.genMother1Id = MuonMother1->pdgId();

	     if( (MuonMother1->numberOfMothers()>=1) ){
	       const reco::Candidate *MuonGrandMother0 = MuonMother1->mother(0);
	       const reco::Candidate *MuonGrandMother1 = 0;

	       if( (MuonMother0->numberOfMothers()>=2) ) MuonGrandMother1 = MuonMother1->mother(1);

	       int gmother0id = MuonGrandMother0->pdgId();
	       int gmother1id = ( MuonGrandMother1!=0 ) ? MuonGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother1id || gmother1id==mother1id) && staytrapped ){
		 if( gmother0id==mother1id && (MuonGrandMother0!=0) ){
		   if( MuonGrandMother0->numberOfMothers()>=1 ){
		     MuonGrandMother0 = MuonGrandMother0->mother(0);
		     gmother0id = MuonGrandMother0->pdgId();
		     gmother1id = -99;
		     if( MuonGrandMother0->numberOfMothers()>=2 ){
		       MuonGrandMother1 = MuonGrandMother0->mother(1);
		       gmother1id = MuonGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother1id && (MuonGrandMother1!=0) ){
		   if( MuonGrandMother1->numberOfMothers()>=1 ){
		     MuonGrandMother1 = MuonGrandMother1->mother(0);
		     gmother1id = MuonGrandMother1->pdgId();
		     gmother0id = -99;
		     if( MuonGrandMother1->numberOfMothers()>=2 ){
		       MuonGrandMother0 = MuonGrandMother1->mother(1);
		       gmother0id = MuonGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyMuon.genGrandMother10Id = gmother0id;
	       MyMuon.genGrandMother11Id = gmother1id;
	     }
	   }
	 }

       }// end check if genLepton

       bnmuons->push_back(MyMuon);
     }
   }


   /////////////////////////////////////////////
   ///////
   ///////   Fill the pfmuon collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNmuonCollection> bnpfmuons(new BNmuonCollection);
   if( producePFMuon ){
     edm::View<pat::Muon> pfmuons = *pfmuonHandle;

     for( edm::View<pat::Muon>::const_iterator pfmuon = pfmuons.begin(); pfmuon!=pfmuons.end(); ++pfmuon ){

       BNmuon MyPfmuon;

       // general kinematic variables
       MyPfmuon.energy = pfmuon->energy();
       MyPfmuon.et = pfmuon->et();
       MyPfmuon.pt = pfmuon->pt();
       MyPfmuon.px = pfmuon->px();
       MyPfmuon.py = pfmuon->py();
       MyPfmuon.pz = pfmuon->pz();
       MyPfmuon.phi = pfmuon->phi();
       MyPfmuon.eta = pfmuon->eta();
       MyPfmuon.theta = pfmuon->theta();


       MyPfmuon.charge = pfmuon->charge();
       MyPfmuon.vx = pfmuon->vx();
       MyPfmuon.vy = pfmuon->vy();
       MyPfmuon.vz = pfmuon->vz();


       MyPfmuon.particleIso = pfmuon->particleIso();     
       MyPfmuon.chargedHadronIso = pfmuon->chargedHadronIso();     
       MyPfmuon.neutralHadronIso = pfmuon->neutralHadronIso();     
       MyPfmuon.photonIso = pfmuon->photonIso();     
       MyPfmuon.puChargedHadronIso = pfmuon->puChargedHadronIso();     

       MyPfmuon.chargedHadronIsoDR03 = pfmuon->userIso(0);     
       MyPfmuon.neutralHadronIsoDR03 = pfmuon->userIso(1);
       MyPfmuon.photonIsoDR03 = pfmuon->userIso(2);
       MyPfmuon.puChargedHadronIsoDR03 = pfmuon->userIso(3);

       MyPfmuon.chargedHadronIsoDR04 = pfmuon->userIso(4);     
       MyPfmuon.neutralHadronIsoDR04 = pfmuon->userIso(5);
       MyPfmuon.photonIsoDR04 = pfmuon->userIso(6);
       MyPfmuon.puChargedHadronIsoDR04 = pfmuon->userIso(7);

       MyPfmuon.rhoPrime = std::max(rho_event, 0.0);
       if( sample_<0 ){
	 MyPfmuon.AEffDr03 = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaAndNeutralHadronIso03, pfmuon->eta(), MuonEffectiveArea::kMuEAData2012);
	 MyPfmuon.AEffDr04 = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaAndNeutralHadronIso04, pfmuon->eta(), MuonEffectiveArea::kMuEAData2012);
       }
       else {
	 MyPfmuon.AEffDr03 = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaAndNeutralHadronIso03, pfmuon->eta(), MuonEffectiveArea::kMuEAFall11MC);
	 MyPfmuon.AEffDr04 = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaAndNeutralHadronIso04, pfmuon->eta(), MuonEffectiveArea::kMuEAFall11MC);
       }

       MyPfmuon.pfIsoR03SumChargedHadronPt = pfmuon->pfIsolationR03().sumChargedHadronPt;
       MyPfmuon.pfIsoR03SumNeutralHadronEt = pfmuon->pfIsolationR03().sumNeutralHadronEt;
       MyPfmuon.pfIsoR03SumPhotonEt = pfmuon->pfIsolationR03().sumPhotonEt;
       MyPfmuon.pfIsoR03SumPUPt = pfmuon->pfIsolationR03().sumPUPt;

       MyPfmuon.pfIsoR04SumChargedHadronPt = pfmuon->pfIsolationR04().sumChargedHadronPt;
       MyPfmuon.pfIsoR04SumNeutralHadronEt = pfmuon->pfIsolationR04().sumNeutralHadronEt;
       MyPfmuon.pfIsoR04SumPhotonEt = pfmuon->pfIsolationR04().sumPhotonEt;
       MyPfmuon.pfIsoR04SumPUPt = pfmuon->pfIsolationR04().sumPUPt;



       MyPfmuon.trackIso = pfmuon->trackIso();
       MyPfmuon.ecalIso = pfmuon->ecalIso();
       MyPfmuon.hcalIso = pfmuon->hcalIso();
       MyPfmuon.caloIso = pfmuon->caloIso();

       MyPfmuon.trackIsoDR03 = pfmuon->isolationR03().sumPt;
       MyPfmuon.ecalIsoDR03 = pfmuon->isolationR03().emEt;
       MyPfmuon.hcalIsoDR03 = pfmuon->isolationR03().hadEt;
       MyPfmuon.caloIsoDR03 = pfmuon->isolationR03().emEt + pfmuon->isolationR03().hadEt;

       MyPfmuon.trackVetoIsoDR03 = pfmuon->isolationR03().trackerVetoPt;
       MyPfmuon.ecalVetoIsoDR03 = pfmuon->isolationR03().emVetoEt;
       MyPfmuon.hcalVetoIsoDR03 = pfmuon->isolationR03().hadVetoEt;
       MyPfmuon.caloVetoIsoDR03 = pfmuon->isolationR03().emVetoEt + pfmuon->isolationR03().hadVetoEt;

       MyPfmuon.trackIsoDR05 = pfmuon->isolationR05().sumPt;
       MyPfmuon.ecalIsoDR05 = pfmuon->isolationR05().emEt;
       MyPfmuon.hcalIsoDR05 = pfmuon->isolationR05().hadEt;
       MyPfmuon.caloIsoDR05 = pfmuon->isolationR05().emEt + pfmuon->isolationR05().hadEt;

       MyPfmuon.trackVetoIsoDR05 = pfmuon->isolationR05().trackerVetoPt;
       MyPfmuon.ecalVetoIsoDR05 = pfmuon->isolationR05().emVetoEt;
       MyPfmuon.hcalVetoIsoDR05 = pfmuon->isolationR05().hadVetoEt;
       MyPfmuon.caloVetoIsoDR05 = pfmuon->isolationR05().emVetoEt + pfmuon->isolationR05().hadVetoEt;

       ////   Removed due to not understood seg faults. Seems ok in bare ROOT.
       ////   Investigate and uncomment in future.
       //MyMuon.hcalE = muon->hcalIsoDeposit()->candEnergy();
       //MyMuon.ecalE = muon->ecalIsoDeposit()->candEnergy();

       MyPfmuon.timeAtIpInOut = pfmuon->time().timeAtIpInOut;
       MyPfmuon.timeAtIpInOutErr = pfmuon->time().timeAtIpInOutErr;
       MyPfmuon.timeAtIpOutIn = pfmuon->time().timeAtIpOutIn;
       MyPfmuon.timeAtIpOutInErr = pfmuon->time().timeAtIpOutInErr;
       MyPfmuon.time_ndof = pfmuon->time().nDof;

       if( pfmuon->isEnergyValid() ){
	 MyPfmuon.ecal_time = pfmuon->calEnergy().ecal_time;
	 MyPfmuon.hcal_time = pfmuon->calEnergy().hcal_time;
	 MyPfmuon.ecal_timeError = pfmuon->calEnergy().ecal_timeError;
	 MyPfmuon.hcal_timeError = pfmuon->calEnergy().hcal_timeError;
	 MyPfmuon.energy_ecal = pfmuon->calEnergy().em;
	 MyPfmuon.energy_hcal = pfmuon->calEnergy().had;
	 MyPfmuon.e3x3_ecal = pfmuon->calEnergy().emS9;
	 MyPfmuon.e3x3_hcal = pfmuon->calEnergy().hadS9;
	 MyPfmuon.energyMax_ecal = pfmuon->calEnergy().emMax;
	 MyPfmuon.energyMax_hcal = pfmuon->calEnergy().hadMax;
       }


       MyPfmuon.IDGMPTight = ( pfmuon->isGood("GlobalMuonPromptTight") ) ? 1 : 0;

       MyPfmuon.isPFMuon = ( pfmuon->isPFMuon() ) ? 1 : 0;
       MyPfmuon.isGoodMuon_1StationTight = ( pfmuon->isGood("TMOneStationTight") )  ? 1 : 0;
       MyPfmuon.isGlobalMuon = ( pfmuon->isGlobalMuon() ) ? 1 : 0;
       MyPfmuon.isTrackerMuon = ( pfmuon->isTrackerMuon() ) ? 1 : 0;
       MyPfmuon.isStandAloneMuon = ( pfmuon->isStandAloneMuon() ) ? 1 : 0;
       MyPfmuon.isGlobalMuonPromptTight = ( pfmuon->isGood("GlobalMuonPromptTight") ) ? 1 : 0;

       MyPfmuon.numberOfMatches = pfmuon->numberOfMatches();
       MyPfmuon.numberOfMatchedStations = pfmuon->numberOfMatchedStations();

       MyPfmuon.dVzPVz = pfmuon->vz() - PVz;
       MyPfmuon.dB = pfmuon->dB();

       MyPfmuon.IP = pfmuon->dB(pat::Muon::PV3D);
       MyPfmuon.IPError = pfmuon->edB(pat::Muon::PV3D);

       if( (pfmuon->globalTrack().isAvailable()) ){
	 MyPfmuon.normalizedChi2 = pfmuon->globalTrack()->normalizedChi2();
	 MyPfmuon.numberOfValidMuonHits = pfmuon->globalTrack()->hitPattern().numberOfValidMuonHits();
	 MyPfmuon.numberOfValidTrackerHits = pfmuon->globalTrack()->hitPattern().numberOfValidTrackerHits();
         MyPfmuon.numberOfLayersWithMeasurement = pfmuon->track()->hitPattern().trackerLayersWithMeasurement();
	 MyPfmuon.ptErr = pfmuon->globalTrack()->ptError();
       }
       if( (pfmuon->innerTrack().isAvailable()) ){
         MyPfmuon.innerTrackNormChi2 = pfmuon->innerTrack()->normalizedChi2();
	 MyPfmuon.numberOfValidTrackerHitsInnerTrack = pfmuon->innerTrack()->numberOfValidHits();
	 MyPfmuon.pixelLayersWithMeasurement = pfmuon->innerTrack()->hitPattern().pixelLayersWithMeasurement();
         MyPfmuon.numberOfValidPixelHits = pfmuon->innerTrack()->hitPattern().numberOfValidPixelHits();

	 MyPfmuon.correctedD0 = pfmuon->innerTrack()->dxy(beamSpotPosition);
         MyPfmuon.correctedD0Vertex = pfmuon->innerTrack()->dxy(vertexPosition);
	 MyPfmuon.correctedDZ = pfmuon->innerTrack()->dz(vertexPosition);
       }

       // Get track pfmuon info
       if( (pfmuon->track().isAvailable()) ){
	 double tkvx = pfmuon->track()->vx();
	 double tkvy = pfmuon->track()->vy();
	 double tkpx = pfmuon->track()->px();
	 double tkpy = pfmuon->track()->py();
	 double tkpt = pfmuon->track()->pt();

	 double ndof = pfmuon->track()->ndof();
	 if( (ndof!=0) ) MyPfmuon.tkNormChi2 = pfmuon->track()->chi2()/ndof;
	 MyPfmuon.tkPT = pfmuon->track()->pt();
	 MyPfmuon.tkEta = pfmuon->track()->eta();
	 MyPfmuon.tkPhi = pfmuon->track()->phi();
	 MyPfmuon.tkDZ = pfmuon->track()->dz();
	 MyPfmuon.tkD0 = pfmuon->track()->d0();
	 MyPfmuon.tkD0bs = (tkvx-BSx)*tkpy/tkpt-(tkvy-BSy)*tkpx/tkpt;
	 MyPfmuon.tkD0err = pfmuon->track()->d0Error();
	 MyPfmuon.tkNumValidHits = pfmuon->track()->numberOfValidHits();
	 MyPfmuon.tkCharge = pfmuon->track()->charge();
       }
       // Get standalone pfmuon info
       if( (pfmuon->standAloneMuon().isAvailable()) ){
	 double samvx = pfmuon->standAloneMuon()->vx();
	 double samvy = pfmuon->standAloneMuon()->vy();
	 double sampx = pfmuon->standAloneMuon()->px();
	 double sampy = pfmuon->standAloneMuon()->py();
	 double sampt = pfmuon->standAloneMuon()->pt();

	 double ndof = pfmuon->standAloneMuon()->ndof();
	 if( (ndof!=0) ) MyPfmuon.samNormChi2 = pfmuon->standAloneMuon()->chi2()/ndof;
	 MyPfmuon.samPT = pfmuon->standAloneMuon()->pt();
	 MyPfmuon.samEta = pfmuon->standAloneMuon()->eta();
	 MyPfmuon.samPhi = pfmuon->standAloneMuon()->phi();
	 MyPfmuon.samDZ = pfmuon->standAloneMuon()->dz();
	 MyPfmuon.samD0 = pfmuon->standAloneMuon()->d0();
	 MyPfmuon.samD0bs = (samvx-BSx)*sampy/sampt-(samvy-BSy)*sampx/sampt;
	 MyPfmuon.samD0err = pfmuon->standAloneMuon()->d0Error();
	 MyPfmuon.samNumValidHits = pfmuon->standAloneMuon()->numberOfValidHits();
	 MyPfmuon.samCharge = pfmuon->standAloneMuon()->charge();
       }
       // Get global pfmuon info
       if( (pfmuon->combinedMuon().isAvailable()) ){
	 double comvx = pfmuon->combinedMuon()->vx();
	 double comvy = pfmuon->combinedMuon()->vy();
	 double compx = pfmuon->combinedMuon()->px();
	 double compy = pfmuon->combinedMuon()->py();
	 double compt = pfmuon->combinedMuon()->pt();

	 double ndof = pfmuon->combinedMuon()->ndof();
	 if( (ndof!=0) ) MyPfmuon.comNormChi2 = pfmuon->combinedMuon()->chi2()/ndof;
	 MyPfmuon.comPT = pfmuon->combinedMuon()->pt();
	 MyPfmuon.comEta = pfmuon->combinedMuon()->eta();
	 MyPfmuon.comPhi = pfmuon->combinedMuon()->phi();
	 MyPfmuon.comDZ = pfmuon->combinedMuon()->dz();
	 MyPfmuon.comD0 = pfmuon->combinedMuon()->d0();
	 MyPfmuon.comD0bs = (comvx-BSx)*compy/compt-(comvy-BSy)*compx/compt;
	 MyPfmuon.comD0err = pfmuon->combinedMuon()->d0Error();
	 MyPfmuon.comNumValidHits = pfmuon->combinedMuon()->numberOfValidHits();
	 MyPfmuon.comCharge = pfmuon->combinedMuon()->charge();
       }

       if( (pfmuon->genLepton()) ){
	 int genId = pfmuon->genLepton()->pdgId();

	 MyPfmuon.genId = pfmuon->genLepton()->pdgId();
	 MyPfmuon.genET = pfmuon->genLepton()->et();
	 MyPfmuon.genPT = pfmuon->genLepton()->pt();
	 MyPfmuon.genPhi = pfmuon->genLepton()->phi();
	 MyPfmuon.genEta = pfmuon->genLepton()->eta();
	 MyPfmuon.genCharge = pfmuon->genLepton()->charge();

	 MyPfmuon.genNumberOfMothers = pfmuon->genLepton()->numberOfMothers();

	 if( (pfmuon->genLepton()->numberOfMothers()>0) ){
	   const reco::Candidate *PfmuonMother = pfmuon->genLepton()->mother();
	   bool staytrapped = true;
	   while( (PfmuonMother->pdgId()==genId && staytrapped) ){
	     if( PfmuonMother->numberOfMothers()>=1 ) PfmuonMother = PfmuonMother->mother();
	     else staytrapped = false;
	   }
       
	   MyPfmuon.genMotherId = PfmuonMother->pdgId();
	   MyPfmuon.genMotherET = PfmuonMother->et();
	   MyPfmuon.genMotherPT = PfmuonMother->pt();
	   MyPfmuon.genMotherPhi = PfmuonMother->phi();
	   MyPfmuon.genMotherEta = PfmuonMother->eta();
	   MyPfmuon.genMotherCharge = PfmuonMother->charge();
	 }

	 if( (pfmuon->genLepton()->numberOfMothers()>=1) ){
	   const reco::Candidate *PfmuonMother0 = pfmuon->genLepton()->mother(0);
	   const reco::Candidate *PfmuonMother1 = 0;

	   if( (pfmuon->genLepton()->numberOfMothers()>=2) ) PfmuonMother1 = pfmuon->genLepton()->mother(1);

	   int mother0id = PfmuonMother0->pdgId();
	   int mother1id = ( PfmuonMother1!=0 ) ? PfmuonMother1->pdgId() : -99;

	   bool staytrapped = true;
	   while( (mother0id==genId || mother1id==genId) && staytrapped ){
	     if( mother0id==genId && (PfmuonMother0!=0) ){
	       if( PfmuonMother0->numberOfMothers()>=1 ){
		 PfmuonMother0 = PfmuonMother0->mother(0);
		 mother0id = PfmuonMother0->pdgId();
		 mother1id = -99;
		 if( PfmuonMother0->numberOfMothers()>=2 ){
		   PfmuonMother1 = PfmuonMother0->mother(1);
		   mother1id = PfmuonMother1->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else if( mother1id==genId && (PfmuonMother1!=0) ){
	       if( PfmuonMother1->numberOfMothers()>=1 ){
		 PfmuonMother1 = PfmuonMother1->mother(0);
		 mother1id = PfmuonMother1->pdgId();
		 mother0id = -99;
		 if( PfmuonMother1->numberOfMothers()>=2 ){
		   PfmuonMother0 = PfmuonMother1->mother(1);
		   mother0id = PfmuonMother0->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else staytrapped = false;
	   }

	   if( mother0id!=-99 ){
	     MyPfmuon.genMother0Id = PfmuonMother0->pdgId();

	     if( (PfmuonMother0->numberOfMothers()>=1) ){
	       const reco::Candidate *PfmuonGrandMother0 = PfmuonMother0->mother(0);
	       const reco::Candidate *PfmuonGrandMother1 = 0;

	       if( (PfmuonMother0->numberOfMothers()>=2) ) PfmuonGrandMother1 = PfmuonMother0->mother(1);

	       int gmother0id = PfmuonGrandMother0->pdgId();
	       int gmother1id = ( PfmuonGrandMother1!=0 ) ? PfmuonGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother0id || gmother1id==mother0id) && staytrapped ){
		 if( gmother0id==mother0id && (PfmuonGrandMother0!=0) ){
		   if( PfmuonGrandMother0->numberOfMothers()>=1 ){
		     PfmuonGrandMother0 = PfmuonGrandMother0->mother(0);
		     gmother0id = PfmuonGrandMother0->pdgId();
		     gmother1id = -99;
		     if( PfmuonGrandMother0->numberOfMothers()>=2 ){
		       PfmuonGrandMother1 = PfmuonGrandMother0->mother(1);
		       gmother1id = PfmuonGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother0id && (PfmuonGrandMother1!=0) ){
		   if( PfmuonGrandMother1->numberOfMothers()>=1 ){
		     PfmuonGrandMother1 = PfmuonGrandMother1->mother(0);
		     gmother1id = PfmuonGrandMother1->pdgId();
		     gmother0id = -99;
		     if( PfmuonGrandMother1->numberOfMothers()>=2 ){
		       PfmuonGrandMother0 = PfmuonGrandMother1->mother(1);
		       gmother0id = PfmuonGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyPfmuon.genGrandMother00Id = gmother0id;
	       MyPfmuon.genGrandMother01Id = gmother1id;
	     }
	   }
	   if( mother1id!=-99 ){
	     MyPfmuon.genMother1Id = PfmuonMother1->pdgId();

	     if( (PfmuonMother1->numberOfMothers()>=1) ){
	       const reco::Candidate *PfmuonGrandMother0 = PfmuonMother1->mother(0);
	       const reco::Candidate *PfmuonGrandMother1 = 0;

	       if( (PfmuonMother0->numberOfMothers()>=2) ) PfmuonGrandMother1 = PfmuonMother1->mother(1);

	       int gmother0id = PfmuonGrandMother0->pdgId();
	       int gmother1id = ( PfmuonGrandMother1!=0 ) ? PfmuonGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother1id || gmother1id==mother1id) && staytrapped ){
		 if( gmother0id==mother1id && (PfmuonGrandMother0!=0) ){
		   if( PfmuonGrandMother0->numberOfMothers()>=1 ){
		     PfmuonGrandMother0 = PfmuonGrandMother0->mother(0);
		     gmother0id = PfmuonGrandMother0->pdgId();
		     gmother1id = -99;
		     if( PfmuonGrandMother0->numberOfMothers()>=2 ){
		       PfmuonGrandMother1 = PfmuonGrandMother0->mother(1);
		       gmother1id = PfmuonGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother1id && (PfmuonGrandMother1!=0) ){
		   if( PfmuonGrandMother1->numberOfMothers()>=1 ){
		     PfmuonGrandMother1 = PfmuonGrandMother1->mother(0);
		     gmother1id = PfmuonGrandMother1->pdgId();
		     gmother0id = -99;
		     if( PfmuonGrandMother1->numberOfMothers()>=2 ){
		       PfmuonGrandMother0 = PfmuonGrandMother1->mother(1);
		       gmother0id = PfmuonGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyPfmuon.genGrandMother10Id = gmother0id;
	       MyPfmuon.genGrandMother11Id = gmother1id;
	     }
	   }
	 }

       }// end check if genLepton

       bnpfmuons->push_back(MyPfmuon);
     }
   }




   /////////////////////////////////////////////
   ///////
   ///////   Fill the pfmuon collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNmuonCollection> bnpfmuonsLoose(new BNmuonCollection);
   if( producePFMuon ){
     edm::View<pat::Muon> pfmuonsLoose = *pfmuonLooseHandle;

     for( edm::View<pat::Muon>::const_iterator pfmuon = pfmuonsLoose.begin(); pfmuon!=pfmuonsLoose.end(); ++pfmuon ){

       BNmuon MyPfmuonLoose;

       // general kinematic variables
       MyPfmuonLoose.energy = pfmuon->energy();
       MyPfmuonLoose.et = pfmuon->et();
       MyPfmuonLoose.pt = pfmuon->pt();
       MyPfmuonLoose.px = pfmuon->px();
       MyPfmuonLoose.py = pfmuon->py();
       MyPfmuonLoose.pz = pfmuon->pz();
       MyPfmuonLoose.phi = pfmuon->phi();
       MyPfmuonLoose.eta = pfmuon->eta();
       MyPfmuonLoose.theta = pfmuon->theta();


       MyPfmuonLoose.charge = pfmuon->charge();
       MyPfmuonLoose.vx = pfmuon->vx();
       MyPfmuonLoose.vy = pfmuon->vy();
       MyPfmuonLoose.vz = pfmuon->vz();


       MyPfmuonLoose.particleIso = pfmuon->particleIso();     
       MyPfmuonLoose.chargedHadronIso = pfmuon->chargedHadronIso();     
       MyPfmuonLoose.neutralHadronIso = pfmuon->neutralHadronIso();     
       MyPfmuonLoose.photonIso = pfmuon->photonIso();     
       MyPfmuonLoose.puChargedHadronIso = pfmuon->puChargedHadronIso();     

       MyPfmuonLoose.chargedHadronIsoDR03 = pfmuon->userIso(0);     
       MyPfmuonLoose.neutralHadronIsoDR03 = pfmuon->userIso(1);
       MyPfmuonLoose.photonIsoDR03 = pfmuon->userIso(2);
       MyPfmuonLoose.puChargedHadronIsoDR03 = pfmuon->userIso(3);

       MyPfmuonLoose.chargedHadronIsoDR04 = pfmuon->userIso(4);     
       MyPfmuonLoose.neutralHadronIsoDR04 = pfmuon->userIso(5);
       MyPfmuonLoose.photonIsoDR04 = pfmuon->userIso(6);
       MyPfmuonLoose.puChargedHadronIsoDR04 = pfmuon->userIso(7);

       MyPfmuonLoose.rhoPrime = std::max(rho_event, 0.0);
       if( sample_<0 ){
	 MyPfmuonLoose.AEffDr03 = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaAndNeutralHadronIso03, pfmuon->eta(), MuonEffectiveArea::kMuEAData2012);
	 MyPfmuonLoose.AEffDr04 = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaAndNeutralHadronIso04, pfmuon->eta(), MuonEffectiveArea::kMuEAData2012);
       }
       else {
	 MyPfmuonLoose.AEffDr03 = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaAndNeutralHadronIso03, pfmuon->eta(), MuonEffectiveArea::kMuEAFall11MC);
	 MyPfmuonLoose.AEffDr04 = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaAndNeutralHadronIso04, pfmuon->eta(), MuonEffectiveArea::kMuEAFall11MC);
       }

       MyPfmuonLoose.pfIsoR03SumChargedHadronPt = pfmuon->pfIsolationR03().sumChargedHadronPt;
       MyPfmuonLoose.pfIsoR03SumNeutralHadronEt = pfmuon->pfIsolationR03().sumNeutralHadronEt;
       MyPfmuonLoose.pfIsoR03SumPhotonEt = pfmuon->pfIsolationR03().sumPhotonEt;
       MyPfmuonLoose.pfIsoR03SumPUPt = pfmuon->pfIsolationR03().sumPUPt;

       MyPfmuonLoose.pfIsoR04SumChargedHadronPt = pfmuon->pfIsolationR04().sumChargedHadronPt;
       MyPfmuonLoose.pfIsoR04SumNeutralHadronEt = pfmuon->pfIsolationR04().sumNeutralHadronEt;
       MyPfmuonLoose.pfIsoR04SumPhotonEt = pfmuon->pfIsolationR04().sumPhotonEt;
       MyPfmuonLoose.pfIsoR04SumPUPt = pfmuon->pfIsolationR04().sumPUPt;



       MyPfmuonLoose.trackIso = pfmuon->trackIso();
       MyPfmuonLoose.ecalIso = pfmuon->ecalIso();
       MyPfmuonLoose.hcalIso = pfmuon->hcalIso();
       MyPfmuonLoose.caloIso = pfmuon->caloIso();

       MyPfmuonLoose.trackIsoDR03 = pfmuon->isolationR03().sumPt;
       MyPfmuonLoose.ecalIsoDR03 = pfmuon->isolationR03().emEt;
       MyPfmuonLoose.hcalIsoDR03 = pfmuon->isolationR03().hadEt;
       MyPfmuonLoose.caloIsoDR03 = pfmuon->isolationR03().emEt + pfmuon->isolationR03().hadEt;

       MyPfmuonLoose.trackVetoIsoDR03 = pfmuon->isolationR03().trackerVetoPt;
       MyPfmuonLoose.ecalVetoIsoDR03 = pfmuon->isolationR03().emVetoEt;
       MyPfmuonLoose.hcalVetoIsoDR03 = pfmuon->isolationR03().hadVetoEt;
       MyPfmuonLoose.caloVetoIsoDR03 = pfmuon->isolationR03().emVetoEt + pfmuon->isolationR03().hadVetoEt;

       MyPfmuonLoose.trackIsoDR05 = pfmuon->isolationR05().sumPt;
       MyPfmuonLoose.ecalIsoDR05 = pfmuon->isolationR05().emEt;
       MyPfmuonLoose.hcalIsoDR05 = pfmuon->isolationR05().hadEt;
       MyPfmuonLoose.caloIsoDR05 = pfmuon->isolationR05().emEt + pfmuon->isolationR05().hadEt;

       MyPfmuonLoose.trackVetoIsoDR05 = pfmuon->isolationR05().trackerVetoPt;
       MyPfmuonLoose.ecalVetoIsoDR05 = pfmuon->isolationR05().emVetoEt;
       MyPfmuonLoose.hcalVetoIsoDR05 = pfmuon->isolationR05().hadVetoEt;
       MyPfmuonLoose.caloVetoIsoDR05 = pfmuon->isolationR05().emVetoEt + pfmuon->isolationR05().hadVetoEt;

       ////   Removed due to not understood seg faults. Seems ok in bare ROOT.
       ////   Investigate and uncomment in future.
       //MyMuon.hcalE = muon->hcalIsoDeposit()->candEnergy();
       //MyMuon.ecalE = muon->ecalIsoDeposit()->candEnergy();

       MyPfmuonLoose.timeAtIpInOut = pfmuon->time().timeAtIpInOut;
       MyPfmuonLoose.timeAtIpInOutErr = pfmuon->time().timeAtIpInOutErr;
       MyPfmuonLoose.timeAtIpOutIn = pfmuon->time().timeAtIpOutIn;
       MyPfmuonLoose.timeAtIpOutInErr = pfmuon->time().timeAtIpOutInErr;
       MyPfmuonLoose.time_ndof = pfmuon->time().nDof;

       if( pfmuon->isEnergyValid() ){
	 MyPfmuonLoose.ecal_time = pfmuon->calEnergy().ecal_time;
	 MyPfmuonLoose.hcal_time = pfmuon->calEnergy().hcal_time;
	 MyPfmuonLoose.ecal_timeError = pfmuon->calEnergy().ecal_timeError;
	 MyPfmuonLoose.hcal_timeError = pfmuon->calEnergy().hcal_timeError;
	 MyPfmuonLoose.energy_ecal = pfmuon->calEnergy().em;
	 MyPfmuonLoose.energy_hcal = pfmuon->calEnergy().had;
	 MyPfmuonLoose.e3x3_ecal = pfmuon->calEnergy().emS9;
	 MyPfmuonLoose.e3x3_hcal = pfmuon->calEnergy().hadS9;
	 MyPfmuonLoose.energyMax_ecal = pfmuon->calEnergy().emMax;
	 MyPfmuonLoose.energyMax_hcal = pfmuon->calEnergy().hadMax;
       }


       MyPfmuonLoose.IDGMPTight = ( pfmuon->isGood("GlobalMuonPromptTight") ) ? 1 : 0;

       MyPfmuonLoose.isPFMuon = ( pfmuon->isPFMuon() ) ? 1 : 0;
       MyPfmuonLoose.isGoodMuon_1StationTight = ( pfmuon->isGood("TMOneStationTight") )  ? 1 : 0;
       MyPfmuonLoose.isGlobalMuon = ( pfmuon->isGlobalMuon() ) ? 1 : 0;
       MyPfmuonLoose.isTrackerMuon = ( pfmuon->isTrackerMuon() ) ? 1 : 0;
       MyPfmuonLoose.isStandAloneMuon = ( pfmuon->isStandAloneMuon() ) ? 1 : 0;
       MyPfmuonLoose.isGlobalMuonPromptTight = ( pfmuon->isGood("GlobalMuonPromptTight") ) ? 1 : 0;

       MyPfmuonLoose.numberOfMatches = pfmuon->numberOfMatches();
       MyPfmuonLoose.numberOfMatchedStations = pfmuon->numberOfMatchedStations();

       MyPfmuonLoose.dVzPVz = pfmuon->vz() - PVz;
       MyPfmuonLoose.dB = pfmuon->dB();

       MyPfmuonLoose.IP = pfmuon->dB(pat::Muon::PV3D);
       MyPfmuonLoose.IPError = pfmuon->edB(pat::Muon::PV3D);

       if( (pfmuon->globalTrack().isAvailable()) ){
	 MyPfmuonLoose.normalizedChi2 = pfmuon->globalTrack()->normalizedChi2();
	 MyPfmuonLoose.numberOfValidMuonHits = pfmuon->globalTrack()->hitPattern().numberOfValidMuonHits();
	 MyPfmuonLoose.numberOfValidTrackerHits = pfmuon->globalTrack()->hitPattern().numberOfValidTrackerHits();
         MyPfmuonLoose.numberOfLayersWithMeasurement = pfmuon->track()->hitPattern().trackerLayersWithMeasurement();
	 MyPfmuonLoose.ptErr = pfmuon->globalTrack()->ptError();
       }
       if( (pfmuon->innerTrack().isAvailable()) ){
         MyPfmuonLoose.innerTrackNormChi2 = pfmuon->innerTrack()->normalizedChi2();
	 MyPfmuonLoose.numberOfValidTrackerHitsInnerTrack = pfmuon->innerTrack()->numberOfValidHits();
	 MyPfmuonLoose.pixelLayersWithMeasurement = pfmuon->innerTrack()->hitPattern().pixelLayersWithMeasurement();
         MyPfmuonLoose.numberOfValidPixelHits = pfmuon->innerTrack()->hitPattern().numberOfValidPixelHits();

	 MyPfmuonLoose.correctedD0 = pfmuon->innerTrack()->dxy(beamSpotPosition);
         MyPfmuonLoose.correctedD0Vertex = pfmuon->innerTrack()->dxy(vertexPosition);
	 MyPfmuonLoose.correctedDZ = pfmuon->innerTrack()->dz(vertexPosition);
       }

       // Get track pfmuon info
       if( (pfmuon->track().isAvailable()) ){
	 double tkvx = pfmuon->track()->vx();
	 double tkvy = pfmuon->track()->vy();
	 double tkpx = pfmuon->track()->px();
	 double tkpy = pfmuon->track()->py();
	 double tkpt = pfmuon->track()->pt();

	 double ndof = pfmuon->track()->ndof();
	 if( (ndof!=0) ) MyPfmuonLoose.tkNormChi2 = pfmuon->track()->chi2()/ndof;
	 MyPfmuonLoose.tkPT = pfmuon->track()->pt();
	 MyPfmuonLoose.tkEta = pfmuon->track()->eta();
	 MyPfmuonLoose.tkPhi = pfmuon->track()->phi();
	 MyPfmuonLoose.tkDZ = pfmuon->track()->dz();
	 MyPfmuonLoose.tkD0 = pfmuon->track()->d0();
	 MyPfmuonLoose.tkD0bs = (tkvx-BSx)*tkpy/tkpt-(tkvy-BSy)*tkpx/tkpt;
	 MyPfmuonLoose.tkD0err = pfmuon->track()->d0Error();
	 MyPfmuonLoose.tkNumValidHits = pfmuon->track()->numberOfValidHits();
	 MyPfmuonLoose.tkCharge = pfmuon->track()->charge();
       }
       // Get standalone pfmuon info
       if( (pfmuon->standAloneMuon().isAvailable()) ){
	 double samvx = pfmuon->standAloneMuon()->vx();
	 double samvy = pfmuon->standAloneMuon()->vy();
	 double sampx = pfmuon->standAloneMuon()->px();
	 double sampy = pfmuon->standAloneMuon()->py();
	 double sampt = pfmuon->standAloneMuon()->pt();

	 double ndof = pfmuon->standAloneMuon()->ndof();
	 if( (ndof!=0) ) MyPfmuonLoose.samNormChi2 = pfmuon->standAloneMuon()->chi2()/ndof;
	 MyPfmuonLoose.samPT = pfmuon->standAloneMuon()->pt();
	 MyPfmuonLoose.samEta = pfmuon->standAloneMuon()->eta();
	 MyPfmuonLoose.samPhi = pfmuon->standAloneMuon()->phi();
	 MyPfmuonLoose.samDZ = pfmuon->standAloneMuon()->dz();
	 MyPfmuonLoose.samD0 = pfmuon->standAloneMuon()->d0();
	 MyPfmuonLoose.samD0bs = (samvx-BSx)*sampy/sampt-(samvy-BSy)*sampx/sampt;
	 MyPfmuonLoose.samD0err = pfmuon->standAloneMuon()->d0Error();
	 MyPfmuonLoose.samNumValidHits = pfmuon->standAloneMuon()->numberOfValidHits();
	 MyPfmuonLoose.samCharge = pfmuon->standAloneMuon()->charge();
       }
       // Get global pfmuon info
       if( (pfmuon->combinedMuon().isAvailable()) ){
	 double comvx = pfmuon->combinedMuon()->vx();
	 double comvy = pfmuon->combinedMuon()->vy();
	 double compx = pfmuon->combinedMuon()->px();
	 double compy = pfmuon->combinedMuon()->py();
	 double compt = pfmuon->combinedMuon()->pt();

	 double ndof = pfmuon->combinedMuon()->ndof();
	 if( (ndof!=0) ) MyPfmuonLoose.comNormChi2 = pfmuon->combinedMuon()->chi2()/ndof;
	 MyPfmuonLoose.comPT = pfmuon->combinedMuon()->pt();
	 MyPfmuonLoose.comEta = pfmuon->combinedMuon()->eta();
	 MyPfmuonLoose.comPhi = pfmuon->combinedMuon()->phi();
	 MyPfmuonLoose.comDZ = pfmuon->combinedMuon()->dz();
	 MyPfmuonLoose.comD0 = pfmuon->combinedMuon()->d0();
	 MyPfmuonLoose.comD0bs = (comvx-BSx)*compy/compt-(comvy-BSy)*compx/compt;
	 MyPfmuonLoose.comD0err = pfmuon->combinedMuon()->d0Error();
	 MyPfmuonLoose.comNumValidHits = pfmuon->combinedMuon()->numberOfValidHits();
	 MyPfmuonLoose.comCharge = pfmuon->combinedMuon()->charge();
       }

       if( (pfmuon->genLepton()) ){
	 int genId = pfmuon->genLepton()->pdgId();

	 MyPfmuonLoose.genId = pfmuon->genLepton()->pdgId();
	 MyPfmuonLoose.genET = pfmuon->genLepton()->et();
	 MyPfmuonLoose.genPT = pfmuon->genLepton()->pt();
	 MyPfmuonLoose.genPhi = pfmuon->genLepton()->phi();
	 MyPfmuonLoose.genEta = pfmuon->genLepton()->eta();
	 MyPfmuonLoose.genCharge = pfmuon->genLepton()->charge();

	 MyPfmuonLoose.genNumberOfMothers = pfmuon->genLepton()->numberOfMothers();

	 if( (pfmuon->genLepton()->numberOfMothers()>0) ){
	   const reco::Candidate *PfmuonMother = pfmuon->genLepton()->mother();
	   bool staytrapped = true;
	   while( (PfmuonMother->pdgId()==genId && staytrapped) ){
	     if( PfmuonMother->numberOfMothers()>=1 ) PfmuonMother = PfmuonMother->mother();
	     else staytrapped = false;
	   }
       
	   MyPfmuonLoose.genMotherId = PfmuonMother->pdgId();
	   MyPfmuonLoose.genMotherET = PfmuonMother->et();
	   MyPfmuonLoose.genMotherPT = PfmuonMother->pt();
	   MyPfmuonLoose.genMotherPhi = PfmuonMother->phi();
	   MyPfmuonLoose.genMotherEta = PfmuonMother->eta();
	   MyPfmuonLoose.genMotherCharge = PfmuonMother->charge();
	 }

	 if( (pfmuon->genLepton()->numberOfMothers()>=1) ){
	   const reco::Candidate *PfmuonMother0 = pfmuon->genLepton()->mother(0);
	   const reco::Candidate *PfmuonMother1 = 0;

	   if( (pfmuon->genLepton()->numberOfMothers()>=2) ) PfmuonMother1 = pfmuon->genLepton()->mother(1);

	   int mother0id = PfmuonMother0->pdgId();
	   int mother1id = ( PfmuonMother1!=0 ) ? PfmuonMother1->pdgId() : -99;

	   bool staytrapped = true;
	   while( (mother0id==genId || mother1id==genId) && staytrapped ){
	     if( mother0id==genId && (PfmuonMother0!=0) ){
	       if( PfmuonMother0->numberOfMothers()>=1 ){
		 PfmuonMother0 = PfmuonMother0->mother(0);
		 mother0id = PfmuonMother0->pdgId();
		 mother1id = -99;
		 if( PfmuonMother0->numberOfMothers()>=2 ){
		   PfmuonMother1 = PfmuonMother0->mother(1);
		   mother1id = PfmuonMother1->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else if( mother1id==genId && (PfmuonMother1!=0) ){
	       if( PfmuonMother1->numberOfMothers()>=1 ){
		 PfmuonMother1 = PfmuonMother1->mother(0);
		 mother1id = PfmuonMother1->pdgId();
		 mother0id = -99;
		 if( PfmuonMother1->numberOfMothers()>=2 ){
		   PfmuonMother0 = PfmuonMother1->mother(1);
		   mother0id = PfmuonMother0->pdgId();
		 }
	       }
	       else staytrapped = false;
	     }
	     else staytrapped = false;
	   }

	   if( mother0id!=-99 ){
	     MyPfmuonLoose.genMother0Id = PfmuonMother0->pdgId();

	     if( (PfmuonMother0->numberOfMothers()>=1) ){
	       const reco::Candidate *PfmuonGrandMother0 = PfmuonMother0->mother(0);
	       const reco::Candidate *PfmuonGrandMother1 = 0;

	       if( (PfmuonMother0->numberOfMothers()>=2) ) PfmuonGrandMother1 = PfmuonMother0->mother(1);

	       int gmother0id = PfmuonGrandMother0->pdgId();
	       int gmother1id = ( PfmuonGrandMother1!=0 ) ? PfmuonGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother0id || gmother1id==mother0id) && staytrapped ){
		 if( gmother0id==mother0id && (PfmuonGrandMother0!=0) ){
		   if( PfmuonGrandMother0->numberOfMothers()>=1 ){
		     PfmuonGrandMother0 = PfmuonGrandMother0->mother(0);
		     gmother0id = PfmuonGrandMother0->pdgId();
		     gmother1id = -99;
		     if( PfmuonGrandMother0->numberOfMothers()>=2 ){
		       PfmuonGrandMother1 = PfmuonGrandMother0->mother(1);
		       gmother1id = PfmuonGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother0id && (PfmuonGrandMother1!=0) ){
		   if( PfmuonGrandMother1->numberOfMothers()>=1 ){
		     PfmuonGrandMother1 = PfmuonGrandMother1->mother(0);
		     gmother1id = PfmuonGrandMother1->pdgId();
		     gmother0id = -99;
		     if( PfmuonGrandMother1->numberOfMothers()>=2 ){
		       PfmuonGrandMother0 = PfmuonGrandMother1->mother(1);
		       gmother0id = PfmuonGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyPfmuonLoose.genGrandMother00Id = gmother0id;
	       MyPfmuonLoose.genGrandMother01Id = gmother1id;
	     }
	   }
	   if( mother1id!=-99 ){
	     MyPfmuonLoose.genMother1Id = PfmuonMother1->pdgId();

	     if( (PfmuonMother1->numberOfMothers()>=1) ){
	       const reco::Candidate *PfmuonGrandMother0 = PfmuonMother1->mother(0);
	       const reco::Candidate *PfmuonGrandMother1 = 0;

	       if( (PfmuonMother0->numberOfMothers()>=2) ) PfmuonGrandMother1 = PfmuonMother1->mother(1);

	       int gmother0id = PfmuonGrandMother0->pdgId();
	       int gmother1id = ( PfmuonGrandMother1!=0 ) ? PfmuonGrandMother1->pdgId() : -99;

	       bool staytrapped = true;
	       while( (gmother0id==mother1id || gmother1id==mother1id) && staytrapped ){
		 if( gmother0id==mother1id && (PfmuonGrandMother0!=0) ){
		   if( PfmuonGrandMother0->numberOfMothers()>=1 ){
		     PfmuonGrandMother0 = PfmuonGrandMother0->mother(0);
		     gmother0id = PfmuonGrandMother0->pdgId();
		     gmother1id = -99;
		     if( PfmuonGrandMother0->numberOfMothers()>=2 ){
		       PfmuonGrandMother1 = PfmuonGrandMother0->mother(1);
		       gmother1id = PfmuonGrandMother1->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else if( gmother1id==mother1id && (PfmuonGrandMother1!=0) ){
		   if( PfmuonGrandMother1->numberOfMothers()>=1 ){
		     PfmuonGrandMother1 = PfmuonGrandMother1->mother(0);
		     gmother1id = PfmuonGrandMother1->pdgId();
		     gmother0id = -99;
		     if( PfmuonGrandMother1->numberOfMothers()>=2 ){
		       PfmuonGrandMother0 = PfmuonGrandMother1->mother(1);
		       gmother0id = PfmuonGrandMother0->pdgId();
		     }
		   }
		   else staytrapped = false;
		 }
		 else staytrapped = false;
	       }

	       MyPfmuonLoose.genGrandMother10Id = gmother0id;
	       MyPfmuonLoose.genGrandMother11Id = gmother1id;
	     }
	   }
	 }

       }// end check if genLepton

       bnpfmuonsLoose->push_back(MyPfmuonLoose);
     }
   }


   /////////////////////////////////////////////
   ///////
   ///////   Fill the cocktail muon collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNmuonCollection> bncocktailmuons(new BNmuonCollection);
   if( produceCocktailMuon ){
     edm::View<reco::Muon> cocktailmuons = *cocktailmuonHandle;

     for( edm::View<reco::Muon>::const_iterator cocktailmuon = cocktailmuons.begin(); cocktailmuon!=cocktailmuons.end(); ++cocktailmuon ){

       BNmuon MyCocktailmuon;

       // general kinematic variables
       MyCocktailmuon.energy = cocktailmuon->energy();
       MyCocktailmuon.et = cocktailmuon->et();
       MyCocktailmuon.pt = cocktailmuon->pt();
       MyCocktailmuon.px = cocktailmuon->px();
       MyCocktailmuon.py = cocktailmuon->py();
       MyCocktailmuon.pz = cocktailmuon->pz();
       MyCocktailmuon.phi = cocktailmuon->phi();
       MyCocktailmuon.eta = cocktailmuon->eta();
       MyCocktailmuon.theta = cocktailmuon->theta();


       MyCocktailmuon.charge = cocktailmuon->charge();
       MyCocktailmuon.vx = cocktailmuon->vx();
       MyCocktailmuon.vy = cocktailmuon->vy();
       MyCocktailmuon.vz = cocktailmuon->vz();

       // MyCocktailmuon.trackIso = cocktailmuon->trackIso();
       // MyCocktailmuon.ecalIso = cocktailmuon->ecalIso();
       // MyCocktailmuon.hcalIso = cocktailmuon->hcalIso();
       // MyCocktailmuon.caloIso = cocktailmuon->caloIso();

       MyCocktailmuon.trackIsoDR03 = cocktailmuon->isolationR03().sumPt;
       MyCocktailmuon.ecalIsoDR03 = cocktailmuon->isolationR03().emEt;
       MyCocktailmuon.hcalIsoDR03 = cocktailmuon->isolationR03().hadEt;
       MyCocktailmuon.caloIsoDR03 = cocktailmuon->isolationR03().emEt + cocktailmuon->isolationR03().hadEt;

       MyCocktailmuon.trackVetoIsoDR03 = cocktailmuon->isolationR03().trackerVetoPt;
       MyCocktailmuon.ecalVetoIsoDR03 = cocktailmuon->isolationR03().emVetoEt;
       MyCocktailmuon.hcalVetoIsoDR03 = cocktailmuon->isolationR03().hadVetoEt;
       MyCocktailmuon.caloVetoIsoDR03 = cocktailmuon->isolationR03().emVetoEt + cocktailmuon->isolationR03().hadVetoEt;

       MyCocktailmuon.trackIsoDR05 = cocktailmuon->isolationR05().sumPt;
       MyCocktailmuon.ecalIsoDR05 = cocktailmuon->isolationR05().emEt;
       MyCocktailmuon.hcalIsoDR05 = cocktailmuon->isolationR05().hadEt;
       MyCocktailmuon.caloIsoDR05 = cocktailmuon->isolationR05().emEt + cocktailmuon->isolationR05().hadEt;

       MyCocktailmuon.trackVetoIsoDR05 = cocktailmuon->isolationR05().trackerVetoPt;
       MyCocktailmuon.ecalVetoIsoDR05 = cocktailmuon->isolationR05().emVetoEt;
       MyCocktailmuon.hcalVetoIsoDR05 = cocktailmuon->isolationR05().hadVetoEt;
       MyCocktailmuon.caloVetoIsoDR05 = cocktailmuon->isolationR05().emVetoEt + cocktailmuon->isolationR05().hadVetoEt;

       ////   Removed due to not understood seg faults. Seems ok in bare ROOT.
       ////   Investigate and uncomment in future.
       //MyCocktailmuon.hcalE = cocktailmuon->hcalIsoDeposit()->candEnergy();
       //MyCocktailmuon.ecalE = cocktailmuon->ecalIsoDeposit()->candEnergy();


       // MyCocktailmuon.IDGMPTight = ( cocktailmuon->isGood("GlobalMuonPromptTight") ) ? 1 : 0;

       MyCocktailmuon.isGlobalMuon = ( cocktailmuon->isGlobalMuon() ) ? 1 : 0;
       MyCocktailmuon.isTrackerMuon = ( cocktailmuon->isTrackerMuon() ) ? 1 : 0;
       MyCocktailmuon.isStandAloneMuon = ( cocktailmuon->isStandAloneMuon() ) ? 1 : 0;
       // MyCocktailmuon.isGlobalMuonPromptTight = ( cocktailmuon->isGood("GlobalMuonPromptTight") ) ? 1 : 0;

       MyCocktailmuon.numberOfMatches = cocktailmuon->numberOfMatches();
       MyCocktailmuon.numberOfMatchedStations = cocktailmuon->numberOfMatchedStations();

       MyCocktailmuon.dVzPVz = cocktailmuon->vz() - PVz;
       //MyCocktailmuon.dB = cocktailmuon->dB();

       if( (cocktailmuon->globalTrack().isAvailable()) ){
	 MyCocktailmuon.normalizedChi2 = cocktailmuon->globalTrack()->normalizedChi2();
	 MyCocktailmuon.numberOfValidMuonHits = cocktailmuon->globalTrack()->hitPattern().numberOfValidMuonHits();
	 MyCocktailmuon.numberOfValidTrackerHits = cocktailmuon->globalTrack()->hitPattern().numberOfValidTrackerHits();
         MyCocktailmuon.numberOfLayersWithMeasurement = cocktailmuon->track()->hitPattern().trackerLayersWithMeasurement();
	 MyCocktailmuon.ptErr = cocktailmuon->globalTrack()->ptError();
       }
       if( (cocktailmuon->innerTrack().isAvailable()) ){
	 MyCocktailmuon.numberOfValidTrackerHitsInnerTrack = cocktailmuon->innerTrack()->numberOfValidHits();
	 MyCocktailmuon.pixelLayersWithMeasurement = cocktailmuon->innerTrack()->hitPattern().pixelLayersWithMeasurement();
         MyCocktailmuon.numberOfValidPixelHits = cocktailmuon->innerTrack()->hitPattern().numberOfValidPixelHits();

	 MyCocktailmuon.correctedD0 = cocktailmuon->innerTrack()->dxy(beamSpotPosition);
         MyCocktailmuon.correctedD0Vertex = cocktailmuon->innerTrack()->dxy(vertexPosition);
	 MyCocktailmuon.correctedDZ = cocktailmuon->innerTrack()->dz(vertexPosition);
       }

       // Get track cocktailmuon info
       if( (cocktailmuon->track().isAvailable()) ){
	 double tkvx = cocktailmuon->track()->vx();
	 double tkvy = cocktailmuon->track()->vy();
	 double tkpx = cocktailmuon->track()->px();
	 double tkpy = cocktailmuon->track()->py();
	 double tkpt = cocktailmuon->track()->pt();

	 double ndof = cocktailmuon->track()->ndof();
	 if( (ndof!=0) ) MyCocktailmuon.tkNormChi2 = cocktailmuon->track()->chi2()/ndof;
	 MyCocktailmuon.tkPT = cocktailmuon->track()->pt();
	 MyCocktailmuon.tkEta = cocktailmuon->track()->eta();
	 MyCocktailmuon.tkPhi = cocktailmuon->track()->phi();
	 MyCocktailmuon.tkDZ = cocktailmuon->track()->dz();
	 MyCocktailmuon.tkD0 = cocktailmuon->track()->d0();
	 MyCocktailmuon.tkD0bs = (tkvx-BSx)*tkpy/tkpt-(tkvy-BSy)*tkpx/tkpt;
	 MyCocktailmuon.tkD0err = cocktailmuon->track()->d0Error();
	 MyCocktailmuon.tkNumValidHits = cocktailmuon->track()->numberOfValidHits();
	 MyCocktailmuon.tkCharge = cocktailmuon->track()->charge();
       }
       // Get standalone cocktailmuon info
       if( (cocktailmuon->standAloneMuon().isAvailable()) ){
	 double samvx = cocktailmuon->standAloneMuon()->vx();
	 double samvy = cocktailmuon->standAloneMuon()->vy();
	 double sampx = cocktailmuon->standAloneMuon()->px();
	 double sampy = cocktailmuon->standAloneMuon()->py();
	 double sampt = cocktailmuon->standAloneMuon()->pt();

	 double ndof = cocktailmuon->standAloneMuon()->ndof();
	 if( (ndof!=0) ) MyCocktailmuon.samNormChi2 = cocktailmuon->standAloneMuon()->chi2()/ndof;
	 MyCocktailmuon.samPT = cocktailmuon->standAloneMuon()->pt();
	 MyCocktailmuon.samEta = cocktailmuon->standAloneMuon()->eta();
	 MyCocktailmuon.samPhi = cocktailmuon->standAloneMuon()->phi();
	 MyCocktailmuon.samDZ = cocktailmuon->standAloneMuon()->dz();
	 MyCocktailmuon.samD0 = cocktailmuon->standAloneMuon()->d0();
	 MyCocktailmuon.samD0bs = (samvx-BSx)*sampy/sampt-(samvy-BSy)*sampx/sampt;
	 MyCocktailmuon.samD0err = cocktailmuon->standAloneMuon()->d0Error();
	 MyCocktailmuon.samNumValidHits = cocktailmuon->standAloneMuon()->numberOfValidHits();
	 MyCocktailmuon.samCharge = cocktailmuon->standAloneMuon()->charge();
       }
       // Get global cocktailmuon info
       if( (cocktailmuon->combinedMuon().isAvailable()) ){
	 double comvx = cocktailmuon->combinedMuon()->vx();
	 double comvy = cocktailmuon->combinedMuon()->vy();
	 double compx = cocktailmuon->combinedMuon()->px();
	 double compy = cocktailmuon->combinedMuon()->py();
	 double compt = cocktailmuon->combinedMuon()->pt();

	 double ndof = cocktailmuon->combinedMuon()->ndof();
	 if( (ndof!=0) ) MyCocktailmuon.comNormChi2 = cocktailmuon->combinedMuon()->chi2()/ndof;
	 MyCocktailmuon.comPT = cocktailmuon->combinedMuon()->pt();
	 MyCocktailmuon.comEta = cocktailmuon->combinedMuon()->eta();
	 MyCocktailmuon.comPhi = cocktailmuon->combinedMuon()->phi();
	 MyCocktailmuon.comDZ = cocktailmuon->combinedMuon()->dz();
	 MyCocktailmuon.comD0 = cocktailmuon->combinedMuon()->d0();
	 MyCocktailmuon.comD0bs = (comvx-BSx)*compy/compt-(comvy-BSy)*compx/compt;
	 MyCocktailmuon.comD0err = cocktailmuon->combinedMuon()->d0Error();
	 MyCocktailmuon.comNumValidHits = cocktailmuon->combinedMuon()->numberOfValidHits();
	 MyCocktailmuon.comCharge = cocktailmuon->combinedMuon()->charge();
       }

       bncocktailmuons->push_back(MyCocktailmuon);
     }
   }



   /////////////////////////////////////////////
   ///////
   ///////   Fill the tau collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNtauCollection> bntaus(new BNtauCollection);
   if( produceTau ){
     edm::View<pat::Tau> taus = *tauHandle;
     for( edm::View<pat::Tau>::const_iterator tau = taus.begin(); tau!=taus.end(); ++tau ){
        
		BNtau MyTau;

		MyTau.px											= tau->px();
		MyTau.py											= tau->py();
		MyTau.pz											= tau->pz();
		MyTau.pt											= tau->pt();
		MyTau.energy										= tau->energy();
		MyTau.et											= tau->et();
		MyTau.eta											= tau->eta();
		MyTau.phi											= tau->phi();
		MyTau.numProngs										= tau->signalPFChargedHadrCands().size();
		MyTau.numSignalGammas								= tau->signalPFGammaCands().size();
		MyTau.numSignalNeutrals								= tau->signalPFNeutrHadrCands().size();
		MyTau.numSignalPiZeros								= tau->signalPiZeroCandidates().size();
		MyTau.decayMode										= tau->decayMode();
		MyTau.emFraction									= tau->emFraction();
		MyTau.inTheCracks									= tauIsInTheCracks(tau->eta());
		MyTau.HPSagainstElectronLoose						= tau->tauID("againstElectronLoose");
		MyTau.HPSagainstElectronMVA							= tau->tauID("againstElectronMVA");
		MyTau.HPSagainstElectronMedium						= tau->tauID("againstElectronMedium");
		MyTau.HPSagainstElectronTight						= tau->tauID("againstElectronTight");
		MyTau.HPSagainstMuonLoose							= tau->tauID("againstMuonLoose");
		MyTau.HPSagainstMuonMedium							= tau->tauID("againstMuonMedium");
		MyTau.HPSagainstMuonTight							= tau->tauID("againstMuonTight");
		MyTau.HPSbyLooseCombinedIsolationDeltaBetaCorr		= tau->tauID("byLooseCombinedIsolationDeltaBetaCorr");
		MyTau.HPSbyMediumCombinedIsolationDeltaBetaCorr		= tau->tauID("byMediumCombinedIsolationDeltaBetaCorr");
		MyTau.HPSbyTightCombinedIsolationDeltaBetaCorr		= tau->tauID("byTightCombinedIsolationDeltaBetaCorr");
		MyTau.HPSbyVLooseCombinedIsolationDeltaBetaCorr		= tau->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
		MyTau.HPSdecayModeFinding							= tau->tauID("decayModeFinding");

		if(tau->leadPFChargedHadrCand().isNonnull()){
			MyTau.leadingTrackPt				= tau->leadPFChargedHadrCand()->pt();
			MyTau.charge						= tau->leadPFChargedHadrCand()->charge();

			if(tau->leadPFChargedHadrCand()->trackRef().isNonnull()){
				MyTau.leadingTrackValid				= 1;
				MyTau.leadingTrackIpVtdxy			= tau->leadPFChargedHadrCand()->trackRef()->dxy(vertexPosition);
				MyTau.leadingTrackIpVtdz			= tau->leadPFChargedHadrCand()->trackRef()->dz(vertexPosition);
				MyTau.leadingTrackIpVtdxyError		= tau->leadPFChargedHadrCand()->trackRef()->dxyError();
				MyTau.leadingTrackIpVtdzError		= tau->leadPFChargedHadrCand()->trackRef()->dzError();
				MyTau.leadingTrackVx				= tau->leadPFChargedHadrCand()->trackRef()->vx();
				MyTau.leadingTrackVy				= tau->leadPFChargedHadrCand()->trackRef()->vy();
				MyTau.leadingTrackVz				= tau->leadPFChargedHadrCand()->trackRef()->vz();
				MyTau.leadingTrackValidHits			= tau->leadPFChargedHadrCand()->trackRef()->numberOfValidHits();
				MyTau.leadingTrackNormChiSqrd		= tau->leadPFChargedHadrCand()->trackRef()->normalizedChi2();
			}else{
				MyTau.leadingTrackValid				= 0;
				MyTau.leadingTrackIpVtdxy			= -99;
				MyTau.leadingTrackIpVtdz			= -99;
				MyTau.leadingTrackIpVtdxyError		= -99;
				MyTau.leadingTrackIpVtdzError		= -99;
				MyTau.leadingTrackVx				= -99;
				MyTau.leadingTrackVy				= -99;
				MyTau.leadingTrackVz				= -99;
				MyTau.leadingTrackValidHits			= -99;
				MyTau.leadingTrackNormChiSqrd		= -99;
			}
		}else{
			MyTau.leadingTrackValid				= 0;
			MyTau.leadingTrackPt				= -99;
			MyTau.charge						= -99;
			MyTau.leadingTrackIpVtdxy			= -99;
			MyTau.leadingTrackIpVtdz			= -99;
			MyTau.leadingTrackIpVtdxyError		= -99;
			MyTau.leadingTrackIpVtdzError		= -99;
			MyTau.leadingTrackVx				= -99;
			MyTau.leadingTrackVy				= -99;
			MyTau.leadingTrackVz				= -99;
			MyTau.leadingTrackValidHits			= -99;
			MyTau.leadingTrackNormChiSqrd		= -99;
		}
        bntaus->push_back(MyTau);
     }
   }

   /////////////////////////////////////////////
   ///////
   ///////   Fill the photon collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNphotonCollection> bnphotons(new BNphotonCollection);
   if( producePhoton ){
     edm::View<pat::Photon> photons = *photonHandle;

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
   }


   int numhighpurity=0;
   reco::TrackBase::TrackQuality _trackQuality = reco::TrackBase::qualityByName("highPurity");

   /////////////////////////////////////////////
   ///////
   ///////   Fill the track collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNtrackCollection> bntracks(new BNtrackCollection);
   int tracksSize = 1;

   if( produceTrack ){
     reco::TrackCollection tracks = *trackHandle;

     for(reco::TrackCollection::const_iterator track = tracks.begin(); track!=tracks.end(); ++track ){

       bool isHighPurity = false;
       if( track->quality(_trackQuality) ){
	 numhighpurity++;
	 isHighPurity = true;
       }

       if( !(track->pt()>minTrackPt_) ) continue;

       BNtrack MyTrack;

       MyTrack.isHighPurity = ( isHighPurity ) ? 1 : 0;

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

       // Now get hit information. 
       // See reference for this class:  http://cmssdt.cern.ch/SDT/doxygen/CMSSW_4_2_5/doc/html/d9/d6d/classTrackingRecHit-members.html 
       if( fillTrackHitInfo_ ){
         for (trackingRecHit_iterator hit = track->recHitsBegin(); hit!=track->recHitsEnd(); hit++) {
           DetId detid = (*hit)->geographicalId();
           const GeomDetUnit *det = m_tracker->idToDetUnit(detid);  
           if (dynamic_cast<const StripGeomDetUnit*>(det)==0 && dynamic_cast<const PixelGeomDetUnit*>(det)==0) {
	     // std::cout << "this detID doesn't seem to belong to the Tracker" 
	     // 	       << "; subdetId = " << int(detid.subdetId())
	     // 	       << "; rawId = " << int(detid.rawId())
	     // 	       << endl;
             continue;
           }
           GlobalPoint center = det->surface().toGlobal(LocalPoint(0,0,0)); //should be the center of this Det. 

           MyTrack.subDetIdHits.push_back(int(detid.subdetId()));
           MyTrack.rawDetIdHits.push_back(int(detid.rawId()));
           MyTrack.isValidHits.push_back((*hit)->isValid());
           MyTrack.modulePerpHits.push_back(center.perp());
           MyTrack.moduleZHits.push_back(center.z());
           MyTrack.modulePhiHits.push_back(center.phi());           
         }  // end loop over hits
       }


       bntracks->push_back(MyTrack);
     }
     
     tracksSize = int(tracks.size());
   }

   bool FilterOutScraping = false;
   double FilterOutScrapingFraction = ( tracksSize!=0 ) ? (double)numhighpurity/(double)tracksSize : 0;

   if( tracksSize>numtrack_ ){
     if( FilterOutScrapingFraction>thresh_ ) FilterOutScraping=true;
   }
   else FilterOutScraping = true;



   /////////////////////////////////////////////
   ///////
   ///////   Fill the supercluster collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNsuperclusterCollection> bnsuperclusters(new BNsuperclusterCollection);
   if( produceSCsEB ){
     reco::SuperClusterCollection EBsuperclusters = *EBsuperclusterHandle;

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
   }
   if( produceSCsEE ){
     reco::SuperClusterCollection EEsuperclusters = *EEsuperclusterHandle;

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
   }



   /////////////////////////////////////////////
   ///////
   ///////   Fill the mcparticles collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNmcparticleCollection> bnmcparticles(new BNmcparticleCollection);
   std::auto_ptr<BNmcparticleCollection> bnmcelectrons(new BNmcparticleCollection);
   std::auto_ptr<BNmcparticleCollection> bnmcmuons(new BNmcparticleCollection);
   std::vector<int> vecZdecay;
   std::vector<int> vecWdecay;
   std::vector<int> vecHdecay;
   vecZdecay.clear();
   vecWdecay.clear();
   vecHdecay.clear();

   if( (genParticles.isValid() && sample_>=0 && produceGenParticle) ){
     for( size_t k = 0; k < genParticles->size(); k++ ){
       const reco::Candidate & mcParticle = (*genParticles)[k];

       int status = mcParticle.status();
       int pdgId  = mcParticle.pdgId();
       int numdgt = mcParticle.numberOfDaughters();

       BNmcparticle MyMCelectron;
       BNmcparticle MyMCmuon;

       if( (status==1) ){
	 if( abs(pdgId)==11 ){
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
	 else if( abs(pdgId)==13 ){
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

       int aId = abs(pdgId);
       bool keep = false;
       if( (status==3) ||
	   (aId==23 || aId==24 || aId==25 || aId==6 || aId==5 || aId==4 || aId==11 || aId==13 || aId==15) ) keep = true;

       if( !keep ) continue;


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
       if( (mcParticle.numberOfMothers()>=1) ){
	 const reco::Candidate *Mother0 = mcParticle.mother(0);
	 const reco::Candidate *Mother1 = 0;

	 if( (mcParticle.numberOfMothers()>=2) ) Mother1 = mcParticle.mother(1);

	 int mother0id = Mother0->pdgId();
	 int mother1id = ( Mother1!=0 ) ? Mother1->pdgId() : -99;

	 bool staytrapped = true;
	 while( (mother0id==pdgId || mother1id==pdgId) && staytrapped ){
	   if( mother0id==pdgId && (Mother0!=0) ){
	     if( Mother0->numberOfMothers()>=1 ){
	       Mother0 = Mother0->mother(0);
	       mother0id = Mother0->pdgId();
	       mother1id = -99;
	       if( Mother0->numberOfMothers()>=2 ){
		 Mother1 = Mother0->mother(1);
		 mother1id = Mother1->pdgId();
	       }
	     }
	     else staytrapped = false;
	   }
	   else if( mother1id==pdgId && (Mother1!=0) ){
	     if( Mother1->numberOfMothers()>=1 ){
	       Mother1 = Mother1->mother(0);
	       mother1id = Mother1->pdgId();
	       mother0id = -99;
	       if( Mother1->numberOfMothers()>=2 ){
		 Mother0 = Mother1->mother(1);
		 mother0id = Mother0->pdgId();
	       }
	     }
	     else staytrapped = false;
	   }
	   else staytrapped = false;
	 }

	 if( mother0id!=-99 ){
	   MyMCparticle.mother0Id = Mother0->pdgId();
	   MyMCparticle.mother0Status = Mother0->status();
	   MyMCparticle.mother0Charge = Mother0->charge();
	   MyMCparticle.mother0ET = Mother0->et();
	   MyMCparticle.mother0PT = Mother0->pt();
	   MyMCparticle.mother0Phi = Mother0->phi();
	   MyMCparticle.mother0Eta = Mother0->eta();

	   if( (Mother0->numberOfMothers()>=1) ){
	     const reco::Candidate *GrandMother0 = Mother0->mother(0);
	     const reco::Candidate *GrandMother1 = 0;

	     if( (Mother0->numberOfMothers()>=2) ) GrandMother1 = Mother0->mother(1);

	     int gmother0id = GrandMother0->pdgId();
	     int gmother1id = ( GrandMother1!=0 ) ? GrandMother1->pdgId() : -99;

	     bool staytrapped = true;
	     while( (gmother0id==mother0id || gmother1id==mother0id) && staytrapped ){
	       if( gmother0id==mother0id && (GrandMother0!=0) ){
		 if( GrandMother0->numberOfMothers()>=1 ){
		   GrandMother0 = GrandMother0->mother(0);
		   gmother0id = GrandMother0->pdgId();
		   gmother1id = -99;
		   if( GrandMother0->numberOfMothers()>=2 ){
		     GrandMother1 = GrandMother0->mother(1);
		     gmother1id = GrandMother1->pdgId();
		   }
		 }
		 else staytrapped = false;
	       }
	       else if( gmother1id==mother0id && (GrandMother1!=0) ){
		 if( GrandMother1->numberOfMothers()>=1 ){
		   GrandMother1 = GrandMother1->mother(0);
		   gmother1id = GrandMother1->pdgId();
		   gmother0id = -99;
		   if( GrandMother1->numberOfMothers()>=2 ){
		     GrandMother0 = GrandMother1->mother(1);
		     gmother0id = GrandMother0->pdgId();
		   }
		 }
		 else staytrapped = false;
	       }
	       else staytrapped = false;
	     }

	     if( gmother0id!=-99 ){
	       MyMCparticle.grandMother00Id = GrandMother0->pdgId();
	       MyMCparticle.grandMother00Status = GrandMother0->status();
	       MyMCparticle.grandMother00Charge = GrandMother0->charge();
	       MyMCparticle.grandMother00ET = GrandMother0->et();
	       MyMCparticle.grandMother00PT = GrandMother0->pt();
	       MyMCparticle.grandMother00Phi = GrandMother0->phi();
	       MyMCparticle.grandMother00Eta = GrandMother0->eta();
	     }
	     if( gmother1id!=-99 ){
	       MyMCparticle.grandMother01Id = GrandMother1->pdgId();
	       MyMCparticle.grandMother01Status = GrandMother1->status();
	       MyMCparticle.grandMother01Charge = GrandMother1->charge();
	       MyMCparticle.grandMother01ET = GrandMother1->et();
	       MyMCparticle.grandMother01PT = GrandMother1->pt();
	       MyMCparticle.grandMother01Phi = GrandMother1->phi();
	       MyMCparticle.grandMother01Eta = GrandMother1->eta();
	     }
	   }
	 }
	 if( mother1id!=-99 ){
	   MyMCparticle.mother1Id = Mother1->pdgId();
	   MyMCparticle.mother1Status = Mother1->status();
	   MyMCparticle.mother1Charge = Mother1->charge();
	   MyMCparticle.mother1ET = Mother1->et();
	   MyMCparticle.mother1PT = Mother1->pt();
	   MyMCparticle.mother1Phi = Mother1->phi();
	   MyMCparticle.mother1Eta = Mother1->eta();

	   if( (Mother1->numberOfMothers()>=1) ){
	     const reco::Candidate *GrandMother0 = Mother1->mother(0);
	     const reco::Candidate *GrandMother1 = 0;

	     if( (Mother0->numberOfMothers()>=2) ) GrandMother1 = Mother1->mother(1);

	     int gmother0id = GrandMother0->pdgId();
	     int gmother1id = ( GrandMother1!=0 ) ? GrandMother1->pdgId() : -99;

	     bool staytrapped = true;
	     while( (gmother0id==mother1id || gmother1id==mother1id) && staytrapped ){
	       if( gmother0id==mother1id && (GrandMother0!=0) ){
		 if( GrandMother0->numberOfMothers()>=1 ){
		   GrandMother0 = GrandMother0->mother(0);
		   gmother0id = GrandMother0->pdgId();
		   gmother1id = -99;
		   if( GrandMother0->numberOfMothers()>=2 ){
		     GrandMother1 = GrandMother0->mother(1);
		     gmother1id = GrandMother1->pdgId();
		   }
		 }
		 else staytrapped = false;
	       }
	       else if( gmother1id==mother1id && (GrandMother1!=0) ){
		 if( GrandMother1->numberOfMothers()>=1 ){
		   GrandMother1 = GrandMother1->mother(0);
		   gmother1id = GrandMother1->pdgId();
		   gmother0id = -99;
		   if( GrandMother1->numberOfMothers()>=2 ){
		     GrandMother0 = GrandMother1->mother(1);
		     gmother0id = GrandMother0->pdgId();
		   }
		 }
		 else staytrapped = false;
	       }
	       else staytrapped = false;
	     }

	     if( gmother0id!=-99 ){
	       MyMCparticle.grandMother10Id = GrandMother0->pdgId();
	       MyMCparticle.grandMother10Status = GrandMother0->status();
	       MyMCparticle.grandMother10Charge = GrandMother0->charge();
	       MyMCparticle.grandMother10ET = GrandMother0->et();
	       MyMCparticle.grandMother10PT = GrandMother0->pt();
	       MyMCparticle.grandMother10Phi = GrandMother0->phi();
	       MyMCparticle.grandMother10Eta = GrandMother0->eta();
	     }
	     if( gmother1id!=-99 ){
	       MyMCparticle.grandMother11Id = GrandMother1->pdgId();
	       MyMCparticle.grandMother11Status = GrandMother1->status();
	       MyMCparticle.grandMother11Charge = GrandMother1->charge();
	       MyMCparticle.grandMother11ET = GrandMother1->et();
	       MyMCparticle.grandMother11PT = GrandMother1->pt();
	       MyMCparticle.grandMother11Phi = GrandMother1->phi();
	       MyMCparticle.grandMother11Eta = GrandMother1->eta();
	     }
	   }
	 }
       }


       if( (mcParticle.numberOfDaughters()>=1) ){
	 MyMCparticle.daughter0Id = mcParticle.daughter(0)->pdgId();
	 MyMCparticle.daughter0Status = mcParticle.daughter(0)->status();
	 MyMCparticle.daughter0Charge = mcParticle.daughter(0)->charge();
	 MyMCparticle.daughter0ET = mcParticle.daughter(0)->et();
	 MyMCparticle.daughter0PT = mcParticle.daughter(0)->pt();
	 MyMCparticle.daughter0Phi = mcParticle.daughter(0)->phi();
	 MyMCparticle.daughter0Eta = mcParticle.daughter(0)->eta();

	 if( (mcParticle.numberOfDaughters()>=2) ){
	   MyMCparticle.daughter1Id = mcParticle.daughter(1)->pdgId();
	   MyMCparticle.daughter1Status = mcParticle.daughter(1)->status();
	   MyMCparticle.daughter1Charge = mcParticle.daughter(1)->charge();
	   MyMCparticle.daughter1ET = mcParticle.daughter(1)->et();
	   MyMCparticle.daughter1PT = mcParticle.daughter(1)->pt();
	   MyMCparticle.daughter1Phi = mcParticle.daughter(1)->phi();
	   MyMCparticle.daughter1Eta = mcParticle.daughter(1)->eta();
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
       else if( (abs(pdgId)==25 && numdgt>=2 && vecHdecay.size()<=2) ){
	 int dId0 = abs(mcParticle.daughter(0)->pdgId());
	 int dId1 = abs(mcParticle.daughter(1)->pdgId());

	 int hdecay0 = 100*dId0 + dId1;
	 int hdecay1 = 100*dId1 + dId0;

	 int hdecay = ( hdecay0 < hdecay1 ) ? hdecay0 : hdecay1;

	 vecHdecay.push_back(hdecay);
       }

     }
   }



   /////////////////////////////////////////////
   ///////
   ///////   Fill the mcparticles collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNgenjetCollection> bngenjets(new BNgenjetCollection);
   if( sample_>=0 && produceGenJet ){
     reco::GenJetCollection genjets = *genjetHandle;

     for(reco::GenJetCollection::const_iterator genjet = genjets.begin(); genjet!=genjets.end(); ++genjet ){

       if( !(genjet->pt()>10. && fabs(genjet->eta())<3.5) ) continue;

       BNgenjet MyGenjet;

       MyGenjet.pt = genjet->pt();
       MyGenjet.px = genjet->px();
       MyGenjet.py = genjet->py();
       MyGenjet.pz = genjet->pz();
       MyGenjet.phi = genjet->phi();
       MyGenjet.eta = genjet->eta();
       MyGenjet.et = genjet->et();
       MyGenjet.energy = genjet->energy();
       MyGenjet.mass = genjet->mass();
       MyGenjet.emEnergy = genjet->emEnergy();
       MyGenjet.hadEnergy = genjet->hadEnergy();
       MyGenjet.invisibleEnergy = genjet->invisibleEnergy();
       MyGenjet.auxiliaryEnergy = genjet->auxiliaryEnergy();
       MyGenjet.charge = genjet->charge();

       bngenjets->push_back(MyGenjet);
     }
   }


   /////////////////////////////////////////////
   ///////
   ///////   Fill the event collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNeventCollection> bnevent(new BNeventCollection);
   BNevent MyEvent;


   MyEvent.weight = weight;
   MyEvent.pthat = pthat;
   MyEvent.qScale = qScale;
   MyEvent.alphaQCD = alphaQCD;
   MyEvent.alphaQED = alphaQED;
   MyEvent.scalePDF = scalePDF;
   MyEvent.id1 = id1;
   MyEvent.id2 = id2;
   MyEvent.x1 = x1;
   MyEvent.x2 = x2;
   MyEvent.xPDF1 = xPDF1;
   MyEvent.xPDF2 = xPDF2;
   MyEvent.BSx = BSx;
   MyEvent.BSy = BSy;
   MyEvent.BSz = BSz;
   MyEvent.PVx = PVx;
   MyEvent.PVy = PVy;
   MyEvent.PVz = PVz;
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

   MyEvent.sumNVtx     = sum_nvtx_gen;
   MyEvent.sumTrueNVtx = sum_nvtx_true;
   MyEvent.numGenPV = npv_gen;
   MyEvent.nm1 = nm1;
   MyEvent.n0  = n0;
   MyEvent.np1 = np1;
   MyEvent.numTruePV = npv_true;
   MyEvent.nm1_true = nm1_true;
   MyEvent.n0_true  = n0_true;
   MyEvent.np1_true = np1_true;

   MyEvent.bField = evt_bField;

   MyEvent.Q2ScaleUpWgt = Q2ScaleUpWgt;
   MyEvent.Q2ScaleDownWgt = Q2ScaleDownWgt;

   MyEvent.rho_kt6PFJets = rho_event;
   MyEvent.rho_kt6PFJetsCentralChargedPileUp = rho_event_CentralChargedPileUp;
   MyEvent.rho_kt6PFJetsCentralNeutral       = rho_event_CentralNeutral;
   MyEvent.rho_kt6PFJetsCentralNeutralTight  = rho_event_CentralNeutralTight;


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
  if( (vecHdecay.size()>0) ){
     MyEvent.H0decay = vecHdecay[0];
     if( (vecHdecay.size()>1) ){
       MyEvent.H1decay = vecHdecay[1];
     }
   }


   bnevent->push_back(MyEvent);




   /////////////////////////////////////////////
   ///////
   ///////   Fill the calomet collection
   ///////
   /////////////////////////////////////////////
   /*
   std::auto_ptr<BNmetCollection> bncalomet(new BNmetCollection);
   BNmet MyCalomet;

   if( produceCaloMET ){
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

     double sigmaX2_calo = (calometHandle->front()).getSignificanceMatrix()(0,0);
     double sigmaY2_calo = (calometHandle->front()).getSignificanceMatrix()(1,1);
     double sigmaXY_calo = (calometHandle->front()).getSignificanceMatrix()(0,1);
     double sigmaYX_calo = (calometHandle->front()).getSignificanceMatrix()(1,0);

     double significance_calo = 99;
     if(sigmaX2_calo<1.e10 && sigmaY2_calo<1.e10) significance_calo = (calometHandle->front()).significance();

     MyCalomet.significance = significance_calo;
     MyCalomet.sigmaX2 = sigmaX2_calo;
     MyCalomet.sigmaY2 = sigmaY2_calo;
     MyCalomet.sigmaXY = sigmaXY_calo;
     MyCalomet.sigmaYX = sigmaYX_calo;

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
   }
   */


   /////////////////////////////////////////////
   ///////
   ///////   Fill the pfmet collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNmetCollection> bnpfmet(new BNmetCollection);
   BNmet MyPfmet;

   if( producePFMET ){
     MyPfmet.pt = pfmetHandle->front().pt();
     MyPfmet.px = pfmetHandle->front().px();
     MyPfmet.py = pfmetHandle->front().py();
     MyPfmet.phi = pfmetHandle->front().phi();
     MyPfmet.sumET = pfmetHandle->front().sumEt();
     MyPfmet.corSumET = pfmetHandle->front().corSumEt();
     MyPfmet.Upt = pfmetHandle->front().uncorrectedPt();
     MyPfmet.Uphi = pfmetHandle->front().uncorrectedPhi();
     if (pfmetHandle->front().isPFMET()) { //Appears to always return false
       MyPfmet.NeutralEMFraction = pfmetHandle->front().NeutralEMFraction();
       MyPfmet.NeutralHadEtFraction = pfmetHandle->front().NeutralHadEtFraction();
       MyPfmet.ChargedEMEtFraction = pfmetHandle->front().ChargedEMEtFraction();
       MyPfmet.ChargedHadEtFraction = pfmetHandle->front().ChargedHadEtFraction();
       MyPfmet.MuonEtFraction = pfmetHandle->front().MuonEtFraction();
       MyPfmet.Type6EtFraction = pfmetHandle->front().Type6EtFraction();
       MyPfmet.Type7EtFraction = pfmetHandle->front().Type7EtFraction();
     }

     double sigmaX2_pf = (pfmetHandle->front()).getSignificanceMatrix()(0,0);
     double sigmaY2_pf = (pfmetHandle->front()).getSignificanceMatrix()(1,1);
     double sigmaXY_pf = (pfmetHandle->front()).getSignificanceMatrix()(0,1);
     double sigmaYX_pf = (pfmetHandle->front()).getSignificanceMatrix()(1,0);

     double significance_pf = 99;
     if(sigmaX2_pf<1.e10 && sigmaY2_pf<1.e10) significance_pf = (pfmetHandle->front()).significance();

     MyPfmet.significance = significance_pf;
     MyPfmet.sigmaX2 = sigmaX2_pf;
     MyPfmet.sigmaY2 = sigmaY2_pf;
     MyPfmet.sigmaXY = sigmaXY_pf;
     MyPfmet.sigmaYX = sigmaYX_pf;

     if( (pfmetHandle->front().genMET()) ){
       MyPfmet.genPT = pfmetHandle->front().genMET()->pt();
       MyPfmet.genPhi = pfmetHandle->front().genMET()->phi();
     }

     bnpfmet->push_back(MyPfmet);
   }


   std::auto_ptr<BNmetCollection> bnpfmet_type1correctedRECO(new BNmetCollection);
   BNmet MyPfmet_type1correctedRECO;

   if( producePFMET_type1correctedRECO ){
     MyPfmet_type1correctedRECO.pt = pfmetHandle_type1correctedRECO->front().pt();
     MyPfmet_type1correctedRECO.px = pfmetHandle_type1correctedRECO->front().px();
     MyPfmet_type1correctedRECO.py = pfmetHandle_type1correctedRECO->front().py();
     MyPfmet_type1correctedRECO.phi = pfmetHandle_type1correctedRECO->front().phi();
     MyPfmet_type1correctedRECO.sumET = pfmetHandle_type1correctedRECO->front().sumEt();

     //pfmet specific quantities
     MyPfmet_type1correctedRECO.NeutralEMFraction = pfmetHandle_type1correctedRECO->front().photonEtFraction();
     MyPfmet_type1correctedRECO.NeutralHadEtFraction = pfmetHandle_type1correctedRECO->front().neutralHadronEtFraction();
     MyPfmet_type1correctedRECO.ChargedEMEtFraction = pfmetHandle_type1correctedRECO->front().electronEtFraction();
     MyPfmet_type1correctedRECO.ChargedHadEtFraction = pfmetHandle_type1correctedRECO->front().chargedHadronEtFraction();
     MyPfmet_type1correctedRECO.MuonEtFraction = pfmetHandle_type1correctedRECO->front().muonEtFraction();
     MyPfmet_type1correctedRECO.Type6EtFraction = pfmetHandle_type1correctedRECO->front().HFHadronEtFraction();
     MyPfmet_type1correctedRECO.Type7EtFraction = pfmetHandle_type1correctedRECO->front().HFEMEtFraction();

     double sigmaX2_pf = (pfmetHandle_type1correctedRECO->front()).getSignificanceMatrix()(0,0);
     double sigmaY2_pf = (pfmetHandle_type1correctedRECO->front()).getSignificanceMatrix()(1,1);
     double sigmaXY_pf = (pfmetHandle_type1correctedRECO->front()).getSignificanceMatrix()(0,1);
     double sigmaYX_pf = (pfmetHandle_type1correctedRECO->front()).getSignificanceMatrix()(1,0);

     double significance_pf = 99;
     if(sigmaX2_pf<1.e10 && sigmaY2_pf<1.e10) significance_pf = (pfmetHandle_type1correctedRECO->front()).significance();

     MyPfmet_type1correctedRECO.significance = significance_pf;
     MyPfmet_type1correctedRECO.sigmaX2 = sigmaX2_pf;
     MyPfmet_type1correctedRECO.sigmaY2 = sigmaY2_pf;
     MyPfmet_type1correctedRECO.sigmaXY = sigmaXY_pf;
     MyPfmet_type1correctedRECO.sigmaYX = sigmaYX_pf;

     bnpfmet_type1correctedRECO->push_back(MyPfmet_type1correctedRECO);
   }



   std::auto_ptr<BNmetCollection> bnpfmet_uncorrectedPF2PAT(new BNmetCollection);
   BNmet MyPfmet_uncorrectedPF2PAT;

   if( producePFMET_uncorrectedPF2PAT ){
     MyPfmet_uncorrectedPF2PAT.pt = pfmetHandle_uncorrectedPF2PAT->front().pt();
     MyPfmet_uncorrectedPF2PAT.px = pfmetHandle_uncorrectedPF2PAT->front().px();
     MyPfmet_uncorrectedPF2PAT.py = pfmetHandle_uncorrectedPF2PAT->front().py();
     MyPfmet_uncorrectedPF2PAT.phi = pfmetHandle_uncorrectedPF2PAT->front().phi();
     MyPfmet_uncorrectedPF2PAT.sumET = pfmetHandle_uncorrectedPF2PAT->front().sumEt();
     MyPfmet_uncorrectedPF2PAT.corSumET = pfmetHandle_uncorrectedPF2PAT->front().corSumEt();
     MyPfmet_uncorrectedPF2PAT.Upt = pfmetHandle_uncorrectedPF2PAT->front().uncorrectedPt();
     MyPfmet_uncorrectedPF2PAT.Uphi = pfmetHandle_uncorrectedPF2PAT->front().uncorrectedPhi();
     if (pfmetHandle_uncorrectedPF2PAT->front().isPFMET()) { //Appears to always return false
       MyPfmet_uncorrectedPF2PAT.NeutralEMFraction = pfmetHandle_uncorrectedPF2PAT->front().NeutralEMFraction();
       MyPfmet_uncorrectedPF2PAT.NeutralHadEtFraction = pfmetHandle_uncorrectedPF2PAT->front().NeutralHadEtFraction();
       MyPfmet_uncorrectedPF2PAT.ChargedEMEtFraction = pfmetHandle_uncorrectedPF2PAT->front().ChargedEMEtFraction();
       MyPfmet_uncorrectedPF2PAT.ChargedHadEtFraction = pfmetHandle_uncorrectedPF2PAT->front().ChargedHadEtFraction();
       MyPfmet_uncorrectedPF2PAT.MuonEtFraction = pfmetHandle_uncorrectedPF2PAT->front().MuonEtFraction();
       MyPfmet_uncorrectedPF2PAT.Type6EtFraction = pfmetHandle_uncorrectedPF2PAT->front().Type6EtFraction();
       MyPfmet_uncorrectedPF2PAT.Type7EtFraction = pfmetHandle_uncorrectedPF2PAT->front().Type7EtFraction();
     }
     double sigmaX2_pf = (pfmetHandle_uncorrectedPF2PAT->front()).getSignificanceMatrix()(0,0);
     double sigmaY2_pf = (pfmetHandle_uncorrectedPF2PAT->front()).getSignificanceMatrix()(1,1);
     double sigmaXY_pf = (pfmetHandle_uncorrectedPF2PAT->front()).getSignificanceMatrix()(0,1);
     double sigmaYX_pf = (pfmetHandle_uncorrectedPF2PAT->front()).getSignificanceMatrix()(1,0);

     double significance_pf = 99;
     if(sigmaX2_pf<1.e10 && sigmaY2_pf<1.e10) significance_pf = (pfmetHandle_uncorrectedPF2PAT->front()).significance();

     MyPfmet_uncorrectedPF2PAT.significance = significance_pf;
     MyPfmet_uncorrectedPF2PAT.sigmaX2 = sigmaX2_pf;
     MyPfmet_uncorrectedPF2PAT.sigmaY2 = sigmaY2_pf;
     MyPfmet_uncorrectedPF2PAT.sigmaXY = sigmaXY_pf;
     MyPfmet_uncorrectedPF2PAT.sigmaYX = sigmaYX_pf;

     if( (pfmetHandle_uncorrectedPF2PAT->front().genMET()) ){
       MyPfmet_uncorrectedPF2PAT.genPT = pfmetHandle_uncorrectedPF2PAT->front().genMET()->pt();
       MyPfmet_uncorrectedPF2PAT.genPhi = pfmetHandle_uncorrectedPF2PAT->front().genMET()->phi();
     }

     bnpfmet_uncorrectedPF2PAT->push_back(MyPfmet_uncorrectedPF2PAT);
   }



   std::auto_ptr<BNmetCollection> bnpfmet_uncorrectedRECO(new BNmetCollection);
   BNmet MyPfmet_uncorrectedRECO;

   if( producePFMET_uncorrectedRECO ){
     MyPfmet_uncorrectedRECO.pt = pfmetHandle_uncorrectedRECO->front().pt();
     MyPfmet_uncorrectedRECO.px = pfmetHandle_uncorrectedRECO->front().px();
     MyPfmet_uncorrectedRECO.py = pfmetHandle_uncorrectedRECO->front().py();
     MyPfmet_uncorrectedRECO.phi = pfmetHandle_uncorrectedRECO->front().phi();
     MyPfmet_uncorrectedRECO.sumET = pfmetHandle_uncorrectedRECO->front().sumEt();

     //pfmet specific quantities
     MyPfmet_uncorrectedRECO.NeutralEMFraction = pfmetHandle_uncorrectedRECO->front().photonEtFraction();
     MyPfmet_uncorrectedRECO.NeutralHadEtFraction = pfmetHandle_uncorrectedRECO->front().neutralHadronEtFraction();
     MyPfmet_uncorrectedRECO.ChargedEMEtFraction = pfmetHandle_uncorrectedRECO->front().electronEtFraction();
     MyPfmet_uncorrectedRECO.ChargedHadEtFraction = pfmetHandle_uncorrectedRECO->front().chargedHadronEtFraction();
     MyPfmet_uncorrectedRECO.MuonEtFraction = pfmetHandle_uncorrectedRECO->front().muonEtFraction();
     MyPfmet_uncorrectedRECO.Type6EtFraction = pfmetHandle_uncorrectedRECO->front().HFHadronEtFraction();
     MyPfmet_uncorrectedRECO.Type7EtFraction = pfmetHandle_uncorrectedRECO->front().HFEMEtFraction();

     double sigmaX2_pf = (pfmetHandle_uncorrectedRECO->front()).getSignificanceMatrix()(0,0);
     double sigmaY2_pf = (pfmetHandle_uncorrectedRECO->front()).getSignificanceMatrix()(1,1);
     double sigmaXY_pf = (pfmetHandle_uncorrectedRECO->front()).getSignificanceMatrix()(0,1);
     double sigmaYX_pf = (pfmetHandle_uncorrectedRECO->front()).getSignificanceMatrix()(1,0);

     double significance_pf = 99;
     if(sigmaX2_pf<1.e10 && sigmaY2_pf<1.e10) significance_pf = (pfmetHandle_uncorrectedRECO->front()).significance();

     MyPfmet_uncorrectedRECO.significance = significance_pf;
     MyPfmet_uncorrectedRECO.sigmaX2 = sigmaX2_pf;
     MyPfmet_uncorrectedRECO.sigmaY2 = sigmaY2_pf;
     MyPfmet_uncorrectedRECO.sigmaXY = sigmaXY_pf;
     MyPfmet_uncorrectedRECO.sigmaYX = sigmaYX_pf;

     bnpfmet_uncorrectedRECO->push_back(MyPfmet_uncorrectedRECO);
   }




   /////////////////////////////////////////////
   ///////
   ///////   Fill the tcmet collection
   ///////
   /////////////////////////////////////////////

   std::auto_ptr<BNmetCollection> bntcmet(new BNmetCollection);
   BNmet MyTcmet;

   if( produceTCMET ){
     MyTcmet.pt = tcmetHandle->front().pt();
     MyTcmet.px = tcmetHandle->front().px();
     MyTcmet.py = tcmetHandle->front().py();
     MyTcmet.phi = tcmetHandle->front().phi();
     MyTcmet.sumET = tcmetHandle->front().sumEt();
     MyTcmet.corSumET = tcmetHandle->front().corSumEt();
     MyTcmet.Upt = tcmetHandle->front().uncorrectedPt();
     MyTcmet.Uphi = tcmetHandle->front().uncorrectedPhi();

     double sigmaX2_tc = (tcmetHandle->front()).getSignificanceMatrix()(0,0);
     double sigmaY2_tc = (tcmetHandle->front()).getSignificanceMatrix()(1,1);
     double sigmaXY_tc = (tcmetHandle->front()).getSignificanceMatrix()(0,1);
     double sigmaYX_tc = (tcmetHandle->front()).getSignificanceMatrix()(1,0);

     double significance_tc = 99;
     if(sigmaX2_tc<1.e10 && sigmaY2_tc<1.e10) significance_tc = (tcmetHandle->front()).significance();

     MyTcmet.significance = significance_tc;
     MyTcmet.sigmaX2 = sigmaX2_tc;
     MyTcmet.sigmaY2 = sigmaY2_tc;
     MyTcmet.sigmaXY = sigmaXY_tc;
     MyTcmet.sigmaYX = sigmaYX_tc;

     if( (tcmetHandle->front().genMET()) ){
       MyTcmet.genPT = tcmetHandle->front().genMET()->pt();
       MyTcmet.genPhi = tcmetHandle->front().genMET()->phi();
     }


     bntcmet->push_back(MyTcmet);
   }

   // std::cout << "   ======> calomet: pt = " << calometHandle->front().pt() << ", phi = " << calometHandle->front().phi() << ", significance = " << significance_calo << ", nCorr = " << calometHandle->front().nCorrections() << ", uncorrected pt = " << calometHandle->front().uncorrectedPt() << std::endl;
   // std::cout << "   ======> pfmet:   pt = " << pfmetHandle->front().pt()   << ", phi = " << pfmetHandle->front().phi()   << ", significance = " << significance_pf << ", nCorr = " << pfmetHandle->front().nCorrections() << ", uncorrected pt = " << pfmetHandle->front().uncorrectedPt() << std::endl;
   // std::cout << "   ======> tcmet:   pt = " << tcmetHandle->front().pt()   << ", phi = " << tcmetHandle->front().phi()   << ", significance = " << significance_tc << ", nCorr = " << tcmetHandle->front().nCorrections() << ", uncorrected pt = " << tcmetHandle->front().uncorrectedPt() << std::endl;


   /////////////////////////////////////////////
   ///////
   ///////   Fill the trigger collection
   ///////
   /////////////////////////////////////////////

   /// HLTrigger/Configuration/python/HLT_8E29_cff.py
   //  hltL1NonIsoHLTNonIsoSingleElectronLWEt15PixelMatchFilter
   edm::Handle< trigger::TriggerEvent > hltHandle;
   iEvent.getByLabel(triggerSummaryTag_, hltHandle);

   std::auto_ptr<BNtrigobjCollection> bnhltobjs(new BNtrigobjCollection);

   std::vector<trigger::size_type> filtKeys;
   std::vector<std::string> filtKeys_string;
   for(unsigned int i=0; i<hltHandle->sizeFilters(); i++) {
     const edm::InputTag filterTag = hltHandle->filterTag(i);
     const std::string filt = filterTag.encode();

     // if( filt.find("Electron")!=std::string::npos ){
     //   std::cout << "  filt string = " << filt << std::endl;
     // }
     //hltSingleMu or hltSingleMuIso
     //if( filt == electronTriggerFilter_.encode() ) {
     if( ( (filt.find("Electron")!=std::string::npos) ||
	   (filt.find("hltSingleMu")!=std::string::npos) ||
	   (filt.find("hltL1IsoL1sMu")!=std::string::npos) ||
	   (filt.find("hltL2IsoL1sMu")!=std::string::npos) ||
	   (filt.find("hltL3IsoL1sMu")!=std::string::npos) ||
	   (filt.find("hltL1fL1sMu")!=std::string::npos) ||
	   (filt.find("hltL2fL1sMu")!=std::string::npos) ||
	   (filt.find("hltL3fL1sMu")!=std::string::npos) ||
	   (filt.find("hltIsoMu17")!=std::string::npos) ||
	   (filt.find("hltMuEta2p1")!=std::string::npos) ||
	   (filt.find("hltL1Mu0HTT")!=std::string::npos) ||
	   (filt.find("hltHT300")!=std::string::npos) ||
	   (filt.find("hlt")!=std::string::npos) ||
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

     MyTrigobj.id  = obj.id();
     MyTrigobj.px  = obj.px();
     MyTrigobj.py  = obj.py();
     MyTrigobj.pz  = obj.pz();
     MyTrigobj.pt  = obj.pt();
     MyTrigobj.eta = obj.eta();
     MyTrigobj.phi = obj.phi();
     MyTrigobj.et  = obj.et();
     MyTrigobj.energy = obj.energy();
     MyTrigobj.filter = filtKeys_string[i];

     //std::cout << " ===> HLT: filter = " << filtKeys_string[i] << ",\t id = " << obj.id() << ",\t pt = " << obj.pt() << ",\t eta = " << obj.eta() << ",\t phi = " << obj.phi() << std::endl;
     bnhltobjs->push_back(MyTrigobj);
   }



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

       //std::cout << " =====>  HLT: path name = " << currentTrigName << ",\t prescale = " << prescale << ",\t pass = " << accept << std::endl; 

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
   //int gtfeBx = -1;

   // get info from GTFE DAQ record
   const L1GtfeWord& gtfeWord = gtReadoutRecord->gtfeWord();
   //gtfeBx = gtfeWord.bxNr();
   int gtfeActiveBoards = gtfeWord.activeBoards();

   ///////////
   unsigned int pfIndexAlgo = 0;
   unsigned long long l1GtPfAlgoCacheID = iSetup.get<L1GtPrescaleFactorsAlgoTrigRcd>().cacheIdentifier();

   if (m_l1GtPfAlgoCacheID != l1GtPfAlgoCacheID) {

     edm::ESHandle<L1GtPrescaleFactors> l1GtPfAlgo;
     iSetup.get<L1GtPrescaleFactorsAlgoTrigRcd>().get(l1GtPfAlgo);
     m_l1GtPfAlgo = l1GtPfAlgo.product();

     m_prescaleFactorsAlgoTrig = & ( m_l1GtPfAlgo->gtPrescaleFactors() );

     m_l1GtPfAlgoCacheID = l1GtPfAlgoCacheID;
    }
    /////////////

   // get info from FDL if active (including decision word)
   if( isActive(gtfeActiveBoards,FDL) ) {
     /// get Global Trigger algo and technical triger bit statistics
     const DecisionWord& gtDecisionWord = gtReadoutRecord->decisionWord();
     const TechnicalTriggerWord& gtTTWord = gtReadoutRecord->technicalTriggerWord();
     //
     const std::vector<int>& prescaleFactorsAlgoTrig = ( *m_prescaleFactorsAlgoTrig ).at(pfIndexAlgo);

     // L1 algos
     for( CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
       BNtrigger MyTriggerL1TAlgo;
       int algoBitNumber = (algo->second).algoBitNumber();

       int prescaleFactor = prescaleFactorsAlgoTrig.at(algoBitNumber);

       MyTriggerL1TAlgo.pass = gtDecisionWord.at(algoBitNumber);
       MyTriggerL1TAlgo.name = (algo->second).algoName();
       MyTriggerL1TAlgo.prescale = prescaleFactor;

       //std::cout << " =====>  L1T algo: path name = " << (algo->second).algoName() << ",\t prescale = " << prescaleFactor << ",\t pass = " << gtDecisionWord.at(algoBitNumber) << std::endl; 

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

   ////////////////////////////////////
   ///////////////////////////////////
   ///////////////////////////////////
   // Isolated EM particles

   std::auto_ptr<BNtrigobjCollection> bnl1isoemobjs(new BNtrigobjCollection);

   Handle< l1extra::L1EmParticleCollection > isoEmColl ;
   iEvent.getByLabel( "l1extraParticles","Isolated", isoEmColl ) ;
   // std::cout << "Number of isolated EM " << isoEmColl->size() << std::endl ;

   for( l1extra::L1EmParticleCollection::const_iterator emItr = isoEmColl->begin(); emItr != isoEmColl->end(); ++emItr ){

     BNtrigobj MyL1obj;

     MyL1obj.px  = emItr->px();
     MyL1obj.py  = emItr->py();
     MyL1obj.pz  = emItr->pz();
     MyL1obj.pt  = emItr->pt();
     MyL1obj.eta = emItr->eta();
     MyL1obj.phi = emItr->phi();
     MyL1obj.et  = emItr->et();
     MyL1obj.energy = emItr->energy();
     MyL1obj.bx = emItr->bx();

     bnl1isoemobjs->push_back(MyL1obj);

     // std::cout << "  p4 (" << emItr->px()
     // 	       << ", " << emItr->py()
     // 	       << ", " << emItr->pz()
     // 	       << ", " << emItr->energy()
     // 	       << ") et " << emItr->et()
     // 	       << " eta " << emItr->eta()
     // 	       << " phi " << emItr->phi()
     // 	       << " pt " << emItr->pt()
     // 	       << " bx " << emItr->bx()
     // 	       << std::endl ;
   }

   // Non-isolated EM particles

   std::auto_ptr<BNtrigobjCollection> bnl1nonisoemobjs(new BNtrigobjCollection);

   Handle< l1extra::L1EmParticleCollection > nonIsoEmColl ;
   iEvent.getByLabel( "l1extraParticles","NonIsolated", nonIsoEmColl ) ;
   // std::cout << "Number of non-isolated EM " << nonIsoEmColl->size() << std::endl ;

   for( l1extra::L1EmParticleCollection::const_iterator emItr = nonIsoEmColl->begin(); emItr != nonIsoEmColl->end(); ++emItr ){

     BNtrigobj MyL1obj;

     MyL1obj.px  = emItr->px();
     MyL1obj.py  = emItr->py();
     MyL1obj.pz  = emItr->pz();
     MyL1obj.pt  = emItr->pt();
     MyL1obj.eta = emItr->eta();
     MyL1obj.phi = emItr->phi();
     MyL1obj.et  = emItr->et();
     MyL1obj.energy = emItr->energy();
     MyL1obj.bx = emItr->bx();

     bnl1nonisoemobjs->push_back(MyL1obj);

     // std::cout << "  p4 (" << emItr->px()
     // 	       << ", " << emItr->py()
     // 	       << ", " << emItr->pz()
     // 	       << ", " << emItr->energy()
     // 	       << ") et " << emItr->et()
     // 	       << " eta " << emItr->eta()
     // 	       << " phi " << emItr->phi()
     // 	       << " pt " << emItr->pt()
     // 	       << " bx " << emItr->bx()
     // 	       << std::endl ;
   }

   // Jet particles

   std::auto_ptr<BNtrigobjCollection> bnl1cenjetobjs(new BNtrigobjCollection);

   Handle< l1extra::L1JetParticleCollection > cenJetColl ;
   iEvent.getByLabel( "l1extraParticles","Central", cenJetColl ) ;
   // std::cout << "Number of central jets " << cenJetColl->size() << std::endl ;

   for( l1extra::L1JetParticleCollection::const_iterator jetItr = cenJetColl->begin(); jetItr != cenJetColl->end(); ++jetItr ){

     BNtrigobj MyL1obj;

     MyL1obj.px  = jetItr->px();
     MyL1obj.py  = jetItr->py();
     MyL1obj.pz  = jetItr->pz();
     MyL1obj.pt  = jetItr->pt();
     MyL1obj.eta = jetItr->eta();
     MyL1obj.phi = jetItr->phi();
     MyL1obj.et  = jetItr->et();
     MyL1obj.energy = jetItr->energy();
     MyL1obj.bx = jetItr->bx();

     bnl1cenjetobjs->push_back(MyL1obj);

     // std::cout << "  p4 (" << jetItr->px()
     // 	       << ", " << jetItr->py()
     // 	       << ", " << jetItr->pz()
     // 	       << ", " << jetItr->energy()
     // 	       << ") et " << jetItr->et()
     // 	       << " eta " << jetItr->eta()
     // 	       << " phi " << jetItr->phi()
     // 	       << " pt " << jetItr->pt()
     // 	       << " bx " << jetItr->bx()
     // 	       << std::endl ;
   }


   std::auto_ptr<BNtrigobjCollection> bnl1forjetobjs(new BNtrigobjCollection);

   Handle< l1extra::L1JetParticleCollection > forJetColl ;
   iEvent.getByLabel( "l1extraParticles","Forward", forJetColl ) ;
   // std::cout << "Number of forward jets " << forJetColl->size() << std::endl ;

   for( l1extra::L1JetParticleCollection::const_iterator jetItr = forJetColl->begin(); jetItr != forJetColl->end(); ++jetItr ){

     BNtrigobj MyL1obj;

     MyL1obj.px  = jetItr->px();
     MyL1obj.py  = jetItr->py();
     MyL1obj.pz  = jetItr->pz();
     MyL1obj.pt  = jetItr->pt();
     MyL1obj.eta = jetItr->eta();
     MyL1obj.phi = jetItr->phi();
     MyL1obj.et  = jetItr->et();
     MyL1obj.energy = jetItr->energy();
     MyL1obj.bx = jetItr->bx();

     bnl1forjetobjs->push_back(MyL1obj);

     // std::cout << "  p4 (" << jetItr->px()
     // 	       << ", " << jetItr->py()
     // 	       << ", " << jetItr->pz()
     // 	       << ", " << jetItr->energy()
     // 	       << ") et " << jetItr->et()
     // 	       << " eta " << jetItr->eta()
     // 	       << " phi " << jetItr->phi()
     // 	       << " pt " << jetItr->pt()
     // 	       << " bx " << jetItr->bx()
     // 	       << std::endl ;
   }


   std::auto_ptr<BNtrigobjCollection> bnl1taujetobjs(new BNtrigobjCollection);

   Handle< l1extra::L1JetParticleCollection > tauColl ;
   iEvent.getByLabel( "l1extraParticles","Tau", tauColl ) ;
   // std::cout << "Number of tau jets " << tauColl->size() << std::endl ;

   for( l1extra::L1JetParticleCollection::const_iterator tauItr = tauColl->begin(); tauItr != tauColl->end(); ++tauItr ){

     BNtrigobj MyL1obj;

     MyL1obj.px  = tauItr->px();
     MyL1obj.py  = tauItr->py();
     MyL1obj.pz  = tauItr->pz();
     MyL1obj.pt  = tauItr->pt();
     MyL1obj.eta = tauItr->eta();
     MyL1obj.phi = tauItr->phi();
     MyL1obj.et  = tauItr->et();
     MyL1obj.energy = tauItr->energy();
     MyL1obj.bx = tauItr->bx();

     bnl1taujetobjs->push_back(MyL1obj);

     // std::cout << "  p4 (" << tauItr->px()
     // 	       << ", " << tauItr->py()
     // 	       << ", " << tauItr->pz()
     // 	       << ", " << tauItr->energy()
     // 	       << ") et " << tauItr->et()
     // 	       << " eta " << tauItr->eta()
     // 	       << " phi " << tauItr->phi()
     // 	       << " pt " << tauItr->pt()
     // 	       << " bx " << tauItr->bx()
     // 	       << std::endl ;
   }

   // Muon particles

   std::auto_ptr<BNtrigobjCollection> bnl1muonobjs(new BNtrigobjCollection);

   Handle< l1extra::L1MuonParticleCollection > muColl ;
   iEvent.getByLabel( "l1extraParticles", muColl ) ;
   // std::cout << "Number of muons " << muColl->size() << std::endl ;

   for( l1extra::L1MuonParticleCollection::const_iterator muItr = muColl->begin(); muItr != muColl->end(); ++muItr ){

     BNtrigobj MyL1obj;

     MyL1obj.px  = muItr->px();
     MyL1obj.py  = muItr->py();
     MyL1obj.pz  = muItr->pz();
     MyL1obj.pt  = muItr->pt();
     MyL1obj.eta = muItr->eta();
     MyL1obj.phi = muItr->phi();
     MyL1obj.et  = muItr->et();
     MyL1obj.energy = muItr->energy();
     MyL1obj.bx = muItr->bx();

     MyL1obj.charge = muItr->charge();
     MyL1obj.isIsolated = muItr->isIsolated();
     MyL1obj.isMip = muItr->isMip();
     MyL1obj.isForward = muItr->isForward();
     MyL1obj.isRPC = muItr->isRPC();

     bnl1muonobjs->push_back(MyL1obj);

     // std::cout << "  q " << muItr->charge()
     // 	       << " p4 (" << muItr->px()
     // 	       << ", " << muItr->py()
     // 	       << ", " << muItr->pz()
     // 	       << ", " << muItr->energy()
     // 	       << ") et " << muItr->et()
     // 	       << " eta " << muItr->eta() << std::endl
     // 	       << "    phi " << muItr->phi()
     // 	       << "    pt " << muItr->pt()
     // 	       << " iso " << muItr->isIsolated()
     // 	       << " mip " << muItr->isMip()
     // 	       << " fwd " << muItr->isForward()
     // 	       << " rpc " << muItr->isRPC()
     // 	       << " bx " << muItr->bx()
     // 	       << std::endl ;
   }

   // MET

   std::auto_ptr<BNtrigobjCollection> bnl1metobjs(new BNtrigobjCollection);

   Handle< l1extra::L1EtMissParticleCollection > etMissColl ;
   iEvent.getByLabel( "l1extraParticles","MET", etMissColl ) ;

   BNtrigobj MyL1METobj;

   MyL1METobj.px  = etMissColl->begin()->px();
   MyL1METobj.py  = etMissColl->begin()->py();
   MyL1METobj.pz  = etMissColl->begin()->pz();
   MyL1METobj.pt  = etMissColl->begin()->pt();
   MyL1METobj.eta = etMissColl->begin()->eta();
   MyL1METobj.phi = etMissColl->begin()->phi();
   MyL1METobj.et  = etMissColl->begin()->et();
   MyL1METobj.energy = etMissColl->begin()->energy();
   MyL1METobj.bx = etMissColl->begin()->bx();

   MyL1METobj.etTotal = etMissColl->begin()->etTotal();

   bnl1metobjs->push_back(MyL1METobj);

   // std::cout << "MET Coll (" << etMissColl->begin()->px()
   // 	     << ", " << etMissColl->begin()->py()
   // 	     << ", " << etMissColl->begin()->pz()
   // 	     << ", " << etMissColl->begin()->energy()
   // 	     << ") phi " << etMissColl->begin()->phi()
   // 	     << " pt " << etMissColl->begin()->pt()
   // 	     << " bx " << etMissColl->begin()->bx()
   // 	     << " EtTot " << etMissColl->begin()->etTotal()
   // 	     << std::endl ;

   // MHT

   std::auto_ptr<BNtrigobjCollection> bnl1mhtobjs(new BNtrigobjCollection);

   Handle< l1extra::L1EtMissParticleCollection > htMissColl ;
   iEvent.getByLabel( "l1extraParticles","MHT", htMissColl ) ;

   BNtrigobj MyL1MHTobj;

   MyL1MHTobj.px  = htMissColl->begin()->px();
   MyL1MHTobj.py  = htMissColl->begin()->py();
   MyL1MHTobj.pz  = htMissColl->begin()->pz();
   MyL1MHTobj.pt  = htMissColl->begin()->pt();
   MyL1MHTobj.eta = htMissColl->begin()->eta();
   MyL1MHTobj.phi = htMissColl->begin()->phi();
   MyL1MHTobj.et  = htMissColl->begin()->et();
   MyL1MHTobj.energy = htMissColl->begin()->energy();
   MyL1MHTobj.bx = htMissColl->begin()->bx();

   MyL1MHTobj.etTotal = htMissColl->begin()->etTotal();

   bnl1mhtobjs->push_back(MyL1MHTobj);

   // std::cout << "MHT Coll (" << htMissColl->begin()->px()
   // 	     << ", " << htMissColl->begin()->py()
   // 	     << ", " << htMissColl->begin()->pz()
   // 	     << ", " << htMissColl->begin()->energy()
   // 	     << ") phi " << htMissColl->begin()->phi()
   // 	     << " pt " << htMissColl->begin()->pt()
   // 	     << " bx " << htMissColl->begin()->bx()
   // 	     << " HtTot " << htMissColl->begin()->etTotal()
   // 	     << std::endl ;

   /*
   // HF Rings
   Handle< l1extra::L1HFRingsCollection > hfRingsColl ;
   iEvent.getByLabel( "l1extraParticles", hfRingsColl ) ;
   std::cout << "HF Rings:" << std::endl ;
   for( int i = 0 ; i < L1HFRings::kNumRings ; ++i ){
     std::cout << "  " << i << ": et sum = "
	       << hfRingsColl->begin()->hfEtSum( (L1HFRings::HFRingLabels) i )
	       << ", bit count = "
	       << hfRingsColl->begin()->hfBitCount( (L1HFRings::HFRingLabels) i )
	       << std::endl ;
   }
   std::cout << std::endl ;
   */

   //////////////////////////////
   ////////////////////
   ///////////
   ////

   if( GoodVertex ) numGoodVertexEvents++;
   if( FilterOutScraping ) numFilterOutScrapingEvents++;


   // Put the collections into the event
   iEvent.put(bnevent);
   if( produceElectron ) iEvent.put(bnelectrons,eleTag_.label());
   if( producePFelectron ) iEvent.put(bnpfelectrons,pfeleTag_.label());
   if( producePFelectronLoose ) iEvent.put(bnpfelectronsLoose,pfeleLooseTag_.label());
   if( produceCaloJet ) iEvent.put(bncalojets,calojetTag_.label());
   if( producePFJet ) iEvent.put(bnpfjets,pfjetTag_.label());
   if( produceGenJet ) iEvent.put(bngenjets,genjetTag_.label());
   //if( produceCaloMET ) iEvent.put(bncalomet,calometTag_.label());
   if( producePFMET ) iEvent.put(bnpfmet,pfmetTag_.label());
   if( producePFMET_type1correctedRECO ) iEvent.put(bnpfmet_type1correctedRECO,std::string(pfmetTag_type1correctedRECO_.label() + "BN"));
   if( producePFMET_uncorrectedPF2PAT )  iEvent.put(bnpfmet_uncorrectedPF2PAT,std::string(pfmetTag_uncorrectedPF2PAT_.label() + "BN"));
   if( producePFMET_uncorrectedRECO )    iEvent.put(bnpfmet_uncorrectedRECO,std::string(pfmetTag_uncorrectedRECO_.label() + "BN"));
   if( produceTCMET ) iEvent.put(bntcmet,tcmetTag_.label());
   if( produceMuon ) iEvent.put(bnmuons,muonTag_.label());
   if( producePFMuon ) iEvent.put(bnpfmuons,pfmuonTag_.label());
   if( producePFMuonLoose ) iEvent.put(bnpfmuonsLoose,pfmuonLooseTag_.label());
   if( produceCocktailMuon ) iEvent.put(bncocktailmuons,cocktailmuonTag_.label());
   if( produceTau ) iEvent.put(bntaus,tauTag_.label());
   if( producePhoton ) iEvent.put(bnphotons,photonTag_.label());
   if( produceSCsEB || produceSCsEE ) iEvent.put(bnsuperclusters,kSC);
   if( produceTrack ) iEvent.put(bntracks,trackTag_.label());
   iEvent.put(bntrigger,kHLT);
   if( produceGenParticle ) iEvent.put(bnmcparticles,kMCpar);
   if( produceGenParticle ) iEvent.put(bnmcelectrons,kMCele);
   if( produceGenParticle ) iEvent.put(bnmcmuons,kMCmu);
   iEvent.put(bntriggerl1talgo,kL1Talgo);
   iEvent.put(bntriggerl1ttech,kL1Ttech);
   iEvent.put(bnhltobjs,kHLTobj);
   iEvent.put(bnl1isoemobjs,kL1EmParticlesIso);
   iEvent.put(bnl1nonisoemobjs,kL1EmParticlesNonIso);
   iEvent.put(bnl1metobjs,kL1EtMissParticlesMET);
   iEvent.put(bnl1mhtobjs,kL1EtMissParticlesMHT);
   iEvent.put(bnl1cenjetobjs,kL1JetParticlesCentral);
   iEvent.put(bnl1forjetobjs,kL1JetParticlesForward);
   iEvent.put(bnl1taujetobjs,kL1JetParticlesTau);
   iEvent.put(bnl1muonobjs,kL1MuonParticles);
   iEvent.put(bnpvs,pvTag_.label());

   iEvent.put(bnbxlumis,kBXlumi);


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

  // tracker geometry used for hit information in BNtrack
  edm::ESHandle<TrackerGeometry> tkGeom;
  c.get<TrackerDigiGeometryRecord>().get( tkGeom );
  m_tracker = tkGeom.product();

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

// === Tau Crack veto === //
bool BEANmaker::tauIsInTheCracks(float etaValue){
	return (fabs(etaValue) < 0.018 ||  
			(fabs(etaValue)>0.423 && fabs(etaValue)<0.461) ||
			(fabs(etaValue)>0.770 && fabs(etaValue)<0.806) ||
			(fabs(etaValue)>1.127 && fabs(etaValue)<1.163) ||
			(fabs(etaValue)>1.460 && fabs(etaValue)<1.558));
}

//define this as a plug-in
DEFINE_FWK_MODULE(BEANmaker);
