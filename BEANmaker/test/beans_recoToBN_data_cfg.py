import FWCore.ParameterSet.Config as cms

process = cms.Process("BEAN")

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

runOnMC = False

process.GlobalTag.globaltag = cms.string('GR_R_42_V21::All') ### global tag for data
if runOnMC == True:
    process.GlobalTag.globaltag = cms.string('MC_42_V13::All') ## global tag for mc


# Get a list of good primary vertices, in 42x, these are DAF vertices
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0), maxRho = cms.double(2.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )


# add pf met
from PhysicsTools.PatAlgos.tools.metTools import *

# tag information not available in AOD (but discriminator there)
process.patJets.addTagInfos = cms.bool(False)

# get the jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *


inputJetCorrLabelPF = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']))
inputJetCorrLabelCalo = ('AK5Calo', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']))

jetMetCorrections = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])

if runOnMC == False:
    jetMetCorrections = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
    inputJetCorrLabelPF = ('AK5PF', jetMetCorrections)
    inputJetCorrLabelCalo = ('AK5Calo', jetMetCorrections)


##################################
#
## load the PAT config
#process.load("PhysicsTools.PatAlgos.patSequences_cff")
#
## Configure PAT to use PF2PAT instead of AOD sources
## this function will modify the PAT sequences. It is currently 
## not possible to run PF2PAT+PAT and standart PAT at the same time
from PhysicsTools.PatAlgos.tools.pfTools import *
#
## An empty postfix means that only PF2PAT is run,
## otherwise both standard PAT and PF2PAT are run. In the latter case PF2PAT
## collections have standard names + postfix (e.g. patElectronPFlow)  
#
postfix = "PF"
pfSuffix = 'PF'
jetAlgo = "AK5"
usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo, runOnMC=runOnMC, postfix=postfix, jetCorrections=('AK5PFchs', jetMetCorrections))


if runOnMC == False:
    # removing MC matching for standard PAT sequence
    # for the PF2PAT+PAT sequence, it is done in the usePF2PAT function
    removeMCMatching( process, ['All'] ) 
    removeMCMatchingPF2PAT( process, '' ) 

#addPfMET(process, 'PF')  ## included in PF2PAT
addTcMET(process,'TC')
process.patMETsAK5Calo = process.patMETs


# addPfMET TypeI
from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5PFJet 
process.metJESCorAK5PFTypeI = metJESCorAK5PFJet.clone( 
    inputUncorJetsLabel = "patJetsPF", 
    metType = "pat",                  
    inputUncorMetLabel = "pfMet",
    )
process.patMETsTypeIPF = process.patMETsPF.clone(
    metSource = cms.InputTag("metJESCorAK5PFTypeI")
    )
# Add to producersLayer1 sequence
process.patDefaultSequencePF.replace(
    process.patMETsPF,
    process.patMETsPF+
    process.metJESCorAK5PFTypeI+
    process.patMETsTypeIPF
    )


# Add AK5 Calo jets
addJetCollection(process,cms.InputTag('ak5CaloJets'),
                 'AK5', 'Calo',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = inputJetCorrLabelCalo,
#                 jetCorrLabel = ('AK5Calo', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
                 doType1MET   = True,
                 doL1Cleaning = True,
                 #doAreaFastjet = True,
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True
                 )

# Add PF jets
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = inputJetCorrLabelPF,
#                 jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
                 doType1MET   = True,
                 doL1Cleaning = True,
                 #doAreaFastjet = True,
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True
                 )

#
## to run second PF2PAT+PAT with differnt postfix uncomment the following lines
## and add it to path
##postfix2 = "PFlow2"
##jetAlgo2="AK7"
##usePF2PAT(process,runPF2PAT=True, jetAlgo=jetAlgo2, runOnMC=True, postfix=postfix2)
#
## to use tau-cleaned jet collection uncomment the following:
##getattr(process,"pfNoTau"+postfix).enable = True 
#
## to switch default tau to HPS tau uncomment the following:
##adaptPFTaus(process,"hpsPFTau",postfix=postfix)
#
## subtract charged hadronic pile-up particles (from wrong PVs)
## effects also JECs
usePFnoPU = True # before any top projection
#
## other switches for PF top projections (default: all 'True')
useNoMuon     = True # before electron top projection
useNoElectron = True # before jet top projection
useNoJet      = True # before tau top projection
#useNoTau      = True # before MET top projection
useNoTau      = False # before MET top projection
#
pfVertices  = 'goodOfflinePrimaryVertices'
pfD0Cut     = 0.2
pfDzCut     = 0.5
#
useMuonCutBasePF = False # use minimal (veto) muon selection cut on top of 'pfMuonSelectionCut'
##pfMuonIsoConeR   = 0.4
##pfMuonCombIsoCut = 0.2
## electrons
##pfElectronSelectionCut  = 'pt > 5. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits < 2'
useElectronCutBasePF = False # use minimal (veto) electron selection cut on top of 'pfElectronSelectionCut'
##pfElectronIsoConeR   = 0.3
##pfElectronCombIsoCut = 0.2
#
#### JEC levels
#
## levels to be accessible from the jets
## jets are corrected to L3Absolute (MC), L2L3Residual (data) automatically, if enabled here
## and remain uncorrected, if none of these levels is enabled here
useL1FastJet    = True  # needs useL1Offset being off, error otherwise
useL1Offset     = False # needs useL1FastJet being off, error otherwise
useL2Relative   = True
useL3Absolute   = True
useL2L3Residual = True  # takes effect only on data
useL5Flavor     = True
useL7Parton     = True
#
applyPostfix( process, 'pfNoPileUp'  , postfix ).enable = usePFnoPU
applyPostfix( process, 'pfNoMuon'    , postfix ).enable = useNoMuon
applyPostfix( process, 'pfNoElectron', postfix ).enable = useNoElectron
applyPostfix( process, 'pfNoJet'     , postfix ).enable = useNoJet
applyPostfix( process, 'pfNoTau'     , postfix ).enable = useNoTau
applyPostfix( process, 'pfPileUp', postfix ).Vertices = cms.InputTag( pfVertices )
if useL1FastJet:
    applyPostfix( process, 'pfPileUp', postfix ).checkClosestZVertex = False
    applyPostfix( process, 'pfJets', postfix ).doAreaFastjet = True
    applyPostfix( process, 'pfJets', postfix ).doRhoFastjet  = False



process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets.doRhoFastjet = True
#    #process.kt6PFJets.Rho_EtaMax = cms.double(3.0)
#
## apply FastJet corrections only if demanded
#    #print jetMetCorrections
process.pfJets.doAreaFastjet = True
#        #process.pfJets.Rho_EtaMax = cms.double(3.0)
#
process.pfJetsPF.doAreaFastjet = True
#
##process.kt6PFJets = kt6PFJets.clone( src          = cms.InputTag( 'particleFlow' )
# #                                    , doRhoFastjet = True
#  #                                   )
#
process.patDefaultSequence.replace(process.patJetCorrFactors,
                                   process.kt6PFJets + process.patJetCorrFactors)


process.patMuons.usePV = cms.bool(False)
process.patElectrons.usePV = cms.bool(False)
#
#
## top projections in PF2PAT:
#getattr(process,"pfNoPileUp"+postfix).enable = True 
#getattr(process,"pfNoMuon"+postfix).enable = True 
#getattr(process,"pfNoElectron"+postfix).enable = True 
#getattr(process,"pfNoTau"+postfix).enable = False 
#getattr(process,"pfNoJet"+postfix).enable = True
#
## verbose flags for the PF2PAT modules
#getattr(process,"pfNoMuon"+postfix).verbose = False
#
## enable delta beta correction for muon selection in PF2PAT? 
#getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = False

##################################


process.skimMuon = cms.EDFilter("SkimMuon",
                                ptCutMu1 = cms.untracked.double(30),
                                ptCutMu2 = cms.untracked.double(5),
                                caloIsoCut = cms.untracked.double(5),
                                trackIsoCut = cms.untracked.double(3),
                                minNumGoodMu = cms.untracked.int32(2)
                                )

# require physics declared
process.load('HLTrigger.special.hltPhysicsDeclared_cfi')
process.hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'

# require scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.2)
                                    )


# switch on PAT trigger
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
switchOnTrigger( process )

process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24), 
                                           maxd0 = cms.double(2) 
                                           )

# Select jets
#process.selectedPatJets.cut = cms.string('pt > 10')
#process.selectedPatJetsAK5Calo.cut = cms.string('pt > 10')
#process.selectedPatJetsAK5PF.cut = cms.string('pt > 10')
#process.selectedPatJetsAK5JPT.cut = cms.string('pt > 10')


process.BNproducer = cms.EDProducer('BEANmaker',
        calometTag = cms.InputTag("patMETsAK5Calo"), 
        pfmetTag = cms.InputTag("patMETsTypeIPF"), 
        tcmetTag = cms.InputTag("patMETsTC"), 
        eleTag = cms.InputTag("selectedPatElectrons"),
        pfeleTag = cms.InputTag("selectedPatElectronsPF"),
        genParticleTag = cms.InputTag("genParticles"),
        calojetTag = cms.InputTag("selectedPatJetsAK5Calo"), 
        pfjetTag = cms.InputTag("selectedPatJetsAK5PF"), 
        jptjetTag = cms.InputTag("selectedPatJetsAK5JPT"), 
        muonTag = cms.InputTag("selectedPatMuons"),
        pfmuonTag = cms.InputTag("selectedPatMuonsPF"),
        cocktailmuonTag = cms.InputTag("refitMuons"),
        photonTag = cms.InputTag("selectedPatPhotons"),
        EBsuperclusterTag = cms.InputTag("correctedHybridSuperClusters"),
        EEsuperclusterTag = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
        reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
        reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
        trackTag = cms.InputTag("generalTracks"),
        triggerResultsTag = cms.InputTag("TriggerResults::HLT"),
        gtSource = cms.InputTag("gtDigis"),
        pvTag = cms.InputTag("offlinePrimaryVertices"),
        triggerSummaryTag = cms.InputTag("hltTriggerSummaryAOD"),
        dcsTag = cms.InputTag("scalersRawToDigi"),
        hltProcessName = cms.string( "HLT"),  
        eventWeight = cms.double(1),
        minSCEt = cms.double(10),
        minPhotonEt = cms.double(10),
        minJetPt = cms.double(10),
        minTrackPt = cms.double(10),
        verbose = cms.bool(True),
        sample = cms.int32(-1)
)


# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [

        'file:/tmp/puigh/EE0408B7-B8F8-E011-A4A9-001D09F244DE.root',
#	'/store/data/Run2011B/SingleMu/RECO/PromptReco-v1/000/178/479/EE0408B7-B8F8-E011-A4A9-001D09F244DE.root',
#	'/store/data/Run2011B/SingleMu/AOD/PromptReco-v1/000/177/515/F6FFF2B0-F0EC-E011-9060-00215AEDFD98.root',

#        'file:/tmp/puigh/8A118A9A-83AC-E011-8BA6-003048C69406.root'
#        'rfio:/castor/cern.ch/user/p/puigh/CUJEM/8A118A9A-83AC-E011-8BA6-003048C69406.root'
#        '/store/mc/Spring11/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S1_START311_V1G1-v1/0011/F0F54048-4A50-E011-9CDE-003048D47792.root',
#         'file:/tmp/puigh/F0F54048-4A50-E011-9CDE-003048D47792.root',
#        '/store/mc/Spring11/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S1_START311_V1G1-v1/0011/106AFBB1-2D50-E011-A529-00E0817918AD.root',
#        '/store/mc/Spring11/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S1_START311_V1G1-v1/0010/88EA5113-B04F-E011-9704-003048D45F48.root',
#        '/store/mc/Spring11/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S1_START311_V1G1-v1/0010/76C441E2-AB4F-E011-A71B-003048D45FE2.root',
#        '/store/mc/Spring11/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S1_START311_V1G1-v1/0010/44C4D221-A94F-E011-AFA2-0025B31E3C00.root',

        ] );
process.source.fileNames = readFiles

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

## HBHE Noise Filter
#process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
#process.p = cms.Path(process.HBHENoiseFilterResultProducer)

process.load('RecoMuon/MuonIdentification/refitMuons_cfi')

process.p = cms.Path(
#    #process.scrapingVeto*
#    #process.primaryVertexFilter*
##    process.HBHENoiseFilter*
    process.goodOfflinePrimaryVertices*
    process.skimMuon*
    process.refitMuons*
    process.HBHENoiseFilterResultProducer*
    (getattr(process,"patPF2PATSequence"+postfix)+
    process.patDefaultSequence)*
    process.BNproducer
    )



# rename output file
process.out.fileName = cms.untracked.string('reco_7TeV_426_BN_v3.root')

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process all the events
process.maxEvents.input = 10000
process.options.wantSummary = True

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent

process.out.outputCommands = [ 'drop *' ]

process.out.outputCommands.extend( [ # BEAN Objects
                                    'keep *_BNproducer_*_*',
				    #'keep *',
				    #'keep recoMuons_refitMuons_*_*',
                                     ] )
