import FWCore.ParameterSet.Config as cms

process = cms.Process("BEAN")

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

### global tag for data
#process.GlobalTag.globaltag = cms.string('GR_R_38X_V7::All')

## global tag for mc
process.GlobalTag.globaltag = cms.string('MC_42_V13::All')

# add pf met
from PhysicsTools.PatAlgos.tools.metTools import *
#removeMCMatching(process, ['All'])
addPfMET(process, 'PF')
addTcMET(process,'TC')
process.patMETsAK5Calo = process.patMETs

# tag information not available in AOD (but discriminator there)
process.patJets.addTagInfos = cms.bool(False)

# get the jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *
#switchJECSet( process, "Spring11")

# Add AK5 Calo jets
addJetCollection(process,cms.InputTag('ak5CaloJets'),
                 'AK5', 'Calo',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5Calo', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
#                 jetCorrLabel = ('AK5Calo', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
                 doType1MET   = True,
                 doL1Cleaning = True,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True
                 )

# Add PF jets
addJetCollection(process,cms.InputTag('ak5PFJets'),
                 'AK5', 'PF',
                 doJTA        = True,
                 doBTagging   = True,
                 jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])),
#                 jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset', 'L2Relative', 'L3Absolute'])),
                 doType1MET   = True,
                 doL1Cleaning = True,                 
                 doL1Counters = False,
                 genJetCollection=cms.InputTag("ak5GenJets"),
                 doJetID      = True
                 )

## Add JPT jets
#addJetCollection(process,cms.InputTag('JetPlusTrackZSPCorJetAntiKt5'),
#                 'AK5', 'JPT',
#                 doJTA        = True,
#                 doBTagging   = True,
#                 jetCorrLabel = ('AK5JPT', cms.vstring(['L1FastJet', 'L1JPTOffset', 'L2Relative', 'L3Absolute'])),
#                 doType1MET   = True,
#                 doL1Cleaning = True,                 
#                 doL1Counters = False,
#                 genJetCollection=cms.InputTag("ak5GenJets"),
#                 doJetID      = True
#                 )



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
        pfmetTag = cms.InputTag("patMETsPF"), 
        tcmetTag = cms.InputTag("patMETsTC"), 
        eleTag = cms.InputTag("selectedPatElectrons"),
        genParticleTag = cms.InputTag("genParticles"),
        calojetTag = cms.InputTag("selectedPatJetsAK5Calo"), 
        pfjetTag = cms.InputTag("selectedPatJetsAK5PF"), 
        jptjetTag = cms.InputTag("selectedPatJetsAK5JPT"), 
        muonTag = cms.InputTag("selectedPatMuons"),
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
        sample = cms.int32(2323)
)


# Add the files 
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

readFiles.extend( [
        'file:/tmp/puigh/8A118A9A-83AC-E011-8BA6-003048C69406.root'
#        'rfio:/castor/cern.ch/user/p/puigh/CUJEM/8A118A9A-83AC-E011-8BA6-003048C69406.root'
#        '/store/mc/Spring11/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S1_START311_V1G1-v1/0011/F0F54048-4A50-E011-9CDE-003048D47792.root',
#        '/store/mc/Spring11/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S1_START311_V1G1-v1/0011/106AFBB1-2D50-E011-A529-00E0817918AD.root',
#        '/store/mc/Spring11/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S1_START311_V1G1-v1/0010/88EA5113-B04F-E011-9704-003048D45F48.root',
#        '/store/mc/Spring11/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S1_START311_V1G1-v1/0010/76C441E2-AB4F-E011-A71B-003048D45FE2.root',
#        '/store/mc/Spring11/WZtoAnything_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S1_START311_V1G1-v1/0010/44C4D221-A94F-E011-AFA2-0025B31E3C00.root',

        ] );
process.source.fileNames = readFiles

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")


process.p = cms.Path(
    #process.scrapingVeto*
    #process.primaryVertexFilter*
    process.patDefaultSequence*
    process.BNproducer
    )

# rename output file
process.out.fileName = cms.untracked.string('reco_7TeV_426_BN.root')

# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process all the events
process.maxEvents.input = -1
process.options.wantSummary = True

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerEventContent

process.out.outputCommands = [ 'drop *' ]

process.out.outputCommands.extend( [ # BEAN Objects
                                    'keep *_BNproducer_*_*',
                                     ] )
