import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
import os


####################################################################
## Global job options

REPORTEVERY = 100
WANTSUMMARY = True



####################################################################
## Define the process

process = cms.Process('BEANs')#"topDileptonNtuple")

#SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )


op_runOnMC = True#
op_runOnAOD = True#
op_globalTag = ''#
op_mode = ''#
op_samplename = 'ttbarz'#
op_inputScript = 'TopAnalysis.Configuration.Summer12.TTH_Inclusive_M_130_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1_cff'#
#op_inputScript = 'TopAnalysis.Configuration.Summer12.WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v2_cff'
op_outputFile = 'ttbarZ.root'#
op_systematicsName = 'Nominal'#
op_json = ''#
op_skipEvents = 0#
op_includePDFWeights = False#
op_maxEvents = 30#
####################################################################
## Set up samplename

if op_samplename == '':
    print 'cannot run without specifying a samplename'
    exit(8888)

if op_samplename == 'data':
    op_runOnMC = False



####################################################################
## Define input

if op_inputScript != '':
    process.load(op_inputScript)
    #inputFiles = cms.untracked.vstring('file:ttH_MC.root')
    #process.source = cms.Source("PoolSource", fileNames = inputFiles, secondaryFileNames = cms.untracked.vstring())
else:
    print 'need an input script'
    exit(8889)

print "max events: ", op_maxEvents
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(op_maxEvents)
)

if op_skipEvents > 0:
    process.source.skipEvents = cms.untracked.uint32(op_skipEvents)

# Limit to json file (if passed as parameter)
if op_json != '':
    import FWCore.PythonUtilities.LumiList as LumiList
    import FWCore.ParameterSet.Types as CfgTypes
    myLumis = LumiList.LumiList(filename = op_json).getCMSSWString().split(',')
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)



####################################################################
## Create output path

if op_outputFile == '':
    fn = op_mode + '_test.root'
else:
    fn = op_outputFile

#process.TFileService = cms.Service("TFileService",
##    fileName = cms.string(fn)
#)



####################################################################
## Configure message logger

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = REPORTEVERY

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(WANTSUMMARY)
)



####################################################################
## Geometry and Detector Conditions

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if op_globalTag != '':
    print "Setting global tag to the command-line value"
    process.GlobalTag.globaltag = cms.string( op_globalTag )
else:
    print "Determine global tag automatically"
    if op_runOnMC:
        process.GlobalTag.globaltag = cms.string('START53_V7G::All')#('START53_V26::All')
    else:
	process.GlobalTag.globaltag = cms.string('FT_53_V21_AN6::All')

print "Using global tag: ", process.GlobalTag.globaltag

process.load("Configuration.StandardSequences.MagneticField_cff")



####################################################################
## HCAL Noise filter

process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
process.HBHENoiseFilter.minIsolatedNoiseSumE = cms.double(999999.)
process.HBHENoiseFilter.minNumIsolatedNoiseChannels = cms.int32(999999)
process.HBHENoiseFilter.minIsolatedNoiseSumEt = cms.double(999999.)



####################################################################
## Beam scraping filter

process.scrapingFilter = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn     = cms.untracked.bool(False),
    numtrack    = cms.untracked.uint32(10),
    thresh      = cms.untracked.double(0.25)
    )



####################################################################
## ECAL laser correction filter

process.load("RecoMET.METFilters.ecalLaserCorrFilter_cfi")



####################################################################
## Primary vertex filtering

from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter(
    "PrimaryVertexObjectFilter",
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) ),
    src=cms.InputTag('offlinePrimaryVertices')
    )



####################################################################
## Trigger filtering

# Get the central diLepton trigger lists, and set up filter
from TopAnalysis.TopFilter.sequences.diLeptonTriggers_cff import *
process.load("TopAnalysis.TopFilter.filters.TriggerFilter_cfi")
process.filterTrigger.TriggerResults = cms.InputTag('TriggerResults','','HLT')
process.filterTrigger.printTriggers = False
if op_mode == 'mumu':
    process.filterTrigger.hltPaths  = mumuTriggers
elif op_mode == 'emu':
    process.filterTrigger.hltPaths  = emuTriggers
elif op_mode == 'ee':
    process.filterTrigger.hltPaths  = eeTriggers
else:
    process.filterTrigger.hltPaths = eeTriggers + emuTriggers + mumuTriggers
    
#print "Printing triggers: ", process.filterTrigger.printTriggers



####################################################################
## Jet corrections

if op_runOnMC:
    jetCorr = ('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute'])
else:
    jetCorr = ('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute', 'L2L3Residual'])



####################################################################
## PF2PAT sequence

# process.load("Configuration.EventContent.EventContent_cff")
# process.out = cms.OutputModule("PoolOutputModule",
#     #fileName = cms.untracked.string(".root"),
#     process.FEVTEventContent,
#     #dataset = cms.untracked.PSet(dataTier = cms.untracked.string('RECO')),
#     fileName = cms.untracked.string("ttbarz.root"),
# )

process.load( "TopQuarkAnalysis.Configuration.patRefSel_outputModule_cff" )
process.out.fileName = cms.untracked.string("merged_bean.root")
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out.outputCommands +=patEventContent
process.out.SelectEvents.SelectEvents = []



process.load("PhysicsTools.PatAlgos.patSequences_cff")





#pfpostfix = "PFlow"
pfpostfix = ""

from PhysicsTools.PatAlgos.tools.pfTools import *

usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=op_runOnMC, postfix=pfpostfix, jetCorrections=jetCorr, pvCollection=cms.InputTag('goodOfflinePrimaryVertices'), typeIMetCorrections=True) 

from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *

# Produce pat trigger content
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff")

#turn off PF2PAT Top Projection....Charlie Projection!

applyPostfix( process, 'pfNoPileUp'  , pfpostfix ).enable = True
applyPostfix( process, 'pfNoMuon'    , pfpostfix ).enable = False
applyPostfix( process, 'pfNoElectron', pfpostfix ).enable = False
applyPostfix( process, 'pfNoJet'     , pfpostfix ).enable = False
applyPostfix( process, 'pfNoTau'     , pfpostfix ).enable = False

####################################################################
## Set up selections for PF2PAT & PAT objects: Electrons

# MVA ID
process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(process.mvaTrigV0)
getattr(process,'patElectrons'+pfpostfix).electronIDSources.mvaTrigV0 = cms.InputTag("mvaTrigV0")
getattr(process, 'patPF2PATSequence'+pfpostfix).replace(getattr(process,'patElectrons'+pfpostfix),
                                                process.eidMVASequence *
                                                getattr(process,'patElectrons'+pfpostfix)
                                                )

process.pfPileUp.checkClosestZVertex = cms.bool(False)

process.pfSelectedElectrons.cut = 'pt > 5.'# && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits < 2'

# Switch isolation cone to 0.3 and set cut to 0.15
process.pfIsolatedElectrons.doDeltaBetaCorrection = True   # not really a 'deltaBeta' correction, but it serves
process.pfIsolatedElectrons.deltaBetaIsolationValueMap = cms.InputTag("elPFIsoValuePU03PFId")
process.pfIsolatedElectrons.isolationValueMapsCharged = cms.VInputTag(cms.InputTag("elPFIsoValueCharged03PFId"))
process.pfIsolatedElectrons.isolationValueMapsNeutral = cms.VInputTag(cms.InputTag("elPFIsoValueNeutral03PFId"), cms.InputTag("elPFIsoValueGamma03PFId"))
process.pfIsolatedElectrons.isolationCut = 99.0#0.2

process.pfElectronsFromVertex.d0Cut = 99.0
process.pfElectronsFromVertex.dzCut = 99.0

process.patElectrons.isolationValues = cms.PSet(
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFId"),
    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFId"),
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFId"),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFId"),
    pfPhotons = cms.InputTag("elPFIsoValueGamma03PFId") )


process.selectedPatElectrons.cut = ''#'electronID("mvaTrigV0") > 0.5 && passConversionVeto'
            #cant do this on python level :-(
            #' && abs(gsfTrack().dxy(vertex_.position())) < 0.04'


#process.selectedPatElectronsAfterScaling = selectedPatElectrons.clone(
#    src = 'scaledJetEnergy:selectedPatElectrons',
#    cut = 'pt > 20 && abs(eta) < 2.5'
#)



####################################################################
## Set up selections for PF2PAT & PAT objects: Muons

process.pfSelectedMuons.cut = 'pt > 5.'


# Switch isolation cone to 0.3 and set cut to 0.15
process.pfIsolatedMuons.doDeltaBetaCorrection = True
process.pfIsolatedMuons.deltaBetaIsolationValueMap = cms.InputTag("muPFIsoValuePU03", "", "")
process.pfIsolatedMuons.isolationValueMapsCharged = [cms.InputTag("muPFIsoValueCharged03")]
process.pfIsolatedMuons.isolationValueMapsNeutral = [cms.InputTag("muPFIsoValueNeutral03"), cms.InputTag("muPFIsoValueGamma03")]
process.pfIsolatedMuons.isolationCut = 99.0#0.2

process.pfMuonsFromVertex.d0Cut = 99.0
process.pfMuonsFromVertex.dzCut = 99.0

process.patMuons.isolationValues = cms.PSet(
        pfNeutralHadrons = cms.InputTag("muPFIsoValueNeutral03"),
        pfChargedAll = cms.InputTag("muPFIsoValueChargedAll03"),
        pfPUChargedHadrons = cms.InputTag("muPFIsoValuePU03"),
        pfPhotons = cms.InputTag("muPFIsoValueGamma03"),
        pfChargedHadrons = cms.InputTag("muPFIsoValueCharged03")
        )


process.selectedPatMuons.cut = ''#'isGlobalMuon && pt > 20 && abs(eta) < 2.5'



####################################################################
## Set up selections for PF2PAT & PAT objects: Jets

#process.selectedPatJets.cut = 'abs(eta)<5.4'




## taus
tauCut                 = 'pt > 5. && abs(eta) < 2.5 && tauID("decayModeFinding")'
process.selectedPatTaus.cut = tauCut


process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")


####################################################################
## Basic debugging analyzer

#process.load("TopAnalysis.TopAnalyzer.CheckDiLeptonAnalyzer_cfi")
#process.analyzeDiLepton.electrons = 'fullySelectedPatElectronsCiC'
#process.analyzeDiLepton.muons = 'fullySelectedPatMuons'



####################################################################
## Set up sample-specific flags for individual treatment in nTuple production

zGenInfo = False
zproducer = False
topfilter = False
signal = False
higgsSignal = False
alsoViaTau = False
ttbarV = False

if op_samplename == 'ttbarsignal':
    topfilter = True
    signal = True
    viaTau = False
elif op_samplename == 'ttbarsignalviatau':
    topfilter = True
    signal = True
    viaTau = True
elif op_samplename == 'ttbarsignalplustau':
    topfilter = True
    signal = True
    viaTau = False
    alsoViaTau = True
elif op_samplename == 'ttbarbg':
    topfilter = True
elif op_samplename == 'dy1050' or op_samplename == 'dy50inf':
    zproducer = True
    zGenInfo = True
elif op_samplename == 'ttbarhiggstobbbar' or op_samplename == 'ttbarhiggsinclusive':
    topfilter = True
    signal = True
    viaTau = False
    alsoViaTau = True
    higgsSignal = True
elif op_samplename == 'gghiggstozzto4l' or op_samplename == 'vbfhiggstozzto4l':
    zGenInfo = True
    higgsSignal = True
elif op_samplename == 'ttbarw' or op_samplename == 'ttbarz':
    topfilter = True
    signal = True
    viaTau = False
    alsoViaTau = True
    ttbarV = True
elif op_samplename in ['data', 'singletop', 'singleantitop','ww',
        'wz','zz','wjets',
        'qcdmu15','qcdem2030','qcdem3080','qcdem80170',
        'qcdbcem2030','qcdbcem3080','qcdbcem80170',
        'zzz','wwz','www','ttww','ttg','wwg']:
    #no special treatment needed, put here to avoid typos
    pass
else:
    print "Error: Unknown samplename!"
    exit(8)



####################################################################
## Define which collections (including which corrections) to be used in nTuple

isolatedMuonCollection = "selectedPatMuons"

isolatedElecCollection = "selectedPatElectrons"
#isolatedElecCollection = "selectedPatElectronsAfterScaling"

jetCollection = "hardJets"

jetForMETCollection = "scaledJetEnergy:selectedPatJets"

metCollection = "scaledJetEnergy:patMETs"

genJetCollection = "ak5GenJetsPlusHadron"

genLevelBJetProducerInput = "produceGenLevelBJets"

genHFHadronMatcherInput = "matchGenHFHadronJets"

# Lepton collection used for kinematic reconstruction (has further selections, and can thus deviate from KinReco in analysis)
finalLeptons = 'filterDiLeptonMassQCDveto'



####################################################################
## Separation of ttbar samples in dileptonic and other decays

if topfilter:
    process.load("TopAnalysis.TopFilter.filters.GeneratorTopFilter_cfi")
    process.generatorTopFilter.rejectNonBottomDecaysOfTops = False
    if higgsSignal or ttbarV:
        process.generatorTopFilter.invert_selection = True
        process.generatorTopFilter.channels = ["none"] #empty array would use some defaults
    else:
        all = ['ElectronElectron', 'ElectronElectronViaTau', 
               'MuonMuon', 'MuonMuonViaTau', 
               'ElectronMuon', 'ElectronMuonViaTau']
        if signal:
                process.generatorTopFilter.invert_selection = False
                if viaTau:
                        process.generatorTopFilter.channels = ['ElectronElectronViaTau', 'MuonMuonViaTau', 'ElectronMuonViaTau']
                elif alsoViaTau:
                        process.generatorTopFilter.channels = all
                else:
                        process.generatorTopFilter.channels = ['ElectronElectron', 'ElectronMuon', 'MuonMuon']
        else:
                process.generatorTopFilter.channels = all
                process.generatorTopFilter.invert_selection = True



####################################################################
## Build Jet Collections

process.load("TopAnalysis.TopUtils.JetEnergyScale_cfi")

process.load("TopAnalysis.TopFilter.filters.JetIdFunctorFilter_cfi")
process.goodIdJets.jets    = cms.InputTag("scaledJetEnergy:selectedPatJets")
process.goodIdJets.jetType = cms.string('PF')
process.goodIdJets.version = cms.string('FIRSTDATA')
process.goodIdJets.quality = cms.string('LOOSE')

process.hardJets = selectedPatJets.clone(src = 'goodIdJets', cut = 'pt > 5 & abs(eta) < 2.4') 

# Additional properties for jets like jet charges
process.load("TopAnalysis.HiggsUtils.producers.JetPropertiesProducer_cfi")
process.jetProperties.src = jetCollection

process.buildJets = cms.Sequence(
            process.scaledJetEnergy *
	    #process.selectedPatElectronsAfterScaling *
            process.goodIdJets * 
            process.hardJets *
            process.jetProperties
            )



####################################################################
## Filter on events containing dilepton system of opposite charge and above m(ll) > 12 GeV

##### WARNING: This tool selects lepton pair based on highest pt-sum of the leptons. Since opposite charge is required, it will form the combination of pt-leading lepton and pt-leading antilepton.
##### WARNING: So it might have strange side effects in selections of events with >=3 leptons ?!
##### WARNING: It will for sure have strange side effects in case of additional lepton selections on nTuple level...
##### WARNING: We should replace it by a simple dilepton selector checking for any acceptable combination, not making any specific choice
from TopAnalysis.TopFilter.filters.DiLeptonFilter_cfi import *
process.filterOppositeCharge = filterLeptonPair.clone(
    electrons    = isolatedElecCollection,
    muons        = isolatedMuonCollection,
    Cut          = (0.,0.),
    filterCharge = -1,
)


##### WARNING: This tool uses the lepton pair from the DiLeptonFilter_cfi, based on highest pt-sum of the leptons.
from PhysicsTools.PatAlgos.selectionLayer1.leptonCountFilter_cfi import *
process.filterChannel = countPatLeptons.clone()
process.filterChannel.electronSource    = 'filterOppositeCharge'
process.filterChannel.muonSource        = 'filterOppositeCharge'
process.filterChannel.minNumber         = 2
process.filterChannel.countTaus         = False
if op_mode == 'ee':
    process.filterChannel.countElectrons    = True
    process.filterChannel.countMuons        = False
elif op_mode == 'mumu':
    process.filterChannel.countElectrons    = False
    process.filterChannel.countMuons        = True
elif op_mode == 'emu':
    process.filterChannel.minNumber         = 1
    process.filterChannel1 = process.filterChannel.clone()
    process.filterChannel2 = process.filterChannel1.clone()
    process.filterChannel1.countElectrons    = True
    process.filterChannel1.countMuons        = False
    process.filterChannel2.countElectrons    = False
    process.filterChannel2.countMuons        = True
    process.filterChannel = cms.Sequence(process.filterChannel1 * process.filterChannel2)
else:
    process.filterChannel.countElectrons    = True
    process.filterChannel.countMuons        = True


##### WARNING: This tool uses the lepton pair from the DiLeptonFilter_cfi, based on highest pt-sum of the leptons.
process.filterDiLeptonMassQCDveto           = filterLeptonPair.clone()
process.filterDiLeptonMassQCDveto.muons     = 'filterOppositeCharge'
process.filterDiLeptonMassQCDveto.electrons = 'filterOppositeCharge'
process.filterDiLeptonMassQCDveto.Cut       = (0.,12.)



####################################################################
## Write Ntuple

from TopAnalysis.TopAnalyzer.NTupleWriter_cfi import writeNTuple
writeNTuple.sampleName = op_samplename
writeNTuple.channelName = op_mode
writeNTuple.systematicsName = op_systematicsName
writeNTuple.isMC = op_runOnMC
writeNTuple.isTtBarSample = signal
writeNTuple.isHiggsSample = higgsSignal
writeNTuple.isZSample = zGenInfo
writeNTuple.includePDFWeights = op_includePDFWeights
writeNTuple.pdfWeights = "pdfWeights:cteq66"
writeNTuple.includeZdecay = zproducer
writeNTuple.saveHadronMothers = False

process.writeNTuple = writeNTuple.clone(
    muons = isolatedMuonCollection,
    elecs = isolatedElecCollection,
    jets = jetCollection,
    met = metCollection,
    genMET = "genMetTrue",
    genJets = genJetCollection,

    BHadJetIndex = cms.InputTag(genLevelBJetProducerInput, "BHadJetIndex"),
    AntiBHadJetIndex = cms.InputTag(genLevelBJetProducerInput, "AntiBHadJetIndex"),
    BHadrons = cms.InputTag(genLevelBJetProducerInput, "BHadrons"),
    AntiBHadrons = cms.InputTag(genLevelBJetProducerInput, "AntiBHadrons"),
    BHadronFromTopB = cms.InputTag(genLevelBJetProducerInput, "BHadronFromTopB"),
    AntiBHadronFromTopB = cms.InputTag(genLevelBJetProducerInput, "AntiBHadronFromTopB"),
    BHadronVsJet = cms.InputTag(genLevelBJetProducerInput, "BHadronVsJet"),
    AntiBHadronVsJet = cms.InputTag(genLevelBJetProducerInput, "AntiBHadronVsJet"),
    genBHadPlusMothers = cms.InputTag(genHFHadronMatcherInput,"genBHadPlusMothers"),
    genBHadPlusMothersIndices = cms.InputTag(genHFHadronMatcherInput,"genBHadPlusMothersIndices"),
    genBHadIndex = cms.InputTag(genHFHadronMatcherInput,"genBHadIndex"),
    genBHadFlavour = cms.InputTag(genHFHadronMatcherInput,"genBHadFlavour"),
    genBHadJetIndex = cms.InputTag(genHFHadronMatcherInput,"genBHadJetIndex"),
)
process.writeNTuple.jetsForMET    = cms.InputTag("scaledJetEnergy:selectedPatJets")
process.writeNTuple.jetsForMETuncorr    = cms.InputTag("selectedPatJets")



####################################################################
## Include PDF weights for systematic signal samples

if op_includePDFWeights:
    if not signal:
        print "PDF variations only supported for the signal"
        exit(5615)
    process.pdfWeights = cms.EDProducer("PdfWeightProducer",
                # Fix POWHEG if buggy (this PDF set will also appear on output,
                # so only two more PDF sets can be added in PdfSetNames if not "")
                #FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
                GenTag = cms.untracked.InputTag("genParticles"),
                PdfInfoTag = cms.untracked.InputTag("generator"),
                PdfSetNames = cms.untracked.vstring(
                        "cteq66.LHgrid"
                        #, "MRST2006nnlo.LHgrid"
                        #, "NNPDF10_100.LHgrid"
                        #"cteq6mE.LHgrid"
                        # ,"cteq6m.LHpdf"
                        #"cteq6m.LHpdf"
                ))
else:
    process.pdfWeights = cms.Sequence()



####################################################################
## Kinematic reconstruction

# std sequence to produce the ttFullLepEvent
process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttFullLepEvtBuilder_cff")
from TopQuarkAnalysis.TopEventProducers.sequences.ttFullLepEvtBuilder_cff import *
if not signal:
    removeTtFullLepHypGenMatch(process)

setForAllTtFullLepHypotheses(process,"muons"    ,finalLeptons)
setForAllTtFullLepHypotheses(process,"electrons",finalLeptons)
# WARNING! The jet.pt > 30 cut is currently hardcoded in the NTupleWriter.cc file
# adding a collections like
#     process.jetsForKinReco = process.hardJets.clone(src = 'hardJets', cut = 'pt > 30')
# will cause problems because the selection of the "best" solution is hardcoded!!!!!
#setForAllTtFullLepHypotheses(process,"jets"     ,'jetsForKinReco')
setForAllTtFullLepHypotheses(process,"jets"     ,jetCollection)
setForAllTtFullLepHypotheses(process,"mets"     ,metCollection)
if op_runOnMC:
    setForAllTtFullLepHypotheses(process,"jetCorrectionLevel","L3Absolute")
    print "L3Absolute"
else:
    setForAllTtFullLepHypotheses(process,"jetCorrectionLevel","L2L3Residual")
    print "L2L3Residual"
setForAllTtFullLepHypotheses(process,"maxNJets",-1)

process.kinSolutionTtFullLepEventHypothesis.maxNComb = -1
process.kinSolutionTtFullLepEventHypothesis.searchWrongCharge = True
process.kinSolutionTtFullLepEventHypothesis.tmassbegin = 100.0
process.kinSolutionTtFullLepEventHypothesis.tmassend   = 300.0
process.kinSolutionTtFullLepEventHypothesis.neutrino_parameters = (30.641, 57.941, 22.344, 57.533, 22.232)
#according to our MC 8 TeV values are:
#nu    mpv 40.567 sigma = 16.876
#nubar mpv 40.639 sigma = 17.021

#process.kinSolutionTtFullLepEventHypothesis.mumuChannel = False
#process.kinSolutionTtFullLepEventHypothesis.eeChannel = False
#process.kinSolutionTtFullLepEventHypothesis.emuChannel = True
process.ttFullLepEvent.decayChannel1 = cms.int32(1)
process.ttFullLepEvent.decayChannel2 = cms.int32(2)

process.kinSolutionTtFullLepEventHypothesis.mumuChannel = True
process.kinSolutionTtFullLepEventHypothesis.emuChannel  = True
process.kinSolutionTtFullLepEventHypothesis.eeChannel = True

#process.ttFullLepEvent.decayChannel1 = cms.int32(1)
#process.ttFullLepEvent.decayChannel2 = cms.int32(2)


#if op_mode == 'mumu':
    #process.kinSolutionTtFullLepEventHypothesis.mumuChannel = True
    #process.ttFullLepEvent.decayChannel1 = cms.int32(2)
    #process.ttFullLepEvent.decayChannel2 = cms.int32(2)
#elif op_mode == 'emu':
    #process.kinSolutionTtFullLepEventHypothesis.emuChannel = True
    #process.ttFullLepEvent.decayChannel1 = cms.int32(1)
    #process.ttFullLepEvent.decayChannel2 = cms.int32(2)
#elif op_mode == 'ee':
    #process.kinSolutionTtFullLepEventHypothesis.eeChannel = True
    #process.ttFullLepEvent.decayChannel1 = cms.int32(1)
    #process.ttFullLepEvent.decayChannel2 = cms.int32(1)



####################################################################
## Sample-specific sequences

if zproducer:
    process.load("TopAnalysis.TopUtils.ZDecayProducer_cfi")
    process.zsequence = cms.Sequence(process.ZDecayProducer)
else:
    process.zsequence = cms.Sequence()

if zGenInfo:
    process.load("TopAnalysis.HiggsUtils.producers.GenZDecay_cfi")
    process.zGenSequence = cms.Sequence(process.genZDecay)
else:
    process.zGenSequence = cms.Sequence()

if topfilter:
    process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
    #process.load("TopAnalysis.TopUtils.HadronLevelBJetProducer_cfi")

    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi") # supplies PDG ID to real name resolution of MC particles, necessary for GenLevelBJetProducer
    process.load("TopAnalysis.TopUtils.GenLevelBJetProducer_cfi")
    process.produceGenLevelBJets.deltaR = 5.0
    process.produceGenLevelBJets.noBBbarResonances = True

    process.load("TopAnalysis.TopUtils.GenHFHadronMatcher_cff")
    process.matchGenHFHadronJets.flavour = 5
    process.matchGenHFHadronJets.noBBbarResonances = True

    process.load("TopAnalysis.TopUtils.sequences.improvedJetHadronQuarkMatching_cff")

    process.decaySubset.fillMode = "kME" # Status3, use kStable for Status2
    if signal:
        process.topsequence = cms.Sequence(
            process.makeGenEvt *
            process.improvedJetHadronQuarkMatchingSequence *
            process.generatorTopFilter *
            process.produceGenLevelBJets *
            process.matchGenHFHadronJets)
    else:
        process.topsequence = cms.Sequence(
            process.makeGenEvt *
            process.generatorTopFilter)

else:
    process.topsequence = cms.Sequence()


if higgsSignal:
    process.load("TopAnalysis.HiggsUtils.filters.GeneratorHiggsFilter_cfi")
    process.generatorHiggsFilter.channels = ["none"]
    process.generatorHiggsFilter.invert_selection = True
    process.load("TopAnalysis.HiggsUtils.sequences.higgsGenEvent_cff")
    process.decaySubsetHiggs.fillMode = "kME" # Status3, use kStable for Status2
    process.higgssequence = cms.Sequence(
        process.makeGenEvtHiggs *
        process.generatorHiggsFilter
    )
else:
    process.higgssequence = cms.Sequence()

if signal or higgsSignal or zGenInfo:
    process.ntupleInRecoSeq = cms.Sequence()
else:
    process.ntupleInRecoSeq = cms.Sequence(process.zsequence * process.writeNTuple) #doesn't affect writing to output file
    


####################################################################
## Remove all the tau stuff

# from PhysicsTools.PatAlgos.tools.coreTools import removeSpecificPATObjects
# removeSpecificPATObjects( process
#                         , names = ['Taus', 'Photons']
#                         , outputModules = []
#                         , postfix = pfpostfix
#                         )
# # Remove the full pftau sequence as it is not needed for us
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTauPFJets08Region'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTauPileUpVertices'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTauTagInfoProducer'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfJetsPiZeros'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfJetsLegacyTaNCPiZeros'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfJetsLegacyHPSPiZeros'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTausBase'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsSelectionDiscriminator'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauProducerSansRefs'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauProducer'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTausBaseDiscriminationByDecayModeFinding'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTausBaseDiscriminationByLooseCombinedIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfTaus'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfNoTau'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByDecayModeFinding'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByVLooseChargedIsolation'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseChargedIsolation'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumChargedIsolation'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightChargedIsolation'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByVLooseIsolation'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseIsolation'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumIsolation'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightIsolation'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByRawChargedIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByRawGammaIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseElectronRejection'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumElectronRejection'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightElectronRejection'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMVAElectronRejection'+pfpostfix))
#getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByLooseMuonRejection'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByMediumMuonRejection'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'hpsPFTauDiscriminationByTightMuonRejection'+pfpostfix))

# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'tauIsoDepositPFCandidates'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'tauIsoDepositPFChargedHadrons'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'tauIsoDepositPFNeutralHadrons'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'tauIsoDepositPFGammas'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'patTaus'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'selectedPatTaus'+pfpostfix))
# getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatTaus'+pfpostfix))

## removal of unnecessary modules
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'patPFParticles'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'patCandidateSummary'+pfpostfix))

getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'selectedPatPFParticles'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'selectedPatCandidateSummary'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatElectrons'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatMuons'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatLeptons'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatJets'+pfpostfix))
getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'countPatPFParticles'+pfpostfix))

getattr(process,'patPF2PATSequence'+pfpostfix).remove(getattr(process,'pfPhotonSequence'+pfpostfix))

massSearchReplaceAnyInputTag(getattr(process,'patPF2PATSequence'+pfpostfix),'pfNoTau'+pfpostfix,'pfJets'+pfpostfix)
##############LETS MAKE SOME BEEEEAAAANNNSSSSSS#####################

process.load("CMGTools.External.pujetidsequence_cff")

# OK, out of laziness, I'm only going to configure the pu jet ID for our PF2PAT jets
#process.puJetIdChs.vertexes = 'goodOfflinePrimaryVertices'
#process.puJetIdChs.jets = 'selectedPatJets'

#process.puJetMvaChs.vertexes = 'goodOfflinePrimaryVertices'
#process.puJetMvaChs.jets = 'selectedPatJets'



from RecoMuon.MuonIsolationProducers.caloExtractorByAssociatorBlocks_cff import *
from RecoMuon.MuonIsolationProducers.trackExtractorBlocks_cff import *
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.BNproducer = cms.EDProducer('BEANmaker',
                                    pfmetTag = cms.InputTag("patMETs"),
                                    pfmetTag_type1correctedRECO = cms.InputTag("pfType1CorrectedMet"),
                                    pfmetTag_uncorrectedPF2PAT  = cms.InputTag("patPFMet"),
                                    pfmetTag_uncorrectedRECO    = cms.InputTag("pfMET"),
                                    eleTag = cms.InputTag("selectedPatElectrons"),
                                    genParticleTag = cms.InputTag("genParticles"),
                                    pfjetTag = cms.InputTag("selectedPatJets"),
                                    genjetTag = cms.InputTag("ak5GenJets"),
                                    muonTag = cms.InputTag("selectedPatMuons"),
                                    photonTag = cms.InputTag("none"),
                                    EBsuperclusterTag = cms.InputTag("correctedHybridSuperClusters"),
                                    EEsuperclusterTag = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
                                    reducedBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                                    reducedEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                                    trackTag = cms.InputTag("generalTracks"),
                                    tauTag = cms.InputTag("selectedPatTaus"),
                                    triggerResultsTag = cms.InputTag("TriggerResults::HLT"),
                                    gtSource = cms.InputTag("gtDigis"),
                                    pvTag = cms.InputTag("offlinePrimaryVertices"),
                                    triggerSummaryTag = cms.InputTag("hltTriggerSummaryAOD"),
                                    dcsTag = cms.InputTag("scalersRawToDigi"),
                                    hltProcessName = cms.string("HLT"),
                                    eventWeight = cms.double(0.0126976903794279),
                                    minSCEt = cms.double(10),
                                    minPhotonEt = cms.double(10),
                                    minJetPt = cms.double(10),
                                    minTrackPt = cms.double(10),
                                    verbose = cms.bool(True),
                                    fillTrackHitInfo = cms.bool(False),
                                    fillTrackIsoInfo = cms.bool(False),
                                    CaloExtractorPSet  = cms.PSet(MIsoCaloExtractorByAssociatorTowersBlock),
                                    TrackExtractorPSet = cms.PSet(MIsoTrackExtractorBlock),
                                    sample = cms.int32(120),
                                    maxAbsZ = cms.untracked.double(24)
                                    )




process.q2weights = cms.EDProducer('Q2Weights')



####################################################################
## Paths, one with preselection, one without for signal samples

process.metseq = cms.Sequence(
    process.pfJetMETcorr *
    process.pfType1CorrectedMet
    )

process.p = cms.Path(
    #process.goodOfflinePrimaryVertices *
    #getattr(process,'patPF2PATSequence'+pfpostfix) *
    #process.buildJets                     *
    #process.filterOppositeCharge          *
    #process.filterChannel                 *
    #     process.filterDiLeptonMassQCDveto     *
    #     process.makeTtFullLepEvent            *
    #process.ntupleInRecoSeq               

 )

if signal or higgsSignal or zGenInfo:
    process.pNtuple = cms.Path(
        process.goodOfflinePrimaryVertices *
        process.q2weights *
        getattr(process,'patPF2PATSequence'+pfpostfix) *
        process.metseq *
        process.recoTauClassicHPSSequence *
        process.BNproducer
        )

#process.out.outputCommands = ['drop *']
process.out.outputCommands.extend( [
    'drop *',
    'keep *_BNproducer_*_*'
    ])


####################################################################
## Prepend PF2PAT

from TopAnalysis.TopAnalyzer.CountEventAnalyzer_cfi import countEvents
process.EventsBeforeSelection = countEvents.clone()
process.EventsBeforeSelection.includePDFWeights = op_includePDFWeights
process.EventsBeforeSelection.pdfWeights = "pdfWeights:cteq66"
    

# pathnames = process.paths_().keys()
# print 'prepending trigger sequence to paths:', pathnames
# for pathname in pathnames:
#     getattr(process, pathname).insert(0, cms.Sequence(
#         process.pdfWeights *
#         process.EventsBeforeSelection * 
#         process.topsequence *
#         process.higgssequence *
#         process.zGenSequence *
#         process.filterTrigger
#         ))
# if signal or higgsSignal or zGenInfo:
#     process.pNtuple.remove(process.filterTrigger)

process.scaledJetEnergy.inputElectrons       = "selectedPatElectrons"
process.scaledJetEnergy.inputJets            = "selectedPatJets"
process.scaledJetEnergy.inputMETs            = "patMETs"
process.scaledJetEnergy.JECUncSrcFile        = cms.FileInPath("TopAnalysis/TopUtils/data/Summer13_V4_DATA_UncertaintySources_AK5PFchs.txt")
process.scaledJetEnergy.scaleType = "abs"   #abs = 1, jes:up, jes:down

if op_runOnMC:
    process.scaledJetEnergy.resolutionEtaRanges  = cms.vdouble(0, 0.5, 0.5, 1.1, 1.1, 1.7, 1.7, 2.3, 2.3, 5.4)
    process.scaledJetEnergy.resolutionFactors    = cms.vdouble(1.052, 1.057, 1.096, 1.134, 1.288) # JER standard

    #please change this on the top where the defaults for the VarParsing are given
    if op_systematicsName == "JES_UP":
        process.scaledJetEnergy.scaleType = "jes:up"
    if op_systematicsName == "JES_DOWN":
        process.scaledJetEnergy.scaleType = "jes:down"
    if op_systematicsName == "JER_UP":
        process.scaledJetEnergy.resolutionFactors = cms.vdouble(1.115, 1.114, 1.161, 1.228, 1.488)
    if op_systematicsName == "JER_DOWN":
        process.scaledJetEnergy.resolutionFactors = cms.vdouble(0.990, 1.001, 1.032, 1.042, 1.089)
else:
    process.scaledJetEnergy.resolutionEtaRanges  = cms.vdouble(0, -1)
    process.scaledJetEnergy.resolutionFactors    = cms.vdouble(1.0) # JER standard
#     for pathname in pathnames:
#         getattr(process, pathname).replace(process.goodOfflinePrimaryVertices,
#                                            process.HBHENoiseFilter * 
#                                            process.scrapingFilter * 
#                                            process.ecalLaserCorrFilter * 
#                                            process.goodOfflinePrimaryVertices)
        


# process.load("Configuration.EventContent.EventContent_cff")
# process.out = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string('ttbarZ.root'),
#     SelectEvents = cms.untracked.PSet(
#     SelectEvents = cms.vstring('p','pNtuple','pNtuple')),
#     dropMetaData = cms.untracked.string('ALL'),
#     outputCommands = cms.untracked.vstring('drop *',
#                                'keep *_BNproducer_*_*',
#                                'keep double_kt6PFJets*_rho_*',
#                                'keep *')
#                                )

####################################################################
## Particle tree drawer

# see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCandidateModules#ParticleTreeDrawer_Utility
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
#                                   src = cms.InputTag("genParticles"),                                                                 
#             #                      printP4 = cms.untracked.bool(False),
#             #                      printPtEtaPhi = cms.untracked.bool(False),
#             #                      printVertex = cms.untracked.bool(False),
#             #                      printStatus = cms.untracked.bool(False),
#             #                      printIndex = cms.untracked.bool(False),
#             #                      status = cms.untracked.vint32( 3 )
#                                   )
#process.p = cms.Path(process.printTree)
#process.pNtuple = cms.Path()



####################################################################
## Signal catcher for more information on errors

process.load("TopAnalysis.TopUtils.SignalCatcher_cfi")

#Dump python config if wished
outfile = open('dumpedConfig.py','w'); print >> outfile,process.dumpPython(); outfile.close()

