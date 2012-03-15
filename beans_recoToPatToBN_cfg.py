
#
# This file contains the NOVA group reference selection work-flow for ele and mu + jets analysis.
# See https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopLeptonPlusJetsRefSel_mu#Selection_Version_SelV4_valid_fr
#

import sys

import FWCore.ParameterSet.Config as cms

process = cms.Process( 'BEAN' )


### ======================================================================== ###
###                                                                          ###
###                                 Constants                                ###
###                            (user job steering)                           ###
###                                                                          ###
### ======================================================================== ###


### Data or MC?
runOnMC = False

### Standard and PF work flow

# Standard
runStandardPAT = True
usePFJets      = True
useCaloJets    = True

# PF2PAT
runPF2PAT = True

### Switch on/off selection steps

# Whether you want to include trigger skimming in this selection
useTrigger  = False

# Whether you want to include muon skimming in this selection
useSkimMuon = False


### Reference selection

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets import *
muonsUsePV         = True
electronsUsePV     = True

photonCut          = 'et > 10. && abs(eta) < 3.5'

#Re-jigger the PF muon and PF electron isolation cuts...
# PF muon: used for loose veto
muonCutBase  =     'pt > 10.'                                                    # transverse momentum
muonCutBase += ' && abs(eta) < 2.5'                                              # pseudo-rapisity range

muonCut    =  muonCutBase
muonCutPF  =  muonCutBase
muonCutPF += ' && (chargedHadronIso+neutralHadronIso+photonIso)/pt < 0.20'     # relative isolation

muonCutLoosePF  =  muonCutBase

# PF electron
electronCutBase  =     'et > 15.'                                                  # transverse energy
electronCutBase += ' && abs(eta) < 2.5'                                            # pseudo-rapisity range

electronCut    = electronCutBase
electronCutPF  = electronCutBase
electronCutPF += ' && (chargedHadronIso+neutralHadronIso+photonIso)/et < 0.20' # relative isolation

electronCutLoosePF  = electronCutBase


# Trigger selection according to run range resp. MC sample:
# lower range limits for data available as suffix;
# available are: 000000, 147196, 160404, 163270 (default)
# sample name for MC available as suffix;
# available are: Summer11 (default)
#triggerSelectionData       = triggerSelection_163270
#triggerObjectSelectionData = triggerObjectSelection_163270
#triggerSelectionMC       = triggerSelection_Summer11
#triggerObjectSelectionMC = triggerObjectSelection_Summer11
triggerSelectionMC   = 'HLT_IsoMu17_v*'
triggerSelectionData = 'HLT_IsoMu17_v*'


### Particle flow
### takes effect only, if 'runPF2PAT' = True

postfix = 'PFlow' # needs to be a non-empty string, if 'runStandardPAT' = True

# subtract charged hadronic pile-up particles (from wrong PVs)
# effects also JECs
usePFnoPU = True # before any top projection

# other switches for PF top projections (default: all 'True')
useNoMuon     = True # before electron top projection
useNoElectron = True # before jet top projection
useNoJet      = True # before tau top projection
useNoTau      = False # before MET top projection

# cuts used in top projections
from TopQuarkAnalysis.Configuration.patRefSel_PF2PAT import *
# vertices
#pfD0Cut   = 0.2
#pfDzCut   = 0.5
#pfVertices = 'goodOfflinePrimaryVertices'
# muons
pfMuonSelectionCut =  'pt > 10. && abs(eta) < 2.5'
pfMuonIsoConeR   = 0.4
pfMuonCombIsoCut = cms.double(0.2) #Loose lepton isolation
# electrons
pfElectronSelectionCut  = 'et > 15. && abs(eta) < 2.5'
pfElectronIsoConeR   = 0.4
pfElectronCombIsoCut = cms.double(0.2)

# What to use for default pfIso cone (default 0.4)
pfMuonIsoConeR03 = False
pfElectronnIsoConeR03 = False

### JEC levels

# levels to be accessible from the jets
# jets are corrected to L3Absolute (MC), L2L3Residual (data) automatically, if enabled here
# and remain uncorrected, if none of these levels is enabled here
useL1FastJet    = True  # needs useL1Offset being off, error otherwise
useL1Offset     = False # needs useL1FastJet being off, error otherwise
useL2Relative   = True
useL3Absolute   = True
useL2L3Residual = True  
useL5Flavor     = False
useL7Parton     = False

useType1Met = True

### Input

# list of input files
useRelVals = False # if 'False', "inputFiles" is used
inputFiles = [
  #  'file:test.root', #ttH120 sample
  #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/FEE3E76C-1EFA-E011-BAF4-002618943864.root',
  #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v2/0000/0040932D-A80F-E111-BBF7-00304867BFAA.root',
  #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v2/0000/00583179-970F-E111-9C30-003048678BB8.root',
  #'/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v2/0000/00859C66-AD0F-E111-A82D-0018F3D096DE.root',

  '/store/data/Run2011B/SingleMu/RECO/PromptReco-v1/000/175/832/10EBFA08-84DB-E011-9D72-003048D2C020.root',
  ] # overwritten, if "useRelVals" is 'True'

# maximum number of events
maxInputEvents = 100 # reduce for testing

### Conditions

# GlobalTags (w/o suffix '::All')
globalTagData = 'GR_R_42_V24' 
globalTagMC   = 'START42_V17' 

### Output

# output file
#outputFile = 'TTJets_TuneZ2_7TeV_madgraph_tauola_Fall11_PU_S6_START42_V14B_v2_BEAN_V05.root'
outputFile = 'SingleMu_2011B_PromptReco_v1_BEAN_V05.root'

# event frequency of Fwk report
fwkReportEvery = 10000 # reduce for testing

# switch for 'TrigReport'/'TimeReport' at job end
wantSummary = True


### ======================================================================== ###
###                                                                          ###
###                              End of constants                            ###
###                                                                          ###
### ======================================================================== ###



###
### Basic configuration
###

keepSC = False
keepTK = False

process.load( "TopQuarkAnalysis.Configuration.patRefSel_basics_cff" )
process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery
process.options.wantSummary = wantSummary
if runOnMC:
  process.GlobalTag.globaltag = globalTagMC   + '::All'
else:
  process.GlobalTag.globaltag = globalTagData + '::All'


###
### Input configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_inputModule_cfi" )
if useRelVals:
  from PhysicsTools.PatAlgos.tools.cmsswVersionTools import pickRelValInputFiles
  if runOnMC:
    inputFiles = pickRelValInputFiles( cmsswVersion  = 'CMSSW_4_2_3'
                                     , relVal        = 'RelValTTbar'
                                     , globalTag     = globalTagMC
                                     , numberOfFiles = -1 # "-1" means "all"
                                     )
  else:
    inputFiles = pickRelValInputFiles( cmsswVersion  = 'CMSSW_4_2_3'
                                     , relVal        = 'Mu'
                                     , dataTier      = 'RECO'
                                     #, globalTag     = globalTagData + '_RelVal_mu2010B'
                                     , globalTag     = globalTagData + '_mu2010B' # wrong naming scheme in CMSSW_4_2_3
                                     , numberOfFiles = -1 # "-1" means "all"
                                     )
process.source.fileNames = inputFiles
process.maxEvents.input  = maxInputEvents


###
### Output configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_outputModule_cff" )
# output file name
process.out.fileName = outputFile
# event content
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out.outputCommands += patEventContent
# clear event selection
process.out.SelectEvents.SelectEvents = []


###
### Cleaning and trigger selection configuration
###

### Event cleaning
process.load( 'TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff' )

### Trigger selection
if runOnMC:
  triggerSelection = triggerSelectionMC
else:
  triggerSelection = triggerSelectionData
from TopQuarkAnalysis.Configuration.patRefSel_triggerSelection_cff import triggerResults
process.step0a = triggerResults.clone(
  triggerConditions = [ triggerSelection ]
)


### Good vertex selection
process.load( "TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi" )


###
### PAT/PF2PAT configuration
###

pfSuffix = 'PF'
caloSuffix = 'Calo'
if runStandardPAT and runPF2PAT:
  if postfix == '':
    sys.exit( 'ERROR: running standard PAT and PF2PAT in parallel requires a defined "postfix" for PF2PAT' )
  if usePFJets:
    if postfix == 'Add' + pfSuffix or postfix == jetAlgo + pfSuffix:
      sys.exit( 'ERROR: running standard PAT with additional PF jets  and PF2PAT in parallel does not allow for the "postfix" %s'%( postfix ) )
if not runStandardPAT and not runPF2PAT:
  sys.exit( 'ERROR: standard PAT and PF2PAT are both switched off' )

process.load( "PhysicsTools.PatAlgos.patSequences_cff" )
from PhysicsTools.PatAlgos.tools.coreTools import *

# Add PAT trigger information to the configuration
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, hltProcess = '*' )


### Check JECs

# JEC set
jecSet        = jecSetBase + 'Calo'
jecSetAddCalo = jecSetBase + caloSuffix
jecSetAddPF   = jecSetBase + pfSuffix
jecSetPF      = jecSetAddPF
if usePFnoPU:
  jecSetPF += 'chs'

# JEC levels
if useL1FastJet and useL1Offset:
  sys.exit( 'ERROR: switch off either "L1FastJet" or "L1Offset"' )
jecLevels = []
if useL1FastJet:
  jecLevels.append( 'L1FastJet' )
if useL1Offset:
  jecLevels.append( 'L1Offset' )
if useL2Relative:
  jecLevels.append( 'L2Relative' )
if useL3Absolute:
  jecLevels.append( 'L3Absolute' )
if useL2L3Residual and not runOnMC:
  jecLevels.append( 'L2L3Residual' )
if useL5Flavor:
  jecLevels.append( 'L5Flavor' )
if useL7Parton:
  jecLevels.append( 'L7Parton' )

### Switch configuration

if runPF2PAT:
  from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
  usePF2PAT( process
           , runPF2PAT      = runPF2PAT
           , runOnMC        = runOnMC
           , jetAlgo        = jetAlgo
           , postfix        = postfix
           , jetCorrections = ( jecSetPF
                              , jecLevels
                              )
           )
  applyPostfix( process, 'pfNoPileUp'  , postfix ).enable = usePFnoPU
  applyPostfix( process, 'pfNoMuon'    , postfix ).enable = useNoMuon
  applyPostfix( process, 'pfNoElectron', postfix ).enable = useNoElectron
  applyPostfix( process, 'pfNoJet'     , postfix ).enable = useNoJet
  applyPostfix( process, 'pfNoTau'     , postfix ).enable = useNoTau
  if useL1FastJet:
    applyPostfix( process, 'pfPileUp', postfix ).Vertices            = cms.InputTag( pfVertices )
    applyPostfix( process, 'pfPileUp', postfix ).checkClosestZVertex = False
    applyPostfix( process, 'pfJets', postfix ).doAreaFastjet = True
    applyPostfix( process, 'pfJets', postfix ).doRhoFastjet  = False
  applyPostfix( process, 'pfMuonsFromVertex'    , postfix ).vertices = cms.InputTag( pfVertices )
  applyPostfix( process, 'pfMuonsFromVertex'    , postfix ).d0Cut    = pfD0Cut
  applyPostfix( process, 'pfMuonsFromVertex'    , postfix ).dzCut    = pfDzCut
  applyPostfix( process, 'pfSelectedMuons'      , postfix ).cut = pfMuonSelectionCut
  getattr(process, 'pfIsolatedMuons'+postfix ).isolationCut = pfMuonCombIsoCut
  getattr(process, 'pfIsolatedMuons'+postfix ).combinedIsolationCut = pfMuonCombIsoCut
  if pfMuonIsoConeR03:
    applyPostfix( process, 'pfIsolatedMuons', postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03' + postfix ) )
    applyPostfix( process, 'pfIsolatedMuons', postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03' + postfix )
    applyPostfix( process, 'pfIsolatedMuons', postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
                                                                                                  , cms.InputTag( 'muPFIsoValueGamma03' + postfix ) )
    applyPostfix( process, 'patMuons', postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
    applyPostfix( process, 'patMuons', postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'muPFIsoValuePU03' + postfix )
    applyPostfix( process, 'patMuons', postfix ).isolationValues.pfPhotons          = cms.InputTag( 'muPFIsoValueGamma03' + postfix )
    applyPostfix( process, 'patMuons', postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'muPFIsoValueCharged03' + postfix )
  applyPostfix( process, 'pfElectronsFromVertex'    , postfix ).vertices = cms.InputTag( pfVertices )
  applyPostfix( process, 'pfElectronsFromVertex'    , postfix ).d0Cut    = pfD0Cut
  applyPostfix( process, 'pfElectronsFromVertex'    , postfix ).dzCut    = pfDzCut
  applyPostfix( process, 'pfSelectedElectrons'    , postfix ).cut = pfElectronSelectionCut
  getattr(process, 'pfIsolatedElectrons'+postfix ).isolationCut = pfElectronCombIsoCut
  getattr(process, 'pfIsolatedElectrons'+postfix ).combinedIsolationCut = pfElectronCombIsoCut
  if pfElectronIsoConeR03:
    applyPostfix( process, 'isoValElectronWithCharged', postfix ).deposits[0].deltaR = 0.3
    applyPostfix( process, 'isoValElectronWithNeutral', postfix ).deposits[0].deltaR = 0.3
    applyPostfix( process, 'isoValElectronWithPhotons', postfix ).deposits[0].deltaR = 0.3

# remove MC matching, object cleaning, objects etc.
if runStandardPAT:
  if not runOnMC:
    runOnData( process )
  if useCaloJets:
    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    addJetCollection( process
                    , cms.InputTag( jetAlgo.lower() + caloSuffix + 'Jets' )
                    , jetAlgo
                    , caloSuffix
                    , doJTA            = True
                    , doBTagging       = True
                    , jetCorrLabel     = ( jecSetAddCalo, jecLevels )
                    , doType1MET       = False
                    , doL1Cleaning     = False
                    , doL1Counters     = True
                    , genJetCollection = cms.InputTag( jetAlgo.lower() + 'GenJets' )
                    , doJetID          = True
                    )
  if usePFJets:
    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    from PhysicsTools.PatAlgos.tools.metTools import addPfMET
    from PhysicsTools.PatAlgos.tools.metTools import addTcMET
    inputTag   = cms.InputTag( jetAlgo.lower() + pfSuffix + 'Jets' )
    inputTagMC = cms.InputTag( jetAlgo.lower() + 'GenJets' )
    addJetCollection( process
                    , cms.InputTag( jetAlgo.lower() + pfSuffix + 'Jets' )
                    , jetAlgo
                    , pfSuffix
                    , doJTA            = True
                    , doBTagging       = True
                    , jetCorrLabel     = ( jecSetAddPF, jecLevels )
                    , doType1MET       = False
                    , doL1Cleaning     = False
                    , doL1Counters     = True
                    , genJetCollection = cms.InputTag( jetAlgo.lower() + 'GenJets' )
                    , doJetID          = True
                    )
    addPfMET( process
            , jetAlgo + pfSuffix
            )
    addTcMET(process, jetAlgo+'TC')
    process.patMETsAK5Calo = process.patMETs
  removeSpecificPATObjects( process
                          , names = [ 'Taus' ]
                          #                          , names = [ 'Photons', 'Taus' ]
                          ) # includes 'removeCleaning'
if runPF2PAT:
  if not runOnMC:
    runOnData( process
             , names = [ 'PFAll' ]
             , postfix = postfix
             )
  removeSpecificPATObjects( process
                          , names = [ 'Photons', 'Taus' ]
                          , postfix = postfix
                          ) # includes 'removeCleaning'


# JetCorrFactorsProducer configuration has to be fixed in standard work flow after a call to 'runOnData()':
if runStandardPAT:
  process.patJetCorrFactors.payload = jecSet
  process.patJetCorrFactors.levels  = jecLevels
# additional event content has to be (re-)added _after_ the call to 'removeCleaning()':
process.out.outputCommands += [ 'keep edmTriggerResults_*_*_*'
                              , 'keep *_hltTriggerSummaryAOD_*_*'
                              , 'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*'
                              # vertices and beam spot
                              , 'keep *_offlineBeamSpot_*_*'
                              , 'keep *_offlinePrimaryVertices*_*_*'
                              , 'keep *_goodOfflinePrimaryVertices*_*_*'
                              , 'keep patTriggerObjects_patTrigger_*_*'
                              , 'keep patTriggerFilters_patTrigger_*_*'
                              , 'keep patTriggerPaths_patTrigger_*_*'
                              , 'keep patTriggerEvent_patTriggerEvent_*_*'
                              , 'keep *_l1extraParticles_*_*'
                                ]

if keepSC:
  process.out.outputCommands += [ 'keep *_correctedHybridSuperClusters_*_*'
                                , 'keep *_correctedMulti5x5SuperClustersWithPreshower_*_*'
                                ]
if keepTK:
  process.highPtTracks = cms.EDFilter (
    "TrackSelector",
    src = cms.InputTag ('generalTracks'),
    cut = cms.string ('pt > 10')
  )
#  process.out.outputCommands += [ 'keep *_generalTracks_*_*'
  process.out.outputCommands += [ 'keep *_highPtTracks_*_*'
                                ]
if runOnMC:
  process.out.outputCommands += [ 'keep GenEventInfoProduct_*_*_*'
                                , 'keep recoGenParticles_*_*_*'
                                , 'keep *_addPileupInfo_*_*'
                                , 'keep LHEEventProduct_*_*_*'
                                ]


###
### Additional configuration
###

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import *

if runStandardPAT:
  ### Jets

  if useL1FastJet:
    process.kt6PFJets = kt6PFJets.clone( src          = cms.InputTag( 'particleFlow' )
                                         , doRhoFastjet = True
                                         )
    
    # compute FastJet rho to correct isolation
    process.kt6PFJetsForIsolation = kt6PFJets.clone()
    process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

    process.load('RecoJets.Configuration.RecoPFJets_cff')
    process.ak5PFJets.doAreaFastjet = True
    setattr(process,'ak5PFJets',ak5PFJets)
    process.patDefaultSequence.replace( process.patJetCorrFactors
                                      , process.kt6PFJets * process.kt6PFJetsForIsolation * process.ak5PFJets * process.patJetCorrFactors
                                      )

    process.out.outputCommands.append( 'keep double_*_*_' + process.name_() )

if runPF2PAT:

  ### Jets

  applyPostfix( process, 'patJetCorrFactors', '' ).useRho = cms.bool(True)
  applyPostfix( process, 'patJetCorrFactors', postfix ).primaryVertices = cms.InputTag( pfVertices )
  if useL1FastJet:
    if usePFnoPU:
      kt6PFJetsPFChs = kt6PFJetsChs.clone( src = cms.InputTag( 'pfNoElectron' + postfix ) )
      #Temporarily set FastJet parameters here.  Need to migrate into top group default
      kt6PFJetsPFChs.voronoiRfact = -0.9  #Should be default 
      setattr( process, 'kt6PFJetsChs' + postfix, kt6PFJetsPFChs )
      
      # compute FastJet rho to correct isolation
      kt6PFJetsPFChsForIsolation = kt6PFJetsPFChs.clone()
      kt6PFJetsPFChsForIsolation.Rho_EtaMax = cms.double(2.5)
      setattr( process, 'kt6PFJetsChsForIsolation' + postfix, kt6PFJetsPFChsForIsolation )
      
      getattr( process, 'patPF2PATSequence' + postfix).replace( getattr( process, 'patJetCorrFactors' + postfix )
                                                              , getattr( process, 'kt6PFJetsChs' + postfix )
                                                                * getattr( process, 'kt6PFJetsChsForIsolation' + postfix )
                                                                * getattr( process, 'patJetCorrFactors' + postfix )
                                                              )
      applyPostfix( process, 'patJetCorrFactors', postfix ).rho = cms.InputTag( 'kt6PFJetsChs' + postfix, 'rho' )
    else:
      kt6PFJetsPF = kt6PFJets.clone( doRhoFastjet = True )
      setattr( process, 'kt6PFJets' + postfix, kt6PFJetsPF )
      getattr( process, 'patPF2PATSequence' + postfix).replace( getattr( process, 'patJetCorrFactors' + postfix )
                                                              , getattr( process, 'kt6PFJets' + postfix ) * getattr( process, 'patJetCorrFactors' + postfix )
                                                              )
      applyPostfix( process, 'patJetCorrFactors', postfix ).rho = cms.InputTag( 'kt6PFJets' + postfix, 'rho' )

    if useType1Met:
      # addPfMET TypeI
      pfMEtCorrector = None
      if runOnMC:
        pfMEtCorrector = "ak5PFL2L3"
        #pfMEtCorrector = "ak5PFL1FastL2L3"  ## crashes
      else:
        pfMEtCorrector = "ak5PFL2L3Residual"
        #pfMEtCorrector = "ak5PFL1FastL2L3Residual"  ## crashes
      from JetMETCorrections.Type1MET.MetType1Corrections_cff import metJESCorAK5PFJet 
      metJESCorAK5PFTypeI = metJESCorAK5PFJet.clone( 
                                                    inputUncorJetsLabel = "patJetsPF", 
                                                    metType = "pat",                  
                                                    inputUncorMetLabel = "pfMet",
      )
      metJESCorAK5PFTypeI.inputUncorJetsLabel = cms.string( 'patJets' + postfix )
      metJESCorAK5PFTypeI.jetPTthreshold = cms.double(10.0)
      metJESCorAK5PFTypeI.corrector = cms.string(pfMEtCorrector)
      setattr(process,'metJESCorAK5PFTypeI'+postfix,metJESCorAK5PFTypeI)
      patMETsTypeI = applyPostfix( process, 'patMETs', postfix ).clone( metSource = cms.InputTag( 'metJESCorAK5PFTypeI' + postfix ) )
      setattr(process,'patMETsTypeI'+postfix,patMETsTypeI)
      # Add to producersLayer1 sequence
      getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'patMETs'+postfix),
                                                                getattr(process, 'patMETs'+postfix)+getattr(process, 'metJESCorAK5PFTypeI'+postfix)+getattr(process, 'patMETsTypeI'+postfix)
                                                            )
      
    process.out.outputCommands.append( 'keep double_*' + postfix + '*_*_' + process.name_() )

###
### Selection configuration
###

if runStandardPAT:

  ### Muons
  process.patMuons.usePV      = muonsUsePV
  process.patMuons.embedTrack = muonEmbedTrack

  process.selectedPatMuons.cut = muonCut

  ### Electrons
  process.patElectrons.electronIDSources = electronIDSources
  process.patElectrons.usePV = electronsUsePV

  process.selectedPatElectrons.cut = electronCut

  ### Photons
  process.selectedPatPhotons.cut = photonCut

  ### Jets
  process.selectedPatJets.cut = 'pt > 10.'
  if usePFJets:
	  getattr(process, 'selectedPatJets'+jetAlgo+pfSuffix).cut = 'pt > 10.'


if runPF2PAT:

  applyPostfix( process, 'patMuons', postfix ).usePV      = muonsUsePV
  applyPostfix( process, 'patMuons', postfix ).embedTrack = muonEmbedTrack

  applyPostfix( process, 'selectedPatMuons', postfix ).cut = muonCutPF

  ### Electrons

  applyPostfix( process, 'patElectrons', postfix ).electronIDSources = electronIDSources
  applyPostfix( process, 'patElectrons', postfix ).usePV = electronsUsePV

  applyPostfix( process, 'selectedPatElectrons', postfix ).cut = electronCutPF

  ### Jets
  applyPostfix( process, 'selectedPatJets', postfix ).cut = 'pt > 10.'

###
### Additional Isolation Values: User Isolation
###
if runPF2PAT:

  #Muon,Charged
  isoValMuonWithChargedDR03 = applyPostfix(process, 'muPFIsoValueCharged03', postfix).clone()
  isoValMuonWithChargedDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValMuonWithChargedDR03'+postfix,isoValMuonWithChargedDR03)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValueCharged03'+postfix),
                                                            getattr(process, 'muPFIsoValueCharged03'+postfix)+getattr(process, 'isoValMuonWithChargedDR03'+postfix)
                                                            )
  #Muon, Neutral
  isoValMuonWithNeutralDR03 = applyPostfix(process, 'muPFIsoValueNeutral03', postfix).clone()
  isoValMuonWithNeutralDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValMuonWithNeutralDR03'+postfix,isoValMuonWithNeutralDR03)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValueNeutral03'+postfix),
                                                            getattr(process, 'muPFIsoValueNeutral03'+postfix)+getattr(process, 'isoValMuonWithNeutralDR03'+postfix)
                                                            )
  #Muon, Photons
  isoValMuonWithPhotonsDR03 = applyPostfix(process, 'muPFIsoValueGamma03', postfix).clone()
  isoValMuonWithPhotonsDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValMuonWithPhotonsDR03'+postfix,isoValMuonWithPhotonsDR03)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValueGamma03'+postfix),
                                                            getattr(process, 'muPFIsoValueGamma03'+postfix)+getattr(process, 'isoValMuonWithPhotonsDR03'+postfix)
                                                            )
  #Muon, PU
  isoValMuonWithPUDR03 = applyPostfix(process, 'muPFIsoValuePU03', postfix).clone()
  isoValMuonWithPUDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValMuonWithPUDR03'+postfix,isoValMuonWithPUDR03)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValuePU03'+postfix),
                                                            getattr(process, 'muPFIsoValuePU03'+postfix)+getattr(process, 'isoValMuonWithPUDR03'+postfix)
                                                            )
  
  #Muon,Charged
  isoValMuonWithChargedDR04 = applyPostfix(process, 'muPFIsoValueCharged04', postfix).clone()
  isoValMuonWithChargedDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValMuonWithChargedDR04'+postfix,isoValMuonWithChargedDR04)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValueCharged04'+postfix),
                                                            getattr(process, 'muPFIsoValueCharged04'+postfix)+getattr(process, 'isoValMuonWithChargedDR04'+postfix)
                                                            )
  #Muon, Neutral
  isoValMuonWithNeutralDR04 = applyPostfix(process, 'muPFIsoValueNeutral04', postfix).clone()
  isoValMuonWithNeutralDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValMuonWithNeutralDR04'+postfix,isoValMuonWithNeutralDR04)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValueNeutral04'+postfix),
                                                            getattr(process, 'muPFIsoValueNeutral04'+postfix)+getattr(process, 'isoValMuonWithNeutralDR04'+postfix)
                                                            )
  #Muon, Photons
  isoValMuonWithPhotonsDR04 = applyPostfix(process, 'muPFIsoValueGamma04', postfix).clone()
  isoValMuonWithPhotonsDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValMuonWithPhotonsDR04'+postfix,isoValMuonWithPhotonsDR04)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValueGamma04'+postfix),
                                                            getattr(process, 'muPFIsoValueGamma04'+postfix)+getattr(process, 'isoValMuonWithPhotonsDR04'+postfix)
                                                            )
  #Muon, PU
  isoValMuonWithPUDR04 = applyPostfix(process, 'muPFIsoValuePU04', postfix).clone()
  isoValMuonWithPUDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValMuonWithPUDR04'+postfix,isoValMuonWithPUDR04)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValuePU04'+postfix),
                                                            getattr(process, 'muPFIsoValuePU04'+postfix)+getattr(process, 'isoValMuonWithPUDR04'+postfix)
                                                            )
  applyPostfix(process,'patMuons',postfix).isolationValues.user = cms.VInputTag("isoValMuonWithChargedDR03"+postfix
                                                                                , "isoValMuonWithNeutralDR03"+postfix
                                                                                , "isoValMuonWithPhotonsDR03"+postfix                                                                                                                                           , "isoValMuonWithPUDR03"+postfix                                                           
                                                                                , "isoValMuonWithChargedDR04"+postfix
                                                                                , "isoValMuonWithNeutralDR04"+postfix
                                                                                , "isoValMuonWithPhotonsDR04"+postfix                                                                                                                                           , "isoValMuonWithPUDR04"+postfix                                                    
                                                                                )

  #Electron, Charged
  isoValElectronWithChargedDR03 = applyPostfix(process, 'elPFIsoValueCharged03', postfix).clone()
  isoValElectronWithChargedDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValElectronWithChargedDR03'+postfix,isoValElectronWithChargedDR03)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueCharged03'+postfix),
                                                            getattr(process, 'elPFIsoValueCharged03'+postfix)+getattr(process, 'isoValElectronWithChargedDR03'+postfix)
                                                            )
  #Electron, Neutral
  isoValElectronWithNeutralDR03 = applyPostfix(process, 'elPFIsoValueNeutral03', postfix).clone()
  isoValElectronWithNeutralDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValElectronWithNeutralDR03'+postfix,isoValElectronWithNeutralDR03)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueNeutral03'+postfix),
                                                            getattr(process, 'elPFIsoValueNeutral03'+postfix)+getattr(process, 'isoValElectronWithNeutralDR03'+postfix)
                                                            )
  #Electron, Photons
  isoValElectronWithPhotonsDR03 = applyPostfix(process, 'elPFIsoValueGamma03', postfix).clone()
  isoValElectronWithPhotonsDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValElectronWithPhotonsDR03'+postfix,isoValElectronWithPhotonsDR03)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueGamma03'+postfix),
                                                            getattr(process, 'elPFIsoValueGamma03'+postfix)+getattr(process, 'isoValElectronWithPhotonsDR03'+postfix)
                                                            )
  #Electron, PU
  isoValElectronWithPUDR03 = applyPostfix(process, 'elPFIsoValuePU03', postfix).clone()
  isoValElectronWithPUDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValElectronWithPUDR03'+postfix,isoValElectronWithPUDR03)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValuePU03'+postfix),
                                                            getattr(process, 'elPFIsoValuePU03'+postfix)+getattr(process, 'isoValElectronWithPUDR03'+postfix)
                                                            )
  
  #Electron, Charged
  isoValElectronWithChargedDR04 = applyPostfix(process, 'elPFIsoValueCharged04', postfix).clone()
  isoValElectronWithChargedDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValElectronWithChargedDR04'+postfix,isoValElectronWithChargedDR04)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueCharged04'+postfix),
                                                            getattr(process, 'elPFIsoValueCharged04'+postfix)+getattr(process, 'isoValElectronWithChargedDR04'+postfix)
                                                            )
  #Electron, Neutral
  isoValElectronWithNeutralDR04 = applyPostfix(process, 'elPFIsoValueNeutral04', postfix).clone()
  isoValElectronWithNeutralDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValElectronWithNeutralDR04'+postfix,isoValElectronWithNeutralDR04)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueNeutral04'+postfix),
                                                            getattr(process, 'elPFIsoValueNeutral04'+postfix)+getattr(process, 'isoValElectronWithNeutralDR04'+postfix)
                                                            )
  #Electron, Photons
  isoValElectronWithPhotonsDR04 = applyPostfix(process, 'elPFIsoValueGamma04', postfix).clone()
  isoValElectronWithPhotonsDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValElectronWithPhotonsDR04'+postfix,isoValElectronWithPhotonsDR04)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueGamma04'+postfix),
                                                            getattr(process, 'elPFIsoValueGamma04'+postfix)+getattr(process, 'isoValElectronWithPhotonsDR04'+postfix)
                                                            )
  #Electron, PU
  isoValElectronWithPUDR04 = applyPostfix(process, 'elPFIsoValuePU04', postfix).clone()
  isoValElectronWithPUDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValElectronWithPUDR04'+postfix,isoValElectronWithPUDR04)
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValuePU04'+postfix),
                                                            getattr(process, 'elPFIsoValuePU04'+postfix)+getattr(process, 'isoValElectronWithPUDR04'+postfix)
                                                            )
 
  applyPostfix(process,'patElectrons',postfix).isolationValues.user = cms.VInputTag("isoValElectronWithChargedDR03"+postfix
                                                                                , "isoValElectronWithNeutralDR03"+postfix
                                                                                , "isoValElectronWithPhotonsDR03"+postfix
                                                                                , "isoValElectronWithPUDR03"+postfix
                                                                                , "isoValElectronWithChargedDR04"+postfix
                                                                                , "isoValElectronWithNeutralDR04"+postfix
                                                                                , "isoValElectronWithPhotonsDR04"+postfix
                                                                                , "isoValElectronWithPUDR04"+postfix
                                                                                )




  from PhysicsTools.PatAlgos.tools.pfTools import adaptPFMuons
  from PhysicsTools.PatAlgos.tools.pfTools import adaptPFElectrons
  
  ## Muons
  pfIsolatedMuonsLoose = applyPostfix(process, 'pfIsolatedMuons', postfix).clone()
  pfIsolatedMuonsLoose.isolationCut = cms.double(999.0)
  pfIsolatedMuonsLoose.combinedIsolationCut = cms.double(999.0)
  setattr(process,'pfIsolatedMuonsLoose'+postfix,pfIsolatedMuonsLoose)

  patMuonsLoose = applyPostfix(process, 'patMuons', postfix).clone()
  patMuonsLoose.addResolutions = cms.bool(False)
  patMuonsLoose.pfMuonSource = cms.InputTag( 'pfIsolatedMuonsLoose' + postfix )
  patMuonsLoose.genParticleMatch = cms.InputTag( 'muonMatchLoose' + postfix )

  setattr(process,'patMuonsLoose'+postfix,patMuonsLoose)

  adaptPFMuons( process, getattr(process, 'patMuonsLoose'+postfix), postfix )

  getattr(process, 'patMuonsLoose'+postfix ).isolationValues.user = cms.VInputTag("isoValMuonWithChargedDR03"+postfix
                                                                                  , "isoValMuonWithNeutralDR03"+postfix
                                                                                  , "isoValMuonWithPhotonsDR03"+postfix                                                                                                                                     , "isoValMuonWithPUDR03"+postfix                                                        
                                                                                  , "isoValMuonWithChargedDR04"+postfix
                                                                                  , "isoValMuonWithNeutralDR04"+postfix
                                                                                  , "isoValMuonWithPhotonsDR04"+postfix                                                                                                                                     , "isoValMuonWithPUDR04"+postfix                                                    
    )
  

  muonMatchLoose = applyPostfix(process, 'muonMatch', postfix).clone()
  muonMatchLoose.src = cms.InputTag( 'pfIsolatedMuonsLoose' + postfix )
  setattr(process,'muonMatchLoose'+postfix,muonMatchLoose)

  applyPostfix( process, 'muonMatch', postfix ).src = cms.InputTag( 'pfIsolatedMuons' + postfix )

  selectedPatMuonsLoose = applyPostfix(process, 'selectedPatMuons', postfix).clone()
  selectedPatMuonsLoose.src = cms.InputTag( 'patMuonsLoose' + postfix )
  selectedPatMuonsLoose.cut = muonCutLoosePF
  setattr(process,'selectedPatMuonsLoose'+postfix,selectedPatMuonsLoose)

  ## Electrons
  pfIsolatedElectronsLoose = applyPostfix(process, 'pfIsolatedElectrons', postfix).clone()
  pfIsolatedElectronsLoose.isolationCut = cms.double(999.0)
  pfIsolatedElectronsLoose.combinedIsolationCut = cms.double(999.0)
  setattr(process,'pfIsolatedElectronsLoose'+postfix,pfIsolatedElectronsLoose)

  patElectronsLoose = applyPostfix(process, 'patElectrons', postfix).clone()
  patElectronsLoose.pfElectronSource = cms.InputTag( 'pfIsolatedElectronsLoose' + postfix )

  setattr(process,'patElectronsLoose'+postfix,patElectronsLoose)

  adaptPFElectrons( process, getattr(process, 'patElectronsLoose'+postfix), postfix )

  getattr(process, 'patElectronsLoose'+postfix ).isolationValues.user = cms.VInputTag("isoValElectronWithChargedDR03"+postfix
                                                                                      , "isoValElectronWithNeutralDR03"+postfix
                                                                                      , "isoValElectronWithPhotonsDR03"+postfix
                                                                                      , "isoValElectronWithPUDR03"+postfix
                                                                                      , "isoValElectronWithChargedDR04"+postfix
                                                                                      , "isoValElectronWithNeutralDR04"+postfix
                                                                                      , "isoValElectronWithPhotonsDR04"+postfix
                                                                                      , "isoValElectronWithPUDR04"+postfix
    )

  selectedPatElectronsLoose = applyPostfix(process, 'selectedPatElectrons', postfix).clone()
  selectedPatElectronsLoose.src = cms.InputTag( 'patElectronsLoose' + postfix )
  selectedPatElectronsLoose.cut = electronCutLoosePF
  setattr(process,'selectedPatElectronsLoose'+postfix,selectedPatElectronsLoose)


  if not runOnMC:
    process.looseLeptonSequence = cms.Sequence(
      getattr(process, 'pfIsolatedMuonsLoose'+postfix) +
      getattr(process, 'patMuonsLoose'+postfix) +
      getattr(process, 'selectedPatMuonsLoose'+postfix) +
      getattr(process, 'pfIsolatedElectronsLoose'+postfix) +
      getattr(process, 'patElectronsLoose'+postfix) +
      getattr(process, 'selectedPatElectronsLoose'+postfix)
      )

  if runOnMC:
    process.looseLeptonSequence = cms.Sequence(
      getattr(process, 'pfIsolatedMuonsLoose'+postfix) +
      getattr(process, 'muonMatchLoose'+postfix) +
      getattr(process, 'patMuonsLoose'+postfix) +
      getattr(process, 'selectedPatMuonsLoose'+postfix) +
      getattr(process, 'pfIsolatedElectronsLoose'+postfix) +
      getattr(process, 'patElectronsLoose'+postfix) +
      getattr(process, 'selectedPatElectronsLoose'+postfix)
      )


###
### Darren Specific Stuff
###

process.BNproducer = cms.EDProducer('BEANmaker',
        calometTag = cms.InputTag("patMETsAK5Calo"), 
        pfmetTag = cms.InputTag("patMETsTypeIPFlow"), 
        tcmetTag = cms.InputTag("patMETsAK5TC"), 
        eleTag = cms.InputTag("selectedPatElectrons"),
        pfeleTag = cms.InputTag("selectedPatElectronsLoosePFlow"),
        genParticleTag = cms.InputTag("genParticles"),
        calojetTag = cms.InputTag("selectedPatJetsAK5Calo"), 
        pfjetTag = cms.InputTag("selectedPatJetsPFlow"), 
        jptjetTag = cms.InputTag("selectedPatJetsAK5JPT"), 
        muonTag = cms.InputTag("selectedPatMuons"),
        pfmuonTag = cms.InputTag("selectedPatMuonsLoosePFlow"),
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
        hltProcessName = cms.string("HLT"),  
        eventWeight = cms.double(0.0126976903794279),
        minSCEt = cms.double(10),
        minPhotonEt = cms.double(10),
        minJetPt = cms.double(10),
        minTrackPt = cms.double(10),
        verbose = cms.bool(True),
        sample = cms.int32(-2500)
)


process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

process.load('RecoMuon/MuonIdentification/refitMuons_cfi')

process.skimMuon = cms.EDFilter("SkimMuon",
                                ptCutMu1 = cms.untracked.double(15),
                                ptCutMu2 = cms.untracked.double(0),
                                caloIsoCut = cms.untracked.double(999999999),
                                trackIsoCut = cms.untracked.double(999999999),
                                minNumGoodMu = cms.untracked.int32(1)
                                )


###
### Scheduling
###

# CiC electron ID

process.load( "RecoEgamma.ElectronIdentification.cutsInCategoriesElectronIdentificationV06_cfi" )
process.eidCiCSequence = cms.Sequence(
  process.eidVeryLooseMC
+ process.eidLooseMC
+ process.eidMediumMC
+ process.eidTightMC
+ process.eidSuperTightMC
+ process.eidHyperTight1MC
+ process.eidHyperTight2MC
+ process.eidHyperTight3MC
+ process.eidHyperTight4MC
)


# Filter for bad LHE events
process.load("GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi")
process.totalKinematicsFilter.tolerance = cms.double(5.0)


# The paths
if runStandardPAT:
  if useCaloJets:
    process.p = cms.Path()
    if useSkimMuon:
      process.p += process.skimMuon
    if not runOnMC:
      process.p += process.eventCleaning
    else :
      process.p += process.totalKinematicsFilter
    if useTrigger:
      process.p += process.step0a
    process.p += process.goodOfflinePrimaryVertices
    process.p += process.eidCiCSequence
    process.p += process.patDefaultSequence
    if usePFJets:
      process.p.remove( getattr( process, 'patJetCorrFactors'                    + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'jetTracksAssociatorAtVertex'          + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'impactParameterTagInfos'              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'secondaryVertexTagInfos'              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'softMuonTagInfos'                     + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'jetBProbabilityBJetTags'              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'jetProbabilityBJetTags'               + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'trackCountingHighPurBJetTags'         + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'trackCountingHighEffBJetTags'         + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'simpleSecondaryVertexHighEffBJetTags' + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'simpleSecondaryVertexHighPurBJetTags' + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'combinedSecondaryVertexBJetTags'      + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'combinedSecondaryVertexMVABJetTags'   + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'softMuonBJetTags'                     + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'softMuonByPtBJetTags'                 + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'softMuonByIP3dBJetTags'               + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJetCharge'                         + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJetPartonMatch'                    + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJetGenJetMatch'                    + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJetPartonAssociation'              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJetFlavourAssociation'             + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJets'                              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patMETs'                              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'selectedPatJets'                      + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'countPatJets'                         + jetAlgo + pfSuffix ) )
    process.out.SelectEvents.SelectEvents.append( 'p' )

  if usePFJets:

    pAddPF = cms.Path()
    if useSkimMuon:
      process.pAddPF += process.skimMuon
    if not runOnMC:
      pAddPF += process.eventCleaning
    else :
      pAddPF += process.totalKinematicsFilter
    if useTrigger:
      pAddPF += process.step0a
    pAddPF += process.goodOfflinePrimaryVertices
    pAddPF += process.eidCiCSequence
    pAddPF += process.patDefaultSequence
    pAddPF.remove( process.patJetCorrFactors )
    pAddPF.remove( process.patJetCharge )
    pAddPF.remove( process.patJetPartonMatch )
    pAddPF.remove( process.patJetGenJetMatch )
    pAddPF.remove( process.patJetPartonAssociation )
    pAddPF.remove( process.patJetFlavourAssociation )
    pAddPF.remove( process.patJets )
    pAddPF.remove( process.patMETs )
    pAddPF.remove( process.selectedPatJets )
    pAddPF.remove( process.countPatJets )
    setattr( process, 'p' + jetAlgo + pfSuffix, pAddPF )
    process.out.SelectEvents.SelectEvents.append( 'p' + jetAlgo + pfSuffix )

if runPF2PAT:
  pPF = cms.Path()
  if useSkimMuon:
    pPF += process.skimMuon
  if not runOnMC:
    pPF += process.eventCleaning
  else:
    pPF += process.totalKinematicsFilter
  if useTrigger:
    pPF += process.step0a
  if keepTK:
    pPF += process.highPtTracks
  pPF += process.goodOfflinePrimaryVertices
  pPF += process.eidCiCSequence
  pPF += getattr( process, 'patPF2PATSequence' + postfix )
  ## new
  pPF += process.looseLeptonSequence
  ## end new
  pPF.remove( getattr( process, 'patTrigger' + postfix ) )
  pPF.remove( getattr( process, 'patTriggerEvent' + postfix ) )
  setattr( process, 'p' + postfix, pPF )
  process.out.SelectEvents.SelectEvents.append( 'p' + postfix )


if runStandardPAT and runPF2PAT:
  pB = cms.Path()
  if useSkimMuon:
    pB += process.skimMuon
  if not runOnMC:
    pB += process.eventCleaning
    pB += process.HBHENoiseFilterResultProducer
  else :
    pB += process.totalKinematicsFilter
  if useTrigger:
    pB += process.step0a
  pB += process.refitMuons
  pB += process.BNproducer
  setattr( process, 'pBeans', pB )
  process.out.SelectEvents.SelectEvents.append( 'pBeans' )
  process.out.outputCommands = [ 'drop *' ]
  process.out.outputCommands.extend( [ # BEAN Objects
                                    'keep *_BNproducer_*_*',
   				    #'keep *',
                                    ] )


#Set the preselection
process.pfSelectedElectrons.cut = cms.string('et > 10. && abs(eta) < 2.5')
process.pfSelectedMuons.cut = cms.string('pt > 10. && abs(eta) < 2.5')

# Additional drop commands
process.out.outputCommands.append( 'drop *_selectedPatPFParticlesPFlow_*_*')
process.out.outputCommands.append( 'drop *_*_caloTowers_*')


#outfile = open('config_pat_test.py','w')
#print >> outfile,process.dumpPython()
#outfile.close()
