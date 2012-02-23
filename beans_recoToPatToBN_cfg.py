
#
# This file contains the Top PAG reference selection work-flow for mu + jets analysis.
# as defined in
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopLeptonPlusJetsRefSel_mu#Selection_Version_SelV4_valid_fr
#

import sys

import FWCore.ParameterSet.Config as cms

process = cms.Process( 'PAT' )


### ======================================================================== ###
###                                                                          ###
###                                 Constants                                ###
###                            (user job steering)                           ###
###                                                                          ###
### ======================================================================== ###


### Data or MC?
runOnMC = True

### Standard and PF work flow

# Standard
runStandardPAT = True
usePFJets      = True
useCaloJets    = True

# PF2PAT
runPF2PAT = True

### Switch on/off selection steps

# Whether you want to include trigger skimming in this selection
useTrigger      = False

# Whether you want to include trigger skimming in this selection
useSkimMuon = False

### Reference selection

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets import *
muonsUsePV             = True
electronsUsePV         = True
#muonEmbedTrack         = True
muonCut                  = 'isGlobalMuon && pt > 10. && abs(eta) < 2.5'
muonCutPF                = 'isGlobalMuon && pt > 10. && abs(eta) < 2.5'
#looseMuonCutPF           = ''
#tightMuonCutPF           = ''
#muonJetsDR             = 0.3
#jetCutPF               = ''
#jetMuonsDRPF           = 0.1

photonCut              = 'et > 10. && abs(eta) < 3.5'

electronCut              = 'et > 15. && abs(eta) < 2.5'
electronCutPF            = 'et > 15. && abs(eta) < 2.5 && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits<2'

#Re-jigger the PF muon and PF electron isolation cuts...
# PF muon: used for loose veto
muonCutPF  =  muonCutBase
muonCutPF += ' && (chargedHadronIso+neutralHadronIso+photonIso)/pt < 0.20'       # relative isolation

# PF electron
electronCutPF  = electronCutBase
electronCutPF += ' && (chargedHadronIso+neutralHadronIso+photonIso)/et < 0.20' # relative isolation


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
pfMuonCombIsoCut = 0.20 #Loose lepton isolation
# electrons
pfElectronSelectionCut  = 'et > 15. && abs(eta) < 2.5'
pfElectronIsoConeR   = 0.4
pfElectronCombIsoCut = 0.20

pfMuonIsoConeR03 = False
pfElectronnIsoConeR03 = False

### JEC levels

# levels to be accessible from the jets
# jets are corrected to L3Absolute (MC), L2L3Residual (data) automatically, if enabled here
# and remain uncorrected, if none of these levels is enabled here
useL1FastJet   = True  # needs useL1Offset being off, error otherwise
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
  '/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/FEE3E76C-1EFA-E011-BAF4-002618943864.root',
  #'/store/data/Run2011B/SingleMu/RECO/PromptReco-v1/000/175/832/10EBFA08-84DB-E011-9D72-003048D2C020.root',
#  '/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/FE7F7567-13FA-E011-86E3-003048678F1C.root',
#  '/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/FCB456E6-3EFA-E011-9020-0018F3D0961A.root',
#  '/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/FC2D57E7-2DFA-E011-B0A0-0018F3D096C0.root',
#  '/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/FA602C8E-3AFA-E011-811C-001A9281174A.root',
  ] # overwritten, if "useRelVals" is 'True'

# maximum number of events
maxInputEvents = 1000 # reduce for testing

### Conditions

# GlobalTags (w/o suffix '::All')
globalTagData = 'GR_R_42_V19' 
globalTagMC   = 'START42_V13' 

### Output

# output file
outputFile = 'pat_42x_fall11_withANDfilter_ttbar.root'

# event frequency of Fwk report
fwkReportEvery = 100 # reduce for testing

# switch for 'TrigReport'/'TimeReport' at job end
wantSummary = True


###                              End of constants                            ###
###                                                                          ###
### ======================================================================== ###


###
###   Try to be tricky and override some settings from a crab config
###

#Detect that this is a CRAB submission:
keepSC = False
keepTK = False
if ('crab.py' in sys.argv[0]) and (sys.argv[1] == '-create'):

  print "This config is being used for CRAB!"
  print "Checking to see if the CRAB config file contains SQWaT parameters"

  index = -1
  try:
    index = sys.argv.index('-cfg') + 1
  except:
    print 'Error cannot find the CRAB config file.  No "-cfg" argument found.  Command line args:'
    print sys.argv

  if (index < 0) or (index >= len(sys.argv)):
    print 'Could not identify the CRAB config file.  Not extracting any special parameters from there.'
  else:

    #Read in the crab config file and process it for any special SQWaT arguments
    import ConfigParser
    cp = ConfigParser.ConfigParser()
    cp.read(sys.argv[index])

    if not cp.has_section('SQWaT') :
      print 'No SQWaT section in %s' % sys.argv[index]
      print 'Not extracting any parameters'
    else:
      if cp.has_option('SQWaT','runOnMC'):
        runOnMC = cp.getboolean('SQWaT','runOnMC')
        print 'Setting runOnMC to %s' % runOnMC

      if cp.has_option('SQWaT','useTrigger'):
        useTrigger = cp.getboolean('SQWaT','useTrigger')
        print 'Setting useTrigger to %s' % useTrigger
    
      if cp.has_option('SQWaT','triggerSelection'):
        triggerSelectionData = cp.get('SQWaT','triggerSelection')
        triggerSelectionMC = triggerSelectionData # Assume that we want this for whichever we're running, if we're overriding it here
        print 'Setting triggerSelection to %s' % triggerSelectionData

      if cp.has_option('SQWaT','keepSuperClusters'):
        keepSC = cp.getboolean('SQWaT','keepSuperClusters')
        print 'Setting keepSC to %s' % keepSC
      if cp.has_option('SQWaT','keeptracks'):
        keepTK = cp.getboolean('SQWaT','keeptracks')
        print 'Setting keepTK to %s' % keepTK
###
### Basic configuration
###

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
  applyPostfix( process, 'pfIsolatedMuons'      , postfix ).isolationCut = pfMuonCombIsoCut
  # applyPostfix( process, 'pfIsolatedMuons'      , postfix ).combinedIsolationCut = pfMuonCombIsoCut
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
  #  applyPostfix( process, 'pfIsolatedElectrons'    , postfix ).combinedIsolationCut = pfElectronCombIsoCut
  applyPostfix( process, 'pfIsolatedElectrons'      , postfix ).isolationCut = pfElectronCombIsoCut
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
    cut = cms.string ('pt > 15')
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
#                                      , process.kt6PFJets * process.kt6PFJetsForIsolation * getattr(process, 'ak5PFJets') * process.patJetCorrFactors
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
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValueNeutral03'+postfix),
                                                            getattr(process, 'muPFIsoValueNeutral03'+postfix)+getattr(process, 'isoValMuonWithNeutralDR03'+postfix)
                                                            )
  #Muon, Photons
  isoValMuonWithPhotonsDR03 = applyPostfix(process, 'muPFIsoValueGamma03', postfix).clone()
  isoValMuonWithPhotonsDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValMuonWithPhotonsDR03'+postfix,isoValMuonWithPhotonsDR03)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValueGamma03'+postfix),
                                                            getattr(process, 'muPFIsoValueGamma03'+postfix)+getattr(process, 'isoValMuonWithPhotonsDR03'+postfix)
                                                            )
  #Muon, PU
  isoValMuonWithPUDR03 = applyPostfix(process, 'muPFIsoValuePU03', postfix).clone()
  isoValMuonWithPUDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValMuonWithPUDR03'+postfix,isoValMuonWithPUDR03)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValuePU03'+postfix),
                                                            getattr(process, 'muPFIsoValuePU03'+postfix)+getattr(process, 'isoValMuonWithPUDR03'+postfix)
                                                            )
  
  #Muon,Charged
  isoValMuonWithChargedDR04 = applyPostfix(process, 'muPFIsoValueCharged04', postfix).clone()
  isoValMuonWithChargedDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValMuonWithChargedDR04'+postfix,isoValMuonWithChargedDR04)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValueCharged04'+postfix),
                                                            getattr(process, 'muPFIsoValueCharged04'+postfix)+getattr(process, 'isoValMuonWithChargedDR04'+postfix)
                                                            )
  #Muon, Neutral
  isoValMuonWithNeutralDR04 = applyPostfix(process, 'muPFIsoValueNeutral04', postfix).clone()
  isoValMuonWithNeutralDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValMuonWithNeutralDR04'+postfix,isoValMuonWithNeutralDR04)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValueNeutral04'+postfix),
                                                            getattr(process, 'muPFIsoValueNeutral04'+postfix)+getattr(process, 'isoValMuonWithNeutralDR04'+postfix)
                                                            )
  #Muon, Photons
  isoValMuonWithPhotonsDR04 = applyPostfix(process, 'muPFIsoValueGamma04', postfix).clone()
  isoValMuonWithPhotonsDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValMuonWithPhotonsDR04'+postfix,isoValMuonWithPhotonsDR04)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'muPFIsoValueGamma04'+postfix),
                                                            getattr(process, 'muPFIsoValueGamma04'+postfix)+getattr(process, 'isoValMuonWithPhotonsDR04'+postfix)
                                                            )
  #Muon, PU
  isoValMuonWithPUDR04 = applyPostfix(process, 'muPFIsoValuePU04', postfix).clone()
  isoValMuonWithPUDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValMuonWithPUDR04'+postfix,isoValMuonWithPUDR04)
  # Install alternative isolation in path
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
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueCharged03'+postfix),
                                                            getattr(process, 'elPFIsoValueCharged03'+postfix)+getattr(process, 'isoValElectronWithChargedDR03'+postfix)
                                                            )
  #Electron, Neutral
  isoValElectronWithNeutralDR03 = applyPostfix(process, 'elPFIsoValueNeutral03', postfix).clone()
  isoValElectronWithNeutralDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValElectronWithNeutralDR03'+postfix,isoValElectronWithNeutralDR03)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueNeutral03'+postfix),
                                                            getattr(process, 'elPFIsoValueNeutral03'+postfix)+getattr(process, 'isoValElectronWithNeutralDR03'+postfix)
                                                            )
  #Electron, Photons
  isoValElectronWithPhotonsDR03 = applyPostfix(process, 'elPFIsoValueGamma03', postfix).clone()
  isoValElectronWithPhotonsDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValElectronWithPhotonsDR03'+postfix,isoValElectronWithPhotonsDR03)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueGamma03'+postfix),
                                                            getattr(process, 'elPFIsoValueGamma03'+postfix)+getattr(process, 'isoValElectronWithPhotonsDR03'+postfix)
                                                            )
  #Electron, PU
  isoValElectronWithPUDR03 = applyPostfix(process, 'elPFIsoValuePU03', postfix).clone()
  isoValElectronWithPUDR03.deposits[0].deltaR = 0.3
  setattr(process,'isoValElectronWithPUDR03'+postfix,isoValElectronWithPUDR03)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValuePU03'+postfix),
                                                            getattr(process, 'elPFIsoValuePU03'+postfix)+getattr(process, 'isoValElectronWithPUDR03'+postfix)
                                                            )
  
  #Electron, Charged
  isoValElectronWithChargedDR04 = applyPostfix(process, 'elPFIsoValueCharged04', postfix).clone()
  isoValElectronWithChargedDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValElectronWithChargedDR04'+postfix,isoValElectronWithChargedDR04)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueCharged04'+postfix),
                                                            getattr(process, 'elPFIsoValueCharged04'+postfix)+getattr(process, 'isoValElectronWithChargedDR04'+postfix)
                                                            )
  #Electron, Neutral
  isoValElectronWithNeutralDR04 = applyPostfix(process, 'elPFIsoValueNeutral04', postfix).clone()
  isoValElectronWithNeutralDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValElectronWithNeutralDR04'+postfix,isoValElectronWithNeutralDR04)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueNeutral04'+postfix),
                                                            getattr(process, 'elPFIsoValueNeutral04'+postfix)+getattr(process, 'isoValElectronWithNeutralDR04'+postfix)
                                                            )
  #Electron, Photons
  isoValElectronWithPhotonsDR04 = applyPostfix(process, 'elPFIsoValueGamma04', postfix).clone()
  isoValElectronWithPhotonsDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValElectronWithPhotonsDR04'+postfix,isoValElectronWithPhotonsDR04)
  # Install alternative isolation in path
  getattr( process, 'patPF2PATSequence' + postfix ).replace(getattr(process, 'elPFIsoValueGamma04'+postfix),
                                                            getattr(process, 'elPFIsoValueGamma04'+postfix)+getattr(process, 'isoValElectronWithPhotonsDR04'+postfix)
                                                            )
  #Electron, PU
  isoValElectronWithPUDR04 = applyPostfix(process, 'elPFIsoValuePU04', postfix).clone()
  isoValElectronWithPUDR04.deposits[0].deltaR = 0.4
  setattr(process,'isoValElectronWithPUDR04'+postfix,isoValElectronWithPUDR04)
  # Install alternative isolation in path
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



###
### Darren Specific Stuff
###

usePatPFMETtype1 = True

#process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
#process.selectedPatJetsForMETtype1p2Corr.src = cms.InputTag('selectedPatJetsPFlow')
#process.selectedPatJetsForMETtype2Corr.src = cms.InputTag('selectedPatJetsPFlow')
## DEFAULT
#process.patPFJetMETtype1p2Corr.type1JetPtThreshold = cms.double(10.0)
## DEFAULT
# skipEMfractionThreshold = cms.double(0.90)
# skipMuonSelection = cms.string("isGlobalMuon | isStandAloneMuon")
#process.patPFJetMETtype1p2Corr.skipEM = cms.bool(False)
#process.patPFJetMETtype1p2Corr.skipMuons = cms.bool(False)
#if not runOnMC:
#  process.patPFJetMETtype1p2Corr.jetCorrLabel = 'L2L3Residual'
#  process.patPFMet.addGenMET = cms.bool(False)

process.BNproducer = cms.EDProducer('BEANmaker',
        calometTag = cms.InputTag("patMETsAK5Calo"), 
        pfmetTag = cms.InputTag("patMETsTypeIPFlow"), 
        #pfmetTag = cms.InputTag("patMETsPFlow"), 
        #        pfmetTag = cms.InputTag("patMETsTypeIPF"), 
        tcmetTag = cms.InputTag("patMETsAK5TC"), 
        eleTag = cms.InputTag("selectedPatElectrons"),
        pfeleTag = cms.InputTag("selectedPatElectronsPFlow"),
        genParticleTag = cms.InputTag("genParticles"),
        calojetTag = cms.InputTag("selectedPatJetsAK5Calo"), 
        pfjetTag = cms.InputTag("selectedPatJetsPFlow"), 
        jptjetTag = cms.InputTag("selectedPatJetsAK5JPT"), 
        muonTag = cms.InputTag("selectedPatMuons"),
        pfmuonTag = cms.InputTag("selectedPatMuonsPFlow"),
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
        eventWeight = cms.double(0.00043132197696737),
        minSCEt = cms.double(10),
        minPhotonEt = cms.double(10),
        minJetPt = cms.double(10),
        minTrackPt = cms.double(10),
        verbose = cms.bool(True),
        sample = cms.int32(2523)
)


process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
#process.p = cms.Path(process.HBHENoiseFilterResultProducer)

process.load('RecoMuon/MuonIdentification/refitMuons_cfi')

process.skimMuon = cms.EDFilter("SkimMuon",
                                ptCutMu1 = cms.untracked.double(30),
                                ptCutMu2 = cms.untracked.double(5),
                                caloIsoCut = cms.untracked.double(999999999),
                                trackIsoCut = cms.untracked.double(999999999),
                                minNumGoodMu = cms.untracked.int32(2)
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

  #  eidTight            = cms.InputTag( 'eidTight' )
  #, eidLoose            = cms.InputTag( 'eidLoose' )
  #, eidRobustTight      = cms.InputTag( 'eidRobustTight' )
  #, eidRobustHighEnergy = cms.InputTag( 'eidRobustHighEnergy' )
  #, eidRobustLoose      = cms.InputTag( 'eidRobustLoose' )


# Milano event counting code
#process.load("PhysicsTools.NtupleUtils.AllPassFilter_cfi")
#process.TFileService = cms.Service("TFileService", fileName = cms.string("histo.root") )
#--------------------------
# Counter1: All read events
#process.AllEvents = process.AllPassFilter.clone()
#------------------------------
# Counter2a: Non-scraping events
#process.NonScrapedEvents = process.AllPassFilter.clone()
# Counter2b: Non-scraping events
#process.TotalKinEvents = process.AllPassFilter.clone()
#------------------------------
# Counter3: Trigger
#process.TriggeredEvents = process.AllPassFilter.clone()

# Filter for bad LHE events
process.load("GeneratorInterface.GenFilters.TotalKinematicsFilter_cfi")


# The paths
if runStandardPAT:
  if useCaloJets:
    process.p = cms.Path()
    #process.p += process.AllEvents
    if useSkimMuon:
      process.p += process.skimMuon
    if not runOnMC:
      process.p += process.eventCleaning
      #process.p += process.NonScrapedEvents
    else :
      process.p += process.totalKinematicsFilter
      #process.p += process.TotalKinEvents
    if useTrigger:
      process.p += process.step0a
      #process.p += process.TriggeredEvents
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
    #pAddPF += process.AllEvents
    if not runOnMC:
      pAddPF += process.eventCleaning
      #pAddPF += process.NonScrapedEvents
    else :
      pAddPF += process.totalKinematicsFilter
      #pAddPF += process.TotalKinEvents
    if useTrigger:
      pAddPF += process.step0a
      #pAddPF += process.TriggeredEvents
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
  #pPF += process.AllEvents
  if useSkimMuon:
    pPF += process.skimMuon
  if not runOnMC:
    pPF += process.eventCleaning
    #pPF += process.NonScrapedEvents
  else:
    pPF += process.totalKinematicsFilter
    #pPF += process.TotalKinEvents
  if useTrigger:
    pPF += process.step0a
    #pPF += process.TriggeredEvents
  if keepTK:
    pPF += process.highPtTracks
  pPF += process.goodOfflinePrimaryVertices
  pPF += process.eidCiCSequence
  pPF += getattr( process, 'patPF2PATSequence' + postfix )
  #pPF.remove( process.countPatJets )
  #pPF.remove( getattr( process, 'countPatJets' + postfix ) )
  pPF.remove( getattr( process, 'patTrigger' + postfix ) )
  pPF.remove( getattr( process, 'patTriggerEvent' + postfix ) )
  #pPF += process.BNproducer
  setattr( process, 'p' + postfix, pPF )
  process.out.SelectEvents.SelectEvents.append( 'p' + postfix )


if runStandardPAT and runPF2PAT:
  pB = cms.Path()
  if useSkimMuon:
    pB += process.skimMuon
  if not runOnMC:
    pB += process.eventCleaning
    pB += process.HBHENoiseFilterResultProducer
    #pB += process.NonScrapedEvents
  else :
    pB += process.totalKinematicsFilter
    #pB += process.TotalKinEvents
  if useTrigger:
    pB += process.step0a
    #pB += process.TriggeredEvents
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
process.pfSelectedElectrons.cut = cms.string('et > 15. && abs(eta) < 2.5')
process.pfSelectedMuons.cut = cms.string('pt > 10. && abs(eta) < 2.5')

# Additional drop commands
process.out.outputCommands.append( 'drop *_selectedPatPFParticlesPFlow_*_*')
process.out.outputCommands.append( 'drop *_*_caloTowers_*')






outfile = open('config_pat_test.py','w')
print >> outfile,process.dumpPython()
outfile.close()

