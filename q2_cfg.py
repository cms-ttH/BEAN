import FWCore.ParameterSet.Config as cms

process = cms.Process("q2weights")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'WARNING'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/Summer12/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S6_START52_V9-v1/0000/004B4749-A1B8-E111-B16C-002618943896.root'
    ),
     skipEvents = cms.untracked.uint32(0)
)

process.q2weights = cms.EDProducer('Q2Weights'
)

#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('tree.root')
#)

process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('cmssw-out.root'),
#        SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
        outputCommands = cms.untracked.vstring('keep *')
      )


process.p = cms.Path(process.q2weights)
process.outpath = cms.EndPath( process.out )

