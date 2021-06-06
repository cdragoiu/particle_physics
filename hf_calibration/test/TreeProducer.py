# running parameters -------------------------------------------------------------------------------
runOnMC = True
inputFile = 'file:/uscmst1b_scratch/lpc1/lpctrig/cdragoiu/test.root'
ntupleName = 'ntuple.root'
nrEvents = 1000
# --------------------------------------------------------------------------------------------------

import FWCore.ParameterSet.Config as cms
process = cms.Process('TP')

# message logger
process.load('FWCore.MessageLogger.MessageLogger_cfi')

# options
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# source
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring(inputFile))

# number of events
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(nrEvents))

# detector geometry and conditions
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
if runOnMC: process.GlobalTag.globaltag = cms.string('POSTLS162_V1::All')
else: process.GlobalTag.globaltag = cms.string('FT_53_V21_AN3::All')

# output module
process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('patTuple.root'),
    outputCommands = cms.untracked.vstring('drop *')
)

# tree producer
process.treeProducer = cms.EDAnalyzer('TreeProducer',
    outFileName = cms.string(ntupleName),
    vertexTag = cms.InputTag('offlinePrimaryVertices'),
    electronTag = cms.InputTag('gedGsfElectrons'),
    conversionTag = cms.InputTag('conversions'),
    beamSpotTag = cms.InputTag('offlineBeamSpot'),
    hfElectronTag = cms.InputTag('hfRecoEcalCandidate'),
    hfClusterMapTag = cms.InputTag('hfEMClusters'),
)

# main path
if runOnMC: process.p = cms.Path(process.treeProducer)
else: process.p = cms.Path(process.treeProducer)

# output path
#process.out.outputCommands = cms.untracked.vstring('keep *',
#                                                   'drop recoSecondaryVertexTagInfos_*_*_*',
#                                                   'drop recoTrackIPTagInfos_*_*_*')
#process.op = cms.EndPath(process.out)
