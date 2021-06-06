import FWCore.ParameterSet.Config as cms

process = cms.Process('SE')

# input files
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        '/store/data/Run2012C/SingleElectron/AOD/PromptReco-v1/000/198/487/FE46A94D-FCCA-E111-92D5-003048CF9B28.root'
    )
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# select events
process.selectEvent = cms.EDFilter('SelectEvent',
    runNrs = cms.vuint32(198487),
    lumiNrs = cms.vuint32(1500),
    eventNrs = cms.vuint32(1495052793)
)

# output module
process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('events.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('sp'))
)

process.sp = cms.Path(process.selectEvent)
process.op = cms.EndPath(process.out)
