import sys

# running parameters -------------------------------------------------------------------------------
runType = 'data' # mc/data
inputFile = '/store/data/Run2012A/SingleElectron/AOD/13Jul2012-v1/0000/' + \
            '001A2EB8-47D4-E111-B527-003048679070.root'
ntupleName = 'ntuple.root'
nrEvents = 30
if len(sys.argv) > 2: runType = sys.argv[2]
if len(sys.argv) > 3: inputFile = sys.argv[3]
if len(sys.argv) > 4: ntupleName = sys.argv[4]
if len(sys.argv) > 5: nrEvents = int(sys.argv[5])
print '\n\n RUNNING PARAMETERS \n   runType = ' + runType + '\n   inputFile = ' + inputFile + \
      '\n   ntupleName = ' + ntupleName + '\n   nrEvents = ' + str(nrEvents) + '\n\n'
runOnMC = False
if 'mc' in runType: runOnMC = True
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
if runOnMC: process.GlobalTag.globaltag = cms.string('START53_V21::All')
else: process.GlobalTag.globaltag = cms.string('FT_53_V21_AN3::All')

# output module
process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('patTuple.root'),
    outputCommands = cms.untracked.vstring('drop *')
)

# standard PAT configuration
process.load('PhysicsTools.PatAlgos.patSequences_cff')

# JEC for calo jets
inputJetCorrLabel = ['L1Offset', 'L2Relative', 'L3Absolute']
if not runOnMC: inputJetCorrLabel.append('L2L3Residual')
process.patJetCorrFactors.levels = cms.vstring(inputJetCorrLabel)

# good offline PVs (for PF jets)
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter('PrimaryVertexObjectFilter',
    filterParams = pvSelector.clone(minNdof = cms.double(4.0), maxZ = cms.double(24.0)),
    src=cms.InputTag('offlinePrimaryVertices')
)

# use PF2PAT
from PhysicsTools.PatAlgos.tools.pfTools import *
inputJetCorrLabel = ['L1FastJet', 'L2Relative', 'L3Absolute']
if not runOnMC: inputJetCorrLabel.append('L2L3Residual')
usePF2PAT(process, runPF2PAT=True, jetAlgo='AK5', runOnMC=runOnMC, postfix='PFlow',
          jetCorrections=('AK5PFchs', inputJetCorrLabel),
          pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
          typeIMetCorrections=True
)
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)

# top projections in PF2PAT
process.pfNoPileUpPFlow.enable = True
process.pfNoMuonPFlow.enable = False
process.pfNoElectronPFlow.enable = False
process.pfNoTauPFlow.enable = False
process.pfNoJetPFlow.enable = True

# use PF isolation for electrons
usePFIso(process)

# change electron isolation cone to 0.3
process.patElectrons.isolationValues = cms.PSet(
    pfChargedHadrons = cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
    pfChargedAll = cms.InputTag('elPFIsoValueChargedAll03PFIdPFIso'),
    pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU03PFIdPFIso'),
    pfNeutralHadrons = cms.InputTag('elPFIsoValueNeutral03PFIdPFIso'),
    pfPhotons = cms.InputTag('elPFIsoValueGamma03PFIdPFIso')
)
process.patElectrons.isolationValuesNoPFId = cms.PSet(
    pfChargedHadrons = cms.InputTag('elPFIsoValueCharged03NoPFIdPFIso'),
    pfChargedAll = cms.InputTag('elPFIsoValueChargedAll03NoPFIdPFIso'),
    pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU03NoPFIdPFIso'),
    pfNeutralHadrons = cms.InputTag('elPFIsoValueNeutral03NoPFIdPFIso'),
    pfPhotons = cms.InputTag('elPFIsoValueGamma03NoPFIdPFIso')
)

# MVA electron ID
process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.mvaID = cms.Sequence(process.mvaTrigV0 + process.mvaNonTrigV0)
process.patElectrons.electronIDSources = cms.PSet(
    mvaTrigV0 = cms.InputTag('mvaTrigV0'),
    mvaNonTrigV0 = cms.InputTag('mvaNonTrigV0')
)

# calibrate electrons
process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('selectedPatElectrons')
process.RandomNumberGeneratorService = cms.Service('RandomNumberGeneratorService',
    calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
    )
)
process.load('EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi')
process.calibratedPatElectrons.isMC = cms.bool(runOnMC)
if runOnMC: process.calibratedPatElectrons.inputDataset = cms.string('Summer12_DR53X_HCP2012')
else: process.calibratedPatElectrons.inputDataset = cms.string('Moriond2013')
process.calibratedPatElectrons.correctionsType = cms.int32(1)
process.calibratedPatElectrons.combinationType = cms.int32(3)
process.calibratedPatElectrons.lumiRatio = cms.double(0.5955)

# remove MC matching
if not runOnMC:
    from PhysicsTools.PatAlgos.tools.coreTools import *
    removeMCMatching(process, ['All'])

# scraping filter
process.scrapingVeto = cms.EDFilter('FilterOutScraping',
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.2)
)

# PV filter
process.primaryVertexFilter = cms.EDFilter('GoodVertexFilter',
    vertexCollection = cms.InputTag('offlinePrimaryVertices'),
    minimumNDOF = cms.uint32(4),
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)

# HB + HE noise filter
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')

# EE bad super crystal filter
process.load('RecoMET.METFilters.eeBadScFilter_cfi')

# HCAL laser filter
process.load('RecoMET.METFilters.hcalLaserEventFilter_cfi')
process.hcalLaserEventFilter.vetoByRunEventNumber=cms.untracked.bool(False)
process.hcalLaserEventFilter.vetoByHBHEOccupancy=cms.untracked.bool(True)

# trigger filter
process.triggerResultsFilter = cms.EDFilter('TriggerResultsFilter',
    hltResults = cms.InputTag('TriggerResults::HLT'),
    l1tResults = cms.InputTag(''),
    l1tIgnoreMask = cms.bool(False),
    l1techIgnorePrescales = cms.bool(False),
    daqPartitions = cms.uint32(0x01),
    throw = cms.bool(False),
    triggerConditions = cms.vstring(
        'HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
        'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
        'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
        'HLT_Ele23_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_HFT30_v*',
        'HLT_Ele27_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_HFT15_v*',
        'HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v*',
        'HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
        'HLT_Ele27_WP80_v*'
    )
)

# FSR weight
process.fsrWeight = cms.EDProducer('FSRWeightProducer',
    GenTag = cms.untracked.InputTag('genParticles')
)

# PDF weights
process.pdfWeights = cms.EDProducer('PdfWeightProducer',
    FixPOWHEG = cms.untracked.string('CT10.LHgrid'),
    GenTag = cms.untracked.InputTag('genParticles'),
    PdfInfoTag = cms.untracked.InputTag('generator'),
    PdfSetNames = cms.untracked.vstring(
        'MSTW2008nlo68cl.LHgrid',
        'NNPDF20_100.LHgrid'
    )
)

# tree producer
process.treeProducer = cms.EDAnalyzer('TreeProducer',
    isMC = cms.bool(runOnMC),
    outFileName = cms.string(ntupleName),
    fsrWeightTag = cms.InputTag('fsrWeight'),
    pdfWeightTag_ct10 = cms.InputTag('pdfWeights:CT10'),
    pdfWeightTag_mstw = cms.InputTag('pdfWeights:MSTW2008nlo68cl'),
    pdfWeightTag_nnpdf = cms.InputTag('pdfWeights:NNPDF20'),
    genParticleTag = cms.InputTag('genParticles'),
    vertexTag = cms.InputTag('offlinePrimaryVertices'),
    triggerResultsTag = cms.InputTag('TriggerResults::HLT'),
    triggerSummaryTag = cms.InputTag('hltTriggerSummaryAOD::HLT'),
    electronTag = cms.InputTag('calibratedPatElectrons'),
    isoRhoTag = cms.InputTag('kt6PFJets', 'rho'),
    conversionTag = cms.InputTag('allConversions'),
    beamSpotTag = cms.InputTag('offlineBeamSpot'),
    hfElectronTag = cms.InputTag('hfRecoEcalCandidate'),
    hfClusterMapTag = cms.InputTag('hfEMClusters'),
    pfJetTag = cms.InputTag('selectedPatJetsPFlow'),
    pfMetTag = cms.InputTag('patMETsPFlow')
)

# main path
if runOnMC:
    process.p = cms.Path(
        process.primaryVertexFilter *
        process.goodOfflinePrimaryVertices *
        process.mvaID *
        process.patDefaultSequence *
        process.eleRegressionEnergy *
        process.calibratedPatElectrons *
        process.patPF2PATSequencePFlow *
        process.fsrWeight *
        process.pdfWeights *
        process.treeProducer
    )
else:
    process.p = cms.Path(
        process.scrapingVeto *
        process.primaryVertexFilter *
        process.HBHENoiseFilter *
        process.eeBadScFilter *
        process.hcalLaserEventFilter *
        process.triggerResultsFilter *
        process.goodOfflinePrimaryVertices *
        process.mvaID *
        process.patDefaultSequence *
        process.eleRegressionEnergy *
        process.calibratedPatElectrons *
        process.patPF2PATSequencePFlow *
        process.treeProducer
    )

# output path
#process.out.outputCommands = cms.untracked.vstring('keep *')
#process.op = cms.EndPath(process.out)
