import FWCore.ParameterSet.Config as cms

process = cms.Process('ttH')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring(
    'file:/eos/uscms/store/user/ztao/ttHToTauTau_M125_13TeV_Spring15_AOD.root'
    )
)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.ParticleListDrawer = cms.EDAnalyzer('ParticleListDrawer',
                                            maxEventsToPrint = cms.untracked.int32(5),
                                            src = cms.InputTag('genParticles'),
                                            printOnlyHardInteraction = cms.untracked.bool(False),
                                            useMessageLogger = cms.untracked.bool(False)
                                            )
process.p = cms.Path(process.ParticleListDrawer)
