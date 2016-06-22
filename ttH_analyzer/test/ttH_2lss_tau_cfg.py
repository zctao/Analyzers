### CMSSW_7_6_3

import FWCore.ParameterSet.Config as cms

process = cms.Process("ttH")

### initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

### Run in unscheduled mode
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

### Switches
isData = False
isSignal = True

### Global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )

if isData:
    process.GlobalTag.globaltag = '76X_dataRun2_v15'
else:
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
)

### Inputs
if isData:
    input_file = ''
elif isSignal:
    # signal sample
    input_file = '/store/mc/RunIIFall15MiniAODv2/ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3/60000/0C6DA13E-38C8-E511-8F6E-00259055220A.root'
else:
    # tt+jet
    input_file = '/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/00DF0A73-17C2-E511-B086-E41D2D08DE30.root'

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(input_file)
)

### JEC
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated #updatedPatJetCorrFactors
process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(  #updatedPatJetCorrFactors.clone(
  src = cms.InputTag("slimmedJets"),
  levels = ['L1FastJet', 
        'L2Relative', 
        'L3Absolute',
        'L2L3Residual'],
  payload = 'AK4PFchs' ) 

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated #updatedPatJets
process.patJetsReapplyJEC = patJetsUpdated.clone( #updatedPatJets.clone(
  jetSource = cms.InputTag("slimmedJets"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
  )

# Get jetCorrector in case needed
#from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
#process.ak4PFCHSL1Fastjet = cms.ESProducer(
#        'L1FastjetCorrectionESProducer',
#        level       = cms.string('L1FastJet'),
#        algorithm   = cms.string('AK4PFchs'),
#        srcRho      = cms.InputTag('fixedGridRhoFastjetAll')
#)
#process.ak4PFchsL2Relative = ak4CaloL2Relative.clone(algorithm = 'AK4PFchs')
#process.ak4PFchsL3Absolute = ak4CaloL3Absolute.clone(algorithm = 'AK4PFchs')
#process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
#        correctors = cms.vstring(
#            'ak4PFCHSL1Fastjet',
#            'ak4PFchsL2Relative',
#            'ak4PFchsL3Absolute')
#)
#
#if isData:
#    process.ak4PFchsResidual  = ak4CaloResidual.clone(algorithm = 'AK4PFchs')
#    process.ak4PFchsL1L2L3.correctors.append('ak4PFchsResidual')

### load the analysis
# electron MVA developed by the EGamma POG
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
# LeptonID producer from ttH Multi-lepton group
process.load("ttH.LeptonID.ttHLeptons_cfi")
# Analyzer
process.load("Analyzers.ttH_analyzer.ttHtaus_cfi")


### Redefine parameter sets
process.ttHLeptons.rhoParam = "fixedGridRhoFastjetCentralNeutral"
process.ttHLeptons.jets = cms.InputTag("patJetsReapplyJEC") # use JEC's from tag
process.ttHtaus.input_tags.electrons = cms.InputTag("ttHLeptons")
process.ttHtaus.input_tags.muons = cms.InputTag("ttHLeptons")
process.ttHtaus.input_tags.taus = cms.InputTag("ttHLeptons")
process.ttHtaus.input_tags.jets = cms.InputTag("patJetsReapplyJEC")
process.ttHtaus.input_tags.rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral")

### Outputs
if isSignal:
    out_file = 'testNtuple.root'
else:
    out_file = 'testNtuple_bkg.root'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(out_file)

)

#Path
process.p = cms.Path(
    process.patJetCorrFactorsReapplyJEC * process.patJetsReapplyJEC *
    process.electronMVAValueMapProducer *
    process.ttHLeptons *
    process.ttHtaus
)
