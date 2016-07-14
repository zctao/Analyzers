### CMSSW_7_6_3

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import HLTrigger.HLTfilters.hltHighLevel_cfi

process = cms.Process("ttH")

### initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.option = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

### Run in unscheduled mode
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

### 
options = VarParsing.VarParsing('analysis')

options.register('doSync', False, 
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Do synchronization or not")
options.register('isData', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run on real data or not")
options.register('doScale', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Scale histogram or not")
options.register('IntLumi', 1.,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Integrated luminosity of the sample")
options.register('SampleName','',  #'ttH_htt'
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Sample name") 
options.register('CrossSection', 1.,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Cross section of the channel")
options.register('doSystematics', True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Include systematics or not")
options.register('OutputDir', '',  #'/uscms/home/ztao/nobackup/'
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "output file directory")

options.maxEvents = -1

# get and parse the command line arguments
options.parseArguments()

### Global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )

if options.isData:
    process.GlobalTag.globaltag = '76X_dataRun2_v15'
else:
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(options.maxEvents)
)

### Filters for running on collision data
if options.isData:
    # primary vertex filter
    process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                               vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
                                               minimumNDOF = cms.uint32(4) ,
                                               maxAbsZ = cms.double(24), 
                                               maxd0 = cms.double(2) 
    )

### Inputs
#'ttH_htt' signal sample
#'/store/mc/RunIIFall15MiniAODv2/ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3/60000/0C6DA13E-38C8-E511-8F6E-00259055220A.root'
#'file:/uscms/home/ztao/nobackup/data/ttH_76X/ttH_htt_sync.root'
#'ttbar'  tt+jet
#'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/00DF0A73-17C2-E511-B086-E41D2D08DE30.root'
#'file:../data/local_tmp/ttbar.root'

options.inputFiles='file:/uscms/home/ztao/nobackup/data/ttH_76X/ttH_htt.root'
    
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(options.inputFiles)
)

### HLT Filter
process.HLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    throw = cms.bool(False),
    HLTPaths = [
        # single lepton trigger
        'HLT_Ele23_WPLoose_Gsf_v*',
        'HLT_IsoMu20_v*',
        'HLT_IsoTkMu20_v*',
        # dilepton trigger
        'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*',
        'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*',
        'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*',
        'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*'
        ]
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

process.ttHtaus.do_systematics = cms.bool(options.doSystematics)
process.ttHtaus.produce_sync_ntuple = cms.bool(options.doSync)
process.ttHtaus.doScale = cms.bool(options.doScale)
process.ttHtaus.sample_xs = cms.double(options.CrossSection)
process.ttHtaus.int_lumi = cms.double(options.IntLumi)
process.ttHtaus.using_real_data = cms.bool(options.isData)

### Outputs
if options.isData:
    out_file = options.OutputDir + 'output_data.root'
else:
    out_file = options.OutputDir + 'output_' + options.SampleName + '.root'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(out_file)
)

#Path
if options.isData:
    process.p = cms.Path(
        process.primaryVertexFilter *
        process.HLTFilter *
        process.patJetCorrFactorsReapplyJEC *
        process.patJetsReapplyJEC *
        process.electronMVAValueMapProducer *
        process.ttHLeptons *
        process.ttHtaus
    )
else:
    process.p = cms.Path(
        process.HLTFilter *
        process.patJetCorrFactorsReapplyJEC * process.patJetsReapplyJEC *
        process.electronMVAValueMapProducer *
        process.ttHLeptons *
        process.ttHtaus
    )


    
