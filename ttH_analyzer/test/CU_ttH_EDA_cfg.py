### CMSSW_8_0_16

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

options.register('evtSelOff', False, 
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Turn off event selection")
options.register('isData', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run on real data or not")
options.register('doLumiScale', False,
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
options.register('OutputDir', '/uscms/home/ztao/nobackup/',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "output file directory")
options.register('SelectionRegion', 'signal_2lss1tau',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Which selection region to apply: signal_2lss1tau, control_2los1tau, control_1lfakeable, control_WZ")
options.register('TurnOffHLTCut', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Turn off HLT path check in event selection")
options.register('JECType', 'NA',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "JEC type: NA, JESUp, JESDown, JERUp, JERDown")
options.register('Debug', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Debug mode")
options.register('reHLT', False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "if reHLT use tag HLT2")
options.register('GridMode', True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "run on grid or interactively")
options.register("doJERSmearing", True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "apply jet energy smearing for MC or not")

options.maxEvents = -1
options.inputFiles='file:/uscms/home/ztao/nobackup/datasample/ttH_80X/ttHnonbb.root'

# get and parse the command line arguments
options.parseArguments()

### Global tag
process.load('Configuration.StandardSequences.Services_cff')
process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )

if options.isData:
    process.GlobalTag.globaltag = '80X_dataRun2_Prompt_ICHEP16JEC_v0'
else:
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'

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
##'ttbar'  tt+jet

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(options.inputFiles)
)

### JEC
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    labelName = 'UpdatedJEC',
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual']), 'None')
)

### MET Uncertainty
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#runMetCorAndUncFromMiniAOD(process,
#                           isData=options.isData,
#                           )

### load the analysis
# electron MVA developed by the EGamma POG
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
# LeptonID producer from ttH Multi-lepton group
process.load("ttH.LeptonID.ttHLeptons_cfi")
# Analyzer
process.load("Analyzers.ttH_analyzer.ttHtaus_cfi")


### Redefine parameter sets
process.ttHLeptons.rhoParam = "fixedGridRhoFastjetCentralNeutral"
process.ttHLeptons.jets = cms.InputTag("updatedPatJetsUpdatedJEC") # use JEC's from tag

process.ttHtaus.input_tags.electrons = cms.InputTag("ttHLeptons")
process.ttHtaus.input_tags.muons = cms.InputTag("ttHLeptons")
process.ttHtaus.input_tags.taus = cms.InputTag("ttHLeptons")
process.ttHtaus.input_tags.jets = cms.InputTag("updatedPatJetsUpdatedJEC")
#process.ttHtaus.input_tags.mets = cms.InputTag("slimmedMETs","","ttH")
process.ttHtaus.input_tags.rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral")

process.ttHtaus.do_systematics = cms.bool(options.doSystematics)
process.ttHtaus.turn_off_event_sel = cms.bool(options.evtSelOff)
process.ttHtaus.doLumiScale = cms.bool(options.doLumiScale)
process.ttHtaus.sampleName = cms.string(options.SampleName)
process.ttHtaus.sample_xs = cms.double(options.CrossSection)
process.ttHtaus.int_lumi = cms.double(options.IntLumi)
process.ttHtaus.JECType = cms.string(options.JECType)
process.ttHtaus.using_real_data = cms.bool(options.isData)
process.ttHtaus.selection_region = cms.string(options.SelectionRegion)
process.ttHtaus.turn_off_HLT_cut = cms.bool(options.TurnOffHLTCut)
process.ttHtaus.debug_mode = cms.bool(options.Debug)
process.ttHtaus.doJERsmear = cms.bool(options.doJERSmearing)
# for reHLT
if options.reHLT:
    process.ttHtaus.HLT_config_tag = cms.string("HLT2")
    process.ttHtaus.filter_config_tag = cms.string("HLT2")

### Outputs
if options.GridMode:
    out_file = 'output_' + options.SampleName + '.root'
else:
    out_file = options.OutputDir + 'output_' + options.SampleName + '.root'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(out_file)
)

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#                                        ignoreTotal = cms.untracked.int32(1)
#                                        )

#Path
if options.isData:
    process.p = cms.Path(
        process.primaryVertexFilter *
        process.patJetCorrFactorsUpdatedJEC *
        process.updatedPatJetsUpdatedJEC *
        process.electronMVAValueMapProducer *
        #process.fullPatMetSequence *
        process.ttHLeptons *
        process.ttHtaus
    )
else:
    process.p = cms.Path(
        process.patJetCorrFactorsUpdatedJEC * process.updatedPatJetsUpdatedJEC *
        process.electronMVAValueMapProducer *
        #process.fullPatMetSequence *
        process.ttHLeptons *
        process.ttHtaus
    )


    
