import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.RecoJets_cff import *
from RecoJets.Configuration.RecoPFJets_cff import *
from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process = cms.Process("MAOD")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
)

process.ak4PFCHSL1Fastjet = cms.ESProducer(
	'L1FastjetCorrectionESProducer',
	level       = cms.string('L1FastJet'),
	algorithm   = cms.string('AK4PFchs'),
	#srcRho      = cms.InputTag( 'fixedGridRhoFastjetAll' )
        srcRho      = cms.InputTag( 'fixedGridRhoFastjetCentralNeutral' ),
        useCondDB   = cms.untracked.bool(True)
)

#process.ak4PFchsL2Relative = ak4CaloL2Relative.clone( algorithm = 'AK4PFchs' )
process.ak4PFchsL2Relative = ak5PFL2Relative.clone( algorithm = 'AK4PFchs' )
#process.ak4PFchsL3Absolute = ak4CaloL3Absolute.clone( algorithm = 'AK4PFchs' )
process.ak4PFchsL3Absolute = ak5PFL3Absolute.clone( algorithm = 'AK4PFchs' )

process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
	correctors = cms.vstring(
		'ak4PFCHSL1Fastjet', 
		'ak4PFchsL2Relative', 
		'ak4PFchsL3Absolute'),
        useCondDB = cms.untracked.bool(True)
)

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
        # signal sample
        '/store/mc/RunIIFall15MiniAODv2/ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3/60000/0C6DA13E-38C8-E511-8F6E-00259055220A.root'
        # tt+jet 
        #'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/00DF0A73-17C2-E511-B086-E41D2D08DE30.root'
	)
)

# LeptonID producer from ttH Multi-lepton group
process.load("ttH.LeptonID.ttHLeptons_cfi")
# new electron MVA developed by the EGamma POG 
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
# load the analysis:
process.load("Analyzers.ttH_analyzer.ttHtaus_cfi")

# re-define parameter sets here if necessary
process.ttHLeptons.rhoParam = "fixedGridRhoFastjetCentralNeutral"
# use leptons from LeptonID producer
process.ttHtaus.input_tags.electrons = cms.InputTag("ttHLeptons")
process.ttHtaus.input_tags.muons = cms.InputTag("ttHLeptons")
process.ttHtaus.input_tags.taus = cms.InputTag("ttHLeptons")

    
process.TFileService = cms.Service("TFileService",
	fileName = cms.string('ttHtausNtuple.root')
        #fileName = cms.string('ttHtausNtuple_ttJets.root')
)

process.p = cms.Path(
    process.electronMVAValueMapProducer
    * process.ttHLeptons
    * process.ttHtaus
)
