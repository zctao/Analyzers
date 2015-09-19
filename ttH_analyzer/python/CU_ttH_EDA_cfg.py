import FWCore.ParameterSet.Config as cms

process = cms.Process("MAOD")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )
process.GlobalTag.globaltag = 'PHYS14_25_V2::All'
#process.GlobalTag.globaltag = 'MCRUN2_74_V7'

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(5000)
)

from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process.ak4PFCHSL1Fastjet = cms.ESProducer(
	'L1FastjetCorrectionESProducer',
	level       = cms.string('L1FastJet'),
	algorithm   = cms.string('AK4PFchs'),
	srcRho      = cms.InputTag( 'fixedGridRhoFastjetAll' )
)

process.ak4PFchsL2Relative = ak4CaloL2Relative.clone( algorithm = 'AK4PFchs' )
process.ak4PFchsL3Absolute = ak4CaloL3Absolute.clone( algorithm = 'AK4PFchs' )

process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
	correctors = cms.vstring(
		'ak4PFCHSL1Fastjet', 
		'ak4PFchsL2Relative', 
		'ak4PFchsL3Absolute'
	)
)

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
        'file:/eos/uscms/store/user/ztao/ttHToTauTau_M125_13TeV_Spring15_miniAOD.root'
        #'file:/eos/uscms/store/user/ztao/ttHToTauTau_M125_13TeV_Spring15_AOD.root'
        
        #'/store/mc/RunIISpring15DR74/ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9_ext3-v1/30000/028D30EF-2530-E511-9C88-002590AC4CEC.root'
        #'root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9_ext1-v1/00000/046ADB1C-9C07-E511-AA40-002590A887F2.root'
	
	#'/store/user/puigh/TTHSync/ttjets_phys14_20bx25_withfatjets_v2.root'
        #'/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/08B36E8F-5E7F-E411-9D5A-002590200AE4.root'
        #'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root'
        #'/store/mc/Spring14miniaod/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/1E4F9BDC-3E1E-E411-A56C-001E67396EAA.root'
        #'/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/004C6DA7-FB03-E411-96BD-0025905A497A.root'
	)
)

process.ttHsyncExercise = cms.EDAnalyzer('CU_ttH_EDA',
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('Outputs/CU_ttH_EDA_output.root')
)


process.p = cms.Path(process.ttHsyncExercise)



