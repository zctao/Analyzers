import FWCore.ParameterSet.Config as cms

process = cms.Process("MAOD")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )
process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_v12'

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
        '/store/mc/RunIIFall15MiniAODv2/ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v3/60000/0C6DA13E-38C8-E511-8F6E-00259055220A.root'
	)
)

process.ttHtaus = cms.EDAnalyzer('CU_ttH_EDA',
        # Analysis type choice: 'lepton+jet', 'dilepton', 'tau_ssleptons', 'ditaus_lepton'
        analysis_type = cms.string("tau_ssleptons"),
        # Generic
        verbosity = cms.bool(False),
        print_HLT_event_path = cms.bool(False),
        HLT_config_tag = cms.string('HLT'),
        filter_config_tag = cms.string('PAT'),
        # Triggers
        collect_trigger_stats = cms.bool(False),
        ## Single lepton triggers:
        HLT_electron_triggers = cms.vstring([
            'HLT_Ele27_eta2p1_WP75_Gsf_v1'
        ]),
        HLT_muon_triggers = cms.vstring([
            'HLT_IsoMu24_eta2p1_v1'
        ]),
        ## Dilepton triggers:
        HLT_electron_electron_triggers = cms.vstring([
            'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v1'
        ]),
        HLT_electron_muon_triggers = cms.vstring([
            'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v1',
            'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v1'
        ]),
        HLT_muon_muon_triggers = cms.vstring([
            'HLT_Mu30_TkMu11_v1',
            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1',
            'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1'
        ]),
        # Cuts
        min_tight_lepton_pT = cms.double(20),
        min_tau_pT = cms.double(20),
        min_jet_pT = cms.double(30),
        min_bjet_pT = cms.double(20),
        max_jet_eta = cms.double(2.5),
        max_bjet_eta = cms.double(2.5),
        min_njets = cms.int32(2),
        min_nbtags = cms.int32(1),
        # Jets
        jet_corrector = cms.string('ak4PFchsL1L2L3'),
        # MiniAODhelper
        using_real_data = cms.bool(False),
        ## available choices '-': none, 'L': loose, 'M': medium, 'T': tight
        b_tag_strength = cms.string('M')
)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('CU_ttH_EDA_output.root')
)

process.p = cms.Path(process.ttHtaus)
