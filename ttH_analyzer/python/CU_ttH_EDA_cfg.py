import FWCore.ParameterSet.Config as cms

process = cms.Process("MAOD")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load( "Configuration.StandardSequences.FrontierConditions_GlobalTag_cff" )
#process.GlobalTag.globaltag = 'PHYS14_25_V2::All'
process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(-1)
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
        #'file:/eos/uscms/store/user/ztao/ttHToTauTau_M125_13TeV_Spring15_miniAOD.root'
        #'file:/eos/uscms/store/user/ztao/ttHToTauTau_M125_13TeV_Spring15_AOD.root'
        'file:/eos/uscms/store/user/ztao/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_miniAOD.root'

        #'/store/mc/RunIISpring15DR74/ttHToTT_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/20000/0249345F-032D-E511-A21D-0025905C95F8.root'
        #/store/mc/RunIISpring15DR74/ttHToTT_M125_13TeV_powheg_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v2/20000/008B5DD7-EB2B-E511-8E73-0025904C63F8.root
        #'/store/mc/RunIISpring15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v3/60000/00181849-176A-E511-8B11-848F69FD4C94.root'

        #'/store/mc/RunIISpring15DR74/ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9_ext3-v1/30000/028D30EF-2530-E511-9C88-002590AC4CEC.root'
        #'root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/ttHJetToTT_M125_13TeV_amcatnloFXFX_madspin_pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9_ext1-v1/00000/046ADB1C-9C07-E511-AA40-002590A887F2.root'
	)
)

process.ttHsyncExercise = cms.EDAnalyzer('CU_ttH_EDA',
        # Analysis type choices: 'lepton+jet', 'dilepton', 'taus_lepton+jet', 'taus_dilepton'
        analysis_type = cms.string("taus_dilepton"),
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
        min_njets = cms.int32(0),
        min_nbtags = cms.int32(0),
        # Jets
        jet_corrector = cms.string('ak4PFchsL1L2L3'),
        # MiniAODhelper
        using_real_data = cms.bool(False),
        ## available choices '-': none, 'L': loose, 'M': medium, 'T': tight
        b_tag_strength = cms.string('M')
)

process.TFileService = cms.Service("TFileService",
	#fileName = cms.string('Outputs/CU_ttH_EDA_output_sig.root')
        #fileName = cms.string('Outputs/CU_ttH_EDA_output_TTJets.root')
        fileName = cms.string('CU_ttH_EDA_output.root')
)


process.p = cms.Path(process.ttHsyncExercise)



