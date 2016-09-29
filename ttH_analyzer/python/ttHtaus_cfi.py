# Configurations for ttH tautau analysis
#
# Zhengcheng Tao
#

import FWCore.ParameterSet.Config as cms

ttHtaus =  cms.EDAnalyzer('CU_ttH_EDA',
        # Analysis type choice: '2lss1tau', '3leptons'
        analysis_type = cms.string("2lss1tau"),
        # Sync ntuple
        produce_sync_ntuple = cms.bool(False),
        # Systematics
        do_systematics = cms.bool(False),
        # Sample parameter
        doLumiScale = cms.bool(False),
        sampleName = cms.string(''),
        sample_xs = cms.double(1.),
        int_lumi = cms.double(1.),
        # Generic
        verbosity = cms.bool(False),
        # Triggers
        print_HLT_event_path = cms.bool(False),
        turn_off_HLT_cut = cms.bool(False),
        HLT_config_tag = cms.string('HLT'),
        filter_config_tag = cms.string('HLT'),
        collect_trigger_stats = cms.bool(False),
        ## Single lepton triggers:
        HLT_electron_triggers = cms.vstring([
            'HLT_Ele27_eta2p1_WPLoose_Gsf_v7'
        ]),
        HLT_muon_triggers = cms.vstring([
            'HLT_IsoMu22_v3',
            'HLT_IsoTkMu22_v4'
        ]),
        ## Dilepton triggers:
        HLT_electron_electron_triggers = cms.vstring([
            'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v8'
        ]),
        HLT_electron_muon_triggers = cms.vstring([
            'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v3',
            'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v3'
        ]),
        HLT_muon_muon_triggers = cms.vstring([
            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v4',
            'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v3'
        ]),
        # Cuts
        min_tight_lepton_pT = cms.double(20),
        min_ele_pT = cms.double(7.),
        min_mu_pT = cms.double(5.),
        min_tau_pT = cms.double(20.),
        min_jet_pT = cms.double(25.),
        min_bjet_pT = cms.double(20.),
        max_jet_eta = cms.double(2.4),
        max_bjet_eta = cms.double(2.5),
        min_njets = cms.int32(2),
        min_nbtags = cms.int32(1),
        # JEC
        #jet_corrector = cms.string('ak4PFchsL1L2L3'),
        using_real_data = cms.bool(False),
        selection_region = cms.string('signal'),
        # InputTags
        input_tags = cms.PSet(
            pv = cms.InputTag("offlineSlimmedPrimaryVertices"),
            sv = cms.InputTag("slimmedSecondaryVertices"),
            pileup = cms.InputTag("addPileupInfo"),
            rho = cms.InputTag("fixedGridRhoFastjetAll"),
            electrons = cms.InputTag("slimmedElectrons"),
            muons = cms.InputTag("slimmedMuons"),
            taus = cms.InputTag("slimmedTaus"),
            jets = cms.InputTag("slimmedJets"),
            mets = cms.InputTag("slimmedMETs"),
            pfcand = cms.InputTag("packedPFCandidates"),
            beamspot = cms.InputTag("offlineBeamSpot"),
            packedgen = cms.InputTag("packedGenParticles"),
            prunedgen = cms.InputTag("prunedGenParticles")
        )
)
