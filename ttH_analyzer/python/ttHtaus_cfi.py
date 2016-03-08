# Configurations for ttH tautau analysis
#
# Zhengcheng Tao
#

import FWCore.ParameterSet.Config as cms

ttHtaus =  cms.EDAnalyzer('CU_ttH_EDA',
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
        b_tag_strength = cms.string('M'),
                          
        # InputTags
        input_tags = cms.PSet(
            pv = cms.InputTag("offlineSlimmedPrimaryVertices"),
            sv = cms.InputTag("slimmedSecondaryVertices"),
            pileup = cms.InputTag("addPileupInfo"),
            rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
            #"fixedGridRhoFastjetAll"
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
