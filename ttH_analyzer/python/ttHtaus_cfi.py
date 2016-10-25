# Configurations for ttH tautau analysis
#
# Zhengcheng Tao
#

import FWCore.ParameterSet.Config as cms

ttHtaus =  cms.EDAnalyzer('CU_ttH_EDA',
        # Analysis type choice: '2lss1tau', '3leptons'
        analysis_type = cms.string("2lss1tau"),
        # Sync ntuple
        turn_off_event_sel = cms.bool(False),
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
            'HLT_Ele27_eta2p1_WPLoose_Gsf_v3'
        ]),
        HLT_muon_triggers = cms.vstring([
            'HLT_IsoMu22_v',
            'HLT_IsoTkMu22_v'
        ]),
        ## Dilepton triggers:
        HLT_electron_electron_triggers = cms.vstring([
            'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v'
        ]),
        HLT_electron_muon_triggers = cms.vstring([
            'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v',
            'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v'
        ]),
        HLT_muon_muon_triggers = cms.vstring([
            'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v',
            'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v'
        ]),
                          
        using_real_data = cms.bool(False),
        # JEC
        JECType = cms.string('JECType'),
        #jet_corrector = cms.string('ak4PFchsL1L2L3'),
        selection_region = cms.string('signal'),
        # debug flag
        debug_mode = cms.bool(False),
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
