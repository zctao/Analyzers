import FWCore.ParameterSet.Config as cms

# Give the process a name
process = cms.Process("PickEvent")

# Tell the process which files to use as the sourdce
process.source = cms.Source ("PoolSource",
          fileNames = cms.untracked.vstring ('/store/mc/RunIISpring15DR74/ttHToTT_M125_13TeV_powheg_pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v2/20000/0249345F-032D-E511-A21D-0025905C95F8.root')
)

# tell the process to only run over 100 events (-1 would mean run over
#  everything
process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32 (100)

)

# Tell the process what filename to use to save the output
process.Out = cms.OutputModule("PoolOutputModule",
         fileName = cms.untracked.string ("/uscms/home/ztao/work/CU_ttH_WD/Outputs/SignalSample100.root")
)

# make sure everything is hooked up
process.end = cms.EndPath(process.Out)
