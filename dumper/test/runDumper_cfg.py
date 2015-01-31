import FWCore.ParameterSet.Config as cms

process = cms.Process("ExREG")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'START53_V10::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    )

process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC23_patch1/RelValQCD_Pt_80_120_14TeV/GEN-SIM-RECO/PH2_1K_FB_V6_UPG23SHNoTaper-v2/00000/3A4509FB-6E9F-E411-A11B-002590494C62.root')
                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/micheli/shashlikPerformance_2/CMSSW_6_2_0_SLHC23_patch1/src/HLTTest/MyCandidates/test/0683A170-619F-E411-A723-002481E9458E.root')
                            )




process.dumper = cms.EDAnalyzer('dumper',
                                   OutputFileName = cms.string("oot_qcd_multifit.root"),
                                   hitEBLabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEB"),
                                   hitEELabel = cms.InputTag("hltEcalRecHit","EcalRecHitsEE"),
                                   isData = cms.bool(False),
                                   saveReco = cms.bool(True),
                                   trgSelection = cms.vstring("HLT_Ele27WP85_Gsf_v1", "HLT_Ele20WP60_Ele8_Mass55_v1","HLT_Ele25WP60_SC4_Mass5_v1",)
                                )




process.p = cms.Path(process.dumper)
