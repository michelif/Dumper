import FWCore.ParameterSet.Config as cms

process = cms.Process("ExREG")
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryExtended2023SHCalReco_cff')

process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = 'PH2_1K_FB_V6::All'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
#    input = cms.untracked.int32(44),
   )

process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC23_patch1/RelValQCD_Pt_80_120_14TeV/GEN-SIM-RECO/PH2_1K_FB_V6_UPG23SHNoTaper-v2/00000/3A4509FB-6E9F-E411-A11B-002590494C62.root')
#                            fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/TP2023SHCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-RECO/SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1/00000/002CD53C-57A7-E411-BB32-002618943901.root')

                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/micheli/shashlikPerformance_2/CMSSW_6_2_0_SLHC23_patch1/src/HLTTest/MyCandidates/test/0683A170-619F-E411-A723-002481E9458E.root')
#                            fileNames = cms.untracked.vstring('file:./ZEE.root')
                            )




process.dumper = cms.EDAnalyzer('dumper',
                                   OutputFileName = cms.string("oot_qcd_multifit.root"),
                                   isData = cms.bool(False),
                                   saveReco = cms.bool(True),
                                   trgSelection = cms.vstring("HLT_Ele27WP85_Gsf_v1", "HLT_Ele20WP60_Ele8_Mass55_v1","HLT_Ele25WP60_SC4_Mass5_v1",)
                                )




process.p = cms.Path(process.dumper)
