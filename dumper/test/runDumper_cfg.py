import FWCore.ParameterSet.Config as cms

# multi5x5 clusters
#from RecoEcal.EgammaClusterProducers.multi5x5ClusteringSequence_cff import *
# preshower sequence for multi5x5 clusters
#from RecoEcal.EgammaClusterProducers.multi5x5PreshowerClusteringSequence_cff import *

process = cms.Process("ExREG")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V4::All', '')

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(-1),
    input = cms.untracked.int32(200),
   )

process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC23_patch1/RelValQCD_Pt_80_120_14TeV/GEN-SIM-RECO/PH2_1K_FB_V6_UPG23SHNoTaper-v2/00000/3A4509FB-6E9F-E411-A11B-002590494C62.root')
#                            fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/relval/CMSSW_6_2_0_SLHC23_patch1/RelValSingleMuPt100Extended/GEN-SIM-RECO/PH2_1K_FB_V6_UPG23SHNoTaper-v2/00000/089F7CA0-609F-E411-BF99-02163E00E791.root')
#                            fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/TP2023SHCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-RECO/SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1/00000/002CD53C-57A7-E411-BB32-002618943901.root')
#                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/m/micheli/shashlikPerformance_2/CMSSW_6_2_0_SLHC23_patch1/src/HLTTest/MyCandidates/test/0683A170-619F-E411-A723-002481E9458E.root')
#                            fileNames = cms.untracked.vstring('/store/relval/CMSSW_6_2_0_SLHC23_patch1/RelValZEE_14TeV/GEN-SIM-RECO/PH2_1K_FB_V6_UPG23SHNoTaper-v2/00000/44D00C34-6B9F-E411-A697-02163E00F518.root')
                            fileNames = cms.untracked.vstring('file:/tmp/micheli/089F7CA0-609F-E411-BF99-02163E00E791.root')
                            )

process.load("RecoEcal.EgammaClusterProducers.multi5x5ClusteringSequence_cff")
process.multi5x5BasicClustersCleaned.endcapHitCollection = cms.string("EcalRecHitsEK")
process.multi5x5BasicClustersUncleaned.endcapHitCollection = cms.string("EcalRecHitsEK")

process.dumper = cms.EDAnalyzer('dumper',
                                   OutputFileName = cms.string("outDumper.root"),
                                   isData = cms.bool(False),
                                   saveReco = cms.bool(True),
                                   recHitNEvts= cms.int32(101),
                                   trgSelection = cms.vstring("HLT_Ele27WP85_Gsf_v1", "HLT_Ele20WP60_Ele8_Mass55_v1","HLT_Ele25WP60_SC4_Mass5_v1",)
                                )

process.p = cms.Path(process.multi5x5ClusteringSequence * process.dumper)
