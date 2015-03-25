# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023SHCal --conditions PH2_1K_FB_V4::All --step RECO --magField 38T_PostLS1 --geometry Extended2023SHCalNoTaper,Extended2023SHCalNoTaperReco --python_filename reReco.py --no_exec -n 1
import FWCore.ParameterSet.Config as cms

process = cms.Process('RERECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
                            fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/TP2023SHCALDR/DYToEE_M-20_TuneZ2star_14TeV-pythia6-tauola/GEN-SIM-RECO/SHCALJan23_NoPU_PH2_1K_FB_V6-v1/10000/0A43CD37-29C6-E411-A16F-0025905B8598.root')
#    fileNames = cms.untracked.vstring('/store/relval/CMSSW_6_2_0_SLHC23_patch1/RelValDYToLL_M_50_TuneZ2star_14TeV/GEN-SIM-RECO/PH2_1K_FB_V6_UPG23SHNoTaper-v1/00000/0046E050-2AA0-E411-9AD1-02163E00E959.root')
#    fileNames = cms.untracked.vstring('/store/relval/CMSSW_6_2_0_SLHC23_patch1/RelValDYToLL_M_50_TuneZ2star_14TeV/GEN-SIM-RECO/PU_PH2_1K_FB_V6_SHNoTapPU140-v1/00000/02D81807-B6A0-E411-9B26-02163E00EB60.root')
)

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('--customise nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition
#process.FEVTDEBUGEventContent.outputCommands.append('drop *_*_*_RECO')
#process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
#    fileName = cms.untracked.string('reRECO.root'),
#    dataset = cms.untracked.PSet(
#        filterName = cms.untracked.string(''),
#        dataTier = cms.untracked.string('')
#    )
#)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '')

from SLHCUpgradeSimulations.Configuration.combinedCustoms import *

process.reconstruction_step = cms.Path()
process=cust_2023Muon(process)
process=customise_shashlikElectronHOverE(process)
process.ecalRecHit.EEuncalibRecHitCollection = cms.InputTag("","")
process.reducedEcalRecHitsSequence.remove(process.reducedEcalRecHitsES)
        #remove the old EE pfrechit producer
del process.particleFlowRecHitECAL.producers[1]
process.particleFlowClusterEBEKMerger = cms.EDProducer('PFClusterCollectionMerger',
                                                       inputs = cms.VInputTag(cms.InputTag('particleFlowClusterECALUncorrected'),
                                                                              cms.InputTag('particleFlowClusterEKUncorrected')
                                                                              )
                                                       )   
process.pfClusteringECAL.remove(process.particleFlowClusterECAL)
#        process.pfClusteringECAL.remove(process.particleFlowClusterECALWithTimeSelected)
process.pfClusteringECAL += process.pfClusteringEK 
process.pfClusteringECAL += process.particleFlowClusterEBEKMerger
process.pfClusteringECAL += process.particleFlowClusterECAL        
process.particleFlowClusterECAL.inputECAL = cms.InputTag('particleFlowClusterEBEKMerger')
process.pfClusteringECAL += process.particleFlowClusterECAL
        #process.particleFlowCluster += process.pfClusteringEK

        #clone photons to mustache photons so we can compare back to old reco
process.mustachePhotonCore = process.photonCore.clone(scHybridBarrelProducer = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel"),scIslandEndcapProducer = cms.InputTag("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower"))
process.mustachePhotons = process.photons.clone(photonCoreProducer = cms.InputTag('mustachePhotonCore'), endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEK"))        
process.photonSequence += process.mustachePhotonCore
process.photonSequence += process.mustachePhotons
        #point particle flow at the right photon collection     
process.particleFlowBlock.elementImporters[2].source = cms.InputTag('mustachePhotons')
process.gedPhotons.endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEK")

process.towerMaker.ecalInputs = cms.VInputTag(cms.InputTag("ecalRecHit","EcalRecHitsEB"), cms.InputTag("ecalRecHit","EcalRecHitsEK"))
process.towerMakerPF.ecalInputs = cms.VInputTag(cms.InputTag("ecalRecHit","EcalRecHitsEB"), cms.InputTag("ecalRecHit","EcalRecHitsEK"))
process.towerMakerWithHO.ecalInputs = cms.VInputTag(cms.InputTag("ecalRecHit","EcalRecHitsEB"), cms.InputTag("ecalRecHit","EcalRecHitsEK"))
process.towerMaker.EESumThreshold = cms.double(0.1)
process.towerMakerPF.EESumThreshold = cms.double(0.1)
process.towerMakerWithHO.EESumThreshold = cms.double(0.1)
process.towerMaker.EEThreshold = cms.double(0.035)
process.towerMakerPF.EEThreshold = cms.double(0.035)
process.towerMakerWithHO.EEThreshold = cms.double(0.035)

# Change all processes to use EcalRecHitsEK instead of EcalRecHitsEE
process.EcalHaloData.EERecHitLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.JPTeidTight.reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.calomuons.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.conversionTrackCandidates.endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.earlyMuons.CaloExtractorPSet.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.earlyMuons.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.earlyMuons.JetExtractorPSet.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.ecalDrivenGsfElectrons.endcapRecHitCollectionTag = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.eidLoose.reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.eidRobustHighEnergy.reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.eidRobustLoose.reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.eidRobustTight.reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.eidTight.reducedEndcapRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.gedGsfElectrons.endcapRecHitCollectionTag = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.gedPhotons.mipVariableSet.endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.gedPhotons.isolationSumsCalculatorSet.endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.gsfElectrons.endcapRecHitCollectionTag = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.interestingEleIsoDetIdEE.recHitsLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.interestingGamIsoDetIdEE.recHitsLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.muonMETValueMapProducer.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.muons1stStep.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.muons1stStep.JetExtractorPSet.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.muons1stStep.CaloExtractorPSet.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.muonsFromCosmics.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.muonsFromCosmics.JetExtractorPSet.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.muonsFromCosmics.CaloExtractorPSet.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.muonsFromCosmics1Leg.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.muonsFromCosmics1Leg.JetExtractorPSet.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.muonsFromCosmics1Leg.JetExtractorPSet.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.muonsFromCosmics1Leg.CaloExtractorPSet.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.particleFlowSuperClusterECAL.regressionConfig.ecalRecHitsEE = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.pfElectronInterestingEcalDetIdEE.recHitsLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.pfPhotonInterestingEcalDetIdEE.recHitsLabel = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.pfPhotonTranslator.endcapEcalHits = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.pfPhotonTranslator.EGPhotons = cms.string("mustachePhotons")
process.photons.mipVariableSet.endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.photons.isolationSumsCalculatorSet.endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.uncleanedOnlyConversionTrackCandidates.endcapEcalRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEK")
process.uncleanedOnlyGsfElectrons.endcapRecHitCollectionTag = cms.InputTag("ecalRecHit","EcalRecHitsEK")

# Path and EndPath definitions
process.reReco=cms.Sequence(process.trackerlocalreco*process.muonGlobalReco*process.trackingGlobalReco*process.particleFlowCluster*process.ecalClusters*process.egammaGlobalReco*process.pfTrackingGlobalReco*process.egammaHighLevelRecoPrePF*process.particleFlowReco*process.egammaHighLevelRecoPostPF)

process.reconstruction_step *= process.reReco
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

process.dumper = cms.EDAnalyzer('dumper',
                                   OutputFileName = cms.string("outDumper.root"),
                                   isData = cms.bool(False),
                                   saveReco = cms.bool(True),
                                   recHitNEvts= cms.int32(101),
                                   trgSelection = cms.vstring("HLT_Ele27WP85_Gsf_v1", "HLT_Ele20WP60_Ele8_Mass55_v1","HLT_Ele25WP60_SC4_Mass5_v1",)
                                )
process.dumper_step = cms.EndPath(process.dumper)



# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.dumper_step)


