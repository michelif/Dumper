// -*- C++ -*-
//
// Package:    Dumper
// Class:      dumper
//
// Original Author:  Francesco Micheli
//         Created:  Sat Jan 31 11:26:52 CET 2015
// $Id$
//
//



#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaHLTAlgos/interface/EgammaHLTHcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"
#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

//#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
//#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/FwdPtr.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h" 
#include "RecoCaloTools/Selectors/interface/CaloConeSelector.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"

//#include "RecoEcal/EgammaCoreTools/interface/ClusterShapeAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "RecoEgamma/Examples/plugins/MCElectronAnalyzer.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruthFinder.h"
#include "RecoEgamma/EgammaMCTools/interface/ElectronMCTruth.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Dumper/dumper/interface/EcalClusterTools.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <string>
#include <iostream>
#define MAXPARTICLESTOSAVE 100
#define MAXPHOTONSTOSAVE 20
#define MAXSCTOSAVE 200
#define MAXBCTOSAVE 200
#define MAXPFCANDTOSAVE 200
#define MAXRECHITTOSAVE 20000

//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array

class dumper : public edm::EDAnalyzer {
public:
  explicit dumper(const edm::ParameterSet&);
  ~dumper();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void clearVector();
  void mcTruth(edm::Handle<reco::GenParticleCollection> genParticleH,std::vector<ElectronMCTruth> MCElectrons);
  void phoReco(edm::Handle<reco::PhotonCollection> photonH, edm::Handle<reco::TrackCollection> traH, const EBRecHitCollection* rhitseb,const EERecHitCollection* rhitsee,edm::Handle<reco::PFCandidateCollection>  PFCandidates);
  void eleReco(edm::Handle<reco::GsfElectronCollection> ElectronHandle,edm::Handle<reco::ConversionCollection> hConversions, edm::Handle<reco::BeamSpot> recoBeamSpotHandle, const EBRecHitCollection* rhitseb,const EERecHitCollection* rhitsee,edm::Handle<reco::PFCandidateCollection>  PFCandidates);
  void scReco(edm::Handle<reco::SuperClusterCollection> superClustersEBHandle, edm::Handle<reco::SuperClusterCollection> superClustersEEHandle);
  void multi5x5scReco(edm::Handle<reco::SuperClusterCollection> multi5x5Handle);
  void hybridscReco(edm::Handle<reco::SuperClusterCollection> hybridHandle);
  void recHitReco(const EERecHitCollection* rhitsee);
private:
  virtual void beginRun(const edm::Run& run,const edm::EventSetup& setup);
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;


  std::string outputFileName;
  TFile* f;
  TTree* t;

  const CaloGeometry* geometry;
  const CaloTopology *topology;
  
  Float_t rho;
  Int_t n, npf, gp_n, gpho_n, gele_n, pho_n,ele_n, pfSC_n, multi5x5SC_n, hybridSC_n, pfcandPho_n, pfcandEle_n, rechit_n;

  Float_t gp_pt[MAXPARTICLESTOSAVE];
  Float_t gp_eta[MAXPARTICLESTOSAVE];
  Float_t gp_phi[MAXPARTICLESTOSAVE];
  Int_t gp_idMC[MAXPARTICLESTOSAVE];
  Int_t gp_statusMC[MAXPARTICLESTOSAVE];
  Int_t gp_motherIdMC[MAXPARTICLESTOSAVE];

  Float_t gpho_pt[MAXPHOTONSTOSAVE];
  Float_t gpho_eta[MAXPHOTONSTOSAVE];
  Float_t gpho_phi[MAXPHOTONSTOSAVE];
  Int_t gpho_index[MAXPHOTONSTOSAVE];

  Float_t gele_pt[MAXPHOTONSTOSAVE];
  Float_t gele_eta[MAXPHOTONSTOSAVE];
  Float_t gele_phi[MAXPHOTONSTOSAVE];
  Float_t gele_fbrem80[MAXPHOTONSTOSAVE];
  Float_t gele_fbrem120[MAXPHOTONSTOSAVE];
  Int_t gele_index[MAXPHOTONSTOSAVE];

  Float_t rechit_e[MAXRECHITTOSAVE];
  Float_t rechit_ix[MAXRECHITTOSAVE];
  Float_t rechit_iy[MAXRECHITTOSAVE];


  Float_t pho_pt[MAXPHOTONSTOSAVE];
  Float_t pho_eta[MAXPHOTONSTOSAVE];
  Float_t pho_phi[MAXPHOTONSTOSAVE];

  Float_t   pho_e[MAXPHOTONSTOSAVE];
  Float_t   pho_esc[MAXPHOTONSTOSAVE];
  Float_t   pho_etasc[MAXPHOTONSTOSAVE];
  Float_t   pho_phisc[MAXPHOTONSTOSAVE];
  
  Float_t   pho_xsc[MAXPHOTONSTOSAVE];
  Float_t   pho_ysc[MAXPHOTONSTOSAVE];
  Float_t   pho_zsc[MAXPHOTONSTOSAVE];
  
  Float_t   pho_xcalo[MAXPHOTONSTOSAVE];
  Float_t   pho_ycalo[MAXPHOTONSTOSAVE];
  Float_t   pho_zcalo[MAXPHOTONSTOSAVE];
  
  Float_t   pho_hasPixelSeed[MAXPHOTONSTOSAVE];
  Float_t   pho_isEB[MAXPHOTONSTOSAVE];
  Float_t   pho_isEE[MAXPHOTONSTOSAVE];
  Float_t   pho_isEBEEGap[MAXPHOTONSTOSAVE];
  Float_t   pho_R9[MAXPHOTONSTOSAVE];
  Float_t   pho_E9[MAXPHOTONSTOSAVE];
  Float_t   pho_E25[MAXPHOTONSTOSAVE];

  Float_t  pho_jurECAL[MAXPHOTONSTOSAVE];
  Float_t  pho_twrHCAL[MAXPHOTONSTOSAVE];
  Float_t  pho_HoverE[MAXPHOTONSTOSAVE];
  Float_t  pho_hlwTrack[MAXPHOTONSTOSAVE];
  Float_t  pho_siEtaiEtaZS[MAXPHOTONSTOSAVE];
  Float_t  pho_siEtaiEtaNoZS[MAXPHOTONSTOSAVE];
  Float_t  pho_siEtaiEtaWrong[MAXPHOTONSTOSAVE];
  Float_t  pho_siPhiiPhiZS[MAXPHOTONSTOSAVE];
  Float_t  pho_siPhiiPhiNoZS[MAXPHOTONSTOSAVE];
  Float_t  pho_siPhiiPhiWrong[MAXPHOTONSTOSAVE];
  Float_t  pho_sMajZS[MAXPHOTONSTOSAVE];
  Float_t  pho_sMajNoZS[MAXPHOTONSTOSAVE];
  Float_t  pho_sMinZS[MAXPHOTONSTOSAVE];
  Float_t  pho_sMinNoZS[MAXPHOTONSTOSAVE];
  Float_t  pho_alphaZS[MAXPHOTONSTOSAVE];
  Float_t  pho_alphaNoZS[MAXPHOTONSTOSAVE];
  Float_t  pho_scetawid[MAXPHOTONSTOSAVE];
  Float_t  pho_scphiwid[MAXPHOTONSTOSAVE];
  Float_t  pho_jurECAL03[MAXPHOTONSTOSAVE];
  Float_t  pho_twrHCAL03[MAXPHOTONSTOSAVE];
  Float_t  pho_hlwTrack03[MAXPHOTONSTOSAVE];


  Float_t  pho_ptiso004[MAXPHOTONSTOSAVE];
  Int_t  pho_ntrkiso004[MAXPHOTONSTOSAVE];
  Float_t  pho_ptiso035[MAXPHOTONSTOSAVE];
  Int_t  pho_ntrkiso035[MAXPHOTONSTOSAVE];
  Float_t  pho_ptiso04[MAXPHOTONSTOSAVE];
  Int_t  pho_ntrkiso04[MAXPHOTONSTOSAVE];
      
  Float_t  pho_pfSumChargedHadronPt[MAXPHOTONSTOSAVE];
  Float_t  pho_pfsumNeutralHadronEt[MAXPHOTONSTOSAVE];
  Float_t  pho_pfsumPhotonEt[MAXPHOTONSTOSAVE];
  Float_t  pho_pfsumPUPt[MAXPHOTONSTOSAVE];

  Float_t pfSC_eta[MAXSCTOSAVE];
  Float_t pfSC_phi[MAXSCTOSAVE];
  Float_t   pfSC_e[MAXSCTOSAVE];
  Int_t pfSC_nBC[MAXSCTOSAVE];
  Int_t pfSC_nXtalsSeed[MAXSCTOSAVE];
  //bc info
  Float_t pfSC_bcEta[MAXSCTOSAVE][MAXBCTOSAVE];//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
  Float_t pfSC_bcPhi[MAXSCTOSAVE][MAXBCTOSAVE];
  Float_t   pfSC_bcE[MAXSCTOSAVE][MAXBCTOSAVE];
  Int_t pfSC_bcNXtals[MAXSCTOSAVE][MAXBCTOSAVE];

  //pfcand
  Float_t pho_pfCandPt  [MAXPHOTONSTOSAVE][MAXPFCANDTOSAVE]  ;//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
  Float_t pho_pfCandEta [MAXPHOTONSTOSAVE][MAXPFCANDTOSAVE] ;
  Float_t pho_pfCandPhi [MAXPHOTONSTOSAVE][MAXPFCANDTOSAVE] ;
  Float_t pho_pfCandVtx [MAXPHOTONSTOSAVE][MAXPFCANDTOSAVE] ;
  Float_t pho_pfCandType[MAXPHOTONSTOSAVE][MAXPFCANDTOSAVE];

  Float_t ele_pfCandPt  [MAXPHOTONSTOSAVE][MAXPFCANDTOSAVE]  ;//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
  Float_t ele_pfCandEta [MAXPHOTONSTOSAVE][MAXPFCANDTOSAVE] ;
  Float_t ele_pfCandPhi [MAXPHOTONSTOSAVE][MAXPFCANDTOSAVE] ;
  Float_t ele_pfCandVtx [MAXPHOTONSTOSAVE][MAXPFCANDTOSAVE] ;
  Float_t ele_pfCandType[MAXPHOTONSTOSAVE][MAXPFCANDTOSAVE];


  Float_t multi5x5SC_eta[MAXSCTOSAVE];
  Float_t multi5x5SC_phi[MAXSCTOSAVE];
  Float_t   multi5x5SC_e[MAXSCTOSAVE];
  Float_t multi5x5SC_nBC[MAXSCTOSAVE];
  Int_t multi5x5SC_nXtalsSeed[MAXSCTOSAVE];
  Float_t multi5x5SC_bcEta[MAXSCTOSAVE][MAXBCTOSAVE];//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
  Float_t multi5x5SC_bcPhi[MAXSCTOSAVE][MAXBCTOSAVE];
  Float_t   multi5x5SC_bcE[MAXSCTOSAVE][MAXBCTOSAVE];
  Int_t   multi5x5SC_bcNXtals[MAXSCTOSAVE][MAXBCTOSAVE];


  Float_t hybridSC_eta[MAXSCTOSAVE];
  Float_t hybridSC_phi[MAXSCTOSAVE];
  Float_t   hybridSC_e[MAXSCTOSAVE];
  Float_t hybridSC_nBC[MAXSCTOSAVE];
  Int_t hybridSC_nXtalsSeed[MAXSCTOSAVE];
  Float_t hybridSC_bcEta[MAXSCTOSAVE][MAXBCTOSAVE];//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
  Float_t hybridSC_bcPhi[MAXSCTOSAVE][MAXBCTOSAVE];
  Float_t   hybridSC_bcE[MAXSCTOSAVE][MAXBCTOSAVE];
  Int_t   hybridSC_bcNXtals[MAXSCTOSAVE][MAXBCTOSAVE];


  Float_t ele_pt[MAXPHOTONSTOSAVE];
  Float_t ele_eta[MAXPHOTONSTOSAVE];
  Float_t ele_phi[MAXPHOTONSTOSAVE];

  Float_t   ele_esc[MAXPHOTONSTOSAVE];
  Float_t   ele_etasc[MAXPHOTONSTOSAVE];
  Float_t   ele_phisc[MAXPHOTONSTOSAVE];
  
  Float_t   ele_xsc[MAXPHOTONSTOSAVE];
  Float_t   ele_ysc[MAXPHOTONSTOSAVE];
  Float_t   ele_zsc[MAXPHOTONSTOSAVE];
  
  Float_t   ele_xcalo[MAXPHOTONSTOSAVE];
  Float_t   ele_ycalo[MAXPHOTONSTOSAVE];
  Float_t   ele_zcalo[MAXPHOTONSTOSAVE];

  Float_t   ele_isEB[MAXPHOTONSTOSAVE];
  Float_t   ele_isEE[MAXPHOTONSTOSAVE];
  Float_t   ele_isEBEEGap[MAXPHOTONSTOSAVE];
  Float_t   ele_E9[MAXPHOTONSTOSAVE];
  Float_t   ele_E25[MAXPHOTONSTOSAVE];



  Float_t    ele_e[MAXPHOTONSTOSAVE];
  Float_t    ele_ecalEnergy[MAXPHOTONSTOSAVE];
  Float_t    ele_trackPatVtx[MAXPHOTONSTOSAVE];
  Float_t    ele_trackAtVtxPt[MAXPHOTONSTOSAVE];
  Float_t    ele_trackAtVtxEta[MAXPHOTONSTOSAVE];
  Float_t    ele_trackAtVtxPhi[MAXPHOTONSTOSAVE];
  Float_t    ele_trackAtCaloPt[MAXPHOTONSTOSAVE];
  Float_t    ele_trackAtCaloEta[MAXPHOTONSTOSAVE];
  Float_t    ele_trackAtCaloPhi[MAXPHOTONSTOSAVE];
  Float_t    ele_px[MAXPHOTONSTOSAVE];
  Float_t    ele_py[MAXPHOTONSTOSAVE];
  Float_t    ele_pz[MAXPHOTONSTOSAVE];
  Float_t    ele_vx[MAXPHOTONSTOSAVE];
  Float_t    ele_vy[MAXPHOTONSTOSAVE];
  Float_t    ele_vz[MAXPHOTONSTOSAVE];
  Float_t    ele_charge[MAXPHOTONSTOSAVE];
  Float_t    ele_fbrem[MAXPHOTONSTOSAVE];
  Float_t    elesc_fbrem[MAXPHOTONSTOSAVE];
  Float_t    elepfsc_fbrem[MAXPHOTONSTOSAVE];
  Float_t    ele_eSeedClusterOverPout[MAXPHOTONSTOSAVE];
  Float_t    ele_EoP[MAXPHOTONSTOSAVE];
  Float_t    ele_OneOverEMinusOneOverP[MAXPHOTONSTOSAVE];
  Float_t    ele_r9[MAXPHOTONSTOSAVE];
  Int_t      ele_misHits[MAXPHOTONSTOSAVE];
  Float_t  ele_dist[MAXPHOTONSTOSAVE];
  Float_t  ele_dcot[MAXPHOTONSTOSAVE];

  Int_t  ele_matchedConv[MAXPHOTONSTOSAVE];
  Int_t  ele_seedType[MAXPHOTONSTOSAVE];
  Int_t  ele_nSubClusters[MAXPHOTONSTOSAVE];
  Float_t  ele_HoE[MAXPHOTONSTOSAVE];
  Float_t  ele_pFlowMVA[MAXPHOTONSTOSAVE];
  Float_t  ele_siEtaiEtaZS[MAXPHOTONSTOSAVE];
  Float_t  ele_siEtaiEtaNoZS[MAXPHOTONSTOSAVE];
  Float_t  ele_siEtaiEtaWrong[MAXPHOTONSTOSAVE];
  Float_t  ele_siPhiiPhiZS[MAXPHOTONSTOSAVE];
  Float_t  ele_siPhiiPhiNoZS[MAXPHOTONSTOSAVE];
  Float_t  ele_siPhiiPhiWrong[MAXPHOTONSTOSAVE];
  Float_t  ele_sMajZS[MAXPHOTONSTOSAVE];
  Float_t  ele_sMajNoZS[MAXPHOTONSTOSAVE];
  Float_t  ele_sMinZS[MAXPHOTONSTOSAVE];
  Float_t  ele_sMinNoZS[MAXPHOTONSTOSAVE];
  Float_t  ele_alphaZS[MAXPHOTONSTOSAVE];
  Float_t  ele_alphaNoZS[MAXPHOTONSTOSAVE];
  Float_t  ele_scetawid[MAXPHOTONSTOSAVE];
  Float_t  ele_scphiwid[MAXPHOTONSTOSAVE];
  Float_t  ele_jurECAL03[MAXPHOTONSTOSAVE];
  Float_t  ele_twrHCAL03[MAXPHOTONSTOSAVE];
  Float_t  ele_hlwTrack03[MAXPHOTONSTOSAVE];
  Float_t  ele_trkIso[MAXPHOTONSTOSAVE];
  Float_t  ele_ecalIso[MAXPHOTONSTOSAVE];
  Float_t  ele_hcalIso[MAXPHOTONSTOSAVE];
  Float_t  ele_trkIso03[MAXPHOTONSTOSAVE];
  Float_t  ele_ecalIso03[MAXPHOTONSTOSAVE];
  Float_t  ele_hcalIso03[MAXPHOTONSTOSAVE];
  Float_t  ele_dEtaIn[MAXPHOTONSTOSAVE];
  Float_t  ele_dPhiIn[MAXPHOTONSTOSAVE];
  Float_t  ele_sc_energy[MAXPHOTONSTOSAVE];
  Float_t  ele_sc_eta[MAXPHOTONSTOSAVE];
  Float_t  ele_sc_phi[MAXPHOTONSTOSAVE];
  Float_t  ele_pfSumChargedHadronPt[MAXPHOTONSTOSAVE];
  Float_t  ele_pfsumNeutralHadronEt[MAXPHOTONSTOSAVE];
  Float_t  ele_pfsumPhotonEt[MAXPHOTONSTOSAVE];
  Float_t  ele_pfsumPUPt[MAXPHOTONSTOSAVE];


  float truePU;
  int bxPU[16];
  Int_t nvtx;
  Float_t vertex_x[200];
  Float_t vertex_y[200];
  Float_t vertex_z[200];

  bool isData;
  bool saveReco;

  int entryNumber;

};

void dumper::beginRun(const edm::Run& run,const edm::EventSetup& setup) {
  std::cout <<"begining run "<<std::endl;
  entryNumber=0;
}


dumper::dumper(const edm::ParameterSet& iConfig) {
  outputFileName  = iConfig.getParameter<std::string>("OutputFileName");
  isData          = iConfig.getParameter<bool>("isData");
  saveReco        = iConfig.getParameter<bool>("saveReco");
}

dumper::~dumper() 
{}


void dumper::clearVector(){

  for(int i=0;i<MAXSCTOSAVE;++i){
    for(int j=0;j<MAXBCTOSAVE;++j){
      pfSC_bcEta[i][j]=0.;
      pfSC_bcPhi[i][j]=0.;
      pfSC_bcE[i][j]=0.;
      pfSC_bcNXtals[i][j]=0;

      multi5x5SC_bcEta[i][j]=0.;
      multi5x5SC_bcPhi[i][j]=0.;
      multi5x5SC_bcE[i][j]=0.;
      multi5x5SC_bcNXtals[i][j]=0;

      hybridSC_bcEta[i][j]=0.;
      hybridSC_bcPhi[i][j]=0.;
      hybridSC_bcE[i][j]=0.;
      hybridSC_bcNXtals[i][j]=0;

    }
  }


  for(int i=0;i<MAXPHOTONSTOSAVE;++i){
    for(int j=0;j<MAXPFCANDTOSAVE;++j){ 
      pho_pfCandPt[i][j]=-1.;//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
      pho_pfCandEta[i][j]=0.;
      pho_pfCandPhi[i][j]=0.;
      pho_pfCandVtx[i][j]=-1.;
      pho_pfCandType[i][j]=-1.;
      
      ele_pfCandPt[i][j]=-1.;//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
      ele_pfCandEta[i][j]=0.;
      ele_pfCandPhi[i][j]=0.;
      ele_pfCandVtx[i][j]=-1.;
      ele_pfCandType[i][j]=-1.;
    }
  }

}


void dumper::recHitReco(const EERecHitCollection* rhitsee){
  
  //save rechits only for first 100 evts
  rechit_n=0;
  
  for (EcalRecHitCollection::const_iterator itRecHit = rhitsee->begin();
       itRecHit != rhitsee->end(); ++itRecHit) {
    if(itRecHit->energy()>0){
      rechit_e[rechit_n]=itRecHit->energy();
      EKDetId ekId(itRecHit->id());
      rechit_ix[rechit_n]=ekId.ix();
      rechit_iy[rechit_n]=ekId.iy();
      rechit_n++;
    }
  }
  
  
}

void dumper::phoReco(edm::Handle<reco::PhotonCollection> phoH, edm::Handle<reco::TrackCollection> tracksH, const EBRecHitCollection* rhitseb,const EERecHitCollection* rhitsee,edm::Handle<reco::PFCandidateCollection>  PFCandidates){

  pho_n=0;

  for(size_t i = 0; i < phoH->size(); ++i) {


    const reco::PhotonRef pho(phoH,i);

    if (pho->pt() > 0. && pho_n<MAXPHOTONSTOSAVE) {
      
      pho_pt[pho_n]  = pho->pt();
      pho_eta[pho_n] = pho->eta();
      pho_phi[pho_n] = pho->phi();

      pho_e[pho_n] = pho->energy();
      pho_esc[pho_n] = pho->superCluster()->energy();
      pho_etasc[pho_n] = pho->superCluster()->position().eta();
      pho_phisc[pho_n] = pho->superCluster()->position().phi();
      
      pho_xsc[pho_n] = pho->superCluster()->position().x();
      pho_ysc[pho_n] = pho->superCluster()->position().y();
      pho_zsc[pho_n] = pho->superCluster()->position().z();
      
      pho_xcalo[pho_n] = pho->caloPosition().x();
      pho_ycalo[pho_n] = pho->caloPosition().y();
      pho_zcalo[pho_n] = pho->caloPosition().z();
      
      pho_hasPixelSeed[pho_n] = pho->hasPixelSeed();
      pho_isEB[pho_n] = pho->isEB();
      pho_isEE[pho_n] = pho->isEE();
      pho_isEBEEGap[pho_n] = pho->isEBEEGap();
      pho_R9[pho_n] = pho->r9();
      pho_E9[pho_n] = pho->e3x3();
      pho_E25[pho_n] = pho->e5x5();

      //default photon Variables
      pho_jurECAL[pho_n] = pho->ecalRecHitSumEtConeDR04();//isolationEcalRecHit
      pho_twrHCAL[pho_n] = pho->hcalTowerSumEtConeDR04();//isolationHcalRecHit
      pho_HoverE[pho_n] = pho->hadronicOverEm();
      pho_hlwTrack[pho_n] = pho->trkSumPtHollowConeDR04();//isolationHollowTrkCone

      const edm::Ptr<reco::CaloCluster> SCseed = pho->superCluster()->seed(); 
      const EBRecHitCollection* rechits = (SCseed->seed().subdetId()==EcalBarrel) ? rhitseb : rhitsee;

      myClusterTools::Cluster2ndMoments momentsNoZS = myClusterTools::noZSEcalClusterTools::cluster2ndMoments(*SCseed, *rechits);
      myClusterTools::Cluster2ndMoments momentsZS = myClusterTools::EcalClusterTools::cluster2ndMoments(*SCseed, *rechits);//this namespace is because in cmssw6 ecalclustertoos was not with ZS, i had to import it from cmmssw7 by ahand and redefine namespaces in Dumper/dumper/interface/EcalClusterTools.h


      pho_sMajNoZS[pho_n]=momentsNoZS.sMaj;
      pho_sMajZS[pho_n]=momentsZS.sMaj;



      pho_sMinNoZS[pho_n]=momentsNoZS.sMin;
      pho_sMinZS[pho_n]=momentsZS.sMin;

      pho_alphaNoZS[pho_n]=momentsNoZS.alpha;
      pho_alphaZS[pho_n]=momentsZS.alpha;


      std::vector<float> etaphimomentsnoZS = myClusterTools::noZSEcalClusterTools::localCovariances(*SCseed, &(*rechits), &(*topology));
      std::vector<float> etaphimomentsZS = myClusterTools::EcalClusterTools::localCovariances(*SCseed, &(*rechits), &(*topology));
      std::vector<float> etaphimomentsWrong = EcalClusterTools::localCovariances(*SCseed, &(*rechits), &(*topology));

      pho_siEtaiEtaNoZS[pho_n] = sqrt(etaphimomentsnoZS[0]);
      pho_siEtaiEtaZS[pho_n] = sqrt(etaphimomentsZS[0]);
      pho_siEtaiEtaWrong[pho_n] = sqrt(etaphimomentsWrong[0]);

      pho_siPhiiPhiNoZS[pho_n] = sqrt(etaphimomentsnoZS[2]);
      pho_siPhiiPhiZS[pho_n] = sqrt(etaphimomentsZS[2]);
      pho_siPhiiPhiWrong[pho_n] = sqrt(etaphimomentsWrong[2]);


      pho_scetawid[pho_n]       = pho->superCluster()->etaWidth();
      pho_scphiwid[pho_n]        = pho->superCluster()->phiWidth();
      pho_jurECAL03[pho_n] = pho->ecalRecHitSumEtConeDR03();//isolationEcalRecHit 
      pho_twrHCAL03[pho_n] = pho->hcalTowerSumEtConeDR03();//isolationHcalRecHit
      pho_hlwTrack03[pho_n] = pho->trkSumPtHollowConeDR03();//isolationHollowTrkCone



      //isolation with different cones
      double ptiso004(0.), ptiso035(0.), ptiso04(0.);
      int ntrkiso004(0), ntrkiso035(0), ntrkiso04(0);

      for (reco::TrackCollection::const_iterator itTrack = tracksH->begin();
	   itTrack != tracksH->end(); ++itTrack) {
	
	double etaTrack = itTrack->eta();
	double phiTrack = itTrack->phi();
	
	double deltaPhi = phiTrack-pho->phi();
	double deltaEta = etaTrack-pho->eta();
	if (deltaPhi > Geom::pi()) deltaPhi -= 2.*Geom::pi();
	if (deltaPhi < -Geom::pi()) deltaPhi += 2.*Geom::pi();
	double deltaR = std::sqrt(deltaEta*deltaEta+deltaPhi*deltaPhi);

	if (deltaR < .04)  {ptiso004  += sqrt(itTrack->pt()); ntrkiso004++; }
	if (deltaR < .35)   {ptiso035 += sqrt(itTrack->pt()); ntrkiso035++;}
	if (deltaR < .4)   {ptiso04 += sqrt(itTrack->pt()); ntrkiso04++;}

      }

      pho_ptiso004[pho_n] = ptiso004;
      pho_ntrkiso004[pho_n] = ntrkiso004;
      pho_ptiso035[pho_n] = ptiso035;
      pho_ntrkiso035[pho_n] = ntrkiso035;
      pho_ptiso04[pho_n] = ptiso04;
      pho_ntrkiso04[pho_n] = ntrkiso04;

      math::XYZTLorentzVectorD const& p4phot = pho->p4();
      TLorentzVector p4pho(p4phot.px(), p4phot.py(), p4phot.pz(), p4phot.energy());

      pfcandPho_n=0;
      for (reco::PFCandidateCollection::const_iterator jt = PFCandidates->begin();
	   jt != PFCandidates->end(); ++jt) {

	reco::PFCandidate::ParticleType id = jt->particleId();
	// Convert particle momentum to normal TLorentzVector, wrong type :(
	math::XYZTLorentzVectorD const& p4t = jt->p4();
	TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());


	if(p4pho.DeltaR(p4)<0.4 && pfcandPho_n<MAXPFCANDTOSAVE){
	  pho_pfCandPt  [pho_n][pfcandPho_n]  =p4t.pt();
	  pho_pfCandEta [pho_n][pfcandPho_n] =p4t.eta();
	  pho_pfCandPhi [pho_n][pfcandPho_n] =p4t.phi();
	  pho_pfCandType[pho_n][pfcandPho_n]=id;

	  //	  std::cout<<"index"<<pfcandPho_n<<" "<<pho_n<< " pt:"<<pho_pfCandPt[pfcandPho_n][pho_n]<<" "<<pho_pfCandEta[pfcandPho_n][pho_n]<<" "<<pho_pfCandPhi[pfcandPho_n][pho_n]<<" "<<pho_pfCandType[pfcandPho_n][pho_n]<<std::endl;
	  // Get the primary vertex coordinates
	  int closestVertex=0;
	  float dist_vtx=999;
	  if(nvtx>0){
	  dist_vtx=fabs(vertex_z[0]-jt->vz());
	  for(int i=0;i<nvtx;++i){
	    float dist=fabs(vertex_z[i]-jt->vz());
	    if(dist<dist_vtx){
	      closestVertex=i;
	      dist_vtx=dist;
	    }
	  }
	  }
	  if(dist_vtx<0.5){
	    pho_pfCandVtx[pho_n][pfcandPho_n]=closestVertex;
	  }
	  pfcandPho_n++;
	}



      }
      

      pho_pfSumChargedHadronPt[pho_n] = pho->chargedHadronIso();
      pho_pfsumNeutralHadronEt[pho_n] = pho->neutralHadronIso();
      pho_pfsumPhotonEt[pho_n] = pho->photonIso();



      pho_n++;
    }

  }

}

void dumper::eleReco(edm::Handle<reco::GsfElectronCollection> ElectronHandle, edm::Handle<reco::ConversionCollection> hConversions,edm::Handle<reco::BeamSpot> recoBeamSpotHandle, const EBRecHitCollection* rhitseb,const EERecHitCollection* rhitsee,edm::Handle<reco::PFCandidateCollection>  PFCandidates){

  ele_n=0;  

  unsigned index_gsf = 0;
  for (reco::GsfElectronCollection::const_iterator itElectron = ElectronHandle->begin();
       itElectron != ElectronHandle->end(); ++itElectron, ++index_gsf) {
     
    reco::GsfElectronRef eleRef = reco::GsfElectronRef(ElectronHandle, index_gsf);

    if (itElectron->pt()<2.5||ele_n>=MAXPHOTONSTOSAVE)
      continue;


    ele_pt[ele_n]          = itElectron->pt();
    ele_e[ele_n]      = itElectron->energy();
    ele_ecalEnergy[ele_n]  = itElectron->ecalEnergy();              
    ele_trackPatVtx[ele_n] = itElectron->trackMomentumAtVtx().R();

    ele_esc[ele_n] = itElectron->superCluster()->energy();
    ele_etasc[ele_n] = itElectron->superCluster()->position().eta();
    ele_phisc[ele_n] = itElectron->superCluster()->position().phi();
    
    ele_xsc[ele_n] = itElectron->superCluster()->position().x();
    ele_ysc[ele_n] = itElectron->superCluster()->position().y();
    ele_zsc[ele_n] = itElectron->superCluster()->position().z();
    
    ele_xcalo[ele_n] = itElectron->caloPosition().x();
    ele_ycalo[ele_n] = itElectron->caloPosition().y();
    ele_zcalo[ele_n] = itElectron->caloPosition().z();
    
    ele_isEB[ele_n] = itElectron->isEB();
    ele_isEE[ele_n] = itElectron->isEE();
    ele_isEBEEGap[ele_n] = itElectron->isEBEEGap();
    ele_r9[ele_n] = itElectron->r9();
    ele_E9[ele_n] = ele_r9[ele_n]*ele_esc[ele_n];
    ele_E25[ele_n] = itElectron->e5x5();



    math::XYZVectorF p3Vtx = itElectron ->trackMomentumAtVtx();
    math::XYZVectorF p3Calo = itElectron ->trackMomentumAtCalo();
    TLorentzVector p4Vtx;  
    TLorentzVector p4Calo;  

    p4Vtx.SetXYZM(p3Vtx.x(),p3Vtx.y(),p3Vtx.z(),0.);
    p4Calo.SetXYZM(p3Calo.x(),p3Calo.y(),p3Calo.z(),0.);

    ele_trackAtVtxPt[ele_n] = p4Vtx.Pt();
    ele_trackAtVtxEta[ele_n] = p4Vtx.Eta();
    ele_trackAtVtxPhi[ele_n] = p4Vtx.Phi();
    ele_trackAtCaloPt[ele_n] = p4Calo.Pt();
    ele_trackAtCaloEta[ele_n] = p4Calo.Eta();
    ele_trackAtCaloPhi[ele_n] = p4Calo.Phi();



    ele_px[ele_n]       = itElectron->px();
    ele_py[ele_n]       = itElectron->py();
    ele_pz[ele_n]       = itElectron->pz();
    ele_vx[ele_n]       = itElectron->vx();
    ele_vy[ele_n]       = itElectron->vy();
    ele_vz[ele_n]       = itElectron->vz();
    ele_phi[ele_n]      = itElectron->phi();
    ele_eta[ele_n]      = itElectron->eta();
    ele_charge[ele_n]       = itElectron->charge();
    ele_fbrem[ele_n]       = itElectron->fbrem();
    elesc_fbrem[ele_n]       = itElectron->superClusterFbrem();
    elepfsc_fbrem[ele_n]       = itElectron->pfSuperClusterFbrem();
    ele_eSeedClusterOverPout[ele_n]       = itElectron->eSeedClusterOverPout();
    ele_EoP[ele_n]       = itElectron->eSuperClusterOverP();
    ele_OneOverEMinusOneOverP[ele_n]       = (1/itElectron->caloEnergy())-(1/itElectron ->trackMomentumAtVtx().R());

    ele_misHits[ele_n]       = itElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits();

    ele_dist[ele_n] = itElectron->convDist();
    ele_dcot[ele_n] = itElectron->convDcot();


    bool matchesConv = ConversionTools::hasMatchedConversion(*eleRef,hConversions,recoBeamSpotHandle->position());
     ele_matchedConv[ele_n] = matchesConv;
     ele_seedType[ele_n]       = itElectron->ecalDrivenSeed()+2*itElectron->trackerDrivenSeed();
     ele_nSubClusters[ele_n]       = itElectron->numberOfBrems();
     ele_HoE[ele_n]          = itElectron->hadronicOverEm();
     ele_pFlowMVA[ele_n]          = itElectron->mvaOutput().mva;



      const edm::Ptr<reco::CaloCluster> SCseed = itElectron->superCluster()->seed(); 
      const EBRecHitCollection* rechits = (SCseed->seed().subdetId()==EcalBarrel) ? rhitseb : rhitsee;

      myClusterTools::Cluster2ndMoments momentsNoZS = myClusterTools::noZSEcalClusterTools::cluster2ndMoments(*SCseed, *rechits);
      myClusterTools::Cluster2ndMoments momentsZS = myClusterTools::EcalClusterTools::cluster2ndMoments(*SCseed, *rechits);//this namespace is because in cmssw6 ecalclustertoos was not with ZS, i had to import it from cmmssw7 by ahand and redefine namespaces in Dumper/dumper/interface/EcalClusterTools.h


      ele_sMajNoZS[ele_n]=momentsNoZS.sMaj;
      ele_sMajZS[ele_n]=momentsZS.sMaj;



      ele_sMinNoZS[ele_n]=momentsNoZS.sMin;
      ele_sMinZS[ele_n]=momentsZS.sMin;

      ele_alphaNoZS[ele_n]=momentsNoZS.alpha;
      ele_alphaZS[ele_n]=momentsZS.alpha;



      std::vector<float> etaphimomentsnoZS = myClusterTools::noZSEcalClusterTools::localCovariances(*SCseed, &(*rechits), &(*topology));
      std::vector<float> etaphimomentsZS = myClusterTools::EcalClusterTools::localCovariances(*SCseed, &(*rechits), &(*topology));
      std::vector<float> etaphimomentsWrong = EcalClusterTools::localCovariances(*SCseed, &(*rechits), &(*topology));

      ele_siEtaiEtaNoZS[ele_n] = sqrt(etaphimomentsnoZS[0]);
      ele_siEtaiEtaZS[ele_n] = sqrt(etaphimomentsZS[0]);
      ele_siEtaiEtaWrong[ele_n] = sqrt(etaphimomentsWrong[0]);

      ele_siPhiiPhiNoZS[ele_n] = sqrt(etaphimomentsnoZS[2]);
      ele_siPhiiPhiZS[ele_n] = sqrt(etaphimomentsZS[2]);
      ele_siPhiiPhiWrong[ele_n] = sqrt(etaphimomentsWrong[2]);


      math::XYZTLorentzVectorD const& p4elet = itElectron->p4();
      TLorentzVector p4ele(p4elet.px(), p4elet.py(), p4elet.pz(), p4elet.energy());

      pfcandEle_n=0;
      for (reco::PFCandidateCollection::const_iterator jt = PFCandidates->begin();
	   jt != PFCandidates->end(); ++jt) {
	
	reco::PFCandidate::ParticleType id = jt->particleId();
	// Convert particle momentum to normal TLorentzVector, wrong type :(
	math::XYZTLorentzVectorD const& p4t = jt->p4();
	TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());


	if(p4ele.DeltaR(p4)<0.4 && pfcandEle_n<MAXPFCANDTOSAVE){
	  ele_pfCandPt  [ele_n][pfcandEle_n]  =p4t.pt();
	  ele_pfCandEta [ele_n][pfcandEle_n] =p4t.eta();
	  ele_pfCandPhi [ele_n][pfcandEle_n] =p4t.phi();
	  ele_pfCandType[ele_n][pfcandEle_n]=id;

	  //	  std::cout<<"index"<<pfcandEle_n<<" "<<ele_n<< " pt:"<<ele_pfCandPt[pfcandEle_n][ele_n]<<" "<<ele_pfCandEta[pfcandEle_n][ele_n]<<" "<<ele_pfCandPhi[pfcandEle_n][ele_n]<<" "<<ele_pfCandType[pfcandEle_n][ele_n]<<"elen"<<ele_n<<std::endl;

	  // Get the primary vertex coordinates
	  int closestVertex=0;
	  float dist_vtx=999;
	  if (nvtx>0){
	    dist_vtx=fabs(vertex_z[0]-jt->vz());
	    for(int i=0;i<nvtx;++i){
	      float dist=fabs(vertex_z[i]-jt->vz());
	      if(dist<dist_vtx){
		closestVertex=i;
		dist_vtx=dist;
	      }
	    }
	  }
	  //	  if(pfcandEle_n==0 && ele_n ==200)	    std::cout<<"closestVtx"<<closestVertex<<" "<<ele_pfCandPt[ele_n][pfcandEle_n]<<"'"<<p4t.pt()<<std::endl;
	  if(dist_vtx<0.5){

	    ele_pfCandVtx[ele_n][pfcandEle_n]=closestVertex;
	  }
	  pfcandEle_n++;
	}



      }


     

     ele_trkIso[ele_n]       = itElectron->dr04TkSumPt() ;
     ele_ecalIso[ele_n]      = itElectron->dr04EcalRecHitSumEt();
     ele_hcalIso[ele_n]      = itElectron->dr04HcalTowerSumEt();
     ele_trkIso03[ele_n]       = itElectron->dr03TkSumPt() ;
     ele_ecalIso03[ele_n]      = itElectron->dr03EcalRecHitSumEt();
     ele_hcalIso03[ele_n]      = itElectron->dr03HcalTowerSumEt();
     ele_dEtaIn[ele_n]       = itElectron->deltaEtaSuperClusterTrackAtVtx();
     ele_dPhiIn[ele_n]       = itElectron->deltaPhiSuperClusterTrackAtVtx();
     ele_sc_energy[ele_n]    = itElectron->superCluster()->energy();
     ele_sc_eta[ele_n]       = itElectron->superCluster()->eta();
     ele_sc_phi[ele_n]       = itElectron->superCluster()->phi();
     reco::GsfElectron::PflowIsolationVariables pfIso = itElectron->pfIsolationVariables();
     ele_pfSumChargedHadronPt[ele_n] = pfIso.chargedHadronIso;
     ele_pfsumNeutralHadronEt[ele_n] = pfIso.neutralHadronIso;
     ele_pfsumPhotonEt[ele_n] = pfIso.photonIso;
     //     ele_pfsumPUPt[ele_n] = pfIso.sumPUPt; //not  in cmssw6


     ele_n++;
  }



}


void dumper::multi5x5scReco(edm::Handle<reco::SuperClusterCollection> multi5x5Handle){
  multi5x5SC_n=0;
  for (reco::SuperClusterCollection::const_iterator itSC = multi5x5Handle->begin();
       itSC != multi5x5Handle->end(); ++itSC) {

    if (itSC->energy() > 0. && multi5x5SC_n<MAXSCTOSAVE) {


      multi5x5SC_eta[multi5x5SC_n] = itSC->eta();
      multi5x5SC_phi[multi5x5SC_n] = itSC->phi();

      multi5x5SC_e[multi5x5SC_n] = itSC->energy();
      multi5x5SC_nBC[multi5x5SC_n] = itSC->clustersSize();
      multi5x5SC_nXtalsSeed[multi5x5SC_n] = itSC->seed()->size();

      //basicClusters
      int nBC=0;
      for (reco::CaloCluster_iterator bclus = (itSC->clustersBegin()); bclus != (itSC->clustersEnd()); ++bclus) {
	if((*bclus)->energy() > 0 && multi5x5SC_nBC[multi5x5SC_n]<MAXBCTOSAVE){
	  multi5x5SC_bcPhi[multi5x5SC_n][nBC]=(*bclus)->phi();
	  multi5x5SC_bcEta[multi5x5SC_n][nBC]=(*bclus)->eta(); 
	  multi5x5SC_bcE[multi5x5SC_n][nBC]=(*bclus)->energy(); 
	  multi5x5SC_bcNXtals[multi5x5SC_n][nBC]=(*bclus)->size(); 
	  //	  std::cout<<"multi 5x5 "<<multi5x5SC_n<<" : "<<nBC<<" : "<< multi5x5SC_nBC[multi5x5SC_n]<< " en: "<<multi5x5SC_bcE[multi5x5SC_n][nBC]<<" eta:"<<multi5x5SC_bcEta[multi5x5SC_n][nBC]<<std::endl;
	  nBC++;
	}

      }


      multi5x5SC_n++;

    }
  

  }

}


void dumper::hybridscReco(edm::Handle<reco::SuperClusterCollection> hybridHandle){
  hybridSC_n=0;
  for (reco::SuperClusterCollection::const_iterator itSC = hybridHandle->begin();
       itSC != hybridHandle->end(); ++itSC) {

    if (itSC->energy() > 0. && hybridSC_n<MAXSCTOSAVE) {


      hybridSC_eta[hybridSC_n] = itSC->eta();
      hybridSC_phi[hybridSC_n] = itSC->phi();

      hybridSC_e[hybridSC_n] = itSC->energy();
      hybridSC_nBC[hybridSC_n] = itSC->clustersSize();
      hybridSC_nXtalsSeed[hybridSC_n] = itSC->seed()->size();

      //basicClusters
      int nBC=0;
      for (reco::CaloCluster_iterator bclus = (itSC->clustersBegin()); bclus != (itSC->clustersEnd()); ++bclus) {
	if((*bclus)->energy() > 0 && hybridSC_nBC[hybridSC_n]<MAXBCTOSAVE){
	  hybridSC_bcPhi[hybridSC_n][nBC]=(*bclus)->phi();
	  hybridSC_bcEta[hybridSC_n][nBC]=(*bclus)->eta(); 
	  hybridSC_bcE[hybridSC_n][nBC]=(*bclus)->energy(); 
	  hybridSC_bcNXtals[hybridSC_n][nBC]=(*bclus)->size(); 
	  //	  std::cout<<hybridSC_n<<" : "<<nBC<<" : "<< hybridSC_nBC[hybridSC_n]<< " en: "<<hybridSC_bcE[hybridSC_n][nBC]<<std::endl;
	  nBC++;
	}

      }



      hybridSC_n++;

    }


  }

}


void dumper::scReco(edm::Handle<reco::SuperClusterCollection> superClustersEBHandle, edm::Handle<reco::SuperClusterCollection> superClustersEEHandle){

  pfSC_n=0;

  for (reco::SuperClusterCollection::const_iterator itSC = superClustersEBHandle->begin();
       itSC != superClustersEBHandle->end(); ++itSC) {

    if (itSC->energy() > 0. && pfSC_n<MAXSCTOSAVE) {


      pfSC_eta[pfSC_n] = itSC->eta();
      pfSC_phi[pfSC_n] = itSC->phi();

      pfSC_e[pfSC_n] = itSC->energy();
      pfSC_nBC[pfSC_n] = itSC->clustersSize();
      pfSC_nXtalsSeed[pfSC_n] = itSC->seed()->size();
      //basicClusters
      int nBC=0;
      for (reco::CaloCluster_iterator bclus = (itSC->clustersBegin()); bclus != (itSC->clustersEnd()); ++bclus) {
	if((*bclus)->energy() > 0 && pfSC_nBC[pfSC_n]<MAXBCTOSAVE){
	  pfSC_bcPhi[pfSC_n][nBC]=(*bclus)->phi();
	  pfSC_bcEta[pfSC_n][nBC]=(*bclus)->eta(); 
	  pfSC_bcE[pfSC_n][nBC]=(*bclus)->energy(); 
	  pfSC_bcNXtals[pfSC_n][nBC]=(*bclus)->size(); 
	  //	  std::cout<<"pfSC "<<pfSC_n<<" : "<<nBC<<" : "<< pfSC_nBC[pfSC_n]<< " en: "<<pfSC_bcE[pfSC_n][nBC]<<std::endl;
	  nBC++;
	}

      }

      pfSC_n++;
    }


  }


  for (reco::SuperClusterCollection::const_iterator itSC = superClustersEEHandle->begin();
       itSC != superClustersEEHandle->end(); ++itSC) {
    if (itSC->energy() > 0. && pfSC_n<MAXSCTOSAVE) {
      

      pfSC_eta[pfSC_n] = itSC->eta();
      pfSC_phi[pfSC_n] = itSC->phi();

      pfSC_e[pfSC_n] = itSC->energy();
      pfSC_nBC[pfSC_n] = itSC->clustersSize();
      pfSC_nXtalsSeed[pfSC_n] = itSC->seed()->size();


      //      std::cout<<" seed E:"<<itSC->seed()->energy()<<" seed Eta:"<<itSC->seed()->eta()<<std::endl;
      //basicClusters
      int nBC=0;
      for (reco::CaloCluster_iterator bclus = (itSC->clustersBegin()); bclus != (itSC->clustersEnd()); ++bclus) {
	if((*bclus)->energy() > 0 && pfSC_nBC[pfSC_n]<MAXBCTOSAVE){
	  pfSC_bcPhi[pfSC_n][nBC]=(*bclus)->phi();
	  pfSC_bcEta[pfSC_n][nBC]=(*bclus)->eta(); 
	  pfSC_bcE[pfSC_n][nBC]=(*bclus)->energy(); 
	  pfSC_bcNXtals[pfSC_n][nBC]=(*bclus)->size(); 
	  //	  std::cout<<"pfsc"<<pfSC_n<<" : "<<nBC<<" : "<< pfSC_nBC[pfSC_n]<< " en: "<<pfSC_bcE[pfSC_n][nBC]<<std::endl;
	  nBC++;
	}

      }

      pfSC_n++;
    }


  }


}




void dumper::mcTruth(edm::Handle<reco::GenParticleCollection> gpH,std::vector<ElectronMCTruth> MCElectrons) {
  
  gp_n = 0;
  gpho_n = 0;
  gele_n = 0;

  for(size_t i = 0; i < gpH->size(); ++i) {

    const reco::GenParticleRef gp(gpH, i);
    //    if(gp->status()==3)    std::cout<<gp->status()<<" "<<gp->pdgId()<<std::endl;    
    //    if ((gp->status() >=20 and gp->status() <=29) and abs(gp->pdgId()) == 11) {
    //    if ((gp->status() ==3) ){
    if (gp->pt() > 3.) {//loose pt cut just to save space
      if(gp_n<MAXPARTICLESTOSAVE){
	gp_pt[gp_n]  = gp->pt();
	gp_eta[gp_n] = gp->eta();
	gp_phi[gp_n] = gp->phi();
	gp_idMC [gp_n] = gp->pdgId();
	gp_statusMC [gp_n] = gp->status();
	gp_motherIdMC [gp_n] = gp->mother()->pdgId();
	gp_n++;

	if((abs(gp->pdgId()) == 22)){
	  if(gpho_n<MAXPHOTONSTOSAVE){
	    gpho_pt[gpho_n]  = gp->pt();
	    gpho_eta[gpho_n] = gp->eta();
	    gpho_phi[gpho_n] = gp->phi();
	    gpho_index [gpho_n] = gp_n-1;
	    gpho_n++;
	  }
	}
	
	if(gele_n<MAXPHOTONSTOSAVE){
	  if((abs(gp->pdgId()) == 11)){
	    gele_pt[gele_n]  = gp->pt();
	    gele_eta[gele_n] = gp->eta();
	    gele_phi[gele_n] = gp->phi();
	    gele_index [gele_n] = gp_n-1;
	    for ( std::vector<ElectronMCTruth>::const_iterator iEl=MCElectrons.begin(); iEl !=MCElectrons.end(); ++iEl ){ 

	      //matching
	      float eleEta=(*iEl).fourMomentum().eta();
	      float elePhi=(*iEl).fourMomentum().phi();

	      if(sqrt((gele_eta[gele_n]-eleEta)*(gele_eta[gele_n]-eleEta)+(gele_phi[gele_n]-elePhi)*(gele_phi[gele_n]-elePhi))>0.5) continue;

	      float totBrem80=0 ;
	      float totBrem120=0 ;
	      unsigned int iBrem ;
	      for ( iBrem=0; iBrem < (*iEl).bremVertices().size(); ++iBrem ) {

		float rBrem= (*iEl).bremVertices()[iBrem].perp();
		if(rBrem < 80){
		  totBrem80 +=  (*iEl).bremMomentum()[iBrem].e();   
		} 
		if ( rBrem < 120 ) {
		  totBrem120 +=  (*iEl).bremMomentum()[iBrem].e();   
		}
	      }
	      gele_fbrem80[gele_n]=totBrem80/((*iEl).fourMomentum().e());
	      gele_fbrem120[gele_n]=totBrem120/((*iEl).fourMomentum().e());
	    }
	    gele_n++;
	  }
	}
	
      }else{
	continue;
      }

    }
  }

}

void dumper::analyze(const edm::Event& event, const edm::EventSetup& iSetup) {

  entryNumber++;

  ///////////////////////Get PU informations
  // rho from fast jet
  edm::Handle<double> rhoH;
  //iEvent.getByLabel(edm::InputTag("kt6PFJets","rho","Iso"),rhoH); 
  if( event.getByLabel(edm::InputTag("kt6PFJets","rho"),rhoH) )
     rho = *rhoH;
  else 
    rho = 0;
  
  
  //PU info    
  edm::Handle<std::vector<PileupSummaryInfo>> puH;
  event.getByLabel("addPileupInfo", puH);
  truePU = (*puH)[0].getTrueNumInteractions();
  for (unsigned int j=0; j<puH->size(); j++)
    bxPU[j] = (*puH)[j].getPU_NumInteractions();
  
  edm::Handle<reco::VertexCollection> VertexHandle;
  event.getByLabel("offlinePrimaryVerticesWithBS", VertexHandle);
  // Get the primary vertex coordinates
  nvtx=0;
  for (reco::VertexCollection::const_iterator it = VertexHandle->begin(); 
       it != VertexHandle->end(); ++it) {
    vertex_x[nvtx] = (it->isValid()) ? it->x() : 999.;
    vertex_y[nvtx] = (it->isValid()) ? it->y() : 999.;
    vertex_z[nvtx] = (it->isValid()) ? it->z() : 999.;

    nvtx++;
  }
  

    
  
  if (!isData) {

    //gen particles
    edm::Handle<reco::GenParticleCollection> gpH;
    event.getByLabel("genParticles", gpH);
    
    //get simtrack info
    std::vector<SimTrack> theSimTracks;
    std::vector<SimVertex> theSimVertices;
    
    edm::Handle<edm::SimTrackContainer> SimTk;
    edm::Handle<edm::SimVertexContainer> SimVtx;
    
    event.getByLabel("g4SimHits",SimTk);
    event.getByLabel("g4SimHits",SimVtx); 
    theSimTracks.insert(theSimTracks.end(),SimTk->begin(),SimTk->end());
    theSimVertices.insert(theSimVertices.end(),SimVtx->begin(),SimVtx->end());
    
    ElectronMCTruthFinder* theElectronMCTruthFinder_ = new ElectronMCTruthFinder();

    std::vector<SimTrack>::const_iterator iFirstSimTk = theSimTracks.begin();
    std::vector<ElectronMCTruth> MCElectrons=theElectronMCTruthFinder_->find (theSimTracks,  theSimVertices);
    
    mcTruth(gpH,MCElectrons);

  }
  
  
  if (saveReco) {
      //photons
      edm::Handle<reco::PhotonCollection>  phoH;
      event.getByLabel("mustachePhotons", phoH);


      edm::Handle<reco::TrackCollection> tracks;
      event.getByLabel("generalTracks",tracks);

      //topology

      // get geometry
      edm::ESHandle<CaloGeometry> geoHandle;
      //   iSetup.get<IdealGeometryRecord>().get(geoHandle);
      iSetup.get<CaloGeometryRecord>().get(geoHandle);
      geometry = geoHandle.product();


      edm::ESHandle<CaloTopology> pTopology;
      iSetup.get<CaloTopologyRecord>().get(pTopology);
      topology = pTopology.product();
      
      edm::Handle<EBRecHitCollection> ecalhitseb;
      const EBRecHitCollection* rhitseb=0;
      event.getByLabel("ecalRecHit","EcalRecHitsEB", ecalhitseb);
      rhitseb = ecalhitseb.product(); // get a ptr to the product

      edm::Handle<EERecHitCollection> ecalhitsee;
      const EERecHitCollection* rhitsee=0;
      event.getByLabel("ecalRecHit","EcalRecHitsEK", ecalhitsee);
      rhitsee = ecalhitsee.product(); // get a ptr to the product

      
      edm::Handle<reco::PFCandidateCollection>  PFCandidates;
      event.getByLabel("particleFlow", PFCandidates);

      clearVector();

      //rechit are too much, just save them for the first 100 events
      if(entryNumber<101)recHitReco(rhitsee);

      //PHOTONS
      phoReco(phoH,tracks,rhitseb,rhitsee,PFCandidates);

      //electrons
      edm::Handle<reco::GsfElectronCollection>  ElectronHandle;
      event.getByLabel("gsfElectrons", ElectronHandle);


      edm::Handle<reco::ConversionCollection> hConversions;
      event.getByLabel("allConversions", hConversions);

      edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
      event.getByLabel("offlineBeamSpot",recoBeamSpotHandle);

 
      //ELECTRONS
      eleReco(ElectronHandle,hConversions,recoBeamSpotHandle,rhitseb,rhitsee,PFCandidates);
     

      //pfsuperclusters
      edm::Handle<reco::SuperClusterCollection> superClustersEBHandle;
      event.getByLabel("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel",superClustersEBHandle);

      edm::Handle<reco::SuperClusterCollection> superClustersEEHandle;
      event.getByLabel("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower",superClustersEEHandle);
      //PFSUPERCLUSTERS
      scReco(superClustersEBHandle,superClustersEEHandle);


      //multi5x5SuperClusters
      edm::Handle<reco::SuperClusterCollection> multi5x5Handle;
      event.getByLabel("multi5x5SuperClustersWithPreshower",multi5x5Handle);

      //PFSUPERCLUSTERS
      multi5x5scReco(multi5x5Handle);


      //hybridSuperClusters
      edm::Handle<reco::SuperClusterCollection> hybridHandle;
      event.getByLabel("hybridSuperClusters",hybridHandle);

      //PFSUPERCLUSTERS
      hybridscReco(hybridHandle);

      
  }
  t->Fill();
}
 
void dumper::beginJob() {
  f = new TFile(outputFileName.c_str(), "recreate");
  t = new TTree("tree", "tree");
  
  
  t->Branch("nvtx", &nvtx, "nvtx/I");
  t->Branch("vertexx", &vertex_x, "vertexx[nvtx]/F");
  t->Branch("vertexy", &vertex_y, "vertexy[nvtx]/F");
  t->Branch("vertexz", &vertex_z, "vertexz[nvtx]/F");
  t->Branch("rho",  &rho, "rho/F");
  

  if (saveReco) {

    t->Branch("phon",   &pho_n,   "phon/I");
    t->Branch("phopt",  &pho_pt,  "phopt[phon]/F");
    t->Branch("phoeta", &pho_eta, "phoeta[phon]/F");
    t->Branch("phophi", &pho_phi, "phophi[phon]/F");
    t->Branch("phoe", &pho_e, "phoe[phon]/F");
    t->Branch("phoesc", &pho_esc, "phoesc[phon]/F");
    t->Branch("phoetasc", &pho_etasc, "phoetasc[phon]/F");
    t->Branch("phophisc", &pho_phisc, "phophisc[phon]/F");

    t->Branch("phoxsc", &pho_xsc, "phoxsc[phon]/F");
    t->Branch("phoysc", &pho_ysc, "phoysc[phon]/F");
    t->Branch("phozsc", &pho_zsc, "phozsc[phon]/F");

    t->Branch("phoxcalo", &pho_xcalo, "phoxcalo[phon]/F");
    t->Branch("phoycalo", &pho_ycalo, "phoycalo[phon]/F");
    t->Branch("phozcalo", &pho_zcalo, "phozcalo[phon]/F");

    t->Branch("phohasPixelSeed", &pho_hasPixelSeed, "phohasPixelSeed[phon]/F");
    t->Branch("phoisEB", &pho_isEB, "phoisEB[phon]/F");
    t->Branch("phoisEE", &pho_isEE, "phoisEE[phon]/F");
    t->Branch("phoisEBEEGap", &pho_isEBEEGap, "phoisEBEEGap[phon]/F");
    t->Branch("phoR9", &pho_R9, "phoR9[phon]/F");
    t->Branch("phoE9", &pho_E9, "phoE9[phon]/F");
    t->Branch("phoE25", &pho_E25, "phoE25[phon]/F");

    t->Branch("phojurECAL", &pho_jurECAL, "phojurECAL[phon]/F");	
    t->Branch("photwrHCAL", &pho_twrHCAL, "photwrHCAL[phon]/F");	
    t->Branch("phoHoverE", &pho_HoverE, "phoHoverE[phon]/F");
    t->Branch("phohlwTrack", &pho_hlwTrack, "phohlwTrack[phon]/F");
    t->Branch("phosiEtaiEtaZS", &pho_siEtaiEtaZS, "phosiEtaiEtaZS[phon]/F");
    t->Branch("phosiEtaiEtaNoZS", &pho_siEtaiEtaNoZS, "phosiEtaiEtaNoZS[phon]/F");
    t->Branch("phosiEtaiEtaWrong", &pho_siEtaiEtaWrong, "phosiEtaiEtaWrong[phon]/F");

    t->Branch("phosiPhiiPhiZS", &pho_siPhiiPhiZS, "phosiPhiiPhiZS[phon]/F");
    t->Branch("phosiPhiiPhiNoZS", &pho_siPhiiPhiNoZS, "phosiPhiiPhiNoZS[phon]/F");
    t->Branch("phosiPhiiPhiWrong", &pho_siPhiiPhiWrong, "phosiPhiiPhiWrong[phon]/F");


    t->Branch("phosMajZS", &pho_sMajZS, "phosMajZS[phon]/F");
    t->Branch("phosMajNoZS", &pho_sMajNoZS, "phosMajNoZS[phon]/F");
    t->Branch("phosMinZS", &pho_sMinZS, "phosMinZS[phon]/F");
    t->Branch("phosMinNoZS", &pho_sMinNoZS, "phosMinNoZS[phon]/F");
    t->Branch("phoalphaZS", &pho_alphaZS, "phoalphaZS[phon]/F");
    t->Branch("phoalphaNoZS", &pho_alphaNoZS, "phoalphaNoZS[phon]/F");


    t->Branch("phoscetawid",&pho_scetawid,"phoscetawid[phon]/F");
    t->Branch("phoscphiwid",&pho_scphiwid,"phoscphiwid[phon]/F");
    t->Branch("phojurECAL03",&pho_jurECAL03,"phojurECAL03[phon]/F");
    t->Branch("photwrHCAL03",&pho_twrHCAL03,"photwrHCAL03[phon]/F");
    t->Branch("phohlwTrack03",&pho_hlwTrack03,"phohlwTrack03[phon]/F");

    t->Branch("phopfCandPt",&pho_pfCandPt,"phopfCandPt[20][200]/F");//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
    t->Branch("phopfCandEta",&pho_pfCandEta,"phopfCandEta[20][200]/F");
    t->Branch("phopfCandPhi",&pho_pfCandPhi,"phopfCandPhi[20][200]/F");
    t->Branch("phopfCandVtx",&pho_pfCandVtx,"phopfCandVtx[20][200]/F");
    t->Branch("phopfCandType",&pho_pfCandType,"phopfCandType[20][200]/F");

    t->Branch("phoptiso004", &pho_ptiso004, "phoptiso004[phon]/F");
    t->Branch("phontrkiso004", &pho_ntrkiso004, "phontrkiso004[phon]/I");
    t->Branch("phoptiso035", &pho_ptiso035, "phoptiso035[phon]/F");
    t->Branch("phontrkiso035", &pho_ntrkiso035, "phontrkiso035[phon]/I");
    t->Branch("phoptiso04", &pho_ptiso04, "phoptiso04[phon]/F");
    t->Branch("phontrkiso04", &pho_ntrkiso04, "phontrkiso04[phon]/I");

    t->Branch("phoPfSumChargedHadronPt",&pho_pfSumChargedHadronPt,"phoPfSumChargedHadronPt[phon]/F");
    t->Branch("phoPfSumNeutralHadronEt",&pho_pfsumNeutralHadronEt,"phoPfSumNeutralHadronEt[phon]/F");
    t->Branch("phoPfSumPhotonEt",&pho_pfsumPhotonEt,"phoPfSumPhotonEt[phon]/F");

    t->Branch("pfSCn",   &pfSC_n,   "pfSCn/I");
    t->Branch("pfSCeta", &pfSC_eta, "pfSCeta[pfSCn]/F");
    t->Branch("pfSCphi", &pfSC_phi, "pfSCphi[pfSCn]/F");
    t->Branch("pfSCe", &pfSC_e, "pfSCe[pfSCn]/F");
    t->Branch("pfSCnBC", &pfSC_nBC, "pfSCnBC[pfSCn]/I");
    t->Branch("pfSCnXtalsSeed", &pfSC_nXtalsSeed, "pfSCnXtalsSeed[pfSCn]/I");

    t->Branch("pfSCbcEta", &pfSC_bcEta, "pfSCbcEta[200][200]/F");//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
    t->Branch("pfSCbcPhi", &pfSC_bcPhi, "pfSCbcPhi[200][200]/F");
    t->Branch("pfSCbcE", &pfSC_bcE, "pfSCbcE[200][200]/F");
    t->Branch("pfSCbcNXtals", &pfSC_bcNXtals, "pfSCbcNXtals[200][200]/I");

    t->Branch("multi5x5SCn",   &multi5x5SC_n,   "multi5x5SCn/I");
    t->Branch("multi5x5SCeta", &multi5x5SC_eta, "multi5x5SCeta[multi5x5SCn]/F");
    t->Branch("multi5x5SCphi", &multi5x5SC_phi, "multi5x5SCphi[multi5x5SCn]/F");
    t->Branch("multi5x5SCe", &multi5x5SC_e, "multi5x5SCe[multi5x5SCn]/F");
    t->Branch("multi5x5SCnBC", &multi5x5SC_nBC, "multi5x5SCnBC[multi5x5SCn]/F");
    t->Branch("multi5x5SCnXtalsSeed", &multi5x5SC_nXtalsSeed, "multi5x5SCnXtalsSeed[multi5x5SCn]/I");

    t->Branch("multi5x5SCbcEta", &multi5x5SC_bcEta, "multi5x5SCbcEta[200][200]/F");//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
    t->Branch("multi5x5SCbcPhi", &multi5x5SC_bcPhi, "multi5x5SCbcPhi[200][200]/F");
    t->Branch("multi5x5SCbcE", &multi5x5SC_bcE, "multi5x5SCbcE[200][200]/F");
    t->Branch("multi5x5SCbcNXtals", &multi5x5SC_bcNXtals, "multi5x5SCbcNXtals[200][200]/I");

    t->Branch("hybridSCn",   &hybridSC_n,   "hybridSCn/I");
    t->Branch("hybridSCeta", &hybridSC_eta, "hybridSCeta[hybridSCn]/F");
    t->Branch("hybridSCphi", &hybridSC_phi, "hybridSCphi[hybridSCn]/F");
    t->Branch("hybridSCe", &hybridSC_e, "hybridSCe[hybridSCn]/F");
    t->Branch("hybridSCnBC", &hybridSC_nBC, "hybridSCnBC[hybridSCn]/F");
    t->Branch("hybridSCnXtalsSeed", &hybridSC_nXtalsSeed, "hybridSCnXtalsSeed[hybridSCn]/I");

    t->Branch("hybridSCbcEta", &hybridSC_bcEta, "hybridSCbcEta[200][200]/F");//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
    t->Branch("hybridSCbcPhi", &hybridSC_bcPhi, "hybridSCbcPhi[200][200]/F");
    t->Branch("hybridSCbcE", &hybridSC_bcE, "hybridSCbcE[200][200]/F");
    t->Branch("hybridSCbcNXtals", &hybridSC_bcNXtals, "hybridSCbcNXtals[200][200]/I");
      
    t->Branch("elen", &ele_n, "elen/I");
    t->Branch("elepx",&ele_px,"elepx[elen]/F");
    t->Branch("elepy",&ele_py,"elepy[elen]/F");
    t->Branch("elepz",&ele_pz,"elepz[elen]/F");
    t->Branch("elevx",&ele_vx,"elevx[elen]/F");
    t->Branch("elevy",&ele_vy,"elevy[elen]/F");
    t->Branch("elevz",&ele_vz,"elevz[elen]/F");
    t->Branch("elept",  &ele_pt,  "elept[elen]/F");
    t->Branch("eleeta",&ele_eta,"eleeta[elen]/F");
    t->Branch("elephi",&ele_phi,"elephi[elen]/F");
    t->Branch("elee",&ele_e,"elee[elen]/F");
    t->Branch("eleecalEnergy",&ele_ecalEnergy,"eleecalEnergy[elen]/F");     
    t->Branch("eletrackPatVtx",&ele_trackPatVtx,"eletrackPatVtx[elen]/F");  
    t->Branch("eletrackAtVtxPt",&ele_trackAtVtxPt,"eletrackAtVtxPt[elen]/F");  
    t->Branch("eletrackAtVtxEta",&ele_trackAtVtxEta,"eletrackAtVtxEta[elen]/F");  
    t->Branch("eletrackAtVtxPhi",&ele_trackAtVtxPhi,"eletrackAtVtxPhi[elen]/F");  
    t->Branch("eletrackAtCaloPt",&ele_trackAtCaloPt,"eletrackAtCaloPt[elen]/F");  
    t->Branch("eletrackAtCaloEta",&ele_trackAtCaloEta,"eletrackAtCaloEta[elen]/F");  
    t->Branch("eletrackAtCaloPhi",&ele_trackAtCaloPhi,"eletrackAtCaloPhi[elen]/F");  

    t->Branch("eleesc", &ele_esc, "eleesc[elen]/F");
    t->Branch("eleetasc", &ele_etasc, "eleetasc[elen]/F");
    t->Branch("elephisc", &ele_phisc, "elephisc[elen]/F");

    t->Branch("elexsc", &ele_xsc, "elexsc[elen]/F");
    t->Branch("eleysc", &ele_ysc, "eleysc[elen]/F");
    t->Branch("elezsc", &ele_zsc, "elezsc[elen]/F");

    t->Branch("elexcalo", &ele_xcalo, "elexcalo[elen]/F");
    t->Branch("eleycalo", &ele_ycalo, "eleycalo[elen]/F");
    t->Branch("elezcalo", &ele_zcalo, "elezcalo[elen]/F");

    t->Branch("eleisEB", &ele_isEB, "eleisEB[elen]/F");
    t->Branch("eleisEE", &ele_isEE, "eleisEE[elen]/F");
    t->Branch("eleisEBEEGap", &ele_isEBEEGap, "eleisEBEEGap[elen]/F");
    t->Branch("eleE9", &ele_E9, "eleE9[elen]/F");
    t->Branch("eleE25", &ele_E25, "eleE25[elen]/F");




    t->Branch("elecharge",&ele_charge,"elecharge[elen]/F");
    t->Branch("eleFbrem",&ele_fbrem,"eleFbrem[elen]/F");
    t->Branch("eleScFbrem",&elesc_fbrem,"eleScfbrem[elen]/F");
    t->Branch("elePfScFbrem",&elepfsc_fbrem,"elePfScFbrem[elen]/F");
    t->Branch("eleESeedClusterOverPout",&ele_eSeedClusterOverPout,"eleESeedClusterOverPout[elen]/F");
    t->Branch("eledist",&ele_dist,"eledist[elen]/F");
    t->Branch("eledcot",&ele_dcot,"eledcot[elen]/F");
    t->Branch("elemisHits",&ele_misHits,"elemisHits[elen]/I");
    t->Branch("elematchedConv",&ele_matchedConv,"elematchedConv[elen]/I");  
    t->Branch("eleseedType",&ele_seedType,"eleseedType[elen]/I");
    t->Branch("eleEoP",&ele_EoP,"eleEoP[elen]/F");
    t->Branch("eleOneOverEMinusOneOverP",&ele_OneOverEMinusOneOverP,"eleOneOverEMinusOneOverP[elen]/F");
    t->Branch("eler9",&ele_r9,"eler9[elen]/F");
    t->Branch("elenSubClusters",&ele_nSubClusters,"elenSubClusters[elen]/I");
    t->Branch("eletrkIso",&ele_trkIso,"eletrkIso[elen]/F");
    t->Branch("eleecalIso",&ele_ecalIso,"eleecalIso[elen]/F");
    t->Branch("elehcalIso",&ele_hcalIso,"elehcalIso[elen]/F");
    t->Branch("eletrkIso03",&ele_trkIso03,"eletrkIso03[elen]/F");
    t->Branch("eleecalIso03",&ele_ecalIso03,"eleecalIso03[elen]/F");
    t->Branch("elehcalIso03",&ele_hcalIso03,"elehcalIso03[elen]/F");
    t->Branch("elesiEtaiEtaZS", &ele_siEtaiEtaZS, "elesiEtaiEtaZS[elen]/F");
    t->Branch("elesiEtaiEtaNoZS", &ele_siEtaiEtaNoZS, "elesiEtaiEtaNoZS[elen]/F");
    t->Branch("elesiEtaiEtaWrong", &ele_siEtaiEtaWrong, "elesiEtaiEtaWrong[elen]/F");

    t->Branch("elesiPhiiPhiZS", &ele_siPhiiPhiZS, "elesiPhiiPhiZS[elen]/F");
    t->Branch("elesiPhiiPhiNoZS", &ele_siPhiiPhiNoZS, "elesiPhiiPhiNoZS[elen]/F");
    t->Branch("elesiPhiiPhiWrong", &ele_siPhiiPhiWrong, "elesiPhiiPhiWrong[elen]/F");


    t->Branch("elesMajZS", &ele_sMajZS, "elesMajZS[elen]/F");
    t->Branch("elesMajNoZS", &ele_sMajNoZS, "elesMajNoZS[elen]/F");
    t->Branch("elesMinZS", &ele_sMinZS, "elesMinZS[elen]/F");
    t->Branch("elesMinNoZS", &ele_sMinNoZS, "elesMinNoZS[elen]/F");
    t->Branch("elealphaZS", &ele_alphaZS, "elealphaZS[elen]/F");
    t->Branch("elealphaNoZS", &ele_alphaNoZS, "elealphaNoZS[elen]/F");

    t->Branch("elepfCandPt",&ele_pfCandPt,"elepfCandPt[20][200]/F");//note:for bidimensional arrays in root the dimension is hardcoded. check if you change a maxdim used by a 2d array
    t->Branch("elepfCandEta",&ele_pfCandEta,"elepfCandEta[20][200]/F");
    t->Branch("elepfCandPhi",&ele_pfCandPhi,"elepfCandPhi[20][200]/F");
    t->Branch("elepfCandVtx",&ele_pfCandVtx,"elepfCandVtx[20][200]/F");
    t->Branch("elepfCandType",&ele_pfCandType,"elepfCandType[20][200]/F");



    t->Branch("eledEtaIn",&ele_dEtaIn,"eledEtaIn[elen]/F");
    t->Branch("eledPhiIn",&ele_dPhiIn,"eledPhiIn[elen]/F");
    t->Branch("eleHoE",&ele_HoE,"eleHoE[elen]/F");
    t->Branch("elepFlowMVA",&ele_pFlowMVA,"elepFlowMVA[elen]/F");
    t->Branch("eleScEnergy",&ele_sc_energy,"eleScEnergy[elen]/F");
    t->Branch("eleScEta",&ele_sc_eta,"eleScEta[elen]/F");
    t->Branch("eleScPhi",&ele_sc_phi,"eleScPhi[elen]/F");
    t->Branch("elePfSumChargedHadronPt",&ele_pfSumChargedHadronPt,"elePfSumChargedHadronPt[elen]/F");
    t->Branch("elePfSumNeutralHadronEt",&ele_pfsumNeutralHadronEt,"elePfSumNeutralHadronEt[elen]/F");
    t->Branch("elePfSumPhotonEt",&ele_pfsumPhotonEt,"elePfSumPhotonEt[elen]/F");

    t->Branch("rechitn",   &rechit_n,   "rechitn/I");
    t->Branch("rechite",  &rechit_e,  "rechite[rechitn]/F");
    t->Branch("rechitix", &rechit_ix, "rechitix[rechitn]/F");
    t->Branch("rechitiy", &rechit_iy, "rechitiy[rechitn]/F");


  }   
      
  if (!isData) {
    t->Branch("gpn",   &gp_n,   "gpn/I");
    t->Branch("gppt",  &gp_pt,  "gppt[gpn]/F");
    t->Branch("gpeta", &gp_eta, "gpeta[gpn]/F");
    t->Branch("gpphi", &gp_phi, "gpphi[gpn]/F");
    t->Branch("gpidMC", &gp_idMC, "gpidMC[gpn]/I");
    t->Branch("gpstatusMC", &gp_statusMC, "gpstatusMC[gpn]/I");
    t->Branch("gpmotherIdMC", &gp_motherIdMC, "gpmotherIdMC[gpn]/I");

    t->Branch("gphon",   &gpho_n,   "gphon/I");
    t->Branch("gphopt",  &gpho_pt,  "gphopt[gphon]/F");
    t->Branch("gphoeta", &gpho_eta, "gphoeta[gphon]/F");
    t->Branch("gphophi", &gpho_phi, "gphophi[gphon]/F");
    t->Branch("gphoindex", &gpho_index, "gphoindex[gphon]/I");

    t->Branch("gelen",   &gele_n,   "gelen/I");
    t->Branch("gelept",  &gele_pt,  "gelept[gelen]/F");
    t->Branch("geleeta", &gele_eta, "geleeta[gelen]/F");
    t->Branch("gelephi", &gele_phi, "gelephi[gelen]/F");
    t->Branch("gelefbrem80", &gele_fbrem80, "gelefbrem80[gelen]/F");
    t->Branch("gelefbrem120", &gele_fbrem120, "gelefbrem120[gelen]/F");
    t->Branch("geleindex", &gele_index, "geleindex[gelen]/I");
    
    t->Branch("truePU", &truePU, "truePU/F");
    t->Branch("bxPU", &bxPU, "bxPU[16]/I");


  }
}

void dumper::endJob() {
  f->cd();
  t->Write();
  f->Close();
}

void dumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(dumper);
