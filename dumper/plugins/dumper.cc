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

#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/FwdPtr.h"

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

class dumper : public edm::EDAnalyzer {
public:
  explicit dumper(const edm::ParameterSet&);
  ~dumper();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void mcTruth(edm::Handle<reco::GenParticleCollection> genParticleH,std::vector<ElectronMCTruth> MCElectrons);
  void phoReco(edm::Handle<reco::PhotonCollection> photonH, edm::Handle<reco::TrackCollection> traH);
  void eleReco(edm::Handle<reco::GsfElectronCollection> ElectronHandle,edm::Handle<reco::ConversionCollection> hConversions, edm::Handle<reco::BeamSpot> recoBeamSpotHandle);
  void scReco(edm::Handle<reco::SuperClusterCollection> superClustersEBHandle, edm::Handle<reco::SuperClusterCollection> superClustersEEHandle);

private:
  virtual void beginRun(const edm::Run& run,const edm::EventSetup& setup);
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::InputTag hitEBLabel_, hitEELabel_;

  std::string outputFileName;
  TFile* f;
  TTree* t;

  const CaloTopology *topology;
  const CaloTopology *geometry;
  
  Float_t rho;
  Int_t n, npf, gp_n, gpho_n, gele_n, pho_n,ele_n, pfSC_n;

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
  Float_t gele_fbrem[MAXPHOTONSTOSAVE];
  Int_t gele_index[MAXPHOTONSTOSAVE];



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
  Float_t   pho_E9[MAXPHOTONSTOSAVE];
  Float_t   pho_E25[MAXPHOTONSTOSAVE];

  Float_t  pho_jurECAL[MAXPHOTONSTOSAVE];
  Float_t  pho_twrHCAL[MAXPHOTONSTOSAVE];
  Float_t  pho_HoverE[MAXPHOTONSTOSAVE];
  Float_t  pho_hlwTrack[MAXPHOTONSTOSAVE];
  Float_t  pho_etawid[MAXPHOTONSTOSAVE];
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
      
  Float_t pfSC_eta[MAXSCTOSAVE];
  Float_t pfSC_phi[MAXSCTOSAVE];
  Float_t   pfSC_e[MAXSCTOSAVE];
  Float_t pfSC_nBC[MAXSCTOSAVE];
  Float_t pfSC_nXtals[MAXSCTOSAVE];


  Float_t ele_pt[MAXPHOTONSTOSAVE];
  Float_t ele_eta[MAXPHOTONSTOSAVE];
  Float_t ele_phi[MAXPHOTONSTOSAVE];

  Float_t    ele_energy[MAXPHOTONSTOSAVE];
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
  Float_t  ele_SigmaIetaIeta[MAXPHOTONSTOSAVE];
  Float_t  ele_SigmaIphiIphi[MAXPHOTONSTOSAVE];
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
  


  float truePU;
  int bxPU[16];
  Int_t nvtx;

  bool isData;
  bool saveReco;


  TH1F* timeEB, *timeEE;
  TH2F* timeEB2D, *timeEE2D;
};

void dumper::beginRun(const edm::Run& run,const edm::EventSetup& setup) {
  std::cout <<"begining run "<<std::endl;
}


dumper::dumper(const edm::ParameterSet& iConfig) {
  outputFileName  = iConfig.getParameter<std::string>("OutputFileName");
  isData          = iConfig.getParameter<bool>("isData");
  saveReco        = iConfig.getParameter<bool>("saveReco");
  hitEBLabel_     = iConfig.getParameter<edm::InputTag>("hitEBLabel");
  hitEELabel_     = iConfig.getParameter<edm::InputTag>("hitEELabel");

  timeEB   = new TH1F("timeEB",   "", 400, -100, 100);
  timeEE   = new TH1F("timeEE",   "", 400, -100, 100);
  timeEB2D = new TH2F("timeEB2D", "", 400, -100, 100, 100, 0, 100);
  timeEE2D = new TH2F("timeEE2D", "", 400, -100, 100, 100, 0, 100);
}

dumper::~dumper() 
{}




void dumper::phoReco(edm::Handle<reco::PhotonCollection> phoH, edm::Handle<reco::TrackCollection> tracksH) {

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
      pho_E9[pho_n] = pho->e3x3();
      pho_E25[pho_n] = pho->e5x5();


      //default photon Variable
      pho_jurECAL[pho_n] = pho->ecalRecHitSumEtConeDR04();//isolationEcalRecHit
      pho_twrHCAL[pho_n] = pho->hcalTowerSumEtConeDR04();//isolationHcalRecHit
      pho_HoverE[pho_n] = pho->hadronicOverEm();
      pho_hlwTrack[pho_n] = pho->trkSumPtHollowConeDR04();//isolationHollowTrkCone
      pho_etawid[pho_n] = pho->sigmaIetaIeta();//sqrt(pho->covEtaEta());
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
      



      pho_n++;
    }

  }

}

void dumper::eleReco(edm::Handle<reco::GsfElectronCollection> ElectronHandle, edm::Handle<reco::ConversionCollection> hConversions,edm::Handle<reco::BeamSpot> recoBeamSpotHandle){

  ele_n=0;  

  unsigned index_gsf = 0;
  for (reco::GsfElectronCollection::const_iterator itElectron = ElectronHandle->begin();
       itElectron != ElectronHandle->end(); ++itElectron, ++index_gsf) {
     
    reco::GsfElectronRef eleRef = reco::GsfElectronRef(ElectronHandle, index_gsf);

    if (itElectron->pt()<2.5||ele_n>=MAXPHOTONSTOSAVE)
      continue;


    ele_pt[ele_n]          = itElectron->pt();
    ele_energy[ele_n]      = itElectron->energy();
    ele_ecalEnergy[ele_n]  = itElectron->ecalEnergy();              
    ele_trackPatVtx[ele_n] = itElectron->trackMomentumAtVtx().R();


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
    ele_r9[ele_n] = itElectron->r9();
    ele_misHits[ele_n]       = itElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits();

    ele_dist[ele_n] = itElectron->convDist();
    ele_dcot[ele_n] = itElectron->convDcot();


    bool matchesConv = ConversionTools::hasMatchedConversion(*eleRef,hConversions,recoBeamSpotHandle->position());
     ele_matchedConv[ele_n] = matchesConv;
     ele_seedType[ele_n]       = itElectron->ecalDrivenSeed()+2*itElectron->trackerDrivenSeed();
     ele_nSubClusters[ele_n]       = itElectron->numberOfBrems();
     ele_HoE[ele_n]          = itElectron->hadronicOverEm();
     ele_pFlowMVA[ele_n]          = itElectron->mvaOutput().mva;
     ele_SigmaIetaIeta[ele_n]= itElectron->sigmaIetaIeta();
     ele_SigmaIphiIphi[ele_n]= itElectron->sigmaIphiIphi();
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

     ele_n++;
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
      pfSC_nXtals[pfSC_n] = itSC->seed()->size();

    }
    pfSC_n++;
  }

  for (reco::SuperClusterCollection::const_iterator itSC = superClustersEEHandle->begin();
       itSC != superClustersEEHandle->end(); ++itSC) {
    if (itSC->energy() > 0. && pfSC_n<MAXSCTOSAVE) {
      

      pfSC_eta[pfSC_n] = itSC->eta();
      pfSC_phi[pfSC_n] = itSC->phi();

      pfSC_e[pfSC_n] = itSC->energy();
      pfSC_nBC[pfSC_n] = itSC->clustersSize();
      pfSC_nXtals[pfSC_n] = itSC->seed()->size();

    }

    pfSC_n++;
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
    if (gp->pt() > 5.) {//loose pt cut just to save space
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

	      float totBrem=0 ;
	      unsigned int iBrem ;
	      for ( iBrem=0; iBrem < (*iEl).bremVertices().size(); ++iBrem ) {

		float rBrem= (*iEl).bremVertices()[iBrem].perp();
		if ( rBrem < 120 ) {
		  totBrem +=  (*iEl).bremMomentum()[iBrem].e();   
		}
	      }
	      gele_fbrem[gele_n]=totBrem/((*iEl).fourMomentum().e());
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


   // rho from fast jet
   edm::Handle<double> rhoH;
   //iEvent.getByLabel(edm::InputTag("kt6PFJets","rho","Iso"),rhoH); 
   if( event.getByLabel(edm::InputTag("kt6PFJets","rho"),rhoH) )
     rho = *rhoH;
   else 
     rho = 0;

  
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
    
    edm::Handle<std::vector<PileupSummaryInfo>> puH;
    event.getByLabel("addPileupInfo", puH);
    truePU = (*puH)[0].getTrueNumInteractions();
    for (unsigned int j=0; j<puH->size(); j++)
      bxPU[j] = (*puH)[j].getPU_NumInteractions();
  }
  
  
  if (saveReco) {
      //photons
      edm::Handle<reco::PhotonCollection>  phoH;
      event.getByLabel("mustachePhotons", phoH);


      edm::Handle<reco::TrackCollection> tracks;
      event.getByLabel("generalTracks",tracks);
      
      //PHOTONS
      phoReco(phoH,tracks);

      //electrons
      edm::Handle<reco::GsfElectronCollection>  ElectronHandle;
      event.getByLabel("gsfElectrons", ElectronHandle);


      edm::Handle<reco::ConversionCollection> hConversions;
      event.getByLabel("allConversions", hConversions);

      edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
      event.getByLabel("offlineBeamSpot",recoBeamSpotHandle);

 
      //ELECTRONS
      eleReco(ElectronHandle,hConversions,recoBeamSpotHandle);
     
      //superclusters
      edm::Handle<reco::SuperClusterCollection> superClustersEBHandle;
      event.getByLabel("particleFlowSuperClusterECAL","particleFlowSuperClusterECALBarrel",superClustersEBHandle);

      edm::Handle<reco::SuperClusterCollection> superClustersEEHandle;
      event.getByLabel("particleFlowSuperClusterECAL","particleFlowSuperClusterECALEndcapWithPreshower",superClustersEEHandle);
      //SUPERCLUSTERS
      scReco(superClustersEBHandle,superClustersEEHandle);


      
      
  }
  t->Fill();
}
 
void dumper::beginJob() {
  f = new TFile(outputFileName.c_str(), "recreate");
  t = new TTree("tree", "tree");
  
  
  t->Branch("nvtx", &nvtx, "nvtx/I");
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
    t->Branch("phoE9", &pho_E9, "phoE9[phon]/F");
    t->Branch("phoE25", &pho_E25, "phoE25[phon]/F");

    t->Branch("phojurECAL", &pho_jurECAL, "phojurECAL[phon]/F");	
    t->Branch("photwrHCAL", &pho_twrHCAL, "photwrHCAL[phon]/F");	
    t->Branch("phoHoverE", &pho_HoverE, "phoHoverE[phon]/F");
    t->Branch("phohlwTrack", &pho_hlwTrack, "phohlwTrack[phon]/F");
    t->Branch("phoetawid", &pho_etawid, "phoetawid[phon]/F");
    t->Branch("phoscetawid",&pho_scetawid,"phoscetawid[phon]/F");
    t->Branch("phoscphiwid",&pho_scphiwid,"phoscphiwid[phon]/F");
    t->Branch("phojurECAL03",&pho_jurECAL03,"phojurECAL03[phon]/F");
    t->Branch("photwrHCAL03",&pho_twrHCAL03,"photwrHCAL03[phon]/F");
    t->Branch("phohlwTrack03",&pho_hlwTrack03,"phohlwTrack03[phon]/F");



    t->Branch("phoptiso004", &pho_ptiso004, "phoptiso004[phon]/F");
    t->Branch("phontrkiso004", &pho_ntrkiso004, "phontrkiso004[phon]/I");
    t->Branch("phoptiso035", &pho_ptiso035, "phoptiso035[phon]/F");
    t->Branch("phontrkiso035", &pho_ntrkiso035, "phontrkiso035[phon]/I");
    t->Branch("phoptiso04", &pho_ptiso04, "phoptiso04[phon]/F");
    t->Branch("phontrkiso04", &pho_ntrkiso04, "phontrkiso04[phon]/I");

    t->Branch("pfSCn",   &pfSC_n,   "pfSCn/I");
    t->Branch("pfSCeta", &pfSC_eta, "pfSCeta[pfSCn]/F");
    t->Branch("pfSCphi", &pfSC_phi, "pfSCphi[pfSCn]/F");
    t->Branch("pfSCe", &pfSC_e, "pfSCe[pfSCn]/F");
    t->Branch("pfSCnBC", &pfSC_nBC, "pfSCnBC[pfSCn]/F");
    t->Branch("pfSCnXtals", &pfSC_nXtals, "pfSCnXtals[pfSCn]/F");
      
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
    t->Branch("eleenergy",&ele_energy,"eleenergy[elen]/F");
    t->Branch("eleecalEnergy",&ele_ecalEnergy,"eleecalEnergy[elen]/F");     
    t->Branch("eletrackPatVtx",&ele_trackPatVtx,"eletrackPatVtx[elen]/F");  
    t->Branch("eletrackAtVtxPt",&ele_trackAtVtxPt,"eletrackAtVtxPt[elen]/F");  
    t->Branch("eletrackAtVtxEta",&ele_trackAtVtxEta,"eletrackAtVtxEta[elen]/F");  
    t->Branch("eletrackAtVtxPhi",&ele_trackAtVtxPhi,"eletrackAtVtxPhi[elen]/F");  
    t->Branch("eletrackAtCaloPt",&ele_trackAtCaloPt,"eletrackAtCaloPt[elen]/F");  
    t->Branch("eletrackAtCaloEta",&ele_trackAtCaloEta,"eletrackAtCaloEta[elen]/F");  
    t->Branch("eletrackAtCaloPhi",&ele_trackAtCaloPhi,"eletrackAtCaloPhi[elen]/F");  


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
    t->Branch("eleSigmaIetaIeta",&ele_SigmaIetaIeta,"eleSigmaIetaIeta[elen]/F");
    t->Branch("eleSigmaIphiIphi",&ele_SigmaIphiIphi,"eleSigmaIphiIphi[elen]/F");
    t->Branch("eledEtaIn",&ele_dEtaIn,"eledEtaIn[elen]/F");
    t->Branch("eledPhiIn",&ele_dPhiIn,"eledPhiIn[elen]/F");
    t->Branch("eleHoE",&ele_HoE,"eleHoE[elen]/F");
    t->Branch("elepFlowMVA",&ele_pFlowMVA,"elepFlowMVA[elen]/F");
    t->Branch("eleScEnergy",&ele_sc_energy,"eleScEnergy[elen]/F");
    t->Branch("eleScEta",&ele_sc_eta,"eleScEta[elen]/F");
    t->Branch("eleScPhi",&ele_sc_phi,"eleScPhi[elen]/F");

      
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
    t->Branch("gelefbrem", &gele_fbrem, "gelefbrem[gelen]/F");
    t->Branch("geleindex", &gele_index, "geleindex[gelen]/I");
    
    t->Branch("truePU", &truePU, "truePU/F");
    t->Branch("bxPU", &bxPU, "bxPU[16]/I");
  }
}

void dumper::endJob() {
  f->cd();
  t->Write();
  timeEB->Write();
  timeEE->Write();
  timeEB2D->Write();
  timeEE2D->Write();
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
