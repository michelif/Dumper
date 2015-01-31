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

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <string>
#include <iostream>
#define MAXPARTICLESTOSAVE 100
#define MAXPHOTONSTOSAVE 20

class dumper : public edm::EDAnalyzer {
public:
  explicit dumper(const edm::ParameterSet&);
  ~dumper();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void mcTruth(edm::Handle<reco::GenParticleCollection> genParticleH);
  void phoReco(edm::Handle<reco::PhotonCollection> photonH, edm::Handle<reco::TrackCollection> traH);
  void eleReco(edm::Handle<reco::GsfElectronCollection> ElectronHandle,const EcalRecHitCollection* rhitseb,const EcalRecHitCollection* rhitsee);

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
  Int_t n, npf, gp_n, gpho_n, gele_n, pho_n,ele_n;

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
      
  
  Float_t ele_pt[MAXPHOTONSTOSAVE];
  Float_t ele_eta[MAXPHOTONSTOSAVE];
  Float_t ele_phi[MAXPHOTONSTOSAVE];

  Float_t    electron_pt[MAXPHOTONSTOSAVE];
  Float_t    electron_energy[MAXPHOTONSTOSAVE];
  Float_t    electron_ecalEnergy[MAXPHOTONSTOSAVE];
  Float_t    electron_trackPatVtx[MAXPHOTONSTOSAVE];
  Float_t    electron_px[MAXPHOTONSTOSAVE];
  Float_t    electron_py[MAXPHOTONSTOSAVE];
  Float_t    electron_pz[MAXPHOTONSTOSAVE];
  Float_t    electron_vx[MAXPHOTONSTOSAVE];
  Float_t    electron_vy[MAXPHOTONSTOSAVE];
  Float_t    electron_vz[MAXPHOTONSTOSAVE];
  Float_t    electron_phi[MAXPHOTONSTOSAVE];
  Float_t    electron_eta[MAXPHOTONSTOSAVE];
  Float_t    electron_charge[MAXPHOTONSTOSAVE];
  Float_t    electron_fBrem[MAXPHOTONSTOSAVE];
  Float_t    electron_EoP[MAXPHOTONSTOSAVE];
  Float_t    electron_OneOverEMinusOneOverP[MAXPHOTONSTOSAVE];
  Float_t    electron_r9[MAXPHOTONSTOSAVE];
  Int_t      electron_misHits[MAXPHOTONSTOSAVE];
  Float_t  electron_dist[MAXPHOTONSTOSAVE];
  Float_t  electron_dcot[MAXPHOTONSTOSAVE];

  Int_t  electron_matchedConv[MAXPHOTONSTOSAVE];
  Int_t  electron_seedType[MAXPHOTONSTOSAVE];
  Int_t  electron_nSubClusters[MAXPHOTONSTOSAVE];
  Float_t  electron_HoE[MAXPHOTONSTOSAVE];
  Float_t  electron_pFlowMVA[MAXPHOTONSTOSAVE];
  Float_t  electron_SigmaIetaIeta[MAXPHOTONSTOSAVE];
  Float_t  electron_SigmaIphiIphi[MAXPHOTONSTOSAVE];
  Float_t  electron_trkIso[MAXPHOTONSTOSAVE];
  Float_t  electron_ecalIso[MAXPHOTONSTOSAVE];
  Float_t  electron_hcalIso[MAXPHOTONSTOSAVE];
  Float_t  electron_trkIso03[MAXPHOTONSTOSAVE];
  Float_t  electron_ecalIso03[MAXPHOTONSTOSAVE];
  Float_t  electron_hcalIso03[MAXPHOTONSTOSAVE];
  Float_t  electron_dEtaIn[MAXPHOTONSTOSAVE];
  Float_t  electron_dPhiIn[MAXPHOTONSTOSAVE];
  Float_t  electron_sc_energy[MAXPHOTONSTOSAVE];
  Float_t  electron_sc_eta[MAXPHOTONSTOSAVE];
  Float_t  electron_sc_phi[MAXPHOTONSTOSAVE];
  


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

    if (pho->pt() > 0.) {
      
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

void dumper::eleReco(edm::Handle<reco::GsfElectronCollection> ElectronHandle,const EcalRecHitCollection* rhitseb,const EcalRecHitCollection* rhitsee){

  ele_n=0;  

  unsigned index_gsf = 0;
  for (reco::GsfElectronCollection::const_iterator itElectron = ElectronHandle->begin();
       itElectron != ElectronHandle->end(); ++itElectron, ++index_gsf) {
     
    reco::GsfElectronRef eleRef = reco::GsfElectronRef(ElectronHandle, index_gsf);

    //metti controllo su pt
    //vedi perche'non va
    //    const EcalRecHitCollection* rechits = ( itElectron->isEB()) ? rhitseb : rhitsee;   
    std::cout<<"itelectron works"<<std::endl;     

    //    const edm::Ptr<reco::CaloCluster> theSeed = itElectron->superCluster()->seed();
    //    float e9=EcalClusterTools::e3x3( *theSeed, &(*rechits), topology );
    //    float e9=EcalClusterLazyTools::e3x3( *theSeed);
    //vedi sotto echiedi se va bene
    //vedi perche'non va
    //    float sigIPhiIPhi=sqrt((EcalClusterTools::localCovariances( *theSeed, &(*rechits), topology ))[2]);

    electron_pt[ele_n]          = itElectron->pt();
    std::cout<<electron_pt[ele_n]<<std::endl;
    electron_energy[ele_n]      = itElectron->energy();
    electron_ecalEnergy[ele_n]  = itElectron->ecalEnergy();                 // chiara
    electron_trackPatVtx[ele_n] = itElectron->trackMomentumAtVtx().R();     // chiara
    electron_px[ele_n]       = itElectron->px();
    electron_py[ele_n]       = itElectron->py();
    electron_pz[ele_n]       = itElectron->pz();
    electron_vx[ele_n]       = itElectron->vx();
    electron_vy[ele_n]       = itElectron->vy();
    electron_vz[ele_n]       = itElectron->vz();
    electron_phi[ele_n]      = itElectron->phi();
    electron_eta[ele_n]      = itElectron->eta();
    electron_charge[ele_n]       = itElectron->charge();
    electron_fBrem[ele_n]       = itElectron->fbrem();
    electron_EoP[ele_n]       = itElectron->eSuperClusterOverP();
    electron_OneOverEMinusOneOverP[ele_n]       = (1/itElectron->caloEnergy())-(1/itElectron->trackMomentumAtVtx().R());
    //    electron_r9[ele_n]       = e9/itElectron->superCluster()->rawEnergy(); ridefinita sotto non funziona ecalclustertools chiedi
    electron_r9[ele_n] = itElectron->r9();
    electron_misHits[ele_n]       = itElectron->gsfTrack()->trackerExpectedHitsInner().numberOfHits();

    ele_n++;
  }

}

void dumper::mcTruth(edm::Handle<reco::GenParticleCollection> gpH) {
  
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

  truePU = 9999.;
  for (int i=0; i<16; i++)
    bxPU[i] = 9999;

  edm::Handle<reco::GenParticleCollection> gpH;
  if (!isData) {
    event.getByLabel("genParticles", gpH);
    mcTruth(gpH);
    
    edm::Handle<std::vector<PileupSummaryInfo>> puH;
    event.getByLabel("addPileupInfo", puH);
    truePU = (*puH)[0].getTrueNumInteractions();
    for (unsigned int j=0; j<puH->size(); j++)
      bxPU[j] = (*puH)[j].getPU_NumInteractions();
  }

  
    if (saveReco) {

      edm::Handle<reco::PhotonCollection>  phoH;
      event.getByLabel("photons", phoH);


      edm::Handle<reco::TrackCollection> tracks;
      event.getByLabel("generalTracks",tracks);
      

      phoReco(phoH,tracks);

      edm::Handle<reco::GsfElectronCollection>  ElectronHandle;
      event.getByLabel("gsfElectrons", ElectronHandle);

      edm::Handle<EcalRecHitCollection> ecalhitseb;
      const EcalRecHitCollection* rhitseb=0;
      event.getByLabel("ecalRecHit","EcalRecHitsEB",ecalhitseb);
      rhitseb = ecalhitseb.product(); // get a ptr to the product

      edm::Handle<EcalRecHitCollection> ecalhitsee;
      const EcalRecHitCollection* rhitsee=0;
      event.getByLabel("ecalRecHit","EcalRecHitsEE",ecalhitsee);
      rhitsee = ecalhitsee.product(); // get a ptr to the product


      // get topology
      //      edm::ESHandle<CaloTopology> pTopology;
      //      iSetup.get<CaloTopologyRecord>().get(pTopology);
      //      topology = pTopology.product();
      //per ora commentto chiedi se va bene

      eleReco(ElectronHandle,rhitseb,rhitsee);
     
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
      
    t->Branch("elen",&ele_n,"elen/I");
    t->Branch("elepx",&electron_px,"elepx[elen]/F");
    t->Branch("elepy",&electron_py,"elepy[elen]/F");
    t->Branch("elepz",&electron_pz,"elepz[elen]/F");
    t->Branch("elevx",&electron_vx,"elevx[elen]/F");
    t->Branch("elevy",&electron_vy,"elevy[elen]/F");
    t->Branch("elevz",&electron_vz,"elevz[elen]/F");
    t->Branch("elept",&electron_pt,"elept[elen]/F");
    t->Branch("eleeta",&electron_eta,"eleeta[elen]/F");
    t->Branch("elephi",&electron_phi,"elephi[elen]/F");
    t->Branch("eleenergy",&electron_energy,"eleenergy[elen]/F");
    t->Branch("eleecalEnergy",&electron_ecalEnergy,"eleecalEnergy[elen]/F");     // chiara
    t->Branch("eletrackPatVtx",&electron_trackPatVtx,"eletrackPatVtx[elen]/F");  // chiara
    t->Branch("elecharge",&electron_charge,"elecharge[elen]/F");
    t->Branch("elefBrem",&electron_fBrem,"elefBrem[elen]/F");
    t->Branch("eledist",&electron_dist,"eledist[elen]/F");
    t->Branch("eledcot",&electron_dcot,"eledcot[elen]/F");
    t->Branch("elemisHits",&electron_misHits,"elemisHits[elen]/I");
    t->Branch("elematchedConv",&electron_matchedConv,"elematchedConv[elen]/I");  // chiara
    t->Branch("eleseedType",&electron_seedType,"eleseedType[elen]/I");
    t->Branch("eleEoP",&electron_EoP,"eleEoP[elen]/F");
    t->Branch("eleOneOverEMinusOneOverP",&electron_OneOverEMinusOneOverP,"eleOneOverEMinusOneOverP[elen]/F");
    t->Branch("eler9",&electron_r9,"eler9[elen]/F");
    t->Branch("elenSubClusters",&electron_nSubClusters,"elenSubClusters[elen]/I");
    t->Branch("eletrkIso",&electron_trkIso,"eletrkIso[elen]/F");
    t->Branch("eleecalIso",&electron_ecalIso,"eleecalIso[elen]/F");
    t->Branch("elehcalIso",&electron_hcalIso,"elehcalIso[elen]/F");
    t->Branch("eletrkIso03",&electron_trkIso03,"eletrkIso03[elen]/F");
    t->Branch("eleecalIso03",&electron_ecalIso03,"eleecalIso03[elen]/F");
    t->Branch("elehcalIso03",&electron_hcalIso03,"elehcalIso03[elen]/F");
    t->Branch("eleSigmaIetaIeta",&electron_SigmaIetaIeta,"eleSigmaIetaIeta[elen]/F");
    t->Branch("eleSigmaIphiIphi",&electron_SigmaIphiIphi,"eleSigmaIphiIphi[elen]/F");
    t->Branch("eledEtaIn",&electron_dEtaIn,"eledEtaIn[elen]/F");
    t->Branch("eledPhiIn",&electron_dPhiIn,"eledPhiIn[elen]/F");
    t->Branch("eleHoE",&electron_HoE,"eleHoE[elen]/F");
    t->Branch("elepFlowMVA",&electron_pFlowMVA,"elepFlowMVA[elen]/F");
    t->Branch("elesc_energy",&electron_sc_energy,"elesc_energy[elen]/F");
    t->Branch("elesc_eta",&electron_sc_eta,"elesc_eta[elen]/F");
    t->Branch("elesc_phi",&electron_sc_phi,"elesc_phi[elen]/F");

      
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
