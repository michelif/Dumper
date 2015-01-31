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

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"

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
  int findEleRef(reco::RecoEcalCandidateRef ref, edm::Handle<reco::ElectronCollection> hltEleH);
  void mcTruth(edm::Handle<reco::GenParticleCollection> genParticleH);
  void phoReco(edm::Handle<reco::PhotonCollection> photonH, edm::Handle<reco::TrackCollection> traH);

private:
  virtual void beginRun(const edm::Run& run,const edm::EventSetup& setup);
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::InputTag hitEBLabel_, hitEELabel_;
  edm::InputTag trigResultsTag_;
  HLTConfigProvider hltConfig_; 
  std::vector<std::string> pathNames_;

  std::string outputFileName;
  TFile* f;
  TTree* t;
  
  Float_t rho;
  Int_t passSel[3];
  Int_t n, npf, gp_n, reco_n, nuns, gpho_n, gele_n, pho_n;
  Float_t etawidth[MAXPHOTONSTOSAVE], phiwidth[MAXPHOTONSTOSAVE];
  Float_t etawidthpf[MAXPHOTONSTOSAVE], phiwidthpf[MAXPHOTONSTOSAVE];
  Float_t eraw[MAXPHOTONSTOSAVE];
  Float_t e[MAXPHOTONSTOSAVE];
  Float_t et[MAXPHOTONSTOSAVE];
  Float_t se[MAXPHOTONSTOSAVE];
  Float_t erawpf[MAXPHOTONSTOSAVE];
  Float_t epf[MAXPHOTONSTOSAVE];  
  Float_t eoppf[MAXPHOTONSTOSAVE];  
  Float_t etpf[MAXPHOTONSTOSAVE];  
  Float_t sepf[MAXPHOTONSTOSAVE];  
  Float_t eta[MAXPHOTONSTOSAVE];
  Float_t etapf[MAXPHOTONSTOSAVE];  
  Float_t euns[MAXPHOTONSTOSAVE];
  Float_t etuns[MAXPHOTONSTOSAVE];
  Float_t etauns[MAXPHOTONSTOSAVE];  
  Float_t seta[MAXPHOTONSTOSAVE];
  Float_t setapf[MAXPHOTONSTOSAVE];  
  Float_t phi[MAXPHOTONSTOSAVE];
  Float_t phipf[MAXPHOTONSTOSAVE];  
  Float_t phiuns[MAXPHOTONSTOSAVE];  
  Float_t sphi[MAXPHOTONSTOSAVE];
  Float_t sphipf[MAXPHOTONSTOSAVE];  
  Float_t ecal[MAXPHOTONSTOSAVE];
  Float_t ecalpf[MAXPHOTONSTOSAVE];
  Float_t sieie[MAXPHOTONSTOSAVE];
  Float_t sieiepf[MAXPHOTONSTOSAVE];
  //Float_t eop[MAXPHOTONSTOSAVE];
  //Float_t eoppf[MAXPHOTONSTOSAVE];
  Float_t deta[MAXPHOTONSTOSAVE];
  Float_t detapf[MAXPHOTONSTOSAVE];
  Float_t dphi[MAXPHOTONSTOSAVE];
  Float_t dphipf[MAXPHOTONSTOSAVE];
  Float_t tkptpf[MAXPHOTONSTOSAVE];
  Float_t tketapf[MAXPHOTONSTOSAVE];
  Float_t tkphipf[MAXPHOTONSTOSAVE];
  Float_t tkpt[MAXPHOTONSTOSAVE];
  Float_t tketa[MAXPHOTONSTOSAVE];
  Float_t tkphi[MAXPHOTONSTOSAVE];
  Float_t hoe[MAXPHOTONSTOSAVE];
  Float_t hcal[MAXPHOTONSTOSAVE];
  Float_t tkiso[MAXPHOTONSTOSAVE];
  Float_t hoepf[MAXPHOTONSTOSAVE];
  Float_t hcalpf[MAXPHOTONSTOSAVE];
  Float_t tkisopf[MAXPHOTONSTOSAVE];
  Float_t chiso[MAXPHOTONSTOSAVE];
  Float_t phiso[MAXPHOTONSTOSAVE];
  Float_t neiso[MAXPHOTONSTOSAVE];
  Float_t chi2pf[MAXPHOTONSTOSAVE];
  Int_t mishitspf[MAXPHOTONSTOSAVE];
  Int_t hitspf[MAXPHOTONSTOSAVE];

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


  Float_t reco_et[MAXPHOTONSTOSAVE];
  Float_t reco_pt[MAXPHOTONSTOSAVE];
  Float_t reco_eta[MAXPHOTONSTOSAVE];
  Float_t reco_phi[MAXPHOTONSTOSAVE];

  float truePU;
  int bxPU[16];
  Int_t nvtx;

  bool isData;
  bool saveReco;
  bool saveUnseeded;
  bool newClustering;
  bool oldClustering;

  TH1F* timeEB, *timeEE;
  TH2F* timeEB2D, *timeEE2D;
};

void dumper::beginRun(const edm::Run& run,const edm::EventSetup& setup) {
  std::cout <<"begining run "<<std::endl;
  bool changed = false;
  hltConfig_.init(run,setup,trigResultsTag_.process(),changed); //as we need the orginal HLT config...
  std::cout <<"table name "<<hltConfig_.tableName()<<std::endl;
}


dumper::dumper(const edm::ParameterSet& iConfig) {
  outputFileName  = iConfig.getParameter<std::string>("OutputFileName");
  isData          = iConfig.getParameter<bool>("isData");
  newClustering   = iConfig.getParameter<bool>("activateNewClustering");
  oldClustering   = iConfig.getParameter<bool>("activateOldClustering");
  saveReco        = iConfig.getParameter<bool>("saveReco");
  saveUnseeded    = iConfig.getParameter<bool>("saveUnseeded");
  pathNames_      = iConfig.getParameter<std::vector<std::string>>("trgSelection");
  trigResultsTag_ = iConfig.getParameter<edm::InputTag>("trgResults");
  hitEBLabel_     = iConfig.getParameter<edm::InputTag>("hitEBLabel");
  hitEELabel_     = iConfig.getParameter<edm::InputTag>("hitEELabel");

  timeEB   = new TH1F("timeEB",   "", 400, -100, 100);
  timeEE   = new TH1F("timeEE",   "", 400, -100, 100);
  timeEB2D = new TH2F("timeEB2D", "", 400, -100, 100, 100, 0, 100);
  timeEE2D = new TH2F("timeEE2D", "", 400, -100, 100, 100, 0, 100);
}

dumper::~dumper() 
{}

int dumper::findEleRef(reco::RecoEcalCandidateRef ref, edm::Handle<reco::ElectronCollection> hltEleH) {

  int index = -1;
  for (unsigned int i=0; i<hltEleH->size(); i++) {
    reco::ElectronRef cand(hltEleH, i);
    if (cand->superCluster() == ref->superCluster()) 
      return i;
  }

  return index;
}


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
//  edm::Handle<edm::TriggerResults> trigResultsHandle;
//  event.getByLabel(trigResultsTag_,trigResultsHandle);
//
//  const edm::TriggerResults& trigResults = *trigResultsHandle;
//  const edm::TriggerNames& trigNames = event.triggerNames(trigResults);  
//
//  for(size_t pathNr=0;pathNr<pathNames_.size();pathNr++){
//    passSel[pathNr] = 0;
//    size_t pathIndex = trigNames.triggerIndex(pathNames_[pathNr]);
//    if(pathIndex<trigResults.size() &&  trigResults.accept(pathIndex)) 
//      passSel[pathNr] = 1;
//fra  }
  
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

      edm::Handle<std::vector<reco::GsfElectron> > elH;
      event.getByLabel("gsfElectrons", elH);
      reco_n = 0;
      for (unsigned int i=0; i<elH->size(); i++) {
	if (reco_n == 9)
	  continue;
	reco_et[reco_n] = 9999.;
	reco_eta[reco_n] = 9999.;
	reco_phi[reco_n] = 9999.;
	reco_pt[reco_n] = 9999.;
      
	reco::GsfElectronRef el(elH, i);
	reco_et[reco_n] = el->superCluster()->energy()*sin(el->superCluster()->position().theta());
	reco_pt[reco_n] = el->pt();
	reco_eta[reco_n] = el->superCluster()->eta();
	reco_phi[reco_n] = el->superCluster()->phi();
      reco_n++;
      }
  }


  t->Fill();
}
 
void dumper::beginJob() {
  f = new TFile(outputFileName.c_str(), "recreate");
  t = new TTree("tree", "tree");
  
  if (oldClustering) {
    t->Branch("n",  &n, "n/I");
    t->Branch("ewidth", &etawidth, "ewidth[n]/F");
    t->Branch("pwidth", &phiwidth, "pwidth[n]/F");
    t->Branch("e", &e, "e[n]/F");
    t->Branch("et", &et, "et[n]/F");
    t->Branch("se", &se, "se[n]/F");
    t->Branch("eraw", &eraw, "eraw[n]/F");
    t->Branch("eta", &eta, "eta[n]/F");
    t->Branch("phi", &phi, "phi[n]/F");
    t->Branch("seta", &seta, "seta[n]/F");
    t->Branch("sphi", &sphi, "sphi[n]/F");
    t->Branch("sieie", &sieie, "sieie[n]/F");
    t->Branch("ecal", &ecal, "ecal[n]/F");
    //t->Branch("eop", &eop, "eop[n]/F");
    t->Branch("dphi", &dphi, "dphi[n]/F");
    t->Branch("deta", &deta, "deta[n]/F");
    t->Branch("tkpt",    &tkpt,  "tkpt[n]/F");
    t->Branch("tketa",   &tketa, "tketa[n]/F");
    t->Branch("tkphi",   &tkphi, "tkphi[n]/F");
    t->Branch("hoe",   &hoe, "hoe[n]/F");
    t->Branch("hcal",   &hcal, "hcal[n]/F");
    t->Branch("tkiso",   &tkiso, "tkiso[n]/F");
  }
  
  t->Branch("passHLT", &passSel, "passHLT[3]/I");
  t->Branch("nvtx", &nvtx, "nvtx/I");
  t->Branch("rho",  &rho, "rho/F");
  
  if (newClustering) {
    t->Branch("npf",  &npf, "npf/I");
    //t->Branch("ewidthpf", &etawidthpf, "ewidthpf[npf]/F");
    //t->Branch("pwidthpf", &phiwidthpf, "pwidthpf[npf]/F");
    t->Branch("epf", &epf, "epf[npf]/F");
    t->Branch("etpf", &etpf, "etpf[npf]/F");
    t->Branch("sepf", &sepf, "sepf[npf]/F");
    //t->Branch("erawpf", &erawpf, "erawpf[npf]/F");
    t->Branch("etapf", &etapf, "etapf[npf]/F");
    t->Branch("phipf", &phipf, "phipf[npf]/F");
    //t->Branch("setapf", &setapf, "setapf[npf]/F");
    //t->Branch("sphipf", &sphipf, "sphipf[npf]/F");
    t->Branch("sieiepf", &sieiepf, "sieiepf[npf]/F");
    t->Branch("ecalpf", &ecalpf, "ecalpf[npf]/F");
    //t->Branch("eop", &eoppf, "eoppf[npf]/F");
    t->Branch("dphipf", &dphipf, "dphipf[npf]/F");
    t->Branch("detapf", &detapf, "detapf[npf]/F");
    t->Branch("tkptpf",  &tkptpf, "tkptpf[npf]/F");
    t->Branch("tketapf", &tketapf, "tketapf[npf]/F");
    t->Branch("tkphipf", &tkphipf, "tkphipf[npf]/F");
    t->Branch("hoepf",   &hoepf, "hoepf[npf]/F");
    t->Branch("hcalpf",   &hcalpf, "hcalpf[npf]/F");
    t->Branch("tkisopf",   &tkisopf, "tkisopf[npf]/F");
    t->Branch("chiso",   &chiso, "chiso[npf]/F");
    t->Branch("phiso",   &phiso, "phiso[npf]/F");
    t->Branch("neiso",   &neiso, "neiso[npf]/F");
    t->Branch("eoppf", &eoppf, "eoppf[npf]/F");
    t->Branch("chi2pf", &chi2pf, "chi2pf[npf]/F");
    t->Branch("mishitspf", &mishitspf, "mishitspf[npf]/I");
    t->Branch("hitspf", &hitspf, "hitspf[npf]/I");
  }
  
  if (saveUnseeded) {   
    t->Branch("nuns",  &nuns, "nuns/I");
    t->Branch("euns", &euns, "euns[nuns]/F");
    t->Branch("etuns", &etuns, "etuns[nuns]/F");
    t->Branch("etauns", &etauns, "etauns[nuns]/F");
    t->Branch("phiuns", &phiuns, "phiuns[nuns]/F");
  }

  if (saveReco) {
    t->Branch("recon",   &reco_n,   "recon/I");
    t->Branch("recoet",  &reco_et,  "recoet[recon]/F");
    t->Branch("recopt",  &reco_pt,  "recopt[recon]/F");
    t->Branch("recoeta", &reco_eta, "recoeta[recon]/F");
    t->Branch("recophi", &reco_phi, "recophi[recon]/F");

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
