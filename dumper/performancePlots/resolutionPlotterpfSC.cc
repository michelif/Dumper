#define resolutionPlotterpfSC_cxx
#include <stdlib.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include "resolutionPlotterpfSC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#define NBINSETA 20
#define NBINSET 20


void resolutionPlotterpfSC::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L resolutionPlotterpfSC.C
//      Root > resolutionPlotterpfSC t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   TFile * outFile=TFile::Open("resPlotterPhotons.root","recreate");

   float etaMin=1.6;
   float etaMax=2.7;
   TH1F* pfSCErecoOverEtrueEta[NBINSETA];
   for(int i=0;i<NBINSETA;i++){
     TString bin;
     bin.Form("%d",i);
     pfSCErecoOverEtrueEta[i]=new TH1F("pfSCErecoOverEtrueEta"+bin,"pfSCErecoOverEtrueEta"+bin,200,0,2);
   }


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(jentry%500==0)std::cout<<"jentry:"<<jentry<<"/"<<nentries<<std::endl;

      for(int i=0;i<NBINSETA;i++){
	if(i!= NBINSETA-1){
	  if(TMath::Abs(pfSCeta)>etaMin+i*(etaMax-etaMin)/NBINSETA && TMath::Abs(pfSCeta)<etaMin+(i+1)*(etaMax-etaMin)){
	    pfSCErecoOverEtrueEta[i]->Fill(pfSCErecoOverEtrue);
	  }
	}else{
	  if(TMath::Abs(pfSCeta)>etaMin+i*(etaMax-etaMin)/NBINSETA && TMath::Abs(pfSCeta)<2.7){
	    pfSCErecoOverEtrueEta[i]->Fill(pfSCErecoOverEtrue);
	  }
	}
      }
   }

   float xEta[NBINSETA];
   float xEtaErr[NBINSETA];
   float yEta[NBINSETA];
   float yEtaErr[NBINSETA];
   float ySigmaEff[NBINSETA];
   for(int i=0;i<NBINSETA;++i){
     xEta[i]=etaMin+i*(etaMax-etaMin)/NBINSETA+(etaMax-etaMin)/(2*NBINSETA);
     xEtaErr[i]=0.;
     yEta[i]=pfSCErecoOverEtrueEta[i]->GetMean();
     yEtaErr[i]=pfSCErecoOverEtrueEta[i]->GetMeanError();
     ySigmaEff[i]=effSigma(pfSCErecoOverEtrueEta[i]);
   }

   TGraphErrors* pfSCMeanEta = new TGraphErrors(NBINSET,xEta,yEta,xEtaErr,yEtaErr);
   pfSCMeanEta->SetName("pfSCErecoOverEtrueEta");
   pfSCMeanEta->SetTitle("pfSCErecoOverEtrueEta");
   pfSCMeanEta->GetXaxis()->SetTitle("|#eta|");
   pfSCMeanEta->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSigmaEffEta = new TGraphErrors(NBINSET,xEta,ySigmaEff);
   pfSCSigmaEffEta->SetName("pfSCSigmaEffEta");
   pfSCSigmaEffEta->SetTitle("pfSCSigmaEffEta");
   pfSCSigmaEffEta->GetXaxis()->SetTitle("|#eta|");
   pfSCSigmaEffEta->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");


   outFile->cd();
   pfSCMeanEta->Write();
   pfSCSigmaEffEta->Write();
   outFile->Write();
   outFile->Close();


}


Double_t resolutionPlotterpfSC::effSigma(TH1 * hist)
{
  using namespace std;
  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return 0.;
  }
 
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  //   if(total < 100.) {
  //     cout << "effsigma: Too few entries " << total << endl;
  //     return 0.;
  //   }
  Int_t ierr=0;
  Int_t ismin=999;
 
  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }  
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
 
  return widmin;
 
}

