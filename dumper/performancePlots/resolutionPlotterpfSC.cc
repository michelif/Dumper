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
#define NBINSET 10
#define NBINSE 10

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

   float etaMin=1.6;
   float etaMax=2.7;

   float etMin=20;
   float etMax=100;

   float eMin=30;
   float eMax=200;


   TH1F* pfSCErecoOverEtrueEta[NBINSETA];
   TH1F* pfSCErecoOverEtrueEtaUnconv[NBINSETA];
   TH1F* pfSCErecoOverEtrueEtaConv[NBINSETA];

   TH1F* pfSCErecoOverEtrueEt[NBINSET];
   TH1F* pfSCErecoOverEtrueEtUnconv[NBINSET];
   TH1F* pfSCErecoOverEtrueEtConv[NBINSET];

   TH1F* pfSCErecoOverEtrueE[NBINSE];
   TH1F* pfSCErecoOverEtrueEUnconv[NBINSE];
   TH1F* pfSCErecoOverEtrueEConv[NBINSE];

   TH1F* pfSCEseedOverEtrueEta[NBINSETA];
   TH1F* pfSCEseedOverEtrueEtaUnconv[NBINSETA];
   TH1F* pfSCEseedOverEtrueEtaConv[NBINSETA];

   TH1F* pfSCEseedOverEtrueEt[NBINSET];
   TH1F* pfSCEseedOverEtrueEtUnconv[NBINSET];
   TH1F* pfSCEseedOverEtrueEtConv[NBINSET];

   TH1F* pfSCEseedOverEtrueE[NBINSE];
   TH1F* pfSCEseedOverEtrueEUnconv[NBINSE];
   TH1F* pfSCEseedOverEtrueEConv[NBINSE];


   for(int i=0;i<NBINSETA;i++){
     TString bin;
     bin.Form("%d",i);
     pfSCErecoOverEtrueEta[i]=new TH1F("pfSCErecoOverEtrueEta"+bin,"pfSCErecoOverEtrueEta"+bin,200,0,2);
     pfSCErecoOverEtrueEtaUnconv[i]=new TH1F("pfSCErecoOverEtrueEtaUnconv"+bin,"pfSCErecoOverEtrueEtaUnconv"+bin,200,0,2);
     pfSCErecoOverEtrueEtaConv[i]=new TH1F("pfSCErecoOverEtrueEtaConv"+bin,"pfSCErecoOverEtrueEtaConv"+bin,200,0,2);

     pfSCEseedOverEtrueEta[i]=new TH1F("pfSCEseedOverEtrueEta"+bin,"pfSCEseedOverEtrueEta"+bin,200,0,2);
     pfSCEseedOverEtrueEtaUnconv[i]=new TH1F("pfSCEseedOverEtrueEtaUnconv"+bin,"pfSCEseedOverEtrueEtaUnconv"+bin,200,0,2);
     pfSCEseedOverEtrueEtaConv[i]=new TH1F("pfSCEseedOverEtrueEtaConv"+bin,"pfSCEseedOverEtrueEtaConv"+bin,200,0,2);

   }

   for(int i=0;i<NBINSET;i++){
     TString bin;
     bin.Form("%d",i);
     pfSCErecoOverEtrueEt[i]=new TH1F("pfSCErecoOverEtrueEt"+bin,"pfSCErecoOverEtrueEt"+bin,200,0,2);
     pfSCErecoOverEtrueEtUnconv[i]=new TH1F("pfSCErecoOverEtrueEtUnconv"+bin,"pfSCErecoOverEtrueEtUnconv"+bin,200,0,2);
     pfSCErecoOverEtrueEtConv[i]=new TH1F("pfSCErecoOverEtrueEtConv"+bin,"pfSCErecoOverEtrueEtConv"+bin,200,0,2);

     pfSCEseedOverEtrueEt[i]=new TH1F("pfSCEseedOverEtrueEt"+bin,"pfSCEseedOverEtrueEt"+bin,200,0,2);
     pfSCEseedOverEtrueEtUnconv[i]=new TH1F("pfSCEseedOverEtrueEtUnconv"+bin,"pfSCEseedOverEtrueEtUnconv"+bin,200,0,2);
     pfSCEseedOverEtrueEtConv[i]=new TH1F("pfSCEseedOverEtrueEtConv"+bin,"pfSCEseedOverEtrueEtConv"+bin,200,0,2);

   }

   for(int i=0;i<NBINSE;i++){
     TString bin;
     bin.Form("%d",i);
     pfSCErecoOverEtrueE[i]=new TH1F("pfSCErecoOverEtrueE"+bin,"pfSCErecoOverEtrueE"+bin,200,0,2);
     pfSCErecoOverEtrueEUnconv[i]=new TH1F("pfSCErecoOverEtrueEUnconv"+bin,"pfSCErecoOverEtrueEUnconv"+bin,200,0,2);
     pfSCErecoOverEtrueEConv[i]=new TH1F("pfSCErecoOverEtrueEConv"+bin,"pfSCErecoOverEtrueEConv"+bin,200,0,2);

     pfSCEseedOverEtrueE[i]=new TH1F("pfSCEseedOverEtrueE"+bin,"pfSCEseedOverEtrueE"+bin,200,0,2);
     pfSCEseedOverEtrueEUnconv[i]=new TH1F("pfSCEseedOverEtrueEUnconv"+bin,"pfSCEseedOverEtrueEUnconv"+bin,200,0,2);
     pfSCEseedOverEtrueEConv[i]=new TH1F("pfSCEseedOverEtrueEConv"+bin,"pfSCEseedOverEtrueEConv"+bin,200,0,2);

   }


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(jentry%500==0)std::cout<<"jentry:"<<jentry<<"/"<<nentries<<std::endl;

      for(int i=0;i<NBINSETA;i++){
	if(i!= NBINSETA-1){
	  if(TMath::Abs(pfSCeta)>etaMin+i*(etaMax-etaMin)/NBINSETA && TMath::Abs(pfSCeta)<etaMin+(i+1)*(etaMax-etaMin)/NBINSETA){
	    //	    std::cout<<pfSCeta<<" "<<etaMin+i*(etaMax-etaMin)/NBINSETA<<" "<<etaMin+(i+1)*(etaMax-etaMin)/NBINSETA<<std::endl;
	    pfSCErecoOverEtrueEta[i]->Fill(pfSCErecoOverEtrue);
	    pfSCEseedOverEtrueEta[i]->Fill(pfSCEseedOverEtrue);
	    if(pfSCisConv){
	      pfSCErecoOverEtrueEtaConv[i]->Fill(pfSCErecoOverEtrue);
	      pfSCEseedOverEtrueEtaConv[i]->Fill(pfSCEseedOverEtrue);
	    }else{
	      pfSCErecoOverEtrueEtaUnconv[i]->Fill(pfSCErecoOverEtrue);
	      pfSCEseedOverEtrueEtaUnconv[i]->Fill(pfSCEseedOverEtrue);
	    }
	  }
	}else{
	  if(TMath::Abs(pfSCeta)>etaMin+i*(etaMax-etaMin)/NBINSETA && TMath::Abs(pfSCeta)<etaMax){
	    pfSCErecoOverEtrueEta[i]->Fill(pfSCErecoOverEtrue);
	    pfSCEseedOverEtrueEta[i]->Fill(pfSCEseedOverEtrue);
	    if(pfSCisConv){
	      pfSCErecoOverEtrueEtaConv[i]->Fill(pfSCErecoOverEtrue);
	      pfSCEseedOverEtrueEtaConv[i]->Fill(pfSCEseedOverEtrue);
	    } else{
	      pfSCErecoOverEtrueEtaUnconv[i]->Fill(pfSCErecoOverEtrue);
	      pfSCEseedOverEtrueEtaUnconv[i]->Fill(pfSCEseedOverEtrue);
	    }
	  }
	}
      }



      for(int i=0;i<NBINSET;i++){
	if(i!= NBINSET-1){
	  if(TMath::Abs(pfSCpt)>etMin+i*(etMax-etMin)/NBINSET && TMath::Abs(pfSCpt)<etMin+(i+1)*(etMax-etMin)/NBINSET){
	    pfSCErecoOverEtrueEt[i]->Fill(pfSCErecoOverEtrue);
	    pfSCEseedOverEtrueEt[i]->Fill(pfSCEseedOverEtrue);
	    if(pfSCisConv){
	      pfSCErecoOverEtrueEtConv[i]->Fill(pfSCErecoOverEtrue);
	      pfSCEseedOverEtrueEtConv[i]->Fill(pfSCEseedOverEtrue);
	    }else {
	      pfSCErecoOverEtrueEtUnconv[i]->Fill(pfSCErecoOverEtrue);
	      pfSCEseedOverEtrueEtUnconv[i]->Fill(pfSCEseedOverEtrue);
	    }
	  }
	}else{
	  if(TMath::Abs(pfSCpt)>etMin+i*(etMax-etMin)/NBINSET && TMath::Abs(pfSCpt)<etMax){
	    pfSCErecoOverEtrueEt[i]->Fill(pfSCErecoOverEtrue);
	    pfSCEseedOverEtrueEt[i]->Fill(pfSCEseedOverEtrue);
	    if(pfSCisConv){
	      pfSCErecoOverEtrueEtConv[i]->Fill(pfSCErecoOverEtrue);
	      pfSCEseedOverEtrueEtConv[i]->Fill(pfSCEseedOverEtrue);
	    } else {
	      pfSCErecoOverEtrueEtUnconv[i]->Fill(pfSCErecoOverEtrue);
	      pfSCEseedOverEtrueEtUnconv[i]->Fill(pfSCEseedOverEtrue);
	    }
	  }
	}
      }

      for(int i=0;i<NBINSE;i++){
	if(i!= NBINSE-1){
	  if(TMath::Abs(pfSCe)>eMin+i*(eMax-eMin)/NBINSE && TMath::Abs(pfSCe)<eMin+(i+1)*(eMax-eMin)/NBINSE){
	    //	    std::cout<<pfSCe<<" "<<eMin+i*(eMax-eMin)/NBINSE<<" "<<eMin+(i+1)*(eMax-eMin)/NBINSE<<std::endl;
	    pfSCErecoOverEtrueE[i]->Fill(pfSCErecoOverEtrue);
	    pfSCEseedOverEtrueE[i]->Fill(pfSCEseedOverEtrue);
	    if(pfSCisConv){
	      pfSCErecoOverEtrueEConv[i]->Fill(pfSCErecoOverEtrue);	      
	      pfSCEseedOverEtrueEConv[i]->Fill(pfSCEseedOverEtrue);
	    }else{
	      pfSCErecoOverEtrueEUnconv[i]->Fill(pfSCErecoOverEtrue);
	      pfSCEseedOverEtrueEUnconv[i]->Fill(pfSCEseedOverEtrue);
	    }
	  }
	}else{
	  if(TMath::Abs(pfSCe)>eMin+i*(eMax-eMin)/NBINSE && TMath::Abs(pfSCe)<eMax){
	    pfSCErecoOverEtrueE[i]->Fill(pfSCErecoOverEtrue);
	    pfSCEseedOverEtrueE[i]->Fill(pfSCEseedOverEtrue);
	    if(pfSCisConv){
	      pfSCErecoOverEtrueEConv[i]->Fill(pfSCErecoOverEtrue);
	      pfSCEseedOverEtrueEConv[i]->Fill(pfSCEseedOverEtrue);
	    } else {
	      pfSCErecoOverEtrueEUnconv[i]->Fill(pfSCErecoOverEtrue);
	      pfSCEseedOverEtrueEUnconv[i]->Fill(pfSCEseedOverEtrue);
	    }

	  }
	}

      }




   }
   //plot vs eta
   float xEta[NBINSETA];
   float xEtaErr[NBINSETA];
   float yEta[NBINSETA];
   float yEtaErr[NBINSETA];
   float ySigmaEffEta[NBINSETA];
   float ySigmaEffEtaErr[NBINSETA];

   float ySeedEta[NBINSETA];
   float ySeedEtaErr[NBINSETA];
   float ySeedSigmaEffEta[NBINSETA];
   float ySeedSigmaEffEtaErr[NBINSETA];

   float yEtaConv[NBINSETA];
   float yEtaConvErr[NBINSETA];
   float ySigmaEffEtaConv[NBINSETA];
   float ySigmaEffEtaConvErr[NBINSETA];
   float yEtaUnconv[NBINSETA];
   float yEtaUnconvErr[NBINSETA];
   float ySigmaEffEtaUnconv[NBINSETA];
   float ySigmaEffEtaUnconvErr[NBINSETA];

   float ySeedEtaConv[NBINSETA];
   float ySeedEtaConvErr[NBINSETA];
   float ySeedSigmaEffEtaConv[NBINSETA];
   float ySeedSigmaEffEtaConvErr[NBINSETA];
   float ySeedEtaUnconv[NBINSETA];
   float ySeedEtaUnconvErr[NBINSETA];
   float ySeedSigmaEffEtaUnconv[NBINSETA];
   float ySeedSigmaEffEtaUnconvErr[NBINSETA];


   for(int i=0;i<NBINSETA;++i){
     xEta[i]=etaMin+i*(etaMax-etaMin)/NBINSETA+(etaMax-etaMin)/(2*NBINSETA);
     xEtaErr[i]=0.;
     yEta[i]=pfSCErecoOverEtrueEta[i]->GetMean();
     yEtaErr[i]=pfSCErecoOverEtrueEta[i]->GetMeanError();
     ySigmaEffEta[i]=effSigma(pfSCErecoOverEtrueEta[i]);
     ySigmaEffEtaErr[i]=pfSCErecoOverEtrueEta[i]->GetRMSError();

     yEtaConv[i]=pfSCErecoOverEtrueEtaConv[i]->GetMean();
     yEtaConvErr[i]=pfSCErecoOverEtrueEtaConv[i]->GetMeanError();
     ySigmaEffEtaConv[i]=effSigma(pfSCErecoOverEtrueEtaConv[i]);
     ySigmaEffEtaConvErr[i]=pfSCErecoOverEtrueEtaConv[i]->GetRMSError();

     yEtaUnconv[i]=pfSCErecoOverEtrueEtaUnconv[i]->GetMean();
     yEtaUnconvErr[i]=pfSCErecoOverEtrueEtaUnconv[i]->GetMeanError();
     ySigmaEffEtaUnconv[i]=effSigma(pfSCErecoOverEtrueEtaUnconv[i]);
     ySigmaEffEtaUnconvErr[i]=pfSCErecoOverEtrueEtaUnconv[i]->GetRMSError();

     ySeedEta[i]=pfSCEseedOverEtrueEta[i]->GetMean();
     ySeedEtaErr[i]=pfSCEseedOverEtrueEta[i]->GetMeanError();
     ySeedSigmaEffEta[i]=effSigma(pfSCEseedOverEtrueEta[i]);
     ySeedSigmaEffEtaErr[i]=pfSCEseedOverEtrueEta[i]->GetRMSError();

     ySeedEtaConv[i]=pfSCEseedOverEtrueEtaConv[i]->GetMean();
     ySeedEtaConvErr[i]=pfSCEseedOverEtrueEtaConv[i]->GetMeanError();
     ySeedSigmaEffEtaConv[i]=effSigma(pfSCEseedOverEtrueEtaConv[i]);
     ySeedSigmaEffEtaConvErr[i]=pfSCEseedOverEtrueEtaConv[i]->GetRMSError();

     ySeedEtaUnconv[i]=pfSCEseedOverEtrueEtaUnconv[i]->GetMean();
     ySeedEtaUnconvErr[i]=pfSCEseedOverEtrueEtaUnconv[i]->GetMeanError();
     ySeedSigmaEffEtaUnconv[i]=effSigma(pfSCEseedOverEtrueEtaUnconv[i]);
     ySeedSigmaEffEtaUnconvErr[i]=pfSCEseedOverEtrueEtaUnconv[i]->GetRMSError();
   }
   
   TGraphErrors* pfSCMeanEta = new TGraphErrors(NBINSETA,xEta,yEta,xEtaErr,yEtaErr);
   pfSCMeanEta->SetName("pfSCErecoOverEtrueEta");
   pfSCMeanEta->SetTitle("pfSCErecoOverEtrueEta");
   pfSCMeanEta->GetXaxis()->SetTitle("|#eta|");
   pfSCMeanEta->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSigmaEffEta = new TGraphErrors(NBINSETA,xEta,ySigmaEffEta,xEtaErr,ySigmaEffEtaErr);
   pfSCSigmaEffEta->SetName("pfSCSigmaEffEta");
   pfSCSigmaEffEta->SetTitle("pfSCSigmaEffEta");
   pfSCSigmaEffEta->GetXaxis()->SetTitle("|#eta|");
   pfSCSigmaEffEta->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCMeanEtaConv = new TGraphErrors(NBINSETA,xEta,yEtaConv,xEtaErr,yEtaConvErr);
   pfSCMeanEtaConv->SetName("pfSCErecoOverEtrueEtaConv");
   pfSCMeanEtaConv->SetTitle("pfSCErecoOverEtrueEtaConv");
   pfSCMeanEtaConv->GetXaxis()->SetTitle("|#eta|");
   pfSCMeanEtaConv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSigmaEffEtaConv = new TGraphErrors(NBINSETA,xEta,ySigmaEffEtaConv,xEtaErr,ySigmaEffEtaConvErr);
   pfSCSigmaEffEtaConv->SetName("pfSCSigmaEffEtaConv");
   pfSCSigmaEffEtaConv->SetTitle("pfSCSigmaEffEtaConv");
   pfSCSigmaEffEtaConv->GetXaxis()->SetTitle("|#eta|");
   pfSCSigmaEffEtaConv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCMeanEtaUnconv = new TGraphErrors(NBINSETA,xEta,yEtaUnconv,xEtaErr,yEtaUnconvErr);
   pfSCMeanEtaUnconv->SetName("pfSCErecoOverEtrueEtaUnconv");
   pfSCMeanEtaUnconv->SetTitle("pfSCErecoOverEtrueEtaUnconv");
   pfSCMeanEtaUnconv->GetXaxis()->SetTitle("|#eta|");
   pfSCMeanEtaUnconv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSigmaEffEtaUnconv = new TGraphErrors(NBINSETA,xEta,ySigmaEffEtaUnconv,xEtaErr,ySigmaEffEtaUnconvErr);
   pfSCSigmaEffEtaUnconv->SetName("pfSCSigmaEffEtaUnconv");
   pfSCSigmaEffEtaUnconv->SetTitle("pfSCSigmaEffEtaUnconv");
   pfSCSigmaEffEtaUnconv->GetXaxis()->SetTitle("|#eta|");
   pfSCSigmaEffEtaUnconv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");


   TGraphErrors* pfSCSeedMeanEta = new TGraphErrors(NBINSETA,xEta,ySeedEta,xEtaErr,ySeedEtaErr);
   pfSCSeedMeanEta->SetName("pfSCEseedOverEtrueEta");
   pfSCSeedMeanEta->SetTitle("pfSCEseedOverEtrueEta");
   pfSCSeedMeanEta->GetXaxis()->SetTitle("|#eta|");
   pfSCSeedMeanEta->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSeedSigmaEffEta = new TGraphErrors(NBINSETA,xEta,ySeedSigmaEffEta,xEtaErr,ySeedSigmaEffEtaErr);
   pfSCSeedSigmaEffEta->SetName("pfSCSeedSigmaEffEta");
   pfSCSeedSigmaEffEta->SetTitle("pfSCSeedSigmaEffEta");
   pfSCSeedSigmaEffEta->GetXaxis()->SetTitle("|#eta|");
   pfSCSeedSigmaEffEta->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCSeedMeanEtaConv = new TGraphErrors(NBINSETA,xEta,ySeedEtaConv,xEtaErr,ySeedEtaConvErr);
   pfSCSeedMeanEtaConv->SetName("pfSCEseedOverEtrueEtaConv");
   pfSCSeedMeanEtaConv->SetTitle("pfSCEseedOverEtrueEtaConv");
   pfSCSeedMeanEtaConv->GetXaxis()->SetTitle("|#eta|");
   pfSCSeedMeanEtaConv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSeedSigmaEffEtaConv = new TGraphErrors(NBINSETA,xEta,ySeedSigmaEffEtaConv,xEtaErr,ySeedSigmaEffEtaConvErr);
   pfSCSeedSigmaEffEtaConv->SetName("pfSCSeedSigmaEffEtaConv");
   pfSCSeedSigmaEffEtaConv->SetTitle("pfSCSeedSigmaEffEtaConv");
   pfSCSeedSigmaEffEtaConv->GetXaxis()->SetTitle("|#eta|");
   pfSCSeedSigmaEffEtaConv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCSeedMeanEtaUnconv = new TGraphErrors(NBINSETA,xEta,ySeedEtaUnconv,xEtaErr,ySeedEtaUnconvErr);
   pfSCSeedMeanEtaUnconv->SetName("pfSCEseedOverEtrueEtaUnconv");
   pfSCSeedMeanEtaUnconv->SetTitle("pfSCEseedOverEtrueEtaUnconv");
   pfSCSeedMeanEtaUnconv->GetXaxis()->SetTitle("|#eta|");
   pfSCSeedMeanEtaUnconv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSeedSigmaEffEtaUnconv = new TGraphErrors(NBINSETA,xEta,ySeedSigmaEffEtaUnconv,xEtaErr,ySeedSigmaEffEtaUnconvErr);
   pfSCSeedSigmaEffEtaUnconv->SetName("pfSCSeedSigmaEffEtaUnconv");
   pfSCSeedSigmaEffEtaUnconv->SetTitle("pfSCSeedSigmaEffEtaUnconv");
   pfSCSeedSigmaEffEtaUnconv->GetXaxis()->SetTitle("|#eta|");
   pfSCSeedSigmaEffEtaUnconv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");


   //plot vs et
   float xEt[NBINSET];
   float xEtErr[NBINSET];
   float yEt[NBINSET];
   float yEtErr[NBINSET];
   float ySigmaEffEt[NBINSET];

   float ySeedEt[NBINSETA];
   float ySeedEtErr[NBINSETA];
   float ySeedSigmaEffEt[NBINSETA];


   float yEtConv[NBINSET];
   float yEtConvErr[NBINSET];
   float ySigmaEffEtConv[NBINSET];
   float yEtUnconv[NBINSET];
   float yEtUnconvErr[NBINSET];
   float ySigmaEffEtUnconv[NBINSET];

   float ySeedEtConv[NBINSETA];
   float ySeedEtConvErr[NBINSETA];
   float ySeedSigmaEffEtConv[NBINSETA];
   float ySeedEtUnconv[NBINSETA];
   float ySeedEtUnconvErr[NBINSETA];
   float ySeedSigmaEffEtUnconv[NBINSETA];


   for(int i=0;i<NBINSET;++i){
     xEt[i]=etMin+i*(etMax-etMin)/NBINSET+(etMax-etMin)/(2*NBINSET);
     xEtErr[i]=0.;
     yEt[i]=pfSCErecoOverEtrueEt[i]->GetMean();
     yEtErr[i]=pfSCErecoOverEtrueEt[i]->GetMeanError();
     ySigmaEffEt[i]=effSigma(pfSCErecoOverEtrueEt[i]);

     yEtConv[i]=pfSCErecoOverEtrueEtConv[i]->GetMean();
     yEtConvErr[i]=pfSCErecoOverEtrueEtConv[i]->GetMeanError();
     ySigmaEffEtConv[i]=effSigma(pfSCErecoOverEtrueEtConv[i]);

     yEtUnconv[i]=pfSCErecoOverEtrueEtUnconv[i]->GetMean();
     yEtUnconvErr[i]=pfSCErecoOverEtrueEtUnconv[i]->GetMeanError();
     ySigmaEffEtUnconv[i]=effSigma(pfSCErecoOverEtrueEtUnconv[i]);

     ySeedEt[i]=pfSCEseedOverEtrueEt[i]->GetMean();
     ySeedEtErr[i]=pfSCEseedOverEtrueEt[i]->GetMeanError();
     ySeedSigmaEffEt[i]=effSigma(pfSCEseedOverEtrueEt[i]);

     ySeedEtConv[i]=pfSCEseedOverEtrueEtConv[i]->GetMean();
     ySeedEtConvErr[i]=pfSCEseedOverEtrueEtConv[i]->GetMeanError();
     ySeedSigmaEffEtConv[i]=effSigma(pfSCEseedOverEtrueEtConv[i]);

     ySeedEtUnconv[i]=pfSCEseedOverEtrueEtUnconv[i]->GetMean();
     ySeedEtUnconvErr[i]=pfSCEseedOverEtrueEtUnconv[i]->GetMeanError();
     ySeedSigmaEffEtUnconv[i]=effSigma(pfSCEseedOverEtrueEtUnconv[i]);


   }

   TGraphErrors* pfSCMeanEt = new TGraphErrors(NBINSET,xEt,yEt,xEtErr,yEtErr);
   pfSCMeanEt->SetName("pfSCErecoOverEtrueEt");
   pfSCMeanEt->SetTitle("pfSCErecoOverEtrueEt");
   pfSCMeanEt->GetXaxis()->SetTitle("Et");
   pfSCMeanEt->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSigmaEffEt = new TGraphErrors(NBINSET,xEt,ySigmaEffEt);
   pfSCSigmaEffEt->SetName("pfSCSeedSigmaEffEt");
   pfSCSigmaEffEt->SetTitle("pfSCSeedSigmaEffEt");
   pfSCSigmaEffEt->GetXaxis()->SetTitle("Et");
   pfSCSigmaEffEt->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCMeanEtConv = new TGraphErrors(NBINSET,xEt,yEtConv,xEtErr,yEtConvErr);
   pfSCMeanEtConv->SetName("pfSCErecoOverEtrueEtConv");
   pfSCMeanEtConv->SetTitle("pfSCErecoOverEtrueEtConv");
   pfSCMeanEtConv->GetXaxis()->SetTitle("Et");
   pfSCMeanEtConv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSigmaEffEtConv = new TGraphErrors(NBINSET,xEt,ySigmaEffEtConv);
   pfSCSigmaEffEtConv->SetName("pfSCSigmaEffEtConv");
   pfSCSigmaEffEtConv->SetTitle("pfSCSigmaEffEtConv");
   pfSCSigmaEffEtConv->GetXaxis()->SetTitle("Et");
   pfSCSigmaEffEtConv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCMeanEtUnconv = new TGraphErrors(NBINSET,xEt,yEtUnconv,xEtErr,yEtUnconvErr);
   pfSCMeanEtUnconv->SetName("pfSCErecoOverEtrueEtUnconv");
   pfSCMeanEtUnconv->SetTitle("pfSCErecoOverEtrueEtUnconv");
   pfSCMeanEtUnconv->GetXaxis()->SetTitle("Et");
   pfSCMeanEtUnconv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSigmaEffEtUnconv = new TGraphErrors(NBINSET,xEt,ySigmaEffEtUnconv);
   pfSCSigmaEffEtUnconv->SetName("pfSCSigmaEffEtUnconv");
   pfSCSigmaEffEtUnconv->SetTitle("pfSCSigmaEffEtUnconv");
   pfSCSigmaEffEtUnconv->GetXaxis()->SetTitle("Et");
   pfSCSigmaEffEtUnconv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCSeedMeanEt = new TGraphErrors(NBINSET,xEt,ySeedEt,xEtErr,ySeedEtErr);
   pfSCSeedMeanEt->SetName("pfSCEseedOverEtrueEt");
   pfSCSeedMeanEt->SetTitle("pfSCEseedOverEtrueEt");
   pfSCSeedMeanEt->GetXaxis()->SetTitle("Et");
   pfSCSeedMeanEt->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSeedSigmaEffEt = new TGraphErrors(NBINSET,xEt,ySeedSigmaEffEt);
   pfSCSeedSigmaEffEt->SetName("pfSCSeedSigmaEffEt");
   pfSCSeedSigmaEffEt->SetTitle("pfSCSeedSigmaEffEt");
   pfSCSeedSigmaEffEt->GetXaxis()->SetTitle("Et");
   pfSCSeedSigmaEffEt->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCSeedMeanEtConv = new TGraphErrors(NBINSET,xEt,ySeedEtConv,xEtErr,ySeedEtConvErr);
   pfSCSeedMeanEtConv->SetName("pfSCEseedOverEtrueEtConv");
   pfSCSeedMeanEtConv->SetTitle("pfSCEseedOverEtrueEtConv");
   pfSCSeedMeanEtConv->GetXaxis()->SetTitle("Et");
   pfSCSeedMeanEtConv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSeedSigmaEffEtConv = new TGraphErrors(NBINSET,xEt,ySeedSigmaEffEtConv);
   pfSCSeedSigmaEffEtConv->SetName("pfSCSeedSigmaEffEtConv");
   pfSCSeedSigmaEffEtConv->SetTitle("pfSCSeedSigmaEffEtConv");
   pfSCSeedSigmaEffEtConv->GetXaxis()->SetTitle("Et");
   pfSCSeedSigmaEffEtConv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCSeedMeanEtUnconv = new TGraphErrors(NBINSET,xEt,ySeedEtUnconv,xEtErr,ySeedEtUnconvErr);
   pfSCSeedMeanEtUnconv->SetName("pfSCEseedOverEtrueEtUnconv");
   pfSCSeedMeanEtUnconv->SetTitle("pfSCEseedOverEtrueEtUnconv");
   pfSCSeedMeanEtUnconv->GetXaxis()->SetTitle("Et");
   pfSCSeedMeanEtUnconv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSeedSigmaEffEtUnconv = new TGraphErrors(NBINSET,xEt,ySeedSigmaEffEtUnconv);
   pfSCSeedSigmaEffEtUnconv->SetName("pfSCSeedSigmaEffEtUnconv");
   pfSCSeedSigmaEffEtUnconv->SetTitle("pfSCSeedSigmaEffEtUnconv");
   pfSCSeedSigmaEffEtUnconv->GetXaxis()->SetTitle("Et");
   pfSCSeedSigmaEffEtUnconv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");



   //plot vs e
   float xE[NBINSE];
   float xEErr[NBINSE];
   float yE[NBINSE];
   float yEErr[NBINSE];
   float ySigmaEffE[NBINSE];

   float ySeedE[NBINSETA];
   float ySeedEErr[NBINSETA];
   float ySeedSigmaEffE[NBINSETA];


   float yEConv[NBINSE];
   float yEConvErr[NBINSE];
   float ySigmaEffEConv[NBINSE];
   float yEUnconv[NBINSE];
   float yEUnconvErr[NBINSE];
   float ySigmaEffEUnconv[NBINSE];

   float ySeedEConv[NBINSETA];
   float ySeedEConvErr[NBINSETA];
   float ySeedSigmaEffEConv[NBINSETA];
   float ySeedEUnconv[NBINSETA];
   float ySeedEUnconvErr[NBINSETA];
   float ySeedSigmaEffEUnconv[NBINSETA];



   for(int i=0;i<NBINSE;++i){
     xE[i]=eMin+i*(eMax-eMin)/NBINSE+(eMax-eMin)/(2*NBINSE);
     xEErr[i]=0.;
     yE[i]=pfSCErecoOverEtrueE[i]->GetMean();
     yEErr[i]=pfSCErecoOverEtrueE[i]->GetMeanError();
     ySigmaEffE[i]=effSigma(pfSCErecoOverEtrueE[i]);

     yEConv[i]=pfSCErecoOverEtrueEConv[i]->GetMean();
     yEConvErr[i]=pfSCErecoOverEtrueEConv[i]->GetMeanError();
     ySigmaEffEConv[i]=effSigma(pfSCErecoOverEtrueEConv[i]);

     yEUnconv[i]=pfSCErecoOverEtrueEUnconv[i]->GetMean();
     yEUnconvErr[i]=pfSCErecoOverEtrueEUnconv[i]->GetMeanError();
     ySigmaEffEUnconv[i]=effSigma(pfSCErecoOverEtrueEUnconv[i]);

     ySeedE[i]=pfSCEseedOverEtrueE[i]->GetMean();
     ySeedEErr[i]=pfSCEseedOverEtrueE[i]->GetMeanError();
     ySeedSigmaEffE[i]=effSigma(pfSCEseedOverEtrueE[i]);

     ySeedEConv[i]=pfSCEseedOverEtrueEConv[i]->GetMean();
     ySeedEConvErr[i]=pfSCEseedOverEtrueEConv[i]->GetMeanError();
     ySeedSigmaEffEConv[i]=effSigma(pfSCEseedOverEtrueEConv[i]);

     ySeedEUnconv[i]=pfSCEseedOverEtrueEUnconv[i]->GetMean();
     ySeedEUnconvErr[i]=pfSCEseedOverEtrueEUnconv[i]->GetMeanError();
     ySeedSigmaEffEUnconv[i]=effSigma(pfSCEseedOverEtrueEUnconv[i]);



   }

   TGraphErrors* pfSCMeanE = new TGraphErrors(NBINSE,xE,yE,xEErr,yEErr);
   pfSCMeanE->SetName("pfSCErecoOverEtrueE");
   pfSCMeanE->SetTitle("pfSCErecoOverEtrueE");
   pfSCMeanE->GetXaxis()->SetTitle("E");
   pfSCMeanE->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSigmaEffE = new TGraphErrors(NBINSE,xE,ySigmaEffE);
   pfSCSigmaEffE->SetName("pfSCSigmaEffE");
   pfSCSigmaEffE->SetTitle("pfSCSigmaEffE");
   pfSCSigmaEffE->GetXaxis()->SetTitle("E");
   pfSCSigmaEffE->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCMeanEConv = new TGraphErrors(NBINSE,xE,yEConv,xEErr,yEConvErr);
   pfSCMeanEConv->SetName("pfSCErecoOverEtrueEConv");
   pfSCMeanEConv->SetTitle("pfSCErecoOverEtrueEConv");
   pfSCMeanEConv->GetXaxis()->SetTitle("E");
   pfSCMeanEConv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSigmaEffEConv = new TGraphErrors(NBINSE,xE,ySigmaEffEConv);
   pfSCSigmaEffEConv->SetName("pfSCSigmaEffEConv");
   pfSCSigmaEffEConv->SetTitle("pfSCSigmaEffEConv");
   pfSCSigmaEffEConv->GetXaxis()->SetTitle("E");
   pfSCSigmaEffEConv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCMeanEUnconv = new TGraphErrors(NBINSE,xE,yEUnconv,xEErr,yEUnconvErr);
   pfSCMeanEUnconv->SetName("pfSCErecoOverEtrueEUnconv");
   pfSCMeanEUnconv->SetTitle("pfSCErecoOverEtrueEUnconv");
   pfSCMeanEUnconv->GetXaxis()->SetTitle("E");
   pfSCMeanEtUnconv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSigmaEffEUnconv = new TGraphErrors(NBINSE,xE,ySigmaEffEUnconv);
   pfSCSigmaEffEUnconv->SetName("pfSCSigmaEffEUnconv");
   pfSCSigmaEffEUnconv->SetTitle("pfSCSigmaEffEUnconv");
   pfSCSigmaEffEUnconv->GetXaxis()->SetTitle("E");
   pfSCSigmaEffEUnconv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");


   TGraphErrors* pfSCSeedMeanE = new TGraphErrors(NBINSE,xE,ySeedE,xEErr,ySeedEErr);
   pfSCSeedMeanE->SetName("pfSCEseedOverEtrueE");
   pfSCSeedMeanE->SetTitle("pfSCEseedOverEtrueE");
   pfSCSeedMeanE->GetXaxis()->SetTitle("Et");
   pfSCSeedMeanE->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSeedSigmaEffE = new TGraphErrors(NBINSE,xE,ySeedSigmaEffE);
   pfSCSeedSigmaEffE->SetName("pfSCSeedSigmaEffE");
   pfSCSeedSigmaEffE->SetTitle("pfSCSeedSigmaEffE");
   pfSCSeedSigmaEffE->GetXaxis()->SetTitle("Et");
   pfSCSeedSigmaEffE->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCSeedMeanEConv = new TGraphErrors(NBINSE,xE,ySeedEConv,xEErr,ySeedEConvErr);
   pfSCSeedMeanEConv->SetName("pfSCEseedOverEtrueEConv");
   pfSCSeedMeanEConv->SetTitle("pfSCEseedOverEtrueEConv");
   pfSCSeedMeanEConv->GetXaxis()->SetTitle("Et");
   pfSCSeedMeanEConv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSeedSigmaEffEConv = new TGraphErrors(NBINSE,xE,ySeedSigmaEffEConv);
   pfSCSeedSigmaEffEConv->SetName("pfSCSeedSigmaEffEConv");
   pfSCSeedSigmaEffEConv->SetTitle("pfSCSeedSigmaEffEConv");
   pfSCSeedSigmaEffEConv->GetXaxis()->SetTitle("Et");
   pfSCSeedSigmaEffEConv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");

   TGraphErrors* pfSCSeedMeanEUnconv = new TGraphErrors(NBINSE,xE,ySeedEUnconv,xEErr,ySeedEUnconvErr);
   pfSCSeedMeanEUnconv->SetName("pfSCEseedOverEtrueEUnconv");
   pfSCSeedMeanEUnconv->SetTitle("pfSCEseedOverEtrueEUnconv");
   pfSCSeedMeanEUnconv->GetXaxis()->SetTitle("Et");
   pfSCSeedMeanEUnconv->GetYaxis()->SetTitle("<E_{reco}/E_{true}>");

   TGraphErrors* pfSCSeedSigmaEffEUnconv = new TGraphErrors(NBINSE,xE,ySeedSigmaEffEUnconv);
   pfSCSeedSigmaEffEUnconv->SetName("pfSCSeedSigmaEffEUnconv");
   pfSCSeedSigmaEffEUnconv->SetTitle("pfSCSeedSigmaEffEUnconv");
   pfSCSeedSigmaEffEUnconv->GetXaxis()->SetTitle("Et");
   pfSCSeedSigmaEffEUnconv->GetYaxis()->SetTitle("#sigma_{eff}^{SC}");



   outFile_->cd();
   pfSCMeanEta->Write();
   pfSCSigmaEffEta->Write();
   pfSCMeanEtaConv->Write();
   pfSCSigmaEffEtaConv->Write();
   pfSCMeanEtaUnconv->Write();
   pfSCSigmaEffEtaUnconv->Write();

   pfSCSeedMeanEta->Write();
   pfSCSeedSigmaEffEta->Write();
   pfSCSeedMeanEtaConv->Write();
   pfSCSeedSigmaEffEtaConv->Write();
   pfSCSeedMeanEtaUnconv->Write();
   pfSCSeedSigmaEffEtaUnconv->Write();


   pfSCMeanEt->Write();
   pfSCSigmaEffEt->Write();
   pfSCMeanEtConv->Write();
   pfSCSigmaEffEtConv->Write();
   pfSCMeanEtUnconv->Write();
   pfSCSigmaEffEtUnconv->Write();

   pfSCMeanE->Write();
   pfSCSigmaEffE->Write();
   pfSCMeanEConv->Write();
   pfSCSigmaEffEConv->Write();
   pfSCMeanEUnconv->Write();
   pfSCSigmaEffEUnconv->Write();



   pfSCMeanE->Write();
   pfSCSigmaEffE->Write();

   pfSCSeedMeanE->Write();
   pfSCSeedSigmaEffE->Write();
   pfSCSeedMeanEConv->Write();
   pfSCSeedSigmaEffEConv->Write();
   pfSCSeedMeanEUnconv->Write();
   pfSCSeedSigmaEffEUnconv->Write();

   pfSCSeedMeanEt->Write();
   pfSCSeedSigmaEffEt->Write();
   pfSCSeedMeanEtConv->Write();
   pfSCSeedSigmaEffEtConv->Write();
   pfSCSeedMeanEtUnconv->Write();
   pfSCSeedSigmaEffEtUnconv->Write();



   outFile_->Write();
   outFile_->Close();


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

