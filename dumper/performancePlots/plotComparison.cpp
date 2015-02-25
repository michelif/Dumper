#include <stdlib.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include "TLegend.h"
#include "createHistos.h"

using namespace std;



int main( int argc, char* argv[] ) {



  if(argc<3 || argc>8) {
    cout << "Usage:  ./tmp/runCreateHistos listfile1  listfile2  listfile3  listfile4 outputfile\n"
	 << "  listfile:    list of root files incusing protocol eg dcap:/// .....\n"
	 << "  outputfile:  name of output root file  eg output.root\n"
	 << endl;
    exit(-1);
  }

  int nLists=argc-2;  
  std::cout<<"nLists:"<<nLists<<std::endl;
  //  Input list 1
  char listName1[500];
  char listName2[500];
  char listName3[500];
  char listName4[500];
  TString nameString[4];

  sprintf(listName1,argv[1]); 
  std::cout<<listName1<<std::endl;

  // Input list 2
  if (nLists >= 2){
  sprintf(listName2,argv[2]);
  }

  // Input list 3
  if (nLists >= 3){
  sprintf(listName3,argv[3]);
  }

  // Input list 4
  if (nLists >= 4){
  sprintf(listName4,argv[4]);
  }

  // Output filename (.root)  
  TString OutputFileName(argv[nLists+1]);
  
  // Name of input tree objects in (.root) files 
  char treeName[100] = "tree";

  TChain* chain[4];

  //creating  TChain
  for (int i=0;i<nLists;i++){
    chain[i] = new TChain(treeName);
    char pName[500];
    char listName[500];
    if(i==0){
      for (int i =0;i<500;++i){
	listName[i]=listName1[i];
      }
    }
    if(i==1){
      for (int i =0;i<500;++i){
	listName[i]=listName2[i];
      }
    }
    if(i==2){
      for (int i =0;i<500;++i){
	listName[i]=listName3[i];
      }
    }
    if(i==3){
      for (int i =0;i<500;++i){
	listName[i]=listName4[i];
      }
    }
    ifstream is(listName);

    if(! is.good()) {
      cout << "int main() >> ERROR : file " << listName << " not read" << endl;
      is.close();
      exit(-1);
    }
    cout << "Reading list : " << listName << " ......." << endl;
  
    while( is.getline(pName, 500, '\n') ) {
      if (pName[0] == '#') continue;
      std::cout<<pName<<std::endl;
      chain[i]->Add(pName); 
    }
    is.close();
  }

    TCanvas* c1=new TCanvas();
    c1->cd();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gPad->SetLogy();
    TLegend* legend=new TLegend(0.7,0.7,0.95,0.95);
    TH1F* histos[4];
    TH1F* histostime[4];
    for (int i=0;i<nLists;++i){
    std::cout<<"processing file:"<<i<<std::endl;
    TString istring;
    istring.Form("%d",i);
    chain[i]->SetLineWidth(2);
    histos[i]=new TH1F("histo"+istring,"histo"+istring,100,0.,20.);
    histostime[i]=new TH1F("histotime"+istring,"histotime"+istring,100,-50.5,49.5);

    chain[i]->Project("histo"+istring,"rechite","rechite>0");
    chain[i]->Project("histotime"+istring,"rechittime","rechite>0");
    histos[i]->Print();
    if(i==0){
      histos[i]->SetLineColor(kBlack);
      histostime[i]->SetLineColor(kBlack);
      legend->AddEntry(histos[i],"Single Mu noPU","l");
    }
    if(i==1){
      histos[i]->SetLineColor(kRed);
      histostime[i]->SetLineColor(kRed);
      legend->AddEntry(histos[i],"Single Mu PU140","l");
    }
    if(i==2){
      histostime[i]->SetLineColor(kBlue);
      legend->AddEntry(histos[i],"DY noPU","l");
    }
    if(i==3){
      histos[i]->SetLineColor(kViolet);
      histostime[i]->SetLineColor(kViolet);
      legend->AddEntry(histos[i],"DY PU140","l");
    }


  }

    legend->SetFillColor(kWhite);

    TFile *outfile=TFile::Open(OutputFileName,"recreate");


    for (int i=0;i<nLists;++i){
      if(i==0){
	histos[i]->DrawNormalized();
      } else{
	histos[i]->DrawNormalized("same");
      }
      histos[i]->Write();
    }
    
    
    legend->Draw("same");    
    c1->SaveAs("comparisonE.png");
    c1->Write();

    TCanvas* c2=new TCanvas();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gPad->SetLogy();
    c2->cd();

    for (int i=0;i<nLists;++i){
      if(i==0){
	histostime[i]->DrawNormalized();
      } else{
	histostime[i]->DrawNormalized("same");
      }
      histostime[i]->Write();
    }
    legend->Draw("same");        

    c2->SaveAs("comparisonTime.png");

    legend->Write();
    outfile->Write();    
    outfile->Close();
    return 0;

}
