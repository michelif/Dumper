#include <stdlib.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <TProfile.h>
#include "createHistos.h"
#include "TKey.h"



int main( int argc, char* argv[] ) {


  if(argc<3 || argc>4) {
    std::cout << "Usage:  ./tmp/plotHistos noPUFile.root PUFile.root output.root"<<std::endl
	 << std::endl;
    exit(-1);
  }


  char noPUFileName[500];
  sprintf(noPUFileName,argv[1]);

  char PUFileName[500];
  sprintf(PUFileName,argv[2]);

  char outputFileName[500];
  sprintf(outputFileName,argv[3]);


  gROOT->ProcessLine(".x ~/rootlogon.C");

  TFile* noPUFile=TFile::Open(noPUFileName);
  TFile* PUFile=TFile::Open(PUFileName);


  std::map<TString,TH1F*> histoNamesnoPU;
  std::map<TString,TH1F*> histoNamesPU;

  std::map<TString,TH2F*> histoNames2DnoPU;
  std::map<TString,TH2F*> histoNames2DPU;

  TList* list = noPUFile->GetListOfKeys() ;
  if (!list) { printf("<E> No keys found in file\n") ; exit(1) ; }
  TIter next(list) ;
  TKey* key ;
  TObject* obj ;
  
  while ( key = (TKey*)next() ) {
    obj = key->ReadObj() ;
    if (    (strcmp(obj->IsA()->GetName(),"TProfile")!=0)
	    && (!obj->InheritsFrom("TH2"))
	    && (!obj->InheritsFrom("TH1")) 
	    ) {
      printf("<W> Object %s is not 1D or 2D histogram : "
	     "will not be converted\n",obj->GetName()) ;
    }
    printf("Histo name:%s title:%s\n",obj->GetName(),obj->GetTitle());
       if((strcmp(obj->IsA()->GetName(),"TH1F"))==0){
      std::cout<<"h1"<<std::endl;
      histoNamesnoPU[obj->GetName()]=(TH1F*)obj;
    }else if ((strcmp(obj->IsA()->GetName(),"TH2F"))==0){
      std::cout<<"h2"<<std::endl;
      histoNames2DnoPU[obj->GetName()]=(TH2F*)obj;
    }
  }

  
  TFile *outfile=TFile::Open(outputFileName,"recreate");

  for(std::map<TString,TH1F*>::const_iterator out=histoNamesnoPU.begin();out!=histoNamesnoPU.end();++out){
    std::cout<<"plotting histo 1d:"<<out->first<<std::endl;
    histoNamesnoPU[out->first]->SetLineWidth(2);
    histoNamesnoPU[out->first]->SetLineColor(kBlack);

    histoNamesPU[out->first]=(TH1F*)PUFile->Get(out->first);
    histoNamesPU[out->first]->SetLineWidth(2);
    histoNamesPU[out->first]->SetLineColor(kRed);


    float max = histoNamesnoPU[out->first]->GetMaximum()/histoNamesnoPU[out->first]->Integral();
    bool isnoPUhigher=true;
    if (histoNamesPU[out->first]->GetMaximum()/histoNamesPU[out->first]->Integral()>max)isnoPUhigher=false;
    //plot 1D PU vs no PU
    outfile->cd();
    TCanvas* c1=new TCanvas();
    c1->cd();
    if(isnoPUhigher){
      histoNamesnoPU[out->first]->DrawNormalized();
      histoNamesPU[out->first]->DrawNormalized("same");
    }else{
      histoNamesPU[out->first]->DrawNormalized();
      histoNamesnoPU[out->first]->DrawNormalized("same");
    }
    c1->SaveAs("plots/h1_PUvsnoPU_"+out->first+".png");
    c1->SaveAs("plots/h1_PUvsnoPU_"+out->first+".pdf");
    c1->Write(out->first);
  }


  for(std::map<TString,TH2F*>::const_iterator out=histoNames2DnoPU.begin();out!=histoNames2DnoPU.end();++out){
    std::cout<<"plotting histo "<<out->first<<std::endl;
    histoNames2DnoPU[out->first]->SetLineWidth(2);
    histoNames2DnoPU[out->first]->SetLineColor(kBlack);
    histoNames2DPU[out->first]=(TH2F*)PUFile->Get(out->first);
    histoNames2DPU[out->first]->SetLineWidth(2);
    histoNames2DPU[out->first]->SetLineColor(kRed);

    //plot 1D PU vs no PU
    TCanvas* c1=new TCanvas();
    c1->cd();
    histoNames2DnoPU[out->first]->Draw("colz");
    TProfile* prof= (TProfile*) histoNames2DnoPU[out->first]->ProfileX();
    prof->Draw("samep");
    c1->SaveAs("plots/h2_noPU_"+out->first+".png");
    c1->SaveAs("plots/h2_noPU_"+out->first+".pdf");
    c1->Write(out->first);
  
    histoNames2DPU[out->first]->Draw("colz");
    TProfile* prof2= (TProfile*) histoNames2DPU[out->first]->ProfileX();
    prof2->SetLineColor(kBlack);
    prof2->Draw("samep");
    c1->SaveAs("plots/h2_PU_"+out->first+".png");
    c1->SaveAs("plots/h2_PU_"+out->first+".pdf");
    c1->Write(out->first);
}

  outfile->Write();
  outfile->Close();

}
