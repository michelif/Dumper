#include <stdlib.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <TProfile.h>
#include "createHistos.h"
#include "TKey.h"
#include "TPaveText.h"

Double_t effSigma(TH1 * hist)
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
      histoNamesnoPU[obj->GetName()]=(TH1F*)obj;
    }else if ((strcmp(obj->IsA()->GetName(),"TH2F"))==0){
      histoNames2DnoPU[obj->GetName()]=(TH2F*)obj;
    }
  }

  
  TFile *outfile=TFile::Open(outputFileName,"recreate");

  for(std::map<TString,TH1F*>::const_iterator out=histoNamesnoPU.begin();out!=histoNamesnoPU.end();++out){
    std::cout<<"plotting histo 1d:"<<out->first<<std::endl;
    TPaveText * pave= new TPaveText(0.7,0.7,0.90,0.85,"NDC");
    pave->SetFillColor(kWhite);
    pave->SetTextSize(0.030);
    pave->SetTextAlign(22);
    pave->SetTextFont(62);

    histoNamesPU[out->first]=(TH1F*)PUFile->Get(out->first);

    if(out->first.Contains("OverETrue")){
      TString meannoPU;
      meannoPU.Form("%.3f",histoNamesnoPU[out->first]->GetMean());
      double effectiveSigmanoPU=effSigma(histoNamesnoPU[out->first]);
      TString sigmanoPU=Form("%.3f",effectiveSigmanoPU);
      pave->AddText("noPU Mean:"+meannoPU+" #sigma_{eff}:"+sigmanoPU);
      
    

      TString meanPU;
      meanPU.Form("%.3f",histoNamesPU[out->first]->GetMean());
      double effectiveSigmaPU=effSigma(histoNamesPU[out->first]);
      TString sigmaPU=Form("%.3f",effectiveSigmaPU);
      pave->AddText("PU Mean:"+meanPU+" #sigma_{eff}:"+sigmaPU);
      pave->SetAllWith("PU Mean:"+meanPU+" #sigma_{eff}:"+sigmaPU,"color",kRed);
      

      TString name(noPUFileName);
      
      if(name.Contains("RelVal")|| name.Contains("dummy")){
//	histoNamesnoPU[out->first]->Rebin(4);
//	histoNamesPU[out->first]->Rebin(4);
      }
      
      histoNamesPU[out->first]->GetXaxis()->SetRangeUser(0.4,1.6);
      histoNamesnoPU[out->first]->GetXaxis()->SetRangeUser(0.4,1.6); 
    }


    histoNamesnoPU[out->first]->SetLineWidth(2);
    histoNamesnoPU[out->first]->SetLineColor(kBlack);


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

    if(out->first.Contains("OverETrue"))    pave->Draw("same");
    if(out->first.Contains("pfSC_RecHits"))c1->SetLogy();
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
