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
    std::cout << "Usage:  ./tmp/plotHistos file.root output.root"<<std::endl
	 << std::endl;
    exit(-1);
  }


  char fileName[500];
  sprintf(fileName,argv[1]);

  char outputFileName[500];
  sprintf(outputFileName,argv[2]);
  TString outName(outputFileName);

  gROOT->ProcessLine(".x ~/rootlogon.C");

  TFile* file=TFile::Open(fileName);


  std::map<TString,TH1F*> histoNames;
  std::map<TString,TH1F*> histoNamesPU;

  std::map<TString,TH2F*> histoNames2D;
  std::map<TString,TH2F*> histoNames2DPU;

  TList* list = file->GetListOfKeys() ;
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
      histoNames[obj->GetName()]=(TH1F*)obj;
    }else if ((strcmp(obj->IsA()->GetName(),"TH2F"))==0){
      histoNames2D[obj->GetName()]=(TH2F*)obj;
    }
  }

  
  TFile *outfile=TFile::Open(outputFileName,"recreate");

  for(std::map<TString,TH1F*>::const_iterator out=histoNames.begin();out!=histoNames.end();++out){

    TPaveText * pave= new TPaveText(0.63,0.83,0.78,0.93,"NDC");
    pave->SetFillColor(kWhite);
    pave->SetTextSize(0.030);
    pave->SetTextAlign(22);
    pave->SetTextFont(62);
    TH1F* histo=out->second;
    TString name(fileName);
    if(!name.Contains("noPU"))histo->SetLineColor(kRed);
    histo->SetLineWidth(2);
    if(out->first.Contains("OverETrue")||out->first.Contains("mass")){
     
      double effectiveSigma=effSigma(histo);

      if(out->first.Contains("OverETrue")){
	histo->GetXaxis()->SetRangeUser(0.4,1.6);
      }else{
	std::cout<<effectiveSigma<<std::endl;
	effectiveSigma=effectiveSigma/histo->GetMean();
	std::cout<<effectiveSigma<<std::endl;
      }
      TString mean=Form("%.3f",histo->GetMean());
      TString sigma=Form("%.3f",effectiveSigma);
      pave->AddText("Mean:"+mean+" #sigma_{eff}:"+sigma);
      
    
      if(name.Contains("RelVal")|| name.Contains("dummy")){
	//	histo->Rebin(5);
      }
    }


      outfile->cd();
      TCanvas* c1=new TCanvas();
      c1->cd();
      histo->Draw();
      pave->Draw("same");
      c1->SaveAs("plots_splitted/"+outName+out->first+".png");

  }
}


