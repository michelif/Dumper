#include <iostream>
#include <string>
#include <stdlib.h>
#include "createHistos.h"


void createHistos::bookHistos(){

  //electrons
  bookHisto("ele_ErecoOverETrue",200,0,2,"E_{reco}/E_{true}",false);
  bookHisto("ele_sMaj",200,0,3,"sMaj",false);
  bookHisto("ele_sMin",200,0,2,"sMin",false);
  bookHisto("ele_alpha",200,-2,2,"alpha",false);

  //sc variables
  std::cout<<"booking"<<std::endl;
  bookHisto("pfSC_ErecoOverETrue",200,0,2,"E_{reco}/E_{true}");
  bookHisto("pfSC_EnergyWrong",150,0,300,"E_{reco}");
  bookHisto("pfSC_EnergyRight",150,0,300,"E_{reco}");
  bookHisto("pfSC_nXtalsTotalWrong",250,-0.5,249.5,"N_{xtals}");
  bookHisto("pfSC_nXtalsTotalRight",250,-0.5,249.5,"N_{xtals}");
  bookHisto("pfSC_nXtalsSeedWrong",250,-0.5,249.5,"N_{xtals}");
  bookHisto("pfSC_nXtalsSeedRight",250,-0.5,249.5,"N_{xtals}");
  bookHisto("pfSC_nXtalsSeed",100,-0.5,99.5,"N_{xtals}^{seed}");
  bookHisto("pfSC_nXtalsTotal",250,-0.5,249.5,"N_{xtals}");
  bookHisto("pfSC_nXtalsBC",100,-0.5,99.5,"N_{xtals}^{BC}");
  bookHisto("pfSC_nBCForSC",25,-0.5,24.5,   "N_{BC} for SC");
  bookHisto("pfSC_maxDistFromSeedinRinSCEle",300,0,0.8  ,"max#Delta R_{BC}^{seed}");
  bookHisto("pfSC_maxDistFromSeedinEtainSCEle",300,0,0.8  ,"max#Delta R_{BC}^{seed}");
  bookHisto("pfSC_maxDistFromSeedinPhiinSCEle",300,0,0.8  ,"max#Delta R_{BC}^{seed}");
  bookHisto("pfSC_bcNXtals",100,-0.5,99.5,"N_{xtals}^{BC}");
  bookHisto2D("pfSC_EBCseedVsDeltaPhiBCSeedEle",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
  bookHisto2D("pfSC_EBCseedVsDeltaEtaBCSeedEle",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
  bookHisto2D("pfSC_EnergyVsnXtalsTotalWrong",250,-0.5,249.5,150,0,300,"N_{xtals}","E_{reco}");
  bookHisto2D("pfSC_EnergyVsnXtalsTotalRight",250,-0.5,249.5,150,0,300,"N_{xtals}","E_{reco}");
  bookHisto2D("pfSC_EnergyVsnXtalsSeedWrong",250,-0.5,249.5,150,0,300,"N_{xtals}^{seed}","E_{reco}");
  bookHisto2D("pfSC_EnergyVsnXtalsSeedRight",250,-0.5,249.5,150,0,300,"N_{xtals}^{seed}","E_{reco}");
  bookHisto2D("pfSC_ErecoMinusEtrueVsEffectiveArea",250,0,400,150,-1.,1,"#rhoxN_{xtals}","#frac{E_{reco}-E_{true}}{E_{reco}}");


//  bookHisto("pfSC_FirstEtaBin_maxDistFromSeedinRinSCPhoUnconv",300,0,0.8,"max#Delta R_{BC}^{seed}");
//  bookHisto("pfSC_FirstEtaBin_maxDistFromSeedinEtainSCPhoUnconv",300,0,0.8,"max#Delta #eta_{BC}^{seed}");
//  bookHisto("pfSC_FirstEtaBin_maxDistFromSeedinPhiinSCPhoUnconv",300,0,0.8,"max#Delta #phi_{BC}^{seed}");
//  bookHisto("pfSC_FirstEtaBin_maxDistFromSeedinRinSCPhoConv",300,0,0.8,"max#Delta R_{BC}^{seed}");
//  bookHisto("pfSC_FirstEtaBin_maxDistFromSeedinEtainSCPhoConv",300,0,0.8,"max#Delta #eta_{BC}^{seed}");
//  bookHisto("pfSC_FirstEtaBin_maxDistFromSeedinPhiinSCPhoConv",300,0,0.8,"max#Delta #phi_{BC}^{seed}");
//

//  bookHisto("pfSC_SecondEtaBin_maxDistFromSeedinRinSCPhoUnconv",300,0,0.8,"max#Delta R_{BC}^{seed}");
//  bookHisto("pfSC_SecondEtaBin_maxDistFromSeedinEtainSCPhoUnconv",300,0,0.8,"max#Delta #eta_{BC}^{seed}");
//  bookHisto("pfSC_SecondEtaBin_maxDistFromSeedinPhiinSCPhoUnconv",300,0,0.8,"max#Delta #phi_{BC}^{seed}");
//  bookHisto("pfSC_SecondEtaBin_maxDistFromSeedinRinSCPhoConv",300,0,0.8,"max#Delta R_{BC}^{seed}");
//  bookHisto("pfSC_SecondEtaBin_maxDistFromSeedinEtainSCPhoConv",300,0,0.8,"max#Delta #eta_{BC}^{seed}");
//  bookHisto("pfSC_SecondEtaBin_maxDistFromSeedinPhiinSCPhoConv",300,0,0.8,"max#Delta #phi_{BC}^{seed}");
//
//  bookHisto("pfSC_ThirdEtaBin_maxDistFromSeedinRinSCPhoUnconv",300,0,0.8,"max#Delta R_{BC}^{seed}");
//  bookHisto("pfSC_ThirdEtaBin_maxDistFromSeedinEtainSCPhoUnconv",300,0,0.8,"max#Delta #eta_{BC}^{seed}");
//  bookHisto("pfSC_ThirdEtaBin_maxDistFromSeedinPhiinSCPhoUnconv",300,0,0.8,"max#Delta #phi_{BC}^{seed}");
//  bookHisto("pfSC_ThirdEtaBin_maxDistFromSeedinRinSCPhoConv",300,0,0.8,"max#Delta R_{BC}^{seed}");
//  bookHisto("pfSC_ThirdEtaBin_maxDistFromSeedinEtainSCPhoConv",300,0,0.8,"max#Delta #eta_{BC}^{seed}");
//  bookHisto("pfSC_ThirdEtaBin_maxDistFromSeedinPhiinSCPhoConv",300,0,0.8,"max#Delta #phi_{BC}^{seed}");
//
//  //multi5x5 bc variables
//  bookHisto("multi5x5SC_ErecoOverETrueFirstEtaBin",200,0,2,"E_{reco}/E_{true}");
//  bookHisto("multi5x5SC_ErecoOverETrueSecondEtaBin",200,0,2,"E_{reco}/E_{true}");
//  bookHisto("multi5x5SC_ErecoOverETrueThirdEtaBin",200,0,2,"E_{reco}/E_{true}");
//
//
//
//
//
//  bookHisto2D("sieieVsPhi",50,-3.,3.,50,0.,0.1,"#phi","#sigma_{i#etai#eta}");
//  bookHisto2D("sieieVsPhiFirstEtaBin",50,-3.,3.,50,0.,0.1,"#phi","#sigma_{i#etai#eta}");
//  bookHisto2D("sieieVsPhiSecondEtaBin",50,-3.,3.,50,0.,0.1,"#phi","#sigma_{i#etai#eta}");
//  bookHisto2D("sieieVsPhiThirdEtaBin",50,-3.,3.,50,0.,0.1,"#phi","#sigma_{i#etai#eta}");
//

//  bookHisto2D("pfSC_FirstEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_FirstEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_FirstEtaBin_EBCseedVsDeltaEtaBCSeedEle",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_FirstEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_FirstEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//
//  bookHisto2D("pfSC_SecondEtaBin_EBCseedVsDeltaPhiBCSeedEle",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_SecondEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_SecondEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_SecondEtaBin_EBCseedVsDeltaEtaBCSeedEle",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_SecondEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_SecondEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//

//  bookHisto2D("pfSC_ThirdEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_ThirdEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_ThirdEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//  bookHisto2D("pfSC_ThirdEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//
//

//  bookHisto2D("multi5x5SC_FirstEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("multi5x5SC_FirstEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("multi5x5SC_FirstEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//  bookHisto2D("multi5x5SC_FirstEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//
//  bookHisto2D("multi5x5SC_SecondEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("multi5x5SC_SecondEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("multi5x5SC_SecondEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//  bookHisto2D("multi5x5SC_SecondEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//
//  bookHisto2D("multi5x5SC_ThirdEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("multi5x5SC_ThirdEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv",150,0.,0.7,150,0.,50.,"#Delta#phi_{BC}^{seed}","E_{BC}");
//  bookHisto2D("multi5x5SC_ThirdEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");
//  bookHisto2D("multi5x5SC_ThirdEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv",150,0.,0.3,150,0.,50.,"#Delta#eta_{BC}^{seed}","E_{BC}");


}


void createHistos::bookHisto(TString name, int nbins, float xLow, float xUp){
  histos_[name]=new TH1F(name, name, nbins,xLow,xUp);

}

void createHistos::bookHisto(TString name, int nbins, float xLow, float xUp,TString xAxisName, bool bookMulti5x5){

    TString histonameIncl=name+"_inclusive";
    histos_[histonameIncl]=new TH1F(histonameIncl, histonameIncl, nbins,xLow,xUp);    
    setAxisTitle(histonameIncl,xAxisName);
    if(bookMulti5x5){
      histonameIncl.Replace(0,2,"multi5x5");
      histos_[histonameIncl]=new TH1F(histonameIncl, histonameIncl, nbins,xLow,xUp);    
      setAxisTitle(histonameIncl,xAxisName);
    }
  
  for (std::map<TString,std::vector<float> >::const_iterator itCatCuts=categoriesAndCuts_.begin();itCatCuts!=categoriesAndCuts_.end();++itCatCuts){


    for(int i=0;i<itCatCuts->second.size();++i){
      TString bin;
      bin.Form("%d",i);
      TString histoname=name+"_"+itCatCuts->first+"bin"+bin;
      histos_[histoname]=new TH1F(histoname, histoname, nbins,xLow,xUp);    
      setAxisTitle(histoname,xAxisName);
      if(bookMulti5x5){
	histoname.Replace(0,2,"multi5x5");
	histos_[histoname]=new TH1F(histoname, histoname, nbins,xLow,xUp);    
	setAxisTitle(histoname,xAxisName);
      }
    }
    
  }
  
}

//farla double float int
void createHistos::fillHisto2D(TString name, double valueX, double valueY, TLorentzVector* p4){
  for (std::map<TString,std::vector<float> >::const_iterator itCatCuts=categoriesAndCuts_.begin();itCatCuts!=categoriesAndCuts_.end();++itCatCuts){

    TString histonameIncl=name+"_inclusive";
    if(p4->Eta() > itCatCuts->second[0]) histos2D_[histonameIncl]->Fill(valueX, valueY);

    if(itCatCuts->first.EqualTo("eta")){
      for(int i=0;i<itCatCuts->second.size();++i){
	TString bin;
	bin.Form("%d",i);
	TString histoname=name+"_"+itCatCuts->first+"bin"+bin;
	if(i<itCatCuts->second.size()-1){
	  if(p4->Eta() > itCatCuts->second[i] && p4->Eta() < itCatCuts->second[i+1] ){
	    histos2D_[histoname]->Fill(valueX, valueY);
	  }
	}else {
	  if(p4->Eta() > itCatCuts->second[i])
	    histos2D_[histoname]->Fill(valueX,valueY);
	}
      }
    }
  }
  
}


//farla double float int
void createHistos::fillHisto(TString name, double value,TLorentzVector* p4){
  for (std::map<TString,std::vector<float> >::const_iterator itCatCuts=categoriesAndCuts_.begin();itCatCuts!=categoriesAndCuts_.end();++itCatCuts){
    TString histonameIncl=name+"_inclusive";
    if(p4->Eta() > itCatCuts->second[0])    histos_[histonameIncl]->Fill(value);
    

    if(itCatCuts->first.EqualTo("eta")){
      for(int i=0;i<itCatCuts->second.size();++i){
	TString bin;
	bin.Form("%d",i);
	TString histoname=name+"_"+itCatCuts->first+"bin"+bin;
	if(i<itCatCuts->second.size()-1){
	  if(p4->Eta() > itCatCuts->second[i] && p4->Eta() < itCatCuts->second[i+1] ){
	    histos_[histoname]->Fill(value);
	  }
	}else {
	  if(p4->Eta() > itCatCuts->second[i])
	    histos_[histoname]->Fill(value);
	}
      }
    }
  }
  
}

void createHistos::bookHisto2D(TString name, int nbins, float xLow, float xUp,int nbinsY, float yLow, float yUp){
  histos2D_[name]=new TH2F(name, name, nbins,xLow,xUp,nbins, yLow,yUp);

}

void createHistos::bookHisto2D(TString name, int nbins, float xLow, float xUp,int nbinsY, float yLow, float yUp,TString xAxisName, TString yAxisName, bool bookMulti5x5){

    TString histonameIncl=name+"_inclusive";
    histos2D_[histonameIncl]=new TH2F(histonameIncl, histonameIncl, nbins,xLow,xUp,nbins, yLow,yUp);
    setAxisTitle(histonameIncl,xAxisName,yAxisName);
     if(bookMulti5x5){
       histonameIncl.Replace(0,2,"multi5x5");
       histos2D_[histonameIncl]=new TH2F(histonameIncl, histonameIncl, nbins,xLow,xUp,nbins, yLow,yUp);
       setAxisTitle(histonameIncl,xAxisName,yAxisName);
     }

  for (std::map<TString,std::vector<float> >::const_iterator itCatCuts=categoriesAndCuts_.begin();itCatCuts!=categoriesAndCuts_.end();++itCatCuts){
    for(int i=0;i<itCatCuts->second.size();++i){
      TString bin;
      bin.Form("%d",i);
      TString histoname=name+"_"+itCatCuts->first+"bin"+bin;
      histos2D_[histoname]=new TH2F(histoname, histoname, nbins,xLow,xUp,nbins, yLow,yUp);
      setAxisTitle(histoname,xAxisName,yAxisName);
      if(bookMulti5x5){
	histoname.Replace(0,2,"multi5x5");
	histos2D_[histoname]=new TH2F(histoname, histoname, nbins,xLow,xUp,nbins, yLow,yUp);
	setAxisTitle(histoname,xAxisName,yAxisName);
      }
    }
  }
}


void createHistos::setAxisTitle(TString name,TString xAxisName){
  histos_[name]->GetXaxis()->SetTitle(xAxisName);
}

void createHistos::setAxisTitle(TString name,TString xAxisName, TString yAxisName){
  histos2D_[name]->GetXaxis()->SetTitle(xAxisName);
  histos2D_[name]->GetYaxis()->SetTitle(yAxisName);
}


void createHistos::writeHistos(){
  for(std::map<TString,TH1F*>::const_iterator out=histos_.begin();out!=histos_.end();++out){
    out->second->Write(out->first,TObject::kOverwrite);
  }

  for(std::map<TString,TH2F*>::const_iterator out=histos2D_.begin();out!=histos2D_.end();++out){
    out->second->Write(out->first,TObject::kOverwrite);
  }

 
}


void createHistos::defineCategories(){

  categoriesAndCuts_["eta"].push_back(1.5);
  categoriesAndCuts_["eta"].push_back(1.8);
  categoriesAndCuts_["eta"].push_back(2.4);

  //  categoriesAndCuts_["fbrem"].push_back(0.2);
}

TLorentzVector* createHistos::createTLorentzVector(float pt, float eta, float phi, float e){
  TLorentzVector* theVector=new TLorentzVector();
  theVector->SetPtEtaPhiE(pt,eta,phi,e);
  return theVector;

}

void createHistos::buildGenEle(){

  for(int j=0;j<gelen;j++){
    if(gpstatusMC[geleindex[j]]!=3)continue;
    if(gelept[j]<0.5)continue;
    
    TLorentzVector* gelep4 = createTLorentzVector(gelept[j],geleeta[j],gelephi[j],gelept[j]*cosh(geleeta[j]));
    theGenElectrons_.push_back(gelep4);	
    //    std::cout<<"build:"<<gelep4->E()<<" "<<gelep4->Eta()<<" "<<gelep4->Phi()<<std::endl;  
  }
}

void createHistos::buildGenPho(){

  for(int j=0;j<gphon;j++){
    if(gpstatusMC[gphoindex[j]]!=3)continue;
    if(gphopt[j]<5)continue;

    TLorentzVector* gphop4 = createTLorentzVector(gphopt[j],gphoeta[j],gphophi[j],gphopt[j]*cosh(gphoeta[j]));
    theGenPhotons_.push_back(gphop4);	
  }
}

int createHistos::matchesGenEle(TLorentzVector* objectToMatch, float DeltaR ){
  bool matches=false;
  int indexMatch=-1;  
  if(objectToMatch->Pt()<0.0000001)return indexMatch;
  
  float drMatch=999;

  for(int i=0;i<theGenElectrons_.size();i++){
    if(theGenElectrons_[i]->Pt()<0.0000001)continue;
    float dr = objectToMatch->DeltaR(*theGenElectrons_[i]);
    if(dr<DeltaR && dr<drMatch){
      drMatch=dr;
      indexMatch=i;
    }
  }

  return indexMatch;
}


int createHistos::matchesGenPho(TLorentzVector* objectToMatch, float DeltaR ){
  bool matches=false;
  int indexMatch=-1;  
  if(objectToMatch->Pt()<0.0000001)return indexMatch;
  
  float drMatch=999;

  for(int i=0;i<theGenPhotons_.size();i++){
    if(theGenPhotons_[i]->Pt()<0.0000001)continue;
    float dr = objectToMatch->DeltaR(*theGenPhotons_[i]);
    if(dr<DeltaR && dr<drMatch){
      drMatch=dr;
      indexMatch=i;
    }
  }
  return indexMatch;
}




void createHistos::Loop2(){

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntries();
  std::cout<<"nentries:"<<nentries<<std::endl;
  Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
  //  for (Long64_t jentry=0; jentry<20;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;           
    
    if(jentry%500==0)std::cout<<"jentry:"<<jentry<<"/"<<nentries<<std::endl;
 
    buildGenEle();
    buildGenPho();

    //loop on ele
    for(int i=0;i<elen;i++){
      if(elept[i]<0.1)continue;

      TLorentzVector* elep4=createTLorentzVector(elept[i],eleeta[i],elephi[i],elee[i]);
      int indexMatchEle=matchesGenEle(elep4);
      if(indexMatchEle<0)continue;
      fillHisto("ele_ErecoOverETrue",elee[i]/theGenElectrons_[indexMatchEle]->E(),elep4);
      fillHisto("ele_sMaj",elesMajZS[i],elep4);
      fillHisto("ele_sMin",elesMinZS[i],elep4);
      fillHisto("ele_alpha",elealphaZS[i],elep4);
   }


    //loop on pfSC
    if(pfSCn>0){
      for (int i=0;i<pfSCn;i++){
	
	TLorentzVector* pfscp4 = createTLorentzVector(pfSCe[i]/cosh(pfSCeta[i]),pfSCeta[i],pfSCphi[i],pfSCe[i]);
	
	int indexMatchEle=matchesGenEle(pfscp4);
	if(indexMatchEle<0)continue;
	//	if(indexMatchEle>-1)    std::cout<<"gen:"<<theGenElectrons_[indexMatchEle]->E()<<" "<<theGenElectrons_[indexMatchEle]->Eta()<<" "<<theGenElectrons_[indexMatchEle]->Phi()<<std::endl;  
	//	std::cout<<"pfsc:"<<pfscp4->E()<<" "<<pfscp4->Eta()<<" "<<pfscp4->Phi()<<std::endl;

	fillHisto("pfSC_ErecoOverETrue",pfSCe[i]/theGenElectrons_[indexMatchEle]->E(),pfscp4);
	if(pfSCe[i]/theGenElectrons_[indexMatchEle]->E()<0.5){
	  fillHisto("pfSC_EnergyWrong",pfSCe[i],pfscp4);
	  fillHisto("pfSC_nXtalsSeedWrong",pfSCnXtalsSeed[i],pfscp4);
	  fillHisto("pfSC_nXtalsTotalWrong",pfSCnXtalsTotal[i],pfscp4);
	  fillHisto2D("pfSC_EnergyVsnXtalsSeedWrong",pfSCnXtalsSeed[i],pfSCe[i],pfscp4);
	  fillHisto2D("pfSC_EnergyVsnXtalsTotalWrong",pfSCnXtalsTotal[i],pfSCe[i],pfscp4);
	}else{
	  fillHisto("pfSC_EnergyRight",pfSCe[i],pfscp4);
	  fillHisto("pfSC_nXtalsSeedWrong",pfSCnXtalsSeed[i],pfscp4);
	  fillHisto("pfSC_nXtalsTotalRight",pfSCnXtalsTotal[i],pfscp4);
	  fillHisto2D("pfSC_EnergyVsnXtalsSeedRight",pfSCnXtalsSeed[i],pfSCe[i],pfscp4);
	  fillHisto2D("pfSC_EnergyVsnXtalsTotalRight",pfSCnXtalsTotal[i],pfSCe[i],pfscp4);
	} 

	fillHisto("pfSC_nXtalsSeed",pfSCnXtalsSeed[i],pfscp4);
	fillHisto("pfSC_nXtalsTotal",pfSCnXtalsTotal[i],pfscp4);
	fillHisto("pfSC_nBCForSC",pfSCnBC[i],pfscp4);
	fillHisto2D("pfSC_ErecoMinusEtrueVsEffectiveArea",rho*pfSCnXtalsTotal[i],(pfSCe[i]-theGenElectrons_[indexMatchEle]->E())/pfSCe[i],pfscp4);

	TLorentzVector seed;
	seed.SetPtEtaPhiE(pfSCbcE[i][0]/cosh(pfSCbcEta[i][0]),pfSCbcEta[i][0],pfSCbcPhi[i][0],pfSCbcE[i][0]);

	float maxDistR=0;
	float maxDistEta=0;
	float maxDistPhi=0;

	for (int j=0;j< pfSCnBC[i];j++){
	  //	  std::cout<<pfSCnBC[i]<< " i,j:"<<i<<","<<j<<" "<<pfSCbcE[i][j]<<std::endl;
	  if(pfSCbcE[i][j]<0.01)continue;
	  TLorentzVector* bc = createTLorentzVector(pfSCbcE[i][j]/cosh(pfSCbcEta[i][j]),pfSCbcEta[i][j],pfSCbcPhi[i][j],pfSCbcE[i][j]);
	  fillHisto("pfSC_bcNXtals",pfSCbcNXtals[i][j],bc);

	  if(j>0){
	    float distR=bc->DeltaR(seed);
	    if(distR>maxDistR)maxDistR=distR;
	    float distPhi=bc->DeltaPhi(seed);
	    if(distPhi>maxDistPhi)maxDistPhi=distPhi;
	    float distEta=sqrt(distR*distR-distPhi*distPhi);
	    if(distEta>maxDistEta)maxDistEta=distEta;
	    fillHisto2D("pfSC_EBCseedVsDeltaPhiBCSeedEle",fabs(distPhi),pfSCbcE[i][j],bc);
	    fillHisto2D("pfSC_EBCseedVsDeltaEtaBCSeedEle",fabs(distEta),pfSCbcE[i][j],bc);
	  }
	}//pfscnBC
	if(maxDistR>0)fillHisto("pfSC_maxDistFromSeedinRinSCEle", maxDistR, pfscp4);
	if(maxDistEta>0)fillHisto("pfSC_maxDistFromSeedinEtainSCEle", maxDistEta, pfscp4);
	if(maxDistPhi>0)fillHisto("pfSC_maxDistFromSeedinPhiinSCEle", maxDistPhi, pfscp4);
      }//pfscn
    }//if pfscn>0
    

    //loop on multi5x5SC
    if(multi5x5SCn>0){
      for (int i=0;i<multi5x5SCn;i++){
	
	TLorentzVector* multi5x5scp4 = createTLorentzVector(multi5x5SCe[i]/cosh(multi5x5SCeta[i]),multi5x5SCeta[i],multi5x5SCphi[i],multi5x5SCe[i]);
	
	int indexMatchEle=matchesGenEle(multi5x5scp4);
	if(indexMatchEle<0)continue;

	fillHisto("multi5x5SC_ErecoOverETrue",multi5x5SCe[i]/theGenElectrons_[indexMatchEle]->E(),multi5x5scp4);

	if(multi5x5SCe[i]/theGenElectrons_[indexMatchEle]->E()<0.5){
	  fillHisto("multi5x5SC_EnergyWrong",multi5x5SCe[i],multi5x5scp4);
	  fillHisto("multi5x5SC_nXtalsSeedWrong",multi5x5SCnXtalsSeed[i],multi5x5scp4);
	  fillHisto("multi5x5SC_nXtalsTotalWrong",multi5x5SCnXtalsTotal[i],multi5x5scp4);
	  fillHisto2D("multi5x5SC_EnergyVsnXtalsSeedWrong",multi5x5SCnXtalsSeed[i],multi5x5SCe[i],multi5x5scp4);
	  fillHisto2D("multi5x5SC_EnergyVsnXtalsTotalWrong",multi5x5SCnXtalsTotal[i],multi5x5SCe[i],multi5x5scp4);
	}else{
	  fillHisto("multi5x5SC_EnergyRight",multi5x5SCe[i],multi5x5scp4);
	  fillHisto("multi5x5SC_nXtalsSeedRight",multi5x5SCnXtalsSeed[i],multi5x5scp4);
	  fillHisto("multi5x5SC_nXtalsTotalRight",multi5x5SCnXtalsTotal[i],multi5x5scp4);
	  fillHisto2D("multi5x5SC_EnergyVsnXtalsSeedRight",multi5x5SCnXtalsSeed[i],multi5x5SCe[i],multi5x5scp4);
	  fillHisto2D("multi5x5SC_EnergyVsnXtalsTotalRight",multi5x5SCnXtalsTotal[i],multi5x5SCe[i],multi5x5scp4);
	}

	fillHisto("multi5x5SC_nXtalsSeed",multi5x5SCnXtalsSeed[i],multi5x5scp4);
	fillHisto("multi5x5SC_nXtalsTotal",multi5x5SCnXtalsTotal[i],multi5x5scp4);
	fillHisto("multi5x5SC_nBCForSC",multi5x5SCnBC[i],multi5x5scp4);
	fillHisto2D("multi5x5SC_ErecoMinusEtrueVsEffectiveArea",rho*multi5x5SCnXtalsTotal[i],(multi5x5SCe[i]-theGenElectrons_[indexMatchEle]->E())/multi5x5SCe[i],multi5x5scp4);

	TLorentzVector seed;
	seed.SetPtEtaPhiE(multi5x5SCbcE[i][0]/cosh(multi5x5SCbcEta[i][0]),multi5x5SCbcEta[i][0],multi5x5SCbcPhi[i][0],multi5x5SCbcE[i][0]);

	float maxDistR=0;
	float maxDistEta=0;
	float maxDistPhi=0;

	for (int j=0;j< multi5x5SCnBC[i];j++){
	  //	  std::cout<<multi5x5SCnBC[i]<< " i,j:"<<i<<","<<j<<" "<<multi5x5SCbcE[i][j]<<std::endl;
	  if(multi5x5SCbcE[i][j]<0.01)continue;
	  

	  if(j>0){
	    TLorentzVector* bc = createTLorentzVector(multi5x5SCbcE[i][j]/cosh(multi5x5SCbcEta[i][j]),multi5x5SCbcEta[i][j],multi5x5SCbcPhi[i][j],multi5x5SCbcE[i][j]);
	    float distR=bc->DeltaR(seed);
	    if(distR>maxDistR)maxDistR=distR;
	    float distPhi=bc->DeltaPhi(seed);
	    if(distPhi>maxDistPhi)maxDistPhi=distPhi;
	    float distEta=sqrt(distR*distR-distPhi*distPhi);
	    if(distEta>maxDistEta)maxDistEta=distEta;
	    fillHisto2D("multi5x5SC_EBCseedVsDeltaPhiBCSeedEle",fabs(distPhi),multi5x5SCbcE[i][j],bc);
	    fillHisto2D("multi5x5SC_EBCseedVsDeltaEtaBCSeedEle",fabs(distEta),multi5x5SCbcE[i][j],bc);
	  }
	}//multi5x5scnBC
	if(maxDistR>0)fillHisto("multi5x5SC_maxDistFromSeedinRinSCEle", maxDistR, multi5x5scp4);
	if(maxDistEta>0)fillHisto("multi5x5SC_maxDistFromSeedinEtainSCEle", maxDistEta, multi5x5scp4);
	if(maxDistPhi>0)fillHisto("multi5x5SC_maxDistFromSeedinPhiinSCEle", maxDistPhi, multi5x5scp4);
      }//multi5x5scn
    }//if multi5x5scn>0


    }
  
}

void createHistos::Loop(){

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  std::cout<<"nentries:"<<nentries<<std::endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;           

    if(jentry%500==0)std::cout<<"jentry:"<<jentry<<"/"<<nentries<<std::endl;

    float endcapBoundaryLow=1.5;
    float firstEtaBinUp=1.8;
    float secondEtaBinDown=1.8;
    float secondEtaBinUp=2.4;
    float thirdEtaBinDown=2.4;

    //loop on electrons
    for(int i=0;i<elen;i++){
      if(elept[i]<0.1)continue;

      //matching with gen ele
      TLorentzVector elep4;
      elep4.SetPtEtaPhiM(elept[i],eleeta[i],elephi[i],0.);
      //      std::cout<<"etaele:"<<elep4.Eta();                                                                                                                                                                                           

      int geleindexMatch=-1;
      float drMatch=999;
      for(int j=0;j<gelen;j++){
	if(gelept[j]<5)continue;

	TLorentzVector gelep4;
	gelep4.SetPtEtaPhiM(gelept[j],geleeta[j],gelephi[j],0.);

	float dr=elep4.DeltaR(gelep4);

	if(dr<0.1 && dr<drMatch){
	  geleindexMatch=j;
	  drMatch=dr;
	}
      }
      //filling histos for electrons
      if(geleindexMatch!=-1 && TMath::Abs(eleeta[i])>endcapBoundaryLow){
	histos2D_["sieieVsPhi"]->Fill(elephi[i],elesiEtaiEtaZS[i]);


	//	std::cout<<elephi[i]<<" "<<elesiEtaiEtaZS[i]<<std::endl;
	if(TMath::Abs(eleeta[i])<firstEtaBinUp && TMath::Abs(eleeta[i])>endcapBoundaryLow){
	  histos_["eleErecoOverETrueFirstEtaBin"]->Fill(elee[i]/(gelept[geleindexMatch]*cosh(geleeta[geleindexMatch])));
	  if(gelefbrem80[geleindexMatch]<0.2)histos_["eleErecoOverETrueFirstEtaBinFbrem02"]->Fill(elee[i]/(gelept[geleindexMatch]*cosh(geleeta[geleindexMatch])));
	  histos2D_["sieieVsPhiFirstEtaBin"]->Fill(elephi[i],elesiEtaiEtaZS[i]);
	}
	else if(TMath::Abs(eleeta[i])>secondEtaBinDown && TMath::Abs(eleeta[i])<secondEtaBinUp){
	  histos_["eleErecoOverETrueSecondEtaBin"]->Fill(elee[i]/(gelept[geleindexMatch]*cosh(geleeta[geleindexMatch])));
	  if(gelefbrem80[geleindexMatch]<0.2)histos_["eleErecoOverETrueSecondEtaBinFbrem02"]->Fill(elee[i]/(gelept[geleindexMatch]*cosh(geleeta[geleindexMatch])));
	  histos2D_["sieieVsPhiSecondEtaBin"]->Fill(elephi[i],elesiEtaiEtaZS[i]);
	}
	else if(TMath::Abs(eleeta[i])>thirdEtaBinDown){
	  histos_["eleErecoOverETrueThirdEtaBin"]->Fill(elee[i]/(gelept[geleindexMatch]*cosh(geleeta[geleindexMatch])));
	  if(gelefbrem80[geleindexMatch]<0.2)histos_["eleErecoOverETrueThirdEtaBinFbrem02"]->Fill(elee[i]/(gelept[geleindexMatch]*cosh(geleeta[geleindexMatch])));
	  histos2D_["sieieVsPhiThirdEtaBin"]->Fill(elephi[i],elesiEtaiEtaZS[i]);
	}
      }
    }//elen


    //loop on photons
    for(int i=0;i<phon;i++){
      if(phopt[i]<0.1)continue;

      //matching with gen pho
      TLorentzVector phop4;
      phop4.SetPtEtaPhiM(phopt[i],phoeta[i],phophi[i],0.);
      //      std::cout<<"etapho:"<<phop4.Eta();                                                                                                                                                                                           

      int gphoindexMatch=-1;
      float drMatch=999;
      for(int j=0;j<gphon;j++){
	if(gphopt[j]<5)continue;

	TLorentzVector gphop4;
	gphop4.SetPtEtaPhiM(gphopt[j],gphoeta[j],gphophi[j],0.);

	float dr=phop4.DeltaR(gphop4);

	if(dr<0.1 && dr<drMatch){
	  gphoindexMatch=j;
	  drMatch=dr;
	}
      }
      //filling histos for photons
      if(gphoindexMatch!=-1){
	if(TMath::Abs(phoeta[i])<firstEtaBinUp && TMath::Abs(phoeta[i])>endcapBoundaryLow){
	  histos_["phoErecoOverETrueFirstEtaBin"]->Fill(phoe[i]/(gphopt[gphoindexMatch]*cosh(gphoeta[gphoindexMatch])));
	}
	else if(TMath::Abs(phoeta[i])>secondEtaBinDown && TMath::Abs(phoeta[i])<secondEtaBinUp){
	  histos_["phoErecoOverETrueSecondEtaBin"]->Fill(phoe[i]/(gphopt[gphoindexMatch]*cosh(gphoeta[gphoindexMatch])));
	}
	else if(TMath::Abs(phoeta[i])>thirdEtaBinDown){
	  histos_["phoErecoOverETrueThirdEtaBin"]->Fill(phoe[i]/(gphopt[gphoindexMatch]*cosh(gphoeta[gphoindexMatch])));
	}
      }
    }//phon



    //loop on pfSC
    if(pfSCn>0){
      for (int i=0;i<pfSCn;i++){
	

	//matching with gen ele
	bool matchesGenEle=false;
	bool matchesGenPhoConv=false;
	bool matchesGenPhoUnconv=false;

	bool matchesRecoEle=false;
	int  recoEleindexMatch=-1;
	float drRecoMatch=999;
	TLorentzVector pfscp4;
	pfscp4.SetPtEtaPhiE(pfSCe[i]/cosh(pfSCeta[i]),pfSCeta[i],pfSCphi[i],pfSCe[i]);

	for(int ii=0;ii<elen;ii++){
	  if(elept[ii]<0.1)continue;
	  
	  //matching with reco ele
	  TLorentzVector elep4;
	  elep4.SetPtEtaPhiM(elept[ii],eleeta[ii],elephi[ii],0.5);
	  //      std::cout<<"etaele:"<<elep4.Eta(); 
	  float drReco=elep4.DeltaR(pfscp4);
	  
	  if(drReco<drRecoMatch && drReco<0.1){
	    recoEleindexMatch=ii;
	    drRecoMatch=drReco;
	    matchesRecoEle=true;
	  }
	}
	


	if(!matchesRecoEle)continue;

	TLorentzVector recoEleMatchp4;
	if(recoEleindexMatch!=-1) recoEleMatchp4.SetPtEtaPhiM(elept[recoEleindexMatch],eleeta[recoEleindexMatch],elephi[recoEleindexMatch],0.5);
	
	int geleindexMatch=-1;
	float drMatch=999;
	float geleMatchedE=-1;
	for(int j=0;j<gelen;j++){
	  if(gelept[j]<5)continue;
	  TLorentzVector gelep4;
	  gelep4.SetPtEtaPhiM(gelept[j],geleeta[j],gelephi[j],0.5);
	  
	  float dr=gelep4.DeltaR(recoEleMatchp4);
	  
	  if(dr<0.1 && dr<drMatch){
	    geleindexMatch=j;
	    drMatch=dr;
	    geleMatchedE=gelep4.E();
	  }
	}

	if(fabs(geleeta[geleindexMatch])<endcapBoundaryLow)continue;

	if(geleindexMatch != -1){
	  matchesGenEle=true;
	  //	  std::cout<<gelept[geleindexMatch]<< " " <<elept[recoEleindexMatch]<<" "<<pfSCe[i]/cosh(pfSCeta[i])<<std::endl;

	  //ereco

	  //	  std::cout<<pfSCe[i]/(gelept[geleindexMatch]*cosh(geleeta[geleindexMatch]))<<std::endl;

	  if(TMath::Abs(pfSCeta[i])<firstEtaBinUp){
	    histos_["pfSC_ErecoOverETrueFirstEtaBin"]->Fill(pfSCe[i]/geleMatchedE);
	    histos_["pfSC_FirstEtaBin_nXtalsSeed"]->Fill(pfSCnXtalsSeed[i]);
	    histos_["pfSC_FirstEtaBin_nBCForSC"]->Fill(pfSCnBC[i]);
	  }else if(TMath::Abs(pfSCeta[i])>secondEtaBinDown && TMath::Abs(pfSCeta[i])<secondEtaBinUp){
	    histos_["pfSC_ErecoOverETrueSecondEtaBin"]->Fill(pfSCe[i]/geleMatchedE);
	    histos_["pfSC_SecondEtaBin_nXtalsSeed"]->Fill(pfSCnXtalsSeed[i]);
	    histos_["pfSC_SecondEtaBin_nBCForSC"]->Fill(pfSCnBC[i]);
	  }else if(TMath::Abs(pfSCeta[i])>thirdEtaBinDown){
	    histos_["pfSC_ErecoOverETrueThirdEtaBin"]->Fill(pfSCe[i]/geleMatchedE);
	    histos_["pfSC_ThirdEtaBin_nXtalsSeed"]->Fill(pfSCnXtalsSeed[i]);
	    histos_["pfSC_ThirdEtaBin_nBCForSC"]->Fill(pfSCnBC[i]);
	  }      


//	  std::cout<<"pfsc: pt,eta,phi "<<pfSCe[i]/cosh(pfSCeta[i])<<","<<pfSCeta[i]<<","<<pfSCphi[i]<<","<<pfSCe[i]<<std::endl;
//	  std::cout<<"gele: pt,eta,phi "<<gelept[geleindexMatch]<<","<<geleeta[geleindexMatch]<<","<<gelephi[geleindexMatch]<<std::endl;
	}

	TLorentzVector matchedGenPhop4;
	int gphoindexMatch=-1;
	float drMatchPho=999;
	for(int j=0;j<gphon;j++){
	  if(gphopt[j]<5)continue;
	  
	  TLorentzVector gphop4;
	  gphop4.SetPtEtaPhiM(gphopt[j],gphoeta[j],gphophi[j],0.);
	  
	  float dr=pfscp4.DeltaR(gphop4);
	  
	  if(dr<0.1 && dr<drMatchPho){
	    gphoindexMatch=j;
	    drMatchPho=dr;
	    matchedGenPhop4=gphop4;
	  }
	}

	if(gphoindexMatch != -1) {
	  //FIX ME: for photons for the moment you don't have the info on conversion so you have to reloop on photons, try to add when producing ntuples
	  int recoPhoindexMatch=-1;
	  float drPhoMatch=999;

	    for(int ii=0;ii<phon;ii++){
	      if(phopt[ii]<0.1)continue;
	      TLorentzVector phop4;
	      phop4.SetPtEtaPhiM(phopt[ii],phoeta[ii],phophi[ii],0.);

	      float dr=matchedGenPhop4.DeltaR(phop4);

	      if(dr<0.1 && dr<drPhoMatch){
		recoPhoindexMatch=ii;
		drPhoMatch=dr;
	      }
	    }

	    if(recoPhoindexMatch>-1){
//	      std::cout<<"pfsc: pt,eta,phi "<<pfSCe[i]/cosh(pfSCeta[i])<<","<<pfSCeta[i]<<","<<pfSCphi[i]<<","<<pfSCe[i]<<std::endl;
//	      std::cout<<"genPho: pt,eta,phi "<<gphopt[gphoindexMatch]<<","<<gphoeta[gphoindexMatch]<<","<<gphophi[gphoindexMatch]<<std::endl;
//	      std::cout<<"recoPho: pt,eta,phi "<<phoe[recoPhoindexMatch]/cosh(phoeta[recoPhoindexMatch])<<","<<phoeta[recoPhoindexMatch]<<","<<phophi[recoPhoindexMatch]<<","<<phoe[recoPhoindexMatch]<<std::endl;
	      if(phoE9[recoPhoindexMatch]/phoe[recoPhoindexMatch]>0.94) matchesGenPhoUnconv=true;
	      else matchesGenPhoConv=true;
	    }
	}


	if(!(matchesGenEle || matchesGenPhoUnconv || matchesGenPhoConv)) continue;

	//max distance from seed
	if(pfSCbcE[i][0]<0.01)continue;
	//	std::cout<<pfSCbcE[i][0]<<std::endl;
	float maxDistR=0;
	float maxDistEta=0;
	float maxDistPhi=0;


	TLorentzVector seed;
	seed.SetPtEtaPhiE(pfSCbcE[i][0]/cosh(pfSCbcEta[i][0]),pfSCbcEta[i][0],pfSCbcPhi[i][0],pfSCbcE[i][0]);



	for (int j=0;j< pfSCnBC[i];j++){
	  //	  std::cout<<pfSCnBC[i]<< " i,j:"<<i<<","<<j<<" "<<pfSCbcE[i][j]<<std::endl;
	  if(pfSCbcE[i][j]<0.01)continue;
	  if(j>0){
	    TLorentzVector bc;
	    bc.SetPtEtaPhiE(pfSCbcE[i][j]/cosh(pfSCbcEta[i][j]),pfSCbcEta[i][j],pfSCbcPhi[i][j],pfSCbcE[i][j]);
	    float distR=seed.DeltaR(bc);
	    if(distR>maxDistR)maxDistR=distR;
	    float distPhi=seed.DeltaPhi(bc);
	    if(distPhi>maxDistPhi)maxDistPhi=distPhi;
	    float distEta=sqrt(distR*distR-distPhi*distPhi);
	    if(distEta>maxDistEta)maxDistEta=distEta;
	    
	    //2d: ebcseed vs deltaphi(bc,seed)
	    if(matchesGenEle){
	      //	      std::cout<<fabs(distPhi)<<","<<pfSCbcE[i][j]<<std::endl;
	      if(TMath::Abs(pfSCbcEta[i][j])<firstEtaBinUp){
		histos2D_["pfSC_FirstEtaBin_EBCseedVsDeltaPhiBCSeedEle"]->Fill(fabs(distPhi),pfSCbcE[i][j]);
		histos2D_["pfSC_FirstEtaBin_EBCseedVsDeltaEtaBCSeedEle"]->Fill(fabs(distEta),pfSCbcE[i][j]);
	      }else if(TMath::Abs(pfSCbcEta[i][j])>secondEtaBinDown && TMath::Abs(pfSCbcEta[i][j])<secondEtaBinUp){
		histos2D_["pfSC_SecondEtaBin_EBCseedVsDeltaPhiBCSeedEle"]->Fill(fabs(distPhi),pfSCbcE[i][j]);
		histos2D_["pfSC_SecondEtaBin_EBCseedVsDeltaEtaBCSeedEle"]->Fill(fabs(distEta),pfSCbcE[i][j]);
	      }else if(TMath::Abs(pfSCbcEta[i][j])>thirdEtaBinDown){
		histos2D_["pfSC_ThirdEtaBin_EBCseedVsDeltaPhiBCSeedEle"]->Fill(fabs(distPhi),pfSCbcE[i][j]);
		histos2D_["pfSC_ThirdEtaBin_EBCseedVsDeltaEtaBCSeedEle"]->Fill(fabs(distEta),pfSCbcE[i][j]);	
	      }
	    }
	    else if(matchesGenPhoUnconv){
	      if(TMath::Abs(pfSCbcEta[i][j])<firstEtaBinUp){
		histos2D_["pfSC_FirstEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv"]->Fill(fabs(distPhi),pfSCbcE[i][j]);
		histos2D_["pfSC_FirstEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv"]->Fill(fabs(distEta),pfSCbcE[i][j]);
	      }else if(TMath::Abs(pfSCbcEta[i][j])>secondEtaBinDown && TMath::Abs(pfSCbcEta[i][j])<secondEtaBinUp){
		histos2D_["pfSC_SecondEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv"]->Fill(fabs(distPhi),pfSCbcE[i][j]);
		histos2D_["pfSC_SecondEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv"]->Fill(fabs(distEta),pfSCbcE[i][j]);
	      }else if(TMath::Abs(pfSCbcEta[i][j])>thirdEtaBinDown){
		histos2D_["pfSC_ThirdEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv"]->Fill(fabs(distPhi),pfSCbcE[i][j]);
		histos2D_["pfSC_ThirdEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv"]->Fill(fabs(distEta),pfSCbcE[i][j]);
	      }
	    }else if(matchesGenPhoConv){
	      if(TMath::Abs(pfSCbcEta[i][j])<firstEtaBinUp){
		histos2D_["pfSC_FirstEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv"]->Fill(fabs(distPhi),pfSCbcE[i][j]);
		histos2D_["pfSC_FirstEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv"]->Fill(fabs(distEta),pfSCbcE[i][j]);
	      }else if(TMath::Abs(pfSCbcEta[i][j])>secondEtaBinDown && TMath::Abs(pfSCbcEta[i][j])<secondEtaBinUp){
		histos2D_["pfSC_SecondEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv"]->Fill(fabs(distPhi),pfSCbcE[i][j]);
		histos2D_["pfSC_SecondEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv"]->Fill(fabs(distEta),pfSCbcE[i][j]);
	      }else if(TMath::Abs(pfSCbcEta[i][j])>thirdEtaBinDown){
		histos2D_["pfSC_ThirdEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv"]->Fill(fabs(distPhi),pfSCbcE[i][j]);
		histos2D_["pfSC_ThirdEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv"]->Fill(fabs(distEta),pfSCbcE[i][j]);
	      }
	    }

	    //	    std::cout<<maxDistR<<std::endl;
	  }
	}//pfSCnBC
	if(matchesGenEle){
	  if(TMath::Abs(pfSCeta[i])<firstEtaBinUp){
	    if(maxDistR>0)histos_[  "pfSC_FirstEtaBin_maxDistFromSeedinRinSCEle"]->Fill(maxDistR);
	    if(maxDistEta>0)histos_["pfSC_FirstEtaBin_maxDistFromSeedinEtainSCEle"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["pfSC_FirstEtaBin_maxDistFromSeedinPhiinSCEle"]->Fill(maxDistPhi);	
	  } else if(TMath::Abs(pfSCeta[i])>secondEtaBinDown && TMath::Abs(pfSCeta[i])<secondEtaBinUp){
	    if(maxDistR>0)histos_[  "pfSC_SecondEtaBin_maxDistFromSeedinRinSCEle"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["pfSC_SecondEtaBin_maxDistFromSeedinEtainSCEle"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["pfSC_SecondEtaBin_maxDistFromSeedinPhiinSCEle"]->Fill(maxDistPhi);
	  }else if(TMath::Abs(pfSCeta[i])>thirdEtaBinDown){
	    if(maxDistR>0)histos_[  "pfSC_ThirdEtaBin_maxDistFromSeedinRinSCEle"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["pfSC_ThirdEtaBin_maxDistFromSeedinEtainSCEle"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["pfSC_ThirdEtaBin_maxDistFromSeedinPhiinSCEle"]->Fill(maxDistPhi);
	  }
	}else if (matchesGenPhoUnconv){
	  if(TMath::Abs(pfSCeta[i])<firstEtaBinUp){
	    if(maxDistR>0)histos_[  "pfSC_FirstEtaBin_maxDistFromSeedinRinSCPhoUnconv"]->Fill(maxDistR);
	    if(maxDistEta>0)histos_["pfSC_FirstEtaBin_maxDistFromSeedinEtainSCPhoUnconv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["pfSC_FirstEtaBin_maxDistFromSeedinPhiinSCPhoUnconv"]->Fill(maxDistPhi);
	  } else if(TMath::Abs(pfSCeta[i])>secondEtaBinDown && TMath::Abs(pfSCeta[i])<secondEtaBinUp){
	    if(maxDistR>0)histos_[  "pfSC_SecondEtaBin_maxDistFromSeedinRinSCPhoUnconv"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["pfSC_SecondEtaBin_maxDistFromSeedinEtainSCPhoUnconv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["pfSC_SecondEtaBin_maxDistFromSeedinPhiinSCPhoUnconv"]->Fill(maxDistPhi);
	  }else if(TMath::Abs(pfSCeta[i])>thirdEtaBinDown){
	    if(maxDistR>0)histos_[  "pfSC_ThirdEtaBin_maxDistFromSeedinRinSCPhoUnconv"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["pfSC_ThirdEtaBin_maxDistFromSeedinEtainSCPhoUnconv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["pfSC_ThirdEtaBin_maxDistFromSeedinPhiinSCPhoUnconv"]->Fill(maxDistPhi);
	  }

	}else if (matchesGenPhoConv){
	  if(TMath::Abs(pfSCeta[i])<firstEtaBinUp){
	    if(maxDistR>0)histos_[  "pfSC_FirstEtaBin_maxDistFromSeedinRinSCPhoConv"]->Fill(maxDistR);
	    if(maxDistEta>0)histos_["pfSC_FirstEtaBin_maxDistFromSeedinEtainSCPhoConv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["pfSC_FirstEtaBin_maxDistFromSeedinPhiinSCPhoConv"]->Fill(maxDistPhi);	
	  } else if(TMath::Abs(pfSCeta[i])>secondEtaBinDown && TMath::Abs(pfSCeta[i])<secondEtaBinUp){
	    if(maxDistR>0)histos_[  "pfSC_SecondEtaBin_maxDistFromSeedinRinSCPhoConv"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["pfSC_SecondEtaBin_maxDistFromSeedinEtainSCPhoConv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["pfSC_SecondEtaBin_maxDistFromSeedinPhiinSCPhoConv"]->Fill(maxDistPhi);
	  }else if(TMath::Abs(pfSCeta[i])>thirdEtaBinDown){
	    if(maxDistR>0)histos_[  "pfSC_ThirdEtaBin_maxDistFromSeedinRinSCPhoConv"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["pfSC_ThirdEtaBin_maxDistFromSeedinEtainSCPhoConv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["pfSC_ThirdEtaBin_maxDistFromSeedinPhiinSCPhoConv"]->Fill(maxDistPhi);
	  }
	}
      }//pfSCn
    }





    //loop on multi5x5SC
    if(multi5x5SCn>0){
      for (int i=0;i<multi5x5SCn;i++){
	

	//matching with gen ele
	bool matchesGenEle=false;
	bool matchesGenPhoConv=false;
	bool matchesGenPhoUnconv=false;

	bool matchesRecoEle=false;
	int  recoEleindexMatch=-1;
	float drRecoMatch=999;
	TLorentzVector pfscp4;
	pfscp4.SetPtEtaPhiE(multi5x5SCe[i]/cosh(multi5x5SCeta[i]),multi5x5SCeta[i],multi5x5SCphi[i],multi5x5SCe[i]);

	for(int ii=0;ii<elen;ii++){
	  if(elept[ii]<0.1)continue;
	  
	  //matching with reco ele
	  TLorentzVector elep4;
	  elep4.SetPtEtaPhiM(elept[ii],eleeta[ii],elephi[ii],0.5);
	  //      std::cout<<"etaele:"<<elep4.Eta(); 
	  float drReco=elep4.DeltaR(pfscp4);
	  
	  if(drReco<drRecoMatch && drReco<0.1){
	    recoEleindexMatch=ii;
	    drRecoMatch=drReco;
	    matchesRecoEle=true;
	  }
	}
	
	if(!matchesRecoEle)continue;

	TLorentzVector recoEleMatchp4;
	if(recoEleindexMatch!=-1) recoEleMatchp4.SetPtEtaPhiM(elept[recoEleindexMatch],eleeta[recoEleindexMatch],elephi[recoEleindexMatch],0.5);
	
	int geleindexMatch=-1;
	float drMatch=999;
	float geleMatchedE=-1;
	for(int j=0;j<gelen;j++){
	  if(gelept[j]<5)continue;
	  TLorentzVector gelep4;
	  gelep4.SetPtEtaPhiM(gelept[j],geleeta[j],gelephi[j],0.5);
	  
	  float dr=gelep4.DeltaR(recoEleMatchp4);
	  
	  if(dr<0.1 && dr<drMatch){
	    geleindexMatch=j;
	    drMatch=dr;
	    geleMatchedE=gelep4.E();
	  }
	}

	if(fabs(geleeta[geleindexMatch])<endcapBoundaryLow)continue;

	if(geleindexMatch != -1){
	  matchesGenEle=true;
	  //	  std::cout<<gelept[geleindexMatch]<< " " <<elept[recoEleindexMatch]<<" "<<multi5x5SCe[i]/cosh(multi5x5SCeta[i])<<std::endl;

	  //ereco

	  //	  std::cout<<multi5x5SCe[i]/(gelept[geleindexMatch]*cosh(geleeta[geleindexMatch]))<<std::endl;

	  if(TMath::Abs(multi5x5SCeta[i])<firstEtaBinUp){
	    histos_["multi5x5SC_ErecoOverETrueFirstEtaBin"]->Fill(multi5x5SCe[i]/geleMatchedE);
	    histos_["multi5x5SC_FirstEtaBin_nXtalsSeed"]->Fill(multi5x5SCnXtalsSeed[i]);
	    histos_["multi5x5SC_FirstEtaBin_nBCForSC"]->Fill(multi5x5SCnBC[i]);
	  }else if(TMath::Abs(multi5x5SCeta[i])>secondEtaBinDown && TMath::Abs(multi5x5SCeta[i])<secondEtaBinUp){
	    histos_["multi5x5SC_ErecoOverETrueSecondEtaBin"]->Fill(multi5x5SCe[i]/geleMatchedE);
	    histos_["multi5x5SC_SecondEtaBin_nXtalsSeed"]->Fill(multi5x5SCnXtalsSeed[i]);
	    histos_["multi5x5SC_SecondEtaBin_nBCForSC"]->Fill(multi5x5SCnBC[i]);
	  }else if(TMath::Abs(multi5x5SCeta[i])>thirdEtaBinDown){
	    histos_["multi5x5SC_ErecoOverETrueThirdEtaBin"]->Fill(multi5x5SCe[i]/geleMatchedE);
	    histos_["multi5x5SC_ThirdEtaBin_nXtalsSeed"]->Fill(multi5x5SCnXtalsSeed[i]);
	    histos_["multi5x5SC_ThirdEtaBin_nBCForSC"]->Fill(multi5x5SCnBC[i]);
	  }      


//	  std::cout<<"pfsc: pt,eta,phi "<<multi5x5SCe[i]/cosh(multi5x5SCeta[i])<<","<<multi5x5SCeta[i]<<","<<multi5x5SCphi[i]<<","<<multi5x5SCe[i]<<std::endl;
//	  std::cout<<"gele: pt,eta,phi "<<gelept[geleindexMatch]<<","<<geleeta[geleindexMatch]<<","<<gelephi[geleindexMatch]<<std::endl;
	}

	TLorentzVector matchedGenPhop4;
	int gphoindexMatch=-1;
	float drMatchPho=999;
	for(int j=0;j<gphon;j++){
	  if(gphopt[j]<5)continue;
	  
	  TLorentzVector gphop4;
	  gphop4.SetPtEtaPhiM(gphopt[j],gphoeta[j],gphophi[j],0.);
	  
	  float dr=pfscp4.DeltaR(gphop4);
	  
	  if(dr<0.1 && dr<drMatchPho){
	    gphoindexMatch=j;
	    drMatchPho=dr;
	    matchedGenPhop4=gphop4;
	  }
	}

	if(gphoindexMatch != -1) {
	  //FIX ME: for photons for the moment you don't have the info on conversion so you have to reloop on photons, try to add when producing ntuples
	  int recoPhoindexMatch=-1;
	  float drPhoMatch=999;

	    for(int ii=0;ii<phon;ii++){
	      if(phopt[ii]<0.1)continue;
	      TLorentzVector phop4;
	      phop4.SetPtEtaPhiM(phopt[ii],phoeta[ii],phophi[ii],0.);

	      float dr=matchedGenPhop4.DeltaR(phop4);

	      if(dr<0.1 && dr<drPhoMatch){
		recoPhoindexMatch=ii;
		drPhoMatch=dr;
	      }
	    }

	    if(recoPhoindexMatch>-1){
//	      std::cout<<"pfsc: pt,eta,phi "<<multi5x5SCe[i]/cosh(multi5x5SCeta[i])<<","<<multi5x5SCeta[i]<<","<<multi5x5SCphi[i]<<","<<multi5x5SCe[i]<<std::endl;
//	      std::cout<<"genPho: pt,eta,phi "<<gphopt[gphoindexMatch]<<","<<gphoeta[gphoindexMatch]<<","<<gphophi[gphoindexMatch]<<std::endl;
//	      std::cout<<"recoPho: pt,eta,phi "<<phoe[recoPhoindexMatch]/cosh(phoeta[recoPhoindexMatch])<<","<<phoeta[recoPhoindexMatch]<<","<<phophi[recoPhoindexMatch]<<","<<phoe[recoPhoindexMatch]<<std::endl;
	      if(phoE9[recoPhoindexMatch]/phoe[recoPhoindexMatch]>0.94) matchesGenPhoUnconv=true;
	      else matchesGenPhoConv=true;
	    }
	}


	if(!(matchesGenEle || matchesGenPhoUnconv || matchesGenPhoConv)) continue;

	//max distance from seed
	if(multi5x5SCbcE[i][0]<0.01)continue;
	//	std::cout<<multi5x5SCbcE[i][0]<<std::endl;
	float maxDistR=0;
	float maxDistEta=0;
	float maxDistPhi=0;


	TLorentzVector seed;
	seed.SetPtEtaPhiE(multi5x5SCbcE[i][0]/cosh(multi5x5SCbcEta[i][0]),multi5x5SCbcEta[i][0],multi5x5SCbcPhi[i][0],multi5x5SCbcE[i][0]);



	for (int j=0;j< multi5x5SCnBC[i];j++){
	  //	  std::cout<<multi5x5SCnBC[i]<< " i,j:"<<i<<","<<j<<" "<<multi5x5SCbcE[i][j]<<std::endl;
	  if(multi5x5SCbcE[i][j]<0.01)continue;
	  if(j>0){
	    TLorentzVector bc;
	    bc.SetPtEtaPhiE(multi5x5SCbcE[i][j]/cosh(multi5x5SCbcEta[i][j]),multi5x5SCbcEta[i][j],multi5x5SCbcPhi[i][j],multi5x5SCbcE[i][j]);
	    float distR=seed.DeltaR(bc);
	    if(distR>maxDistR)maxDistR=distR;
	    float distPhi=seed.DeltaPhi(bc);
	    if(distPhi>maxDistPhi)maxDistPhi=distPhi;
	    float distEta=sqrt(distR*distR-distPhi*distPhi);
	    if(distEta>maxDistEta)maxDistEta=distEta;
	    
	    //2d: ebcseed vs deltaphi(bc,seed)
	    if(matchesGenEle){
	      //	      std::cout<<fabs(distPhi)<<","<<multi5x5SCbcE[i][j]<<std::endl;
	      if(TMath::Abs(multi5x5SCbcEta[i][j])<firstEtaBinUp){
		histos2D_["multi5x5SC_FirstEtaBin_EBCseedVsDeltaPhiBCSeedEle"]->Fill(fabs(distPhi),multi5x5SCbcE[i][j]);
		histos2D_["multi5x5SC_FirstEtaBin_EBCseedVsDeltaEtaBCSeedEle"]->Fill(fabs(distEta),multi5x5SCbcE[i][j]);
	      }else if(TMath::Abs(multi5x5SCbcEta[i][j])>secondEtaBinDown && TMath::Abs(multi5x5SCbcEta[i][j])<secondEtaBinUp){
		histos2D_["multi5x5SC_SecondEtaBin_EBCseedVsDeltaPhiBCSeedEle"]->Fill(fabs(distPhi),multi5x5SCbcE[i][j]);
		histos2D_["multi5x5SC_SecondEtaBin_EBCseedVsDeltaEtaBCSeedEle"]->Fill(fabs(distEta),multi5x5SCbcE[i][j]);
	      }else if(TMath::Abs(multi5x5SCbcEta[i][j])>thirdEtaBinDown){
		histos2D_["multi5x5SC_ThirdEtaBin_EBCseedVsDeltaPhiBCSeedEle"]->Fill(fabs(distPhi),multi5x5SCbcE[i][j]);
		histos2D_["multi5x5SC_ThirdEtaBin_EBCseedVsDeltaEtaBCSeedEle"]->Fill(fabs(distEta),multi5x5SCbcE[i][j]);	
	      }
	    }
	    else if(matchesGenPhoUnconv){
	      if(TMath::Abs(multi5x5SCbcEta[i][j])<firstEtaBinUp){
		histos2D_["multi5x5SC_FirstEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv"]->Fill(fabs(distPhi),multi5x5SCbcE[i][j]);
		histos2D_["multi5x5SC_FirstEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv"]->Fill(fabs(distEta),multi5x5SCbcE[i][j]);
	      }else if(TMath::Abs(multi5x5SCbcEta[i][j])>secondEtaBinDown && TMath::Abs(multi5x5SCbcEta[i][j])<secondEtaBinUp){
		histos2D_["multi5x5SC_SecondEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv"]->Fill(fabs(distPhi),multi5x5SCbcE[i][j]);
		histos2D_["multi5x5SC_SecondEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv"]->Fill(fabs(distEta),multi5x5SCbcE[i][j]);
	      }else if(TMath::Abs(multi5x5SCbcEta[i][j])>thirdEtaBinDown){
		histos2D_["multi5x5SC_ThirdEtaBin_EBCseedVsDeltaPhiBCSeedPhoUnconv"]->Fill(fabs(distPhi),multi5x5SCbcE[i][j]);
		histos2D_["multi5x5SC_ThirdEtaBin_EBCseedVsDeltaEtaBCSeedPhoUnconv"]->Fill(fabs(distEta),multi5x5SCbcE[i][j]);
	      }
	    }else if(matchesGenPhoConv){
	      if(TMath::Abs(multi5x5SCbcEta[i][j])<firstEtaBinUp){
		histos2D_["multi5x5SC_FirstEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv"]->Fill(fabs(distPhi),multi5x5SCbcE[i][j]);
		histos2D_["multi5x5SC_FirstEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv"]->Fill(fabs(distEta),multi5x5SCbcE[i][j]);
	      }else if(TMath::Abs(multi5x5SCbcEta[i][j])>secondEtaBinDown && TMath::Abs(multi5x5SCbcEta[i][j])<secondEtaBinUp){
		histos2D_["multi5x5SC_SecondEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv"]->Fill(fabs(distPhi),multi5x5SCbcE[i][j]);
		histos2D_["multi5x5SC_SecondEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv"]->Fill(fabs(distEta),multi5x5SCbcE[i][j]);
	      }else if(TMath::Abs(multi5x5SCbcEta[i][j])>thirdEtaBinDown){
		histos2D_["multi5x5SC_ThirdEtaBin_EBCseedVsDeltaPhiBCSeedPhoConv"]->Fill(fabs(distPhi),multi5x5SCbcE[i][j]);
		histos2D_["multi5x5SC_ThirdEtaBin_EBCseedVsDeltaEtaBCSeedPhoConv"]->Fill(fabs(distEta),multi5x5SCbcE[i][j]);
	      }
	    }

	    //	    std::cout<<maxDistR<<std::endl;
	  }
	}//multi5x5SCnBC
	if(matchesGenEle){
	  if(TMath::Abs(multi5x5SCeta[i])<firstEtaBinUp){
	    if(maxDistR>0)histos_[  "multi5x5SC_FirstEtaBin_maxDistFromSeedinRinSCEle"]->Fill(maxDistR);
	    if(maxDistEta>0)histos_["multi5x5SC_FirstEtaBin_maxDistFromSeedinEtainSCEle"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["multi5x5SC_FirstEtaBin_maxDistFromSeedinPhiinSCEle"]->Fill(maxDistPhi);	
	  } else if(TMath::Abs(multi5x5SCeta[i])>secondEtaBinDown && TMath::Abs(multi5x5SCeta[i])<secondEtaBinUp){
	    if(maxDistR>0)histos_[  "multi5x5SC_SecondEtaBin_maxDistFromSeedinRinSCEle"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["multi5x5SC_SecondEtaBin_maxDistFromSeedinEtainSCEle"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["multi5x5SC_SecondEtaBin_maxDistFromSeedinPhiinSCEle"]->Fill(maxDistPhi);
	  }else if(TMath::Abs(multi5x5SCeta[i])>thirdEtaBinDown){
	    if(maxDistR>0)histos_[  "multi5x5SC_ThirdEtaBin_maxDistFromSeedinRinSCEle"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["multi5x5SC_ThirdEtaBin_maxDistFromSeedinEtainSCEle"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["multi5x5SC_ThirdEtaBin_maxDistFromSeedinPhiinSCEle"]->Fill(maxDistPhi);
	  }
	}else if (matchesGenPhoUnconv){
	  if(TMath::Abs(multi5x5SCeta[i])<firstEtaBinUp){
	    if(maxDistR>0)histos_[  "multi5x5SC_FirstEtaBin_maxDistFromSeedinRinSCPhoUnconv"]->Fill(maxDistR);
	    if(maxDistEta>0)histos_["multi5x5SC_FirstEtaBin_maxDistFromSeedinEtainSCPhoUnconv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["multi5x5SC_FirstEtaBin_maxDistFromSeedinPhiinSCPhoUnconv"]->Fill(maxDistPhi);
	  } else if(TMath::Abs(multi5x5SCeta[i])>secondEtaBinDown && TMath::Abs(multi5x5SCeta[i])<secondEtaBinUp){
	    if(maxDistR>0)histos_[  "multi5x5SC_SecondEtaBin_maxDistFromSeedinRinSCPhoUnconv"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["multi5x5SC_SecondEtaBin_maxDistFromSeedinEtainSCPhoUnconv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["multi5x5SC_SecondEtaBin_maxDistFromSeedinPhiinSCPhoUnconv"]->Fill(maxDistPhi);
	  }else if(TMath::Abs(multi5x5SCeta[i])>thirdEtaBinDown){
	    if(maxDistR>0)histos_[  "multi5x5SC_ThirdEtaBin_maxDistFromSeedinRinSCPhoUnconv"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["multi5x5SC_ThirdEtaBin_maxDistFromSeedinEtainSCPhoUnconv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["multi5x5SC_ThirdEtaBin_maxDistFromSeedinPhiinSCPhoUnconv"]->Fill(maxDistPhi);
	  }

	}else if (matchesGenPhoConv){
	  if(TMath::Abs(multi5x5SCeta[i])<firstEtaBinUp){
	    if(maxDistR>0)histos_[  "multi5x5SC_FirstEtaBin_maxDistFromSeedinRinSCPhoConv"]->Fill(maxDistR);
	    if(maxDistEta>0)histos_["multi5x5SC_FirstEtaBin_maxDistFromSeedinEtainSCPhoConv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["multi5x5SC_FirstEtaBin_maxDistFromSeedinPhiinSCPhoConv"]->Fill(maxDistPhi);	
	  } else if(TMath::Abs(multi5x5SCeta[i])>secondEtaBinDown && TMath::Abs(multi5x5SCeta[i])<secondEtaBinUp){
	    if(maxDistR>0)histos_[  "multi5x5SC_SecondEtaBin_maxDistFromSeedinRinSCPhoConv"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["multi5x5SC_SecondEtaBin_maxDistFromSeedinEtainSCPhoConv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["multi5x5SC_SecondEtaBin_maxDistFromSeedinPhiinSCPhoConv"]->Fill(maxDistPhi);
	  }else if(TMath::Abs(multi5x5SCeta[i])>thirdEtaBinDown){
	    if(maxDistR>0)histos_[  "multi5x5SC_ThirdEtaBin_maxDistFromSeedinRinSCPhoConv"]->Fill(maxDistR);    
	    if(maxDistEta>0)histos_["multi5x5SC_ThirdEtaBin_maxDistFromSeedinEtainSCPhoConv"]->Fill(maxDistEta);
	    if(maxDistPhi>0)histos_["multi5x5SC_ThirdEtaBin_maxDistFromSeedinPhiinSCPhoConv"]->Fill(maxDistPhi);
	  }
	}
      }//multi5x5SCn
    }

	

  }//jentry

  
}

//void createHistos::setInputFile(TString inputFileString){
//  inputFile_=TFile::Open(inputFileString);
//  Init((TTree*)inputFile_->Get("tree"));
//
//}

void createHistos::setOutputFile(TString outputFileString){
  outputFile_=TFile::Open(outputFileString,"recreate");
}

createHistos::createHistos(TChain* chain){
  Init(chain);
}

createHistos::~createHistos()
{
  outputFile_->Write();
  outputFile_->Close();
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t createHistos::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t createHistos::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void createHistos::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("vertexx", vertexx, &b_vertexx);
   fChain->SetBranchAddress("vertexy", vertexy, &b_vertexy);
   fChain->SetBranchAddress("vertexz", vertexz, &b_vertexz);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("phon", &phon, &b_phon);
   fChain->SetBranchAddress("phopt", phopt, &b_phopt);
   fChain->SetBranchAddress("phoeta", phoeta, &b_phoeta);
   fChain->SetBranchAddress("phophi", phophi, &b_phophi);
   fChain->SetBranchAddress("phoe", phoe, &b_phoe);
   fChain->SetBranchAddress("phoesc", phoesc, &b_phoesc);
   fChain->SetBranchAddress("phoetasc", phoetasc, &b_phoetasc);
   fChain->SetBranchAddress("phophisc", phophisc, &b_phophisc);
   fChain->SetBranchAddress("phoxsc", phoxsc, &b_phoxsc);
   fChain->SetBranchAddress("phoysc", phoysc, &b_phoysc);
   fChain->SetBranchAddress("phozsc", phozsc, &b_phozsc);
   fChain->SetBranchAddress("phoxcalo", phoxcalo, &b_phoxcalo);
   fChain->SetBranchAddress("phoycalo", phoycalo, &b_phoycalo);
   fChain->SetBranchAddress("phozcalo", phozcalo, &b_phozcalo);
   fChain->SetBranchAddress("phohasPixelSeed", phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoisEB", phoisEB, &b_phoisEB);
   fChain->SetBranchAddress("phoisEE", phoisEE, &b_phoisEE);
   fChain->SetBranchAddress("phoisEBEEGap", phoisEBEEGap, &b_phoisEBEEGap);
   fChain->SetBranchAddress("phoE9", phoE9, &b_phoE9);
   fChain->SetBranchAddress("phoE25", phoE25, &b_phoE25);
   fChain->SetBranchAddress("phojurECAL", phojurECAL, &b_phojurECAL);
   fChain->SetBranchAddress("photwrHCAL", photwrHCAL, &b_photwrHCAL);
   fChain->SetBranchAddress("phoHoverE", phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phohlwTrack", phohlwTrack, &b_phohlwTrack);
   fChain->SetBranchAddress("phosiEtaiEtaZS", phosiEtaiEtaZS, &b_phosiEtaiEtaZS);
   fChain->SetBranchAddress("phosiEtaiEtaNoZS", phosiEtaiEtaNoZS, &b_phosiEtaiEtaNoZS);
   fChain->SetBranchAddress("phosiEtaiEtaWrong", phosiEtaiEtaWrong, &b_phosiEtaiEtaWrong);
   fChain->SetBranchAddress("phosiPhiiPhiZS", phosiPhiiPhiZS, &b_phosiPhiiPhiZS);
   fChain->SetBranchAddress("phosiPhiiPhiNoZS", phosiPhiiPhiNoZS, &b_phosiPhiiPhiNoZS);
   fChain->SetBranchAddress("phosiPhiiPhiWrong", phosiPhiiPhiWrong, &b_phosiPhiiPhiWrong);
   fChain->SetBranchAddress("phosMajZS", phosMajZS, &b_phosMajZS);
   fChain->SetBranchAddress("phosMajNoZS", phosMajNoZS, &b_phosMajNoZS);
   fChain->SetBranchAddress("phosMinZS", phosMinZS, &b_phosMinZS);
   fChain->SetBranchAddress("phosMinNoZS", phosMinNoZS, &b_phosMinNoZS);
   fChain->SetBranchAddress("phoalphaZS", phoalphaZS, &b_phoalphaZS);
   fChain->SetBranchAddress("phoalphaNoZS", phoalphaNoZS, &b_phoalphaNoZS);
   fChain->SetBranchAddress("phoscetawid", phoscetawid, &b_phoscetawid);
   fChain->SetBranchAddress("phoscphiwid", phoscphiwid, &b_phoscphiwid);
   fChain->SetBranchAddress("phojurECAL03", phojurECAL03, &b_phojurECAL03);
   fChain->SetBranchAddress("photwrHCAL03", photwrHCAL03, &b_photwrHCAL03);
   fChain->SetBranchAddress("phohlwTrack03", phohlwTrack03, &b_phohlwTrack03);
   fChain->SetBranchAddress("phopfCandPt", phopfCandPt, &b_phopfCandPt);
   fChain->SetBranchAddress("phopfCandEta", phopfCandEta, &b_phopfCandEta);
   fChain->SetBranchAddress("phopfCandPhi", phopfCandPhi, &b_phopfCandPhi);
   fChain->SetBranchAddress("phopfCandVtx", phopfCandVtx, &b_phopfCandVtx);
   fChain->SetBranchAddress("phopfCandType", phopfCandType, &b_phopfCandType);
   fChain->SetBranchAddress("phoptiso004", phoptiso004, &b_phoptiso004);
   fChain->SetBranchAddress("phontrkiso004", phontrkiso004, &b_phontrkiso004);
   fChain->SetBranchAddress("phoptiso035", phoptiso035, &b_phoptiso035);
   fChain->SetBranchAddress("phontrkiso035", phontrkiso035, &b_phontrkiso035);
   fChain->SetBranchAddress("phoptiso04", phoptiso04, &b_phoptiso04);
   fChain->SetBranchAddress("phontrkiso04", phontrkiso04, &b_phontrkiso04);
   fChain->SetBranchAddress("phoPfSumChargedHadronPt", phoPfSumChargedHadronPt, &b_phoPfSumChargedHadronPt);
   fChain->SetBranchAddress("phoPfSumNeutralHadronEt", phoPfSumNeutralHadronEt, &b_phoPfSumNeutralHadronEt);
   fChain->SetBranchAddress("phoPfSumPhotonEt", phoPfSumPhotonEt, &b_phoPfSumPhotonEt);
   fChain->SetBranchAddress("pfSCn", &pfSCn, &b_pfSCn);
   fChain->SetBranchAddress("pfSCeta", pfSCeta, &b_pfSCeta);
   fChain->SetBranchAddress("pfSCphi", pfSCphi, &b_pfSCphi);
   fChain->SetBranchAddress("pfSCe", pfSCe, &b_pfSCe);
   fChain->SetBranchAddress("pfSCnBC", pfSCnBC, &b_pfSCnBC);
   fChain->SetBranchAddress("pfSCnXtalsSeed", pfSCnXtalsSeed, &b_pfSCnXtalsSeed);
   fChain->SetBranchAddress("pfSCnXtalsTotal", pfSCnXtalsTotal, &b_pfSCnXtalsTotal);
   fChain->SetBranchAddress("pfSCbcEta", pfSCbcEta, &b_pfSCbcEta);
   fChain->SetBranchAddress("pfSCbcPhi", pfSCbcPhi, &b_pfSCbcPhi);
   fChain->SetBranchAddress("pfSCbcE", pfSCbcE, &b_pfSCbcE);
   fChain->SetBranchAddress("pfSCbcNXtals", pfSCbcNXtals, &b_pfSCbcNXtals);
   fChain->SetBranchAddress("multi5x5SCn", &multi5x5SCn, &b_multi5x5SCn);
   fChain->SetBranchAddress("multi5x5SCeta", &multi5x5SCeta, &b_multi5x5SCeta);
   fChain->SetBranchAddress("multi5x5SCphi", &multi5x5SCphi, &b_multi5x5SCphi);
   fChain->SetBranchAddress("multi5x5SCe", &multi5x5SCe, &b_multi5x5SCe);
   fChain->SetBranchAddress("multi5x5SCnBC", &multi5x5SCnBC, &b_multi5x5SCnBC);
   fChain->SetBranchAddress("multi5x5SCnXtalsSeed", &multi5x5SCnXtalsSeed, &b_multi5x5SCnXtalsSeed);
   fChain->SetBranchAddress("multi5x5SCnXtalsTotal", &multi5x5SCnXtalsTotal, &b_multi5x5SCnXtalsTotal);
   fChain->SetBranchAddress("multi5x5SCbcEta", multi5x5SCbcEta, &b_multi5x5SCbcEta);
   fChain->SetBranchAddress("multi5x5SCbcPhi", multi5x5SCbcPhi, &b_multi5x5SCbcPhi);
   fChain->SetBranchAddress("multi5x5SCbcE", multi5x5SCbcE, &b_multi5x5SCbcE);
   fChain->SetBranchAddress("multi5x5SCbcNXtals", multi5x5SCbcNXtals, &b_multi5x5SCbcNXtals);
   fChain->SetBranchAddress("hybridSCn", &hybridSCn, &b_hybridSCn);
   fChain->SetBranchAddress("hybridSCeta", hybridSCeta, &b_hybridSCeta);
   fChain->SetBranchAddress("hybridSCphi", hybridSCphi, &b_hybridSCphi);
   fChain->SetBranchAddress("hybridSCe", hybridSCe, &b_hybridSCe);
   fChain->SetBranchAddress("hybridSCnBC", hybridSCnBC, &b_hybridSCnBC);
   fChain->SetBranchAddress("hybridSCnXtalsSeed", hybridSCnXtalsSeed, &b_hybridSCnXtalsSeed);
   fChain->SetBranchAddress("hybridSCnXtalsTotal", hybridSCnXtalsTotal, &b_hybridSCnXtalsTotal);
   fChain->SetBranchAddress("hybridSCbcEta", hybridSCbcEta, &b_hybridSCbcEta);
   fChain->SetBranchAddress("hybridSCbcPhi", hybridSCbcPhi, &b_hybridSCbcPhi);
   fChain->SetBranchAddress("hybridSCbcE", hybridSCbcE, &b_hybridSCbcE);
   fChain->SetBranchAddress("hybridSCbcNXtals", hybridSCbcNXtals, &b_hybridSCbcNXtals);
   fChain->SetBranchAddress("elen", &elen, &b_elen);
   fChain->SetBranchAddress("elepx", elepx, &b_elepx);
   fChain->SetBranchAddress("elepy", elepy, &b_elepy);
   fChain->SetBranchAddress("elepz", elepz, &b_elepz);
   fChain->SetBranchAddress("elevx", elevx, &b_elevx);
   fChain->SetBranchAddress("elevy", elevy, &b_elevy);
   fChain->SetBranchAddress("elevz", elevz, &b_elevz);
   fChain->SetBranchAddress("elept", elept, &b_elept);
   fChain->SetBranchAddress("eleeta", eleeta, &b_eleeta);
   fChain->SetBranchAddress("elephi", elephi, &b_elephi);
   fChain->SetBranchAddress("elee", elee, &b_elee);
   fChain->SetBranchAddress("eleecalEnergy", eleecalEnergy, &b_eleecalEnergy);
   fChain->SetBranchAddress("eletrackPatVtx", eletrackPatVtx, &b_eletrackPatVtx);
   fChain->SetBranchAddress("eletrackAtVtxPt", eletrackAtVtxPt, &b_eletrackAtVtxPt);
   fChain->SetBranchAddress("eletrackAtVtxEta", eletrackAtVtxEta, &b_eletrackAtVtxEta);
   fChain->SetBranchAddress("eletrackAtVtxPhi", eletrackAtVtxPhi, &b_eletrackAtVtxPhi);
   fChain->SetBranchAddress("eletrackAtCaloPt", eletrackAtCaloPt, &b_eletrackAtCaloPt);
   fChain->SetBranchAddress("eletrackAtCaloEta", eletrackAtCaloEta, &b_eletrackAtCaloEta);
   fChain->SetBranchAddress("eletrackAtCaloPhi", eletrackAtCaloPhi, &b_eletrackAtCaloPhi);
   fChain->SetBranchAddress("eleesc", eleesc, &b_eleesc);
   fChain->SetBranchAddress("eleetasc", eleetasc, &b_eleetasc);
   fChain->SetBranchAddress("elephisc", elephisc, &b_elephisc);
   fChain->SetBranchAddress("elexsc", elexsc, &b_elexsc);
   fChain->SetBranchAddress("eleysc", eleysc, &b_eleysc);
   fChain->SetBranchAddress("elezsc", elezsc, &b_elezsc);
   fChain->SetBranchAddress("elexcalo", elexcalo, &b_elexcalo);
   fChain->SetBranchAddress("eleycalo", eleycalo, &b_eleycalo);
   fChain->SetBranchAddress("elezcalo", elezcalo, &b_elezcalo);
   fChain->SetBranchAddress("eleisEB", eleisEB, &b_eleisEB);
   fChain->SetBranchAddress("eleisEE", eleisEE, &b_eleisEE);
   fChain->SetBranchAddress("eleisEBEEGap", eleisEBEEGap, &b_eleisEBEEGap);
   fChain->SetBranchAddress("eleE9", eleE9, &b_eleE9);
   fChain->SetBranchAddress("eleE25", eleE25, &b_eleE25);
   fChain->SetBranchAddress("elecharge", elecharge, &b_elecharge);
   fChain->SetBranchAddress("eleFbrem", eleFbrem, &b_eleFbrem);
   fChain->SetBranchAddress("eleScFbrem", eleScFbrem, &b_eleScFbrem);
   fChain->SetBranchAddress("elePfScFbrem", elePfScFbrem, &b_elePfScFbrem);
   fChain->SetBranchAddress("eleESeedClusterOverPout", eleESeedClusterOverPout, &b_eleESeedClusterOverPout);
   fChain->SetBranchAddress("eledist", eledist, &b_eledist);
   fChain->SetBranchAddress("eledcot", eledcot, &b_eledcot);
   fChain->SetBranchAddress("elemisHits", elemisHits, &b_elemisHits);
   fChain->SetBranchAddress("elematchedConv", elematchedConv, &b_elematchedConv);
   fChain->SetBranchAddress("eleseedType", eleseedType, &b_eleseedType);
   fChain->SetBranchAddress("eleEoP", eleEoP, &b_eleEoP);
   fChain->SetBranchAddress("eleOneOverEMinusOneOverP", eleOneOverEMinusOneOverP, &b_eleOneOverEMinusOneOverP);
   fChain->SetBranchAddress("eler9", eler9, &b_eler9);
   fChain->SetBranchAddress("elenSubClusters", elenSubClusters, &b_elenSubClusters);
   fChain->SetBranchAddress("eletrkIso", eletrkIso, &b_eletrkIso);
   fChain->SetBranchAddress("eleecalIso", eleecalIso, &b_eleecalIso);
   fChain->SetBranchAddress("elehcalIso", elehcalIso, &b_elehcalIso);
   fChain->SetBranchAddress("eletrkIso03", eletrkIso03, &b_eletrkIso03);
   fChain->SetBranchAddress("eleecalIso03", eleecalIso03, &b_eleecalIso03);
   fChain->SetBranchAddress("elehcalIso03", elehcalIso03, &b_elehcalIso03);
   fChain->SetBranchAddress("elesiEtaiEtaZS", elesiEtaiEtaZS, &b_elesiEtaiEtaZS);
   fChain->SetBranchAddress("elesiEtaiEtaNoZS", elesiEtaiEtaNoZS, &b_elesiEtaiEtaNoZS);
   fChain->SetBranchAddress("elesiEtaiEtaWrong", elesiEtaiEtaWrong, &b_elesiEtaiEtaWrong);
   fChain->SetBranchAddress("elesiPhiiPhiZS", elesiPhiiPhiZS, &b_elesiPhiiPhiZS);
   fChain->SetBranchAddress("elesiPhiiPhiNoZS", elesiPhiiPhiNoZS, &b_elesiPhiiPhiNoZS);
   fChain->SetBranchAddress("elesiPhiiPhiWrong", elesiPhiiPhiWrong, &b_elesiPhiiPhiWrong);
   fChain->SetBranchAddress("elesMajZS", elesMajZS, &b_elesMajZS);
   fChain->SetBranchAddress("elesMajNoZS", elesMajNoZS, &b_elesMajNoZS);
   fChain->SetBranchAddress("elesMinZS", elesMinZS, &b_elesMinZS);
   fChain->SetBranchAddress("elesMinNoZS", elesMinNoZS, &b_elesMinNoZS);
   fChain->SetBranchAddress("elealphaZS", elealphaZS, &b_elealphaZS);
   fChain->SetBranchAddress("elealphaNoZS", elealphaNoZS, &b_elealphaNoZS);
   fChain->SetBranchAddress("elepfCandPt", elepfCandPt, &b_elepfCandPt);
   fChain->SetBranchAddress("elepfCandEta", elepfCandEta, &b_elepfCandEta);
   fChain->SetBranchAddress("elepfCandPhi", elepfCandPhi, &b_elepfCandPhi);
   fChain->SetBranchAddress("elepfCandVtx", elepfCandVtx, &b_elepfCandVtx);
   fChain->SetBranchAddress("elepfCandType", elepfCandType, &b_elepfCandType);
   fChain->SetBranchAddress("eledEtaIn", eledEtaIn, &b_eledEtaIn);
   fChain->SetBranchAddress("eledPhiIn", eledPhiIn, &b_eledPhiIn);
   fChain->SetBranchAddress("eleHoE", eleHoE, &b_eleHoE);
   fChain->SetBranchAddress("elepFlowMVA", elepFlowMVA, &b_elepFlowMVA);
   fChain->SetBranchAddress("eleScEnergy", eleScEnergy, &b_eleScEnergy);
   fChain->SetBranchAddress("eleScEta", eleScEta, &b_eleScEta);
   fChain->SetBranchAddress("eleScPhi", eleScPhi, &b_eleScPhi);
   fChain->SetBranchAddress("elePfSumChargedHadronPt", elePfSumChargedHadronPt, &b_elePfSumChargedHadronPt);
   fChain->SetBranchAddress("elePfSumNeutralHadronEt", elePfSumNeutralHadronEt, &b_elePfSumNeutralHadronEt);
   fChain->SetBranchAddress("elePfSumPhotonEt", elePfSumPhotonEt, &b_elePfSumPhotonEt);
   fChain->SetBranchAddress("gpn", &gpn, &b_gpn);
   fChain->SetBranchAddress("gppt", gppt, &b_gppt);
   fChain->SetBranchAddress("gpeta", gpeta, &b_gpeta);
   fChain->SetBranchAddress("gpphi", gpphi, &b_gpphi);
   fChain->SetBranchAddress("gpidMC", gpidMC, &b_gpidMC);
   fChain->SetBranchAddress("gpstatusMC", gpstatusMC, &b_gpstatusMC);
   fChain->SetBranchAddress("gpmotherIdMC", gpmotherIdMC, &b_gpmotherIdMC);
   fChain->SetBranchAddress("gphon", &gphon, &b_gphon);
   fChain->SetBranchAddress("gphopt", gphopt, &b_gphopt);
   fChain->SetBranchAddress("gphoeta", gphoeta, &b_gphoeta);
   fChain->SetBranchAddress("gphophi", gphophi, &b_gphophi);
   fChain->SetBranchAddress("gphoindex", gphoindex, &b_gphoindex);
   fChain->SetBranchAddress("gelen", &gelen, &b_gelen);
   fChain->SetBranchAddress("gelept", gelept, &b_gelept);
   fChain->SetBranchAddress("geleeta", geleeta, &b_geleeta);
   fChain->SetBranchAddress("gelephi", gelephi, &b_gelephi);
   fChain->SetBranchAddress("gelefbrem80", gelefbrem80, &b_gelefbrem80);
   fChain->SetBranchAddress("gelefbrem120", gelefbrem120, &b_gelefbrem120);
   fChain->SetBranchAddress("geleindex", geleindex, &b_geleindex);
   fChain->SetBranchAddress("truePU", &truePU, &b_truePU);
   fChain->SetBranchAddress("bxPU", bxPU, &b_bxPU);
   Notify();
}

Bool_t createHistos::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void createHistos::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t createHistos::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
