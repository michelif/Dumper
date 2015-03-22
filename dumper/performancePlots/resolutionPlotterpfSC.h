//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Mar 22 16:04:57 2015 by ROOT version 5.34/07
// from TTree outTreepfSC/outTreepfSC
// found on file: dummyPho.root
//////////////////////////////////////////////////////////

#ifndef resolutionPlotterpfSC_h
#define resolutionPlotterpfSC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class resolutionPlotterpfSC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   TFile* outFile_;


   // Declaration of leaf types
   Float_t         pfSCpt;
   Float_t         pfSCe;
   Float_t         pfSCeta;
   Float_t         pfSCphi;
   Float_t         pfSCErecoOverEtrue;
   Float_t         pfSC_nXtalsSeed;
   Float_t         pfSC_nXtalsTotal;
   Float_t         pfSC_nBCforSC;
   Float_t         pfSCEseedOverEtrue;
   Float_t         pfSCEtrue;
   Float_t         pfSCptTrue;
   Float_t         pfSCisEle;
   Float_t         pfSCisConv;
   Float_t         pfSCfBrem;
   Float_t         pfSCR9;

   // List of branches
   TBranch        *b_pfSCpt;   //!
   TBranch        *b_pfSCe;   //!
   TBranch        *b_pfSCeta;   //!
   TBranch        *b_pfSCphi;   //!
   TBranch        *b_pfSCErecoOverEtrue;   //!
   TBranch        *b_pfSC_nXtalsSeed;   //!
   TBranch        *b_pfSC_nXtalsTotal;   //!
   TBranch        *b_pfSC_nBCforSC;   //!
   TBranch        *b_pfSCEseedOverEtrue;   //!
   TBranch        *b_pfSCEtrue;   //!
   TBranch        *b_pfSCptTrue;   //!
   TBranch        *b_pfSCisEle;   //!
   TBranch        *b_pfSCisConv;   //!
   TBranch        *b_pfSCfBrem;   //!
   TBranch        *b_pfSCR9;   //!

   resolutionPlotterpfSC(TTree *tree=0);
   virtual ~resolutionPlotterpfSC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   Double_t effSigma(TH1 * hist);
   void setOutFile(TString outFileName);


};

#endif

#ifdef resolutionPlotterpfSC_cxx
resolutionPlotterpfSC::resolutionPlotterpfSC(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("dummyPho.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("dummyPho.root");
      }
      f->GetObject("outTreepfSC",tree);

   }
   Init(tree);
}

resolutionPlotterpfSC::~resolutionPlotterpfSC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t resolutionPlotterpfSC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t resolutionPlotterpfSC::LoadTree(Long64_t entry)
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

void resolutionPlotterpfSC::Init(TTree *tree)
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
 
   fChain->SetBranchAddress("pfSCpt", &pfSCpt, &b_pfSCpt);
   fChain->SetBranchAddress("pfSCe", &pfSCe, &b_pfSCe);
   fChain->SetBranchAddress("pfSCeta", &pfSCeta, &b_pfSCeta);
   fChain->SetBranchAddress("pfSCphi", &pfSCphi, &b_pfSCphi);
   fChain->SetBranchAddress("pfSCErecoOverEtrue", &pfSCErecoOverEtrue, &b_pfSCErecoOverEtrue);
   fChain->SetBranchAddress("pfSC_nXtalsSeed", &pfSC_nXtalsSeed, &b_pfSC_nXtalsSeed);
   fChain->SetBranchAddress("pfSC_nXtalsTotal", &pfSC_nXtalsTotal, &b_pfSC_nXtalsTotal);
   fChain->SetBranchAddress("pfSC_nBCforSC", &pfSC_nBCforSC, &b_pfSC_nBCforSC);
   fChain->SetBranchAddress("pfSCEseedOverEtrue", &pfSCEseedOverEtrue, &b_pfSCEseedOverEtrue);
   fChain->SetBranchAddress("pfSCEtrue", &pfSCEtrue, &b_pfSCEtrue);
   fChain->SetBranchAddress("pfSCptTrue", &pfSCptTrue, &b_pfSCptTrue);
   fChain->SetBranchAddress("pfSCisEle", &pfSCisEle, &b_pfSCisEle);
   fChain->SetBranchAddress("pfSCisConv", &pfSCisConv, &b_pfSCisConv);
   fChain->SetBranchAddress("pfSCfBrem", &pfSCfBrem, &b_pfSCfBrem);
   fChain->SetBranchAddress("pfSCR9", &pfSCR9, &b_pfSCR9);
   
   Notify();
}

Bool_t resolutionPlotterpfSC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void resolutionPlotterpfSC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t resolutionPlotterpfSC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void resolutionPlotterpfSC::setOutFile(TString outFileName){
  outFile_=TFile::Open(outFileName,"recreate");
}

#endif // #ifdef resolutionPlotterpfSC_cxx
