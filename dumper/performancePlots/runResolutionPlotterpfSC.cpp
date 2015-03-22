#include <stdlib.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include "resolutionPlotterpfSC.h"

using namespace std;



int main( int argc, char* argv[] ) {



  if(argc<2 || argc>3) {
    cout << "Usage:  ./tmp/runResolutionPlotterpfSC inputFile outputfile\n"
	 << "  listfile:    list of root files incusing protocol eg dcap:/// .....\n"
	 << "  outputfile:  name of output root file  eg output.root\n"
	 << endl;
    exit(-1);
  }


  TString inputFileName(argv[1]);
  TFile * inputFile=TFile::Open(inputFileName);

  // Output filename (.root)  
  TString OutputFileName(argv[2]);

    
  //creating  TChain
  TTree *chain = (TTree*)inputFile->Get("outTreepfSC");

 
  resolutionPlotterpfSC t(chain);
  t.setOutFile(OutputFileName);
  t.Loop();

}
