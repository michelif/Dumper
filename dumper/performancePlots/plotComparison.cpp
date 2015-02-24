#include <stdlib.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
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

  int nLists=argc-1;  

  //  Input list 1
  char listName1[500];
  char listName2[500];
  char listName3[500];
  char listName4[500];
  sprintf(listName1,argv[1]); 

  // Input list 2
  if (nLists == 2){
  sprintf(listName2,argv[2]);
  }

  // Input list 3
  if (nLists == 3){
  sprintf(listName3,argv[3]);
  }

  // Input list 4
  if (nLists == 4){
  sprintf(listName4,argv[4]);
  }

  // Output filename (.root)  
  TString OutputFileName(argv[nLists+1]);
  
  // Name of input tree objects in (.root) files 
  char treeName[100] = "tree";

  //creating  TChain
  TChain *chain = new TChain(treeName);
  char pName[500];
  ifstream is(listName);
  if(! is.good()) {
    cout << "int main() >> ERROR : file " << listName << " not read" << endl;
    is.close();
    exit(-1);
  }
  cout << "Reading list : " << listName << " ......." << endl;
  
  while( is.getline(pName, 500, '\n') ) {
    if (pName[0] == '#') continue;
    chain->Add(pName); 
  }
  is.close();

  return 0;

}
