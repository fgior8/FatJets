#include "Analyzer.h"
#include "ChainMaker.h"
#include "mtrand.h"

void display_usage() {
  cout << "simple instrunctions" << endl;
}

int main (int argc, const char* argv[]) {

  const char * _output   = 0;
  const char * _input    = 0;
  const char * _dir      = "/uscms_data/d2/fgior8/SubStructure/CMSSW_5_3_14_patch2/src/AWhitbeck/DataSetList/";
  
  cout << " Hello " << endl;
  // Arguments used
  //std::set<int> usedargs;
  //Parsing input options
  if(argc == 1){
    display_usage();
    return -1;
  }
  else{
    //Argument 1 must be a valid input fileName
    for (int i = 1; i < argc; i++){
      if( strcmp(argv[i],"-i") == 0 ){
        _input = argv[i+1];
        i++;
      }
      if( strcmp(argv[i],"-o") == 0 ){
        _output= argv[i+1];
        i++;
      }
      if( strcmp(argv[i],"-d") == 0 ){
        _dir = argv[i+1];
        i++;
      }
      if( strcmp(argv[i],"-h") == 0 ||
          strcmp(argv[i],"--help") == 0 ){
        display_usage();
        return 0;
      }
    }
  }//else
  if( _input ==0){
    std::cerr << "\033[1;31mskimfile ERROR:\033[1;m The '-i' option is mandatory!"
              << std::endl;
    display_usage();
    return -1;
  }

cout << "ma " << endl;
  // reassigning
  TString fname(_input);
  TString hname(_output);
  TString fdir(_dir);

  Analyzer Pippo;
cout << "fname " << fname << endl;
  fdir.Append(fname);
  cout << "Running on " << fdir << endl;
  TChain* chain = ChainMaker(fdir);
  Pippo.Init(chain); 
  Pippo.SetName(hname,2);
  Pippo.SetWeight(1, 1);
  cout << "Saved in " << hname << endl;  
  Pippo.Loop();

  if (0) {
    Analyzer Pippo;  
    TChain* chain = ChainMaker("/uscms_data/d2/fgior8/SubStructure/CMSSW_5_3_14_patch2/src/AWhitbeck/DataSetList/ZJets.txt");
    Pippo.Init(chain); 
    Pippo.SetName("ZJets",1);
    Pippo.SetWeight(1, 1);
    cout << "ZJets" << endl;  
    Pippo.Loop();
    cout<<"This example illustrates a few ways kernelsmoother can be used"<<endl;
  }
  if (0) {
    Analyzer Pippo;
    TChain* chain = ChainMaker("/uscms_data/d2/fgior8/SubStructure/CMSSW_5_3_14_patch2/src/AWhitbeck/DataSetList/WJets.txt");
    Pippo.Init(chain);
    Pippo.SetName("WJets",1);
    Pippo.SetWeight(1, 1);
    cout << "WJets" << endl;
    Pippo.Loop();
  }
  if (0) {
    Analyzer Pippo;
    TChain* chain = ChainMaker("/uscms_data/d2/fgior8/SubStructure/CMSSW_5_3_14_patch2/src/AWhitbeck/DataSetList/TTJets.txt");
    Pippo.Init(chain);
    Pippo.SetName("TTbarJets",1);
    Pippo.SetWeight(1, 1);
    cout << "ttbar" << endl;
    Pippo.Loop();
  }
  if (0) {
    Analyzer Pippo;
    TChain* chain = ChainMaker("/uscms_data/d2/fgior8/SubStructure/CMSSW_5_3_14_patch2/src/AWhitbeck/DataSetList/QCD_1000.txt");
    Pippo.Init(chain);
    Pippo.SetName("QCD_1000",1);
    Pippo.SetWeight(1, 1);
    cout << "QCD_1000" << endl;
    Pippo.Loop();
  } cout << endl << endl;
  if (0) {
    Analyzer Pippo;
    TChain* chain = ChainMaker("/uscms_data/d2/fgior8/SubStructure/CMSSW_5_3_14_patch2/src/AWhitbeck/DataSetList/QCD_500HT1000.txt");
    Pippo.Init(chain);
    Pippo.SetName("QCD_500HT1000",1);
    Pippo.SetWeight(1, 1);
    cout << "QCD_500HT1000" << endl;
    Pippo.Loop();
  }
  if (0) {
    Analyzer Pippo;
    TChain* chain = ChainMaker("/uscms_data/d2/fgior8/SubStructure/CMSSW_5_3_14_patch2/src/AWhitbeck/DataSetList/SMS_T1tttt_mGo1400_mLSP25.txt");
    Pippo.Init(chain);
    Pippo.SetName("SMS_T1tttt_mGo1400_mLSP25",1);
    Pippo.SetWeight(1, 1);
    cout << "SMS_T1tttt_mGo1400_mLSP25" << endl;
    Pippo.Loop();
  }  
cout << endl << endl;
}
