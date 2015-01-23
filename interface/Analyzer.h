#ifndef Analyzer_h
#define Analyzer_h

#include <set>
#include "Data.h"
#include "TH2F.h"
#include "Reweight.cc"
#include "BTagSFUtil.h"
#include "SignalPlots.h"
//#include "OtherFunctions.h"

#include <iostream>
#include <cmath>
//random number generator library
#include "mtrand.h"

//make_list have functions to make vectors easily
#include "make_vector.h"

#include "dataset.h"
#include "template.h"


using namespace std;

class Analyzer : public Data {


  static const Bool_t debug = false; 
  static const Double_t integratedlumi = 19762.501;
//  static const Double_t integratedlumi = 1.0; //for Fakes
  static const Double_t Mass_Z = 91.1876;
  static const Double_t Mass_W = 80.398;

  TString completename;
  TFile *outfile;
  Long64_t entrieslimit;
  ReweightPU *reweightPU;
  BTagSFUtil *fBTagSF;


  Double_t MCweight, weight;

  SignalPlots *h_control, *h_smoothing, *h_dressed, *h_signal;
  TH1F *h_prova;
  TH1F *h_onevar, *h_onevarSmall, *h_onevarSmooth; 
  TH2F *h_twovar, *h_twovarSmall, *h_twovarSmooth, *h_twovarDress; 

 public:
  static const Bool_t MC_pu = false; 


  Analyzer();
  ~Analyzer();
  void Loop();
  void SetWeight(Double_t CrossSection, Double_t nevents);
  void SetName(TString name, Int_t version);
  void SetEvtN(Long64_t events);

};
#endif
