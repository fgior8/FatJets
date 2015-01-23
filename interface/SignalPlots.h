#ifndef SignalPlots_h
#define SignalPlots_h

#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <iostream>
#include <cmath>

#include "template.h"
#include "make_vector.h"

class SignalPlots {
  
  TString internalName;
 public:
  TH1F *h_leadingJetPt, *h_leadingJetPtLog, *h_leadingJetMass, *h_leadingJetMassLog;
  TH2F *h_leadingJetPtvsMassLog;
 
  SignalPlots(TString name);
  ~SignalPlots();
  void Fill(Double_t pt, Double_t mass, Double_t weight=1., Bool_t nolog=true);
  void Fill(Template MyPdf);
  void Write();

};


#endif
