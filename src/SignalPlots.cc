#include "SignalPlots.h"

SignalPlots::SignalPlots(TString name) {
  h_leadingJetPt      =   new TH1F("h_leadingJetPt_"           + name,"leading jet pt",200,0,2000);
  h_leadingJetPtLog   =   new TH1F("h_leadingJetPtLog_"        + name,"leading jet pt/200 log",300,0,3);
  h_leadingJetMass    =   new TH1F("h_leadingJetMass_"         + name,"leading jet mass",150,0,1500);
  h_leadingJetMassLog =   new TH1F("h_leadingJetMassLog_"      + name,"leading jet mas/pt log",300,0,3);

  h_leadingJetPtvsMassLog =   new TH2F("h_leadingJetPtvsMassLog_"  + name,"leading jet pt vs mass",300,0,3,300,0,3);
  internalName = name;
}

SignalPlots::~SignalPlots() {
  delete h_leadingJetPt;
  delete h_leadingJetPtLog;
  delete h_leadingJetMass;
  delete h_leadingJetMassLog;
  delete h_leadingJetPtvsMassLog;
}

void SignalPlots::Fill(Double_t pt, Double_t mass, Double_t weight, Bool_t nolog) {
  if (nolog) {
    h_leadingJetPt->Fill(pt, weight);
    h_leadingJetPtLog->Fill(log(pt/200.), weight);
    h_leadingJetMass->Fill(mass, weight);
    h_leadingJetMassLog->Fill(-log10(mass/pt), weight);
      
    h_leadingJetPtvsMassLog->Fill(-log10(mass/pt), log(pt/200.), weight);
  }
  else {
    h_leadingJetPtLog->Fill(pt, weight);
    h_leadingJetPt->Fill( exp(pt)*200., weight);
    h_leadingJetMassLog->Fill(mass, weight);
    h_leadingJetMass->Fill( pow(10,mass)*pt, weight);
      
    h_leadingJetPtvsMassLog->Fill( mass, pt, weight);
  }
}

void SignalPlots::Fill(Template MyPdf) {
  for (UInt_t i=1; i<=h_leadingJetPtvsMassLog->GetNbinsX(); i++)
    for (UInt_t j=1; j<=h_leadingJetPtvsMassLog->GetNbinsY(); j++) {
      h_leadingJetPtvsMassLog->SetBinContent(i, j,  MyPdf( make_vector<double>(h_leadingJetPtvsMassLog->GetXaxis()->GetBinCenter(i), h_leadingJetPtvsMassLog->GetYaxis()->GetBinCenter(j) ) ) ); 
  }
  h_leadingJetPtLog=(TH1F*)h_leadingJetPtvsMassLog->ProjectionY("h_leadingJetPtLog_"+internalName);
  h_leadingJetMass =(TH1F*)h_leadingJetPtvsMassLog->ProjectionX("h_leadingJetMass_" +internalName);
  
}

void SignalPlots::Write() {
  h_leadingJetPt->Write();
  h_leadingJetPtLog->Write();
  h_leadingJetMass->Write();
  h_leadingJetMassLog->Write();
  h_leadingJetPtvsMassLog->Write();
}
