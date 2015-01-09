#include "Analyzer.h"

Analyzer::Analyzer() {

  if (debug) cout<<"inizio"<<endl;

  h_control   = new SignalPlots("control");
  h_smoothing = new SignalPlots("smoothing");
  h_dressed   = new SignalPlots("dressed");
  h_signal    = new SignalPlots("signal");

  
  h_onevar = new TH1F("h_onevar","m / p_T",300,0,3);
  h_onevar->SetDefaultSumw2(true);
  h_onevarSmall = new TH1F("h_onevarSmall","m / p_T",300,0,3);
  h_onevarSmooth = new TH1F("h_onevarSmooth","m / p_T",300,0,3);

  h_twovar = new TH2F("h_twovar","m / p_T, p/200",30,0,3,30,0,3);
  h_twovarSmall = new TH2F("h_twovarSmall","m / p_T, p/200",30,0,3,30,0,3);
  h_twovarSmooth = new TH2F("h_twovarSmooth","m / p_T, p/200",30,0,3,30,0,3);
  h_twovarDress = new TH2F("h_twovarDress","m / p_T, p/200",30,0,3,30,0,3);

  h_prova = new TH1F("h_prova", "prova", 30,0,3);

  if (debug) cout<<"fine"<<endl;
}

Analyzer::~Analyzer() { }

void Analyzer::SetName(TString name, Int_t version) {
  completename = name + "_";
  completename += version;
  completename += ".root";
  outfile = new TFile(completename,"RECREATE");
}

void Analyzer::SetWeight(Double_t CrossSection, Double_t nevents) {

  MCweight = integratedlumi * CrossSection / nevents;
// lumi *  cs(pb) * gen filter efficiency / MCevents
  cout<<"MCweight = "<<MCweight<<endl;
 
}

void Analyzer::SetEvtN(Long64_t events) {
  events ? entrieslimit=events :  entrieslimit=-1;
  cout<<"events "<<events<<endl<<"entrieslimit "<<entrieslimit<<endl;
}

void Analyzer::Loop() {
  DataSet mydata(make_vector<double>(0.0), make_vector<double>(3.0), make_vector<double>(0.01));
  DataSet mydataTwo(make_vector<double>(0.0, 0.0), make_vector<double>(3.0, 3.0), make_vector<double>(0.01, 0.01));

  vector<Double_t> jetpt;
  cout << "total number of entries " <<nentries<<endl;

  if (debug) cout<< "loop begins" <<endl;

  fBTagSF = new BTagSFUtil("CSVM");

//  reweightPU = new ReweightPU("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/MyDataPileupHistogram_69400.root");

  if (debug) cout<< "PU histos loaded" <<endl;


  if(!MCweight) MCweight=1; 

  weight=MCweight;

  if (fChain == 0) 
    cout << "Ciao!" << endl;

//  cout << "Do you want limited events?" <<endl;
//  cin >> entrieslimit;
//  if (entrieslimit != -1)
//    nentries=entrieslimit;
//  entrieslimit = 1000000;

  if (debug) cout<< "at the loop" <<endl;
  std::set<int> runs;
  for (Long64_t jentry = 0; jentry < nentries; jentry++ ) {

    //    watch_getentry.Start(false);
    if (debug) cout<< "Event number " <<jentry<<endl;
    if (debug) cout<<"begin loop"<<endl;
    if (!(jentry % 10000)) 
      cout << jentry << endl;

    if (!fChain) cout<<"porcaccia"<<endl;
    fChain->GetEntry(jentry);
    //control region
    if (NJets_ak1p2Jets>=2 && NJets_ak1p2Jets<=3) {      
      if (pt_ak1p2Jets->at(0) > 200 && mass_ak1p2Jets->at(0)>1 && pt_ak1p2Jets->at(1) > 100) {

	h_onevarSmall->Fill( -log10( mass_ak1p2Jets->at(0)/pt_ak1p2Jets->at(0) ) );
	h_twovarSmall->Fill( -log10( mass_ak1p2Jets->at(0)/pt_ak1p2Jets->at(0) ), log( pt_ak1p2Jets->at(0)/200. ) );

	h_control->Fill( pt_ak1p2Jets->at(0), mass_ak1p2Jets->at(0) );

	mydata.Fill( make_vector<double>( -log10( mass_ak1p2Jets->at(0)/pt_ak1p2Jets->at(0) ) ) );
	mydataTwo.Fill( make_vector<double>( -log10( mass_ak1p2Jets->at(0)/pt_ak1p2Jets->at(0) ), log( pt_ak1p2Jets->at(0)/200. ) ) );

      }
    }
    //signal region
    if (NJets_ak1p2Jets>3) {
      if( pt_ak1p2Jets->at(0) > 200 && mass_ak1p2Jets->at(0)>1 && pt_ak1p2Jets->at(1) > 100) {
	h_onevar->Fill( -log10( mass_ak1p2Jets->at(0)/pt_ak1p2Jets->at(0) ) );
	h_twovar->Fill( -log10( mass_ak1p2Jets->at(0)/pt_ak1p2Jets->at(0) ), log( pt_ak1p2Jets->at(0)/200. ) );

	h_signal->Fill( pt_ak1p2Jets->at(0), mass_ak1p2Jets->at(0) );
	
	jetpt.push_back(pt_ak1p2Jets->at(0));
      }
    } 

  }
  if(debug) cout<< "out of the loop" <<endl;
  
  cout<<"Example 1: kernel smoothing with N=1000"<<endl;

  //let's try a kernel smoothing for n=1000 events
  //declare a dataset first (we use make_vector in the make_vector.h file)
  
  

  //generate some data
  // for(int i=0; i<100; ++i) 
  //mydata.Fill(make_vector<double>( i ), h_prova_small->GetBinContent(i+1)/h_prova_small->Integral() );
		//cout << "pro " << generate_test_func() << endl;
  //  cout << "checking " << mydata.GetBin(make_vector<double>(0.1)) << endl;
  
  //make a pdf
  //rule_of_thumb > 1.0 = bigger bandwidth = oversmoothing
  double rule_of_thumb=10.;
  Template mypdf=mydata.ComputeTemplate(rule_of_thumb);

  Template mypdfTwo=mydataTwo.ComputeTemplate(rule_of_thumb);


  //compare the results
  cout<<"x, rho(x), rho_hat(x)"<<endl;
  for(int i=1; i<=300; i++) {
    //  cout<<i<<", "<<h_onevarSmall->GetBinContent(i+1)/h_onevarSmall->Integral()<<", "<<mypdf(make_vector<double>( h_onevarSmall->GetBinCenter(i) ))<<endl;
    h_onevarSmooth->SetBinContent(i,mypdf(make_vector<double>( h_onevarSmall->GetBinCenter(i) )));
  }

  for(int i=1; i<=30; i++) 
    for (int j=1; j<=30; j++){
      //   cout<<i<<", "<< j << ", "<<h_twovarSmall->GetBinContent(i+1,j+1)/h_twovarSmall->Integral()<<", "<<mypdfTwo(make_vector<double>( h_twovarSmall->GetXaxis()->GetBinCenter(i), h_twovarSmall->GetYaxis()->GetBinCenter(j) ))<<endl;
      h_twovarSmooth->SetBinContent(i, j, mypdfTwo(make_vector<double>( h_twovarSmall->GetXaxis()->GetBinCenter(i), h_twovarSmall->GetYaxis()->GetBinCenter(j) )));

  }
/*
  for(int i=0; i<300; i++) 
    for (int j=0; j<300; j++){
      
      h_smoothing->Fill( j/100., i/100., mypdfTwo(make_vector<double>( i/100., j/100. )), false);

    }
  */
  h_smoothing->Fill( mypdfTwo );
   //example 3
  cout<<endl;
  cout<<"Example 3: full analysis computing efficiencies for x < 0.5"<<endl;

  //first compute the true efficiency
  double true_eff=h_twovar->Integral(10,30,0,30)/h_twovar->Integral();


  cout<<"true eff: "<<true_eff<<endl;

  //then compute the efficiency estimate, bias corrected
  //first declare a MC generator, specifying with is input/output
  //in our simple case, all the entries are output
  Template::Dresser mygenTwo = mypdfTwo.Generator(make_vector<Template::Flag>(Template::OUTPUT, Template::INPUT));

  //generate enough events for MC integration
  //enough events mean the variance from MC integration is smaller
  //than the smoothing variance
  vector<Template::DressedEvent> myevtTwo;
  double myeff=0;
  double myeff_normalization=0;
  //keep track of uncorrected efficiency
  double myeff_un=0;
  double myeff_normalization_un=0;
  double bias_star = 0;
  cout << " jetpt.size() " << jetpt.size() << endl;
  for (Long64_t jentry = 0; jentry < jetpt.size(); jentry++ ) {
    if (jetpt[jentry]>200) {
      myevtTwo = mygenTwo.Generate(make_vector<double>( log(jetpt[jentry]/200.) ), 100);
      h_prova->Fill(log(jetpt[jentry]/200.));
      // myevt = mygen.Generate(1e4);
      for(int i=0; i<myevtTwo.size(); ++i)
	{ 
	  //	  if (myevtTwo[i].rho_star() > 0.02 && log(jetpt[jentry]/200.)>1.) cout << "AZZ " << myevtTwo[i].rho_star() << " and "  << log(jetpt[jentry]/200.) << " at jentry " << jentry << endl;
	  Double_t cazzo = myevtTwo[i].rho_star();
	  if (cazzo <= 0.) {
	    // cout << "CAZZO " << cazzo << " so the weight is " << myevtTwo[i].rho() << " and star " << myevtTwo[i].rho_star() << " for pt " << jetpt[jentry] << endl;
	    continue;
	  }
	  h_twovarDress->Fill( myevtTwo[i].value()[0], log(jetpt[jentry]/200.), myevtTwo[i].rho_star() );
	  h_dressed->Fill(jetpt[jentry], pow(10,-(myevtTwo[i].value()[0]))*jetpt[jentry], myevtTwo[i].rho_star() );
	  bias_star +=  myevtTwo[i].bias();
	  if(myevtTwo[i].value()[0] > h_twovar->GetXaxis()->GetBinLowEdge(10) )
	    myeff += myevtTwo[i].rho_star();
	  myeff_normalization += myevtTwo[i].rho_star();	  
	  //keep track of uncorrected eff, for bias correction computation
	  if(myevtTwo[i].value()[0] > h_twovar->GetXaxis()->GetBinLowEdge(10) )
	    myeff_un += myevtTwo[i].rho();
	  myeff_normalization_un += myevtTwo[i].rho();	
	}
    }
  }
  h_twovarDress->Scale(1./myeff_normalization);
  myeff/=myeff_normalization;
  myeff_un /= myeff_normalization_un;
  cout<<"eff hat: "<<myeff<<endl;
  cout<<"bias correction: "<<myeff_un - myeff<<endl;

  //now do the variance computation
  //repeat the analysis 100 times with smeared dataset
  double s=0;
  double s2=0;
 /* 
  for(int i=0; i<10; ++i)
    {
      DataSet mydata_smearedTwo = mydataTwo.GenerateDataSet();
      Template mypdf_smearedTwo = mydata_smearedTwo.ComputeTemplate();
      
      Template::Dresser mygen_smearedTwo = mypdf_smearedTwo.Generator(make_vector<Template::Flag>(Template::OUTPUT, Template::INPUT));

      //generate lots of events for MC integration
      // vector<Template::DressedEvent> myevt_smeared = mygen_smeared.Generate(1e4);
      vector<Template::DressedEvent> myevt_smearedTwo;
      //compute the efficiency
      double eff_temp=0;
      double eff_norm_temp=0;
      for (Long64_t jentry = 0; jentry < jetpt.size(); jentry++ ) {
	if (jetpt[jentry]>200) {
	  myevt_smearedTwo = mygen_smearedTwo.Generate(make_vector<double>( log(jetpt[jentry]/200.) ),100);
	  
	  for(int ev=0; ev<myevt_smearedTwo.size(); ++ev)
	    {
	      if(myevt_smearedTwo[ev].value()[0] > h_twovar->GetXaxis()->GetBinLowEdge(10))
		eff_temp += myevt_smearedTwo[ev].rho_star();
	      eff_norm_temp += myevt_smearedTwo[ev].rho_star();
	    }
	}
      }
      if (i%(100/100) == 0) cout << "Computing variation " << (i+1.)/100.*100. << "%" <<endl;
      eff_temp /= eff_norm_temp;
      s += eff_temp;
      s2 += eff_temp*eff_temp;
    }
  s/=10;
  s2/=10;
  cout << " s = " << s << "  s2 = " << s2 << " s2-s*s = " << s2-s*s << "  and (s2-s*s)*100./99. = " << (s2-s*s)*100./99. << endl;
  var=sqrt(fabs(s2-s*s)*10./9.);
  cout<<"var_star: "<<var<<endl;
  cout<<endl;
*/
  
  cout<<"program ends."<<endl;

  h_control->Write();
  h_smoothing->Write();
  h_dressed->Write();
  h_signal->Write();
  
  h_onevar->Write();
  h_onevarSmall->Write();
  h_onevarSmooth->Write();

  h_twovar->Write();
  h_twovarSmall->Write();
  h_twovarSmooth->Write();
  h_twovarDress->Write();

  h_prova->Write();
  outfile->Close();

}

