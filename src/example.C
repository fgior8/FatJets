#include <iostream>
#include <cmath>
//random number generator library
#include "mtrand.h"

//make_list have functions to make vectors easily
#include "make_vector.h"

#include "dataset.h"
#include "template.h"

using namespace std;

//initialize random number generator with seed 123
MTRand_int32 irand(123);

//we need a test function for the kernel smoothing
double test_func(double x)
{
  //function is (1-x)^2 * x^4 / (2-x)
  // x from 0 to 1, normalized to one
  double normalization = 142.418582618187;
  
  return normalization*(1-x)*(1-x)*x*x*x*x/(2-x);
}

//generate monte carlo for test function
double generate_test_func()
{
  //max of function
  double max=2.36866988642;
  
  double gen_x, gen_fx;

  do {
    gen_x=irand.rand_double();
    gen_fx=irand.rand_double()*max;
  }  while(test_func(gen_x) < gen_fx);

  return gen_x;
}

int main()
{
  cout<<"This example illustrates a few ways kernelsmoother can be used"<<endl;


  cout<<endl;
  cout<<"Example 1: kernel smoothing with N=1000"<<endl;
  //let's try a kernel smoothing for n=1000 events
  //declare a dataset first (we use make_vector in the make_vector.h file)
  DataSet mydata(make_vector<double>(0.0),
		 make_vector<double>(1.0),
		 make_vector<double>(0.001));
  
  //generate some data
  for(int i=0; i<1000; ++i)
    mydata.Fill(make_vector<double>(generate_test_func()));


  //make a pdf
  //rule_of_thumb > 1.0 = bigger bandwidth = oversmoothing
  double rule_of_thumb=1.0;
  Template mypdf=mydata.ComputeTemplate(rule_of_thumb, true);

  //compare the results
  cout<<"x, rho(x), rho_hat(x)"<<endl;
  for(double i=0; i<=1.0; i+=0.1)
    cout<<i<<", "<<test_func(i)<<", "<<mypdf(make_vector<double>(i))<<endl;

  
  cout<<endl;
  cout<<"Example 2: kernel smoothing with N=1000, with errors"<<endl;
  cout<<"for x=0.5, compute the bias corrected rho_hat, and variance"<<endl;

  //here we need to store the kernel smoothing estimate for 500 results
  //to compute the fluctuation at a given x
  double bias=0;
  double rho_hat=mypdf(make_vector<double>(0.5), &bias);
  double rho_star=rho_hat-bias;

  //generate smeared pdf for variance computation
  double s=0;
  double s2=0;
  for(int i=0; i<1000; ++i)
    {
      DataSet mydata_smeared = mydata.GenerateDataSet();
      Template mypdf_smeared = mydata_smeared.ComputeTemplate();
      double bias_smeared = 0;
      double rho_hat_smeared = mypdf_smeared(make_vector<double>(0.5), 
					     &bias_smeared);
      double rho_star_smeared = rho_hat_smeared - bias_smeared;
      
      s+= rho_star_smeared;
      s2+= rho_star_smeared*rho_star_smeared;
    }

  s/=1000;
  s2/=1000;
  double var= sqrt((s2-s*s)*1000/999);
  
  cout<<"rho(0.5) = "<<test_func(0.5)<<endl;
  cout<<"rho_hat_star(0.5) = "<<rho_star<<endl;
  cout<<"bias_hat(0.5) = "<<bias<<endl;
  cout<<"var_star(0.5) = "<<var<<endl;


  //example 3
  cout<<endl;
  cout<<"Example 3: full analysis computing efficiencies for x < 0.5"<<endl;

  //first compute the true efficiency
  double true_eff=0;
  for(int i=0; i<1e6; ++i)
    if(generate_test_func() < 0.5)
      true_eff+= 1e-6;

  cout<<"true eff: "<<true_eff<<endl;

  //then compute the efficiency estimate, bias corrected
  //first declare a MC generator, specifying with is input/output
  //in our simple case, all the entries are output
  Template::Dresser mygen = mypdf.Generator
    (make_vector<Template::Flag>(Template::OUTPUT));

  //generate enough events for MC integration
  //enough events mean the variance from MC integration is smaller
  //than the smoothing variance
  vector<Template::DressedEvent> myevt = mygen.Generate(1e4);

  double myeff=0;
  double myeff_normalization=0;

  //keep track of uncorrected efficiency
  double myeff_un=0;
  double myeff_normalization_un=0;  

  for(int i=0; i<myevt.size(); ++i)
    {
      if(myevt[i].value()[0] < 0.5)
	myeff += myevt[i].rho_star();
      myeff_normalization += myevt[i].rho_star();

      //keep track of uncorrected eff, for bias correction computation
      if(myevt[i].value()[0] < 0.5)
	myeff_un += myevt[i].rho();
      myeff_normalization_un += myevt[i].rho();

    }

  myeff/=myeff_normalization;
  myeff_un /= myeff_normalization_un;

  cout<<"eff hat: "<<myeff<<endl;
  cout<<"bias correction: "<<myeff_un - myeff<<endl;

  //now do the variance computation
  //repeat the analysis 100 times with smeared dataset
  s=0;
  s2=0;
  for(int i=0; i<100; ++i)
    {
      DataSet mydata_smeared = mydata.GenerateDataSet();
      Template mypdf_smeared = mydata_smeared.ComputeTemplate();
      
      Template::Dresser mygen_smeared = mypdf_smeared.Generator
	(make_vector<Template::Flag>(Template::OUTPUT));

      //generate lots of events for MC integration
      vector<Template::DressedEvent> myevt_smeared = mygen_smeared.Generate(1e4);

      //compute the efficiency
      double eff_temp=0;
      double eff_norm_temp=0;
      
      for(int i=0; i<myevt_smeared.size(); ++i)
	{
	  if(myevt_smeared[i].value()[0] < 0.5)
	    eff_temp += myevt_smeared[i].rho_star();
	  eff_norm_temp += myevt_smeared[i].rho_star();
	}
      eff_temp /= eff_norm_temp;
      s += eff_temp;
      s2 += eff_temp*eff_temp;
    }

  s/=100;
  s2/=100;

  var=sqrt((s2-s*s)*100/99);
  cout<<"var_star: "<<var<<endl;
  cout<<endl;


  cout<<"program ends."<<endl;
  return 1;
}
