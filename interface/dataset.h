#ifndef __DATASET__
#define __DATASET__

#include "covariance.h"
#include "tensor.h"
#include "template.h"
#include "mtrand.h"
#include "time.h"
#include <vector>
#include <cmath>
#include <iostream>

/*! 
  \brief a class to store multi-dimensional data
  

  A class for multi-dimension histogram
  included fucntions to compute Multivariate Kernel smoothing
*/

class DataSet{
 public:

  DataSet(){}
  //!< empty constructor

  DataSet(const vector<double>& min, const vector<double>& max,
	  const vector<double>& precision);
  //!< constructor using vectors
  //!< supply min & max range

  int Rank() const
  {
    return data.rank;
  }
  //!< return rank of the dataset

  inline vector<int> GetBin(const vector<double>& ntuple) const;
  //!< return bin for a given ntuple
  //!< returns the index in 1D vector

  void Fill(const vector<double>& ntuple, double weight=1.0);
  //!< fill dataset by a data sample

  Template ComputeTemplate(double c=1.0, bool verbose_mode=false);
  //!< core code to generate a Template based on multivariate Kernel Smoothing
  //!< fast algorithm implementing discrete fourier transform

  DataSet GenerateDataSet(double scale_factor=1.0) const;
  //!< a class to generate new dataset by fluctuate each bin
  //!< extra scale_factor may be supplied

  static void SetRandomSeed(int seed)
  {
    irand.seed(seed);
  }
  //!< function to set a fixed random seed

  friend ostream& operator<<(ostream& os, const DataSet&);
  //!< streaming operator for debugging purposes

  vector<double> one_over_width;
  //!< vector to store width info

  vector<double> width;
  //!< vector to store width info

  vector<double> min;
  //!< vector to store min values

  vector<double> max;
  //!< vector to store max values

  vector<int> nbins;
  //!< vector to store number of bins

  vector<double> precision;
  //!< vector to store precision of variables

  static bool b_bias_only;
  //!< setting to enable bias only in computation of Template

  static double Random_Gaus()
  {
    //we've already generated another random number
    if(b_random_gaus)
      {
	b_random_gaus=false;
	return next_gaus;
      }

    //else generate two new normal random number
    b_random_gaus=true;
    
    double s,u,v;
    double num1=0;
    double num2=0;

    while(true)
      {
	u=irand.rand_double()*2-1;
	v=irand.rand_double()*2-1;
	
	s=u*u+v*v;
	
	if(s<1.0)
	  {
	    //num1=num2=0
	    if(s>0)
	      {
		double part=sqrt(-2*log(s)/s);
		//convert into integers
		num1=static_cast<int>(u*part+0.5);
		num2=static_cast<int>(v*part+0.5);
	      }
	    break;
	  }
      }
    
    //num1 and num2 are independ normally distribution random number
    next_gaus=num2;
    return num1;      
  }
  //!< output a normally distributed random variable

  
  static double Random_Poisson(double lambda)
  {
    if(lambda > 30)
      {
	double result=-1;
	double mean=lambda;
	double var=sqrt(lambda);
	
	//make sure things are non-negative
	while(result<0)
	  {
	    result= Random_Gaus()*var + mean;
	  }
	return result;
      }

    //if lambda is small, do poisson
    //by counting in second, how many decays are in
    double time=1.0;
    int count=0;
    double exp_lambda=exp(-lambda);

    while(time > exp_lambda)
      {
	//generate an Exponential number
	double random=irand.rand_double();
	if(random==0) continue;
	time*=random;
	if(time > exp_lambda)
	  count++;
      }
    return count;
  }
  //!< output a Poisson distributed random number (integer)
  //!< if lambda > 20, output a normally distributed number 
  //(non-integer and non-negative)


 private:
  Tensor<double> data;
  //!< to store multidimensional info

  Tensor<double> data2;
  //!< to store multidimensional info (store value^2)
  //!< needed to compute error
  
  Covariance co;
  //!< to store covariance matrix
  //!< and bandwidth info for Kernel Smoothing

  void Initialize(int nvar, const double* min, const double* max, 
		  const double* precision_=NULL);
  //!< unsafe way to initialize a dataset
  //!< supply number of variables, min & max range
  //!< precision (default to 0.01 for all variables)

  void Initialize(const vector<double>& min, const vector<double>& max, 
		  const vector<double>& precision);
  //!< unsafe way to initialize a dataset
  //!< supply number of variables, min & max range
  //!< precision (default to 0.01 for all variables)

  static MTRand_int32 irand;
  //!< random number generator

  static bool b_random_gaus;
  //!< boolean to keep track of Gaussian random variables

  static double next_gaus;
  //!< the next gaussian random variable
  //!< they are calculated two at a time

  double fix_eff_n;
  //!< variable used to produced smeared dataset
  //!< with a fixed number of events


  //private inline functions
  inline double GetBinCenter(int index, int bin) const;
  //!< compute value centered at bin (a particular dimension)
  
  inline vector<double> GetBinCenter(const vector<int>& index) const;
  //!< compute value array centered at an vectorial index 
  //!< i.e. pt, m, t21...etc

  inline int GetBin_helper(double ntuple, int bin) const;
  //!< return bin for a given ntuple (in multidimension)
  //!< helper function for GetBin
  //!< UNSAEF function, assume length of ntuple is same as rank  

  void SetN(double eff_n);
  //!< fixed the number of events
  //!< used for producing smeared Template
};


inline double DataSet::GetBinCenter(int index, int bin) const
{
  //DEBUG
  /*
  cout<<"GetBin small input : "<<index<<","<<bin<<endl;
  cout<<"GetBin small output : "<<(index+0.5)*(max[bin]-min[bin])<<","<<min[bin]<<endl;
  */

  return (index+0.5)*width[bin]+min[bin];
}

inline vector<double> DataSet::GetBinCenter(const vector<int>& index) const
{
  assert(index.size()==data.rank);
  vector<double> result(data.rank);
  
  //grab central value of bin in each dimension
  for(unsigned int i=0; i<result.size(); i++)
    {
      result[i]=GetBinCenter(index[i], i);
    }

  return result;
}

inline int DataSet::GetBin_helper(double value, int bin) const
{
  double temp_bin= (value-min[bin])*one_over_width[bin];
  //if value is in underflowed, go to the zeroth bin
  if(temp_bin<0) return 0;
  
  //if value is overflowed, go to the last bin
  if(temp_bin>=nbins[bin]) return nbins[bin]-1;
  return static_cast<int>(temp_bin);
}

inline vector<int> DataSet::GetBin(const vector<double>& ntuple) const
{
  assert(ntuple.size()==data.rank);
  vector<int> result_bin(data.rank);
  
  //convert multimensional bin to 1D binning
  for(int i=0; i<data.rank; i++)
    {
      result_bin[i]=GetBin_helper(ntuple[i], i);
    }
  return result_bin;
}


#endif
