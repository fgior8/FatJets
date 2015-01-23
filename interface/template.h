#ifndef __Template__
#define __Template__

#include "tensor.h"
#include "mtrand.h"
#include <iostream>
#include <vector>

using namespace std;

/*! 
  \brief Template functions


  this class contains a table of templates
  values are interpolated between bins.
*/
  
class Template{
 public:
  Template(){valid=false;}
  //!< Template must be initialize externally, a trivial
  //!< constructor is presented but Template will be invalid
  //!< until constructed by DataSet class or read from a file

  double operator()(const vector<double>& ntuple, double* err=NULL) const;
  //!< evaluate the Template at a location
  //!< using linear interpolation
  //!< also obtain error if a pointer is specified

  bool IsValid() const
  {return valid;}
  //!< returns true if Template is valid

  int Rank() const
  {return data.Rank();}

  friend ostream& operator<<(ostream&, const Template&);
  //!< streaming operator for debugging

  friend class DataSet;
  //!< allow DataSet to read private members

  Tensor<double> data;
  //!< store Template value at a multidimensional grid

  /*! enum for setting which is input and out in Dresser
   */
  
  enum Flag{
    INPUT, OUTPUT
  };

  /*! \brief Monte Carlo Event


    the Dresser generates DressedEvents,
    DressedEvents store weights and errors for each event
    the erros ARE NOT INDEPENDENT, treat them as another 
    variable and compute them by averaging, not adding
    in quadrature
   */

  class DressedEvent;

  /*! \brief Monte Carlo Integrator


    this is used to generate monte carlo events based
    on the given Template. Slices can be taken by the 
    appropriate settings
   */

  class Dresser;

  Dresser Generator(const vector<Flag>& flag_) const;
  //!< returns a Dresser given some flags


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

 private:
  
  Tensor<double> data_bias;
  //!< Template value bias

  bool valid;
  //!< to keep track of whether not the Template is valid

  void Initialize(int nvar, const double* min, const double* max,
		  const double* precision_);
  //!< initialize the Template
  //!< Template could still invalid!

  //private inline functions
  inline double GetBinCenter(int index, int bin) const;
  //!< compute value centered at bin (a particular dimension)
  
  inline vector<double> GetBinCenter(const vector<int>& index) const;
  //!< compute value array centered at an vectorial index 
  //!< i.e. pt, m, t21...etc

  inline int GetBin_helper(double value, int bin) const;
  //!< return bin given value
};


//more detailed implementation of the Dresser class
class Template::Dresser{
 public:
  
  Dresser(const vector<Flag>& flag_, const Template* template_);
  //!< default constructor, fix the flags and Template
  
  vector<DressedEvent> Generate(const vector<double>& input, int nevt=1000) const;
  //!< generate a number of Monte Carlo events

  vector<DressedEvent> Generate(int nevt=1000) const
    {
      return Generate(vector<double>(), nevt);
    }
  //!< generate a number of Monte Carlo events (without input)

  vector<DressedEvent> GenerateFast(const vector<double>& input, int nevt=1000) 
    const;
  //!< generate a number of Monte Carlo events, much faster but binned

  vector<DressedEvent> GenerateFast(int nevt=1000) const
    {
      return Generate(vector<double>(), nevt);
    }
  //!< generate a number of Monte Carlo events (without input)

  vector<DressedEvent> GenerateFull(const vector<double>& input) const;
  //!< this does not give MC events, but rather simply dump the whole 
  //!< Template into DressedEvents

			 
 private:
  vector<Flag> flag;
  //!< settings for which variables to generate MC
  
  const Template* ttemplate;
  //!< the corresponding Template

  int input_size;
  //!< size of the input
  
  vector<int> index;
  //!< which index to generate Monte Carlo

  static MTRand_int32 irand;
  //!< random number generator
  
};


//more detailed implementation of the DressedEvent class
class Template::DressedEvent{
 public:
  
  vector<double> extra;
  //!< extra variable

  DressedEvent(const vector<double>& value__, 
	   double bias__, double rho__):
  value_(value__), bias_(bias__), rho_(rho__) {}
  //!< default constructor

  vector<double>& value()
  {return value_;}

  double rho() const
  {return rho_;}

  double rho_star() const
  {return rho_-bias_;}

  double bias() const
  {return bias_;}



  void SetExtra(const vector<double>& extra_)
  {
    extra=extra_;
  }
  
  vector<double>& GetExtra()
  {
    return extra;
  }

  double operator[](int i) const
  {
    assert(i>=0 && i<value_.size());
    return value_[i];
  }
  //!< access operator (checks for bounds)

 private:
  
  double bias_;
  //!< bias for the event

  double rho_;
  //!< weight for the event

  vector<double> value_;
  //!< value for the monte carlo



  friend class Template;
  

};


inline double Template::GetBinCenter(int index, int bin) const
{
  //DEBUG
  /*
  cout<<"GetBin small input : "<<index<<","<<bin<<endl;
  cout<<"GetBin small output : "<<(index+0.5)*(max[bin]-min[bin])<<","<<min[bin]<<endl;
  */

  return (index+0.5)*width[bin]+min[bin];
}

inline vector<double> Template::GetBinCenter(const vector<int>& index) const
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


inline int Template::GetBin_helper(double value, int bin) const
{
  double temp_bin= (value-min[bin])*one_over_width[bin];
  //if value is in underflowed, go to the zeroth bin
  if(temp_bin<0) return 0;
  
  //if value is overflowed, go to the last bin
  if(temp_bin>=nbins[bin]) return nbins[bin]-1;
  return static_cast<int>(temp_bin);
}


#endif
