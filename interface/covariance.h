#ifndef __COVARIANCE__
#define __COVARIANCE__

#include <vector>
#include <iostream>

using namespace std;

/*!
  \brief class for Covariance matrix


  Implementation of 2D symmetric square matrix via vector,
  used to store and compute covariance matrix of multiple variables
*/

class Covariance
{
 public:

  vector<double> bandwidth; 
  //!< bandwidth matrix^2 computed from ary2
  //!< done using Silverman's Rule of Thumb

  vector<double> bandwidth_inverse; 
  //!< bandwidth inverse^2
  //!< done using Silverman's Rule of Thumb
  //!< in computation, bandwidth_inverse is used instead of bandwidth

  
  Covariance(int rank_=1);
  //!< default constructor, supply the rank of matrix
 
  int Rank() const
  {return rank;}
  //!< return rank of tensor

  void ComputeBandwidth(double min_size=-1.0, double extra_factor=1.0,
			double fix_eff_n=-1.0);
  //!< Compute Bandwidth Matrix using Silverman's Rule
  //!< argument determines the minimum smoothing size
  //!< this will automatically modify error computations

  double BandwidthNorm(const vector<double>& band, int power=-1);
  //!< compute double^T*bandwidth^(2*power)*double*0.5
  //!< using valarry to enforce safety

  double Correlation(int row, int col) const;
  //!< return correlation for variable at row and col

  double operator()(int row, int col) const;
  //!< return covariance matrix elements
  //!< this DOES NOT access elements of ary

  void Fill(const vector<double>& ntuple, double weight=1.0);
  //!< fill covariance with entry and weight
  //!< vector is used to enforce safety
  
  friend ostream& operator<<(ostream&, const Covariance&);
  //!< streaming operator for printing matrices

  friend class DataSet;

  static double fix_width;
  //!< adhoc external width for kernel smoothing

 private:
  vector<double> ary2; 
  //!< vector implementation for 2D matrix
  //!< stores sum of xy
  //!< this is NOT the same as the covariance matrix

  vector<double> ary;
  //!< stores sum of x

  bool compute_bandwidth;
  //!< used to keep track of whether bandwidth is computed
  
  int rank; //!< rank of matrix (not the same as tensor rank)
  double sum; //!< sum of weights or number of entries in covariance matrix
  double sum2; //!< sum of weights^2 or (number of entries)^2 in covariance matrix
  //!< coefficient for variance computation
  //!< computed in ComputeBandwidth

  double second_derivative_factor;
  //!< extra factor used for computing second derivative

};



#endif
