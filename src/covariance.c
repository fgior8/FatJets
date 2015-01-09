#include "covariance.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <cstdio>
#include "matrix_inverse.h"
#include "vector_tools.h"

double Covariance::fix_width=-1.0;
//disable to fix bandwidth by default

//invert matrix and return determinant
double Invert_Matrix(double* matrix, int rank)
{
  MatrixClass<double> M(rank, rank);
  
  //dump matrix entries to external class
  for(int i=0; i<rank; i++)
  for(int j=0; j<rank; j++)
    M.setvalue(i,j,matrix[i*rank+j]);
  
  //invert matrix
  double matrix_determinant=M.invert();
  assert(matrix_determinant!=0);
  
  bool state;
  //dump values back to matrix
  for(int i=0; i<rank; i++)
  for(int j=0; j<rank; j++)
    {
      M.getvalue(i,j,matrix[i*rank+j], state);
      //make sure there is no error
    }

  return matrix_determinant;
}



using namespace std;

Covariance::Covariance(int rank_): rank(rank_),
  sum(0), sum2(0), second_derivative_factor(0)
{
  //demand nonzero rank
  assert(rank_>0);

  //initialize matrix
  ary2.resize(rank*rank);
  bandwidth_inverse.resize(rank*rank);
  bandwidth.resize(rank*rank);
  ary.resize(rank);
  compute_bandwidth=false;
}

//default value = -1, then custom_coefficient not used
void Covariance::ComputeBandwidth(double min_size, double extra_factor, 
				  double fix_eff_n)
{
  assert(sum2>0);
  //dump values to a double*
  vector<double> matrix(rank*rank, 0);
  
  double eff_n=sum*sum/sum2;

  //if n is fixed, we will fix eff_n instead
  if(fix_eff_n > 0)
    {
      /*
      cout<<"Fixing eff_n "<<fix_eff_n<<" in bandwidth computation"<<endl;
      cout<<"Compared to original: "<<eff_n<<endl;
      cout<<"Extra factor: "<<extra_factor<<endl;
      */
      eff_n=fix_eff_n;
    }

  else
    {
      eff_n*=extra_factor;
    }

  //Silverman's Rule of Thumb
  //H = n^(-1/(d+4))*sqrt(covariance)
  //H = sqrt(covariance/coefficient)
  //we are computing H^-2
  double coefficient=pow(eff_n, 2.0/(rank+4));
  
  second_derivative_factor=pow(eff_n, -4.0/((rank+4)*(rank+6)));
  
  double one_over_coefficient=1.0/coefficient;

  //fill matrix with covariance
  for(int i=0; i<rank; i++)
  for(int j=i; j<rank; j++)
    {
      //(*this)(i,j) computes elements in covariance matrix
      matrix[i*rank+j]=(*this)(i, j);
      //enforce symmetric matrix
      matrix[j*rank+i]=matrix[i*rank+j];
      
    }

  //fill the bandwidth matrix
  //for storage, the actual computation uses bandwidth_inverse
  for(int i=0; i<rank; i++)
  for(int j=i; j<rank; j++)
    {
      bandwidth[i*rank+j]= one_over_coefficient*matrix[i*rank+j];
      //enforce symmetry
      bandwidth[j*rank+i]= bandwidth[i*rank+j];
      //make sure bandwidth_inverse is finite
      assert(isfinite(bandwidth[j*rank+i]));
    }

  //now invert matrix
  //and get determinant
  double determinant=Invert_Matrix(&matrix[0], rank);
  
  //matrix should be positive definite
  //since our matrix really is covariance^2
  determinant=sqrt(determinant);

  //now check and see if the current smoothing size is too small
  //smoothing_size = det(bandwidth)
  //bandwidth_inverse = sqrt(matrix/coefficient)
  double smoothing_size=determinant/pow(coefficient, rank/2.0);

  //if smoothing area is too small
  //rescale coefficient
  if(smoothing_size < min_size)
    {
      cout<<"WARNING: approaching minimum bandwidth size, "<<endl
	  <<"consider decreasing precision"<<endl;

      cout<<"Kernel smoothing size: "<<smoothing_size<<endl;
      cout<<"Bin size: "<<min_size<<endl;
      double coefficient_new=pow(determinant/min_size, 2.0/rank);
      
      //rescale the bandwidth matrix
      bandwidth*= (coefficient/coefficient_new);

      //use the new coefficient
      coefficient=coefficient_new;
    }

  
  //populate bandwidth_inverse
  //use for calculations
  //i.e. ComputeBandwidthNorm
  for(int i=0; i<rank; i++)
  for(int j=i; j<rank; j++)
    {
      bandwidth_inverse[i*rank+j]= coefficient*matrix[i*rank+j];
      //enforce symmetry
      bandwidth_inverse[j*rank+i]= bandwidth_inverse[i*rank+j];
      //make sure bandwidth_inverse is finite
      assert(isfinite(bandwidth_inverse[j*rank+i]));
    }

  //if fix_width > 0, fix kernel to be identity * fix_width
  if(fix_width > 0)
    {
      for(int i=0; i<rank; i++)
	for(int j=i; j<rank; j++)
	  {
	    bandwidth_inverse[i*rank+j]= 1.0/(fix_width*fix_width);
	    //enforce symmetry
	    bandwidth_inverse[j*rank+i]= bandwidth_inverse[i*rank+j];
	  }
    }


  
  //output the bandwidth matrix^2
  //this should be small!
  //DEBUG

  
  //cout<<"Bandwidth Matrix should be reasonably small): "<<endl;

  /*
  for(int i=0; i<rank; i++)
    {
      for(int j=0; j<rank; j++)
	{
	  cout<<bandwidth[i*rank+j]<<" ";
	}
      cout<<endl;
    }
  */

  //turn on flag to indicate bandwidth is computed
  compute_bandwidth=true;
}

double Covariance::BandwidthNorm(const vector<double>& vec, int power)
{
  assert(power==1 || power==-1);

  //see if bandwidth is already computed
  if(!compute_bandwidth)
    ComputeBandwidth();

  //make sure vector length is same as rank
  assert(vec.size()==rank);
  
  double result=0;
 
  //if negative power
  if(power==-1)
    for(int i=0; i<rank; i++)
    for(int j=0; j<rank; j++)
      result+= vec[i]*bandwidth_inverse[i*rank+j]*vec[j];

  else if(power==1)
    for(int i=0; i<rank; i++)
    for(int j=0; j<rank; j++)
      result+= vec[i]*bandwidth[i*rank+j]*vec[j];

  return result;
}

double Covariance::Correlation(int r, int c) const
{
  //return covariance(r,c)/(rms(r)*rms(c))
  double covarx=(*this)(r,r);
  double covary=(*this)(c,c);
  assert(covarx>0);
  assert(covary>0);

  return (*this)(r, c)/sqrt(covarx*covary);
}

double Covariance::operator()(int row, int col) const
{
  //demand indices to be within range
  assert(row>=0 && row<rank);
  assert(col>=0 && col<rank);
  //make sure entry isn't zero
  assert(sum2>0);

  //effective number of entries
  double eff_n=sum*sum/sum2;

  //apply N/(N-1) to get an unbiased estimate
  //if eff_n < 2, just ignore such factor
  double unbias_factor=1.0;
  if(eff_n>2) unbias_factor=eff_n/(eff_n-1);

  double meanxy=ary2[row*rank+col]/sum;
  double meanx=ary[row]/sum;
  double meany=ary[col]/sum; 
  
  return unbias_factor*(meanxy-meanx*meany);
}

void Covariance::Fill(const vector<double>& ntuple, double weight)
{
  //make sure ntuple is of right size
  assert(ntuple.size()==rank);
  //bandwidth needs to be changed
  compute_bandwidth=false;
  //increment counts
  sum2+=weight*weight;
  sum+=weight;
  
  //loop over ary2 and ary and increment x*y and x
  for(int i=0; i<rank; i++)
    {
      ary[i]+=(ntuple[i]*weight);
      for(int j=i; j<rank; j++)
	{
	  ary2[i*rank+j]+= (ntuple[i]*ntuple[j]*weight);
	  //symmetrize results
	  ary2[j*rank+i]=ary2[i*rank+j];
	}
    }
}


ostream& operator<<(ostream& os, const Covariance& co)
{
  //cout covariance matrix
  os<<endl;
  for(int i=0; i<co.rank; i++)
    {
      os<<"[ ";
	for(int j=0; j<co.rank; j++)
	  {
	    os<<setw(5)<<setprecision(3)<<co(i,j)<<" ";
	  }
      os<<"]"<<endl;
    }

  os<<endl;

  //cout bandwidth * covariance
  //to check and see if it is proportional to unity
  /*
  for(int i=0; i<co.rank; i++)
    {
      os<<"[ ";
	for(int j=0; j<co.rank; j++)
	  {
	    double value=0;
	    for(int r=0; r<co.rank; r++)
	      {
		value+=co(i,r)*co.bandwidth[r*co.rank+j];
	      }
	    os<<setw(5)<<setprecision(3)<<value<<" ";

	  }
      os<<"]"<<endl;
    }
  */
  return os;
}
