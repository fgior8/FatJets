#include "covariance.h"
#include "dataset.h"
#include "tensor.h"
#include "time.h"
#include "vector_tools.h"
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <cstring>
#include <sstream>
#include <fftw3.h>

//initialize random number generator for DataSet
MTRand_int32 DataSet::irand(time(0));
bool DataSet::b_random_gaus=false;
double DataSet::next_gaus=0;

//remove power ambiguity
int int_pow(double a, double b)
{
  return static_cast<int>(pow(a,b));
}


//constructor
DataSet::DataSet(const vector<double>& min, const vector<double>& max,
		 const vector<double>& precision)
{
  if( min.size()!= max.size() || min.size()!= precision.size())
    throw std::runtime_error
 ("Trying to initialize a DataSet with different sizes of vectors as input");
  
  //now all is good, proceed
  //Initialize(min.size(), &min[0], &max[0], &precision[0]);
  Initialize(min, max, precision);

}

void DataSet::Initialize(const vector<double>& min_, 
			 const vector<double>& max_, 
			 const vector<double>& precision_)
{
  //make sure inputs are well defined
  int nvar=min_.size();
  assert(min_.size()==max_.size());
  assert(min_.size()==precision_.size());
  
  co=Covariance(nvar);
  fix_eff_n=-1.0;
  min.resize(nvar);
  max.resize(nvar);
  one_over_width.resize(nvar);
  width.resize(nvar);
  nbins.resize(nvar);
  precision.resize(nvar);

  //debug, see if vector constructor works
  assert(min.size() == max.size());

  //size for covariance matrix and data
  int datasize=1;

  for(int i=0; i<nvar; i++)
    {
      //max must be greater than min
      assert(max_[i]>min_[i]);

      //set min and max values
      min[i]=min_[i];
      max[i]=max_[i];

      assert(precision_[i]>0);
      nbins[i]=ceil((max[i]-min[i])/precision_[i]);
      precision[i]=precision_[i];
      
      //the size of each dimension must be greater than zero
      assert(nbins[i]>0);

      datasize*=nbins[i];
      //now we may compute the width
      one_over_width[i]=static_cast<double>(nbins[i])/(max[i]-min[i]);
      width[i]=(max[i]-min[i])/static_cast<double>(nbins[i]);
    }

  //initialize the data grid
  data=Tensor<double>(nbins);
  data2=Tensor<double>(nbins);
}


void DataSet::SetN(double eff_n)
{
  fix_eff_n=eff_n;
}

void DataSet::Fill(const vector<double>& ntuple, double weight)
{
  //fill covariance matrix for bandwidth computation
  co.Fill(ntuple, weight);
  
  //increment the given bin
  vector<int> bin=GetBin(ntuple);

  data[bin]+=weight;
  data2[bin]+= (weight*weight);
  
}

Template DataSet::ComputeTemplate(double c, bool verbose_mode)
{
  //We proceed in the following way
  //1. bin the kernel and data
  //2. add padding zeroes to prevent looping
  //3. proceed to fftw routines (fourier transform)
  //4. do convolutions and what not
  //5. estimate errors via convolution
  

  ///////////////////////////////////
  /// Initialize Template
  ///////////////////////////////////

  //first compute the bandwidth matrix
  double smoothing_size=1.0;
  for(unsigned int i=0; i<width.size(); ++i)
    smoothing_size*=width[i];

  co.ComputeBandwidth(smoothing_size, 1.0, fix_eff_n);

  if(verbose_mode)
    for(unsigned int i=0; i<width.size(); ++i)
      {
	//output the width
	cout<<"width size: "<<width[i]<<endl;
	//output the smoothing size
	cout<<"smoothing size: "
	    <<1.0/sqrt(co.bandwidth_inverse[i*Rank()+i])<<endl;
      }
  
  ///////////////////////////////////////////////////
  //// Computing Kernel+padding and initialize Template
  ///////////////////////////////////////////////////
  
  //next we get a grid representing the kernel function
  //we do not allow smoothing size that is too big
  //i.e. across half of the pdf
  Tensor<double> kernel(data.length);

  //over-smoothed kernel for computing bias
  Tensor<double> kernel_oversmooth(data.length);

  //we will only need the kernel in the positive "cube"
  //i.e. all coordinates are >= 0
  //the fftw will take care of symmetric reflections
  //we also keep a list of maximum smoothing distance
  //this will be useful for adding padding zeroes
  vector<int> max_size(data.rank);
  
  //initial point
  //bin center at bin 0
  vector<double> initial_pt = GetBinCenter(vector<int>(data.rank));  

  //now proceed to construct the kernel function
  int smoothing_num=0;
  for(unsigned int i=0; i<kernel.size(); ++i)
    { 
      //get the n dimensional representation of the bin
      vector<int> index = kernel.IndexVec(i);
      
      vector<double> dist = 
	(GetBinCenter(index) - initial_pt);

      //get actual 1/bandwidth^2*distance^2
      //before exponentiating
      double dist2=co.BandwidthNorm(dist);
      
      //over-smoothed distance for computing bias
      //can be different
      double dist2_oversmooth=dist2 ;


      //don't do anything if further away than 3 sigma
      if(dist2 > 12.) continue;
      
      //find max_grid_distance
      for(unsigned int j=0; j<index.size(); ++j)
	if(index[j]>max_size[j])
	  max_size[j]=index[j];
      
      //since this bin has a non-zero smoothing
      //store the amount of smoothing
      //c is the extra factor in the silverman's rule of thumb
      kernel[i]= exp(-0.5*dist2/(c*c));
      kernel_oversmooth[i]= exp(-0.5*dist2_oversmooth/(c*c));
      ++smoothing_num;
    }

  //initialize a Template
  Template pdf;

  pdf.Initialize(data.rank, &min[0], &max[0], &precision[0]);
  Tensor<double> pdf_oversmooth(data.length);
  Tensor<double> pdf_doublesmooth(data.length);

  //initialize timer to see how much time is taken
  clock_t start=clock();
  clock_t diff;


  ////////////////////////////////////
  //// the following is strictly for performances
  ////////////////////////////////////

  
  //now we need to determine padding sizes 
  //to carry out fftw
  //large composite numbers are fastests, so things like
  //2^n, or 2^n*3^m where n,m are large are ideal
  int max_length= maximum(data.length);
  max_length=log2(max_length);
  
  vector<int> list_of_sizes;
  list_of_sizes.reserve(100);
  //generate a list of good composite numbers
  //using small primes
  //we also demand it to be even
  //for easier manipulation
  for(int i=0; i<max_length; ++i)
  for(int j=0; j<max_length; ++j)
  for(int k=0; k<max_length; ++k)
    list_of_sizes.push_back(int_pow(2,i+1)*int_pow(3,j)*int_pow(5,k));

  //get the list
  sort(list_of_sizes.begin(), list_of_sizes.end());
  
  //now we'll figure out all the sizes
  vector<int> sizes= data.length + max_size*2;
  vector<int> kernel_sizes(data.rank);

  for(unsigned int i=0; i<sizes.size(); ++i)
    for(unsigned int j=0; j<list_of_sizes.size(); ++j)
      if(list_of_sizes[j] > sizes[i])
	{
	  //we've found a good size!
	  sizes[i]=list_of_sizes[j];
	  kernel_sizes[i]=list_of_sizes[j]/2+1;
	  break;
	}

  /////// padding sizes determined
  /////////////////////////////////
  /////////////////////////////////

  //now get the padded sizes for data and kernel
  //and helper array to convert 1d index to multi-d index
  int data_padded_size=1;
  vector<int> data_padded_helper(data.rank);

  int data_complex_size=1;
  vector<int> data_complex_helper(data.rank);
  
  int kernel_padded_size=1;
  vector<int> kernel_padded_helper(data.rank);
  
  for(int i=data.rank-1; i>=0; --i)
    {
      data_padded_helper[i]=data_padded_size;
      data_complex_helper[i]=data_complex_size;
      kernel_padded_helper[i]=kernel_padded_size;
      data_padded_size*= sizes[i];

      if(i==data.rank-1)
	data_complex_size*= (sizes[i]/2+1);
      else
	data_complex_size*= sizes[i];

      kernel_padded_size*= kernel_sizes[i];
    }

  ///////////////////////////////////
  /// Initialize FFTW
  ///////////////////////////////////

  //now initialize fftw routines
  fftw_complex* out_data_complex;
  double* in_data;
  double* in_kernel, *out_kernel;

  //we need to use fftw routines
  //to allocate memeory for different arrays

  out_data_complex= (fftw_complex*) 
    fftw_malloc(sizeof(fftw_complex)* data_complex_size);

  in_data= (double*) 
    fftw_malloc(sizeof(double)* data_padded_size);

  in_kernel = (double *)
    fftw_malloc(sizeof(double)* kernel_padded_size);

  out_kernel = (double *)
    fftw_malloc(sizeof(double)* kernel_padded_size);

  //now we want to set all the arrays to zero first
  memset(out_data_complex, 0, sizeof(fftw_complex)* data_complex_size);
  memset(in_data, 0, sizeof(double)* data_padded_size);
  memset(in_kernel, 0, sizeof(double)* kernel_padded_size);
  memset(out_kernel, 0, sizeof(double)* kernel_padded_size);
  
  //initialize plans for doing fftw
  fftw_plan fftw_data, fftw_kernel, fftw_final;
  
  // plans for fourier transforming data
  fftw_data= fftw_plan_dft_r2c(data.rank, &sizes[0],
			       in_data, out_data_complex, FFTW_ESTIMATE);

  // plans for fourier transforming kernel
  // the plans is different since we are exploiting
  // specify boundary conditions 
  // i.e. even for all reflections
  vector<fftw_r2r_kind> kind(data.rank, FFTW_REDFT00);
  fftw_kernel= fftw_plan_r2r(data.rank, &kernel_sizes[0],
			     in_kernel, out_kernel,
			     &kind[0], FFTW_ESTIMATE);

  // plans to fourier transforming back the convolution
  fftw_final= fftw_plan_dft_c2r(data.rank, &sizes[0],
				out_data_complex, in_data, FFTW_ESTIMATE);



  //fftw initialization for bias
  //need to declare new variables for doubly-smoothed distribution
  double* in_data_double_kernel;
  in_data_double_kernel = (double *)
    fftw_malloc(sizeof(double)* data_padded_size);
  
  fftw_complex* out_data_double_kernel;
  out_data_double_kernel= (fftw_complex*) 
    fftw_malloc(sizeof(fftw_complex)* data_complex_size);
  
  memset(in_data_double_kernel, 0, sizeof(double)* data_padded_size);
  memset(out_data_double_kernel, 0, sizeof(fftw_complex)* data_complex_size);
  
  // plans to fourier transforming back the convolution
  fftw_plan fftw_final_double_kernel;
  fftw_final_double_kernel = fftw_plan_dft_c2r
    (data.rank, &sizes[0], out_data_double_kernel, in_data_double_kernel,
     FFTW_ESTIMATE);



  ///////////////////////////////////
  /// Initialize data input
  ///////////////////////////////////
  
  // populate data to the input to fftw include padding
  for(unsigned int i=0; i<data.size(); i++)
    {
      //get the multidimensional index
      vector<int> index= data.IndexVec(i);
      
      //populate the data input
      in_data[ (index*data_padded_helper) ] = data[i];
      //the rest is zero
    }

  // populate kernel to the input to fftw include padding
  for(unsigned int i=0; i<kernel.size(); i++)
    if(kernel[i]>0)
      {
	//get the multidimensional index
	vector<int> index= kernel.IndexVec(i);
	
	//populate the kernel input
	in_kernel[ (index*kernel_padded_helper) ] = kernel[i];
	
	//the rest is zero
      }

  //see how much time it has taken
  if(verbose_mode)
    {
      diff=clock();
      cout<<"ComputeTemplate FFTW Initialization done..."
	  <<Time(diff-start)<<" elapsed"<<endl;
    }

  //perform fourier transform for data and kernel
  fftw_execute(fftw_data);
  fftw_execute(fftw_kernel);

  //compute fourier transform of convolution
  for(unsigned int i=0; i<data_complex_size; ++i)
    {
      //since out_kernel have different dimensions
      //when comparing to out_complex
      //we need to figure out what is multiplying what
      vector<int> index_vec(data.rank);
      int i_temp=i;
      for(unsigned int j=0; j<data.rank; ++j)
	{
	  index_vec[j]=i_temp/data_complex_helper[j];
	  //now if the index is too large
	  //wrap around
	  if(index_vec[j] >= kernel_sizes[j])
	    {
	      index_vec[j] = sizes[j] - index_vec[j];
	    }

	  i_temp%=data_complex_helper[j];
	}
      
      //now we know which kernel coefficient to multiply
      int index_kernel= (index_vec*kernel_padded_helper);

      //multiply both real and complex part
      out_data_complex[i][0]*=out_kernel[index_kernel];
      out_data_complex[i][1]*=out_kernel[index_kernel];

      //keep the values in out_data_double_kernel
      out_data_double_kernel[i][0] = out_data_complex[i][0];
      out_data_double_kernel[i][1] = out_data_complex[i][1];
    }

  //convolute back
  fftw_execute(fftw_final);

  //populate the output
  for(unsigned int i=0; i<data.size(); i++)
    {
      //get the multidimensional index
      vector<int> index= data.IndexVec(i);
      
      //populate the data input
      pdf.data[i] = in_data[ (index*data_padded_helper) ];      

      //it may happen that the pdf result is negative
      //due to rounding errors
      //reset to zero
      if(pdf.data[i] < 0.)
	pdf.data[i]=0;
    }

  //see how much time it has taken
  if(verbose_mode)
    {
      diff=clock();
      cout<<"Template Convolution done..."
	  <<Time(diff-start)<<" elapsed"<<endl;
    }
    
  //now rescale the pdf
  //get a norm and compute integral
  //double sum=summation(pdf.data.ary);
  double sum=summation(pdf.data.ary);

  //get volume of the base space
  double volume=1.0;
  for(unsigned int i=0; i<width.size(); i++)
    volume*= width[i];

  if(sum <= 0)
    throw std::runtime_error
	("Template total integral is identically 0, please make sure the dataset is non-empty");
  

  //approximate integral normalization
  double scale= 1.0/(sum*volume);
  
  //now rescale the pdf and compute full errors
  //valarry does vector multiply
  pdf.data.ary*=scale;
  
  ////////////////////////////
  /// bias computation
  ////////////////////////////

  //now we want to set all the arrays to zero first
  memset(out_data_complex, 0, sizeof(fftw_complex)* data_complex_size);
  memset(in_data, 0, sizeof(double)* data_padded_size);
  memset(in_kernel, 0, sizeof(double)* kernel_padded_size);
  memset(out_kernel, 0, sizeof(double)* kernel_padded_size);
  
  //now we need to estimate the bias
  //first we obtain another smoothed dataset with a larger bandwidth
  // populate data to the input to fftw include padding
  for(unsigned int i=0; i<data.size(); i++)
    {
      //get the multidimensional index
      vector<int> index= data.IndexVec(i);
      
      //populate the data input
      in_data[ (index*data_padded_helper) ] = data[i];
      //the rest is zero
    }

  // populate kernel to the input to fftw include padding
  for(unsigned int i=0; i<kernel_oversmooth.size(); i++)
    if(kernel_oversmooth[i]>0)
      {
	//get the multidimensional index
	vector<int> index= kernel_oversmooth.IndexVec(i);
	
	//populate the kernel input
	in_kernel[ (index*kernel_padded_helper) ] = kernel_oversmooth[i];
	
	//the rest is zero
      }

  //perform fourier transform for data and kernel
  fftw_execute(fftw_data);
  fftw_execute(fftw_kernel);

  //compute fourier transform of convolution
  for(unsigned int i=0; i<data_complex_size; ++i)
    {
      //since out_kernel have different dimensions
      //when comparing to out_complex
      //we need to figure out what is multiplying what
      vector<int> index_vec(data.rank);
      int i_temp=i;
      for(unsigned int j=0; j<data.rank; ++j)
	{
	  index_vec[j]=i_temp/data_complex_helper[j];
	  //now if the index is too large
	  //wrap around
	  if(index_vec[j] >= kernel_sizes[j])
	    {
	      index_vec[j] = sizes[j] - index_vec[j];
	    }

	  i_temp%=data_complex_helper[j];
	}
      
      //now we know which kernel coefficient to multiply
      int index_kernel= (index_vec*kernel_padded_helper);

      //multiply both real and complex part
      out_data_complex[i][0]*=out_kernel[index_kernel];
      out_data_complex[i][1]*=out_kernel[index_kernel];

      //now get the doubly smoothed pdf
      out_data_double_kernel[i][0]*=out_kernel[index_kernel];
      out_data_double_kernel[i][1]*=out_kernel[index_kernel];
    }

  //convolute back
  fftw_execute(fftw_final);

  //also for the double kernel
  fftw_execute(fftw_final_double_kernel);
  
  //populate the over_smoothed output
  for(unsigned int i=0; i<pdf_oversmooth.size(); i++)
    {
      //get the multidimensional index
      vector<int> index= data.IndexVec(i);
      
      //populate the data input
      int temp_index=(index*data_padded_helper);
      pdf_oversmooth[i] = in_data[ temp_index ];      
      pdf_doublesmooth[i] = in_data_double_kernel[ temp_index ];

      //it may happen that the pdf result is negative
      //due to rounding errors
      //reset to zero
      if(pdf_oversmooth[i] < 0.)
	pdf_oversmooth[i]=0;

      if(pdf_doublesmooth[i] < 0.)
	pdf_doublesmooth[i]=0;
    }


  //scale the pdf to integrate to one
  double scale_oversmooth =
    1.0/(summation(pdf_oversmooth.ary)*volume);
  pdf_oversmooth.ary*=scale_oversmooth;

  double scale_doublesmooth = 
    1.0/(summation(pdf_doublesmooth.ary)*volume);
  pdf_doublesmooth.ary*=scale_doublesmooth;

  //bias will be obtained from the difference of the two;
  pdf.data_bias.ary= pdf_doublesmooth.ary - pdf_oversmooth.ary;

  //Finally, our pdf is computed!
  //set status to valid
  pdf.valid=true;  

  
  //clean up after ourselves
  //free up memeory used in fftw
  fftw_destroy_plan(fftw_data);
  fftw_destroy_plan(fftw_kernel);
  fftw_destroy_plan(fftw_final);
  fftw_destroy_plan(fftw_final_double_kernel);

  fftw_free(out_data_complex);
  fftw_free(in_data);
  fftw_free(in_kernel);
  fftw_free(out_kernel);
  fftw_free(in_data_double_kernel);
  fftw_free(out_data_double_kernel);

  //see how much time it has taken
  if(verbose_mode)
    {
      diff=clock();
      cout<<"Template computed ! ...total "
	  <<Time(diff-start)<<" elapsed"<<endl;
    }


  return pdf;
}

//Generate new DataSet based on fluctuation
//input = scale factor 
DataSet DataSet::GenerateDataSet(double scale_factor) const
{
  //make sure scale_factor is greater than zero
  assert(scale_factor > 0);

  //first declare a new DataSet
  //using copy constructor
  DataSet result(*this);

  double eff_n_total=0;

  //we now modify all the values in the dataset
  //the generated dataset is unweighted and events are integers,
  //i.e. data = data2
  
  //the covariance matrix on the other hand is unchanged
  
  for(unsigned int i=0; i<data.size(); ++i)
    if(data[i]>0)
    {
      double eff_n=data[i]*scale_factor;
      double random_poisson=Random_Poisson(eff_n);
      result.data[i]=random_poisson;
      result.data2[i]=random_poisson;	  
      eff_n_total+=random_poisson;
    }

  //used for normalizing covariance matrix to eff_n_total
  result.SetN(eff_n_total);

  return result;
}
