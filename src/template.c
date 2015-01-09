#include "template.h"
#include "tensor.h"
#include "vector_tools.h"
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <vector>

//return sign of double
inline int sign(double input)
{
  if(input>0) return 1;
  else return -1;
}

//core function to do the interpolation
double Template::operator()(const vector<double>& ntuple, double* bias) const
{
  //make sure input has the same dimension
  assert(ntuple.size()==data.rank);

  //grid coordinate = relative coordinate with respect
  //to the bin in which the extrapolation is carried out
  vector<double> grid_coordinate(data.rank);
  vector<int> grid_index(data.rank);

  //see which direction to compute derivative
  //+ or -
  vector<int> bin_direction(data.rank);

  //find the grid coordinate
  for(unsigned int i=0; i<ntuple.size(); ++i)
    {
      //compute bin_coordinate
      double temp_coordinate=(ntuple[i]-min[i])*one_over_width[i];

      //compute which bin this falls in
      int temp_bin=static_cast<int>(temp_coordinate);

      //shift grid coordinate to -0.5
      //since the first bin stores the value at bin coordinate 0.5
      //not 0
      temp_coordinate-=0.5;
      
      //minimum is -0.5, 0 is the middle of the first bin
      if(temp_coordinate < -0.5) temp_coordinate=-0.5;

      //maximum is nbins[i]-0.5, since nbins[i]-1 
      //stores values not at the maximum but middle of last bin
      else if(temp_coordinate > nbins[i]-0.5) 
	temp_coordinate=nbins[i]-0.5;

      //minimum bin is 0
      if(temp_bin<0) temp_bin=0;

      //maximum bin is nbins[i]-1
      else if(temp_bin>=nbins[i])
	temp_bin=nbins[i]-1;

      //get the bin direction
      //this is useful mainly for points near the boundary
      int temp_bin_dir=sign(temp_coordinate-temp_bin);

      //if bin is at lower boundary
      //and evaluation is done below the middle of the boundary bound
      //this will demand special treatment
      //set direction to zero for now
      //since there is no more data below
      if(temp_bin==0 && temp_bin_dir<0)
	temp_bin_dir=0;	

      //same thing for maximum
      //there is no more bin above
      else if(temp_bin==nbins[i]-1 && temp_bin_dir>0)
	temp_bin_dir=0;

      //now populate the vectors
      grid_coordinate[i]=temp_coordinate-temp_bin;
      grid_index[i]=temp_bin;
      bin_direction[i]=temp_bin_dir;
    }

  //now we know which bin to look for info and
  //what bin to look for derivative
  //fist compute the gradient
  vector<double> gradient(data.rank);

  //do the same for bias
  vector<double> gradient_bias;

  //only initialize err if it is present
  if(bias)
    gradient_bias.resize(data.rank);

  double bin_value=data[grid_index];
  double bias_value=data_bias[grid_index];

  for(int i=0; i<gradient.size(); ++i)
    {
      vector<int> current_index=grid_index;

      //special treatment for boundary
      if(bin_direction[i]==0)
	{
	  bin_direction[i]=sign(grid_coordinate[i]);
	  //simply extrapolate gradient assuming zero at boundary
	  gradient[i]=bin_direction[i]*(0.0-bin_value);

	  if(bias)
	    gradient_bias[i]=bin_direction[i]*(0.0-bias_value);
	}

      //for other cases
      //gradient is easily calculated
      else
	{
	//in the ith direction
	//we want to look in either the next
	//or previous cell
	current_index[i]+=bin_direction[i];
	
	//Tensor class takes input indices
	//and grab their values      
	gradient[i]= bin_direction[i]*(data[current_index]-bin_value);

	if(bias)
	  gradient_bias[i]= bin_direction[i]*(data_bias[current_index]-bias_value);
	}
    }

  //DEBUG
  /*
  cout<<"Template: gradient "<<gradient<<endl;
  cout<<"Template: bin value "<<bin_value<<endl;
  */

  //store error if pointer is specified
  //no need to do complicated interpolation for error
  double result=(gradient*grid_coordinate) + bin_value;

  if(bias)
    {
      //do a little linear interpolation
      (*bias)=(gradient_bias*grid_coordinate) + bias_value;
    }
  //now gradient is computed
  //we can interpolate!
  //slick use of vector! dotproduct
  return result;
}


void Template::Initialize(int nvar, const double* min_, const double* max_, 
		     const double* precision_)
{
  //make sure inputs are well defined
  assert(nvar>0);
  assert(min_!=NULL);
  assert(max_!=NULL);
  
  min.resize(nvar);
  max.resize(nvar);
  one_over_width.resize(nvar);
  width.resize(nvar);
  nbins.resize(nvar);
  
  //size for covariance matrix and data
  int datasize=1;

  for(int i=0; i<nvar; i++)
    {
      //set min and max values
      min[i]=min_[i];
      max[i]=max_[i];

      //set the number of bins
      if(precision_==NULL)
	{
	  nbins[i]=ceil((max[i]-min[i])/0.01);
	}
      else
	{
	  assert(precision_[i]>0);
	  nbins[i]=ceil((max[i]-min[i])/precision_[i]);
	}
      
      //the size of each dimension must be greater than zero
      assert(nbins[i]>0);

      datasize*=nbins[i];
      //now we may compute the width
      one_over_width[i]=static_cast<double>(nbins[i])/(max[i]-min[i]);
      width[i]=(max[i]-min[i])/static_cast<double>(nbins[i]);
    }

  //initialize the data grid
  data=Tensor<double>(nbins);
  data_bias=Tensor<double>(nbins);
}

/*
ostream& operator<<(ostream& os, const Template& template)
{
  os<<endl<<"Template info: (size="<<template.data.size()<<")"<<endl;
  for(int i=0; i<template.data.size(); ++i)
    os<<template.data[i]<<endl;

  return os;
}
*/

//code for Dresser sub class in Template
//initialize random number generator at run time
MTRand_int32 Template::Dresser::irand(time(0));

Template::Dresser::Dresser(const vector<Flag>& flag_, const Template* template_):
flag(flag_), ttemplate(template_)
{
  if(ttemplate==NULL)
    throw std::runtime_error
      ("Template is null when creating Dresser");

  if(flag.size() != ttemplate->Rank())
    throw std::runtime_error
      ("Dresser input does not match rank of Template");
  
  input_size=0;

  for(unsigned i=0; i<flag.size(); ++i)
    {
      if(flag[i]==Template::INPUT)
	++input_size;
      else
	index.push_back(i);
    }

  if(input_size == ttemplate->Rank())
    throw std::runtime_error
      ("All Dresser flags are input, there is nothing to generate!");
}

//generate event
vector<Template::DressedEvent> Template::Dresser::Generate(const vector<double>& input, 
					  int nevt) const
{

  //the number of input must match the input_size
  assert(input.size()==input_size);

  //make sure the Template Dresser refers to exists
  assert(ttemplate!=NULL);

  //initialize result
  vector<Template::DressedEvent> result;
  result.reserve(nevt);

  //initialize an input and output vector
  vector<double> base_input(ttemplate->Rank());
  vector<double> base_output(ttemplate->Rank() - input_size);

  //populate input values

  //index for the input vector
  int input_index=0;

  //index to see if current i is an input index
  vector<int>::const_iterator current_index_check=index.begin();
  for(unsigned int i=0; i<base_input.size(); ++i)
    {
      //if the index i is an input, ignore
      if((*current_index_check)==i)
	{
	  if(current_index_check!=index.end())
	    ++current_index_check;
	  continue;
	}

      //else it's an input, simply grab it
      else
	{
	  base_input[i]=input[input_index];
	  //on to the next input
	  //potential problem if input_index goes out of bound??
	  ++input_index;
	}
    }

  //first loop over however many events
  for(int i=0; i<nevt; i++)
    {

      //loop over all the output that needs Monte Carlo
      for(unsigned j=0; j<index.size(); ++j)
	{
	  base_input[index[j]]=
	    ttemplate->min[index[j]] + 
	    (ttemplate->max[index[j]]-ttemplate->min[index[j]])*irand.rand_double();
	  
	  base_output[j]=base_input[index[j]];
	}      
      //evaluate the template at that point
      double output_err;
      double output_value=(*ttemplate)(base_input, &output_err);

      //we have weight, error, time to push_back the new event!
      result.push_back(Template::DressedEvent(base_output, output_err, output_value));
    }

  //now we have to rescale all the weight so they all sum to 1
  //first get the sum of all weight
  double sum_of_weight=0;
  for(unsigned int i=0; i<result.size(); ++i)
    sum_of_weight+=result[i].rho();

  //see if the sum of weights are zero (something wrong??)
  if(sum_of_weight==0)
    cout<<"ERROR: template_values from Dresser are all zero, results may not make sense"<<endl;
  
  if(sum_of_weight>0)
    {
      //get 1/sum first to make calculation more efficient
      sum_of_weight=1.0/sum_of_weight;
      
      for(unsigned int i=0; i<result.size(); ++i)
	{
	  result[i].rho_*=sum_of_weight;
	  result[i].bias_*=sum_of_weight;
	}
    }

  //and that's it! return the list of events
  return result;
}


//generate event quickly
vector<Template::DressedEvent> Template::Dresser::GenerateFast
(const vector<double>& input, int nevt) const
{
  //the number of input must match the input_size
  assert(input.size()==input_size);

  //make sure the Template Dresser refers to exists
  assert(ttemplate!=NULL);

  //initialize result
  vector<Template::DressedEvent> result;
  result.reserve(nevt);

  //get a list of indices for MC generation
  vector<int> list_of_indices;
  vector<vector<double> > list_of_output;

  list_of_indices.reserve(nevt);
  
  //first get the initial index the input corresponds to
  int initial_index=0;

  int input_index=0;
  for(unsigned int i=0; i<flag.size(); ++i)
    if(flag[i]==INPUT)
      {
	//compute bin_coordinate
	double temp_coordinate=(input[input_index]-ttemplate->min[i])*
	  ttemplate->one_over_width[i];

	//compute which bin this falls in
	int temp_bin=ttemplate->GetBin_helper(input[input_index], i);

	//add the corrdinate index
	initial_index+= temp_bin* (ttemplate->data.index_helper[i]);
	++input_index;
      }
  list_of_indices.push_back(initial_index);
  
  //now populate all the indices
  for(int i=flag.size()-1; i>=0; --i)
    if(flag[i]==OUTPUT)
      {
	//for all indices already in the list
	//add the possible index + n*helper
	int current_size=list_of_indices.size();
	for(int k=0; k<(ttemplate->nbins[i]); ++k)
	  for(int j=0; j<current_size; ++j)
	    {
	      list_of_indices.push_back(list_of_indices[j]+ 
					k*ttemplate->data.index_helper[i]);
	    }
	
      }

  //get a list of nonzero indices
  vector<int> list_of_indices_nonzero;  

  //now grab a bunch of output vectors
  for(unsigned int i=0; i<list_of_indices.size(); ++i)
    {
      //output list does not contain the input
      vector<double> temp_output(ttemplate->Rank() - input_size);
      vector<int> bin_vec=ttemplate->data.IndexVec(list_of_indices[i]);
     
      int temp_output_index=0;
      for(unsigned j=0; j<bin_vec.size(); ++j)
	{
	  if(flag[j]==OUTPUT)
	    temp_output[temp_output_index]=
	      ttemplate->GetBinCenter(bin_vec[j],j);
	  ++temp_output_index;
	}

      //check and see if template is nonzero
      
      if(ttemplate->data[list_of_indices[i]] != 0 || 
	 ttemplate->data_bias[list_of_indices[i]] != 0)
	{
	  list_of_output.push_back(temp_output);
	  list_of_indices_nonzero.push_back(list_of_indices[i]);
	}
    }

  //now we got a full list of indices
  //time to do some MC!
  
  for(int i=0; i<nevt; ++i)
    {
      //generate an int
      int random_int = irand() % list_of_indices_nonzero.size();
      
      //push back event and that's it!
      result.push_back
	(Template::DressedEvent(list_of_output[random_int], 
		       ttemplate->data_bias[list_of_indices_nonzero[random_int]],
		       ttemplate->data[list_of_indices_nonzero[random_int]]));
    }


  //now we have to rescale all the weight so they all sum to 1
  //first get the sum of all weight
  double sum_of_weight=0;
  for(unsigned int i=0; i<result.size(); ++i)
    sum_of_weight+=result[i].rho();

  //see if the sum of weights are zero (something wrong??)
  if(sum_of_weight==0)
    cout<<"ERROR: template_values from Dresser are all zero, results may not make sense"<<endl;
  
  if(sum_of_weight>0)
    {
      //get 1/sum first to make calculation more efficient
      sum_of_weight=1.0/sum_of_weight;
      
      for(unsigned int i=0; i<result.size(); ++i)
	{
	  result[i].rho_*=sum_of_weight;
	  result[i].bias_*=sum_of_weight;
	}
    }
  
  //and that's it
  return result;
}

vector<Template::DressedEvent> Template::Dresser::GenerateFull
(const vector<double>& input) const
{
  //the number of input must match the input_size
  assert(input.size()==input_size);

  //make sure the Template Dresser refers to exists
  assert(ttemplate!=NULL);

  //initialize result
  vector<Template::DressedEvent> result;

  //get a list of indices for MC generation
  vector<int> list_of_indices;
  vector<vector<double> > list_of_output;

  
  //first get the initial index the input corresponds to
  int initial_index=0;

  int input_index=0;
  for(unsigned int i=0; i<flag.size(); ++i)
    if(flag[i]==INPUT)
      {
	//compute bin_coordinate
	double temp_coordinate=(input[input_index]-ttemplate->min[i])*
	  ttemplate->one_over_width[i];

	//compute which bin this falls in
	int temp_bin=ttemplate->GetBin_helper(input[input_index], i);

	//add the corrdinate index
	initial_index+= temp_bin* (ttemplate->data.index_helper[i]);
	++input_index;

      }
  list_of_indices.push_back(initial_index);
  
  //now populate all the indices
  for(int i=flag.size()-1; i>=0; --i)
    if(flag[i]==OUTPUT)
      {
	//for all indices already in the list
	//add the possible index + n*helper
	int current_size=list_of_indices.size();
	for(int k=0; k<(ttemplate->nbins[i]); ++k)
	  for(int j=0; j<current_size; ++j)
	    {
	      list_of_indices.push_back(list_of_indices[j]+ 
					k*ttemplate->data.index_helper[i]);
	    }

      }

  //now grab a bunch of output vectors
  for(unsigned int i=0; i<list_of_indices.size(); ++i)
    {
      //output list does not contain the input
      vector<double> temp_output(ttemplate->Rank() - input_size);
      vector<int> bin_vec=ttemplate->data.IndexVec(list_of_indices[i]);
     
      int temp_output_index=0;
      for(unsigned j=0; j<bin_vec.size(); ++j)
	{
	  if(flag[j]==OUTPUT)
	    temp_output[temp_output_index]=
	      ttemplate->GetBinCenter(bin_vec[j],j);
	  ++temp_output_index;
	}
      list_of_output.push_back(temp_output);
    }

  //now we got a full list of indices
  //since this is GenerateFull, we simply dump all the values out
  
  for(int i=0; i<list_of_output.size(); ++i)
    {     
      //push back event and that's it!
      result.push_back
	(Template::DressedEvent(list_of_output[i], 
		       ttemplate->data_bias[list_of_indices[i]],
		       ttemplate->data[list_of_indices[i]]));
    }


  //rescale the weights if desired
  /*
  
  //now we have to rescale all the weight so they all sum to 1
  //first get the sum of all weight
  double sum_of_weight=0;
  for(unsigned int i=0; i<result.size(); ++i)
    sum_of_weight+=result[i].weight;

  //see if the sum of weights are zero (something wrong??)
  if(sum_of_weight==0)
    cout<<"ERROR: template_values from Dresser are all zero, results may not make sense"<<endl;
  
  if(sum_of_weight>0)
    {
      //get 1/sum first to make calculation more efficient
      sum_of_weight=1.0/sum_of_weight;
      
      for(unsigned int i=0; i<result.size(); ++i)
	{
	  result[i].weight*=sum_of_weight;
	  result[i].err*=sum_of_weight;
	}
    }
  
  */

  //and that's it
  return result;
}



Template::Dresser Template::Generator(const vector<Flag>& flag_) const
{
  return Template::Dresser(flag_, this);
}
