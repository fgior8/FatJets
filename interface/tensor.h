#ifndef __TENSOR__
#define __TENSOR__

#include <vector>
#include <iostream>
#include <cassert>
#include "vector_tools.h"

using namespace std;

//streamer method for vector
template < typename T >
inline ostream& operator<<(ostream& os, const vector<T>& v)
{
  os<<"vector: ";
  for(int i=0; i<v.size(); i++)
    os<<v[i]<<",";

  return os;
}

/*!
  \brief class for manipulating tensors


  Implementation of multi-dimensional tensor via vector,
  need not be squre. This is used as a multi-dimensional histogram
*/

template < typename T >
class Tensor
{
 public:

  Tensor();
  //!< default constructor
  
  /*
  Tensor(int rank_, const int* length_);
  //!< constructor
  */

  Tensor(const vector<int>& length_);
  //!< safer constructor
  
  int Length(int i) const
  {
    assert(i>=0);
    assert(i<rank);
    return length[i];
  }
  //!< access length of a particular dimension

  int Rank() const
  {return rank;}
  //!< return rank of tensor

  T& operator[](const vector<int>& ntuple);
  //!< quick access elements of the tensor
  //!< might return wrong elements if out of bound

  T operator[](const vector<int>& ntuple) const;
  //!< quick const access elements of the tensor
  //!< might return wrong elements if out of bound

  T& operator[](int);
  //!< quick access elements of the tensor
  //!< might return wrong elements if out of bound

  T operator[](int) const;
  //!< quick const access elements of the tensor
  //!< might return wrong elements if out of bound

  unsigned int size() const;
  //!< return size of array

  vector<int> IndexVec(int index) const;
  //!< return a vector from an index

  int Index(const vector<int>& index) const;
  //!< return an index from a vector

  bool ValidIndex(int) const;
  //!< see if index is within range

  int GridDistance(int index1, int index2) const;
  //!< a function to figure out the "grid" norm between
  //!< two index, i.e. max(abs(xi-y-)) between two vector of indices
  //!< function is unsafe (assuming index1 and index2) both make sense


  friend class DataSet;
  friend class Template;

 private:
  vector<T> ary;
  //!< vector to store data

  int rank; //!< rank/dimension of tensor
  vector<int> length; //!< array to store the length of each dimension
  vector<int> index_helper; 
  //!< array to help get tensor element
  //!< its elements are
  //!< [lengthn*length(n-1)*...*length1, length(n-1)*...*length1, ..., length1]
  //!< to grab an object in tensor with index=[n1, n2, n3..., nn]
  //!< simply take the index dotproduct index_helper

  void initialize_index_helper();
  //!< function to initialize index_helper
  //!< used in access operator[]

};


template < typename T >
Tensor<T>::Tensor() : rank(0)
{
}

/*
template < typename T >
Tensor<T>::Tensor(int rank_, const int* length_) : rank(rank_)
{
  int size=1;
  for(int i=0; i<rank; i++)
    {
      //length must be nonzero integer
      assert(length_[i]>0);
      size*=length_[i];
    }

  ary.resize(size);
  length=vector<int>(length_, rank_);

  //now initialize index_helper
  initialize_index_herlper();
}
*/

template < typename T >
Tensor<T>::Tensor(const vector<int>& length_) : rank(length_.size())
{
  int size=1;
  for(int i=0; i<rank; i++)
    {
      //length must be nonzero integer
      assert(length_[i]>0);
      size*=length_[i];
    }

  ary.resize(size);
  length.resize(length_.size());
  length=length_;

  //now initialize index_helper
  initialize_index_helper();
			       
}


template < typename T >
T& Tensor<T>::operator[](const vector<int>& ntuple)
{
  //make sure index is not out of bound
  for(int i=0; i<rank; i++)
    {
      /*
      cout<<"ntuple: "<<ntuple[i]<<endl;
      cout<<"length: "<<length[i]<<endl;
      cout<<"index_helper: "<<index_helper[i]<<endl;
      */

      assert(ntuple[i]<length[i]);
    }

  /*
  cout<<"computed size: "<<ntuple*index_helper<<endl;
  cout<<"max size: "<<ary.size()<<endl;
  */

  return ary[ntuple*index_helper];
}

template < typename T >
T Tensor<T>::operator[](const vector<int>& ntuple) const
{
  //make sure index is not out of bound
  for(int i=0; i<rank; i++)
    assert(ntuple[i]<length[i]);
  
  return ary[ntuple*index_helper];
}

template < typename T >
T& Tensor<T>::operator[](int index)
{
  //make sure index is not out of bound
  assert(ValidIndex(index));
  
  return ary[index];
}

template < typename T >
T Tensor<T>::operator[](int index) const
{
  //make sure index is not out of bound
  assert(ValidIndex(index));
  
  return ary[index];
}

template < typename T >
unsigned int Tensor<T>::size() const
{
  return ary.size();
}


template < typename T >
vector<int> Tensor<T>::IndexVec(int index) const
{
  //for a given index
  //we want to write
  //index= n1*(index_helper[0]) + n2*(index_helper[1]) + ...
  //and extract the resulting vector of ni

  //define result
  vector<int> result(rank);

  //grab ni
  for(unsigned int i=0; i<index_helper.size(); i++)
    {
      //extract ni
      result[i]=index/index_helper[i];

      //substract away ni*index_helper[i]
      //using module operator
      index%=index_helper[i];
    }
  return result;
}

template < typename T >
int Tensor<T>::Index(const vector<int>& index) const
{
  return index*index_helper;
}

template < typename T >
bool Tensor<T>::ValidIndex(int index) const
{
  return (index>=0 && index< ary.size());
}

template < typename T >
int Tensor<T>::GridDistance(int index1, int index2) const
{
  //assume index1,2 are both written in the form
  // n1*(index_helper[0]) + n2*(index_helper[1]) + ...
  // same for
  // m1*(index_helper[0]) + m2*(index_helper[1]) + ...
  // we wnat to find max(ni-mi)

  return maximum(abs(IndexVec(index1)-IndexVec(index2)));
}

template < typename T >
void Tensor<T>::initialize_index_helper()
{
  //the elements of index_helpers are
  //[lengthn*length(n-1)*...*length1, length(n-1)*...*length1, ..., length1]
  index_helper.resize(rank);
  
  int size=1;
  //iterate backward
  for(int i=rank-1; i>=0; --i)
    {
      index_helper[i]=size;
      size*=length[i];
    }
}



#endif
