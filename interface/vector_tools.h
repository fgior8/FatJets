#ifndef __VECTORTOOLS__
#define __VECTORTOOLS__

#include <vector> 
#include <cmath>

using namespace std;

template<class T>
T operator*(const vector<T>& a, const vector<T>& b)
{
  if(a.size() == 0)
    return T();

  assert(a.size() == b.size());
  T result=a[0]*b[0];

  for(unsigned int i=1; i<a.size(); ++i)
    result+= a[i]*b[i];

  return result;
}

template<class T>
T maximum(const vector<T>& a)
{
  if(a.size() == 0)
    return T();
  T result=a[0];

  for(unsigned int i=1; i<a.size(); ++i)
    if(a[i] > result)
      result=a[i];

  return result;
}

template<class T>
T summation(const vector<T>& a)
{
  if(a.size()==0)
    return T();
  T result=a[0];

  for(unsigned int i=1; i<a.size(); ++i)
    result+=a[i];

  return result;
}

template<class T>
T average(const vector<T>& a)
{
  if(a.size()==0)
    return T();
  T result=a[0];

  for(unsigned int i=1; i<a.size(); ++i)
    result+= ( (a[i]-result)/((T) i+1) );

  return result;
}

template<class T>
vector<T> abs(const vector<T>& a)
{
  vector<T> result(a.size());

  for(unsigned int i=0; i<a.size(); ++i)
    result[i]=abs(a[i]);

  return result;
}

template<class T>
vector<T> sqrt(const vector<T>& a)
{
  vector<T> result(a.size());

  for(unsigned int i=0; i<a.size(); ++i)
    result[i]=sqrt(a[i]);

  return result;
}

template<class T>
vector<T> prod(const vector<T>& a, const vector<T>& b)
{
  assert(a.size() == b.size());
  vector<T> result(a.size());

  for(unsigned int i=0; i<a.size(); ++i)
    result[i]=a[i]*b[i];

  return result;
}

template<class T>
vector<T> operator-(const vector<T>& a, const vector<T>& b)
{
  assert(a.size() == b.size());
  vector<T> result(a.size());

  for(unsigned int i=0; i<a.size(); ++i)
    result[i]=a[i]-b[i];

  return result;
}

template<class T>
vector<T> operator+(const vector<T>& a, const vector<T>& b)
{
  assert(a.size() == b.size());
  vector<T> result(a.size());

  for(unsigned int i=0; i<a.size(); ++i)
    result[i]=a[i]+b[i];

  return result;
}


template<class T>
vector<T> operator/(const vector<T>& a, const T& b)
{
  vector<T> result=a;

  for(unsigned int i=0; i<a.size(); ++i)
    result/=b;

  return result;
}

template<class T>
vector<T> operator*(const vector<T>& a, const T& b)
{
  vector<T> result=a;

  for(unsigned int i=0; i<a.size(); ++i)
    result*=b;

  return result;
}

template<class T>
void operator/=(vector<T>& a, const T& b)
{
  for(unsigned int i=0; i<a.size(); ++i)
    a[i]/=b;
}

template<class T>
void operator*=(vector<T>& a, const T& b)
{
  for(unsigned int i=0; i<a.size(); ++i)
    a[i]*=b;
}


#endif
