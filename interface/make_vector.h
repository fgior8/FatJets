#ifndef __MAKEVECTOR__
#define __MAKEVECTOR__

#include <vector>

using namespace std;

template<class T>
vector<T> make_vector()
{
  vector<T> result(0);
  return result;
}

template<class T>
vector<T> make_vector(const T& a)
{
  vector<T> result(1);
  result[0]=a;
  return result;
}

template<class T>
vector<T> make_vector(const T& a1, const T& a2)
{
  vector<T> result(2);
  result[0]=a1;
  result[1]=a2;
  return result;
}

template<class T>
vector<T> make_vector(const T& a1, const T& a2, const T& a3)
{
  vector<T> result(3);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  return result;
}

template<class T>
vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4)
{
  vector<T> result(4);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  return result;
}

template<class T>
vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4, const T& a5)
{
  vector<T> result(5);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  result[4]=a5;
  return result;
}

template<class T>
vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4, const T& a5, const T& a6)
{
  vector<T> result(6);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  result[4]=a5;
  result[5]=a6;
  return result;
}

template<class T>
  vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4, const T& a5, const T& a6, const T& a7)
{
  vector<T> result(7);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  result[4]=a5;
  result[5]=a6;
  result[6]=a7;
  return result;
}

template<class T>
  vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4, const T& a5, const T& a6, const T& a7, const T& a8)
{
  vector<T> result(8);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  result[4]=a5;
  result[5]=a6;
  result[6]=a7;
  result[7]=a8;
  return result;
}

template<class T>
  vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4, const T& a5, const T& a6, const T& a7, const T& a8, const T& a9)
{
  vector<T> result(9);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  result[4]=a5;
  result[5]=a6;
  result[6]=a7;
  result[7]=a8;
  result[8]=a9;
  return result;
}


template<class T>
  vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4, const T& a5, const T& a6, const T& a7, 
			const T& a8, const T& a9, const T& a10)
{
  vector<T> result(10);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  result[4]=a5;
  result[5]=a6;
  result[6]=a7;
  result[7]=a8;
  result[8]=a9;
  result[9]=a10;
  return result;
}

template<class T>
  vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4, const T& a5, const T& a6, const T& a7, 
			const T& a8, const T& a9, const T& a10, const T& a11)
{
  vector<T> result(11);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  result[4]=a5;
  result[5]=a6;
  result[6]=a7;
  result[7]=a8;
  result[8]=a9;
  result[9]=a10;
  result[10]=a11;
  return result;
}

template<class T>
  vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4, const T& a5, const T& a6, const T& a7, 
			const T& a8, const T& a9, const T& a10, const T& a11, const T& a12)
{
  vector<T> result(12);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  result[4]=a5;
  result[5]=a6;
  result[6]=a7;
  result[7]=a8;
  result[8]=a9;
  result[9]=a10;
  result[10]=a11;
  result[11]=a12;
  return result;
}

template<class T>
  vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4, const T& a5, const T& a6, const T& a7, 
			const T& a8, const T& a9, const T& a10, const T& a11, const T& a12, const T& a13)
{
  vector<T> result(13);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  result[4]=a5;
  result[5]=a6;
  result[6]=a7;
  result[7]=a8;
  result[8]=a9;
  result[9]=a10;
  result[10]=a11;
  result[11]=a12;
  result[12]=a13;
  return result;
}

template<class T>
  vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4, const T& a5, const T& a6, const T& a7, 
			const T& a8, const T& a9, const T& a10, const T& a11, const T& a12, const T& a13, const T& a14)
{
  vector<T> result(14);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  result[4]=a5;
  result[5]=a6;
  result[6]=a7;
  result[7]=a8;
  result[8]=a9;
  result[9]=a10;
  result[10]=a11;
  result[11]=a12;
  result[12]=a13;
  result[13]=a14;
  return result;
}

template<class T>
  vector<T> make_vector(const T& a1, const T& a2, const T& a3, const T& a4, const T& a5, const T& a6, const T& a7, 
			const T& a8, const T& a9, const T& a10, const T& a11, const T& a12, const T& a13, const T& a14, const T& a15)
{
  vector<T> result(15);
  result[0]=a1;
  result[1]=a2;
  result[2]=a3;
  result[3]=a4;
  result[4]=a5;
  result[5]=a6;
  result[6]=a7;
  result[7]=a8;
  result[8]=a9;
  result[9]=a10;
  result[10]=a11;
  result[11]=a12;
  result[12]=a13;
  result[13]=a14;
  result[14]=a15;
  return result;
}


//convert vectors of one kind to another
template<class A, class B>
  vector<B> convert(const vector<A> a)
{
  vector<B> result(a.size());
  for(unsigned int i=0; i<a.size(); ++i)
    {
      result[i]=a[i];
    }
  return result;
}

#endif
