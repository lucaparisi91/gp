#ifndef TOOLSGP_H
#define TOOLSGP_H


#include <vector>
#include <cstdio>
#include <ostream>
#include <istream>
#include <iostream>
using namespace std;


template<class T,int N >
class shortVector
{
public:
  T& operator[](int i){return data[i];}
  const T& operator[](int i) const {return data[i];}
  void operator*=(T& factor)
  {
    for(int i=0;i<N;i++)
      {
	data[i]*=factor;
      }
  }

  template<class T1,int N1>
  friend ostream& operator<<(ostream& out,const shortVector<T1,N1> & vec);
  
  template<class T1,int N1>
  friend istream& operator>>(istream& in, shortVector<T1,N1> & vec);
  
  shortVector<T,N>& operator=(T& factor)
  {
    for(int i=0;i<N;i++)
      {
	data[i]=factor;
      }
  }

  
  shortVector<T,N>& operator/=(const T& factor)
  {
    for(int i=0;i<N;i++)
      {
	data[i]/=(T)(factor);
      }
    return this;
  }
  
  template<class T2>
  shortVector<T,N>& operator=(const shortVector<T2,N>& vec2)
  {
    for(int i=0;i<N;i++)
      {
	data[i]=vec2[i];
      }
    return *this;
  }
  
  template<class T2>
  shortVector<T,N>& operator/=(const shortVector<T2,N>& vec2)
  {
    for(int i=0;i<N;i++)
      {
	data[i]/=vec2[i];
      }
    return *this;
  }

  template<class T2>
  shortVector<T,N>& operator*=(const shortVector<T2,N>& vec2)
  {
    for(int i=0;i<N;i++)
      {
	data[i]*=vec2[i];
      }
    return *this;
  }
  
  template<class T2>
  shortVector<T,N>& operator*=(const T2&  factor)
  {
    for(int i=0;i<N;i++)
      {
	data[i]*=factor;
      }
    return *this;
  }
  
  template<class T2>
  shortVector<T,N>& operator+=(const shortVector<T2,N>& vec2)
  {
    for(int i=0;i<N;i++)
      {
	data[i]+=vec2[i];
      }
    return *this;
  }
  
  shortVector<T,N>& operator+=(T& factor)
  {
    for(int i=0;i<N;i++)
      {
	data[i]+=factor;
      }
    return *this;
  }
  
  
private:
  T data[N];
  
};



template<class T,class T1,int N>
shortVector<T1,N> operator/(const shortVector<T,N> & vec ,const shortVector<T1,N> & vec1)
  {
    shortVector<T1,N> vec2;
    for(int i=0;i<N;i++)
      {
	vec2[i]=vec[i]/vec1[i];
      }
    return vec2;
  }

template<class T,class T1,int N>
shortVector<T,N> operator+(const shortVector<T,N> & vec ,const shortVector<T1,N> & vec1)
  {
    shortVector<T,N> vec2;
    for(int i=0;i<N;i++)
      {
	vec2[i]=vec[i]+vec1[i];
      }
    return vec2;
  }

template<class T,class T1,int N>
shortVector<T,N> operator-(const shortVector<T,N> & vec ,const shortVector<T1,N> & vec1)
  {
    shortVector<T,N> vec2;
    for(int i=0;i<N;i++)
      {
	vec2[i]=vec[i]-vec1[i];
      }
    return vec2;
  }

template<class T,int N>
ostream& operator<<(ostream& out,const shortVector<T,N> & vec)
{
  for(int i=0;i<N;i++)
    {
      out<<vec[i]<<" ";
    }
  return out;
  
}

template<class T,int N>
istream& operator<<(istream& in,shortVector<T,N> & vec)
{
  for(int i=0;i<N;i++)
    {
      in>>vec[i]<<" ";
    }
  
  return in;
}

double  norm(double x);

double norm(const shortVector<double,1> & vec);
double norm(const shortVector<double,2> & vec);
double norm(const shortVector<double,3> & vec);

namespace indexTools
{
  
  int flatten(const shortVector<int,3> & vec,const shortVector<int,3> & shape);
  
  int flatten(const shortVector<int,1> & vec,const shortVector<int,1> & shape);
  
  shortVector<int,3> flattenInverse(int index,const shortVector<int,3> & shape);
  shortVector<int,1> flattenInverse(int index,const shortVector<int,1> & shape);

  int getRank(const shortVector<int,3> & vec);
  int getRank(const shortVector<int,1> & vec);

}


class index1D
{
public:
  typedef int shape_t;
  
  index1D(int size_){sizeIndex=size_;value=0;}
  index1D(){sizeIndex=0;value=0;}
  void resize(int size_){sizeIndex=size_;}
  inline index1D& operator()(int j){value=j;return(*this);}
  inline index1D& operator++(){value++;return (*this);}
  
  int flatten() const {return value;}
  friend ostream& operator<<(ostream& out, const index1D& index )
  {
    out<<index.value;
    return out;
  }
  
  friend istream& operator>>(istream& in, index1D& index )
  {
    in>>index.value;
    return in;
    
  }
  
  int size() const {return sizeIndex;};
  bool isEnd(){return (value>=(size()));}
  bool isValid(){return ((value>=0) and (value<size()));}
  bool operator<(int i){return value<i;}
  index1D& operator=(int i){value=i;return *this;}
  bool operator<=(int i){return i<=value;}
  index1D& operator++ (int){value++;return *this;}
  int shape(){return size();}
  
private:
  int sizeIndex;
  int value;
};

template<int N>
class indexND
{
public:
  typedef shortVector<int,N> shape_t;
  
  indexND(const shape_t & shape){value=0;reshape(shape);}
  
  void reshape(const shortVector<int,N> & shape_ ){currentShape=shape_;resize(indexTools::getRank(shape_));}
  
  inline indexND<N>& operator()(int j){value=j;return(*this);}
  
  inline indexND<N>& operator++(){value++;return (*this);}
  
  inline indexND<N>& operator()(const shortVector<int,N> & in){value=indexTools::flatten(in);return(*this);}
  
  shortVector<int,N> operator()() const {return indexTools::flattenInverse(value,currentShape);} const

  int flatten() const {return value;}
  
  friend ostream& operator<<(ostream& out, const indexND<N>& index )
  {
    out<<index.value;
    return out;
  }
  
  friend istream& operator>>(istream& in, indexND<N>& index )
  {
    in>>index.value;
    return in;
  }
  
  int size() const {return sizeIndex;};
  bool isEnd(){return (value>=(size()));}
  bool isValid(){return ((value>=0) and (value<size()));}
  bool operator<(int i){return value<i;}
  indexND<N>& operator=(int i){value=i;return *this;}
  bool operator<=(int i){return i<=value;}
  indexND<N>& operator++ (int){value++;return *this;}
  const shortVector<int,N> & shape() const {return currentShape;}
private:

  void resize(int n){sizeIndex=n;}
  int sizeIndex;
  shortVector<int,N> currentShape;
  int value;
};


template<class T,class indexType>
class tensor
{
public:
  typedef T pos_t;
  typedef indexType index_t;
  
  typedef typename index_t::shape_t shape_t;
  
  tensor(const shape_t &shape_ )
  {
    reshape(shape_);
  }
  
  tensor(const shape_t &shape_,const T& T0 )
  {
    reshape(shape_,T0);
  }
  
  void reshape(const shape_t &shape_)
  {
    currentShape=shape_;
    resize(indexTools::getRank(currentShape));
   }
  
  void reshape(const shape_t & shape_,const T& T0)
  {
    currentShape=shape_;
    resize(indexTools::getRank(currentShape),T0);
   }
  
  const shape_t & shape() const {return currentShape;}
  
  
  
  
  T& operator()(const indexType &index)
  {
    return data[index.flatten()];
  }
  
  T& operator()(int index)
  {
    return data[index];
  }
  
  const T& operator()(int index) const
  {
    return data[index];
  }
  const T& operator()(const indexType &index) const
  {
    return data[index.flatten()];
  }
    
    void resize(int size)
    {
      data.resize(size);
    };
    
    void resize(int size,const T& valueInit)
    {
      data.resize(size,valueInit);
    };

  template<class T1,class index_t1>
  friend ostream& operator<<(ostream& out, const tensor<T1,index_t1> &tensor1);

  
  T* first(){return &data[0];}
    
  size_t size() const {return data.size();}
  
  tensor& operator+=(tensor& t1)
    {
      assert(t1.size()==size());
      for(int i=0;i<t1.size();i++)
	{
	  data[i]+=t1;
	}   
    }
  
  tensor& operator/=(double t1)
    {
      
      for(int i=0;i<size();i++)
	{
	  data[i]/=t1;
	}
      return *this;
    }
  
  
private:
  vector<T> data;
  shape_t currentShape;
  
};



template<class T1,class index_t1>
 ostream& operator<<(ostream& out, const tensor<T1,index_t1> &tensor1)
  {
    for(int i=0;i<tensor1.data.size();i++)
      {
	out<<i<<" "<<tensor1.data[i]<<endl;
      }
    return out;
  }


#endif
