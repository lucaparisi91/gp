#ifndef FOURIER_H
#define FOURIER_H

#include <complex>
#include "fftw3.h"
#include "traits.h"
#include "tools.h"

int shiftNegativeFreq(int index,int size);

shortVector<int, 3> shiftNegativeFreqIndex( const shortVector<int, 3> &vec,const shortVector<int,3> & shape);

shortVector<int, 1> shiftNegativeFreqIndex( const shortVector<int, 1> vec,const shortVector<int,1> & shape);
template<class dimensions_t,class T>
class fourierStrategy
{
  
};

template<>
class fourierStrategy<D1_t,complex<double> >
{
public:
  
  void setPlan(fftw_plan& plan, complex<double>* in,complex<double> *out,int N,int direction)
  {
    plan=fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(in), reinterpret_cast<fftw_complex*>(out), direction, FFTW_ESTIMATE);
  }

  
  template<class tensor_t>
  void getFrequencies(tensor_t& frequencies,const typename tensor_t::pos_t & boxDimensions)
  {
    
    typename tensor_t::index_t index(frequencies.shape());
    
    for (index=0;index<index.size();index++)
      {
	frequencies(index)=shiftNegativeFreqIndex(index(),index.shape())/(boxDimensions);
	frequencies(index)*=2*M_PI;
	
      }
  }
  
};

#endif
