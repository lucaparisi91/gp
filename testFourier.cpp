
#include <complex>
#include <vector>
#include "fftw3.h"
#include <iostream>
using namespace std;

int main(int argc,char ** argv)
{
  int size;
  size=1000;
  
  complex<double> *x=new complex<double>[size];
  complex<double> *y=new complex<double>[size];
  
  double pos;
  
  fftw_plan plan;
  plan=fftw_plan_dft_1d(size, reinterpret_cast<fftw_complex*>(&x[0]), reinterpret_cast<fftw_complex*>(&y[0]), +1, FFTW_ESTIMATE);

  for(int i=0;i<size;i++)
    {
      pos=i/1000.*20-10.;
      x[i]=exp(-abs(pos*pos));
      y[i]=0;
    }
  
  
  fftw_execute(plan);

  for(int i=0;i<size;i++)
    {
      pos=i/1000.*20-10.;
      cout << pos <<" "<< abs(y[i].real())<< " "<< y[i].imag()<<endl;
    }
  
}
