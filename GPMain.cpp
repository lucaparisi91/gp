#include "GP.h"
#include "./tools.h"
#include <vector>
#include <iostream>
#include <cmath>
#include "fftw3.h"
#include <complex>

using namespace std;

int main(int argc,char** argv)
{
  shortVector<int,1> bins;
  shortVector<double,1> max;
  shortVector<double,1> min;
  int N=100;
  double ratio=0.6;
  double g=2;
  double g_tilde=g*ratio;
  max[0]=40.;
  min[0]=-40.;
  bins[0]=1000;
  
  GPSolver<gPEquationLY1D> solver(bins,min,max);
  solver.setParameter("deltag",(g-g_tilde)*N/2.);
  
  double c=pow( (pow( (g+g_tilde)/2.,3/2.) + pow((g-g_tilde)/2.,3/2.)),2/3.);
  //double c=g*pow(3./2.,2./3);
  //solver.setParameter("omega",1e-5);
  solver.setParameter("g",c*pow(N/2.,1./3));
  solver.setTimeStep(0.01);
  solver.run(1000,10);
  
  
}


