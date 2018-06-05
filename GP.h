#ifndef GP_H
#define GP_H

#include "tools.h"
#include <complex>
#include "traits.h"
#include <iostream>
#include "fourier.h"
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;
template<class T,class T1>
class waveGeneral
{
public:
  typedef T pos_t;
  typedef T1 index_t;
  typedef typename index_t::shape_t shape_t;
  typedef tensor<complex<double> ,index_t> values_t;
  typedef tensor<pos_t ,index_t> positions_t;
  typedef complex<double> value_t;
  
  waveGeneral(shape_t &bins) : positions(bins),waveValue(bins,0),index(bins)
  {
    currentShape=positions.shape();
    
  };
  
  inline size_t size() const {return positions.size();}
  tensor<pos_t,index_t>& grid(){return positions;};
  tensor<complex<double>,index_t>& value(){return waveValue;};
  const tensor<complex<double> ,index_t>& value() const {return waveValue;};
  const tensor<pos_t,index_t>& grid() const {return positions;};
  
  

  void setMin(const pos_t &min_){min=min_;}
  void setMax(const pos_t &max_){max=max_;}
  
  pos_t getMin() const {return min;}
  pos_t getMax() const {return max;}
  const shape_t & shape() const {return currentShape;}
  pos_t getLength() const {return (max-min);}

  double getArea()
  {
    pos_t length;
    length=getLength();
    length/=shape();
    return area(length);
  }

  double sumNorms()
  {
    double waveNorm;
    waveNorm=0;
    for(index=0;index<index.size();index++)
      {
	waveNorm+=(std::norm(waveValue(index)));
      }
    return waveNorm;
  }
  
  double norm()
  {
    return sumNorms()*getArea();
  }
  
  void normalize()
  {
    double normWave;
    normWave=norm();
    for(int i=0;i<size();i++)
      {
	value()(i)/=sqrt(normWave);
      }
  }
  
private:
  tensor<complex<double> ,index_t> waveValue;
  tensor<pos_t,index_t> positions;
  pos_t min;
  pos_t max;
  shape_t currentShape;
  index_t index;
  template<class OT,class OT1>
  friend ostream& operator<<(ostream& out, const waveGeneral<OT,OT1>& index );
  
};

double area(shortVector<double,1> & vec){return vec[0];}
double area(shortVector<double,3> & vec){return vec[0]*vec[1]*vec[2];}

class wave1D : public waveGeneral<shortVector<double,1>,indexND<1> >
{
public:

  wave1D(shape_t &bins) : waveGeneral<shortVector<double,1>,indexND<1> >(bins){};
  
public:
  typedef D1_t dimensions_t;
  
  
};


class harmonicPotential1D
{
public:
  typedef wave1D::index_t index_t;
  typedef wave1D::values_t values_t;
  typedef wave1D::positions_t positions_t;
  
  harmonicPotential1D(){omega=1;g=0;}
  void setParameter(string name,double parameter)
  {
    if (name == "omega")
      {
	omega=parameter;
      }
    else if (name == "g")
      {
	g=parameter;
      }
  }

  
  double getParameter(string name) const
  {
    if (name == "omega")
      {
	return omega;
      }
    else if (name == "g")
      {
	return g;
      }
    
  }
  
  
  void evaluate(const values_t  &valuesIn,const positions_t &positions,values_t &valuesOut)
  {
    for(int i=0;i<valuesIn.size();i++)
      {
	valuesOut(i)=0.5*omega*norm(positions(i)) + g*norm(valuesIn(i));
      }
  }
  
 private:
  double omega;
  double g;
};


class LYPotential
{
public:
  typedef wave1D::index_t index_t;
  typedef wave1D::values_t values_t;
  typedef wave1D::positions_t positions_t;
  
  LYPotential(){deltag=0;g=1;}
  
  void evaluate(const values_t  &valuesIn,const positions_t &positions,values_t &valuesOut)
  {
    double C=sqrt(2)/M_PI*pow(g,3/2.);
    double tmp=0;
    for(int i=0;i<valuesIn.size();i++)
      {
	tmp=norm(valuesIn(i));
	
	valuesOut(i)=deltag*tmp -C*sqrt(tmp);
      }
  }

  void setParameter(string name,double parameter)
  {
    if (name == "deltag")
      {
	deltag=parameter;
      }
    else if (name == "g")
      {
	g=parameter;
      }
    else
      {
	cout << "Unkown parameter"<<endl;
	exit(2);
      }
  }
  
  double getParameter(string name) const
  {
    if (name == "deltag")
      {
	return deltag;
      }
    else if (name == "g")
      {
	return g;
      }
    else
      {
	cout << "Unkown parameter"<<endl;
	exit(2);
      }
  }

  
 private:
  
  double g;
  
  double deltag;
  
};

class chemicalPotentialModelFitted
{
public:
  typedef wave1D::index_t index_t;
  typedef wave1D::values_t values_t;
  typedef wave1D::positions_t positions_t;
  
  chemicalPotentialModelFitted()
  {
    deltag=0;g=1;
    a=-0.232606;
    b=-0.208116;
    c=0.175507;
    d=0.658604;
    e=0.38154;
    
  }

  void evaluate(const values_t  &valuesIn,const positions_t &positions,values_t &valuesOut)
  {
    double C=sqrt(2)/M_PI*pow(g,3/2.);
    double x=0;
    for(int i=0;i<valuesIn.size();i++)
      {
	x=norm(valuesIn(i))*100;	
	valuesOut(i)=a + 2*b*x+2*c*x*atan(e*pow(x,d)) + c*e*d*pow(x,d+1)/(1+e*e*pow(x,2*d));
	
      }
    
    
  }
  
  void setParameter(string name,double parameter)
  {
    if (name == "deltag")
      {
	deltag=parameter;
      }
    else if (name == "g")
      {
	g=parameter;
      }
    else
      {
	cout << "Unkown parameter"<<endl;
	exit(2);
      }
  }
  
  double getParameter(string name) const
  {
    if (name == "deltag")
      {
	return deltag;
      }
    else if (name == "g")
      {
	return g;
      }
    else
      {
	cout << "Unkown parameter"<<endl;
	exit(2);
      }
  }

  
 private:
  
  double g;
  
  double deltag;

  //

  double a;
  double b;
  double c;
  double d;
  double e;
  
};





class LYPotentialHarmonicTrapping
{
public:
  typedef wave1D::index_t index_t;
  typedef wave1D::values_t values_t;
  typedef wave1D::positions_t positions_t;
  
  LYPotentialHarmonicTrapping(){deltag=0;g=1;omega=1;}
  
  void evaluate(const values_t  &valuesIn,const positions_t &positions,values_t &valuesOut)
  {
    double C=sqrt(2)/M_PI*pow(g,3/2.);
    double tmp=0;
    for(int i=0;i<valuesIn.size();i++)
      {
	tmp=norm(valuesIn(i));
	
	valuesOut(i)=deltag*tmp -C*sqrt(tmp) + 0.5*omega*norm(positions(i));
      }
  }

  void setParameter(string name,double parameter)
  {
    if (name == "deltag")
      {
	deltag=parameter;
      }
    else if (name == "g")
      {
	g=parameter;
      }
    else if (name=="omega")
      {
	omega=parameter;
      }
    else
      {
	cout << "Unkown parameter"<<endl;
	exit(2);
      }
  }
  
  double getParameter(string name) const
  {
    if (name == "deltag")
      {
	return deltag;
      }
    else if (name == "g")
      {
	return g;
      }
    else if (name=="omega")
      {
	return omega;
      }
    else
      {
	cout << "Unkown parameter"<<endl;
	exit(2);
      }
  }

  
 private:
  
  double g;
  double omega;
  double deltag;
  
};


template<class wave_t,class potential_t>
class splitFourierStepper
{
public:
  typedef typename wave_t::pos_t pos_t;
  typedef typename wave_t::index_t index_t;
  typedef typename index_t::shape_t shape_t;
  
  splitFourierStepper(shape_t & shape_ ): waveTmp(shape_),index(shape_){unity=1;}
  
  void attachWavefunction(wave_t & wave_)
  {
    wave=&wave_;
    waveTmp=(*wave);
    fS.setPlan(fourierPlanForward,wave->value().first(),waveTmp.value().first(),wave->size(),FFTW_FORWARD);
    
    fS.setPlan(fourierPlanBackward,waveTmp.value().first(),wave->value().first(),wave->size(),FFTW_BACKWARD);
    
    fS.getFrequencies(waveTmp.grid(),wave->getLength() );
    
  }
  
  void attachPotential(potential_t& potential_)
  {
    potential=&potential_;
  }
  
  void stepSplitFourier(const complex<double> &tau)
  {
    
    fftw_execute(fourierPlanForward);
    
    
    for(int i=0;i<wave->size();i++)
      {
	waveTmp.value()(i)*=exp(-tau/2.*norm(waveTmp.grid()(i)));
      }
    
    /*
    for(int i=0;i<wave->size();i++)
      {
	//waveTmp.value()(i)*=-norm(waveTmp.grid()(i));
	waveTmp.value()(i)*=waveTmp.grid()(i);
      }  
    */
    
    fftw_execute(fourierPlanBackward);
    wave->value()/=(double)wave->size();
  }

  
  void stepPotential(const complex<double> &delta)
  {
    
    potential->evaluate(wave->value(),wave->grid(),waveTmp.value());
    
    for(index=0;index<index.size();index++)
      {
	wave->value()(index)*=exp(-delta*waveTmp.value()(index));
      }
    
  }

  
protected:
  wave_t* getWave(){return wave;}
private:
  index_t index;
  wave_t *wave;
  wave_t waveTmp;
  potential_t *potential;
  complex<double> unity;
  vector<pos_t> frequencies;
  fftw_plan fourierPlanForward;
  fftw_plan fourierPlanBackward;
  
  fourierStrategy<typename wave_t::dimensions_t,typename wave_t::value_t> fS;
};


template<class wave_t,class potential_t>
class splitFourierStepper1Order : public splitFourierStepper<wave_t,potential_t>
{
public:
  typedef typename splitFourierStepper<wave_t,potential_t>::shape_t shape_t;
  splitFourierStepper1Order(shape_t & shape_) :  splitFourierStepper<wave_t,potential_t>(shape_){};
  
  void stepRealTime(double tau)
  {
    
    this->stepPotential(tau);
    this->stepSplitFourier(tau);
    this->getWave()->normalize();
  }
  private:
};


template<class stepper_t>
class stepperMixture
{
public:
  typedef typename stepper_t::wave_t wave_t;
  typedef typename stepper_t::potential_t potential_t;
  typedef typename wave_t::pos_t pos_t;
  typedef typename wave_t::index_t index_t;
  typedef typename index_t::shape_t shape_t;

  stepperMixture(shape_t &shape1,shape_t &shape2) : stepper1(shape1),stepper2(shape2){};
  
  void stepRealTime(double tau)
  {
    stepper2.setPotentialArg(stepper1.getWave());
    stepper1.stepRealTime(tau);
    stepper2.stepRealTime(tau);
    
  }
  
private:
  
  stepper_t stepper1;
  stepper_t stepper2;
  
};

class trappedGPEquation1D
{
public:
  typedef wave1D wave_t;
  typedef harmonicPotential1D potential_t;
  
};

class gPEquationLY1D
{
public:
  typedef wave1D wave_t;
  typedef LYPotential potential_t;
  
};

class gPEquationLY1DTrapped
{
public:
  typedef wave1D wave_t;
  typedef LYPotentialHarmonicTrapping potential_t;
  
};

class gP1DModel
{
public:
  typedef wave1D wave_t;
  typedef chemicalPotentialModelFitted potential_t;
  
};



template<class strategy_t>
class GPSolver
{
public:
  typedef typename strategy_t::wave_t wave_t;
  typedef typename strategy_t::potential_t potential_t;
  
  typedef typename wave_t::index_t index_t;
  typedef typename wave_t::positions_t positions_t;
  
  typedef splitFourierStepper1Order <wave_t,potential_t> stepper_t;
  typedef typename wave_t::values_t values_t;
  typedef typename index_t::shape_t shape_t;
  typedef typename wave_t::pos_t pos_t;
  
  GPSolver(shape_t bins,pos_t xMin,pos_t xMax) : index(bins),waveObject(bins),potential(),stepper(bins)
  {
    positions_t &positions=waveObject.grid();
    values_t& waveValues=waveObject.value();;
    
    for(index=0;index<index.size();index++)
      {
	positions(index)=index();
	positions(index)*=(xMax-xMin);
	positions(index)/=(bins);
	positions(index)+= xMin;
	waveValues(index)=exp(-norm(positions(index)));
      }
    
    waveObject.setMin(xMin);
    waveObject.setMax(xMax);
    
    stepper.attachWavefunction(waveObject);
    stepper.attachPotential(potential);
    
    
    filename="wave.dat";
    tau=0.01;
    
  }
  
  void setTimeStep(double tau_){tau=tau_;}
  void setParameter(string name,double param){potential.setParameter(name,param);}
  
  void setOutputFile(string filename_){filename=filename_;}
  
  void run(double time,double timeBlock)
  {
    int nBlocks=(int)(time/timeBlock);
    int stepsPerBlock=timeBlock/tau;
    for(int i=0;i<nBlocks;i++)
      {
	for(int j=0;j<stepsPerBlock;j++)
	  {
	    stepper.stepRealTime(tau);
	  }
	
	ofile.open(filename.c_str());
	ofile<<waveObject<<endl;
	ofile.close();
      }
  }
  
private:
  double tau;
  ofstream ofile;
  string filename;
  int blockSize;
  wave_t waveObject;
  potential_t potential;
  index_t index;
  stepper_t stepper;
  
};

template<class OT,class OT1>
  ostream& operator<<(ostream& out, const waveGeneral<OT,OT1>& wave)
{
  typedef typename waveGeneral<OT,OT1>::index_t index_t;
  index_t index(wave.shape());
  
  for(index=0;index<index.size();index++)
    {
      out<<wave.grid()(index)<<" "<<(wave.value()(index)).real()<<" "<<(wave.value()(index)).imag()<<endl;
    }
  
  return out;
  
}

#endif
