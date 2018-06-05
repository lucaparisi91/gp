#include "tools.h"

double  norm(double x){return x*x;};
double norm(const shortVector<double,2> & vec){return (norm(vec[0]) + norm(vec[1])) ;}

double norm(const shortVector<double,1> & vec){return vec[0]*vec[0];}

double  norm(const shortVector<double,3> & vec){return (norm(vec[0]) + norm(vec[1]) + norm(vec[2])) ;}

namespace indexTools
{
  int getRank(const shortVector<int,3> & vec){return vec[0]*vec[1]*vec[2];}
  int getRank(const shortVector<int,1> & vec){return vec[0];}
  
  int flatten(const shortVector<int,3> & vec,const shortVector<int,3> & shape){return vec[2] + shape[2]*( vec[1] + shape[1]* vec[0]); }
  int flatten(const shortVector<int,1> & vec,const shortVector<int,1> & shape){return vec[0]; }
  shortVector<int,3> flattenInverse(int index,const shortVector<int,3> & shape)
{
  shortVector<int,3> multIndex;
  multIndex[2]=index%shape[2];
  index=index/shape[2];
  multIndex[1]=index%shape[1];
  index=index/shape[1];
  multIndex[0]=index%shape[0];
  index=index/shape[0];
  
  return multIndex;
}
 
 shortVector<int,1> flattenInverse(int index,const shortVector<int,1> & shape){shortVector<int,1> ret;ret[0]=index;return ret;}
};
