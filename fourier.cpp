#include "fourier.h"

shortVector<int, 3> shiftNegativeFreqIndex( const shortVector<int, 3> &vec,const shortVector<int,3> & shape)
{ 
  shortVector<int, 3> vecRet;
  vecRet[0]=shiftNegativeFreq(vec[0],shape[0]);
  vecRet[1]=shiftNegativeFreq(vec[1],shape[1]);
  vecRet[2]=shiftNegativeFreq(vec[2],shape[2]);
  
  return vecRet;
  
}

shortVector<int, 1> shiftNegativeFreqIndex( const shortVector<int, 1> vec,const shortVector<int,1> & shape)
{
  
  shortVector<int, 1> vecRet;
  vecRet[0]=shiftNegativeFreq(vec[0],shape[0]);
  
  return vecRet;
  
}

int shiftNegativeFreq(int index,int size)
{
  if (index <= size/2)
    {
      return index;
    }
  else
    {
      return -(size-index);
    }
  
}
