#include <iostream> //FOR file IO
#include <fstream>//For file IO
//#include <cstdlib> //
//#include <cctype>//
//#include <cstring>
#include <string>//to use string
#include <sstream>//to use sstream
#include <cmath>//For sin cos functions
#include <iomanip>//For  setiosflags
#include <ctime>//For timer
#include <omp.h>//For openmp
#include <gmm/gmm.h>//for gmm library
int neq=4;
//int RealComp(int ,int );
using namespace std;
int RealComp(int i,int j)
{
  int c=0;
     if(i==j) c=i;
     else
       if(j>i) c=(2*neq-1-i)*i/2+(j-i)+neq-1;
       else c=(2*neq-1-j)*j/2+(i-j)+neq-1;

  return c;

}

int ImagComp(int i,int j)
{
  int c=0;
     if(i==j) c=0;
     else
       if(j>i) c=(2*neq-1-i)*i/2+(j-i)+neq-1+neq*(neq-1)/2;
       else c=(2*neq-1-j)*j/2+(i-j)+neq-1+neq*(neq-1)/2;

  return c;

}

int main()
{
  gmm::row_matrix< gmm::wsvector<long double> > Trans_E(neq*neq,neq*neq),Trans_I(neq*neq,neq*neq),Trans_B(neq*neq,neq*neq),Trans_C(neq*neq,neq*neq);
  for(int i=0;i<10;i++)
  gmm::mult(Trans_I,Trans_E,Trans_C);



}


