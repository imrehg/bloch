#define GMM_USES_LAPACK
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
#include "atom.cpp" //for atom struct
using namespace gmm;
using namespace std;

typedef double doub;
#include <ctime>
#include <time.h>
int neq=32;

int Real(int i,int j)
{
  int c=0;
     if(i==j) c=i;
     else
       if(j>i) c=(2*neq-1-i)*i/2+(j-i)+neq-1;
       else c=(2*neq-1-j)*j/2+(i-j)+neq-1;

  return c;

}

int main(){

    dense_matrix<doub>  A(1000,1000);
    dense_matrix<doub>  B(1000,1000);
    dense_matrix<doub>  C(1000,1000);

    A((Real(1,1)),Real(1,1))=0;
    cout<<Real(1,1);
    fill_random(A,0.5);
    fill_random(B,0.5);

    time_t start=clock();

    for(int i=0;i<1;i++)
        mult(A,B,C);
// cout<<A;

    cout<<(clock()-start)*1.0/CLOCKS_PER_SEC<<endl;





}
