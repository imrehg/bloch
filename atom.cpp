#include <math.h>
#include <algorithm>
#include <iostream> //FOR file IO

struct Atom {

  int factorial (int num){
   if (num==1||num==0)
    return 1;
    return factorial(num-1)*num; // recursive call
  }

  int bad_values (int j1,int j2,int j3,int m1,int m2,int m3){
        if (j1<fabs(j2-j3)||j1>(j2+j3))
            return 1;
        if (fabs(m1)>j1 || fabs(m2)>j2 || fabs(m3)>j3)
            return 1;
        if (m1+m2+m3 !=0)
            return 1;
        return 0;
  }


  long double threej (int j1,int j2,int j3,int m1,int m2,int m3){

   if (bad_values (j1,j2,j3,m1,m2,m3)){
       std::cout<<"Bad Values";
       return 0;
   }

   int jphase = (-1)^(j1-j2-m3);
   int fac[10];

   fac[0] = factorial(j1+j2-j3);
   fac[1] = factorial(j1-j2+j3);
   fac[2] = factorial(-j1+j2+j3);
   fac[3] = factorial(j1+m1);
   fac[4] = factorial(j1-m1);
   fac[5] = factorial(j2+m2);
   fac[6] = factorial(j2-m2);
   fac[7] = factorial(j3+m3);
   fac[8] = factorial(j3-m3);
   fac[9] = factorial(j1+j2+j3+1);


   long double jprodfac = sqrt(1.0/fac[9]);

   for  (int i=0 ; i<9 ; i++)
      jprodfac=jprodfac*sqrt(1.0*fac[i]);

    int kmax = (std::min( std::min( (j1+j2-j3) , (j1-m1) ) , (j2+m2) ) );
    int kmin = (std::max( std::max( 0 , -(j3-j2+m1)) , -(j3-j1-m2)));
    long double jsum=0;

    for (int k=kmin;k<kmax+1;k++)
        jsum += ((-1)^k)*1.0 / (factorial(k)*factorial(j1+j2-j3-k)*factorial(j1-m1-k)*factorial(j2+m2-k)*factorial(j3-j2+m1+k)*factorial(j3-j1-m2+k));


    return jphase*jprodfac*jsum;

  }




  };



