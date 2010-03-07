#include <cmath>

struct Atomic {

  int factorial (int num){
   if (num==1)
    return 1;
    return factorial(num-1)*num; // recursive call
  }

  int bad_values (int j1,int j2,int j3,int m1,int m2,int m3){
        if (j1<abs(j2-j3)||j1>(j2+j3))
            return 1;
        if (abs(m1)>j1 || abs(m2)>j2 || abs(m3)>j3)
            return 1;
        if (m1+m2+m3 !=0)
            return 1;
        return 0;
  }


  long double 3j (int j1,int j2,int j3,int m1,int m2,int m3){

   if (bad_values (j1,j2,j3,m1,m2,m3)){
       std::cout<<"Bad Values";
       return 0;
   }

   int jphase = (-1)^(j1-j2-m3);
   int fac[10]=0;

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


   long double jprodfac = sqrt(1/fac[9]);

   for  (int i=0 ; i<9 ; i++)
      jprodfac=fprodfac*sqrt(fac[i]);

    int kmax = int(min( min( (j1+j2-j3) , (j1-m1) ) , (j2+m2) ) );
    int kmin = int(max(max( 0 , -(j3-j2+m1)) , -(j3-j1-m2)));
    int jsum=0;

    for (int k=kmin;k<kmax+1;k++)
        jsum += (-1)^k / (factorial(k)*factorial(j1+j2-j3-k)*factorial(j1-m1-k)*factorial(j2+m2-k)*factorial(j3-j2+m1+k)*factorial(j3-j1-m2+k));


    return jphase*jprodfac*jsum

  }







  }














}
