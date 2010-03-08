#include <cmath>

struct Atomic {

  int factorial (int num){
   if (num==1)
    return 1;
    return factorial(num-1)*num; // recursive call
  }




  long double threej (int j1,int j2,int j3,int m1,int m2,int m3){

   int bad_values (int j1,int j2,int j3,int m1,int m2,int m3){
        if (j1<abs(j2-j3)||j1>(j2+j3))
            return 1;
        if (abs(m1)>j1 || abs(m2)>j2 || abs(m3)>j3)
            return 1;
        if (m1+m2+m3 !=0)
            return 1;
        return 0;
  }

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

  long double sixj(int j1, int j2, int j3, int l1, int l2, int l3){

    int bad_values(j1,j2,j3,l1,l2,l3){
        if (j1<(abs(j2-j3)) || j1>(j2+j3)):
            return 1;
        if (j1<(abs(l2-l3)) || j1>(l2+l3)):
            return 1;
        if (l1<(abs(j2-l3)) || l1>(j2+l3)):
            return 1;
        if (l1<(abs(l2-j3)) || l1>(l2+j3)):
            return 1;
        return 0;
        }


    long double delta(a,b,c){
        return sqrt(1.0*((factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c))/factorial(a+b+c+1)));
    }

    if bad_values(j1,j2,j3,l1,l2,l3):
        return 0

    int jphase=(-1)^(j1+j2+l1+l2);
    long double proddelt=delta(j1,j2,j3)*delta(l1,l2,j3)*delta(l1,j2,l3)*delta(j1,l2,l3);


    int kmax = min(min(j1+j2+l1+l2+1,j1+j2-j3),min(min(l1+l2-j3,j1+l2-l3),l1+j2-l3))
    int kmin = max(0, -j1-l1+j3+l3, -j2-l2+j3+l3)
    long double jsum = 0

    for (int k=kmin; k<kmax+1;k++){
        long double jsfac[8];
        jsfac[0] = factorial(j1+j2+l1+l2+1-k);
        jsfac[1] = factorial(k);
        jsfac[2] = factorial(j1+j2-j3-k);
        jsfac[3] = factorial(l1+l2-j3-k);
        jsfac[4] = factorial(j1+l2-l3-k);
        jsfac[5] = factorial(l1+j2-l3-k);
        jsfac[6] = factorial(-j1-l1+j3+l3+k);
        jsfac[7] = factorial(-j2-l2+j3+l3+k);
        long double product=1;

        for (int j=1; j<8; j++)
          product=product/jsfac[j];

          jsum += (-1)^k * jsfac[0]*product;

    }

    return jphase*proddelt*jsum






  }














}
