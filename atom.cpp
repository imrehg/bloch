#include <math.h>
#include <algorithm>
#include <iostream> //FOR file IO

struct Atom {

 int factorial (int num ){
   if (num==1.0||num==0.0)
    return 1;
    return factorial(num-1)*num; // recursive call
  }

  int bad_values (double j1,double j2,double j3,double m1,double m2,double m3){
        if (j1<fabs(j2-j3)||j1>(j2+j3))
            return 1;
        if (fabs(m1)>j1 || fabs(m2)>j2 || fabs(m3)>j3)
            return 2;
        if (m1+m2+m3 !=0)
            return 3;
        return 0;
  }


  long double threej (double j1,double j2,double j3,double m1,double m2,double m3){

   if (bad_values (j1,j2,j3,m1,m2,m3)){
       return 0;
   }

   int jphase = pow((-1.0),(j1-j2-m3));
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


   for  (int i=0 ; i<9 ; i++){
      jprodfac=jprodfac*sqrt(fac[i]);
      }


    int kmax = (std::min( std::min( (j1+j2-j3) , (j1-m1) ) , (j2+m2) ) );
    int kmin = (std::max( std::max( 0.0 , -(j3-j2+m1)) , -(j3-j1-m2)));

    long double jsum=0;

    for (int k=kmin;k<kmax+1;k++)
        jsum += pow((-1),k)*1.0 / (factorial(k)*factorial(j3-j2+m1+k)*factorial(j3-j1-m2+k)*factorial(j1+j2-j3-k)*factorial(j1-k-m1)*factorial(j2+m2-k));

    return jphase*jprodfac*jsum;

  }


  int bad_values6j(double j1,double j2,double j3,double l1,double l2,double l3){
        if (j1<(abs(j2-j3)) || j1>(j2+j3))
            return 1;
        if (j1<(abs(l2-l3)) || j1>(l2+l3))
            return 1;
        if (l1<(abs(j2-l3)) || l1>(j2+l3))
            return 1;
        if (l1<(abs(l2-j3)) || l1>(l2+j3))
            return 1;
        return 0;
        }

  long double delta(double a,double b,double c){
        return sqrt((1.0*(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c))/factorial(a+b+c+1)));
    }

  long double sixj(double j1, double j2, double j3, double l1, double l2, double l3){

    if (bad_values6j(j1,j2,j3,l1,l2,l3))
        return 0;

    int jphase=pow((-1.0),(j1+j2+l1+l2));
    long double proddelt=delta(j1,j2,j3)*delta(l1,l2,j3)*delta(l1,j2,l3)*delta(j1,l2,l3);


    int kmax = std::min(std::min(j1+j2+l1+l2+1,j1+j2-j3),std::min(std::min(l1+l2-j3,j1+l2-l3),l1+j2-l3));
    int kmin = std::max(std::max(0.0, -(-j1-l1+j3+l3)), -(-j2-l2+j3+l3));
    long double jsum = 0;

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

          jsum += pow(-1,k )* jsfac[0]*product;

    }
    return (jphase*proddelt*jsum);


  }

  long double coef(double q,double F1,double F2,double mf1, double mf2,double J1, double J2,double I){
    return pow(-1,F2-1+mf1)*sqrt(2*F1+1)*threej(F2,1,F1,mf2,q,-mf1)*pow(-1,F2+J1+1+I)*sqrt((2*F2+1)*(2*J1+1))*sixj(J1,J2,1,F2,F1,I);
  }


};



