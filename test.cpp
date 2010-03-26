#include "atom.cpp"
#include "gmm/gmm.h"
#include <iomanip>
const int neq=4;
using namespace gmm;
int D1_coef (int L,int F,int mf){
     if(F==3)
      return 32-(L*16+(mf+4));
     if(F==4)
      return 32-(L*16+7+(mf+5));
 }
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

int main(){

   col_matrix< wsvector<long double> > M(32,32);
   row_matrix< rsvector<long double> > A(32*32,32*32);
   rsvector<long double> EnergyDiff(32-1);

//  for(int i=0; i<2;i++)
//     for(int j=3; j<5;j++)
//        for(int t=-j;t<j+1;t++)
//          for(int l=0; l<2;l++)
//            for(int m=3; m<5;m++)
//               for(int n=-m;n<m+1;n++)
   Atom atom;
//   cout<<atom.coef(-1,1,0,4,4,-4,-3,0.5,0.5,3.5)<<endl;
//   cout<<atom.coef(1,0,1,4,4,-3,-4,0.5,0.5,3.5);
//
//   for(int i=0; i<2;i++)
//     for(int j=3; j<5;j++)
//        for(int k=-j;k<j+1;k++)
//          for(int l=0; l<2;l++)
//            for(int m=3; m<5;m++)
//               for(int n=-m;n<m+1;n++)
//                     M(D1_coef(i,j,k),D1_coef(l,m,n))+=atom.coef(+1,i,l,j,m,k,n,0.5,0.5,3.5)+atom.coef(-1,i,l,j,m,k,n,0.5,0.5,3.5);
//
// long double sum=0;
//      for(int j=3; j<=5;j++)
//         for(int k=-j;k<j+1;k++){
//                    std::cout<<j<<"  "<<k<<" ";
//                    std::cout<<atom.coef(0,0,1,4,j,k,k,0.5,1.5,3.5)<<std::endl;}
//             for(int m=3; m<5;m++)
//                for(int n=-m;n<m+1;n++)
//                  for(int q=-1;q<2;q++){
//                     A(D1_coef(0,j,k),D1_coef(1,m,n))+=pow(atom.coef(q,0,1,j,m,k,n,1.5,0.5,3.5),2);
//                     sum+=pow(atom.coef(q,0,1,j,m,k,n,1.5,0.5,3.5),2);

//                  }
//                  std::cout<<sum;
//
//    for(int i=0; i<2;i++)
//     for(int j=3; j<5;j++)
//        for(int k=-j;k<j+1;k++)
//          for(int l=0; l<2;l++)
//            for(int m=3; m<5;m++)
//               for(int n=-m;n<m+1;n++)
//                     M(D1_coef(i,j,k),k*neq+D1_coef(l,m,n))+=(atom.coef(+1,i,l,j,m,k,n,0.5,0.5,3.5)+atom.coef(-1,i,l,j,m,k,n,0.5,0.5,3.5));
//
//   EnergyDiff[8]=10;
//   EnergyDiff[24]=20;

//  std::cout << std::setprecision(10);
// std::cout<<RealComp(15,15);
 cout<<4%2;
 return 0;




}
