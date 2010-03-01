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
using namespace std;
const long double pi=3.14159265358979323846264338327950288419716939937511;
const int npulse=500000000;
int ninterval_1=50,ninterval_2=500;//npulse = number of pulse; interval_1 =steps in interval 1 ..
long double period0=10.87827848197104833208251261254802895928271242476719;
long double frequency=0,peakO=1.34163815218652164669542605053/2,FWHM=0.0007; //about 150uW/cm2 about 1ps
const int  neq=4,ninterval=npulse*(ninterval_1+ninterval_2); // neq= nuber of equations, nexp= terms of expansion, ninterval= iteration terms
int nexp=12;
long double r[neq]={0.0052227*2*pi,0.0052227*2*pi,0,0};//total decay constant
long double R[neq]={0.0052227*2*pi,0.0052227*2*pi,0,0};//relaxation rate
long double Rc[neq][neq]={{0,0,0,0},
                          {0,0,0,0},
                          {0,0,0,0.0000005*2*pi},
                          {0,0,0.0000005*2*pi,0}};//coherence relaxation rate
long double A[neq][neq]={{0,0,0,0},{0,0,0,0},{0.0052227*2*pi/2,0.0052227*2*pi/2,0,0},{0.0052227*2*pi/2,0.0052227*2*pi/2,0,0}};//Einstein A coefficient
long double R_L[neq]={0,0,0,0};//laser line width
long double lasDe = 0;
long double d0[neq][neq]={{0,-0.20124*2*pi,-0.20124*2*pi+lasDe,-0.20124*2*pi-9.192631*2*pi+lasDe},
                         {+0.20124*2*pi,0,+lasDe,-9.192631*2*pi+lasDe},
                         {0.20124*2*pi-lasDe,0-lasDe,0,-9.192631*2*pi},
                         {+9.192631*2*pi+0.20124*2*pi-lasDe,+9.192631*2*pi-lasDe,9.192631*2*pi,0}};//laser */
long double d[neq][neq]={{0,-0.20124*2*pi,-0.20124*2*pi+lasDe,-0.20124*2*pi-9.192631*2*pi+lasDe},
                         {+0.20124*2*pi,0,+lasDe,-9.192631*2*pi+lasDe},
                         {0.20124*2*pi-lasDe,0-lasDe,0,-9.192631*2*pi},
                         {+9.192631*2*pi+0.20124*2*pi-lasDe,+9.192631*2*pi-lasDe,9.192631*2*pi,0}};//laser */
long double y0I[neq][neq]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};//initial condition
long double y0R[neq][neq]={{0,0,0,0},{0,0,0,0},{0,0,0.5,0},{0,0,0,0.5}};//initial condiion
double phase=0;
int pulse_average=100;
//void fun(long double ***,long double ***,long double **,int );//聯立方程式
//void solve(long double ***,long double ***,long double ***,long double ***,long double ***,long double*,int);//演算法
//long double ReRabi(long double,long double,long double );//脈衝包絡線函數(實部)
//long double ImRabi(long double,long double,long double );//脈衝包絡線函數(虛部)
//int factorial (int);
//void fun_Matrix(long double *****,long double **);
//void solve_Martix(long double ***,long double ****,long double ****,long double *);
//void Matrix_Multiply(long double ****,long double ****);
//int sweep (int,int,long double,long double,int,int,long double);

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


int factorial (int num)
{
 if (num==1)
  return 1;
 return factorial(num-1)*num; // recursive call
}

long double ReRabi(long double x,long double period,long double peak)//脈衝包絡線函數(實部)，高斯函數*Re[e^{-i*phase}]
{
  long double value=0,time=0,factor=0;
  int i=0;
  if(x<0.5*period)
    value=exp(-pow(x/FWHM,2))*cos(-i*phase);
  else
  {
  i=int ((x-0.5*period)/period+1);
  value=exp(-pow((x-i*period)/FWHM,2))*cos(-i*phase);
  }
  return peak*value;
}

long double ImRabi(long double x,long double period,long double peak)//脈衝包絡線函數(虛部)，高斯函數*Im[e^{-i*phase}]
{
  long double value=0,time=0,factor=0;
  int i=0;
  if(x<0.5*period)
    value=exp(-pow(x/FWHM,2))*sin(-i*phase);
  else
  {
  i=int ((x-0.5*period)/period+1);
  value=exp(-pow((x-i*period)/FWHM,2))*sin(-i*phase);
  }

  return peak*value;
}
//      long double ****C= new long double***[neq*2];
//
//      for(int i=0;i<neq;i++){
//         C[i]=new long double**[neq];
//         C[i+neq]=new long double**[neq];
//           for(int j=0;j<neq;j++){
//              C[i][j]=new long double*[neq*2];
//              C[i+neq][j]=new long double*[neq*2];
//                for(int k=0;k<neq*2;k++){
//                   C[i][j][k]=new long double[neq];
//                   C[i+neq][j][k]=new long double[neq];
//             }
//          }
//       }
//
//      for(int a=0;a<neq;a++){
//           for(int b=0;b<=a;b++)
//             for(int c=0;c<neq;c++)
//                 for(int d=0;d<neq;d++){
//                     C[a+neq][b][c][d]=0;
//                     C[a+neq][b][c+neq][d]=0;
//                     C[a][b][c+neq][d]=0;
//                     C[a][b][c][d]=0;
//                          for(int e=0;e<neq;e++)
//                                 for(int f=0;f<neq;f++){
//                                        C[a+neq][b][c][d]+=B[a+neq][b][e+neq][f]*A[e+neq][f][c][d]+B[a+neq][b][e][f]*A[e][f][c][d];
//                                        C[a+neq][b][c+neq][d]+=B[a+neq][b][e+neq][f]*A[e+neq][f][c+neq][d]+B[a+neq][b][e][f]*A[e][f][c+neq][d];
//                                        C[a][b][c+neq][d]+=B[a][b][e+neq][f]*A[e+neq][f][c+neq][d]+B[a][b][e][f]*A[e][f][c+neq][d];
//                                        C[a][b][c][d]+=B[a][b][e+neq][f]*A[e+neq][f][c][d]+B[a][b][e][f]*A[e][f][c][d];
//                                 }
//
//                                 C[b+neq][a][c][d]=-C[a+neq][b][c][d];
//                                 C[b+neq][a][c+neq][d]=-C[a+neq][b][c+neq][d];
//                                 C[b][a][c+neq][d]=C[a][b][c+neq][d];
//                                 C[b][a][c][d]=C[a][b][c][d];
//
//                 }}
//
//       for(int a=0;a<neq;a++)
//           for(int b=0;b<=a;b++)
//             for(int c=0;c<neq;c++)
//                 for(int d=0;d<neq;d++){
//                            A[a+neq][b][c+neq][d]=C[a+neq][b][c+neq][d];
//                            A[b+neq][a][c+neq][d]=-C[a+neq][b][c+neq][d];
//                            A[a+neq][b][c][d]=C[a+neq][b][c][d];
//                            A[b+neq][a][c][d]=-C[a+neq][b][c][d];
//                            A[a][b][c+neq][d]=C[a][b][c+neq][d];
//                            A[b][a][c+neq][d]=C[a][b][c+neq][d];
//                            A[a][b][c][d]=C[a][b][c][d];
//                            A[b][a][c][d]=C[a][b][c][d];
//                            }//M'=M(t)*M
//
//    for(int i=0;i<neq;i++)
//        for(int j=0;j<neq;j++)
//            for(int k=0;k<neq;k++){
//                  delete[] C[i+neq][j][k];
//                  delete[] C[i][j][k];
//                  delete[] C[i+neq][j][k+neq];
//                  delete[] C[i][j][k+neq];
//            }
//
//     for(int i=0;i<neq;i++)
//         for(int j=0;j<neq;j++){
//               delete[] C[i][j];
//               delete[] C[i+neq][j];
//         }
//
//     for(int i=0;i<neq;i++){
//               delete[] C[i];
//               delete[] C[i+neq];
//     }
//
//      delete[] C;
//
//
//
//}

void fun_Matrix(gmm::row_matrix< gmm::wsvector<long double> > &Trans,long double **H)//B imaginary part ; C real part; H rabi requence; w time
{

  for (int i=0;i<neq;i++){
     Trans(RealComp(i,i),RealComp(i,i))=-r[i];
     for (int j=0;j<neq;j++){
          if(i<j){
          Trans(RealComp(i,i),ImagComp(i,j))+=2*H[j][i];
          } else Trans(RealComp(i,i),ImagComp(i,j))+=-2*H[j][i];

          Trans(RealComp(i,i),RealComp(j,j))+=A[i][j];
      }
 } // Bloch eq for population part

 for (int i=1;i<neq;i++){
    for (int j=0;j<i;j++){
       Trans(RealComp(i,j),RealComp(i,j))=-(R[i]+R[j]+R_L[i])/2-Rc[i][j];
       Trans(RealComp(i,j),ImagComp(i,j))=d[i][j];
       for (int l=0;l<neq;l++){
           if(l<j){
             Trans(RealComp(i,j),ImagComp(l,j))+=-H[i][l];
           } else{
             Trans(RealComp(i,j),ImagComp(l,j))+=H[i][l];
           }
           if(i<l){
             Trans(RealComp(i,j),ImagComp(i,l))+=H[l][j];
           } else{
             Trans(RealComp(i,j),ImagComp(i,l))+=-H[l][j];
           }
         }
   }
 }   // Bloch eq for off diagonal real part


 for (int i=0;i<neq;i++){
   for (int j=0;j<i;j++){
       Trans(ImagComp(i,j),ImagComp(i,j))=-(R[i]+R[j]+R_L[i])/2-Rc[i][j];
       Trans(ImagComp(i,j),RealComp(i,j))=-d[i][j];
      for (int l=0;l<neq;l++){
            Trans(ImagComp(i,j),RealComp(l,j))+=H[i][l];
            Trans(ImagComp(i,j),RealComp(i,l))+=-H[l][j];
      }
    }
  }    //  Bloch eq for off diagonal imaginary part


}

void solve_Martix(long double ***M, gmm::row_matrix< gmm::wsvector<long double> >&Trans, gmm::row_matrix< gmm::wsvector<long double> >&Trans_Ave,long double *T)// solve(presultI,presultR,M,k)
{



   cout<<Trans<<endl;


  for(int t=1;t<(ninterval_1+ninterval_2)+1;t++){

      gmm::row_matrix< gmm::wsvector<long double> > Trans_B(neq*neq,neq*neq),Trans_I(neq*neq,neq*neq),Trans_C(neq*neq,neq*neq),Trans_E(neq*neq,neq*neq);
       for(int i=0;i<neq*neq;i++)
             Trans_I(i,i)=1;
      fun_Matrix(Trans_E,M[t-1]);

      for(int j=1;j<=nexp;j++){
        gmm::mult(Trans_I,Trans_E,Trans_C);
        gmm::copy(Trans_C,Trans_I);

         for(int a=0;a<neq*neq;a++)
             for(int b=0;b<neq*neq;b++)
                      Trans_B(a,b)+=Trans_I(a,b)*pow((T[t]-T[t-1]),j)/factorial(j);
      }
        gmm::add(Trans_B,Trans_Ave);
        gmm::mult(Trans_B,Trans,Trans_C);
        gmm::copy(Trans_C,Trans);
        cout<<Trans_C<<endl;
      } cout<<"End of Solve"<<endl;
  }





void adj_detune(long double detune)
{

 d[0][2]=d0[0][2]+detune;
 d[0][3]=d0[0][3]+detune;
 d[1][2]=d0[1][2]+detune;
 d[1][3]=d0[1][3]+detune;
 d[2][0]=d0[2][0]-detune;
 d[3][0]=d0[3][0]-detune;
 d[2][1]=d0[2][1]-detune;
 d[3][1]=d0[3][1]-detune;

}




int sweep(int steps,int total_steps,long double PeakPower,long double convergence,int conS,int expN,int n1, int n2,long double detune)
{


  long double phase=0;
  ninterval_1 =n1;
  ninterval_2 =n2;
  fstream file1,file2;//file1:紀錄輸入的參數。file2://紀錄計算結果
  peakO = PeakPower/150*1.34163815218652164669542605053/2;
  nexp=expN;
  stringstream strstream;
  string filename;
  strstream<<PeakPower<<"uWcm2_"<<convergence<<"_conS_"<<conS<<"_O="<<nexp<<"_N1_"<<n1<<"_N2_"<<n2<<"_D_"<<detune/2/pi*1000<<"MHz.txt";
  strstream>>filename;
  cout<<filename.c_str()<<endl;
  file2.open(filename.c_str(),ios::out | ios::trunc);
  file1.open("inputMP.txt", ios::out | ios::trunc);
  file2.precision(15);
  adj_detune(detune);
  pulse_average=conS;

#pragma omp num_threads(1)
#pragma omp parallel for
for(int thread=0;thread<2;thread++)
{
   long double *Time= new long double[ninterval_1+ninterval_2+1];
   long double ***M= new long double**[ninterval_1+ninterval_2+1];//Rabifrequence*2

      for(int i=0;i<(ninterval_1+ninterval_2+1);i++){
          M[i]=new long double*[neq];
          for(int j=0;j<neq;j++){
              M[i][j]=new long double[neq];
           }
     }




for(int m=0;m<=steps;m++)
{

   cout<<m<<endl;
   long double De=m*(pow(-1,omp_get_thread_num()))*1.0/total_steps;
   long double period=10.878278481971048332082512612548/100*(100+De);
   long double peak=peakO*(100+De)/100;
   long double interval_1=FWHM*10,interval_2=period-interval_1;
   long double dt_1=interval_1/ninterval_1,dt_2=interval_2/ninterval_2;
   interval_2=period-interval_1;
   dt_2=interval_2/ninterval_2;


    gmm::row_matrix< gmm::wsvector<long double> > Trans(neq*neq,neq*neq),Trans_AVE(neq*neq,neq*neq);
    gmm::col_matrix< gmm::wsvector<long double> > Result(neq*neq,pulse_average+1);

    for (int i=0;i<neq;i++){ //initailizing
    for (int j=0;j<neq;j++){
            Result(ImagComp(i,j),0)=y0I[i][j];
            Result(RealComp(i,j),0)=y0R[i][j];
            }
    }

    for(int i=0;i<neq*neq;i++){
        Trans(i,i)=1.0;
        Trans_AVE(i,i)=1.0;
    }

    for(int k=0;k<(ninterval_1+ninterval_2+1);k++){

            long double buffer=0;

            if( k%(ninterval_2+ninterval_1)>=ninterval_1/2 && (ninterval_2+ninterval_1/2)>k%(ninterval_2+ninterval_1) )
               buffer=dt_2;
            else
                buffer=dt_1;

           if(k==0)
              Time[k]=0;
            else
              Time[k]=Time[k-1]+buffer;

              buffer= Time[k]-period*int(Time[k]/period);

              M[k][0][0]=0;
              M[k][0][1]=0;
              M[k][0][2]=0.38188130791298666722*ReRabi(buffer,period,peak);
              M[k][0][3]=0.27277236279499047658*ReRabi(buffer,period,peak);
              M[k][1][0]=0;
              M[k][1][1]=0;
              M[k][1][2]=0.14433756729740644113*ReRabi(buffer,period,peak);
              M[k][1][3]=0.43301270189221932338*ReRabi(buffer,period,peak);
              M[k][2][0]=0.38188130791298666722*ReRabi(buffer,period,peak);
              M[k][2][1]=0.14433756729740644113*ReRabi(buffer,period,peak);
              M[k][2][2]=0;
              M[k][2][3]=0;
              M[k][3][0]=0.27277236279499047658*ReRabi(buffer,period,peak);
              M[k][3][1]=0.43301270189221932338*ReRabi(buffer,period,peak);
              M[k][3][2]=0;
              M[k][3][3]=0;
           }

solve_Martix(M,Trans,Trans_AVE,Time);
cout<<Trans<<endl;
int k=0,flag=0;
double diff=0,diff0=0;

while(flag<pulse_average){


     for(int a=0;a<neq;a++)
           for(int b=0;b<neq;b++){
           Result(RealComp(a,b),(k+1)%(pulse_average+1))=0;
           Result(ImagComp(a,b),(k+1)%(pulse_average+1))=0;}

           gmm::mult(Trans,gmm::mat_col(Result,(k)%(pulse_average+1)),gmm::mat_col(Result,(k+1)%(pulse_average+1)));
//           cout<<gmm::mat_col(Result,(k+1)%(pulse_average+1));
           k+=1;

      if(k>pulse_average){
          diff0=0;
          for(int c=0;c<neq;c++){
            for(int d=0;d<neq;d++)
                diff0=diff0+abs((1-Result(RealComp(c,d),(k%(pulse_average+1)))/Result(RealComp(c,d),(k-pulse_average)%(pulse_average+1))))/(neq*neq);

          }
        if(diff0<diff){
           if(abs(diff)<(convergence))
              flag+=1;
           else
             flag=0;}
         else
            flag=0;

         diff=diff0;

      }

     if(k==(npulse-1))
       flag=pulse_average+1;

}

long double buffer=0,bufferC=0;

                 for(int d=0;d<neq*neq;d++){
                  buffer+=Trans_AVE(0,d)*Result(d,k%(pulse_average+1));
                  buffer+=Trans_AVE(1,d)*Result(d,k%(pulse_average+1));
                 }
            for(int c=0;c<neq;c++)
                 bufferC+=Result(c,k%(pulse_average+1));

buffer=buffer/(ninterval_1+ninterval_2+1);

       file2<<setiosflags(ios::left)<<setw(30)<<1/period;
       file2<<setiosflags(ios::left)<<setw(30)<<buffer;
       file2<<setiosflags(ios::left)<<setw(30)<<bufferC;
       file2<<setiosflags(ios::left)<<setw(30)<<k;
       file2<<setiosflags(ios::left)<<setw(30)<<m<<endl;

}

///////////////////////////////End of Sweeping//////////////////////////////////




       for(int i=0;i<ninterval_1+ninterval_2+1;i++)
         for(int j=0;j<neq;j++){
               delete[] M[i][j];
       }

       for(int i=0;i<ninterval_1+ninterval_2+1;i++)
               delete[] M[i];

      delete[] Time;
      delete[] M;



}


  return 0;




}


