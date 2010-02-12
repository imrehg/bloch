#include <iostream> //FOR file IO
#include <fstream>//For file IO
//#include <cstdlib> //
//#include <cctype>
//#include <cstring>
#include <string>//to use string
#include <sstream>//to use sstream
#include <cmath>//For sin cos functions
#include <iomanip>//For  setiosflags
#include <ctime>//For timer
#include <omp.h>//For openmp
using namespace std;
const long double pi=3.14159265358979323846264338327950288419716939937511;
const int npulse=50000000;
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
void fun(long double ***,long double ***,long double **,int );//聯立方程式
void solve(long double ***,long double ***,long double ***,long double ***,long double ***,long double*,int);//演算法
long double ReRabi(long double,long double,long double );//脈衝包絡線函數(實部)
long double ImRabi(long double,long double,long double );//脈衝包絡線函數(虛部)
int factorial (int);
void fun_Matrix(long double *****,long double **);
void solve_Martix(long double ***,long double ****,long double ****,long double *);
void Matrix_Multiply(long double ****,long double ****);
double phase=0;
 int pulse_average=100;
int sweep (int,int,long double,long double,int,long double);


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

void Matrix_Multiply(long double ****A,long double ****B)//A=B*A ;C:Buffer
{
      long double ****C= new long double***[neq*2];

      for(int i=0;i<neq;i++){
         C[i]=new long double**[neq];
         C[i+neq]=new long double**[neq];
           for(int j=0;j<neq;j++){
              C[i][j]=new long double*[neq*2];
              C[i+neq][j]=new long double*[neq*2];
                for(int k=0;k<neq*2;k++){
                   C[i][j][k]=new long double[neq];
                   C[i+neq][j][k]=new long double[neq];
             }
          }
       }

      for(int a=0;a<neq;a++){
           for(int b=0;b<=a;b++)
             for(int c=0;c<neq;c++)
                 for(int d=0;d<neq;d++){
                     C[a+neq][b][c][d]=0;
                     C[a+neq][b][c+neq][d]=0;
                     C[a][b][c+neq][d]=0;
                     C[a][b][c][d]=0;
                          for(int e=0;e<neq;e++)
                                 for(int f=0;f<neq;f++){
                                        C[a+neq][b][c][d]+=B[a+neq][b][e+neq][f]*A[e+neq][f][c][d]+B[a+neq][b][e][f]*A[e][f][c][d];
                                        C[a+neq][b][c+neq][d]+=B[a+neq][b][e+neq][f]*A[e+neq][f][c+neq][d]+B[a+neq][b][e][f]*A[e][f][c+neq][d];
                                        C[a][b][c+neq][d]+=B[a][b][e+neq][f]*A[e+neq][f][c+neq][d]+B[a][b][e][f]*A[e][f][c+neq][d];
                                        C[a][b][c][d]+=B[a][b][e+neq][f]*A[e+neq][f][c][d]+B[a][b][e][f]*A[e][f][c][d];
                                 }

                                 C[b+neq][a][c][d]=-C[a+neq][b][c][d];
                                 C[b+neq][a][c+neq][d]=-C[a+neq][b][c+neq][d];
                                 C[b][a][c+neq][d]=C[a][b][c+neq][d];
                                 C[b][a][c][d]=C[a][b][c][d];

                 }}

       for(int a=0;a<neq;a++)
           for(int b=0;b<=a;b++)
             for(int c=0;c<neq;c++)
                 for(int d=0;d<neq;d++){
                            A[a+neq][b][c+neq][d]=C[a+neq][b][c+neq][d];
                            A[b+neq][a][c+neq][d]=-C[a+neq][b][c+neq][d];
                            A[a+neq][b][c][d]=C[a+neq][b][c][d];
                            A[b+neq][a][c][d]=-C[a+neq][b][c][d];
                            A[a][b][c+neq][d]=C[a][b][c+neq][d];
                            A[b][a][c+neq][d]=C[a][b][c+neq][d];
                            A[a][b][c][d]=C[a][b][c][d];
                            A[b][a][c][d]=C[a][b][c][d];
                            }//M'=M(t)*M

    for(int i=0;i<neq;i++)
        for(int j=0;j<neq;j++)
            for(int k=0;k<neq;k++){
                  delete[] C[i+neq][j][k];
                  delete[] C[i][j][k];
                  delete[] C[i+neq][j][k+neq];
                  delete[] C[i][j][k+neq];
            }

     for(int i=0;i<neq;i++)
         for(int j=0;j<neq;j++){
               delete[] C[i][j];
               delete[] C[i+neq][j];
         }

     for(int i=0;i<neq;i++){
               delete[] C[i];
               delete[] C[i+neq];
     }

      delete[] C;



}

void fun_Matrix(long double ****Trans,long double **H)//B imaginary part ; C real part; H rabi requence; w time
{

  for (int i=0;i<neq;i++){
     Trans[i][i][i][i]=-r[i];
     for (int j=0;j<neq;j++){
          Trans[i][i][i+neq][j]+=2*H[j][i];
          Trans[i][i][j][j]+=A[i][j];
      }
 } // Bloch eq for population part

 for (int i=1;i<neq;i++){
    for (int j=0;j<i;j++){
       Trans[i][j][i][j]=-(R[i]+R[j]+R_L[i])/2-Rc[i][j];
       Trans[i][j][i+neq][j]=d[i][j];
       for (int l=0;l<neq;l++){
             Trans[i][j][l+neq][j]+=-H[i][l];
             Trans[i][j][i+neq][l]+=H[l][j];
         }
       for(int a=0;a<2*neq;a++)
         for(int b=0;b<neq;b++)
            Trans[j][i][a][b]=Trans[i][j][a][b];
   }
 }   // Bloch eq for off diagonal real part


 for (int i=1;i<neq;i++){
   for (int j=0;j<i;j++){
       Trans[i+neq][j][i+neq][j]=-(R[i]+R[j]+R_L[i])/2-Rc[i][j];
       Trans[i+neq][j][i][j]=-d[i][j];
      for (int l=0;l<neq;l++){
            Trans[i+neq][j][l][j]+=H[i][l];
            Trans[i+neq][j][i][l]+=-H[l][j];
      }
      for(int a=0;a<2*neq;a++)
         for(int b=0;b<neq;b++)
            Trans[j+neq][i][a][b]=-Trans[i+neq][j][a][b];
    }
  }    //  Bloch eq for off diagonal imaginary part
}

void solve_Martix(long double ***M,long double ****Trans,long double ****Trans_Ave,long double *T)// solve(presultI,presultR,M,k)
{
     long double ****Trans_B= new long double***[neq*2],****Trans_E=new long double***[neq*2],****Trans_I=new long double***[neq*2];

      for(int i=0;i<neq*2;i++){
          Trans_B[i]=new long double**[neq];
          Trans_E[i]=new long double**[neq];
          Trans_I[i]=new long double**[neq];
           for(int j=0;j<neq;j++){
              Trans_B[i][j]=new long double*[neq*2];
              Trans_E[i][j]=new long double*[neq*2];
              Trans_I[i][j]=new long double*[neq*2];
                for(int k=0;k<neq*2;k++){
                    Trans_B[i][j][k]=new long double[neq];
                    Trans_E[i][j][k]=new long double[neq];
                    Trans_I[i][j][k]=new long double[neq];
           }
     }
  }

  for(int t=1;t<(ninterval_1+ninterval_2)+1;t++){

    for(int n=0;n<2*neq;n++)
     for(int i=0;i<neq;i++)
         for(int j=0;j<2*neq;j++)
             for(int k=0;k<neq;k++){
                   Trans_E[n][i][j][k]=0;
                if(n==j&&i==k){
                   Trans_B[n][i][j][k]=1;
                   Trans_I[n][i][j][k]=1;
                }
                 else{
                   Trans_B[n][i][j][k]=0;
                   Trans_I[n][i][j][k]=0;
                 }
             }//initialing

      fun_Matrix(Trans_E,M[t-1]);

      for(int j=1;j<=nexp;j++){
        Matrix_Multiply(Trans_I,Trans_E);
         for(int a=0;a<2*neq;a++)
           for(int b=0;b<neq;b++)
             for(int c=0;c<2*neq;c++)
                 for(int d=0;d<neq;d++)
                   {
                      Trans_B[a][b][c][d]+=Trans_I[a][b][c][d]*pow((T[t]-T[t-1]),j)/factorial(j);
                    }}

       for(int a=0;a<2*neq;a++)
           for(int b=0;b<neq;b++)
             for(int c=0;c<2*neq;c++)
                 for(int d=0;d<neq;d++)
                      Trans_Ave[a][b][c][d]+=Trans_B[a][b][c][d];

      Matrix_Multiply(Trans,Trans_B);

      }



     for(int i=0;i<neq*2;i++)
         for(int j=0;j<neq;j++)
             for(int k=0;k<neq*2;k++){
               delete[] Trans_B[i][j][k];
               delete[] Trans_E[i][j][k];
               delete[] Trans_I[i][j][k];
       }

     for(int i=0;i<neq*2;i++){
         for(int j=0;j<neq;j++){
               delete[] Trans_B[i][j];
               delete[] Trans_E[i][j];
               delete[] Trans_I[i][j];
       }}

     for(int i=0;i<neq*2;i++){
               delete[] Trans_B[i];
               delete[] Trans_E[i];
               delete[] Trans_I[i];
     }

      delete[] Trans_B;
      delete[] Trans_E;
      delete[] Trans_I;


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

int sweep(int steps,int total_steps,long double PeakPower,long double convergence,int expN,int n1, int n2,long double detune)
{


  long double phase=0;
  ninterval_1 =n1;
  ninterval_2 =n2;
  fstream file1,file2;//file1:紀錄輸入的參數。file2://紀錄計算結果
  peakO = PeakPower/150*1.34163815218652164669542605053/2;
  nexp=expN;
  stringstream strstream;
  string filename;
  strstream<<"Car_"<<PeakPower<<"uWcm2_"<<convergence<<"_O="<<nexp<<"_N1_"<<n1<<"_N2_"<<n2<<".txt";
  strstream>>filename;
  cout<<filename.c_str()<<endl;
  file2.open(filename.c_str(),ios::out | ios::trunc);
  file1.open("inputMP.txt", ios::out | ios::trunc);
  file2.precision(15);


   long double *Time= new long double[ninterval_1+ninterval_2+1];
   long double ***M= new long double**[ninterval_1+ninterval_2+1];//Rabifrequence*2

      for(int i=0;i<(ninterval_1+ninterval_2+1);i++){
          M[i]=new long double*[neq];
          for(int j=0;j<neq;j++){
              M[i][j]=new long double[neq];
           }
     }

 long double ****Trans= new long double***[neq*2],****Trans_AVE= new long double***[neq*2];
      for(int i=0;i<neq*2;i++){
          Trans[i]=new long double**[neq];
          Trans_AVE[i]=new long double**[neq];
           for(int j=0;j<neq;j++){
              Trans[i][j]=new long double*[neq*2];
              Trans_AVE[i][j]=new long double*[neq*2];
                for(int k=0;k<neq*2;k++){
                    Trans[i][j][k]=new long double[neq];
                    Trans_AVE[i][j][k]=new long double[neq];
           }
     }
    }
    long double ***presultR= new long double**[pulse_average+1];//所有時間點的數值存於此指標(real)
      for(int i=0;i<pulse_average+1;i++){
          presultR[i]=new long double*[neq];
          for(int j=0;j<neq;j++){
              presultR[i][j]=new long double[neq];
           }
     }
   long double ***presultI= new long double**[pulse_average+1];//所有時間點的數值存於此指標(imaginary) presultI[][o][o]為時間參數
      for(int i=0;i<pulse_average+1;i++){
          presultI[i]=new long double*[neq];
          for(int j=0;j<neq;j++){
              presultI[i][j]=new long double[neq];
           }
     }
/////////////////////////////Sweeping//////////////////////////////////

for(int m=-steps;m<=steps;m++)
{

   cout<<m<<endl;
   long double De=0;
   adj_detune(detune*m*(pow(-1,omp_get_thread_num()))*1.0/total_steps);
   long double period=10.878278481971048332082512612548/100*(100+De);
   long double peak=peakO*(100+De)/100;
   long double interval_1=FWHM*10,interval_2=period-interval_1;
   long double dt_1=interval_1/ninterval_1,dt_2=interval_2/ninterval_2;
   interval_2=period-interval_1;
   dt_2=interval_2/ninterval_2;

  for (int i=0;i<(pulse_average+1);i++){
        for (int k=0;k<neq;k++){
            for (int l=0;l<neq;l++){
                    presultI[i][k][l]=0;
                    presultR[i][k][l]=0;
                    }}}

    for (int i=0;i<neq;i++){ //initailizing
        for (int j=0;j<neq;j++){
            presultI[0][i][j]=y0I[i][j];
            presultR[0][i][j]=y0R[i][j];
            }
    }

    for(int n=0;n<2*neq;n++)
     for(int i=0;i<neq;i++)
         for(int j=0;j<2*neq;j++)
             for(int k=0;k<neq;k++){
                  if(n==j&&i==k){
                   Trans[n][i][j][k]=1;
                   Trans_AVE[n][i][j][k]=1;}
                 else{
                   Trans[n][i][j][k]=0;
                   Trans_AVE[n][i][j][k]=0;}
             }//initialing


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

 int k=0,flag=0;
 double diff=0;


while(flag<pulse_average){







     for(int a=0;a<neq;a++)
           for(int b=0;b<neq;b++){
           presultI[(k+1)%(pulse_average+1)][a][b]=0;
           presultR[(k+1)%(pulse_average+1)][a][b]=0;
             for(int c=0;c<neq;c++)
                 for(int d=0;d<neq;d++){
                  presultI[(k+1)%(pulse_average+1)][a][b]+=Trans[a+neq][b][c][d]*presultR[k%(pulse_average+1)][c][d]+Trans[a+neq][b][c+neq][d]*presultI[k%(pulse_average+1)][c][d];
                  presultR[(k+1)%(pulse_average+1)][a][b]+=Trans[a][b][c][d]*presultR[k%(pulse_average+1)][c][d]+Trans[a][b][c+neq][d]*presultI[k%(pulse_average+1)][c][d];
                 }}


    k+=1;

      if(k>pulse_average){
          diff=presultR[k%(pulse_average+1)][0][0]-presultR[(k-pulse_average)%(pulse_average+1)][0][0];
        if(abs(diff)<convergence)
         flag+=1;
        else
         flag=0;
      }

     if(k==(npulse-1))
       flag=pulse_average+1;

}

long double buffer=0;

             for(int c=0;c<neq;c++)
                 for(int d=0;d<neq;d++){
                  buffer+=Trans_AVE[0][0][c][d]*presultR[k%(pulse_average+1)][c][d]+Trans_AVE[0][0][c+neq][d]*presultI[k%(pulse_average+1)][c][d];
                  buffer+=Trans_AVE[1][1][c][d]*presultR[k%(pulse_average+1)][c][d]+Trans_AVE[1][1][c+neq][d]*presultI[k%(pulse_average+1)][c][d];
                 }

buffer=buffer/(ninterval_1+ninterval_2+1);

       file2<<setiosflags(ios::left)<<setw(30)<<detune*m*(pow(-1,omp_get_thread_num()))*1.0/total_steps/2/pi*1000;
       file2<<setiosflags(ios::left)<<setw(30)<<buffer;
       file2<<setiosflags(ios::left)<<setw(30)<<k;
       file2<<setiosflags(ios::left)<<setw(30)<<m<<endl;

}

///////////////////////////////End of Sweeping//////////////////////////////////



     for(int i=0;i<neq*2;i++)
         for(int j=0;j<neq;j++){
             for(int k=0;k<neq*2;k++){
               delete[] Trans[i][j][k];
               delete[] Trans_AVE[i][j][k];
       }}

     for(int i=0;i<neq*2;i++){
         for(int j=0;j<neq;j++){
               delete[] Trans[i][j];
               delete[] Trans_AVE[i][j];
       }}

     for(int i=0;i<neq*2;i++){
               delete[] Trans[i];
               delete[] Trans_AVE[i];
     }
      delete[] Trans;
      delete[] Trans_AVE;


      for(int i=0;i<pulse_average+1;i++)
         for(int j=0;j<neq;j++){
               delete[] presultI[i][j];
               delete[] presultR[i][j];
       }

        for(int i=0;i<pulse_average+1;i++){
               delete[] presultI[i];
               delete[] presultR[i];
       }

       delete[] presultI;
       delete[] presultR;

       for(int i=0;i<ninterval_1+ninterval_2+1;i++)
         for(int j=0;j<neq;j++){
               delete[] M[i][j];
       }

       for(int i=0;i<ninterval_1+ninterval_2+1;i++)
               delete[] M[i];

      delete[] Time;
      delete[] M;





  return 0;




}


