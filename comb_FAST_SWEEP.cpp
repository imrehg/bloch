#include <iostream>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <time.h>
using namespace std;
const long double pi=3.141592654;
const int npulse=5000,ninterval_1=50,ninterval_2=500;//npulse = number of pulse; interval_1 =steps in interval 1 ..
long double De=0,period=10.878278481971048332082512612548/100*(100+De);
long double frequency=0,peakO=37.6834589/2,peak=37.6834589/2,FWHM=0.001,interval_1=FWHM*10,interval_2=period-interval_1; //frequency:載波角頻率。peroid：脈衝周期。FWHM：脈衝半高寬。peak：拉比頻率最大值
long double dt_1=interval_1/ninterval_1,dt_2=interval_2/ninterval_2;
const int neq=4,nexp=12,ninterval=npulse*(ninterval_1+ninterval_2); // neq= nuber of equations, nexp= terms of expansion, ninterval= iteration terms
long double r[neq]={0.0052227*2*pi,0.0052227*2*pi,0,0};//total decay constant
long double R[neq]={0.0052227*2*pi,0.0052227*2*pi,0,0};//relaxation rate
long double A[neq][neq]={{0,0,0,0},{0,0,0,0},{0.0052227*2*pi/2,0.0052227*2*pi/2,0,0},{0.0052227*2*pi/2,0.0052227*2*pi/2,0,0}};//Einstein A coefficient
long double R_L[neq]={0,0,0,0};//laser line width
long double d[neq][neq]={{0,0,0,-0.20124*2*pi-9.192631*2*pi},
                         {0,0,0,-9.192631*2*pi},
                         {0,0,0,-9.192631*2*pi},
                         {+9.192631*2*pi,9.192631*2*pi,9.192631*2*pi,0}};//laser */
long double y0I[neq][neq]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};//initial condition
long double y0R[neq][neq]={{0,0,0,0},{0,0,0,0},{0,0,0.5,0},{0,0,0,0.5}};//initial condiion

void fun(long double ***,long double ***,long double **,int );//聯立方程式
void solve(long double ***,long double ***,long double ***,long double ***,long double ***,long double*,int);//演算法
long double ReRabi(long double );//脈衝包絡線函數(實部)
long double ImRabi(long double );//脈衝包絡線函數(虛部)
int factorial (int);
void fun_Matrix(long double *****,long double **);
void solve_Martix(long double ***,long double ****,long double ****,long double *);
void Matrix_Multiply(long double ****,long double ****);
double phase=0;


int factorial (int num)
{
 if (num==1)
  return 1;
 return factorial(num-1)*num; // recursive call
}

long double ReRabi(long double x)//脈衝包絡線函數(實部)，高斯函數*Re[e^{-i*phase}]
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

long double ImRabi(long double x)//脈衝包絡線函數(虛部)，高斯函數*Im[e^{-i*phase}]
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


      for(int a=0;a<neq;a++)
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

                 }


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
       Trans[i][j][i][j]=-(R[i]+R[j]+R_L[i])/2;
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
       Trans[i+neq][j][i+neq][j]=-(R[i]+R[j]+R_L[i])/2;
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

int main()
{
  int start=clock();
  long double phase=0;
  fstream file1,file2;//file1:紀錄輸入的參數。file2://紀錄計算結果
  file1.open("input.txt", ios::out | ios::trunc);
  file2.open("data.txt", ios::out | ios::trunc);
  file2.precision(10);

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
    long double ***presultR= new long double**[npulse];//所有時間點的數值存於此指標(real)
      for(int i=0;i<npulse;i++){
          presultR[i]=new long double*[neq];
          for(int j=0;j<neq;j++){
              presultR[i][j]=new long double[neq];
           }
     }
   long double ***presultI= new long double**[npulse];//所有時間點的數值存於此指標(imaginary) presultI[][o][o]為時間參數
      for(int i=0;i<npulse;i++){
          presultI[i]=new long double*[neq];
          for(int j=0;j<neq;j++){
              presultI[i][j]=new long double[neq];
           }
     }






/////////////////////////////Sweeping//////////////////////////////////

for(int m=-50;m<=50;m++){

cout<<m<<endl;

   De=m/200.0;
   period=10.878278481971048332082512612548/100*(100+De);
   peak=peakO*(100+De)/100;
   interval_2=period-interval_1;
   dt_2=interval_2/ninterval_2;

  for (int i=0;i<npulse;i++){
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
              M[k][0][2]=-2.269*ReRabi(buffer);
              M[k][0][3]=1.906*ReRabi(buffer);
              M[k][1][0]=0;
              M[k][1][1]=0;
              M[k][1][2]=2.269*ReRabi(buffer);
              M[k][1][3]=1.906*ReRabi(buffer);
              M[k][2][0]=-2.269*ReRabi(buffer);
              M[k][2][1]=2.269*ReRabi(buffer);
              M[k][2][2]=0;
              M[k][2][3]=0;
              M[k][3][0]=1.906*ReRabi(buffer);
              M[k][3][1]=1.906*ReRabi(buffer);
              M[k][3][2]=0;
              M[k][3][3]=0;
           }

solve_Martix(M,Trans,Trans_AVE,Time);

 int k=0,flag=0;
 double diff=0;
 int pulse_average=100;

while(flag<pulse_average){

     for(int a=0;a<neq;a++)
           for(int b=0;b<neq;b++)
             for(int c=0;c<neq;c++)
                 for(int d=0;d<neq;d++){
                  presultI[k+1][a][b]+=Trans[a+neq][b][c][d]*presultR[k][c][d]+Trans[a+neq][b][c+neq][d]*presultI[k][c][d];
                  presultR[k+1][a][b]+=Trans[a][b][c][d]*presultR[k][c][d]+Trans[a][b][c+neq][d]*presultI[k][c][d];
                 }

    k+=1;

      if(k>pulse_average){
          diff=presultR[k][0][0]-presultR[k-pulse_average][0][0];
        if(abs(diff)<0.000001)
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
                  buffer+=Trans_AVE[0][0][c][d]*presultR[k][c][d]+Trans_AVE[0][0][c+neq][d]*presultI[k][c][d];
                  buffer+=Trans_AVE[1][1][c][d]*presultR[k][c][d]+Trans_AVE[1][1][c+neq][d]*presultI[k][c][d];
                 }

buffer=buffer/(ninterval_1+ninterval_2+1);

       file2<<setiosflags(ios::left)<<setw(20)<<1/period;
       file2<<setiosflags(ios::left)<<setw(20)<<buffer;
       file2<<setiosflags(ios::left)<<setw(20)<<k<<endl;

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


  cout<<"time spent:"<<(clock()-start)/1000.0<<"sec";

  return 0;




}


