#define GMM_USES_LAPACK
#include "comb.h"
doub phase=0;
int pulse_average=100;
const doub pi=3.14159265358979323846264338327950288419716939937511;
const int npulse=100000;
int ninterval_1=50,ninterval_2=500;//npulse = number of pulse; interval_1 =steps in interval 1 ..
doub period0=10.87827848197104833208;
doub frequency=0,peakO=1.34163815218652164669542605053/2,FWHM=0.0007; //about 150uW/cm2 about 1ps
const int  neq=32,ninterval=npulse*(ninterval_1+ninterval_2); // neq= nuber of equations, nexp= terms of expansion, ninterval= iteration terms
int nexp=12;
vector<doub> r(neq);
//total decay constant
vector<doub> R(neq);
//relaxation rate
col_matrix< vector<doub> > Rc(neq,neq);
col_matrix< vector<doub> > A(neq,neq);
vector<doub> R_L(neq);
doub lasDe = 0;
vector<doub> EnergyDiff(neq-1);


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

doub ReRabi(doub x,doub period,doub peak)//脈衝包絡線函數(實部)，高斯函數*Re[e^{-i*phase}]
{
  doub value=0;

  value=exp(-pow((x-0.5*period)/FWHM,2));

  return peak*value;
}

doub ImRabi(doub x,doub period,doub peak)//脈衝包絡線函數(虛部)，高斯函數*Im[e^{-i*phase}]
{
  doub value=0,time=0,factor=0;
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


void fun_Matrix(col_matrix< vector<doub> > &Trans, col_matrix< vector<doub> >&H,col_matrix< vector<doub> >&D)//B imaginary part ; C real part; H rabi requence; w time
{

  for (int i=0;i<neq;i++){
     Trans(RealComp(i,i),RealComp(i,i))=-r[i];
     for (int j=0;j<neq;j++){
          if(i<j){
          Trans(RealComp(i,i),ImagComp(i,j))+=2*H(j,i);
          } else{
              if(i>j)
              Trans(RealComp(i,i),ImagComp(i,j))+=-2*H(j,i);
          }
          Trans(RealComp(i,i),RealComp(j,j))+=A(i,j);
      }
 } // Bloch eq for population part

 for (int j=1;j<neq;j++){
    for (int i=0;i<j;i++){
       Trans(RealComp(i,j),RealComp(i,j))=-(R[i]+R[j]+R_L[i])/2-Rc(i,j);
       Trans(RealComp(i,j),ImagComp(i,j))=D(i,j);
       for (int l=0;l<neq;l++){
           if(l<j){
             Trans(RealComp(i,j),ImagComp(l,j))+=-H(i,l);
           } else{
              if(l>j)
                Trans(RealComp(i,j),ImagComp(l,j))+=H(i,l);
           }
           if(i<l){
             Trans(RealComp(i,j),ImagComp(i,l))+=H(l,j);
           } else{
              if(l>i)
              Trans(RealComp(i,j),ImagComp(i,l))+=-H(l,j);
           }
         }
   }
 }   // Bloch eq for off diagonal real part

 for (int j=1;j<neq;j++){
   for (int i=0;i<j;i++){
       Trans(ImagComp(i,j),ImagComp(i,j))=-(R[i]+R[j]+R_L[i])/2-Rc(i,j);
       Trans(ImagComp(i,j),RealComp(i,j))=-D(i,j);
      for (int l=0;l<neq;l++){
            Trans(ImagComp(i,j),RealComp(l,j))+=H(i,l);
            Trans(ImagComp(i,j),RealComp(i,l))+=-H(l,j);
      }
    }
  }    //  Bloch eq for off diagonal imaginary part

}

void solve_Martix(col_matrix< vector<doub> >&M, col_matrix<vector<doub> > &Trans, col_matrix< vector<doub> >&Trans_Ave,doub *T, col_matrix< vector<doub> >&D)// solve(presultI,presultR,M,k)
{
  col_matrix< vector<doub> > Trans_I(neq*neq,neq*neq),Trans_C(neq*neq,neq*neq),Trans_E(neq*neq,neq*neq);
  col_matrix< vector<doub> > Msub(neq,neq);
  col_matrix<vector<doub> > Trans_D(neq*neq,neq*neq),Trans_B(neq*neq,neq*neq);

  for(int t=1;t<(ninterval_1+ninterval_2)+1;t++){


      clear(Trans_B);
      clear(Trans_D);
      clear(Trans_E);
      clear(Trans_I);
      clear(Trans_C);
      cout<<t<<endl;

       for(int i=0;i<neq*neq;i++){
             Trans_I(i,i)=1;
             Trans_B(i,i)=1;
       }

      copy(sub_matrix(M,sub_interval(0,neq),sub_interval((t-1)*neq,neq)),Msub);
      fun_Matrix(Trans_E,Msub,D);
//    copy(Trans_E,Trans_E_R);

      for(int j=1;j<=nexp;j++){
          mult(Trans_E,Trans_I,Trans_C);
          copy(Trans_C,Trans_I);
          add(scaled(Trans_I,pow((T[t]-T[t-1]),j)/factorial(j)),Trans_B);
//          cout<<"nnzE="<<nnz(Trans_E)<<endl;
//          cout<<"nnzI="<<nnz(Trans_I)<<endl;
//          cout<<"nnzC="<<nnz(Trans_C)<<endl;
      }
//        time_t start=clock();
        add(Trans_B,Trans_Ave);
        mult(Trans_B,Trans,Trans_D);
        copy(Trans_D,Trans);
//        cout<<"nnz="<<nnz(Trans_E)<<endl;
//        cout<<(clock()-start)*1.0/CLOCKS_PER_SEC<<endl;

  }
}


int D1_coef (int L,int F,int mf){
     if(F==3)
      return 32-(L*16+(mf+4));
     if(F==4)
      return 32-(L*16+7+(mf+5));
 }


int sweep(int steps,int total_steps,doub PeakPower,doub convergence,int conS,int expN,int n1, int n2,doub detune)
{

  clear(A);
  clear(R);
  clear(r);
  clear(Rc);
  clear(R_L);
  clear(EnergyDiff);
  doub phase=0;
  ninterval_1 =n1;
  ninterval_2 =n2;
  fstream file1,file2;//file1:紀錄輸入的參數。file2://紀錄計算結果
  peakO = PeakPower/150*1.34163815218652164669542605053/2;
  nexp=expN;
  stringstream strstream,strstream2;
  string filename,filename2;
  strstream<<PeakPower<<"uWcm2_"<<convergence<<"_conS_"<<conS<<"_O="<<nexp<<"_N1_"<<n1<<"_N2_"<<n2<<"_D_"<<detune/2/pi*1000<<"MHz_S.txt";
  strstream>>filename;
  cout<<filename.c_str()<<endl;
  file2.open(filename.c_str(),ios::out | ios::trunc);
  strstream2<<"FS_"<<PeakPower<<"uWcm2_"<<convergence<<"_conS_"<<conS<<"_O="<<nexp<<"_N1_"<<n1<<"_N2_"<<n2<<"_D_"<<detune/2/pi*1000<<"MHz_S.txt";
  strstream2>>filename2;
  file1.open(filename2.c_str(), ios::out | ios::trunc);
  file2.precision(15);
  pulse_average=conS;
  col_matrix< vector<doub> > y0I(neq,neq);
  //initial condition
  col_matrix< vector<doub> > y0R(neq,neq);
  //initial condiion
  col_matrix< vector<doub> > EnerDet(neq,neq);
  Atom atom;

   EnergyDiff[8]=0.2012871*2*pi;
   EnergyDiff[24]=9.192631*2*pi;

    for (int i=0; i<neq; i++)
      for (int j=i+1; j<neq; j++)
        for (int k=i; k<j; k++){
           EnerDet(i,j)-=EnergyDiff[k];
           EnerDet(j,i)+=EnergyDiff[k];
        }
  //initialzing the energy level difference for D2 line


      for(int j=3; j<5;j++)
         for(int k=-j;k<j+1;k++)
             for(int m=3; m<5;m++)
                for(int n=-m;n<m+1;n++)
                  for(int q=-1;q<2;q++)
                     A(D1_coef(0,j,k),D1_coef(1,m,n))+=pow(atom.coef(q,0,1,j,m,k,n,0.5,0.5,3.5),2)*0.0052227*2*pi;

  //initialzing the A coefficients



#pragma omp num_threads(2)
#pragma omp parallel for
for(int thread=0;thread<2;thread++)
{
   doub *Time= new doub[ninterval_1+ninterval_2+1];
   col_matrix< vector<doub> > M(neq,(ninterval_1+ninterval_2+1)*neq);//Rabifrequence*2
   int ninterval_m = (npulse-1);


for(int m=0;m<=steps;m++)
{

   cout<<m<<endl;
   doub De=1.0*m*(pow(-1,omp_get_thread_num()))*1.0/total_steps;
   doub period=period0/100*(100+De);
   doub peak=peakO*(100.0+De)/100;
   doub interval_1=FWHM*5,interval_2=period-interval_1;
   doub dt_1=interval_1/ninterval_1,dt_2=interval_2/ninterval_2;
   interval_2=period-interval_1;
   dt_2=interval_2/ninterval_2;

    col_matrix< vector<doub> > Result(neq*neq,pulse_average+1);
    col_matrix<vector<doub> > Trans(neq*neq,neq*neq);
    // IMPORTANT!!!!!!!
    // Remember to change it back to dense_matrix<doub> when dealing with left plus right circular polarization, since the matrix would be larger than.
    // IMPORTANT!!!!!!!!
    col_matrix< vector<doub> > Trans_AVE(neq*neq,neq*neq);



   for(int i=16; i<32;i++)
     y0R(i,i)=1.0/16;

    for (int i=0;i<neq;i++){
    for (int j=0;j<neq;j++){

            Result(RealComp(i,j),0)=y0R(i,j);
            if(i!=j) Result(ImagComp(i,j),0)=y0I(i,j);

            }
    }
 //initailizing for density matrices


 for(int i=0;i<16;i++){
   r[i]=0.0052227*2*pi;
   R[i]=0.0052227*2*pi;
 }

 //initailizing for relaxation rate

    for(int i=0;i<neq*neq;i++){
        Trans(i,i)=1.0;
        Trans_AVE(i,i)=1.0;
    }


    for(int k=0;k<(ninterval_1+ninterval_2+1);k++){

      doub buffer=0;

      if( k>=ninterval_2/2 && k<(ninterval_2/2+ninterval_1) )
	buffer=dt_1;
      else
	buffer=dt_2;

      if(k==0)
	Time[k]=dt_2;
      else
	Time[k]=Time[k-1]+buffer;



      for(int j=3; j<5;j++)
        for(int t=-j;t<j+1;t++)
	  for(int m=3; m<5;m++)
	    for(int n=-m;n<m+1;n++){
	      M(D1_coef(1,j,t),k*neq+D1_coef(0,m,n))=(atom.coef(+0,1,0,j,m,t,n,0.5,0.5,3.5)+atom.coef(-0,1,0,j,m,t,n,0.5,0.5,3.5))/2*ReRabi(Time[k],period,peak);
	      M(D1_coef(0,m,n),k*neq+D1_coef(1,j,t))=M(D1_coef(1,j,t),k*neq+D1_coef(0,m,n));
               }
//initailizing for M matrices
//The reason to set the matrix this way(the second equation) is that actually calulated transition would be pure imaginary, but we set is to real(multiply a phase).
//If we directly set M and run through the parameter, we will get a extra munus sign in the symmetric terms, which can't be used in the formalism applied in fun_matrix.
           }


solve_Martix(M,Trans,Trans_AVE,Time,EnerDet);


cout<<"end of solve"<<endl;
int k=0,flag=0;
doub diff=0;

if(m==0){
  while(flag<pulse_average){

    for(int a=0;a<neq;a++)
      for(int b=0;b<neq;b++){
	Result(RealComp(a,b),(k+1)%(pulse_average+1))=0;
	Result(ImagComp(a,b),(k+1)%(pulse_average+1))=0;}

    mult(Trans,mat_col(Result,(k)%(pulse_average+1)),mat_col(Result,(k+1)%(pulse_average+1)));
//    cout<<Trans<<endl;
    k+=1;

    if(k>pulse_average){
      diff=0;

  for(int c=0;c<neq;c++){
	  for(int d=0;d<neq;d++){
	      if(Result(RealComp(c,d),(k%(pulse_average+1)))==0)
	         diff+=1;
	       else if(abs(1.0-Result(RealComp(c,d),(k%(pulse_average+1)))/(Result(RealComp(c,d),(k-pulse_average)%(pulse_average+1))))<convergence)
             diff+=1;
	  }
    }
          cout<<"diff="<<diff<<endl;
      if(diff==neq*neq)
	flag+=1;
      else
	flag=0;
    }
    if(k==ninterval_m)
      flag=pulse_average+1;
    if(omp_get_thread_num()==0){
      doub data_sum=0;
     for(int c=0;c<16;c++)
      data_sum+=Result(c,k%(pulse_average+1));
     file1<<setiosflags(ios::left)<<setw(30)<<k*period<<setiosflags(ios::left)<<setw(30)<<data_sum<<endl;
   }

  }

  ninterval_m = k;

}else{

  while(flag<pulse_average){

    for(int a=0;a<neq;a++)
      for(int b=0;b<neq;b++){
	Result(RealComp(a,b),(k+1)%(pulse_average+1))=0;
	Result(ImagComp(a,b),(k+1)%(pulse_average+1))=0;}

    mult(Trans,mat_col(Result,(k)%(pulse_average+1)),mat_col(Result,(k+1)%(pulse_average+1)));
    k+=1;

    if(k==ninterval_m)
      flag=pulse_average+1;
  }

}

cout<<"npulse="<<k<<endl;

doub buffer=0,buffer2=0,bufferC=0;

 for (int j=0;j<16;j++)
   for(int d=0;d<neq*neq;d++){
     buffer+=Trans_AVE(j,d)*Result(d,k%(pulse_average+1));
     //                  buffer2+=Trans_AVE(0,d)*Result(d,k%(pulse_average+1));
   }
 for(int c=0;c<neq;c++)
   bufferC+=Result(c,k%(pulse_average+1));

 buffer=buffer/(ninterval_1+ninterval_2+1);
 buffer2=buffer2/(ninterval_1+ninterval_2+1);

 file2<<setiosflags(ios::left)<<setw(30)<<1/period;
 file2<<setiosflags(ios::left)<<setw(30)<<buffer;
 //       file2<<setiosflags(ios::left)<<setw(30)<<buffer2;
 file2<<setiosflags(ios::left)<<setw(30)<<bufferC;
 file2<<setiosflags(ios::left)<<setw(30)<<k;
 file2<<setiosflags(ios::left)<<setw(30)<<m<<endl;


}

///////////////////////////////End of Sweeping//////////////////////////////////



 delete[] Time;



}


 return 0;




}


