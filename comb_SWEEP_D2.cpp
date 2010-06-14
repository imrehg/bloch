#define GMM_USES_LAPACK
#include "comb.h"
doub phase=0;
int pulse_average=100;
const int npulse=10000000;
int ninterval_1=50,ninterval_2=500,ninterval_b=0;//npulse = number of pulse; interval_1 =steps in interval 1 ;ninterval_b= ninterval in a free iteration
doub period0=10.87827848197104833208;
doub frequency=0,peakO=0.02428235615859994948/2,FWHM=1; //about 100uW/cm2 A = 1ns
const int  neq=32,neq_gr=16,ninterval=npulse*(ninterval_1+ninterval_2); // neq= nuber of equations, nexp= terms of expansion, ninterval= iteration terms
int nexp=12;
vector<doub> r(neq);
//total decay constant
vector<doub> R(neq);
//relaxation rate
vector<doub> R_gr(neq);
//relazation rate of ground state
col_matrix< vector<doub> > Rc(neq,neq);
col_matrix< vector<doub> > A(neq,neq);
vector<doub> R_L(neq);
doub lasDe = 0;
vector<doub> EnergyDiff(neq-1);
doub LineWidth=0.7*2*pi;//0.0052227*2*pi;//


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

doub Gaussian(doub x,peroid,doub A)
{
    Return exp(-pow((x-0.5*period)/FWHM,2)/2);
}

doub ReRabi(doub x,doub period,doub peak)//脈衝包絡線函數(實部)，高斯函數*Re[e^{-i*phase}]
{
  doub value=0;

  value=exp(-pow((x-0.5*period)/FWHM,2)/2);

  return peak*value;
}

doub ImRabi(doub x,doub period,doub peak)//脈衝包絡線函數(虛部)，高斯函數*Im[e^{-i*phase}]
{
  doub value=0,time=0,factor=0;
  int i=0;
  if(x<0.5*period)
    value=exp(-pow(x/FWHM,2)/2)*sin(-i*phase);
  else
  {
  i=int ((x-0.5*period)/period+1);
  value=exp(-pow((x-i*period)/FWHM,2)/2)*sin(-i*phase);
  }

  return peak*value;
}


void fun_Matrix(col_matrix< vector<doub> > &Trans, col_matrix< vector<doub> >&H,col_matrix< vector<doub> >&D)//B imaginary part ; C real part; H rabi requence; w time
{

  for (int i=0;i<neq;i++){
     Trans(RealComp(i,i),RealComp(i,i))=-r[i]-R_gr[i];
     for (int j=0;j<neq;j++){
          if(i<j){
          Trans(RealComp(i,i),ImagComp(i,j))+=2*H(j,i);
          } else{
              if(i>j)
              Trans(RealComp(i,i),ImagComp(i,j))+=-2*H(j,i);
          }
          Trans(RealComp(i,i),RealComp(j,j))+=A(i,j);
      }
 }
 // Bloch eq for population part

 for(int i=0;i<neq_gr;i++){
  for(int j=0;j<neq_gr;j++){
      if(i!=j)
        Trans(RealComp(i+neq-neq_gr,i+neq-neq_gr),RealComp(j+neq-neq_gr,j+neq-neq_gr))+=R_gr[j+neq-neq_gr]/(neq_gr-1);
  }
 }

 //  Additional setup for Bloch eq of ground state ## This term should be combine to A matrix later##

 for (int j=1;j<neq;j++){
    for (int i=0;i<j;i++){
       Trans(RealComp(i,j),RealComp(i,j))=-(R[i]+R[j]+R_L[i]+R_gr[i]+R_gr[j])/2;
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
       Trans(ImagComp(i,j),ImagComp(i,j))=-(R[i]+R[j]+R_L[i]+R_gr[i]+R_gr[j])/2;
       Trans(ImagComp(i,j),RealComp(i,j))=-D(i,j);
      for (int l=0;l<neq;l++){
            Trans(ImagComp(i,j),RealComp(l,j))+=H(i,l);
            Trans(ImagComp(i,j),RealComp(i,l))+=-H(l,j);
      }
    }
  }    //  Bloch eq for off diagonal imaginary part

}
void solve_Martix_new(col_matrix< vector<doub> > &M, col_matrix<vector<doub> > &Trans,col_matrix<vector<doub> > &Trans_Blank, col_matrix< vector<doub> >&D,doub dt1,doub dt2)// solve(presultI,presultR,M,k)
{
  col_matrix< vector<doub> > Trans_I(neq*neq,neq*neq),Trans_C(neq*neq,neq*neq),Trans_E(neq*neq,neq*neq),Trans_2B(neq*neq,neq*neq),Trans_Ave_B(neq*neq,neq*neq);
  col_matrix< vector<doub> > Msub(neq,neq);
  col_matrix<vector<doub> > Trans_D(neq*neq,neq*neq),Trans_B(neq*neq,neq*neq);


/* Build up for a free iteration interval*/

      clear(Trans_B);
      clear(Trans_E);
      clear(Trans_I);
      clear(Trans_C);

    for(int i=0;i<neq*neq;i++)
             Trans_B(i,i)=1,Trans_I(i,i)=1;

      copy(sub_matrix(M,sub_interval(0,neq),sub_interval(0,neq)),Msub);
      fun_Matrix(Trans_E,Msub,D);

      for(int j=1;j<=nexp;j++){
          if((j%2)==1){
          mult(Trans_E,Trans_I,Trans_C);
          add(scaled(Trans_C,pow((dt2),j)/factorial(j)),Trans_B);
          }else{
            mult(Trans_E,Trans_C,Trans_I);
            add(scaled(Trans_I,pow((dt2),j)/factorial(j)),Trans_B);
          }
         }

        copy(Trans_B,Trans_2B);

for(int j=0;j<(log(ninterval_b/2)/log(2));j++){
   mult(Trans_B,Trans_2B,Trans_I);
   copy(Trans_I,Trans_B);
   copy(Trans_I,Trans_2B);
}

/* End of building*/


/*  bulid up iteration matrix for zero field part  */


      clear(Trans_B);
      clear(Trans_E);
      clear(Trans_I);
      clear(Trans_C);

    for(int i=0;i<neq*neq;i++)
             Trans_B(i,i)=1;

       for(int i=0;i<neq*neq;i++)
             Trans_I(i,i)=1;

      copy(sub_matrix(M,sub_interval(0,neq),sub_interval(0,neq)),Msub);
      fun_Matrix(Trans_E,Msub,D);

      for(int j=1;j<=nexp;j++){
          if((j%2)==1){
          mult(Trans_E,Trans_I,Trans_C);
          add(scaled(Trans_C,pow((dt2),j)/factorial(j)),Trans_B);
          }else{
            mult(Trans_E,Trans_C,Trans_I);
            add(scaled(Trans_I,pow((dt2),j)/factorial(j)),Trans_B);
          }
         }

        copy(Trans_B,Trans_2B);

for(int j=0;j<(log(ninterval_2/2)/log(2));j++){
   mult(Trans_B,Trans_2B,Trans_I);
   copy(Trans_I,Trans_B);
   copy(Trans_I,Trans_2B);
}

/*  End of building  */


/*  bulid up iteration matrix for pulse field part  */

  for(int t=ninterval_2/2;t<(ninterval_2/2+ninterval_1);t++){
  cout<<t<<endl;
      clear(Trans_E);
      clear(Trans_I);
      clear(Trans_C);
      clear(Trans_B);

       for(int i=0;i<neq*neq;i++)
             Trans_I(i,i)=1,Trans_B(i,i)=1;

      copy(sub_matrix(M,sub_interval(0,neq),sub_interval((t-1)*neq,neq)),Msub);
      fun_Matrix(Trans_E,Msub,D);

      for(int j=1;j<=nexp;j++){
          if((j%2)==1){
          mult(Trans_E,Trans_I,Trans_C);
          add(scaled(Trans_C,pow((dt1),j)/factorial(j)),Trans_B);
          }else{
            mult(Trans_E,Trans_C,Trans_I);
            add(scaled(Trans_I,pow((dt1),j)/factorial(j)),Trans_B);
          }
         }

     if((t%2)==((ninterval_2/2)%2)){
        mult(Trans_B,Trans,Trans_D);
     }else{
        mult(Trans_B,Trans_D,Trans);
     }

  }

  if((ninterval_1%2)==1)
      copy(Trans_D,Trans);

/*  End of building  */

   mult(Trans,Trans_2B,Trans_D);
   mult(Trans_2B,Trans_D,Trans);

}


void solve_Martix(col_matrix< vector<doub> >&M, col_matrix<vector<doub> > &Trans, col_matrix< vector<doub> >&Trans_Ave, col_matrix< vector<doub> >&D,doub dt1,doub dt2)// solve(presultI,presultR,M,k)
{
  col_matrix< vector<doub> > Trans_I(neq*neq,neq*neq),Trans_C(neq*neq,neq*neq),Trans_E(neq*neq,neq*neq),Trans_2B(neq*neq,neq*neq),Trans_Ave_B(neq*neq,neq*neq);
  col_matrix< vector<doub> > Msub(neq,neq);
  col_matrix<vector<doub> > Trans_D(neq*neq,neq*neq),Trans_B(neq*neq,neq*neq);

    for(int i=0;i<neq*neq;i++)
             Trans_B(i,i)=1;


  for(int t=1;t<(ninterval_2/2+ninterval_1+2);t++){

      clear(Trans_E);
      clear(Trans_I);
      clear(Trans_C);
//      cout<<t<<endl;

       for(int i=0;i<neq*neq;i++)
             Trans_I(i,i)=1;

      copy(sub_matrix(M,sub_interval(0,neq),sub_interval((t-1)*neq,neq)),Msub);
      fun_Matrix(Trans_E,Msub,D);
//time_t start=clock();

    if(t<=ninterval_2/2||t>(ninterval_2/2+ninterval_1)){
     if(t==1){
      for(int j=1;j<=nexp;j++){
          if((j%2)==1){
          mult(Trans_E,Trans_I,Trans_C);
          add(scaled(Trans_C,pow((dt2),j)/factorial(j)),Trans_B);
          }else{
            mult(Trans_E,Trans_C,Trans_I);
            add(scaled(Trans_I,pow((dt2),j)/factorial(j)),Trans_B);
          }
         }
        copy(Trans_B,Trans_2B);
      }else{
         if(t==(ninterval_2/2+ninterval_1+1))
            copy(Trans_2B,Trans_B);
      }
    }else{
      clear(Trans_B);
     for(int i=0;i<neq*neq;i++)
             Trans_B(i,i)=1;
      for(int j=1;j<=nexp;j++){
          if((j%2)==1){
          mult(Trans_E,Trans_I,Trans_C);
          add(scaled(Trans_C,pow((dt1),j)/factorial(j)),Trans_B);
          }else{
            mult(Trans_E,Trans_C,Trans_I);
            add(scaled(Trans_I,pow((dt1),j)/factorial(j)),Trans_B);
          }
         }
    }

//cout<<(clock()-start)*1.0/CLOCKS_PER_SEC<<endl;


     if(t%2==1){
        mult(Trans_B,Trans,Trans_D);
        if(t!=(ninterval_2/2+ninterval_1+1))
         add(Trans_D,Trans_Ave);
     }else{
        mult(Trans_B,Trans_D,Trans);
        if(t!=(ninterval_2/2+ninterval_1+1))
         add(Trans,Trans_Ave);
     }

     if(t==ninterval_2/2){
           copy(Trans,Trans_2B);
           copy(Trans_Ave,Trans_Ave_B);
     }

//        copy(Trans_D,Trans);
//        cout<<"nnz="<<nnz(Trans_E)<<endl;

  }

  mult(Trans_Ave_B,Trans,Trans_2B);
  add(Trans_2B,Trans_Ave);
  copy(Trans_D,Trans);

}


int D1_coef (int L,int F,int mf){
     if(F==3)
      return 32-(L*16+(mf+4));
     if(F==4)
      return 32-(L*16+7+(mf+5));
 }


int sweep(doub g2,doub LineW,int steps,int total_steps,doub PeakPower,doub convergence,doub convergence_threshold,int conS,int expN,int n1, int n2,int Msteps,doub detune,doub A_Factor)
{

  clear(A);
  clear(R);
  clear(r);
  clear(Rc);
  clear(R_L);
  clear(EnergyDiff);
  doub phase=0;
  ninterval_1 =n1-(n1%4);
  ninterval_2 = pow(2,int(log(n2)/log(2)));
  fstream file1,file2;//file1:紀錄輸入的參數。file2://紀錄計算結果
  FWHM = A_Factor;//0.00175
  peakO = sqrt(PeakPower/100/A_Factor)*0.02428235615859994948/2;
  nexp=expN;
  stringstream strstream,strstream2;
  string filename,filename2;
  strstream<<"./Data/comb_g2_"<<g2/2/pi<<"_LW_"<<LineW/2/pi<<"GHz_"<<PeakPower<<"uWcm2_"<<convergence<<"_conS_"<<conS<<"_O="<<nexp<<"_N1_"<<n1<<"_N2_"<<n2<<"_D_"<<detune/2/pi*1000<<"MHz_A_"<<A_Factor<<"_S.txt";
  strstream>>filename;
  cout<<filename.c_str()<<endl;
  file2.open(filename.c_str(),ios::out | ios::trunc);
  strstream2<<"./Data/comb_FS_"<<"g2_"<<g2/2/pi<<"_LW_"<<LineW/2/pi<<"GHz_"<<PeakPower<<"uWcm2_"<<convergence<<"_conS_"<<conS<<"_O="<<nexp<<"_N1_"<<n1<<"_N2_"<<n2<<"_D_"<<detune/2/pi*1000<<"MHz_A_"<<A_Factor<<"_S.txt";
  strstream2>>filename2;
  file1.open(filename2.c_str(), ios::out | ios::trunc);
  file2.precision(15);
  col_matrix< vector<doub> > y0I(neq,neq);
  //initial condition
  col_matrix< vector<doub> > y0R(neq,neq);
  //initial condiion
  col_matrix< vector<doub> > EnerDet(neq,neq);
  Atom atom;
  LineWidth=LineW;
  doub gamma2=g2/1E9;

   EnergyDiff[8]=0.2012871*2*pi;
   EnergyDiff[24]=9.192631*2*pi;

    for (int i=0; i<neq; i++)
      for (int j=i+1; j<neq; j++)
        for (int k=i; k<j; k++){
           EnerDet(i,j)+=EnergyDiff[k];
           EnerDet(j,i)-=EnergyDiff[k];
        }
  //initialzing the energy level difference for D2 line

    for(int j=3; j<5;j++)
         for(int k=-j;k<j+1;k++)
             for(int m=3; m<5;m++)
                for(int n=-m;n<m+1;n++){
                     EnerDet(D1_coef(0,j,k),D1_coef(1,m,n))+=detune;
                     EnerDet(D1_coef(1,m,n),D1_coef(0,j,k))-=detune;
                }
  //set detuning of laser field positve represent blue detuning.

      for(int j=3; j<5;j++)
         for(int k=-j;k<j+1;k++)
             for(int m=3; m<5;m++)
                for(int n=-m;n<m+1;n++)
                  for(int q=-1;q<2;q++)
                     A(D1_coef(0,j,k),D1_coef(1,m,n))+=pow(atom.coef(q,1,0,m,j,n,k,1.5,0.5,3.5),2)*LineWidth;

  //initialzing the A coefficients

 int Matrix_Step = pow(2,(Msteps-1));
 pulse_average=(conS/Matrix_Step+1);

 int num_thread = 8;

#pragma omp num_threads(num_thread)
#pragma omp parallel for
for(int thread=0;thread<num_thread;thread++)
{
   doub *Time= new doub[ninterval_1+ninterval_2+1];
   col_matrix< vector<doub> > M(neq,(ninterval_1+ninterval_2+1)*neq);//Rabifrequence*2
   int ninterval_m = (npulse-1);


for(int m= int(omp_get_thread_num()/2)*steps;m<=steps*(int(omp_get_thread_num()/2)+1);m++)
{
   cout<<"T_"<<omp_get_thread_num()<<endl;
   cout<<m<<endl;
   doub De=1.0*m*(pow(-1,(omp_get_thread_num()%2+1)))*1.0/total_steps;
   doub period=period0/100*(100+De);
   doub peak=peakO*sqrt((100.0+De)/100);
   doub interval_1=FWHM*5,interval_2=period-interval_1;
   doub dt_1=interval_1/ninterval_1,dt_2=interval_2/ninterval_2;
//   ninterval_b= int((phase_shift/2/pi*period)/dt_2);
   interval_2=period-interval_1;
   dt_2=interval_2/ninterval_2;

    col_matrix< vector<doub> > Result(neq*neq,pulse_average+1);
    col_matrix<vector<doub> > Trans(neq*neq,neq*neq);
    col_matrix<vector<doub> > Trans_Blank(neq*neq,neq*neq);
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


 for(int i=0;i<(neq-neq_gr);i++){
   r[i]=LineWidth;
   R[i]=LineWidth;
 }

for(int i=0;i<neq_gr;i++)
   R_gr[i+neq-neq_gr]=gamma2;

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

if( k>=ninterval_2/2 && k<(ninterval_2/2+ninterval_1) ){
  for(int j=3; j<5;j++)
    for(int t=-j;t<j+1;t++)
	  for(int m=3; m<5;m++)
	    for(int n=-m;n<m+1;n++){
	      M(D1_coef(1,j,t),k*neq+D1_coef(0,m,n))=(atom.coef(1,1,0,j,m,t,n,1.5,0.5,3.5)+atom.coef(1,1,0,j,m,t,n,1.5,0.5,3.5))/2*ReRabi(Time[k],period,peak);
	      M(D1_coef(0,m,n),k*neq+D1_coef(1,j,t))=M(D1_coef(1,j,t),k*neq+D1_coef(0,m,n));}
}
//initailizing for M matrices
//The reason to set the matrix this way(the second equation) is that actually calulated transition would be pure imaginary, but we set is to real(multiply a phase).
//If we directly set M and run through the parameter, we will get a extra munus sign in the symmetric terms, which can't be used in the formalism applied in fun_matrix.
           }


solve_Martix(M,Trans,Trans_AVE,EnerDet,dt_1,dt_2);

dense_matrix < doub > Trans_1(neq*neq,neq*neq),Trans_2(neq*neq,neq*neq),Trans_3(neq*neq,neq*neq);

copy(Trans,Trans_1);
copy(Trans,Trans_2);

for(int j=0;j<(Msteps-1);j++){
   mult(Trans_1,Trans_2,Trans_3);
   copy(Trans_3,Trans_1);
   copy(Trans_3,Trans_2);
   cout<<"S_"<<j<<endl;
}

copy(Trans_1,Trans);



cout<<"end of solve"<<endl;
int k=0,flag=0;
doub diff=0;
//
//if(m==0){

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
	      if(Result(RealComp(c,d),(k%(pulse_average+1)))<=convergence_threshold)
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
//     for(int c=0;c<16;c++)
//      data_sum+=Result(c,k%(pulse_average+1));

 for(int l=-3;l<4;l++)
  for(int n=-4;n<5;n++){
    data_sum+=(pow(Result(RealComp(D1_coef(0,3,l),D1_coef(0,4,n)),k%(pulse_average+1)),2)+pow(Result(ImagComp(D1_coef(0,3,l),D1_coef(0,4,n)),k%(pulse_average+1)),2));
    }
    data_sum=data_sum/63;
    if(m==0)
     file1<<k*period*Matrix_Step<<"\t"<<data_sum<<endl;
   }

  }

cout<<"npulse="<<k<<endl;

doub buffer=0,buffer2=0,bufferP=0,bufferC=0;

 for (int j=-3;j<4;j++)
   for(int d=0;d<neq*neq;d++){
     buffer+=Trans_AVE(RealComp(D1_coef(1,3,j),D1_coef(1,3,j)),d)*Result(d,k%(pulse_average+1));
   }
 for (int j=-3;j<5;j++)
   for(int d=0;d<neq*neq;d++){
     buffer2+=Trans_AVE(RealComp(D1_coef(1,4,j),D1_coef(1,4,j)),d)*Result(d,k%(pulse_average+1));
   }

 for(int c=0;c<neq;c++)
   bufferP+=Result(c,k%(pulse_average+1));

 for(int l=-3;l<4;l++)
 for(int n=-4;n<5;n++){
   bufferC+=(pow(Result(RealComp(D1_coef(0,3,l),D1_coef(0,4,n)),k%(pulse_average+1)),2)+pow(Result(ImagComp(D1_coef(0,3,l),D1_coef(0,4,n)),k%(pulse_average+1)),2));
   }
 bufferC=bufferC/63;

 buffer=buffer/(ninterval_1+ninterval_2+1);
 buffer2=buffer2/(ninterval_1+ninterval_2+1);

 file2<<1/period<<"\t";
 file2<<buffer+buffer2<<"\t";
 file2<<buffer<<"\t";
 file2<<buffer2<<"\t";
 file2<<bufferC<<"\t";
 file2<<bufferP<<"\t";
 file2<<k*Matrix_Step*period<<"\t";
 file2<<m<<endl;
}

///////////////////////////////End of Sweeping//////////////////////////////////

 delete[] Time;

}
 return 0;

}
