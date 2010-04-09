#define GMM_USES_LAPACK
#include "CW.h"
doub phase=0;
int pulse_average=100;
const int npulse=100000;
int ninterval=50;//npulse = number of pulse; interval_1 =steps in interval 1 ..
doub dt;
doub period0=10.87827848197104833208;
doub frequency=0,peakO=0.00989116604990337187/2; //about 100uW/cm2 about 5ps  peak0=1.34163815218652164669542605053 for 2ps
const int  neq=32,neq_gr=16; // neq= nuber of equations, nexp= terms of expansion, ninterval= iteration terms
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
doub LineWidth=0.0052227*2*pi;//0.7*2*pi;//0.0052227*2*pi



int factorial (int num)
{
 if (num==1)
  return 1;
 return factorial(num-1)*num; // recursive call
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

void solve_Martix(col_matrix< vector<doub> >&M, col_matrix<vector<doub> > &Trans, col_matrix< vector<doub> >&D,doub dt)// solve(presultI,presultR,M,k)
{
  col_matrix< vector<doub> > Trans_I(neq*neq,neq*neq),Trans_E(neq*neq,neq*neq),Trans_Ave_B(neq*neq,neq*neq);
  col_matrix<vector<doub> > Trans_D(neq*neq,neq*neq),Trans_B(neq*neq,neq*neq);

    for(int i=0;i<neq*neq;i++){
             Trans_B(i,i)=1;
              Trans_I(i,i)=1;
    }

    fun_Matrix(Trans_E,M,D);

    for(int j=1;j<=nexp;j++){
          if((j%2)==1){
          mult(Trans_E,Trans_I,Trans_D);
          add(scaled(Trans_D,pow((dt),j)/factorial(j)),Trans_B);
          }else{
            mult(Trans_E,Trans_D,Trans_I);
            add(scaled(Trans_I,pow((dt),j)/factorial(j)),Trans_B);
          }
         }

  for(int t=1;t<(ninterval+1);t++){

      cout<<t<<endl;
     if((t%2)==1){
      mult(Trans_B,Trans,Trans_D);
     }else{
      mult(Trans_B,Trans_D,Trans);
     }
  }

  if((ninterval%2)==1)
  copy(Trans_D,Trans);

}


int D1_coef (int L,int F,int mf){
     if(F==3)
      return 32-(L*16+(mf+4));
     if(F==4)
      return 32-(L*16+7+(mf+5));
 }


int sweep(doub period,int period_steps,int sweep_steps,doub Max_detune,doub peak_1,doub peak_2,doub convergence,doub convergence_threshold,int conS,int expN,int Msteps)
{

  clear(A);
  clear(R);
  clear(r);
  clear(Rc);
  clear(R_L);
  doub phase=0;
  ninterval=period_steps-(period_steps%2);
  dt=period/ninterval;
  fstream file1,file2;//file1:紀錄輸入的參數。file2://紀錄計算結果
  nexp=expN;
  stringstream strstream,strstream2;
  string filename,filename2;
  strstream<<"CW"<<"dt_"<<period/period_steps<<peak_1<<"_uWcm2_"<<convergence<<"_conS_"<<conS<<"_O="<<nexp<<".txt";
  strstream>>filename;
  cout<<filename.c_str()<<endl;
  file2.open(filename.c_str(),ios::out | ios::trunc);
  strstream2<<"FS_"<<"CW"<<"dt_"<<period/period_steps<<peak_1<<"_uWcm2_"<<convergence<<"_conS_"<<conS<<"_O="<<nexp<<".txt";
  strstream2>>filename2;
  file1.open(filename2.c_str(), ios::out | ios::trunc);
  file2.precision(15);
  col_matrix< vector<doub> > y0I(neq,neq);
  //initial condition
  col_matrix< vector<doub> > y0R(neq,neq);
  //initial condiion
  Atom atom;

      for(int j=3; j<5;j++)
         for(int k=-j;k<j+1;k++)
             for(int m=3; m<5;m++)
                for(int n=-m;n<m+1;n++)
                  for(int q=-1;q<2;q++)
                     A(D1_coef(0,j,k),D1_coef(1,m,n))+=pow(atom.coef(q,1,0,m,j,n,k,1.5,0.5,3.5),2)*LineWidth;

  //initialzing the A coefficients

 long Matrix_Step = pow(2,(Msteps-1));
 pulse_average=(conS/Matrix_Step+1);

#pragma omp num_threads(2)
#pragma omp parallel for
for(int thread=0;thread<2;thread++)
{

   col_matrix< vector<doub> > M(neq,neq);//Rabifrequence*2
   int ninterval_m = (npulse-1);
   col_matrix< vector<doub> > EnerDet(neq,neq);
   vector<doub> EnergyDiff(neq-1);

for(int m=0;m<=sweep_steps;m++)
{

   cout<<m<<endl;
   doub De=1.0*m*(pow(-1,omp_get_thread_num()))*1.0/sweep_steps;
   doub detune_1=Max_detune*De,detune_2=0;
   dt=period/period_steps;

    col_matrix< vector<doub> > Result(neq*neq,pulse_average+1);
    col_matrix<vector<doub> > Trans(neq*neq,neq*neq);
    // IMPORTANT!!!!!!!
    // Remember to change it back to dense_matrix<doub> when dealing with left plus right circular polarization, since the matrix would be larger than.
    // IMPORTANT!!!!!!!!
    col_matrix< vector<doub> > Trans_AVE(neq*neq,neq*neq);

   clear(EnergyDiff);
   clear(EnerDet);
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
                    if(m==3){
                    EnerDet(D1_coef(0,j,k),D1_coef(1,m,n))=0;
                    EnerDet(D1_coef(1,m,n),D1_coef(0,j,k))=0;
                    }else if(j==3){
                       EnerDet(D1_coef(0,j,k),D1_coef(1,m,n))+=EnergyDiff[24];
                       EnerDet(D1_coef(1,m,n),D1_coef(0,j,k))=-EnerDet(D1_coef(0,j,k),D1_coef(1,m,n));
                    }
                }

clear(EnerDet);
    for(int j=3; j<5;j++)
         for(int k=-j;k<j+1;k++)
             for(int m=3; m<5;m++)
                for(int n=-m;n<m+1;n++){
                if(j==4){
                     EnerDet(D1_coef(0,j,k),D1_coef(1,m,n))+=detune_1;
                     EnerDet(D1_coef(1,m,n),D1_coef(0,j,k))-=detune_1;
                  }else{
                    EnerDet(D1_coef(0,j,k),D1_coef(1,m,n))+=detune_2;
                    EnerDet(D1_coef(1,m,n),D1_coef(0,j,k))-=detune_2;
                  }
                }

    for(int j=3; j<5;j++)
         for(int k=-j;k<j+1;k++)
             for(int m=3; m<5;m++)
                for(int n=-m;n<m+1;n++){
                if(j==4&&m==3){
                     EnerDet(D1_coef(0,j,k),D1_coef(0,m,n))+=detune_1;
                     EnerDet(D1_coef(0,m,n),D1_coef(0,j,k))-=detune_1;
                }
                if(j==3&&m==4){
                     EnerDet(D1_coef(1,j,k),D1_coef(1,m,n))+=detune_1;
                     EnerDet(D1_coef(1,m,n),D1_coef(1,j,k))-=detune_1;
                }

                }

  //set detuning of laser field positve represent blue detuning.



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
   R_gr[i+neq-neq_gr]=0.000001*2*pi;

 //initailizing for relaxation rate

    for(int i=0;i<neq*neq;i++){
        Trans(i,i)=1.0;
    }



  for(int j=3; j<5;j++)
    for(int t=-j;t<j+1;t++)
	  for(int m=3; m<5;m++)
	    for(int n=-m;n<m+1;n++){
         if(m==4){
	      M(D1_coef(1,j,t),D1_coef(0,m,n))=(atom.coef(1,1,0,j,m,t,n,1.5,0.5,3.5)+atom.coef(1,1,0,j,m,t,n,1.5,0.5,3.5))/2*sqrt(peak_1/100)*peakO;
	      }else{
	      M(D1_coef(1,j,t),D1_coef(0,m,n))=(atom.coef(1,1,0,j,m,t,n,1.5,0.5,3.5)+atom.coef(1,1,0,j,m,t,n,1.5,0.5,3.5))/2*sqrt(peak_2/100)*peakO;
	      }
	      M(D1_coef(0,m,n),D1_coef(1,j,t))=M(D1_coef(1,j,t),D1_coef(0,m,n));}
//initailizing for M matrices
//The reason to set the matrix this way(the second equation) is that actually calulated transition would be pure imaginary, but we set is to real(multiply a phase).
//If we directly set M and run through the parameter, we will get a extra munus sign in the symmetric terms, which can't be used in the formalism applied in fun_matrix.


solve_Martix(M,Trans,EnerDet,dt);

dense_matrix < doub > Trans_1(neq*neq,neq*neq),Trans_2(neq*neq,neq*neq),Trans_3(neq*neq,neq*neq);

cout<<"end of solve"<<endl;

copy(Trans,Trans_1);
copy(Trans,Trans_2);

for(long j=0;j<(Msteps-1);j++){
   mult(Trans_1,Trans_2,Trans_3);
   copy(Trans_3,Trans_1);
   copy(Trans_3,Trans_2);
   cout<<j<<endl;
}

copy(Trans_1,Trans);

cout<<"aa"<<endl;
int k=0,flag=0;
doub diff=0;


while(flag<pulse_average){

    for(int a=0;a<neq;a++)
      for(int b=0;b<neq;b++){
	Result(RealComp(a,b),(k+1)%(pulse_average+1))=0;
	Result(ImagComp(a,b),(k+1)%(pulse_average+1))=0;}

    mult(Trans,mat_col(Result,(k)%(pulse_average+1)),mat_col(Result,(k+1)%(pulse_average+1)));
   cout<<k<<endl;
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

     file1<<setiosflags(ios::left)<<setw(30)<<k*period*Matrix_Step<<setiosflags(ios::left)<<setw(30)<<data_sum<<endl;
   }

  }

cout<<"n_period="<<k<<endl;

doub buffer=0,bufferC=0;

 for (int j=0;j<16;j++)
     buffer+=Result(j,k%(pulse_average+1));



 for(int l=-3;l<4;l++)
 for(int n=-4;n<5;n++){
   bufferC+=(pow(Result(RealComp(D1_coef(0,3,l),D1_coef(0,4,n)),k%(pulse_average+1)),2)+pow(Result(ImagComp(D1_coef(0,3,l),D1_coef(0,4,n)),k%(pulse_average+1)),2));
   }

bufferC=bufferC/63;

 file2<<setiosflags(ios::left)<<setw(30)<<detune_1/2/pi;
 file2<<setiosflags(ios::left)<<setw(30)<<buffer;
 file2<<setiosflags(ios::left)<<setw(30)<<bufferC;
 file2<<setiosflags(ios::left)<<setw(30)<<k*Matrix_Step;
 file2<<setiosflags(ios::left)<<setw(30)<<m<<endl;
}

///////////////////////////////End of Sweeping//////////////////////////////////


}
 return 0;

}

