#include "comb.h"
#include <ctime>

int main()
{

sweep(20*2*pi,0.7*2*pi,20,30000000,100,1/pow(10,10),0,10,13,50,600,15,0,0.001,"Sinc");
sweep(20*2*pi,0.7*2*pi,20,30000000,100,1/pow(10,10),0,10,13,30,600,15,0,0.001,"Gaussian");

//gamma2(Hz),line width(GHz),sweep steps, total steps, power(uW/cm2),convergence condition,convergence threshold,convergency steps,ADM order,n1 interval,n2   interval,Matrix self multiplication,laser detune,A Factor
// The last argument is the function type "Sinc" for Sinc Function other for Gaussian function

//0.7*2*pi;//0.0052227*2*pi;//

system("./mail.py");

return 0;

}

//  fstream file1;
//  file1.open("./Data/result_A.txt",ios::out | ios::trunc);
//  file1<<"Function\t"<<"A"<<endl;
//  doub A;
//
//   for(int j=2;j<5;j++){
//      time_t T1 = time(NULL),T2;
//      A = 1/pow(10,j);
//      sweep(20*2*pi,0.7*2*pi,5,8000000,100,1/pow(10,10),0,10,13,50,600,15,0,A,"Sinc");
//      T2 = time(NULL);
//      file1<<"Sinc\t"<<A<<"\t"<<float(T2-T1)<<endl;
//   }
//
//    for(int j=2;j<5;j++){
//      time_t T1 = time(NULL),T2;
//      A = 1/pow(10,j);
//      sweep(20*2*pi,0.7*2*pi,5,8000000,100,1/pow(10,10),0,10,13,30,600,15,0,A,"Gaussian");
//      T2 = time(NULL);
//      file1<<"Gaussian\t"<<A<<"\t"<<float(T2-T1)<<endl;
//   }
