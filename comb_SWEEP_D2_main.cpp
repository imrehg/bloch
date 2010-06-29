#include "comb.h"
#include <ctime>

int main()
{

  fstream file1;
  file1.open("./Data/result.txt",ios::out | ios::trunc);
  file1<<"Function\t"<<"Order\t"<<"Steps\t"<<endl;

  for(int i=1;i<6;i++){
   for(int j=1;j<4;j++){
    time_t T1 = time(NULL),T2;
    sweep(20*2*pi,0.7*2*pi,1,400000,1000,1/pow(10,10),0,10,12-i,10*j,600,15,0,0.001,"Gaussian");
    T2 = time(NULL);
    file1<<"Gaussian\t"<<12-i<<"\t"<<10*j<<"\t"<<(T2-T1)<<endl;
    sweep(20*2*pi,0.7*2*pi,1,400000,1000,1/pow(10,10),0,10,12-i,10*j,600,15,0,0.001,"Sinc");
    T1 = time(NULL);
    file1<<"sinc\t"<<12-i<<"\t"<<10*j<<"\t"<<(T1-T2)<<endl;
   }
  }

      //gamma2(Hz),line width(GHz),sweep steps, total steps, power(uW/cm2),convergence condition,convergence threshold,convergency steps,ADM order,n1 interval,n2   interval,Matrix self multiplication,laser detune,A Factor
      // The last argument is the function type "Sinc" for Sinc Function other for Gaussian function

//0.7*2*pi;//0.0052227*2*pi;//

system("./mail.py");

return 0;

}
