#include "comb.h"
#include <ctime>

int main()
{

for(int i=-3;i<2;i++)
  sweep(20*2*pi,0.7*2*pi,1,200000,1000,1/pow(10,10),0,10,12,10,600,15,0,0.001*pow(10,i),"Gaussian");

  //  sweep(20*2*pi,0.7*2*pi,20,30000000,100,1/pow(10,10),0,10,13,50,600,15,0,0.001,"Sinc");
  //sweep(20*2*pi,0.7*2*pi,20,30000000,100,1/pow(10,10),0,10,13,30,600,15,0,0.001,"Gaussian");

//gamma2(Hz),line width(GHz),sweep steps, total steps, power(uW/cm2),convergence condition,convergence threshold,convergency steps,ADM order,n1 interval,n2   interval,Matrix self multiplication,laser detune,A Factor
// The last argument is the function type "Sinc" for Sinc Function other for Gaussian function


    system("./mail.py");

    return 0;

}


