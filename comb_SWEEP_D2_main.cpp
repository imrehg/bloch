#include "comb.h"


int main()
{

    sweep(20*2*pi,0.7*2*pi,50,20000000,1000,1/pow(10,10),0,10,12,100,600,15,0,0.001,"Gaussian");
      //gamma2(Hz),line width(GHz),sweep steps, total steps, power(uW/cm2),convergence condition,convergence threshold,convergency steps,ADM order,n1 interval,n2 interval,Matrix self multiplication,laser detune,A Factor
      // The last argument is the function type "Sinc" for Sinc Function other for Gaussian function

//0.7*2*pi;//0.0052227*2*pi;//

system("./mail.py");
return 0;

}
