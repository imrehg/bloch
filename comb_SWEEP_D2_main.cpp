#include "comb.h"


int main()
{
//    for(int i=0;i<20;i++)
    sweep(20*2*pi,0.7*2*pi,10,10000000,140,1/pow(10,10),0,10,12,10,600,15,0);
      //gamma2(Hz),line width(GHz),sweep steps, total steps, power(uW/cm2),convergence condition,convergence threshold,convergency steps,ADM order,n1 interval,n2 interval,Matrix self multiplication,laser detune

//0.7*2*pi;//0.0052227*2*pi;//

system("./mail.py");
return 0;

}
