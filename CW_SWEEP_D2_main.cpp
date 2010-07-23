#include "CW.h"


int main()
{  
  for(int j =1;j<11;j++)
       sweep(20*2*pi,0.7*2*pi,0.1,10,40,0.00000002*2*pi,50*(11-j),50*(11-j),1/pow(10,10),0,10,13,25);
      //gamma2(Hz),Line Width,period,steps in period,sweep_steps,Max_detune, power of laser 1(uW/cm2), power of laser 2(uW/cm2),convergence condition,convergence
      //threshold,convergency steps,ADM order,Matrix self multiplication

return 0;

}
