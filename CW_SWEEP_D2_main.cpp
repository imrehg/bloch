#include "CW.h"


int main()
{
     for(int i=1;i<15;i++)
     sweep(50*2*pi,0.7*2*pi,0.1,10,200,0.000040*2*pi,50*i,50*i,5/pow(10,7),0,10,12,22);
      //gamma2(Hz),Line Width,period,steps in period,sweep_steps,Max_detune, power of laser 1(uW/cm2), power of laser 2(uW/cm2),convergence condition,convergence
      //threshold,convergency steps,ADM order,Matrix self multiplication



return 0;

}
