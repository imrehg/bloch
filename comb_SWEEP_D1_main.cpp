#include "comb.h"


int main()
{

    for(int i=0;i<10;i++)
     sweep(10,1000000,200,1/pow(10,4),0,10,14-i,200,600,15,0);
      //sweep steps, total steps, power(uW/cm2),convergence condition,convergence threshold,convergency steps,ADM order,n1 interval,n2 interval,Matrix self multiplication,laser detune

return 0;

}
