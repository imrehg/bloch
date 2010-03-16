#include "comb.h"


int main()
{

    for(int i=1;i<5;i++)
      sweep(20,50000,2000,5.0/pow(10,3)/i,30,8,200,600,0);
      //sweep steps, total steps, power(uW/cm2),convergence condition,convergency steps,ADM order,n1 interval,n2 interval,laser detune

return 0;

}
