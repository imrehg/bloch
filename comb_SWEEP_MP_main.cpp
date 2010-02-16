#include </media/DATA/C/master/comb_SWEEP_MP.cpp>
extern int sweep();

int main()
{

//    for(int i=0;i<10;i++)
      sweep(40,500000,200,5/pow(10,4),500,12,200,600,0);
      //sweep steps, total steps, power(uW/cm2),convergence condition,convergency steps,ADM order,n1 interval,n2 interval,laser detune

return 0;

}
