#include </media/ProgramData/C/Multilevel/comb_SWEEP_MP_Carrier.cpp>
extern int sweep();

int main()
{

      sweep(40,10,2000,1/pow(10,14),20,200,600,0.005*2*pi*2);
      //sweep steps, total steps, power(uW/cm2),convergence condition,ADM order,n1 interval,n2 interval,laser detune

return 0;

}
