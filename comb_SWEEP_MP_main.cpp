#include </media/ProgramData/C/Multilevel/comb_SWEEP_MP_test.cpp>
extern int sweep();

int main()
{

    for(int i=0;i<10;i++)
      sweep(40,500000,200,1/pow(10,16),i+10,200,600,0);
      //sweep steps, total steps, power(uW/cm2),convergence condition,ADM order,n1 interval,n2 interval,laser detune

return 0;

}
