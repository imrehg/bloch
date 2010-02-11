#include </media/ProgramData/C/Multilevel/comb_SWEEP_MP.cpp>
extern int sweep();

int main()
{

    for(int i=-5;i<5;i++)
      sweep(50,100000,2000,1/pow(10,12),20,200,600,0.0087*i/5*2*pi);
      //sweep steps, total steps, power(uW/cm2),convergence condition,ADM order,n1 interval,n2 interval,laser detune

return 0;

}
