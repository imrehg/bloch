#include </media/ProgramData/C/Multilevel/comb_SWEEP_MP.cpp>
extern int sweep();

int main()
{
//
//    for(int i=0;i<5;i++)
      sweep(100,1000000,150,1/pow(10,16),20,200,600);//sweep steps, total steps, power(uW/cm2),convergence condition,ADM order,n1 interval,n2 interval

return 0;

}
