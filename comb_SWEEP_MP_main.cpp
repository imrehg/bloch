#include </media/ProgramData/C/Multilevel/comb_SWEEP_MP.cpp>
extern int sweep();

int main()
{
    for(int i=-5;i<5;i++)
      sweep(200,1000,10000+i*500);//last parameter is power uW/cm2

return 0;

}
