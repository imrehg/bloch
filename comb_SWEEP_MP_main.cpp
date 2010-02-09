#include </media/ProgramData/C/Multilevel/comb_SWEEP_MP.cpp>
extern int sweep();

int main()
{
    for(int i=0;i<10;i++)
      sweep(200,1000,150*pow(2,i));

return 0;

}
