#include </media/ProgramData/C/FourLevelComb/comb_FAST_SWEEP_MPtest.cpp>
extern int sweep();

int main()
{
 #pragma omp parallel for
for(int thread=0;thread<2;thread++)
      sweep();
return 0;
}
