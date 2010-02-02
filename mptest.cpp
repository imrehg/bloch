#include <omp.h>
#include <iostream>
using namespace std;
int main()
{

   #pragma omp parallel for
    for(int i=0;i<100;i++)
       cout<<i<<endl;
    return 0;
}
