#include "comb.h"


int main()
{

//    for(int i=0;i<41;i+=4)
//     sweep(15,1000000,200,1/pow(10,4),0,10,12,10,600,15,0.001*2*pi*i);
//      //sweep steps, total steps, power(uW/cm2),convergence condition,convergence threshold,convergency steps,ADM order,n1 interval,n2 interval,Matrix self multiplication,laser detune


fstream file;
stringstream strstream;
string filename;
strstream<<"CenterSweep.txt";
strstream>>filename;
cout<<filename.c_str()<<endl;
file.open(filename.c_str(),ios::out | ios::trunc);

    for(int i=-40;i<41;i+=1){
     file<<setiosflags(ios::left)<<setw(30)<<0.001*i;
     file<<setiosflags(ios::left)<<setw(30)<<sweep_single(10.87827848197104833208,200,1/pow(10,4),0,10,12,10,600,15,0.001*2*pi*i);
    }
    //period, power(uW/cm2),convergence condition,convergence threshold,convergency steps,ADM order,n1 interval,n2 interval,Matrix self multiplication,laser detune


return 0;

}
