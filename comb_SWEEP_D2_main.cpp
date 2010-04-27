#include "comb.h"


int main()
{
//    for(int i=1;i<11;i++)
     sweep(0.7*2*pi,80,2000000,100000,1/pow(10,6),0,10,12,10,600,10,0);
      //line width,sweep steps, total steps, power(uW/cm2),convergence condition,convergence threshold,convergency steps,ADM order,n1 interval,n2 interval,Matrix self multiplication,laser detune

//0.7*2*pi;//0.0052227*2*pi;//
//fstream file;
//stringstream strstream;
//string filename;
//
//for(int j=2;j<3;j++){
//
//filename.clear();
//strstream.clear();
//strstream<<"CenterSweep_"<<j*100000<<"uW.txt";
//strstream>>filename;
//cout<<filename.c_str()<<endl;
//file.open(filename.c_str(),ios::out | ios::trunc);
//
//    for(int i=-200;i<201;i+=1){
//     file<<setiosflags(ios::left)<<setw(30)<<0.0002*i;
//     file<<setiosflags(ios::left)<<setw(30)<<sweep_single(10.87827848197104833208,100000*j,1/pow(10,4),0,10,12,10,600,15,0.0002*2*pi*i)<<endl;
//    }
//    //period, power(uW/cm2),convergence condition,convergence threshold,convergency steps,ADM order,n1 interval,n2 interval,Matrix self multiplication,laser detune
//file.close();
//
//}

return 0;

}
