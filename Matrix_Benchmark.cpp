#define GMM_USES_LAPACK
#include <iostream> //FOR file IO
#include <fstream>//For file IO
//#include <cstdlib> //
//#include <cctype>//
//#include <cstring>
#include <string>//to use string
#include <sstream>//to use sstream
#include <cmath>//For sin cos functions
#include <iomanip>//For  setiosflags
#include <ctime>//For timer
#include <omp.h>//For openmp
#include <gmm/gmm.h>//for gmm library
#include "atom.cpp" //for atom struct
using namespace gmm;
using namespace std;
typedef double doub;
#include <ctime>
#include <time.h>


int main(){
      cout<<"start!";
    fstream file1;
    file1.open("Matrix_Benchmark.txt", ios::out | ios::trunc);

    int n=500;
    double con=0.05;
    int time=1;

    file1<<setiosflags(ios::left)<<setw(30)<<"Matrix connectivity"<<setiosflags(ios::left)<<setw(30)<<"Lapack Matrix"<<setiosflags(ios::left)<<setw(30)<<"Gmm dense Matrix"<<setiosflags(ios::left)<<setw(30)<<"Gmm read sparse col Matrix"<<setiosflags(ios::left)<<setw(30)<<"Gmm write sparse col Matrix"<<endl;
    dense_matrix<doub>  A(n,n);
    dense_matrix<doub>  B(n,n);
    dense_matrix<doub>  C(n,n);

    col_matrix< vector<doub> >  A2(n,n);
    col_matrix< vector<doub> >  B2(n,n);
    col_matrix< vector<doub> >  C2(n,n);

    col_matrix< rsvector<doub> > A3(n,n);
    col_matrix< rsvector<doub> > B3(n,n);
    col_matrix< rsvector<doub> > C3(n,n);

    col_matrix< wsvector<doub> > A4(n,n);
    col_matrix< wsvector<doub> > B4(n,n);
    col_matrix< wsvector<doub> > C4(n,n);


for (int i=1;i<11;i++){
    clear(A);
    clear(B);
    cout<<i*0.1<<endl;
    file1<<setiosflags(ios::left)<<setw(30)<<i*0.1;
    con=i*0.1;
    fill_random(A,con);
    fill_random(B,con);
//    fill_random(A4,0.5);
//    fill_random(B4,0.5);
//    fill_random(A2,0.5);
//    fill_random(B2,0.5);
//    fill_random(A3,0.5);
//    fill_random(B3,0.5);

    copy(A,A2);
    copy(A,A3);
    copy(A,A4);
    copy(B,B2);
    copy(B,B3);
    copy(B,B4);

    cout<<"start!";

    time_t T1=clock(),T2;

    for(int i=0;i<time;i++){
        mult(A,B,C);
        copy(scaled(C,0.5),A);
    }

    T2=clock();

    file1<<setiosflags(ios::left)<<setw(30)<<(T2-T1)*1.0/CLOCKS_PER_SEC/time;

    T1=clock();

    for(int i=0;i<time;i++){
        mult(A2,B2,C2);
        copy(scaled(C2,0.5),A2);
    }
    T2=clock();

    file1<<setiosflags(ios::left)<<setw(30)<<(T2-T1)*1.0/CLOCKS_PER_SEC/time;

    T1=clock();

    for(int i=0;i<time;i++){
        mult(A3,B3,C3);
        copy(scaled(C3,0.5),A3);
    }
    T2=clock();

    file1<<setiosflags(ios::left)<<setw(30)<<(T2-T1)*1.0/CLOCKS_PER_SEC/time;

    T1=clock();

    for(int i=0;i<time;i++){
        mult(A4,B4,C4);
        copy(scaled(C4,0.5),A4);
    }
    T2=clock();

    file1<<setiosflags(ios::left)<<setw(30)<<(T2-T1)*1.0/CLOCKS_PER_SEC/time<<endl;


}







}
