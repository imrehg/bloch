#ifndef COMB_H_INCLUDED
#define COMB_H_INCLUDED
//#define GMM_USES_LAPACK
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
#include <time.h>
#include <omp.h>//For openmp
#include <gmm/gmm.h>//for gmm library
#include "atom.cpp" //for atom struct
using namespace gmm;
using namespace std;
typedef long double doub;
extern int sweep(doub g2, doub LineW,int steps,int total_steps,doub PeakPower,doub convergence,doub convergence_threshold,int conS,int expN,int n1, int n2,int Msteps,doub detune,doub A_Factor,string Func);
extern doub sweep_single(doub period_set,doub PeakPower,doub convergence,doub convergence_threshold,int conS,int expN,int n1, int n2,int Msteps,doub detune);
const doub pi=3.14159265358979323846264338327950288419716939937511;
#endif // COMB_H_INCLUDED
