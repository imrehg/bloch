#ifndef COMB_H_INCLUDED
#define COMB_H_INCLUDED
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
extern int sweep(int steps,int total_steps,long double PeakPower,long double convergence,int conS,int expN,int n1, int n2,long double detune);
#endif // COMB_H_INCLUDED
