#define GMM_USES_LAPACK
#include "comb.h"

#ifndef COMB
#define COMB

class comb {


////*Laser parameters*////

doub phase;
//Phase of the comb laser, which set to zero for all current case

int pulse_average;
//number of pulses to average over for stop condition

int npulse;
//number for the maxium pulses applied

int ninterval_1,ninterval_2,ninterval_b,ninterval;
//interval_1 =steps in interval 1 ;ninterval_b= ninterval in a free iteration,ninterval= iteration terms

doub period0;
//the center rep rate

doub lasDe,peakO,FWHM;
//frequency is the detune of the carrier frequency,peak0 is the peak for the gaussian for about 150uW in FWHM 5ps
//peak0=1.34163815218652164669542605053 for 2ps

////*Laser parameters setup end*////

////*Atom parameters*////

doub LineWidth;
//The line width of upper level,0.0052227*2*pi is natural line width;

const int  neq,neq_gr;
// neq = nuber of equations(states), neq_gr = number of the ground states

int nexp;
//nexp= terms of ADM expansion,

vector<doub> r;
//total decay constant

vector<doub> R;
//relaxation rate

vector<doub> R_gr;
//relazation rate of ground state

col_matrix< vector<doub> > A;
//A matrix

vector<doub> EnergyDiff;
//Level differency Vector



////*Atom parameters end*////

comb();
int RealComp(int i,int j);
int ImagComp(int i,int j);
int factorial (int num);
doub ReRabi(doub x,doub period,doub peak);
doub ImRabi(doub x,doub period,doub peak);
void fun_Matrix(col_matrix< vector<doub> > &Trans, col_matrix< vector<doub> >&H,col_matrix< vector<doub> >&D);
void solve_Martix_new(col_matrix< vector<doub> > &M, col_matrix<vector<doub> > &Trans,col_matrix<vector<doub> > &Trans_Blank, col_matrix< vector<doub> >&D,doub dt1,doub dt2);
void solve_Martix(col_matrix< vector<doub> >&M, col_matrix<vector<doub> > &Trans, col_matrix< vector<doub> >&Trans_Ave, col_matrix< vector<doub> >&D,doub dt1,doub dt2);// solve(presultI,presultR,M,k)
int D1_coef (int L,int F,int mf);
int sweep(doub g2,doub LineW,int steps,int total_steps,doub PeakPower,doub convergence,doub convergence_threshold,int conS,int expN,int n1, int n2,int Msteps,doub detune);

};
#endif
