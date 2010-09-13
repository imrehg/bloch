#include "CW.h"


int main()
{

    doub gamma_2 = 1000*2*pi;                   // Decoherence Rate.(Hz)
    doub upperLevelLineWidth = 0.7*2*pi;        // Upper Leve Line Width.(GHz)
    int stepsInPeriod = 10;                     // Steps in one Matrxi unit.
    doub period = 0.1;                          // period in one Matrix unit(ns).
    doub maxDetune = 0.0001*2*pi;               // Maxium Detuen of laser angular frequency.(GHz)
    int sweepSteps = 40;                        // Steps in sweep of max detune.
    doub laserPower1 = 50;                      // Laser Power 1.(uW/cm2).
    doub laserPower2 = 50;                      // Laser Power 1.(uW/cm2).
    doub convergenceCondition = 1/pow(10,10);   // If the reletive difference of one density matrix element in subsequent steps is smaller than this value, the element is considered as converge.
    doub convergenceThreshold = 0;              // If any density matrix element smaller than this value is considered as alreay converge.
    int convergenceSteps = 10;                  // If all elements of density matrix converge in the all last pulses steps of this number, the density matrix considered as converge.
    int orderOfADM = 13;                        // The expanding order of ADM method.
    int matrixSelfMult = 25;                    // The self multiplication of the final pulse evolution matrix.
    doub polarizationOfLight = 1;               // The Polarization of ligh +-1 represent cirular polarizaed ligh; 0 for linear polarized light. Default is 1;
    int dLine = 2;                              // D line. 2 represent D2 line 1 for D1 line. Defualt is D2 line.
    int numOfThread = 2;                        // Number of CPU cores.

    sweep(gamma_2,upperLevelLineWidth,period,stepsInPeriod,sweepSteps,maxDetune,laserPower1,laserPower2,convergenceCondition,convergenceThreshold,convergenceSteps,orderOfADM,matrixSelfMult,polarizationOfLight,dLine,numOfThread);
    //gamma2(Hz),Line Width,period,steps in period,sweep_steps,Max_detune, power of laser 1(uW/cm2), power of laser 2(uW/cm2),convergence condition,convergence
    //threshold,convergency steps,ADM order,Matrix self multiplication

    return 0;

}
