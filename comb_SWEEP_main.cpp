#include "comb.h"
#include <ctime>

int main()
{

    doub gamma_2 = 1000*2*pi;                   // Decoherence Rate.(Hz)
    doub upperLevelLineWidth = 0.7*2*pi;        // Upper Leve Line Width.(GHz)
    int numOfThread = 8;                         // Number of CPU cores.
    int sweepStep = 1;                          // The Actually sweep steps for each core.
    int totalSweepStep = 20000;                 // This step is set to be the total step in the rep rate = clock transition /100 ~ clock transtion/101. About 9.1016KHz.
    doub laserPower = 800;                      // Average Laser Power(uW).
    doub convergenceCondition = 1/pow(10,10);   // If the reletive difference of one density matrix element in subsequent steps is smaller than this value, the element is considered as converge.
    doub convergenceThreshold = 0;              // If any density matrix element smaller than this value is considered as alreay converge.
    int convergenceSteps = 10;                  // If all elements of density matrix converge in the all last pulses steps of this number, the density matrix considered as converge.
    int orderOfADM = 13;                        // The expanding order of ADM method.
    int n1Interval = 50;                        // The time interval in the first interval of pulse - pulse region. 50 should be set for sinc 30 for gaussina under 1000uW.
    int n2Interval = 600;                       // The time interval in the second interval of pulse - free evolusion region.
    int matrixSelfMult = 15;                    // The self multiplication of the final pulse evolution matrix.
    doub laserDetune = 0;                       // The laser detuning of carrier frequency.(GHz)
    doub pulseAFactor = 0.001;                  // The A factor of the pulse.(ns)
    string pulseFunction = "Sinc";              // The type of the pulse shape is either "Gaussian" or "Sinc". Anything else is set to Gaussina.
    doub polarizationOfLight = 1;               // The Polarization of ligh +-1 represent cirular polarizaed ligh; 0 for linear polarized light. Default is 1;
    int dLine = 2;                              // D line. 2 represent D2 line 1 for D1 line. Defualt is D2 line.

    sweep(gamma_2,upperLevelLineWidth,sweepStep,totalSweepStep,laserPower,convergenceCondition,convergenceThreshold,convergenceSteps,orderOfADM,n1Interval,n2Interval,matrixSelfMult,laserDetune,pulseAFactor,pulseFunction, numOfThread, polarizationOfLight, dLine);

    //gamma2(Hz),line width(GHz),sweep steps, total steps, power(uW/cm2),convergence condition,convergence threshold,convergency steps,ADM order,n1 interval,n2   interval,Matrix self multiplication,laser detune,A Factor
    // The last argument is the function type "Sinc" for Sinc Function other for Gaussian function


    //system("./mail.py");

    return 0;

}


