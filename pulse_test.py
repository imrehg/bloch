#!/usr/bin/env python

from pylab import *
from cmath import exp
import sys
import ConfigParser

ion()

try:
    configfile = 'pulse_test.config'
    config = ConfigParser.ConfigParser()
    config.readfp(open(configfile))
except:
    print "Cannot find configuration file."
    sys.exit(1)


def shift(X):
    Y = copy(X)
    Y[:Points*Number_Pulses/2], Y[Points*Number_Pulses/2:] = X[Points*Number_Pulses/2:], X[:Points*Number_Pulses/2]
    return Y

Gaussian = lambda width,x: exp(-(x/width)**2)

Sech_2 = lambda width,x: (2/(exp(-x/width)+exp(x/width)))**2

Window_Function = lambda x,start_x,end_x: 1 if start_x<=x<=end_x else 0

Period = config.getfloat("Other","period")
Number_Pulses =  config.getint("Other","number_pulses")
Points = config.getint("Other","points")

def get_pulse_function():
    
    Function = Gaussian
    Width = config.getfloat("Pulse_Function","w_0")
    wstart = config.getfloat("Pulse_Function","window_s")
    wend = config.getfloat("Pulse_Function","window_e")
    Freq = [ 1/(Period*Number_Pulses)*a for a in range( -Points*Number_Pulses/2, Points*Number_Pulses/2 ) ]
    x = [ Period*(a)/Points for a in range( -Points/2, Points/2) ]
    X = [ Period*(a)/Points for a in range( -Points*Number_Pulses/2, Points*Number_Pulses/2 ) ]
    y = [ Function( Width,a) for a in x ] * Number_Pulses
    Y=  shift(fft(y))
    W = [ Window_Function( j, wstart, wend) for j in Freq ]
    Y = ifft( shift( [ Y[i]*W[i] for i in range(len(Y)) ]) )
    
    return X, Y.real, y


if  __name__ == '__main__':

    ######################## Try to Choose the Pulse Function #############################################

    Input = 0

    while( not Input or Input not in (1,2,3)):
      try: Input = int (raw_input("Enter file func to use (1): Gaussian (2): Sech_2 (3)All:"))
      except: Input = 0

    Function = [] 
    Function_Name = []
     
    if  Input == 1: Function.append(Gaussian); Function_Name.append('Gaussian')
    elif Input == 2: Function.append(Sech_2); Function_Name.append('Sech_2')
    else: Function.append(Gaussian); Function.append(Sech_2); Function_Name.append('Gaussian'); Function_Name.append('Sech_2')

    ####################### End of Section ############################################################

    flag = ''

    while (flag != 'exit'):

        Width = ['','']

        for i in range(len(Function)):
            while( not Width[i] ):     
                try: 
                    Width[i] = float( raw_input("Input width for %s (%f ns):"%(Function_Name[i], config.getfloat("Pulse_Function","w_%d"%i) ) ) ) 
                    config.set("Pulse_Function","w_%d"%i, "%f"%Width[i])
                except: print "Use default value"; Width[i] = config.getfloat("Pulse_Function","w_%d"%i)
        

        try:
            wstart = float( raw_input("Input start point of window function(%f Ghz):"%config.getfloat("Pulse_Function","window_s")))
            config.set("Pulse_Function","window_s", "%f"%wstart)
        except: print  "Use default value"; wstart = config.getfloat("Pulse_Function","window_s")

        try: 
            wend = float( raw_input("Input end point of window function(%f Ghz):"%config.getfloat("Pulse_Function","window_e") )) 
            config.set("Pulse_Function","Window_E", "%f"%wend)
        except: print  "Use default value"; wend = config.getfloat("Pulse_Function","window_e")

        x = [ Period*(a)/Points for a in range( -Points/2, Points/2) ]
        X = [ Period*(a)/Points for a in range( -Points*Number_Pulses/2, Points*Number_Pulses/2 ) ]
        y = [ [ Function[i]( Width[i],a) for a in x ] * Number_Pulses  for i in range( len( Function ) ) ]
        Y = [ shift( ( fft(spam) ) ) for spam in y ]
        Freq = [ 1/(Period*Number_Pulses)*a for a in range( -Points*Number_Pulses/2, Points*Number_Pulses/2 ) ]
        W = [ Window_Function( j, wstart, wend) for j in Freq ]
        max_value= max([max(i) for i in Y])
        W_mod = [ i* max_value for i in W]        
        
        subplot(311)

        for a in range(len(Y)):
         plot( Freq, abs(Y[a]), '.' if Function_Name[a] == 'Gaussian'  else 'r-' ) 
         xlabel('Spectrum')
        
        plot( Freq, W_mod, 'r-' )

        subplot(312)

        for a in range(len(y)):
         plot(X, y[a], '.' if Function_Name[a] == 'Gaussian' else 'r-' ) 
         xlabel('Pulse Shape')
        
    ##    xlim( -Width[0]*20, Width[0]*20 )

        subplot(313)
      
        Y = [ ifft(shift( [ a[i]*W[i] for i in range(len(a)) ]))  for a in Y ]

        for a in range(len(Y)):
             plot(X, Y[a], 'r.' if Function_Name[a] == 'Gaussian'  else '-' ) 
             xlabel('Pulse after windowed')
      
        
    ##    xlim( -Width[0]*20, Width[0]*20 )

        subplots_adjust(hspace=0.3)

        flag = raw_input("Type exit to exti or anything to continue:")
        clf()
        
    config.write(open(configfile,'w'))



