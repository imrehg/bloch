#!/usr/bin/env python
from pylab import *
from scipy import stats

def sc( rang, points, a):
    
    x = [ i *rang/ points for i in range( -points/2, points/2 ) ]
    y = [ sinc( i / a)**2 for i in x]
    return sum(y)*rang/points

if __name__ == "__main__":

    a =  range(1,10)
    data = [ sc(i*20.0,50000,i) for i in a ]
    gradient, intercept, r_value, p_value, std_err = stats.linregress(a,data)
    plot( a, data )
    show()
    print  gradient, intercept, r_value, p_value, std_err
    
