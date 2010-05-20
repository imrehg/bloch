#!/usr/bin/env python

from numpy import loadtxt
from pylab import *
from scipy import *
from scipy.optimize import leastsq
import sys

#filename = 'comb_g2_50_LW_0.7GHz_140uWcm2_5e-07_conS_10_O=12_N1_10_N2_600_D_0MHz_S.txt'
filename = sys.argv[1]
data = loadtxt(filename)
x = data[:,0]
y = data[:,1]
x = (x - 0.09192631)*1e10
y = y*1e8
fit_func = lambda p,x : p[0]+p[3]*(p[2]/2/pi)/((x-p[1])**2+(p[2]/2)**2) 
error_func = lambda p,x,y : fit_func(p,x)-y
p0 = [ 1.60 , 0 , 10 , -1.2]

p1,sucess = leastsq(error_func,p0[:],args=(x,y), maxfev=2000)
print p1
xfit = linspace(min(x),max(x), 1001)
plot(x,y,'.')
#plot(x,error_func(p1,x,y),'.')
plot(xfit,fit_func(p1,xfit),'-')
show()
