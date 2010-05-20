#!/usr/bin/env python

from numpy import loadtxt
from pylab import *
from scipy import *
from scipy.optimize import leastsq
#import sys

g2 = ['%g' %(20.0/2**(i)) for i in range(9)]

def fit(line):
  filename = 'CW_g2_%s_LW_0.7_dt_0.01_5000_uWcm2_1e-10_conS_10_O=12.txt' %line
  print filename
  #filename = sys.argv[1]
  data = loadtxt(filename,skiprows=1)
  x = data[:,0]
  y = data[:,1]
  x = x*1e7
  y = y*1e10
  fit_func = lambda p,x : p[0]+p[3]*(p[2]/2/pi)/((x-p[1])**2+(p[2]/2)**2) 
  error_func = lambda p,x,y : fit_func(p,x)-y
  p0 = [ 5 , 0 , 1 , -10]
  p1,sucess = leastsq(error_func,p0[:],args=(x,y), maxfev=10000)
  xfit = linspace(min(x),max(x), 1001)
  plot(x,y,'.')
  plot(xfit,fit_func(p1,xfit),'r-')
  savefig('CW_g2_%s_LW_0.7_dt_0.01_5000_uWcm2_1e-10_conS_10_O=12.png'%line) 
  clf()
  pfinal=r_[float(line),p1]
  return pfinal

a = [fit(gg) for gg in g2]
#ll = array(['detune','y_offset','x_offset','width','peak'])
#a.insert(0,ll)
print a
savetxt('data.txt', a)
