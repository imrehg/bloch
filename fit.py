#!/usr/bin/env python

from numpy import loadtxt
from pylab import *
from scipy import *
from scipy.optimize import leastsq
#import sys
<<<<<<< HEAD
filename = 'data.txt'

g2 = ['%g' %(0.5*200*(i+1)) for i in range(5)]
print g2

def fit(line):
  filename = 'CW_g2_20_LW_0.7_dt_0.01_%s_uWcm2_1e-10_conS_10_O=12.txt' %line
=======

g2 = ['%g' %(20.0/2**(i)) for i in range(9)]

def fit(line):
  filename = 'CW_g2_%s_LW_0.7_dt_0.01_5000_uWcm2_1e-10_conS_10_O=12.txt' %line
>>>>>>> c2708de87bb33bfe2255409564cfce30b0d0d6cf
  print filename
  #filename = sys.argv[1]
  data = loadtxt(filename,skiprows=1)
  x = data[:,0]
  y = data[:,1]
  x = x*1e7
  y = y*1e10
<<<<<<< HEAD
  fit_func = lambda p,x : p[0]+p[3]*(p[2]/2/pi)/((x-p[1])**2+(p[2]/2)**2)+p[4]*x
  error_func = lambda p,x,y : fit_func(p,x)-y
  p0 = [ 1000 , 0 , 1 , -50, 0]
  p1,sucess = leastsq(error_func,p0[:],args=(x,y), maxfev=50000)
  xfit = linspace(min(x),max(x), 1001)
  plot(x,y,'.')
  plot(xfit,fit_func(p1,xfit),'r-')
  savefig('CW_g2_20_LW_0.7_dt_0.01_%s_uWcm2_1e-10_conS_10_O=12.png'%line) 
  clf()
  p1[0]=p1[0]/1e10
  p1[1]=p1[1]/1e7
  p1[2]=p1[2]/1e7
  p1[3]=p1[3]/1e10
  pfinal=r_[float(line)*2,p1]
  return pfinal


a = [fit(gg) for gg in g2]
data = open( filename, 'w')
ll = ['detune(uW)\t','y_offset\t','x_offset(Hz)\t','width(Hz)\t','peak\t','slope\t']

for name in ll:
 data.write(name)

data.write('\n')

for element in a:

 for value in element:
  data.write('%g\t' %value)   
 data.write('\n')

data.close()
#data.write('')
#savetxt('data.txt', a)
=======
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
>>>>>>> c2708de87bb33bfe2255409564cfce30b0d0d6cf
