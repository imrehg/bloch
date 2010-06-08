#!/usr/bin/env python

from numpy import loadtxt
from pylab import *
from scipy import *
from scipy.optimize import leastsq
#import sys
filename = 'data_comb.txt'

g2 = ['%g' %(1e-6*10**i) for i in range(5)]
print g2

def fit(line):
  filename = 'comb_g2_20_LW_0.7GHz_1000uWcm2_1e-10_conS_10_O=12_N1_10_N2_600_D_0MHz_A_%s_S.txt' %line
  print filename
  #filename = sys.argv[1]
  data = loadtxt(filename,skiprows=1)
  x = data[:,0]
  y = data[:,1]
  x = (x-0.09192631)*1e7
  y = y*1e10
  fit_func = lambda p,x : p[0]+p[3]*(p[2]/2/pi)/((x-p[1])**2+(p[2]/2)**2)
  error_func = lambda p,x,y : fit_func(p,x)-y
  p0 = [ 500 , 0 , 1 , -1000]
  p1,sucess = leastsq(error_func,p0[:],args=(x,y), maxfev=50000)
  xfit = linspace(min(x),max(x), 1001)
  plot(x,y,'.')
  plot(xfit,fit_func(p1,xfit),'r-')
  savefig('comb_g2_20_LW_0.7GHz_1000uWcm2_1e-10_conS_10_O=12_N1_10_N2_600_D_0MHz_A_%s_S.png'%line) 
  clf()
  p1[0]=p1[0]/1e10
  p1[1]=p1[1]/1e7*1e11
  p1[2]=p1[2]/1e7*1e11
  p1[3]=p1[3]/1e10
  pfinal=r_[float(line),p1]
  return pfinal

a = [fit(gg) for gg in g2]

data = open( filename, 'w')
data.write('detune(uW)\ty_offset\tx_offset(Hz)\twidth(Hz)\tpeak\n')


for element in a:

 for index,value in enumerate(element):
  data.write('%g' %value)
  if (index+1 < len(element)):
   data.write("\t")   
 
 data.write('\n')

data.close()


