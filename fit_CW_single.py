#!/usr/bin/env python

from numpy import loadtxt
from pylab import *
from scipy import *
from scipy.optimize import leastsq
import sys


filename = sys.argv[1]
data = loadtxt(filename,skiprows=1)
x = data[:,0]
y = data[:,1]
x = x*1e7
y = y*1e10
fit_func = lambda p,x : p[0]+p[3]*(p[2]/2/pi)/((x-p[1])**2+(p[2]/2)**2)+p[4]*x
error_func = lambda p,x,y : fit_func(p,x)-y
p0 = [ 1000 , 0 , 1 , -50, 0]
p1,sucess = leastsq(error_func,p0[:],args=(x,y), maxfev=50000)
xfit = linspace(min(x),max(x), 1001)
plot(x,y,'.')
plot(xfit,fit_func(p1,xfit),'r-')
savefig(filename+".png") 
clf()
p1[0]=p1[0]/1e10
p1[1]=p1[1]/1e7*1e9
p1[3]=p1[3]/1e10/p1[2]*2/pi
p1[2]=p1[2]/1e7*1e9
pfinal=r_[p1,abs(p1[3]/p1[0])]

data = open( "data_"+filename, 'w')
data.write('y_offset\tx_offset(Hz)\twidth(Hz)\tpeak\tslope\tsignal_ratio\n')

for element in pfinal:
  data.write('%g\t' %element)

data.close()

