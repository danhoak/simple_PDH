#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to fit the excess noise due to charging (?) of the L1 ITMY test mass 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np

import matplotlib
matplotlib.use("Agg")
import pylab
matplotlib.rcParams.update({'savefig.dpi':250})

import os, sys
import scipy.optimize as opt
from uncertainties import ufloat


# power law function to fit
def plaw(f,A,n):
    return A/(1+f**n)


# load the data

"""
# data files were created on the LIGO CIT cluster from LLO science data in S5
#seg_start = 830599214 # May  2 2006 10:00
#seg_start = 831186014 # May  9 2006 05:00
seg_start = 824622494 # Feb 22 2006 05:48


# the 'coarse' data files used only ~300 seconds of data (I forget exactly how much),
# the other data files have 1200 seconds of data and so have less noise
duration = 1200
Fs = 16384
chans = 'L1:LSC-STRAIN'
frame_type = 'L1_RDS_C03_L2'

# The FFT parameters to calculate the PSDs were:
stride = 6.0   # FFT stride in seconds
overlap = 3.0  # overlap in seconds (50%)
"""



#data1 = np.genfromtxt('L1_S5_ASD_2May2006_coarse.txt')
data1 = np.genfromtxt('L1_S5_ASD_2May2006.txt')
freq = data1[:,0]
May2 = data1[:,1]

#data2 = np.genfromtxt('L1_S5_ASD_9May2006_coarse.txt')
data2 = np.genfromtxt('L1_S5_ASD_9May2006.txt')
May9 = data2[:,1]

#data3 = np.genfromtxt('L1_S5_ASD_22Feb2006_coarse.txt')
data3 = np.genfromtxt('L1_S5_ASD_22Feb2006.txt')
Feb22 = data3[:,1]


# calculate the subtraction, pre-vent minus post-vent
# label of 1 is Feb22-May9, label of 2 is May2-May9
x1 = np.sqrt(Feb22**2-May9**2)
x2 = np.sqrt(May2**2-May9**2)



# data cleaning; drop the NaN values, ignore the 60Hz region
idx1 = np.argwhere(np.logical_or(np.isnan(x1), (freq>57) * (freq<63)))
idx2 = np.argwhere(np.logical_or(np.isnan(x2), (freq>57) * (freq<63)))

fit_freq1 = np.delete(freq, idx1)
fit_freq2 = np.delete(freq, idx2)
fit_x1 = np.delete(x1,idx1)
fit_x2 = np.delete(x2,idx2)


# get the indices of the frequency range to fit over
a1 = np.argmin(np.abs(fit_freq1-50))
b1 = np.argmin(np.abs(fit_freq1-102))
a2 = np.argmin(np.abs(fit_freq2-47))
b2 = np.argmin(np.abs(fit_freq2-115))


# guess the starting parameters for the fit
v0 = np.array([1e-18, 2.])

pfit1, fitCov1 = opt.curve_fit(plaw, fit_freq1[a1:b1], fit_x1[a1:b1], p0=v0)
perr1 = np.sqrt(np.diag(fitCov1))
n1 = ufloat(pfit1[1],perr1[1])

pfit2, fitCov2 = opt.curve_fit(plaw, fit_freq2[a2:b2], fit_x2[a2:b2], p0=v0)
perr2 = np.sqrt(np.diag(fitCov2))
n2 = ufloat(pfit2[1],perr2[1])

print 'Feb22-May9:', n1
print 'May2-May9:', n2

yfit1 = plaw(freq, pfit1[0], pfit1[1])
yfit2 = plaw(freq, pfit2[0], pfit2[1])


# by-hand guesses
y1 = 5.3e-19/(1+freq**2)
y2 = 2.7e-17/(1+freq**3)
y3 = 3.9e-18/(1+freq**2.5)




fignum=0

fignum+=1
pylab.figure(fignum)

ax1 = pylab.subplot(2,1,1)

pylab.loglog(freq,Feb22,'darkorange',linestyle='-',linewidth=1.0,label='L1 Before Vent (22Feb2006')
pylab.loglog(freq,May2,'b-',linewidth=1.0,label='L1 Before Vent (2May2006)')
pylab.loglog(freq,May9,'r-',linewidth=1.0,label='L1 After Vent (9May2006),\nViton EQ stop retracted')

#pylab.loglog(freq,4.5e-18/(1+freq**2.5),'k--',linewidth=1.4,label='1/f^2.5')
#pylab.loglog(freq,4.5e-18/(1+freq**2.5),'k--',linewidth=1.4,label='1/f^2.5')

pylab.xticks(visible=False)
pylab.ylim(2e-23,3e-21)
pylab.grid(True, which='both', linestyle=':',alpha=0.8)
pylab.ylabel('Strain ASD [1/rt(Hz)]',fontsize=12)
pylab.xlim(30,300)
pylab.yticks(fontsize=10)
ax1.get_yaxis().set_label_coords(-0.08,0.5)
pylab.legend(loc=1,fancybox=True,prop={'size':8},bbox_to_anchor=[1.08,1])
pylab.title('Charge migration noise at L1 - May 2006',fontsize=12)

ax2 = pylab.subplot(2,1,2)

#pylab.plot(fit_freq1,fit_x1,'c-',linewidth=0.8,label='Subtraction (feb22-May9)')
#pylab.plot(fit_freq1[a1:b1],fit_x1[a1:b1],'m.',markersize=4,label='Subtraction - data for fit')
#pylab.loglog(freq, yfit1, 'k--',linewidth=1.4,label='fit, index='+str(n1))

pylab.plot(fit_freq2,fit_x2,'c-',linewidth=0.8,label='Subtraction (May2-May9)')
pylab.plot(fit_freq2[a2:b2],fit_x2[a2:b2],'m.',markersize=4,label='Subtraction - data for fit')
pylab.loglog(freq, yfit2, 'k--',linewidth=1.4,label='powerlaw fit, index='+str(n2))

#pylab.loglog(freq, y1, 'b--',linewidth=1.4,label='1/f^2')
#pylab.loglog(freq, y3, 'k--',linewidth=1.4,label='1/f^2.5')
#pylab.loglog(freq, y2, 'r--',linewidth=1.4,label='1/f^3')

pylab.grid(True, which='both', linestyle=':',alpha=0.8)
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)
pylab.xlim(30,300)
pylab.ylim(0.5e-23,1e-21)
pylab.xlabel('Frequency [Hz]',fontsize=12)
pylab.ylabel('Strain ASD [1/rt(Hz)]',fontsize=12)
pylab.legend(loc=1,fancybox=True,prop={'size':8})

ax2.get_yaxis().set_label_coords(-0.08,0.5)

pylab.savefig('S5_L1_charging_noise.png',bbox_inches='tight')
pylab.close


