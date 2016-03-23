#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script to plot PDH signals and cavity buildup as a function of detuning 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
from matplotlib import pyplot
import pylab

# Functions for cavity reflectivity and transmissivity
# r_i, r_e are the *amplitude* reflectivities of the input and end mirror
# The detuning parameter 'w' is in units of f / f_FSR
def refl_coeff(r_i, r_e, w):
    return (r_e * exp(2j*pi*w) - r_i) / (1.0 - r_i*r_e*exp(2j*pi*w))

def trans_coeff(t_i, t_e, r_i, r_e, w):
    return (t_i * t_e * exp(1j*pi*w)) / (1.0 - r_i*r_e*exp(2j*pi*w))



c = 299792458.0
lam = 1064e-9

x = arange(-0.9,0.9,0.00005)

E_in = 1.0

r_i = sqrt(0.95165)
r_e = sqrt(0.25)
#r_e = sqrt(0.98624)

t_i = sqrt(1-r_i**2)
t_e = sqrt(1-r_e**2)

finesse = pi * sqrt(r_i*r_e) / (1-r_i*r_e)
print 'The finesse of this cavity is %.2f' % finesse

L = 11.952
FSR = c / (2*L)

E_circ = E_in * t_i *exp(2j*pi*x)/ (1.0 - r_i*r_e*exp(2j*pi*x))
P_circ = real(E_circ * conj(E_circ))

fmod = 6270676.17
#fmod = 56436085.57
#fmod = 8360901.56
#fmod = 131684199.67

W = fmod/FSR

# These are the equations for the PDH signals in reflection
PD_I = imag( refl_coeff(r_i, r_e, x) * conj(refl_coeff(r_i, r_e, x+W)) - conj(refl_coeff(r_i, r_e, x)) * refl_coeff(r_i, r_e, x-W) )
PD_Q = real( refl_coeff(r_i, r_e, x) * conj(refl_coeff(r_i, r_e, x+W)) - conj(refl_coeff(r_i, r_e, x)) * refl_coeff(r_i, r_e, x-W) )

y = arange(-1.12, 1.5, 0.00001)

E_trans = trans_coeff(t_i, t_e, r_i, r_e, y)
trans_power = abs(E_trans * conj(E_trans))
trans_phase = angle(E_trans, deg=True)

E_refl = refl_coeff(r_i, r_e, y)
refl_power = abs(E_refl * conj(E_refl))
refl_phase = angle(-1*E_refl, deg=True)


fignum=0

matplotlib.rcParams.update({'savefig.dpi':250})

fignum=fignum+1
pylab.figure(fignum)

ax1 = pylab.subplot(211)

pylab.plot(2*pi*x,P_circ,'g-',linewidth=1.4,label='Circulating Power')
pylab.legend(loc=1,prop={'size':10},fancybox=True)
pylab.grid(True,which='both')
pylab.xlim(-pi,pi)
pylab.ylabel(r'$P_{circ}$ (arb)')
pylab.ylabel(r'$P_{circ}/P_0$')
pylab.xticks(visible=False)
pylab.yticks(fontsize=10)
ax1.get_yaxis().set_label_coords(-0.08,0.5)

ax2 = pylab.subplot(212)

pylab.plot(2*pi*x,PD_I,'r-',linewidth=1.2,label='I-phase')
pylab.plot(2*pi*x,PD_Q,'b-',linewidth=1.2,label='Q-phase')
pylab.legend(loc=1,prop={'size':10},fancybox=True)
pylab.grid(True,which='both')
pylab.xlim(-pi,pi)
pylab.xlabel(r'$\phi_f$ (rad)')
pylab.ylabel(r'$S_{refl}/P_0$ (arb)')
pylab.xticks(fontsize=10)
pylab.yticks(fontsize=10)
ax2.get_yaxis().set_label_coords(-0.08,0.5)


######################

# Now begins an extremely hacky way to make two x-axes on the bottom of the plot
# (Making another x-axis at the top of the plot is not so hard, lots of answers via googling)
pyplot.subplots_adjust(bottom=0.25)

# Set up a function to draw another bottom spline
def make_second_bottom_spine(ax=None, label=None, offset=0, labeloffset=25):
    if ax is None:
        ax = pyplot.gca()
    second_bottom = matplotlib.spines.Spine(ax, 'bottom', ax.spines['bottom']._path)
    second_bottom.set_position(('outward', offset))
    ax.spines['second_bottom'] = second_bottom

    if label is not None:
        # Make a new xlabel
        ax.annotate(label, 
                xy=(0.5, 0), xycoords='axes fraction', 
                xytext=(0, -labeloffset), textcoords='offset points', 
                verticalalignment='top', horizontalalignment='center')

# Move the bottom x-axis down, and draw a second spline with the label for the second axis
ax2.spines['bottom'].set_position(('outward', 55))
make_second_bottom_spine(label='cavity displacement [nm]')
pylab.xticks(fontsize=10)


# Make a copy of the x-axis, move the tick marks to the bottom, move the labels to the bottom
ax3 = ax2.twiny()
ax3.xaxis.tick_bottom()
ax3.xaxis.set_label_position('bottom')

# Position and label the tickmarks for the new axis
# note the map() routine - the labels (string values) are given as a list, but we want to do some math on the tick postions which needs an array
x_tick_locations = array([-500,-400,-300,-200,-100,0,100,200,300,400,500])
pylab.xticks(2*pi*x_tick_locations*1e-9/lam, map(str,x_tick_locations.tolist()), fontsize=10)

##### End hacky twin-x-axis method

pylab.savefig('cavity.png',bbox_inches='tight')





fignum=fignum+1
pylab.figure(fignum)

ax1 = pylab.subplot(211)

pylab.plot(y*2*pi,trans_power,'b-',label='Transmitted')
pylab.plot(y*2*pi,refl_power,'r-',label='Reflected')

pylab.grid(True)
pylab.xlim(-1,8)
pylab.ylim(0,1)
pylab.ylabel(r'Power, $P/P_0$')
pylab.legend(loc=1,prop={'size':10},fancybox=True,bbox_to_anchor=(0.65,0.95))
ax1.get_yaxis().set_label_coords(-0.08,0.5)


ax2 = pylab.subplot(212)
pylab.plot(y*2*pi,trans_phase,'b-',label='Transmitted')
pylab.plot(y*2*pi,refl_phase,'r-',label='Reflected')
pylab.grid(True)
pylab.xlim(-1,8)
pylab.ylabel('Phase [deg]')
pylab.xlabel(r'$\phi_f$ (rad)')
ax2.get_yaxis().set_label_coords(-0.08,0.5)

pylab.savefig('cavity_trans.png',bbox_inches='tight')
