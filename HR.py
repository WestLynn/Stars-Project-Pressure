import pickle
import matplotlib.pyplot as plt
import math as m
import numpy as np
from scipy.interpolate import spline
lum = []
temp = []
mbol = []
mass = []
radius =[]
coeff = [1.07e-4, 1.07e-5, 1.07e-6, 1.07e-7, 1.07e-8, 1.07e-9, 1.07e-10]
L_sun = 3.828e26
mbol_sun = 4.74
M_sun = 1.989e30
R_sun = 7e8
with open('stars_ppcoeff.pkl', 'rb') as of:
    d = pickle.load(of)
    
def plotHR(c): 
    lum = []
    temp = []
    mbol = []
    mass = []
    radius =[]
    for index in range(len(d[c])):
        if d[c][index] is not None:
            lum.append((d[c][index].L[-1]) / L_sun)
            temp.append(d[c][index].T[-1])
            mass.append(d[c][index].M[-1] / M_sun)
            radius.append(d[c][index].r[-1] / R_sun)
            mbol.append(-2.5 * m.log((d[c][index].L[-1] / L_sun), 10) + mbol_sun)
    
    
    return lum, temp, mbol, mass, radius
    
full = []
for c in coeff:
    full.append(plotHR(c))
#HR Diagram    
    
fig, ax1 = plt.subplots()
ax1.set_ylabel('$L/L_{\odot}$')
ax1.set_xlabel('Temperature')
ax1.set_xlim(10000, 1000)
ax1.set_ylim(0.001, 600)
ax1.yaxis.set_label_position('right')
ax1.xaxis.set_label_position('top')
#xnew = np.linspace(min(full[0][1]), max(full[0][1]), 300)
#spl = make_interp_spline(full[0][1], full[0][0], k=3)
#power_smooth = spline(full[0][1], full[0][0], xnew)
#ax1.plot(xnew,power_smooth)
ax1.plot(full[0][1],full[0][0], 'o', label = '1.07e-4')
ax1.plot(full[1][1],full[1][0], 'o', label = '1.07e-5')
ax1.plot(full[2][1],full[2][0], 'o', label = '1.07e-6')
ax1.plot(full[3][1],full[3][0], 'o', label = '1.07e-7')
ax1.plot(full[4][1],full[4][0], 'o', label = '1.07e-8')
ax1.plot(full[5][1],full[5][0], 'o', label = '1.07e-9')
ax1.plot(full[6][1],full[6][0], 'o', label = '1.07e-10')
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend(loc = 'best')
 
ax2 = ax1.twiny()
ax2.set_xlabel('Spectral Class')
ax2.xaxis.set_label_position('bottom')
ax2.xaxis.tick_bottom()
ax2.set_xscale('log')
ax2.tick_params(axis=u'both', which=u'both',length=0)
ax2.set_xticks([40000, 20000, 10000, 7500, 5500, 4000, 3000, 2000, 1300])
ax2.set_xlim(10000, 1000)
ax2.set_xticklabels(['O','B','A','F','G','K','M','L','T'])
ax1.xaxis.tick_top()
ax3 = ax1.twinx()
ax3.set_ylabel('$M_{bol}$')
ax3.yaxis.set_label_position('left')
ax3.yaxis.tick_left()
#ax3.plot(temp,mbol, 'ok')
ax3.set_ylim(18, -9)
ax3.set_xlim(10000, 1000)
ax1.yaxis.tick_right()
#L/Lsun versus M/Msun
fig, ax1 = plt.subplots()
ax1.set_ylabel('$L/L_{\odot}$')
ax1.set_xlabel('$M/M_{\odot}$')
ax1.plot(full[0][3],full[0][0], 'o', label = '1.07e-4')
ax1.plot(full[1][3],full[1][0], 'o', label = '1.07e-5')
ax1.plot(full[2][3],full[2][0], 'o', label = '1.07e-6')
ax1.plot(full[3][3],full[3][0], 'o', label = '1.07e-7')
ax1.plot(full[4][3],full[4][0], 'o', label = '1.07e-8')
ax1.plot(full[5][3],full[5][0], 'o', label = '1.07e-9')
ax1.plot(full[6][3],full[6][0], 'o', label = '1.07e-10')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(0.02, 12)
ax1.set_ylim(1e-3, 7e4)
ax1.legend(loc = 'best')
#R/Rsun versus M/Msun
    
fig, ax1 = plt.subplots()
ax1.set_ylabel('$R/R_{\odot}$')
ax1.set_xlabel('$M/M_{\odot}$')
ax1.plot(full[0][3],full[0][4], 'o', label = '1.07e-4')
ax1.plot(full[1][3],full[1][4], 'o', label = '1.07e-5')
ax1.plot(full[2][3],full[2][4], 'o', label = '1.07e-6')
ax1.plot(full[3][3],full[3][4], 'o', label = '1.07e-7')
ax1.plot(full[4][3],full[4][4], 'o', label = '1.07e-8')
ax1.plot(full[5][3],full[5][4], 'o', label = '1.07e-9')
ax1.plot(full[6][3],full[6][4], 'o', label = '1.07e-10')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(0.02, 12)
ax1.set_ylim(1e-1, 10)
ax1.legend(loc = 'best')

