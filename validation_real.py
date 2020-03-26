import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad
from glob import glob
from constants import *

files = sorted(set(glob('./PP.CNO_0/*')))

test_data  = np.genfromtxt('hipparcus.dat', skip_header=1)
test_data2 = np.genfromtxt('data.dat')

def func_R(M):
    if M <= 1.66*M_sun:
        return 1.06*(M/M_sun)**(0.945)*R_sun
    else:
        return 1.33*(M/M_sun)**(0.555)*R_sun
    
def func_L(M):
    if M <= 0.7*M_sun:
        return 0.35*(M/M_sun)**(2.62)*L_sun
    else:
        return 1.02*(M/M_sun)**(3.92)*L_sun

def rho(r):
    return rhos - rhos/(Rf-RR)*r

def integrate(r):
    return 4*pi*r**2.*rho(r)

def func_M_f(R_f,R,M,rho_s):
    M_output = []
    for i in range(len(R)):
        global rhos, Rf,RR
        
        rhos = rho_s[i]
        Mm   = M[i]
        Rf   = R_f[i]
        RR   = R[i]
        
        M_to_add = quad(integrate,Rf,RR)[0]
        
        M_output.append(Mm+M_to_add)
        
    return np.array(M_output)
    
test_L_hp = 10**((15-test_data[:,1]-5*np.log10(test_data[:,4]))/2.5)*L_sun
test_T_hp   = 9000/(0.93+test_data[:,8])
#test_T_hp   = 7090/(0.71+test_data[:,8])


prob   = test_data2[:,-1]
test_L = []
test_T = []
for i,p in enumerate(prob):
    test_L.append(10**test_data2[i,7]*L_sun)
    test_T.append(10**test_data2[i,6])

test_L = np.array(test_L)
test_T = np.array(test_T)

fig1,ax1 = plt.subplots(figsize=(8,6))
fig2,ax2 = plt.subplots(figsize=(8,6))
fig3,ax3 = plt.subplots(figsize=(8,6))

#ax1.loglog(test_T,test_L/L_sun,'ko',alpha=0.3,label='Watko, 1997')
#ax1.loglog(test_T_hp,test_L_hp/L_sun,'ro',alpha=0.3, label='ESA, 1997')

for file_name in files:
    contents_MS = np.loadtxt(file_name,skiprows=1,delimiter=',')

    exp = file_name.split('Tpp')[1].split('_')[0]

    R = contents_MS[:,0]
    L = contents_MS[:,1]
    M = contents_MS[:,2]
    T = contents_MS[:,3]
    rho_s = contents_MS[:,4]

    T_f = (L/(4*pi*R**2.*sigma))**(1/4.)
    T_a = (T_f + T)/2.
    R_f = np.sqrt(L/(4*pi*sigma*T_a**4.))
    
    M_f = func_M_f(R_f,R,M,rho_s)

    #M_f = M*(R_f/R)**3.

    mass_range = np.linspace(min(M),max(M),10000)

    #ax1.loglog(T_f,L/L_sun,'o',label='$PP \propto T^{%s}$' % exp)
    ax1.loglog(T_f,L/L_sun,'o',label=r'$\rho_{c,tol}: %s$' % exp)
    #ax1.loglog(T,L,'o')
    
    

    L_range = np.array([func_L(mass) for mass in mass_range])
    R_range = np.array([func_R(mass) for mass in mass_range])

    emp_L = np.array([func_L(m) for m in M]) # empirical shit
    emp_R = np.array([func_R(m) for m in M])
    
    rel_err_L = abs(emp_L - L)/emp_L
    rel_err_R = abs(emp_R - R)/emp_R
    
    rel_err_avg_L = np.mean(rel_err_L)
    rel_err_avg_R = np.mean(rel_err_R)
    
    #ax2.loglog(M/M_sun,L/L_sun,'o',label=r"$PP \propto T^{%s}, \bar\epsilon$:%s%s" %(exp, int(rel_err_avg_L*1000)/10.,'%') )
    ax2.loglog(M/M_sun,L/L_sun,'o',label=r"$\rho_{c,tol}:%s, \bar\epsilon$:%s%s" %(exp, int(rel_err_avg_L*1000)/10.,'%') )
    ax2.loglog(mass_range/M_sun,L_range/L_sun,'k-')
    ax2.loglog([0.7,0.7],[min(L_range)/L_sun,max(L_range)/L_sun],'b-')

    
    #ax3.loglog(M/M_sun,R/R_sun,'o',label=r"$PP \propto T^{%s}, \bar\epsilon$:%s%s" %(exp, int(rel_err_avg_R*1000)/10.,'%') )
    ax3.loglog(M/M_sun,R/R_sun,'o',label=r"$\rho_{c,tol}:%s, \bar\epsilon$:%s%s" %(exp, int(rel_err_avg_R*1000)/10.,'%') )
    ax3.loglog(mass_range/M_sun,R_range/R_sun,'k-')
    ax3.loglog([1.66,1.66],[min(R_range)/R_sun,max(R_range)/R_sun],'b-')


csfont = {'fontname':'Times New Roman'}
#plt.rcParams["font.family"] = "Times New Roman"

ax1.set_xlabel(r'Temperature in K', size=20,**csfont)
ax1.set_ylabel(r'Luminosity in $L_\odot$', size=20,**csfont)
ax1.set_xlim([1610,42625])
ax1.set_ylim([3e-4,2e5])
ax1.invert_xaxis()
ax1.tick_params('both', labelsize=16)
ax1.grid()
fig1.tight_layout()

ax2.set_xlabel(r'Mass in $M_\odot$', size=20,**csfont)
ax2.set_ylabel(r'Luminosity in $L_\odot$', size=20,**csfont)
ax2.tick_params('both', labelsize=16)
ax2.grid()
fig2.tight_layout()


ax3.set_xlabel(r'Mass in $M_\odot$', size=20,**csfont)
ax3.set_ylabel(r'Radius in $R_\odot$', size=20,**csfont)
ax3.tick_params('both', labelsize=16)
ax3.grid()
fig3.tight_layout()

fig1.savefig('HR_PP_validation.png',dpi=200)
fig2.savefig('ML_PP_validation.png',dpi=200)
fig3.savefig('MR_PP_validation.png',dpi=200)

plt.show()

