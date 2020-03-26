import numpy as np, matplotlib.pyplot as plt, time, math, sys
import scipy.integrate as scipy

######################## Constants for Stars Project ###############################
# all constants are in basic SI units

sigma = 5.670373e-8      # Stefan-Boltzmann constant
k = 1.3806488e-23         # Boltzmann constant
c = 2.99792458e8           # speed of light
me = 9.10938291e-31        # mass of electron
mp = 1.67262178e-27        # mass of proton
hbar = 1.05457173e-34      # h/2pi
G = 6.67384e-11         # Gravitation constant
a = (4.0*sigma)/c     # radiation constant
X = 0.70              # Hydrogen mass fraction (Sun)
Y = 0.28              # Helium mass fraction (Sun)
Z = 0.02              # Metals mass fraction (Sun)
mu = 0.617283951         # mean molecular weight (100% ionized, Sun's metallicity)
gamma = 5.0/3.0       # adiabatic constant
Kappa_es = 0.02*(1.0 + X) # [m^2/kg]  Opacity from electron scattering

#Constants from the example star
rho_c = 58560.0
T_c = 8.23e6
M_sun = 1.989e30
R_sun = 6.955e8
L_sun = 3.846e26
M_star = 0.673*M_sun
R_star = 0.865*R_sun
L_star = 5.86e-2*L_sun

class trialStar: #Class to be used to create one star

    def __init__(self,dr,central_density,central_temp,p):
        """
        dr = step size to use for numerical integration
        p = a boolean array indicating which pressures we want on or off
        """
        self.individual_pressures = np.zeros(3)
        self.dP_drho_values = np.zeros(3)
        self.dP_dT_values = np.zeros(3)

        self.dr = dr # step size
        #self.accuracy = 1.0e-7*np.array([central_density/1.0e8,central_temp/1.0e8,1.0e29,1.0e24])
        self.accuracy = 1.5e-16#*np.array([4.38e-4,2.73e-2,5.338e21,1.572e18])
        self.a = 0
        
        self.dr_values = [self.dr]
        #A bunch of lists, with values for each attribute at points between the center and surface of the star
        self.control = np.array(p)
        
        self.radius = [10.0]
        self.density = [central_density]
        self.temp = [central_temp]
        self.mass = [(4.0/3.0)*np.pi*self.radius[0]**3.0*self.density[0]]
        self.luminosity = [(4.0/3.0)*np.pi*self.radius[0]**3.0*self.density[0]*self.epsilon(self.density[0],self.temp[0])]
        self.opacity = [self.Kappa(self.density[0], self.temp[0])]
        self.tau = [self.Kappa(self.density[0], self.temp[0])*self.density[0]]
        self.pressure = [self.P(self.density[0], self.temp[0])]
        self.dL_dr_values = [0.0]
        
        self.hminus = [self.Kappa_Hminus(self.density[0],self.temp[0])]
        self.ff = [self.Kappa_ff(self.density[0],self.temp[0])]
        
        self.convective = [self.dT_dr_convective(self.mass[0],self.density[0],self.radius[0],self.temp[0],self.luminosity[0])]
        self.radiative = [self.dT_dr_radiative(self.mass[0],self.density[0],self.radius[0],self.temp[0],self.luminosity[0])]
        
        self.dTau = [self.Kappa(self.density[0], self.temp[0]) * self.density[0]**2 / abs(self.drho_dr(self.mass[0],self.density[0],self.radius[0],self.temp[0],self.luminosity[0]))]
        
        # David: I'm adding two new lists to hold dT_dr for convective and radiative in order to determine
        # which one to use
        self.derivative_radiative = [self.dT_dr_radiative(self.mass[0],self.density[0],self.radius[0],self.temp[0],self.luminosity[0])]
        self.derivative_convective = [self.dT_dr_convective(self.mass[0],self.density[0],self.radius[0],self.temp[0],self.luminosity[0])]
        
        self.Generate()
        self.a = self.surfaceRad()
        
####Runge-Kutta Witch-craft happens here: #########################################

    def Generate(self):
        """
        Loops through the Runge-Kutta until opacity condition is met.
        """
        i = 0
        stop_Integration = False # always run through the loop at least once
        while stop_Integration == False:
            rk4_input = np.array([self.density[-1], self.temp[-1], self.mass[-1], self.luminosity[-1], self.tau[-1]],float)
            #print 'rk4 call', self.dr
            rk4_output, new_radius = self.AdaptiveRK4(rk4_input, self.radius[-1], self.dr, self.f)
            
            # adding new values
            self.UpdateValues(rk4_output,new_radius)
            # check termination condition
            stop_Integration = self.opacityCheck()
            i += 1
            #if i > 100:
                #break
            if self.mass[-1] >= 100.0*M_sun:
                break

    def AdaptiveRK4(self,s,r,h,f):
        """
        Implements adaptive step-sizes for our Runge-Kutta method.
        """
        rho_error = 0.
        increase_h = True
        hmax = 7.0e7 # 0.1% of the suns radius
        i = 0
        hmin = 100.0

        while rho_error <= 1.:
            
            #print rho_error
            xsmall_1, new_radius = self.rk4(s,r,self.dr,f)
            xsmall_2, new_radius2 = self.rk4(xsmall_1,r+self.dr,self.dr,f)
            xbig, notused_radius = self.rk4(s,r,2.*self.dr, f)

            """
            Be careful: Probably should check if dtau is nan and recalculate steps if it is!
            """
            # accuracy for the surface
            #accuracy2 = 1.0e-11
            top = 30. * self.dr * self.accuracy
            #bottom = np.sqrt(np.sum(np.power((xbig - xsmall_2),2)))
            bottom = max(abs(xbig-xsmall_2)/xsmall_2) + 1.0e-15
            #print bottom
            #print bottom
            #print self.dr
            #print top
            
            #if min(top/bottom[:4]) == 0:
            #if min(bottom[:4]) == 0:
            #    #bottom[:4] = max(bottom[:4])
            #    rho_error = 10.0
            #else:
            rho_error = top / bottom
            if math.isnan(rho_error):
                rho_error = 0.9
            #print rho_error
            #print self.dr
            # If temp is getting close to the surface, change hmax to 10km
            #if self.temp[-1]/self.temp[0] < 0.008*(self.temp[0]/8.23e6):
            #    hmax = 1.0e4*(self.temp[0]/8.23e6)
            #    if self.dr > hmax:
            #        print self.radius[-1]/R_star
            #        print 'maxed'
            #        self.dr = hmax
            #        increase_h = False
            #    
            #    top = 30.*self.dr*accuracy2
            #    bottom = abs(xbig[1]-xsmall_2[1])/xsmall_2[1]
            #    rho_error = top/bottom
            #    #print rho_error
                #print self.dr
                
                #i = 1
            
            # i is a holder for a condition to override the rho_error. i = 1 when self.dr gets larger than hmax
            if rho_error > 1. or i == 1:
                if increase_h == True: self.dr *= 2.0 #Increase step-size only by a factor of two (just to be careful)!
                if self.dr > hmax:
                    self.dr = hmax
                    increase_h = False
                    i = 1
                    continue
                    
                #Keep the values from the itty-bitty steps (to save computation time):
                if rho_error < 2.0 or i == 1:
                    #self.UpdateValues(xsmall_1,new_radius)
                    #print self.dr
                    break
                else: 
                    #print 'hello'
                    continue
            else:
                self.dr *= 1.0 / 2.0
                if self.dr < hmin:
                    self.dr = hmin
                    i = 1
                #print "decrease", rho_error, "top", top, "bottom", bottom
                    increase_h = False

        #print "Step size:", self.dr, "  Iteration:", self.iterations, "Radius:", self.radius[-1]
        return xsmall_2, new_radius2

    def rk4(self,s,r,h,f):
        """
        Suggestion:
        s = self.RK_variables
        r = self.radius[-1]
        h = self.dr
        """
        """
        s       - is the vector of dependent variable(s)
        e.g s = [0,1], s[0] represents x in this example, s[1] represents y
        r       - is the independent variable (radius for the final project)
        h       - is the step size
        f       - is the function that the program calls to runge-kutta
        e.g the function will be of form f(s,t).
        """
        """
        With vector form of s, multiple variables can be hidden and everything
        can be solved step-wise at the same time
        """
        # fourth order runge-kutta
        k1 = h*f(s, r)
        k2 = h*f(s + 0.5*k1, r + 0.5*h)
        k3 = h*f(s + 0.5*k2, r + 0.5*h)
        k4 = h*f(s + k3, r + h)
        return s + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0, r+h

    def f(self,RK_variables,r):
        # Variables for function calls
        density = RK_variables[0]
        temp = RK_variables[1] 
        mass = RK_variables[2]
        luminosity = RK_variables[3]
        #tau = RK_variables[4]
        radius = r
        #pressure = self.P(density,temp)
        
        frho = self.drho_dr(mass, density, radius, temp, luminosity)
        fT = self.dT_dr(mass, density, radius, temp, luminosity)
        fM = self.dM_dr(density, radius)
        fL = self.dL_dr(density, radius, temp)
        ftau = self.dtau_dr(density, temp)
        return np.array([frho, fT, fM, fL, ftau],float)

    def UpdateValues(self,rk4_output,new_radius):
        #Update our arrays for our star's properties.
        self.dL_dr_values.append(self.dL_dr(rk4_output[0],new_radius,rk4_output[1]))
        
        self.density.append(rk4_output[0])
        self.temp.append(rk4_output[1])
        self.mass.append(rk4_output[2])
        self.luminosity.append(rk4_output[3])
        
        self.dr_values.append(self.dr)
        #if self.dTau[-1] < 100.0:
        self.tau.append(rk4_output[4])
        self.ff.append(self.Kappa_ff(self.density[-1],self.temp[-1]))
        #else:
            #self.tau.append(0)
        
        self.hminus.append(self.Kappa_Hminus(self.density[-1],self.temp[-1]))
        self.convective.append(self.dT_dr_convective(self.mass[-1],self.density[-1],self.radius[-1],self.temp[-1],self.luminosity[-1]))
        self.radiative.append(self.dT_dr_radiative(self.mass[-1],self.density[-1],self.radius[-1],self.temp[-1],self.luminosity[-1]))
        
        self.pressure.append(self.P(self.density[-1], self.temp[-1]))
        self.opacity.append(self.Kappa(self.density[-1],self.temp[-1]))
        self.radius.append(new_radius)
        
        # David: Updating lists of derivative of dT_dr for convective and radiative (must be after new lists have new values)
        # Would appending be quicker then resetting self.derivative_radiative[0] to a different value?
        self.derivative_radiative.append(abs(self.dT_dr_radiative(self.mass[-2],self.density[-2],self.radius[-2],self.temp[-2],self.luminosity[-2])-self.dT_dr_radiative(self.mass[-1],self.density[-1],self.radius[-1],self.temp[-1],self.luminosity[-1])))
        self.derivative_convective.append(abs(self.dT_dr_convective(self.mass[-2],self.density[-2],self.radius[-2],self.temp[-2],self.luminosity[-2])-self.dT_dr_convective(self.mass[-1],self.density[-1],self.radius[-1],self.temp[-1],self.luminosity[-1])))
    
    #code that checks boundary conditions####################
    def opacityCheck(self):
        #This returns a True or False statement about whether we should continue to integrate outward or not.
        dtau = self.Kappa(self.density[-1], self.temp[-1]) * self.density[-1]**2 / abs(self.drho_dr(self.mass[-1],self.density[-1],self.radius[-1],self.temp[-1],self.luminosity[-1]))
        #print dtau
        self.dTau.append(dtau)
        
        if dtau < 0.001:
            opacityCondition = True

        else:
            opacityCondition = False
        
        if self.radius[-1] > 1.E10: opacityCondition = True #Just in case it wants to integrate forever, we'll cut it off.
        
        return opacityCondition

    #finds the radius that corresponds to tau_surface
    def surfaceRad(self):
        # integrate dtau
        #self.tau = np.array(self.dr_values) * np.cumsum(np.array(self.density)*np.array(self.opacity))
        self.a = np.argmin(abs(self.tau[-1]-np.array(self.tau) - (2.0/3.0)))
        if self.tau[self.a] == 0:
            self.a = len(self.tau)-1
        print self.a
        return self.a
    #########################################################

#####Pressure Functions: ###########################################################
    def P(self,density,temp):
        if self.control[0] != 0:
            self.individual_pressures[0] = self.pressure1(density)

        if self.control[1] != 0:
            self.individual_pressures[1] = self.pressure2(density,temp)

        if self.control[2] != 0:
            self.individual_pressures[2] = self.pressure3(temp)

        return np.sum(self.control*self.individual_pressures)

    def dP_drho(self,density,temp):
        if self.control[0] != 0:
            self.dP_drho_values[0] = self.dP_drho_1(density)

        if self.control[1] != 0:
            self.dP_drho_values[1] = self.dP_drho_2(temp)

        return np.sum(self.control*self.dP_drho_values)

    def dP_dT(self,density,temp):
        if self.control[1] != 0:
            self.dP_dT_values[1] = self.dP_dT_2(density)

        if self.control[2] != 0:
            self.dP_dT_values[2] = self.dP_dT_3(temp)

        return np.sum(self.control*self.dP_dT_values)


#####Individual Pressure and Partial Pressure Functions ###########################

    # Degeneracy pressure
    def pressure1(self,density):
        #David Stevens
        return ((3.0*(np.pi)**2.0)**(2.0/3.0)*(hbar**2.0)*(density/mp)**(5.0/3.0))/(5.0*me)

    # Degeneracy part of differential
    def dP_drho_1(self,density):
        #David Stevens
        return ((3.0*(np.pi)**2.0)**(2.0/3.0)*(hbar**2.0)*(density/mp)**(2.0/3.0))/(3.0*me*mp)

    # Ideal gas part of differential
    def dP_drho_2(self,temp):
        #David Stevens
        return (k*temp)/(mu*mp)
        
    # Ideal gas pressure
    def pressure2(self,density,temp):
        #David Stevens
        return (density*k*temp)/(mu*mp)

    # Radiative pressure
    def pressure3(self,temp):
        #David Stevens
        return (a*(temp)**4.0)/(3.0)

    # Ideal gas part of differential
    def dP_dT_2(self,density):
        #David Stevens
        return (density*k)/(mu*mp)

    # Radiative part of differential
    def dP_dT_3(self,temp):
        #David Stevens
        return (4.0*a*temp**3.0)/(3.0)


#####Functions for the five equations of stellar structure: ########################
    def drho_dr(self, mass, density, radius, temp, luminosity):
        #David Stevens
        ########### need appropriate dT_dr ################
        return -((G*mass*density)/(radius**2.0)+self.dP_dT(density,temp)*self.dT_dr(mass,density,radius,temp,luminosity))/self.dP_drho(density,temp)

    # Returns the current rate of change of temperature (this is chosen by which has the minimum derivative)
    def dT_dr(self, mass, density, radius, temp, luminosity):
        #Joey & Rufus
        #if min(abs(self.derivative_radiative[-1]),abs(self.derivative_convective[-1])) is self.derivative_radiative[-1]:
        #    return self.dT_dr_radiative(mass, density, radius, temp, luminosity)
        #else:
        #    return self.dT_dr_convective(mass, density, radius, temp, luminosity)
        return -min(abs(self.dT_dr_radiative(mass, density, radius, temp, luminosity)),abs(self.dT_dr_convective(mass,density,radius,temp,luminosity)))
    
    # Defining new dT_dr functions that just return values for radiative and convective
    #######################################################################################
    def dT_dr_radiative(self, mass, density, radius, temp, luminosity):
        T_radiative = (3.0*self.Kappa(density,temp)*density*luminosity)/(16.0*np.pi*4.0*sigma*temp**3.0 *radius**2.0)
        return -T_radiative
        
    def dT_dr_convective(self, mass, density, radius, temp, luminosity):
        pressure = self.P(density, temp)
        T_convective = (1.0 - (1.0 / gamma)) * (temp / pressure) * ((G * mass * density) / (radius**2.0))
        return -T_convective
    #######################################################################################
    
    def dM_dr(self, density, radius):
        #Adam G.
        return 4.0 * np.pi * (radius**2.0) * density

    def dL_dr(self,density,radius,temp):
        #Jamie
        return 4.0 * np.pi * (radius**2.0) * density * self.epsilon(density,temp)

    def dtau_dr(self,density,temp):
        #Adam H.
        return self.Kappa(density,temp) * density


####Kappa-specific functions (note: K_es is a pre-defined constant) ################

    def Kappa(self,density,temp):
        #Adam H.
        #if self.temp[-1] > 1.0e4:
        Kappa = ((1.0 / self.Kappa_Hminus(density, temp)) + (1.0 / max(Kappa_es, self.Kappa_ff(density,temp))))**-1.0
        #else:
            #Kappa = ((1.0/self.Kappa_Hminus(density,temp)) + (1.0/min(self.Kappa_ff(density,temp), Kappa_es)))**-1.0
        return Kappa

    def Kappa_ff(self,density,temp):
        #Adam H.
        return (1.0e24) * (1.0+X) * (Z + 0.0001)* ((density/1.0e3)**0.7) * (temp**(-3.5))

    def Kappa_Hminus(self,density,temp):
        #Adam H.
        return (2.5e-32) * (Z / 0.02) * ((density/1.0e3)**0.5) * ((temp)**9.0)


#####Energy-production equations (epsilon): #######################################
    def epsilon(self,density,temp):
        #Jamie
        return self.epsilon_PP(density,temp) + self.epsilon_CNO(density,temp)

    def epsilon_PP(self,density,temp):
        # Jamie
        return 1.07e-7*(density/1.0e5) * (X**2.0) * ((temp/1.0e6)**4.0)

    def epsilon_CNO(self,density,temp):
        #Jamie
        return 8.24e-26 * (density/1.0e5) * 0.03 * (X**2.0) * ((temp / (1.0e6))**19.9)
        
class Star:
    
    def __init__(self,dr,central_temp,p):
        self.dr = dr
        self.central_temp = central_temp
        self.p = p

        self.star_a = trialStar(1000.0e3,0.3e3,central_temp,p)
        self.star_b = trialStar(1000.0e3,500.0e3,central_temp,p)
        self.star_c = trialStar(1000.0e3,(0.3e3+500.0e3)/2.0,central_temp,p)

        
        self.final_star = self.bisection(self.star_a, self.star_c, self.star_b, 1.0e-9)
    
    ####Step-3 Bisection Method: ######################################################
    ## Bisection method
    ## def bisection(a,b,tol)
    ## Usage: applies the bisection method to find the root of a 
    ## previously defined function f (which we will define as the 
    ## function of the central pressure). This method takes two 
    ## real values as the interval where the root exists and a tolerance
    ## which determines how accurate we want to find the root
    
    def bisection(self,star_a,star_c,star_b,tol):
        #Jamie
        i = 0
        while abs((star_b.density[0] - star_a.density[0])) > tol:
            # you might want f(c) < threshold, where threshold is a small number
            # instead of checking equality to zero
            i += 1
            print "delta rho", star_b.density[0] - star_a.density[0] ,"A rho", star_a.density[0], "B rho", star_b.density[0]
            #print star_a.temp[star_a.a]
            #print star_b.temp[star_b.a]
            #print i
            #if self.f_bisect(self.star_c) == 0:
            #    return self.star_c
            if self.f_bisect(star_a) * self.f_bisect(star_c) < 0:
                star_b = star_c
            else :
                star_a = star_c
            if i > 80:
                break
            star_c_density = (star_a.density[0] + star_b.density[0]) / 2.0
            star_c = trialStar(self.dr,star_c_density,self.central_temp,self.p)
        return star_c
        
    def f_bisect(self,trial_star):
        #Jamie
        # Will write the rest of this function once integrating for L_star = L(R_star)
        # and T_star = T(R_star) is done.
        a = trial_star.a
        
        top = trial_star.luminosity[a]-(4.0*np.pi*sigma*(trial_star.radius[a])**2.0*(trial_star.temp[a])**4.0)
        bottom = np.sqrt(4.0*np.pi*sigma*(trial_star.radius[a])**2.0*(trial_star.temp[a])**4.0*trial_star.luminosity[a])
        #print top
        return (top/bottom)
    
    

#test = trialStar(100.0e3,103748.9,1.7e7,[1.0,1.0,1.0])
test_star = Star(100.0e3,2.3e7,[1.0,1.0,1.0])
test = test_star.final_star


###### Looping for a bunch of stars ############
#stars = []
#central_temperatures = []
#surface_temperatures = []
#stars_luminosities = []
#star_number = 0
#for i in np.arange(6.5, 7.3, 0.2):
#        star_number += 1
#        print "Star: ", star_number, "/", int((7.05-6.8)/0.05)+1
#        t = float(10**(i))
#        central_temperatures.append(t)
#        star = Star(100.0e3,t,[1.0,1.0,1.0]).final_star
#        stars.append(star)
#        a = star.a
#        if star.luminosity[a] > 0.0:
#            stars_luminosities.append(star.luminosity[a])
#        else:
#            stars_luminosities.append(0.001)
#
#        surface_temp = star.temp[a]
#        surface_temperatures.append(star.temp[a])
#        print "a", a, "temp", surface_temp, "mass:", star.mass[a]
#
#
#fig, ax = plt.subplots()
#ax.set_xscale('log')
#ax.set_yscale('log')
#plt.plot(surface_temperatures, stars_luminosities)
#plt.gca().invert_xaxis()
#plt.show()


#central_temperatures = []
#surface_temperatures = []
#stars_luminosities = []
#for i in np.arange(6.6, 7.5, 0.1):
#        print "Star: ", int((i-6.6)/0.1)+1, "/", int((7.5-6.6)/0.1)
#        t = 10**(i)
#        central_temperatures.append(t)
#        star = Star(100.0e3,t,[1.0,1.0,1.0]).final_star
#        a = star.a
#        print "a", a, "temp", t
#        if star.luminosity[a] > 0.0:
#            stars_luminosities.append(star.luminosity[a])
#        else:
#            stars_luminosities.append(0.001)
#        surface_temperatures.append(star.temp[a])
#
#
#fig, ax = plt.subplots()
#ax.set_xscale('log')
#ax.set_yscale('log')
#plt.plot(surface_temperatures, stars_luminosities)
#plt.gca().invert_xaxis()
plt.show()

def Plots():
    plt.figure('Radius')
    plt.plot(test.radius, 'o')
    plt.figure('Temperature')
    plt.plot(np.array(test.radius)/R_sun, np.array(test.temp)/test.temp[0],'r-')
    #plt.figure('Density')
    plt.plot(np.array(test.radius)/R_sun, np.array(test.density)/test.density[0],'k-')
    #plt.figure('Mass')
    #plt.plot(np.array(test.radius)/R_sun, np.array(test.mass)/M_sun,'g-')
    #plt.figure('Luminosity')
    #plt.plot(np.array(test.radius)/R_sun, np.array(test.luminosity)/L_sun,'b-')
    plt.figure('Pressure')
    plt.plot(np.array(test.radius)/R_sun, np.array(test.pressure)/test.pressure[0])
    #plt.figure('dL/dr')
    #plt.plot(np.array(test.radius)/R_sun, np.array(test.dL_dr_values)*(R_star/L_star))
    plt.figure('Opacity')
    plt.plot(np.array(test.radius)/R_star, np.log10(test.opacity))
    
    #plt.figure('rho_changes')
    #plt.plot(np.array(test.radius)/R_star,np.array(test.drho_P),'r-')
    #plt.plot(np.array(test.radius)/R_star,np.array(test.dt_P),'b-')
    #
    #plt.figure('drho')
    #plt.plot(np.array(test.radius)/R_star,np.array(test.drho))
    #plt.ylim([-0.001,0])
    #plt.figure('dT_dr')
    #plt.plot(np.array(test.radius)/R_star,test.convective,'b-')
    #plt.plot(np.array(test.radius)/R_star,test.radiative,'k-')
    #plt.ylim([-0.01,0.01])
    #plt.xlim([0,0.15])
    plt.show()

Plots() #unhashtag Plots() to see the plots.
#print (np.array(test.ff)**-1.0 + np.array(test.hminus)**-1.0)**-1.0
print "Surface Temperature", test.temp[test.a]
print "Recalculated Temperature", (test.luminosity[test.a]/(4.0 * np.pi * sigma * test.radius[test.a]**2))**(1.0/4.0)
print "Luminosity", test.luminosity[test.a]
print "Temperature^4", 4.0 * np.pi * sigma * test.radius[test.a]**2.0 * test.temp[test.a]**4.0
print "radius", test.radius[test.a]
print "mass", test.mass[test.a]
print "density", test.density[test.a]
