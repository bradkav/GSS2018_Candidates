#The Hubble Constant
from __future__ import division
import numpy as np
import scipy
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import Symbol
from scipy.integrate import odeint
from scipy.optimize import newton


M_pl = 2.435e18 #GeV
G_F = 1.166e-5

class Freezeout(object):
    def __init__(self,  mass = 1):
        self.mass = mass
        self.Tlist = np.logspace(4, -2, 200)

    def H(self, T, g_star = None):
        if g_star is None:
            g_star = self.g_star(T)
        return np.pi/3 * (g_star/10)**(0.5) * T**2/M_pl

    def g_star_int(self, T):
        '''Takes temperature in GeV, returns g_eff, numbers found from https://arxiv.org/pdf/1609.04979.pdf'''
        if T > 170:
            return 106.25
        if T <= 170 and T > 125: #annihilation of  t, tbar
            return 96.25
        if T <= 125  and T > 91: #annihilation of H0
            return 95.25
        if T <= 91  and T > 80: #annihilation of Z0
            return 92.25
        if T <= 80  and T > 4: #annihilation of W+,W-
            return 86.25
        if T <= 4 and T > 1.8: #annihilation of tau+, tau -
            return 75.75
        if T <= 1.8 and T > 1.2: #annihilation of b bbar
            return 72.25
        if T <= 1.2 and T > 0.2: #annihilation of c, cbar
            return 61.75
        if T <= 0.2 and T > 0.139: #QCD transition
            return 17.25
        if T <= 0.139 and T > 0.131: #annihilation of pi+, pi-
            return 15.25
        if T <= 0.131 and T > 0.106: #annihilation of pi0
            return 14.25
        if T <= 0.106 and T > 0.0008: #annihilation of mu+, mu-
            return 10.75
        if T <= 0.0008 and T > 0.0005: #neutrino decoupling
            return 6.863
        if T <= 0.0005:
            return 3.363

    def g_star(self, Tlist):
        '''Takes temperature in GeV, returns g_eff, numbers found from https://arxiv.org/pdf/1609.04979.pdf'''
        gs = []
        for T in Tlist:
            if T > 170:
                gs.append(106.25)
            if T <= 170 and T > 125: #annihilation of  t, tbar
                gs.append(96.25)
            if T <= 125  and T > 91: #annihilation of H0
                gs.append(95.25)
            if T <= 91  and T > 80: #annihilation of Z0
                gs.append(92.25)
            if T <= 80  and T > 4: #annihilation of W+,W-
                gs.append(86.25)
            if T <= 4 and T > 1.8: #annihilation of tau+, tau -
                gs.append(75.75)
            if T <= 1.8 and T > 1.2: #annihilation of b bbar
                gs.append(72.25)
            if T <= 1.2 and T > 0.2: #annihilation of c, cbar
                gs.append(61.75)
            if T <= 0.2 and T > 0.139: #QCD transition
                gs.append(17.25)
            if T <= 0.139 and T > 0.131: #annihilation of pi+, pi-
                gs.append(15.25)
            if T <= 0.131 and T > 0.106: #annihilation of pi0
                gs.append(14.25)
            if T <= 0.106 and T > 0.0008: #annihilation of mu+, mu-
                gs.append(10.75)
            if T <= 0.0008 and T > 0.0005: #neutrino decoupling
                gs.append(6.863)
            if T <= 0.0005:
                gs.append(3.363)
        return gs

    def g_starDM(self, mass, dof, Tlist):
        '''takes a descending Tlist and returns geff including extra dof'''
        geff = []
        l = 0
        for T in self.Tlist:
            if T > mass:
                l +=1
        for n in range(0, l):
            geff = self.g_star(T[n]) + dof
        for m in range(l, len(T)):
            geff = self.g_star(T[n])
        return geff

    def number_density(self, mass, T, limit = 'rel', g_star = None):
        if g_star is None:
            g_star = self.g_star(T)
        '''Number density at equilibrium'''
        # non-relativistic m_i >> T
        if limit != 'rel':
            return g_star * (mass * T / (2*np.pi))**(3/2) * np.exp(-mass/T)
        # relativistic m_i << T
        if limit == 'rel':
            return (g_star * (T)**3)/np.pi**2

    def dYdx(self, Y, x, sigv):
        '''The differential equation for dY/dx, Y = nx/T^3, x = m/T'''
        return ((self.mass)**3 * sigv / self.H(self.mass, g_star = self.g_star_int(self.mass)) / x**2) * (Y**2 - (self.number_density(self.mass/x, g_star=self.g_star_int(self.mass/x))/(self.mass/x)**3)**2)

    def getY(self, xlist = None):
        '''Solves dYdx'''
        if xlist is None:
            xlist = self.mass / self.Tlist
        T0 = 10000
        Y0 = self.number_density(T0, g_star = self.g_star_int(T0)) / T0**3
        Y = odeint(self.dYdx, Y0, xlist, args = (1e-20,))[:,0]
        return lambda x: np.interp(x, xlist, Y)

    def freezeNonRel(self,x, mass,mz):
        return np.sqrt(x)*np.exp(-x)  - 1/(mass * M_pl * self.crossSecNonRel(mass, mz))  #* self.g_star_int(self.mass/x)*10 * (3 * self.mass * M_pl*sigv)**2 / ((np.pi**(5)) * (2**(3)))

    def xFreezeRel(self, mass):
        x = []
        for m in mass:
            T = (np.pi**3 / 3 * np.sqrt(1/10*106.25) * 1/(self.crossSecRel(m)*(np.sqrt(3 *8.61e-14/m))*M_pl))**(2/3)
            x.append(m / T)
        return x

    def xFreezeNonRel(self, mass, mz):
        x0 = 10
        x = []
        for m in mass:
            x.append(newton(self.freezeNonRel, x0, args = (m,mz), maxiter = 10000))
        return x

    def RelicDensityNonRel(self, mass, mz):
        pc = 1.055e-5 /(5.06e13**3) #GeV
        T0 = 2.75e-13 #GeV
        xf = self.xFreezeNonRel(mass, mz)
        for x in xf:
            r = (x * T0**3) / (M_pl * self.crossSecNonRel(mass, mz) * pc)
        return r

    def RelicDensityRel(self, mass):
        pc = 1.055e-5 /(5.06e13**3) #GeV
        xf = self.xFreezeRel(mass)
        for x in xf:
            r = (mass * self.number_density(mass, mass/x)) / (pc)
        return r

    def crossSecNonRel(self, mass, mz):
        return mass**2 / (((mass+mz)**2 - mz**2)**2 + mz**4)

    def crossSecRel(self, mass):
        return mass**(-2)

F2 = Freezeout(mass = 100)
mass = np.logspace(-1, 5, 100)
plt.loglog(mass, F2.RelicDensityNonRel(mass, 10))
plt.loglog(mass, F2.RelicDensityNonRel(mass, 100))
plt.loglog(mass, F2.RelicDensityNonRel(mass, 1000))
mass = np.logspace(-5.5, -3.5, 100)
plt.loglog(mass, F2.RelicDensityRel(mass))
#plt.xlim()
#mDM = np.logspace(-7, 3)
#plt.semilogx(mDM, F2.RelicDensity(F2.crossSecNonRel(mDM, 10)))
#plt.semilogx(F2.mass/F2.Tlist, F2.getY()(F2.mass/F2.Tlist), label = 'm = 10')
# plt.xlabel("$x = m/T$/ GeV")
# plt.ylabel("$Y = n_x / T^3$")
# plt.legend()
plt.show()
