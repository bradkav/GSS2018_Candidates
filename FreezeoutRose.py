#The Hubble Constant
from __future__ import division
import numpy as np
import scipy
import matplotlib.pyplot as plt
from sympy.solvers import solve
from sympy import Symbol
from scipy.integrate import odeint


M_pl = 2.435e18 #GeV/c^2


class Freezeout(object):
    def __init__(self,  mass = 1):
        self.mass = mass
        self.Tlist = np.linspace(200, 1, 200)

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

    def number_density(self, T, limit = 'rel', g_star = None):
        if g_star is None:
            g_star = self.g_star(T)
        '''Number density at equilibrium'''
        # non-relativistic m_i >> T
        if limit != 'rel':
            return g_star * (self.mass * T / (2*np.pi))**(3/2) * np.exp(-self.mass/T)
        # relativistic m_i << T
        if limit == 'rel':
            return (g_star * (T)**3)/np.pi**2

    def dYdx(self, Y, x, sigv):
        '''The differential equation for dY/dx, Y = nx/T^3, x = m/T'''
        return ((self.mass)**3 * sigv / self.H(self.mass, g_star = self.g_star_int(self.mass)) / x**2) * (Y - self.number_density(self.mass/x, g_star = self.g_star_int(self.mass/x)) / (self.mass/x)**3)


    def getY(self, xlist = None):
        '''Solves dYdx'''
        if xlist is None:
            xlist = self.mass / self.Tlist
        T0 = 1000
        Y0 = self.number_density(T0, g_star = self.g_star_int(T0)) / T0**3
        print(Y0)
        Y = odeint(self.dYdx, Y0, xlist, args = (1e-20,))[:,0]
        #return lambda x: np.interp(x, xlist, Y)
        return Y

    def xfreeze(self, sigv):
        Yinf = 1 #Y at very late times,
        return ((self.mass)**3 * sigv / H(self.mass)) * Yinf

    def Tfreeze(self, sigv):
        T = Symbol('T')
        return solve(self.H(T, 106.25) == self.number_density(T, g_star = 106.25) * sigv, T)

F = Freezeout()
plt.plot(F.mass/ F.Tlist, F.getY())
plt.show()
