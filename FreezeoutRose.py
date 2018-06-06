#The Hubble Constant

import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt


M_pl = 2.435e18 #GeV/c^2


class Freezeout(object):
    def __init__(self,  Tlist, mass = 1):
        self.mass = mass
        self.Tlist = Tlist

    def Ht(self, T):
        return np.pi/3 * (g_star(T)/10)**(0.5) * T**2/M_pl

    def g_star(self, T):
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

    def g_starDM(self, mass, dof, Tlist):
        '''takes a descending Tlist and returns geff including extra dof'''
        geff = []
        l = 0
        for T in self.Tlist:
            if T > mass
            l +=1
        for n in range(0, l):
            geff = g_star(T[n]) + dof
        for m in range(l, len(T)):
            geff = g_star(T[n])
        return geff

    def number_density(self, T, limit = 'rel'):
        '''Number density at equilibrium'''
        # non-relativistic m_i >> T
        if limit == 'nonrel' :
            n_i = g_i * (m_i * T / (2*np.pi))**(3/2) * np.exp(-m_i/T)

        # relativistic m_i << T
        if limit == 'rel'
            n_i = g_i * T**3/np.pi**2

        return n_i

    def dYdx(self, Y, x):
        '''The differential equation for dY/dx, Y = nx/T^3, x = m/T'''
        sigv = 1e-30
        Yeq = number_density(self.Tlist) / (self.Tlist)**3
        lambd = (self.mass)**3 * sigv / H(self.mass)
        return (lambd / x**2) * (Y - Yeq)

    def xfreeze(self, sigv):
        Yinf = 1 #Y at very late times,
        return ((self.mass)**3 * sigv / H(self.mass)) * Yinf
