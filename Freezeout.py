#The Hubble Constant

import numpy as np
import scipy

M_pl = 2.435e18 #GeV/c^2


class Freezeout(object):
    def __init__(mass = 1):


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
        for T in Tlist:
            if T > mass
            l +=1
        for n in range(0, l):
            geff = g_star(T[n]) + dof
        for m in range(l, len(T)):
            geff = g_star(T[n])
        return geff


    def dYdx(self, Yn, x):
        '''The '''
        lambda = self.mass**3 sigv / H(m)
