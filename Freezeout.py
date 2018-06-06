#The Hubble Constant

import numpy as np
import scipy

M_pl = 2.435e18 #GeV/c^2



def Ht(T):
    return np.pi/3 * (g_star(T)/10)**(0.5) * T**2/M_pl

def g_star(T):
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
    if T <= 0.2 and T > 0.13: #QCD transition
        return 17.25
