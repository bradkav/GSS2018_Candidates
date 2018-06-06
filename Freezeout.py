#The Hubble Constant

import numpy as np
import scipy
import matplotlib.pyplot as plt

M_pl = 2.435e18 #GeV/c^2

def Ht(T):
    return np.pi/3 * (g_star(T)/10)**(0.5) * T**2/M_pl

def g_star(T):
    return 10.75

def number_density(T, limit='rel'):

    # non-relativistic m_i >> T
    if limit == 'nonrel' :
        n_i = g_i * (m_i * T / (2*np.pi))**(3/2) * np.exp(-m_i/T)

    # relativistic m_i << T
    if limit == 'rel'
        n_i = g_i * T**3/np.pi**2

    return n_i
