#The Hubble Constant

import numpy as np
import scipy

M_pl = 2.435e18 #GeV/c^2

def Ht(T):
    return np.pi/3 * (g_star(T)/10)**(0.5) * T**2/M_pl

def g_star(T):
    return 10.75
