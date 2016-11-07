import scipy
import scipy.special as scs
import numpy as np

def sombrero():
    a0 = 3.
    S  = lambda r: a0*scs.j1(a0*r)/r 
    H  = lambda r: r<a0 
    return S,H

def exp():
    a0  = 1.5
    S   = lambda r: np.exp(-a0*r)
    H   = lambda r: a0/(a0*a0+r*r)**1.5
    return S,H

def Gauss():
    a0 = 1./np.sqrt(np.pi)
    S  = lambda r: np.exp(-r*r/(4*np.pi)) 
    H  = lambda r: np.exp(-np.pi*r*r)*2*np.pi     
    return S,H

def Bessel():
    a0 = 20.
    S  = lambda r: scs.j0(a0*r)
    H  = lambda r: np.pi*2./r if r==a0 else 0
    return S,H

