""" FILE: test_dFBT.py

Unittest module implementing unit tests for modul discreteFourierBesselTrafo.py
which implements the discrete Fourier Bessel transform of order zero as 
detailed in Ref. [1].

The general operational rules for nth order transform are thoroughly discussed
in Ref. [2]. Here, they are geared towards the case n=0, allowing to perform
test on the Fourier Bessel transform or zeroth order, also referred to as a
Hankel transform.

Refs:
    [1] An Improved Method for Computing a Discrete Hankel Transform
        H. Fisk Johnson
        Comp. Phys. Commun. 43 (1987) 181-202

    [2] Theory and operational rules for the discrete Hankel transform
        N. Baddour, U. Chouinard
        J. Opt. Soc. Am. A 32 (2015) 611

"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2017, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import sys; sys.path.append("../src/")
import unittest
import scipy.special as scs
import numpy as np
import discreteFourierBesselTrafo as dFBT


def Gauss():
    a0 = 1./np.sqrt(np.pi)
    S  = lambda r: np.exp(-r*r/a0/a0/(2*np.pi)**2) 
    H  = lambda r: np.exp(-r*r/a0/a0)*2*np.pi     
    return S,H


def eRMS(Fn,Fx):
    """Compute root mean square error for input arrays.

    Implements root mean square error (eRMS) according to Eq (25) in Ref. [1]

    Args:
        Fn (numpy array, ndim=1): numeric transform of test function
        Fx (numpy array, ndim=1): exact transform of test function

    Refs:
        [1] Algorithms to Numerically Evaluate the Hankel Transform
            Cree, M. J. and Bones P. J.
            Computers Math. Applic. 26 (1993) 1-12 

    """
    return np.sqrt(((Fn-Fx)**2).mean()/(Fx*Fx).mean())


class FourierBesselTransformTestCase(unittest.TestCase):
    """Unit test for discrete Fourier Bessel transform
    
    Implements three unit tests for discrete fourier bessel transform
    """

    def setUp(self):
        """Sets up unit test prerequisits
        
        Attributes:
            T (float): truncation threshold for objective function
            N (int): maximal order of Bessel function zeros to consider
            r: numpy array containing mesh points for test function
            f: test function
            Fx: exact transform of test function
        """
        rMin, rMax, Nr = 0.0, 20.0, 1000
        self.T, self.N = 18.0, 20
        self.r = np.linspace(rMin,rMax,Nr,retstep=False,endpoint=False)
        self.f, self.Fx = Gauss()

    def tearDown(self):
        """Deletes all attributes set by method setUp()"""
        del self.N
        del self.T
        del self.r
        del self.f
        del self.Fx

    def test_dFBT_fourierPair(self):
        """Perform unit test on Fourier Bessel pair"""
        # FISK JOHNSON FWD TRAFO - SEE EQ. (12), REF. [1]
        rhoFJC,F0FJC,T = dFBT.FiskJohnsonContinuousFuncFWD(self.r,self.f,self.T,self.N)
        # EXTRAPOLATION TO DESIRED SAMPLE POINTS - SEE EQ. (9), REF. [1] 
        F0FJCEx = dFBT.FiskJohnsonDiscreteFuncExtrapolate(self.r,F0FJC,self.T)
        print eRMS(F0FJCEx,self.Fx(self.r))
        self.assertLessEqual(eRMS(F0FJCEx,self.Fx(self.r)), 1e-6)

    def test_dFBT_selfReciprocality(self): 
        """Perform unit test checking self reciprocality"""
        # FISK JOHNSON FWD TRAFO - SEE EQ. (12), REF. [1]
        rhoFJC,F0FJC,T = dFBT.FiskJohnsonContinuousFuncFWD(self.r,self.f,self.T,self.N)
        # FISK JOHNSON BCKWD TRAFO - SEE EQ. (10), REF. [1]
        fFJC = dFBT.FiskJohnsonDiscreteFuncBCKWD(self.r,F0FJC,self.T)
        print eRMS(fFJC,self.f(self.r))
        self.assertLessEqual(eRMS(fFJC,self.f(self.r)), 1e-6)

    def test_dFBT_generalizedParsevalTheorem(self):
        """Perform unit test based on generalized Parseval theorem"""
        j = scs.jn_zeros(0,self.N)
        Fm = np.zeros(self.N)

        # DFT KERNEL - SEE EQ. (19), REF. [2] 
        Y = lambda m,k: 2.0*scs.j0(j[m]*j[k]/j[-1])/j[-1]/scs.j1(j[k])**2
        # FISK JOHNSON FWD TRAFO - SEE EQ. (12), REF. [1]
        rhoFJC,F0FJC,T = dFBT.FiskJohnsonContinuousFuncFWD(self.r,self.f,self.T,self.N)
        # FOURIER-BESSEL COEFFICIENTS - SEE EQ. (8), REF. [2] 
        fk  = F0FJC*2./(self.T*self.T*scs.j1(j)**2)

        # FORWARD TRANSFORM VIA DFT KERNEL Y - SEE EQ. (33), REF. [2] 
        for im in range(self.N):
           for ik in range(self.N):
                 Fm[im] += Y(im,ik)*fk[ik]

        # SCALED COEFFICIENTS - SEE EQ. (46), REF. [2] 
        fkScaled = fk /scs.j1(j)
        FmScaled = Fm /scs.j1(j)

        self.assertAlmostEqual(np.sum(FmScaled**2), np.sum(fkScaled**2), 6)
            

if __name__ == "__main__":
        unittest.main()

# EOF: test_dFBT.py
