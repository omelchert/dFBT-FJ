import unittest
import dFBT
import scipy.special as scs
import numpy as np


def FBPair():
    f  = lambda r: np.exp(-r*r/4/np.pi) 
    F0 = lambda r: np.exp(-r*r*np.pi)*2*np.pi     
    return f, F0


def eRMS(Fn,Fx):
    return np.sqrt(((Fn-Fx)**2).mean()/(Fx*Fx).mean())


class FourierBesselTransformTestCase(unittest.TestCase):

    def setUp(self):
        self.r = np.linspace(0.0, 20.0, 1000)
        self.f, self.Fx = FBPair()
        self.T, self.N = 18.0, 20

    def tearDown(self):
        del self.N, self.T, self.r, self.f, self.Fx

    def test_dFBT_fourierPair(self):
        rhoFJC,F0FJC,T = dFBT.fwdTrafo(self.r,self.f,self.T,self.N)
        F0FJCEx = dFBT.extrapolate(self.r,F0FJC,self.T)
        print eRMS(F0FJCEx,self.Fx(self.r))
        self.assertLessEqual(eRMS(F0FJCEx,self.Fx(self.r)), 1e-6)

    def test_dFBT_selfReciprocality(self): 
        rhoFJC,F0FJC,T = dFBT.fwdTrafo(self.r,self.f,self.T,self.N)
        fFJC = dFBT.bckwdTrafo(self.r,F0FJC,self.T)
        print eRMS(fFJC,self.f(self.r))
        self.assertLessEqual(eRMS(fFJC,self.f(self.r)), 1e-6)

    def test_dFBT_generalizedParsevalTheorem(self):
        j = scs.jn_zeros(0,self.N)
        Fm = np.zeros(self.N)
        Y = lambda m,k: 2.0*scs.j0(j[m]*j[k]/j[-1])/j[-1]/scs.j1(j[k])**2

        rhoFJC,F0FJC,T = dFBT.fwdTrafo(self.r,self.f,self.T,self.N)
        fk  = F0FJC*2./(self.T*self.T*scs.j1(j)**2)

        for im in range(self.N):
           for ik in range(self.N):
                 Fm[im] += Y(im,ik)*fk[ik]

        fkScaled = fk /scs.j1(j)
        FmScaled = Fm /scs.j1(j)

        print np.sum(FmScaled**2) -  np.sum(fkScaled**2)

        self.assertAlmostEqual(np.sum(FmScaled**2), np.sum(fkScaled**2), 6)
            

if __name__ == "__main__":
    unittest.main()
