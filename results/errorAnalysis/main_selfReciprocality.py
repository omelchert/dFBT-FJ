import sys; sys.path.append('../../src/')
import scipy
import scipy.special as scs
import FourierBesselPairs as FBP
import numpy as np 
import discreteFourierBesselTrafo as dFBT 


def eRMS(Fn,Fx):
    """Compute root mean square error for input arrays.

    Implements root mean square error (eRMS) according to Eq (25) in Ref

    Algorithms to Numerically Evaluate the Hankel Transform
    Cree, M. J. and Bones P. J.
    Computers Math. Applic. 26 (1993) 1-12 

    Args:
        Fn: numpy array containing values of Hankel transform of test function.
        Fx: numpy array with same length as Fn containing values of exact
            Hankel transform of test function.

    """
    return np.sqrt(((Fn-Fx)**2).mean()/(Fx*Fx).mean())


def main():
    # SIMULATION PARAMETERS -----------------------------------------------
    #T = 10. 
    N = 5
    rMin = 0.01
    rMax = 10.
    Nr = 1000

    for T in np.linspace(1.,rMax, 20.):
        # SET OBJECTIVE FUNCTION -----------------------------------------------
        r = np.linspace(rMin,rMax,Nr,retstep=False,endpoint=False)
        f, F = FBP.sombrero()

        print T, N,

        # FISK JOHNSON dFBT FOR CONTINUOUS FUNCTION ---------------------------
        rhoFJC,F0FJC,T = dFBT.FiskJohnsonContinuousFuncFWD(r,f,T,N)
        fFJC = dFBT.FiskJohnsonDiscreteFuncBCKWD(r,F0FJC,T)
        print eRMS(f(r),fFJC),

        # FISK JOHNSON dFBT FOR DISCRETE FUNCTION -----------------------------
        rhoFJD,F0FJD,T = dFBT.FiskJohnsonDiscreteFuncFWD(r,f(r),T,N)
        fFJD = dFBT.FiskJohnsonDiscreteFuncBCKWD(r,F0FJD,T)
        print eRMS(f(r),fFJD),

        # CREE BONES dFBT FOR CONTINUOUS FUNCTION -----------------------------
        rhoCB,F0CB  = dFBT.CreeBonesDiscreteFunc(r,f(r))
        rCB,fCB =  dFBT.CreeBonesDiscreteFunc(rhoCB,F0CB)
        print eRMS(f(r),fCB)


main()
# EOF: main_selfReciprocality.py
