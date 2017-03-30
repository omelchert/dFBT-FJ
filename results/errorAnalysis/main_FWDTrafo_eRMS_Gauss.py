import sys; sys.path.append('../../src/')
import numpy as np 
import FourierBesselPairs as FBP
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
    # SIMULATION PARAMETERS ---------------------------------------------------
    #T = 10. 
    N = int(sys.argv[1])
    rMin = 0.01
    rMax = 10.
    Nr = 1000
    # SET OBJECTIVE FUNCTION --------------------------------------------------
    r = np.linspace(rMin,rMax,Nr,retstep=False,endpoint=False)
    f, F = FBP.Gauss()

    for T in np.linspace(20*rMin,rMax, 100.):

        # CREE BONES dFBT FOR CONTINUOUS FUNCTION -----------------------------
        rhoCB,F0CB  = dFBT.CreeBonesDiscreteFunc(r[r<T],f(r[r<T]))

        # FISK JOHNSON dFBT FOR CONTINUOUS FUNCTION ---------------------------
        rhoFJC,F0FJC,T = dFBT.FiskJohnsonContinuousFuncFWD(r,f,T,N)

        # FISK JOHNSON EXTRAPOLATION TO DESIRED SAMPLE POINTS -----------------
        F0FJCEx = dFBT.FiskJohnsonDiscreteFuncExtrapolate(rhoCB,F0FJC,T)

        print T,N, eRMS(F(rhoCB),F0CB), eRMS(F(rhoCB),F0FJCEx)


main()
# EOF: main_FWDTrafo_eRMS_Gauss.py
