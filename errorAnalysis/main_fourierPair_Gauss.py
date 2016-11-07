import sys; sys.path.append('../src/')
import numpy as np 
import FourierBesselPairs as FBP
import discreteFourierBesselTrafo as dFBT 


def main():

        # SIMULATION PARAMETERS -----------------------------------------------
        T = float(sys.argv[1]) #18. 
        N = int(sys.argv[2]) #50
        
        rMin = 0.0
        rMax = 10.
        Nr = 1000

        r = np.linspace(rMin,rMax,Nr,retstep=False,endpoint=False)

        # OBJECTIVE FUNCTION --------------------------------------------------
        f, F = FBP.Gauss()
        for i in range(r.size):
                print "O", r[i],float(F(r[i]))

        # QUADRATURE USING TRAPEZOIDAL RULE ----------------------------------- 
        rhoCB,F0CB  = dFBT.CreeBonesDiscreteFunc(r,f(r))
        for i in range(rhoCB.size):
                print "CB ", rhoCB[i], F0CB[i] 

        # FISK JOHNSON FWD TRAFO FOR CONTINUOUS FUNCTION ----------------------
        rhoFJC,F0FJC,T = dFBT.FiskJohnsonContinuousFuncFWD(r,f,T,N)
        for i in range(rhoFJC.size):
                print "FJC ", rhoFJC[i], F0FJC[i] 
        
        # FISK JOHNSON FWD TRAFO FOR DISCRETE FUNCTION ------------------------
        rhoFJD,F0FJD,T = dFBT.FiskJohnsonDiscreteFuncFWD(r,f(r),T,N)
        for i in range(rhoFJD.size):
                print "FJD ", rhoFJD[i], F0FJD[i] 

        # FISK JOHNSON EXTRAPOLATION TO DESIRED SAMPLE POINTS -----------------
        F0FJCEx = dFBT.FiskJohnsonDiscreteFuncExtrapolate(rhoCB,F0FJD,T)
        for i in range(rhoCB.size):
                print "FJCEx ", rhoCB[i], F0FJCEx[i] 


main()
# EOF: main_fourierTransform_sombrero.py
