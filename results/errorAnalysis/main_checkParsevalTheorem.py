import sys; sys.path.append('../../src/')
import numpy as np 
import FourierBesselPairs as FBP
import discreteFourierBesselTrafo as dFBT 

import scipy
import scipy.special as scs

def main():

        # SIMULATION PARAMETERS -----------------------------------------------
        T = float(sys.argv[1]) #18. 
        N = int(sys.argv[2]) #50
        
        rMin = 0.01
        rMax = 20.
        Nr = 1000

        r = np.linspace(rMin,rMax,Nr,retstep=False,endpoint=False)

        # OBJECTIVE FUNCTION --------------------------------------------------
        #f, F = FBP.sombrero()
        f, F = FBP.Gauss()
        for i in range(r.size):
                print "O", r[i],float(F(r[i]))


        # FISK JOHNSON FWD TRAFO FOR CONTINUOUS FUNCTION ----------------------
        rhoFJC,F0FJC,T = dFBT.FiskJohnsonContinuousFuncFWD(r,f,T,N)
        for i in range(rhoFJC.size):
                print "FJC ", rhoFJC[i], F0FJC[i] 

        def Y(m,k):
                j = scs.jn_zeros(0,N)
                jN = j[-1]
                J0 = lambda x: scs.j0(x)
                J1 = lambda x: scs.j1(x)
                return 2.0*J0(j[m]*j[k]/jN)/jN/J1(j[k])**2

        jm = scs.jn_zeros(0,N)

        fk  = F0FJC*2./(T*T*scs.j1(jm)**2)

        Fm = np.zeros(N)
        for im in range(N):
           for ik in range(N):
                 Fm[im] += Y(im,ik)*fk[ik]

        fks = fk /scs.j1(jm)
        Fms = Fm /scs.j1(jm)

        print fks
        print Fms
        print np.sum(fks[:-1]**2)
        print np.sum(Fms[:-1]**2)



main()
# EOF: main_checkParsevalTheorem.py
