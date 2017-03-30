import scipy.special as scs
import numpy as np 

def fwdTrafo(r,f,T=None,N=10):
        # SET TRUNCATION THRESHOLD FOR OBJECTIVE FUNCTION
        T = max(r) if T==None else min(T,max(r))

        jm = scs.jn_zeros(0,N)
        jp = jm[:-1,np.newaxis]
        jN = jm[-1]
        x = jp/jN

        # DISCRETE TRANSFORMATION EQUATION RELATING PARTICULAR VALUES F(j[m]/T) 
        # m=0...N-1 OF TRANSFORMED FUNCTION TO PARTICULAR VALUES f(x[p] T) 
        # p=0...N-2 OF OBJECTIVE FUNCTION. SEE EQ. (12) OF REF. [2].
        F0 = 2.0*(T/jN)**2 * np.sum(f(x*T)*scs.j0(x*jm)/scs.j1(jp)**2 ,axis=0)

        return jm/T, F0, T 


def bckwdTrafo(r,F0,T):
        f = np.zeros(r.size)
        jm = scs.jn_zeros(0,F0.size)
        
        # REVERSE TRANSFORM YIELDING ARBITRARY FUNCTION VALUES f(xT) FROM ITS 
        # FOURIER BESSEL TRANSFORM F(j[m]/T) m=0...N-1 AT SCALED BESSEL ZEROS 
        # j[m]/T. SEE EQ. (10) OF REF. [2].
        x = r/T
        f[x<1] = 2.0/T**2*np.sum(
            F0*scs.j0(jm*x[x<1,np.newaxis])/scs.j1(jm)**2, 
            axis=1)

        return f


def extrapolate(rho,F0,T):
        jFull = scs.jn_zeros(0,F0.size)
        jm = jFull[:-1]
        jN = jFull[-1] 

        # SET SEQUENCE OF SCALED TARGET COORDINATES
        r  = rho*T/jN 
        
        # EXTRAPOLATION FORMULA EQ. (9) OF REF. [2].
        F0Ex = 2.0*np.sum(
            F0[:-1]*scs.j0(r[:,np.newaxis]*jN)*jm/
            scs.j1(jm)/
            (jm**2-r[:,np.newaxis]**2*jN*jN),
            axis=1)

        return F0Ex

