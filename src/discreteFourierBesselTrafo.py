""" FILE: discreteFourierBesselTrafo.py

Module implementing functions to compute Fourier-Bessel function (zero order 
Hankel transform) in terms of a rapidly convergent Fourier-Bessel expansion
of the function of interest. 

The algorithmic procedure follows the method detailed in Ref. [1]. The theory
and operational rules for the general nth order Hankel transform are 
thoroughly discussed in Ref. [2].

Refs:
    [1] An Improved Method for Computing a Discrete Hankel Transform
        H. Fisk Johnson
        Comp. Phys. Commun. 43 (1987) 181-202

    [2] Theory and operational rules for the discrete Hankel transform
        N. Baddour, U. Chouinard
        J. Opt. Soc. Am. A 32 (2015) 611
        
"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import scipy
import scipy.special as scs
import numpy as np 


def FiskJohnsonContinuousFuncFWD(r,f,T=None,N=10):
        """Compute Fourier-Bessel transformation of continuous function via 
            Fisk Johnson procedure.

        Compute Fourier-Bessel transform (i.e. 0th order Hankel transform) 
        using a rapidly convergent summation of a Fourier-Bessel expansion 
        following the metod introduced in Ref. [1] and further detailed in 
        Ref. [2].


        Args:
            r (numpy array, ndim=1): equispaced 1D grid of coordinates.
            f (numpy array, ndim=1): objective function.
            T (float): truncation threshold for objective function. 
            N (int): maximal order of Bessel function zeros to consider.

        Returns:
            rho (numpy array, ndim=1): sampling points for transformed grid. 
            F0 (numpy array, ndim=1): Fourier Bessel transform of objective 
                function.
            T (float): truncation radius for objective function used during
                forward Fourier-Bessel transformation.

        Notes:
            - Fisk Johnson procedure for continuous function f(r)
            - Note that the notation follows Ref. [1]. The computational 
                efficiency is O(N^2).
            - Yields significant overall reduction of computation time if 
                follow up back transformation is needed.
            - Implements Eq. (12) of Ref. [1].

        Refs:
            [1] An Improved Method for Computing a Discrete Hankel Transform
                H. Fisk Johnson
                Comp. Phys. Commun. 43 (1987) 181-202
                
            [2] Theory and operational rules for the discrete Hankel transform
                N. Baddour, U. Chouinard
                J. Opt. Soc. Am. A 32 (2015) 611
        """
        # SET TRUNCATION THRESHOLD FOR OBJECTIVE FUNCTION
        T = max(r) if T==None else min(T,max(r))

        # COMPUTE FIRST N ZEROS OF 0TH ORDER BESSEL FUNCTION IN ASCENDING ORDER
        jm = scs.jn_zeros(0,N)

        # CONVENIENT ABBREVIATIONS
        jp = jm[:-1,np.newaxis]
        jN = jm[-1]
        x = jp/jN

        # DISCRETE TRANSFORMATION EQUATION RELATING PARTICULAR VALUES F(j[m]/T) 
        # m=0...N-1 OF TRANSFORMED FUNCTION TO PARTICULAR VALUES f(x[p] T) 
        # p=0...N-2 OF OBJECTIVE FUNCTION. SEE EQ. [12] OF REF. [1].
        F0 = 2.0*(T/jN)**2 * np.sum(f(x*T)*scs.j0(x*jm)/scs.j1(jp)**2 ,axis=0)

        return jm/T, F0, T 


def FiskJohnsonDiscreteFuncFWD(r,fr,T=None,N=10):
        """Compute Fourier-Bessel transformation of discrete function via 
            Fisk Johnson procedure.

        Compute Fourier-Bessel transform (i.e. 0th order Hankel transform) 
        using a rapidly convergent summation of a Fourier-Bessel expansion 
        following the metod introduced in Ref. [1] and further detailed in 
        Ref. [2].


        Args:
            r (numpy array, ndim=1): equispaced 1D grid of coordinates.
            fr (numpy array, ndim=1): objective function at discrete coordinate
                values of r.
            T (float): truncation threshold for objective function. 
            N (int): maximal order of Bessel function zeros to consider.

        Returns:
            rho (numpy array, ndim=1): sampling points for transformed grid. 
            F0 (numpy array, ndim=1): Fourier Bessel transform of objective 
                function.
            T (float): truncation radius for objective function used during
                forward Fourier-Bessel transformation.

        Notes:
            - Fisk Johnson procedure for function fr known at discrete 
                coordinate values of r, only.
            - Note that the notation follows Ref. [1]. The computational 
                efficiency is O(N^2).
            - Considerable reduction of computation time if follow up back
                transformation is needed.
            - Yields significant overall reduction of computation time if 
                follow up back transformation is needed.
            - Implements Eqs. (7), (8) of Ref. [1].

        Refs:
            [1] An Improved Method for Computing a Discrete Hankel Transform
                H. Fisk Johnson
                Comp. Phys. Commun. 43 (1987) 181-202
                
            [2] Theory and operational rules for the discrete Hankel transform
                N. Baddour, U. Chouinard
                J. Opt. Soc. Am. A 32 (2015) 611
        """
        # SET TRUNCATION THRESHOLD FOR OBJECTIVE FUNCTION AND TRUNCATE
        # COORDINATE AXIS AND OBJECTIVE FUNCTION ARRAY CORRESPONDINGLY
        T = max(r) if T==None else min(T,max(r))
        rTrunc = r[r<T]
        fTrunc = fr[r<T]
        
        # COMPUTE FIRST N ZEROS OF 0TH ORDER BESSEL FUNCTION IN ASCENDING ORDER
        jm = scs.jn_zeros(0,N)

        # COMPUTATION OF FIRST m FOURIER BESSEL COEFFICIENTS FOLLOWING EQ. (7)
        # OF REF. [1].
        x = rTrunc/T
        C = np.trapz( x*fTrunc*scs.j0(jm[:,np.newaxis]*x),dx=x[1]-x[0],axis=1) 
        
        return jm/T, T*T*C, T


def FiskJohnsonDiscreteFuncBCKWD(r,F0,T):
        """Compute reverse Fourier-Bessel transformation via Fisk Johnson 
            procedure.

        Compute reverse Fourier-Bessel transform (i.e. 0th order reverse Hankel 
        transform) using a rapidly convergent summation of a Fourier-Bessel 
        expansion following the metod introduced in Ref. [1] and further
        detailed in Ref. [2].


        Args:
            r (numpy array, ndim=1): equispaced 1D grid of target coordinates.
            F0 (numpy array, ndim=1): Fourier-Bessel transformed function 
                at discrete coordinates given by its scaled bessel zeros. 
            T (float): truncation threshold for objective function. 

        Returns:
            f (numpy array, ndim=1): reverse Fourier-Bessel transform of input
                function.

        Notes:
            - Fisk Johnson procedure for reverse Fourier-Bessel transformation. 
            - Implements Eq. (10) of Ref. [1].
            - above truncation threshold it holds that f(r>T) = 0.
            - on input F0 = F0[jm/T] for m = 0...N-1 where jm are the first
                N zeros of the 0th order Bessel function in ascending order.

        Refs:
            [1] An Improved Method for Computing a Discrete Hankel Transform
                H. Fisk Johnson
                Comp. Phys. Commun. 43 (1987) 181-202
                
            [2] Theory and operational rules for the discrete Hankel transform
                N. Baddour, U. Chouinard
                J. Opt. Soc. Am. A 32 (2015) 611
        """
        # INITIALIZE EMPTY ARRAY FOR REVESE TRANSFORM
        f = np.zeros(r.size)

        # COMPUTE FIRST N ZEROS OF 0TH ORDER BESSEL FUNCTION IN ASCENDING ORDER
        jm = scs.jn_zeros(0,F0.size)
        
        # REVERSE TRANSFORM YIELDING ARBITRARY FUNCTION VALUES f(xT) FROM ITS 
        # FOURIER BESSEL TRANSFORM F(j[m]/T) m=0...N-1 AT SCALED BESSEL ZEROS 
        # j[m]/T. SEE EQ. (10) OF REF. [1].
        x = r/T
        f[x<1] = 2.0/T**2*np.sum(
            F0*scs.j0(jm*x[x<1,np.newaxis])/scs.j1(jm)**2, 
            axis=1)

        return f


def FiskJohnsonDiscreteFuncExtrapolate(rho,F0,T):
        """Extrapolate samples of transform to sequence of target coordinates.

        Extrapolate sequence of (few) Fourier-Bessel transform samples to 
        sequence of target coordinates following Eq. (9) of Ref. [1].

        Args:
            rho (numpy array, ndim=1): equispaced 1D grid of target coordinates.
            F0 (numpy array, ndim=1): Fourier-Bessel transformed function 
                at discrete coordinates given by its scaled bessel zeros. 
            T (float): truncation threshold for objective function. 

        Returns:
            F0Ex (numpy array, ndim=1): Fourier-Bessel transform extrapolated
                to sequence of target coordinates.

        Notes:
            - Implements Eq. (9) of Ref. [1].
            - Extrapolation for 0 < r < \infty.
            - on input F0 = F0[jm/T] for m = 0...N-1 where jm are the first
                N zeros of the 0th order Bessel function in ascending order.

        Refs:
            [1] An Improved Method for Computing a Discrete Hankel Transform
                H. Fisk Johnson
                Comp. Phys. Commun. 43 (1987) 181-202
        """
        # COMPUTE FIRST N ZEROS OF 0TH ORDER BESSEL FUNCTION IN ASCENDING ORDER
        jFull = scs.jn_zeros(0,F0.size)
        jm = jFull[:-1]
        jN = jFull[-1] 

        # SET SEQUENCE OF SCALED TARGET COORDINATES
        r  = rho*T/jN 
        
        # EXTRAPOLATION FORMULA EQ. (9) OF REF. [1].
        F0Ex = 2.0*np.sum(
            F0[:-1]*scs.j0(r[:,np.newaxis]*jN)*jm/
            scs.j1(jm)/
            (jm**2-r[:,np.newaxis]**2*jN*jN),
            axis=1)

        return F0Ex


def CreeBonesDiscreteFunc(r,fr):
        """Compute zero order Hankel transform.

        Compute zeroth order hankel transform using trapezoidal quadrature 
        rule for equally spaced samples as explained in section 4 of Ref. [2].

        Args:
            r (numpy array, ndim=1): equispaced 1D grid for radial coordinate.
            fr (numpy array, ndim=1): objective function.

        Returns:
            rho (numpy array, ndim=1): equi-spaced complementary grid.
            F0 (numpy array, ndim=1): zeroth order Hankel transform of 
                objective function.

        Notes:
            Note that the notation follows Refs. [1], [3]. The computational 
            efficiency is O(N^2).

        Refs:
            [1] Operational and convolution properties of two-dimensional 
                Fourier transforms in polar coordinates 
                Baddour, N.
                J. Opt. Soc. Am. A 26 (2009) 1767-1777 
          
            [2] Algorithms to Numerically Evaluate the Hankel Transform
                Cree, M. J. and Bones, P. J.
                Computers Math. Applic., 26 (1993) 1-12
          
            [3] Fast computation of zero order Hankel transform
                Gopalan, K. and Chen, C. S. 
                Journal of the Franklin Institute, 316 (1983) 317-326
        """
        dr    = r[1]-r[0]
        rho   = np.linspace(0,1./(2.*dr),r.size,endpoint=False)
        rfrJ0 = r*fr*scs.j0(r*rho[:,np.newaxis])
        return rho,np.trapz(rfrJ0,x=r,axis=1) 

# EOF: discreteFourierBesselTrafo.py
