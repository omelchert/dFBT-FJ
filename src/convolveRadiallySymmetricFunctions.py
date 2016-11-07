""" FILE: convolveRadiallySymmetricFunctions.py

Module implementing 2D convolution of two radially symmetric functions 
based on two differnt approaches to compute Fourier-Bessel functions 
as explained in Refs. [2] and [3].

Refs:
    [1] Operational and convolution properties of two-dimensional Fourier 
        transforms in polar coordinates 
        Baddour, N.
        J. Opt. Soc. Am. A 26 (2009) 1767-1777 

    [2] Algorithms to Numerically Evaluate the Hankel Transform
        Cree, M. J. and Bones, P. J.
        Computers Math. Applic., 26 (1993) 1-12

    [3] An Improved Method for Computing a Discrete Hankel Transform
        H. Fisk Johnson
        Comp. Phys. Commun. 43 (1987) 181-202
"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import scipy
import scipy.special as scs
import numpy as np 
import irradiationSourceProfile as isp
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


def CreeBonesConvolvedResponse(r,z,grz,f,P=1.0):
        """Compute response to spatially extended irradiation source profile.

        Computes full 2D convolution of pencil beam response with extended
        beam profile at each z-coordinate of region of interest (ROI).

        Args:
            r (numpy array, ndim=1): equispaced 1D grid for radial coordinate.
            z (numpy array, ndim=1): equispaced 1D grid for depth coordinate.
            grz (numpy array, ndim=2): response to infinitely narrow beam.
            isp (lambda function): radial beam profile.
            P (float): beam intensity (default: 1.0).

        Returns:
            Wrz (numpy array, ndim=2): response to spatially extended 
                irradiation source profile.

        See Also:
            convolve: 2D convolution of two radially symmetric functions.

        Notes: 
            Due to the spatial invariance of the region of interest, the 
            material response to an extended photon beam with given radial 
            profile is obtained as convolution of two radially symmetric 
            functions for given z-coordinate, see Eq. (84) of Ref. [1]. For the
            integration, a numpy native trapezoidal rule is used.

            The 2D convolution integral given by Eq. (84) of see Ref. [1] 
            relies on the zeroth order coefficients G0(rho)=2\pi H0(g(r)) and 
            F0(rho)=2\pi H0(f(r)) of the Fourier Transforms of f and g, NOT 
            their plain Hankel transforms. For radially symmetric functions 
            only this zeroth order Fourier coefficient is nonzero. Thus, the 
            integrahd formulated in terms of their Hankel transforms introduces 
            an additional factor of (2\pi)^2 in Eq. (84).

        Refs:
            [1] Operational and convolution properties of two-dimensional 
                Fourier transforms in polar coordinates 
                Baddour, N.
                J. Opt. Soc. Am. A 26 (2009) 1767-1777 
        """
        Wrz  = np.zeros(grz.shape)

        fr = f(r)
        # Normalize beam profile using function from isp module
        f0 = isp.beamProfNormalization(r,fr,P)
        # Hankel transform of beam profile will be used repeatedly so 
        # precompute it here for time-efficiency 
        rho,F0 = dFBT.CreeBonesDiscreteFunc(r,fr)

        for iz in range(z.size):
            #print iz
            # Hankel transform of response to pencil beam
            rho,G0 = dFBT.CreeBonesDiscreteFunc(r,grz[:,iz])
            hr = 2.*np.pi*dFBT.CreeBonesDiscreteFunc(rho,G0*F0)[1]
            # Gibbs-phenomenon workaround: remove all negative valued 
            # contributions to volumetric energy density at given depth
            hr = np.abs(hr)
            # Normalize contribution to volumetric energy density
            Wrz[:,iz] = f0*hr
        return Wrz


def FiskJohnsonConvolvedResponse(r,z,grz,f,(T,N),P=1.0,):
        """Compute response to spatially extended irradiation source profile.

        Computes full 2D convolution of pencil beam response with extended
        beam profile at each z-coordinate of region of interest (ROI) using
        a polar convolution algorithm based on the procedure reported in 
        Refs.~[1], [2].

        Args:
            r (numpy array, ndim=1): equispaced 1D grid for radial coordinate.
            z (numpy array, ndim=1): equispaced 1D grid for depth coordinate.
            grz (numpy array, ndim=2): response to infinitely narrow beam.
            isp (lambda function): radial beam profile.
            T (float): truncation radius for beam profile
            N (int): number of transformed samples to consider.
            P (float): beam intensity (default: 1.0).

        Returns:
            Wrz (numpy array, ndim=2): response to spatially extended 
                irradiation source profile.
            recErr (float): reconstruction error.

        Notes: 
            Due to the spatial invariance of the region of interest, the 
            material response to an extended photon beam with given radial 
            profile is obtained as convolution of two radially symmetric 
            functions for given z-coordinate, see Eq. (84) of Ref. [3]. For the
            integration, a numpy native trapezoidal rule is used.

            The 2D convolution integral given by Eq. (84) of see Ref. [3] 
            relies on the zeroth order coefficients G0(rho)=2\pi H0(g(r)) and 
            F0(rho)=2\pi H0(f(r)) of the Fourier Transforms of f and g, NOT 
            their plain Hankel transforms. For radially symmetric functions 
            only this zeroth order Fourier coefficient is nonzero. Thus, the 
            integrand formulated in terms of their Hankel transforms introduces 
            an additional factor of (2\pi)^2 in Eq. (84).

        Refs:
            [1] An Improved Method for Computing a Discrete Hankel Transform
                H. Fisk Johnson
                Comp. Phys. Commun. 43 (1987) 181-202
            [2] Theory and operational rules for the discrete Hankel transform
                N. Baddour, U. Chouinard
                J. Opt. Soc. Am. A 32 (2015) 611
            [3] Operational and convolution properties of two-dimensional 
                Fourier transforms in polar coordinates 
                Baddour, N.
                J. Opt. Soc. Am. A 26 (2009) 1767-1777 
        """
        Wrz  = np.zeros(grz.shape)

        # Normalize beam profile using function from isp module
        f0 = isp.beamProfNormalization(r,f(r),P)
        # Hankel transform of beam profile will be used repeatedly so 
        # precompute it here for time-efficiency 
        rho,F0,T = dFBT.FiskJohnsonContinuousFuncFWD(r,f,T,N)
        fRec = dFBT.FiskJohnsonDiscreteFuncBCKWD(r,F0,T)
        recErr = eRMS(f(r),fRec)

        for iz in range(z.size):
            # print iz
            # Hankel transform of response to pencil beam
            rho,G0,T = dFBT.FiskJohnsonDiscreteFuncFWD(r,grz[:,iz],T,N)
            hr = dFBT.FiskJohnsonDiscreteFuncBCKWD(r,2.*np.pi*G0*F0,T)
            # Gibbs-phenomenon workaround: remove all negative valued 
            # contributions to volumetric energy density at given depth
            hr = np.abs(hr)
            # Normalize contribution to volumetric energy density
            Wrz[:,iz] = f0*hr

        return Wrz, recErr

# EOF: convolveRadiallySymmetricFunctions.py
