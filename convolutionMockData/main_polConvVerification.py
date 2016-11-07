import sys; sys.path.append('../src/')
import numpy as np 
import irradiationSourceProfile as isp
import discreteFourierBesselTrafo as dFBT 

def flatTop(r, (R0,D0)=(0.5,0.05)):
        condList = [r<R0, r>=R0]
        funcList = [lambda r: 1, 
                    lambda r: np.exp(-(r-R0)**2/D0**2)]
        return np.piecewise(r,condList,funcList)

def deltaG(r,E0):
        return np.exp(-0.5*r*r/E0/E0)/(2*np.pi*E0*E0)

def eRMS(Fn,Fx):
    return np.sqrt(((Fn-Fx)**2).mean()/(Fx*Fx).mean())

def norm(r,f):
        return 1./(2*np.pi*np.trapz(r*f,r))


def main():
        T = 1.0 #float(sys.argv[1]) # truncation thresholf for ISP
        N = int(sys.argv[1]) # truncation parameter for FB expansion
        R0 = 0.3 #float(sys.argv[3]) # width of the top-hat part 
        A0 = float(sys.argv[2]) # decay width of Gaussian 
        E0 = float(sys.argv[3])

        # SEQUENCE OF SAMPLING POINTS FOR OBJECTIVE FUNCTIONS
        r = np.linspace(0,1.,1000)

        # SET CUSTOM PROBE FUNCTIONS 
        f = lambda x: flatTop(x, (R0,A0))
        g = lambda x: deltaG(x, E0)
        
        f0 = norm(r,f(r))
        g0 = norm(r,g(r))

        # FORWARD TRANSFORM FLAT-TOP
        rho,F0,T = dFBT.FiskJohnsonContinuousFuncFWD(r,f,T,N)
        # FORWARD TRANSFORM DELTA-APPROXIMATION 
        rho,G0,T = dFBT.FiskJohnsonContinuousFuncFWD(r,g,T,N)
        # BACKWARD TRANSFORM POINTWISE PRODUCT
        hr = dFBT.FiskJohnsonDiscreteFuncBCKWD(r,2.*np.pi*G0*F0,T)

        for ir in range(rho.size):
                print "FWD ", rho[ir],F0[ir],G0[ir]

        fr =f(r)
        # NORMALIZATION ACCOUNTS FOR DEVIATION OF ORIGINAL NORM IF 
        # WIDTH FALLS FAR BELOW LENGTHSCALE INDUCED BY MESH WIDTH
        for ir in range(r.size):
                print "BCKWD ", r[ir], fr[ir] , hr[ir]*fr[0]/hr[0]


def main3():
        T = 1.0 #float(sys.argv[1]) # truncation thresholf for ISP
        N = int(sys.argv[1]) # truncation parameter for FB expansion
        R0 = 0.3 #float(sys.argv[3]) # width of the top-hat part 
        A0 = float(sys.argv[2]) # decay width of Gaussian 
        E0 = float(sys.argv[3])

        # SEQUENCE OF SAMPLING POINTS FOR OBJECTIVE FUNCTIONS
        r = np.linspace(0,1.,1000)

        # SET CUSTOM PROBE FUNCTIONS 
        f = lambda x: flatTop(x, (R0,A0))
        g = lambda x: deltaG(x, E0)
        
        f0 = norm(r,f(r))
        g0 = norm(r,g(r))

        # FORWARD TRANSFORM FLAT-TOP
        rho,F0,T = dFBT.FiskJohnsonContinuousFuncFWD(r,f,T,N)
        # FORWARD TRANSFORM DELTA-APPROXIMATION 
        rho,G0,T = dFBT.FiskJohnsonContinuousFuncFWD(r,g,T,N)
        # BACKWARD TRANSFORM POINTWISE PRODUCT
        hr = dFBT.FiskJohnsonDiscreteFuncBCKWD(r,2.*np.pi*G0*F0,T)

        for ir in range(rho.size):
                print "FWD ", rho[ir],F0[ir],G0[ir]

        fr =f(r)
        # NORMALIZATION ACCOUNTS FOR DEVIATION OF ORIGINAL NORM IF 
        # WIDTH FALLS FAR BELOW LENGTHSCALE INDUCED BY MESH WIDTH
        for ir in range(r.size):
                print "BCKWD ", r[ir], fr[ir] , hr[ir]*fr[0]/hr[0]


def main2():
        T = 1.0 #float(sys.argv[1]) # truncation thresholf for ISP
        R0 = 0.4 #float(sys.argv[3]) # width of the top-hat part 
        A0 = 0.1 #float(sys.argv[4]) # decay width of Gaussian 
        E0 = float(sys.argv[1])

        # SEQUENCE OF SAMPLING POINTS FOR OBJECTIVE FUNCTIONS
        r = np.linspace(0,1.,1000)

        # SET CUSTOM PROBE FUNCTIONS 
        f = lambda x: flatTop(x, (R0,A0))
        g = lambda x: deltaG(x, E0)


        for N in range(4,50):
                # FORWARD TRANSFORM FLAT-TOP
                rho,F0,T = dFBT.FiskJohnsonContinuousFuncFWD(r,f,T,N)
                # FORWARD TRANSFORM DELTA-APPROXIMATION 
                rho,G0,T = dFBT.FiskJohnsonContinuousFuncFWD(r,g,T,N)
                # BACKWARD TRANSFORM POINTWISE PRODUCT
                hr = dFBT.FiskJohnsonDiscreteFuncBCKWD(r,2.*np.pi*G0*F0,T)
                print E0, N,  eRMS(f(r),hr)


        
main()
# EOF
