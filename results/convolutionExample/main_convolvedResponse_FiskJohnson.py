import sys; sys.path.append('../../src/')
import numpy as np
import dataIO as io
import convolveRadiallySymmetricFunctions as conv                                     

def donut(r, (R0,R1,D0,D1)=(0.5,0.7,0.1,0.1)):
        condList = [r<R0, (r>=R0) & (r<=R1), r>R1]
        funcList = [lambda r: np.exp(-(r-R0)**2/D0**2), 
                    lambda r: 1., 
                    lambda r: np.exp(-(r-R1)**2/D1**2)]
        return np.piecewise(r,condList,funcList)

def main():
        inFileName = sys.argv[1] # input .mco file
        outFileName = sys.argv[8] # output .prof file
        T = float(sys.argv[2]) # truncation thresholf for ISP
        N = int(sys.argv[3]) # truncation parameter for FB expansion
        R0 = float(sys.argv[4]) # width of the top-hat part 
        A0 = float(sys.argv[5]) # decay width of Gaussian 
        R1 = float(sys.argv[6]) # width of the top-hat part 
        A1 = float(sys.argv[7]) # decay width of Gaussian 

        # FETCH IMPULSE RESPONSE FROM mco. FILE
        r,z,Arz = io.fetchRawData(inFileName)

        # CORRECT AXES TO REFER TO BIN-CENTERS: REF [1], EQS. (4.11), (4.19)
        r = r + 0.5*(r[1]-r[0])
        z = z + 0.5*(z[1]-z[0])
        # CORRECT FOR ACCUMULUATED VOID WALKERS: REF [2], EQ. (22) 
        Arz[-1:,:]=0

        # SET CUSTOM IRRADIATION SOURCE PROFILE FROM OSG LIBRARY
        myIsp = lambda x: donut(x, (R0,R1,A0,A1))
        # CONVOLVE IMPULSE RESPONSE USING CUSTOM IRRADIATION SOURCE PROFILE
        Wrz,recErr = conv.FiskJohnsonConvolvedResponse(r, z, Arz, myIsp, (T,N))

        # SAVE DATA IN GNUPLOT FORMAT
        io.writeGpl(r, z, Wrz, outFileName)
        f = open(outFileName,'a')
        f.write("# ISPRecErr: %lf\n"%(recErr))
        f.write("# maxW: %lf\n"%(Wrz.max()))
        f.close()
        
main()
# EOF
