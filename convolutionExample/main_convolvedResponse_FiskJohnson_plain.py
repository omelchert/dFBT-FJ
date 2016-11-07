import sys; sys.path.append('../src/')
import dataIO as io
import convolveRadiallySymmetricFunctions as conv                                     
import irradiationSourceProfile as isp

def main():
        inFileName = sys.argv[1]
        outFileName = "testFJ"

        T = None # Truncation radius for beam profile
        N = 32 # number of samples to consider during FB transform 

        # FETCH IMPULSE RESPONSE FROM mco. FILE
        r,z,Arz = io.fetchRawData(inFileName)
        # ACCOUNT FOR VOID-WALKER
        Arz[-1,:] = 0.

        # SET CUSTOM IRRADIATION SOURCE PROFILE FROM OSG LIBRARY
        myIsp = lambda x: isp.flatTop(x, (0.05,0.05))
        # CONVOLVE IMPULSE RESPONSE USING CUSTOM IRRADIATION SOURCE PROFILE
        Wrz = conv.FiskJohnsonConvolvedResponse(r, z, Arz, myIsp, (T,N))

        # SAVE DATA IN GNUPLOT FORMAT
        io.writeGpl(r, z, Wrz, outFileName)
        
main()
# EOF
