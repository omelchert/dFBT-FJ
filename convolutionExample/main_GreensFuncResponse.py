import sys; sys.path.append('../src/')
import dataIO as io
import convolveRadiallySymmetricFunctions as conv                                     
import irradiationSourceProfile as isp

def main():
        inFileName = sys.argv[1] # input .mco file
        outFileName = sys.argv[2] # output .prof file

        # FETCH IMPULSE RESPONSE FROM mco. FILE
        r,z,Arz = io.fetchRawData(inFileName)
        # ACCOUNT FOR VOID-WALKER
        Arz[-1,:] = 0.

        # SAVE DATA IN GNUPLOT FORMAT
        io.writeGpl(r, z, Arz, outFileName)
        
main()
# EOF
