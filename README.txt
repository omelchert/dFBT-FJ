1. INTRODUCTION

Efficient polar convolution based on the discrete Fourier-Bessel transform for
application in computational biophotonics

We illustrate efficient algorithms for the accurate forward and reverse
evaluation of the discrete Fourier-Bessel transform (dFBT) as numerical tools
to assist in the $2$D polar convolution of two radially symmetric functions,
relevant, e.g., to applications in computational biophotonics.
  In our survey of the numerical procedure we account for the circumstance that
the objective function might result from a more complex measurement process and
is, in the worst case, known on a finite sequence of coordinate values, only.
  We contrast the performance of the resulting algorithms with a procedure
based on a straight forward numerical quadrature of the underlying integral
transform and asses its efficienty for two benchmark Fourier-Bessel pairs. 
  An application to the problem of finite-size beam-shape convolution in polar
coordinates, relevant in the context of tissue optics and optoacoustics, is
used to illustrate the versatility and computational efficiency of the
numerical procedure.


2. DEPENDENCIES

The software was developed and tested under OS X Yosemite (Version: 10.10.3)
but should be able to run on any system that features the necessary version 
of Python and has all dependency modules available.

Python -- Version 2.7.6 or higher
numpy  -- Version 1.8.0rc1 or higher
scipy  -- Version 0.13.0b1 or higher 


3. CONTENT

README.txt              -- this readme document
LICENSE.txt             -- BSD 3-Clause License file

dFBT/
     convolutionExample/
         convolvedResponse_semiInf
         getData_FJ.sh
         getData_GreensFuncResponse.sh
         getData_profiling.sh
         GP/
             beamShapeConvolution.gpi
             FIGS/
                 beamShapeConvolution.eps       -- FIG 4
                 timeEfficiency.eps             -- FIG 6
                 WrzSlices.eps                  -- FIG 5
             magma.pal
             timeEfficiency.gpi
             WrzSlices.gpi
         impulseResponse_semiInf
         main_convolvedResponse_CreeBones.py
         main_convolvedResponse_FiskJohnson.py
         main_convolvedResponse_FiskJohnson_plain.py
         main_GreensFuncResponse.py
         main_pstat.py
         main_WrzSlice.py
         NoteOnTiming.txt
         profilingResults/
         WrzSlices/
     convolutionMockData/
         data/
         getData.sh
         GP/
             examples_polConv.gpi
             FIGS/
                 examples_polConv.eps           -- FIG 3
         main_polConvVerification.py
         main_polConvVerification_varyN.py
     errorAnalysis/
         data_Gauss/
         data_sombrero/
         FourierBesselPairs.py
         FourierBesselPairs.pyc
         getData_Gauss.sh
         getData_sombrero.sh
         GP/
             FIGS/
                 forwardTrafo.eps               -- FIG 1
                 reverseTrafo.eps               -- FIG 2
             forwardTrafo.gpi
             reverseTrafo.gpi
         main_checkParsevalTheorem.py
         main_fourierPair_Gauss.py
         main_fourierPair_sombrero.py
         main_FWDTrafo.py
         main_FWDTrafo_eRMS_Gauss.py
         main_FWDTrafo_eRMS_sombrero.py
         main_FWDTrafo_eRMS_sombrero_varyN.py
         main_selfReciprocality.py
     src/
         convolveRadiallySymmetricFunctions.py
         dataIO.py
         discreteFourierBesselTrafo.py
         irradiationSourceProfile.py


4. LICENSE

BSD 3-Clause License


5. ACKNOWLEDGEMENTS

O. Melchert acknowledges support from the VolkswagenStiftung within the
Nieders\"achsisches Vorab program in the framework of the project Hybrid
Numerical Optics (Grant ZN 3061). 

