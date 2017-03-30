function singleRun {
N=$1
A0=$2
E0=$3
python main_polConvVerification.py $N $A0 $E0 > ./data/polConvTest_N${N}_A0${A0}_E0${E0}.dat
}

#A0=0.01
#N=100
#for E0 in 0.002 0.004 0.008 0.01 0.02 0.04 0.08 0.1 0.2 0.4 0.8; 
#do
#        singleRun $N $A0 $E0
#done

#A0=0.2
#N=40
#for E0 in 0.002 0.004 0.008 0.01 0.02 0.04 0.08 0.1 0.2 0.4 0.8; 
#do
#        singleRun $N $A0 $E0
#done

#A0=0.2
#E0=0.04
#python main_polConvVerification_varyN.py $A0 $E0 > ./data/polConvTest_varyN_A0${A0}_E0${E0}.dat

A0=0.01
E0=0.01
python main_polConvVerification_varyN.py $A0 $E0 > ./data/polConvTest_varyN_A0${A0}_E0${E0}.dat

