function singleRun {
  T=$1
  N=$2
  echo $T $N
  time python main_fourierPair_sombrero.py $T $N > \
          ./data_sombrero/dFBT_fourierPair_sombrero_T${T}_N${N}.dat
}


N=20

for T in 4.0 6.0 8.0 10.0 12.0 14.0 16.0 18.0 20.0 ; 
do
  #singleRun $T $N
  time python main_FWDTrafo_eRMS_sombrero_varyN.py $T \
          > ./data_sombrero/dFBT_FWDTrafo_eRMS_sombrero_varyN_T${T}.dat
done

T=10.0
for N in 5 8 10 20 40 50 100;
do
  singleRun $T $N
  time python main_FWDTrafo_eRMS_sombrero.py $N \
          > ./data_sombrero/dFBT_FWDTrafo_eRMS_sombrero_${N}.dat
done









