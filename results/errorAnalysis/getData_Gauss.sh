function singleRun {
  T=$1
  N=$2
  echo $T $N
  time python main_fourierPair_Gauss.py $T $N > \
          ./data_Gauss/dFBT_fourierPair_Gauss_T${T}_N${N}.dat
}


N=20

for T in 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 ; 
do
  singleRun $T $N
done

for N in 5 10 20 40 50 100;
do
  time python main_FWDTrafo_eRMS_Gauss.py $N \
          > ./data_Gauss/dFBT_FWDTrafo_eRMS_Gauss_${N}.dat
done
