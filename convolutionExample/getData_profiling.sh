
T=4.0
N=60
R0=0.4
A0=0.1

for N in 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160;
do
  echo $N
  python -m cProfile -o profFJ_N${N}.dat main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95.mco \
        $T $N 0.0 0.0 $R0 $A0 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95_FJ_T${T}_N${N}_FT_R0${R0}_A0${A0}.prof 

  python main_pstat.py profFJ_N${N}.dat > profFJ_N${N}.time 
done     

python -m cProfile -o "profCB.dat" main_convolvedResponse_CreeBones.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95.mco \
        0.0 0.0 $R0 $A0 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95_CB_FT_R0${R0}_A0${A0}.prof 

