
function singleRun {

python main_GreensFuncResponse.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95.mco \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95.prof 

python main_GreensFuncResponse.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.70.mco \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.70.prof 

python main_GreensFuncResponse.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.10.mco \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.10.prof 

}

singleRun
