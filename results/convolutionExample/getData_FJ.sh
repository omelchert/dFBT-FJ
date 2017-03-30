
function Gaussian {
T=4.0
N=40
A0=0.25

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95.mco \
        $T $N 0.0 0.0 0.0 $A0 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95_FJ_T${T}_N${N}_G_A0${A0}.prof 

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.90.mco \
        $T $N 0.0 0.0 0.0 $A0 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.90_FJ_T${T}_N${N}_G_A0${A0}.prof 

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.70.mco \
        $T $N 0.0 0.0 0.0 $A0 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.70_FJ_T${T}_N${N}_G_A0${A0}.prof 

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.10.mco \
        $T $N 0.0 0.0 0.0 $A0 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.10_FJ_T${T}_N${N}_G_A0${A0}.prof 
}


function flatTop {
T=4.0
N=80
R0=0.4
A0=0.1

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95.mco \
        $T $N 0.0 0.0 $R0 $A0 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95_FJ_T${T}_N${N}_FT_R0${R0}_A0${A0}.prof 

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.90.mco \
        $T $N 0.0 0.0 $R0 $A0 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.90_FJ_T${T}_N${N}_FT_R0${R0}_A0${A0}.prof 

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.70.mco \
        $T $N 0.0 0.0 $R0 $A0 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.70_FJ_T${T}_N${N}_FT_R0${R0}_A0${A0}.prof 

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.10.mco \
        $T $N 0.0 0.0 $R0 $A0 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.10_FJ_T${T}_N${N}_FT_R0${R0}_A0${A0}.prof 
}


function donut {
T=4.0
N=150
R0=0.35
A0=0.07
R1=0.7
A1=0.07

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95.mco \
        $T $N $R0 $A0 $R1 $A1 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95_FJ_T${T}_N${N}_D_R0${R0}_A0${A0}_R1${R1}_A1${A1}.prof 

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.90.mco \
        $T $N $R0 $A0 $R1 $A1 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.90_FJ_T${T}_N${N}_D_R0${R0}_A0${A0}_R1${R1}_A1${A1}.prof 

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.70.mco \
        $T $N $R0 $A0 $R1 $A1 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.70_FJ_T${T}_N${N}_D_R0${R0}_A0${A0}_R1${R1}_A1${A1}.prof 

python main_convolvedResponse_FiskJohnson.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.10.mco \
        $T $N $R0 $A0 $R1 $A1 \
        ./convolvedResponse_semiInf/semiInf_ma0.1_ms10.0_g0.10_FJ_T${T}_N${N}_D_R0${R0}_A0${A0}_R1${R1}_A1${A1}.prof 
}

#Gaussian
#flatTop
#donut

function GaussianSlices {
T=4.0
N=40
A0=0.25

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95.mco \
        $T $N 0.0 0.0 0.0 $A0 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.95_FJ_T${T}_N${N}_G_A0${A0}.prof 

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.90.mco \
        $T $N 0.0 0.0 0.0 $A0 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.90_FJ_T${T}_N${N}_G_A0${A0}.prof 

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.70.mco \
        $T $N 0.0 0.0 0.0 $A0 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.70_FJ_T${T}_N${N}_G_A0${A0}.prof 

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.10.mco \
        $T $N 0.0 0.0 0.0 $A0 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.10_FJ_T${T}_N${N}_G_A0${A0}.prof 
}



function flatTopSlices {
T=4.0
N=80
R0=0.4
A0=0.1

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95.mco \
        $T $N 0.0 0.0 $R0 $A0 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.95_FJ_T${T}_N${N}_FT_R0${R0}_A0${A0}.prof 

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.90.mco \
        $T $N 0.0 0.0 $R0 $A0 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.90_FJ_T${T}_N${N}_FT_R0${R0}_A0${A0}.prof 

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.70.mco \
        $T $N 0.0 0.0 $R0 $A0 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.70_FJ_T${T}_N${N}_FT_R0${R0}_A0${A0}.prof 

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.10.mco \
        $T $N 0.0 0.0 $R0 $A0 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.10_FJ_T${T}_N${N}_FT_R0${R0}_A0${A0}.prof 
}



function donutSlices {
T=4.0
N=150
R0=0.35
A0=0.07
R1=0.7
A1=0.07

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.95.mco \
        $T $N $R0 $A0 $R1 $A1 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.95_FJ_T${T}_N${N}_D_R0${R0}_A0${A0}_R1${R1}_A1${A1}.prof 

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.90.mco \
        $T $N $R0 $A0 $R1 $A1 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.90_FJ_T${T}_N${N}_D_R0${R0}_A0${A0}_R1${R1}_A1${A1}.prof 

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.70.mco \
        $T $N $R0 $A0 $R1 $A1 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.70_FJ_T${T}_N${N}_D_R0${R0}_A0${A0}_R1${R1}_A1${A1}.prof 

python main_WrzSlice.py \
        ./impulseResponse_semiInf/semiInf_ma0.1_ms10.0_g0.10.mco \
        $T $N $R0 $A0 $R1 $A1 \
        > ./WrzSlices/semiInf_ma0.1_ms10.0_g0.10_FJ_T${T}_N${N}_D_R0${R0}_A0${A0}_R1${R1}_A1${A1}.prof 
}

GaussianSlices
flatTopSlices
donutSlices
