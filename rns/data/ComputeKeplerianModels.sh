#!/bin/bash
# ########################################################
# Run rns to compute keplerian sequences
# ########################################################

# ########################################################
# Exe and options
# ########################################################
#RUN=../rns.v1.1d/rns -t kepler -p 2 -d 0 -q tab
RUN=../rns.v2.0/kepler.x -p 2 -q tab

# Set EOS for -f option
EOS=("./eos/SFHo_25-Sept-2017.rns" "./eos/LS220B0v2.rns" "./eos/SLy4.rns" "./eos/2B.rns")

# Set output
OUT=("SFHo_kepler.dat" "LS220B0v2_kepler.dat" "SLy4_kepler.dat" "2B_kepler.dat")

# Set first/last energy for -e -l options and number of models in sequence
ECF="5e14" # 1e15 ~ 2Mo
ECL="3e15" #
N="10"

# ########################################################
# Run for each EOS
# ########################################################
NEOS=${#EOS}
for((i=0;i<$numA;i++)); do
    echo ">>>" ${EOS[${i}]}
    $RUN -f ${EOS[${i}]} -e $ECF -l $ECL -n $N >  ${OUT[${i}]}
done

exit 0
