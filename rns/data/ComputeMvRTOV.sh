#!/bin/bash
# ########################################################
# BASH SCRIPT TO COMPUTE M Vs R DIGRAM FOR NON-ROTATING NS
# ########################################################
# Uses RNS v1
#

# ########################################################
# CHOOSE EXE
# ########################################################

EXE=rns.v1.1d/rns

# ########################################################
# CHOOSE EOS
# ########################################################

EOS="./eos/eosFPS ./eos/eosSly ./eos/eosAPR ./eos/eosFPSfit ./eos/eosSlyfit ./eos/eosAPRfit"

# max central energy for each EOS
maxe=(3.8e15 3.8e15 3.2e15 3.8e15 3.8e15 3.2e15)

# ########################################################
# COMPUTE SEQUENCE 
# ########################################################

echo $EOS

c=0

for e in $EOS
do
    n=${e##*/} 
    echo 'sequence for eos ' $n ' -- max energy ' ${maxe[c]}
    
    #
    $EXE -q tab -f $e -t static -e 3e14 -l ${maxe[c]} -n 51 -p 2 -d 0 > MvR_$n.dat
    
    awk 'NR>6' MvR_$n.dat | awk '{print $4 " " $2}' > mr$c
    awk 'NR>6' MvR_$n.dat | awk '{print $1}'        > ec$c	
    let c++	
    
done

echo "# $(date)" > MvR.dat
echo "# R vs M for $EOS" >> MvR.dat
echo "# data computed with RNS code" >> MvR.dat
paste mr0 mr1 mr2 mr3 mr4 mr5 >> MvR.dat
rm mr0 mr1 mr2 mr3 mr4 mr5 

echo "# $(date)" > ec.dat
echo "# ec for $EOS" >> ec.dat
echo "# data computed with RNS code" >> ec.dat
paste ec0 ec1 ec2 ec3 ec4 ec5 >> ec.dat
rm ec0 ec1 ec2 ec3 ec4 ec5 

# ########################################################
# GNUPLOT
# ########################################################

#echo "set terminal post eps color enhanced 'Times' 30 dl 3 lw 1.0"
#echo "p 'MvR.dat' u 1:2 w l lt 2 lw 2 lc rgb 'red' t 'FPS'\" > MvR.gp
#echo "'MvR.dat' u 3:4 w l lt 2 lw 2 lc rgb 'blue' t 'Sly4'\" >> MvR.gp
#echo "'MvR.dat' u 5:6 w l lt 2 lw 2 lc rgb 'green' t 'APR'\" >> MvR.gp
#echo "'MvR.dat' u 7:8 w l lt 2 lw 2 lc rgb 'red' t 'FPSfit'\" >> MvR.gp
#echo "'MvR.dat' u 9:10 w l lt 2 lw 2 lc rgb 'blue' t 'Sly4fit'\" >> MvR.gp
#echo "'MvR.dat' u 11:12 w l lt 2 lw 2 lc rgb 'green' t 'APRfit'" >> MvR.gp

exit 0

