#!/bin/bash
# ###########################################################
# RUN RNS TO COMPUTE INITIAL DATA MODEL IN AU AND BU SEQUENCE
# ###########################################################
# Uses RNS version 2/rotstar
# 
# Resolution is:  
#                MDIVxSDIV = 301x601
# and must be set in the makefile
# Accuracy can be set in rotstar.c with the parameter 
#                accuracy = 1e-6
# Other parameters related to accuracy are LMAX and MMAX
# (see in const.h)  
# 
# Standard dimensionless units are used
# It takes about 15' on standard desk machine.
# 

# ###########################################################
# CHOOSE EXE
# ###########################################################

EXE=rotstar/rotstar.x

# ###########################################################
# COMPUTE MODELS IN SEQUENCES AU 
# ###########################################################

Aec=("1.4440e-03" "1.3000e-03" "1.1870e-03" "1.0740e-03" "9.6100e-04" "8.6300e-04")
Arr=("1.0000e+00" "9.1900e-01" "8.5200e-01" "7.8000e-01" "6.9800e-01" "5.7500e-01")
nameA=("AU0.ini" "AU1.ini" "AU2.ini" "AU3.ini" "AU4.ini" "AU5.ini")

numA=${#Aec}
for((i=0;i<$numA;i++)); do
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" ${nameA[${i}]} 
    $EXE -e ${Aec[${i}]} -a ${Arr[${i}]} -s yes -O ${nameA[${i}]}
done

# ###########################################################
# COMPUTE MODELS IN SEQUENCE BU
# ###########################################################

Bec=("1.4440e-03" "1.4440e-03" "1.4440e-03" "1.4440e-03" "1.4440e-03" "1.4440e-03" "1.4440e-03" "1.4440e-03" "1.4440e-03" "1.4440e-03")
Brr=("1.0000e+00" "9.5000e-01" "9.0000e-01" "8.5000e-01" "8.0000e-01" "7.5000e-01" "7.0000e-01" "6.5000e-01" "6.0000e-01" "5.8000e-01")
nameB=("BU0.ini" "BU1.ini" "BU2.ini" "BU3.ini" "BU4.ini" "BU5.ini" "BU6.ini" "BU7.ini" "BU8.ini" "BU9.ini")

numB=${#Bec}
for((i=0;i<$numB;i++)); do
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" ${nameB[${i}]} 
    $EXE -e ${Bec[${i}]} -a ${Arr[${i}]} -s yes -O ${nameB[${i}]}
done

exit 0

