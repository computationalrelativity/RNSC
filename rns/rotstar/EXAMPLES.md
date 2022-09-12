# EXAMPLES

## EXAMPLE 1

```
./rotstar.x
```

Output:

```
>>>     RNS Code - start
SETUP :: MDIVxSDIV      =       301x601
SETUP :: EOS            =       poly
SETUP :: EOS_N          =       1.000000
SETUP :: EOS_K          =       100.000000
SETUP :: e_c [c=G=Ms=1] =       0.00144
SETUP :: r_p/r_e        =       1
SETUP :: e_c [poly]     =       0.144
SETUP :: rho_c [poly]   =       0.127694
SETUP :: p_c [poly]     =       0.0163058
SETUP :: h_c [poly]     =       0.227445
>>>     Spherical star guessed
>>>     Iterating equilibrium model...
accuracy = 1.000e+00
accuracy = 1.741e-06
accuracy = 9.512e-06
accuracy = 7.708e-06
accuracy = 5.546e-06
accuracy = 2.348e-06
>>>     ...Done
>>>     Equilibrium quantites
EQUIL :: ratio          e_15            M               M_0             r_star          spin            Omega_K         I               J/M^2
EQUIL :: 1.000000e+00   1.440000e-03    1.398988e+00    1.504785e+00    9.590188e+00    0.000000e+00    3.982622e-02    0.000000e+00    0.000000e+00
>>>     RNS Code - now exiting
```

## EXAMPLE 2

```
./rotstar.x -e 1.1870e-03 -a 8.5200e-01
```

Output:

```
>>>     RNS Code - start
SETUP :: MDIVxSDIV      =       301x601
SETUP :: EOS            =       poly
SETUP :: EOS_N          =       1.000000
SETUP :: EOS_K          =       100.000000
SETUP :: e_c [c=G=Ms=1] =       0.001187
SETUP :: r_p/r_e        =       0.852
SETUP :: e_c [poly]     =       0.1187
SETUP :: rho_c [poly]   =       0.107207
SETUP :: p_c [poly]     =       0.0114933
SETUP :: h_c [poly]     =       0.194261
>>>     Spherical star guessed
>>>     Iterating equilibrium model...
>>>     Ratio=0.900000
accuracy = 1.000e+00
accuracy = 4.454e-02
accuracy = 2.682e-02
accuracy = 2.365e-02
accuracy = 6.425e-03
accuracy = 2.481e-03
accuracy = 8.956e-04
accuracy = 3.355e-04
accuracy = 1.291e-04
accuracy = 5.073e-05
accuracy = 2.039e-05
accuracy = 8.360e-06
accuracy = 3.489e-06
accuracy = 1.478e-06
accuracy = 1.000e+00
accuracy = 2.331e-02
accuracy = 1.070e-02
accuracy = 1.167e-02
accuracy = 3.157e-03
accuracy = 1.203e-03
accuracy = 4.247e-04
accuracy = 1.594e-04
accuracy = 6.233e-05
accuracy = 2.526e-05
accuracy = 1.060e-05
accuracy = 4.573e-06
accuracy = 2.020e-06
>>>     ...Done
>>>     Equilibrium quantites
EQUIL :: ratio          e_15            M               M_0             r_star          spin            Omega_K         I               J/M^2
EQUIL :: 8.520000e-01   1.187000e-03    1.407580e+00    1.506335e+00    1.078968e+01    1.656386e-03    3.358128e-02    4.820630e+02    3.519024e-01
>>>     RNS Code - now exiting
```

## EXAMPLE 3

```
./rotstar.x -q tab -f ../eos/eosC -e 2e15 -a 5.68311e-01
```

Output

```
>>>     RNS Code - start
SETUP :: MDIVxSDIV      =       301x601
SETUP :: EOS            =       tab
SETUP :: EOS_tab        =       ../eos/eosC
SETUP :: e_c [gr/cm^3]  =       2e+15
SETUP :: r_p/r_e        =       0.568311
SETUP :: e_c [poly]     =       2
SETUP :: rho_c [poly]   =       1.59926
SETUP :: p_c [poly]     =       0.586518
SETUP :: h_c [poly]     =       0.480771
>>>     Spherical star guessed
>>>     Iterating equilibrium model...
>>>     Ratio=0.900000
accuracy = 1.000e+00
accuracy = 4.085e-02
accuracy = 6.791e-02
accuracy = 3.879e-02
accuracy = 7.563e-03
accuracy = 7.231e-03
accuracy = 3.007e-03
accuracy = 1.524e-03
accuracy = 8.325e-04
accuracy = 4.056e-04
accuracy = 2.297e-04
accuracy = 1.201e-04
accuracy = 6.671e-05
accuracy = 3.658e-05
accuracy = 2.019e-05
accuracy = 1.125e-05
accuracy = 6.236e-06
accuracy = 3.484e-06
accuracy = 1.941e-06
accuracy = 1.085e-06
>>>     Ratio=0.800000
accuracy = 1.000e+00
accuracy = 4.738e-02
accuracy = 6.372e-02
accuracy = 3.698e-02
accuracy = 7.368e-03
accuracy = 8.382e-03
accuracy = 3.258e-03
accuracy = 1.785e-03
accuracy = 1.008e-03
accuracy = 5.309e-04
accuracy = 3.340e-04
accuracy = 1.912e-04
accuracy = 1.185e-04
accuracy = 7.139e-05
accuracy = 4.361e-05
accuracy = 2.669e-05
accuracy = 1.628e-05
accuracy = 9.980e-06
accuracy = 6.100e-06
accuracy = 3.738e-06
accuracy = 2.288e-06
accuracy = 1.402e-06
>>>     Ratio=0.700000
accuracy = 1.000e+00
accuracy = 5.607e-02
accuracy = 7.490e-02
accuracy = 4.744e-02
accuracy = 7.887e-03
accuracy = 1.287e-02
accuracy = 4.372e-03
accuracy = 3.110e-03
accuracy = 1.786e-03
accuracy = 1.148e-03
accuracy = 7.858e-04
accuracy = 5.118e-04
accuracy = 3.497e-04
accuracy = 2.323e-04
accuracy = 1.566e-04
accuracy = 1.049e-04
accuracy = 7.044e-05
accuracy = 4.731e-05
accuracy = 3.174e-05
accuracy = 2.132e-05
accuracy = 1.431e-05
accuracy = 9.609e-06
accuracy = 6.452e-06
accuracy = 4.332e-06
accuracy = 2.908e-06
accuracy = 1.953e-06
accuracy = 1.311e-06
>>>     Ratio=0.600000
accuracy = 1.000e+00
accuracy = 6.725e-02
accuracy = 8.923e-02
accuracy = 6.401e-02
accuracy = 1.047e-02
accuracy = 2.091e-02
accuracy = 6.786e-03
accuracy = 6.291e-03
accuracy = 3.721e-03
accuracy = 2.949e-03
accuracy = 2.105e-03
accuracy = 1.571e-03
accuracy = 1.148e-03
accuracy = 8.449e-04
accuracy = 6.192e-04
accuracy = 4.549e-04
accuracy = 3.336e-04
accuracy = 2.449e-04
accuracy = 1.797e-04
accuracy = 1.319e-04
accuracy = 9.673e-05
accuracy = 7.097e-05
accuracy = 5.206e-05
accuracy = 3.819e-05
accuracy = 2.802e-05
accuracy = 2.055e-05
accuracy = 1.507e-05
accuracy = 1.106e-05
accuracy = 8.111e-06
accuracy = 5.950e-06
accuracy = 4.364e-06
accuracy = 3.201e-06
accuracy = 2.348e-06
accuracy = 1.723e-06
accuracy = 1.264e-06
accuracy = 1.000e+00
accuracy = 2.309e-02
accuracy = 2.296e-02
accuracy = 2.236e-02
accuracy = 6.413e-03
accuracy = 7.553e-03
accuracy = 3.188e-03
accuracy = 2.780e-03
accuracy = 1.778e-03
accuracy = 1.435e-03
accuracy = 1.037e-03
accuracy = 7.967e-04
accuracy = 5.912e-04
accuracy = 4.467e-04
accuracy = 3.342e-04
accuracy = 2.513e-04
accuracy = 1.884e-04
accuracy = 1.415e-04
accuracy = 1.062e-04
accuracy = 7.971e-05
accuracy = 5.982e-05
accuracy = 4.490e-05
accuracy = 3.369e-05
accuracy = 2.528e-05
accuracy = 1.897e-05
accuracy = 1.424e-05
accuracy = 1.069e-05
accuracy = 8.019e-06
accuracy = 6.017e-06
accuracy = 4.516e-06
accuracy = 3.389e-06
accuracy = 2.543e-06
accuracy = 1.908e-06
accuracy = 1.432e-06
accuracy = 1.075e-06
>>>     ...Done
>>>     Equilibrium quantites
WARNING :: Unit Conversion for tab EOS models not yet implemented !!!
WARNING :: Units for equilibrium quantities are given below,
WARNING :: Units for the fields should be the polytropic dimensionless ones...
EQUIL :: ratio  e_15    M       M_0     r_star  spin    Omega_K I       J/M^2
EQUIL ::        [g/cm^3] [sun]  [sun]   [km]    [s-1]   [s-1]   [g cm^2]
EQUIL :: 0.5683 2.0000  2.1340  2.4350  14.3283 9637.3859       9632.2914       2.6537  0.6390
>>>     RNS Code - now exiting
```
