# Test output parsing

import os, re
import numpy as np

def grep_for(expr, text):
  out = []
  for line in text.splitlines():
    match = re.match(expr, line)
    if match is not None:
      out.append(match.group(1))
  return out

def grep(expr, text):
  out = []
  for line in text.splitlines():
    if re.search(expr, line):
        out.append(line)
  return out

# output samples

OUTPUT = {}

# v1

OUTPUT['v1a']="""
    ../rns.v1.1d/eos/eosC,  MDIVxSDIV=151x301,  accuracy=1e-05
  ------------------------------------------
  2.00000e+00  e_c           (10^15 gr/cm^3)
  2.13269e+00  M             (M_sun)
  2.43374e+00  M_0           (M_sun)
  1.39567e+01  R_e           (km)
  9.61439e-01  Omega         (10^4 s^-1)
  1.00250e+00  Omega_p       (10^4 s^-1)
  1.09776e-01  T/W
  2.89532e+00  cJ/GM_sun^2
  2.64660e+00  I             (10^45 gr cm^2)
  2.75872e+02  Phi_2         (10^42 gr cm^2)
  1.14118e-01  h+            (km)
  1.19727e+01  h-            (km)
  5.90997e-01  Z_p       
 -2.96127e-01  Z_f        
  1.66219e+00  Z_b       
  7.22500e-01  omega_c/Omega
  1.04187e+01  r_e           (km)
  5.90000e-01  r_p/r_e
  ------------------------------------------
"""

OUTPUT['v1b']="""
-------------------------------------------------------------------------------
../rns.v1.1d/eos/eosC,  MDIVxSDIV=151x301,  accuracy=1e-05
-------------------------------------------------------------------------------
   rho_c        M         M_0         R         Omega      Omega_p       T/W           C*J/GM_s^2      I        h_plus       h_minus       Z_p      Z_b        Z_f       omega_c/Omega   r_e     r_ratio     Omega_pa    Omega+    u_phi
-------------------------------------------------------------------------------
6.00000e+14 8.09439e-01 8.47592e-01 1.26934e+01 0.00000e+00 7.24437e-01 0.00000e+00 0.00000e+00    ---     0.00000e+00 0.00000e+00 1.09853e-01 1.09853e-01 1.09853e-01 0.00000e+00 3.12495e-01 1.00000e+00 7.24437e-01 7.24437e-01 3.11499e+00
2.00000e+15 1.79388e+00 2.05282e+00 1.07711e+01 0.00000e+00 1.37974e+00 0.00000e+00 0.00000e+00    ---     5.11020e+00 5.11020e+00 4.02250e-01 4.02250e-01 4.02250e-01 0.00000e+00 2.15342e-01 1.00000e+00 1.28656e+00 1.27752e+00 6.21434e+00
"""

# kepler

OUTPUT['kepler']="""
../rns.v1.1d/eos/eosA,  MDIVxSDIV=151x301
ratio	e_15	M	M_0	r_star	spin	Omega_K	I	J/M^2
	g/cm^3	sun	sun	km	s-1	s-1	g cm^2	

1.000 	 1.0 	0.807 	0.859 	9.936 	 0.0 	10451.7 	0.00 	0.000 
0.950 	 1.0 	0.831 	0.885 	10.209 	2771.5 	10157.3 	0.52 	0.238 
0.900 	 1.0 	0.854 	0.910 	10.501 	3897.4 	9902.0 	0.56 	0.341 
0.850 	 1.0 	0.879 	0.937 	10.831 	4736.4 	9625.8 	0.60 	0.422 
0.800 	 1.0 	0.905 	0.966 	11.207 	5413.2 	9315.8 	0.65 	0.491 
0.750 	 1.0 	0.934 	0.997 	11.643 	5970.2 	8958.7 	0.71 	0.551 
0.700 	 1.0 	0.962 	1.028 	12.159 	6420.4 	8535.7 	0.76 	0.603 
0.650 	 1.0 	0.987 	1.055 	12.790 	6757.0 	8018.2 	0.82 	0.644 
0.600 	 1.0 	1.004 	1.073 	13.596 	6951.0 	7365.2 	0.85 	0.668 
0.550 	 1.0 	1.003 	1.073 	14.681 	6947.6 	6539.1 	0.85 	0.668 
0.575 	 1.0 	1.006 	1.076 	14.091 	6977.7 	6977.7 	0.86 	0.672"""

# simple test

OUTPUT['rotstar1']="""
# ratio		e_15		M		M_0		r_star		spin		Omega_K		I		J/M^2
1.000000e+00	1.440000e-03	1.398988e+00	1.504785e+00	9.590188e+00	0.000000e+00	3.982622e-02	0.000000e+00	0.000000e+00
2.274450e-01
"""

# rotstar

OUTPUT['rotstar']="""
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
"""


# ---------------------------------
print()
print(OUTPUT['v1a'])

t = [s for s in OUTPUT['v1a'].splitlines() if s] # rm empty lines
t = t[2:-1] # rm comments
t = [s.strip() for s in t] # rm leading/ending spaces
t = [("".join(s.split()[1:]),float(s.split()[0])) for s in t] # split every line in value/key (list of tuple)
t = dict(t)
#print(t)

def parse1(out):
    out = [s for s in out.splitlines() if s] 
    out = [s.strip() for s in out[2:-1]]
    return dict([("".join(s.split()[1:]),float(s.split()[0])) for s in out])

print(parse1(OUTPUT['v1a']))


# ---------------------------------
print()
print(OUTPUT['v1b'])

t = [s for s in OUTPUT['v1b'].splitlines() if s] # rm empty lines
keys = t[3].split() # keys
t = t[5:] # data
t = [s.split() for s in t]
t = list(np.concatenate(t))
t = [np.nan  if s == '---' else s for s in t]
#t = np.array([float(s) for s in t]).reshape(2,20) #(len(t)/len(keys),len(keys)))
t = np.array([float(s) for s in t]).reshape( (len(t)//len(keys),len(keys)) )
t= [(k, t[:,i]) for i,k in enumerate(keys)]
#print(dict(t))

def parse2(out, headerlines=5, keysline=3):
    out = [s for s in out.splitlines() if s] 
    keys = out[keysline].split()
    out = list(np.concatenate([s.split() for s in out[headerlines:]]))
    out = [np.nan  if s == '---' else s for s in out]
    out = np.array([float(s) for s in out]).reshape((len(out)//len(keys),len(keys)))
    return dict([(k, out[:,i]) for i,k in enumerate(keys)])

print(parse2(OUTPUT['v1b']))


# ---------------------------------
print()
print(OUTPUT['rotstar1'])

def parse3(out,headerlines=1,keysline=0):
    # tabulated EOS can have 2 rows of headers
    out = [s for s in out.splitlines() if s] 
    keys = out[keysline].split()[1:] # rm #
    #keys = out[keysline].split()
    vals = [float(s) for s in out[headerlines].split()]
    return dict(zip(keys, vals))

print(parse3(OUTPUT['rotstar1']))


# ---------------------------------
print()
print(OUTPUT['rotstar'])

o = grep("EQUIL ::*",OUTPUT['rotstar'])
o = [s.replace("EQUIL ::","").strip() for s in o]
print(o) # -> form ready for parse3 
