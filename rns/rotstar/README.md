# ROTSTAR

This program uses the routines of RNS (v2) to compute equilibrium models of
uniformly rotating neutron stars described either by polytropic or
realistic EOS (tables). The inputs are the the central energy density
and the ratio between the polar the and equatorial radius. For
polytropes also the polytropic index and the constant must be given. 
              
The RNS routines use polytropic units: c=G=K=1, see
Cook, Shapiro & Teukolsky ApJ 398:203-223 (1992).
In the original code, realistic EOS models computation required input
in CGS units, while polytropes models in polytropic units. Ouput is
consistent.

This program uses dimensionless units c=G=Msun=1 for I/O of polytropes
while I/O for realistic EOS still follows the original (to be fixed).
     
Quick help to run

```
 ./rotstar -h
```

Main quantities are output at screen.
Output 2D arrays can be output in a ASCII file                        

 