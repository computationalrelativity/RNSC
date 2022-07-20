// Minimal working example for making use of this library

#include "RNS.h"

#ifndef PRINTDATA
#define PRINTDATA (0)
#endif

#ifndef TESTINTERP
#define TESTINTERP (0)
#endif

int main(int argc, char* argv[]) {
  
  char * inputfile = NULL;
  if(argc == 2){
    inputfile = argv[1];
    printf("Input file: %s \n", inputfile);
    RNS_params_set_inputfile(inputfile);
  } else {
    RNS_params_set_default();
  }
  
  /*
    Sans parameter injection the following is equivalent to prior all NULL call
  */
  ini_data *data = RNS_make_initial_data();
  
#if(PRINTDATA)
  RNS_data_print(data);
#endif
  
#if (TESTINTERP)
  
  // example of interpolating lapse and rho to a few nodes
  int imin[3] = {0, 0, 0};
  int imax[3] = {2, 2, 2};
  int n[3] = {2, 2, 2};

  int sz = (n[0]) * (n[1]) * (n[2]);

  double *lapse = (double *) malloc(sz * sizeof(double));
  double *rho = (double *) malloc(sz * sizeof(double));
  double *tmp = (double *) malloc(sz * sizeof(double));
  
  double x[2] = {0., -5.0};
  double y[2] = {0., +5.0};
  double z[2] = {0., -5.0};

  RNS_Cartesian_interpolation
    (data,
     imin,
     imax,
     n,
     x,   
     y,
     z,
     lapse,          // lapse
     tmp, tmp, tmp,  // betax, betay, betaz,
     tmp, tmp, tmp,  // adm_gxx, adm_gxy, adm_gxz,
     tmp, tmp, tmp,  // adm_gyy, adm_gyz, adm_gzz,
     tmp, tmp, tmp,  // adm_Kxx, adm_Kxy, adm_Kxz,
     tmp, tmp, tmp,  // adm_Kyy, adm_Kyz, adm_Kzz,
     rho,  
     tmp,            // grhd_epsl,
     tmp, tmp, tmp,  // grhd_vx, grhd_vy, grhd_vz,
     tmp, tmp, tmp,  // grhd_ux, grhd_uy, grhd_uz,
     tmp             //grhd_p
     );

  printf("Interpolation resuts:\n");
  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      for (int i = imin[0]; i < imax[0]; ++i) {
	int flat_ix = i + n[0] * j + n[0]*n[1]*k;
	printf("(x,y,z) = (%4.8f,%4.8f,%4.8f) flat_index = %d\n",
	       x[i],y[j],z[k], flat_ix);
	printf("  lapse = %4.8f\n",lapse[flat_ix]);
	printf("  rho   = %4.8f\n",rho[flat_ix]);
      }
    }
  }
  
  free(lapse);
  free(rho);
  free(tmp);
  
#endif
  
  // make sure to take care of any internal memory that must be freed!
  RNS_finalise(data);
  

  return 0;
}
