/***************************************************************************
 * consts.h 
 *
 * This file contains definitions of constants used in the procedures.
 *
 **************************************************************************/

#define TABP 200001               /* MAX POINTS IN TAB e_vs_P_filename */

#define DM (1.0/(MDIV-1.0))          /* spacing in mu direction */ 
#define RDIV 900                     /* grid point in RK integration */ 
#define SMAX 0.9999                  /* maximum value of s-coordinate */  
#define DS (SMAX/(SDIV-1.0))         /* spacing in s-direction */
#define LMAX 24                      /* max. term in Legendre poly. */
#define C 2.9979e10                  /* speed of light in vacuum [cm s^-1]*/
#define G 6.6732e-8                  /* gravitational constant [g^-1 cm^3 s^-2] */ 
#define KAPPA (1.0e-15*C*C/G)        /* scaling factor [cm^-1] (assuming 1e-15 are [g]) */
#define KSCALE (KAPPA*G/(C*C*C*C))   /* another scaling factor */
#define MSUN 1.987e33                /* Mass of Sun [g] */
#define SQ(x) ((x)*(x))              /* square macro */
#define MB 1.66e-24                  /* baryon mass [g] */
#define RMIN 1.0e-15                 /* use approximate TOV equations when
                                        computing spherical star and r<RMIN */
#define UNUSED (-1.11e30)            /* Used in the Ridder zero-finder */ 
#ifndef PI
#define PI 3.1415926535              /* what else */
#endif

#define M_MAX 12   /* Maximum l=m mode considered */
#define BOUND 1    /* Excludes last grid point, since we set everything to zero there */
#define MAXIT 500  /* Maximum number of iterations in rtsec_G */ 
#include <float.h>
#include <limits.h>

#ifndef VERBOSE
#define VERBOSE (0)
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-15
#endif
