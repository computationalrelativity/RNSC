/* RNS.h */

#ifndef RNS_HEADER_H
#define RNS_HEADER_H

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "stdbool.h"

/** Stuff from numerical recipes */

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

/* Moved NR macros IMIN, IMAX, SIGN where needed.
 * Commented out unused macros:
 */
/* static float sqrarg; */
/* #define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg) */
/* static float maxarg1,maxarg2; */
/* #define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?	\ */
/* 		   (maxarg1) : (maxarg2)) */
/* static float minarg1,minarg2; */
/* #define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?	\ */
/* 		   (minarg1) : (minarg2)) */
/* static long lmaxarg1,lmaxarg2; */
/* #define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ? \ */
/* (lmaxarg1) : (lmaxarg2)) */
/* static long lminarg1,lminarg2; */
/* #define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?	\ */
/* 		   (lminarg1) : (lminarg2)) */
/* static double dmaxarg1,dmaxarg2; */
/* #define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\ */
/* 		   (dmaxarg1) : (dmaxarg2)) */
/* static double dminarg1,dminarg2; */
/* #define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?	\ */
/* 		   (dminarg1) : (dminarg2)) */
/* static double dsqrarg; */
/* #define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg) */

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI)
/* ANSI */
void nrerror(char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);
#else /* ANSI */
/* traditional - K&R */
void nrerror();
float *vector();
float **matrix();
float **submatrix();
float **convert_matrix();
float ***f3tensor();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
unsigned char *cvector();
unsigned long *lvector();
void free_vector();
void free_dvector();
void free_ivector();
void free_cvector();
void free_lvector();
void free_matrix();
void free_submatrix();
void free_convert_matrix();
void free_dmatrix();
void free_imatrix();
void free_f3tensor();
#endif /* ANSI */
#endif /* _NR_UTILS_H_ */


/** Constants & parameters */

/* Following grid parameters are now parameters/local variables */
/* #define MDIV 301 */
/* #define SDIV 601 */
/* #define DM (1.0/(MDIV-1.0)) */    /* spacing in mu direction */
/* #define DS (SMAX/(SDIV-1.0)) */   /* spacing in s-direction */
/* #define SMAX (0.9999) */          /* maximum value of s-coordinate */

/* moved where needed (src/RNS_equil.c) */
/* #define RDIV (900) */             /* grid point in RK integration */
/* #define RMIN (1.0e-15) */         /* use approximate TOV equations when
                                        computing spherical star and r<RMIN */

#define C_SPEED (2.9979e10)          /* speed of light in vacuum */
#define G (6.6732e-8)                /* gravitational constant */
#define KAPPA (1.0e-15*C_SPEED*C_SPEED/G)        /* scaling factor */
#define KSCALE (KAPPA*G/(C_SPEED*C_SPEED*C_SPEED*C_SPEED))   /* another scaling factor */
#define MSUN (1.987e33)              /* Mass of Sun */
#define SQ(x) ((x)*(x))              /* square macro */
#define MB (1.66e-24)                /* baryon mass */

/* constants to convert from cgs to cactus units*/
#define cactusM (0.5027e-33)
#define cactusL (0.6770e-5)
#define cactusT (2.0295e5)
#define cactusV (1.0/(cactusL*cactusL*cactusL))

#define TABPTS (2001) /* Warning: fixed max size of EOS tables */

/* moved where needed (src/RNS_equil_util.c) */
/* #define MAXIT (1000) */ /* max iterations for secant */
/* #define ITMAX (100) */ /* max iterations for zbrent */

#ifndef PI
#define PI (3.1415926535)              /* what else */
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-15
#endif
#define EPS DBL_EPSILON


/* verbose levels */
enum{
  print_nothing,
  print_globalprop1, // print global quantities
  print_info,        // print various info
  verbose_levels
};


/** RNS_equil.c */

void rns_make_grid(double *s_gp,
		   double *mu);
void load_eos(char eos_file[],
	      double log_e_table[TABPTS],
	      double log_p_table[TABPTS],
	      double log_h_table[TABPTS],
	      double log_n0_table[TABPTS],
	      int *n_tab);
double e_of_rho0(double rho0, double Gamma_P, double eos_k);
double e_at_p(double pp,
              double log_e_tab[TABPTS],
              double log_p_tab[TABPTS],
              int    n_tab,
              int    *n_nearest_pt,
              char eos_type[],
	      double eos_k,
              double Gamma_P);
double p_at_e(double ee,
              double log_p_tab[TABPTS],
              double log_e_tab[TABPTS],
              int    n_tab,
              int    *n_nearest_pt);
double p_at_h(double hh,
              double log_p_tab[TABPTS],
              double log_h_tab[TABPTS],
              int    n_tab,
              int    *n_nearest_pt);
double h_at_p(double pp,
              double log_h_tab[TABPTS],
              double log_p_tab[TABPTS],
              int    n_tab,
              int    *n_nearest_pt);
double n0_at_e(double ee,
               double log_n0_tab[TABPTS],
               double log_e_tab[TABPTS],
               int    n_tab,
               int    *n_nearest_pt);
void make_center(double e_center,
                 double log_e_tab[TABPTS],
                 double log_p_tab[TABPTS],
                 double log_h_tab[TABPTS],
                 int n_tab,
                 char eos_type[],
		 double eos_k,
                 double Gamma_P,
                 double *p_center,
                 double *h_center);
void comp_values(double *s_gp,
                 double *mu,
                 double r_ratio,
                 double e_surface,
                 double r_e,
                 char   eos_type[],
                 double log_e_tab[TABPTS],
                 double log_n0_tab[TABPTS],
                 int    n_tab,
                 double Omega,
                 double **rho,
                 double **gama,
                 double **alpha,
                 double **omega,
                 double **energy,
                 double **pressure,
                 double **enthalpy,
                 double **velocity_sq,
                 double *Mass,
                 double *Mass_0,
                 double *T,
                 double *W,
                 double *Omega_K,
                 double *R_e,
                 char   rotation_type[],
		 double **Omega_diff,
		 double *J);
double dm_dr_is(double r_is,
                double r,
                double m,
                double p,
                double e_center,
                double p_surface,
                double log_e_tab[TABPTS],
                double log_p_tab[TABPTS],
                int    n_tab,
                int    *n_nearest_pt,
                char   eos_type[],
		double eos_k,
                double Gamma_P);
double dp_dr_is(double r_is,
                double r,
                double m,
                double p,
                double e_center,
                double p_surface,
                double log_e_tab[TABPTS],
                double log_p_tab[TABPTS],
                int    n_tab,
                int    *n_nearest_pt,
                char eos_type[],
		double eos_k,
                double Gamma_P);
double dr_dr_is(double r_is, double r, double m);
void integrate(int    i_check,
               char   eos_type[],
	       double eos_k,
               double e_center,
               double p_center,
               double p_surface,
               double e_surface,
               double Gamma_P,
               double log_e_tab[TABPTS],
               double log_p_tab[TABPTS],
               double log_h_tab[TABPTS],
               int    n_tab,
               double r_is_gp[], //[RDIV+1],
               double lambda_gp[], //[RDIV+1],
               double nu_gp[], //[RDIV+1],
               double *r_is_final,
               double *r_final,
               double *m_final);
void guess(double *s_gp,
           char   eos_type[],
	   double eos_k,
           double e_center,
           double p_center,
           double p_surface,
           double e_surface,
           double Gamma_P,
           double log_e_tab[TABPTS],
           double log_p_tab[TABPTS],
           double log_h_tab[TABPTS],
           int    n_tab,
           double **rho,
           double **gama,
           double **alpha,
           double **omega,
           double *r_e);
void iterate(double *s_gp,
             double *mu,
             char   eos_type[],
	     double eos_k,
             double log_e_tab[TABPTS],
             double log_p_tab[TABPTS],
             double log_h_tab[TABPTS],
             int    n_tab,
             double Gamma_P,
             double r_ratio,
             double h_center,
             double enthalpy_min,
             int    *a_check,
             double accuracy,
             int    print_dif,
             double cf,
             double *r_e_new,
             double **rho,
             double **gama,
             double **alpha,
             double **omega,
             double **energy,
             double **pressure,
             double **enthalpy,
             double **velocity_sq,
             double *Omega,
             char   rotation_type[],
	     double A_diff,
	     double *Omega_e,
	     double **Omega_diff,
             int RNS_lmax);

/** RNS_util.c */
double ***d3tensor(long nrl,
                   long nrh,
                   long ncl,
                   long nch,
                   long ndl,
                   long ndh);
void free_d3tensor(double ***t,
                   long nrl,
                   long nrh,
                   long ncl,
                   long nch,
               	   long ndl,
                   long ndh);
double interp_4(double xp[5],
                double yp[5],
                int    np ,
                double xb);
void grid_interp(double **old,
                 double *s_gp,
                 double *mu,
                 double r_e,
                 int nx,
                 int ny,
                 int nz,
                 double *x_grid,
                 double *y_grid,
                 double *z_grid,
                 int i,
                 int j,
                 int k,
                 double *new_,
                 int sign);
void grid_interp_all( double *s_gp,
                       double *mu,
                      double  r_e,
                      int     nx,
                      int     ny,
                      int     nz,
                      double *x_grid,
                      double *y_grid,
                      double *z_grid,
                      int i,
                      int j,
                      int k,
                      double **nu,
                      double **B,
                      double **alpha,
                      double **omega,
                      double **nu_dr,
                      double **B_dr,
                      double **alpha_dr,
                      double **omega_dr,
                      double **nu_dth,
                      double **B_dth,
                      double **alpha_dth,
                      double **omega_dth,
                      double **rho_0,
                      double **energy,
                      double **pressure,
                      double *nu_c,
                      double *B_c,
                      double *alpha_c,
                      double *omega_c,
                      double *nu_dr_c,
                      double *B_dr_c,
                      double *alpha_dr_c,
                      double *omega_dr_c,
                      double *nu_dth_c,
                      double *B_dth_c,
                      double *alpha_dth_c,
                      double *omega_dth_c,
                      double *rho_0_c,
                      double *energy_c,
                      double *pressure_c,
                      double *distance_c,
		      double **Omega_diff,
		      double *Omega_diff_c);
void grid_interp_new( int k_s,
		      int k_m,
		      int sign,
		      double s_c,
		      double mu_c,
		      double zp,
		      double *s_4,
		      double *mu_4,
		      double **old,
		      double *new_ );
void grid_interp_all_new (double *s_gp,
			  double *mu,
			  double r_e,
			  double xp,  double yp, double zp,
			  double **nu, double **B, double **alpha, double **omega,
			  double **nu_dr, double **B_dr, double **alpha_dr, double **omega_dr,
			  double **nu_dth, double **B_dth, double **alpha_dth, double **omega_dth,
			   double **rho_0, double **energy, double **pressure,
			  double **Omega_diff,
			  double *nu_c, double *B_c, double *alpha_c, double *omega_c,
			  double *nu_dr_c, double *B_dr_c, double *alpha_dr_c, double *omega_dr_c,
			  double *nu_dth_c, double *B_dth_c, double *alpha_dth_c, double *omega_dth_c,
			  double *rho_0_c, double *energy_c, double *pressure_c,
			  double *Omega_diff_c);
void transform_units(char   eos_type[],
		     double  n_P,
		     double  eos_k,
		     double *rho0_center,
		     double *e_center,
		     double *p_center,
		     double *r_e,
		     double **omega,
		     double **energy,
		     double **pressure,
		     double *Mass,
		     double *Mass_0,
		     double *T,
		     double *W,
		     double *Omega,
		     double *Omega_K,
		     double *R_e,
		     double *Omega_e,
		     double **Omega_diff,
		     double *J);


/** RNS_equil_util.c */

void nr_hunt(double xx[], int n, double x, int *jlo);
double interp(double xp[],
              double yp[],
              int    np ,
              double xb,
              int    *n_nearest_pt);
double deriv_s(double **f,int s, int m);
double deriv_ss(double **f,int s, int m);
double deriv_m(double **f,int s, int m);
double deriv_mm(double **f,int s, int m);
double deriv_sm(double **f,int s, int m);
double nr_legendre( int n, double x );
double nr_plgndr(int l, int m, double x);
double rtsec_G( double (*func)(double, double, double),
                double Gamma_P,
                double x1,
                double x2,
                double xacc,
                double ee,
		double eos_k);
double zbrent_diff(double (*func)(double, double, double, double, double,
                                  double, double, double),
                   double r_e,
                   double rho_equator_h,
                   double gama_equator_h,
                   double omega_equator_h,
                   double rho_pole_h,
                   double gama_pole_h,
                   double x1,
                   double x2,
                   double tol,
		   double A_diff);
double zbrent_rot(double (*func)(double, double, double, double, double,
                                  double, double, double),
                   double r_e,
                   double rhogp,
                   double omegagp,
                   double sgp,
                   double mugp,
                   double Omega_c,
                   double x1,
                   double x2,
    		   double tol,
		   double A_diff);
double diff_rotation(double x,
                     double r_e,
                     double rho_equator_h,
                     double gama_equator_h,
                     double omega_equator_h,
                     double rho_pole_h,
                     double gama_pole_h,
		     double A_diff);
double rotation_law(double x,
                     double r_e,
                     double rhogp,
                     double omegagp,
                     double sgp,
                     double mugp,
                     double Omega_c,
		     double A_diff);


/** RNS_rnsid_util.c */

double ***d3tensor(long nrl,
                   long nrh,
                   long ncl,
                   long nch,
                   long ndl,
                   long ndh);
void free_d3tensor(double ***t,
                   long nrl,
                   long nrh,
                   long ncl,
                   long nch,
               	   long ndl,
                   long ndh);
double interp_4(double xp[5],
                double yp[5],
                int    np ,
                double xb);
void grid_interp(double **old,
                 double *s_gp,
                 double *mu,
                 double r_e,
                 int nx,
                 int ny,
                 int nz,
                 double *x_grid,
                 double *y_grid,
                 double *z_grid,
                 int i,
                 int j,
                 int k,
                 double *new_,
                 int sign);
void grid_interp_all(double *s_gp,
		     double *mu,
		     double  r_e,
		     int     nx,
		     int     ny,
		     int     nz,
		     double *x_grid,
		     double *y_grid,
		     double *z_grid,
		     int i,
		     int j,
		     int k,
		     double **nu,
		     double **B,
		     double **alpha,
		     double **omega,
		     double **nu_dr,
		     double **B_dr,
		     double **alpha_dr,
		     double **omega_dr,
		     double **nu_dth,
		     double **B_dth,
		     double **alpha_dth,
		     double **omega_dth,
		     double **rho_0,
		     double **energy,
		     double **pressure,
		     double *nu_c,
		     double *B_c,
		     double *alpha_c,
		     double *omega_c,
		     double *nu_dr_c,
		     double *B_dr_c,
		     double *alpha_dr_c,
		     double *omega_dr_c,
		     double *nu_dth_c,
		     double *B_dth_c,
		     double *alpha_dth_c,
		     double *omega_dth_c,
		     double *rho_0_c,
		     double *energy_c,
		     double *pressure_c,
		     double *distance_c,
		     double **Omega_diff,
		     double *Omega_diff_c);
void print_arrays_check(
                      int     nx,
                      int     ny,
                      double *x_grid,
                      double *y_grid,
                      int     z_print,
                      double nu_c,
                      double B_c,
                      double alpha_c,
                      double omega_c,
                      double nu_dr_c,
                      double B_dr_c,
                      double alpha_dr_c,
                      double omega_dr_c,
                      double nu_dth_c,
                      double B_dth_c,
                      double alpha_dth_c,
                      double omega_dth_c,
                      double rho_0_c,
                      double energy_c,
                      double pressure_c);
/*
  void print_id(
  int     rdiv,
  int     thdiv,
  double *rc_grid,
  double *thc_grid,
  double **nu_c,
  double **B_c,
  double **alpha_c,
  double **omega_c,
  double **nu_dr,
  double **B_dr,
  double **alpha_dr,
  double **omega_dr,
  double **nu_dth,
  double **B_dth,
  double **alpha_dth,
  double **omega_dth,
  double **rho_0_c,
  double **e_int_c,
  double **ut_c,
  double **uphi_c);
*/
void transform_units(char   eos_type[],
		     double  n_P,
		     double  eos_k,
		     double *rho0_center,
		     double *e_center,
		     double *p_center,
		     double *r_e,
		     double **omega,
		     double **energy,
		     double **pressure,
		     double *Mass,
		     double *Mass_0,
		     double *T,
		     double *W,
		     double *Omega,
		     double *Omega_K,
		     double *R_e,
		     double *Omega_e,
		     double **Omega_diff,
		     double *J);

/** RNS_extra.c */

/* stuff for params */
#define NPARAMS (100) /* Warning: max pars is fixed */
#define STRLEN (256)
#define STREQL(s1,s2) ((strcmp((s1),(s2))==0))
#define ERROR(s) {printf(s); exit(0);}
#define DBGSTOP(s) {printf("DEBUGSTOP: %s",s); exit(0);}
enum{
  INTEGER,
  REAL,
  STRING,
  N_PAR_TYPES,
};
static const char* str_par_type[N_PAR_TYPES] = {
  "INTEGER", "REAL", "STRING",
};
typedef struct {
  int type[NPARAMS];
  char key[NPARAMS][STRLEN];
  char val[NPARAMS][STRLEN];
  int n;
  //parameters *next;
} parameters;

void params_alloc();
void params_free();
void params_read(char *fname);
void params_write(char * fname);
double params_get_real(char * key);
int params_get_int(char * key);
char *params_get_str(char * key);
void params_setadd(char * key, int type, char *val, int addpar);
void params_set_int(char * key, int val);
void params_set_bool(char * key, bool val);
void params_set_real(char * key, double val);
void params_set_str(char * key, char* val);
void params_add_int(char * key, int val);
void params_add_real(char * key, double val);
void params_add_str(char * key, char *val);
void make_output_dir();

void print_global_quantities(char * eos_type,
			     char * rotation_type,
			     double rho_center,
			     double a_ratio,
			     double e_center,
			     double Mass,
			     double Mass_0,
			     double R_e,
			     double J,
			     double T,
			     double W,
			     double Omega,
			     double Omega_K,
			     double Omega_e);

void contract_vect_v2(double vlx, double vly, double vlz,
		      double gxx, double gxy, double gxz,
		      double gyy, double gyz, double gzz,
		      double *vx, double *vy, double *vz,
		      double *v2);


/** RNS.c */

// Expose for C++ intefacing
#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
    int SDIV;
    int MDIV;

    double *s_gp;
    double *mu;
    double r_e;

    double **rho_0;    /* rest mass density */
    double **energy;
    double **pressure;
    double **Omega_diff; /* Diff. ang. vel. */

    double **nu;       /* potential nu */
    double **B;        /* potential B */
    double **alpha;    /* potential alpha */
    double **omega;    /* potential omega */

    double **nu_dr;    /* r-der. in s-coord. of nu */
    double **B_dr;     /* r-der. in s-coord. of B */
    double **alpha_dr; /* r-der. in s-coord. of alpha */
    double **omega_dr; /* r-der. in s-coord. of omega */

    double **nu_dth;   /* theta-der. in mu-coord. of nu */
    double **B_dth;    /* theta-der. in mu-coord. of B */
    double **alpha_dth;/* theta-der. in mu-coord. of alpha */
    double **omega_dth;/* theta-der. in mu-coord. of omega */

  } ini_data;

  // set default parameters
  void RNS_params_set_default();

  // set based on input file
  void RNS_params_set_inputfile(char *inputfile);
  void RNS_params_set_Real(char *key, double value);
  void RNS_params_set_Int(char *key, int value);
  void RNS_params_set_Boolean(char *key, bool value);

  ini_data * RNS_make_initial_data();

  void RNS_Cartesian_interpolation
  (ini_data *data,    // struct containing the previously calculated solution
   int *imin,         // min, max idxs of Cartesian Grid in the three directions
   int *imax,         // in the three dirs
   int *nxyz,
   double *x,         // Cartesian coordinates
   double *y,
   double *z,
   double *lapse,             // ADM metric
   double *betax, double *betay, double *betaz,
   double *adm_gxx, double *adm_gxy, double *adm_gxz,
   double *adm_gyy, double *adm_gyz, double *adm_gzz,
   double *adm_Kxx, double *adm_Kxy, double *adm_Kxz,
   double *adm_Kyy, double *adm_Kyz, double *adm_Kzz,
   double *grhd_rho,          // hydro
   double *grhd_epsl,
   double *grhd_vx, double *grhd_vy, double *grhd_vz,
   double *grhd_ux, double *grhd_uy, double *grhd_uz,
   double *grhd_p);

  // for cleanup purposes
  void RNS_finalise(ini_data *data);

  // output
  void RNS_data_print(ini_data *data);
  void RNS_data_tofile(ini_data *data);

#ifdef __cplusplus
}
#endif

#endif /* RNS_HEADER_H */
