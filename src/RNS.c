/* RNS.c */

#include "RNS.h"

/* -------------------------------------------------------------------------*/
/* PARAMETERS                                                               */
/* -------------------------------------------------------------------------*/

/* Manage memory internally */
static int allocated_params = 0;

static void _alloc_params_mem_if_req(){
  if(!allocated_params){
    params_alloc();
    allocated_params = 1;
  }
}

static void _dealloc_params_mem_if_req(){
  if(allocated_params){
    params_free();
    allocated_params = 0;
  }
}

/* Control interaction with internal parameters in the following */ //

void RNS_params_set_Real(char *key, double value){
  /*
    Set parameters according to input.
  */
  params_set_real(key, value);       // replace, don't append new
}

void RNS_params_set_Int(char *key, int value){
  /*
    Set parameters according to input.
  */
  params_set_int(key, value);   // replace, don't append new
}

void RNS_params_set_Boolean(char *key, bool value){
  /*
    Set parameters according to input.
  */
  params_set_bool(key, value);   // replace, don't append new
}

void RNS_params_set_String(char *key, char * value){
  /*
    Set parameters according to input.
  */
  params_set_str(key, value);   // replace, don't append new
}

void RNS_params_set_inputfile(char* inputfile){
  /*
    Set parameters from input file
  */

  // allocate and seed with defaults such that we have a guaranteee of params
  _alloc_params_mem_if_req();
  RNS_params_set_default();

  if(inputfile!=NULL){
    params_read(inputfile);

    // Set outputdir from parfile
    //SB: this should be improved.
    /*
    char od[STRLEN];
    int n = strlen(inputfile);
    if (STREQL(inputfile+(n-4),".par")) {
      strncpy (od,inputfile,n-4);
      od[n-4] = '\0';
      //printf("%s\n",s);
      params_set_str("outputdir",od);
    }
    */
  }
}

void RNS_params_set_default(){
  /*
    Set default parameters.
  */
  _alloc_params_mem_if_req();

  /* Add parameters here */
  params_add_int("rns_verbose",0); // verb. levels = 0,1,2 see RNS.h
  params_add_int("rns_SDIV",601);  // number of points in radial dir
  params_add_int("rns_MDIV",301);  // number of points in angular dir
  params_add_int("rns_lmax",10);
  params_add_int("rns_save_2ddata",0);

  params_add_real("rns_SMAX",0.9999);
  params_add_real("rns_accuracy",1e-10);
  params_add_real("rns_conv_factor",1.);
  params_add_real("rns_axes_ratio",1.0);
  params_add_real("rns_rhoc",1.28e-3); // central density
  params_add_real("rns_A_diff",1.0);   // differential rotation's A parameters
  params_add_real("rns_eos_K",100.);   // polytropic EOS K
  params_add_real("rns_eos_Gamma",2.); // polytropic EOS Gamma
  params_add_real("rns_star_radius_e",0.); // internally set!
  params_add_real("rns_star_radius_p",0.); // internally set!
  params_add_real("atm_level_rho",1e-10); // rho level of atmosphere
  params_add_real("atm_level_e",0.); // internally reset!
  params_add_real("atm_level_p",0.); // internally reset!

  params_add_str("rns_eos_type","poly"); // {"poly","tab"}
  params_add_str("rns_eos_tab",""); // EOS tab filename
  params_add_str("rns_rotation_type","uniform"); // {"uniform","diff"}
  params_add_str("outputdir","./");

}

/* -------------------------------------------------------------------------*/
/* MAIN ROUTINE TO BUILD INITIAL DATA                                       */
/* -------------------------------------------------------------------------*/

/* the folloing array contains treshold for the iteration of the initial guess
   can be changed by hand if needed */
#define ARSTEPS (11)
double axes_ratio_steps[ARSTEPS] = {.8, .6, .5, .4, .3, .2, .1, .05, .0025, .00125, .00625};

ini_data* RNS_make_initial_data() {

  char rotation_type[STRLEN];
  strcpy(rotation_type,params_get_str("rns_rotation_type"));
  const double accuracy = params_get_real("rns_accuracy");
  const double A_diff   = params_get_real("rns_A_diff");
  const double cf       = params_get_real("rns_conv_factor");
  const int RNS_lmax    = params_get_int("rns_lmax");
  const int SDIV        = params_get_int("rns_SDIV");
  const int MDIV        = params_get_int("rns_MDIV");
  char eos_type[STRLEN];
  strcpy(eos_type,params_get_str("rns_eos_type"));
  if( (strcmp(eos_type,"poly")!=0) && (strcmp(eos_type,"tab")!=0) ){
    ERROR("unknown eos for rns\n");
  }

  const double axes_ratio = params_get_real("rns_axes_ratio");
  double rho0_center      = params_get_real("rns_rhoc");

  const double atm_factor = 1e-4;
  const double rho0_atm   =  atm_factor * params_get_real("atm_level_rho");

  const double eos_ideal_fluid_gamma = params_get_real("rns_eos_Gamma");
  const double eos_k                 = params_get_real("rns_eos_K");

  const int print_dif = params_get_int("rns_verbose");


  /* EQUILIBRIUM VARIABLES */

  int n_tab,                      /* Number of points in EOS file */
    a_check=0,                   /* if =200, iteration diverges */
    n_near,
    n_nearest;

  double log_e_tab[TABPTS],       /* energy dens./c^2 in tab. EOS */
    log_p_tab[TABPTS],            /* pressure in tabulated EOS */
    log_h_tab[TABPTS],            /* enthalpy in EOS file */
    log_n0_tab[TABPTS],           /* number density in EOS file */
    Gamma_tab[TABPTS],            /* Gamma in tab. EOS file */
    p_center,                     /* central pressure */
    h_center,                     /* central enthalpy */
    e_surface,                    /* surface en. density */
    p_surface,                    /* surface pressure */
    enthalpy_min,                 /* minimum enthalpy in EOS */
    n0;

  double *s_gp,                   /* s grid points */
    *mu,                          /* \mu grid points */
    **rho_potential,              /* potential \rho_potential */
    **gama,                       /* potential \gamma */
    **omega,                      /* potential \omega */
    **alpha,                      /* potential \alpha */
    **energy,                     /* energy density \epsilon */
    **pressure,                   /* pressure */
    **enthalpy,                   /* enthalpy */
    **velocity_sq,                /* square of velocity */
    **Omega_diff;                 /* Diff. ang. vel. */

  double R_e,                     /* Circumferential radius */
    Mass,                         /* Gravitational mass */
    Mass_0,                       /* Baryon Mass */
    T,                            /* Rotational kinetic energy */
    W,                            /* Gravitational binding energy */
    Omega,                        /* Angular velocity */
    Omega_K,                      /* Ang. vel. of part. in orbit at eq.*/
    r_e,                          /* coord. radius at equator */
    Omega_e,                      /* Ang. vel. at equator, when difrot. */
    J;                            /* Angular momentun */

  if (print_dif>print_globalprop1) {
    printf("Set rotating star initial data\n");
    printf("Parameters:\n");
    printf("  rotation  = %s\n", rotation_type);
    printf("  axesratio = %e\n", axes_ratio);
    printf("  rhoc      = %e\n", rho0_center);
    printf("  eos_type  = %s\n",eos_type);
    printf("  eos_tab   = %s\n",params_get_str("rns_eos_tab"));
    printf("  eos_gamma = %e\n",eos_ideal_fluid_gamma);
    printf("  eos_K     = %e\n",eos_k);
    printf("Solving rns ... \n");
  }

  /* COMPUTE POLYTROPIC INDEX AND CENTRAL ENERGY DENSITY */
  double n_P = 1.0/(eos_ideal_fluid_gamma-1.0);
  double e_center = (eos_k*pow(rho0_center,eos_ideal_fluid_gamma)/(eos_ideal_fluid_gamma-1.0)+rho0_center);

  /* these can be also 0, atm set later on */
  double e_atm = rho0_atm;
  double p_atm = eos_k*pow(rho0_atm,eos_ideal_fluid_gamma);

  /* TABULATED EOS OPTION */
  if(strcmp(eos_type,"tab")==0) {

    load_eos(params_get_str("rns_eos_tab"), log_e_tab, log_p_tab, log_h_tab, log_n0_tab, &n_tab );
    n_nearest = 50;
    n0 = rho0_center/(MB*cactusM);
    e_center = pow(10.0,interp(log_n0_tab, log_e_tab, n_tab,log10(n0), &n_nearest));

    e_atm = 0.;//e_at_n0(rho_0_atm, log_e_tab, log_n0_tab, n_tab, &n_near);
    p_atm = 0.;//p_at_e( e_atm, log_p_tab, log_e_tab, n_tab, &n_near);

  }

  params_set_real("atm_level_e",e_atm);
  params_set_real("atm_level_p", p_atm);
  /* SET UP GRID */
  s_gp = malloc((SDIV+1)*sizeof(double));
  mu = malloc((MDIV+1)*sizeof(double));
  rns_make_grid(s_gp, mu);

  /* ALLLOCATE MEMORY */
  rho_potential = dmatrix(1,SDIV,1,MDIV);
  gama = dmatrix(1,SDIV,1,MDIV);
  alpha = dmatrix(1,SDIV,1,MDIV);
  omega = dmatrix(1,SDIV,1,MDIV);
  energy = dmatrix(1,SDIV,1,MDIV);
  pressure = dmatrix(1,SDIV,1,MDIV);
  enthalpy = dmatrix(1,SDIV,1,MDIV);
  velocity_sq = dmatrix(1,SDIV,1,MDIV);
  Omega_diff = dmatrix(1,SDIV,1,MDIV);

  /* INITIALIZE VARIABLES WITH ZERO */
  for(int s=1;s<=SDIV;s++)
    for(int m=1;m<=MDIV;m++) {
      rho_potential[s][m] = 0.0e0;
      gama[s][m] = 0.0e0;
      alpha[s][m] = 0.0e0;
      omega[s][m] = 0.0e0;
      energy[s][m] = 0.0e0;
      pressure[s][m] = 0.0e0;
      enthalpy[s][m] = 0.0e0;
      velocity_sq[s][m] = 0.0e0;
      Omega_diff[s][m] = 0.0e0;
    }

  /* SET DEFAULT EQUILIBRIUM PARAMETERS */
  e_surface=0.0;
  p_surface=0.0;
  enthalpy_min=0.0;
  Omega_e=0.0; /* initialize ang. vel. at equator for diff. rot. */

  // if(strcmp(eos_type,"tab")==0) {
  //   e_surface=7.8*C_SPEED*C_SPEED*KSCALE;
  //   p_surface=1.01e8*KSCALE;
  //   enthalpy_min=1.0/(C_SPEED*C_SPEED);
  //
  //   /* MAKE e_center DIMENSIONLESS FOR TAB. EOS */
  //   e_center *= (C_SPEED*C_SPEED*KSCALE);
  //
  // }

  /* COMPUTE DIMENSIONLESS CENTRAL PRESSURE AND ENTHALPY */
  make_center( e_center, log_e_tab, log_p_tab, log_h_tab, n_tab,
	       eos_type, eos_k,eos_ideal_fluid_gamma, &p_center, &h_center);
  rho0_center = (e_center+p_center)*exp(-h_center);
  /* COMPUTE THE MODEL */

  /* COMPUTE A SPHERICAL STAR AS A FIRST GUESS */

  guess( s_gp, eos_type, eos_k,e_center, p_center, p_surface, e_surface,
	 eos_ideal_fluid_gamma, log_e_tab, log_p_tab, log_h_tab, n_tab, rho_potential, gama,
	 alpha, omega, &r_e );
  /* ITERATE AND SPIN UP THE GUESS */
  for(int i=0;i<ARSTEPS;i++) {

    if (axes_ratio>axes_ratio_steps[i]) break;
    if(print_dif>print_globalprop1)
      printf("iter %d (%d) a=%g (%g)\n",i,ARSTEPS,axes_ratio_steps[i],axes_ratio);

    iterate( s_gp, mu, eos_type, eos_k, log_e_tab, log_p_tab, log_h_tab,
	     n_tab, eos_ideal_fluid_gamma, axes_ratio_steps[i],
	     h_center, enthalpy_min, &a_check,
	     accuracy,print_dif,cf, &r_e, rho_potential, gama, alpha, omega,
	     energy, pressure, enthalpy, velocity_sq, &Omega,
	     rotation_type,A_diff,&Omega_e, Omega_diff,RNS_lmax);


  }

  /* last one */
  iterate( s_gp, mu, eos_type, eos_k, log_e_tab, log_p_tab, log_h_tab,
	   n_tab, eos_ideal_fluid_gamma, axes_ratio,
	   h_center, enthalpy_min, &a_check,
	   accuracy,print_dif,cf, &r_e, rho_potential, gama, alpha, omega,
	   energy, pressure, enthalpy, velocity_sq, &Omega,
	   rotation_type,A_diff,&Omega_e, Omega_diff,RNS_lmax);

  if (a_check) ERROR("rns did not converged, nans during iterations");

  /* COMPUTE EQUILIBRIUM QUANTITIES (Mass, Radius, T/W etc.) */
  comp_values( s_gp, mu, axes_ratio, e_surface, r_e, eos_type, log_e_tab,
	       log_n0_tab, n_tab, Omega, rho_potential, gama, alpha, omega,
	       energy, pressure, enthalpy, velocity_sq, &Mass,
	       &Mass_0, &T, &W, &Omega_K, &R_e, rotation_type,Omega_diff,
	       &J);

  params_set_real("rns_star_radius_e", R_e);
  params_set_real("rns_star_radius_p", R_e*axes_ratio);

  /* PRINT-OUT SOME EQUILIBRIUM QUANTITIES */
  if(print_dif){
    print_global_quantities(eos_type,
			    rotation_type,
			    rho0_center,
			    axes_ratio,
			    e_center,
			    Mass,
			    Mass_0,
			    R_e,
			    J,
			    T,
			     W,
			    Omega,
			    Omega_K,
			    Omega_e);
  }

  /* PREPARE OUTPUT: CONSTRUCT ARRAYS WITH NEEDED POLAR QUANTITIES */
  ini_data *data;
  data = (ini_data *) calloc(1, sizeof(ini_data));
  if (data == NULL) ERROR("Out of memory");

  data->MDIV = MDIV;
  data->SDIV = SDIV;
  data->r_e = r_e;
  data->mass = Mass;
  data->mass_0 = Mass_0;

  data->s_gp = malloc((SDIV+1)*sizeof(double));
  data->mu = malloc((MDIV+1)*sizeof(double));
  rns_make_grid(data->s_gp, data->mu);

  data->rho_0 = dmatrix(1,SDIV,1,MDIV);
  data->energy = dmatrix(1,SDIV,1,MDIV);
  data->pressure = dmatrix(1,SDIV,1,MDIV);
  data->Omega_diff = dmatrix(1,SDIV,1,MDIV);

  data->nu = dmatrix(1,SDIV,1,MDIV);
  data->B = dmatrix(1,SDIV,1,MDIV);
  data->alpha = dmatrix(1,SDIV,1,MDIV);
  data->omega = dmatrix(1,SDIV,1,MDIV);

  data->nu_dr = dmatrix(1,SDIV,1,MDIV);
  data->B_dr = dmatrix(1,SDIV,1,MDIV);
  data->alpha_dr = dmatrix(1,SDIV,1,MDIV);
  data->omega_dr = dmatrix(1,SDIV,1,MDIV);

  data->nu_dth = dmatrix(1,SDIV,1,MDIV);
  data->B_dth = dmatrix(1,SDIV,1,MDIV);
  data->alpha_dth = dmatrix(1,SDIV,1,MDIV);
  data->omega_dth = dmatrix(1,SDIV,1,MDIV);

  for(int m=1;m<=MDIV;m++) {
    for(int s=1;s<=SDIV;s++) {
      data->rho_0[s][m] = (energy[s][m]+pressure[s][m])*exp(-enthalpy[s][m]);
      data->energy[s][m] = energy[s][m];
      data->pressure[s][m] = pressure[s][m];
      data->Omega_diff[s][m] = Omega_diff[s][m];
    }
  }

  for(int m=1;m<=MDIV;m++) {
    for(int s=1;s<=SDIV;s++) {
      data->nu[s][m] = (gama[s][m]+rho_potential[s][m])/2.0;
      data->B[s][m] = exp(gama[s][m]);
      data->alpha[s][m] = alpha[s][m];
      data->omega[s][m] = omega[s][m];
    }
  }

  for(int m=1;m<=MDIV;m++)
    for(int s=1;s<=SDIV;s++) {
      data->nu_dr[s][m] = deriv_s(data->nu,s,m)*SQ(1.0-s_gp[s])/r_e;
      data->B_dr[s][m] = deriv_s(data->B,s,m)*SQ(1.0-s_gp[s])/r_e;
      data->alpha_dr[s][m] = deriv_s(data->alpha,s,m)*SQ(1.0-s_gp[s])/r_e;
      data->omega_dr[s][m] = deriv_s(data->omega,s,m)*SQ(1.0-s_gp[s])/r_e;
      data->nu_dth[s][m] = deriv_m(data->nu,s,m)*(-sqrt(1.0-SQ(mu[m])));
      data->B_dth[s][m] = deriv_m(data->B,s,m)*(-sqrt(1.0-SQ(mu[m])));
      data->alpha_dth[s][m] = deriv_m(data->alpha,s,m)*(-sqrt(1.0-SQ(mu[m])));
      data->omega_dth[s][m] = deriv_m(data->omega,s,m)*(-sqrt(1.0-SQ(mu[m])));
    }

  /* FREE MEMORY */
  free(s_gp);
  free(mu);
  free_dmatrix(energy,1,SDIV,1,MDIV);
  free_dmatrix(pressure,1,SDIV,1,MDIV);
  free_dmatrix(Omega_diff,1,SDIV,1,MDIV);
  free_dmatrix(enthalpy,1,SDIV,1,MDIV);
  free_dmatrix(velocity_sq,1,SDIV,1,MDIV);

  free_dmatrix(rho_potential,1,SDIV,1,MDIV);
  free_dmatrix(gama,1,SDIV,1,MDIV);
  free_dmatrix(alpha,1,SDIV,1,MDIV);
  free_dmatrix(omega,1,SDIV,1,MDIV);

  /* SAVE 2D FILE (todo) */
  if (params_get_int("rns_save_2ddata")) {
    make_output_dir();
    RNS_data_tofile(data);
    if(print_dif) {
      printf("Data saved in: %s\n",params_get_str("outputdir"));
    }
  }

  return data;

}

void RNS_finalise(ini_data *data){
  /*
    Clean up all internally allocated objects.
  */

  const int SDIV = params_get_int("rns_SDIV");
  const int MDIV = params_get_int("rns_MDIV");

  if (data) {

    if(data->s_gp) free(data->s_gp);
    if(data->mu) free(data->mu);

    if (data->rho_0) free_dmatrix(data->rho_0,1,SDIV,1,MDIV);
    if (data->energy) free_dmatrix(data->energy,1,SDIV,1,MDIV);
    if (data->pressure) free_dmatrix(data->pressure,1,SDIV,1,MDIV);
    if (data->Omega_diff) free_dmatrix(data->Omega_diff,1,SDIV,1,MDIV);

    if (data->nu) free_dmatrix(data->nu,1,SDIV,1,MDIV);
    if (data->B) free_dmatrix(data->B,1,SDIV,1,MDIV);
    if (data->alpha) free_dmatrix(data->alpha,1,SDIV,1,MDIV);
    if (data->omega) free_dmatrix(data->omega,1,SDIV,1,MDIV);

    if (data->nu_dr) free_dmatrix(data->nu_dr,1,SDIV,1,MDIV);
    if (data->B_dr) free_dmatrix(data->B_dr,1,SDIV,1,MDIV);
    if (data->alpha_dr) free_dmatrix(data->alpha_dr,1,SDIV,1,MDIV);
    if (data->omega_dr) free_dmatrix(data->omega_dr,1,SDIV,1,MDIV);

    if (data->nu_dth) free_dmatrix(data->nu_dth,1,SDIV,1,MDIV);
    if (data->B_dth) free_dmatrix(data->B_dth,1,SDIV,1,MDIV);
    if (data->alpha_dth) free_dmatrix(data->alpha_dth,1,SDIV,1,MDIV);
    if (data->omega_dth) free_dmatrix(data->omega_dth,1,SDIV,1,MDIV);

    free(data);
  }

  _dealloc_params_mem_if_req();
}


/* -------------------------------------------------------------------------*/
/* CARTESIAN GRID INTERPOLATION                                             */
/* -------------------------------------------------------------------------*/

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
 double *grhd_vx, double *grhd_vy, double *grhd_vz, // v^i
 double *grhd_ux, double *grhd_uy, double *grhd_uz, // u^i = W v^i
 double *grhd_p
 ) {

  const int verbose = params_get_int("rns_verbose");
  if (verbose>print_globalprop1){
    printf("Interpolating on Cartesian grid...\n");
  }

  if (data == NULL) ERROR("No data to interpolate");

  const double rho0_atm = params_get_real("atm_level_rho");
  const double e_atm = params_get_real("atm_level_e");
  const double p_atm = params_get_real("atm_level_p");

  double nu_ijk, B_ijk, alpha_ijk, omega_ijk;
  double nu_dr_ijk, B_dr_ijk, alpha_dr_ijk, omega_dr_ijk;
  double nu_dtheta_ijk, B_dtheta_ijk, alpha_dtheta_ijk, omega_dtheta_ijk;
  double rho_0_ijk, energy_ijk, pressure_ijk, Omega_ijk;

  const int xxx = nxyz[0];
  const int xxxyyy = nxyz[0]*nxyz[1];

  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
      for (int i = imin[0]; i < imax[0]; ++i) {

	const int ijk = i + xxx * j + xxxyyy* k;

	double x_i = x[i];
        double y_j = y[j];
        double z_k = z[k];

	grid_interp_all_new (data->s_gp, data->mu, data->r_e,
			     x_i, y_j,z_k,
			     data->nu, data->B, data->alpha, data->omega,
			     data->nu_dr, data->B_dr, data->alpha_dr, data->omega_dr,
			     data->nu_dth, data->B_dth, data->alpha_dth, data->omega_dth,
			     data->rho_0, data->energy, data->pressure, data->Omega_diff,
			     &nu_ijk, &B_ijk, &alpha_ijk, &omega_ijk,
			     &nu_dr_ijk, &B_dr_ijk, &alpha_dr_ijk, &omega_dr_ijk,
			     &nu_dtheta_ijk, &B_dtheta_ijk, &alpha_dtheta_ijk, &omega_dtheta_ijk,
			     &rho_0_ijk, &energy_ijk, &pressure_ijk,
			     &Omega_ijk);

	double exp_nu_ijk = exp(nu_ijk);
	double exp_alpha_ijk = exp(alpha_ijk);

	double r_ijk = sqrt(SQ(x_i)+SQ(y_j)+SQ(z_k));
	double r_bar_ijk = sqrt(SQ(x_i)+SQ(y_j));

	double dr_dx = x_i / r_ijk;
	double dr_dy = y_j / r_ijk;
	double dr_dz = z_k / r_ijk;

	double dtheta_dx = x_i*z_k/(SQ(r_ijk)*r_bar_ijk);
	double dtheta_dy = y_j*z_k/(SQ(r_ijk)*r_bar_ijk);
	double dtheta_dz = - r_bar_ijk/SQ(r_ijk);

	double nu_dx = dr_dx*nu_dr_ijk + dtheta_dx*nu_dtheta_ijk;
	double nu_dy = dr_dy*nu_dr_ijk + dtheta_dy*nu_dtheta_ijk;

	double B_dx = dr_dx*B_dr_ijk + dtheta_dx*B_dtheta_ijk;
	double B_dy = dr_dy*B_dr_ijk + dtheta_dy*B_dtheta_ijk;

	double alpha_dx = dr_dx*alpha_dr_ijk + dtheta_dx*alpha_dtheta_ijk;
	double alpha_dy = dr_dy*alpha_dr_ijk + dtheta_dy*alpha_dtheta_ijk;

	double omega_dx = dr_dx*omega_dr_ijk + dtheta_dx*omega_dtheta_ijk;
	double omega_dy = dr_dy*omega_dr_ijk + dtheta_dy*omega_dtheta_ijk;

	/* enforce omega_dz=0 at z=0 (it is slightly nonzero due
	   to O(h) forwards formula in computing derivative) */
	double omega_dz;
	if(z_k==0.0)
	  omega_dz = 0.0;
	else
	  omega_dz = dr_dz*omega_dr_ijk + dtheta_dz*omega_dtheta_ijk;


	double gxx = ( SQ(B_ijk*y_j/exp_nu_ijk)
		       +SQ(exp_alpha_ijk*x_i) ) /
	  (SQ(x_i)+SQ(y_j));

	double gxy = ( SQ(exp_alpha_ijk)
		       -SQ(B_ijk/exp_nu_ijk) ) *
	  x_i*y_j/(SQ(x_i)+SQ(y_j));

	double gxz = 0.0;

	double gyy = ( SQ(B_ijk*x_i/exp_nu_ijk)
		       +SQ(exp_alpha_ijk*y_j) ) /
	  (SQ(x_i)+SQ(y_j));

	double gyz = 0.0;

	double gzz = SQ(exp_alpha_ijk);


	double kxx = (  ( SQ(r_bar_ijk)*y_j*omega_dx+
			  (x_i*nu_dy-y_j*nu_dx)*SQ(y_j)
			  *omega_ijk)*SQ(B_ijk)
			+(y_j*B_dx-x_i*B_dy)*omega_ijk
			*SQ(y_j)*B_ijk
			+(y_j*alpha_dx-x_i*alpha_dy)
			*omega_ijk*SQ(x_i*exp_alpha_ijk
				      *exp_nu_ijk))/(SQ(r_bar_ijk
							*exp_nu_ijk)*exp_nu_ijk);

	double kxy = ( ( 0.5*SQ(r_bar_ijk)*
			 (y_j*omega_dy - x_i*omega_dx) +
			 (y_j*nu_dx-x_i*nu_dy)*x_i*y_j*
			 omega_ijk )*SQ(B_ijk)
		       +(-y_j*B_dx+x_i*B_dy)*omega_ijk
		       *x_i*y_j*B_ijk
		       +(y_j*alpha_dx-x_i*alpha_dy)
		       *omega_ijk*x_i*y_j
		       *SQ(exp_alpha_ijk*exp_nu_ijk))/
	  (SQ(r_bar_ijk*exp_nu_ijk)
	   *exp_nu_ijk);

	double kxz = 0.5*SQ(B_ijk)*y_j*omega_dz/
	  ( SQ(exp_nu_ijk)*exp_nu_ijk );

	double kyy = ( ( -SQ(r_bar_ijk)*x_i*omega_dy+
			 (x_i*nu_dy-y_j*nu_dx)*SQ(x_i)*
			 omega_ijk )*SQ(B_ijk)
		       +(y_j*B_dx-x_i*B_dy)*omega_ijk
		       *SQ(x_i)*B_ijk
		       +(y_j*alpha_dx-x_i*alpha_dy)
		       *omega_ijk*SQ(y_j*exp_alpha_ijk
				     *exp_nu_ijk))/(SQ(r_bar_ijk
						       *exp_nu_ijk)*exp_nu_ijk);

	double kyz = -0.5*SQ(B_ijk)*x_i*omega_dz/
	  ( SQ(exp_nu_ijk)*exp_nu_ijk );

	double kzz = (y_j*alpha_dx-x_i*alpha_dy)*
	  omega_ijk*SQ(exp_alpha_ijk)/
	  exp_nu_ijk;


	if(x_i==0.0 && y_j==0.0) {

	  gxx = SQ(exp_alpha_ijk);
	  gyy = SQ(exp_alpha_ijk);
	  gzz = SQ(exp_alpha_ijk);

	  gxy = 0.0;
	  gxz = 0.0;
	  gyz = 0.0;

	  kxx = 0.0;
	  kyy = 0.0;
	  kzz = 0.0;
	  kxy = 0.0;
	  kxz = 0.0;
	  kyz = 0.0;

	}

	// populate Cartesian tensors
	if (lapse) lapse[ijk] = exp_nu_ijk;
      	if (betax) betax[ijk] = omega_ijk*y_j;
	if (betay) betay[ijk] = -omega_ijk*x_i;
	if (betaz) betaz[ijk] = 0.0;

	if (adm_gxx) adm_gxx[ijk] = gxx;
	if (adm_gxy) adm_gxy[ijk] = gxy;
	if (adm_gxz) adm_gxz[ijk] = gxz;
	if (adm_gyy) adm_gyy[ijk] = gyy;
	if (adm_gyz) adm_gyz[ijk] = gyz;
	if (adm_gzz) adm_gzz[ijk] = gzz;

	if (adm_Kxx) adm_Kxx[ijk] = kxx;
	if (adm_Kxy) adm_Kxy[ijk] = kxy;
	if (adm_Kxz) adm_Kxz[ijk] = kxz;
	if (adm_Kyy) adm_Kyy[ijk] = kyy;
	if (adm_Kyz) adm_Kyz[ijk] = kyz;
	if (adm_Kzz) adm_Kzz[ijk] = kzz;

	if( (rho_0_ijk<=0.0) || (energy_ijk<=0.0) || (pressure_ijk<=0.0) ) {
	  rho_0_ijk    = rho0_atm;
	  energy_ijk   = e_atm;
	  pressure_ijk = p_atm;
	}

	double vx = (omega_ijk-Omega_ijk)*y_j/exp_nu_ijk;
	double vy = -(omega_ijk-Omega_ijk)*x_i/exp_nu_ijk;
	double vz = 0.;
	double vlx = 0, vly = 0, vlz = 0, v2 = 0;
	contract_vect_v2(vx, vy, vz,
			 gxx, gxy, gxz,
			 gyy, gyz, gzz,
			 &vlx, &vly, &vlz,
			 &v2);

	if (fabs(v2) < 1e-20) v2 = 0.;
	double W = 1.0/sqrt(1-v2);

	if (grhd_rho) grhd_rho[ijk]  = rho_0_ijk;
	if (grhd_epsl) grhd_epsl[ijk] = energy_ijk/rho_0_ijk-1.0;
	if (grhd_p) grhd_p[ijk]  = pressure_ijk;

	if (grhd_vx) grhd_vx[ijk] = vx;
	if (grhd_vy) grhd_vy[ijk] = vy;
	if (grhd_vz) grhd_vz[ijk] = vz;

	if (grhd_ux) grhd_ux[ijk] = W * vx;
	if (grhd_uy) grhd_uy[ijk] = W * vy;
	if (grhd_uz) grhd_uz[ijk] = W * vz;

      }
    }
  }

  if (verbose>print_globalprop1) printf("Done\n");

  return;

}

/* -------------------------------------------------------------------------*/
/* OUTPUT STUFF                                                             */
/* -------------------------------------------------------------------------*/

void RNS_data_print(ini_data *data)
{
  const int SDIV = params_get_int("rns_SDIV");
  const int MDIV = params_get_int("rns_MDIV");

  printf("# MDIV %d\n# SDIV %d\n", data->MDIV, data->SDIV);
  printf("# rho_center = %e\n# axes_ratio = %e\n",
	  params_get_real("rns_rhoc"),
	  params_get_real("rns_axes_ratio"));
  printf("# rotation = %s\n",params_get_str("rns_rotation_type"));
  printf("# EOS = %s\n",params_get_str("rns_eos_type"));
  if (strcmp(params_get_str("rns_eos_type"),"poly")==0) {
    printf("# EOS_gamma = %e EOS_K = %e\n",
	    params_get_real("rns_eos_Gamma"),
	    params_get_real("rns_eos_K"));
  } else {
    printf("# EOS_file = %s\n",params_get_str("rns_eos_tab"));
  }
  printf("# \n");
  printf("#data: s_gp:1 mu:2 nu:3 nu_dr:4 nu_dth:5 B:6 B_dr:7 B_dth:8 omega:9 omega_dr:10 omega_dth:11 alpha:12 alpha_dr:13 alpha_dth:14 rho0:15 energy:16 pressure:17 Omega_diff:18\n");
  for(int m=1;m<=MDIV;m++) {
    for(int s=1;s<=SDIV;s++) {
      printf("%le %le"
	     "%le %le %le %le %le %le "
	     "%le %le %le %le %le %le "
	     "%le %le %le %le \n",
	      data->s_gp[s], data->mu[m],
	     data->nu[s][m], data->nu_dr[s][m], data->nu_dth[s][m],
	     data->B[s][m], data->B_dr[s][m],data->B_dth[s][m],
	     data->omega[s][m],data->omega_dr[s][m],data->omega_dth[s][m],
	     data->alpha[s][m],data->alpha_dr[s][m],data->alpha_dth[s][m],
	     data->rho_0[s][m],
	     data->energy[s][m],
	     data->pressure[s][m],
	     data->Omega_diff[s][m]);
    }
  }

}

void RNS_data_tofile(ini_data *data)
{
  const int SDIV = params_get_int("rns_SDIV");
  const int MDIV = params_get_int("rns_MDIV");

  FILE *fp;
  char file[STRLEN*10];
  sprintf(file,"%s/RNSdata_2d.dat",params_get_str("outputdir"));

  fp = fopen(file, "w");
  fprintf(fp,"# MDIV %d\n# SDIV %d\n", data->MDIV, data->SDIV);
  fprintf(fp,"# rho_center = %e\n# axes_ratio = %e\n",
	  params_get_real("rns_rhoc"),
	  params_get_real("rns_axes_ratio"));
  fprintf(fp,"# rotation = %s\n",params_get_str("rns_rotation_type"));
  fprintf(fp,"# EOS = %s\n",params_get_str("rns_eos_type"));
  if (strcmp(params_get_str("rns_eos_type"),"poly")==0) {
    fprintf(fp,"# EOS_gamma = %e EOS_K = %e\n",
	    params_get_real("rns_eos_Gamma"),
	    params_get_real("rns_eos_K"));
  } else {
    fprintf(fp,"# EOS_file = %s\n",params_get_str("rns_eos_tab"));
  }
  fprintf(fp,"# r_e = %e\n",data->r_e);
  fprintf(fp,"# mass = %e\n",data->mass);
  fprintf(fp,"# mass_0 = %e\n",data->mass_0);
  fprintf(fp,"# \n");
  fprintf(fp,"#data: s_gp:1 mu:2 nu:3 nu_dr:4 nu_dth:5 B:6 B_dr:7 B_dth:8 omega:9 omega_dr:10 omega_dth:11 alpha:12 alpha_dr:13 alpha_dth:14 rho0:15 energy:16 pressure:17 Omega_diff:18\n");
  for(int m=1;m<=MDIV;m++) {
    for(int s=1;s<=SDIV;s++) {
      fprintf(fp,"%le %le "
	     "%le  %le  %le %le %le %le "
	     "%le  %le  %le %le %le %le "
	     "%le  %le  %le %le \n",
	      data->s_gp[s], data->mu[m],
	     data->nu[s][m], data->nu_dr[s][m], data->nu_dth[s][m],
	     data->B[s][m], data->B_dr[s][m],data->B_dth[s][m],
	     data->omega[s][m],data->omega_dr[s][m],data->omega_dth[s][m],
	     data->alpha[s][m],data->alpha_dr[s][m],data->alpha_dth[s][m],
	     data->rho_0[s][m],
	     data->energy[s][m],
	     data->pressure[s][m],
	     data->Omega_diff[s][m]);
    }
  }
  fclose(fp);

}
