/************************************************************************** 
 *                            ROTSTAR.C                                    * 
 *                                                                         *
 * This program uses the routines of RNS to compute equilibrium models of  *    
 * uniformly rotating neutron stars described either by polytropic or      *
 * realistic EOS (tables)                                                  *
 * The inputs are the the central energy density and the ratio between the *
 * polar the and equatorial radius                                         *
 * For polytropes also the polytropic index and the constant must be given.*
 *                                                                         *
 * The RNS routines use (always? some kind of) polytropic units:           *
 *                            c=G=K=1 , see:                               *
 *       Cook, Shapiro & Teukolsky ApJ 398:203-223 (1992)                  *
 * Hystorically realistic EOS models computation required input in CGS     *
 * units, while polytropes models in polytropic units.                     *
 * Ouput was consistent                                                    *
 *                                                                         *
 * Since it is generically better to use the standard dimensionless units: *
 *                            c=G=Msun=1 ,                                 *
 * in this program I/O for polytropes is in this standard units made with  *
 * conversions inside.                                                     *
 * I/O for realistic EOS still follows the original and will be fixed soon *
 *                                                                         *
 * Quick help to run it                                                    *

 ./rotstar -h

 * Main quantities are output at screen                                    *
 * Relevant 2D arrays can be output in a ASCII file                        *
 * S.Bernuzzi                                                              *
 **************************************************************************/

#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "consts.h"
#include "nrutil.h"
#include "equil.h"

/***************************************************************************
 * printing routine for realistic star
 ***************************************************************************/

void printtab(
    double r_ratio,
    double e_center, double Mass, double Mass_0, double R_e,
    double Omega, double Omega_K, double J
    )
{
    
    double I_45;
    
    if( Omega == 0.0) I_45 = 0.0;
    else I_45 = J/(Omega*1.0e45);

    //printf("EQUIL :: ratio\te_15\tM\tM_0\tr_star\tspin\tOmega_K\tI\tJ/M^2\n");
    //printf("EQUIL :: \t[g/cm^3] [sun]\t[sun]\t[km]\t[s-1]\t[s-1]\t[g cm^2]\t\nEQUIL :: ");  
    printf("# ratio\te_15\tM\tM_0\tr_star\tspin\tOmega_K\tI\tJ/M^2\n");
    printf("#\t[g/cm^3] [sun]\t[sun]\t[km]\t[s-1]\t[s-1]\t[g cm^2]\t\n");  
    
    printf(
	/* "%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\n", */
	"%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\n", 
	/* "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", */
	r_ratio,
	e_center,
	Mass/MSUN,
	Mass_0/MSUN,
	R_e/1.0e5,
	Omega,
	Omega_K,
	I_45,
	( C*J/(G*Mass*Mass))
	);
}

/***************************************************************************
 * printing routine for polytropic stars
 ***************************************************************************/

void printpoly(
    double r_ratio,
    double e_center, double Mass, double Mass_0, double R_e,
    double Omega, double Omega_K, double J
    )
{

    double I;
    
    if( Omega == 0.0) I = 0.0;
    else I = J/(Omega);
    
    /* printf("EQUIL :: ratio\te_15\tM\tM_0\tr_star\tspin\tOmega_K\tI\tJ/M^2\nEQUIL :: "); */
    //printf("EQUIL :: ratio\t\te_15\t\tM\t\tM_0\t\tr_star\t\tspin\t\tOmega_K\t\tI\t\tJ/M^2\nEQUIL :: "); 
    printf("# ratio\t\te_15\t\tM\t\tM_0\t\tr_star\t\tspin\t\tOmega_K\t\tI\t\tJ/M^2\n"); 

    printf(
	/* "%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\n", */
	/* "%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\n", */
	"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
	r_ratio,
	e_center,
	Mass,
	Mass_0,
	R_e,
	Omega,
	Omega_K,
	I,
	J/(Mass_0*Mass_0)
	);
}

/**************************************************************************/                                                                                  
/* Transform units from polytropic dimensionless to c=G=M_sun=1           */                                 
/**************************************************************************/

void transform_units(
    char    eos_type[],
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
/*     double **Omega_diff, */ /* Differential angular velocity : NOT USED HERE */
    double *J
    )

{
    
    int i,j; /* counter */
    
    if(strcmp(eos_type,"poly")==0) {    
	
      /* Polytropic units: the length scale is L_poly = K^{N/2} */

	for(i=1;i<=SDIV;i++)
	    for(j=1;j<=MDIV;j++) {
		omega[i][j]    *= ( 1.0/pow(eos_k, n_P/2.0) );
		energy[i][j]   *= ( 1.0/pow(eos_k, n_P) );
		pressure[i][j] *= ( 1.0/pow(eos_k, n_P) );
		/* Omega_diff[i][j] *= ( 1.0/pow(eos_k, n_P/2.0) ); */
	    }
	
	(*rho0_center) *= ( 1.0/pow(eos_k, n_P) );
	(*e_center)    *= ( 1.0/pow(eos_k, n_P) );
	(*p_center)    *= ( 1.0/pow(eos_k, n_P) );
	(*r_e)         *= pow(eos_k, n_P/2.0);
	(*Mass)        *= pow(eos_k, n_P/2.0);
	(*Mass_0)      *= pow(eos_k, n_P/2.0);
 	(*T)           *= pow(eos_k, n_P/2.0); 
 	(*W)           *= pow(eos_k, n_P/2.0); 
	(*Omega)       *= ( 1.0/pow(eos_k, n_P/2.0) );
	(*Omega_K)     *= ( 1.0/pow(eos_k, n_P/2.0) );
	(*R_e)         *= pow(eos_k, n_P/2.0);
	(*Omega_e)     *= ( 1.0/pow(eos_k, n_P/2.0) );
	(*J)           *= pow(eos_k, n_P);

    } else if(strcmp(eos_type,"tab")==0) { 

      if (VERBOSE) {
	printf("WARNING :: Unit Conversion for tab EOS models not yet completely implemented !!!\n");
	printf("WARNING :: Units for equilibrium quantities are given below (mass_radius in equil.c)\n");	     
	printf("WARNING :: Units for the fields should be the following dimensionless units:\n");
	printf("WARNING ::           c=G=KAPPA=1   where:   KAPPA=(1.0e-15 g)xc^2/G=%g [cm^-1]\n",KAPPA);	
      }
      // function mass_radius in equil.c scales the equilibrium quantities
      // here we should scale the fields...
      
      /*
	for(i=1;i<=SDIV;i++)
	for(j=1;j<=MDIV;j++) {
	omega[i][j]    *= KAPPA;
	energy[i][j]   *= 1./KAPPA*KAPPA;
	pressure[i][j] *= 1./KAPPA*KAPPA;
	}
	(*r_e) *= 1.0/KAPPA;
      */
      
    } else {
	
	printf("ERROR :: Wrong value for 'eos_type' in 'transform_units' routine.\n");
	exit(1);
    }
    
}


/*************************************************************************/
/* Main program.                                                         */
/*************************************************************************/

int main(
    int argc,                    /* Number of command line arguments */ 
    char **argv                 /* Command line arguments */
    )
{

    /* EQUILIBRIUM VARIABLES */
    
    int    n_tab;                     /* Number of points in EOS file */
    
    int i,s,m;
    
    double log_e_tab[TABP],            /* energy density/c^2 in tabulated EOS */
	log_p_tab[TABP],               /* pressure in tabulated EOS */
	log_h_tab[TABP],               /* enthalpy in EOS file */
	log_n0_tab[TABP],              /* number density in EOS file */  
	e_center,                     /* central en. density */
	p_center,                     /* central pressure */
	h_center,                     /* central enthalpy */
	rho0_center,                  /* central rest-mass density */
	n_P,                          /* Polytropic index N */
	Gamma_P,                      /* Gamma for polytropic EOS */     
	k_P,                          /* K for polytropic EOS */
	r_ratio,                      /* axis ratio */
	s_gp[SDIV+1],                 /* s grid points */
	mu[MDIV+1],                   /* \mu grid points */
	**rho,                        /* potential \rho */ 
	**gama,                       /* potential \gamma */ 
	**omega,                      /* potential \omega */ 
	**alpha,                      /* potential \alpha */ 
	Mass,                         /* Gravitational mass */
	e_surface,                    /* surface en. density */ 
	p_surface,                    /* surface pressure */
	enthalpy_min,                 /* minimum enthalpy in EOS */
	**energy,                     /* energy density \epsilon */
	**pressure,                   /* pressure */ 
	**enthalpy,                   /* enthalpy */
	**velocity_sq,                /* square of velocity */ 
	Mass_0,                       /* Baryon Mass */
	T,                            /* Kin rot energy */
	W,                            /* Grav bin energy */
	Omega,			      /* Angular Velocity */
	J,			      /* Angular Momentum */
	R_e,                          /* Circumferential radius at equator */
	*v_plus,		      /* vel. of co-rot. particle wrt ZAMO */
	*v_minus,		      /* vel. of counter-rot. ... */
	Omega_K,                      /* Keplerian velocity of particle orbiting at equator */
	r_e                           /* coord. radius at equator 	*/
	;
    
    double
	cf=1, /* convergence factor */
	accuracy,
	xacc,
	a_check
	;
    
    FILE *file2D; 
    
    double tmprr;
    
    /* SOME DEFAULT VALUES */
    
    char eos_file[80] = "no EOS file specified";   /* EOS file name */
    char eos_type[80] = "poly";                    /* EOS type (poly or tab) */
    char save_2Dmodel[10] = "no";
    char output_2Dfile[80] = "no output file specified"; 

    r_ratio = 1.0;
    n_P = 1.0;
    k_P = 100.0;
    Gamma_P = 1.0+1.0/n_P;
    e_center = 1.44e-3;          /* Standard dimensionless units !!! */
    
    /* READ IN THE COMMAND LINE OPTIONS */
    
    for(i=1;i<argc;i++) 
	if(argv[i][0]=='-'){
	    switch(argv[i][1]){
		
		case 'q':
		    /* CHOOSE THE EOS TYPE: EITHER "tab" or "poly"
		       (default is poly) */
		    sscanf(argv[i+1],"%s",eos_type);
		    break;
		    
		case 'N':
		    /* IF A POLYTROPIC EOS WAS CHOSEN, CHOOSE THE 
		       POLYTROPIC INDEX "N" */
		    sscanf(argv[i+1],"%lf",&n_P);		   
		    Gamma_P=1.0+1.0/n_P;
		    break;               
		    
		case 'K':
		    /* IF A POLYTROPIC EOS WAS CHOSEN, CHOOSE THE 
		       POLYTROPIC CONSTANT "K" */
		    sscanf(argv[i+1],"%lf",&k_P);
		    break;               
		    
		case 'f':
		    /* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
		       NAME OF THE FILE */
		    sscanf(argv[i+1],"%s",eos_file);
		    break;
		    
		case 'e':
		    /* CHOOSE THE CENTRAL ENERGY DENSITY OF THE NS
		       (IN g/cm^3 for TAB, IN G=c=Msun=1 units for POLY) */
		    sscanf(argv[i+1],"%lf",&e_center);
		    break;
		    
		case 'a':
		    /* CHOOSE R_RATIO
		       (default is 1) */
		    sscanf(argv[i+1],"%lf",&r_ratio); 
		    break;
		    
		case 's':
		    /* SAVE 2D OUTPUT ? */
		    sscanf(argv[i+1],"%s",save_2Dmodel);
		    break;
		    
		case 'O':
		    /* IF SAVE FILE OPTION WAS CHOSEN, CHOOSE THE
		       NAME OF THE OUTPUT FILE */
		    sscanf(argv[i+1],"%s",output_2Dfile);
		    break;
		    
		case 'h': 
		    fprintf(stderr,"\n");
		    fprintf(stderr,"Quick help:\n");
		    fprintf(stderr,"\n");
		    fprintf(stderr,"  -q EOS type {poly}\n"); 
		    fprintf(stderr,"     tab  : tabulated \n");
		    fprintf(stderr,"     poly : analytic polytropic \n");           
		    fprintf(stderr,"  -N polytropic index (P=K*e^(1+1/N)) {1}\n");  
		    fprintf(stderr,"  -K polytropic constant (P=K*e^(1+1/N)) {100}\n");  
		    fprintf(stderr,"  -f EOS file {none}\n");
		    fprintf(stderr,"  -e central energy density\n");
		    fprintf(stderr,"     [CGS ]       tab EOS (gr/cm^3)\n");
		    fprintf(stderr,"     [c=G=Msun=1] poly EOS \n");
		    fprintf(stderr,"  -a axies ratio {1}\n");
		    fprintf(stderr,"  -s save output file? yes/no {no} \n");
		    fprintf(stderr,"  -O output 2D filename {none}\n");
		    fprintf(stderr,"  -h this menu\n");
		    fprintf(stderr,"\n");
		    exit(1);
		    break;
		    
	    }
	}
    
    /* PRINT THE HEADER */

    if (VERBOSE) {

      printf(">>>\tRNS Code - start\n");    
      printf("SETUP :: MDIVxSDIV\t=\t%dx%d\n",MDIV,SDIV);
      printf("SETUP :: EOS\t\t=\t%s\n",eos_type);
      if(strcmp(eos_type,"tab")==0) {
	printf("SETUP :: EOS_tab\t=\t%s\n",eos_file);
	printf("SETUP :: e_c [gr/cm^3]\t=\t%g\n",e_center);
      }
      else {
	printf("SETUP :: EOS_N\t\t=\t%f\n",n_P);
	printf("SETUP :: EOS_K\t\t=\t%f\n",k_P);
	printf("SETUP :: e_c [c=G=Ms=1]\t=\t%g\n",e_center);
      }
      printf("SETUP :: r_p/r_e\t=\t%g\n",r_ratio);
    
    }

    /* LOAD TABULATED EOS */ 
    
    if(strcmp(eos_type,"tab")==0) 
	load_eos( eos_file, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, &n_tab );
    
    /* SET UP GRID */
    
    make_grid(s_gp, mu);
    
    /* ALLLOCATE MEMORY */
    
    rho = dmatrix(1,SDIV,1,MDIV);
    gama = dmatrix(1,SDIV,1,MDIV);
    alpha = dmatrix(1,SDIV,1,MDIV);
    omega = dmatrix(1,SDIV,1,MDIV);
    
    energy = dmatrix(1,SDIV,1,MDIV);
    pressure = dmatrix(1,SDIV,1,MDIV);
    enthalpy = dmatrix(1,SDIV,1,MDIV);
    velocity_sq = dmatrix(1,SDIV,1,MDIV);
    
    v_plus = dvector(1,SDIV);
    v_minus = dvector(1,SDIV);
    
    /* SET PROGRAM DEFAULTS AND 
       SWITCH UNITS TO POLYTROPIC DIMENSIONLESS UNITS */
    
    cf=1.0;
    accuracy=1e-6;    
    xacc = 1e-5;  
    
    if(strcmp(eos_type,"tab")==0) {
	
	e_center *= C*C*KSCALE;
	e_surface=7.8*C*C*KSCALE;
	p_surface=1.01e8*KSCALE;
	enthalpy_min=1.0/(C*C);

    } else {
	
	e_center /= ( 1.0/pow(k_P, n_P) );
	e_surface=0.0;
	p_surface=0.0;
	enthalpy_min=0.0;
	
    }
    
    /* CALCULATE THE PRESSURE, ENTHALPY AND RHO AT THE CENTRE OF THE STAR*/
    
    make_center(eos_file, log_e_tab, log_p_tab, 
		log_h_tab, log_n0_tab, n_tab,eos_type, Gamma_P, 
		e_center, &p_center, &h_center);
    
    rho0_center =  (e_center+p_center)*exp(-h_center);
    
    if (VERBOSE) {
      printf("SETUP :: e_c [poly]\t=\t%g\n",e_center);
      printf("SETUP :: rho_c [poly]\t=\t%g\n",rho0_center);
      printf("SETUP :: p_c [poly]\t=\t%g\n",p_center);
      printf("SETUP :: h_c [poly]\t=\t%g\n",h_center);
    }

    /* COMPUTE A SPHERICAL STAR AS A FIRST GUESS FOR THE ROTATING STAR */
    /* PRINT OUT INFORMATION ABOUT THE GUESS */
    
    sphere( s_gp, log_e_tab, log_p_tab, 
	    log_h_tab, log_n0_tab, n_tab, eos_type, Gamma_P, 
	    e_center, p_center, h_center, p_surface,e_surface,
	    rho, gama, alpha, omega, &r_e );
    
    if (VERBOSE) printf(">>>\tSpherical star guessed\n");
    
    /* COMPUTE A STAR WITH THE SPECIFID VALUE OF r_ratio	*/
    
    /* THE PROCEDURE SPIN() WILL COMPUTE THE METRIC OF A STAR WITH
       GIVEN OBLATENESS. THE OBLATENESS IS SPECIFIED BY GIVING 
       THE RATIO OF THE LENGTH OF THE AXIS CONNECTING THE CENTRE OF THE STAR 
       TO ONE OF THE POLES TO THE RADIUS OF THE STAR'S EQUATOR. 
       THIS RATIO IS NAMED r_ratio.
       WHEN r_ratio = 1.0, THE STAR IS SPHERICAL */
    
    /* THE METRIC FUNCTIONS ARE STORED IN THE FUNCTIONS 
       alpha, rho, gama, omega (see user's manual for the definition
       of the metric */
    
    /* If the axis ratio is less than 0.8 (0.6), one needs to compute 
       the model with r_ratio=0.8 (0.6) first and then compute the 
       desired model (otherwise the iteration may not converge). The 
       subroutine "iterate" starts with the current guess and 
       converges to the desired rotating model. The guess is either a 
       spherical star or a previously computed (slower) rotating 
       star. 
    */
    
    if (VERBOSE) printf(">>>\tIterating equilibrium model...\n");
    
    tmprr = r_ratio;
    
    if(r_ratio<0.9) {
	
	tmprr = 0.9;
	
	if (VERBOSE) printf(">>>\tRatio=%f\n",tmprr);
	
	spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	      n_tab, eos_type, Gamma_P, 
	      h_center, enthalpy_min,
	      rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	      a_check, accuracy, cf,
	      tmprr, &r_e, &Omega, VERBOSE );
	
    }
    
    if(r_ratio<0.8) {
	
	tmprr = 0.8;
	
	if (VERBOSE) printf(">>>\tRatio=%f\n",tmprr);
	
	spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	      n_tab, eos_type, Gamma_P, 
	      h_center, enthalpy_min,
	      rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	      a_check, accuracy, cf,
	      tmprr, &r_e, &Omega, VERBOSE );
	
    }      
    
    if(r_ratio<0.7) {
	
      tmprr = 0.7;
	
      if (VERBOSE)	printf(">>>\tRatio=%f\n",tmprr);
	
	spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	      n_tab, eos_type, Gamma_P, 
	      h_center, enthalpy_min,
	      rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	      a_check, accuracy, cf,
	      tmprr, &r_e, &Omega, VERBOSE); 
	
    }      
    
    if(r_ratio<0.6) {
   
	tmprr = 0.6;
 
	if (VERBOSE) printf(">>>\tRatio=%f\n",tmprr);
	
 	spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	      n_tab, eos_type, Gamma_P, 
	      h_center, enthalpy_min,
	      rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	      a_check, accuracy, cf,
	      tmprr, &r_e, &Omega, VERBOSE); 
	
    }      
 
    if(r_ratio<0.5) {
   
	tmprr = 0.5;
	
	if (VERBOSE) printf(">>>\tRatio=%f\n",tmprr);
	
	spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	      n_tab, eos_type, Gamma_P, 
	      h_center, enthalpy_min,
	      rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	      a_check, accuracy, cf,
	      tmprr, &r_e, &Omega, VERBOSE); 
	
    }      
    
    if(r_ratio<0.4) {
	
	tmprr = 0.4;
	
	if (VERBOSE) printf(">>>\tRatio=%f\n",tmprr);
	
	spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	      n_tab, eos_type, Gamma_P, 
	      h_center, enthalpy_min,
	      rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	      a_check, accuracy, cf,
	      tmprr, &r_e, &Omega, VERBOSE); 
	
    }      
    
    if(r_ratio<0.3) {
	
	tmprr = 0.3;
	
	if (VERBOSE) printf(">>>\tRatio=%f\n",tmprr);
	
	spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	      n_tab, eos_type, Gamma_P, 
	      h_center, enthalpy_min,
	      rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	      a_check, accuracy, cf,
	      tmprr, &r_e, &Omega, VERBOSE); 
	
    }      
    
    if(r_ratio<0.2) {
	
	tmprr = 0.2;
	
	if (VERBOSE) printf(">>>\tRatio=%f\n",tmprr);
	
	spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	      n_tab, eos_type, Gamma_P, 
	      h_center, enthalpy_min,
	      rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	      a_check, accuracy, cf,
	      tmprr, &r_e, &Omega, VERBOSE); 
	
    }      
    
    if(r_ratio<0.1) {
	
	tmprr = 0.1;
	
	if (VERBOSE) printf(">>>\tRatio=%f\n",tmprr);
	
	spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	      n_tab, eos_type, Gamma_P, 
	      h_center, enthalpy_min,
	      rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	      a_check, accuracy, cf,
	      tmprr, &r_e, &Omega, VERBOSE); 
	
    }      
    
    if(r_ratio<0.05) {
	
	tmprr = 0.05;
	
	if (VERBOSE) printf(">>>\tRatio=%f\n",tmprr);
	
	spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	      n_tab, eos_type, Gamma_P, 
	      h_center, enthalpy_min,
	      rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	      a_check, accuracy, cf,
	      tmprr, &r_e, &Omega, VERBOSE); 
	
  }      
    
    spin( s_gp, mu, log_e_tab, log_p_tab, log_h_tab, log_n0_tab, 
	  n_tab, eos_type, Gamma_P, 
	  h_center, enthalpy_min,
	  rho, gama, alpha, omega, energy, pressure, enthalpy, velocity_sq,
	  a_check, accuracy, cf,
	  r_ratio, &r_e, &Omega, VERBOSE); 
    
    if (VERBOSE)  printf(">>>\t...Done\n");
    
    /* COMPUTE THE VALUES OF VARIOUS EQUILIBRIUM QUANTITIES, SUCH
       AS MASS (Mass), RADIUS (R_e), BARYON MASS(Mass_0), 
       ANGULAR MOMENTUM (J), 
       KEPLERIAN ANGULAR VELOCITY OF PARTICLE ORBITING AT 
       THE EQUATOR,
       VELOCITIES OF CO-ROTATING PARTICLES (v_plus),
       AND COUNTER-ROTATING PARTICLES (v_minus) */
    
    if (VERBOSE) printf(">>>\tEquilibrium quantites\n");
    
    mass_radius( s_gp, mu, log_e_tab, log_p_tab, 
		 log_h_tab, log_n0_tab, n_tab, eos_type, Gamma_P, 
		 rho, gama, alpha, omega, 
		 energy, pressure, enthalpy, velocity_sq,
		 r_ratio, e_surface, r_e, Omega,
		 &Mass, &Mass_0, &T, &W, &J, &R_e, v_plus, v_minus, &Omega_K );
    
    /* TRANSFORM UNITS FROM POLYTROPIC TO STANDARD DIMENSIONLESS c=G=M_sun=1 */
    
    transform_units( eos_type, n_P, k_P, &rho0_center, &e_center,
		     &p_center, &r_e, omega, energy, pressure, &Mass,
		     &Mass_0, &T, &W, &Omega, &Omega_K, &R_e,
		     &Omega, &J );
    
    /* PRINT OUT INFORMATION ABOUT THE STELLAR MODEL */
    
    if(strcmp(eos_type,"tab")==0) 
	printtab(r_ratio,e_center, Mass, Mass_0, R_e, Omega, Omega_K, J);
    else 
	printpoly(r_ratio,e_center, Mass, Mass_0, R_e, Omega, Omega_K, J);
    
    
    if( strcmp(save_2Dmodel, "yes")==0) {
	
	/* SAVE 2D FILE */ 
	
	if((file2D=fopen(output_2Dfile,"w")) == NULL ) {    
	    printf("ERROR Cannot open file for saving 2D data:  %s\n",output_2Dfile); 
	    fflush(stdout);
	    exit(0);
	}
	
	fprintf(file2D, "#RNS initial data MDIVxSDIV=%dx%d\n",MDIV,SDIV); 
	fprintf(file2D, "#cols(10): s,mu,alpha,omega,gama,rho,e,p,h,v^2\n"); 
	fprintf(file2D, "#e_c=%13.12e r_ratio=%13.12e\n", e_center, r_ratio); 
	fprintf(file2D, "#r_e=%13.12e Omega=%13.12e\n", r_e, Omega); 
	fprintf(file2D, "#M=%13.12e R=%13.12e\n", Mass, R_e); 
	
	for(s=1;s<=SDIV;s++) 
	    for(m=1;m<=MDIV;m++) {
		
		fprintf(file2D, "%13.12e %13.12e %13.12e %13.12e %13.12e %13.12e %13.12e %13.12e %13.12e %13.12e\n",
			s_gp[s], mu[m], alpha[s][m], omega[s][m], gama[s][m], rho[s][m],
			energy[s][m], pressure[s][m], enthalpy[s][m],
			velocity_sq[s][m]);
	    }
	
	
	printf(">>>\tDone saving 2D file.\n");
	
    }  /* END SAVE 2D FILE */
    

    printf("%e\n",h_center);
    
    if (VERBOSE) printf(">>>\tRNS Code - now exiting\n");  
    fflush(stdout);

    return 0;
    
}








