/* RNS_rnsid_util.c */
#include "RNS.h"

#define NR_END 1
#define FREE_ARG char*
 
/*************************************************************************
* Allocate memory for double 3tensor (adapted from num. rec. f3tensor
*************************************************************************/
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;

	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in d3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}



/*************************************************************************
* Free memory of double 3tensor (adapted from num. rec. f3tensor
*************************************************************************/
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a double d3tensor allocated by d3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}



/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */  
/*************************************************************************/
double interp_4(double xp[5], 
                double yp[5], 
                int    np ,
                double xb)
{ 
 int k=1;      /* index of 1st point */
 
 double y;     /* intermediate value */


 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) 
    xb += DBL_EPSILON;

 /*
 if( yp[k]==yp[k+1] &&  yp[k]==yp[k+2] && yp[k]==yp[k+3]) 
    y=yp[k];
 else{ */ 
 y= (xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
       ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
       ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
       ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));
 /* }*/
 return (y);
 
}



/*************************************************************************
* Interpolation between two different grids
*************************************************************************/
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
                 double *new,
                 int sign) {

  int MDIV,SDIV;
  SDIV        = params_get_int("rns_SDIV");
  MDIV        = params_get_int("rns_MDIV");

int s,
    m,
    s_nearest,
    m_nearest,                /* nearest points in interpolation */
    k_s,                      /* first s point in interpolation */
    k_m;                      /* first s point in interpolation */


double r_c,                   /* r of cartesian x,y,z point */
       s_c,                   /* s of cartesian x,y,z point */
       mu_c,                  /* mu of cartesian x,y,z point */
       s_4[5],                /* s of the 4 nearest points */
       mu_4[5],               /* mu of the 4 nearest points */
       old_s[5],              /* old at 4 nearest constant s points */
       old_m[5];              /* old at 4 nearest constant mu points */  


	    
               r_c = sqrt(   SQ(x_grid[i-1+nx*(j-1+ny*(k-1))])   
                                + SQ(y_grid[i-1+nx*(j-1+ny*(k-1))]) 
                                + SQ(z_grid[i-1+nx*(j-1+ny*(k-1))]) ); 
               s_c = r_c/(r_e+r_c);
      
               if(r_c==0.0) 
                 mu_c = 0.0;
               else
                 mu_c = fabs(z_grid[i-1+nx*(j-1+ny*(k-1))])/r_c;
 
               s_nearest = 0; m_nearest = 0;

               nr_hunt(s_gp, SDIV, s_c, &s_nearest);
               nr_hunt(mu, MDIV, mu_c, &m_nearest);
  
               k_s = IMIN(IMAX((s_nearest)-(4-1)/2,1),SDIV+1-4);
               k_m = IMIN(IMAX((m_nearest)-(4-1)/2,1),MDIV+1-4);
	       	       
               for(s=1;s<=4;s++) 
                     s_4[s] = s_gp[k_s-1+s];

               for(m=1;m<=4;m++) 
                     mu_4[m] = mu[k_m-1+m];
 
               for(s=1;s<=4;s++) {
                     for(m=1;m<=4;m++) {
                           old_s[m] = old[k_s-1+s][k_m-1+m];
	         }
                     old_m[s] = interp_4(mu_4, old_s, 4, mu_c);  
               }

         if(z_grid[i-1+nx*(j-1+ny*(k-1))]<0.0)                               
              (*new) = (1.0*sign)*interp_4(s_4, old_m, 4, s_c); 
         else 
              (*new) = interp_4(s_4, old_m, 4, s_c); 
}



/*************************************************************************
* Interpolation between two different grids - all
*************************************************************************/
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
		      double *Omega_diff_c)
{  
  
      grid_interp( nu, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, nu_c, 1);
      
      grid_interp( B, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, B_c,  1);  

      grid_interp( alpha, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, alpha_c,  1);

      grid_interp( omega, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, omega_c,  1);

      grid_interp( nu_dr, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, nu_dr_c, 1);

      grid_interp( B_dr, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, B_dr_c, 1);

      grid_interp( alpha_dr, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, alpha_dr_c, 1);

      grid_interp( omega_dr, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, omega_dr_c, 1);
		   
      grid_interp( nu_dth, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, nu_dth_c, -1);
      
      grid_interp( B_dth, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, B_dth_c, -1);

      grid_interp( alpha_dth, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, alpha_dth_c, -1);

      grid_interp( omega_dth, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, omega_dth_c, -1);
		   	   
      grid_interp( rho_0, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, rho_0_c, 1);
      
      grid_interp( energy, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, energy_c, 1);

      grid_interp( pressure, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, pressure_c, 1);       	       

      grid_interp( Omega_diff, s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid,
                   z_grid, i, j, k, Omega_diff_c, 1);

}


/*************************************************************************
* Interpolation between two different grids
*************************************************************************/
void grid_interp_new( int k_s, int k_m, int sign, 
		      double s_c, double mu_c, double zp, 
		      double *s_4, double *mu_4, 
		      double **old, 
		      double *new) 
{
  
  int s,m;
  double old_s[5],old_m[5];   
  
  for(s=1;s<=4;s++) {
    for(m=1;m<=4;m++) {
      old_s[m] = old[k_s-1+s][k_m-1+m];
    }
    old_m[s] = interp_4(mu_4, old_s, 4, mu_c);  
  }
  
  if(zp<0.)                               
    (*new) = (1.0*sign)*interp_4(s_4, old_m, 4, s_c); 
  else 
    (*new) = interp_4(s_4, old_m, 4, s_c); 
  
}



/*************************************************************************
* Interpolation between two different grids - all
*************************************************************************/
void grid_interp_all_new ( double *s_gp, 
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
			   double *Omega_diff_c)
{  
 

  int MDIV,SDIV;
  SDIV        = params_get_int("rns_SDIV");
  MDIV        = params_get_int("rns_MDIV");

 
  int s,
    m,
    s_nearest,            /* nearest points in interpolation */
    m_nearest,               
    k_s,                  /* first s point in interpolation */
    k_m;                  /* first m point in interpolation */
    
  double r_c,              /* r of cartesian x,y,z point */
    s_c,                   /* s of cartesian x,y,z point */
    mu_c,                  /* mu of cartesian x,y,z point */
    s_4[5],                /* s of the 4 nearest points */
    mu_4[5];               /* mu of the 4 nearest points */

  const int sign = 1;



  /* FIND NN */

  r_c = sqrt( SQ(xp) + SQ(yp) + SQ(zp) ); 
  s_c = r_c/(r_e + r_c);
  
  if(r_c==0.) mu_c = 0.;
  else        mu_c = fabs(zp)/r_c;
  
  s_nearest = 0; m_nearest = 0;
 
  nr_hunt(s_gp, SDIV, s_c, &s_nearest);
  nr_hunt(mu, MDIV, mu_c, &m_nearest);
 
  k_s = IMIN(IMAX((s_nearest)-(4-1)/2,1),SDIV+1-4);
  k_m = IMIN(IMAX((m_nearest)-(4-1)/2,1),MDIV+1-4);
 
  for(s=1;s<=4;s++) 
    s_4[s] = s_gp[k_s-1+s];
  
  for(m=1;m<=4;m++) 
    mu_4[m] = mu[k_m-1+m];
  

  /* INTERPOLATE EACH VAR, 4 pts interp */

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   nu, nu_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   B, B_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   alpha, alpha_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   omega, omega_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   nu_dr, nu_dr_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   B_dr, B_dr_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   alpha_dr, alpha_dr_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   omega_dr, omega_dr_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   nu_dth, nu_dth_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   B_dth, B_dth_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   alpha_dth, alpha_dth_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   omega_dth, omega_dth_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   rho_0, rho_0_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   energy , energy_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   pressure, pressure_c );

  grid_interp_new( k_s, k_m, sign, 
		   s_c, mu_c, zp, 
		   s_4, mu_4, 
		   Omega_diff, Omega_diff_c );

}


/*************************************************************************
* Transform units to c=G=M_sun=1
*************************************************************************/
void transform_units( 
                 char   eos_type[],
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
		 double *J)
{ 
int i,                                       /* counter */
      j;                                      /* counter */
 

  int MDIV,SDIV;
  SDIV        = params_get_int("rns_SDIV");
  MDIV        = params_get_int("rns_MDIV");

   for(i=1;i<=SDIV;i++)
        for(j=1;j<=MDIV;j++) {
             omega[i][j] *= ( 1.0/pow(eos_k, n_P/2.0) );
             energy[i][j] *= ( 1.0/pow(eos_k, n_P) );
             pressure[i][j] *= ( 1.0/pow(eos_k, n_P) );
             Omega_diff[i][j] *= ( 1.0/pow(eos_k, n_P/2.0) );
        }

   (*rho0_center) *= ( 1.0/pow(eos_k, n_P) );
   (*e_center) *= ( 1.0/pow(eos_k, n_P) );
   (*p_center) *= ( 1.0/pow(eos_k, n_P) );
   (*r_e) *= pow(eos_k, n_P/2.0);
   (*Mass) *= pow(eos_k, n_P/2.0);
   (*Mass_0) *= pow(eos_k, n_P/2.0);
   (*T) *= pow(eos_k, n_P/2.0);
   (*W) *= pow(eos_k, n_P/2.0);
   (*Omega) *= ( 1.0/pow(eos_k, n_P/2.0) );
   (*Omega_K) *= ( 1.0/pow(eos_k, n_P/2.0) );
   (*R_e) *= pow(eos_k, n_P/2.0);
   (*Omega_e) *= ( 1.0/pow(eos_k, n_P/2.0) );
   (*J) *= pow(eos_k, n_P);
} 
