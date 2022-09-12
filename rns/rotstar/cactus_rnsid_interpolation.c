/* ===================================================================== */
/* ===================================================================== */
/* ===================================================================== */
/* CACTUS/WHISKY INTERPOLATION OF THE STAR ON THE CARTESIAN GRID */
/* THIS CODE DOESN'T WORK: IT IS JUST A CUT&PASTE FOR REFERENCE  */
/* ===================================================================== */
/* ===================================================================== */
/* ===================================================================== */


void rnsid_cactus() 

{

    /* ===================================================================== */
    /* ===================================================================== */
    /* ===================================================================== */
    /* ABOVE IS EXACLTY AS ROTSTAR ....THEN INTERPOLATION ON 3D GRID:        */
    /* ===================================================================== */
    /* ===================================================================== */
    /* ===================================================================== */

      /* CONSTRUCT ARRAYS WITH NEEDED POLAR QUANTITIES */


      nu = dmatrix(1,SDIV,1,MDIV);
      B = dmatrix(1,SDIV,1,MDIV);
      rho_0 = dmatrix(1,SDIV,1,MDIV);      

      for(m=1;m<=MDIV;m++) 
         for(s=1;s<=SDIV;s++) {
 	      nu[s][m] = (gama[s][m]+rho_potential[s][m])/2.0; 
	      B[s][m] = exp(gama[s][m]);	       
              rho_0[s][m] = (energy[s][m]+pressure[s][m])
                             *exp(-enthalpy[s][m]);
         } 
      
      free_dmatrix(rho_potential,1,SDIV,1,MDIV);
      free_dmatrix(gama,1,SDIV,1,MDIV);
      free_dmatrix(enthalpy,1,SDIV,1,MDIV);
      free_dmatrix(velocity_sq,1,SDIV,1,MDIV);
      
      nu_dr = dmatrix(1,SDIV,1,MDIV);
      B_dr = dmatrix(1,SDIV,1,MDIV);
      alpha_dr = dmatrix(1,SDIV,1,MDIV);
      omega_dr = dmatrix(1,SDIV,1,MDIV);
      nu_dth = dmatrix(1,SDIV,1,MDIV);
      B_dth = dmatrix(1,SDIV,1,MDIV);
      alpha_dth = dmatrix(1,SDIV,1,MDIV);
      omega_dth = dmatrix(1,SDIV,1,MDIV);

      for(m=1;m<=MDIV;m++) 
         for(s=1;s<=SDIV;s++) {
            nu_dr[s][m] = deriv_s(nu,s,m)*SQ(1.0-s_gp[s])/r_e;
            B_dr[s][m] = deriv_s(B,s,m)*SQ(1.0-s_gp[s])/r_e;
            alpha_dr[s][m] = deriv_s(alpha,s,m)*SQ(1.0-s_gp[s])/r_e;
            omega_dr[s][m] = deriv_s(omega,s,m)*SQ(1.0-s_gp[s])/r_e;
            nu_dth[s][m] = deriv_m(nu,s,m)*(-sqrt(1.0-SQ(mu[m])));
            B_dth[s][m] = deriv_m(B,s,m)*(-sqrt(1.0-SQ(mu[m])));
            alpha_dth[s][m] = deriv_m(alpha,s,m)*(-sqrt(1.0-SQ(mu[m])));
            omega_dth[s][m] = deriv_m(omega,s,m)*(-sqrt(1.0-SQ(mu[m])));
         }  

      /* COMPUTE INITIAL DATA */

      rho_0_atm = rnsid_rho_min; /* rename the constant for historical reasons */

      e_atm = rho_0_atm;
      p_atm = eos_k*pow(rho_0_atm,eos_ideal_fluid_gamma);

      
      for(i=1;i<=nx;i++)
         for(j=1;j<=ny;j++) 
           for(k=1;k<=nz;k++) {

               x_i = x_grid[i-1+nx*(j-1+ny*(k-1))];
               y_j = y_grid[i-1+nx*(j-1+ny*(k-1))];
               z_k = z_grid[i-1+nx*(j-1+ny*(k-1))];

               grid_interp_all( s_gp, mu, r_e, nx, ny, nz, 
                                x_grid, y_grid, z_grid,
                                i, j, k,
                                nu, B, alpha, omega, 
                                nu_dr, B_dr, alpha_dr, omega_dr, 
                                nu_dth, B_dth, alpha_dth, omega_dth,
                                rho_0, energy, pressure,
                                &nu_ijk, &B_ijk, &alpha_ijk, &omega_ijk, 
                                &nu_dr_ijk, &B_dr_ijk, &alpha_dr_ijk, 
                                &omega_dr_ijk,             
                                &nu_dtheta_ijk, &B_dtheta_ijk, 
                                &alpha_dtheta_ijk, &omega_dtheta_ijk,
                                &rho_0_ijk, &energy_ijk, &pressure_ijk,
                                &distance_ijk,
				Omega_diff, &Omega_ijk);

              if( (rho_0_ijk<=0.0) || (energy_ijk<=0.0) || 
                  (pressure_ijk<=0.0) ) {
 
                rho_0_ijk= rho_0_atm;
                energy_ijk= e_atm + 1.e-20;
                pressure_ijk= p_atm;
                
              }
	      
              exp_nu_ijk = exp(nu_ijk);
              exp_alpha_ijk = exp(alpha_ijk);
 
              r_ijk = sqrt(SQ(x_i)+SQ(y_j)+SQ(z_k));
              r_bar_ijk = sqrt(SQ(x_i)+SQ(y_j));


              alp[i-1+nx*(j-1+ny*(k-1))] = exp_nu_ijk;


              if(x_i==0.0 && y_j==0.0) {
  
                gxx[i-1+nx*(j-1+ny*(k-1))] = SQ(exp_alpha_ijk);
                gyy[i-1+nx*(j-1+ny*(k-1))] = SQ(exp_alpha_ijk);
                gzz[i-1+nx*(j-1+ny*(k-1))] = SQ(exp_alpha_ijk);

                gxy[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                gxz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                gyz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
  
                kxx[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                kyy[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                kzz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                kxy[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                kxz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                kyz[i-1+nx*(j-1+ny*(k-1))] = 0.0;

              }else{
 
                dr_dx = x_i / r_ijk;
                dr_dy = y_j / r_ijk;
                dr_dz = z_k / r_ijk;

                dtheta_dx = x_i*z_k/(SQ(r_ijk)*r_bar_ijk);
                dtheta_dy = y_j*z_k/(SQ(r_ijk)*r_bar_ijk);
                dtheta_dz = - r_bar_ijk/SQ(r_ijk);
 
                nu_dx = dr_dx*nu_dr_ijk + dtheta_dx*nu_dtheta_ijk;
                nu_dy = dr_dy*nu_dr_ijk + dtheta_dy*nu_dtheta_ijk;

                B_dx = dr_dx*B_dr_ijk + dtheta_dx*B_dtheta_ijk;
                B_dy = dr_dy*B_dr_ijk + dtheta_dy*B_dtheta_ijk;
 
                alpha_dx = dr_dx*alpha_dr_ijk + dtheta_dx*alpha_dtheta_ijk;
                alpha_dy = dr_dy*alpha_dr_ijk + dtheta_dy*alpha_dtheta_ijk;
  
                omega_dx = dr_dx*omega_dr_ijk + dtheta_dx*omega_dtheta_ijk;
                omega_dy = dr_dy*omega_dr_ijk + dtheta_dy*omega_dtheta_ijk;

		/* enforce omega_dz=0 at z=0 (it is slightly nonzero due
		   to O(h) forwards formula in computing derivative) */

                if(z_k==0.0) 
                  omega_dz = 0.0;
                else
                  omega_dz = dr_dz*omega_dr_ijk + dtheta_dz*omega_dtheta_ijk;

 
                gxx[i-1+nx*(j-1+ny*(k-1))] = ( SQ(B_ijk*y_j/exp_nu_ijk) 
                                             +SQ(exp_alpha_ijk*x_i) ) /
                                             (SQ(x_i)+SQ(y_j));

                gxy[i-1+nx*(j-1+ny*(k-1))] = ( SQ(exp_alpha_ijk) 
                                            -SQ(B_ijk/exp_nu_ijk) ) *
                                            x_i*y_j/(SQ(x_i)+SQ(y_j));

                gxz[i-1+nx*(j-1+ny*(k-1))] = 0.0;

                gyy[i-1+nx*(j-1+ny*(k-1))] = ( SQ(B_ijk*x_i/exp_nu_ijk) 
                                             +SQ(exp_alpha_ijk*y_j) ) /
                                             (SQ(x_i)+SQ(y_j));
  
                gyz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
  
                gzz[i-1+nx*(j-1+ny*(k-1))] = SQ(exp_alpha_ijk);
 

                kxx[i-1+nx*(j-1+ny*(k-1))] = (  ( SQ(r_bar_ijk)*y_j*omega_dx+
                                             (x_i*nu_dy-y_j*nu_dx)*SQ(y_j)
                                             *omega_ijk)*SQ(B_ijk)         
                                             +(y_j*B_dx-x_i*B_dy)*omega_ijk
                                             *SQ(y_j)*B_ijk
                                             +(y_j*alpha_dx-x_i*alpha_dy)
                                             *omega_ijk*SQ(x_i*exp_alpha_ijk
                                             *exp_nu_ijk))/(SQ(r_bar_ijk
                                             *exp_nu_ijk)*exp_nu_ijk);
		
                kxy[i-1+nx*(j-1+ny*(k-1))] = ( ( 0.5*SQ(r_bar_ijk)*
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
	        
                kxz[i-1+nx*(j-1+ny*(k-1))] = 0.5*SQ(B_ijk)*y_j*omega_dz/
                                             ( SQ(exp_nu_ijk)*exp_nu_ijk );
   
                kyy[i-1+nx*(j-1+ny*(k-1))] = ( ( -SQ(r_bar_ijk)*x_i*omega_dy+
                                             (x_i*nu_dy-y_j*nu_dx)*SQ(x_i)* 
                                             omega_ijk )*SQ(B_ijk) 
                                             +(y_j*B_dx-x_i*B_dy)*omega_ijk
                                             *SQ(x_i)*B_ijk
                                             +(y_j*alpha_dx-x_i*alpha_dy)
                                             *omega_ijk*SQ(y_j*exp_alpha_ijk
                                             *exp_nu_ijk))/(SQ(r_bar_ijk
                                             *exp_nu_ijk)*exp_nu_ijk);
	       
                kyz[i-1+nx*(j-1+ny*(k-1))] = -0.5*SQ(B_ijk)*x_i*omega_dz/
                                             ( SQ(exp_nu_ijk)*exp_nu_ijk );
 
                kzz[i-1+nx*(j-1+ny*(k-1))] = (y_j*alpha_dx-x_i*alpha_dy)*
                                             omega_ijk*SQ(exp_alpha_ijk)/
                                             exp_nu_ijk;
	      }
 
	       
              betax[i-1+nx*(j-1+ny*(k-1))] = omega_ijk*y_j;
              betay[i-1+nx*(j-1+ny*(k-1))] = -omega_ijk*x_i;
	       
              betaz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
	       
              rho0[i-1+nx*(j-1+ny*(k-1))] = rho_0_ijk;
	       
              eps[i-1+nx*(j-1+ny*(k-1))] = energy_ijk/rho_0_ijk-1.0;
              h_ijk = (energy_ijk+pressure_ijk)/rho_0_ijk;
             
              gamma_ijk = SQ(exp_alpha_ijk)*SQ(B_ijk*exp_alpha_ijk/
                          exp_nu_ijk); 

              W_ijk = 1.0/sqrt(1.0-SQ((omega_ijk-Omega_ijk)*B_ijk*
                      r_bar_ijk/SQ(exp_nu_ijk)));

              dens[i-1+nx*(j-1+ny*(k-1))] = sqrt(gamma_ijk)*W_ijk*rho_0_ijk;

              tau[i-1+nx*(j-1+ny*(k-1))] = sqrt(gamma_ijk)*( rho_0_ijk*h_ijk
                                           *SQ(W_ijk)-pressure_ijk - 
                                           W_ijk*rho_0_ijk );
              
/*
	      sx[i-1+nx*(j-1+ny*(k-1))] = sqrt(gamma_ijk)*rho_0_ijk*h_ijk
                                          *SQ(W_ijk)*
                                          (omega_ijk-Omega)*y_j/exp_nu_ijk;

              sy[i-1+nx*(j-1+ny*(k-1))] = - sqrt(gamma_ijk)*rho_0_ijk*h_ijk
                                          *SQ(W_ijk)*
                                          (omega_ijk-Omega)*x_i/exp_nu_ijk;

              sz[i-1+nx*(j-1+ny*(k-1))] = 0.0;

	      velx[i-1+nx*(j-1+ny*(k-1))] = (omega_ijk-Omega)*y_j/exp_nu_ijk;

              vely[i-1+nx*(j-1+ny*(k-1))] = -(omega_ijk-Omega)*x_i/exp_nu_ijk;

              velz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
*/

	      sx[i-1+nx*(j-1+ny*(k-1))] = sqrt(gamma_ijk)*rho_0_ijk*h_ijk
                                          *SQ(W_ijk*B_ijk/exp_nu_ijk)*
                                          (omega_ijk-Omega_ijk)*y_j/exp_nu_ijk;

              sy[i-1+nx*(j-1+ny*(k-1))] = - sqrt(gamma_ijk)*rho_0_ijk*h_ijk
                                          *SQ(W_ijk*B_ijk/exp_nu_ijk)*
                                          (omega_ijk-Omega_ijk)*x_i/exp_nu_ijk;

              sz[i-1+nx*(j-1+ny*(k-1))] = 0.0;

	      velx[i-1+nx*(j-1+ny*(k-1))] = (omega_ijk-Omega_ijk)
		                              *y_j/exp_nu_ijk;

              vely[i-1+nx*(j-1+ny*(k-1))] = -(omega_ijk-Omega_ijk)
		                              *x_i/exp_nu_ijk;

              velz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
	      
	      dens_atm = sqrt(gamma_ijk)*rho_0_atm;
              tau_atm = sqrt(gamma_ijk)*eos_k*pow(rho_0_atm,eos_ideal_fluid_gamma) / 
                (eos_ideal_fluid_gamma - 1.0);
	      if ( (dens[i-1+nx*(j-1+ny*(k-1))] < dens_atm)||
		   (tau[i-1+nx*(j-1+ny*(k-1))] < tau_atm)||
                   (rho0[i-1+nx*(j-1+ny*(k-1))] < (1.0 + whisky_atmo_tolerance) * rho_0_atm) ) {
		dens[i-1+nx*(j-1+ny*(k-1))] = dens_atm;
		tau[i-1+nx*(j-1+ny*(k-1))] = tau_atm;
		sx[i-1+nx*(j-1+ny*(k-1))] = 0.0;
		sy[i-1+nx*(j-1+ny*(k-1))] = 0.0;
		sz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
		velx[i-1+nx*(j-1+ny*(k-1))] = 0.0;
		vely[i-1+nx*(j-1+ny*(k-1))] = 0.0;
		velz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
	      } 

              if (tracer)
              {
                tracer[i-1+nx*(j-1+ny*(k-1))] = distance_ijk;
              }

	   } /* END FOR LOOP OF LINE 525 */




      if( strcmp(zero_shift, "yes")==0) {

        /* SET SHIFT TO ZERO */
 
        for(i=1;i<=nx;i++)
           for(j=1;j<=ny;j++)
              for(k=1;k<=nz;k++) {
                 betax[i-1+nx*(j-1+ny*(k-1))] = 0.0;
                 betay[i-1+nx*(j-1+ny*(k-1))] = 0.0;     
                 betaz[i-1+nx*(j-1+ny*(k-1))] = 0.0;
              }
      }

      /* compute central value of 3-determinant of metric */

      *gamma_center = SQ(B[1][1])*exp( 4.0*alpha[1][1] - 2.0*nu[1][1]);

      /* FREE MEMORY */
	      
      free_dmatrix(alpha,1,SDIV,1,MDIV);
      free_dmatrix(omega,1,SDIV,1,MDIV);      
      free_dmatrix(rho_0,1,SDIV,1,MDIV);      
      free_dmatrix(energy,1,SDIV,1,MDIV);
      free_dmatrix(pressure,1,SDIV,1,MDIV);

      free_dmatrix(nu,1,SDIV,1,MDIV);
      free_dmatrix(B,1,SDIV,1,MDIV);

      free_dmatrix(nu_dr,1,SDIV,1,MDIV);
      free_dmatrix(B_dr,1,SDIV,1,MDIV);
      free_dmatrix(alpha_dr,1,SDIV,1,MDIV);
      free_dmatrix(omega_dr,1,SDIV,1,MDIV);
      free_dmatrix(nu_dth,1,SDIV,1,MDIV);
      free_dmatrix(B_dth,1,SDIV,1,MDIV);
      free_dmatrix(alpha_dth,1,SDIV,1,MDIV);
      free_dmatrix(omega_dth,1,SDIV,1,MDIV);

      free(s_gp);
      free(mu);
    
    

      return;

} /* END RNSID */

/* ===================================================================== */
/* ===================================================================== */
/* ===================================================================== */
/* HERE BELOW THE ROUTINES USED (THE IMPORTANT ONES ARE ALSO IN ROTSTAR) */
/* ===================================================================== */
/* ===================================================================== */
/* ===================================================================== */

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

 np = 0; /* dummy assignement: unused variable*/

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

 DECLARE_CCTK_PARAMETERS;
	    
 nz = 0.0; /* dummy assignement: unused variable*/

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
  
  distance_c = 0; /* dummy assignement: unused variable*/

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

      /*
      find_surface( s_gp, mu, r_e, nx, ny, nz, x_grid, y_grid, z_grid,
                    i, j, k, distance);
      */
}
