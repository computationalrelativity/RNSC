
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
                 double *new,
                 int sign);

void hunt(double xx[], int n, double x, int *jlo);

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

double legendre( int n, double x );

double plgndr(int l, int m, double x);
 
double rtsec_G( double (*func)(double, double), 
                double Gamma_P,
                double x1, 
                double x2, 
                double xacc,
                double ee);

