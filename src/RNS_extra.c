#include "RNS.h"

/* Code not present in the original RNS code */

/*
 * Dummy/minimal stuff for managing parameters
 * - Params are stored as double and can be retrived only as int or real
 * - Input file parser assumes input file is correctly written with lines
 * 'key=val\n'
 * or
 * '# comment'
 */

parameters *params;

/* alloc/free  */
void params_alloc()
{
  params = calloc(sizeof(parameters), 1);
  params->n=0;
}
void params_free()
{
  free(params);
}

/* parse parameter file */
void params_read(char *fname) {

#define verbose (0)

  FILE *fp;
  char line[STRLEN], key[STRLEN], val[STRLEN];
  char *s;

  fp = fopen(fname, "r");
  if (fp == NULL)
    ERROR("[params_read] Failed to open input file");
  while ( fgets(line,sizeof line,fp) != NULL ) {
    /* TODO: better parsing */
    if (line[0]=='#') continue;
    s = strtok(line, "=");
    strcpy(key,s);
    /* s = strtok(NULL, s); */
    s = strtok(0, "=");
    strcpy(val,s);
    strtok(val, "\n");
    if (verbose && 0) printf("|%s|%s|\n",key,val);
    // Do not add parameters from input file, make sure they exist already
    for (int i =0; i < params->n; i++) {
      if (STREQL(key,params->key[i])) {
	strcpy(params->val[i],val);
	if (verbose) printf("[params_read] read: %s = %s\n",key,val);
	break; // goto next line
      }
    }

  }
  fclose(fp);
}

/* dump parameters to file */
void params_write(char * fname)
{
  FILE *fp;
  fp = fopen(fname, "w");
  if (fp == NULL)
    ERROR("[params_write] Failed to open file\n");
  for (int i =0; i < params->n; i++) {
    fprintf(fp,"%s=%s\n",params->key[i],params->val[i]);
  }
  fclose(fp);
}

/* get double */
double params_get_real(char * key)
{
  for (int i =0; i < params->n; i++)  {
    if (STREQL(key,params->key[i])) {
      if (params->type[i]==REAL) {
	return atof(params->val[i]);
      }  else
	ERROR("[params_get] parameter is not of type REAL\n");
    }
  }
  ERROR("[params_get] parameter not found");
}

/* get int */
int params_get_int(char * key) 
{
  for (int i =0; i < params->n; i++)  {
    if (STREQL(key,params->key[i])) {
      if (params->type[i]==INTEGER) 
	return atoi(params->val[i]);
      else
	ERROR("[params_get] parameter is not of type INTEGER\n");
    }
  }
  ERROR("[params_get] parameter not found");
}

/* get string */
char * params_get_str(char * key)
{
  for (int i =0; i < params->n; i++)  {
    if (STREQL(key,params->key[i])) {
      if (params->type[i]==STRING) 
	return params->val[i];
      else
	ERROR("[params_get] parameter is not of type STRING\n");
    }
  }
  ERROR("[params_get] parameter not found");
}

/* set parameter */
void params_setadd(char * key, int type, char* val, int addpar)
{
#define check_unique (0) // addpar=1 assumes parameter is not there,
			 // activate this if you need to enforce this.
  
  if (addpar) {
#if (check_unique)
    //Do nothing if present
    for (int i =0; i < params->n; i++)
      if (strcmp(params->key[i],key) == 0) return;
#endif
    int i = params->n;
    if (i>=NPARAMS) ERROR("Increase NPARAMS.");
    if (type>=N_PAR_TYPES) ERROR("Unknown parameter type.");
    strcpy(params->key[i],key);
    strcpy(params->val[i],val);    
    params->type[i] = type;
    params->n++;
    if (0) printf ("Add parameter: %s TYPE=%s VALUE=%s (%d)\n",
		   key,str_par_type[type],val,params->n);
    return;
  }

  for (int i =0; i < params->n; i++)  {
    if (STREQL(key,params->key[i])) {
      strcpy(params->val[i],val);
      return;
    }
  }

  ERROR("[params_set] parameter not found\n");
}

/* 'set' wrappers 
   Replace values, do not append new (assume exists)
   Note: addpar=0 does not set the type
*/

void params_set_int(char * key, int val) {
  char cval[STRLEN];
  sprintf(cval, "%d", val);
  params_setadd(key, INTEGER, cval, 0); 
}
 
void params_set_bool(char * key, bool val) {
  char cval[STRLEN];
  sprintf(cval, "%d", val?1:0);
  params_setadd(key, INTEGER, cval, 0); 
}

void params_set_real(char * key, double val) {
  char cval[STRLEN];
  sprintf(cval, "%.16e", val);
  params_setadd(key, REAL, cval, 0); 
}

void params_set_str(char * key, char *val) {
  params_setadd(key, STRING, val, 0);
}

/* 'add' wrappers 
   New parameter, append
 */

void params_add_int(char * key, int val) {
  char cval[STRLEN];
  sprintf(cval, "%d", val);
  params_setadd(key, INTEGER, cval, 1);
}

void params_add_bool(char * key, bool val) {
  char cval[STRLEN];
  sprintf(cval, "%d", val?1:0);
  params_setadd(key, INTEGER, cval, 1); 
}

void params_add_real(char * key, double val) {
  char cval[STRLEN];
  sprintf(cval, "%.16e", val);
  params_setadd(key, REAL, cval, 1);
}

void params_add_str(char * key, char *val) {
  params_setadd(key, STRING, val, 1);
}

void make_output_dir()
/* note this is unsafe */
{
  char s[STRLEN*4];
  sprintf(s,"mkdir -p %s",params_get_str("outputdir"));
  if (system(s)==-1) ERROR("System call failed to mkdir");
  strcpy(s,params_get_str("outputdir"));
  strcat (s,"/params.par");
  params_write(s);
}


/*
 * Formatted output routines
 * TODO: units, adapt output to rotation and eos type
 */

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
			     double Omega_e)
{
  const double I = (Omega != 0.0)? J/Omega : 0.0 ;
  const double I_45 = I/(1.0e45);
  const double J_units = C/(G*Mass*Mass);
  const double mass_units = MSUN;
  
  printf("# 1:rho_center");
  printf(" 2:axes_ratio");
  printf(" 3:e_center");
  printf(" 4:Mass");
  printf(" 5:Mass_0");
  printf(" 6:R_e");
  printf(" 7:J");
  printf(" 8:J/M^2");
  printf(" 9:T");
  printf(" 10:W");
  printf(" 11:T/W");
  printf(" 12:Omega");
  printf(" 13:Omega_K");
  printf(" 14:Omega_e");
  printf(" 15:I");
  printf("\n");
  
  printf("%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n", 
	 rho_center,
	 a_ratio,
	 e_center,
	 Mass, 
	 Mass_0, 
	 R_e,
	 (Omega > 0.0) ? J : 0.0,
	 (Omega > 0.0) ? J/(Mass*Mass) : 0.0,
	 T,
	 W,
	 T/W,
	 Omega,
	 Omega_K,
	 Omega_e,
	 I);
  
}


/*
 * miscellanea
 */


/* compute v^i = g^{ij} v_j (or v_i = g_{ij} v^j ) 
   and v_i v^i */
void contract_vect_v2(double vlx, double vly, double vlz,
		      double gxx, double gxy, double gxz,
		      double gyy, double gyz, double gzz,
		      double *vx, double *vy, double *vz,
		      double *v2)
{
  *vx = gxx*vlx + gxy*vly + gxz*vlz;
  *vy = gxy*vlx + gyy*vly + gyz*vlz;
  *vz = gxz*vlx + gyz*vly + gzz*vlz; 
  *v2 = (*vx)*vlx + (*vy)*vly + (*vz)*vlz;
}

