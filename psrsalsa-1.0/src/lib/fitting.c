/*
  The double precision functions are in fitting.c, and the other
  precision code (fitting_ld are long doubles) is derived from
  fitting.c. So always change the code in fitting.c and re-create the
  other precisions code. 

  For other precisions: Compile without math.h to see if all functions
  are in correct precision.

  double                              -> long double
  fit_levmar2                         -> fit_levmar2_ld
  linalg_solve_matrix_eq_gauss_jordan -> linalg_solve_matrix_eq_gauss_jordan_ld
  fabs                                -> fabsl
  sqrt                                -> sqrtl
  %e                                  -> %Le
  1.2e-16                             -> 1.1e-19 (long double) 1.2e-16 (double) 6.0e-8 (float)

 */

//START REGION DEVELOP
//START REGION RELEASE
#include <math.h>           
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_version.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include "psrsalsa.h"

//START REGION DEVELOP
//START REGION RELEASE

/* Information used by the different fitters */
typedef struct {
  size_t n;
  double *x;
  double *y;
  double *sigma;
  fitfunc_collection_type *fitfunction;
  size_t iter;   // Total number of itterations performed to obtain solution
}levmar_internal_fitter_data_def;

// Needs to be globally defined, since NR doesn't pass on a void pointer to a user defined information struct
static levmar_internal_fitter_data_def levmar_internal_fitter_data;


/* These are gsl specific parameters that need to be passed on between functions */
typedef struct {
  gsl_vector_view func_params;
  const gsl_multifit_fdfsolver_type *solver_type;
  gsl_multifit_function_fdf user_itteration_funcs;
  gsl_multifit_fdfsolver *solver;
  gsl_matrix *covar;
}levmar_internal_gsl_info;

/* These are PSRSALSA specific parameters that need to be passed on between functions */
typedef struct {
  double *covar;         // This will hold the covariance matrix (nrfitparams X nrfitparams)
  double *alpha;         // This will be the alpha matrix, which multiplied with the solution vector gives beta
  double *beta;          // The beta vector
  double *alpha_tmp;     // Storage space for a copy of alpha
  double *beta_tmp;      // Storage space for a copy of beta
  double *func_params_laststep;   // The parameters, as found in the last step
  double *derivatives;   // A vector containing the derivatives of the fit function to the different fit parameters
  double alambda;        // A key numerical factor used in the Levenberg-Marquart algorithm
  double chisq;          // The computer chi-square
  long nrfitparams;      // The number of fit parameters
  levmar_internal_fitter_data_def *data;
}levmar_internal_psrsalsa_info;

//START REGION DEVELOP

#ifdef NRAVAIL
/* These are NR specific parameters that need to be passed on between functions */
typedef struct {
  int *ia;
  double **covar;
  double **alpha;
  double *func_params_laststep;
  double alambda;
  double chisq;
}levmar_internal_nr_info;
#endif

//START REGION RELEASE
void print_gsl_version_used(FILE *stream)
{
  //  fprintf(stream, "%s", gsl_version);
  int major, minor;
  major = GSL_VERSION_NUMBER/100.0;
  minor = GSL_VERSION_NUMBER-100*major;
  fprintf(stream, "%s (header) %d.%d (specified during compilation)", GSL_VERSION, major, minor);
}

//START REGION DEVELOP
//START REGION RELEASE

/* Calculates collection of function at value x. */
double evaluate_fitfunc_collection(fitfunc_collection_type *function, double x, verbose_definition verbose)
{
  int i;
  double y;
  y = 0;
  for(i = 0; i < function->nrfuncs; i++) {
    if(function->func[i].type == FUNC_POLYNOMAL) {    /* y = a[1] * x ** p1 */
      y += function->func[i].value[0]*pow(x, function->func[i].param[0]);
//START REGION DEVELOP
    }else if(function->func[i].type == FUNC_ATAN) {   /* y = a[1]*atan(a[2]*(x+a[3])) */
      y += function->func[i].value[0]*atan(function->func[i].value[1]*(x + function->func[i].value[2]));
    }else if(function->func[i].type == FUNC_SIN) {   /* y = a[1]*sin(a[2]*(x+a[3])) */
      y += function->func[i].value[0]*sin(function->func[i].value[1]*(x + function->func[i].value[2]));
    }else if(function->func[i].type == FUNC_COS) {   /* y = a[1]*cos(a[2]*(x+a[3])) */
      y += function->func[i].value[0]*cos(function->func[i].value[1]*(x + function->func[i].value[2]));
    }else if(function->func[i].type == FUNC_NORMAL) {   /* y = a[1]*exp(-a[2]*(x-a[3])^2) */
      double tmp;
      tmp = x - function->func[i].value[2];
      y += function->func[i].value[0]*exp(-function->func[i].value[1]*tmp*tmp);
//START REGION RELEASE
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR evaluate_fitfunction_collection: Unknown funtional type in specified function.");
      exit(0);
    }
  }
  return y;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Calculates the derivative (to specified parameter number) of function functionnr and parameter paramnr at value x. */
double evaluate_fitfunc_collection_deriv_param(fitfunc_collection_type *function, int functionnr, int paramnr, double x, verbose_definition verbose)
{
//START REGION DEVELOP
  double tmp;
//START REGION RELEASE
  if(function->func[functionnr].type == FUNC_POLYNOMAL) {    /* y = a[1] * x ** p1 */
    if(paramnr == 0) {
      return pow(x, function->func[functionnr].param[0]);    /*  x ** p1 */
    }
//START REGION DEVELOP
  }else if(function->func[functionnr].type == FUNC_ATAN) {   /* y = a[1]*atan(a[2]*(x+a[3])) */
    if(paramnr == 0) {                                           // atan(a[2]*(x+a[3]))
      return atan(function->func[functionnr].value[1]*(x + function->func[functionnr].value[2]));
    }else if(paramnr == 1) {
      tmp = (x+function->func[functionnr].value[2]);
      return function->func[functionnr].value[0]*tmp/(1.0+function->func[functionnr].value[1]*function->func[functionnr].value[1]*tmp*tmp);
    }else if(paramnr == 2) {                                 // a[1]*a[2]/(1.0+(a[2]*(x+a[3])*a[2]*(x+a[3])))
      tmp = function->func[functionnr].value[1]*(x+function->func[functionnr].value[2]);
      return function->func[functionnr].value[0]*function->func[functionnr].value[1]/(1.0+tmp*tmp);
    }
  }else if(function->func[functionnr].type == FUNC_SIN) {   /* y = a[1]*sin(a[2]*(x+a[3])) */
    if(paramnr == 0) {                                           // sin(a[2]*(x+a[3]))
      return sin(function->func[functionnr].value[1]*(x+function->func[functionnr].value[2]));
    }else if(paramnr == 1) {                                  // a1*cos(a2*(x+a3))*(x+a3);
      tmp = x + function->func[functionnr].value[2];
      return function->func[functionnr].value[0]*cos(function->func[functionnr].value[1]*tmp)*tmp;
    }else if(paramnr == 2) {                                  // a1*cos(a2*(x+a3))*a2
      return function->func[functionnr].value[0]*cos(function->func[functionnr].value[1]*(x+function->func[functionnr].value[2]))*function->func[functionnr].value[1];
    }
  }else if(function->func[functionnr].type == FUNC_COS) {   /* y = a[1]*cos(a[2]*(x+a[3])) */
    if(paramnr == 0) {                                           // cos(a[2]*(x+a[3]))
      return cos(function->func[functionnr].value[1]*(x+function->func[functionnr].value[2]));
    }else if(paramnr == 1) {                                  // -a1*sin(a2*(x+a3))*(x+a3);
      tmp = x + function->func[functionnr].value[2];
      return -function->func[functionnr].value[0]*sin(function->func[functionnr].value[1]*tmp)*tmp;
    }else if(paramnr == 2) {                                  // -a1*sin(a2*(x+a3))*a2
      return -function->func[functionnr].value[0]*sin(function->func[functionnr].value[1]*(x+function->func[functionnr].value[2]))*function->func[functionnr].value[1];
    }
  }else if(function->func[functionnr].type == FUNC_NORMAL) {   /* y = a[1]*exp(-a[2]*(x-a[3])^2) */
    if(paramnr == 0) { 
      /* y = a1*exp(-a2*(x-a3)^2) */
      /* y = exp(-a2*(x-a3)^2) */
      tmp = x-function->func[functionnr].value[2];
      return exp(-function->func[functionnr].value[1]*tmp*tmp);
    }else if(paramnr == 1) { 
      /* y = a1*exp(-a2*(x-a3)^2) */
      /* y = a1*exp(-a2*(x-a3)^2)*(-(x-a3)^2) */
      tmp = x-function->func[functionnr].value[2];
      tmp *= tmp;
      return function->func[functionnr].value[0]*exp(-function->func[functionnr].value[1]*tmp)*(-tmp);
    }else if(paramnr == 2) { 
      /* y = a1*exp(-a2*(x-a3)^2) */
      /* y = a1*exp(-a2*(x-a3)^2)*(-a2)*2*(x-a3)*(-1) */
      /* y = 2*a1*a2*exp(-a2*(x-a3)^2)*(x-a3) */
    //    a1 = function->func[functionnr].value[0];
    //    a2 = function->func[functionnr].value[1];
    //    a3 = function->func[functionnr].value[2];
      tmp = x-function->func[functionnr].value[2];
      return 2.0*function->func[functionnr].value[0]*function->func[functionnr].value[1]*exp(-function->func[functionnr].value[1]*tmp*tmp)*tmp;
    }
//START REGION RELEASE
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR evaluate_fitfunc_collection_deriv_param: Unknown funtional type in specified function.");
    exit(0);
  }

  fflush(stdout);
  printerror(verbose.debug, "ERROR evaluate_fitfunc_collection_deriv_param: Unknown paramer number specified for given funtional type.");
  exit(0);
}

//START REGION DEVELOP

void show_fitfunctions_commandline_options(FILE *stream, char *cmd_flag, int nrspaces)
{
  int i;
  fprintf(stream, "%s", cmd_flag);
  for(i = 0; i < nrspaces-strlen(cmd_flag); i++) fprintf(stream, " ");
  fprintf(stream, "\"type parameters start_values fit_flags\"\n");
  for(i = 0; i < nrspaces; i++) fprintf(stream, " ");
  fprintf(stream, "if fit_flag is set, than this value is fitted for\n");
  for(i = 0; i < nrspaces; i++) fprintf(stream, " ");
  fprintf(stream, "type = %d: a1*x**p1\n", FUNC_POLYNOMAL); 
  for(i = 0; i < nrspaces; i++) fprintf(stream, " ");
  fprintf(stream, "type = %d: a1*atan(a2*(x+a3))\n", FUNC_ATAN); 
  for(i = 0; i < nrspaces; i++) fprintf(stream, " ");
  fprintf(stream, "type = %d: a1*sin(a2*(x+a3))\n", FUNC_SIN); 
  for(i = 0; i < nrspaces; i++) fprintf(stream, " ");
  fprintf(stream, "type = %d: a1*cos(a2*(x+a3))\n", FUNC_COS); 
  for(i = 0; i < nrspaces; i++) fprintf(stream, " ");
  fprintf(stream, "type = %d: a1*exp(-a2*(x-a3)**2)\n", FUNC_NORMAL); 
}

int parse_commandline_fitfunctions(int argc, char **argv, char *cmd_flag, fitfunc_collection_type *function, verbose_definition verbose)
{
  int i, j, k, nrparams, nrvalues;
  char txt[1000], *txt_ptr;
  function->nrfuncs = 0;
  for(i = 1; i < argc-1; i++) {
    if(strcmp(argv[i], "-f") == 0) {
      strcpy(txt, argv[i+1]);
      txt_ptr = strtok(txt, " ");
      if(txt_ptr == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR parse_commandline_fitfunctions: Cannot parse '%s' for first argument", txt);
	return 0;
      }
      k = sscanf(txt_ptr, "%d", &function->func[function->nrfuncs].type);
      if(k != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR parse_commandline_fitfunctions: Cannot parse '%s' for first argument", txt);
	return 0;
      }
      if(verbose.verbose)
	printf("t=%d ", function->func[function->nrfuncs].type);

      if(function->func[function->nrfuncs].type == FUNC_POLYNOMAL) {
	nrparams = 1;
	nrvalues = 1;
      }else if(function->func[function->nrfuncs].type == FUNC_ATAN) {
	nrparams = 0;
	nrvalues = 3;
      }else if(function->func[function->nrfuncs].type == FUNC_SIN) {
	nrparams = 0;
	nrvalues = 3;
      }else if(function->func[function->nrfuncs].type == FUNC_COS) {
	nrparams = 0;
	nrvalues = 3;
      }else if(function->func[function->nrfuncs].type == FUNC_NORMAL) {
	nrparams = 0;
	nrvalues = 3;
      }else {
	fflush(stdout);
	printerror(verbose.debug, "ERROR parse_commandline_fitfunctions: Unrecognized function type.");
	return 0;
      }

      for(j = 0; j < nrparams; j++) {
	txt_ptr = strtok(NULL, " ");
	if(txt_ptr == NULL) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR parse_commandline_fitfunctions: Cannot parse '%s' for %d arguments", argv[i+1], nrparams+2*nrvalues+1);
	  return 0;
	}
	if(verbose.verbose)
	  printf("p%d = %s ", j+1, txt_ptr);
	k = sscanf(txt_ptr, "%lf", &function->func[function->nrfuncs].param[j]);
	if(k != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR parse_commandline_fitfunctions: Cannot parse %s", argv[i+1]);
	  return 0;
	}
      }

      for(j = 0; j < nrvalues; j++) {
	txt_ptr = strtok(NULL, " ");
	if(txt_ptr == NULL) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR parse_commandline_fitfunctions: Cannot parse '%s' for %d arguments", argv[i+1], nrparams+2*nrvalues+1);
	  return 0;
	}
	if(verbose.verbose)
	  printf("a%d = %s ", j+1, txt_ptr);
	k = sscanf(txt_ptr, "%lf", &function->func[function->nrfuncs].start[j]);
	if(k != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR parse_commandline_fitfunctions: Cannot parse %s", argv[i+1]);
	  return 0;
	}
	function->func[function->nrfuncs].value[j] = function->func[function->nrfuncs].start[j];
      }

      for(j = 0; j < nrvalues; j++) {
	txt_ptr = strtok(NULL, " ");
	if(txt_ptr == NULL) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR parse_commandline_fitfunctions: Cannot parse '%s' for %d arguments", argv[i+1], nrparams+2*nrvalues+1);
	  return 0;
	}
	k = sscanf(txt_ptr, "%d", &function->func[function->nrfuncs].fit_flag[j]);
	if(function->func[function->nrfuncs].fit_flag[j])
	  function->func[function->nrfuncs].fit_flag[j] = 1;
	if(k != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR parse_commandline_fitfunctions: Cannot parse %s", argv[i+1]);
	  return 0;
	}
	if(verbose.verbose)
	  printf("f%d = %d ", j+1, function->func[function->nrfuncs].fit_flag[j]);
      }
      if(verbose.verbose)
	printf("\n");
      function->nrfuncs  += 1;
      i++;
    }
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* 
This function will show the specified function. If novalue is set,
the functional form is shown, but not the actual values of the
parameters. If showerror is set the error is shown as well as the
value.  A negative index means show all parameters, otherwise the
function index of the function collection is shown.  
*/
int print_fitfunctions(fitfunc_collection_type *function, int novalue, int showerror, int index, verbose_definition verbose)
{
  int i, j, n;

  if(showerror)
    showerror = 1;

  for(n = 0; n <= showerror; n++) {
    if(n == 0)
      printf("y = ");
    j = 1;
    for(i = 0; i < function->nrfuncs; i++) {
      if(index == i || index < 0) {
	if(function->func[i].type == FUNC_POLYNOMAL) {
	  if(novalue || showerror) {
	    if(showerror && n == 1) {
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j, function->func[i].value[0], function->func[i].error[0], function->func[i].value[0], function->func[i].error[0]);
	    }else {
	      printf("a%d*x**%lf", j, function->func[i].param[0]);
	    }
	  }else {
	    printf("%lf*x**%lf", function->func[i].value[0], function->func[i].param[0]);
	  }
	  j += 1;
//START REGION DEVELOP
	}else if(function->func[i].type == FUNC_ATAN) {
	  if(novalue || showerror) {
	    if(showerror && n == 1) {
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j,   function->func[i].value[0], function->func[i].error[0],   function->func[i].value[0], function->func[i].error[0]);
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j+1, function->func[i].value[1], function->func[i].error[1], function->func[i].value[1], function->func[i].error[1]);
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j+2, function->func[i].value[2], function->func[i].error[2], function->func[i].value[2], function->func[i].error[2]);
	    }else {
	      printf("a%d*atan(a%d*(x+a%d))", j, j+1, j+2);
	    }
	  }else {
	    printf("%lf*atan(%lf*(x+%lf))", function->func[i].value[0], function->func[i].value[1], function->func[i].value[2]);
	  }
	  j += 3;
	}else if(function->func[i].type == FUNC_SIN) {
	  if(novalue || showerror) {
	    if(showerror && n == 1) {
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j,   function->func[i].value[0], function->func[i].error[0],   function->func[i].value[0], function->func[i].error[0]);
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j+1, function->func[i].value[1], function->func[i].error[1], function->func[i].value[1], function->func[i].error[1]);
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j+2, function->func[i].value[2], function->func[i].error[2], function->func[i].value[2], function->func[i].error[2]);
	    }else {
	      printf("a%d*sin(a%d*(x+a%d))", j, j+1, j+2);
	    }
	  }else {
	    printf("%lf*sin(%lf*(x+%lf))", function->func[i].value[0], function->func[i].value[1], function->func[i].value[2]);
	  }
	  j += 3;
	}else if(function->func[i].type == FUNC_COS) {
	  if(novalue || showerror) {
	    if(showerror && n == 1) {
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j,   function->func[i].value[0], function->func[i].error[0],   function->func[i].value[0], function->func[i].error[0]);
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j+1, function->func[i].value[1], function->func[i].error[1], function->func[i].value[1], function->func[i].error[1]);
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j+2, function->func[i].value[2], function->func[i].error[2], function->func[i].value[2], function->func[i].error[2]);
	    }else {
	      printf("a%d*cos(a%d*(x+a%d))", j, j+1, j+2);
	    }
	  }else {
	    printf("%lf*cos(%lf*(x+%lf))", function->func[i].value[0], function->func[i].value[1], function->func[i].value[2]);
	  }
	  j += 3;
	}else if(function->func[i].type == FUNC_NORMAL) {
	  if(novalue || showerror) {
	    if(showerror && n == 1) {
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j,   function->func[i].value[0], function->func[i].error[0],   function->func[i].value[0], function->func[i].error[0]);
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j+1, function->func[i].value[1], function->func[i].error[1], function->func[i].value[1], function->func[i].error[1]);
	      printf("    a%d = %13e +- %13e   (%lf +- %lf)\n", j+2, function->func[i].value[2], function->func[i].error[2], function->func[i].value[2], function->func[i].error[2]);
	    }else {
	      printf("a%d*exp(-a%d*(x-a%d)**2)", j, j+1, j+2);
	    }
	  }else {
	    printf("%lf*exp(-%lf*(x-%lf)**2)", function->func[i].value[0], function->func[i].value[1], function->func[i].value[2]);
	  }
	  j += 3;
//START REGION RELEASE
	}else {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR print_fitfunctions: Unrecognized function type.");
	  return 0;
	}
	if(i == function->nrfuncs - 1)
	  printf("\n");
	else if(n == 0)
	  printf(" + ");
      }
    }
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* 
   function will set nrfitparameters
 */
void countnrparameters_fitfunction(fitfunc_collection_type *function, int *nrfitparameters)
{
  int n, fnr, nrvalues;
  *nrfitparameters = 0;
  for(fnr = 0; fnr < function->nrfuncs; fnr++) {
    if(function->func[fnr].type == FUNC_POLYNOMAL) {    /* y = a[1] * x ** p1 */
      nrvalues = 1;
//START REGION DEVELOP
    }else if(function->func[fnr].type == FUNC_ATAN) {   /* y = a[1]*atan(a[2]*(x+a[3])) */
      nrvalues = 3;
    }else if(function->func[fnr].type == FUNC_SIN) {   /* y = a1*sin(a2*(x+a3)) */
      nrvalues = 3;
    }else if(function->func[fnr].type == FUNC_COS) {   /* y = a1*cos(a2*(x+a3)) */
      nrvalues = 3;
    }else if(function->func[fnr].type == FUNC_NORMAL) {   /* y = a1*exp(-a2*(x-a3)^2) */
      nrvalues = 3;
//START REGION RELEASE
    }else {
      fflush(stdout);
      printerror(0, "ERROR countnrparameters_fitfunction: Unknown funtional type.");
      exit(0);
    }
    for(n = 0; n < nrvalues; n++) {
      if(function->func[fnr].fit_flag[n]) {
	*nrfitparameters += 1;
      }
    }
  }
}

//START REGION DEVELOP
//START REGION RELEASE

/*
 Given the array fitparameters, update the value parameters in
 function. Note that here it is assumed that the length of the
 fitparameters array is the length of parameters which have fit_flag
 != 0

 Returns 0 on error
*/
int set_fitted_parameters_fitfunc_collection(fitfunc_collection_type *function, double *fitparameters, verbose_definition verbose)
{
  int fnr, nrvalues, valnr, fitparamnr;
  fitparamnr = 0;
  for(fnr = 0; fnr < function->nrfuncs; fnr++) {
    if(function->func[fnr].type == FUNC_POLYNOMAL) {    /* y = a[1] * x ** p1 */
      nrvalues = 1;
//START REGION DEVELOP
    }else if(function->func[fnr].type == FUNC_ATAN) {   /* y = a[1]*atan(a[2]*(x+a[3])) */
      nrvalues = 3;
    }else if(function->func[fnr].type == FUNC_SIN) {   /* y = a1*sin(a2*(x+a3)) */
      nrvalues = 3;
    }else if(function->func[fnr].type == FUNC_COS) {   /* y = a1*cos(a2*(x+a3)) */
      nrvalues = 3;
    }else if(function->func[fnr].type == FUNC_NORMAL) {   /* y = a1*exp(-a2*(x-a3)^2) */
      nrvalues = 3;
//START REGION RELEASE
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR set_fitted_parameters_fitfunc_collection: Unknown funtional type.");
      return 0;
    }
    for(valnr = 0; valnr < nrvalues; valnr++) {
      if(function->func[fnr].fit_flag[valnr]) {
	//	fprintf(stderr, "XXXXXXX Reading index %d\n", fitparamnr);
	function->func[fnr].value[valnr] = fitparameters[fitparamnr++];
      }
    }
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/*
  Using the current values of the fitparams, fill the array f (with a
  length equal to the number of data points that are fitted) with
  (f(x[n])-y[n])/sigma[n]. This function is called by the gsl fitter
  routine.
*/
int levmar_itteration_calc_function_internal_gsl(const gsl_vector *fitparams, void *data, gsl_vector *f)
{
  int n;
  size_t npts   = ((levmar_internal_fitter_data_def *)data)->n;
  double *xdata = ((levmar_internal_fitter_data_def *)data)->x;
  double *ydata = ((levmar_internal_fitter_data_def *)data)->y;
  double *sigma = ((levmar_internal_fitter_data_def *)data)->sigma;
  fitfunc_collection_type *fitfunction_internal = ((levmar_internal_fitter_data_def *)data)->fitfunction;
  double value;
  verbose_definition noverbose;

  cleanVerboseState(&noverbose);

  if(fitparams->stride != 1) {
    printerror(0, "ERROR levmar_itteration_calc_function_internal_gsl: It is assumed the gsl vectors have a stride of 1.");
    exit(0);
  }
  // Make a new fitfunction with the fit parameters of the current itteration
  fitfunc_collection_type newfunction;
  memcpy(&newfunction, fitfunction_internal, sizeof(fitfunc_collection_type));
  if(set_fitted_parameters_fitfunc_collection(&newfunction, fitparams->data, noverbose) == 0) {
    printerror(0, "ERROR levmar_itteration_calc_function_internal_gsl: Copying new fit parameters failed.");
    exit(0);
  }

  for(n = 0; n < npts; n++) {
    value = evaluate_fitfunc_collection(&newfunction, xdata[n], noverbose);
    gsl_vector_set(f, n, (value - ydata[n])/sigma[n]);
  }
  return GSL_SUCCESS;
}

//START REGION DEVELOP
//START REGION RELEASE

/*
  Using the current values of the fitparams, fill the matrix J with
  the Jacobian. This function is called by the gsl fitter routine.  
*/
int levmar_itteration_calc_deriv_internal_gsl(const gsl_vector *fitparams, void *data, gsl_matrix *J)
{
  int fnr, j, n, nrvalues, paramnr;
  size_t npts   = ((levmar_internal_fitter_data_def *)data)->n;
  double *xdata = ((levmar_internal_fitter_data_def *)data)->x;
  double *sigma = ((levmar_internal_fitter_data_def *)data)->sigma;
  fitfunc_collection_type *fitfunction_internal = ((levmar_internal_fitter_data_def *)data)->fitfunction;
  double deriv;
  verbose_definition noverbose;

  cleanVerboseState(&noverbose);

  if(fitparams->stride != 1) {
    printerror(0, "ERROR levmar_itteration_calc_function_internal_gsl: It is assumed the gsl vectors have a stride of 1.");
    exit(0);
  }
  // Make a new fitfunction with the fit parameters of the current itteration
  fitfunc_collection_type newfunction;
  memcpy(&newfunction, fitfunction_internal, sizeof(fitfunc_collection_type));
  if(set_fitted_parameters_fitfunc_collection(&newfunction, fitparams->data, noverbose) == 0) {
    printerror(0, "ERROR levmar_itteration_calc_function_internal_gsl: Copying new fit parameters failed.");
    exit(0);
  }


  // Loop over all parameters
  j = 0;
  for(fnr = 0; fnr < newfunction.nrfuncs; fnr++) {
    if(newfunction.func[fnr].type == FUNC_POLYNOMAL) {    /* y = a[1] * x ** p1 */
      nrvalues = 1;
//START REGION DEVELOP
   }else if(newfunction.func[fnr].type == FUNC_ATAN) {   /* y = a[1]*atan(a[2]*(x+a[3])) */
      nrvalues = 3;
    }else if(newfunction.func[fnr].type == FUNC_SIN) {   /* y = a1*sin(a2*(x+a3)) */
      nrvalues = 3;
    }else if(newfunction.func[fnr].type == FUNC_COS) {   /* y = a1*cos(a2*(x+a3)) */
      nrvalues = 3;
    }else if(newfunction.func[fnr].type == FUNC_NORMAL) {   /* y = a1*exp(-a2*(x-a3)^2) */
      nrvalues = 3;
//START REGION RELEASE
    }else {
      fflush(stdout);
      printerror(noverbose.debug, "ERROR levmar_itteration_calc_deriv_internal_gsl: Unknown funtional type.");
      return 0;
    }
    for(paramnr = 0; paramnr < nrvalues; paramnr++) {
      if(newfunction.func[fnr].fit_flag[paramnr]) {
	for(n = 0; n < npts; n++) {
	  deriv = evaluate_fitfunc_collection_deriv_param(&newfunction, fnr, paramnr, xdata[n], noverbose);
	  gsl_matrix_set(J, n, j, deriv/sigma[n]);
	}
	j++;
      }
    }
  }
  return GSL_SUCCESS;
}

//START REGION DEVELOP
//START REGION RELEASE

/*
  Using the current values of the fitparams, fill the matrix J
(Jacobian) and find vector f with the difference between the function
and the data. This function is called by the gsl fitter routine.  */
int levmar_itteration_calc_func_and_deriv_internal_gsl(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J)
{
  levmar_itteration_calc_function_internal_gsl(x, data, f);
  levmar_itteration_calc_deriv_internal_gsl(x, data, J);
  return GSL_SUCCESS;
}

//START REGION DEVELOP

#ifdef NRAVAIL
/* Evaluates the fitting function
yfit, and its derivatives dyda[1..ma] with respect to the fitting parameters a at x.
*/
void levmar_itteration_calc_func_and_deriv_internal_nr(double x, double a[], double *yfit, double dyda[], int ma) 
{
  verbose_definition noverbose;

  cleanVerboseState(&noverbose);

  //  fprintf(stderr, "XXXX ma = %d\n", ma);
  //  fprintf(stderr, "XXXX a[1] = %lf\n", a[1]);

  // Make a new fitfunction with the fit parameters of the current itteration
  fitfunc_collection_type newfunction;
  memcpy(&newfunction, levmar_internal_fitter_data.fitfunction, sizeof(fitfunc_collection_type));
  if(set_fitted_parameters_fitfunc_collection(&newfunction, a+1, noverbose) == 0) {
    printerror(0, "ERROR levmar_itteration_calc_func_and_deriv_internal_nr: Copying new fit parameters failed.");
    exit(0);
  }

  *yfit = evaluate_fitfunc_collection(&newfunction, x, noverbose);

  int j, fnr, nrvalues, paramnr;
  double deriv;

  j = 0;
  for(fnr = 0; fnr < newfunction.nrfuncs; fnr++) {
    if(newfunction.func[fnr].type == FUNC_POLYNOMAL) {    /* y = a[1] * x ** p1 */
      nrvalues = 1;
    }else if(newfunction.func[fnr].type == FUNC_ATAN) {   /* y = a[1]*atan(a[2]*(x+a[3])) */
      nrvalues = 3;
    }else if(newfunction.func[fnr].type == FUNC_SIN) {   /* y = a1*sin(a2*(x+a3)) */
      nrvalues = 3;
    }else if(newfunction.func[fnr].type == FUNC_COS) {   /* y = a1*cos(a2*(x+a3)) */
      nrvalues = 3;
    }else if(newfunction.func[fnr].type == FUNC_NORMAL) {   /* y = a1*exp(-a2*(x-a3)^2) */
      nrvalues = 3;
    }else {
      fflush(stdout);
      printerror(noverbose.debug, "ERROR levmar_itteration_calc_func_and_deriv_internal_nr: Unknown funtional type.");
      exit(0);
    }
    for(paramnr = 0; paramnr < nrvalues; paramnr++) {
      if(newfunction.func[fnr].fit_flag[paramnr]) {
	deriv = evaluate_fitfunc_collection_deriv_param(&newfunction, fnr, paramnr, x, noverbose);
	dyda[j+1] = deriv;
	j++;
      }
    }
  }
}
#endif

//START REGION DEVELOP
//START REGION RELEASE
/* Evaluates the vector beta and the matrix alpha (NR eq. 15.5.8 on pg
   682 (really pg 706) for the current vector of solution parameters
   params.

   The vector beta can be written as:
   beta_k = sum over datapoints i{ (df/d param_k)*(y_i -f)/variance_i }

   The matrix alpha corresponds to the so-called Hessian
   matrix can be calculates as:

   alpha_jk = sum over datapoints i{ (df/d param_j)*(df/d param_k)/variance_i }


*/
void levmar_itteration_calc_alpha_beta_chi2_internal_psrsalsa(double *params, void *info) 
{
  long curDataPoint, nrDataPoints, curparam, curparam2, nrfitparams, fnr, nrvalues, paramnr;
  double *alpha, *beta, *derivatives, *x, *y, *sigma, chi2, ypred, tmp1, tmp2, tmp3, tmp4;
  verbose_definition noverbose;
  cleanVerboseState(&noverbose);

  // The pointer info points to a levmar_internal_psrsalsa_info
  // struct. We need to extract some information from this.
  nrDataPoints = ((levmar_internal_psrsalsa_info *)info)->data->n;
  x = ((levmar_internal_psrsalsa_info *)info)->data->x;
  y = ((levmar_internal_psrsalsa_info *)info)->data->y;
  sigma = ((levmar_internal_psrsalsa_info *)info)->data->sigma;
  alpha = ((levmar_internal_psrsalsa_info *)info)->alpha;
  beta = ((levmar_internal_psrsalsa_info *)info)->beta;
  derivatives = ((levmar_internal_psrsalsa_info *)info)->derivatives;
  nrfitparams = ((levmar_internal_psrsalsa_info *)info)->nrfitparams;

  //  for(curparam = 0; curparam < nrfitparams; curparam++) {
  //    printf("XXXXX param=%ld -> %lf\n", curparam, params[curparam]);
  //  }

  // Make a new fitfunction with the fit parameters of the current itteration
  fitfunc_collection_type newfunction;
  memcpy(&newfunction, levmar_internal_fitter_data.fitfunction, sizeof(fitfunc_collection_type));
  if(set_fitted_parameters_fitfunc_collection(&newfunction, params, noverbose) == 0) {
    printerror(0, "ERROR levmar_itteration_calc_func_and_deriv_internal_nr: Copying new fit parameters failed.");
    exit(0);
  }

  // Clear the beta vector and the alpha matrix
  // Probably a memset is quicker, as it requires less loops and calculations
  for(curparam = 0; curparam < nrfitparams; curparam++) {
    beta[curparam] = 0;
    // Only clear diagonal and the top-right corner of the alpha matrix. 
    // Since the matrix is symmetric, the bottom-right corner will be copied in later.
    for(curparam2 = curparam; curparam2 < nrfitparams; curparam2++) {
      alpha[curparam*nrfitparams+curparam2] = 0;
    }
  }
  chi2 = 0;

  for(curDataPoint = 0; curDataPoint < nrDataPoints; curDataPoint++) {
    // Get the predicted value of the fit-function at the location of the current data-point
    ypred = evaluate_fitfunc_collection(&newfunction, x[curDataPoint], noverbose);

    // Get the derivatives
    curparam = 0;
    for(fnr = 0; fnr < newfunction.nrfuncs; fnr++) {
      if(newfunction.func[fnr].type == FUNC_POLYNOMAL) {    /* y = a[1] * x ** p1 */
	nrvalues = 1;
//START REGION DEVELOP
      }else if(newfunction.func[fnr].type == FUNC_ATAN) {   /* y = a[1]*atan(a[2]*(x+a[3])) */
	nrvalues = 3;
      }else if(newfunction.func[fnr].type == FUNC_SIN) {   /* y = a1*sin(a2*(x+a3)) */
	nrvalues = 3;
      }else if(newfunction.func[fnr].type == FUNC_COS) {   /* y = a1*cos(a2*(x+a3)) */
	nrvalues = 3;
      }else if(newfunction.func[fnr].type == FUNC_NORMAL) {   /* y = a1*exp(-a2*(x-a3)^2) */
	nrvalues = 3;
//START REGION RELEASE
      }else {
	fflush(stdout);
	printerror(noverbose.debug, "ERROR levmar_itteration_calc_alpha_beta_chi2_internal_psrsalsa: Unknown funtional type.");
	exit(0);
      }
      for(paramnr = 0; paramnr < nrvalues; paramnr++) {
	if(newfunction.func[fnr].fit_flag[paramnr]) {
	  derivatives[curparam++] = evaluate_fitfunc_collection_deriv_param(&newfunction, fnr, paramnr, x[curDataPoint], noverbose);
	}
      }
    }

    tmp1 = 1.0/(sigma[curDataPoint]*sigma[curDataPoint]);
    tmp2 = (y[curDataPoint] - ypred);
    tmp3 = tmp1*tmp2;
    chi2 += tmp2*tmp3;
    // Fill the vectors beta and the matrix alpha
    for(curparam = 0; curparam < nrfitparams; curparam++) {
      tmp4 = derivatives[curparam]*tmp3;
      beta[curparam] += tmp4;
      for(curparam2 = curparam; curparam2 < nrfitparams; curparam2++) {
	alpha[curparam*nrfitparams+curparam2] += derivatives[curparam]*derivatives[curparam2]*tmp1;
      }
    }
  }

  // Fill in a copy of the upper right corner of alpha in the bottom-left corner to complete the matrix
  for(curparam = 0; curparam < nrfitparams; curparam++) {
    for(curparam2 = curparam+1; curparam2 < nrfitparams; curparam2++) {
      alpha[curparam2*nrfitparams+curparam] = alpha[curparam*nrfitparams+curparam2];
    }
  }

  ((levmar_internal_psrsalsa_info *)info)->chisq = chi2;
}

//START REGION DEVELOP
//START REGION RELEASE

// Returns 1 = ok, 0=error
int initialize_levmar_check_start_and_current_are_the_same(fitfunc_collection_type *function, verbose_definition verbose)
{
  int fnr, nrvalues, paramnr;
  for(fnr = 0; fnr < function->nrfuncs; fnr++) {
    if(function->func[fnr].type == FUNC_POLYNOMAL) {    /* y = a[1] * x ** p1 */
      nrvalues = 1;
//START REGION DEVELOP
    }else if(function->func[fnr].type == FUNC_ATAN) {   /* y = a[1]*atan(a[2]*(x+a[3])) */
      nrvalues = 3;
    }else if(function->func[fnr].type == FUNC_SIN) {   /* y = a1*sin(a2*(x+a3)) */
      nrvalues = 3;
    }else if(function->func[fnr].type == FUNC_COS) {   /* y = a1*cos(a2*(x+a3)) */
      nrvalues = 3;
    }else if(function->func[fnr].type == FUNC_NORMAL) {   /* y = a1*exp(-a2*(x-a3)^2) */
      nrvalues = 3;
//START REGION RELEASE
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR levmar_itteration_calc_func_and_deriv_internal_nr: Unknown funtional type.");
      return 0;
    }
    for(paramnr = 0; paramnr < nrvalues; paramnr++) {
      if(function->func[fnr].fit_flag[paramnr]) {
	if(function->func[fnr].value[paramnr] != function->func[fnr].start[paramnr]) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR initialize_levmar_check_start_and_current_are_the_same: Start and current value appear to be different (%lf != %lf)", function->func[fnr].value[paramnr], function->func[fnr].start[paramnr]);
	  return 0;
	}
      }
    }
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/*

mode = 0 = set nrfitparameters only
mode = 1 = do all allocations (func_params, func_params_err, solver, covar)
mode = 2 = free things up and store results in function

if showcovariance is set, the covariance matrix is reported before it is destroyed in mode 2

returns 0 if error
*/
int levmar_initialize(fitfunc_collection_type *function, double **func_params, double **func_params_err, int *nrfitparameters, double *datax, double *datay, double *datasigma, long nrdatapoints, levmar_internal_gsl_info *gsl_params, levmar_internal_psrsalsa_info *psrsalsa_params, 
#ifdef NRAVAIL
levmar_internal_nr_info *nr_params, 
#endif
int algorithm, int mode, int showcovariance, verbose_definition verbose)
{
  int fnr, nrvalues, paramnr, n, n2;
  if(mode == 0) {
    countnrparameters_fitfunction(function, nrfitparameters);
    if(verbose.verbose) printf("  Initializing fitter with %d fit parameters\n", *nrfitparameters);
  }else if(mode == 1) {
    // Some information the function value & derivative calculate
    // functions during the itterations need to know about.
    levmar_internal_fitter_data.n = nrdatapoints;
    levmar_internal_fitter_data.x = datax;
    levmar_internal_fitter_data.y = datay;
    levmar_internal_fitter_data.sigma = datasigma;
    levmar_internal_fitter_data.fitfunction = function;

    // Allocate memory that will fit the found parameters + errors
    *func_params = malloc(*nrfitparameters*sizeof(double));
    *func_params_err = malloc(*nrfitparameters*sizeof(double));
    if(*func_params == NULL || *func_params_err == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR levmar_initialize: Cannot allocate memory");
      return 0;
    }

    // Now fill the func_params struct with the initial values
    paramnr = 0;
    for(fnr = 0; fnr < function->nrfuncs; fnr++) {
      if(function->func[fnr].type == FUNC_POLYNOMAL) {    /* y = a[1] * x ** p1 */
	nrvalues = 1;
//START REGION DEVELOP
      }else if(function->func[fnr].type == FUNC_ATAN) {   /* y = a[1]*atan(a[2]*(x+a[3])) */
	nrvalues = 3;
      }else if(function->func[fnr].type == FUNC_SIN) {   /* y = a1*sin(a2*(x+a3)) */
	nrvalues = 3;
      }else if(function->func[fnr].type == FUNC_COS) {   /* y = a1*cos(a2*(x+a3)) */
	nrvalues = 3;
      }else if(function->func[fnr].type == FUNC_NORMAL) {   /* y = a1*exp(-a2*(x-a3)^2) */
	nrvalues = 3;
//START REGION RELEASE
      }else {
	fflush(stdout);
	printerror(verbose.debug, "ERROR levmar_initialize: Unknown funtional type.");
	return 0;
      }
      for(n = 0; n < nrvalues; n++) {
	if(function->func[fnr].fit_flag[n]) {
	  (*func_params)[paramnr++] = function->func[fnr].value[n];
	}
	function->func[fnr].error[n] = sqrt(-1);    // Set all errors to undefined
      }
    }

    //    int param;
    //    fprintf(stderr, "XXXX Mode = %d\n", mode);
    //    fprintf(stderr, "XXXX Should have set %d initial variables, expected %d\n", j, *nrfitparameters);
    //    for(param = 0; param < *nrfitparameters; param++) {
    //      fprintf(stderr, "XXXX param %d = %lf\n", param+1, (*func_params)[param]);
    //    }
    if(algorithm == 1) {  // Initialize some GSL specific stuff
      gsl_params->func_params = gsl_vector_view_array(*func_params, *nrfitparameters);

      gsl_params->covar = gsl_matrix_alloc(*nrfitparameters, *nrfitparameters);
      gsl_params->user_itteration_funcs.f = &levmar_itteration_calc_function_internal_gsl;
      gsl_params->user_itteration_funcs.df = &levmar_itteration_calc_deriv_internal_gsl;
      gsl_params->user_itteration_funcs.fdf = &levmar_itteration_calc_func_and_deriv_internal_gsl;
      gsl_params->user_itteration_funcs.n = nrdatapoints;
      gsl_params->user_itteration_funcs.p = *nrfitparameters;
      gsl_params->user_itteration_funcs.params = &levmar_internal_fitter_data;

      gsl_params->solver_type = gsl_multifit_fdfsolver_lmsder;
      gsl_params->solver = gsl_multifit_fdfsolver_alloc (gsl_params->solver_type, nrdatapoints, *nrfitparameters);
      gsl_multifit_fdfsolver_set(gsl_params->solver, &(gsl_params->user_itteration_funcs), &(gsl_params->func_params.vector));
    }else if(algorithm == 2) {
      psrsalsa_params->func_params_laststep = malloc(sizeof(double)*(*nrfitparameters));
      psrsalsa_params->beta = malloc(sizeof(double)*(*nrfitparameters));
      psrsalsa_params->derivatives = malloc(sizeof(double)*(*nrfitparameters));
      psrsalsa_params->covar = malloc(sizeof(double)*(*nrfitparameters)*(*nrfitparameters));
      psrsalsa_params->alpha = malloc(sizeof(double)*(*nrfitparameters)*(*nrfitparameters));
      psrsalsa_params->alpha_tmp = malloc(sizeof(double)*(*nrfitparameters)*(*nrfitparameters));
      psrsalsa_params->beta_tmp = malloc(sizeof(double)*(*nrfitparameters));
      if(psrsalsa_params->func_params_laststep == NULL || psrsalsa_params->beta == NULL || psrsalsa_params->derivatives == NULL || psrsalsa_params->alpha == NULL || psrsalsa_params->covar == NULL || psrsalsa_params->beta_tmp == NULL || psrsalsa_params->alpha_tmp == NULL) {
	printerror(verbose.debug, "ERROR levmar_initialize: Memory allocation error");
	return 0;
      }
      psrsalsa_params->nrfitparams = *nrfitparameters;
      psrsalsa_params->data = &levmar_internal_fitter_data;
      psrsalsa_params->alambda = 0.001;    // Initial value of lambda to use
      // Obtain initial alpha, beta and chi2 values
      levmar_itteration_calc_alpha_beta_chi2_internal_psrsalsa(*func_params, psrsalsa_params);
    }else if(algorithm == 3) {
#ifdef NRAVAIL
      nr_params->ia = malloc(sizeof(int)*(*nrfitparameters));
      nr_params->func_params_laststep = malloc(sizeof(double)*(*nrfitparameters));
      if(nr_params->ia == NULL || nr_params->func_params_laststep == NULL) {
	printerror(verbose.debug, "ERROR levmar_initialize: Memory allocation error");
	return 0;
      }
      for(n = 0; n < *nrfitparameters; n++)
	nr_params->ia[n] = 1;  // All parameters should be fitted for, per definition, as non-fitted parameters has already been dealt with. 
      double **matrix_nr_d(long nrl, long nrh, long ncl, long nch);
      nr_params->covar = matrix_nr_d(1, *nrfitparameters, 1, *nrfitparameters);
      nr_params->alpha = matrix_nr_d(1, *nrfitparameters, 1, *nrfitparameters);
      void mrqmin_nr_d(double x[], double y[], double sig[], int ndata, double a[], int ia[],
		   int ma, double **covar, double **alpha, double *chisq,
		   void (*funcs)(double, double [], double *, double [], int), double *alamda);
      // Initialise the NR routine by setting alambda to a negative value. It will changed to a positive value so it can be used for the actual itterations
      nr_params->alambda = -1;
      mrqmin_nr_d(levmar_internal_fitter_data.x-1, levmar_internal_fitter_data.y-1, levmar_internal_fitter_data.sigma-1, nrdatapoints, (*func_params)-1, nr_params->ia-1, *nrfitparameters, nr_params->covar, nr_params->alpha, &nr_params->chisq, levmar_itteration_calc_func_and_deriv_internal_nr, &(nr_params->alambda));
#else
      printerror(verbose.debug, "ERROR levmar_initialize: Code is not compiled with NR support");
      return 0;
#endif
    }
  }else if(mode == 2) {

    // Copy found parameters in the function struct
    paramnr = 0;
    for(fnr = 0; fnr < function->nrfuncs; fnr++) {
      if(function->func[fnr].type == FUNC_POLYNOMAL) {    /* y = a[1] * x ** p1 */
	nrvalues = 1;
//START REGION DEVELOP
      }else if(function->func[fnr].type == FUNC_ATAN) {   /* y = a[1]*atan(a[2]*(x+a[3])) */
	nrvalues = 3;
      }else if(function->func[fnr].type == FUNC_SIN) {   /* y = a1*sin(a2*(x+a3)) */
	nrvalues = 3;
      }else if(function->func[fnr].type == FUNC_COS) {   /* y = a1*cos(a2*(x+a3)) */
	nrvalues = 3;
      }else if(function->func[fnr].type == FUNC_NORMAL) {   /* y = a1*exp(-a2*(x-a3)^2) */
	nrvalues = 3;
//START REGION RELEASE
      }else {
	fflush(stdout);
	printerror(0, "ERROR levmar_initialize: Unknown funtional type.");
	return 0;
      }
      for(n = 0; n < nrvalues; n++) {
	if(function->func[fnr].fit_flag[n]) {
	  function->func[fnr].value[n] = (*func_params)[paramnr];
	  function->func[fnr].error[n] = (*func_params_err)[paramnr];
	  paramnr++;
	}
      }
    }

    if(algorithm == 1) {  // Free some GSL specific stuff
      if(showcovariance) {
	printf("Covariance matrix: \n");
	for(n = 0; n < *nrfitparameters; n++) {
	  for(n2 = 0; n2 < *nrfitparameters; n2++) {
	    printf("%lf ", gsl_matrix_get(gsl_params->covar,n,n2));
	  }
	  printf("\n");
	}
      }
      gsl_multifit_fdfsolver_free(gsl_params->solver);
      gsl_matrix_free(gsl_params->covar);
    }else if(algorithm == 2) {  // Free some psrsalsa specific stuff
      if(showcovariance) {
	printf("Covariance matrix: \n");
	for(n = 0; n < *nrfitparameters; n++) {
	  for(n2 = 0; n2 < *nrfitparameters; n2++) {
	    printf("%lf ", psrsalsa_params->covar[n*(*nrfitparameters)+n2]);
	  }
	  printf("\n");
	}
      }
      free(psrsalsa_params->beta);
      free(psrsalsa_params->derivatives);
      free(psrsalsa_params->alpha);
      free(psrsalsa_params->covar);
      free(psrsalsa_params->func_params_laststep);
      free(psrsalsa_params->alpha_tmp);
      free(psrsalsa_params->beta_tmp);
    }else if(algorithm == 3) {  // Free some NR specific stuff
#ifdef NRAVAIL
      if(showcovariance) {
	printf("Covariance matrix: \n");
	for(n = 0; n < *nrfitparameters; n++) {
	  for(n2 = 0; n2 < *nrfitparameters; n2++) {
	    printf("%lf ", nr_params->covar[n+1][n2+1]);
	  }
	  printf("\n");
	}
      }
      void free_matrix_nr_d(double **m, long nrl, long nrh, long ncl, long nch);
      free_matrix_nr_d(nr_params->covar, 1, *nrfitparameters, 1, *nrfitparameters);
      free_matrix_nr_d(nr_params->alpha, 1, *nrfitparameters, 1, *nrfitparameters);
      free(nr_params->ia);
      free(nr_params->func_params_laststep);
#else
      printerror(verbose.debug, "ERROR levmar_initialize: Code is not compiled with NR support");
      return 0;
#endif
    }
    free(*func_params);
    free(*func_params_err);
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE


/*
Return values:
  0 = Converged
  2 = Max number of itterations reached
  3 = Machine precision limit reached
  4 = Cannot determine suitable trial step.
  9 = Critical error during execution, such as a memory allocation error
 10 = Shouldn't happen, not documented error code of gsl

algorithm = 1 = GSL  (gsl_solver and gsl_covar are only necessary/meaningful for this type. Set to NULL otherwise)
algorithm = 2 = psrsalsa version
algorithm = 3 = NR double-version hack
*/
int fit_levmar_internal(int algorithm, int nrdatapoints, int nrfitparameters, double *func_params, double *func_params_err, levmar_internal_gsl_info *gsl_params, levmar_internal_psrsalsa_info *psrsalsa_params, 
#ifdef NRAVAIL
			levmar_internal_nr_info *nr_params, 
#endif
			double epsabs, double epsrel, int maxiter, verbose_definition verbose)
{
  int status, i, converged;
  long iter;
  double chisq_last_nr;
  
#ifdef NRAVAIL
  void mrqmin_nr_d(double x[], double y[], double sig[], int ndata, double a[], int ia[],
		   int ma, double **covar, double **alpha, double *chisq,
		   void (*funcs)(double, double [], double *, double [], int), double *alamda);
#else
  if(algorithm == 3) {
    printerror(verbose.debug, "ERROR fit_levmar_internal: Code is not compiled with NR support");
    return 0;
  }
#endif

  // Initialisation of some things
  if(algorithm == 2) {
    // Make a copy of the initial parameters
    memcpy(psrsalsa_params->func_params_laststep, func_params, sizeof(double)*nrfitparameters);
    // Remember the initial chi square
    chisq_last_nr = psrsalsa_params->chisq;
#ifdef NRAVAIL
  }else if(algorithm == 3) {
    // Make a copy of the initial parameters
    memcpy(nr_params->func_params_laststep, func_params, sizeof(double)*nrfitparameters);
    // Remember the initial chi square
    chisq_last_nr = nr_params->chisq;
#endif
  }

  iter = 0;
  do {
    iter++;  // Number of iterations in this current call to this function
    levmar_internal_fitter_data.iter += 1;    // This is the total number of itterations, including multiple calls to this function
    if(algorithm == 1) {
      // Do an itteration
      status = gsl_multifit_fdfsolver_iterate(gsl_params->solver);
      if(status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG) { // Machine precision reached
	int status2;
	status2 = gsl_multifit_test_delta(gsl_params->solver->dx, gsl_params->solver->x, epsabs, epsrel);
	if(status2 == GSL_SUCCESS) {
	  status = 0;                // If we reached tollerance, clear the error that we reached the machine precision.
	}else {
	  if(verbose.debug == 0) {
	    printwarning(verbose.debug, "Machine precision reached, but required precision tolerance conditions are not met.");
	  }else {
	    printwarning(verbose.debug, "Machine precision reached, but required precision tolerance conditions are not met:");
	    for(i = 0; i < nrfitparameters; i++) {
	      printwarning(verbose.debug, "parameter %d: value=%e stepsize=%e", i+1, gsl_vector_get(gsl_params->solver->x, i), gsl_vector_get(gsl_params->solver->dx, i));
	    }
	  }
	}
      }
      if(verbose.debug) {
	printf("fit_levmar_internal: itteration %ld - status = %s\n", iter, gsl_strerror(status));
      }
      if(status) {  // Fit failed, so quit further processing
	break;
      }else {
	// Find out if convergence was reached
	status = gsl_multifit_test_delta(gsl_params->solver->dx, gsl_params->solver->x, epsabs, epsrel);
      }
    }else if(algorithm == 2) {
      // Make some copies in case the fit did not improve after next itteration
      for(i = 0; i < nrfitparameters*nrfitparameters; i++)
	psrsalsa_params->alpha_tmp[i] = psrsalsa_params->alpha[i];
      for(i = 0; i < nrfitparameters; i++)
	psrsalsa_params->beta_tmp[i] = psrsalsa_params->beta[i];

      // A key trick of the LevMar algorithm: multiply elements on the diagonal of alpha matrix with a numerical factor
      for (i = 0; i < nrfitparameters; i++) {
	psrsalsa_params->alpha[i*nrfitparameters+i] *= 1.0+psrsalsa_params->alambda;
      }
      // Think I should use LU decomp for this. It is 3X faster when the inverse is not required, and about the same to calculate the inverse. See pg 72 and earlier
      // alpha is replaced by its inverse, da (beta) is replaced by the solution vector
      int linalg_solve_matrix_eq_gauss_jordan(double *matrixa, double *matrixb, int n, int m, verbose_definition verbose);
      if(linalg_solve_matrix_eq_gauss_jordan(psrsalsa_params->alpha, psrsalsa_params->beta, nrfitparameters, 1, verbose) != 0) {
	printerror(verbose.debug, "ERROR fit_levmar_internal: Solving matrix equation failed.");
	return 9;
      }
      // Update solution vector
      for (i = 0; i < nrfitparameters; i++)
	func_params[i] += psrsalsa_params->beta[i];

      // Get new alpha and beta and chisq
      levmar_itteration_calc_alpha_beta_chi2_internal_psrsalsa(func_params, psrsalsa_params);
      // Need to check if new solution is acceptable

      status = GSL_SUCCESS;
      if(chisq_last_nr < psrsalsa_params->chisq) { // Chi square did not improve, so presumably at machine precission
	status = GSL_ETOLF;
	// Restore the previous solution
	memcpy(func_params, psrsalsa_params->func_params_laststep, sizeof(double)*nrfitparameters);
	for(i = 0; i < nrfitparameters*nrfitparameters; i++)
	  psrsalsa_params->alpha[i] = psrsalsa_params->alpha_tmp[i];
	for(i = 0; i < nrfitparameters; i++)
	  psrsalsa_params->beta[i] = psrsalsa_params->beta_tmp[i];
	psrsalsa_params->alambda *= 10.1;   // Increase the value of alambda. Not by the same factor as the decrease to possibly avoid ping-ponging back-forth between solutions. Maybe this is a non-issue.
      }

      converged = 1;
      // Do a convergence test
      for(i = 0; i < nrfitparameters; i++) {
	//	  fprintf(stderr, "XXXXXXX %lf %lf\n", func_params[i], psrsalsa_params->func_params_laststep[i]);
	if(fabs(func_params[i] - psrsalsa_params->func_params_laststep[i]) >= 1.2e-16+epsabs+epsrel*fabs(psrsalsa_params->func_params_laststep[i])) {
	  converged = 0;
	  break;
	}
      }
      if(converged && status == GSL_ETOLF)  // Clear machine precision error if converged
	status = GSL_SUCCESS;
      if(status == GSL_ETOLF) {
      	printwarning(verbose.debug, "Machine precision reached, but required precision tolerance conditions are not met.");
      }else if(converged == 0) {
	status = GSL_CONTINUE;
      }

      if(chisq_last_nr > psrsalsa_params->chisq) { // Chi square did improve
	//	fprintf(stderr, "fit_levmar_internal: itteration %ld - status = %s\n", iter, gsl_strerror(status));
	// Keep a copy of the last parameters to allow future checks of convergence
	memcpy(psrsalsa_params->func_params_laststep, func_params, sizeof(double)*nrfitparameters);
	chisq_last_nr = psrsalsa_params->chisq;
	psrsalsa_params->alambda *= 0.1;   // Decrease the alambda parameter if step was successful
      }
      if(verbose.debug) {
	printf("fit_levmar_internal: itteration %ld - status = %s\n", iter, gsl_strerror(status));
      }
#ifdef NRAVAIL
    }else if(algorithm == 3) {
      // Do another itteration
      mrqmin_nr_d(levmar_internal_fitter_data.x-1, levmar_internal_fitter_data.y-1, levmar_internal_fitter_data.sigma-1, nrdatapoints, func_params-1, nr_params->ia-1, nrfitparameters, nr_params->covar, nr_params->alpha, &(nr_params->chisq), levmar_itteration_calc_func_and_deriv_internal_nr, &(nr_params->alambda));
      // Do a convergence test
      status = GSL_SUCCESS;
      for(i = 0; i < nrfitparameters; i++) {
	if(fabs(func_params[i] - nr_params->func_params_laststep[i]) >= epsabs+epsrel*fabs(nr_params->func_params_laststep[i])) {
	  status = GSL_CONTINUE;
	  break;
	}
      }
      if(verbose.debug) {
	printf("fit_levmar_internal: itteration %ld - status = %s\n", iter, gsl_strerror(status));
      }
      if(status == GSL_CONTINUE) {
	if(chisq_last_nr <= nr_params->chisq) { // Chi square did not improve, so presumably at machine precission
	  status = GSL_ETOLF;
	  printwarning(verbose.debug, "Machine precision reached, but required precision tolerance conditions are not met.");
	  // Note that the solution is not changed in the last itteration because of a check in mrqmin_nr_d
	}else {
	  // Keep a copy of the last parameters to allow further check of convergence
	  memcpy(nr_params->func_params_laststep, func_params, sizeof(double)*nrfitparameters);
	  chisq_last_nr = nr_params->chisq;
	}
      }
#endif
    }
  }while(status == GSL_CONTINUE && iter < maxiter);


  // Done with fitting. Now find errors
  if(algorithm == 1) {
#if GSL_VERSION_NUMBER >= 200
    if(verbose.debug) {
      printf("fit_levmar_internal: Entering gsl >= 2.0 specific part of the code\n");
    }
    gsl_matrix *J = gsl_matrix_alloc(nrdatapoints, nrfitparameters);
    // This might not work in version 2.2 anymore. See pg 504 and compare with 489 in manual of version 2.1
    gsl_multifit_fdfsolver_jac(gsl_params->solver, J);
    gsl_multifit_covar(J, 0.0, gsl_params->covar);
#else
    if(verbose.debug) {
      printf("fit_levmar_internal: Entering gsl < 2.0 specific part of the code\n");
    }
    gsl_multifit_covar(gsl_params->solver->J, 0.0, gsl_params->covar);
#endif
    levmar_internal_fitter_data.fitfunction->chi2 = gsl_blas_dnrm2(gsl_params->solver->f);
    levmar_internal_fitter_data.fitfunction->chi2_red = (levmar_internal_fitter_data.fitfunction->chi2)/sqrt(nrdatapoints - nrfitparameters);
    // Not sure what this was supposed to be doing, as it is not used
    // double c;
    //    c = GSL_MAX_DBL(1, levmar_internal_fitter_data.fitfunction->chi2_red);
    levmar_internal_fitter_data.fitfunction->chi2_red *= levmar_internal_fitter_data.fitfunction->chi2_red;
    levmar_internal_fitter_data.fitfunction->chi2 *= levmar_internal_fitter_data.fitfunction->chi2;
  }else if(algorithm == 2) {
    // Obtain the covariance matrix by inverting the current alpha matrix
    int linalg_solve_matrix_eq_gauss_jordan(double *matrixa, double *matrixb, int n, int m, verbose_definition verbose);
    if(linalg_solve_matrix_eq_gauss_jordan(psrsalsa_params->alpha, psrsalsa_params->beta, nrfitparameters, 1, verbose) != 0) {
      printerror(verbose.debug, "ERROR fit_levmar_internal: Solving matrix equation failed.");
      return 9;
    }
    // Copy the inverse in the covar matrix
    for(i = 0; i < nrfitparameters*nrfitparameters; i++)
      psrsalsa_params->covar[i] = psrsalsa_params->alpha[i];

    //    psrsalsa_params->alambda = 0;
    //    mrqmin_psrsalsa_d(levmar_internal_fitter_data.x, levmar_internal_fitter_data.y, levmar_internal_fitter_data.sigma, nrdatapoints, func_params, nrfitparameters, psrsalsa_params->covar, psrsalsa_params->alpha, psrsalsa_params->beta-1, &psrsalsa_params->chisq, levmar_itteration_calc_alpha_beta_chi2_internal_psrsalsa, psrsalsa_params, &psrsalsa_params->alambda, verbose);

    levmar_internal_fitter_data.fitfunction->chi2 = psrsalsa_params->chisq;
    levmar_internal_fitter_data.fitfunction->chi2_red = psrsalsa_params->chisq/(double)(nrdatapoints - nrfitparameters);
#ifdef NRAVAIL
  }else if(algorithm == 3) {
    // Do a final call with alambda = 0 to obtain the covariance matrix
    nr_params->alambda = 0;
    mrqmin_nr_d(levmar_internal_fitter_data.x-1, levmar_internal_fitter_data.y-1, levmar_internal_fitter_data.sigma-1, nrdatapoints, func_params-1, nr_params->ia-1, nrfitparameters, nr_params->covar, nr_params->alpha, &nr_params->chisq, levmar_itteration_calc_func_and_deriv_internal_nr, &nr_params->alambda);

    levmar_internal_fitter_data.fitfunction->chi2 = nr_params->chisq;
    levmar_internal_fitter_data.fitfunction->chi2_red = nr_params->chisq/(double)(nrdatapoints - nrfitparameters);
#endif
  }
  if(verbose.debug) {
    printf("chisq/dof = %g after %ld iterations\n", levmar_internal_fitter_data.fitfunction->chi2_red, iter);
  }
  if(algorithm == 1) {
    for(i = 0; i < nrfitparameters; i++) {
      func_params[i] = gsl_vector_get(gsl_params->solver->x, i);
      //    func_params_err[i] = c*sqrt(gsl_matrix_get(gsl_params->covar,i,i));
      func_params_err[i] = sqrt(gsl_matrix_get(gsl_params->covar,i,i));
    }
  }else if(algorithm == 2) {
    for(i = 0; i < nrfitparameters; i++) {
      func_params_err[i] = sqrt(psrsalsa_params->covar[i*nrfitparameters+i]);
    }
#ifdef NRAVAIL
  }else if(algorithm == 3) {
    for(i = 0; i < nrfitparameters; i++) {
      func_params_err[i] = sqrt(nr_params->covar[i+1][i+1]);
    }
#endif
  }
  if(levmar_internal_fitter_data.iter == maxiter)
    return 2;
  if(algorithm == 1) {
    if(status == 0)
      return 0;
    if(status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG) {
      return 3;
    }
    if(status == GSL_ENOPROG) {
      fflush(stdout);
      printerror(verbose.debug, "fit_levmar_internal: Cannot determine suitable trial step. Continue itterating might help");
      return 4;
    }
  }else if(algorithm == 2 || algorithm == 3) {
    return 0;
  }
  fflush(stdout);
  printerror(verbose.debug, "fit_levmar_internal: UNKNOWN RETURN CODE");
  return 10;
}


//START REGION DEVELOP
//START REGION RELEASE

/*

This function is similar to fit_levmar2(), so make changes in both if
required.

Use the Levenberg-Marquardt algorithm to fit a function specified by
function to ndata datapoints (set data_sigma to NULL to assume equal
errors, which implies force_chi2_1). If oneatatime is speciefied one
parameter is fit at a time in a loop. This makes the fitting slower,
but potentially more robust. If force_chi2_1 is set, the errorbars of
the datapoints are rescaled such that the reduced chi2 is equal to
unity and the errorbars on the fit parameters are probably more
sensible. You can set either epsabs or epsrel (set one to zero),
otherwise the stringents condition of the two will be used. If the
step size |dx_i| is smaller than epsabs or epsrel|x|, then the
function is taken to be converged. maxiter sets the maximum number of
iterations. If showresults is set, or if verbose is set, the fitted
function + chi2 is reported. If showcovariance is set, the covariance
matrix is reported as well. The return value is the same value as
returned via status.

Return values:
  0 = Converged
  1 = Memory allocation error or problem with the initial parameters
  2 = Max number of itterations reached
  3 = Machine precision limit reached
  4 = Cannot determine suitable trial step.
 10 = Shouldn't happen, not documented error code of gsl

algorithm = 1 = GSL
algorithm = 2 = PSRSALSA version
algorithm = 3 = NR double-version hack

 */
int fit_levmar(int algorithm, fitfunc_collection_type *function, double *data_x, double *data_y, double *data_sigma, long ndata, int oneatatime, int force_chi2_1, double epsabs, double epsrel, int maxiter, int *status, int showresults, int showcovariance, verbose_definition verbose)
{
  int loopnr, improved, i, j, k, fnr, valuenr, nrfitparameters_one, nrfitparameters, free_data_sigma, nrvalues;
  double old_chi2, sig_scale, *func_params, *func_params_err;
  verbose_definition noverbose;
  levmar_internal_gsl_info gsl_params;
  levmar_internal_psrsalsa_info psrsalsa_params;

#ifndef NRAVAIL
  if(algorithm == 3) {
    printerror(verbose.debug, "ERROR fit_levmar: Code is not compiled with NR support");
    return 0;
  }
#else
  levmar_internal_nr_info nr_params;
#endif


  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;

  if(verbose.verbose) {
    printf("Levenberg-Marquardt algorithm:\n");
    if(verbose.debug) {
      printf("  Fit function: ");
      print_fitfunctions(function, 1, 0, -1, verbose);
    }
  }

  // Make sure the start and current variable values are set the same. I think this turns out to be required somewhere in the following code.
  if(initialize_levmar_check_start_and_current_are_the_same(function, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fit_levmar: Start and current value appear to be different");
    return 1;
  }

  // If no errors are assigned, take equal error bars
  free_data_sigma = 0;
  if(data_sigma == NULL) {
    data_sigma = malloc(ndata*sizeof(double));
    for(i = 0; i < ndata; i++)
      data_sigma[i] = 1;
    free_data_sigma = 1;
    force_chi2_1 = 1;
  }

  // This initialises nrfitparameters
  if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
		       &nr_params, 
#endif
		       algorithm, 0, showcovariance, verbose) == 0) {
    return 1;
  }

  levmar_internal_fitter_data.iter = 0;    // Keep track of the number of itterations performed
  if(oneatatime == 0) {
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
			 &nr_params, 
#endif
			 algorithm, 1, showcovariance, verbose) == 0) {
      return 1;
    }
    *status = fit_levmar_internal(algorithm, ndata, nrfitparameters, func_params, func_params_err, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
				  &nr_params, 
#endif
				  epsabs, epsrel, maxiter, verbose);
  }else {
    // Keep a copy of the original function definition, so we can fiddle with the settings of fit parameters to fit one parameter at the time
    fitfunc_collection_type *original_function;
    original_function = malloc(sizeof(fitfunc_collection_type));
    if(original_function == NULL) {
      printerror(verbose.debug, "ERROR fit_levmar: Memory allocation error.");
      return 1;
    }
    memcpy(original_function, function, sizeof(fitfunc_collection_type));
    old_chi2 = -1;
    loopnr = 0;
    do {
      improved = 0;
      i = 0;
      for(fnr = 0; fnr < function->nrfuncs; fnr++) {
	if(function->func[fnr].type == FUNC_POLYNOMAL) {    /* y = a[1] * x ** p1 */
	  nrvalues = 1;
//START REGION DEVELOP
	}else if(function->func[fnr].type == FUNC_ATAN) {   /* y = a[1]*atan(a[2]*(x+a[3])) */
	  nrvalues = 3;
	}else if(function->func[fnr].type == FUNC_SIN) {   /* y = a1*sin(a2*(x+a3)) */
	  nrvalues = 3;
	}else if(function->func[fnr].type == FUNC_COS) {   /* y = a1*cos(a2*(x+a3)) */
	  nrvalues = 3;
	}else if(function->func[fnr].type == FUNC_NORMAL) {   /* y = a1*exp(-a2*(x-a3)^2) */
	  nrvalues = 3;
//START REGION RELEASE
	}else {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR fit_levmar: Unknown funtional type.");
	  return 1;
	}
	for(valuenr = 0; valuenr < nrvalues; valuenr++) {
	  if(original_function->func[fnr].fit_flag[valuenr]) {
	    // Reset all fit flags, except parameter i
	    for(j = 0; j < function->nrfuncs; j++) {
	      for(k = 0; k < MaxNrFitParameters; k++)
		function->func[j].fit_flag[k] = 0;
	    }
	    function->func[fnr].fit_flag[valuenr] = 1;
	    if(verbose.debug) {
	      printf("  Loop %d: chi2 = %f", loopnr, old_chi2);
	      printf("\n");
	    }
	    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters_one, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
				 &nr_params, 
#endif
				 algorithm, 0, showcovariance, noverbose) == 0) {
	      free(original_function);
	      return 1;
	    }
	    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters_one, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
				 &nr_params, 
#endif
				 algorithm, 1, showcovariance, verbose) == 0) {
	      free(original_function);
	      return 1;
	    }
	    *status = fit_levmar_internal(algorithm, ndata, nrfitparameters_one, func_params, func_params_err, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
					  &nr_params, 
#endif
					  epsabs, epsrel, maxiter, verbose);
	    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters_one, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
				 &nr_params, 
#endif
				 algorithm, 2, showcovariance, verbose) == 0) {
	      free(original_function);
	      return 1;
	    }
	    if(function->chi2 < old_chi2 || old_chi2 < 0) {
	      old_chi2 = function->chi2;
	      improved = 1;
	    }
	    loopnr++;
	    i++;  // This is the fit parameter (sequential)
	  }
	}
      }
    }while(improved);

    /* Run fit one more time with all fit parameters on simultaneously */
    for(j = 0; j < function->nrfuncs; j++) {
      for(k = 0; k < MaxNrFitParameters; k++)
	function->func[j].fit_flag[k] = original_function->func[j].fit_flag[k];
    }
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
			 &nr_params, 
#endif
			 algorithm, 0, showcovariance, verbose) == 0) {
      free(original_function);
      return 1;
    }
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
			 &nr_params, 
#endif
			 algorithm, 1, showcovariance, verbose) == 0) {
      free(original_function);
      return 1;
    }
    *status = fit_levmar_internal(algorithm, ndata, nrfitparameters, func_params, func_params_err, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
				  &nr_params, 
#endif
				  epsabs, epsrel, maxiter, verbose);
    free(original_function);
  }

  if(verbose.verbose) {
    printf("  After %ld iterations found chisq/dof = %g\n", levmar_internal_fitter_data.iter, levmar_internal_fitter_data.fitfunction->chi2_red);
  }

  /* Free up memory etc */
  if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
		       &nr_params, 
#endif
		       algorithm, 2, showcovariance, verbose) == 0) {
    return 1;
  }

  /* Rescale errorbars to ensure reduced chi2 will be 1 */
  if(force_chi2_1) {
    sig_scale = sqrt(function->chi2_red);
    for(i = 0; i < ndata; i++) {
      data_sigma[i] *= sig_scale;
    }
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
			 &nr_params, 
#endif
			 algorithm, 0, showcovariance, verbose) == 0) {
      return 1;
    }
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
			 &nr_params, 
#endif
			 algorithm, 1, showcovariance, verbose) == 0) {
      return 1;
    }
    *status = fit_levmar_internal(algorithm, ndata, nrfitparameters, func_params, func_params_err, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
				  &nr_params, 
#endif
				  epsabs, epsrel, maxiter, verbose);
    if(levmar_initialize(function, &func_params, &func_params_err, &nrfitparameters, data_x, data_y, data_sigma, ndata, &gsl_params, &psrsalsa_params, 
#ifdef NRAVAIL
			 &nr_params, 
#endif
			 algorithm, 2, showcovariance, verbose) == 0) {
      return 1;
    }
  }

  if(verbose.verbose || showresults) {
    printf("  Fit function: ");
    print_fitfunctions(function, 0, 1, -1, verbose);
    printf("  ");
    print_fitfunctions(function, 0, 0, -1, verbose);
    printf("\n  chi2 = %f\n  reduced chi2 = %f (%ld points - %d params)\n", function->chi2, function->chi2_red, ndata, nrfitparameters);
  }
  if(free_data_sigma) {
    free(data_sigma);
    data_sigma = NULL;
  }
  return *status;
}

//START REGION DEVELOP

// See, for instance, the general least squares discussion in NR
// Input:
// - The solution = the optimum values for the parameters for which the chi2 is minimised (must be done elsewhere). These are nrparams parameters passed to the function in the array optimum_params.
// - dydparam() is a function that takes for argument 1 the parameter number and for argument 2 the optimum_params array. The return value is the derivative of the fit function w.r.t. the parameter number (argument 1) at x=argument number 3
// - The chi2 of the solution, so NOT the reduced chi2
// - If rescale_errors is non-zero, the error-bars will be rescaled such that the REDUCED chi2=1, hence to incorporate a estimation for systematic errors. The reduced chi2 = chi2 / (nrdatapoints - nrfitparameters).
// - Information about the points that are fitted:
//     - Nr of data points = nrdatapoints
//     - The x values
//     - The y values are NOT required
//     - The y-errorbars = sigma
//     - If sigma = NULL, the errorbars are assumed to be equal and will be choosen such that the reduced chi2=1. The chi2 should be calculated by effectively taking all sigma values to be 1. This implies that rescale_errors is effectively always enabled in this case.
// - The covariance matrix is returned as the array covariance_matrix (nrparams*nrparams values). Memory should already be allocated
//
// NOTE: the error's on the parameters (to be procise the variance, so take the the sqrt if you want sigma, the std. dev.) can be found on the diagonal of the covariance matrix.
//
// Return values
//    0 = ok
//    1 = Error (such as memory allocation error)
//    2 = Not enough data points (nrdatapoints <= nrparams)
//    3 = Covariance matrix cannot be computed (matrix cannot be inverted, probably because there are fully covariant parameters?)
int compute_covariance_matrix(long nrparams, double *optimum_params, double (*dydparam)(long, double [], double), double chi2, int rescale_errors, long nrdatapoints, double *x, double *sigma, double *covariance_matrix, verbose_definition verbose_state)
{
  long datapoint, paramnr, paramnr2;
  double variance_scaling, variance;
  /*
  int nrtoasincluded, j, k;
  int *fitparam_list;
  long double chi2w, rms;
  double sigma;
  parfile_def parfile_tmp;
  double *mrq_alpha, double *mrq_alpha_inv, double *mrq_deriv;
  */

  if(nrdatapoints <= nrparams) {
    printerror(verbose_state.debug, "ERROR compute_covariance_matrix: Nr of data points should exceed the number of parameters.");
    return 2;
  }
  
  double *alpha, *derivatives;
  alpha = malloc(nrparams*nrparams*sizeof(double));
  derivatives = malloc(nrparams*sizeof(double));
  if(alpha == NULL || derivatives == NULL) {
    printerror(verbose_state.debug, "ERROR compute_covariance_matrix: Memory allocation error.");
    return 1;
  }
  for(paramnr = 0; paramnr < nrparams*nrparams; paramnr++) {
    alpha[paramnr] = 0;
  }

  if(sigma == NULL || rescale_errors) {  // No errorbars
    // chi2 = Sum[ (y_i - fit_i)**2 / sigma**2
    // If want reduced chi2 to be 1, then chi2 = dof = nrdatapoints - nrparams
    // nrdatapoints - nrparams = Sum[ (y_i - fit_i)**2 / sigma**2 ]
    // So Sum is too large by factor chi2/(nrdatapoints - nrparams)
    // So elements in sum needs to be divided by this factor
    // So sigma**2 -> sigma**2 * chi2/(nrdatapoints - nrparams) should do the trick
    variance_scaling = chi2/(double)(nrdatapoints - nrparams);
  }else {
    variance_scaling = 1;
  }

  for(datapoint = 0; datapoint < nrdatapoints; datapoint++) {
    for(paramnr = 0; paramnr < nrparams; paramnr++) {
      derivatives[paramnr] = (*dydparam)(paramnr, optimum_params, x[datapoint]);
    }
    if(sigma == NULL) {
      variance = variance_scaling;
    }else {
      variance = sigma[datapoint]*sigma[datapoint]*variance_scaling;
    }
    for(paramnr = 0; paramnr < nrparams; paramnr++) {
      for(paramnr2 = 0; paramnr2 < nrparams; paramnr2++) {
	alpha[nrparams*paramnr+paramnr2] += derivatives[paramnr]*derivatives[paramnr2]/variance;
      }
    }
  }

  if(linalg_solve_matrix_eq(alpha, paramnr, 0, NULL, 0, NULL, 0, covariance_matrix, 1, 1, verbose_state) != 0) {  // Failed
    printwarning(0, "WARNING tempo3 compute_covariance_matrix: Covariance matrix could not be determined. Two or more parameters might be 100%% covariant.");
    return 3;
  }

  free(derivatives);
  free(alpha);

  return 0;
}

