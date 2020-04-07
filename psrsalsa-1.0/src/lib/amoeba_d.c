//START REGION DEVELOP
/* This is an exact copy of amoeba.c (float version) to make the 
   double version. This file should be regeneratable by applying the
   following substitutions:

   float     ->  double 
   _psrsalsa  -> _psrsalsa_d
   _function -> _function_d
   _amoeba   -> _amoeba_d
   _nmsimplex -> _nmsimplex_d
   doAmoeba  -> doAmoeba_d
   %f        -> %lf
   vector    -> vector_d
   matrix    -> matrix_d
   amoeba_nr -> amoeba_nr_d
*/

//START REGION RELEASE
#include <stdlib.h>
#include <math.h>
#include "psrsalsa.h"
//START REGION DEVELOP
#ifdef NRAVAIL
  #include "nr.h"
  #include "nrutil.h"
#endif

//START REGION RELEASE
static int nrparams_internal_psrsalsa_d;   /* The number of parameters, including any fixed (not fitted for) parameters */
static int *fixed_internal_psrsalsa_d;     /* Global variable indicating which parameter numbers are fixed */
static double *xstart_internal_psrsalsa_d;  /* The start position of the searched parameters */
static double *x_internal_psrsalsa_d;       /* Temporary array containing the current probed parameter location */
static double (*funk_remember_user_function_d)(double []);  /* This points to the function that needs to be minimised */
static int algorithm_internal_psrsalsa_d;  /* Keeps track of which algorithm is used */

//START REGION DEVELOP
//START REGION RELEASE

/* Internal function that allows user to supply a function to be
   minimised with fixed parameters, functionality which is not
   supported by the actual algorithms used to minimise those
   functions. */
double funk_internal_psrsalsa_d(double x[])
{
  int i, j;
  if(algorithm_internal_psrsalsa_d == 1)
    j = 1;                /* NR starts counting from 1. */
  else
    j = 0;
  for(i = 0; i < nrparams_internal_psrsalsa_d; i++) {
    if(fixed_internal_psrsalsa_d[i] == 0) {
      x_internal_psrsalsa_d[i+1] = x[j++];
    }else {
      x_internal_psrsalsa_d[i+1] = xstart_internal_psrsalsa_d[i];
    }
  }
  return funk_remember_user_function_d(x_internal_psrsalsa_d+1);
}

//START REGION DEVELOP
//START REGION RELEASE

/* Does an amoeba search to minimize function funk, which has nrparams
   and is called with the array x[0]..x[nrparams-1]. The start point
   is given by xstart and the inital search step by dx. The nonzero
   values of fixed indicate if the parameter is fixed and not fitted
   for. xfit is filled with the position at which funk has a minimum
   with value yfit. All arrays start at entry 0, not 1. ftol (a small
   number) indicates the precision of the fit and nfunc is set to the
   number of iteration steps it took to converge. If finderrors is set
   it will search the parameterspace for all nonfixed parameters to
   find the sigma errors (assuming funk is a chi^2 surface, dplus &
   dmin are the two errorbars). GSL is not implemented. I cannot
   remember exactly why, but probably because it only supports double
   format.

   algorithm speciefies which code to use
     0: Modified code originally written by Michael F. Hutt
     1: Modified code originally from Numerical Recipes
     2: GSL (unsupported)

   Returns 0 on success.
   Returns 1 on maximum number of itterations exceeded error (try smaller ftol)
   Returns 2 on memory error
   Returns 3 if less than 2 fit parameter are not fixed
   Returns 4 if algorithm is not available
*/
int doAmoeba_d(int algorithm, double *xstart, double *dx, int *fixed, double *xfit, double *yfit, int nrparams, double (*funk)(double []), double ftol, int *nfunk, int verbose, int finderrors, double sigma, double *dplus, double *dmin)
{
  int i, j, nfitparameters, ret;
//START REGION DEVELOP
//START REGION RELEASE
  extern double amoeba_nmsimplex_d(double (*objfunc)(double[]), double start[], double dx[], int n, double EPSILON, int *nritterations, double *reachedEpsilon, int verbose);
//START REGION DEVELOP
#ifdef NRAVAIL
  int n;
  double **p_nr, *y_nr;
  extern int amoeba_nr_d(double **p, double y[], int ndim, double ftol, double (*funk)(double []), int *nfunk);
#endif
//START REGION RELEASE

  if(algorithm < 0 || algorithm > 1) {
    fprintf(stderr, "ERROR doAmoeba_d: Unknown algorithm requested.\n");
    return 4;
  }


//START REGION DEVELOP
#ifndef NRAVAIL
  if(algorithm == 1) {
    fprintf(stderr, "ERROR doAmoeba_d: Numerical Recipies requested, but not included during compilation.\n");
    return 4;
  }
#endif
//START REGION RELEASE

  /* Determine the number total number of parameters for which we are
     fitting. */
  nfitparameters = 0;
  for(i = 0; i < nrparams; i++) {
    if(fixed[i] == 0)
      nfitparameters++;
  }
  if(nfitparameters < 2) {
    fprintf(stderr, "ERROR doAmoeba_d: Cannot fit for less than 2 parameters (now have %d).\n", nfitparameters);
    return 3;
  }

  /* Store some global parameters */
  fixed_internal_psrsalsa_d = fixed;
  xstart_internal_psrsalsa_d = xstart;
  nrparams_internal_psrsalsa_d = nrparams;
  funk_remember_user_function_d = funk;
  algorithm_internal_psrsalsa_d = algorithm;

  if(algorithm == 0) {
    //    int nritterations;
    double reachedEpsilon, *xstart_nmsimplex_d, *dx_nmsimplex_d;
    xstart_nmsimplex_d = malloc(nfitparameters*sizeof(double));
    dx_nmsimplex_d = malloc(nfitparameters*sizeof(double));
    x_internal_psrsalsa_d = malloc((nrparams+1)*sizeof(double));
    if(xstart_nmsimplex_d == NULL || dx_nmsimplex_d == NULL || x_internal_psrsalsa_d == NULL) {
      fprintf(stderr, "ERROR doAmoeba_d: Memory allocation error.\n");
      return 2;
    }

    /* Construct start and dx vectors.*/
    j = 0;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
	xstart_nmsimplex_d[j] = xstart[i];
	dx_nmsimplex_d[j] = dx[i];
	j++;
      }
    }
    /* Do the amoeba search. */
    *yfit = amoeba_nmsimplex_d(funk_internal_psrsalsa_d, xstart_nmsimplex_d, dx_nmsimplex_d, nfitparameters, ftol, nfunk, &reachedEpsilon, 0);
    /* Copy the end point at minimum.*/
    j = 0;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
	xfit[i] = xstart_nmsimplex_d[j];
	j++;
      }else {
	xfit[i] = xstart[i];
      }
    }
   
    /* Release memory */
    free(xstart_nmsimplex_d); 
    free(dx_nmsimplex_d);
    free(x_internal_psrsalsa_d);

    /* Check if converged */
    if(reachedEpsilon > ftol)
      return 1;
  }

//START REGION DEVELOP
#ifdef NRAVAIL
  if(algorithm == 1) {
    /* Allocate simplex and a vector with y values of function. */
    p_nr = matrix_d(1, nfitparameters+1, 1, nfitparameters);
    y_nr = vector_d(1, nfitparameters+1);
    x_internal_psrsalsa_d = vector_d(1, nrparams+1);
    if(p_nr == NULL || y_nr == NULL || x_internal_psrsalsa_d == NULL) {
      fprintf(stderr, "ERROR doAmoeba_d: Memory allocation error.\n");
      return 2;
    }
    
    /* Fill matrix_d with starting points */
    j = 1;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
	for(n = 1; n <= nfitparameters+1; n++)
	  p_nr[n][j] = xstart[i];
	j++;
      }
    }

    /* matrix_d except first column has xstart+dx on diagonal */
    j = 1;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
	p_nr[j+1][j] += dx[i];
	j++;
      }
    }
    
    /* Initialize the array y, which contains the y values corresponding
       with funk(x). */
    for(i = 0; i < nfitparameters+1; i++) {
      y_nr[i+1] = funk_internal_psrsalsa_d(p_nr[i+1]);
    }
    
    if(verbose) {
      for(j = 1; j <= nfitparameters; j++) {
	if(j == 1)
	  fprintf(stderr, "matrix_d p: ");
	else
	  fprintf(stderr, "          ");
	for(n = 1; n <= nfitparameters+1; n++) {
	  fprintf(stderr, "%lf ", p_nr[n][j]);
	}
	fprintf(stderr, "\n");
      }
    }
    
    /* Do the amoeba search. */
    ret = amoeba_nr_d(p_nr, y_nr, nfitparameters, ftol, funk_internal_psrsalsa_d, nfunk);
    if(ret != 0)
      return ret;
    /* Remember the y value of funk at minimum. */
    *yfit = funk_internal_psrsalsa_d(p_nr[1]);
    
    /* Copy the end point at minimum.*/
    j = 1;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
	xfit[i] = p_nr[1][j];
	j++;
      }else {
	xfit[i] = xstart[i];
      }
    }
    
    /* Release memory */
    free_matrix_d(p_nr,1, nfitparameters+1, 1, nfitparameters);
    free_vector_d(y_nr, 1, nfitparameters+1);
    free_vector_d(x_internal_psrsalsa_d, 1, nfitparameters+1);
  }
#endif   /* End of NR initialization */
//START REGION RELEASE


  /* Reset errorbars if requested */
  for(i = 0; i < nrparams; i++) {
    if(dplus != NULL)
      dplus[i] = 0;
    if(dmin != NULL)
      dmin[i] = 0;
  }


  if(finderrors && dplus != NULL && dmin != NULL) {
    if(nfitparameters < 3) {
      fprintf(stderr, "ERROR doAmoeba_d: Cannot estimate errors if the number of fit parameters is less than 3.\n");
    }else {
      for(i = 0; i < nrparams; i++) {
	ret = find_errors_amoeba_d(algorithm, dx, fixed, xfit, *yfit, nrparams, funk, ftol, i, &dplus[i], &dmin[i], sigma);
	if(ret != 0) {
	  fprintf(stderr, "ERROR doAmoeba_d: find_errors_amoeba_d failed with error code %d\n", ret);
	  return ret;
	}
	/*
	if(isnan(dplus[i]) || isnan(dmin[i])) {
	  fprintf(stderr, "doAmoeba_d: find_errors_amoeba_d didn't converge for parameter %d with values %lf %lf\n", i, dplus[i], dmin[i]);
	}
	*/
      }
    }
  }

  return 0;
}

//START REGION DEVELOP
//START REGION RELEASE


/* Copy parameters of doAmoeba_d, but also give paramnr for which to
   find the errorbar at sigma times the minimum. 
   Returns 0 on success.
   Returns 1 on NMAX exceed error (smaller ftol doesn't work)
   Returns 2 on memory error
   Returns 3 if less than 2 fit parameter are not fixed
*/
int find_errors_amoeba_d(int algorithm, double *dx, int *fixed, double *xfit, double yfit, int nrparams, double (*funk)(double []), double ftol, int paramnr, double *dplus, double *dmin, double sigma)
{
  double *xstartnew, *xfitnew, *dxnew, yfitnew, x0, x1, yold, step, fsign;
  int *fixednew, j, nfunknew, ret, sign;
  long n, n2;
  fixednew = malloc(nrparams*sizeof(int));
  xstartnew = malloc(nrparams*sizeof(double));
  dxnew = malloc(nrparams*sizeof(double));
  xfitnew = malloc(nrparams*sizeof(double));
  if(fixednew == NULL || xstartnew == NULL || xfitnew == NULL) {
    fprintf(stderr, "ERROR find_point_amoeba_d: Memory allocation error.\n");
    return 2;
  }
  if(fixed[paramnr] == 0) {
    /* Set fixed for the parameter we want to know the errorbar of. */
    for(j = 0; j < nrparams; j++) {
      fixednew[j] = fixed[j];
      if(j == paramnr)
	fixednew[j] = 1;
    }
    for(sign = 0; sign < 2; sign++) {
      if(sign == 0)
	fsign = 1;
      else
	fsign = -1;
      /* Start at minimum */
      for(j = 0; j < nrparams; j++) {
	xstartnew[j] = xfit[j];
	dxnew[j] = dx[j];
      }
      x0 = xstartnew[paramnr];
      step = fsign*dx[paramnr];
      yold = yfit;
      n = 0;
      n2 = 0;
      do {
	n++;
	x1 = x0 + step;
	xstartnew[paramnr] = x1;
	dxnew[paramnr] = step;
	do {
	  ret = doAmoeba_d(algorithm, xstartnew, dxnew, fixednew, xfitnew, &yfitnew, nrparams, funk, ftol, &nfunknew, 0, 0, sigma, NULL, NULL);
	  if(ret == 3) {
	    free(fixednew);
	    free(xstartnew);
	    free(dxnew);
	    free(xfitnew);
	    return 3;
	  }
	  if(ret == 1) {
	    fprintf(stderr, "WARNING find_point_amoeba_d: Adjusting amoeba tollerance to try to converge (for errorbar estimation).\n");
	    ftol *= 10;
	  }
	  if(ret == 0) {
	    break;
	  }
	}while(ftol < 0.01);
	if(ret == 1) {
	  free(fixednew);
	  free(xstartnew);
	  free(dxnew);
	  free(xfitnew);
	  return 1;
	}
	/*	fprintf(stderr, "pre  paramnr=%d sign=%.0f n=%ld x0=%lf x1=%lf dx=%lf y=%lf (%lf) steps=%d\n", paramnr, fsign, n, x0, x1, step, yfitnew, yfit*(sigma+1.0), nfunknew); */
	if(yfitnew < yfit*(sigma+1.0) && yold < yfit*(sigma+1.0)) {
	  /* Step step was not enough to cross point, do same step again. */
	  x0 = x1;
	}else if(yfitnew > yfit*(sigma+1.0) && yold < yfit*(sigma+1.0)) {
	  /* Crossed the point, so take half a step. */
	  step = 0.5*step;
	}else if(yfitnew < yfit*(sigma+1.0) && yold > yfit*(sigma+1.0)) {
	  /* Crossed the point, so take half a step. */
	  step = 0.5*step;
	}else {
	  /* only other option: Both >, so do same step again.*/
	  x0 = x1;
	}
	/* Do a restart so hopefully we get a converging solution */
	if(n % 100 == 0) {
	  /* Quit loop if within 5% of answer, hopefully that is (within the error on the error) good enough.  */
	  if(fabs((yfitnew-yfit*(sigma+1.0))/yfit*(sigma+1.0)) < 0.05) {
	    fprintf(stderr, "ERROR find_errors_amoeba_d: error parameter %d didn't quite converge, but possibly it is close enough: chi %lf != %lf\n", paramnr, yfitnew, yfit*(sigma+1.0));
	    yfitnew = yfit*(sigma+1.0);
	  }else {
	    //	  x0 = x0 + (x0-x1)*1.234;
	    /* Make the new start point inbetween the best fit solution and the current place where the fitting failed */
	    x0 = xfit[paramnr] + (x0-xfit[paramnr])*0.987654321;
	    step = fsign*dx[paramnr];
	    n2 ++;
	  }
	}
	/*	fprintf(stderr, "post paramnr=%d sign=%.0f n=%ld x0=%lf x1=%lf dx=%lf y=%lf (%lf) steps=%d\n", paramnr, fsign, n, x0, x1, step, yfitnew, yfit*(sigma+1.0), nfunknew); */
      }while(fabs((yfitnew-yfit*(sigma+1.0))/yfit*(sigma+1.0)) > 0.001 && n2 < 10);
      if(n2 == 10) {
	fprintf(stderr, "ERROR find_errors_amoeba_d: converging error parameter %d: chi %lf != %lf\n", paramnr, yfitnew, yfit*(sigma+1.0));
	x0 = sqrt(-1);
	/*      }else {
		fprintf(stderr, "find_errors_amoeba_d: found x0=%lf (xstart=%lf)\n", x0, xfit[paramnr]);	*/
      }
      if(sign == 0)
	*dplus = x0 - xfit[paramnr];
      else
	*dmin = x0 - xfit[paramnr];
    }
  }else {   /* Fixed, so no errorbar */
    *dplus = 0;
    *dmin = 0;
  }
  free(fixednew);
  free(xstartnew);
  free(dxnew);
  free(xfitnew);
  return 0;
}

//START REGION DEVELOP
