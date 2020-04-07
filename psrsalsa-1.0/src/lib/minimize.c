//START REGION RELEASE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include "gsl/gsl_roots.h"
#include "psrsalsa.h"
//START REGION DEVELOP


//START REGION RELEASE

// Returns 
//  0 on success
//  2 = Did not find root in specified range
//  4 = Memory allocation error
int minimize_1D_double_refine_borders(int findroot, double (*funk)(double, void *), void *params, int gridsearch, int investigateLocalMinima, double *x_lower, double *x_upper, int debug_verbose)
{
  double x, y, ymin[3], ymin_new[3], dy[3], ymax, x_lower_new, x_upper_new, sign_old, gridsearch_margin, imin_new[3], i_originalgrid, imin1, imin2; //, limit; //, *grid;
  int imin[3], unset[3], i, minima;  // , imax

  /*
  grid = malloc(gridsearch*sizeof(double));
  if(grid == NULL) {
    fflush(stdout);
    fprintf(stderr, "minimize_1D_double_refine_borders: Memory allocation error\n");
    return 4;
  }
  */

  unset[0] = unset[1] = unset[2] = 1;

  for(i = 0; i < gridsearch; i++) {
    x = *x_lower + i*(*x_upper - *x_lower)/(double)(gridsearch-1);
    y = funk(x, params);
    //    grid[i] = y;
    if(debug_verbose) {
      fflush(stdout); fprintf(stderr, "x=%f y=%f\n", x, y);
    }
    if(findroot) {
      if(i == 0) {
	sign_old = y;
	imin[0] = 0;
      }else {
	if(y < 0 && sign_old > 0) {
	  imin[0] = i;
	}else if(y > 0 && sign_old < 0) {
	  imin[0] = i;
	}
	if(imin[0] != 0)
	  break;
      }
    }else {
      if(unset[0]) {
	imin[0] = i;
	ymin[0] = y;
	unset[0] = 0;
      }else if(unset[1]) {
	imin[1] = i;
	ymin[1] = y;
	unset[1] = 0;
      }else if(unset[2]) {
	imin[2] = i;
	ymin[2] = y;
	unset[2] = 0;
      }else {
	for(minima = 0; minima < 3; minima++) {
	  dy[minima] = ymin[minima] - y;
	}
	if(dy[0] > dy[1] && dy[0] > dy[2]) {
	  if(dy[0] > 0) {
	    imin[0] = i;
	    ymin[0] = y;	
	  }  
	}else if(dy[1] > dy[0] && dy[1] > dy[2]) {
	  if(dy[1] > 0) {
	    imin[1] = i;
	    ymin[1] = y;	  
	  }
	}else {
	  if(dy[2] > 0) {
	    imin[2] = i;
	    ymin[2] = y;	  
	  }
	}
      }
      if(i == 0 || y > ymax) {
	//	imax = i;
	ymax = y;
      }
    }
    //    if(debug_verbose) {
    //      fflush(stdout); fprintf(stderr, "imin=%d ymin=%f imin=%d ymin=%f imin=%d ymin=%f\n", imin[0], ymin[0], imin[1], ymin[1], imin[2], ymin[2]);
    //    }
  }

  // Now we know the minimum in the grid, make sure to include any
  // (possibly unresolved) minima within a certain margin, to avoid homing
  // in on a wrong local minimum.
  /*
  if(findroot == 0) {
    limit = ymin + (ymax-ymin)*0.1;
    imin1 = 0;
    for(i = 1; i < gridsearch; i++) {
      if(grid[i] < limit) {
	imin1 = i;
	break;
      }
    }
    imin2 = gridsearch-1;
    for(i = gridsearch-1; i >= 0; i--) {
      if(grid[i] < limit) {
	imin2 = i;
	break;
      }
    }
    if(imin1 != imin2 && debug_verbose) {
      fflush(stdout); fprintf(stderr, "Found more than one possible minima\n");
    }
  }
  */

  // Now we know the minimum in the grid, make sure to include other (possibly unresolved) minima.
  if(findroot == 0) {
    imin1 = imin[0];
    imin2 = imin[0];
    if(imin[1] < imin1)
      imin1 = imin[1];
    if(imin[2] < imin1)
      imin1 = imin[2];
    if(imin[1] > imin2)
      imin2 = imin[1];
    if(imin[2] > imin2)
      imin2 = imin[2];
    
    if(imin2-imin1 != 2 && debug_verbose) {
      fflush(stdout); fprintf(stderr, "Found more than one possible minima (range = %.2f .. %.2f)\n", imin1, imin2);
    }
    if(imin2-imin1 > 0.3*gridsearch) {
      if(investigateLocalMinima == 0) {
	imin1 = imin[0];
	imin2 = imin[0];
	if(ymin[1] < ymin[0]) {
	  imin1 = imin[1];
	  imin2 = imin[1];
	}else if(ymin[2] < ymin[0] && ymin[2] < ymin[1]) {
	  imin1 = imin[2];
	  imin2 = imin[2];
	}
      }else {
	if(debug_verbose) {
	  fflush(stdout); fprintf(stderr, "Refining gridsearch around local minima\n");
	}
	
	unset[0] = unset[1] = unset[2] = 1;
	for(minima = 0; minima < 3; minima++) {
	  for(x = *x_lower + (imin[minima]-1.5)*(*x_upper - *x_lower)/(double)(gridsearch-1); x < *x_lower + (imin[minima]+1.5)*(*x_upper - *x_lower)/(double)(gridsearch-1); x += 0.09*(*x_upper - *x_lower)/(double)(gridsearch-1)) {
	    y = funk(x, params);
	    if(debug_verbose) {
	      fflush(stdout); fprintf(stderr, "x=%f y=%f\n", x, y);
	    }
	    //	  x = *x_lower + i*(*x_upper - *x_lower)/(double)(gridsearch-1)
	    //	    (x - *x_lower)*(gridsearch-1) = i*(*x_upper - *x_lower)
	    i_originalgrid = (x - *x_lower)*(double)(gridsearch-1)/(*x_upper - *x_lower);
	    if(unset[0] || y < ymin_new[0]) {
	      imin_new[0] = i_originalgrid;
	      ymin_new[0] = y;
	      unset[0] = 0;
	    }else if(unset[1] || y < ymin_new[1]) {
	      imin_new[1] = i_originalgrid;
	      ymin_new[1] = y;
	      unset[1] = 0;
	    }else if(unset[2] || y < ymin_new[2]) {
	      imin_new[2] = i_originalgrid;
	      ymin_new[2] = y;
	      unset[2] = 0;
	    }
	  }
	}
      
	imin1 = imin_new[0];
	imin2 = imin_new[0];
	if(imin_new[1] < imin1)
	  imin1 = imin_new[1];
	if(imin_new[2] < imin1)
	  imin1 = imin_new[2];
	if(imin_new[1] > imin2)
	  imin2 = imin_new[1];
	if(imin_new[2] > imin2)
	  imin2 = imin_new[2];
	
	if(imin2-imin1 != 2 && debug_verbose) {
	  fflush(stdout); fprintf(stderr, "Found more than one possible minima (range = %.2f .. %.2f)\n", imin1, imin2);
	}
	if(imin2-imin1 > 0.3*gridsearch) {
	  fflush(stdout); fprintf(stderr, "minimize_1D_double_refine_borders: Refining gridsearch around local minima failed\n");
	  exit(0);
	}
      }
    }
  }

  gridsearch_margin = 0.1*(double)gridsearch;
  if(gridsearch_margin < 1)
    gridsearch_margin = 1;
  if(findroot) {
    if(imin[0] == 0) {
      return 2;
    }
    x_lower_new = *x_lower + (imin[0]-gridsearch_margin)*(*x_upper - *x_lower)/(double)(gridsearch-1);
    x_upper_new = *x_lower + (imin[0]+gridsearch_margin-1)*(*x_upper - *x_lower)/(double)(gridsearch-1);
    *x_lower = x_lower_new;
    *x_upper = x_upper_new;
  }else {
    x_lower_new = *x_lower + (imin1-gridsearch_margin)*(*x_upper - *x_lower)/(double)(gridsearch-1);
    x_upper_new = *x_lower + (imin2+gridsearch_margin)*(*x_upper - *x_lower)/(double)(gridsearch-1);
    *x_lower = x_lower_new;
    *x_upper = x_upper_new;
  }
  //  free(grid);
  return 0;
}

//START REGION DEVELOP
//START REGION RELEASE

/*
  Given the function: double funk(double x, void *params)

  Minimize the function as function of x (if findroot == 0), otherwise
  find the root of the function. The variable params can be used to pass
  on other (fixed) parameters. The minimization is done using gsl with
  the Brent method.
  
  As input, the search range x_lower and x_upper must be known. The
  result is returned as x_minumum.  epsabs sets the minimum required
  uncertainty in x_minimum, while epsrel set the minimum required
  uncertainty as a fraction of x_minimum. It is allowed to set either
  epsabs or epsrel to zero. max_iter is the maximum nr of iterrations
  (maybe 100).
  
  If gridsearch is > 1, the specified x range in searched over in
  gridsearch points to find the rough location of the minimum first,
  before bracketing the minimum further. This initial gridsearch is
  repeated nested + 1 times. If not set correctly, the minimum can be
  missed. If investigateLocalMinima is set, 3 local minima are
  explored in a higher resolution in order to potentially discriminate
  between them and zoom in on what is hopefully the global
  minimum. This is probably not a good idea if the fit parameter is
  circular (like an angle), as the best solution can be at the lower
  and upper border of the parameter range simultaneously, which
  results in the routine not being able to distinguish between the two
  minima, resulting in an error.
  
  If set, verbose determines the nr of spaces before the output. If
  debug_verbose is set, function values are outputted when doing the
  nested gridsearch.
  
  Return values:
  0 = success, converged
  1 = Maximum nr of itterations reached, did not converge fully
  2 = Did not find root in specified range
  3 = Lower and upper limit do not bracket a root
  4 = memory allocation error
  5 = Other unspecified error
*/
int minimize_1D_double(int findroot, double (*funk)(double, void *), void *params, double x_lower, double x_upper, int gridsearch, int investigateLocalMinima, int nested, double *x_minimum, int max_iter, double epsabs, double epsrel, int verbose, int debug_verbose)
{
  const gsl_min_fminimizer_type *T_minimizer;
  gsl_min_fminimizer *s_minimizer;
  const gsl_root_fsolver_type *T_root;
  gsl_root_fsolver *s_root;
  int status;
  int iter, i, ret, nest;
  gsl_function F;
  F.function = funk;
  F.params = params;
  iter = 0;

#if GSL_VERSION_NUMBER < 102
  printerror(0, "ERROR minimize_1D_double: Not supported for GSL < 1.2");  
  exit(0);
#endif

  if(verbose) {
    for(i = 0; i < verbose - 1; i++)
      printf(" ");
    if(findroot == 0)
      printf("Minimising function between %f and %f\n", x_lower, x_upper);
    else
      printf("Finding root of function between %f and %f\n", x_lower, x_upper);
  }

  if(gridsearch > 1) {
    if(verbose) {
      for(i = 0; i < verbose - 1; i++)
	printf(" ");
      printf("  Refining boundaries with a %d point grid search\n", gridsearch);
    }

    for(nest = 0; nest < 1+nested; nest++) {
      ret = minimize_1D_double_refine_borders(findroot, funk, params, gridsearch, investigateLocalMinima, &x_lower, &x_upper, debug_verbose);
      if(ret != 0) {
	if(ret == 1) {
	  if(verbose) {
	    for(i = 0; i < verbose - 1; i++)
	      printf(" ");
	    printf("  Maximum itterations (%d) exceeded while refining borders\n", max_iter);
	  }
	}
	return ret;
      }
      if(verbose) {
	for(i = 0; i < verbose - 1; i++)
	  printf(" ");
	printf("  Refined boundaries are %f and %f\n", x_lower, x_upper);
      }
    }
  }

  *x_minimum = 0.5*(x_upper+x_lower);

  //  fprintf(stderr, "Start init\n");
  if(findroot == 0) {
    T_minimizer = gsl_min_fminimizer_brent;
    s_minimizer = gsl_min_fminimizer_alloc(T_minimizer);
    fflush(stdout);
    //    fprintf(stderr, "Start gsl_min_fminimizer_set: %f %f %f\n", x_lower, *x_minimum, x_upper);
    //    double val1, val2, val3;
    //    val1 = funk(x_lower, params);
    //    val2 = funk(*x_minimum, params);
    //    val3 = funk(x_upper, params);
    //    fprintf(stderr, "y = %f %f %f\n", val1, val2, val3);
    gsl_min_fminimizer_set(s_minimizer, &F, *x_minimum, x_lower, x_upper);
    //    fprintf(stderr, "Stop gsl_min_fminimizer_set\n");
  }else {
    T_root = gsl_root_fsolver_brent;
    s_root = gsl_root_fsolver_alloc(T_root);
    if(gridsearch == 0) {
      double val1, val2;
      val1 = funk(x_lower, params);
      val2 = funk(x_upper, params);
      if((val1 > 0 && val2 > 0) || (val1 < 0 && val2 < 0)) {
	return 3;
      }
    }
    //    fprintf(stderr, "Start gsl_root_fsolver_set: %f %f\n", x_lower, x_upper);
    gsl_root_fsolver_set(s_root, &F, x_lower, x_upper);
    //    fprintf(stderr, "Stop gsl_root_fsolver_set\n");
  }
  //  fprintf(stderr, "Stop init\n");
  /*
    if(findroot) {
    printf ("using %s method\n",
    gsl_min_fminimizer_name (s_minimizer));
}else {
    printf ("using %s method\n",
    gsl_root_fsolver_name(s_root));
}
  printf ("%5s [%9s, %9s] %9s %10s %9s\n",
	  "iter", "lower", "upper", "min",
	  "err", "err(est)");
  printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
	  iter, x_lower, x_upper,
	  *x_minimum, x_upper - x_lower);
	  */
  
  //#include "gsl/gsl_errno.h"
			//			gsl_error_handler_t *oldhandler;
			//			oldhandler = gsl_set_error_handler_off();
			//			gsl_set_error_handler(oldhandler); 
  do {
    iter++;
    //    fprintf(stderr, "Start ittr\n");
    if(findroot == 0) {
      status = gsl_min_fminimizer_iterate(s_minimizer);
#if GSL_VERSION_NUMBER >= 102
      *x_minimum = gsl_min_fminimizer_x_minimum (s_minimizer);
#endif
      x_lower = gsl_min_fminimizer_x_lower(s_minimizer);
      x_upper = gsl_min_fminimizer_x_upper(s_minimizer);
    }else {
      status = gsl_root_fsolver_iterate(s_root);
      *x_minimum = gsl_root_fsolver_root(s_root);
      x_lower = gsl_root_fsolver_x_lower(s_root);
      x_upper = gsl_root_fsolver_x_upper(s_root);
    }
    //    fprintf(stderr, "Stop ittr %d\n", status);
    //    fprintf(stderr, "Test start\n");
    status = gsl_min_test_interval (x_lower, x_upper, epsabs, epsrel);
    //    fprintf(stderr, "Test stop %d\n", status);
    /*
    if (status == GSL_SUCCESS)
      printf ("Converged:\n");
    printf ("%5d [%.7f, %.7f] "
	    "%.7f %.7f\n",
	    iter, x_lower, x_upper,
	    *x_minimum, x_upper - x_lower);
	    */
  }while(status == GSL_CONTINUE && iter < max_iter);
  //  fprintf(stderr, "Start free\n");
  if(findroot)
    gsl_root_fsolver_free(s_root);
  else
    gsl_min_fminimizer_free(s_minimizer);
  //  fprintf(stderr, "Stop free %d\n", status);

  if(status == GSL_SUCCESS) {
    if(verbose) {
      for(i = 0; i < verbose - 1; i++)
	printf(" ");
      if(findroot)
	printf("  Root found at %f with error %f\n", *x_minimum, fabs(x_upper - x_lower));
      else
	printf("  Minimum found at %f with error %f\n", *x_minimum, fabs(x_upper - x_lower));
    }
    return 0;
  }else if(iter == max_iter) {
    if(verbose) {
      for(i = 0; i < verbose - 1; i++)
	printf(" ");
      printf("  Maximum itterations (%d) exceeded. Value confined to [%f, %f]\n", max_iter, x_lower, x_upper);
    }
    return 1;
  }else {
    if(verbose) {
      for(i = 0; i < verbose - 1; i++)
	printf(" ");
      printf("  Unknown error during minimization\n");
    }
    return 5;
  }

  // XXXXX work out return values
  return status;
}

//START REGION DEVELOP
//START REGION RELEASE

double (*internal_funk_provided_by_user)(double *, void *);
int internal_paramnr;
double *internal_xminimum;
double internal_desired_chi2;

double internal_find_1D_error_funk(double x, void *params)
{
  double chi2;
  internal_xminimum[internal_paramnr] = x;
  chi2 = internal_funk_provided_by_user(internal_xminimum, params);
  //    printf("XXXX chi2 = %f\n", chi2);
  chi2 -= internal_desired_chi2;
  //  printf("XXXX chi2 offset = %f\n", chi2);
  return chi2;  
}

//START REGION DEVELOP
//START REGION RELEASE

/*
  You provide the function: double funk(double *x, void *params),
  which returns the chi2 as function of an array of parameters
  x[nrparameters] (set by xminimum) and potenital additional (fixed)
  parameter params. Of these parameters x[], parameter number paramnr
  (counting from zero) is varied to find an errorbar on that
  parameter, which is defined where the chi2 is (1+sigma) times higher
  than the minimum value chi2min. The chi2 at parameters xminimum
  should be within this limit. Set sigma to a positive value to find
  the errorbar in the positive direction, or to a negative value to
  find the errorbar in the negative direction. The size of the
  errorbar is returned as errorbar.

  The minimization is done using gsl with the Brent method.
  
  As input, xminimum is provided, which is the x value for which funk
  has a minimum (potentially found with function
  minimize_1D_double). Also dx must be provided, which is the stepsize
  used to find the specified sigma point. If separation exceeds dxmax
  (set to negative value to disable this feature), the errorbar is set
  to this value and no further searching is done.

  epsabs sets the minimum required uncertainty in the determined errorbar, while
  epsrel set the minimum required uncertainty as a fraction of
  the value. It is allowed to set either epsabs or epsrel to
  zero. max_iter is the maximum nr of iterrations (maybe 100).
  
  If set, verbose determines the nr of spaces before the output.
  
  Return values:
  0 = success, converged
  1 = Maximum nr of itterations reached, did not converge fully
  2 = Did not find root in specified range
  3 = Lower and upper limit do not bracket a root
  4 = Memory allocation error
  5 = Other unspecified error
*/
int find_1D_error(double (*funk)(double *, void *), double *xminimum, int paramnr, int nrparameters, double dx, double dxmax, void *params, double sigma, double chi2min, int max_itr, double epsabs, double epsrel, double *errorbar, int verbose)
{
  int ittr, i, ret, debug_verbose;
  double x_lower, x_upper, diff, xval, *xminimum_fiddle;

  debug_verbose = 0;

  dx = fabs(dx);
  if(verbose) {
    for(i = 0; i < verbose - 1; i++)
      printf(" ");
    printf("Finding error at value %f with stepsize %f\n", xminimum[paramnr], dx);
  }

  xminimum_fiddle = malloc(nrparameters*sizeof(double));
  if(xminimum_fiddle == NULL) {
    fflush(stdout);
    fprintf(stderr, "ERROR find_1D_error: Memory allocation error\n");
    return 4;
  }
  memcpy(xminimum_fiddle, xminimum, nrparameters*sizeof(double));

  internal_paramnr = paramnr;
  internal_xminimum = xminimum_fiddle;
  internal_funk_provided_by_user = funk;

  internal_desired_chi2 = chi2min*(1+fabs(sigma));
  if(verbose) {
    for(i = 0; i < verbose - 1; i++)
      printf(" ");
    printf("  Mimimum chi2 = %f, so %f sigma point corresponds to %f\n", chi2min, fabs(sigma), internal_desired_chi2);
  }

  x_lower = xminimum[paramnr];
  x_upper = xminimum[paramnr];
  ittr = 0;
  diff = internal_find_1D_error_funk(x_upper, params);
  if(diff > 0) {
    fflush(stdout);
    fprintf(stderr, "ERROR find_1D_error: Function called with initial parameters outside specified sigma limit: chi2 = %f higher than sigma border (%f).\n", diff, internal_desired_chi2);
    exit(0);
  }
  do {
    if(sigma >= 0) {
      x_upper += dx;
      if(dxmax >= 0) {
	if(fabs(x_upper - xminimum[paramnr]) > dxmax) {
	  *errorbar = dxmax;
	  free(xminimum_fiddle);
	  return 0;
	}
      }
      diff = internal_find_1D_error_funk(x_upper, params);
    }else {
      x_lower -= dx;
      if(dxmax >= 0) {
	if(fabs(x_lower - xminimum[paramnr]) > dxmax) {
	  *errorbar = dxmax;
	  free(xminimum_fiddle);
	  return 0;
	}
      }
      diff = internal_find_1D_error_funk(x_lower, params);
    }
    ittr++;
    //    printf("Trial %d: [%f %f] -> %f\n", ittr, x_lower, x_upper, diff);
    if(ittr == max_itr) {
      free(xminimum_fiddle);
      return 1;
    }
  }while(diff < 0);   // Continue until found a bracket for which chi2 is larger than requested value

  if(verbose) {
    for(i = 0; i < verbose - 1; i++)
      printf(" ");
    printf("  Found brackets: [%f %f]\n", x_lower, x_upper);
    verbose += 2;
  }

  ret = minimize_1D_double(1, internal_find_1D_error_funk, params, x_lower, x_upper, 0, 0, 0, &xval, max_itr, epsabs, epsrel, verbose, debug_verbose);

  *errorbar = fabs(xval - xminimum[paramnr]);
  if(verbose) {
    verbose -= 2;
    for(i = 0; i < verbose - 1; i++)
      printf(" ");
    if(ret == 0)
      printf("  Found errorbar: %f\n", *errorbar);
    else if(ret == 1)
      printf("  Maximum nr of itterations reached, did not converge fully\n");
    else if(ret == 2)
      printf("  Did not find root in specified range\n");
    else if(ret == 3)
      printf("  Lower and upper limit do not bracket a root\n");
    else {
      ret = 5;
      printf("  Other unspecified error\n");
    }
  }
  free(xminimum_fiddle);
  return ret;
}

//START REGION DEVELOP
