#include <math.h>           
#include <string.h>
#include "psrsalsa.h"

/* This function is similar to fit_levmar2_ld(), so make changes in both
   if required. 

   This function is designed to look like doAmoeba(), allowing the
   user to switch easily between the two optimisation strategies. Only
   the psrsalsa fitter is being supported at the moment (not the gsl
   or NR equivalent). Switching between algorithms would be
   complicated because the user specified function needs to be quite
   different to be used by the different algorithms.

   This function fits a user-defined function to a data-set. This
   function has nrparams parameters and the initial values are defined
   by the array xstart[0]..xstart[nrparams-1]. The nonzero values of
   fixed indicate if the parameter is fixed and not fitted for. xfit
   is filled with the found solution at which the chi2 has a minimum
   with value called "yfit" to use the same language as the doAmoeba()
   function. ftol (a small number) indicates the precision of the fit
   (in this function it corresponds to the maximum required fractional
   change in the fit parameters before calling the solution converged)
   and nfunc is set to the number of iteration steps it took to
   converge. If finderrors is set it will use the determined
   covariance matrix to find the sigma-level errors (as set by sigma)
   for all nonfixed parameters.  This is stored in dplus and/or dmin
   if != NULL. The values in dplus and dmin will be identical per
   definition, but both can be set to make the function more
   comparable with doAmoeba().

   A function func needs to be provided which is of the form:

   void func(func_params, alpha, beta, &chisq, info);

   So for a set of func_params (nrparams long), is should calculate
   alpha (nxn matrix), beta (n vector) and the chisq. Here n is the
   number of parameters which are fitted for. info is a (void *),
   which can point to anything that might be relevant for that
   function. It is directly passed on from the argument list of
   fit_levmar2_ld().



   Returns 0 on success.
   Returns 1 on maximum number of itterations exceeded error (try smaller ftol)
   Returns 2 on memory error
   Returns 3 if less than 1 fit parameter are not fixed
   Returns 4 if algorithm is not available  NOT RELEVANT
   Returns 5 Other fit fail error (such as singular matrix)
*/
int fit_levmar2_ld(long double *xstart, int *fixed, long double *xfit, long double *yfit, int nrparams, 
		void (*func)(long double *, long double *, long double *, long double *, void *), void *info, long double ftol, int *nfunk, int finderrors, long double sigma, long double *dplus, long double *dmin, verbose_definition verbose)
{
  int i, nrfitparameters, paramnr, maxiter, status, converged;
  long iter;
  long double *func_params, *func_params_laststep, *alpha, *beta, *alpha_tmp, *beta_tmp, chisq, alambda, tmpvalue;
  long double chisq_last_nr;


  maxiter = 5000;

  /* Determine the number total number of parameters for which we are
     fitting. */
  nrfitparameters = 0;
  for(i = 0; i < nrparams; i++) {
    if(fixed[i] == 0)
      nrfitparameters++;
  }
  if(nrfitparameters == 0) {
    printerror(verbose.debug, "ERROR fit_levmar2_ld: No parameters are set to be fitted.");
    return 3;
  }

  // Allocate memory that will fit the found parameters + errors
  func_params = malloc(nrparams*sizeof(long double));
  func_params_laststep = malloc(nrparams*sizeof(long double));
  alpha = malloc(nrfitparameters*nrfitparameters*sizeof(long double));
  beta = malloc(nrfitparameters*sizeof(long double));
  alpha_tmp = malloc(nrfitparameters*nrfitparameters*sizeof(long double));
  beta_tmp = malloc(nrfitparameters*sizeof(long double));
  if(func_params == NULL || func_params_laststep == NULL || alpha == NULL || beta == NULL || alpha_tmp == NULL || beta_tmp == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fit_levmar2_ld: Cannot allocate memory");
    return 2;
  }

  // Now fill the func_params struct with the initial values
  for(i = 0; i < nrparams; i++) {
    func_params[i] = xstart[i];
    func_params_laststep[i] = xstart[i];
  }

  // Obtain initial alpha, beta and chi2 values
  func(func_params, alpha, beta, &chisq, info);
  // Remember the initial chi square
  chisq_last_nr = chisq;

  iter = 0;
  alambda = 0.001;    // Initial value of lambda to use
  status = 0;
  do {
    iter++;  // Number of iterations in this current call to this function

    // Make some copies in case the fit did not improve after next itteration
    for (i = 0; i < nrfitparameters*nrfitparameters; i++) {
      alpha_tmp[i] = alpha[i];
    }
    for (i = 0; i < nrfitparameters; i++) {
      beta_tmp[i] = beta[i];
    }

    // A key trick of the LevMar algorithm: multiply elements on the diagonal of alpha matrix with a numerical factor
    for (i = 1; i < nrfitparameters; i++) {
      alpha[i*nrfitparameters+i] *= 1.0 + alambda;
    }

    // Think I should use LU decomp for this. It is 3X faster when the inverse is not required, and about the same to calculate the inverse. See pg 72 and earlier
    // alpha is replaced by its inverse, da (beta) is replaced by the solution vector
    if(linalg_solve_matrix_eq_gauss_jordan_ld(alpha, beta, nrfitparameters, 1, verbose) != 0) {
      printerror(verbose.debug, "ERROR fit_levmar2_ld: Solving matrix equation failed.");
      return 5;
    }

    // Update solution vector
    paramnr = 0;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
	//	printf("XXXX updating param %d from %.15Le to %.15Le\n", paramnr, func_params[i], func_params[i]+beta[paramnr]);
	func_params[i] += beta[paramnr++];
      }
    }

    // Get new alpha and beta and chisq
    func(func_params, alpha, beta, &chisq, info);

    // Need to check if new solution is acceptable
    //    printf("XXXX new chi2 = %Le (instead of %Le)\n", chisq, chisq_last_nr);
    if(chisq > chisq_last_nr) { // Chi square did not improve, so presumably at machine precission
      //      printf("XXXX new chi2 did not improve\n");
      status = -1;
      // Restore the previous solution
      memcpy(func_params, func_params_laststep, sizeof(long double)*nrfitparameters);
      for (i = 0; i < nrfitparameters*nrfitparameters; i++) {
	alpha[i] = alpha_tmp[i];
      }
      for (i = 0; i < nrfitparameters; i++) {
	beta[i] = beta_tmp[i];
      }

      alambda *= 10.1;   // Increase the value of alambda. Not by the same factor as the decrease to possibly avoid ping-ponging back-forth between solutions. Maybe this is a non-issue.
      //    }else {
      //      printf("XXXX new chi2 improved\n");
    }

    converged = 1;
    // Do a convergence test, even if there is a machine precision warning
    for(i = 0; i < nrparams; i++) {
      //	  fprintf(stderr, "XXXXXXX %lf %lf\n", func_params[i], psrsalsa_params->func_params_laststep[i]);
      if(fixed[i] == 0) {
	if(fabsl(func_params[i] - func_params_laststep[i]) >= 1.1e-19 + ftol*fabsl(func_params_laststep[i])) {
	  converged = 0;
	  break;
	}
      }
    }
    if(converged && status == -1)
      status = 0;

    //    if(converged)
    //      printf("XXXX converged, status =  %d\n", status);
    //    else
    //      printf("XXXX not converged, status = %d\n", status);

    if(status == -1) {
      printwarning(verbose.debug, "Machine precision reached, but required precision tolerance conditions are not met.");
    }

    if(chisq < chisq_last_nr) { // Chi square did improve
      // Keep a copy of the last parameters to allow future checks of convergence
      for(i = 0; i < nrparams; i++) {
	func_params_laststep[i] = func_params[i];
      }

      chisq_last_nr = chisq;
      alambda *= 0.1;   // Decrease the alambda parameter if step was successful
    }
    if(verbose.debug) {
      printf("fit_levmar2_ld: itteration %ld - status = %d\n", iter, status);
    }
  }while(status == 0 && converged == 0 && iter < maxiter);


  // Done with fitting. Now find errors
  if(finderrors && dplus != NULL && dmin != NULL) {
    for(i = 0; i < nrparams; i++) {
      if(dplus != NULL)
	dplus[i] = 0;
      if(dmin != NULL)
	dmin[i] = 0;
    }
    // Obtain the covariance matrix by inverting the current alpha matrix
    if(linalg_solve_matrix_eq_gauss_jordan_ld(alpha, beta, nrfitparameters, 1, verbose) != 0) {
      printerror(verbose.debug, "ERROR fit_levmar2_ld: Solving matrix equation failed.");
      return 5;
    }


    paramnr = 0;
    for(i = 0; i < nrparams; i++) {
      if(fixed[i] == 0) {
	tmpvalue = sqrtl(alpha[paramnr*nrfitparameters+paramnr])*sigma;
	if(dplus != NULL) {
	  dplus[i] = tmpvalue;
	}
	if(dmin != NULL) {
	  dmin[i] = tmpvalue;
	}
      }
    }

  }

  // Copy solution
  for(i = 0; i < nrparams; i++) {
    xfit[i] = func_params[i];
  }

  *yfit = chisq;
  *nfunk = iter;

  if(iter == maxiter)
    return 1;
  if(status == -1) 
    return 5;

  return 0;
}
