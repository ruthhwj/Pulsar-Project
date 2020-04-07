#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "psrsalsa.h"

// Test function designed to test the function compute_covariance_matrix(), which is not used any further in this program.
double test_dydparam(long param, double *optimum_params, double x);

int main(int argc, char **argv)
{
  int Read_Sigma, Use_Sigma, dolog, forceChi2, i, oneatatime, status, algorithm, showcovariance;
  double sig_scale, tolerance;
  long ndata;
  double *data_x, *data_y, *data_sig;
  fitfunc_collection_type function;
  verbose_definition verbose;
  int test_compute_covariance_matrix;  // If set, compute_covariance_matrix() is tested, which is not used by this program, but it allows to test if this function behaves correctly

  cleanVerboseState(&verbose);

  /* int j, k, l; */

  Read_Sigma = 0;
  sig_scale = 1;
  Use_Sigma = 1;
  dolog = 0;
  forceChi2 = 0;
  oneatatime = 1;
  test_compute_covariance_matrix = 0;
  algorithm = 1;
  showcovariance = 0;
  tolerance = 1e-4;

  if(argc < 2) {
    fprintf(stderr, "Usage: fit_mrq [options] file\n\n");
    fprintf(stderr, "-s     Read three columns from file: x,y and sigma. Default is only x,y.\n");
    fprintf(stderr, "-i     Read three columns from file: x,y and sigma, but ignore take all sigma's equal to one.\n");
    fprintf(stderr, "-a     Scale errors by this fraction.\n");
    show_fitfunctions_commandline_options(stderr, "-f", 6);
    fprintf(stderr, "-l     Take log of x and y data.\n");
    fprintf(stderr, "-1     Force reduced chi2 to be one.\n");
    fprintf(stderr, "-v     verbose.\n");
    fprintf(stderr, "-all   Do fitting for all parameter simulteneously instead of one at a time. This makes the fitting faster, but less robust.\n");
    fprintf(stderr, "-algorithm  1=GSL, 2=NR.\n");
    fprintf(stderr, "-tol  Stop itterating when changes in fit parameters are less than this fraction of its value [%.2e].\n", tolerance);
    fprintf(stderr, "-covar\n");
    return 0;
  }else if(argc > 2) {
    for(i = 1; i < argc-1; i++) {
      if(strcmp(argv[i], "-s") == 0) {
	Read_Sigma = 1;
      }else if(strcmp(argv[i], "-f") == 0) {
	/* Ignore parsing fitfunction, it will be dealt with later */
	i++;
      }else if(strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "-A") == 0) {
	sig_scale = atof(argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-algorithm") == 0) {
	algorithm = atoi(argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-tolerance") == 0 || strcmp(argv[i], "-tol") == 0) {
	tolerance = atof(argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "-I") == 0) {
	Read_Sigma = 1;
	Use_Sigma = 0;
      }else if(strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "-L") == 0) {
	dolog = 1;
      }else if(strcmp(argv[i], "-all") == 0) {
	oneatatime = 0;
      }else if(strcmp(argv[i], "-v") == 0) {
	verbose.verbose = 1;
      }else if(strcmp(argv[i], "-1") == 0) {
	forceChi2 = 1;
      }else if(strcmp(argv[i], "-test_compute_covariance_matrix") == 0) {
	test_compute_covariance_matrix = 1;
      }else if(strcmp(argv[i], "-covar") == 0) {
	showcovariance = 1;
      }
    }
  }
  if(parse_commandline_fitfunctions(argc, argv, "-f", &function, verbose) == 0)
    return 0;
  //  if(read_ascii_statdata_double(argv[argc-1], &ndata, &data_x, &data_y, &data_sig, Read_Sigma, sig_scale, dolog, dolog, 1) == 0)
  //    return 0;

  verbose_definition verbose2;
  copyVerboseState(verbose, &verbose2);
  verbose2.verbose = 1;

  int nrExpectedColumns;

  if(Read_Sigma) {
    nrExpectedColumns = 3;
  }else {
    nrExpectedColumns = 2;
  }

  if(read_ascii_column_double(argv[argc-1], 0, '#', nrExpectedColumns, 0, &ndata, 1, 1.0, dolog, &data_x, NULL, NULL, NULL, verbose2, 0) == 0)
    return 0;
  if(read_ascii_column_double(argv[argc-1], 0, '#', nrExpectedColumns, 0, &ndata, 2, 1.0, dolog, &data_y, NULL, NULL, NULL, verbose2, 0) == 0)
    return 0;
  if(Read_Sigma) {
    if(read_ascii_column_double(argv[argc-1], 0, '#', nrExpectedColumns, 0, &ndata, 3, sig_scale, dolog, &data_sig, NULL, NULL, NULL, verbose2, 0) == 0)
      return 0;
    if(Use_Sigma == 0) {
      for(i = 0; i < ndata; i++)
	data_sig[i] = 1;
    }
  }else {
    data_sig = NULL;
  }

  /*
  printf("Showing %ld points\n", ndata);
  for(i = 0; i < ndata; i++) {
    printf("%e %e %e\n", data_x[i], data_y[i], data_sig[i]);
  }
  */

  if(fit_levmar(algorithm, &function, data_x, data_y, data_sig, ndata, oneatatime, forceChi2, 0, tolerance, 500, &status, 1, showcovariance, verbose) != 0) {
    if(status) {
      printerror(verbose.debug, "ERROR fit_mrq: Fit did not converge.");
      return 0;
    }
    printwarning(verbose.debug, "WARNING fit_mrq: A warning was generated in fit_levmar().");
  }


  if(test_compute_covariance_matrix) {
    // Fit should be a + b + c*x**2
    int ok;
    ok = 1;
    if(function.nrfuncs != 3) {
      ok = 0;
      if(function.func[0].type != FUNC_POLYNOMAL || function.func[1].type != FUNC_POLYNOMAL || function.func[2].type != FUNC_POLYNOMAL)
	ok = 0;
      if(function.func[0].param[0] != 0.0 || function.func[1].param[0] != 1.0 || function.func[2].param[0] != 2.0)
	ok = 0;
    }
    if(ok == 0) {
      printerror(verbose.debug, "ERROR fit_mrq: Fit should be a + b + c*x**2 to do this test.");
      return 0;
    }
    fflush(stdout);
    printf("Fit function of correct type\n");
    double optimum_params[3];
    optimum_params[0] = function.func[0].value[0];
    optimum_params[1] = function.func[0].value[1];
    optimum_params[2] = function.func[0].value[2];
    double covariance_matrix[3*3];
    if(compute_covariance_matrix(3, optimum_params, test_dydparam, function.chi2, forceChi2, ndata, data_x, data_sig, covariance_matrix, verbose) != 0) {
      printerror(verbose.debug, "ERROR fit_mrq: Covariance matrix couldn't be determined.");
      return 0;
    }
    printf("According to compute_covariance_matrix() the errors should be (for chi2=%e):\n", function.chi2);
    printf("  %e\n", sqrt(covariance_matrix[0*3+0]));
    printf("  %e\n", sqrt(covariance_matrix[1*3+1]));
    printf("  %e\n", sqrt(covariance_matrix[2*3+2]));
  }

  free(data_x);
  free(data_y);
  if(Read_Sigma)
    free(data_sig);
  return 0;
}

double test_dydparam(long param, double *optimum_params, double x)
{
  if(param == 0) { // a -> 1.0
    return 1.0;
  }else if(param == 1) { // a*x -> x
    return x;
  }else if(param == 2) { // a*x^2 -> x^2
    return x*x;
  }
  printf("Bug!!!!\n");
  exit(0);
  return 0;
}
