#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "psrsalsa.h"


double **data;
long nr_p3_bins = -1;
long long_bin_nr;
int pol_eq_nr;
double sigma1, sigma2;

double chi2_funk(double *x)
{
  double chi2;
  double denomenator;
  chi2 = 0;


  //  printf("Try out parameter: %e %e", x[0], x[1]);

  if(pol_eq_nr == 0) {
    denomenator = sigma1*sigma1*(1+x[0]*x[0]) + x[1]*x[1]*sigma2*sigma2;
  }else {
    denomenator = x[0]*x[0]*sigma1*sigma1 + sigma2*sigma2*(1+x[1]*x[1]);
  }


  int i;
  for(i = 0; i < nr_p3_bins; i++) {
    double term;
    double Ol1, Ol2, Ot1, Ot2;
    Ol1 = (data[long_bin_nr])[i+0*nr_p3_bins];
    Ol2 = (data[long_bin_nr])[i+2*nr_p3_bins];
    Ot1 = (data[long_bin_nr])[i+1*nr_p3_bins];
    Ot2 = (data[long_bin_nr])[i+3*nr_p3_bins];
    if(pol_eq_nr == 0) {
      term = Ol1-x[0]*Ot1-x[1]*Ot2;
    }else {
      term = Ol2-x[0]*Ot1-x[1]*Ot2;
    }
    chi2 += term*term/denomenator;
  }

  //  printf("  - got chi2 = %e\n", chi2);

  return chi2;
}

int main(int argc, char **argv)
{
  long i;
  psrsalsaApplication application;

  sigma1 = sigma2 = -1;

  initApplication(&application, "fit_asym", "[options] inputfile1 inputfile2");
  application.switch_verbose = 1;
  application.switch_debug = 1;

  if(argc < 2) {
    printf("Program to fit the Asymmetry matrix, to do Crispin's 0031 stuff.\n");
    printApplicationHelp(&application);
    printf("Input options:\n\n");
    printf("-sigma   \"sigma1 sigma2\"    Define the two RMSes of the two polarizations.\n");
    printf("-np3bins n                    Set the expected the number of P3 bins to be n.\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {

    for(i = 1; i < argc; i++) {
      int index;
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-sigma") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &sigma1, &sigma2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fit_asym: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-np3bins") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &nr_p3_bins, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fit_asym: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else {
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "ERROR fit_asym: Unknown option: %s\n\nRun fit_asym without command line arguments to show help", argv[i]);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }

  if(sigma1 <= 0 || sigma2 <= 0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Please use the -sigma option to provide two positive numbers.");
    return 0;
  }
  if(nr_p3_bins <= 0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Please use the -np3bins option to provide a positive number.");
    return 0;
  }
  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) != 2) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Two input files were expected at the end of the command line.");
    return 0;
  }

  long nr_longitude_bins;
  double *guess1, *guess2, *guess3, *guess4;
  char *filename_ptr;


  filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state);
  if(application.verbose_state.verbose) {
    fprintf(stdout, "Loading initial guesses from ascii file %s\n", filename_ptr);
  }
  // Read first column
  if(read_ascii_column_double(filename_ptr, 0, '#', 4, 0, &nr_longitude_bins, 1, 1.0, 0, &guess1, NULL, NULL, NULL, application.verbose_state, 1) == 0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: cannot load file.\n");
    return 0;
  }
  if(read_ascii_column_double(filename_ptr, 0, '#', 4, 0, &i, 2, 1.0, 0, &guess2, NULL, NULL, NULL, application.verbose_state, 1) == 0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: cannot load file.\n");
    return 0;
  }
  if(i != nr_longitude_bins) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Each column in the input data should have an equal number of (non-empty) rows.\n");
    return 0;
  }
  if(read_ascii_column_double(filename_ptr, 0, '#', 4, 0, &i, 3, 1.0, 0, &guess3, NULL, NULL, NULL, application.verbose_state, 1) == 0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: cannot load file.\n");
    return 0;
  }
  if(i != nr_longitude_bins) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Each column in the input data should have an equal number of (non-empty) rows.\n");
    return 0;
  }
  if(read_ascii_column_double(filename_ptr, 0, '#', 4, 0, &i, 4, 1.0, 0, &guess4, NULL, NULL, NULL, application.verbose_state, 1) == 0) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: cannot load file.\n");
    return 0;
  }
  if(i != nr_longitude_bins) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Each column in the input data should have an equal number of (non-empty) rows.\n");
    return 0;
  }

  if(application.verbose_state.verbose) {
    fprintf(stdout, "  Found %ld longitude bins with initial guesses in %s\n", nr_longitude_bins, filename_ptr);
  }

  /*
  for(i = 0; i < nr_longitude_bins; i++) {
    printf("%e %e %e %e\n", guess1[i], guess2[i], guess3[i], guess4[i]);
  }
  return 0;
  */


  data = malloc(nr_longitude_bins*sizeof(double *));
  if(data == NULL) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: Memory allocation error.\n");
    return 0;
  }
  for(long_bin_nr = 0; long_bin_nr < nr_longitude_bins; long_bin_nr++) {
    data[long_bin_nr] = malloc(4*nr_p3_bins * sizeof(double));
    if(data[long_bin_nr] == NULL) {
      printerror(application.verbose_state.debug, "ERROR fit_asym: Memory allocation error.\n");
      return 0;
    }
  }
  filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state);
  if(application.verbose_state.verbose) {
    fprintf(stdout, "Loading data from ascii file %s\n", filename_ptr);
  }

  FILE *fin;
  fin = fopen(filename_ptr, "r");
  if(fin == NULL) {
    printerror(application.verbose_state.debug, "ERROR fit_asym: cannot open file %s.\n", filename_ptr);
    return 0;
  }
  long row_nr;
  long colnr;
  for(row_nr = 0; row_nr < 4*nr_p3_bins; row_nr++) {
    for(colnr = 0; colnr < nr_longitude_bins; colnr++) {
      int ret = fscanf(fin, "%lf", &((data[colnr])[row_nr]));
      if(ret != 1) {
	printerror(application.verbose_state.debug, "ERROR fit_asym: cannot read file %s. Failed when reading in column %ld, row %ld\n", filename_ptr, colnr, row_nr);
	return 0;
      }
    }
  }
  fclose(fin);

  /*
  for(row_nr = 0; row_nr < 4*nr_p3_bins; row_nr++) {
    for(i = 0; i < nr_longitude_bins; i++) {
      printf("%e ", (data[i])[row_nr]);
    }
    printf("\n");
  }
  return 0;
  */



  for(long_bin_nr = 0; long_bin_nr < nr_longitude_bins; long_bin_nr++) {
    double  chi2, chi2_1;
    chi2_1 = 0;
    for(pol_eq_nr = 0; pol_eq_nr < 2; pol_eq_nr++) {

      double xstart[2];
      double dx[2];
      int fixed[2];
      double solution[2];
      int nritt;

      if(pol_eq_nr == 0) {
	xstart[0] = guess1[long_bin_nr];
	xstart[1] = guess2[long_bin_nr];
      }else {
	xstart[0] = guess3[long_bin_nr];
	xstart[1] = guess4[long_bin_nr];
      }
      dx[0] = 0.1;
      dx[1] = 0.1;
      fixed[0] = 0;
      fixed[1] = 0;
      
      int ret;
      ret = doAmoeba_d(0, xstart, dx, fixed, solution, &chi2, 2, &chi2_funk, 1e-6, &nritt, application.verbose_state.verbose, 0, 0.0, NULL, NULL);
      if(ret == 1) {
	printerror(application.verbose_state.debug, "Error fit_asym: Downhill-Simplex method did not converge. You can try lowering the tolerance with -ftol.");
	return 0;
      }else if(ret != 0) {
	printerror(application.verbose_state.debug, "Error fit_asym: Downhill-Simplex method failed.");
	return 0;
      }

      if(application.verbose_state.verbose) {
	printf("After %d steps the down-hill simplex found the best solution with a chi2 = %e\n\n", nritt, chi2);
      }

      printf("%e %e ", solution[0], solution[1]);
      if(pol_eq_nr == 0) {
	chi2_1 = chi2;
      }
    }
    printf("%e %e\n", chi2_1, chi2);
  }



  return 0;
}
