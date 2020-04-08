//START REGION RELEASE
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "psrsalsa.h"

/*
#include <stdlib.h>
#include <malloc.h>
#include "nr.h"
#include "nrutil.h"
*/


//START REGION DEVELOP
void SHOWREVISIONINFO_prog() {
#include "pdist.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);
//START REGION RELEASE



/*
float makeenergydist(int ft, float x1, float x2, long npoints, char *outputname, int appendflag) 
{
  char cmdline[1000];
  char cmdline2[1000];
  if(ft == 1) {
    printf("m=%f s=%f\n", x1, x2);
    // Prevent amouba to fit negative sigma 
    if(x2 < 0)
      return 1e10;
    sprintf(cmdline, "fake_distribution -N %ld %s -l \"%f %f\" -S measured_off.en -s 1  -r \"%f %f %ld\" ", npoints, extraoptions, x1, x2, min_e, max_e, nrebins);
  }else if(ft == 2) {
    printf("p=%f c=%f\n", x1, x2);
    sprintf(cmdline, "fake_distribution -N %ld %s -p \"1 %f 0\" -S measured_off.en -s 1  -r \"%f %f %ld\" ", npoints, extraoptions, x1, x2, max_e, nrebins);
  }else if(ft == 4) {
    printf("s=%f\n", x1);
    sprintf(cmdline, "fake_distribution -N %ld %s -R \"%f\" -S measured_off.en -s 1  -r \"%f %f %ld\" ", npoints, extraoptions, x1, min_e, max_e, nrebins);
  }else if(ft == 5) {
    printf("k=%f theta=%f\n", x1, x2);
    // Prevent amouba to fit negative sigma
    if(x2 < 0)
      return 1e10;
    sprintf(cmdline, "fake_distribution -N %ld %s -G \"%f %f\" -S measured_off.en -s 1  -r \"%f %f %ld\" ", npoints, extraoptions, x1, x2, min_e, max_e, nrebins);
  }else if(ft == 11) {
    printf("m=%f s=%f\n", x1, x2);
    // Prevent amouba to fit negative sigma
    if(x2 < 0)
      return 1e10;
    sprintf(cmdline, "fake_distribution -N %ld %s -l \"%f %f\" -S measured_off.en -s 1  -r \"%f %f %ld\" ", npoints, extraoptions, x1, x2, min_e, max_e, nrebins);
  }else if(ft == 13) {
    printf("m=%f s=%f\n", x1, x2);
    // Prevent amouba to fit negative sigma
    if(x2 < 0)
      return 1e10;
    sprintf(cmdline, "fake_distribution -N %ld %s -S measured_off.en -s 1 -n \"%f %f\" -r \"%f %f %ld\" ", npoints, extraoptions, x1, x2, min_e, max_e, nrebins);
  }
  if(verbose)
    strcat(cmdline, " -v ");
  if(appendflag)
    strcat(cmdline, " >> ");
  else
    strcat(cmdline, " > ");
  strcat(cmdline, outputname);
  strcat(cmdline, "\n");
  if(verbose) printf("  %s", cmdline);
  if(!appendflag) {
    sprintf(cmdline2, "rm -f %s", outputname);
    system(cmdline2);
  }
  system(cmdline);
  return -1;
}
*/


struct {
  int fittype[2];
  char *cmdline, *txt, *measurement_file, *threshold_values, *pdist_options, *noisefile_id;
  long nr_model_points, seednr;
  int debug, fixseed, colspecified, file_column1, nocounters;
  int method, nrfunctions, frac_paramnr, null_distr_specified;
  double dx, sigma;
}fitter_info;

//return 1 = ok, 0 means parameters make no sense
int make_fakeDist_cmd(double *x, char *filename)
{
  long int n, funcnr, paramnr;
  n = fitter_info.nr_model_points;
  if(fitter_info.nrfunctions == 2)
    n = round((double)n*(1.0-x[fitter_info.frac_paramnr]));    // Calculate nr of points in first distribution
  sprintf(fitter_info.cmdline, "fakeDist -N %ld ", n);
  n = fitter_info.nr_model_points - n;                         // Calculate nr of points in second distribution
  if(filename != NULL) {
    sprintf(fitter_info.txt, "-output model_trial_tmp.dist ");
    strcat(fitter_info.cmdline, fitter_info.txt);
  }
  if(fitter_info.fixseed) {
    strcat(fitter_info.cmdline, "-fixseed ");
  }else {
    sprintf(fitter_info.txt, "-seed %ld ", fitter_info.seednr);
    strcat(fitter_info.cmdline, fitter_info.txt);
  }
  if(fitter_info.sigma > 0) {
    sprintf(fitter_info.txt, "-sigma %e ", fitter_info.sigma);
    strcat(fitter_info.cmdline, fitter_info.txt);
  }
  if(fitter_info.noisefile_id != NULL) {
    sprintf(fitter_info.txt, "-noisefile %s ", fitter_info.noisefile_id);
    strcat(fitter_info.cmdline, fitter_info.txt);
  }
  paramnr = 0;
  for(funcnr = 0; funcnr < fitter_info.nrfunctions; funcnr++) {
    if(funcnr == 1) {
      sprintf(fitter_info.txt, " -N2 %ld ", n);
      strcat(fitter_info.cmdline, fitter_info.txt);
      paramnr++;
    }
    if(fitter_info.fittype[funcnr] == 1) {  // Gamma distribution
      sprintf(fitter_info.txt, "-gamma '%e %e'", x[paramnr], x[paramnr+1]);
      if(x[paramnr] <= 0)
	return 0;
      if(x[paramnr+1] <= 0)
	return 0;
      paramnr += 2;
    }else if(fitter_info.fittype[funcnr] == 2) {  // Flat distribution
      sprintf(fitter_info.txt, "-flat '%e %e'", x[paramnr], x[paramnr+1]);
      if(x[paramnr] >= x[paramnr+1])
	return 0;
      paramnr += 2;
    }else if(fitter_info.fittype[funcnr] == 3) {  // Normal distribution
      sprintf(fitter_info.txt, "-norm '%e %e'", x[paramnr], x[paramnr+1]);
      //      printf("XXXX Normal distr: %e %e\n", x[paramnr], x[paramnr+1]);
      if(x[paramnr+1] <= 0)
	return 0;
      paramnr += 2;
    }else if(fitter_info.fittype[funcnr] == 4) {  // Lognormal distribution
      sprintf(fitter_info.txt, "-lognorm '%e %e'", x[paramnr], x[paramnr+1]);
      if(x[paramnr+1] <= 0)
	return 0;
      paramnr += 2;
    }else if(fitter_info.fittype[funcnr] == 5) {  // Powerlaw distribution
      sprintf(fitter_info.txt, "-pwrlaw '%e %e'", x[paramnr], x[paramnr+1]);
      if(x[paramnr] > 0)
	return 0;
      paramnr += 2;
    }else if(fitter_info.fittype[funcnr] == 6) {  // Rayleigh distribution
      sprintf(fitter_info.txt, "-Rayleigh '%e'", x[paramnr]);
      //      printf("XXXX Rayleigh: %e\n", x[paramnr]);
      if(x[paramnr] <= 0)
	return 0;
      paramnr++;
    }else {
      printerror(fitter_info.debug, "ERROR pdistFit: Bug.");
      exit(0);
    }
    strcat(fitter_info.cmdline, fitter_info.txt);
  }
  if(fitter_info.null_distr_specified) {
    sprintf(fitter_info.txt, " -null %e -quiet", x[paramnr]);
    paramnr++;
    strcat(fitter_info.cmdline, fitter_info.txt);
  }

  return 1;
}

//return 1 = ok, 0 means no binning was required
//If output_filename != NULL, it is used rather than the default extension
// If trial is set, the data is always assumed to be in the first column
int make_pdist_cmd(char *filename, char *output_filename, int trial)
{
  if(fitter_info.method != 2)
    return 0;                    // No binning required
  sprintf(fitter_info.cmdline, "pdist -dx %e -frac -sigma ", fitter_info.dx);
  if(fitter_info.colspecified && trial == 0) {
    sprintf(fitter_info.txt, "-col %d ", fitter_info.file_column1);
    strcat(fitter_info.cmdline, fitter_info.txt);
  }
  if(fitter_info.pdist_options != NULL) {
    strcat(fitter_info.cmdline, fitter_info.pdist_options);
    strcat(fitter_info.cmdline, " ");
  }
  //  if(fitter_info.debug) {
  //    sprintf(fitter_info.txt, "-v ");
  //    strcat(fitter_info.cmdline, fitter_info.txt);
  //  }
  if(output_filename != NULL) {
    sprintf(fitter_info.txt, "-output %s ", output_filename);
    strcat(fitter_info.cmdline, fitter_info.txt);
  }
  strcat(fitter_info.cmdline, filename);
  return 1;
}

double funk(double *x)
{
  FILE *fin;
  static long trial_nr = 0;
  double teststat;
  static double lastteststat = 0;
  int ret, found, rebinning;

  // See if temporary files are not there
  found = 0;
  rebinning = 0;
  if(access("model_trial_tmp.dist", F_OK) == 0) {
    printerror(fitter_info.debug, "ERROR pdistFit: File model_trial_tmp.dist already exist. Remove or move this file first, as it will be overwritten otherwise.");
    found = 1;
  }
  if(make_pdist_cmd("model_trial_tmp.dist", NULL, 1) == 1) { // See if rebinning is required
    rebinning = 1;
    if(access("model_trial_tmp.dist.hist", F_OK) == 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: File model_trial_tmp.dist.hist already exist. Remove or move this file first, as it will be overwritten otherwise.");
      found = 1;
    }
  }
  if(access("teststat_trial_tmp.txt", F_OK) == 0) {
    printerror(fitter_info.debug, "ERROR pdistFit: File teststat_trial_tmp.txt already exist. Remove or move this file first, as it will be overwritten otherwise.");
    found = 1;
  }
  if(found) {
    exit(0);
  }

  if(fitter_info.debug) {
    printf("This is trial %ld\n", trial_nr+1);
  }
  // Generate command line to give modelled distribution
  if(make_fakeDist_cmd(x, "model_trial_tmp.dist") == 0) { // If parameters make no sense
    if(trial_nr == 0) {
      printf("  Rejecting input parameters.\n");
      return 1e10;
    }else {
      printf("  Rejecting input parameters.\n");
      return 1e10*lastteststat;
    }
  }
  if(fitter_info.debug) {
    printf("  Executing: %s\n", fitter_info.cmdline);
  }
  fflush(stdout);
  system(fitter_info.cmdline);
  trial_nr++;

  if(make_pdist_cmd("model_trial_tmp.dist", NULL, 1) == 1) { // See if rebinning is required
    if(fitter_info.debug == 0)
      strcat(fitter_info.cmdline, " > /dev/null");
    if(fitter_info.debug) {
      printf("  Executing: %s\n", fitter_info.cmdline);
    }
    fflush(stdout);
    system(fitter_info.cmdline);
  }

  // Get test statistic
  if(fitter_info.method == 0) {
    sprintf(fitter_info.cmdline, "pstat -chi2cdf ");
  }else if(fitter_info.method == 1) {
    sprintf(fitter_info.cmdline, "pstat -ks -v "); // verbose is needed, because difference in cdf is used rather than probability
  }else if(fitter_info.method == 2) {
    sprintf(fitter_info.cmdline, "pstat -chi2hist "); 
    if(fitter_info.threshold_values == NULL) {
      strcat(fitter_info.cmdline, "-1 ");
    }else {
      sprintf(fitter_info.txt, "\"%s\" ", fitter_info.threshold_values);
      strcat(fitter_info.cmdline, fitter_info.txt);
    }
  }else {
    printerror(fitter_info.debug, "ERROR pdistFit: Bug.");
    exit(0);
  }
  if(rebinning) { // See if rebinning is required
    sprintf(fitter_info.txt, "-col1 '1 2 3' -col2 '1 2 3' measurement_tmp.hist model_trial_tmp.dist.hist ");
  }else {
    if(fitter_info.colspecified) {
      sprintf(fitter_info.txt, "-col1 %d ", fitter_info.file_column1);
      strcat(fitter_info.cmdline, fitter_info.txt);
      //    }else if(fitter_info.polspecified) {
      //      sprintf(fitter_info.txt, "-pol1 %d ", fitter_info.file_column1);
      //      strcat(fitter_info.cmdline, fitter_info.txt);
    }
    sprintf(fitter_info.txt, "-col2 1 %s model_trial_tmp.dist ", fitter_info.measurement_file);
  }
  strcat(fitter_info.cmdline, fitter_info.txt);
  if(fitter_info.method == 1) {
    sprintf(fitter_info.txt, "| grep 'KS-test statistic max_diff' "); // For KS-test, filter out line with KS statistic
    strcat(fitter_info.cmdline, fitter_info.txt);
  }
  sprintf(fitter_info.txt, "> teststat_trial_tmp.txt");
  strcat(fitter_info.cmdline, fitter_info.txt);
  if(fitter_info.debug) {
    printf("  Executing: %s\n", fitter_info.cmdline);
  }
  fflush(stdout);
  system(fitter_info.cmdline);
  if(fitter_info.debug) {
    printf("    ");
    fflush(stdout);
    system("cat teststat_trial_tmp.txt");
  }

  // Read in test statistic
  fin = fopen("teststat_trial_tmp.txt", "r");
  if(fin == NULL) {
    printerror(fitter_info.debug, "ERROR pdistFit: A file teststat_trial_tmp.txt should be generated, but it cannot be opened.");
    exit(0);
  }
  if(fitter_info.method == 0) {
    //Non-weighted total chi square = 417.378262 = 4.173783e+02
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "Non-weighted") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "total") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "chi") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "square") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "=") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    fscanf(fin, "%s", fitter_info.txt);
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "=") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    ret = fscanf(fin, "%lf", &teststat);
    if(ret != 1) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
  }else if(fitter_info.method == 1) {
    //    KS-test statistic max_diff: 0.443609 = 4.436092e-01
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "KS-test") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "statistic") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "max_diff:") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    fscanf(fin, "%s", fitter_info.txt);
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "=") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    ret = fscanf(fin, "%lf", &teststat);
    if(ret != 1) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
  }else if(fitter_info.method == 2) {
    //    Reduced chi square: 38.604024 = 3.860402e+01
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "Reduced") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "chi") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "square:") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    fscanf(fin, "%s", fitter_info.txt);
    fscanf(fin, "%s", fitter_info.txt);
    if(strcmp(fitter_info.txt, "=") != 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
    ret = fscanf(fin, "%lf", &teststat);
    if(ret != 1) {
      printerror(fitter_info.debug, "ERROR pdistFit: Cannot interpret the file teststat_trial_tmp.txt, so something is going wrong.");
      exit(0);
    }
  }else {
    printerror(fitter_info.debug, "ERROR pdistFit: Bug.");
    exit(0);
  }
  if(fitter_info.debug) {
    printf("  Test statistic: %e\n", teststat);
  }

  // Clean up temporary files
  sprintf(fitter_info.cmdline, "rm model_trial_tmp.dist teststat_trial_tmp.txt");
  if(rebinning) { // See if rebinning is required
    strcat(fitter_info.cmdline, " model_trial_tmp.dist.hist");
  }
  if(fitter_info.debug) {
    printf("  Executing: %s\n", fitter_info.cmdline);
  }
  fflush(stdout);
  system(fitter_info.cmdline);

  //  exit(0);

  lastteststat = teststat; // Remember for future use
  if(fitter_info.nocounters == 0) {
    printf("\rThis is trial %ld: teststat = %e              ", trial_nr, teststat);
  }
  return teststat;
}

int main(int argc, char **argv)
{
  int plotcdf, second_distr_specified, polspecified, remove_tmp_dist;
  long i, j;
  double function_param[2][4], ftol, second_distr_frac, second_distr_dfrac, null_distr_av, null_distr_dav;
  int function_param_fixed[4][4];
  psrsalsaApplication application;

//START REGION RELEASE
  initApplication(&application, "pdistFit", "[options] inputfile");
//START REGION DEVELOP
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
//START REGION RELEASE

  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_iformat = 1;
  application.switch_formatlist = 1;
  application.switch_fixseed = 1;
  application.switch_nocounters = 1;
  application.switch_libversions = 1;


  fitter_info.colspecified = 0;
  fitter_info.file_column1 = 1; // Default columns to read in (in ascii mode)
  ftol = 1e-3;
  plotcdf = 0;
  second_distr_specified = 0;
  polspecified = 0;
  remove_tmp_dist = 0;
  fitter_info.nr_model_points = 10000;
  fitter_info.debug = 0;
  fitter_info.fixseed = 0;
  fitter_info.method = 0;
  fitter_info.nocounters = 0;
  fitter_info.threshold_values = NULL;
  fitter_info.pdist_options = NULL;
  fitter_info.nrfunctions = 0;
  fitter_info.frac_paramnr = 0;
  fitter_info.sigma = -1;
  fitter_info.noisefile_id = NULL;
  fitter_info.null_distr_specified = 0;
  for(i = 0; i < 4; i++) {
    for(j = 0; j < 4; j++) {
      function_param_fixed[j][i] = 0;
    }
  }

  if(argc < 2) {
    printf("Program to fit a measured distribution with a model distribution. The input is\n");
    printf("an unbinned list of values. The parameters of the model distribution (of which a\n");
    printf("realisation is generated with the program fakeDist) are optimised using a\n");
    printf("down-hill simplex search for a minimum in a test statistic of which different\n");
    printf("options are available. The default test statistic is obtained by a system call\n");
    printf("to 'pstat -chi2cdf'. Usage:\n\n");
    printApplicationHelp(&application);
    printf("Input options:\n\n");
    printf("-col nr          Specify the column number (counting from 1) which contains\n");
    printf("                 the (not binned) values of the distribution to be fitted.\n");
    printf("                 This option implies the input files are ascii files with each\n");
    printf("                 line having an equal number of columns. Lines starting with a #\n");
    printf("                 will be ignored. The default is 1.\n");
    printf("-pol nr          This option implies that the input file should be recognized\n");
    printf("                 as a pulsar format. Polarization nr, counting from zero, is\n");
    printf("                 used (all bins, freqs and subints). The default is 0, or the\n");
    printf("                 integrated onpulse energy for penergy output files in mode 1.\n");
    printf("Test statistic:\n\n");
    printf("-chi2hist dx     Make a binned histogram with bin width dx and the chi square\n");
    printf("                 obtained via pstat -chi2hist is minimised rather than the\n");
    printf("                 default test statistic. This is a weighted chi square test.\n");
    printf("-chi2hist_thresh \"t1 t2 t3\"  The up to three threshold values to be passed on\n");
    printf("                 to pstat -chi2hist (see help pstat for a description).\n");
    printf("                 The default is \"-1\". Note that the histograms are not\n");
    printf("                 increasing in steps of 1, but fractions 1/N.\n");
    printf("-chi2hist_opt    \".....\"  Extra options to be passed on to pdist for\n");
    printf("                 generation of histograms (for instance -zero).\n");
    printf("-ks              Use pdist -ks to obtain the test statistic rather than the\n");
    printf("                 default test statistic.\n");
    printf("Precision:\n\n");
    printf("-N nr            Generate nr model distribution values for each trial.\n");
    printf("                 Default is %ld.\n", fitter_info.nr_model_points);
    printf("-tol             Set tolerance for fitting (default is %.2e)\n", ftol); 
    printf("\nOther actions:\n\n");
    printf("-plotcdf         Plot the resulting cdf of the fit and data.\n"); 
    printf("\nFit functions (run fakeDist with -v to obtain a mathematical description):\n\n");
    printf("-2nd             \"frac dfrac\". The distribution is the combination of that\n");
    printf("                 specified before and after this command-line option. The\n");
    printf("                 values are drawn either from one or the other distribution\n");
    printf("                 with a probability frac that the second distribution is used.\n");
    printf("                 This fraction is a fit parameter and dfrac is the initial\n");
    printf("                 step-size used. See below for an example using -2nd.\n");
    printf("-flat            \"min dmin max dmax\", specify a flat/uniform distribution\n");
    printf("                 between min and max. dmin and dmax are the initial step sizes\n");
    printf("                 of optimisation.\n");
    printf("-gamma           \"k dk theta dtheta\", specify gamma function gamma(k,theta),\n");
    printf("                 and dk,dtheta are the initial step sizes of optimisation.\n");
    printf("-norm            \"mu dmu sigma dsigma\", specify normal (Gaussian) function\n");
    printf("                 with mean mu and standard devitation sigma. dmu and dsigma\n");
    printf("                 are the initial step sizes to be used in optimisation.\n");
    printf("-lognorm         \"mu dmu sigma dsigma\", specify lognormal function dmu and\n");
    printf("                 dsigma are initial step sizes to be used in optimisation.\n");
    printf("-null            \"av dav\". Add zeroes to the distribution to make average equal\n");
    printf("                 to av. Hence the nr of generated values is larger than what is\n");
    printf("                 specied with -N. av is fitted for with dav being the initial\n");
    printf("                 step size. Since the average can be measured, in general you\n");
    printf("                 want to fix that parameter by adding a third parameter: see\n");
    printf("                 below for .\n");
    printf("-pwrlaw          \"idx didx min dmix\", specify function f(x)=x**idx for x>=min,\n");
    printf("                 where idx should be a negative number. didx and dmix are\n");
    printf("                 the initial step sizes.\n");
    printf("-Rayleigh        \"sigma dsigma\", specify Rayleigh function, where dsigma is the\n");
    printf("                 initial step size.\n");
    printf("-sigma s         Specified distribution(s) are convolved with a Gaussian with\n");
    printf("                 sigma s. This is not a fit-parameter but a given number.\n");
    printf("-noisefile       Use specified file (list of values) instead of Gaussian noise.\n");
    printf("                 The values are read in from the first column.\n");
    printf("\n");
    printf("Example of a combined distribution fit:\n");
    printf("  pdistFit -norm \"10 1 3 1\" -2nd \"0.5 0.1\" -norm \"20 1 3 1\" inputfile\n");
    printf("Example of how to fix fit parameters:\n");
    printf("  pdistFit -norm \"10.0 1.0 0 3.0 0.0 1\" -sigma 0.1 -null \"1.0 0.1\" inputfile\n");
    printf("\nIn the example above, mu is fitted for (3rd param = 0 = not fixed), while sigma\n");
    printf("is fixed for (6th param = 1 = fixed). For all fit parameters (including -2nd and\n");
    printf("-null, it is optional to add these extra 'fixed' flags. Note that there should\n");
    printf("be at least two free parameters for the fitting to work, so in this example the\n");
    printf("average in the -null option is kept as a free parameter.\n");
    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about fitting distributions (in the context of pulse energies) can be found in:\n");
    printf(" - Weltevrede et al. 2006, A&A, 458, 269\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;

    /*
    fprintf(stdout, "Usage: fit_dist options\n\n");
    fprintf(stdout, "-f        type of fit\n");
    fprintf(stdout, "            1 = lognormal  (set m and s)\n");
    fprintf(stdout, "            2 = powerlaw   (set p and s)\n");
    fprintf(stdout, "            4 = Rayleigh   (set s)\n");
    fprintf(stdout, "            5 = Gamma      (set m (k) and s (theta))\n");
    fprintf(stdout, "            11 = lognormal (set m and s, fitting uses cdf)\n");
    fprintf(stdout, "            13 = normal    (set m and s, fitting uses cdf)\n");
    fprintf(stdout, "-f2       type of second fit, also specify relative strength and search range. Use -m2 etc. to define function. Do not combine cdf/non-cdf functions.\n");
    fprintf(stdout, "-fix      Fix value of this parameter value (can use multiple times).\n");
    fprintf(stdout, "-m        start value and range of m (lognormal).\n");
    fprintf(stdout, "-s        start value and range of s (lognormal).\n"); 
    fprintf(stdout, "-p        start value and range of p (powerlaw index).\n");
    fprintf(stdout, "-c        start value and range of c (powerlaw cutoff).\n");  
    fprintf(stdout, "-dx       binning dx (default is %f), not used in cdf mode.\n", binning_dx); 
    fprintf(stdout, "-emax_mod Maximum energy of model distribution [%f].\n", max_e);
    fprintf(stdout, "-emin_mod Minumum energy of model distribution [%f].\n", min_e);
    fprintf(stdout, "-ne_mod   Number energy bins of model distribution [%ld].\n", nrebins);
    fprintf(stdout, "-null     enable additional null distribution which makes the average energy this value.\n"); 
    fprintf(stdout, "-t        set threshold value for counts in measurement bins to accept for fitting (default is %lf). Not used in cdf mode\n", threshold); 
    fprintf(stdout, "-v        enable verbose.\n"); 
    fprintf(stdout, "\n\nYou should provide a measured_on.en and a measured_off.en (just a list of numbers, not the penergy output directly)\n\n");
    return 0;
    */
  }else {
    for(i = 1; i < argc; i++) {
      int index;
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-col") == 0 || strcmp(argv[i], "-col1") == 0 || strcmp(argv[i], "-pol") == 0) {
	if(strcmp(argv[i], "-pol") == 0)
	  polspecified = 1;
	else
	  fitter_info.colspecified = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &fitter_info.file_column1, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-gamma") == 0 || strcmp(argv[i], "-flat") == 0 || strcmp(argv[i], "-norm") == 0 || strcmp(argv[i], "-lognorm") == 0 || strcmp(argv[i], "-pwrlaw") == 0) {
	if(fitter_info.nrfunctions == 1 && second_distr_specified == 0) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: No more than one fit function can be specified without the -2nd option separating the two distributions functions.");
	  return 0;
	}
	if(fitter_info.nrfunctions == 2) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: No more than two fit function can be specified.");
	  return 0;
	}
	// Check if 6 parameters can be read in
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 1, -1, "%lf %lf %d %lf %lf %d", &function_param[fitter_info.nrfunctions][0], &function_param[fitter_info.nrfunctions][1], &function_param_fixed[fitter_info.nrfunctions][0], &function_param[fitter_info.nrfunctions][2], &function_param[fitter_info.nrfunctions][3], &function_param_fixed[fitter_info.nrfunctions][1], NULL) == 6) {
	  // Read in 6 parameters
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %d %lf %lf %d", &function_param[fitter_info.nrfunctions][0], &function_param[fitter_info.nrfunctions][1], &function_param_fixed[fitter_info.nrfunctions][0], &function_param[fitter_info.nrfunctions][2], &function_param[fitter_info.nrfunctions][3], &function_param_fixed[fitter_info.nrfunctions][1], NULL) < 6) {
	    printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	// Check if 4 parameters can be read in
	}else if(parse_command_string(application.verbose_state, argc, argv, i+1, 1, -1, "%lf %lf %lf %lf", &function_param[fitter_info.nrfunctions][0], &function_param[fitter_info.nrfunctions][1], &function_param[fitter_info.nrfunctions][2], &function_param[fitter_info.nrfunctions][3], NULL) == 4) {
	  function_param_fixed[fitter_info.nrfunctions][0] = 0;
	  function_param_fixed[fitter_info.nrfunctions][1] = 0;
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %lf %lf", &function_param[fitter_info.nrfunctions][0], &function_param[fitter_info.nrfunctions][1], &function_param[fitter_info.nrfunctions][2], &function_param[fitter_info.nrfunctions][3], NULL) < 4) {
	    printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	}else {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse the %s option. Needs to have four or six values.\n", argv[i]);
	  return 0;
	}
	if(strcmp(argv[i], "-gamma") == 0) {
	  fitter_info.fittype[fitter_info.nrfunctions] = 1;
	}else if(strcmp(argv[i], "-flat") == 0) {
	  fitter_info.fittype[fitter_info.nrfunctions] = 2;
	}else if(strcmp(argv[i], "-norm") == 0) {
	  fitter_info.fittype[fitter_info.nrfunctions] = 3;
	}else if(strcmp(argv[i], "-lognorm") == 0) {
	  fitter_info.fittype[fitter_info.nrfunctions] = 4;
	}else if(strcmp(argv[i], "-pwrlaw") == 0) {
	  fitter_info.fittype[fitter_info.nrfunctions] = 5;
	}else {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Bug.");
	  return 0;
	}
	fitter_info.nrfunctions++;
        i++;
      }else if(strcasecmp(argv[i], "-Rayleigh") == 0) {
	if(fitter_info.nrfunctions == 1 && second_distr_specified == 0) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: No more than one fit function can be specified without the -2nd option separating the two distributions functions.");
	  return 0;
	}
	if(fitter_info.nrfunctions == 2) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: No more than two fit function can be specified.");
	  return 0;
	}
	// Check if 3 parameters can be read in
	int ret;
	ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 2, "%lf %lf %d", &function_param[fitter_info.nrfunctions][0], &function_param[fitter_info.nrfunctions][1], &function_param_fixed[fitter_info.nrfunctions][0], NULL);
	if(ret == 2) {
	  function_param_fixed[fitter_info.nrfunctions][0] = 0;
	}else if(ret != 3) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse the %s option. Needs to have two or three values.\n", argv[i]);
	  terminateApplication(&application);
	  return 0;
	}
	/*
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 1, -1, "%lf %lf %d", &function_param[fitter_info.nrfunctions][0], &function_param[fitter_info.nrfunctions][1], &function_param_fixed[fitter_info.nrfunctions][0], NULL) == 0) {
	  // Read in 3 parameters
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %d", &function_param[fitter_info.nrfunctions][0], &function_param[fitter_info.nrfunctions][1], &function_param_fixed[fitter_info.nrfunctions][0], NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	// Check if 2 parameters can be read in
	}else if(parse_command_string(application.verbose_state, argc, argv, i+1, 1, -1, "%lf %lf", &function_param[fitter_info.nrfunctions][0], &function_param[fitter_info.nrfunctions][1], NULL) == 0) {
	  function_param_fixed[fitter_info.nrfunctions][0] = 0;
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &function_param[fitter_info.nrfunctions][0], &function_param[fitter_info.nrfunctions][1], NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	}else {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse the %s option. Needs to have two or three values.\n", argv[i]);
	  return 0;
	}
	*/
	if(strcasecmp(argv[i], "-Rayleigh") == 0) {
	  fitter_info.fittype[fitter_info.nrfunctions] = 6;
	}else {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Bug.");
	  return 0;
	}
	fitter_info.nrfunctions++;
        i++;
      }else if(strcmp(argv[i], "-tol") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &ftol, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-N") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &fitter_info.nr_model_points, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-ks") == 0) {
	fitter_info.method = 1;
      }else if(strcmp(argv[i], "-chi2hist") == 0) {
	fitter_info.method = 2;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &fitter_info.dx, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;	
      }else if(strcmp(argv[i], "-chi2hist_thresh") == 0) {
	fitter_info.threshold_values = argv[++i];
      }else if(strcmp(argv[i], "-chi2hist_opt") == 0) {
	fitter_info.pdist_options = argv[++i];
      }else if(strcmp(argv[i], "-plotcdf") == 0) {
	plotcdf = 1;
      }else if(strcmp(argv[i], "-2nd") == 0) {
	second_distr_specified = 1;
	int ret;
	ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 2, "%lf %lf %d", &second_distr_frac, &second_distr_dfrac, &function_param_fixed[2][0], NULL);
	if(ret == 2) {
	  function_param_fixed[2][0] = 0;
	}else if(ret != 3) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse the %s option. Needs to have two or three values.\n", argv[i]);
	  return 0;
	}

	// Check if 3 parameters can be read in
	/*
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 1, -1, "%lf %lf %d", &second_distr_frac, &second_distr_dfrac, &function_param_fixed[2][0], NULL) == 0) {
	  // Read in 3 parameters
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %d", &second_distr_frac, &second_distr_dfrac, &function_param_fixed[2][0], NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	// Check if 2 parameters can be read in
	}else if(parse_command_string(application.verbose_state, argc, argv, i+1, 1, -1, "%lf %lf", &second_distr_frac, &second_distr_dfrac, NULL) == 0) {
	  function_param_fixed[2][0] = 0;
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &second_distr_frac, &second_distr_dfrac, NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	}else {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse the %s option. Needs to have two or three values.\n", argv[i]);
	  return 0;
	}
	*/
	i++;
      }else if(strcmp(argv[i], "-null") == 0) {
	fitter_info.null_distr_specified = 1;
	int ret;
	ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 2, "%lf %lf %d", &null_distr_av, &null_distr_dav, &function_param_fixed[3][0], NULL);
	if(ret == 2) {
	  function_param_fixed[3][0] = 0;
	}else if(ret != 3) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse the %s option. Needs to have two or three values.\n", argv[i]);
	  return 0;
	}
	/*
	// Check if 3 parameters can be read in
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 1, -1, "%lf %lf %d", &null_distr_av, &null_distr_dav, &function_param_fixed[3][0], NULL) == 0) {
	  // Read in 3 parameters
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %d", &null_distr_av, &null_distr_dav, &function_param_fixed[3][0], NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	// Check if 2 parameters can be read in
	}else if(parse_command_string(application.verbose_state, argc, argv, i+1, 1, -1, "%lf %lf", &null_distr_av, &null_distr_dav, NULL) == 0) {
	  function_param_fixed[3][0] = 0;
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &null_distr_av, &null_distr_dav, NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	}else {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse the %s option. Needs to have two or three values.\n", argv[i]);
	  return 0;
	}
	*/
	i++;
      }else if(strcmp(argv[i], "-sigma") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &fitter_info.sigma, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pdistFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-noisefile") == 0) {
	fitter_info.noisefile_id = argv[i+1];
        i++;
      }else {
	/* If the option is not recognized, assume it is a filename */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "pdistFit: Unknown option: %s\n\nRun pdistFit without command line arguments to show help", argv[i]);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    /*
      if(strcmp(argv[i], "-f") == 0) {
	j = sscanf(argv[i+1], "%d", &fittype);
	if(j != 1) {
	  fprintf(stderr, "Cannot parse -f option, need one values.\n");
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-f2") == 0) {
	j = sscanf(argv[i+1], "%d %f %f", &fittype2, &fstart, &df);
	if(j != 3) {
	  fprintf(stderr, "Cannot parse -f option, need three values.\n");
	  return 0;
	}
	nfitparameters = 5;
	i++;
      }else if(strcmp(argv[i], "-p") == 0) {
	j = sscanf(argv[i+1], "%f %f", &pstart, &dp);
	if(j != 2) {
	  fprintf(stderr, "Cannot parse -p option, need two values.\n");
	  return 0;
	}
	pset = 1;
	i++;
      }else if(strcmp(argv[i], "-p2") == 0) {
	j = sscanf(argv[i+1], "%f %f", &pstart2, &dp2);
	if(j != 2) {
	  fprintf(stderr, "Cannot parse -p2 option, need two values.\n");
	  return 0;
	}
	pset2 = 1;
	i++;
      }else if(strcmp(argv[i], "-c") == 0) {
	j = sscanf(argv[i+1], "%f %f", &cstart, &dc);
	if(j != 2) {
	  fprintf(stderr, "Cannot parse -c option, need two values.\n");
	  return 0;
	}
	cset = 1;
	i++;
      }else if(strcmp(argv[i], "-c2") == 0) {
	j = sscanf(argv[i+1], "%f %f", &cstart2, &dc2);
	if(j != 2) {
	  fprintf(stderr, "Cannot parse -c2 option, need two values.\n");
	  return 0;
	}
	cset2 = 1;
	i++;
      }else if(strcmp(argv[i], "-m") == 0) {
	j = sscanf(argv[i+1], "%f %f", &mstart, &dm);
	if(j != 2) {
	  fprintf(stderr, "Cannot parse -m option, need two values.\n");
	  return 0;
	}
	mset = 1;
	i++;
      }else if(strcmp(argv[i], "-m2") == 0) {
	j = sscanf(argv[i+1], "%f %f", &mstart2, &dm2);
	if(j != 2) {
	  fprintf(stderr, "Cannot parse -m2 option, need two values.\n");
	  return 0;
	}
	mset2 = 1;
	i++;
      }else if(strcmp(argv[i], "-s") == 0) {
	j = sscanf(argv[i+1], "%f %f", &sstart, &ds);
	if(j != 2) {
	  fprintf(stderr, "Cannot parse -s option, need two values.\n");
	  return 0;
	}
	sset = 1;
	i++;
      }else if(strcmp(argv[i], "-s2") == 0) {
	j = sscanf(argv[i+1], "%f %f", &sstart2, &ds2);
	if(j != 2) {
	  fprintf(stderr, "Cannot parse -s2 option, need two values.\n");
	  return 0;
	}
	sset2 = 1;
	i++;
      }else if(strcmp(argv[i], "-dx") == 0) {
	j = sscanf(argv[i+1], "%f", &binning_dx);
	if(j != 1) {
	  fprintf(stderr, "Cannot parse -dx option, need one value.\n");
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-fix") == 0) {
	j = sscanf(argv[i+1], "%d", &fixparameters[nrfixparameters++]);
	if(j != 1) {
	  fprintf(stderr, "Cannot parse -fix option, need one value.\n");
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-emax_mod") == 0) {
	j = sscanf(argv[i+1], "%f", &max_e);
	if(j != 1) {
	  fprintf(stderr, "Cannot parse -emax_mod option, need one value.\n");
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-emin_mod") == 0) {
	j = sscanf(argv[i+1], "%f", &min_e);
	if(j != 1) {
	  fprintf(stderr, "Cannot parse -emin_mod option, need one value.\n");
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-ne_mod") == 0) {
	j = sscanf(argv[i+1], "%ld", &nrebins);
	if(j != 1) {
	  fprintf(stderr, "Cannot parse -ne_mod option, need one value.\n");
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-t") == 0) {
	j = sscanf(argv[i+1], "%f", &threshold);
	if(j != 1) {
	  fprintf(stderr, "Cannot parse -t option, need one value.\n");
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-nmod") == 0) {
	j = sscanf(argv[i+1], "%ld", &nmodelpoints);
	if(j != 1) {
	  fprintf(stderr, "Cannot parse -nmod option, need one value.\n");
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-null") == 0) {
	j = sscanf(argv[i+1], "%f", &avE);
	if(j != 1) {
	  fprintf(stderr, "Cannot parse -null option, need one value.\n");
	  return 0;
	}
	sprintf(extraoptions, " -A %f", avE);
	i++;
      }else if(strcmp(argv[i], "-v") == 0) {
	verbose = 1;
      }else {
	fprintf(stderr, "Unknown option \"%s\".\n", argv[i]);
	return 0;
      }
    */
    }
  }


  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) != 1) {
    printerror(application.verbose_state.debug, "ERROR pdistFit: One input file was expected at the end of the command line with the measured distribution to be fitted.");
    return 0;
  }

  if(fitter_info.nrfunctions == 0) {
    printerror(application.verbose_state.debug, "ERROR pdistFit: No fit function specified.");
    return 0;
  }
  if((fitter_info.threshold_values != NULL || fitter_info.pdist_options != NULL) && fitter_info.method != 2) {
    printerror(application.verbose_state.debug, "ERROR pdistFit: The -chi2hist_thresh and -chi2hist_opt options should be used in combination with -chi2hist.");
    return 0;
  }
  if(fitter_info.null_distr_specified && fitter_info.sigma <= 0 && fitter_info.noisefile_id == NULL) {
    printwarning(application.verbose_state.debug, "WARNING pdistFit: The -null option does causes a discontinuity in the cumulative distribution, which can result in unrealistic fits. Inclusion of the -sigma or -noisefile option will get rid of the discontinuity.");
  }

  double xstart[8], dx[8], xfit[8], teststatistic;
  int fixed[8], nrparams, nritt;
  nrparams = 0;
  for(i = 0; i < fitter_info.nrfunctions; i++) {
    if(fitter_info.fittype[i] == 1 || fitter_info.fittype[i] == 2 || fitter_info.fittype[i] == 3 || fitter_info.fittype[i] == 4 || fitter_info.fittype[i] == 5) {
      xstart[nrparams] = function_param[i][0];
      dx[nrparams]     = function_param[i][1];
      fixed[nrparams]  = function_param_fixed[i][0];
      nrparams++;
      xstart[nrparams] = function_param[i][2];
      dx[nrparams]     = function_param[i][3];
      fixed[nrparams]  = function_param_fixed[i][1];
      nrparams++;
    }else if(fitter_info.fittype[i] == 6) {
      xstart[nrparams] = function_param[i][0];
      dx[nrparams]     = function_param[i][1];
      fixed[nrparams]  = function_param_fixed[i][0];
      nrparams++;
    }else {
      printerror(application.verbose_state.debug, "ERROR pdistFit: Bug.");
      return 0;
    }
    if(i == 0 && fitter_info.nrfunctions == 2) {
      xstart[nrparams] = second_distr_frac;
      dx[nrparams]     = second_distr_dfrac;
      fixed[nrparams]  = function_param_fixed[2][0];
      fitter_info.frac_paramnr = nrparams;
      nrparams++;
    }
  }
  if(fitter_info.null_distr_specified) {
    xstart[nrparams] = null_distr_av;
    dx[nrparams]     = null_distr_dav;
    fixed[nrparams]  = function_param_fixed[3][0];
    nrparams++;
  }
  if(nrparams < 2) {
    printerror(application.verbose_state.debug, "ERROR pdistFit: The minimisation algorithm requires more than 1 fit parameter. The current fit only has %d free parameters, hence this type of fit is not supported.", nrparams);
    terminateApplication(&application);
    return 0;
  }
  fitter_info.cmdline = malloc(10000);
  fitter_info.txt = malloc(10000);
  if(fitter_info.cmdline == NULL || fitter_info.txt == NULL) {
    printerror(application.verbose_state.debug, "ERROR pdistFit: Memory allocation error.");
    return 0;
  }

  if(application.verbose_state.debug)
    fitter_info.debug = 1;
  if(application.fixseed) {
    fitter_info.fixseed = 1;
  }else {
    fitter_info.fixseed = 0;
    randomize_idnum(&fitter_info.seednr);
  }
  if(application.verbose_state.nocounters)
    fitter_info.nocounters = 1;

  fitter_info.measurement_file = getNextFilenameFromList(&application, argv, application.verbose_state);
  if(fitter_info.measurement_file == NULL) {
    printerror(fitter_info.debug, "ERROR pdistFit: Bug.");
    return 0;
  }
  if(polspecified || fitter_info.colspecified == 0) { // If -pol is specified, or if -col is not specified, it is not clear if the input is an ascii list of values, or a recognized pulsar format, like a penergy file.
    // If iformat has not been specified on command line, try to guess if the file is in a known format
    if(application.iformat <= 0) { // If -col is used, ascii file is implied 
      application.iformat = guessPSRData_format(fitter_info.measurement_file, 1, application.verbose_state);
      if(application.iformat == -2 || application.iformat == -3)
	return 0;
    }
    
    if(polspecified && application.iformat <= 0) {
      printerror(application.verbose_state.debug, "ERROR pdistFit: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized. Data is not recognized as pulsar data, but -pol was used. Use -col for ascii files.\n", fitter_info.measurement_file);
      return 0;	
    }
  }
  if(fitter_info.colspecified && application.iformat > 0) {
    printerror(application.verbose_state.debug, "ERROR pdistFit: Data is recognized as pulsar data, but -col was used. Use -pol to specify which input polarization channel you want to use in input, or use -col to force file to be read as a simple ascii file.\n");
    return 0;	
  }
  // If input format is recognized, read in the data and dump requested column in a new ascii file
  if(application.iformat > 0) {
    verbose_definition verbose;
    cleanVerboseState(&verbose);
    copyVerboseState(application.verbose_state, &verbose);
    if(application.iformat == PSRCHIVE_ASCII_format) {
      if(verbose.debug == 0) // Do not show header etc. for penergy files, unless in debug mode
	verbose.verbose = 0;
    }
    datafile_definition datain;
    //    cleanPSRData(&datain, application.verbose_state); 
    i = openPSRData(&datain, fitter_info.measurement_file, application.iformat, 0, 1, 1, verbose);
    if(i == 0) {
      printerror(application.verbose_state.debug, "ERROR pdistFit: Error opening data");
      return 0;
    }
    if(polspecified == 0) {
      fitter_info.file_column1 = 0;  // Default is first polarization channel
    }
    if(application.iformat == PSRCHIVE_ASCII_format && datain.gentype == GENTYPE_PENERGY) {
      if(application.verbose_state.verbose)
	printf("The file is a penergy file generated in mode 1.\n");
      if(polspecified == 0) {
	printwarning(application.verbose_state.debug, "WARNING pdist: -pol not specified, but input is penergy output. Using the the integrated on-pulse energy (-pol 2) by default.");
	fitter_info.file_column1 = 2;
      }
    }
    if(application.verbose_state.debug) {
      printf("Loading %ld points from input file and writing it to measurement_tmp.dist\n", datain.NrSubints*datain.NrFreqChan*datain.NrBins);
    }
    if(access("measurement_tmp.dist", F_OK) == 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: File measurement_tmp.dist already exist. Remove or move this file first, as it will be overwritten otherwise.");
      return 0;
    }
    FILE *fout;
    fout = fopen("measurement_tmp.dist", "w");
    if(fout == NULL) {
      printerror(fitter_info.debug, "ERROR pdistFit: Error opening temporary file measurement_tmp.dist.");
      return 0;
    }
    long subintnr, freqnr, binnr;
    float dummy_float;
    for(subintnr = 0; subintnr < datain.NrSubints; subintnr++) {
      for(freqnr = 0; freqnr < datain.NrFreqChan; freqnr++) {
	for(binnr = 0; binnr < datain.NrBins; binnr++) {
	  if(readPulsePSRData(&datain, subintnr, fitter_info.file_column1, freqnr, binnr, 1, &dummy_float, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pdist: Read error, shouldn't happen.\n");
	    return 0;
	  }
	  fprintf(fout, "%e\n", dummy_float);
	}
      }
    }
    fclose(fout);
    // Now specify new input file
    fitter_info.colspecified = 1;
    fitter_info.file_column1 = 1;
    polspecified = 0;
    fitter_info.measurement_file = "measurement_tmp.dist";
    remove_tmp_dist = 1;
    closePSRData(&datain, 0, application.verbose_state);
  } // End of recognized file format
  if(fitter_info.colspecified == 0) {
    fitter_info.colspecified = 1;
    fitter_info.file_column1 = 1;
    printerror(fitter_info.debug, "WARNING pdistFit: Going to use the default: -col 1.");
  }
  // Rebin input if required
  if(make_pdist_cmd(fitter_info.measurement_file, "measurement_tmp.hist", 0)) {
    if(access("measurement_tmp.hist", F_OK) == 0) {
      printerror(fitter_info.debug, "ERROR pdistFit: File measurement_tmp.hist already exist. Remove or move this file first, as it will be overwritten otherwise.");
      return 0;
    }
    if(application.verbose_state.debug) {
      printf("Rebinning input data using:\n  %s\n", fitter_info.cmdline);
    }else {
      printf("Rebinning input data:\n");
    }
    fflush(stdout);
    system(fitter_info.cmdline);
    printf("\n");
  }

  if(application.verbose_state.verbose) {
    printf("Fitting process started\n");
    fflush(stdout);
  }

  int ret;
  ret = doAmoeba_d(0, xstart, dx, fixed, xfit, &teststatistic, nrparams, &funk, ftol, &nritt, application.verbose_state.verbose, 0, 0.0, NULL, NULL);
  if(ret == 1) {
    printerror(application.verbose_state.debug, "Error pdistFit: Downhill-Simplex method did not converge. You can try lowering the tolerance with -ftol.");
    return 0;
  }else if(ret != 0) {
    printerror(application.verbose_state.debug, "Error pdistFit: Downhill-Simplex method failed.");
    return 0;
  }

  printf("\n\n");
  if(application.verbose_state.verbose) {
    printf("After %d steps the down-hill simplex found the best solution with a test statistic = %e\n\n", nritt, teststatistic);
    make_fakeDist_cmd(xfit, NULL);
    printf("The fitted distribution can be generated with:\n  %s\n", fitter_info.cmdline);
    if(make_pdist_cmd(fitter_info.measurement_file, NULL, 0)) {
      printf("The input data can be binned with:\n  %s\n", fitter_info.cmdline);
      printf("A similar command can be used for the produced fitted distribution.\n");
    }
    //    if(application.verbose_state.debug == 0) {
    //      fitter_info.debug = 1;
    //      funk(xfit);
    //      fitter_info.debug = 0;
    //    }
    printf("\n");
  }

  j = 0;
  for(i = 0; i < fitter_info.nrfunctions; i++) {
    if(fitter_info.nrfunctions == 2) {
      printf("Distribution %ld:\n", i+1);
    }
    if(fitter_info.fittype[i] == 1) {
      printf("k     = %e\n", xfit[j++]);
      printf("theta = %e\n", xfit[j++]);
    }else if(fitter_info.fittype[i] == 2) {
      printf("min   = %e\n", xfit[j++]);
      printf("max   = %e\n", xfit[j++]);
    }else if(fitter_info.fittype[i] == 3 || fitter_info.fittype[i] == 4) {
      printf("mu    = %e\n", xfit[j++]);
      printf("sigma = %e\n", xfit[j++]);
    }else if(fitter_info.fittype[i] == 5) {
      printf("idx   = %e\n", xfit[j++]);
      printf("min   = %e\n", xfit[j++]);
    }else if(fitter_info.fittype[i] == 6) {
      printf("sigma = %e\n", xfit[j++]);
    }else {
      printerror(application.verbose_state.debug, "ERROR pdistFit: Bug.");
      return 0;
    }
    if(i == 0 && fitter_info.nrfunctions == 2) {
      j++;
    }
  }
  if(fitter_info.nrfunctions == 2) {
    printf("The second distribution has a fraction of occurance of %e\n", xfit[fitter_info.frac_paramnr]);
  }

  if(fitter_info.method == 2) {
    fflush(stdout);
    system("rm measurement_tmp.hist");
  }

  if(application.verbose_state.verbose || plotcdf) {
    printf("\n");
    make_fakeDist_cmd(xfit, "model_trial_tmp.dist");
    fflush(stdout);
    system(fitter_info.cmdline);

    printf("Running KS-test:\n  ");
    sprintf(fitter_info.cmdline, "pstat -ks -col1 1 ");
    if(fitter_info.colspecified) {
      sprintf(fitter_info.txt, "-col2 %d ", fitter_info.file_column1);
      strcat(fitter_info.cmdline, fitter_info.txt);
    }
    sprintf(fitter_info.txt, "model_trial_tmp.dist %s", fitter_info.measurement_file);
    strcat(fitter_info.cmdline, fitter_info.txt);
    fflush(stdout);
    system(fitter_info.cmdline);

    if(plotcdf) {
      if(fitter_info.colspecified && fitter_info.file_column1 != 1) {
	printwarning(application.verbose_state.debug, "WARNING pdistFit: The -plotcdf function only works properly when the -col option is set to 1. Showing model distribution and data seperately.");
	sprintf(fitter_info.cmdline, "pdist -cdf -plot model_trial_tmp.dist");
	system(fitter_info.cmdline);
	sprintf(fitter_info.cmdline, "pdist -cdf -plot -col %d %s", fitter_info.file_column1, fitter_info.measurement_file);
	system(fitter_info.cmdline);
      }else {
	sprintf(fitter_info.cmdline, "pdist -cdf -overplot model_trial_tmp.dist %s", fitter_info.measurement_file);
	system(fitter_info.cmdline);
      }
    }
    fflush(stdout);
    system("rm model_trial_tmp.dist");
  }

  if(remove_tmp_dist)
    system("rm measurement_tmp.dist");

  free(fitter_info.cmdline);
  free(fitter_info.txt);

  terminateApplication(&application);
  return 0;
}