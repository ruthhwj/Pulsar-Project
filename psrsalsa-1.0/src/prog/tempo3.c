/* 
First do:

gcc -Wall -I$PGPLOT_DIR -c ppgplot.c
gcc -Wall -I$PGPLOT_DIR -c pgplot.c
gcc -Wall -I$PGPLOT_DIR -c myio.c
gcc -Wall -I$PGPLOT_DIR -c amoeba_ld.c
gcc -Wall -I$PGPLOT_DIR -c nr_ld.c
gcc -Wall -I$PGPLOT_DIR -c nrutil.c

To compile on wulfgeat as neil:
source /psr/PSR.cshrc
source /psr/PSR.login

gcc ppgplot.o pgplot.o myio.o nr_ld.o amoeba_ld.o nrutil.o -Wall -L$PGPLOT_DIR -lcpgplot -lpgplot -lm -L/usr/X11R6/lib -lX11 -L/usr/lib/gcc/i386-redhat-linux5E/4.1.2/ -lg2c -o tempo3  `pkg-config --cflags --libs gtk+-3.0` tempo3.c


To compile on windmill as patrick:

gcc ppgplot.o pgplot.o myio.o nr_ld.o amoeba_ld.o nrutil.o -Wall -L$PGPLOT_DIR -lcpgplot -lpgplot -lg2c -lm -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/ -lX11 -lpng -o tempo3 tempo3.c



Really needs to be fixed:

- Want to be able to zoom in frequency and then fit that part of the data?

- Number NRTOA is not written out correctly when zoomed in.
  Look for HaveALookHereIfFixingChi2BugZoomed and for nrtoasincluded

- If change nr waves, wave_om changes and solution explodes. Probably
  need to fix wave_om. If not fit wave, wave_om shouldn't change.

More challenging changes to be made:

- Probably most parameters have to be made a fit parameter, even if it
  is actually never fitted for. This should be done for all parameters
  which are read in.

- put covariancs matrix stuff in amoeba, possibly quit difficult
  because this find_dx must be implemented as well. At the other hand,
  that will be an useful option to have in general as well.

- In the gtk window a button could be added to ignore a
  parameter. When tracking pulse numbers is on this allows the user to
  see what effect a parameter has.

Checks that needs to be done: 

- measurenu/nudot still working properly, also after deleting points?


*/

// Comment the following line out if you do not want to make use of certain development psrsalsa libraries (less functionality, but less dependent code). It is using NR code: matrix_ld(), nrerror_ld(), gaussj_ld(), ivector(), free_ivector(), free_matrix_ld()
#define EnablePSRSALSAdevelop   1

#define MaxNrStack              2

#define PARSEC                    3.08567758149137e16 // m

#define MaxNrGTKLines 200
int nrActiveGTKLines;
extern int GTKinitialised;

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h> 
#include <time.h>
#include <pthread.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include <gsl/gsl_sort.h>
#include "psrsalsa_tempo3.h"
#ifdef TEMPO3_EnableGTK
  #include <gtk/gtk.h>
  #include <signal.h>
#else
  #include <pthread.h>
  #warning Make sure to enable -pthread in CFLAGS in the Makefile
#endif
#include "psrsalsa.h"

#ifdef EnablePSRSALSAdevelop
  void SHOWREVISIONINFO_prog() {
  #include "tempo3.c.svninfo"
  }
#else
  #include "tempo3_nodevel.h"
#endif


//#define HARDCODED_NR_PARAMETERS 2149  // The expected nr of parameters. Very ugly, but this number needs to be known to save time re-writing the code. The program will check and report if this number is not set correctly.

extern parfile_description_def tempo3_parfile_descr;

int repeatable;                             // If set, no lines containing random variables are produced

// Non-zero if the user requested to quit via a button or something
int user_requested_quit = 0; 

typedef struct {
  long double *resid;                       /* The corresponding residual (in turns) */
  float *residy;                            /* Space for resids stored stored as floats. Generate when used (e.g. plotting) */
  long long *nturn;                         /* The number of rotations of the star. */
  int *nturn_fixed;                          /* This flag is 1 if the turn number is read in from tim file*/
  long double mjdmin, mjdmax;               /* The data span which is shown (and for which the residuals are calculated) */
  long double freqmin, freqmax;             /* The data span which is shown (and for which the residuals are calculated) */
  long double variance;                     /* The variance of the solution */
  long double chi2;                         /* The (weighed) chi**2 value of the solution (not reduced). */
}residual_def;

void free_residuals(residual_def *resid)
{
  free(resid->resid);
  free(resid->residy);
  free(resid->nturn);
  free(resid->nturn_fixed);
}


// Copy everything in tempo3_parameters_def, but not the errors for example
void copy_parfile(tempo3_parameters_def *copy, tempo3_parameters_def *orig)
{
  int i;
  long double *orig_dplus_pointer = copy->dplus;
  long double *orig_dmin_pointer = copy->dmin;
  long double *orig_param_pointer = copy->parameter;
  int *orig_companion_defined_pointer = copy->companion_defined;
  int *orig_jumpmeaning_pointer = copy->jumpmeaning;
  char **orig_jumpflag_pointer = copy->jumpflag;

  memcpy(copy, orig, sizeof(tempo3_parameters_def));  // Make a copy of everything

  copy->dplus = orig_dplus_pointer;  // Except the location (and contents) of the errorbars
  copy->dmin = orig_dmin_pointer;

  copy->parameter = orig_param_pointer;   // except the location of the array of parameters
  memcpy(copy->parameter, orig->parameter, TEMPO3_PARAM_TOTAL*sizeof(long double)); // Now make a copy of the parameters

  copy->companion_defined = orig_companion_defined_pointer;   // and the location of the array of companion_defined
  memcpy(copy->companion_defined, orig->companion_defined, TEMPO3_MaxNrCompanions*sizeof(int)); // Now make a copy of the array

  copy->jumpflag = orig_jumpflag_pointer;   // and the location of the array of jumpflag
  for(i = 0; i < orig->nrjumps; i++) {
    memcpy(copy->jumpflag[i], orig->jumpflag[i], TEMPO3_MaxFlagsLength); // Now make a copy of the array
  }

  copy->jumpmeaning = orig_jumpmeaning_pointer;   // and the location of the array of jumpmeaning
  for(i = 0; i < orig->nrjumps; i++) {
    copy->jumpmeaning[i] = orig->jumpmeaning[i]; // Now make a copy of the array
  }
}


/* Global variables only to be used by the fitting function funk. */
struct {
  toa_data_def *toas;                 /* Obviously TOAs are needed by the fitting function */
  residual_def *residuals;            /* Could be defined in funk, but then need to do a lot of mallocs/frees every time */
  tempo3_parameters_def *parfile;               /* Makes sure that non-fit parameters are known by the fit function */
  tempo3_parameters_def parfile_funk;           /* Temporary parfile stucture which is filled with the fitted parameters */
  tempo3_wraps_def *wraps;
  int superverbose;
  long double *best_xfit, best_chi2;   /* Remember best solution in case there was no convergence */
  int *paramlinked;
  int *paramset;
  int numthreads;
  int nosort;
}additional_funk_info;

typedef struct {
  tempo3_parameters_def parfiles[MaxNrStack];
  tempo3_wraps_def wraps[MaxNrStack];
  float previous_xstart[MaxNrStack];
  float previous_xend[MaxNrStack];
  int *fixed[MaxNrStack];
  int *paramset[MaxNrStack];
  int nrinstack;
  /* Only store some of the toa information, as rest does not change or is not important */
  int nrtoas[MaxNrStack];        
  int nrtoasdeleted[MaxNrStack]; 
  int *deleted[MaxNrStack]; 
}tempo3stack;


#define PLOTTYPE_NBREAK        30
#define PLOTTYPE_ORBITALPHASE  31
#define PLOTTYPE_FREQUENCY    100

void calcChi2(residual_def *residuals, toa_data_def toas, long double *variance, long double spinfreq, long double *chi2w);
void draw_residual(toa_data_def toas, int toanr, float x, float y, char *highlightedFlag, float markersize, int plottype, int mode, long double spinfreq);
void calcphases(tempo3_parameters_def *parfile, toa_data_def toas, int showWaves, residual_def *residuals, int *paramlinked, int *paramset, int ignorebinary, int numthreads);
/*void calcresids(tempo3_parameters_def parfile, toa_data_def toas, residual_def *residuals, tempo3_wraps_def wraps);*/
void subtractTurns(toa_data_def toas, residual_def *residuals, tempo3_wraps_def *wraps, int nosort, int fixnturn);
void perturb_parfile(tempo3_parameters_def *parfile, int *paramset, int *paramlinked, long double fraction);
long double pickvalue(int param, tempo3_parameters_def *parfile);
void setvaluevalue(int param, long double value, tempo3_parameters_def *parfile);
void print_ephemeris_ExtraParams(FILE *stream, tempo3_parameters_def *parfile, int *paramset, int *paramlinked, int showerrors, long double *dplus, long double distance, long double distanceerror, int hierarchical);
void fill_xstart(tempo3_parameters_def *parfile, long double *xstart);
long double get_initial_step_size(int param, tempo3_parameters_def *parfile, toa_data_def toas, residual_def residuals, tempo3_wraps_def *wraps, int superverbose, int *paramlinked, int *paramset, int nosort, int numthreads);
void fill_dx(tempo3_parameters_def *parfile, toa_data_def toas, residual_def residuals, tempo3_wraps_def *wraps, int *fixed, long double *dx, int superverbose, int *paramlinked, int *paramset, int nosort, int numthreads);
void use_xfit(tempo3_parameters_def *parfile, long double *xfit);
void calcAlternativePlots(float *ygraph, float *ygraph2, char *ylabel, int plottype, int showdata, toa_data_def *toas, residual_def residuals, tempo3_parameters_def *parfile, int *paramlinked, int *paramset, int nummericalApproxNuNudot);
int doplot(char *device, float *ygraph, float *ygraph2, int nophaselines, int plottype, int showdata, int showWaves, toa_data_def toas, residual_def *residuals, ephemeris_def *eph, int superverbose, int writeoutpieces, int plottitle, float *cur_x1, float *cur_y1, float *cur_x2, float *cur_y2, float previous_xstart, float previous_xend, toa_data_def toas_fiducial, residual_def residuals_fiducial, int showfiducial, int showpepoch, char *highlightedFlag, int nummericalApproxNuNudot, char *title, int amoeba_algorithm, float markersize, float stridefit_nrdays, float stridefit_mintoarange, float stridefit_minstride, int stridefit_minnroftoas, int stridefit_includeF2, int stridefit_includeF3, char *stridefit_scriptfilename, int nosort, int numthreads, verbose_definition verbose);
void dofit(int param, tempo3_parameters_def *parfile, toa_data_def toas, residual_def *residuals, tempo3_wraps_def *wraps, int *paramlinked, int *paramset, int nosort, int numthreads);
int measureNuNudot(long double *nu, long double *nudot, long double *nu_mjd, float stridefit_nrdays, float stridefit_mintoarange, float stridefit_minstride, int minnrtoas, float *ygraph, float *ygraph2, toa_data_def *toas, residual_def residuals, ephemeris_def *eph, int superverbose, int writeoutpieces, int amoeba_algorithm, int stridefit_includeF2, int stridefit_includeF3, char *stridefit_scriptfilename, int nosort, int numthreads, verbose_definition verbose);
void storeInStack(tempo3stack *stack, tempo3_parameters_def *parfile, toa_data_def timfile, tempo3_wraps_def *wraps, float previous_xstart, float previous_xend, int *fixed, int *paramset);
void popstack(tempo3stack *stack, tempo3_parameters_def *parfile, toa_data_def *timfile, tempo3_wraps_def *wraps, float *previous_xstart, float *previous_xend, int *fixed, int *paramset);
float timediff(struct timeval start, struct timeval stop);
void derive_orbital_parameters(tempo3_parameters_def *parfile, int *paramset, int companion, long double mjd, long double inclination, long double innermass, int doprint, long double *xpsini, long double *ypsini, long double *xc, long double *yc, long double *companionmass, int showerrors, long double *dplus);
void plotOrbit(tempo3_parameters_def *parfile, int *paramset, toa_data_def *toas, long double t, long double tmin, long double tmax, long double inclination, long double mpsr, int hierarchical, int markersize, char *device, char *orbit_dump_prefix, int finderrors, long double *dplus);
void finderrors_tempo3(long double *dplus, long double *dmin, tempo3_parameters_def *parfile, toa_data_def toas, residual_def residuals, int *fixed, int *paramset, int *paramlinked, long double *xstart, long double *dx, long double *xfit, double *mrq_alpha, double *mrq_alpha_inv, long max_nr_fitted_params, long double *mrq_deriv, int dumpcovar, verbose_definition verbose_state);
long double calc_derivative(int parameternr, long double *parfile_params, tempo3_parameters_def *parfile_tmp, long double dx, long double mjd, long double freq, long double ssb1, long double ssb2, long double ssb3, char *site, char *flags, int *paramlinked, int *paramset);

int allocate_ephemeris_par_only(tempo3_parameters_def *parfile, int allowerrors, verbose_definition verbose); // Is in library, but should be avoided to be used by program which should use initialise_ephemeris() instead
void free_ephemeris_par_only(tempo3_parameters_def *parfile); // Is in library, but should be avoided to be used by program which should use free_ephemeris() instead
void initialise_ephemeris_par_only(tempo3_parameters_def *parfile);  // Is in library, but should be avoided to be used by program which should use initialise_ephemeris() instead
void ephemeris_setLinkedValues_par_version(tempo3_parameters_def *parfile, int *paramlinked);  // Is in library, but should be avoided to be used by program which should not directly refer to tempo3_parameters_def
int ephemeris_check_if_parameter_is_linked(int param_id, int *paramlinked); // Is in library, but should be avoided to fiddle with linked values directly from program
long double ephemeris_estimateNu_numerically_avoid_using(long double mjd, long double freq, long double ssb1, long double ssb2, long double ssb3, char *site, char *flags, tempo3_parameters_def *parfile, int *paramset); // Is in library, but should be avoided. It probably should be completely abandoned if possible. Note it is also used in library.
long double ephemeris_estimateNudot_numerically_avoid_using(long double mjd, long double freq, long double ssb1, long double ssb2, long double ssb3, char *site, char *flags, tempo3_parameters_def parfile, int *paramset); // Is in library, but should be avoided. It probably should be completely abandoned if possible. Note it is also used in library.
long double fitwave_function(tempo3_parameters_def *parfile, long double mjd, int mode, int showWaves); // Is in library, but should be avoided to be used by program which should use evaluate ephemeris functions. Could easily make a ephemeris version as a wrapper function.
long double ephemeris_get_torb(long double mjd, tempo3_parameters_def *parfile, int *paramset, int calc_deriv, int whichcompanion); // Is in library, but should be avoided to be used by program which should use evaluate ephemeris functions. Could easily make a ephemeris version as a wrapper function.
long double evaluate_ephemeris(long double mjd, long double freq, long double ssb1, long double ssb2, long double ssb3, int showWaves, int noNu, char *site, char *flags, tempo3_parameters_def *parfile, long double dphase, int *paramset, int mode, int ignore_binary); // Is in library, but should be avoided to be used by program which should use evaluate ephemeris functions. 
int filter_parfile_for_tempo2(char *inputfile, char *outputfile, tempo3_parameters_def *parfile, verbose_definition verbose); // Is in library, but should be avoided to be used by program which should use loadtimfile() instead
int filter_tempo2_general2_output(char *filename_in, char *filename_out, int nosort, unsigned long **indx, int repeatable, verbose_definition verbose); // Is in library, but should be avoided to be used by program which should use loadtimfile() instead
int readtimfile_custom(char *filename, toa_data_def *toas, int shortformat, verbose_definition verbose);
 // Is in library, but should be avoided to be used by program which should use loadtimfile() instead


void toggleFixed(int ident, int *fixed, int *paramlinked)
{
  if(ident != TEMPO3_PARAM_PHASE0 && ident != TEMPO3_PARAM_PSRJNAME && ident != TEMPO3_PARAM_WAVEOM && ident != TEMPO3_PARAM_PEPOCH && ident != TEMPO3_PARAM_WAVEEPOCH && ident != TEMPO3_PARAM_TRES && ident != TEMPO3_PARAM_BINARY && (ident < TEMPO3_PARAM_GLEPMIN || ident >= TEMPO3_PARAM_GLEPMIN+TEMPO3_MaxNrGlitches) && (ident < TEMPO3_PARAM_GLEPMAX || ident >= TEMPO3_PARAM_GLEPMAX + TEMPO3_MaxNrGlitches) && ident != TEMPO3_PARAM_MODE && ident != TEMPO3_PARAM_TRACK && ident != TEMPO3_PARAM_TZRMJD && ident != TEMPO3_PARAM_TZRFRQ && ident != TEMPO3_PARAM_TZRSITE && ident != TEMPO3_PARAM_POSEPOCH && ident != TEMPO3_PARAM_DMEPOCH) {
    if(fixed[ident])
      fixed[ident] = 0;
    else
      fixed[ident] = 1;
    if(ephemeris_check_if_parameter_is_linked(ident, paramlinked) != 0)
      fixed[ident] = 1;      
    // Fit for WAVECOS as well if WAVESIN was pressed
    if(ident >= TEMPO3_PARAM_WAVESIN && ident < TEMPO3_PARAM_WAVESIN + TEMPO3_MaxNrWaves) {
      ident = ident - TEMPO3_PARAM_WAVESIN + TEMPO3_PARAM_WAVECOS;
      if(fixed[ident])
	fixed[ident] = 0;
      else
	fixed[ident] = 1;
    }
  }
}



/* int count = 0; */

long double funk(long double params[])
{
  long double chi2, chi2w;
  int i, ignore_binary;
  verbose_definition noverbose;
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;
  /*  if(count++ > 10)
    exit(0);
  */


  // Much is similar to the downhill simplex version

  copy_parfile(&(additional_funk_info.parfile_funk), additional_funk_info.parfile);
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    additional_funk_info.parfile_funk.parameter[i] = params[i];
  } 

  // Check if binary parameters will cause a problem, if so the chi2 is set to a very high number
  ignore_binary = 0;
  for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
    if(additional_funk_info.parfile_funk.companion_defined[i] != 0) {
      if(additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_ECC+i] < 0.0) {
	ignore_binary = 1;
      }else if(additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_ECC+i] > 1.0 && additional_funk_info.parfile_funk.binarymodel != 3) {
	ignore_binary = 1;
      }
    }
  }

  /*  fprintf(stderr, "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxx start %p\n", additional_funk_info.paramlinked); */
  calcphases(&additional_funk_info.parfile_funk, *(additional_funk_info.toas), 0, additional_funk_info.residuals, additional_funk_info.paramlinked, additional_funk_info.paramset, ignore_binary, additional_funk_info.numthreads);
  /*  fprintf(stderr, "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxx end\n"); */
  subtractTurns(*(additional_funk_info.toas), additional_funk_info.residuals, additional_funk_info.wraps, additional_funk_info.nosort, 1);

  calcChi2(additional_funk_info.residuals, *(additional_funk_info.toas), &chi2, additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_F0], &chi2w);
  /* If MODE is set to 1, take weighted chi2 rather than the rms. */
  if(additional_funk_info.parfile_funk.mode != 0) {
    chi2 = chi2w;
  }

  for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
    if(additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEPMIN+i] > 0 && additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEP+i] > 0 && (additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEP+i] < additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEPMIN+i])) {
      if(chi2 > 0)
	chi2 *= 100;
    }
    if(additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEPMAX+i] > 0 && additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEP+i] > 0 && (additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEP+i] > additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEPMAX+i])) {
      if(chi2 > 0)
	chi2 *= 100;
    }
  }

  if(ignore_binary)
    chi2 *= 100;
 

  //  printf("funk: %Le - rms=%Lf\n", additional_funk_info.parfile_funk.nbreak, sqrtl(chi2/(long double)(((*(additional_funk_info.toas)).nrtoas))));

  if(additional_funk_info.superverbose != -1) {

    /*    printf("funk: %Le (%Le) %Le chi2=%Lf\n", additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_F0], additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_F0]-2.45120394490856228316e+00, additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_PHASE0], chi2); */
    printf("funk: rms=%Lf                \r", sqrtl(chi2/(long double)(((*(additional_funk_info.toas)).nrtoas))));

    pgplot_options_definition pgplot_options;
    pgplot_clear_options(&pgplot_options);
    pgplot_options.viewport.noclear = 1;
    pgplot_options.viewport.dontopen = 1;
    pgplot_options.viewport.dontclose = 1;
    strcpy(pgplot_options.box.xlabel, "MJD");
    strcpy(pgplot_options.box.ylabel, "Resid [Phase]");
    strcpy(pgplot_options.box.title, "");

    for(i = 0; i < additional_funk_info.toas->nrtoas + additional_funk_info.toas->nrtoasdeleted; i++) {
      additional_funk_info.residuals->residy[i] = additional_funk_info.residuals->resid[i];
    }
    pgplotGraph1(&pgplot_options, (*(additional_funk_info.residuals)).residy, (*(additional_funk_info.toas)).mjdx, NULL, (*(additional_funk_info.toas)).nrtoas, -1, -1, 0, (*(additional_funk_info.residuals)).mjdmin, (*(additional_funk_info.residuals)).mjdmax, 0, 0, 0, 0, 1, 1, 17, 1, 1, NULL, -1, noverbose); 
    if(additional_funk_info.superverbose) {
      printf("Press a key to continue\n");
      fflush(stdout);
      pgetch();  
    }
  }

  /*
  if(parfile_funk.perturb > 0)
    chi2 += parfile_funk.perturb;
  else
    chi2 -= parfile_funk.perturb;
  */

  if(chi2 <= additional_funk_info.best_chi2 || additional_funk_info.best_chi2 < 0) {
    memcpy(additional_funk_info.best_xfit, params, TEMPO3_PARAM_TOTAL*sizeof(long double));
    additional_funk_info.best_chi2 = chi2;
    //    fprintf(stderr, "\nXXXXXX new best chi2=%Le %Le %Le\n", chi2, additional_funk_info.best_xfit[1], params[1]);
    //    fprintf(stderr, "XXXXXX new best chi2=%Le %Le %Le %Le\n", chi2, params[7], params[8], params[17]);
    //  }else {
    //    fprintf(stderr, "XXXXXX ignoring chi2=%Le %Le %Le %Le\n", chi2, params[7], params[8], params[17]);
  }

  //  printf("XXXXX chi2=%Le\n", chi2);
  return chi2;
}

typedef struct {
  int *fixed;
  int nrfitparameters;
}funk_levmar_info;

/* Evaluates the vector beta and the matrix alpha (NR eq. 15.5.8 on pg
   682 (really pg 706) for the current vector of solution parameters
   params.

   The vector beta can be written as:
   beta_k = sum over datapoints i{ (df/d param_k)*(y_i -f)/variance_i }

   The matrix alpha corresponds to the so-called Hessian
   matrix can be calculates as:

   alpha_jk = sum over datapoints i{ (df/d param_j)*(df/d param_k)/variance_i }


*/
void funk_levmar(long double *params, long double *alpha, long double *beta, long double *chi2, void *info)
{
  int i, param, fitparam, ignore_binary, *fixed, nrfitparameters;
  long nrdatapoints;
  long double *deriv;

  fixed = ((funk_levmar_info *)info)->fixed;
  nrfitparameters = ((funk_levmar_info *)info)->nrfitparameters;
  // Much is similar to the downhill simplex version

  deriv = malloc(nrfitparameters*sizeof(long double));
  if(deriv == NULL) {
    printerror(0, "ERROR tempo3: Memory allocation error.");
    exit(0);
  }

  copy_parfile(&(additional_funk_info.parfile_funk), additional_funk_info.parfile);
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    additional_funk_info.parfile_funk.parameter[i] = params[i];
  } 

  // Check if binary parameters will cause a problem, if so the chi2 is set to a very high number
  ignore_binary = 0;
  for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
    if(additional_funk_info.parfile_funk.companion_defined[i] != 0) {
      if(additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_ECC+i] < 0.0) {
	ignore_binary = 1;
      }else if(additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_ECC+i] > 1.0 && additional_funk_info.parfile_funk.binarymodel != 3) {
	ignore_binary = 1;
      }
    }
  }

  // Calculate the predicted phases 
  calcphases(&additional_funk_info.parfile_funk, *(additional_funk_info.toas), 0, additional_funk_info.residuals, additional_funk_info.paramlinked, additional_funk_info.paramset, ignore_binary, additional_funk_info.numthreads);
  subtractTurns(*(additional_funk_info.toas), additional_funk_info.residuals, additional_funk_info.wraps, additional_funk_info.nosort, 1);

  for(i = 0; i < nrfitparameters; i++) {
    beta[i] = 0;
  }
  for(i = 0; i < nrfitparameters*nrfitparameters; i++) {
    alpha[i] = 0;
  }
  // Calculate chi2
  *chi2 = 0;
  nrdatapoints = 0;
  for(i = 0; i < additional_funk_info.toas->nrtoas + additional_funk_info.toas->nrtoasdeleted; i++) {
    if(additional_funk_info.toas->deleted[i] == 0) {
      if(additional_funk_info.toas->mjd[i] >= additional_funk_info.residuals->mjdmin && additional_funk_info.toas->mjd[i] <= additional_funk_info.residuals->mjdmax) {
	nrdatapoints++;
	/* If MODE is set to 1, take weighted chi2 rather than the rms. */
	if(additional_funk_info.parfile_funk.mode == 0) {
	  *chi2 += additional_funk_info.residuals->resid[i]*additional_funk_info.residuals->resid[i]; // This is in turns
	}else {
	  *chi2 += additional_funk_info.residuals->resid[i]*additional_funk_info.residuals->resid[i]/(long double)(1e-12*additional_funk_info.toas->err[i]*additional_funk_info.toas->err[i]*additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_F0]*additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_F0]);  // Error is defined in micro-sec
	}
	fitparam = 0;
	// Now loop over all the parameters, and calculate the derivatives
	for(param = 0; param < TEMPO3_PARAM_TOTAL; param++) {
	  if(fixed[param] == 0) {
	    long double freq;
	    if(additional_funk_info.toas->freqSSB[i] > 0) {
	      freq = additional_funk_info.toas->freqSSB[i];
	    }else {
	      freq = additional_funk_info.toas->freqSite[i];
	    }
	    deriv[fitparam] = calc_derivative(param, params, &additional_funk_info.parfile_funk, -1.0, additional_funk_info.toas->mjd[i], freq, additional_funk_info.toas->ssb1[i], additional_funk_info.toas->ssb2[i], additional_funk_info.toas->ssb3[i], additional_funk_info.toas->site[i], additional_funk_info.toas->flags[i], additional_funk_info.paramlinked, additional_funk_info.paramset);  // In phase per unit of the parameter
	    deriv[fitparam] *= -1.0;
	    fitparam++;
	  }
	}

	// Now loop over all the parameters, identify which ones we fit for, and calculate alpha and beta
	fitparam = 0;
	for(param = 0; param < TEMPO3_PARAM_TOTAL; param++) {
	  if(fixed[param] == 0) {
	    if(additional_funk_info.parfile_funk.mode == 0) {
	      beta[fitparam] += deriv[fitparam]*additional_funk_info.residuals->resid[i];
	    }else {
	      beta[fitparam] += deriv[fitparam]*additional_funk_info.residuals->resid[i]/(long double)(1e-12*additional_funk_info.toas->err[i]*additional_funk_info.toas->err[i]*additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_F0]*additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_F0]);
	    }
	    int param2, fitparam2;
	    fitparam2 = 0;
	    for(param2 = 0; param2 < TEMPO3_PARAM_TOTAL; param2++) {
	      if(fixed[param2] == 0) {
		if(additional_funk_info.parfile_funk.mode == 0) {
		  alpha[fitparam*nrfitparameters+fitparam2] += deriv[fitparam]*deriv[fitparam2];
		}else {
		  alpha[fitparam*nrfitparameters+fitparam2] += deriv[fitparam]*deriv[fitparam2]/(long double)(1e-12*additional_funk_info.toas->err[i]*additional_funk_info.toas->err[i]*additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_F0]*additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_F0]);
		}
		fitparam2++;
	      }
	    }
	    fitparam++;
	  }
	}
      }
    }
  }


  // Check if glitches out of range
  for(i = 0; i < TEMPO3_MaxNrGlitches; i++) {
    if(additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEPMIN+i] > 0 && additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEP+i] > 0 && (additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEP+i] < additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEPMIN+i])) {
      if(*chi2 > 0)
	*chi2 *= 100;
    }
    if(additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEPMAX+i] > 0 && additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEP+i] > 0 && (additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEP+i] > additional_funk_info.parfile_funk.parameter[TEMPO3_PARAM_GLEPMAX+i])) {
      if(*chi2 > 0)
	*chi2 *= 100;
    }
  }

  if(ignore_binary)
    *chi2 *= 100;
 

  //  printf("funk: %Le - rms=%Lf\n", additional_funk_info.parfile_funk.nbreak, sqrtl(chi2/(long double)(((*(additional_funk_info.toas)).nrtoas))));


  if(*chi2 <= additional_funk_info.best_chi2 || additional_funk_info.best_chi2 < 0) {
    memcpy(additional_funk_info.best_xfit, params, TEMPO3_PARAM_TOTAL*sizeof(long double));
    additional_funk_info.best_chi2 = *chi2;
    /*    fprintf(stderr, "\nOLAOLAOLA %Le %Le %Le\n", chi2, additional_funk_info.best_xfit[1], params[1]); */
  }

  //  printf("XXXXX chi2=%Le %ld data points\n", *chi2, nrdatapoints);

  free(deriv);
  return;
}



// Returns the derivative dphase/dparameter for the specified parameter. This is either found analytically, or numerically.
// parfile_tmp should be the parfile used as the position at which to calculate the derivative, although the values are used from the array parfile_params are the parameter values (this array is not modified). This array is modified when a numerical derivative was used. mjd is the mjd at which the derivative is determined. dx is the change in the specified parameter to use to find the numerical derivative if required.
long double calc_derivative(int parameternr, long double *parfile_params, tempo3_parameters_def *parfile_tmp, long double dx, long double mjd, long double freq, long double ssb1, long double ssb2, long double ssb3, char *site, char *flags, int *paramlinked, int *paramset)
{
  long double dy, tmp_value, t_minus_tepoch, t_minus_tepoch2, t_minus_tepoch4, t_minus_tepoch8;
  long double deriv, deriv2;

  deriv = 0.0;

  if(parameternr == TEMPO3_PARAM_PHASE0) {
    return 1.0;
  }
  if(parfile_tmp->parameter[TEMPO3_PARAM_NBRAKE] == 0) {
    if(parameternr == TEMPO3_PARAM_F0) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      return t_minus_tepoch;
    }else if(parameternr == TEMPO3_PARAM_F0 + 1) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      return 0.5*t_minus_tepoch*t_minus_tepoch;
    }
  }else { // NBRAKE is specified
    /*
    long double expr1, expr2;
    long double f0_pwr_n = powl(parfile.parameter[TEMPO3_PARAM_F0], parfile.parameter[TEMPO3_PARAM_NBRAKE]);
    expr2 = (parfile.parameter[TEMPO3_PARAM_F0]-parfile.parameter[TEMPO3_PARAM_F0+1]*(parfile.parameter[TEMPO3_PARAM_NBRAKE]-1.0)*t_minus_tepoch)/f0_pwr_n;
    expr1 = f0_pwr_n/((2.0-parfile.parameter[TEMPO3_PARAM_NBRAKE])*parfile.parameter[TEMPO3_PARAM_F0+1]);
    phase += expr1*powl(expr2,(2.0-parfile.parameter[TEMPO3_PARAM_NBRAKE])/(1.0-parfile.parameter[TEMPO3_PARAM_NBRAKE]));
    phase -= parfile.parameter[TEMPO3_PARAM_F0]*parfile.parameter[TEMPO3_PARAM_F0]/(parfile.parameter[TEMPO3_PARAM_F0+1]*(2.0-parfile.parameter[TEMPO3_PARAM_NBRAKE]));

    expr2 = (f0-f1*(n-1.0)*(t-t0))/f0^n;
    expr1 = f0^n/((2-n)*f1);
    phase += expr1*(expr2^((2-n)/(1-n)));
    phase -= f0^2/(f1*(2-n));

    expr2 = (f0-f1*(n-1.0)*(t-t0))/f0^n;
    phase += +(f0^n/((2-n)*f1))*expr2^((2-n)/(1-n)) - f0^2/(f1*(2-n))

    DERIV f0
    expr2 = (f0-f1*(n-1.0)*(t-t0))/f0^n = f0^(1-n)-f1*(n-1.0)*(t-t0)/f0^n

    +(n*f0^(n-1)/((2-n)*f1))*expr2^((2-n)/(1-n)) +(f0^n/((2-n)*f1))*((2-n)/(1-n))expr2^((2-n)/(1-n) -1)diff_expr2 - 2*f0/(f1*(2-n))

    diff_expr2 = (1-n)f0^(-n)+n*f1*(n-1.0)*(t-t0)/f0^(1+n)

    +(n*f0^(n-1)/((2-n)*f1))*expr2^((2-n)/(1-n)) +(f0^n/((2-n)*f1))*((2-n)/(1-n))expr2^((2-n)/(1-n) -1)*((1-n)f0^(-n)+n*f1*(n-1.0)*(t-t0)/f0^(1+n)) - 2*f0/(f1*(2-n))
    +(n*f0^(n-1)/((2-n)*f1))*
expr2^((2-n)/(1-n)) +(f0^n/((2-n)*f1))*((2-n)/(1-n))expr2^((2-n)/(1-n) -1)*((1-n)f0^(-n)+n*f1*(n-1.0)*(t-t0)/f0^(1+n)) - 2*f0/(f1*(2-n))

    Need to write it out properly before implementing this
    */
  }

  if(parameternr >= TEMPO3_PARAM_F0 + 2 && parameternr <= TEMPO3_PARAM_F0 + TEMPO3_MaxNrFderivatives) {  // The parameter is a spin-frequency derivative


    if(parameternr == TEMPO3_PARAM_F0 + 2) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      return t_minus_tepoch*t_minus_tepoch*t_minus_tepoch/6.0;
    }else if(parameternr == TEMPO3_PARAM_F0 + 3) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
      return t_minus_tepoch2*t_minus_tepoch2/24.0;
    }else if(parameternr == TEMPO3_PARAM_F0 + 4) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
      return t_minus_tepoch*t_minus_tepoch2*t_minus_tepoch2/120.0;
    }else if(parameternr == TEMPO3_PARAM_F0 + 5) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
      t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
      return t_minus_tepoch2*t_minus_tepoch4/720.0;
    }else if(parameternr == TEMPO3_PARAM_F0 + 6) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
      t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
      return t_minus_tepoch*t_minus_tepoch2*t_minus_tepoch4/5040.0;
    }else if(parameternr == TEMPO3_PARAM_F0 + 7) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
      t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
      t_minus_tepoch8 = t_minus_tepoch4*t_minus_tepoch4;
      return t_minus_tepoch8/40320.0;
    }else if(parameternr == TEMPO3_PARAM_F0 + 8) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
      t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
      t_minus_tepoch8 = t_minus_tepoch4*t_minus_tepoch4;
      return t_minus_tepoch*t_minus_tepoch8/362880.0;
    }else if(parameternr == TEMPO3_PARAM_F0 + 9) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
      t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
      t_minus_tepoch8 = t_minus_tepoch4*t_minus_tepoch4;
      return t_minus_tepoch2*t_minus_tepoch8/3628800.0;
    }else if(parameternr == TEMPO3_PARAM_F0 + 10) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
      t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
      t_minus_tepoch8 = t_minus_tepoch4*t_minus_tepoch4;
      return t_minus_tepoch*t_minus_tepoch2*t_minus_tepoch8/39916800.0;
    }else if(parameternr == TEMPO3_PARAM_F0 + 11) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
      t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
      t_minus_tepoch8 = t_minus_tepoch4*t_minus_tepoch4;
      return t_minus_tepoch4*t_minus_tepoch8/479001600.0;
    }else if(parameternr == TEMPO3_PARAM_F0 + 12) {
      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
      t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
      t_minus_tepoch8 = t_minus_tepoch4*t_minus_tepoch4;
      return t_minus_tepoch*t_minus_tepoch4*t_minus_tepoch8/6227020800.0;
    }else {
      //  Following is not tested:
      //      int order = parameternr-TEMPO3_PARAM_F0;
      //      t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
      //      return powl(t_minus_tepoch, order+1, order)/(long double)tgamma(order+1.0);
      printerror(0, "ERROR tempo3: the %dth frequency derivative is not (yet?) implemented in calc_derivative()\n", parameternr-TEMPO3_PARAM_F0);
      exit(0);
    }
  }
  if(parameternr >= TEMPO3_PARAM_WAVESIN && parameternr < TEMPO3_PARAM_WAVESIN + TEMPO3_MaxNrWaves) {  // The parameter is wave parameter
    deriv = sinl(parfile_tmp->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(parameternr-TEMPO3_PARAM_WAVESIN+1)*(mjd-parfile_tmp->parameter[TEMPO3_PARAM_WAVEEPOCH]));
    deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
    return deriv;
  }
  if(parameternr >= TEMPO3_PARAM_WAVECOS && parameternr < TEMPO3_PARAM_WAVECOS + TEMPO3_MaxNrWaves) {  // The parameter is wave parameter
    deriv = cosl(parfile_tmp->parameter[TEMPO3_PARAM_WAVEOM]*(long double)(parameternr-TEMPO3_PARAM_WAVECOS+1)*(mjd-parfile_tmp->parameter[TEMPO3_PARAM_WAVEEPOCH]));
    deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
    return deriv;
  }

  if(parameternr >= TEMPO3_PARAM_GLPH && parameternr < TEMPO3_PARAM_GLPH + TEMPO3_MaxNrGlitches) {
    int glitchnr = parameternr - TEMPO3_PARAM_GLPH;
    long double glep = parfile_tmp->parameter[TEMPO3_PARAM_GLEP+glitchnr];
    if(mjd >= glep) {
       return 1.0;
     }else {
       return 0.0;
     }
  }
  if(parameternr >= TEMPO3_PARAM_GLF0 && parameternr < TEMPO3_PARAM_GLF0 + TEMPO3_MaxNrGlitches) {   // GLF0
    int glitchnr = parameternr - TEMPO3_PARAM_GLF0;
    long double t_minus_tglep = mjd-parfile_tmp->parameter[TEMPO3_PARAM_GLEP+glitchnr];
    if(t_minus_tglep >= 0) {
       return t_minus_tglep*TEMPO3_SecPerDay;
     }else {
       return 0.0;
     }
  }
  if(parameternr >= TEMPO3_PARAM_GLF0+1*TEMPO3_MaxNrGlitches && parameternr < TEMPO3_PARAM_GLF0+2*TEMPO3_MaxNrGlitches) {   // GLF1
    int glitchnr = parameternr - TEMPO3_PARAM_GLF0 - 1*TEMPO3_MaxNrGlitches;
    long double t_minus_tglep = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_GLEP+glitchnr])*TEMPO3_SecPerDay;
    long double t_minus_tglep2 = t_minus_tglep*t_minus_tglep;
    if(t_minus_tglep >= 0) {
       return t_minus_tglep2/2.0;
     }else {
       return 0.0;
     }
  }
  if(parameternr >= TEMPO3_PARAM_GLF0+2*TEMPO3_MaxNrGlitches && parameternr < TEMPO3_PARAM_GLF0+3*TEMPO3_MaxNrGlitches) {   // GLF2
    int glitchnr = parameternr - TEMPO3_PARAM_GLF0 - 2*TEMPO3_MaxNrGlitches;
    long double t_minus_tglep = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_GLEP+glitchnr])*TEMPO3_SecPerDay;
    long double t_minus_tglep2 = t_minus_tglep*t_minus_tglep;
    long double t_minus_tglep3 = t_minus_tglep2*t_minus_tglep;
    if(t_minus_tglep >= 0) {
       return t_minus_tglep3/6.0;
     }else {
       return 0.0;
     }
  }
  if(parameternr >= TEMPO3_PARAM_GLF0+3*TEMPO3_MaxNrGlitches && parameternr < TEMPO3_PARAM_GLF0+4*TEMPO3_MaxNrGlitches) {   // GLF3
    int glitchnr = parameternr - TEMPO3_PARAM_GLF0 - 3*TEMPO3_MaxNrGlitches;
    long double t_minus_tglep = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_GLEP+glitchnr])*TEMPO3_SecPerDay;
    long double t_minus_tglep2 = t_minus_tglep*t_minus_tglep;
    long double t_minus_tglep4 = t_minus_tglep2*t_minus_tglep2;
    if(t_minus_tglep >= 0) {
       return t_minus_tglep4/24.0;
     }else {
       return 0.0;
     }
  }
  if(parameternr >= TEMPO3_PARAM_GLF0+4*TEMPO3_MaxNrGlitches && parameternr < TEMPO3_PARAM_GLF0+5*TEMPO3_MaxNrGlitches) {   // GLF4
    int glitchnr = parameternr - TEMPO3_PARAM_GLF0 - 4*TEMPO3_MaxNrGlitches;
    long double t_minus_tglep = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_GLEP+glitchnr])*TEMPO3_SecPerDay;
    long double t_minus_tglep2 = t_minus_tglep*t_minus_tglep;
    long double t_minus_tglep4 = t_minus_tglep2*t_minus_tglep2;
    if(t_minus_tglep >= 0) {
       return t_minus_tglep4*t_minus_tglep/120.0;
     }else {
       return 0.0;
     }
  }
  if(parameternr >= TEMPO3_PARAM_GLF0+5*TEMPO3_MaxNrGlitches && parameternr < TEMPO3_PARAM_GLF0+6*TEMPO3_MaxNrGlitches) {   // GLF5
    int glitchnr = parameternr - TEMPO3_PARAM_GLF0 - 5*TEMPO3_MaxNrGlitches;
    long double t_minus_tglep = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_GLEP+glitchnr])*TEMPO3_SecPerDay;
    long double t_minus_tglep2 = t_minus_tglep*t_minus_tglep;
    long double t_minus_tglep4 = t_minus_tglep2*t_minus_tglep2;
    if(t_minus_tglep >= 0) {
       return t_minus_tglep4*t_minus_tglep2/720.0;
     }else {
       return 0.0;
     }
  }
  if(parameternr >= TEMPO3_PARAM_GLF0D && parameternr < TEMPO3_PARAM_GLF0D+TEMPO3_MaxNrGlitches) {   // GLF0D
    int glitchnr = parameternr - TEMPO3_PARAM_GLF0D;
    long double gltd = parfile_tmp->parameter[TEMPO3_PARAM_GLTD+glitchnr]*TEMPO3_SecPerDay;
    long double t_minus_tglep = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_GLEP+glitchnr])*TEMPO3_SecPerDay;
    if(t_minus_tglep >= 0) {
       return gltd*(1.0-expl(-t_minus_tglep/gltd));
     }else {
       return 0.0;
     }
  }
  if(parameternr >= TEMPO3_PARAM_GLTD && parameternr < TEMPO3_PARAM_GLTD+TEMPO3_MaxNrGlitches) {   // GLTD
    int glitchnr = parameternr - TEMPO3_PARAM_GLTD;
    long double gltd = parfile_tmp->parameter[TEMPO3_PARAM_GLTD+glitchnr];
    long double t_minus_tglep = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_GLEP+glitchnr]);
    long double ratio = t_minus_tglep/gltd;
    long double exponent = expl(-ratio);
    if(t_minus_tglep >= 0) {
       return parfile_tmp->parameter[TEMPO3_PARAM_GLF0D+glitchnr]*(1.0-exponent-ratio*exponent)*TEMPO3_SecPerDay;
     }else {
       return 0.0;
     }
  }
  if(parameternr >= TEMPO3_PARAM_GLEP && parameternr < TEMPO3_PARAM_GLEP+TEMPO3_MaxNrGlitches) {   // GLEP
    int glitchnr = parameternr - TEMPO3_PARAM_GLEP;
    long double gltd = parfile_tmp->parameter[TEMPO3_PARAM_GLTD+glitchnr]*TEMPO3_SecPerDay;
    long double t_minus_tglep = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_GLEP+glitchnr])*TEMPO3_SecPerDay;
    long double t_minus_tglep2 = t_minus_tglep*t_minus_tglep;
    long double t_minus_tglep4 = t_minus_tglep2*t_minus_tglep2;
    long double ratio = t_minus_tglep/gltd;
    long double exponent = expl(-ratio);
    if(t_minus_tglep >= 0) {

      //      phase = parfile_tmp->parameter[TEMPO3_PARAM_GLF0+glitchnr]*t_minus_tglep;  // GLF0
      deriv = -parfile_tmp->parameter[TEMPO3_PARAM_GLF0+glitchnr];
      //      phase += parfile_tmp->parameter[TEMPO3_PARAM_GLF0+1*TEMPO3_MaxNrGlitches+glitchnr]*t_minus_tglep2/2.0; //GLF1
      deriv -= parfile_tmp->parameter[TEMPO3_PARAM_GLF0+1*TEMPO3_MaxNrGlitches+glitchnr]*t_minus_tglep;
      //      phase += parfile_tmp->parameter[TEMPO3_PARAM_GLF0+2*TEMPO3_MaxNrGlitches+glitchnr]*t_minus_tglep3/6.0; //GLF2
      deriv -= parfile_tmp->parameter[TEMPO3_PARAM_GLF0+2*TEMPO3_MaxNrGlitches+glitchnr]*t_minus_tglep2/2.0; //GLF2
      //      phase += parfile_tmp->parameter[TEMPO3_PARAM_GLF0+3*TEMPO3_MaxNrGlitches+glitchnr]*t_minus_tglep3*t_minus_tglep/24.0; //GLF3
      deriv -= parfile_tmp->parameter[TEMPO3_PARAM_GLF0+3*TEMPO3_MaxNrGlitches+glitchnr]*t_minus_tglep2*t_minus_tglep/6.0; 
      //      phase += parfile_tmp->parameter[TEMPO3_PARAM_GLF0+4*TEMPO3_MaxNrGlitches+glitchnr]*t_minus_tglep3*t_minus_tglep2/120.0; // GLF4
      deriv -= parfile_tmp->parameter[TEMPO3_PARAM_GLF0+4*TEMPO3_MaxNrGlitches+glitchnr]*t_minus_tglep4/24.0; // GLF4
      //      phase += parfile_tmp->parameter[TEMPO3_PARAM_GLF0+5*TEMPO3_MaxNrGlitches+glitchnr]*t_minus_tglep3*t_minus_tglep3/720.0; // GLF5
      deriv -= parfile_tmp->parameter[TEMPO3_PARAM_GLF0+5*TEMPO3_MaxNrGlitches+glitchnr]*t_minus_tglep4*t_minus_tglep/120.0; // GLF5
      //      phase += parfile_tmp->parameter[TEMPO3_PARAM_GLF0D+glitchnr]*gltd*(1.0-expl(-t_minus_tglep/gltd)); // GLF0D
      deriv -= parfile_tmp->parameter[TEMPO3_PARAM_GLF0D+glitchnr]*exponent; // GLF0D
      deriv *= TEMPO3_SecPerDay;
      return deriv;
     }else {
       return 0.0;
     }
  }

  if(parameternr >= TEMPO3_PARAM_JUMPS && parameternr < TEMPO3_PARAM_JUMPS+TEMPO3_MaxNrJumps) {   // JUMP
    int jumpnr = parameternr - TEMPO3_PARAM_JUMPS;
    if(parfile_tmp->jumpmeaning[jumpnr] == 0) {  // Search for flag
      if(strstr(flags, parfile_tmp->jumpflag[jumpnr]) != NULL) {
	deriv = 1;
	return deriv*parfile_tmp->parameter[TEMPO3_PARAM_F0];
      }else {
	return 0.0;
      }
    }else if(parfile_tmp->jumpmeaning[jumpnr] == 1) {  // Search for site
      if(strcmp(site, parfile_tmp->jumpflag[jumpnr]) == 0) {
	deriv -= 1.0;  // CONFUSINGLY, THE TELESCOPE JUMP IS SUBTRACTED, BUT THE FLAG JUMP IS ADDED?????
	return deriv*parfile_tmp->parameter[TEMPO3_PARAM_F0];
      }else {
	return 0.0;
      }
    }
  }


  if(parameternr == TEMPO3_PARAM_RAJ || parameternr == TEMPO3_PARAM_DECJ || parameternr == TEMPO3_PARAM_PMRA || parameternr == TEMPO3_PARAM_PMDEC) {
    if(fabsl(ssb1) > 1 || fabsl(ssb2) > 1 || fabsl(ssb3) > 1) {
      long double x[3], ra, dec, constant;

      // This bit only needs to be once: undo the barycentre/proper motion bit. It is time dependent, so it needs to be done for each toa. In principle the mjd is changed by jumps etc. I assume those effects can be neglected for any reasonable situation.

      ra = parfile_tmp->raj;
      dec = parfile_tmp->decj;
      if(parfile_tmp->parameter[TEMPO3_PARAM_PMRA] != 0.0 || parfile_tmp->parameter[TEMPO3_PARAM_PMDEC] != 0.0 || parfile_tmp->pmra != 0.0 || parfile_tmp->pmdec != 0.0) {
	// pmra [mas/yr] * dt [yr] * (1e-3 / 3600.0) [deg/mas]
	//	constant = ((mjd/TEMPO3_SecPerDay - parfile_tmp->parameter[TEMPO3_PARAM_POSEPOCH])/365.25 )*(1e-3/3600.0)*(M_PI/180.0);
	constant = ((mjd - parfile_tmp->parameter[TEMPO3_PARAM_POSEPOCH])/365.25 )*(1e-3/3600.0)*(M_PI/180.0);
      	ra  += parfile_tmp->pmra  * constant;
      	dec += parfile_tmp->pmdec * constant;
      }

      //      x[0] = cosl(dec)*cosl(ra);
      //      x[1] = cosl(dec)*sinl(ra);
      //      x[2] = sinl(dec);
      //      dt = ssb1*x[0] + ssb2*x[1] + ssb3*x[2];

      ra  += parfile_tmp->parameter[TEMPO3_PARAM_RAJ];
      dec += parfile_tmp->parameter[TEMPO3_PARAM_DECJ];
      if(parfile_tmp->parameter[TEMPO3_PARAM_PMRA] != 0.0 || parfile_tmp->parameter[TEMPO3_PARAM_PMDEC] != 0.0) {
      	ra  += parfile_tmp->parameter[TEMPO3_PARAM_PMRA]  * constant;
      	dec += parfile_tmp->parameter[TEMPO3_PARAM_PMDEC] * constant;
      }
      //      x[0] = cosl(dec)*cosl(ra);
      //      x[1] = cosl(dec)*sinl(ra);
      //      x[2] = sinl(dec);
      //      dt -= ssb1*x[0] + ssb2*x[1] + ssb3*x[2];

      if(parameternr == TEMPO3_PARAM_RAJ) {
	x[0] = -cosl(dec)*sinl(ra);
	x[1] = cosl(dec)*cosl(ra);
	deriv = ssb1*x[0] + ssb2*x[1];
      }else if(parameternr == TEMPO3_PARAM_DECJ) {
	x[0] = -sinl(dec)*cosl(ra);
	x[1] = -sinl(dec)*sinl(ra);
	x[2] = cosl(dec);
	deriv = ssb1*x[0] + ssb2*x[1] + ssb3*x[2];
      }else if(parameternr == TEMPO3_PARAM_PMRA) {
	x[0] = -cosl(dec)*sinl(ra);
	x[1] = cosl(dec)*cosl(ra);
	deriv = ssb1*x[0] + ssb2*x[1];
	deriv *= constant;
      }else if(parameternr == TEMPO3_PARAM_PMDEC) {
	x[0] = -sinl(dec)*cosl(ra);
	x[1] = -sinl(dec)*sinl(ra);
	x[2] = cosl(dec);
	deriv = ssb1*x[0] + ssb2*x[1] + ssb3*x[2];
	deriv *= constant;
      }
    }else {
      printwarning(0, "WARNING: Position Earth w.r.t. barycenter not known, cannot fit for position/proper motion.");
    }
    deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
    return deriv;
  }

  if(parameternr >= TEMPO3_PARAM_DM && parameternr <= TEMPO3_PARAM_DM+TEMPO3_MaxNrDMderivatives) {
    deriv = -DM_CONST/(freq*freq);
    if(parameternr == TEMPO3_PARAM_DM) {
      deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
    }else if(parameternr == TEMPO3_PARAM_DM+1) {
      if(parfile_tmp->parameter[TEMPO3_PARAM_DM+1] != 0) {
	long double t_minus_tepoch;
	t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_DMEPOCH])/365.25;
	deriv *= t_minus_tepoch;
	deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      }
    }else if(parameternr == TEMPO3_PARAM_DM+2) {
      if(parfile_tmp->parameter[TEMPO3_PARAM_DM+2] != 0) {
	long double t_minus_tepoch;
	t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_DMEPOCH])/365.25;
	deriv *= t_minus_tepoch*t_minus_tepoch/2.0;    // 2!
	deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      }
    }else if(parameternr == TEMPO3_PARAM_DM+3) {
      if(parfile_tmp->parameter[TEMPO3_PARAM_DM+3] != 0) {
	long double t_minus_tepoch, t_minus_tepoch2;
	t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_DMEPOCH])/365.25;
	t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
	deriv *= t_minus_tepoch2*t_minus_tepoch/6.0;    // 3!
	deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      }
    }else if(parameternr == TEMPO3_PARAM_DM+4) {
      if(parfile_tmp->parameter[TEMPO3_PARAM_DM+4] != 0) {
	long double t_minus_tepoch, t_minus_tepoch2;
	t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_DMEPOCH])/365.25;
	t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
	deriv *= t_minus_tepoch2*t_minus_tepoch2/24.0;    // 4!
	deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      }
    }else if(parameternr == TEMPO3_PARAM_DM+5) {
      if(parfile_tmp->parameter[TEMPO3_PARAM_DM+5] != 0) {
	long double t_minus_tepoch, t_minus_tepoch2, t_minus_tepoch4;
	t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_DMEPOCH])/365.25;
	t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
	t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
	deriv *= t_minus_tepoch4*t_minus_tepoch/120.0;    // 5!
	deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      }
    }else if(parameternr == TEMPO3_PARAM_DM+6) {
      if(parfile_tmp->parameter[TEMPO3_PARAM_DM+6] != 0) {
	long double t_minus_tepoch, t_minus_tepoch2, t_minus_tepoch4, t_minus_tepoch6;
	t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_DMEPOCH])/365.25;
	t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
	t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
	t_minus_tepoch6 = t_minus_tepoch4*t_minus_tepoch2;
	deriv *= t_minus_tepoch6/720.0;    // 6!
	deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      }
    }else if(parameternr == TEMPO3_PARAM_DM+7) {
      if(parfile_tmp->parameter[TEMPO3_PARAM_DM+7] != 0) {
	long double t_minus_tepoch, t_minus_tepoch2, t_minus_tepoch4, t_minus_tepoch6;
	t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_DMEPOCH])/365.25;
	t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
	t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
	t_minus_tepoch6 = t_minus_tepoch4*t_minus_tepoch2;
	deriv *= t_minus_tepoch6*t_minus_tepoch/5040.0;    // 7!
	deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      }
    }else if(parameternr == TEMPO3_PARAM_DM+8) {
      if(parfile_tmp->parameter[TEMPO3_PARAM_DM+8] != 0) {
	long double t_minus_tepoch, t_minus_tepoch2, t_minus_tepoch4, t_minus_tepoch6;
	t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_DMEPOCH])/365.25;
	t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
	t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
	t_minus_tepoch6 = t_minus_tepoch4*t_minus_tepoch2;
	deriv *= t_minus_tepoch6*t_minus_tepoch2/40320.0;    // 8!
	deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      }
    }else if(parameternr == TEMPO3_PARAM_DM+9) {
      if(parfile_tmp->parameter[TEMPO3_PARAM_DM+9] != 0) {
	long double t_minus_tepoch, t_minus_tepoch2, t_minus_tepoch4, t_minus_tepoch6;
	t_minus_tepoch = (mjd-parfile_tmp->parameter[TEMPO3_PARAM_DMEPOCH])/365.25;
	t_minus_tepoch2 = t_minus_tepoch*t_minus_tepoch;
	t_minus_tepoch4 = t_minus_tepoch2*t_minus_tepoch2;
	t_minus_tepoch6 = t_minus_tepoch4*t_minus_tepoch2;
	deriv *= t_minus_tepoch6*t_minus_tepoch2*t_minus_tepoch/362880.0;    // 9!
	deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      }
    }else {
      printerror(0, "ERROR tempo3: Analytic derivative to fit parameter is not implemented.\n");
      exit(0);
    }
    return deriv;
  }

  if(parfile_tmp->binarymodel == 3) {
    if(parameternr >= TEMPO3_PARAM_A1 && parameternr < TEMPO3_PARAM_A1 + TEMPO3_MaxNrCompanions) {
      int companionnr;
      companionnr = parameternr - TEMPO3_PARAM_A1;
      //    printf("XXXX: mjd=%Le\n", mjd);
      deriv = ephemeris_get_torb(mjd*TEMPO3_SecPerDay, parfile_tmp, paramset, 1, companionnr);
      deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      return deriv;
    }else if(parameternr >= TEMPO3_PARAM_PB && parameternr < TEMPO3_PARAM_PB + TEMPO3_MaxNrCompanions) {
      int companionnr;
      companionnr = parameternr - TEMPO3_PARAM_PB;
      deriv = ephemeris_get_torb(mjd*TEMPO3_SecPerDay, parfile_tmp, paramset, 2, companionnr);
      deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      return deriv;
    }else if(parameternr >= TEMPO3_PARAM_T0 && parameternr < TEMPO3_PARAM_T0 + TEMPO3_MaxNrCompanions) {
      int companionnr;
      companionnr = parameternr - TEMPO3_PARAM_T0;
      deriv = ephemeris_get_torb(mjd*TEMPO3_SecPerDay, parfile_tmp, paramset, 3, companionnr);
      deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      return deriv;
    }else if(parameternr >= TEMPO3_PARAM_OM && parameternr < TEMPO3_PARAM_OM + TEMPO3_MaxNrCompanions) {
      int companionnr;
      companionnr = parameternr - TEMPO3_PARAM_OM;
      deriv = ephemeris_get_torb(mjd*TEMPO3_SecPerDay, parfile_tmp, paramset, 4, companionnr);
      deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      return deriv;
    }else if(parameternr >= TEMPO3_PARAM_ECC && parameternr < TEMPO3_PARAM_ECC + TEMPO3_MaxNrCompanions) {
      int companionnr;
      companionnr = parameternr - TEMPO3_PARAM_ECC;
      deriv = ephemeris_get_torb(mjd*TEMPO3_SecPerDay, parfile_tmp, paramset, 5, companionnr);
      deriv *= parfile_tmp->parameter[TEMPO3_PARAM_F0];
      return deriv;
    }
  }



  if(dx <= 0) {
    printerror(0, "ERROR tempo3: Analytic derivative to fit parameter is not implemented.\n");
    exit(0);
  }

  use_xfit(parfile_tmp, parfile_params); 
  ephemeris_setLinkedValues_par_version(parfile_tmp, paramlinked);
  dy = -evaluate_ephemeris(mjd, freq, ssb1, ssb2, ssb3, 0, 0, site, flags, parfile_tmp, 0, paramset, 0, 0);
  tmp_value = parfile_params[parameternr];
  /*  dx = 1e-3; */
  //  dx = 1e-3;
  //  printf("XXXX increasing parameter %Le with %Le\n", parfile_params[parameternr], dx);
  parfile_params[parameternr] += dx;
  use_xfit(parfile_tmp, parfile_params);
  ephemeris_setLinkedValues_par_version(parfile_tmp, paramlinked);
  parfile_params[parameternr] = tmp_value;
  dy += evaluate_ephemeris(mjd, freq, ssb1, ssb2, ssb3, 0, 0, site, flags, parfile_tmp, 0, paramset, 0, 0);
  deriv2 = dy/dx;

  // Leave out exit(0) after the error message and do a deriv=.... rather than a return ... to check analytic derivative with nummerical equivalent
  //  if(parameternr != TEMPO3_PARAM_PHASE0)
  //    printf("Check of analytic derivative %Le which should be compared with the numerical derivative %Le\nThe difference is %Lf%%\n", deriv, deriv2, 100.0*(deriv-deriv2)/deriv2);
  return deriv2;
}




#ifdef TEMPO3_EnableGTK
  GtkWidget *gtkLine_param[MaxNrGTKLines];
  GtkWidget *gtkLine_value[MaxNrGTKLines];
  GtkWidget *gtkLine_fit[MaxNrGTKLines];
  GtkWidget *gtkLine_diff[MaxNrGTKLines];
  GtkWidget *gtkLine_err[MaxNrGTKLines];
  int gtkLine_ident[MaxNrGTKLines];
  int enableGTK_flag;
#ifdef HAVEPGPLOT2 
  int enablePGPLOT2;
  int deviceID_pgplot2;
//  GtkWidget *gtk_drawing_area;
//  GtkWidget *gtk_menuBar;
#endif
  int *gtk_fixed;
  tempo3_parameters_def gtk_parfile;
  tempo3_parameters_def *gtk_parfile2;
  toa_data_def *gtk_toas;
  tempo3_wraps_def *gtk_wraps;
  int * gtk_paramset, gtk_finderrors, *gtk_paramlinked;
  long double *gtk_dplus, *gtk_dmin;

// Keeps track if we're in the process of generating gtk function
// calls, which means the wakeup call shouldn't interupt with
// processing pending get operations which could disrupt
// things. Likewise, when threading happens, no gtk calls should be
// made in the alarm_wakeup function.
  int gtk_iteration_blocked = 0;

int alarm_active = 0;

/* ident indicates the number of spaces before the lines. If parfile2
   is set, the difference will be quoted as well. If otherparams is
   set the unrecognized lines will be written out as well. */
void print_ephemeris_gtk(FILE *stream, tempo3_parameters_def parfile, tempo3_parameters_def *parfile2, int *fixed, int ident, int otherparams, toa_data_def *toas, tempo3_wraps_def *wraps, char **unused_parfile_lines, int nrunused_parfile_lines, int *paramset, int finderrors, long double *dplus, long double *dmin, int *paramlinked)
{
  char textBuffer[1000];
  int i, linenr, curparam;

  gtk_iteration_blocked++;
  while(gtk_iteration_blocked > 1) {
    fprintf(stderr, "XXXX waiting until free to process gtk events %d\n", alarm_active);
    sleep(1);
  }
  while(gtk_events_pending()) {
    fprintf(stderr, "XXXX process gtk events at start of print_ephemeris\n");
    gtk_main_iteration();
  }

  ephemeris_setLinkedValues_par_version(&parfile, paramlinked);

  linenr = 0;
  /* Set the identities of all buttons to unknown, as their order can change */
  for(i=0; i < nrActiveGTKLines; i++) 
    gtkLine_ident[i] = -1;

  // Display the parameters
  for(curparam = 0; curparam < TEMPO3_PARAM_TOTAL; curparam++) {
    if(paramset[curparam] && 
       (curparam < TEMPO3_PARAM_WAVECOS || curparam >= TEMPO3_PARAM_WAVECOS + TEMPO3_MaxNrWaves)  // Ignore cos terms, as they were already dealt with when outputting the sin terms
       && curparam != TEMPO3_PARAM_TRES && curparam != TEMPO3_PARAM_BINARY && curparam != TEMPO3_PARAM_MODE && curparam != TEMPO3_PARAM_TRACK && curparam != TEMPO3_PARAM_TZRMJD && curparam != TEMPO3_PARAM_TZRFRQ && curparam != TEMPO3_PARAM_TZRSITE   // Don't output all non-fit parameters in gtk window to avoid cluttering
      ) {
      if(enableGTK_flag && linenr < nrActiveGTKLines) {
	gtkLine_ident[linenr] = curparam;
	if(curparam >= TEMPO3_PARAM_JUMPS && curparam < TEMPO3_PARAM_JUMPS + TEMPO3_MaxNrJumps) {
	  if(parfile.jumpmeaning[curparam - TEMPO3_PARAM_JUMPS] == 0) { // Search for flag
	    sprintf(textBuffer, "%s %s", tempo3_parfile_descr.identifiers[curparam][0], parfile.jumpflag[curparam - TEMPO3_PARAM_JUMPS]);
	  }else if(parfile.jumpmeaning[curparam - TEMPO3_PARAM_JUMPS] == 1) { // Search for site
	    sprintf(textBuffer, "%s TEL %s", tempo3_parfile_descr.identifiers[curparam][0], parfile.jumpflag[curparam - TEMPO3_PARAM_JUMPS]);
	  }else {
	    printerror(0, "ERROR: The jumpmeaning=%d is not defined. This is a bug.", parfile.jumpmeaning[curparam - TEMPO3_PARAM_JUMPS]);
	    exit(0);
	  }
	  gtk_entry_set_text(GTK_ENTRY(gtkLine_param[linenr]), textBuffer);
	}else {
	  gtk_entry_set_text(GTK_ENTRY(gtkLine_param[linenr]), tempo3_parfile_descr.identifiers[curparam][0]);
	}
	int value_id = ephemeris_check_if_parameter_is_linked(curparam, paramlinked);
	if(curparam == TEMPO3_PARAM_PSRJNAME) {
	  gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), parfile.psrjname);
	}else if(curparam == TEMPO3_PARAM_RAJ) {
	  converthms_string(textBuffer, (parfile.parameter[curparam]+parfile.raj)*12.0/M_PI, 19, 1);
	  gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), textBuffer);
	}else if(curparam == TEMPO3_PARAM_DECJ) {
	  converthms_string(textBuffer, (parfile.parameter[curparam]+parfile.decj)*180.0/M_PI, 19, 1);
	  gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), textBuffer);
	}else if(curparam == TEMPO3_PARAM_PMRA) {
	  sprintf(textBuffer, tempo3_parfile_descr.print_format[curparam], parfile.pmra + parfile.parameter[curparam]);
	  gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), textBuffer);
	}else if(curparam == TEMPO3_PARAM_PMDEC) {
	  sprintf(textBuffer, tempo3_parfile_descr.print_format[curparam], parfile.pmdec + parfile.parameter[curparam]);
	  gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), textBuffer);
	}else if(curparam == TEMPO3_PARAM_TZRSITE) {
	  gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), parfile.tzrsite);
	}else if(curparam == TEMPO3_PARAM_BINARY) {
	  if(parfile.binarymodel == 1) {
	    gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), "BT");
	  }else if(parfile.binarymodel == 2) {
	    gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), "T2");
	  }else if(parfile.binarymodel == 3) {
	    gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), "HYPER");
	  }else {
	    gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), "NOT SUPPORTED");
	  }
	}else if(curparam >= TEMPO3_PARAM_DM && curparam <= TEMPO3_PARAM_DM+TEMPO3_MaxNrDMderivatives) {
	  sprintf(textBuffer, tempo3_parfile_descr.print_format[curparam], parfile.dm[curparam - TEMPO3_PARAM_DM] + parfile.parameter[curparam]);
	  gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), textBuffer);
	}else if(curparam >= TEMPO3_PARAM_WAVESIN && curparam < TEMPO3_PARAM_WAVESIN + TEMPO3_MaxNrWaves) {
	  sprintf(textBuffer, "%-12.8Le %-12.8Le", parfile.parameter[curparam], parfile.parameter[curparam-TEMPO3_PARAM_WAVESIN+TEMPO3_PARAM_WAVECOS]);  // Note, didn't use the set print_format, as this one is shorter than that in the parfile output
	  gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), textBuffer);
	}else {
	  if(value_id == 0) {
	    sprintf(textBuffer, tempo3_parfile_descr.print_format[curparam], parfile.parameter[curparam]);
	  }else {
	    if(value_id > 0) {
	      sprintf(textBuffer, "VALUE%d", value_id);
	    }else {
	      sprintf(textBuffer, "-VALUE%d", value_id);
	    }
	  }
	  gtk_entry_set_text(GTK_ENTRY(gtkLine_value[linenr]), textBuffer);
	}
	if(value_id == 0 && curparam != TEMPO3_PARAM_PSRJNAME && curparam != TEMPO3_PARAM_WAVEOM && curparam != TEMPO3_PARAM_PEPOCH && curparam != TEMPO3_PARAM_WAVEEPOCH && curparam != TEMPO3_PARAM_TRES && curparam != TEMPO3_PARAM_BINARY && (curparam < TEMPO3_PARAM_GLEPMIN || curparam >= TEMPO3_PARAM_GLEPMIN+TEMPO3_MaxNrGlitches) && (curparam < TEMPO3_PARAM_GLEPMAX || curparam >= TEMPO3_PARAM_GLEPMAX+TEMPO3_MaxNrGlitches) && curparam != TEMPO3_PARAM_MODE && curparam != TEMPO3_PARAM_TRACK && curparam != TEMPO3_PARAM_TZRMJD && curparam != TEMPO3_PARAM_TZRFRQ && curparam != TEMPO3_PARAM_TZRSITE && curparam != TEMPO3_PARAM_POSEPOCH && curparam != TEMPO3_PARAM_DMEPOCH) {
	  if(fixed != NULL) {
	    if(fixed[curparam]) 
	      gtk_entry_set_text(GTK_ENTRY(gtkLine_fit[linenr]), "fixed");
	    else
	      gtk_entry_set_text(GTK_ENTRY(gtkLine_fit[linenr]), "FIT");
	  }
	  if(parfile2 != NULL) {
	    if(curparam >= TEMPO3_PARAM_DM && curparam <= TEMPO3_PARAM_DM+TEMPO3_MaxNrDMderivatives) {
	      sprintf(textBuffer, "%11.1Le", parfile.dm[curparam - TEMPO3_PARAM_DM] + parfile.parameter[curparam] - parfile2->dm[curparam - TEMPO3_PARAM_DM] - parfile2->parameter[curparam]);
	    }else if(curparam == TEMPO3_PARAM_PMRA) {
	      sprintf(textBuffer, "%11.1Le", parfile.pmra + parfile.parameter[curparam] - parfile2->pmra - parfile2->parameter[curparam]);
	    }else if(curparam == TEMPO3_PARAM_PMDEC) {
	      sprintf(textBuffer, "%11.1Le", parfile.pmdec + parfile.parameter[curparam] - parfile2->pmdec - parfile2->parameter[curparam]);
	    }else if(curparam >= TEMPO3_PARAM_WAVESIN && curparam < TEMPO3_PARAM_WAVESIN + TEMPO3_MaxNrWaves) {
	      sprintf(textBuffer, "%11.1Le %11.1Le", parfile.parameter[curparam] - parfile2->parameter[curparam], parfile.parameter[curparam-TEMPO3_PARAM_WAVESIN+TEMPO3_PARAM_WAVECOS] - parfile2->parameter[curparam-TEMPO3_PARAM_WAVESIN+TEMPO3_PARAM_WAVECOS]);
	    }else {
	      sprintf(textBuffer, "%11.1Le", parfile.parameter[curparam]-parfile2->parameter[curparam]);
	    }
	    gtk_entry_set_text(GTK_ENTRY(gtkLine_diff[linenr]), textBuffer);
	  }
	  if(finderrors) {
	    if(curparam >= TEMPO3_PARAM_DM && curparam <= TEMPO3_PARAM_DM+TEMPO3_MaxNrDMderivatives) {
	      sprintf(textBuffer, "%12.2Le (%.2Le)", dplus[curparam], fabsl((parfile.dm[curparam - TEMPO3_PARAM_DM] + parfile.parameter[curparam])/dplus[curparam]));
	    }else if(curparam == TEMPO3_PARAM_PMRA) {
	      sprintf(textBuffer, "%12.2Le (%.2Le)", dplus[curparam], fabsl((parfile.pmra + parfile.parameter[curparam])/dplus[curparam]));
	    }else if(curparam == TEMPO3_PARAM_PMDEC) {
	      sprintf(textBuffer, "%12.2Le (%.2Le)", dplus[curparam], fabsl((parfile.pmdec + parfile.parameter[curparam])/dplus[curparam]));
	    }else if(curparam >= TEMPO3_PARAM_WAVESIN && curparam < TEMPO3_PARAM_WAVESIN + TEMPO3_MaxNrWaves) {
	      sprintf(textBuffer, "%12.2Le %12.2Le)", dplus[curparam], dplus[curparam-TEMPO3_PARAM_WAVESIN+TEMPO3_PARAM_WAVECOS]);
	    }else if(curparam == TEMPO3_PARAM_RAJ) {
	      sprintf(textBuffer, "%12.2Le (sec)", dplus[curparam]*12.0*60.0*60.0/M_PI);
	    }else if(curparam == TEMPO3_PARAM_DECJ) {
	      sprintf(textBuffer, "%12.2Le (arcsec)", dplus[curparam]*180.0*60.0*60.0/M_PI);
	    }else {
	      sprintf(textBuffer, "%12.2Le (%.2Le)", dplus[curparam], fabsl(parfile.parameter[curparam]/dplus[curparam]));
	    }
	    gtk_entry_set_text(GTK_ENTRY(gtkLine_err[linenr]), textBuffer);
	  }
	}
	linenr++;
      }
    }
  }

  while(gtk_iteration_blocked > 1) {
    fprintf(stderr, "XXXX waiting until free to process gtk events\n");
    sleep(1);
  }
  while(gtk_events_pending()) {
    fprintf(stderr, "XXXX process gtk events at end of print_ephemeris\n");
    gtk_main_iteration();
  }
//  while(gtk_events_pending()) {
//    gtk_main_iteration();
//  }

  gtk_iteration_blocked--;
}

// Callback:  The data passed to this function is printed to stdout.
static void callback1( GtkWidget *widget,
                       int *   data )
//                       gpointer   data )
{
  gtk_iteration_blocked++;
  /*  printf ("Hello again - button %d was pressed\n", *data); */
  toggleFixed(*data, gtk_fixed, gtk_paramlinked);

  //  print_ephemeris_gtk(NULL, tempo3_parameters_def parfile, tempo3_parameters_def *parfile2, int *fixed, int ident, int otherparams, toa_data_def *toas, tempo3_wraps_def *wraps, char unused_parfile_lines[][TEMPO3_MaxNrParLineLength], int nrunused_parfile_lines, int *paramset, int finderrors, long double *dplus, long double *dmin, int *paramlinked);

  //  if(gtk_iteration_blocked == 0) {
    print_ephemeris_gtk(NULL, gtk_parfile, gtk_parfile2, gtk_fixed, 0, 0, gtk_toas, gtk_wraps, NULL, 0, gtk_paramset, gtk_finderrors, gtk_dplus, gtk_dmin, gtk_paramlinked);
    //  }

    //    gtk_widget_hide(greyBox);
    //gtk_widget_hide(redBox);
    //gtk_widget_show(greenBox);
    gtk_iteration_blocked--;
}

void createGTKtextbox(GtkWidget *hbox1, GtkWidget **text_entry_Widget, int text_width)
{
  gtk_iteration_blocked++;
  //  GdkColor colorBlack;
  //  GdkColor colorGreen;
  // black
  //  colorBlack.red=0;
  //  colorBlack.green=0;
  //  colorBlack.blue=0;
  
    // Green
  //  colorGreen.red=0;
  //  colorGreen.green=65535;
  //  colorGreen.blue=0;

  *text_entry_Widget = gtk_entry_new();
  gtk_widget_set_name(*text_entry_Widget, "text_entry");   // Set name to be used in css file

  // Replaced with css code (see gtk_css_provider_load_from_data)
  //  gtk_widget_modify_text(*text_entry_Widget, GTK_STATE_NORMAL, &colorGreen);
  //  gtk_widget_modify_base(*text_entry_Widget, GTK_STATE_NORMAL, &colorBlack);

  //  GtkStyle *style = gtk_widget_get_style(*text_entry_Widget);
  /* PANGO_WEIGHT_BOLD */
  //  pango_font_description_set_weight(style->font_desc, PANGO_WEIGHT_NORMAL);
  //  pango_font_description_set_size(style->font_desc, 9 * PANGO_SCALE);
  //  gtk_widget_modify_font(*text_entry_Widget, style->font_desc);

  gtk_entry_set_width_chars(GTK_ENTRY(*text_entry_Widget), text_width);

  gtk_editable_set_editable(GTK_EDITABLE(*text_entry_Widget), FALSE);


  //  GTK_WIDGET_UNSET_FLAGS(*text_entry_Widget, GTK_CAN_FOCUS);
  // GTK_WIDGET_UNSET_FLAGS has been deprecated since version 2.22 and should not be used in newly-written code.
  //Use the proper function instead. 
  //gtk_widget_set_app_paintable(), gtk_widget_set_can_default(), gtk_widget_set_can_focus(), gtk_widget_set_double_buffered(), gtk_widget_set_has_window(), gtk_widget_set_mapped(), gtk_widget_set_no_show_all(), gtk_widget_set_realized(), gtk_widget_set_receives_default(), gtk_widget_set_sensitive() or gtk_widget_set_visible().
  gtk_widget_set_can_focus(*text_entry_Widget, FALSE);


  char *txtBuffer = "";
  gtk_entry_set_text(GTK_ENTRY(*text_entry_Widget), txtBuffer);

  gtk_box_pack_start (GTK_BOX (hbox1), *text_entry_Widget, TRUE, TRUE, 0);
  gtk_widget_show (*text_entry_Widget);
  gtk_iteration_blocked--;
}

void createGTKline(GtkWidget *vBox_Widget, GtkWidget **text_entry_Widget_param, GtkWidget **text_entry_Widget_value, GtkWidget **text_entry_Widget_fit, GtkWidget **text_entry_Widget_diff, GtkWidget **text_entry_Widget_err, int linenr)
{
  gtk_iteration_blocked++;
  GtkWidget *button1;
  GtkWidget *hbox1;
  
  //  hbox1 = gtk_hbox_new (FALSE, 0);
  // Should use grid instead, but this is the quick gtk2 -> gtk3 fix
  hbox1 = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);

  gtk_container_add(GTK_CONTAINER (vBox_Widget), hbox1);
  gtk_widget_show(hbox1);

  /*
  GtkWidget *boxLabel_Widget;
  GtkWidget *textLabel_Widget;
  char       textLabel[64];
  char *textLabel0="Display Label";
  sprintf(textLabel, "<small>%s</small>",textLabel0);
  textLabel_Widget = gtk_label_new (textLabel);
  gtk_label_set_markup(GTK_LABEL (textLabel_Widget),textLabel);
  gtk_label_set_justify(GTK_LABEL (textLabel_Widget),GTK_JUSTIFY_LEFT);
  gtk_misc_set_alignment (GTK_MISC (textLabel_Widget), 0, 0.5);
  gtk_widget_show (textLabel_Widget);
  gtk_box_pack_start( GTK_BOX (hbox1),textLabel_Widget,TRUE,TRUE, 0);
  */
    // --------------------------------------------------------------------
    // Create text entry widget and modify.
    // --------------------------------------------------------------------


  /*
    Position both images in the same location
    gtk_misc_set_alignment(GTK_MISC(greyBox), 1, 0.5);
  */

  button1 = gtk_button_new_with_label("  Fit  ");
  gtk_widget_set_can_focus(button1, FALSE);  // This means that after clicking the button, pressing space (do fit) will not toggle the button instead/as well. Especially important for the pgplot2 local mode.
  gtkLine_ident[linenr] = -1;
  g_signal_connect(G_OBJECT (button1), "clicked", G_CALLBACK (callback1), &(gtkLine_ident[linenr]));
  //  gtk_box_pack_start(GTK_BOX(hbox1), button1, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(hbox1), button1, FALSE, FALSE, 0);   // Do not make buttons/text boxes to expand
  gtk_widget_show(button1);
  
  
  createGTKtextbox(hbox1, text_entry_Widget_param, 10);
  createGTKtextbox(hbox1, text_entry_Widget_value, 30);
  createGTKtextbox(hbox1, text_entry_Widget_fit, 4);
  createGTKtextbox(hbox1, text_entry_Widget_diff, 9);
  createGTKtextbox(hbox1, text_entry_Widget_err, 19);
  gtk_iteration_blocked--;
}

// interval is is microseconds
void set_alarm_interval(long interval)
{
  struct itimerval it_val;
  it_val.it_value.tv_sec = 0;
  it_val.it_value.tv_usec = 0;
  it_val.it_interval = it_val.it_value;
  if(setitimer(ITIMER_REAL, &it_val, NULL) == -1) {
    fprintf(stderr, "Cannot set alarm interval.\n");
    perror("error calling setitimer");
    exit(0);
  }
}


void alarm_wakeup(int i)
{
  alarm_active++;

  gtk_iteration_blocked++;
  fprintf(stderr, "Wakeup: block=%d!!!\n", gtk_iteration_blocked-1);

  int weshouldsetalarm = 1;
#ifdef HAVEPGPLOT2 // If in pgplot2 mode, let pgplot deal with setting the alarm function
  if(enablePGPLOT2 == 0) {  // If not in local mode, let pgplot2 know that our timer should be called as well.
    weshouldsetalarm = 0;
  }
#endif
  //  weshouldsetalarm = 1;
  //  int remove_line_above = 1;
  if(weshouldsetalarm) {
    fprintf(stderr, "tempo3: disabling timer until completion\n");
    set_alarm_interval(0);
  }

  int allow_gtk_proccessing = 1;

#ifndef TEMPO3_EnableGTK
  allow_gtk_proccessing = 0;          // If gtk is not enabled during runtime, there is nothing to do
#endif

#ifdef TEMPO3_EnableGTK
  if(enableGTK_flag == 0) {
    allow_gtk_proccessing = 0;        // No gtk commands are executed, so there is nothing to do, and finding out if there are pending gtk events might be dodgy
  }
#endif

  while(allow_gtk_proccessing == 1) {   // Keep processing gtk command while there are no things that prevents us from continuing
    if(gtk_iteration_blocked != 1) {
      allow_gtk_proccessing = 0;        // A function which is making gtk calls is interrupted, so possibly not healthy to continue
      break;
    }
    fprintf(stderr, "Wakeup pos 1!!!\n");

#ifdef HAVEPGPLOT2 
    if(pgplot2_gtk_local_mode_callback_active()) {    // Possibly it is not healthy to process gtk event halfway an update
      allow_gtk_proccessing = 0;
      break;
    }
#endif
    fprintf(stderr, "Wakeup pos 2!!!\n");

#ifdef TEMPO3_EnableGTK
    if(gtk_events_pending() == FALSE) {
      fprintf(stderr, "No pending gtk events\n");
      allow_gtk_proccessing = 0;
      break;
    }
#endif
    fprintf(stderr, "Wakeup pos 3!!!\n");

#ifdef TEMPO3_EnableGTK
    GdkEvent *event;
    event = gtk_get_current_event();
    if(event != NULL) {
      allow_gtk_proccessing = 0;
      fprintf(stderr, "Preventing processing gtk events, as event is active\n");
      gdk_event_free(event);
      break;
    }
#endif
    fprintf(stderr, "Wakeup pos 4!!!\n");

#ifdef TEMPO3_EnableGTK
    if(gtk_iteration_blocked == 1) { // Final check for gtk activity
#ifdef HAVEPGPLOT2 
      if(pgplot2_gtk_local_mode_callback_active() == 0) {
#endif
	fprintf(stderr, "Process pending gtk events block=%d\n", gtk_iteration_blocked);
	gtk_main_iteration_do(FALSE);   // Do not block if no events are pending. Sometimes it seems possible that although gtk_events_pending() claims there are pending iterations, gtk_main_iteration() is stuck.
	fprintf(stderr, "Process pending gtk events done\n");
#ifdef HAVEPGPLOT2 
      }
#endif
    }
#endif

    fprintf(stderr, "Wakeup pos 5!!!\n");
  }

  gtk_iteration_blocked--;

  if(weshouldsetalarm) {
    fprintf(stderr, "tempo3 restoring timer\n");
    set_alarm_interval(100000);
  }

  alarm_active = 0;
  fprintf(stderr, "Wakeup finished: block=%d!!!\n", gtk_iteration_blocked);
}

// Closing the gtk window of a plotting device is requested
gboolean delete_event_window(GtkWidget *widget, GdkEvent *event, gpointer user_data)
{
  gtk_iteration_blocked++;
  //  gtk_window_def *window = user_data;

#ifdef HAVEPGPLOT2 
  // Window gets destroyed, which means the program should terminate
  // and not wait for a never occuring keypress.
  if(enablePGPLOT2) {
    user_requested_quit = 1;
  }
#endif

  gtk_iteration_blocked--;
  return FALSE; // Propagate event further -> Window gets destroyed.
}

void initGTKwindow(int argc, char **argv)
{
  GtkWidget *window;
  GtkWidget *scrollable_window;
  GtkWidget *vbox1;
  int i;

  /* This is called in all GTK applications. Arguments are parsed
     from the command line and are returned to the application. */
  gtk_init(&argc, &argv);

  GtkCssProvider *css_provider;
  GdkDisplay *display;
  GdkScreen *screen;
  css_provider = gtk_css_provider_new();  // Set-up the css style file interface 
  display = gdk_display_get_default();
  screen = gdk_display_get_default_screen(display);
  gtk_style_context_add_provider_for_screen(screen, GTK_STYLE_PROVIDER(css_provider), GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
  GError *gdk_error;
  gdk_error = 0;
  //  gtk_css_provider_load_from_path(css_provider, "/local/scratch/wltvrede/puma1soft/trunk/src/patrickSoft/prog/tempo3.css", &gdk_error);
  //  char *css_text1 = "#text_entry {\n    font: 15px Sans;\n    color: #00ff00;\n    background-color: #000000;\n}\n";
  //  char *css_text2 = "#text_entry {\n    color: #00ff00;\n    background-color: #000000;\n}\n";
  //  char *css_text = css_text1;
  char *css_text = "#text_entry {\n    font: 15px Sans;\n    color: #00ff00;\n    background-color: #000000;\n}\n";
#ifdef HAVEPGPLOT2 
  // Have a smaller font to improve window layout
  if(enablePGPLOT2) {
    css_text = "#text_entry {\n    color: #00ff00;\n    background-color: #000000;\n}\n";
  }
#endif
  gtk_css_provider_load_from_data(css_provider, css_text, strlen(css_text), &gdk_error);
  if(gdk_error != 0) {
    fprintf(stderr, "Loading css failed with message '%s'\n", gdk_error->message);
    exit(0);
  }

  /* Create window */
  window = gtk_window_new(GTK_WINDOW_TOPLEVEL); 
  gtk_window_set_title(GTK_WINDOW (window), "Tempo3 parfile");
  //  gtk_widget_set_size_request(window, 730, 500);

  /* Amount of space around buttons etc */
  //  gtk_container_set_border_width(GTK_CONTAINER(window), 10);
  gtk_widget_show(window);  

  GtkWidget *window_grid;
  window_grid = gtk_grid_new();
  gtk_container_add(GTK_CONTAINER(window), window_grid); 

  /* Here we just set a handler for delete_event that immediately exits GTK.
     Ignore for now. Could sent quit signal to program
     g_signal_connect (G_OBJECT (window), "delete_event", G_CALLBACK (delete_event), NULL);
  */

  g_signal_connect(G_OBJECT(window), "delete-event", G_CALLBACK(delete_event_window), NULL);


  /* Make a box with a vertical scrollbar */
  scrollable_window = gtk_scrolled_window_new(NULL, NULL);
  gtk_widget_set_hexpand(scrollable_window, TRUE);
  gtk_widget_set_vexpand(scrollable_window, TRUE);
  gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollable_window), GTK_POLICY_AUTOMATIC, GTK_POLICY_ALWAYS);
  gtk_widget_show(scrollable_window);
  //  gtk_container_add(GTK_CONTAINER (window), scrollable_window); 
  gtk_grid_attach(GTK_GRID(window_grid), scrollable_window, 0, 0, 1, 1);
  gtk_widget_set_size_request(scrollable_window, 730, 500);

#ifdef HAVEPGPLOT2 
  if(enablePGPLOT2) {
    GtkWidget *gtk_drawing_area = gtk_drawing_area_new();
    gtk_grid_attach(GTK_GRID(window_grid), gtk_drawing_area, 1, 0, 1, 1);
    gtk_widget_set_size_request(gtk_drawing_area, 800, 600);
    GtkWidget *gtk_menuBar = gtk_menu_bar_new();
    if(pgplot2_gtk_exclusively_local() == 0) {
      printerror(0, "tempo3: Cannot set pgplot2 in exclusive local mode.\n");
      exit(0);
    }
    deviceID_pgplot2 = pgplot2_gtk_open_device_in_local_mode(gtk_drawing_area, NULL, NULL, NULL, 800, 600, 0, -1, -1);
    pgplot2_gtk_local_add_pgplot_menu_options(gtk_menuBar);
    gtk_grid_attach(GTK_GRID(window_grid), gtk_menuBar, 0, -1, 2, 1);
  }
#endif




  GtkWidget *hbox1;
  hbox1 = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
  gtk_container_add(GTK_CONTAINER(scrollable_window), hbox1);

  // Needed to avoid the lines of parfile to expand to fill window if resized
  GtkWidget *emptylabel;
  emptylabel = gtk_label_new(NULL);
  gtk_widget_set_hexpand(emptylabel, TRUE);
  //  gtk_grid_attach(GTK_GRID(window_grid), scrollable_window, 1, 0, 1, 1);

  /* Make a box in which each next widget will be placed in the vertical direction */
  // Should use grid instead, but this is the quick gtk2 -> gtk3 fix
  //  vbox1 = gtk_vbox_new (FALSE, 0);
  vbox1 = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  //  gtk_scrolled_window_add_with_viewport (GTK_SCROLLED_WINDOW (scrollable_window), vbox1);
  //gtk_container_add(GTK_CONTAINER(scrollable_window), vbox1);
  gtk_container_add(GTK_CONTAINER(hbox1), vbox1);
  gtk_widget_show(vbox1);

  /* Add parfile lines to the vertical box*/
  for(i = 0; i < nrActiveGTKLines; i++) {
    createGTKline(vbox1, &(gtkLine_param[i]), &(gtkLine_value[i]), &(gtkLine_fit[i]), &(gtkLine_diff[i]), &(gtkLine_err[i]), i);
  }

  gtk_widget_show_all(window);

  /* Process all pending gtk operations (window will be drawn on the screen at this point) */
  while(gtk_events_pending()) {
    gtk_main_iteration();
  }

  /* Rest in gtk_main and wait for user input.
          gtk_main ();
  */

}

#endif

/* A global variable used by random number generation code */
gsl_rng *rand_num_gen;




int main(int argc, char **argv)
{
  float xpos, ypos, deltax, deltax_best, seconds, *ygraph, *ygraph2, previous_xstart, previous_xend;
  float cur_x1, cur_y1, cur_x2, cur_y2, markersize, stridefit_nrdays, stridefit_minstride, stridefit_mintoarange;
  char ch, filename[MaxFilenameLength], txt[1000], macroname[MaxFilenameLength], tmpfilename[MaxFilenameLength], title[1000], stridefit_scriptfilename[MaxFilenameLength], orbit_dump_prefix[MaxFilenameLength];
  // , tmpfilename2[MaxFilenameLength], tmpfilename4[MaxFilenameLength], tmpfilename3[MaxFilenameLength]
  int i, j, ibest, param, key, glitchnr, nrwaves_old, plottype, showdata, finderrors, ok, autoquiet, numthreads;
  int do_show_params; // Flag to indicate that the table with parameters needs to be displayed/updated
  int doauto, automode, showWaves, macro_open, ctrlspace, nrSuperFit, ppgplot, stridefit_minnroftoas, stridefit_includeF2, stridefit_includeF3, deviceDefined, nosort, randstart, randstart_itt, randstart_ittmax, randstart_nrparams, readhex, writehex, dumpcovar, hierarchical, fit_algorithm;
  double randstart_sigma;
  long double *xstart, *dx, *xfit, rms, ftol = 1e-6, distance, distanceerror;
  int year, month, day, hours, minutes, showfiducial, showpepoch, loadtimfile_flag;
  int nfunk, runtempo2, plottitle, deviceID, dofake, do_perturb_parfile, amoeba_algorithm, timfile_shortformat;
  residual_def residuals, residuals_fiducial;
  toa_data_def toas, toas_fiducial;
  //  tempo3_wraps_def wraps;
  FILE *fin, *fmacro;
  //  char unused_parfile_lines[TEMPO3_MaxNrParLines][TEMPO3_MaxNrParLineLength]
  char highlightedFlag[1000];
  //  int nrunused_parfile_lines;
  //  int *fixed, *paramset, *paramlinked;
  double *mrq_alpha, *mrq_alpha_inv; // *mrq_covar, 
  long max_nr_fitted_params;  // Set to -1 if all defined parameters are allowed to be fitted for. If this number is set to be smaller than the actual numbers of parameters that can in principle be defined, memory usage could reduce drastically since in principle a nrparams * nrparams matrix should be defined for error calculations. 
  long double *mrq_deriv, perturb_fraction; 
  int writeoutpieces, superfit, dumpphases, nummericalApproxNuNudot, force_colsite, baryssbdumpfile; // nrparamsused, 
  long idnum;
  tempo3stack *stackptr;
  pgplot_options_definition pgplot_options;
  const gsl_rng_type *rand_num_gen_type;
  unsigned long *indx_input_toas;
  verbose_definition verbose_state;
  max_nr_fitted_params = -1;
  cleanVerboseState(&verbose_state);
  verbose_state.nocounters = 1;

#ifdef EnablePSRSALSAdevelop
#else
  printwarning(0, "WARNING: Code is compiled with a possibly old versions of library functions, so functionality is not guaranteed.");
#endif

  // Deal with command-line options first that change the number of recognized parameters
  for(i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-maxglitchnr") == 0 || strcmp(argv[i], "-maxnrglitches") == 0) {
      if(sscanf(argv[i+1], "%d", &TEMPO3_MaxNrGlitches) != 1) {   // Note: This variable should not change anymore after this point!!!!
	printerror(0, "ERROR: Cannot parse \"%s\" as an integer.", argv[i+1]);
	return 0;
      }
      if(TEMPO3_MaxNrGlitches < 1) {
	printerror(0, "ERROR: At least one glitch is required to be supported.");
	return 0;
      }
    }else if(strcmp(argv[i], "-maxwavenr") == 0 || strcmp(argv[i], "-maxnrwaves") == 0 || strcmp(argv[i], "-maxwaves") == 0 || strcmp(argv[i], "-numwaves") == 0 || strcmp(argv[i], "-nrwaves") == 0) {
      if(sscanf(argv[i+1], "%d", &TEMPO3_MaxNrWaves) != 1) {   // Note: This variable should not change anymore after this point!!!!
	printerror(0, "ERROR: Cannot parse \"%s\" as an integer.", argv[i+1]);
	return 0;
      }
      if(TEMPO3_MaxNrWaves < 1) {
	printerror(0, "ERROR: At least one set of wave parameters is required to be supported.");
	return 0;
      }
    }else if(strcmp(argv[i], "-maxcompanionnr") == 0 || strcmp(argv[i], "-maxnrcompanions") == 0) {
      if(sscanf(argv[i+1], "%d", &TEMPO3_MaxNrCompanions) != 1) {   // Note: This variable should not change anymore after this point!!!!
	printerror(0, "ERROR: Cannot parse \"%s\" as an integer.", argv[i+1]);
	return 0;
      }
      if(TEMPO3_MaxNrCompanions < 1) {
	printerror(0, "ERROR: At least one set of binary parameters is required to be supported.");
	return 0;
      }
    }else if(strcmp(argv[i], "-maxjumpnr") == 0 || strcmp(argv[i], "-maxnrjumpss") == 0) {
      if(sscanf(argv[i+1], "%d", &TEMPO3_MaxNrJumps) != 1) {   // Note: This variable should not change anymore after this point!!!!
	printerror(0, "ERROR: Cannot parse \"%s\" as an integer.", argv[i+1]);
	return 0;
      }
      if(TEMPO3_MaxNrJumps < 1) {
	printerror(0, "ERROR: At least one jump parameters is required to be supported.");
	return 0;
      }
    }else if(strcasecmp(argv[i], "-maxnfit") == 0 || strcasecmp(argv[i], "-maxnrfit") == 0) {
      if(sscanf(argv[i+1], "%ld", &max_nr_fitted_params) != 1) {   // Note: This variable should not change anymore after this point!!!!
	printerror(0, "ERROR: Cannot parse \"%s\" as an integer.", argv[i+1]);
	return 0;
      }
      if(max_nr_fitted_params < 2) {
	printerror(0, "ERROR: At least two fit paramters should be allowed to be fitted simultaneously.");
	return 0;
      }
    }else if(strcmp(argv[i], "-v") == 0) {
      verbose_state.verbose = 1;
    }else if(strcmp(argv[i], "-debug") == 0) {
      verbose_state.debug = 1;
    }
  }

  initialise_tempo3_lib(verbose_state);

  ephemeris_def parfile, parfile_old;
  if(initialise_ephemeris(&parfile, 1, verbose_state) == 0) {
    printerror(verbose_state.debug, "ERROR tempo3: Initialising ephemeris failed.");
    return 0;
  }
  if(initialise_ephemeris(&parfile_old, 0, verbose_state) == 0) {
    printerror(verbose_state.debug, "ERROR tempo3: Initialising ephemeris failed.");
    return 0;
  }
  
#ifdef TEMPO3_EnableGTK
  if(allocate_ephemeris_par_only(&gtk_parfile, 0, verbose_state) == 0) {
    printerror(verbose_state.debug, "ERROR tempo3: Initialising ephemeris failed.");
    return 0;
  }
#endif
  //  allocate_ephemeris_par_only(&parfile, verbose_state.debug);
  if(allocate_ephemeris_par_only(&(additional_funk_info.parfile_funk), 0, verbose_state) == 0) {
    printerror(verbose_state.debug, "ERROR tempo3: Initialising ephemeris failed.");
    return 0;
  }



  /* Allocate some memory with malloc, rather than take it from the stack memory. */

  xstart = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
  dx = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
  xfit = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
  //  dplus = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
  //  dmin = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
  mrq_deriv = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
  //  fixed = malloc(TEMPO3_PARAM_TOTAL*sizeof(int));
  //  paramset = malloc(TEMPO3_PARAM_TOTAL*sizeof(int));
  //  paramlinked = malloc((1+2*TEMPO3_PARAM_TOTAL)*sizeof(int));
  additional_funk_info.best_xfit = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
  stackptr = malloc(sizeof(tempo3stack));
  // fixed == NULL || paramset == NULL || paramlinked == NULL || 
  if(xstart == NULL || dx == NULL || xfit == NULL || mrq_deriv == NULL || additional_funk_info.best_xfit == NULL || stackptr == NULL) {
    printerror(verbose_state.debug, "ERROR tempo3: Cannot allocate memory\n");
    return 0;
  }


  stackptr->nrinstack = 0;
  for(i = 0; i < MaxNrStack; i++) {
    if(allocate_ephemeris_par_only(&(stackptr->parfiles[i]), 0, verbose_state) == 0) {
      printerror(verbose_state.debug, "ERROR tempo3: Initialising ephemeris failed.");
      return 0;
    }
    stackptr->fixed[i] = malloc(TEMPO3_PARAM_TOTAL*sizeof(int));
    stackptr->paramset[i] = malloc(TEMPO3_PARAM_TOTAL*sizeof(int));
    if(stackptr->fixed[i] == NULL || stackptr->paramset[i] == NULL) {
      printerror(verbose_state.debug, "ERROR tempo3: Cannot allocate memory\n");
      return 0;
    }
  }

  /*
  long double mjd = 54786.999999;

  mjd2date_old(mjd, &year, &month, &day, &hours, &minutes, &seconds);
  mjd2dateString(mjd, txt, 5, 1, " ");
  printf("%s\n", txt);
  mjd2dateString(mjd, txt, 0, 1, " ");
  printf("%s\n", txt);
  printf("%d-%02d-%02d %02d:%02d:%f\n", year, month, day, hours, minutes, seconds);
  return 0;
  */

  size_t maxallowednrparameters;
  maxallowednrparameters = floor(sqrt((double)SIZE_MAX/sizeof(double)));
  if(TEMPO3_PARAM_TOTAL > maxallowednrparameters) {
    printerror(verbose_state.debug, "ERROR tempo3: There are %ld parameters defined which exceed the maximum possible of %ld before array overflows will occur.\n", TEMPO3_PARAM_TOTAL, maxallowednrparameters);
    return 0;
  }
  /* Some matrices required by error calculation method */
  //  mrq_covar = malloc(TEMPO3_PARAM_TOTAL*TEMPO3_PARAM_TOTAL*sizeof(double));
  size_t allocsize;
  if(max_nr_fitted_params <= 0) 
    allocsize = (size_t)TEMPO3_PARAM_TOTAL*(size_t)TEMPO3_PARAM_TOTAL*(size_t)sizeof(double);
  else
    allocsize = (size_t)max_nr_fitted_params*(size_t)max_nr_fitted_params*(size_t)sizeof(double);
  mrq_alpha = malloc(allocsize);
  mrq_alpha_inv = malloc(allocsize);
  if(mrq_alpha == NULL || mrq_alpha_inv == NULL) { // mrq_covar == NULL || 
    printerror(verbose_state.debug, "ERROR tempo3: Cannot allocate memory for %ld * %ld * %ld = %ld bytes. Maybe using the -maxnfit option would be useful.\n", TEMPO3_PARAM_TOTAL, TEMPO3_PARAM_TOTAL, sizeof(double), allocsize);
    return 0;
  }

  /* By default call tempo2 to make the bary list */
  runtempo2 = 1;
  writeoutpieces = 0;
  /* Clean some variables */
  //  wraps.nrwraps = 0;
  //  paramlinked[0] = 0;  // Means: no parfile parameters are linked to VALUE parameters
  //  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
  //    fixed[i] = 1;
  //    paramset[i] = 0;
  //  }
  //  fixed[TEMPO3_PARAM_PHASE0] = 0;
  //  paramset[TEMPO3_PARAM_PHASE0] = 1;         // Could leave this out, and it will not show up in the output
  additional_funk_info.superverbose = -1;
  //  nrunused_parfile_lines = 0;
  plottype = 0;
  showdata = 0;
  doauto = 0;
  plottitle = 1;
  showWaves = 0;
  dofake = 0;
  do_perturb_parfile = 0;
  showfiducial = 0;
  showpepoch = 0;
  highlightedFlag[0] = 1;  highlightedFlag[1] = 2;  highlightedFlag[2] = 3;  highlightedFlag[3] = 0;   /* Some random non-ascii text*/
  superfit = 0;
  dumpphases = 0;
  macroname[0] = 0;
  macro_open = 0;
  nummericalApproxNuNudot = 0;
  pgplot_clear_options(&pgplot_options);
  xpos = 0;
  ypos = 0;
  ch = 0;
  deviceDefined = 0;
  nrSuperFit = 10;   // By default, ctrl-space fits 10 times
  repeatable = 0;    // If set, lines that contain random numbers etc are not produced
  nosort = 0;    // By default, sort the toa's in mjd order
#ifdef TEMPO3_EnableGTK
  enableGTK_flag = 1;  /* Enable GTK stuff by default if compiled */
#ifdef HAVEPGPLOT2 
  enablePGPLOT2 = 0;
#endif
  gtk_fixed = parfile.fixed;   /* button callback function needs to be able to find fixed array */
  gtk_parfile2 = &(parfile_old.par);
  gtk_toas = &toas;
  gtk_wraps = parfile.wraps;
  gtk_paramset = parfile.paramset;
  gtk_dplus = parfile.par.dplus;
  gtk_dmin = parfile.par.dmin;
  gtk_paramlinked = parfile.paramlinked;
#endif
  previous_xstart = previous_xend = 0;
  ctrlspace = 0;
  title[0] = 0;
  amoeba_algorithm = 0;
  ppgplot = 0;
  markersize = 1;
  stridefit_nrdays = 150;
  stridefit_mintoarange = 0;
  stridefit_minstride = 0;
  stridefit_minnroftoas = 4;
  stridefit_includeF2 = 0;
  stridefit_includeF3 = 0;
  stridefit_scriptfilename[0] = 0;
  timfile_shortformat = 0;
  autoquiet = 0;
  distance = -1;  // distance in kpc, by default neg = disable
  numthreads = 1;
  force_colsite = 0;
  baryssbdumpfile = 0;
  randstart = 0;
  randstart_itt = 0;
  randstart_ittmax = 2;
  randstart_nrparams = 3;
  randstart_sigma = 3.0;
  readhex = 0;
  writehex = 0;
  dumpcovar = 0;
  hierarchical = 0;
  fit_algorithm = 0;

  for(i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-listparameters") == 0 || strcmp(argv[i], "-listparams") == 0 || strcmp(argv[i], "-paramlist") == 0 || strcmp(argv[i], "-paramslist") == 0) {
      printf("The supported parameters are (%d in total taking %ld bytes of memory):\n\n", TEMPO3_PARAM_TOTAL, sizeof(tempo3_parameters_def) + sizeof(long double)*TEMPO3_PARAM_TOTAL + sizeof(int)*TEMPO3_MaxNrCompanions + sizeof(char *)*TEMPO3_MaxNrJumps + TEMPO3_MaxNrJumps*TEMPO3_MaxFlagsLength);
      for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	if(i == TEMPO3_PARAM_VALUE+1 || i == TEMPO3_PARAM_F0+2 || i == TEMPO3_PARAM_GLEP+1 || i == TEMPO3_PARAM_GLF0+1 || i == TEMPO3_PARAM_WAVESIN+1 || i == TEMPO3_PARAM_DM+2) { //  || i == TEMPO3_PARAM_GLF1D+1
	  printf("...\n");
	}else if( (i > TEMPO3_PARAM_VALUE+1 && i < TEMPO3_PARAM_VALUE+TEMPO3_MaxNrValues-1)
		  || (i > TEMPO3_PARAM_F0+2 && i < TEMPO3_PARAM_F0+TEMPO3_MaxNrFderivatives)
		  || (i > TEMPO3_PARAM_DM+2 && i < TEMPO3_PARAM_DM+TEMPO3_MaxNrDMderivatives)
		  || (i > TEMPO3_PARAM_GLEP+1 && i < TEMPO3_PARAM_GLEP+TEMPO3_MaxNrGlitches-1)
		  || (i >= TEMPO3_PARAM_GLEPMIN+1 && i <= TEMPO3_PARAM_GLEPMIN+TEMPO3_MaxNrGlitches-1)
		  || (i >= TEMPO3_PARAM_GLEPMAX+1 && i <= TEMPO3_PARAM_GLEPMAX+TEMPO3_MaxNrGlitches-1)
		  || (i >= TEMPO3_PARAM_GLPH+1 && i <= TEMPO3_PARAM_GLPH+TEMPO3_MaxNrGlitches-1)
		  || (i > TEMPO3_PARAM_GLF0+1 && i < TEMPO3_PARAM_GLF0+TEMPO3_MaxNrGlitches-1+TEMPO3_MaxNrGlitchFderivatives*TEMPO3_MaxNrGlitches)
		  || (i >= TEMPO3_PARAM_GLNBRAKE+1 && i <= TEMPO3_PARAM_GLNBRAKE+TEMPO3_MaxNrGlitches-1)
		  || (i >= TEMPO3_PARAM_GLTD+1 && i <= TEMPO3_PARAM_GLTD+TEMPO3_MaxNrGlitches-1)
		  || (i >= TEMPO3_PARAM_GLF0D+1 && i <= TEMPO3_PARAM_GLF0D+TEMPO3_MaxNrGlitches-1)
		  //		    || (i > TEMPO3_PARAM_GLF1D+1 && i < TEMPO3_PARAM_GLF1D+TEMPO3_MaxNrGlitches-1)
		  || (i >= TEMPO3_PARAM_JUMPS+1 && i < TEMPO3_PARAM_JUMPS+TEMPO3_MaxNrJumps)  // Ignore all but first jump parameter
		|| (i >= TEMPO3_PARAM_WAVECOS && i < TEMPO3_PARAM_WAVECOS+TEMPO3_MaxNrWaves) // Ignore all cosine terms
		  || (i > TEMPO3_PARAM_WAVESIN+1 && i < TEMPO3_PARAM_WAVESIN+TEMPO3_MaxNrWaves-1)
		  || (i >= TEMPO3_PARAM_PB+1 && i < TEMPO3_PARAM_PB+TEMPO3_MaxNrCompanions)
		  || (i >= TEMPO3_PARAM_PBDOT+1 && i < TEMPO3_PARAM_PBDOT+TEMPO3_MaxNrCompanions)
		  || (i >= TEMPO3_PARAM_A1+1 && i < TEMPO3_PARAM_A1+TEMPO3_MaxNrCompanions)
		  || (i >= TEMPO3_PARAM_A1DOT+1 && i < TEMPO3_PARAM_A1DOT+TEMPO3_MaxNrCompanions)
		  || (i >= TEMPO3_PARAM_OM+1 && i < TEMPO3_PARAM_OM+TEMPO3_MaxNrCompanions)
		  || (i >= TEMPO3_PARAM_OMDOT+1 && i < TEMPO3_PARAM_OMDOT+TEMPO3_MaxNrCompanions)
		  || (i >= TEMPO3_PARAM_ECC+1 && i < TEMPO3_PARAM_ECC+TEMPO3_MaxNrCompanions)
		  || (i >= TEMPO3_PARAM_GAMMA+1 && i < TEMPO3_PARAM_GAMMA+TEMPO3_MaxNrCompanions)
		  ) {
	}else {
	  if(i == TEMPO3_PARAM_JUMPS) {
	    printf("\nThe JUMP parameters can be used assign offsets to certain groups of toa's. Up to %d jumps can be defined:\n", TEMPO3_MaxNrJumps);
	  }else if(i == TEMPO3_PARAM_GLEP) {
	    printf("\nUp to %d glitches can be defined (see GLEP parameter for the how to specify glitch numbers):\n", TEMPO3_MaxNrGlitches);
	  }else if(i == TEMPO3_PARAM_WAVEEPOCH) {
	    printf("\nThe WAVE parameters can be used to whiten the data. %d sin/cosine parirs can be defined for subsequent harmonics:\n", TEMPO3_MaxNrWaves);
	  }else if(i == TEMPO3_PARAM_BINARY) {
	    printf("\nFor the binary parameters up to %d companions can be defined (see T0 parameter for the how to specify companion numbers):\n", TEMPO3_MaxNrCompanions);
	  }else if(i == TEMPO3_PARAM_TZRMJD) {
	    printf("\nOther parameters:\n");
	  }
	  
	  printf("%-15s", tempo3_parfile_descr.identifiers[i][0]);
	  if(tempo3_parfile_descr.identifiers[i][1] != NULL) {
	    printf(" or %s", tempo3_parfile_descr.identifiers[i][1]);
	  }
	  printf("%-15s", tempo3_parfile_descr.units[i]);
	  printf("%-10s", tempo3_parfile_descr.description[i]);
	  printf("\n");
	}
      }

      free(mrq_alpha); 
      free(mrq_alpha_inv); 
      free(xstart);
      free(dx);
      free(xfit);
      //      free(dmin);
      free(mrq_deriv);
      //      free(fixed);
      //      free(paramset);
      //      free(paramlinked);
      free(additional_funk_info.best_xfit);
      for(i = 0; i < MaxNrStack; i++) {
	free_ephemeris_par_only(&(stackptr->parfiles[i]));
	free(stackptr->fixed[i]);
	free(stackptr->paramset[i]);
      }
      free(stackptr);
#ifdef TEMPO3_EnableGTK
      free_ephemeris_par_only(&gtk_parfile);
#endif
      //      free_ephemeris_par_only(&parfile_old);
      free_ephemeris_par_only(&(additional_funk_info.parfile_funk));
      free_ephemeris(&parfile, verbose_state);
      free_ephemeris(&parfile_old, verbose_state);
      if(cleanup_tempo3_lib(verbose_state) == 0) {
	printerror(verbose_state.debug, "ERROR tempo3: Releasing tempo3 library related memory failed.");
	return 0;
      }

      return 0;
    }
  }

  if(argc < 3) {
    printf("Usage: tempo3 [options] parfile toafile\n\n");
    printf("Options related to fitting:\n");
    printf("  -auto              Fit automatically for all parameters, except DM.\n");
    printf("  -autoN             Fit automatically for all parameters, except DM and NBREAK.\n");
    printf("  -num               Use a numerical derivatives of the par file for nu/nudot,\n");
    printf("                     rather than the analytic derivatives of the par file, to\n");
    printf("                     determine nu and nudot as function of time.\n");
    printf("  -perturb           Perturb parameters with this fraction (except F0 and F1).\n");
    printf("  -randstart \"n m s\" Every nth iteration of the fitting process starts at a\n");
    printf("                     somewhat random start position, rather than the best-fit\n");
    printf("                     solution so far, by perturbing m parameters by an random\n");
    printf("                     uniform offset up to s times the choosen initial stepsize.\n");
    printf("                     Using this option might help in reaching the global best\n");
    printf("                     solution, although might not work well in practice.\n");
    printf("  -threads N         Specify the number of threads to be used\n");
    printf("  -tol               Set fit tolerance, default is %Le.\n", ftol);
    printf("\nOptions related to what appears in terminal:\n");
    printf("  -autoquiet       Used with -auto or -autoN: output is minimal, which can be\n");
    printf("                   useful in scripts.\n");
    printf("  -debug           Enable debugging messages\n");
    printf("  -dist \"X dX\"     The distance to the pulsar is X kpc, enabling distance\n");
    printf("                     dependent derived parameters (i.e. Shklovskii effect).\n");
    printf("  -dumpcovar       Show the covariance matrix.\n");
    printf("  -dumpcovarnorm   Show the covariance matrix (with the units of the parameters\n");
    printf("                   normalised such that diagonal is 1).\n");
    printf("  -hierarchical    The orbits are assumed to be hierarchical, i.e. a\n");
    printf("                   companions orbit is a pure binary orbit where the\n");
    printf("                   other mass in the centre of mass of the inner orbits\n");
    printf("                   The companions should be defined in order of increasing\n");
    printf("                   orbital size. Only implemented for triple systems.\n");
    printf("  -v               Enable verbose\n");
    printf("\nOptions related to input/output files:\n");
    printf("  -bary            TOAs are provided as a list of barycentred TOAs, errorbars in\n");
    printf("                   microsec, freqs and labels. Optional: Flags after a | symbol.\n");
    printf("                   Using this option avoids tempo2 system calls.\n");
    printf("  -baryssb         Like -bary, but the TOA error is followed by the barycenter\n");
    printf("                   coordinates ssb1, ssb2 and ssb3.\n");
    printf("  -baryssbdump XXX After calling tempo2, the information is copied in file XXX\n");
    printf("                   for later use with the -baryssb option.\n");
    printf("  -dump            Dump the residuals (in phase, not time).\n");
    printf("  -macro FILE      Execute keypressed from this macro (a ^ prefix indicates the\n");
    printf("                   use of the ctrl key, the word RETURN is the return key).\n");
    printf("  -nosort          Do not sort the input toa's in mjd order.\n");
    printf("  -readhex         Parameter file is read in hex format.\n");
    printf("  -sitecolumn i    Assume site codes are in specified column i in toafile\n");
    printf("  -writehex        Parameter file is written in hex format.\n");
    printf("\nOptions related to the available parameters:\n");
    printf("  -listparameters  List all the supported parameters in a parameter file\n");
    printf("  -maxcompanionnr  Change the maximum number of binary companions to be\n");
    printf("                   supported. Default is %d and avoid using 0.\n", TEMPO3_MaxNrCompanions);
    printf("  -maxglitchnr     Change the maximum number of glitches supported.\n");
    printf("                   Default is %d\n", TEMPO3_MaxNrGlitches);
    printf("  -maxjumpnr       Change the maximum number of JUMP parameters supported.\n");
    printf("                   Default is %d\n", TEMPO3_MaxNrJumps);
    printf("  -maxnfit         Change the maximum number of parameters that can fit simul-\n");
    printf("                   taneously. Reducing this number could drastically limit the\n");
    printf("                   memory usage. By default enough memory is allocated to allow\n");
    printf("                   simultanous fitting of all parameters that in priciple can be\n");
    printf("                   defined.\n");
    printf("  -maxwavenr       Change the maximum number of WAVE parameters supported.\n");
    printf("                   Default is %d\n", TEMPO3_MaxNrWaves);
    printf("\nOptions related to graphics:\n");
    printf("  -device          Specify pgplot device to use (same as -dev).\n");
#ifdef TEMPO3_EnableGTK
    printf("  -nogtk           Disable the gtk parfile window.\n");
#ifdef HAVEPGPLOT2 
    printf("  -pgplot2         Use the integrated pgplot/gtk gui.\n");
#endif
#endif
    printf("  -ppgplot FILE    Output pgplot commands to FILE rather than executing them.\n");
    printf("  -title           Set a title for the plot.\n");
    printf("\nOther options:\n");
    printf("  -f               Ignored (in case you type it by accident).\n");
    printf("  -fake            Generate fake toas consistent with the parfile. THE TOAFILE\n");
    printf("                   AT THE END OF COMMAND LINE WILL BE OVERWRITTEN IF IT EXISTS.\n");
    printf("  -gr ...          Ignored (in case you type it by accident).\n");
    free(mrq_alpha); 
    free(mrq_alpha_inv); 
    free(xstart);
    free(dx);
    free(xfit);
    //    free(dmin);
    free(mrq_deriv);
    //    free(fixed);
    //    free(paramset);
    //    free(paramlinked);
    free(additional_funk_info.best_xfit);
    for(i = 0; i < MaxNrStack; i++) {
      free_ephemeris_par_only(&(stackptr->parfiles[i]));
      free(stackptr->fixed[i]);
      free(stackptr->paramset[i]);
    }
    free(stackptr);
#ifdef TEMPO3_EnableGTK
    free_ephemeris_par_only(&gtk_parfile);
#endif
    //    free_ephemeris_par_only(&parfile_old);
    free_ephemeris_par_only(&(additional_funk_info.parfile_funk));
    free_ephemeris(&parfile, verbose_state);
    free_ephemeris(&parfile_old, verbose_state);
    if(cleanup_tempo3_lib(verbose_state) == 0) {
      printerror(verbose_state.debug, "ERROR tempo3: Releasing tempo3 library related memory failed.");
      return 0;
    }
    return 0;
  }      

  if(argc >= 4) {
    for(i = 1; i < argc-2; i++) {
      if(strcmp(argv[i], "-maxglitchnr") == 0 || strcmp(argv[i], "-maxnrglitches") == 0) {   // Already dealt with above
	i++;
      }else if(strcmp(argv[i], "-maxwavenr") == 0 || strcmp(argv[i], "-maxnrwaves") == 0 || strcmp(argv[i], "-maxwaves") == 0 || strcmp(argv[i], "-numwaves") == 0 || strcmp(argv[i], "-nrwaves") == 0) {   // Already dealt with above
	i++;
      }else if(strcmp(argv[i], "-maxcompanionnr") == 0 || strcmp(argv[i], "-maxnrcompanions") == 0) {   // Already dealt with above
	i++;
      }else if(strcmp(argv[i], "-maxjumpnr") == 0 || strcmp(argv[i], "-maxnrjumpss") == 0) {   // Already dealt with above
	i++;
      }else if(strcasecmp(argv[i], "-maxnfit") == 0 || strcasecmp(argv[i], "-maxnrfit") == 0) {    // Already dealt with above
	i++;
      }else if(strcmp(argv[i], "-bary") == 0) {
	runtempo2 = 0;
	timfile_shortformat = 1;
      }else if(strcmp(argv[i], "-baryssb") == 0) {
	runtempo2 = 0;
	timfile_shortformat = 0;
      }else if(strcmp(argv[i], "-baryssbdump") == 0) {
	baryssbdumpfile = i+1;
	i++;
      }else if(strcmp(argv[i], "-gr") == 0) {
	i++;
      }else if(strcmp(argv[i], "-randstart") == 0) {
	j = sscanf(argv[++i], "%d %d %lf", &randstart_ittmax, &randstart_nrparams, &randstart_sigma);
	if(j != 3) {
	  printerror(0, "ERROR tempo3: Expected 3 arguments for option '%s'\n", argv[i-1]);
	  return 0;
	}
	randstart = 1;
      }else if(strcmp(argv[i], "-readhex") == 0) {
	readhex = 1;
      }else if(strcmp(argv[i], "-writehex") == 0) {
	writehex = 1;
      }else if(strcmp(argv[i], "-dumpcovar") == 0) {
	dumpcovar = 1;
      }else if(strcmp(argv[i], "-dumpcovarnorm") == 0) {
	dumpcovar = 2;
      }else if(strcmp(argv[i], "-hierarchical") == 0) {
	hierarchical = 1;
      }else if(strcmp(argv[i], "-f") == 0) {
      }else if(strcmp(argv[i], "-v") == 0) {   // Already dealt with above
      }else if(strcmp(argv[i], "-debug") == 0) {   // Already dealt with above
      }else if(strcmp(argv[i], "-auto") == 0) {
	doauto = 1;
	automode = 1;
      }else if(strcmp(argv[i], "-autoN") == 0) {
	doauto = 1;
	automode = 2;
      }else if(strcmp(argv[i], "-autoquiet") == 0) {
	autoquiet = 1;
	repeatable = 1;  // Gets rid of some messages automatically
      }else if(strcmp(argv[i], "-fake") == 0) {
	dofake = 1;
      }else if(strcmp(argv[i], "-dump") == 0) {
	dumpphases = 1;
      }else if(strcmp(argv[i], "-num") == 0) {
	nummericalApproxNuNudot = 1;
      }else if(strcmp(argv[i], "-nosort") == 0) {
	nosort = 1;
      }else if(strcmp(argv[i], "-repeatable") == 0) {
	repeatable = 1;
#ifdef EnablePSRSALSAdevelop
      }else if(strcmp(argv[i], "-rev") == 0) {
	printf("Revision info:\n");
	SHOWREVISIONINFO_prog();
	printf("\n");
#endif
      }else if(strcmp(argv[i], "-ppgplot") == 0) {
	ppgplot = i+1;
	i++;
      }else if(strcmp(argv[i], "-dev") == 0 || strcmp(argv[i], "-device") == 0) {
	deviceDefined = i+1;
	i++;
      }else if(strcmp(argv[i], "-macro") == 0) {
	sscanf(argv[++i], "%s", macroname);
      }else if(strcmp(argv[i], "-sitecolumn") == 0) {
	sscanf(argv[++i], "%d", &force_colsite);
      }else if(strcmp(argv[i], "-dist") == 0) {
	sscanf(argv[++i], "%Lf %Lf", &distance, &distanceerror);
      }else if(strcmp(argv[i], "-perturb") == 0) {
	sscanf(argv[++i], "%Lf", &perturb_fraction);
	do_perturb_parfile = 1;
      }else if(strcmp(argv[i], "-threads") == 0 || strcmp(argv[i], "-thread") == 0) {
	sscanf(argv[++i], "%d", &numthreads);
      }else if(strcmp(argv[i], "-tol") == 0) {
	sscanf(argv[++i], "%Lf", &ftol);
      }else if(strcmp(argv[i], "-title") == 0) {
	strcpy(title, argv[++i]);
#ifdef TEMPO3_EnableGTK
      }else if(strcmp(argv[i], "-nogtk") == 0) {
	enableGTK_flag = 0;
#ifdef HAVEPGPLOT2 
      }else if(strcmp(argv[i], "-pgplot2") == 0) {
	enablePGPLOT2 = 1;
#endif
#endif
      }else {
	printerror(0, "Unknown option '%s'\n", argv[i]);
	return 0;
      }
    }
  }

#ifdef HAVEPGPLOT2 
  if(enablePGPLOT2 && deviceDefined) {
    printerror(0, "tempo3: Cannot set use the -pgplot2 and -device options together.\n");
    return 0;
  }
#endif

  if(ppgplot) {
    if(macroname[0] == 0) {
      printerror(0, "tempo3: Use -ppgplot together with -macro\n");
      return 0;
    }
#ifdef EnablePSRSALSAdevelop
    if(pgopenoutputfile(argv[ppgplot]) == 0) {
      printerror(0, "tempo3: Cannot open %s\n", argv[ppgplot]);
      return 0;
    }
#endif
  }

  gsl_rng_env_setup();    /* Set the default generators, can be influenced by environment variables */
  rand_num_gen_type = gsl_rng_default;
  rand_num_gen = gsl_rng_alloc(rand_num_gen_type);
  if(repeatable == 0)
    randomize_idnum(&idnum);
  else
    idnum = 123;
  gsl_rng_set(rand_num_gen, idnum);

  //  initialise_ephemeris_par_only(&parfile);

  if(read_ephemeris(argv[argc-2], &parfile, readhex, 0, verbose_state) == 0) { // , autoquiet
    printerror(verbose_state.debug, "ERROR tempo3: Reading ephemeris from '%s' failed.", argv[argc-2]);
    return 0;
  }

#ifdef TEMPO3_EnableGTK
  nrActiveGTKLines = 0;
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    if(parfile.paramset[i]) {
      // MORE PARAMETERS ARE NOT PRINTED, AND SHOULD HERE BE EXCLUDED. MAKE A GENERAL ACCEPT PARAMETER FOR GTK PRINTING FUNCTION. NEED TO BE ABLE TO EXTEND TABLE LATER IF FIT PARAMETERS INCREASES, AS NOW THE PAR FILE WILL NOT FIT IN TABLE IF TOO MANY PARAMETERS ARE ADDED.
      if(i < TEMPO3_PARAM_WAVECOS || i > TEMPO3_PARAM_WAVECOS + TEMPO3_MaxNrWaves)
	nrActiveGTKLines++;
    }
  }
  nrActiveGTKLines += 1;  /* Keep extra space for pulsar name */
  if(nrActiveGTKLines < 15)
    nrActiveGTKLines = 15;
  if(nrActiveGTKLines > MaxNrGTKLines)
    nrActiveGTKLines = MaxNrGTKLines;
  if(dofake != 0 || dumpphases != 0) {
    enableGTK_flag = 0;
#ifdef HAVEPGPLOT2 
    enablePGPLOT2 = 0;
#endif
  }
  // Moved to later to avoid freezing gtk window. Appears to work fine without printparams failing.
  //  if(enableGTK_flag)
  //    initGTKwindow(argc, argv);
#endif
  if(do_perturb_parfile) {
    perturb_parfile(&(parfile.par), parfile.paramset, parfile.paramlinked, perturb_fraction);
  }
  if(verbose_state.verbose && autoquiet == 0) {
    print_ephemeris_ExtraParams(stdout, &(parfile.par), parfile.paramset, parfile.paramlinked, 0, NULL, distance, distanceerror, hierarchical);
    verbose_state.indent = 2;
    if(print_ephemeris(stdout, &parfile, NULL, NULL, 1, 0, 1, 1, autoquiet, 0, verbose_state) == 0) {
      printerror(verbose_state.debug, "ERROR tempo3: Showing ephemeris failed.");
      return 0;
    }
    verbose_state.indent = 0;
  }

  copy_parfile(&(parfile_old.par), &(parfile.par));

  /* Load the fiducial point as a barycentric toa by generating a tim file and running tempo2 */
  if(parfile.paramset[TEMPO3_PARAM_TZRFRQ] && parfile.paramset[TEMPO3_PARAM_TZRMJD] && parfile.paramset[TEMPO3_PARAM_TZRSITE] && runtempo2 != 0) {
    if(verbose_state.verbose) {
      printf("Determine phase of fiducial point\n");
      fflush(stdout);
    }
    char *username;
    if(getUsername(&username, verbose_state) == 0) {
      username = malloc(8);
      if(username == NULL) {
	printerror(verbose_state.debug, "ERROR tempo3: Memory allocation error");
	return 0;
      }
      sprintf(username, "Unknown");
    }
    sprintf(tmpfilename,  "/tmp/junk_%s_%ld.tim", username, randomUnsignedInt());
    free(username);

    if(initialise_toas(&toas_fiducial, 1, 0, 0, 0, 0, verbose_state) == 0) {
      printerror(verbose_state.debug, "ERROR 42ft: Cannot initialise TOAs");
      return 0;
    }
    if(add_toa(&toas_fiducial, 0, parfile.par.parameter[TEMPO3_PARAM_TZRMJD], 1.0, parfile.par.parameter[TEMPO3_PARAM_TZRFRQ], sqrt(-1), 0, "FiducialPoint", NULL, parfile.par.tzrsite, NULL, verbose_state) == 0) {
      printerror(verbose_state.debug, "ERROR 42ft: Adding TOAs failed");
      return 0;
    }
    int oldverbose;
    oldverbose = verbose_state.verbose;
    if(verbose_state.debug == 0 || repeatable) {
      verbose_state.verbose = 0;
    }
    if(writetimfile(tmpfilename, toas_fiducial, NULL, NULL, 0, verbose_state) == 0) {
      printerror(verbose_state.debug, "ERROR 42ft: Writing TOAs failed");
      return 0;
    }
    free_toas(&toas_fiducial);
    if(loadtimfile(tmpfilename, 0, 0, argv[argc-2], parfile, &toas_fiducial, 1, &indx_input_toas, repeatable, NULL, verbose_state) == 0) {
      printerror(verbose_state.debug, "ERROR tempo3: Reading fiducial point TOA failed");
      return 0;
    }
    verbose_state.verbose = oldverbose;
    free(indx_input_toas); // Sorting not relevant with only one TOA
    remove(tmpfilename);

    residuals_fiducial.residy = malloc((toas_fiducial.nrtoas + toas_fiducial.nrtoasdeleted)*sizeof(float));
    residuals_fiducial.resid =  malloc((toas_fiducial.nrtoas + toas_fiducial.nrtoasdeleted)*sizeof(long double));
    residuals_fiducial.nturn =  malloc((toas_fiducial.nrtoas + toas_fiducial.nrtoasdeleted)*sizeof(long long));
    residuals_fiducial.nturn_fixed =  calloc((toas_fiducial.nrtoas + toas_fiducial.nrtoasdeleted), sizeof(int));
    // Already seems to be allocated
    if(toas_fiducial.mjdx != NULL) {
      free(toas_fiducial.mjdx);
      toas_fiducial.mjdx = NULL;
    }
    if(toas_fiducial.freqx != NULL) {
      free(toas_fiducial.freqx);
      toas_fiducial.freqx = NULL;
    }
    toas_fiducial.mjdx  = malloc((toas_fiducial.nrtoas + toas_fiducial.nrtoasdeleted)*sizeof(float));
    toas_fiducial.freqx  = malloc((toas_fiducial.nrtoas + toas_fiducial.nrtoasdeleted)*sizeof(float));
    if(residuals_fiducial.residy == NULL || toas_fiducial.mjdx == NULL || toas_fiducial.freqx == NULL || residuals_fiducial.nturn == NULL || residuals_fiducial.resid == NULL || residuals_fiducial.nturn_fixed == NULL) {
      printerror(verbose_state.debug, "Cannot allocate memory\n");
      return 0;
    }
    residuals_fiducial.mjdmin = toas_fiducial.mjd[0]-1; // Range just needs to be big enough to include the toa so calculations are happening with the fiducial point
    residuals_fiducial.mjdmax = toas_fiducial.mjd[0]+1;
    if(verbose_state.verbose) {
      printf("Determine phase of fiducial point done\n");
      fflush(stdout);
    }
  }else {
    toas_fiducial.nrtoas = 0;
  }

  char fakeInputFilename[1000];
  loadtimfile_flag = 1;
  if(dofake) {
    printf("\nNOTE THAT %s WILL BE OVERWRITTEN. IF THIS IS NOT A GOOD IDEA, QUIT WITH CTRL-C\n\n", argv[argc-1]);
    printf("Name of tim file to get sampling (press return if not used): ");
    fflush(stdout);
    fgets(fakeInputFilename, 1000, stdin);
    if(fakeInputFilename[0] != '\n') {
      for(i = 0; i < 1000; i++) {
	if(fakeInputFilename[i] == '\n')
	  fakeInputFilename[i] = 0;
      }
      FILE *test;
      test = fopen(fakeInputFilename, "rb");
      if(test == NULL) {
	printerror(0, "Cannot open '%s'\n", fakeInputFilename);
	return 0;
      }else {
	fclose(test);
      }
    }else {
      loadtimfile_flag = 0;
    }
  }

  if(loadtimfile_flag) {
    if(runtempo2) {
      char *timfilename;
      if(dofake) {
	timfilename = fakeInputFilename;
      }else {
	timfilename = argv[argc-1];
	printf("Reading in %s\n", argv[argc-1]);
      }
      char *baryssbdumpfilename;
      if(baryssbdumpfile > 0) {
	baryssbdumpfilename = argv[baryssbdumpfile];
      }else {
	baryssbdumpfilename = NULL;
      }
      if(loadtimfile(timfilename, 1, force_colsite, argv[argc-2], parfile, &toas, nosort, &indx_input_toas, repeatable, baryssbdumpfilename, verbose_state) == 0) {
	printerror(verbose_state.debug, "ERROR tempo3: No TOA's loaded!\n");
	return 0;
      }
      residuals.mjdmin = toas.data_mjdmin;
      residuals.mjdmax = toas.data_mjdmax;      
    }else {
      int oldverbose;
      oldverbose = verbose_state.verbose;
      if(autoquiet != 0 && verbose_state.debug == 0) {
	verbose_state.verbose = 0;
      }
      toas.nrtoas = readtimfile_custom(argv[argc-1], &toas, timfile_shortformat, verbose_state);
      //      exit(0);
      residuals.mjdmin = toas.data_mjdmin;
      residuals.mjdmax = toas.data_mjdmax;
      verbose_state.verbose = oldverbose;
    }
    if(toas.nrtoas == 0) {
      printerror(0, "Error, no TOA's loaded!\n");
      return 0;
    }
  }

  if(dofake) {
    long double startmjd, endmjd, freq, err, mjd, phase, tzphase, phase0, timejitter;
    int nrtoas, iformat, loaderrors, randommjds;
    char profileName[1000];
    char flags[2], site[2];
    long long nturn;
    datafile_definition probdist;
    float *cumdist, y;

    flags[0] = 0;
    site[0] = 0;
    loaderrors = 0;
    randommjds = 1;

    if(loadtimfile_flag == 0) {
      printf("Start MJD: ");
      fflush(stdout);
      scanf("%Lf", &startmjd);
      printf("End MJD: ");
      fflush(stdout);
      scanf("%Lf", &endmjd);
    }
    printf("Frequency [MHz]: ");
    fflush(stdout);
    scanf("%Lf", &freq);
    if(loadtimfile_flag) {
      printf("Use errorbars from tim file (yes/no): ");
      fflush(stdout);
      scanf("%s", txt);
      if(strcasecmp(txt, "yes") == 0) {
	loaderrors = 1;
	printf("Take errorbars from tim file\n");
      }else {
	printf("Don't take errorbars from tim file\n");
      }
    }
    if(loaderrors == 0) {
      printf("Errorbars [microseconds]: ");
      fflush(stdout);
      scanf("%Lf", &err);
    }
    if(loadtimfile_flag == 0) {
      printf("Number of TOAs: ");
      fflush(stdout);
      scanf("%d", &nrtoas);
      printf("Generate equally seperated TOA's (yes/no): ");
      fflush(stdout);
      scanf("%s", txt);
      if(strcasecmp(txt, "yes") == 0) {
	randommjds = 0;
	printf("Generate equally seperated TOA's\n");
      }else {
	randommjds = 1;
	printf("Generate randomized MJD's\n");
      }
    }else {
      nrtoas = toas.nrtoas;
      printf("Amount of jitter on observing dates (in days): ");
      fflush(stdout);
      scanf("%Lf", &timejitter);
    }
    printf("Filename (used as a distribution of toa's across pulse phase, empty name/nonexisting file results in use of delta function): ");
    fflush(stdout);
    fgets(profileName, 1000, stdin);
    /* Ignore return of last input. Might be system dependent if this is necessary. */
    if(profileName[0] == 10)
      fgets(profileName, 1000, stdin);
    
    if(profileName[strlen(profileName)-1] == 10)
      profileName[strlen(profileName)-1] = 0;
    iformat = -1;
    if(iformat <= 0) {
      iformat = guessPSRData_format(profileName, 0, verbose_state);
      if(iformat == -2 || iformat == -3)
	return 0;
    }
    if(isValidPSRDATA_format(iformat) == 0) {
      printwarning(0, "ERROR tempo3: Input file cannot be opened. Please check if file %s exists. Otherwise it probably is in an unknown data format. Input file is ignored and a delta function is assumed.\nContinuing\n\n", profileName);
      printwarning(0, "");
      probdist.NrBins = 0;
    }else {
      /* Open input file and read header */
      if(openPSRData(&probdist, profileName, iformat, 0, 1, 0, verbose_state) == 0) {
	printwarning(0, "Assuming a delta function for the probability function.\nContinuing");
	probdist.NrBins = 0;
      }else {
	if(probdist.NrFreqChan > 1)
	  printwarning(0, "WARNING: TAKE FIRST FREQUENCY CHANNEL");
	if(probdist.NrSubints > 1)
	  printwarning(0, "WARNING: TAKE FIRST SUBINT");

	cumdist = (float *)malloc((probdist.NrBins+1)*sizeof(float));
	/* Make cummulative distribution */
	cumdist[0] = 0;
	for(i = 0; i < probdist.NrBins; i++) {
	  cumdist[i+1] = probdist.data[i];
	  if(cumdist[i] < 0) {
	    printwarning(0, "WARNING: All data points should be positive. Setting negative values to zero.");
	    cumdist[i] = 0;
	  }
	}
	for(i = 1; i < probdist.NrBins+1; i++) {
	  cumdist[i] += cumdist[i-1];
	}
	for(i = 0; i <= probdist.NrBins; i++) {
	  cumdist[i] /= cumdist[probdist.NrBins];
	  /*	printf("%d %f\n", i, cumdist[i]); */
	}
      }
    }

  
    fin = fopen(argv[argc-1], "w");
    if(fin == NULL) {
      printerror(0, "Cannot open %s\n", argv[argc-1]);
      return 0;
    }


    ephemeris_setLinkedValues_par_version(&(parfile.par), parfile.paramlinked);
    if(toas_fiducial.nrtoas == 1) {
      long double freq;
      if(toas_fiducial.freqSSB[0] > 0) {
	freq = toas_fiducial.freqSSB[0];
      }else {
	freq = toas_fiducial.freqSite[0];
      }
      tzphase = evaluate_ephemeris(toas_fiducial.mjd[0], freq, 0, 0, 0, 0, 0, site, flags, &(parfile.par), 0, parfile.paramset, 0, 0); 
    }else {
      tzphase = 0;
    }
    fprintf(fin, "FORMAT 1\n"); 
    for(i = 0; i < nrtoas; i++) {
      /* Choose a random TOA within the specified mjd range */
      if(loadtimfile_flag == 0) {
	if(randommjds) {
	  mjd = startmjd + gsl_rng_uniform(rand_num_gen)*(endmjd-startmjd);
	}else {
	  if(nrtoas > 1)
	    mjd = startmjd + i*(endmjd-startmjd)/(long double)(nrtoas-1);
	  else
	    mjd = (endmjd+startmjd)*0.5;
	}
      }else {
	mjd = (gsl_rng_uniform(rand_num_gen)-0.5)*((1.0/parfile.par.parameter[TEMPO3_PARAM_F0])/TEMPO3_SecPerDay + timejitter);
	mjd += toas.mjd[i];
	if(loaderrors) {
	  err = toas.err[i];
	}
      }
      /* Choose a phase according to probability function */
      if(probdist.NrBins != 0) {
	y = gsl_rng_uniform(rand_num_gen);
	for(j = 0; j < probdist.NrBins+1; j++) {
	  if(y < cumdist[j])
	    break;
	}
	/*	y = (x - j + 1)(cumdist[j]-cumdist[j-1]) + cumdist[j-1]; */
	phase0 = (y - cumdist[j-1])/(cumdist[j]-cumdist[j-1]) + j - 1;
	phase0 /= (long double)probdist.NrBins;
	/*	printf("XXXX %Lf (%d %f %f)\n", phase0, j, y, cumdist[j]); */
      }else
	phase0 = 0;
      /* Now refine mjd in order to make the TOA correspond to the fiducial phase. */
      for(j = 0; j < 3; j++) {
	phase = evaluate_ephemeris(mjd, freq, 0, 0, 0, 0, 0, site, flags, &(parfile.par), 0, parfile.paramset, 0, 0);
	phase -= tzphase + phase0;
	nturn = phase;
	phase -= nturn;
	if(phase > 0.5)
	  phase -= 1;
	if(phase < -0.5)
	  phase += 1;
	/*	printf("%Lf %13.2Le | ", mjd, phase); */
	mjd -= phase/(TEMPO3_SecPerDay*parfile.par.parameter[TEMPO3_PARAM_F0]);
      }
      /* Add the errorbar uncertainty to the TOA */
      mjd += gsl_ran_gaussian(rand_num_gen, err*1e-6/TEMPO3_SecPerDay);
      /*      printf("\n");*/
      /* Write out the TOA */
      /*      fprintf(fin, " FakeTempo3TOA%06d 0 0  %.3Lf  %.13Lf    0.00    %.2Lf        @\n", i+1, freq, mjd, err); */
      fprintf(fin, " FakeTempo3TOA%06d %15.8Lf  %.13Lf   %.2Lf  @\n", i+1, freq, mjd, err); 
    }
    
    fclose(fin);
    printf("Writing to %s is done\n", argv[argc-1]);
    free(mrq_alpha); 
    free(mrq_alpha_inv); 
    free(xstart);
    free(dx);
    free(xfit);
    //    free(dmin);
    free(mrq_deriv);
    //    free(fixed);
    //    free(paramset);
    //    free(paramlinked);
    free(additional_funk_info.best_xfit);
    for(i = 0; i < MaxNrStack; i++) {
      free_ephemeris_par_only(&(stackptr->parfiles[i]));
      free(stackptr->fixed[i]);
      free(stackptr->paramset[i]);
    }
    free(stackptr);
#ifdef TEMPO3_EnableGTK
    free_ephemeris_par_only(&gtk_parfile);
#endif
    //    free_ephemeris_par_only(&parfile_old);
    free_ephemeris_par_only(&(additional_funk_info.parfile_funk));
    if(toas_fiducial.nrtoas > 0) {
      free_toas(&toas_fiducial);
      free_residuals(&residuals_fiducial);
    }
    if(probdist.NrBins > 0)
      closePSRData(&probdist, 0, verbose_state);
    free_ephemeris(&parfile, verbose_state);
    free_ephemeris(&parfile_old, verbose_state);
    if(cleanup_tempo3_lib(verbose_state) == 0) {
      printerror(verbose_state.debug, "ERROR tempo3: Releasing tempo3 library related memory failed.");
      return 0;
    }
    return 0;
  }                        /* End of dofake */


  if(parfile.par.parameter[TEMPO3_PARAM_PEPOCH] < toas.data_mjdmin || parfile.par.parameter[TEMPO3_PARAM_PEPOCH] > toas.data_mjdmax) {
    printwarning(0, "WARNING, PEPOCH OUTSIDE DATA RANGE. PROBLEMS WITH CONVERGENCE OF SOLUTION MIGHT BE SOLVED BY FIXING PEPOCH.");
  }


  residuals.residy = malloc((toas.nrtoas + toas.nrtoasdeleted)*sizeof(float));
  residuals.resid =  malloc((toas.nrtoas + toas.nrtoasdeleted)*sizeof(long double));
  residuals.nturn =  malloc((toas.nrtoas + toas.nrtoasdeleted)*sizeof(long long));
  residuals.nturn_fixed =  calloc((toas.nrtoas + toas.nrtoasdeleted), sizeof(int));
  toas.nu = malloc((toas.nrtoas + toas.nrtoasdeleted)*sizeof(long double));
  toas.nudot = malloc((toas.nrtoas + toas.nrtoasdeleted)*sizeof(long double));
  toas.nu_mjd = malloc((toas.nrtoas + toas.nrtoasdeleted)*sizeof(long double));

  ygraph = malloc((toas.nrtoas + toas.nrtoasdeleted)*sizeof(float));
  ygraph2 = malloc((toas.nrtoas + toas.nrtoasdeleted)*sizeof(float));
  // At least when toas loaded via tempo2 these are already allocated
  if(loadtimfile_flag) {
    if(runtempo2) {
      if(toas.mjdx != NULL) {
	free(toas.mjdx);
	toas.mjdx = NULL;
      }
      if(toas.freqx != NULL) {
	free(toas.freqx);
	toas.freqx = NULL;
      }
    }
  }
  toas.mjdx  = malloc((toas.nrtoas + toas.nrtoasdeleted)*sizeof(float));
  toas.freqx  = malloc((toas.nrtoas + toas.nrtoasdeleted)*sizeof(float));
  if(residuals.residy == NULL || toas.mjdx == NULL || toas.freqx == NULL || residuals.nturn == NULL || residuals.nturn_fixed == NULL || residuals.resid == NULL || ygraph == NULL || ygraph2 == NULL || toas.nu == NULL || toas.nudot == NULL || toas.nu_mjd == NULL) {
    printerror(0, "Cannot allocate memory\n");
    return 0;
  }

  // Set the pulse numbers if flag written out with pulse number
  long double testphase;
  int some_toas_no_pulse_number = 0;
  for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
    char *txtptr;
    txtptr = strstr(toas.flags[i], "-pn");
    if(txtptr != NULL) {
      sscanf(txtptr, "-pn %lld", &(residuals.nturn[i]));   // Used to be %Ld
      residuals.nturn_fixed[i] = 1;
    }else {
      residuals.nturn[i] = 0;
      residuals.nturn_fixed[i] = 0;
      some_toas_no_pulse_number = 1;
    }
  }
  if(some_toas_no_pulse_number && parfile.par.track != 0) {
    printwarning(0, "WARNING: For one or more TOA's no pulse number was recorded in the .tim file. The TRACK parameter is set to 0.");
    parfile.par.track = 0;
  }

  if(parfile.par.track) {
    long long expectedturn;
    long double tmp_ld, expectedturn_total;
    int stop;
    expectedturn_total = 0;
    stop = 0;
    ephemeris_setLinkedValues_par_version(&(parfile.par), parfile.paramlinked);
    for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
      if(toas.deleted[i] == 0) {
	testphase = evaluate_ephemeris(toas.mjd[i], toas.freqSSB[i], toas.ssb1[i], toas.ssb2[i], toas.ssb3[i], 0, 0, toas.site[i], toas.flags[i], &(parfile.par), 0, parfile.paramset, 0, 0);
	tmp_ld = testphase - residuals.nturn[i];
	//	tmp_ld = testphase - 0.5;
	//	tmp_ld -= residuals.nturn[i];
	//      if(verbose)
	//	printf("  Expectation off by = %Ld\n", expectedturn);
	expectedturn_total += tmp_ld/(long double)(toas.nrtoas + toas.nrtoasdeleted);
	//	printf("XXXXX expectedturn_total=%Lf\n", expectedturn_total);
	if(verbose_state.verbose && stop == 0) {
	  stop = 1;
	  printf("  The calculated phase of TOA = %Lf and turn number = %lld\n", testphase, residuals.nturn[i]);  // Used to be %Ld
	}
      }
    }
    if(expectedturn_total > 0)
      expectedturn = expectedturn_total+0.5;
    else if(expectedturn_total == 0)
      expectedturn = expectedturn_total;
    else
      expectedturn = expectedturn_total-0.5;
    printf("  Expectation from pulse numbering different by = %Lf = %lld\n", expectedturn_total, expectedturn);   // Used to be %Ld
    if(llabs(expectedturn) > 3) {
      printf("    Adjusting turn numbers\n");
      for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
	if(residuals.nturn_fixed[i]) {
	  residuals.nturn[i] += expectedturn;
	}
      }
    }else {
      printf("    Decided to not adjust the turn numbers\n");
    }
  }

  for(i = 0; i < MaxNrStack; i++) {
    stackptr->deleted[i] = malloc((toas.nrtoas + toas.nrtoasdeleted)*sizeof(int));
    if(stackptr->deleted[i] == NULL) {
      printerror(0, "Cannot allocate memory\n");
      return 0;
    }
  }


  for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
    toas.mjdx[i] = toas.mjd[i];  // Should already be done in add_toa in readtimfile_custom, but somehow gives difference? Because of sorting?
    toas.freqx[i] = toas.freqSSB[i];
    //    toas.err[i] *= parfile.par.parameter[TEMPO3_PARAM_F0];    // Convert the errors from micro seconds to micro phase
  }

#ifdef TEMPO3_EnableGTK
  if(enableGTK_flag) {
    initGTKwindow(argc, argv);
    GTKinitialised = 1;
  }
  int weshouldsetalarm = 1;
#ifdef HAVEPGPLOT2 // If in pgplot2 mode, let pgplot deal with setting the alarm function
  if(enablePGPLOT2 == 0) {  // If not in local mode, let pgplot2 know that our timer should be called as well.
    printf("XXXXXXX pgplot2 deals with alarm\n");
    pgplot2_gtk_server_user_alarm_function(alarm_wakeup);
    weshouldsetalarm = 0;
  }
#endif
  //  weshouldsetalarm = 1;
  //  int remove_line_above = 1;
  if(weshouldsetalarm) {
    /* Start the timer that will periodically process the pending gtk operations. */
    //  printf("Starting timer\n");
    printf("XXXXXXX tempo3 deals with alarm\n");
    if(signal(SIGALRM, (void (*)(int))alarm_wakeup) == SIG_ERR) {
      printerror(0, "tempo3: Cannot set wakeup timer signal.\n");
      perror("error calling SIGALRM signal");
      return 0;
    }
    struct itimerval tout_val;
    tout_val.it_value.tv_sec = 0;
    tout_val.it_value.tv_usec = 100000;   /* function executed 10 times per second */
    tout_val.it_interval = tout_val.it_value;
    if(setitimer(ITIMER_REAL, &tout_val, NULL) == -1) {
      printerror(0, "tempo3: Cannot set wakeup timer signal.\n");
      perror("error calling setitimer");
      return 0;
    }
  }

  //  printf("Starting timer completed\n");

  //  for(i = 0; i < 1000; i++) {
  //    sleep(1);
  //  }

#endif


  /*
  //  int pgplot_initialised = 0;
#ifdef HAVEPGPLOT2 
  if(enablePGPLOT2 && enableGTK_flag) {
    if(pgplot2_gtk_exclusively_local() == 0) {
      printerror(0, "tempo3: Cannot set pgplot2 in exclusive local mode.\n");
      return 0;
    }
    if(deviceDefined) {
      printerror(0, "tempo3: Cannot set use the -pgplot2 and -device options together.\n");
      return 0;
    }
    deviceID = pgplot2_gtk_open_device_in_local_mode(gtk_drawing_area, 800, 600, 1);
    pgplot2_gtk_local_add_pgplot_menu_options(gtk_menuBar);
    // Really should open device in initgtk!!!!
    gtk_widget_show_all(gtk_widget_get_toplevel(gtk_menuBar));
    pgplot_initialised = 1;
  }
#endif
  */


  finderrors = 1;
  int pgplot_initialised = 0;
#ifdef HAVEPGPLOT2 
  if(enablePGPLOT2 && enableGTK_flag) {
    deviceID = deviceID_pgplot2;
    pgplot_initialised = 1;
    if(deviceID <= 0) {
      printerror(0, "tempo3: Cannot open plotting device!\n");
      return 0;
    }
  }
#endif
  if(pgplot_initialised == 0) {
    if(deviceDefined) {
      deviceID = ppgopen(argv[deviceDefined]);
    }else {
#ifdef HAVEPGPLOT2 
      deviceID = ppgopen("11196127/pw");
      if(deviceID <= 0) {
	deviceID = ppgopen("/pw");
      }
#else
      deviceID = ppgopen("11196127/xs");
      if(deviceID <= 0) {
	deviceID = ppgopen("/xs");
      }
#endif
    }
    if(deviceID <= 0) {
      printerror(0, "tempo3: Cannot open plotting device!\n");
      return 0;
    }
  }

  do_show_params = 1;
  /* Not necessary, but saves a warning */
  //  nrparamsused = 0;
  /* Store the start situation in the stack */
  storeInStack(stackptr, &(parfile.par), toas, parfile.wraps, previous_xstart, previous_xend, parfile.fixed, parfile.paramset);
  do {
    if(do_show_params) {
      if(verbose_state.debug) {
	printf("Updating shown parameters\n");
      }
      if(autoquiet == 0) {
	printf("\n");
	print_ephemeris_ExtraParams(stdout, &(parfile.par), parfile.paramset, parfile.paramlinked, finderrors, parfile.par.dplus, distance, distanceerror, hierarchical);
	printf("PARAMETER   VALUE                                        FIT        DIFF       ERROR    SIGMA\n");
      }
      if(doauto != 1 || autoquiet == 0)  // For first pass of this point finderrors is set to 1, making output to be generated even if autoquiet is set. So when in doauto mode, ignore printing of the first step if in quiet mode
	if(print_ephemeris(stdout, &parfile, &parfile_old, &toas, 1, 1, 0, !finderrors, autoquiet, 0, verbose_state) == 0) {
	  printerror(verbose_state.debug, "ERROR tempo3: Showing ephemeris failed.");
	  return 0;
	}
      if(autoquiet == 0) {
	/* Don't copy new parfile yet, because it is needed when button is clicked during pgband */
	/*      do_show_params = 0;
		memcpy(&parfile_old, &parfile, sizeof(tempo3_parameters_def));*/
	printf("\n");
	printf("Press ? in pgplot window for help\n");
      }
      if(verbose_state.debug) {
	printf("Updating shown parameters done\n");
      }
    }

    if(verbose_state.debug) {
      printf("Updating plot\n");
    }
    if(doplot("?", ygraph, ygraph2, 0, plottype, showdata, showWaves, toas, &residuals, &parfile, additional_funk_info.superverbose, writeoutpieces, plottitle, &cur_x1, &cur_y1, &cur_x2, &cur_y2, previous_xstart, previous_xend, toas_fiducial, residuals_fiducial, showfiducial, showpepoch, highlightedFlag, nummericalApproxNuNudot, title, amoeba_algorithm, markersize, stridefit_nrdays, stridefit_mintoarange, stridefit_minstride, stridefit_minnroftoas, stridefit_includeF2, stridefit_includeF3, stridefit_scriptfilename, nosort, numthreads, verbose_state) == 0)
      return 0;
    showdata = 0;
    writeoutpieces = 0;
    if(verbose_state.debug) {
      printf("Updating plot done\n");
    }

    /*    printf("click to see result after fixing.\n");
    ppgcurs(&xpos, &ypos, &ch);
    subtractTurns(toas, &residuals, parfile.wraps, 1);
    pgplot_options.viewport.dontopen = 1;
    pgplot_options.viewport.dontopen = 1;
    for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
      residuals.residy[i] = residuals.resid[i];
    }
    if(pgplotGraph1(&pgplot_options, residuals.residy, toas.mjdx, NULL, toas.nrtoas, -1, -1, residuals.mjdmin, residuals.mjdmax, 0, "MJD", "Resid", "", 1, 1, 17, 3, 1, NULL, -1, noverbose) == 0) {
      return 0;
    }
    */
  
    if(dumpphases != 0) {
      if(verbose_state.debug) {
	printf("Dumping toa's\n");
      }
      for(i = 0; i < toas.nrtoas; i++) {
	printf("DUMPTOA%d %.13Lf %Lf\n", i+1, toas.mjd[i], residuals.resid[i]);
      }
      if(verbose_state.debug) {
	printf("Dumping toa's done\n");
      }
      return 0;
    }

    ppgsci(1);
    if(doauto == 0) {
      /*      ppgcurs(&xpos, &ypos, &ch); */
#ifdef TEMPO3_EnableGTK
      copy_parfile(&gtk_parfile, &(parfile.par));
      //      memcpy(&gtk_parfile, &parfile, sizeof(tempo3_parameters_def));
      gtk_finderrors = finderrors;
      //      if(enableGTK_flag) {
      //	while(gtk_events_pending()) {
      //	  gtk_main_iteration();
      //	}
      //      }
#endif
      if(superfit == 0) {
	ch = 0;
	if(macroname[0] != 0) {
	  if(macro_open == 0) {
	    fmacro = fopen(macroname, "r");
	    if(fmacro == NULL) {
	      printerror(0, "tempo3: Cannot open macro '%s'\n", macroname);
	      return 0;
	    }
	    macro_open = 1;
	    printf("Opening macro %s\n", macroname); 
	  }
	  do {
	    ch = fgetc(fmacro);
	  }while(ch == '\n' || ch == '\r');
	  if(ch == '^') {  /* Need to interpret next character as being ctrl-character */
	    ch = fgetc(fmacro);
	    switch(ch) {
	    case ' ': ch = 0; ctrlspace = 1; break;
	    case 'd': ch = 4; break;
	    case 'f': ch = 6; break;
	    case 'n': ch = 14; break;
	    case 'p': ch = 16; break;
	    case 'r': ch = 18; break;
	    case 't': ch = 20; break;
	    case 'w': ch = 23; break;
	    case 'x': ch = 24; break;
	    default: printerror(0, "tempo3: Don't understand ^%c (%d) in macro.\n", ch, ch); return 0;
	    }
	  }
	  if(ch == 'R') {  // See if it spells out the word RETURN
	    int ok;
	    long pos;
	    pos = ftell(fmacro);
	    //	    printf("POSITION: %ld\n", pos);
	    ok = 1;
	    if(fgetc(fmacro) != 'E') {
	      ok = 0;
	    }else {
	      if(fgetc(fmacro) != 'T') {
		ok = 0;
	      }else {
		if(fgetc(fmacro) != 'U') {
		  ok = 0;
		}else {
		  if(fgetc(fmacro) != 'R') {
		    ok = 0;
		  }else {
		    if(fgetc(fmacro) != 'N') {
		      ok = 0;
		    }
		  }
		}
	      }
	    }
	    if(ok) {
	      ch = 13;
	      printf("Found word RETURN in macro\n");
	      //	      printf("POSITION: %ld\n",  ftell(fmacro));
	    }else {
	      fseek(fmacro, pos, SEEK_SET);
	      ch = fgetc(fmacro);
	    }
	  }
	  /*	  fprintf(stderr, "XXXXXX Getting %c (%d) from macro\n", ch, ch); */
	  if(ch == EOF) {
	    ch = 0;
	    printf("Reached EOF in macro\n");
	  }else {
	    fflush(stdout);
	    if(ch == 13)
	      printf("Executing RETURN (%d) from macro\n", ch); 
	    else if(ch == 0)
	      printf("Executing ctrl-space (%d) from macro\n", ch); 
	    else
	      printf("Executing %c (%d) from macro\n", ch, ch); 
	    fflush(stdout);
	  }
	  if(feof(fmacro)) {  /* Stopping executing macro */
	    fclose(fmacro);
	    ch = 0;
	    macro_open = 0;
	    macroname[0] = 0;
	    printf("Closing macro file\n");
	  }
	}
	if(ch == 0 && ctrlspace == 0)
	  ppgband(7, 0, 0, 0, &xpos, &ypos, &ch);
	ctrlspace = 0;
      }else {
	ch = ' ';
	superfit--;
	if(superfit)
	  printf("\nThere are %d fit iterations left to do\n\n", superfit);
      }
    }else if(doauto == 1) {   // Switch off error calculation for 1st fit
      ch = 'e';
      doauto++;
    }else if(doauto == 2) {   // Fit for phase only first
      ch = 13;
      doauto++;
    }else if(doauto == 3) {   // Enable all parameters
      ch = 'a';
      doauto++;
    }else if(doauto == 4) {   // Do first fit after fixing certain parameters
      if(automode >= 1 && automode <= 2) {
	parfile.fixed[TEMPO3_PARAM_DM] = 1;
      }
      if(automode == 2) {
	parfile.fixed[TEMPO3_PARAM_NBRAKE] = 1;
      }
      ch = ' ';
      doauto++;
    }else if(doauto == 5) {   // Switch on error calculation for 2nd fit
      ch = 'e';
      doauto++;
    }else if(doauto == 6) {   // Do 2nd fit
      ch = ' ';
      doauto++;
    }else if(doauto == 7) {  // Quit
      ch = 'q';
      doauto++;
    }
    if(do_show_params) {
      do_show_params = 0;
      copy_parfile(&(parfile_old.par), &(parfile.par));
      //      memcpy(&parfile_old, &parfile, sizeof(tempo3_parameters_def));
    }
    param = ch;
    if(param == 'q' || param == 27) {
    }else if(param == 'e') {
      if(finderrors) {
	if(autoquiet == 0)
	  printf("Error calculation is switched off.\n");
	finderrors = 0;
      }else {
	if(autoquiet == 0)
	  printf("Error calculation is switched on.\n");
	finderrors = 1;
      }
      /*      additional_funk_info.superverbose = -1; */
    }else if(param == '?' || param == 'h') {
      printf("Commands (in pgplot window):\n");
      printf("  Change value of, or toggle fitting for fit parameters:\n");
      printf("    c             = Change the value of a parameter (fit parameters/options)\n");
      printf("    t             = Toggle a given fit parameter (especially useful if not defined below)\n");
      printf("    0 .. 9        = toggle F0 .. F9 fitting\n");
      //      printf("    )             = toggle F10... fitting\n");
      printf("    b             = toggle various binary parameters fitting\n");
      printf("    p             = toggle various position parameters fitting\n");
      //      printf("    N             = toggle NBRAKE fitting\n");
      printf("    g             = toggle glitch fit params\n"); 
      printf("    w             = toggle wave fitting\n");
      printf("    d             = toggle DM fitting\n");
      printf("    j             = toggle jump fitting\n");
      printf("    m             = toggle MODE (error bars)\n"); 
      printf("    v             = toggle value fitting (new VALUEX parameter)\n"); 
      printf("    n/a           = toggle none or all fit parameters\n");
      printf("  Fitting commands:\n");
      printf("    SPACE/RETURN  = do the fit. By pressing return the parameters are fit one by one.\n");
      printf("    ctrl-SPACE    = do the fit %d times.\n", nrSuperFit);
      printf("    ctrl-n        = change the number of times ctrl-SPACE fits.\n");
      printf("    ctrl-f        = toggle between fit algorithm to be used.\n");
      printf("    e             = toggle determination of errors on fit parameters\n"); 
      printf("    F             = Refine TZR params to corresponds to residual=0\n"); 
      printf("    r             = restore/pop stack (undo change made in par file)\n"); 
      //      printf("    t             = set fit tolerance (%.2Le)\n", ftol);
      printf("    V             = toggle verbose\n");
      printf("  Residuals:\n");
      printf("    +/-/BACKSPACE = add/subtract/remove phase wraps\n");
      printf("    s/f/u         = set start,finish mjd, or unzoom\n");
      printf("    x/right click = delete point\n");
      printf("    ctrl-x        = delete all visuable points\n");
      printf("  Plotting:\n");
      printf("    H             = highlight points with a certain flag value\n");
      printf("    ctrl-p        = try to plot profile (assumed to be the name of the toa)\n");
      printf("    S             = show different style plots\n");
      printf("    ctrl-t        = toggle various plotting options\n");
      printf("  Input/output:\n");
      printf("    D             = show data as well as model (in combination with nu/nudot plots)\n");
      printf("    o             = write copy graphics to other pgplot device\n");
      printf("    O             = write copy graphics to script\n");
      printf("    P             = write par file\n");
      printf("    q             = end\n");
      printf("    T             = write tim file\n");
      printf("    W             = show fitwaves contribution\n");
      printf("    ctrl-w        = show fitwaves contribution from certain order\n");
      printf("    ctrl-d        = Write out the pieces of data (D option)\n");
      printf("    ctrl-r        = Write out residuals\n");
    }else if(param == 'D') {   /* D */
      showdata = 1;
    }else if(param == 20) {  /* Ctrl-t, toggle various parameters */
      printf("Commands to toggle:\n");
      printf("  f = Toggle showing fiducial point (TZ parameters)\n");
      printf("  p = Toggle showing PEPOCH\n");
      printf("  t = Toggle showing title\n");
      printf("  m = Change marker size\n");
      if(macroname[0] != 0) {
	int ret;
	do {
	  ret = fscanf(fmacro, "%c", &ch);
	  //	    printf("XXXXX %d '%c'\n", ch, ch);
	}while(ret == 1 && (ch == '\n' || ch == '\r'));
      }else {
	ppgcurs(&xpos, &ypos, &ch);
      }
      key = ch;
      if(key == 't') {
	if(plottitle)
	  plottitle = 0;
	else
	  plottitle = 1;
	do_show_params = 1;
      }else if(key == 'f') {
        if(toas_fiducial.nrtoas != 1) {
	  printwarning(0, "NO FIDUCIAL POINT LOADED!!");
	}else {
	  do_show_params = 1;
	}
	if(showfiducial == 0) 
	  showfiducial = 1;
	else
	  showfiducial = 0;
      }else if(key == 'p') {
	do_show_params = 1;
	if(showpepoch == 0) 
	  showpepoch = 1;
	else
	  showpepoch = 0;
      }else if(key == 'm') {
	printf("Change marker size (currently %f): ", markersize);
	fflush(stdout);
	if(macroname[0] != 0) {
	  fscanf(fmacro, "%f", &markersize);
	}else {
	  scanf("%f", &markersize);
	}
	printf("Change marker size into %f\n", markersize);
      }
    }else if(param == 4) {  /* ctrl-d */
      printf("\nThe spin parameters will be measured using TOAs in short windows of time. Each window is initially centered on one of the TOAs (the \"central TOA\") and includes other nearby TOAs, although the actual epoch of each measurement is taken to be halfway the first and last included TOA. Each TOA might be used in more than one window, hence making the output measurements dependent of each other. A window is always truncated at a glitch epoch.\n");
      printf("\nFull window length (currently set to %.2f days). Only TOAs which are up to half this value (%.2f) away from the central TOA will be considered: ", stridefit_nrdays, 0.5*stridefit_nrdays);
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%f", &stridefit_nrdays);
      }else {
	scanf("%f", &stridefit_nrdays);
	printf("%f\n", stridefit_nrdays);
      }

      printf("\nMinimum nr of days covered by TOAs (largest TOA - smallest TOA within the window, currently set to %.2f days): ", stridefit_mintoarange);
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%f", &stridefit_mintoarange);
      }else {
	scanf("%f", &stridefit_mintoarange);
	printf("%f\n", stridefit_mintoarange);
      }

      printf("\nAny window which is offset (difference between central TOAs) by less than the \"minimum stride\" from the previous window will be discarded. The TOAs are assumed to be ordered in MJD, which should be fine as long as the -nosort option wasn't used when calling tempo3. The minimum stride is currently set to %.2f days: ", stridefit_minstride);
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%f", &stridefit_minstride);
      }else {
	scanf("%f", &stridefit_minstride);
	printf("%f\n", stridefit_minstride);
      }

      printf("\nDiscard all windows containing less than this minimum nr of TOAs (currently set to %d): ", stridefit_minnroftoas);
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%d", &stridefit_minnroftoas);
      }else {
	scanf("%d", &stridefit_minnroftoas);
	printf("%d\n", stridefit_minnroftoas);
      }

      printf("\nInclude F2 in fit (1 = yes, 0 = no, 2 = instead use NBREAK as fit paramter): ");
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%d", &stridefit_includeF2);
      }else {
	scanf("%d", &stridefit_includeF2);
      }
      if(stridefit_includeF2 == 1) {
	printf("INCLUDE F2\n");
      }else if(stridefit_includeF2 == 2) {
	printf("INCLUDE NBREAK\n");
      }else {
	stridefit_includeF2 = 0;
	printf("NOT INCLUDING F2 OR NBREAK\n");
      }

      printf("\nInclude F3 in fit (1 = yes, 0 = no): ");
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%d", &stridefit_includeF3);
      }else {
	scanf("%d", &stridefit_includeF3);
      }
      if(stridefit_includeF3 == 1) {
	stridefit_includeF3 = 1;
	printf("YES\n");
      }else {
	stridefit_includeF3 = 0;
	printf("NO\n");
      }

      printf("\nFile name of script to be generated to fit pieces (- means skip): ");
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%s", stridefit_scriptfilename);
      }else {
	scanf("%s", stridefit_scriptfilename);
      }
      if(strcmp(stridefit_scriptfilename, "-") == 0) {
	stridefit_scriptfilename[0] = 0;
      }else {
	printf("%s\n", stridefit_scriptfilename);
      }

      writeoutpieces = 1;
      showdata = 1;
    }else if(param == 'H') {
      printf("Type flag to highlight: ");
      fflush(stdout);
      if(macroname[0] != 0) {
	fgets(highlightedFlag, 1000, fmacro);
	if(strlen(highlightedFlag) == 0 || highlightedFlag[0] == '\r' || highlightedFlag[0] == '\n') {  // Check if no flag was entered on same line, then probably on next one
	  fgets(highlightedFlag, 1000, fmacro);
	}
      }else {
	fgets(highlightedFlag, 1000, stdin);
      }
      if(highlightedFlag[strlen(highlightedFlag)-1] == 10)
	highlightedFlag[strlen(highlightedFlag)-1] = 0;
      printf("Highlighting flag = '%s'\n", highlightedFlag);
    }else if(param == 'W') {          /* W */
      if(showWaves == 0) {
	showWaves = 1;  /* showWaves is the first fitwave to be excluded, unless it is zero */
	printf("All fitwave parameters are excluded from timing model\n");
      }else {
	showWaves = 0;
	printf("All fitwave parameters are included in timing model\n");
      }
    }else if(param == 23) {          /* ctrl-w */
      printf("First order of fitwave parameter to exclude from timing model: ");
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%d", &showWaves);
      }else {
	scanf("%d", &showWaves);
      }
      printf("First fitwave parameter excluded from timing model is WAVE%d\n", showWaves);
    }else if(param == 's') {
      if(plottype == PLOTTYPE_FREQUENCY)
	residuals.freqmin = xpos;
      else
	residuals.mjdmin = xpos;
    }else if(param == 'f') {
      if(plottype == PLOTTYPE_FREQUENCY)
	residuals.freqmax = xpos;
      else
	residuals.mjdmax = xpos;
    }else if(param == 'u') {
      if(plottype == PLOTTYPE_FREQUENCY) {
	j = 0;
	for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
	  if(toas.deleted[i] == 0) {
	    if(j == 0) {
	      residuals.freqmin = toas.freqSSB[i];
	      residuals.freqmax = toas.freqSSB[i];
	      j = 1;
	    }else {
	      if(toas.freqSSB[i] > residuals.freqmax)
		residuals.freqmax = toas.freqSSB[i];
	      if(toas.freqSSB[i] < residuals.freqmin)
		residuals.freqmin = toas.freqSSB[i];
	    }
	  }
	}
      }else {
	residuals.mjdmax = toas.data_mjdmax;
	residuals.mjdmin = toas.data_mjdmin;
      }
    }else if(param == '+') {
      storeInStack(stackptr, &(parfile.par), toas, parfile.wraps, previous_xstart, previous_xend, parfile.fixed, parfile.paramset);
      parfile.wraps->phase[parfile.wraps->nrwraps] = +1;
      parfile.wraps->mjd[parfile.wraps->nrwraps++] = xpos;
    }else if(param == '-') {
      storeInStack(stackptr, &(parfile.par), toas, parfile.wraps, previous_xstart, previous_xend, parfile.fixed, parfile.paramset);
      parfile.wraps->phase[parfile.wraps->nrwraps] = -1;
      parfile.wraps->mjd[parfile.wraps->nrwraps++] = xpos;
    }else if(param == 8) {
      storeInStack(stackptr, &(parfile.par), toas, parfile.wraps, previous_xstart, previous_xend, parfile.fixed, parfile.paramset);
      parfile.wraps->nrwraps = 0;
    }else if(param == 'V') {   /* V */
      if(additional_funk_info.superverbose == 0) {
	printf("Verbose is off\n");
	additional_funk_info.superverbose = -1;
      }else if(additional_funk_info.superverbose == -1) {
	printf("Verbose is maximal\n");
	additional_funk_info.superverbose = 1;
      }else {
	printf("Verbose is normal\n");
	additional_funk_info.superverbose = 0;
      }
    }else if(param == 6) {   /* ctrl-f */
      if(fit_algorithm == 0) {
	printf("Switching to the Levenberg-Marquart algorithm\n");
	fit_algorithm = 1;
      }else {
	printf("Switching to the downhill simplex algorithm\n");
	fit_algorithm = 0;
      }
    }else if(param == 'c' || param == 't') {
      if(param == 'c') {
	printf("User requested the change of a parameter. The paramters that can changed include:\n");
	printf("  F0 etc., i.e. parameters as can be defined in parfile\n");
	printf("  TOL     - Tolerance as used during fitting\n");
	printf("Type in terminal the parameter to change, i.e. F0: ");
      }else {
	printf("User requested the toggling the fitting of a parameter.\n");
	printf("Type in terminal the parameter to change, i.e. F0: ");
      }
      fflush(stdout);
      char textBuffer[1000];
      textBuffer[0] = 0;
      if(macroname[0] != 0) {
	fscanf(fmacro, "%s", textBuffer);
      }else {
	scanf("%s", textBuffer);
      }
      if(strcasecmp(textBuffer, "TOL") == 0 && param == 'c') {
	printf("New tolerace (type in xterm, currently set to %Le): ", ftol);
	fflush(stdout);
	if(macroname[0] != 0) {
	  fscanf(fmacro, "%Lf", &ftol);
	}else {
	  scanf("%Lf", &ftol);
	}
	printf("Tolerance is now set to %Le\n", ftol);
      }else {
	int found = 0;
	long curparam;
	int altid;
	for(curparam = 0; curparam < TEMPO3_PARAM_TOTAL; curparam++) {
	  for(altid = 0; altid < 2; altid++) {
	    if(tempo3_parfile_descr.identifiers[curparam][altid] != NULL) {
	      if(strcasecmp(textBuffer, tempo3_parfile_descr.identifiers[curparam][altid]) == 0) {
		found = 1;
		// If it is a string or integer, or a selection of other parameters which are more complicated to implement, then don't allow the user to change it.
		if(strcmp(tempo3_parfile_descr.print_format[curparam], "%s") == 0 || strcmp(tempo3_parfile_descr.print_format[curparam], "%d") == 0 || (curparam >= TEMPO3_PARAM_VALUE && curparam < TEMPO3_PARAM_VALUE+TEMPO3_MaxNrValues) || curparam == TEMPO3_PARAM_TRES) {
		  if(param == 'c') {
		    printerror(0, "Parameter \"%s\" cannot be modified. Please update the ephemeris manually and run tempo3 again.\n", tempo3_parfile_descr.identifiers[curparam][altid]);
		  }else {
		    printerror(0, "Parameter \"%s\" cannot be fitted for.\n", tempo3_parfile_descr.identifiers[curparam][altid]);
		  }
		}else {
		  if(param == 'c') {
		    printf("New %s (type in xterm, currently set to %Le): ", tempo3_parfile_descr.identifiers[curparam][altid], parfile.par.parameter[curparam]);
		    fflush(stdout);
		    if(macroname[0] != 0) {
		      fscanf(fmacro, "%Lf", &(parfile.par.parameter[curparam]));
		    }else {
		      scanf("%Lf", &(parfile.par.parameter[curparam]));
		    }
		    printf("%s is set to %Le\n", tempo3_parfile_descr.identifiers[curparam][altid], parfile.par.parameter[curparam]);
		    parfile.paramset[curparam] = 1;
		    do_show_params = 1;
		  }else {
		    if(parfile.paramset[curparam] == 0) {
		      printerror(0, "Parameter \"%s\" is not initialised. Use the 'c' option to initialise it first.\n", tempo3_parfile_descr.identifiers[curparam][altid]);
		    }else {
		      toggleFixed(curparam, parfile.fixed, parfile.paramlinked);
		      do_show_params = 1;
		    }
		  }
		}
		break;
	      }
	    }
	  }
	  if(found)
	    break;
	}
	if(found == 0) {
	  printerror(0, "Parameter \"%s\" is not recognized. Run \"tempo3 -paramlist\" to obtain a list of supported parameters.\n", textBuffer);
	}
      }
    }else if(param == 'P') {
      printf("Output filename: ");
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%s", filename);
      }else {
	scanf("%s", filename);
      }
      fin = fopen(filename, "w");
      if(fin == NULL) {
	printerror(0, "Cannot open %s\n", filename);
	return 0;
      }
      if(print_ephemeris(fin, &parfile, NULL, &toas, 0, 0, 1, 1, autoquiet, writehex, verbose_state) == 0) {
	printerror(verbose_state.debug, "ERROR tempo3: Showing ephemeris failed.");
	return 0;
      }
      fclose(fin);
      printf("Writing done to %s\n", filename);
    }else if(param == 'T') {
      printf("Output filename (xterm): ");
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%s", filename);
      }else {
	scanf("%s", filename);
      }
      if(writetimfile(filename, toas, residuals.nturn, parfile.wraps, 1, verbose_state) == 0) {
	printerror(verbose_state.debug, "Writing timfile failed");
	return 0;
      }
    }else if(param == 14) {                        /* ctrl-n change nr times ctrl-SPACE fits */
      printf("Type the number of time (it is now %d) you want ctrl-SPACE to fit (in xterm): ", nrSuperFit);
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%d", &nrSuperFit);
      }else {
	scanf("%d", &nrSuperFit);
      }
      printf("Changed to %d\n", nrSuperFit);
    }else if(param == 18) {                        /* ctrl-r write out residuals */
      printf("Output filename (xterm): ");
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%s", filename);
      }else {
	scanf("%s", filename);
      }
      fin = fopen(filename, "w");
      if(fin == NULL) {
	printerror(0, "Cannot open %s\n", filename);
	return 0;
      }
      for(i = 0; i < toas.nrtoas; i++) {
	fprintf(fin, "%.13Lf %.13Le %.13Le\n", toas.mjd[i], residuals.resid[i], toas.err[i]*parfile.par.parameter[TEMPO3_PARAM_F0]);
      }
      fclose(fin);
      printf("Writing done to %s\n", filename);

      /*    }else if(param == '\\') {
      if(nummericalApproxNuNudot) {
	printf("Now using the analytic function\n");
	nummericalApproxNuNudot = 0;
      }else {
	printf("Now using the nummerical approximation of analytic function\n");
	nummericalApproxNuNudot = 1;
	}*/
    }else if(param == 'g') {
      do_show_params = 1;
      param = 1000;
      if(parfile.par.nrglitches > 1) {
	printf("Glitch nr (type in pgplot window, press 'a' for all glitches):\n");
	if(macroname[0] != 0) {
	  int ret;
	  do {
	    ret = fscanf(fmacro, "%c", &ch);
	    //	    printf("XXXXX %d '%c'\n", ch, ch);
	  }while(ret == 1 && (ch == '\n' || ch == '\r'));
	}else {
	  ppgcurs(&xpos, &ypos, &ch);
	}
	key = ch;
      }else {
	key = '1';
      }
      if((key >= '1' && key <= '9') || key == 'a') {
	glitchnr = key-'0';
	if(key == 'a')
	  glitchnr = -1;
	param += 100*glitchnr;
	printf("Commands:\n");
	printf("  0 = GLF0_%c     (perminent F0 change)\n", key);
	printf("  1 = GLF1_%c     (perminent F1 change)\n", key);
	printf("  2 = GLF2_%c     (perminent F2 change)\n", key);
	printf("  3 = GLF3_%c     (perminent F3 change)\n", key);
	printf("  4 = GLF4_%c     (perminent F4 change)\n", key);
	printf("  5 = GLF5_%c     (perminent F5 change)\n", key);
	printf("  a = GLF0D_%c    (non-perminent F0 change)\n", key);
	//	printf("  b = GLF1D_%c    (non-perminent F1 change)\n", key);
	printf("  e = GLEP_%c     (toggle fitting for glitch epoch)\n", key);
	printf("  E = GLEP_%c     (set glitch epoch)\n", key);
	printf("  p = GLPH_%c     (phase step)\n", key);
	printf("  t = GLTD_%c     (timescale of changes)\n", key);
	printf("  N = GLNBRAKE_%c (perminent change in braking index)\n", key);
	printf("  q = back\n");

	if(macroname[0] != 0) {
	  int ret;
	  do {
	    ret = fscanf(fmacro, "%c", &ch);
	    //	    printf("XXXXX %d '%c'\n", ch, ch);
	  }while(ret == 1 && (ch == '\n' || ch == '\r'));
	}else {
	  ppgcurs(&xpos, &ypos, &ch);
	}
	key = ch;
	param += (key-'0');
	if(key >= '0' && key <= '5') {
	  //	  printf("Toggling GLF%d_%d\n", key-'0', glitchnr);
	  if(glitchnr >= 0) {
	    toggleFixed(TEMPO3_PARAM_GLF0+(key-'0')*TEMPO3_MaxNrGlitches+glitchnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(glitchnr = 1; glitchnr <= parfile.par.nrglitches; glitchnr++) {
	      toggleFixed(TEMPO3_PARAM_GLF0+(key-'0')*TEMPO3_MaxNrGlitches+glitchnr-1, parfile.fixed, parfile.paramlinked);
	    }
	    glitchnr = -1;
	  }
	}else if(key == 'p') {
	  if(glitchnr >= 0) {
	    toggleFixed(TEMPO3_PARAM_GLPH+glitchnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(glitchnr = 1; glitchnr <= parfile.par.nrglitches; glitchnr++) {
	      toggleFixed(TEMPO3_PARAM_GLPH+glitchnr-1, parfile.fixed, parfile.paramlinked);
	    }
	    glitchnr = -1;
	  }
	}else if(key == 'e') {
	  if(glitchnr >= 0) {
	    toggleFixed(TEMPO3_PARAM_GLEP+glitchnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(glitchnr = 1; glitchnr <= parfile.par.nrglitches; glitchnr++) {
	      toggleFixed(TEMPO3_PARAM_GLEP+glitchnr-1, parfile.fixed, parfile.paramlinked);
	    }
	    glitchnr = -1;
	  }
	}else if(key == 'a') {
	  if(glitchnr >= 0) {
	    toggleFixed(TEMPO3_PARAM_GLF0D+glitchnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(glitchnr = 1; glitchnr <= parfile.par.nrglitches; glitchnr++) {
	      toggleFixed(TEMPO3_PARAM_GLF0D+glitchnr-1, parfile.fixed, parfile.paramlinked);
	    }
	    glitchnr = -1;
	  }
	  /*
	}else if(key == 'b') {
	  if(glitchnr >= 0) {
	    toggleFixed(TEMPO3_PARAM_GLF1D+glitchnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(glitchnr = 1; glitchnr <= parfile.par.nrglitches; glitchnr++) {
	      toggleFixed(TEMPO3_PARAM_GLF1D+glitchnr-1, parfile.fixed, parfile.paramlinked);
	    }
	    glitchnr = -1;
	    } */
	}else if(key == 't') {
	  if(glitchnr >= 0) {
	    toggleFixed(TEMPO3_PARAM_GLTD+glitchnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(glitchnr = 1; glitchnr <= parfile.par.nrglitches; glitchnr++) {
	      toggleFixed(TEMPO3_PARAM_GLTD+glitchnr-1, parfile.fixed, parfile.paramlinked);
	    }
	    glitchnr = -1;
	  }
	}else if(key == 'N') {
	  if(glitchnr >= 0) {
	    toggleFixed(TEMPO3_PARAM_GLNBRAKE+glitchnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(glitchnr = 1; glitchnr <= parfile.par.nrglitches; glitchnr++) {
	      toggleFixed(TEMPO3_PARAM_GLNBRAKE+glitchnr-1, parfile.fixed, parfile.paramlinked);
	    }
	    glitchnr = -1;
	  }
	}else if(key == 'E') { 
	  if(glitchnr == parfile.par.nrglitches+1) {
	    printf("Adding glitch\n");
	    parfile.par.nrglitches += 1;
	  }
	  printf("Set glitch epoch to: \n");
	  fflush(stdout);
	  scanf("%Lf", &parfile.par.parameter[TEMPO3_PARAM_GLEP+glitchnr-1]);
	}
      }
    }else if(param == 'b') {
      int companionnr;
      int dotimederivs;
      dotimederivs = 0;
      do_show_params = 1;
      param = 1000;
      if(parfile.par.companion_defined[1] != 0) {
	printf("Companion nr (type in pgplot window, press 'a' for all companions):\n");
	if(macroname[0] != 0) {
	  int ret;
	  do {
	    ret = fscanf(fmacro, "%c", &ch);
	    //	    printf("XXXXX %d '%c'\n", ch, ch);
	  }while(ret == 1 && (ch == '\n' || ch == '\r'));
	}else {
	  ppgcurs(&xpos, &ypos, &ch);
	}
	key = ch;
      }else {
	key = '1';
      }
      if((key >= '1' && key <= '9') || key == 'a') {
	companionnr = key-'0';
	if(key == 'a')
	  companionnr = -1;
	do {
	  key = '0' + companionnr;
	  if(companionnr == 'a')
	    key = 'a';
	  if(dotimederivs == 0) {
	    printf("Commands:\n");
	    printf("  . = Switch to time derivatives\n");
	    printf("  a = A1_%c                  \n", key);
	    printf("  g = GAMMA_%c               \n", key);
	    printf("  e = ECC_%c or EPS1_%c/EPS2_%c\n", key, key, key);
	    printf("  o = OM_%c                  \n", key);
	    printf("  p = PB_%c                  \n", key);
	    printf("  t = T0_%c or T0_%c          \n", key, key);
	    printf("  q = back\n");
	  }else {
	    printf("  . = Switch back to non-time derivatives\n");
	    printf("  a = A1DOT_%c      \n", key);
	    printf("  o = OMDOT_%c      \n", key);
	    printf("  p = PBDOT_%c      \n", key);
	    printf("  q = back\n");
	  }
	  if(macroname[0] != 0) {
	    int ret;
	    do {
	      ret = fscanf(fmacro, "%c", &ch);
	      //	    printf("XXXXX %d '%c'\n", ch, ch);
	    }while(ret == 1 && (ch == '\n' || ch == '\r'));
	  }else {
	    ppgcurs(&xpos, &ypos, &ch);
	  }
	  key = ch;
	  if(key == '.') {
	    if(dotimederivs)
	      dotimederivs = 0;
	    else
	      dotimederivs = 1;
	  }
	}while(key == '.');
	//	param += (key-'0');
	if(key == 't') {
	  if(companionnr >= 0) {
	    if(parfile.paramset[TEMPO3_PARAM_T0+companionnr-1])
	      toggleFixed(TEMPO3_PARAM_T0+companionnr-1, parfile.fixed, parfile.paramlinked);
	    else
	      toggleFixed(TEMPO3_PARAM_TASC+companionnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(companionnr = 1; companionnr <= TEMPO3_MaxNrCompanions; companionnr++) {
	      if(parfile.par.companion_defined[companionnr-1]) {
		if(parfile.paramset[TEMPO3_PARAM_T0+companionnr-1])
		  toggleFixed(TEMPO3_PARAM_T0+companionnr-1, parfile.fixed, parfile.paramlinked);
		else
		  toggleFixed(TEMPO3_PARAM_TASC+companionnr-1, parfile.fixed, parfile.paramlinked);
	      }
	    }
	    companionnr = -1;
	  }
	}else if(key == 'p') {
	  if(companionnr >= 0) {
	    if(dotimederivs)
	      toggleFixed(TEMPO3_PARAM_PBDOT+companionnr-1, parfile.fixed, parfile.paramlinked);
	    else
	      toggleFixed(TEMPO3_PARAM_PB+companionnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(companionnr = 1; companionnr <= TEMPO3_MaxNrCompanions; companionnr++) {
	      if(parfile.par.companion_defined[companionnr-1]) {
		if(dotimederivs)
		  toggleFixed(TEMPO3_PARAM_PBDOT+companionnr-1, parfile.fixed, parfile.paramlinked);
		else
		  toggleFixed(TEMPO3_PARAM_PB+companionnr-1, parfile.fixed, parfile.paramlinked);
	      }
	    }
	    companionnr = -1;
	  }
	}else if(key == 'a') {
	  if(companionnr >= 0) {
	    if(dotimederivs)
	      toggleFixed(TEMPO3_PARAM_A1DOT+companionnr-1, parfile.fixed, parfile.paramlinked);
	    else
	      toggleFixed(TEMPO3_PARAM_A1+companionnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(companionnr = 1; companionnr <= TEMPO3_MaxNrCompanions; companionnr++) {
	      if(parfile.par.companion_defined[companionnr-1]) {
		if(dotimederivs)
		  toggleFixed(TEMPO3_PARAM_A1DOT+companionnr-1, parfile.fixed, parfile.paramlinked);
		else
		  toggleFixed(TEMPO3_PARAM_A1+companionnr-1, parfile.fixed, parfile.paramlinked);
	      }
	    }
	    companionnr = -1;
	  }
	}else if(key == 'o') {
	  if(companionnr >= 0) {
	    if(dotimederivs)
	      toggleFixed(TEMPO3_PARAM_OMDOT+companionnr-1, parfile.fixed, parfile.paramlinked);
	    else
	      toggleFixed(TEMPO3_PARAM_OM+companionnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(companionnr = 1; companionnr <= TEMPO3_MaxNrCompanions; companionnr++) {
	      if(parfile.par.companion_defined[companionnr-1]) {
		if(dotimederivs)
		  toggleFixed(TEMPO3_PARAM_OMDOT+companionnr-1, parfile.fixed, parfile.paramlinked);
		else
		  toggleFixed(TEMPO3_PARAM_OM+companionnr-1, parfile.fixed, parfile.paramlinked);
	      }
	    }
	    companionnr = -1;
	  }
	}else if(key == 'e') {
	  if(companionnr >= 0) {
	    if(parfile.paramset[TEMPO3_PARAM_ECC +companionnr-1]) {
	      toggleFixed(TEMPO3_PARAM_ECC +companionnr-1, parfile.fixed, parfile.paramlinked);
	    }else {
	      toggleFixed(TEMPO3_PARAM_EPS1+companionnr-1, parfile.fixed, parfile.paramlinked);
	      toggleFixed(TEMPO3_PARAM_EPS2+companionnr-1, parfile.fixed, parfile.paramlinked);
	    }
	  }else {
	    for(companionnr = 1; companionnr <= TEMPO3_MaxNrCompanions; companionnr++) {
	      if(parfile.par.companion_defined[companionnr-1]) {
		if(parfile.paramset[TEMPO3_PARAM_ECC +companionnr-1]) {
		  toggleFixed(TEMPO3_PARAM_ECC +companionnr-1, parfile.fixed, parfile.paramlinked);
		}else {
		  toggleFixed(TEMPO3_PARAM_EPS1+companionnr-1, parfile.fixed, parfile.paramlinked);
		  toggleFixed(TEMPO3_PARAM_EPS2+companionnr-1, parfile.fixed, parfile.paramlinked);
		}
	      }
	    }
	    companionnr = -1;
	  }
	}else if(key == 'g') {
	  if(companionnr >= 0) {
	    toggleFixed(TEMPO3_PARAM_GAMMA+companionnr-1, parfile.fixed, parfile.paramlinked);
	  }else {
	    for(companionnr = 1; companionnr <= TEMPO3_MaxNrCompanions; companionnr++) {
	      if(parfile.par.companion_defined[companionnr-1])
		toggleFixed(TEMPO3_PARAM_GAMMA+companionnr-1, parfile.fixed, parfile.paramlinked);
	    }
	    companionnr = -1;
	  }
	}
      }
    }else if(param == 'p') {
      do_show_params = 1;
      param = 1000;
      printf("Commands:\n");
      printf("  r = RAJ\n");
      printf("  d = DECJ\n");
      printf("  R = PMRA\n");
      printf("  D = PMDEC\n");
      if(macroname[0] != 0) {
	int ret;
	do {
	  ret = fscanf(fmacro, "%c", &ch);
	  //	    printf("XXXXX %d '%c'\n", ch, ch);
	}while(ret == 1 && (ch == '\n' || ch == '\r'));
      }else {
	ppgcurs(&xpos, &ypos, &ch);
      }
      key = ch;
      if(key == 'r') {
	toggleFixed(TEMPO3_PARAM_RAJ, parfile.fixed, parfile.paramlinked);
      }else if(key == 'd') {
	toggleFixed(TEMPO3_PARAM_DECJ, parfile.fixed, parfile.paramlinked);
      }else if(key == 'R') {
	toggleFixed(TEMPO3_PARAM_PMRA, parfile.fixed, parfile.paramlinked);
      }else if(key == 'D') {
	toggleFixed(TEMPO3_PARAM_PMDEC, parfile.fixed, parfile.paramlinked);
      }
    }else if(param == 'v') {
      do_show_params = 1;
      param = 1000;
      printf("Value nr:\n");
      int valueint;
      if(macroname[0] != 0) {
	fscanf(fmacro, "%d", &valueint);
	printf("%d\n", valueint);
	ch = '0' + valueint;
      }else {
	ppgcurs(&xpos, &ypos, &ch);
      }
      key = ch;
      param += (key-'0');
      if(key > '0' && key <= '9') {
	toggleFixed(TEMPO3_PARAM_VALUE+(key-'0')-1, parfile.fixed, parfile.paramlinked);
      }
    }else if(param >= '0' && param <= '9') {    // F0 - F9
      do_show_params = 1;
      toggleFixed(TEMPO3_PARAM_F0+param-'0', parfile.fixed, parfile.paramlinked);
      /*
    }else if(param == ')') {    // Something beyond F9
      printf("Type in terminal the derivitive number to toggle, i.e. 10 for F10:\n");
      int derivnr = -1;
      if(macroname[0] != 0) {
	fscanf(fmacro, "%d", &derivnr);
	//	    printf("XXXXX %d '%c'\n", ch, ch);
      }else {
	scanf("%d", &derivnr);
      }
      //      printf("XXXXXXX deriv number = %d\n", derivnr);
      if(derivnr < 0 || derivnr > TEMPO3_MaxNrFderivatives) {
	printerror(0, "Invalid derivative number provided by user.\n");
      }else {
	do_show_params = 1;
	toggleFixed(TEMPO3_PARAM_F0+derivnr, parfile.fixed, parfile.paramlinked);
      }
      */
    }else if(param == 'm') {
      if(parfile.par.mode == 0)
	parfile.par.mode = 1;
      else
	parfile.par.mode = 0;
      do_show_params = 1;
    }else if(param == 'j' && parfile.par.nrjumps > 0) {
      do_show_params = 1;

      if(parfile.par.nrjumps > 1) {
	printf("Jump nr:\n");
	fflush(stdout);
	if(macroname[0] != 0) {
	  fscanf(fmacro, "%d", &key);
	  key += '0';
	}else {
	  ppgcurs(&xpos, &ypos, &ch);
	  key = ch;
	}
      }else
	key = '1';
      key -= '1';
      printf("Jump parameter number = %d\n", key+1);
      if(key >= 0 && key < TEMPO3_MaxNrJumps) {
	toggleFixed(TEMPO3_PARAM_JUMPS+key, parfile.fixed, parfile.paramlinked);
      }
    }else if(param == 'd') {
      do_show_params = 1;
      toggleFixed(TEMPO3_PARAM_DM, parfile.fixed, parfile.paramlinked);
      /*
    }else if(param == 'N') {
      do_show_params = 1;
      toggleFixed(TEMPO3_PARAM_NBRAKE, parfile.fixed, parfile.paramlinked);
      */
    }else if(param == 'w') {
      do_show_params = 1;
      printf("How many waves? (type in pgplot window, type 'm' for more than 9 waves)\n");
      if(macroname[0] != 0) {
	fscanf(fmacro, "%c", &ch);
      }else {
	ppgcurs(&xpos, &ypos, &ch);
      }
      nrwaves_old = parfile.par.nrwaves;
      if((ch >= '1' && ch <= '9') || ch == 'm') {
	if(ch == 'm') {
	  printf("How many waves (type in xterm): ");
	  fflush(stdout);
	  scanf("%d", &i);
	  ch = '0'+i;
	}
	if(ch - '0' > TEMPO3_MaxNrWaves) {
	  printwarning(0, "WARNING: Nr or wave parameters is truncated at %d. Use the -maxwavenr option to increase this limit.", TEMPO3_MaxNrWaves);
	  ch = '0' + TEMPO3_MaxNrWaves;
	}
	for(i = 0; i <= ch - '1'; i++) {
	  toggleFixed(TEMPO3_PARAM_WAVESIN+i, parfile.fixed, parfile.paramlinked);
	  if(parfile.paramset[TEMPO3_PARAM_WAVESIN+i] == 0) {
	    parfile.par.parameter[TEMPO3_PARAM_WAVESIN+i] = 0;
	    parfile.par.parameter[TEMPO3_PARAM_WAVECOS+i] = 0;
	    parfile.paramset[TEMPO3_PARAM_WAVESIN+i] = 1;
	    parfile.paramset[TEMPO3_PARAM_WAVECOS+i] = 1;
	    if(parfile.par.nrwaves < i + 1)
	      parfile.par.nrwaves = i + 1;
	  }
	}
	/* If nr waves changed, wave_om changes hence should find new solution */
	if(parfile.par.nrwaves != nrwaves_old) {
	  for(i = 0; i < TEMPO3_MaxNrWaves; i++) {
	    parfile.par.parameter[TEMPO3_PARAM_WAVESIN+i] = 0;
	    parfile.par.parameter[TEMPO3_PARAM_WAVECOS+i] = 0;	    
	  }
	}
	parfile.paramset[TEMPO3_PARAM_WAVEOM] = 1;
	parfile.paramset[TEMPO3_PARAM_WAVEEPOCH] = 1;
	parfile.par.parameter[TEMPO3_PARAM_WAVEOM] = 2.0*M_PI/(toas.data_mjdmax-toas.data_mjdmin)/
	    (1.0+4.0/(long double)(parfile.par.nrwaves));
      }
    }else if(param == 'a') {
      do_show_params = 1;
      for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	if(parfile.paramset[i] && ephemeris_check_if_parameter_is_linked(i, parfile.paramlinked) == 0)
	  parfile.fixed[i] = 0;
      }
      parfile.fixed[TEMPO3_PARAM_PSRJNAME] = 1; 
      parfile.fixed[TEMPO3_PARAM_WAVEOM] = 1; 
      parfile.fixed[TEMPO3_PARAM_PEPOCH] = 1; 
      parfile.fixed[TEMPO3_PARAM_DMEPOCH] = 1; 
      parfile.fixed[TEMPO3_PARAM_WAVEEPOCH] = 1; 
      parfile.fixed[TEMPO3_PARAM_TRES] = 1; 
      parfile.fixed[TEMPO3_PARAM_BINARY] = 1;
      parfile.fixed[TEMPO3_PARAM_MODE] = 1;
      parfile.fixed[TEMPO3_PARAM_TRACK] = 1;
      parfile.fixed[TEMPO3_PARAM_TZRMJD] = 1; 
      parfile.fixed[TEMPO3_PARAM_TZRFRQ] = 1;
      parfile.fixed[TEMPO3_PARAM_TZRSITE] = 1;
      parfile.fixed[TEMPO3_PARAM_POSEPOCH] = 1; 
      /* When pressing 'a', do not fit for glitch epochs */
      for(i = 0; i < parfile.par.nrglitches; i++) {
	if(parfile.paramset[TEMPO3_PARAM_GLEP+i]) {
	  parfile.fixed[TEMPO3_PARAM_GLEP+i] = 1;
	}
	parfile.fixed[TEMPO3_PARAM_GLEPMIN+i] = 1;
	parfile.fixed[TEMPO3_PARAM_GLEPMAX+i] = 1;
      }
      if(parfile.par.nrglitches > 1) {
	for(i = 1; i < parfile.par.nrglitches; i++) {
	  if(parfile.paramset[TEMPO3_PARAM_GLTD+i]) {
	    parfile.fixed[TEMPO3_PARAM_GLTD + i] = 0;
	  }
	}
      }
    }else if(param == 'n') {
      do_show_params = 1;
      for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	parfile.fixed[i] = 1;
      }
      parfile.fixed[TEMPO3_PARAM_PHASE0] = 0;
    }else if(param == 65 || param == 88 || param == 'x' || param == 16) {              /* Mouse click (left = show toa, right = delete) ctrl-p to plot*/
      ibest = 0;
      j = 0;
      printf("You clicked on MJD=%f Resid=%f\n", xpos, ypos);
      for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
	if(toas.deleted[i] == 0) {
	  /*	  deltax = fabs(xpos-toas.mjdx[i]); */
	  if(plottype == PLOTTYPE_FREQUENCY) {
	    deltax = (xpos-toas.freqx[i])*(xpos-toas.freqx[i])/((cur_x2-cur_x1)*(cur_x2-cur_x1)) + (ypos-ygraph[i])*(ypos-ygraph[i])/((cur_y2-cur_y1)*(cur_y2-cur_y1));
	  }else {
	    deltax = (xpos-toas.mjdx[i])*(xpos-toas.mjdx[i])/((cur_x2-cur_x1)*(cur_x2-cur_x1)) + (ypos-ygraph[i])*(ypos-ygraph[i])/((cur_y2-cur_y1)*(cur_y2-cur_y1));
	  }
	  if(j == 0 || deltax < deltax_best) {
	    deltax_best = deltax;
	    ibest = i;
	    j = 1;
	  }
	}
      }
      if(param == 65 || param == 16) {
	printf("------------------------------\n");
	printf("TOA = %d (name=%s site='%s' flags='%s')\n", ibest+1, toas.filenames[ibest], toas.site[ibest], toas.flags[ibest]);
	printf("FREQ = %Lf (%Lf at SSB) MHz\n", toas.freqSite[ibest], toas.freqSSB[ibest]);
	printf("MJD = %f (BAT) Error = %Lf (micro-phase)\n", toas.mjdx[ibest], toas.err[ibest]*parfile.par.parameter[TEMPO3_PARAM_F0]);
	mjd2date(toas.mjd[ibest], &year, &month, &day, &hours, &minutes, &seconds);
	printf("UT  = %d %02d %02d %02d:%02d:%02.0f (BAT)\n", year, month, day, hours, minutes, seconds);
	printf("SSB = (%Lf,%Lf,%Lf)\n", toas.ssb1[ibest], toas.ssb2[ibest], toas.ssb3[ibest]);
	printf("Derived turn number: %lld\n", residuals.nturn[ibest]);   // Used to be %Ld
	char textBuffer[1000];
	long double ra, dec, constant;
	ra = parfile.par.raj+parfile.par.parameter[TEMPO3_PARAM_RAJ];
	dec = parfile.par.decj+parfile.par.parameter[TEMPO3_PARAM_DECJ];
	// pmra [mas/yr] * dt [yr] * (1e-3 / 3600.0) [deg/mas]
	constant = ((toas.mjd[ibest] - parfile.par.parameter[TEMPO3_PARAM_POSEPOCH])/365.25 )*(1e-3/3600.0)*(M_PI/180.0);
      	ra  += parfile.par.pmra  * constant;
      	dec += parfile.par.pmdec * constant;
      	ra  += parfile.par.parameter[TEMPO3_PARAM_PMRA]  * constant;
      	dec += parfile.par.parameter[TEMPO3_PARAM_PMDEC] * constant;
	converthms_string(textBuffer, ra*12.0/M_PI, 20, 1);
	printf("Derived RAJ: %s\n", textBuffer);
	converthms_string(textBuffer, dec*180.0/M_PI, 20, 1);
	printf("Derived DECJ: %s\n", textBuffer);
	printf("------------------------------\n");
	if(param == 16) {
	  sprintf(txt, "pplot -TSCR %s", toas.filenames[ibest]);
	  system(txt);
	}
      }else if(param == 88 || param == 'x') {
	printf("------------------------------\nDeleting following TOA:\n");
	printf("TOA = %d (name=%s site='%s' flags='%s')\n", ibest+1, toas.filenames[ibest], toas.site[ibest], toas.flags[ibest]);
	printf("FREQ = %Lf (%Lf at SSB) MHz\n", toas.freqSite[ibest], toas.freqSSB[ibest]);
	printf("MJD = %f (BAT) Error = %Lf (micro-phase)\n", toas.mjdx[ibest], toas.err[ibest]*parfile.par.parameter[TEMPO3_PARAM_F0]);
	mjd2date(toas.mjd[ibest], &year, &month, &day, &hours, &minutes, &seconds);
	printf("UT  = %d %02d %02d %02d:%02d:%02.0f (BAT)\n", year, month, day, hours, minutes, seconds);
	printf("------------------------------\n");
	if(toas.nrtoas > 1) {
	  storeInStack(stackptr, &(parfile.par), toas, parfile.wraps, previous_xstart, previous_xend, parfile.fixed, parfile.paramset);
	  toas.deleted[ibest] = 1;
	  toas.nrtoasdeleted += 1;
	  toas.nrtoas -= 1;
	}
	/* Recalculate data range */
	j = 0;
	for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
	  if(toas.deleted[i] == 0) {
	    if(j == 0) {
	      toas.data_mjdmin = toas.mjd[i];
	      toas.data_mjdmax = toas.mjd[i];
	      j = 1;
	    }else {
	      if(toas.mjd[i] > toas.data_mjdmax)
		toas.data_mjdmax = toas.mjd[i];
	      if(toas.mjd[i] < toas.data_mjdmin)
		toas.data_mjdmin = toas.mjd[i];
	    }
	  }
	}
      }
    }else if(param == 24) {              /* ctrl-x, delete all points */
      storeInStack(stackptr, &(parfile.par), toas, parfile.wraps, previous_xstart, previous_xend, parfile.fixed, parfile.paramset);
      for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
	if(plottype != PLOTTYPE_FREQUENCY) {
	  if((toas.mjdx[i] >= cur_x1) && (toas.mjdx[i] <= cur_x2)) {
	    if(toas.deleted[i] == 0 && toas.nrtoas > 1) {
	      toas.deleted[i] = 1;
	      toas.nrtoasdeleted += 1;
	      toas.nrtoas -= 1;
	    }
	  }
	}else {
	  if((toas.freqx[i] >= cur_x1) && (toas.freqx[i] <= cur_x2)) {
	    if(toas.deleted[i] == 0 && toas.nrtoas > 1) {
	      toas.deleted[i] = 1;
	      toas.nrtoasdeleted += 1;
	      toas.nrtoas -= 1;
	    }
	  }
	}
      }
      /* Recalculate data range */
      j = 0;
      for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
	if(toas.deleted[i] == 0) {
	  if(j == 0) {
	    toas.data_mjdmin = toas.mjd[i];
	    toas.data_mjdmax = toas.mjd[i];
	    j = 1;
	  }else {
	    if(toas.mjd[i] > toas.data_mjdmax)
	      toas.data_mjdmax = toas.mjd[i];
	    if(toas.mjd[i] < toas.data_mjdmin)
	      toas.data_mjdmin = toas.mjd[i];
	  }
	}
      }
#ifdef EnablePSRSALSAdevelop
    }else if(param == 'o' || param == 'O') {
      ppgend();
      if(param == 'O') {
	printf("Output name script (type in xterm): ");
	fflush(stdout);
	if(macroname[0] != 0) {
	  fscanf(fmacro, "%s", txt);
	}else {
	  scanf("%s", txt);
	}
	if(pgopenoutputfile(txt) != 1) {
	  printerror(0, "Cannot open file\n");
	  return 0;
	}
      }
      printf("Output pgplot device (type in xterm): ");
      fflush(stdout);
      if(macroname[0] != 0) {
	fscanf(fmacro, "%s", txt);
      }else {
	scanf("%s", txt);
      }
      ppgopen(txt);
      pgplot_setWindowsize(1000, 500, -1);
      if(doplot("?", ygraph, ygraph2, 1, plottype, showdata, showWaves, toas, &residuals, &parfile, additional_funk_info.superverbose, 0, plottitle, &cur_x1, &cur_y1, &cur_x2, &cur_y2, 0, 0, toas_fiducial, residuals_fiducial, showfiducial, showpepoch, highlightedFlag, nummericalApproxNuNudot, title, amoeba_algorithm, markersize, stridefit_nrdays, stridefit_mintoarange, stridefit_minstride, stridefit_minnroftoas, stridefit_includeF2, stridefit_includeF3, stridefit_scriptfilename, nosort, numthreads, verbose_state) == 0)
	return 0;
      ppgend();
      if(param == 'O') {
	pgcloseoutputfile();
      }
#ifdef HAVEPGPLOT2 
      ppgopen("/pw");
#else
      ppgopen("/xs");
#endif
#endif
    }else if(param == 'S') {
      printf("Available plot types:\n");
      printf("  1 = residual versus mjd (default)\n");
      printf("  2 = spin freq versus mjd\n");
      printf("  3 = spin freq (minus F0 contribution) versus mjd\n");
      printf("  4 = spin freq (minus F0 & F1 contribution) versus mjd\n");
      printf("  5 = spin freq (minus F0, F1 & F2 contribution) versus mjd\n");
      printf("  6 = spin freq derivative versus mjd\n");
      printf("  7 = spin freq derivative (minus F1 contribution) versus mjd\n");
      printf("  8 = spin freq derivative (minus F1 & F2 contribution) versus mjd\n");
      printf("  9 = seconde spin freq derivative versus mjd\n");
      printf("  f = residual versus sky frequency\n");
      printf("  n = braking index versus mjd\n");
      printf("  o = residual versus \"orbital phase\" for first companion\n");
      printf("  O = Various types of orbital plots\n");
      printf("\n  Type choice in pgplot window \n");
      if(macroname[0] == 0) {
	ppgcurs(&xpos, &ypos, &ch);
      }else {
	do {
	    ch = fgetc(fmacro);
	  }while(ch == '\n' || ch == '\r');
      }
      key = ch;
      if(key >= '1' && key <= '9') {
	plottype = key - '1';
      }else if(key == 'n') {
	plottype = PLOTTYPE_NBREAK;
      }else if(key == 'f') {
	plottype = PLOTTYPE_FREQUENCY;
	j = 0;
	for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
	  if(toas.deleted[i] == 0) {
	    if(j == 0) {
	      residuals.freqmin = toas.freqSSB[i];
	      residuals.freqmax = toas.freqSSB[i];
	      j = 1;
	    }else {
	      if(toas.freqSSB[i] > residuals.freqmax)
		residuals.freqmax = toas.freqSSB[i];
	      if(toas.freqSSB[i] < residuals.freqmin)
		residuals.freqmin = toas.freqSSB[i];
	    }
	  }
	}
      }else if(key == 'o') {
	plottype = PLOTTYPE_ORBITALPHASE;
	printf("Making plot of residual versus \"orbital phase\" = (t-t0)/pb. Variations in the pb etc are ignored.\n");
      }else if(key == 'O') {
	long double mjd_orbit, inclination, tmin, tmax;
	printf("Specify MJD of the pulsar/companion position in orbit plots: ");
	fflush(stdout);
	if(macroname[0] != 0) {
	  fscanf(fmacro, "%Lf", &mjd_orbit);
	}else {
	  scanf("%Lf", &mjd_orbit);
	}
	printf("Specify inclination of the orbits (degrees): ");
	fflush(stdout);
	if(macroname[0] != 0) {
	  fscanf(fmacro, "%Lf", &inclination);
	}else {
	  scanf("%Lf", &inclination);
	}
	printf("Specify start time of orbits (mjd, -1 is full orbit, -2 is full toa range): ");
	fflush(stdout);
	if(macroname[0] != 0) {
	  fscanf(fmacro, "%Lf", &tmin);
	}else {
	  scanf("%Lf", &tmin);
	}
	if(tmin < 0) {
	  tmax = tmin;
	}else {
	  printf("Specify end time of orbits (mjd, -1 is full orbit, -2 is full toa range): ");
	  fflush(stdout);
	  if(macroname[0] != 0) {
	    fscanf(fmacro, "%Lf", &tmax);
	  }else {
	    scanf("%Lf", &tmax);
	  }
	}
	printf("Output pgplot device (type in xterm): ");
	fflush(stdout);
	if(macroname[0] != 0) {
	  fscanf(fmacro, "%s", txt);
	}else {
	  scanf("%s", txt);
	}
	printf("If you don't want the data from the graphs to be written to a file, then press return. Otherwise, specify the prefix for the generated files: ");
	fflush(stdout);
	if(macroname[0] != 0) {
	  fgets(orbit_dump_prefix, 1000, fmacro);
	}else {
	  fgets(orbit_dump_prefix, 1000, stdin);
	}
	/* Ignore return of last input. Might be system dependent if this is necessary. */
	if(orbit_dump_prefix[0] == 10) {
	  if(macroname[0] != 0) {
	    fgets(orbit_dump_prefix, 1000, fmacro);
	  }else {
	    fgets(orbit_dump_prefix, 1000, stdin);
	  }
	}
	if(orbit_dump_prefix[strlen(orbit_dump_prefix)-1] == 10)
	  orbit_dump_prefix[strlen(orbit_dump_prefix)-1] = 0;

	plotOrbit(&(parfile.par), parfile.paramset, &toas, mjd_orbit, tmin, tmax, inclination, 1.35, hierarchical, markersize, txt, orbit_dump_prefix, finderrors, parfile.par.dplus);
      }else {
	printf("Unknown option ('%c', %d)\n", key, key);
      }
      /*    }else if(param == 'l') { */
      /*      void mrq_func_internal_ld(long double x, long double *a, long double *y, long double *dyda, int ma);

      fill_xstart(&parfile, xstart);
      mrqmin_ld(toas.mjd, toas.); */

    }else if(param == 'r') {
      /* Pop stack */
      popstack(stackptr, &(parfile.par), &toas, parfile.wraps, &previous_xstart, &previous_xend, parfile.fixed, parfile.paramset);
      do_show_params = 1;
    }else if(param == 'F') {
      if(toas_fiducial.nrtoas == 1) {
      	calcphases(&(parfile.par), toas_fiducial, showWaves, &residuals_fiducial, parfile.paramlinked, parfile.paramset, 0, numthreads);
	subtractTurns(toas_fiducial, &residuals_fiducial, parfile.wraps, nosort, 0);
	ypos = residuals_fiducial.resid[0];
	sprintf(parfile.par.tzrsite, "@");
	parfile.par.parameter[TEMPO3_PARAM_TZRFRQ] = 0;
	parfile.par.parameter[TEMPO3_PARAM_TZRMJD] = toas_fiducial.mjd[0]-(long double)ypos/(parfile.par.parameter[TEMPO3_PARAM_F0]*24.0*3600.0);
	toas_fiducial.mjd[0] = parfile.par.parameter[TEMPO3_PARAM_TZRMJD];
      }
      do_show_params = 1;
    }else if(param == 0) {  /* ctrl-space */
      superfit = nrSuperFit;
    }else if(param == ' ') {
      storeInStack(stackptr, &(parfile.par), toas, parfile.wraps, previous_xstart, previous_xend, parfile.fixed, parfile.paramset);

      /* Remember what x-range was of last fit */
      previous_xstart = residuals.mjdmin;
      previous_xend   = residuals.mjdmax;

      do_show_params = 1;
      if(autoquiet == 0) {
	fprintf(stderr, "Determining start parameters downhill simplex method\n");
      }
      fill_xstart(&(parfile.par), xstart);
      struct timeval starttime, endtime;
      gettimeofday(&starttime,0x0);
      long double *dx_ptr = NULL;    // The array dx is not always calculated, hence available
      //      if(finderrors || randstart || fit_algorithm != 1) {
      if(randstart || fit_algorithm != 1) {
	fill_dx(&(parfile.par), toas, residuals, parfile.wraps, parfile.fixed, dx, additional_funk_info.superverbose, parfile.paramlinked, parfile.paramset, nosort, numthreads);
	dx_ptr = dx;
      }
      gettimeofday(&endtime,0x0);
      /*      printf("Finding start offsets took %f seconds\n", timediff(starttime, endtime)); */
      additional_funk_info.toas = &toas;
      additional_funk_info.residuals = &residuals;
      additional_funk_info.parfile = &(parfile.par);
      additional_funk_info.wraps = parfile.wraps;
      additional_funk_info.paramlinked = parfile.paramlinked;
      additional_funk_info.paramset = parfile.paramset;
      additional_funk_info.numthreads = numthreads;
      additional_funk_info.nosort = nosort;
     
      /*      if(finderrors == 2) {
	finderrors = 1;
      }else if(finderrors == 1) {
	finderrors = 0;
	additional_funk_info.superverbose = 0;
	}*/
      if(finderrors) {
	for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	  parfile.par.dplus[i] = 0;
	  parfile.par.dmin[i] = 0;
	}
	/*	printf("Measuring\n"); */
      }
      
      if(randstart) {
	if(randstart_itt == randstart_ittmax) {
	  long nrfitparams = 0, perturbparam;
	  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	    if(parfile.fixed[i] == 0) {
	      nrfitparams++;
	    }
	  }
	  for(randstart_itt = 0; randstart_itt < randstart_nrparams; randstart_itt++) { // Change 2 parameters
	    perturbparam = gsl_rng_uniform_int(rand_num_gen, nrfitparams);  // Get an int in range 0 to nrfitparams-1
	    long counter = 0;
	    for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	      if(parfile.fixed[i] == 0) {
		if(counter == perturbparam)
		  break;
		counter++;
	      }
	    }
	    long double offset;
	    offset = (2.0*randstart_sigma*gsl_rng_uniform(rand_num_gen)-randstart_sigma)*dx_ptr[i];
	    xstart[i] += offset;  // Offset start positions by some significant amount
	    printf("APPLYING OFFSET ON PARAMETER %s: %Le\n", tempo3_parfile_descr.identifiers[i][0], offset);
	  }
	  randstart_itt = 0;
	}else {
	  randstart_itt++;
	}
      }
      /* Use covariance matrix below to estimate errors */
      ok = 0;
      additional_funk_info.best_chi2 = -1;      
      if(autoquiet == 0)
	fprintf(stderr, "Starting fitting\n");
      if(fit_algorithm == 1) {
	funk_levmar_info info;
	info.fixed = parfile.fixed;
	info.nrfitparameters = 0;
	for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	  if(parfile.fixed[i] == 0) {
	    info.nrfitparameters++;
	  }
	}
	i = fit_levmar2_ld(xstart, parfile.fixed, xfit, &rms, TEMPO3_PARAM_TOTAL, funk_levmar, &info, ftol, &nfunk, 0, 1.0, NULL, NULL, verbose_state);
      }else {
	i = doAmoeba_ld(amoeba_algorithm, xstart, dx_ptr, parfile.fixed, xfit, &rms, TEMPO3_PARAM_TOTAL, &funk, ftol, &nfunk, 0, 0*finderrors, 1.0, parfile.par.dplus, parfile.par.dmin);
      }
      if(i != 0) {
	if(autoquiet == 0)
	  printf("\n");
	fflush(stdout);
	if(i == 1 || i == 5) {
	  printf("\n\n\nERROR, did not converge!\n");
	  if(additional_funk_info.best_chi2 <= residuals.chi2) {
	    printf("Found better solution.\n");
	    memcpy(xfit, additional_funk_info.best_xfit, TEMPO3_PARAM_TOTAL*sizeof(long double));
	    use_xfit(&(parfile.par), additional_funk_info.best_xfit);
	    rms = additional_funk_info.best_chi2;
	    ok = 1;
	  }else {
	    printf("Ignoring solution.\n");
	    memcpy(xfit, xstart, TEMPO3_PARAM_TOTAL*sizeof(long double));
	    use_xfit(&(parfile.par), xfit);
	    rms = residuals.chi2;
	    ok = 1;
	  }
	}else if (i == 3) {
	  printf("\n\n\nERROR, there should be at least two fit-parameters!\n\n\n");
	  ok = 0;
	}else {
	  printerror(0, "\nERROR amoeba (%d)\n", i);
	  return 0;
	}
      }else {
	if(autoquiet == 0)
	  printf("\n");
	fflush(stdout);
	if(rms <= residuals.chi2) {
	  //residuals->chi2/(long double)(toas.nrtoas-nrparamsused), toas.nrtoas-nrparamsused);
	  if(parfile.par.mode) {
	    int nrparamsused = 0;
	    for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	      if(parfile.fixed[i] == 0)
		nrparamsused++;
	    }
	    if(autoquiet == 0)
	      printf("Accepting parameters with chi2=%Lf with %d dof. giving a reduced-chi2=%Lf\n", rms, toas.nrtoas-nrparamsused, rms/(long double)(toas.nrtoas-nrparamsused));
	  }else {
	    if(autoquiet == 0)
	      printf("Accepting parameters with rms=%Le phase\n", rms);
	  }
	  use_xfit(&(parfile.par), xfit);
	  //	  printf("nbreak -> %Le (%Le) parameter %d\n", parfile.par.parameter[TEMPO3_PARAM_NBRAKE], xfit[TEMPO3_PARAM_NBRAKE], TEMPO3_PARAM_NBRAKE);
	  /* Make sure all parameters are flagged as being set which were fitted for */
	  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	    if(parfile.fixed[i] == 0) {
	      parfile.paramset[i] = 1;
	      /*	      printf("Using %d\n", i);*/
	    }
	  }
	  ok = 1;
	}else {
	  ok = 0;
	  printf("\nIgnoring solution solution (%Lf < %Lf).\n", rms,  residuals.chi2);
	}
	/*	printf("Press a key to continue\n");
	  pgetch(); */
      }
      /* Removed  && ok, should always work */
      if(finderrors && ok != 0) {
	clock_t tstart, tend;
	tstart = 0;   // Should be uneccesary to do this initialisation, but solves a compiler warning
	if(verbose_state.debug) {
	  tstart = clock();
	}
	if(autoquiet == 0 && verbose_state.debug == 0)
	  fprintf(stderr, "Determining error bars\n");

	finderrors_tempo3(parfile.par.dplus, parfile.par.dmin, &(parfile.par), toas, residuals, parfile.fixed, parfile.paramset, parfile.paramlinked, xstart, dx_ptr, xfit, mrq_alpha, mrq_alpha_inv, max_nr_fitted_params, mrq_deriv, dumpcovar, verbose_state);

	double ttotal;
	tend = clock();
	ttotal = (tend-tstart)/(double)CLOCKS_PER_SEC;
	if(verbose_state.debug) {
	  printf("Finding errorbars took %.2lf seconds = %.2lf minutes\n", ttotal, ttotal/60.0);
	}
      }
      gettimeofday(&endtime,0x0);
      if(repeatable == 0)
	printf("Finding solution took %f seconds\n", timediff(starttime, endtime));
    }else if(param == 13) {
      storeInStack(stackptr, &(parfile.par), toas, parfile.wraps, previous_xstart, previous_xend, parfile.fixed, parfile.paramset);
      /* Remember what x-range was of last fit */
      previous_xstart = residuals.mjdmin;
      previous_xend   = residuals.mjdmax;

      if(autoquiet == 0)
	fprintf(stderr, "Start fitting\n");
      do_show_params = 1;
      for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	if(parfile.fixed[i] == 0) {
	  dofit(i, &(parfile.par), toas, &residuals, parfile.wraps, parfile.paramlinked, parfile.paramset, nosort, numthreads);
	}
      }
    }else {
      printf("Unknown command: %c (%d).\n", param, param);
    }
  }while(param != 'q' && param != 27 && user_requested_quit == 0);

  // Stop processing gtk events in wakeup alarm call, as the tempo3
  // would occasionally freeze after quiting. This is possibly because
  // gtk_parfile gets released? But seems to happen even if -nogtk is
  // enabled.
  gtk_iteration_blocked++;

  ppgend();
  free(residuals.residy);
  free(residuals.resid);
  free(residuals.nturn); 
  free(residuals.nturn_fixed);
  free_toas(&toas);
  //  free(toas.mjd);
  //  free(toas.ssb1);
  //  free(toas.ssb2);
  //  free(toas.ssb3);
  //  free(toas.err);
  //  free(toas.freq);
  //  for(i = 0; i < toas.nrtoas; i++) {
  //    free(toas.filenames[i]);
  //    free(toas.flags[i]);
  //    free(toas.site[i]);
  //  }
  //  free(toas.filenames);
  //  free(toas.flags);
  //  free(toas.site);
  //  free(toas.deleted);
  //  free(toas.mjdx);
  //  free(toas.freqx);
  free(ygraph);
  free(ygraph2);
  //  free(toas.nu);
  //  free(toas.nudot); 
  //  free(toas.nu_mjd);
  if(loadtimfile_flag) {
    if(runtempo2) {
      free(indx_input_toas);
    }
  }
 
  //  free(mrq_covar); 
  free(mrq_alpha); 
  free(mrq_alpha_inv); 
  free(xstart);
  free(dx);
  free(xfit);
  //  free(dmin);
  free(mrq_deriv);
  //  free(fixed);
  //  free(paramset);
  //  free(paramlinked);
  free(additional_funk_info.best_xfit);
  for(i = 0; i < MaxNrStack; i++) {
    free(stackptr->deleted[i]);
    free_ephemeris_par_only(&(stackptr->parfiles[i]));
    free(stackptr->fixed[i]);
    free(stackptr->paramset[i]);
  }
  free(stackptr);

  if(macroname[0] != 0) {
    fclose(fmacro);
  }

#ifdef TEMPO3_EnableGTK
  free_ephemeris_par_only(&gtk_parfile);
#endif
  //  free_ephemeris_par_only(&parfile);
  //  free_ephemeris_par_only(&parfile_old);
  free_ephemeris_par_only(&(additional_funk_info.parfile_funk));
  if(toas_fiducial.nrtoas > 0) {
    free_toas(&toas_fiducial);
    free_residuals(&residuals_fiducial);
  }

  gsl_rng_free(rand_num_gen);
  free_ephemeris(&parfile, verbose_state);
  free_ephemeris(&parfile_old, verbose_state);
  if(cleanup_tempo3_lib(verbose_state) == 0) {
    printerror(verbose_state.debug, "ERROR tempo3: Releasing tempo3 library related memory failed.");
    return 0;
  }
  return 0;
}             // End of main




// spinfreq used to conert TOA errors in time to phase
void calcChi2(residual_def *residuals, toa_data_def toas, long double *variance, long double spinfreq, long double *chi2w)
{
  int i;
  long double chi;
  *variance = 0;
  *chi2w = 0;
  for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
    if(toas.deleted[i] == 0) {
      if(toas.mjd[i] >= residuals->mjdmin && toas.mjd[i] <= residuals->mjdmax) {
	chi = residuals->resid[i]; 
	*chi2w += chi*chi/(long double)(1e-12*toas.err[i]*toas.err[i]*spinfreq*spinfreq);
	*variance += chi*chi;
      }
    }
  }
}

typedef struct _thread_data_t {
  int nthreads;  // Number of threads used
  int thread_id; // Counting from 0
  tempo3_parameters_def *parfile;
  toa_data_def *toas;
  int showWaves;
  residual_def *residuals;
  int *paramset;
  int ignorebinary;
}thread_data_calcphases;

/* thread function */
void *calcphases_thr_func(void *arg)
{
  thread_data_calcphases *data = (thread_data_calcphases *)arg;

  int i;

  for(i = 0; i < data->toas->nrtoas + data->toas->nrtoasdeleted; i++) {
    if(i % data->nthreads == data->thread_id) {  // Only calculate residual assigned to this specific thread
      data->residuals->resid[i] = evaluate_ephemeris(data->toas->mjd[i], data->toas->freqSSB[i], data->toas->ssb1[i], data->toas->ssb2[i], data->toas->ssb3[i], data->showWaves, 0, data->toas->site[i], data->toas->flags[i], (data->parfile), 0, data->paramset, 0, data->ignorebinary);
    }
  }
 
  pthread_exit(NULL);
}

// SEE ABOVE FOR THE THREADED EQUIVALENT
void calcphases(tempo3_parameters_def *parfile, toa_data_def toas, int showWaves, residual_def *residuals, int *paramlinked, int *paramset, int ignorebinary, int numthreads)
{
  int i;

  // Do this once, rather for each toa
  ephemeris_setLinkedValues_par_version(parfile, paramlinked);
  paramlinked = NULL;

  if(numthreads == 1) {
    /*  float offset, mjdoff; 
	int n; */
    
    /* Find offset to make first resid appear at zero (more stable for
       phase wrap problems). First toa is put at 0, furthest away from
       wrap. */
    /*  offset = -100;
	mjdoff = 1e10;*/
    for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
      /*    fprintf(stderr, "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxx start2 %p\n", paramlinked); */
      // PATRICK: CAN USE freqSite RATHER THAN freqSSB TO REVERT TO OLD BEHAVIOUR
      residuals->resid[i] = evaluate_ephemeris(toas.mjd[i], toas.freqSSB[i], toas.ssb1[i], toas.ssb2[i], toas.ssb3[i], showWaves, 0, toas.site[i], toas.flags[i], parfile, 0, paramset, 0, ignorebinary);
      /*    fprintf(stderr, "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxx end2\n"); */
      
      //    if(toas.mjdx[i] >= residuals->mjdmin && toas.mjdx[i] <= residuals->mjdmax) {
      /*      if(toas.mjdx[i] < mjdoff) {
	      offset = residuals->resid[i];*/
      /*	printf("OFFSET = %f (%f)\n", offset, toas.mjdx[i]); */
      /*	mjdoff = toas.mjdx[i];
		}*/
      //    }
    }
    
  /* Ignore for now, because offset is a fit parameter */
  /*  
      float av;
int nractivetoas;
offset  = 0;
  av = 0; 
  nractivetoas = 0;
  for(i = 0; i < toas.nrtoas; i++) {
    residuals->resid[i] = evaluate_ephemeris(toas.mjd[i], parfile, showWaves, -offset, 0);
    if(toas.mjdx[i] >= residuals->mjdmin && toas.mjdx[i] <= residuals->mjdmax) {
      av += residuals->resid[i];
      nractivetoas++;
      }*/
    /*    printf("%f %f\n", mjdx[i], resid[i]); */
  /*  }
      av /= (float)nractivetoas;*/
  /* Don't subtract average, it is a fit parameter */
  /*
  for(i = 0; i < toas.nrtoas; i++) {
    residuals->resid[i] -= av;
    }*/
  }else {
    gtk_iteration_blocked++;
    int rc;
    pthread_t *thr;
    thread_data_calcphases *thr_data;
    thr = malloc(numthreads*sizeof(pthread_t));
    thr_data = malloc(numthreads*sizeof(thread_data_calcphases));
    if(thr == NULL || thr_data == NULL) {
      printerror(0, "Memory allocation error while creating threads.");
      exit(0);
    }
    for(i = 0; i < numthreads; i++) {
      thr_data[i].nthreads = numthreads;
      thr_data[i].thread_id = i;
      thr_data[i].parfile = parfile;
      thr_data[i].toas = &toas;
      thr_data[i].showWaves = showWaves;
      thr_data[i].residuals = residuals;
      thr_data[i].paramset = paramset;
      thr_data[i].ignorebinary = ignorebinary;
      if ((rc = pthread_create(&thr[i], NULL, calcphases_thr_func, &thr_data[i]))) {
	fprintf(stderr, "error: pthread_create failed with return code %d\n", rc);
	exit(0);
      }
    }
    // After launching the threads, wait until they are finished before proceeding.
    for (i = 0; i < numthreads; i++) {
      pthread_join(thr[i], NULL);
    }
    free(thr);
    free(thr_data);
    gtk_iteration_blocked--;
  }
}

// If nosort is set, the toa's are not assumed to be in toa order
void subtractTurns(toa_data_def toas, residual_def *residuals, tempo3_wraps_def *wraps, int nosort, int fixnturn)
{
  int i, n, firstpoint;
  long double mjdmin;
  long long nturn, nturn2;
  long double phase;

  if(fixnturn) {  // In this case, the code is more simple/faster, so implement seperately
    for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
      residuals->resid[i] -= residuals->nturn[i];
    }
  }else {
    firstpoint = 1;
    nturn2 = 0;
    mjdmin = 0;
    /* Subtract amount of turns of first point*/
    for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
      if(toas.mjd[i] >= residuals->mjdmin && toas.mjd[i] <= residuals->mjdmax) {  
	if(toas.mjd[i] < mjdmin || firstpoint) {
	  firstpoint = 0;
	  phase = residuals->resid[i];
	  phase += 0.5;
	  nturn = phase;
	  nturn2 = nturn-10;
	  phase -= nturn2;
	  /* Remember how much subtracted, avoid wrap problem. */
	  nturn = phase;
	  nturn2 += nturn;
	  /* Make residual between -0.5 and 0.5 */
	  /*	nturn = residuals->resid[i]+0.5;*/
	  if(nosort == 0) {
	    break;   // No need to search further, all following toa's will be later
	  }
	  mjdmin = toas.mjd[i];
	}
      }
    }
    
    /* Make residuals above 0.5, which is easier later on. */
    nturn2--;
    
    /* Determine nturns of each toa */
    for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
      // 18 sept 2016 - this used to be fixnturn == 0 ||  residuals->nturn_fixed[i] == 0
      // This not implied that if pulse numbering is specified, it is always used, while before fixnturn is not set (after fitting process) the phase can change.
      if(fixnturn == 0 &&  residuals->nturn_fixed[i] == 0) {
	/* Make residuals above 0.5, which is easier later on. */
	phase = residuals->resid[i] - (long double)nturn2 + 0.5;
	/* phase is now at least 1 (so it cannot be negative) because of step above */
	nturn = phase;
	residuals->nturn[i] = nturn + nturn2;
	/* After subtracting nturn, phase should now between 0 and 1, but make residual between -0.5 and 0.5 */
	residuals->resid[i] = phase-nturn-0.5;
      }else {                     /* Only subtract the turns */
	residuals->resid[i] -= residuals->nturn[i];
	/*      if(toas.mjdx[i] >= residuals->mjdmin && toas.mjdx[i] <= residuals->mjdmax)
		printf("XXXX %f\n", residuals->resid[i]);*/
      }
    }
  }

  if(wraps->nrwraps > 0) {
    for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
      for(n = 0; n < wraps->nrwraps; n++) {
	if(toas.mjdx[i] >= wraps->mjd[n]) {
	  residuals->resid[i] += wraps->phase[n];
	}
      }
    }
  }

  // ONLY DONE WHEN USED, TO SPEED UP FITTING
  //  for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
  //    residuals->residy[i] = residuals->resid[i];
  //  }

}





void perturb(long double *value, long double fraction)
{
  long double dx;
  dx = fraction*(gsl_rng_uniform(rand_num_gen)-0.5)*(*value);
  *value += dx;
}


/* Perturb most parameters in the timing solution. All epochs are
   unchanged, as well as F0 and F1. I cannot remember the reason for not perturbing f0 and f1 */
void perturb_parfile(tempo3_parameters_def *parfile, int *paramset, int *paramlinked, long double fraction)
{
  int i;

  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    if(i != TEMPO3_PARAM_F0 && i != paramset[TEMPO3_PARAM_F0 + 1] &&                                    // Ignore F0 & F1
       i != TEMPO3_PARAM_WAVEOM && (i < TEMPO3_PARAM_GLEP || i >= TEMPO3_PARAM_GLEP+TEMPO3_MaxNrGlitches) && i != TEMPO3_PARAM_PEPOCH && i != TEMPO3_PARAM_WAVEEPOCH && (i < TEMPO3_PARAM_GLEPMIN || i >= TEMPO3_PARAM_GLEPMIN+TEMPO3_MaxNrGlitches) && (i < TEMPO3_PARAM_GLEPMAX || i >= TEMPO3_PARAM_GLEPMAX+TEMPO3_MaxNrGlitches) && i != TEMPO3_PARAM_TZRMJD && TEMPO3_PARAM_TZRFRQ && i != TEMPO3_PARAM_TZRSITE && i != TEMPO3_PARAM_POSEPOCH && i != TEMPO3_PARAM_DMEPOCH) {   // And epochs
      if(paramset[i] == 1) {
	if(ephemeris_check_if_parameter_is_linked(i, paramlinked) == 0) {
	  perturb(&(parfile->parameter[i]), fraction);
	}
      }
    }
  }
}

void print_ephemeris_ExtraParams(FILE *stream, tempo3_parameters_def *parfile, int *paramset, int *paramlinked, int showerrors, long double *dplus, long double distance, long double distanceerror, int hierarchical)
{
  int i;
  long double p, perr, pdot, edot, n, tau, bs, blc, rlc, err, pm, pmerr;
  fprintf(stream, "DERIVED PARAMETERS:\n");

  ephemeris_setLinkedValues_par_version(parfile, paramlinked);

  p = 1.0/parfile->parameter[TEMPO3_PARAM_F0];
  fprintf(stream, "  P     = %13Lf s", p);
  if(showerrors) {
    perr = fabsl(p)*sqrtl(powl(dplus[TEMPO3_PARAM_F0]/parfile->parameter[TEMPO3_PARAM_F0], 2.0));
    fprintf(stream, " (%Le)", perr);
  }

  pdot = -parfile->parameter[TEMPO3_PARAM_F0+1]/(parfile->parameter[TEMPO3_PARAM_F0]*parfile->parameter[TEMPO3_PARAM_F0]);
  fprintf(stream, "      Pdot = %13Le s/s", pdot);
  if(showerrors) {
    err = fabsl(pdot)*sqrtl(powl(dplus[TEMPO3_PARAM_F0+1]/parfile->parameter[TEMPO3_PARAM_F0+1], 2.0)+ powl(2.0*dplus[TEMPO3_PARAM_F0]/parfile->parameter[TEMPO3_PARAM_F0], 2.0));
    fprintf(stream, " (%Le)", err);
  }

  n = parfile->parameter[TEMPO3_PARAM_F0]*parfile->parameter[TEMPO3_PARAM_F0+2]/(parfile->parameter[TEMPO3_PARAM_F0+1]*parfile->parameter[TEMPO3_PARAM_F0+1]);
  fprintf(stream, "   n   = %13Lf", n);
  if(showerrors) {
    err = fabsl(n)*sqrtl(powl(dplus[TEMPO3_PARAM_F0]/parfile->parameter[TEMPO3_PARAM_F0], 2.0) + powl(dplus[TEMPO3_PARAM_F0+2]/parfile->parameter[TEMPO3_PARAM_F0+2], 2.0) + powl(2.0*dplus[TEMPO3_PARAM_F0+1]/parfile->parameter[TEMPO3_PARAM_F0+1], 2.0));
    fprintf(stream, " (%Lf)", err);
  }
  fprintf(stream, "\n");


  edot = 3.95e31*(pdot/1e-15)/(p*p*p);
  tau = p/(2.0*pdot);
  tau /= TEMPO3_SecPerDay*365.25*1e3;
  bs = 3.2e19*sqrtl(p*pdot);
  blc = 9.2*powl(p, -5.0/2.0)*sqrtl(pdot/1e-15);
  rlc = 4.77e4*p;
  fprintf(stream, "  Edot  = %13Le erg/s  tau  = %13Lf kyr\n", edot, tau);
  fprintf(stream, "  Bsurf = %13Le G      BLC  = %13Le G     RLC = %13Le km\n", bs, blc, rlc);
  if(paramset[TEMPO3_PARAM_PMRA] && paramset[TEMPO3_PARAM_PMDEC]) {
    long double pmra, pmdec;
    pmra = parfile->pmra + parfile->parameter[TEMPO3_PARAM_PMRA];
    pmdec = parfile->pmdec + parfile->parameter[TEMPO3_PARAM_PMDEC];
    pm = sqrtl(pmra*pmra+pmdec*pmdec);
    fprintf(stream, "  proper motion = %Lf", pm);
    if(showerrors) {
      pmerr = pmra*pmra*dplus[TEMPO3_PARAM_PMRA]*dplus[TEMPO3_PARAM_PMRA]+pmdec*pmdec*dplus[TEMPO3_PARAM_PMDEC]*dplus[TEMPO3_PARAM_PMDEC];
      pmerr = sqrt(pmerr/(pm*pm));
      fprintf(stream, " (%Lf)", pmerr);
    }
    fprintf(stream, " mas/yr\n");
  }

  int nrcompanions = 0;
  for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
    if(parfile->companion_defined[i]) {
      nrcompanions++;
    }
  }
  if(nrcompanions > 1) {
    if(hierarchical) {
      printwarning(0, "WARNING: There are multiple orbits defined and they are treated as a hierarchinal system: all interactions between companions are ignored.");
    }else {
      printwarning(0, "WARNING: There are multiple orbits defined. By default the reported values relating to the companions ignore the presence of the other companions, which will be wrong. It is in approximation only correct for the most inner companion in the limit of the system being an hierarchical system. Consider using the -hierarchical option, if applicable.");
    }
  }

  long double inner_mass90, inner_mass60, inner_mass26;
  long double m_psr = 1.35;
  inner_mass90 = m_psr;
  inner_mass60 = m_psr;
  inner_mass26 = m_psr;
  for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
    if(parfile->companion_defined[i]) {
      if(i == 0) {
	fprintf(stream, "\n  The min/mid/upper companion masses assume m_psr=%.3Lf Msun and correspond to i=90/60/25.842 deg.\n", m_psr);
	fprintf(stream, "  The upper mass is at a 90%% confidence level assuming a random inclination angle.\n");
      }
      if(paramset[TEMPO3_PARAM_A1+i] && paramset[TEMPO3_PARAM_PB+i]) {
	/*
	  P*K_1^3/(2pi*G) = P*(2*pi*a*sin(i)/P)^3/(2pi G) = 4*pi^2(a*sin(i))^3/(G*P^2) = fn
             fn = (m2*sin(i))^3/(m1+m2)**2 < m2
             fn*m1^2 + fn*m2^2 = (m2*sin(i))^3
             fn*m1^2 + fn*m2^2 = (m2*sin(i))^3
	 */
	long double GM          = 1.3271243999e20;      /* Gravitational constant * mass sun */
	long double a1err = 0;
	long double a1 = parfile->parameter[TEMPO3_PARAM_A1+i]*TEMPO3_SPEEDOFLIGHT; // So this is p for a parabolic orbit
	if(showerrors)
	  a1err = dplus[TEMPO3_PARAM_A1+i]*TEMPO3_SPEEDOFLIGHT;
	long double pberr = 0;
	long double pb = parfile->parameter[TEMPO3_PARAM_PB+i]*TEMPO3_SecPerDay;
	if(showerrors)
	  pberr = dplus[TEMPO3_PARAM_PB+i]*TEMPO3_SecPerDay;
	long double fn  = 4*M_PI*M_PI*a1*a1*a1/(GM*pb*pb);	  
	long double fnerr;
	if(showerrors) {
	  fnerr = sqrtl(fn*fn*(9.0*a1err*a1err/(a1*a1) + 4.0*pberr*pberr/(pb*pb)));
	  printf("  Mass function for companion %d: %Le +- %Le Msun = %Le +- %Le Mearth\n", i+1, fn, fnerr, fn*1.989e6/5.972, fnerr*1.989e6/5.972);
	}else {
	  printf("  Mass function for companion %d: %Le Msun = %Le Mearth\n", i+1, fn, fn*1.989e6/5.972);
	}
	// 2nd argument is sin(i)
	long double m2_min = solve_mass_function_m2(fn, 1.0, inner_mass90);
	long double m2_min_err = 0;
	// Probability to see inclination < i
	//   p = 1 - cos(i)
	//   cos(i) = 1 - p
	//   p = 100%   -> cos(i) = 0.0   -> i = 90 deg
	//   p =  50%   -> cos(i) = 0.5   -> i = 60 deg
	//   p =  10%   -> cos(i) = 0.9   -> i = 25.842 deg
	long double m2_med = solve_mass_function_m2(fn, sin(acos(1.0-0.5)), inner_mass60);  // i = 60 deg
	long double m2_med_err = 0;
	long double m2_max = solve_mass_function_m2(fn, sin(acos(1.0-0.1)), inner_mass26);  // i = 26 deg  
	long double m2_max_err = 0;
	if(showerrors) {
	  m2_min_err = m2_min*(inner_mass90+m2_min)*fnerr/(fn*(3.0*inner_mass90+m2_min));
	  m2_med_err = m2_med*(inner_mass60+m2_med)*fnerr/(fn*(3.0*inner_mass60+m2_med));
	  m2_max_err = m2_max*(inner_mass26+m2_max)*fnerr/(fn*(3.0*inner_mass26+m2_max));
	  m2_min_err = sqrtl(m2_min_err*m2_min_err);
	  m2_med_err = sqrtl(m2_med_err*m2_med_err);
	  m2_max_err = sqrtl(m2_max_err*m2_max_err);
	}
	if(showerrors) {
	  printf("    min/median/max mass: %Lf +- %Lf / %Lf +- %Lf / %Lf +- %Lf Msun\n", m2_min, m2_min_err, m2_med, m2_med_err, m2_max, m2_max_err);
	  printf("    min/median/max mass: %Lf +- %Lf / %Lf +- %Lf / %Lf +- %Lf Mearth\n", m2_min*1.989e6/5.972, m2_min_err*1.989e6/5.972, m2_med*1.989e6/5.972, m2_med_err*1.989e6/5.972, m2_max*1.989e6/5.972, m2_max_err*1.989e6/5.972);
	}else {
	  printf("    min/median/max mass: %Lf/%Lf/%Lf Msun\n", m2_min, m2_med, m2_max);
	  printf("    min/median/max mass: %Lf/%Lf/%Lf Mearth\n", m2_min*1.989e6/5.972, m2_med*1.989e6/5.972, m2_max*1.989e6/5.972);
	}
	if(hierarchical) {
	  inner_mass90 += m2_min;
	  inner_mass60 += m2_med;
	  inner_mass26 += m2_max;
	}
      }
    }
  }

  if(distance > 0) {
    fprintf(stream, "\nDISTANT DEPENDENT DERIVED PARAMETERS:\n");
    fprintf(stream, "  Distance                        = %Lf (%Lf) kpc\n", distance, distanceerror);
    if(paramset[TEMPO3_PARAM_PMRA] && paramset[TEMPO3_PARAM_PMDEC]) {
      long double vt, vterr;
      //      vt = (1e-3*(M_PI/180.0)*(1.0/(60.0*60.0))*(1.0/(365.25*24.0*60.0*60.0)) )*(1e3*3.08567758e16)*1e-3*pm*distance;
      vt = (M_PI*PARSEC*1e-3/(365.25*TEMPO3_SecPerDay*3600.0*180.0) )*pm*distance;
      fprintf(stream, "  Transverse velocity             = %Lf", vt);
      if(showerrors && pmerr > 0) {
	//	vterr = sqrt(vt*vt*(pmerr*pmerr/(pm*pm) + disterr*disterr/(dist*dist)));
	vterr = sqrtl(vt*vt*(pmerr*pmerr/(pm*pm)  + distanceerror*distanceerror/(distance*distance) ));
	fprintf(stream, " (%Lf)", vterr);
      }     
      fprintf(stream, " km/s\n");

      //Shklovskii
      long double skl, sklerr;
      skl = (vt*1e3)*(vt*1e3)*p/(TEMPO3_SPEEDOFLIGHT*distance*1e3*PARSEC);
      fprintf(stream, "  Shklovskii contribution Pdot    = %Le", skl);
      if(showerrors && pmerr > 0) {
	//	sklerr = skl*sqrt(4.0*vterr*vterr/(vt*vt)+perr*perr/(p*p)+distanceerror*distanceerror/(distance*distance));
	// Problem with the above: vt already depends on distance, so distance and vt are not independent variables
	// Solution: get the error on vt without the distance contribution
	// The skl formula using (vtnew=vt/d): skl = vtnew**2*P*d/c
	long double vterr_new_frac;
	vterr_new_frac = pmerr/pm;
	sklerr = skl*sqrtl(4.0*vterr_new_frac*vterr_new_frac+perr*perr/(p*p)+distanceerror*distanceerror/(distance*distance));
	fprintf(stream, " (%Le)", sklerr);
      }
      fprintf(stream, "  = %Lf%%\n", 100.0*skl/pdot);
      fprintf(stream, "  Intrinsic Pdot                  = %Le s/s\n", pdot - skl);
      for(i = 0; i < TEMPO3_MaxNrCompanions; i++) {
	if(parfile->companion_defined[i]) {
	  if(paramset[TEMPO3_PARAM_PB+i]) {
	    skl = (vt*1e3)*(vt*1e3)*parfile->parameter[TEMPO3_PARAM_PB+i]*TEMPO3_SecPerDay/(TEMPO3_SPEEDOFLIGHT*distance*1e3*PARSEC);
	    fprintf(stream, "  Shklovskii contribution PBDOT_%d = %Le", i+1, skl);
	    if(showerrors && pmerr > 0) {
	      // Similar to above
	      long double vterr_new_frac;
	      vterr_new_frac = pmerr/pm;
	      sklerr = skl*sqrtl(4.0*vterr_new_frac*vterr_new_frac+dplus[TEMPO3_PARAM_PB+i]*dplus[TEMPO3_PARAM_PB+i]/(parfile->parameter[TEMPO3_PARAM_PB+i]*parfile->parameter[TEMPO3_PARAM_PB+i])+distanceerror*distanceerror/(distance*distance));
	      fprintf(stream, " (%Le)", sklerr);
	    }
	    fprintf(stream, "  = %Lf%%\n", 100.0*skl/parfile->parameter[TEMPO3_PARAM_PBDOT+i]);
	    fprintf(stream, "  Intrinsic PBDOT_%d               = %Le s/s\n", i+1, parfile->parameter[TEMPO3_PARAM_PBDOT+i] - skl);
	  }
	}
      }

    }
  }

  fprintf(stream, "\n");
}





long double pickvalue(int param, tempo3_parameters_def *parfile)
{
  return parfile->parameter[param];
}


void setvaluevalue(int param, long double value, tempo3_parameters_def *parfile)
{
  parfile->parameter[param] = value; 
}




void fill_xstart(tempo3_parameters_def *parfile, long double *xstart)
{
  int i;
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    xstart[i] = parfile->parameter[i];
  }
}

void use_xfit(tempo3_parameters_def *parfile, long double *xfit)
{
  int i;
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    parfile->parameter[i] = xfit[i]; 
  }
}


/* Find decent initial step size */
long double get_initial_step_size(int param, tempo3_parameters_def *parfile, toa_data_def toas, residual_def residuals, tempo3_wraps_def *wraps, int superverbose, int *paramlinked, int *paramset, int nosort, int numthreads)
{
  long double dx, xbest, chi2_old, chi2, chi2w;
  xbest = pickvalue(param, parfile);
  setvaluevalue(param, xbest, parfile);
  calcphases(parfile, toas, 0, &residuals, paramlinked, paramset, 0, numthreads);
  subtractTurns(toas, &residuals, wraps, nosort, 1);
  calcChi2(&residuals, toas, &chi2_old, parfile->parameter[TEMPO3_PARAM_F0], &chi2w);
  if(parfile->mode != 0)
    chi2_old = chi2w;

  /* First get a very rough idea what stepsize would be good */
  dx = 1e-80;
  while(dx < 1) {
    dx *= 1000.0;
    setvaluevalue(param, xbest + dx, parfile);
    calcphases(parfile, toas, 0, &residuals, paramlinked, paramset, 0, numthreads);
    subtractTurns(toas, &residuals, wraps, nosort, 1);
    calcChi2(&residuals, toas, &chi2, parfile->parameter[TEMPO3_PARAM_F0], &chi2w);
    if(parfile->mode != 0)
      chi2 = chi2w;
    if(fabsl((chi2-chi2_old)/chi2_old) > 0.01)
      break;
  }

  /* Now refine */
  dx /= 1000.0;
  while(dx < 1) {
    dx *= 45.0;
    setvaluevalue(param, xbest + dx, parfile);
    calcphases(parfile, toas, 0, &residuals, paramlinked, paramset, 0, numthreads);
    subtractTurns(toas, &residuals, wraps, nosort, 1);
    calcChi2(&residuals, toas, &chi2, parfile->parameter[TEMPO3_PARAM_F0], &chi2w);
    if(parfile->mode != 0)
      chi2 = chi2w;
    if(fabsl((chi2-chi2_old)/chi2_old) > 0.01)
      break;
  }

  /* Now refine more */
  dx /= 45.0;
  while(dx < 1) {
    dx *= 2.0;
    setvaluevalue(param, xbest + dx, parfile);
    calcphases(parfile, toas, 0, &residuals, paramlinked, paramset, 0, numthreads);
    subtractTurns(toas, &residuals, wraps, nosort, 1);
    calcChi2(&residuals, toas, &chi2, parfile->parameter[TEMPO3_PARAM_F0], &chi2w);
    if(parfile->mode != 0)
      chi2 = chi2w;
    if(fabsl((chi2-chi2_old)/chi2_old) > 0.01)
      break;
  }

  setvaluevalue(param, xbest, parfile);
  if(superverbose != -1)
    printf("Initial step size param %d = %Le\n", param, dx);
  return dx;
}

void fill_dx(tempo3_parameters_def *parfile, toa_data_def toas, residual_def residuals, tempo3_wraps_def *wraps, int *fixed, long double *dx, int superverbose, int *paramlinked, int *paramset, int nosort, int numthreads)
{
  int i;
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    if(fixed[i] == 0) {
      dx[i] = get_initial_step_size(i, parfile, toas, residuals, wraps, superverbose, paramlinked, paramset, nosort, numthreads);
      //      printf("XXXXXXX dx[%d] = %Le\n", i, dx[i]);
    }
  }
}




/* if nummericalApproxNuNudot is set, an nummerical derivative of the
   spin-phase calculation is used rather than the analytic
   calculation */
void calcAlternativePlots(float *ygraph, float *ygraph2, char *ylabel, int plottype, int showdata, toa_data_def *toas, residual_def residuals, tempo3_parameters_def *parfile, int *paramlinked, int *paramset, int nummericalApproxNuNudot)
{
  long double phase2; /*phase, , dt;*/
  int i;

  ephemeris_setLinkedValues_par_version(parfile, paramlinked);

  /* First plottype is normal residual plot, so nothing needs to be done */
  if(plottype == 0 || plottype == PLOTTYPE_FREQUENCY || plottype == PLOTTYPE_ORBITALPHASE) {
    for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
      ygraph[i] = residuals.resid[i];
    }
    sprintf(ylabel, "\\gD\\gf [phase]");
  }else if(plottype >= 1 && plottype <= 4) {
    /* Calculate nu versus mjd */
    for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
      /*
      phase = evaluate_ephemeris(toas->mjd[i], toas->freq[i], 0, parfile, 0, paramset, 0, 0);
      dt = 0.5/parfile->parameter[TEMPO3_PARAM_F0];
      phase2 = evaluate_ephemeris(toas->mjd[i]+dt/TEMPO3_SecPerDay, toas->freq[i], 0, parfile, 0, paramset, 0);
      phase2 -= phase; */
      if(nummericalApproxNuNudot)
	phase2 = ephemeris_estimateNu_numerically_avoid_using(toas->mjd[i], toas->freqSSB[i], toas->ssb1[i], toas->ssb2[i], toas->ssb3[i], toas->site[i], toas->flags[i], parfile, paramset);
      else 
	phase2 = evaluate_ephemeris(toas->mjd[i], toas->freqSSB[i], toas->ssb1[i], toas->ssb2[i], toas->ssb3[i], 0, 0, toas->site[i], toas->flags[i], parfile, 0, paramset, 1, 0);
      residuals.resid[i] = phase2; 
      ygraph[i] = residuals.resid[i];
      ygraph2[i] = toas->nu[i];
    }
    sprintf(ylabel, "\\gn [Hz]");

    if(plottype >= 2) {
      /* Calculate nu-F0 versus mjd */
      for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
	/*	phase2 = parfile->f1*(toas->mjd[i]-parfile->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay; */
	residuals.resid[i] -= parfile->parameter[TEMPO3_PARAM_F0];
	toas->nu[i] -= parfile->parameter[TEMPO3_PARAM_F0];
	ygraph[i] = residuals.resid[i];
	ygraph2[i] = toas->nu[i];
      }
      sprintf(ylabel, "\\gn-F\\d0\\u [Hz]");

      /* Calculate nu minus F0 and F1 contribution versus mjd */
      if(plottype >= 3) {
	for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
	  phase2 = parfile->parameter[TEMPO3_PARAM_F0+1]*(toas->mjd[i]-parfile->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
	  residuals.resid[i] -= phase2;
	  ygraph[i] = residuals.resid[i];

	  phase2 = parfile->parameter[TEMPO3_PARAM_F0+1]*(toas->nu_mjd[i]-parfile->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
	  toas->nu[i] -= phase2;
	  ygraph2[i] = toas->nu[i];
	}
	sprintf(ylabel, "\\gn-F\\d0\\u-F\\d1\\u\\gDt [Hz]");

      /* Calculate nu minus F0, F1 and F2 contribution versus mjd */
	if(plottype >= 4) {
	  for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
	    phase2 = 0.5*parfile->parameter[TEMPO3_PARAM_F0+2]*(toas->mjd[i]-parfile->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay*(toas->mjd[i]-parfile->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
	    residuals.resid[i] -= phase2;
	    ygraph[i] = residuals.resid[i];

	    phase2 = 0.5*parfile->parameter[TEMPO3_PARAM_F0+2]*(toas->nu_mjd[i]-parfile->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay*(toas->nu_mjd[i]-parfile->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
	    toas->nu[i] -= phase2;
	    ygraph2[i] = toas->nu[i];
	  }
	  sprintf(ylabel, "\\gn-F\\d0\\u-F\\d1\\u\\gDt-F\\d2\\u\\gDt\\u2\\d [Hz]");
	}
      }
    }
  }else if(plottype >= 5 && plottype <= 7) {
    /* Calculate nudot versus mjd */
    for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
      if(nummericalApproxNuNudot)
	residuals.resid[i] = ephemeris_estimateNudot_numerically_avoid_using(toas->mjd[i], toas->freqSSB[i], toas->ssb1[i], toas->ssb2[i], toas->ssb3[i], toas->site[i], toas->flags[i], *parfile, paramset);
      else 
	residuals.resid[i] = evaluate_ephemeris(toas->mjd[i], toas->freqSSB[i], toas->ssb1[i], toas->ssb2[i], toas->ssb3[i], 0, 0, toas->site[i], toas->flags[i], parfile, 0, paramset, 2, 0);

      /*
      phase = evaluate_ephemeris(toas->mjd[i]+9.5, toas->freqSSB[i], 0, parfile, 0, paramset, 0, 0);
      dt = 0.5/parfile->parameter[TEMPO3_PARAM_F0];
      phase2 = evaluate_ephemeris(toas->mjd[i]+9.5+dt/TEMPO3_SecPerDay, toas->freqSSB[i], 0, parfile, 0, paramset, 0, 0);
      phase2 -= phase;
      phase2 = phase2/dt;
      residuals.resid[i] = phase2;
      phase = evaluate_ephemeris(toas->mjd[i]-9.5, toas->freqSSB[i], 0, parfile, 0, paramset, 0, 0);
      dt = 0.5/parfile->parameter[TEMPO3_PARAM_F0];
      phase2 = evaluate_ephemeris(toas->mjd[i]-9.5+dt/TEMPO3_SecPerDay, toas->freqSSB[i], 0, parfile, 0, paramset, 0, 0);
      phase2 -= phase;
      phase2 = phase2/dt;
      residuals.resid[i] -= phase2;
      residuals.resid[i] /= 19.0*TEMPO3_SecPerDay;
      */
      
      ygraph[i] = residuals.resid[i];
      ygraph2[i] = toas->nudot[i];
    }
    sprintf(ylabel, "\\gn\\b\\u\\u \\(2198)\\(2210)\\d\\d [Hz/s]");

    if(plottype >= 6) {
      /* Calculate nudot minus F1 contribution versus mjd */
      for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
	residuals.resid[i] -= parfile->parameter[TEMPO3_PARAM_F0+1];
	toas->nudot[i]      -= parfile->parameter[TEMPO3_PARAM_F0+1];
	ygraph[i] = residuals.resid[i];
	ygraph2[i] = toas->nudot[i];
      }
      sprintf(ylabel, "\\gn\\b\\u\\u \\(2198)\\(2210)\\d\\d - F\\d1\\u [Hz/s]");

      if(plottype >= 7) {
	/* Calculate nudot minus F1 and F2 contribution versus mjd */
	for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
	  phase2 = parfile->parameter[TEMPO3_PARAM_F0+2]*(toas->mjd[i]-parfile->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
	  residuals.resid[i] -= phase2;
	  ygraph[i] = residuals.resid[i];

	  phase2 = parfile->parameter[TEMPO3_PARAM_F0+2]*(toas->nu_mjd[i]-parfile->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
	  toas->nudot[i] -= phase2;
	  ygraph2[i] = toas->nudot[i];
	}
	sprintf(ylabel, "\\gn\\b\\u\\u \\(2198)\\(2210)\\d\\d - F\\d1\\u - F\\d2\\u\\gDt [Hz/s]");
      }
    }
  }else if(plottype >= 5 && plottype <= 7) {
    /* Calculate nudot versus mjd */
    for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
      if(nummericalApproxNuNudot)
	residuals.resid[i] = ephemeris_estimateNudot_numerically_avoid_using(toas->mjd[i], toas->freqSSB[i], toas->ssb1[i], toas->ssb2[i], toas->ssb3[i], toas->site[i], toas->flags[i], *parfile, paramset);
      else 
	residuals.resid[i] = evaluate_ephemeris(toas->mjd[i], toas->freqSSB[i], toas->ssb1[i], toas->ssb2[i], toas->ssb3[i], 0, 0, toas->site[i], toas->flags[i], parfile, 0, paramset, 2, 0);

      /*
      phase = evaluate_ephemeris(toas->mjd[i]+9.5, toas->freqSSB[i], 0, parfile, 0, paramset, 0, 0);
      dt = 0.5/parfile->parameter[TEMPO3_PARAM_F0];
      phase2 = evaluate_ephemeris(toas->mjd[i]+9.5+dt/TEMPO3_SecPerDay, toas->freqSSB[i], 0, parfile, 0, paramset, 0, 0);
      phase2 -= phase;
      phase2 = phase2/dt;
      residuals.resid[i] = phase2;
      phase = evaluate_ephemeris(toas->mjd[i]-9.5, toas->freqSSB[i], 0, parfile, 0, paramset, 0, 0);
      dt = 0.5/parfile->parameter[TEMPO3_PARAM_F0];
      phase2 = evaluate_ephemeris(toas->mjd[i]-9.5+dt/TEMPO3_SecPerDay, toas->freqSSB[i], 0, parfile, 0, paramset, 0, 0);
      phase2 -= phase;
      phase2 = phase2/dt;
      residuals.resid[i] -= phase2;
      residuals.resid[i] /= 19.0*TEMPO3_SecPerDay;
      */
      
      ygraph[i] = residuals.resid[i];
      ygraph2[i] = toas->nudot[i];
    }
    sprintf(ylabel, "\\gn\\b\\u\\u \\(2198)\\(2210)\\d\\d [Hz/s]");

    if(plottype >= 6) {
      /* Calculate nudot minus F1 contribution versus mjd */
      for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
	residuals.resid[i] -= parfile->parameter[TEMPO3_PARAM_F0+1];
	toas->nudot[i]      -= parfile->parameter[TEMPO3_PARAM_F0+1];
	ygraph[i] = residuals.resid[i];
	ygraph2[i] = toas->nudot[i];
      }
      sprintf(ylabel, "\\gn\\b\\u\\u \\(2198)\\(2210)\\d\\d - F\\d1\\u [Hz/s]");

      if(plottype >= 7) {
	/* Calculate nudot minus F1 and F2 contribution versus mjd */
	for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
	  phase2 = parfile->parameter[TEMPO3_PARAM_F0+2]*(toas->mjd[i]-parfile->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
	  residuals.resid[i] -= phase2;
	  ygraph[i] = residuals.resid[i];

	  phase2 = parfile->parameter[TEMPO3_PARAM_F0+2]*(toas->nu_mjd[i]-parfile->parameter[TEMPO3_PARAM_PEPOCH])*TEMPO3_SecPerDay;
	  toas->nudot[i] -= phase2;
	  ygraph2[i] = toas->nudot[i];
	}
	sprintf(ylabel, "\\gn\\b\\u\\u \\(2198)\\(2210)\\d\\d - F\\d1\\u - F\\d2\\u\\gDt [Hz/s]");
      }
    }
  }else if(plottype == 8) {
    /* Calculate nuddot versus mjd */
    for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
      if(nummericalApproxNuNudot) {
	printerror(0, "Second numerical derivative is not (yet?) implemented\n");
	exit(0);
      }else 
	residuals.resid[i] = evaluate_ephemeris(toas->mjd[i], toas->freqSSB[i], toas->ssb1[i], toas->ssb2[i], toas->ssb3[i], 0, 0, toas->site[i], toas->flags[i], parfile, 0, paramset, 3, 0);

      ygraph[i] = residuals.resid[i];
      ygraph2[i] = residuals.resid[i];   // No measured value yet, should be filled in!!!!
      //      printf("XXXX: %Le\n", residuals.resid[i]);
      if(i == 0)
	printwarning(0, "MEASURED VALUE IS NOT FILLED IN");
    }
    sprintf(ylabel, "\\gn\\b\\u\\u\\(2198)\\(2210)\\(2210)\\d\\d [Hz/s/s]");

  }else if(plottype == PLOTTYPE_NBREAK) {
    /* Calculate braking index versus mjd */
    for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
      if(nummericalApproxNuNudot) {
	printerror(0, "Second numerical derivative is not (yet?) implemented\n");
	exit(0);
      }else {
	long double nu, nudot, nuddot;
	nu     = evaluate_ephemeris(toas->mjd[i], toas->freqSSB[i], toas->ssb1[i], toas->ssb2[i], toas->ssb3[i], 0, 0, toas->site[i], toas->flags[i], parfile, 0, paramset, 1, 0);
	nudot  = evaluate_ephemeris(toas->mjd[i], toas->freqSSB[i], toas->ssb1[i], toas->ssb2[i], toas->ssb3[i], 0, 0, toas->site[i], toas->flags[i], parfile, 0, paramset, 2, 0);
	nuddot = evaluate_ephemeris(toas->mjd[i], toas->freqSSB[i], toas->ssb1[i], toas->ssb2[i], toas->ssb3[i], 0, 0, toas->site[i], toas->flags[i], parfile, 0, paramset, 3, 0);
	residuals.resid[i] = nu*nuddot/(nudot*nudot);
      }
      ygraph[i] = residuals.resid[i];
      ygraph2[i] = residuals.resid[i];   // No measured value yet, should be filled in!!!!
      //      printf("XXXX: %Le\n", residuals.resid[i]);
      if(i == 0)
	printwarning(0, "MEASURED VALUE IS NOT FILLED IN");
    }
    sprintf(ylabel, "n");

  }
}


/*
     - - - - -
      C L D J
     - - - - -

  Gregorian Calendar to Modified Julian Date

  Given:
     IY,IM,ID     int    year, month, day in Gregorian calendar

  Returned:
     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs
     J            int    status:
                           0 = OK
                           1 = bad year   (MJD not computed)
                           2 = bad month  (MJD not computed)
                           3 = bad day    (MJD computed)

  The year must be -4699 (i.e. 4700BC) or later.

  The algorithm is derived from that of Hatcher 1984
  (QJRAS 25, 53-55).

  P.T.Wallace   Starlink   December 1985
-
*/
void sla_CLDJ (int IY, int IM, int ID, double *DJM, int *J) 
{

  /*  Month lengths in days */
  int MTAB[12] = {31,28,31,30,31,30,31,31,30,31,30,31};


  /*  Preset status */
  *J=0;

  /*  Validate year */
  if (IY<-4699)
    *J=1;
  else {
    /*     Validate month */
    if (IM>=1 && IM<=12) {
      /*        Allow for leap year */
      if ((IY % 4) == 0) 
	MTAB[1]=29;
      else
	MTAB[1]=28;
      if ((IY % 100) == 0 && (IY % 400) != 0)
	MTAB[1]=28;

      /*        Validate day */
      if (ID < 1 || ID > MTAB[IM-1]) 
	*J=3;

      /*        Modified Julian Date */
      *DJM=((1461*(IY-(12-IM)/10+4712))/4
	    +(306*((IM+9)%12)+5)/10
	    -(3*((IY-(12-IM)/10+4900)/100))/4
	    +ID-2399904);

      /*        Bad month */
    }else {
      *J=2;
    }
  }  /* End else year */
}

/*
     - - - - - -
      C A L D J
     - - - - - -

  Gregorian Calendar to Modified Julian Date

  (Includes century default feature:  use sla_CLDJ for years
   before 100AD.)

  Given:
     IY,IM,ID     int    year, month, day in Gregorian calendar

  Returned:
     DJM          dp     modified Julian Date (JD-2400000.5) for 0 hrs
     J            int    status:
                           0 = OK
                           1 = bad year   (MJD not computed)
                           2 = bad month  (MJD not computed)
                           3 = bad day    (MJD computed)

  Acceptable years are 00-49, interpreted as 2000-2049,
                       50-99,     "       "  1950-1999,
                       100 upwards, interpreted literally.

  Called:  sla_CLDJ

  P.T.Wallace   Starlink   November 1985
*/

void sla_CALDJ(int IY, int IM, int ID, double *DJM, int *J)
{
  int NY;
  /*  Default century if appropriate */
  if(IY >= 0 && IY < 49) {
    NY=IY+2000;
  }else if(IY>=50 && IY <99) {
    NY=IY+1900;
  }else {
    NY=IY;
  }
  /*  Modified Julian Date */
  sla_CLDJ(NY,IM,ID,DJM,J);
}

/* device is doing nothing */
int doplot(char *device, float *ygraph, float *ygraph2, int nophaselines, int plottype, int showdata, int showWaves, toa_data_def toas, residual_def *residuals, ephemeris_def *eph, int superverbose, int writeoutpieces, int plottitle, float *cur_x1, float *cur_y1, float *cur_x2, float *cur_y2, float previous_xstart, float previous_xend, toa_data_def toas_fiducial, residual_def residuals_fiducial, int showfiducial, int showpepoch, char *highlightedFlag, int nummericalApproxNuNudot, char *title, int amoeba_algorithm, float markersize, float stridefit_nrdays, float stridefit_mintoarange, float stridefit_minstride, int stridefit_minnroftoas, int stridefit_includeF2, int stridefit_includeF3, char *stridefit_scriptfilename, int nosort, int numthreads, verbose_definition verbose)
{
  char txt[1000], txt2[1000];
  float min, max, x, y, yold;
  /*, freqmin, freqmax;*/
  int i, j, notset, nrtoasincluded, nrparamsused;
  FILE *fout;
  long nturn;
  int plotmonths, plotyears, day, yr, mt;
  double djm;
  float ticklength;
  char months[12][10];

  if(((plottype >= 1 && plottype <= 8) || plottype == PLOTTYPE_NBREAK) && showdata) {
    if(measureNuNudot(toas.nu, toas.nudot, toas.nu_mjd, stridefit_nrdays, stridefit_mintoarange, stridefit_minstride, stridefit_minnroftoas, ygraph, ygraph2, &toas, *residuals, eph, superverbose, writeoutpieces, amoeba_algorithm, stridefit_includeF2, stridefit_includeF3, stridefit_scriptfilename, nosort, numthreads, verbose) == 0)
      return 0;
    /*    for(i = 0; i < residuals.nrtoas; i++) {
      ygraph2[i] = toas.nu[i];
    }
    pgplot_options.viewport.dontopen = 1;
    pgplot_options.viewport.dontopen = 1;
    strcpy(pgplot_options.viewport.plotDevice, device);
    if(pgplotGraph1(&pgplot_options, ygraph2, toas.mjdx, NULL, toas.nrtoas, -1, -1, residuals.mjdmin, residuals.mjdmax, 0, "MJD", txt2, txt, 1, 1, 17, 3, 1, NULL, -1, noverbose) == 0) {
      return 0;
    }
    pgetch();
    */
    /*
    if(writeoutpieces) {
      sprintf(txt, "plot%d.txt", plottype);
      fout = fopen(txt, "w");
      if(fout == NULL) {
	printerror(0, "doplot: Cannot open file %s\n", txt);
	return 0;
      }
      printf("Writen file '%s'\n", txt);
      for(i = 0; i < toas.nrtoas; i++) {
	printf("%Le\t%Le\t%Le\n", toas.nu_mjd[i], toas.nu[i], toas.nudot[i]);
      }
      fclose(fout);
    }
    */
  }

  calcphases(&(eph->par), toas, showWaves, residuals, eph->paramlinked, eph->paramset, 0, numthreads);
  subtractTurns(toas, residuals, eph->wraps, nosort, eph->par.track);
  calcChi2(residuals, toas, &(residuals->variance), eph->par.parameter[TEMPO3_PARAM_F0], &(residuals->chi2));
  nrtoasincluded = 0;
  for(j = 0; j < toas.nrtoas+toas.nrtoasdeleted; j++) {
    if(toas.deleted[j] == 0) 
      if(toas.mjd[j] >= residuals->mjdmin && toas.mjd[j] <= residuals->mjdmax) 
	nrtoasincluded += 1;
  }
  eph->par.parameter[TEMPO3_PARAM_TRES] = 1000000*sqrtl(residuals->variance/(long double)nrtoasincluded)/eph->par.parameter[TEMPO3_PARAM_F0];

  if(plottitle) {
    if(eph->par.mode == 0)
      sprintf(txt, "%s: RMS = %.2Lf milli-phases (%Lf \\gms)", eph->par.psrjname, 1000*sqrtl(residuals->variance/(long double)toas.nrtoas), eph->par.parameter[TEMPO3_PARAM_TRES]); 
    else {
      nrparamsused = 0;
      for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
	if(eph->fixed[i] == 0)
	  nrparamsused++;
      }
      sprintf(txt, "%s: \\gx\\u2\\d = %.2Lf  (reduced \\gx\\u2\\d = %.2Lf with %d dof.)", eph->par.psrjname, residuals->chi2, residuals->chi2/(long double)(toas.nrtoas-nrparamsused), toas.nrtoas-nrparamsused);
    }
  }else {
    txt[0] = 0;
  }
  calcAlternativePlots(ygraph, ygraph2, txt2, plottype, showdata, &toas, *residuals, &(eph->par), eph->paramlinked, eph->paramset, nummericalApproxNuNudot);
  if(writeoutpieces) {
    sprintf(txt, "model%d.txt", plottype);
    fout = fopen(txt, "w");
    if(fout == NULL) {
      printerror(0, "doplot: Cannot open file %s\n", txt);
      return 0;
    }
    printf("Writen file '%s'\n", txt);
    for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
      if(toas.deleted[i] == 0) {
	if(toas.nu_mjd[i] > 0) {
	  fprintf(fout, "%Le\t%e\n", toas.nu_mjd[i], ygraph[i]);
	}
      }
    }
    fclose(fout);
  }
  if((plottype >= 1 && plottype <= 8) || plottype == PLOTTYPE_NBREAK)
    sprintf(txt, "Model calculation");



  /*
  pgplot_options.viewport.dontopen = 1;
  pgplot_options.viewport.dontopen = 1;
  strcpy(pgplot_options.viewport.plotDevice, device);
  if(pgplotGraph1(&pgplot_options, ygraph, toas.mjdx, NULL, toas.nrtoas, -1, -1, residuals->mjdmin, residuals->mjdmax, 0, "MJD", txt2, txt, 1, 1, 17, 3, 1, NULL, -1, noverbose) == 0) {
    return 0;
  }
 */
  ppgask(0);
  ppgsch(1);
  ppgslw(1);
  ppgpage();
  ppgsvp(0.1, 0.9, 0.1, 0.9);

  notset = 1;
  min = max = 0;
  for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
    if(toas.deleted[i] == 0) {
      if(toas.mjd[i] >= residuals->mjdmin && toas.mjd[i] <= residuals->mjdmax) {
	if(notset) {
	  min = max = ygraph[i];
	  /*	  freqmin = freqmax = toas.freqx[i]; */
	  notset = 0;
	}
	if(ygraph[i] > max)
	  max = ygraph[i];
	if(ygraph[i] < min)
	  min = ygraph[i];
	/*	if(toas.freqx[i] > freqmax)
	  freqmax = toas.freqx[i];
		if(toas.freqx[i] < freqmin)
	  freqmin = toas.freqx[i]; */
	if(showdata) {
	  /* If mjd is negative it means that the point should be excluded */
	  if(toas.nu_mjd[i] > 0) {
	    if(ygraph2[i] > max)
	      max = ygraph2[i];
	    if(ygraph2[i] < min)
	      min = ygraph2[i];
	  }
	}
      }
    }
  }

  /* Calculate the residual corresponding to the fiducial point if requested */
  if(plottype == 0 && toas_fiducial.nrtoas == 1 && showfiducial != 0) {
    calcphases(&(eph->par), toas_fiducial, showWaves, &residuals_fiducial, eph->paramlinked, eph->paramset, 0, numthreads);
    subtractTurns(toas_fiducial, &residuals_fiducial, eph->wraps, nosort, 0);
    y = residuals_fiducial.resid[0];
    if(y > max)
      max = y;
    if(y < min)
      min = y;
  }

  ppgsci(1);
  ppgmtxt("T", 3.5, 0.5, 0.5, title);
  if(plottype != PLOTTYPE_FREQUENCY && plottype != PLOTTYPE_ORBITALPHASE) {
    ppgswin(residuals->mjdmin, residuals->mjdmax, min, max);
    *cur_x1 = residuals->mjdmin;
    *cur_x2 = residuals->mjdmax;
    *cur_y1 = min;
    *cur_y2 = max;
    ppglab("MJD", "", txt);
    /* Removed 'c' from first params (top) and first (right) */
    if((max-min != 0)) {
    /* Normally plot bottom mjd axis myself, unless too much zoomed in, causing pgaxis to fail */
      ppgbox("b1",0.0,0,"bntsi",0.0,0);
      // Depending on how much zoomed in, try different pgaxis commands. 
      if(floor(residuals->mjdmax) - floor(residuals->mjdmin) > 2) { // If more than 2 days span, pgplot does probably something sensible
	ppgaxis("n1", (float)residuals->mjdmin, min, (float)residuals->mjdmax, min, (float)residuals->mjdmin, (float)residuals->mjdmax, 0.0, 0, 0, 0.7, 0.3, 1.0, 0.0);
      }else if(floor(residuals->mjdmax) - floor(residuals->mjdmin ) > 1) { // If more than 1 days span, force 1 day tickmarks, as otherwise pgaxis tends to come up with corrupt numbers (too many steps away from zero?)
	ppgaxis("n1", (float)residuals->mjdmin, min, (float)residuals->mjdmax, min, (float)residuals->mjdmin, (float)residuals->mjdmax, 1.0, -1, 0, 0.7, 0.3, 1.0, 0.0);
      }else { // If too much zoomed in, do not show axis labels, but with exponential notation, which appears to work better.
	ppgaxis("n2", (float)residuals->mjdmin, min, (float)residuals->mjdmax, min, (float)residuals->mjdmin, (float)residuals->mjdmax, 1.0, -1, 0, 0.7, 0.3, 1.0, 0.0);
      }
      //      printf("XXXXX %f %f\n", (float)residuals->mjdmin, (float)residuals->mjdmax);
    }else {
      if(floor(residuals->mjdmax) - floor(residuals->mjdmin ) > 1) {
	ppgbox("bnsti",0.0,0,"bntsi",0.0,0);
      }else {
	ppgbox("bsti",0.0,0,"bntsi",0.0,0);
      }
    }
    /* Right time axis */
    ppgaxis("n", residuals->mjdmax, min, residuals->mjdmax, max, min*1000.0/eph->par.parameter[TEMPO3_PARAM_F0], max*1000.0/eph->par.parameter[TEMPO3_PARAM_F0], 0.0, 0, 0.0, 0.7, 0.3, 1.0, 0.0);
    /* Now place the y-axis tick marks for the dates */
    ppgaxis("", residuals->mjdmin, max, residuals->mjdmax, max,residuals->mjdmin, residuals->mjdmax, 100000.0, 0, 0, 0.7, 0.3, 1.0, 0.0); 
    sprintf(months[0], "Jan");
    sprintf(months[1], "Feb");
    sprintf(months[2], "Mar");
    sprintf(months[3], "Apr");
    sprintf(months[4], "May");
    sprintf(months[5], "Jun");
    sprintf(months[6], "Jul");
    sprintf(months[7], "Aug");
    sprintf(months[8], "Sep");
    sprintf(months[9], "Oct");
    sprintf(months[10], "Nov");
    sprintf(months[11], "Dec");
    if(residuals->mjdmax-residuals->mjdmin < 365*0.5)
      plotmonths = 1;
    else if(residuals->mjdmax-residuals->mjdmin < 365)
      plotmonths = 4;
    else if(residuals->mjdmax-residuals->mjdmin < 2*365)
      plotmonths = 6;
    else
      plotmonths = 12;
    plotyears = 1;
    if(residuals->mjdmax-residuals->mjdmin > 5*365)
      plotyears = 2;
    if(residuals->mjdmax-residuals->mjdmin > 10*365)
      plotyears = 4;
    day = 1;
    for(yr = 1967; yr < 2050; yr++) {
      for(mt = 1; mt <= 12; mt++) {
	sla_CALDJ(yr, mt, day, &djm, &j);
	if(j != 0) {
	  printf("SLA ERROR %d (%d-%d-%d)\n", j, yr,mt,day);
	  return 0;
	}
	if(djm >= residuals->mjdmin && djm <= residuals->mjdmax) {
	  if((mt-1)%plotmonths == 0 && yr%plotyears == 0) {
	    sprintf(txt, "%s %d", months[mt-1], yr);
	    ticklength = 0.7;
	  }else if(mt == 1) {
	    txt[0] = 0;
	    ticklength = 0.7;
	  }else {
	    txt[0] = 0;
	    ticklength = 0.3;
	  }
	  ppgtick(djm-0.5, max, djm+0.5, max, 0.5, 0.0, ticklength, -0.2, 0, txt);
	}
      }
    }
 
    if(previous_xstart != 0 || previous_xend != 0) {
      ppgarro(previous_xstart, min + 1.00*(max-min), previous_xstart, min + 0.96*(max-min));
      ppgarro(previous_xend, min + 1.00*(max-min), previous_xend, min + 0.96*(max-min));
    }

    /* Plot the fitwaves function */
    ppgsci(6);
    int nwrapsplus, nwrapsmin, wrapnr;
    nwrapsplus = floor(max+0.5);
    nwrapsmin = -floor(-min+0.5);
    if(showWaves && eph->par.nrwaves > 0) {
      for(wrapnr = nwrapsmin; wrapnr <= nwrapsplus; wrapnr++) {
	yold = 1000;
	for(i = 0; i < toas.nrtoas+toas.nrtoasdeleted-1; i++) {
	  for(j = 0; j <= 100; j++) {
	    x = toas.mjdx[i] + (toas.mjdx[i+1]-toas.mjdx[i])*j*0.01;
	    y = -fitwave_function(&(eph->par), x, 0, -showWaves);  /* Calculate shape fitwaves function which was not taken into account in residual calculation */
	    /* Hope F0 is more or less good to convert time into phase */
	    y *= eph->par.parameter[TEMPO3_PARAM_F0];
	    nturn = y+0.5;
	    y -= nturn;
	    if(y < -0.5)
	      y += 1.0;
	    if(y > 0.5)
	      y -= 1.0;
	    /*	printf("x=%f y=%f\n", x, y); */
	    if(fabs(yold-y) > 0.1)
	      ppgmove(x, y+wrapnr);
	    else
	      ppgdraw(x, y+wrapnr);
	    yold = y;
	  }
	}
      }
    }
    ppgsci(1);
  }else if(plottype == PLOTTYPE_FREQUENCY) {
    ppgswin(residuals->freqmin, residuals->freqmax, min, max);
    *cur_x1 = residuals->freqmin;
    *cur_x2 = residuals->freqmax;
    *cur_y1 = min;
    *cur_y2 = max;
    ppglab("Frequency", "", txt);
    ppgbox("bnstic",0.0,0,"bntsi",0.0,0);
    /* Right time axis */
    ppgaxis("n", residuals->freqmax, min, residuals->freqmax, max, min*1000.0/eph->par.parameter[TEMPO3_PARAM_F0], max*1000.0/eph->par.parameter[TEMPO3_PARAM_F0], 0.0, 0, 0.0, 0.7, 0.3, 1.0, 0.0);
  }else if(plottype == PLOTTYPE_ORBITALPHASE) {  // Orbital phase plot
    ppgswin(0, 1, min, max);
    *cur_x1 = 0.0;
    *cur_x2 = 1.0;
    *cur_y1 = min;
    *cur_y2 = max;
    ppglab("(t-t0)/pb", "", txt);
    ppgbox("bnstic",0.0,0,"bntsi",0.0,0);
    /* Right time axis */
    ppgaxis("n", 1.0, min, 1.0, max, min*1000.0/eph->par.parameter[TEMPO3_PARAM_F0], max*1000.0/eph->par.parameter[TEMPO3_PARAM_F0], 0.0, 0, 0.0, 0.7, 0.3, 1.0, 0.0);
  }
  ppgmtxt("L", 3.5, 0.5, 0.5, txt2);
  ppgmtxt("R", 3.5, 0.5, 0.5, "\\gD\\gf [ms]");


  if(showpepoch && plottype != PLOTTYPE_ORBITALPHASE && plottype != PLOTTYPE_FREQUENCY) {
    ppgsci(2);
    ppgslw(1);
    ppgsls(2);
    ppgmove(eph->par.parameter[TEMPO3_PARAM_PEPOCH], *cur_y1);
    ppgdraw(eph->par.parameter[TEMPO3_PARAM_PEPOCH], *cur_y2);
    ppgsls(1);
    ppgsci(1);
  }


  for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
    if(toas.deleted[i] == 0) {
      if((plottype >= 0 && plottype <= 8) || plottype == PLOTTYPE_NBREAK) {
	if(showdata) {
	  y = ygraph2[i];
	  x = toas.nu_mjd[i];
	}else {
	  y = ygraph[i];
	  x = toas.mjdx[i];
	}
      }else if(plottype == PLOTTYPE_FREQUENCY) {
	y = ygraph[i];
	x = toas.freqx[i];
      }else if(plottype == PLOTTYPE_ORBITALPHASE) {
	y = ygraph[i];
	if(eph->paramset[TEMPO3_PARAM_PB] == 0 || eph->paramset[TEMPO3_PARAM_T0] == 0) {
	  printwarning(0, "No binary period or value of t0 defined.");
	  break;
	}
	x = (toas.mjdx[i]-eph->par.parameter[TEMPO3_PARAM_T0])/eph->par.parameter[TEMPO3_PARAM_PB];
	x = derotate_deg_double(x*360.0)/360.0;
      }
      ppgsci(1);

      draw_residual(toas, i, x, y, highlightedFlag, markersize, plottype, eph->par.mode, eph->par.parameter[TEMPO3_PARAM_F0]);
    }
  }
  if(showdata) {
    ppgsci(1);
    for(i = 0; i < toas.nrtoas + toas.nrtoasdeleted; i++) {
      if(toas.deleted[i] == 0) {
	ppgpt1(toas.mjdx[i], ygraph[i], -2);      
      }
    }
  }

  // To ensure glitch epochs are set correctly if VALUE parameters are used
  ephemeris_setLinkedValues_par_version(&(eph->par), eph->paramlinked);
  ppgsci(6);
  for(i = 0; i < eph->par.nrglitches; i++) {
    //    ppgmove(eph->par.parameter[TEMPO3_PARAM_GLEP+i], -100000);
    //    ppgdraw(eph->par.parameter[TEMPO3_PARAM_GLEP+i], 100000);
    ppgmove(eph->par.parameter[TEMPO3_PARAM_GLEP+i], *cur_y1);
    ppgdraw(eph->par.parameter[TEMPO3_PARAM_GLEP+i], *cur_y2);
  }
  if(nophaselines == 0 && plottype == 0) {
    ppgsci(7);
    for(i = 0; i < eph->wraps->nrwraps; i++) {
      ppgmove(eph->wraps->mjd[i], *cur_y1);
      ppgdraw(eph->wraps->mjd[i], *cur_y2);
    }
  }

  /* Show the fiducial point if defined */
  if(plottype == 0) {
    if(toas_fiducial.nrtoas == 1 && showfiducial != 0) {
      y = residuals_fiducial.resid[0];
      x = toas_fiducial.mjd[0];
      ppgsci(8);
      ppgsch(2);
      ppgslw(5);
      ppgpt1(x, y, 12);      
      ppgslw(1);
      ppgsch(1);
    }
  }

  ppgsci(1);
  return 1;
}


void dofit(int param, tempo3_parameters_def *parfile, toa_data_def toas, residual_def *residuals, tempo3_wraps_def *wraps, int *paramlinked, int *paramset, int nosort, int numthreads)
{
  int i;
  long double sign, xbest, xcur, chi2, chi2w, chi2_old;
  long double dx;


  /* Find decent initial step size */
  dx = 1e-80;
  xbest = pickvalue(param, parfile);
  sign = 1;
  while(dx < 1) {
    setvaluevalue(param, xbest, parfile);
    calcphases(parfile, toas, 0, residuals, paramlinked, paramset, 0, numthreads);
    subtractTurns(toas, residuals, wraps, nosort, 1);
    calcChi2(residuals, toas, &chi2_old, parfile->parameter[TEMPO3_PARAM_F0], &chi2w);
    if(parfile->mode != 0)
      chi2_old = chi2w;

    dx *= 2.0;
    setvaluevalue(param, xbest + sign*dx, parfile);
    calcphases(parfile, toas, 0, residuals, paramlinked, paramset, 0, numthreads);
    subtractTurns(toas, residuals, wraps, nosort, 1);
    calcChi2(residuals, toas, &chi2, parfile->parameter[TEMPO3_PARAM_F0], &chi2w);
    if(parfile->mode != 0)
      chi2 = chi2w;
    if(fabsl((chi2-chi2_old)/chi2_old) > 0.01)
      break;
  }
  setvaluevalue(param, xbest, parfile);
  if(additional_funk_info.superverbose != -1)
    printf("Initial step size param %d = %Le\n", param, dx);

  sign = 1;
  xbest = pickvalue(param, parfile);
  xcur = xbest;
  for(i = 0; i < 1000; i++) {
    calcphases(parfile, toas, 0, residuals, paramlinked, paramset, 0, numthreads);
    subtractTurns(toas, residuals, wraps, nosort, 1);
    calcChi2(residuals, toas, &chi2_old, parfile->parameter[TEMPO3_PARAM_F0], &chi2w);
    if(parfile->mode != 0)
      chi2_old = chi2w;
    xcur += sign*dx;
    setvaluevalue(param, xcur, parfile);
    calcphases(parfile, toas, 0, residuals, paramlinked, paramset, 0, numthreads);
    subtractTurns(toas, residuals, wraps, nosort, 1);
    calcChi2(residuals, toas, &chi2, parfile->parameter[TEMPO3_PARAM_F0], &chi2w);
    if(parfile->mode != 0)
      chi2 = chi2w;
    /*	printf("%Le %Le %Lf\n", xstart[param]-f0, sign*dx, chi2); */
    if(chi2 < chi2_old) {
      xbest = pickvalue(param, parfile);
    }else {
      dx = 0.5*dx;
      if(sign > 0)
	sign = -1.0;
      else
	sign = 1.0;
    }
    /*	  
pgplot_options.viewport.noclear = 1;
    pgplot_options.viewport.dontopen = 1;
    pgplot_options.viewport.dontopen = 1;
if(pgplotGraph1(&pgplot_options, resid, mjdx, NULL, nrtoas, -1, -1, -1, -1, 0, "MJD", "Resid", "", 1, 1, 1, 1, 1, NULL, -1, noverbose) == 0) {
      return 0;
      }*/
  }
  setvaluevalue(param, xbest, parfile);
}

// Return 0 on error
/* The _tmp arrays are not that useful anymore now they are no longer backing up global varables. */
int measureNuNudot(long double *nu, long double *nudot, long double *nu_mjd, float stridefit_nrdays, float stridefit_mintoarange, float stridefit_minstride, int minnrtoas, float *ygraph, float *ygraph2, toa_data_def *toas, residual_def residuals, ephemeris_def *eph, int superverbose, int writeoutpieces, int amoeba_algorithm, int stridefit_includeF2, int stridefit_includeF3, char *stridefit_scriptfilename, int nosort, int numthreads, verbose_definition verbose)
{
  int i, j, k, nfunk, ok, enabletracking;
  float prev_toa, prev_center;
  char txt[1000];
  toa_data_def toas_tmp;
  residual_def residuals_tmp;
  tempo3_parameters_def parfile_tmp;
  if(allocate_ephemeris_par_only(&parfile_tmp, 0, verbose) == 0) {
    printerror(verbose.debug, "ERROR measureNuNudot: Initialising ephemeris failed.");
    return 0;
  }
  tempo3_wraps_def wraps_tmp;
  long double *xstart, *dx, *xfit, rms, ftol = 1e-6;
  int *fixed2, *paramset2;
  FILE *fout, *fout_script;
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);
  verbose_definition noverbose;
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;

  enabletracking = eph->par.track;

  xstart = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
  dx = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
  xfit = malloc(TEMPO3_PARAM_TOTAL*sizeof(long double));
  fixed2 = malloc(TEMPO3_PARAM_TOTAL*sizeof(int));
  paramset2 = malloc(TEMPO3_PARAM_TOTAL*sizeof(int));

  if(xstart == NULL || dx == NULL || xfit == NULL || fixed2 == NULL || paramset2 == NULL) {
    printerror(verbose.debug, "ERROR tempo3: Cannot allocate memory\n");
    return 0;
  }


  /*
  if(superverbose == 0) {
    additional_funk_info.superverbose = -1;
  }else {
    additional_funk_info.superverbose = 0;
  }
  */

  /*  printf("XXXXX %d\n", superverbose); */

  /* Estimate F0 and F1 from current model */
  calcAlternativePlots(ygraph, ygraph2, txt, 1, 0, toas, residuals, &(eph->par), eph->paramlinked, eph->paramset, 1);
  for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
    nu[i] = residuals.resid[i];
  }
  calcAlternativePlots(ygraph, ygraph2, txt, 5, 0, toas, residuals, &(eph->par), eph->paramlinked, eph->paramset, 1);
  for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
    nudot[i] = residuals.resid[i];
  }

  /* Make a backup of original stuctures used by funk */
  memcpy(&residuals_tmp, &residuals, sizeof(residual_def));
  memcpy(&toas_tmp, toas, sizeof(toa_data_def));
  copy_parfile(&parfile_tmp, &(eph->par));
  //  memcpy(&parfile_tmp, parfile, sizeof(tempo3_parameters_def));
  memcpy(&wraps_tmp, eph->wraps, sizeof(tempo3_wraps_def));

  ephemeris_setLinkedValues_par_version(&parfile_tmp, eph->paramlinked);



  residuals.residy = malloc((toas_tmp.nrtoas+toas_tmp.nrtoasdeleted)*sizeof(float));
  residuals.resid = malloc((toas_tmp.nrtoas+toas_tmp.nrtoasdeleted)*sizeof(long double));
  residuals.nturn = malloc((toas_tmp.nrtoas+toas_tmp.nrtoasdeleted)*sizeof(long long));
  residuals.nturn_fixed = calloc((toas_tmp.nrtoas+toas_tmp.nrtoasdeleted), sizeof(int));
  toas->mjdx  = malloc((toas_tmp.nrtoas+toas_tmp.nrtoasdeleted)*sizeof(float));
  toas->mjd = malloc((toas_tmp.nrtoas+toas_tmp.nrtoasdeleted)*sizeof(long double));
  toas->freqSite = malloc((toas_tmp.nrtoas+toas_tmp.nrtoasdeleted)*sizeof(long double));
  toas->freqSSB = malloc((toas_tmp.nrtoas+toas_tmp.nrtoasdeleted)*sizeof(long double));
  toas->filenames = malloc((toas_tmp.nrtoas+toas_tmp.nrtoasdeleted)*sizeof(char *));
  toas->flags = malloc((toas_tmp.nrtoas+toas_tmp.nrtoasdeleted)*sizeof(char *));
  if(residuals.residy == NULL || toas->mjdx == NULL || residuals.nturn == NULL || residuals.nturn_fixed == NULL || residuals.resid == NULL || toas->mjd == NULL || toas->freqSite == NULL || toas->freqSSB == NULL || toas->filenames == NULL || toas->flags == NULL) {
    printerror(verbose.debug, "measureNuNudot: Cannot allocate memory (%d)\n", toas_tmp.nrtoas+toas_tmp.nrtoasdeleted);
    free(xstart);
    free(dx);
    free(xfit);
    free(fixed2);
    free(paramset2);
    return 0;
  }

  eph->wraps->nrwraps = 0;
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++)
    fixed2[i] = 1;
  fixed2[TEMPO3_PARAM_PHASE0] = 0;
  /* Only allow fitting for F0 and F1 */
  fixed2[TEMPO3_PARAM_F0] = 0;
  fixed2[TEMPO3_PARAM_F0+1] = 0;
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++)
    paramset2[i] = 0;
  paramset2[TEMPO3_PARAM_PEPOCH] = 1;
  paramset2[TEMPO3_PARAM_PHASE0] = 1;
  paramset2[TEMPO3_PARAM_F0] = 1;
  paramset2[TEMPO3_PARAM_F0+1] = 1;
  if(stridefit_includeF2 == 1)
    paramset2[TEMPO3_PARAM_F0+2] = 1;
  if(stridefit_includeF2 == 2)
    paramset2[TEMPO3_PARAM_NBRAKE] = 1;
  if(stridefit_includeF3)
    paramset2[TEMPO3_PARAM_F0+3] = 1;

  if(stridefit_scriptfilename[0] != 0) {
    fout_script = fopen(stridefit_scriptfilename, "w");
    if(fout_script == NULL) {
      printerror(verbose.debug, "measureNuNudot: Cannot open file %s\n", stridefit_scriptfilename);
      return 0;
    }
  }

  /* Now have to measure it */
  prev_toa = -1;
  prev_center = -1;
  for(i = 0; i < toas_tmp.nrtoas+toas_tmp.nrtoasdeleted; i++) {
    printf("Measuring nu and nudot toa %d/%d             \r", i+1, toas_tmp.nrtoas+toas_tmp.nrtoasdeleted);
    fflush(stdout);
    nu_mjd[i] = -1;
    // Check if this toa isn't too close to previous central toa to be considered
    if(toas_tmp.deleted[i] == 0 && (prev_toa < 0 || fabsl(toas_tmp.mjd[i]- prev_toa) >= stridefit_minstride)) {
      prev_toa = toas_tmp.mjd[i];
	/* Generate new toa list */
      toas->nrtoas = 0;
      toas->nrtoasdeleted = 0;
      toas->data_mjdmax = toas_tmp.data_mjdmin;   // NOTE: set max to min, so always gets updated with first TOA to be added
      toas->data_mjdmin = toas_tmp.data_mjdmax;
      //      printf("XXXXXX %Lf %Lf\n", toas->data_mjdmin, toas->data_mjdmax);
      for(j = 0; j < toas_tmp.nrtoas+toas_tmp.nrtoasdeleted; j++) {
	if(toas_tmp.deleted[j] == 0) {
	  if(toas_tmp.mjd[j] >= toas_tmp.mjd[i] - 0.5*stridefit_nrdays && toas_tmp.mjd[j] <= toas_tmp.mjd[i] + 0.5*stridefit_nrdays) {
	    ok = 1;
	    if(parfile_tmp.nrglitches >= 1) {
	      int glitchnr;
	      for(glitchnr = 0; glitchnr < parfile_tmp.nrglitches; glitchnr++) {
		if(parfile_tmp.parameter[TEMPO3_PARAM_GLEP+glitchnr] < toas_tmp.mjd[j] && parfile_tmp.parameter[TEMPO3_PARAM_GLEP+glitchnr] > toas_tmp.mjd[i])
		  ok = 0;
		if(parfile_tmp.parameter[TEMPO3_PARAM_GLEP+glitchnr] > toas_tmp.mjd[j] && parfile_tmp.parameter[TEMPO3_PARAM_GLEP+glitchnr] < toas_tmp.mjd[i])
		  ok = 0;
	      }
	    }
	    if(ok) {
	      toas->mjd[toas->nrtoas] = toas_tmp.mjd[j];
	      toas->freqSite[toas->nrtoas] = toas_tmp.freqSite[j];
	      toas->freqSSB[toas->nrtoas] = toas_tmp.freqSSB[j];
	      toas->mjdx[toas->nrtoas] = toas_tmp.mjdx[j];
	      toas->filenames[toas->nrtoas] = toas_tmp.filenames[j];
	      toas->flags[toas->nrtoas] = toas_tmp.flags[j];
	      toas->nrtoas++;
	      //	      printf("XXXXXX UPDATE %Lf %Lf with %Lf\n", toas->data_mjdmin, toas->data_mjdmax, toas_tmp.mjd[j]);
	      if(toas_tmp.mjd[j] > toas->data_mjdmax)
		toas->data_mjdmax = toas_tmp.mjd[j];
	      if(toas_tmp.mjd[j] < toas->data_mjdmin)
		toas->data_mjdmin = toas_tmp.mjd[j];
	      //	      printf("XXXXXX UPDATE %Lf %Lf\n", toas->data_mjdmin, toas->data_mjdmax);
	    }else {
	      /*	  printf("XXXX Excluded %Lf (%Lf) %d\n", toas_tmp.mjd[j], toas_tmp.mjd[i], parfile_tmp.nrglitches);*/
	    }
	  }
	}
      }
      residuals.mjdmin = toas->data_mjdmin;
      residuals.mjdmax = toas->data_mjdmax;

      /* Generate initial par file */
      initialise_ephemeris_par_only(&(eph->par));
      eph->par.track = enabletracking;

      /*    eph->par.phase0 = parfile_tmp.phase0; */
      eph->par.parameter[TEMPO3_PARAM_PHASE0] = 0;
      eph->par.parameter[TEMPO3_PARAM_F0] = nu[i];
      eph->par.parameter[TEMPO3_PARAM_F0+1] = nudot[i];
      if(verbose.debug) {
	fprintf(stderr, "\n  F0=%Lf, F1=%Le (initial guess)\n", eph->par.parameter[TEMPO3_PARAM_F0], eph->par.parameter[TEMPO3_PARAM_F0+1]);
      }
      if(stridefit_includeF2 == 1) {
	eph->par.parameter[TEMPO3_PARAM_F0+2] = 0;
      }
      if(stridefit_includeF2 == 2) {
	eph->par.parameter[TEMPO3_PARAM_NBRAKE] = 3;
      }
      if(stridefit_includeF3)
	eph->par.parameter[TEMPO3_PARAM_F0+3] = 0;
      eph->par.parameter[TEMPO3_PARAM_PEPOCH] = toas_tmp.mjd[i];
      /* Make sure turns will be subtracted correctly */
      
      /* I should probably centre pepoch */
      /* Note that nu and nudot are measured at toa, not at midpoint toa-range. Should also correct pdot, but whatever*/
      eph->par.parameter[TEMPO3_PARAM_PEPOCH] = 0.5*(toas->data_mjdmax + toas->data_mjdmin);
      eph->par.parameter[TEMPO3_PARAM_F0] = nu[i] + nudot[i]*(eph->par.parameter[TEMPO3_PARAM_PEPOCH] - toas_tmp.mjd[i])*TEMPO3_SecPerDay;
      eph->par.parameter[TEMPO3_PARAM_F0] = ephemeris_estimateNu_numerically_avoid_using(eph->par.parameter[TEMPO3_PARAM_PEPOCH], 1400, 0, 0, 0, "", "", &parfile_tmp, eph->paramset);
      eph->par.parameter[TEMPO3_PARAM_F0+1] = ephemeris_estimateNudot_numerically_avoid_using(eph->par.parameter[TEMPO3_PARAM_PEPOCH], 1400, 0, 0, 0, "", "", parfile_tmp, eph->paramset);
      if(verbose.debug) {
	printf("  F0=%Lf, F1=%Le (extrapolate pepoch)\n", eph->par.parameter[TEMPO3_PARAM_F0], eph->par.parameter[TEMPO3_PARAM_F0+1]);
      }

      /*
	if(eph->par.parameter[TEMPO3_PARAM_PEPOCH] > 54299 && eph->par.parameter[TEMPO3_PARAM_PEPOCH] < 54311) {
	fprintf(stderr, "XXXXXXXXXX\n");
	additional_funk_info.superverbose = 1;
	}else {
	fprintf(stderr, "YYYYYYYYYYYYY %Lf\n", eph->par.parameter[TEMPO3_PARAM_PEPOCH]);
	additional_funk_info.superverbose = -1;
	}*/

      if(toas->nrtoas < minnrtoas)
	eph->par.parameter[TEMPO3_PARAM_PEPOCH] = -1;
      if(toas->data_mjdmax - toas->data_mjdmin < stridefit_mintoarange)
	eph->par.parameter[TEMPO3_PARAM_PEPOCH] = -1;
      if(fabsl(eph->par.parameter[TEMPO3_PARAM_PEPOCH]- prev_center) < stridefit_minstride) {
	prev_center = eph->par.parameter[TEMPO3_PARAM_PEPOCH];
	eph->par.parameter[TEMPO3_PARAM_PEPOCH] = -1;
      }else {
	prev_center = eph->par.parameter[TEMPO3_PARAM_PEPOCH];
      }
      //      printf("MJDMAX=%Lf, MJDMIN=%Lf: diff = %Lf\n", toas->data_mjdmax, toas->data_mjdmin, toas->data_mjdmax - toas->data_mjdmin);

      if(additional_funk_info.superverbose != -1) {
	printf("Parfile for this piece of data:\n");
	verbose.indent = 10;

	ephemeris_def eph2;
	memcpy(&eph2, eph, sizeof(ephemeris_def));
	eph2.fixed = fixed2;    // So eph2 should be the same as eph, but pointing to fixed2 instead

	if(print_ephemeris(stdout, &eph2, NULL, toas, 0, 1, 0, 1, 0, 0, verbose) == 0) {
	  printerror(verbose.debug, "measureNuNudot: Showing ephemeris failed.");
	  return 0;
	}
	verbose.indent = 0;
      }


      if(eph->par.parameter[TEMPO3_PARAM_PEPOCH] > 0) {
	if(writeoutpieces) {
	  //	  printf("ACCEPTED PIECE WITH %d toas, which is more than %d\n", toas->nrtoas, minnrtoas);
	  sprintf(txt, "par%05d.par", i);
	  fout = fopen(txt, "w");
	  if(fout == NULL) {
	    printerror(verbose.debug, "measureNuNudot: Cannot open file %s\n", txt);
	    free(xstart);
	    free(dx);
	    free(xfit);
	    free(fixed2);
	    free(paramset2);
	    return 0;
	  }

	  ephemeris_def eph2;
	  memcpy(&eph2, eph, sizeof(ephemeris_def));
	  eph2.paramset = paramset2;    // So eph2 should be the same as eph, but pointing to paramset2 instead
	  if(print_ephemeris(fout, &eph2, NULL, toas, 0, 1, 1, 1, 0, 0, verbose) == 0) {
	    printerror(verbose.debug, "measureNuNudot: Showing ephemeris failed.");
	    return 0;
	  }

	  fclose(fout);
	  
	  sprintf(txt, "tim%05d.tim", i);
	  fout = fopen(txt, "w");
	  if(fout == NULL) {
	    printerror(verbose.debug, "measureNuNudot: Cannot open file %s\n", txt);
	    free(xstart);
	    free(dx);
	    free(xfit);
	    free(fixed2);
	    free(paramset2);
	    return 0;
	  }
	  for(k = 0; k < toas->nrtoas; k++) {
	    /* fprintf(fout, "%.30Le\t%Le\t%s\n", toas->mjd[k], toas->freq[k], toas->filenames[k]); */
	    if(strlen(toas->flags[k]) == 0) {
	      fprintf(fout, "%.30Le\t%Lf\t%Le\t%s\n", toas->mjd[k], toas->err[k]*eph->par.parameter[TEMPO3_PARAM_F0], toas->freqSSB[k], toas->filenames[k]);
	    }else {
	      fprintf(fout, "%.30Le\t%Lf\t%Le\t%s\t|%s\n", toas->mjd[k], toas->err[k]*eph->par.parameter[TEMPO3_PARAM_F0], toas->freqSSB[k], toas->filenames[k], toas->flags[k]);
	    }
	  }
	  fclose(fout);
	  if(stridefit_scriptfilename[0] != 0) {
	    // In some cases want to use -autoN option??? To check what get if N=3 is forced for example?
	    fprintf(fout_script, "tempo3 -nogtk -tol 1e-10 -auto -bary -autoquiet par%05d.par tim%05d.tim\n", i, i);
	  }
	}else {
	  /*    calcphases(parfile, *toas, 0, &residuals, paramlinked, paramset, 0);
		subtractTurns(*toas, &residuals, wraps, nosort, 0); 
		eph->par.parameter[TEMPO3_PARAM_PHASE0] -= residuals.resid[0];
		eph->par.parameter[TEMPO3_PARAM_PHASE0] = 0; */
	  calcphases(&(eph->par), *toas, 0, &residuals, eph->paramlinked, eph->paramset, 0, numthreads);
	  subtractTurns(*toas, &residuals, eph->wraps, nosort, 0);
	  dofit(0, &(eph->par), *toas, &residuals, eph->wraps, eph->paramlinked, eph->paramset, nosort, numthreads);
	  calcphases(&(eph->par), *toas, 0, &residuals, eph->paramlinked, eph->paramset, 0, numthreads);
	  subtractTurns(*toas, &residuals, eph->wraps, nosort, 0); 
	  
	  fill_xstart(&(eph->par), xstart);
	  fill_dx(&(eph->par), *toas, residuals, eph->wraps, fixed2, dx, additional_funk_info.superverbose, eph->paramlinked, eph->paramset, nosort, numthreads);
	  /* blaat */ 
	  if(additional_funk_info.superverbose != -1) {
	    pgplot_options.viewport.dontopen = 1;
	    pgplot_options.viewport.dontopen = 1;
	    strcpy(pgplot_options.box.xlabel, "MJD");
	    strcpy(pgplot_options.box.ylabel, "Resid [Phase]");
	    strcpy(pgplot_options.box.title, "");
	    int residctr;
	    for(residctr = 0; residctr < toas->nrtoas + toas->nrtoasdeleted; residctr++) {
	      residuals.residy[residctr] = residuals.resid[residctr];
	    }
	    pgplotGraph1(&pgplot_options, residuals.residy, toas->mjdx, NULL, toas->nrtoas, -1, -1, 0, residuals.mjdmin, residuals.mjdmax, 0, 0, 0, 0, 1, 1, 17, 1, 1, NULL, -1, noverbose);
	    printf("Press a key to continue\n");
	    fflush(stdout);
	    pgetch(); 
	  }
	  /*    if(doplot("?", ygraph, ygraph2, 0, 0, 0, 0, toas, &residuals, parfile, wraps, additional_funk_info.superverbose, paramset, 0, 1, paramlinked, title) == 0)
		return 0;
		pgetch();*/
	  additional_funk_info.toas = toas;
	  additional_funk_info.residuals = &residuals;
	  additional_funk_info.parfile = &(eph->par);
	  additional_funk_info.wraps = eph->wraps;
	  additional_funk_info.paramlinked = NULL;   /* We don't have any linked parameters in this simple parfile */
	  additional_funk_info.paramset    = NULL;   /* We don't use set information in this simple parfile (as at moment is only used for binary information) */
	  additional_funk_info.numthreads = numthreads;
	  additional_funk_info.nosort = nosort;
	  
	  doAmoeba_ld(amoeba_algorithm, xstart, dx, fixed2, xfit, &rms, TEMPO3_PARAM_TOTAL, &funk, ftol, &nfunk, 0, 0, 3, NULL, NULL);
	  
	  /* Refine */
	  eph->par.parameter[TEMPO3_PARAM_PHASE0] = xfit[TEMPO3_PARAM_PHASE0];
	  eph->par.parameter[TEMPO3_PARAM_F0] = xfit[TEMPO3_PARAM_F0];
	  eph->par.parameter[TEMPO3_PARAM_F0+1] = xfit[TEMPO3_PARAM_F0 + 1];
	  calcphases(&(eph->par), *toas, 0, &residuals, eph->paramlinked, eph->paramset, 0, numthreads);
	  subtractTurns(*toas, &residuals, eph->wraps, nosort, 0);
	  if(additional_funk_info.superverbose != -1) {    
	    pgplot_options.viewport.dontopen = 1;
	    pgplot_options.viewport.dontopen = 1;
	    strcpy(pgplot_options.box.xlabel, "MJD");
	    strcpy(pgplot_options.box.ylabel, "Resid [Phase]");
	    strcpy(pgplot_options.box.title, "");
	    int residctr;
	    for(residctr = 0; residctr < toas->nrtoas + toas->nrtoasdeleted; residctr++) {
	      residuals.residy[residctr] = residuals.resid[residctr];
	    }
	    pgplotGraph1(&pgplot_options, residuals.residy, toas->mjdx, NULL, toas->nrtoas, -1, -1, 0, residuals.mjdmin, residuals.mjdmax, 0, 0, 0, 0, 1, 1, 17, 1, 1, NULL, -1, noverbose);
	    pgetch(); 
	  }

	  /*
	    fill_xstart(parfile, xstart);
	    fill_dx(parfile, toas, residuals, wraps, fixed2, dx, additional_funk_info.superverbose, paramlinked, paramset, nosort, numthreads);
	    doAmoeba_ld(amoeba_algorithm, xstart, dx, fixed2, xfit, &rms, TEMPO3_PARAM_TOTAL, &funk, ftol*1e-6, &nfunk, 0, 0, 3, NULL, NULL);
	  */
	  
	  /*
	    if(doplot("?", ygraph, ygraph2, 0, 0, 0, 0, toas, &residuals, parfile, wraps, additional_funk_info.superverbose, paramset, 0, 1, paramlinked, fixed, 1, title, amoeba_algorithm) == 0)
	    return 0;
	    pgetch(); */
	  
	  nu[i] = xfit[TEMPO3_PARAM_F0];
	  nudot[i] = xfit[TEMPO3_PARAM_F0 + 1];
	  if(toas->nrtoas >= minnrtoas)
	    nu_mjd[i] = eph->par.parameter[TEMPO3_PARAM_PEPOCH];
	  else
	    nu_mjd[i] = -1;
	}
      }else {
	nu_mjd[i] = -1;
	//	printf("IGNORED PIECE\n");
      }
    }else {
      nu_mjd[i] = -1;
	//	printf("ALSO, OR DIFFERENTLY IGNORED PIECE?\n");
    }
  }

  if(stridefit_scriptfilename[0] != 0) {
    fclose(fout_script);
  }

  free(residuals.residy);
  free(residuals.resid);
  free(residuals.nturn);
  free(residuals.nturn_fixed);
  free(toas->mjdx);
  free(toas->mjd);
  free(toas->freqSite);
  free(toas->freqSSB);
  free(toas->flags);
  free(toas->filenames);

  memcpy(&residuals, &residuals_tmp, sizeof(residual_def));
  memcpy(toas, &toas_tmp, sizeof(toa_data_def));
  copy_parfile(&(eph->par), &parfile_tmp);
  //  memcpy(&parfile, &parfile_tmp, sizeof(tempo3_parameters_def));
  memcpy(eph->wraps, &wraps_tmp, sizeof(tempo3_wraps_def));

  printf("done                       \n");
  additional_funk_info.superverbose = 0;

  free_ephemeris_par_only(&parfile_tmp);
  free(xstart);
  free(dx);
  free(xfit);
  free(fixed2);
  free(paramset2);
  return 1;
}


void storeInStack(tempo3stack *stack, tempo3_parameters_def *parfile, toa_data_def timfile, tempo3_wraps_def *wraps, float previous_xstart, float previous_xend, int *fixed, int *paramset) 
{
  int i;
  /* Store last state in stack */
  if(stack->nrinstack == MaxNrStack) {
    for(i = 1; i < MaxNrStack; i++) {
      copy_parfile(&stack->parfiles[i-1], &stack->parfiles[i]);
      //      memcpy(&stack->parfiles[i-1], &stack->parfiles[i], sizeof(tempo3_parameters_def));
      memcpy(&stack->wraps[i-1], &stack->wraps[i], sizeof(tempo3_wraps_def));
      stack->previous_xstart[i-1] = stack->previous_xstart[i];
      stack->previous_xend[i-1] = stack->previous_xend[i];
      memcpy(stack->fixed[i-1], stack->fixed[i], TEMPO3_PARAM_TOTAL*sizeof(int));      
      memcpy(stack->paramset[i-1], stack->paramset[i], TEMPO3_PARAM_TOTAL*sizeof(int));      
      stack->nrtoas[i-1] = stack->nrtoas[i];
      stack->nrtoasdeleted[i-1] = stack->nrtoasdeleted[i];
      memcpy(stack->deleted[i-1], stack->deleted[i], stack->nrtoas[i] + stack->nrtoasdeleted[i]);      
    }
    stack->nrinstack--;
  }
  copy_parfile(&stack->parfiles[stack->nrinstack], parfile);
  //  memcpy(&stack->parfiles[stack->nrinstack], &parfile, sizeof(tempo3_parameters_def));
  memcpy(&stack->wraps[stack->nrinstack], wraps, sizeof(tempo3_wraps_def));
  stack->previous_xstart[stack->nrinstack] = previous_xstart;
  stack->previous_xend[stack->nrinstack] = previous_xend;
  memcpy(stack->fixed[stack->nrinstack], fixed, TEMPO3_PARAM_TOTAL*sizeof(int));
  memcpy(stack->paramset[stack->nrinstack], paramset, TEMPO3_PARAM_TOTAL*sizeof(int));
  stack->nrtoas[stack->nrinstack] = timfile.nrtoas;
  stack->nrtoasdeleted[stack->nrinstack] = timfile.nrtoasdeleted;
  memcpy(stack->deleted[stack->nrinstack], timfile.deleted, (timfile.nrtoas + timfile.nrtoasdeleted)*sizeof(int));      
  stack->nrinstack++;
}

void popstack(tempo3stack *stack, tempo3_parameters_def *parfile, toa_data_def *toas, tempo3_wraps_def *wraps, float *previous_xstart, float *previous_xend, int *fixed, int *paramset)
{
  int i, j;
  if(stack->nrinstack > 0) {
    copy_parfile(parfile, &stack->parfiles[stack->nrinstack-1]);
    //    memcpy(parfile, &stack->parfiles[stack->nrinstack-1], sizeof(tempo3_parameters_def));
    memcpy(wraps, &stack->wraps[stack->nrinstack-1], sizeof(tempo3_wraps_def));
    *previous_xstart = stack->previous_xstart[stack->nrinstack-1];
    *previous_xend = stack->previous_xend[stack->nrinstack-1];
    memcpy(fixed, stack->fixed[stack->nrinstack-1], TEMPO3_PARAM_TOTAL*sizeof(int));
    memcpy(paramset, stack->paramset[stack->nrinstack-1], TEMPO3_PARAM_TOTAL*sizeof(int));
    toas->nrtoas = stack->nrtoas[stack->nrinstack-1];
    toas->nrtoasdeleted = stack->nrtoasdeleted[stack->nrinstack-1];
    memcpy(toas->deleted, stack->deleted[stack->nrinstack-1], (toas->nrtoas + toas->nrtoasdeleted)*sizeof(int));      
    j = 0;
    for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
      if(toas->deleted[i] == 0) {
	if(j == 0) {
	  toas->data_mjdmin = toas->mjd[i];
	  toas->data_mjdmax = toas->mjd[i];
	  j = 1;
	}else {
	  if(toas->mjd[i] > toas->data_mjdmax)
	    toas->data_mjdmax = toas->mjd[i];
	  if(toas->mjd[i] < toas->data_mjdmin)
	    toas->data_mjdmin = toas->mjd[i];
	}
      }
    }
    /*
    printf("Restored: ");
    int i;
    for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
      printf("%d ", toas->deleted[i]);
    }
    printf("\n");
    toas->deleted[0] = 1;
    toas->deleted[1] = 1;
    toas->deleted[2] = 1;
    toas->deleted[3] = 1;
    toas->deleted[4] = 1;
    toas->deleted[5] = 1;
    toas->deleted[6] = 1;
    toas->deleted[7] = 1;
    toas->deleted[8] = 1;
    toas->deleted[9] = 1;*/

    /* Leave start situation in stack, so you can never overwrite that one */
    if(stack->nrinstack > 1)
      stack->nrinstack--;
  }else {
    printf("Stack is empty.\n");
  }
}

float timediff (struct timeval start, struct timeval stop)
{
  float diff;
  diff = stop.tv_sec - start.tv_sec;
  diff += 1e-6*(stop.tv_usec - start.tv_usec);
  return diff;
}



// spinfreq used for calculation of error in phase
void draw_residual(toa_data_def toas, int toanr, float x, float y, char *highlightedFlag, float markersize, int plottype, int mode, long double spinfreq)
{
  //      fprintf(stderr, "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXx DO HIGHLIGHTING '%s' == '%s'?\n", highlightedFlag, toas.flags[toanr]);
  //      fprintf(stderr, "  XXXXXXXXXXXXXXXXXXXXXXXXXXXx DO HIGHLIGHTING '%s'\n", highlightedFlag);
  //      fprintf(stderr, "  XXXXXXXXXXXXXXXXXXXXXXXXXXXx DO HIGHLIGHTING '%s'\n", toas.flags[toanr]);
  char *txtptr;
  txtptr = strstr(toas.flags[toanr], highlightedFlag);
  if(txtptr != NULL) {      // Given search string is part of set of flags
    if(txtptr[strlen(highlightedFlag)] == ' ' || txtptr[strlen(highlightedFlag)] == 0) { // Should be followed by a space (other flag is following) or 0 (end of string). Otherwise "-t TEL" will also find "-t TELESCOPE"
      ppgsch(3.5*markersize);
      ppgslw(10*markersize);
      ppgpt1(x, y, 17);
      ppgslw(1);
    }
  }else {
    ppgsch(1.5*markersize);
    ppgpt1(x, y, 17);
  }
  if(toas.freqSite[toanr] <= 0.01)	                                /* White is infinite frequency*/
    ppgsci(1);
  else if(toas.freqSite[toanr] < 1000 && toas.freqSite[toanr] > 0.01)            /* Red is low frequency */
    ppgsci(2);
  else if(toas.freqSite[toanr] > 2000)                              /* Blue is high frequency */
    ppgsci(4);
  else
    ppgsci(3);                                              /* Green is mid frequencies */
  ppgsch(1.0*markersize);
  ppgpt1(x, y, 17);
  if((plottype == 0 || plottype == PLOTTYPE_ORBITALPHASE || plottype == PLOTTYPE_FREQUENCY) && mode != 0) {
    ppgerr1(6, x, y, toas.err[toanr]*1e-6*spinfreq, 1);
  }
}


// inclination is in deg
// mpsr = mass pulsar is solar mass units
// If tmin = tmax = -1: one or two full cycles are plotted
// If tmin = tmax = -2: at least one (or two) full cycle is plotted, but also includes full toa range
void plotOrbit(tempo3_parameters_def *parfile, int *paramset, toa_data_def *toas, long double t, long double tmin, long double tmax, long double inclination, long double mpsr, int hierarchical, int markersize, char *device, char *orbit_dump_prefix, int finderrors, long double *dplus)
{
  long double dummy1, dummy2, xp, yp, xc, yc, inner_mass[TEMPO3_MaxNrCompanions], companionmass;
  float *xgraph, *ygraph, *timearray, range_min, range_max, range_min2, range_max2;
  int companion, dofullcycle, i, idx, deviceID, makeplot, c1, c2, plottype, curc;
  int nrcompanions;
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);
  verbose_definition verbose;
  cleanVerboseState(&verbose);

  int nrpoints_line = 30000;        // Nr of points calculated to make a line drawing
  int larger_points_in_orbit = 100;  // Every so many points in the orbit are larger to see how fast the movement is


  printf("\nCalculate binary position for MJD %Lf assuming an orbital inclination or %Lf deg\n", t, inclination);
  printf("\nIn the residual plots: the green points mark the times for which there is a TOA defined. The red dot indicates the position at the requested time. \n");
  printf("\nIn orbit plots: The origin of the coordinate system is the centre of mass. The black dots are equally spaced in time. The red dot indicates the position at the specified time. The blue cross markes periastron and the blue plus symbol the ascending node (object moves in the positive y direction at that point). If periastron and the ascending node overlap you see a star, which should happen for circular orbits.\n");
  printf("\nThe \"projected on l.o.s.\" plots show the shape of the orbit, but rescaled such that the actually measured radial component is shown. The \"in orbital plane\" plots are identical (since the shape is the same in the orbital plane), but the y-axis is re-labeled. By dividing by sin(i) the actual orbit is obtained.\n");
  nrcompanions = 0;
  for(companion = 0; companion < TEMPO3_MaxNrCompanions; companion++) {
    if(parfile->companion_defined[companion]) {
      nrcompanions++;
    }
  }
  if(nrcompanions == 0) {
    printwarning(0, "WARNING plotOrbit: No binary orbit appears to be defined");
    return;
  }else if(nrcompanions > 1) {
    if(hierarchical == 0) {
      printwarning(0, "WARNING plotOrbit: There are multiple companions. Each companion is treated independently as if the others do not exist. This is completely wrong for all properties of the companion(s), except for the most inner companion in the limit of the system being an hierarchical system. Consider using the -hierarchical option, if applicable. In addition, it is assumed that all orbits are coplanar (identical inclination angles).");
    }else {
      printwarning(0, "WARNING: There are multiple orbits defined and they are treated as a hierarchinal system: all interactions between companions are ignored. In addition, it is assumed that all orbits are coplanar (identical inclination angles).");
    }
  }

  ppgqid(&deviceID);

  i = nrpoints_line;
  if(toas != NULL) {
    if(toas->nrtoas > i)
      i = toas->nrtoas;
  }
  xgraph = malloc(i*sizeof(float));
  ygraph = malloc(i*sizeof(float));
  timearray = malloc(i*sizeof(float));
  if(xgraph == NULL || ygraph == NULL || timearray == NULL) {
    printerror(0, "ERROR plotOrbit(): Cannot allocate memory.");
    exit(0);
  }

  inclination *= M_PI/180.0;
  if(tmin < 0 && tmax < 0) {
    if(fabsl(tmin) > 0.9 && fabsl(tmin) < 1.1) {
      dofullcycle = 1;
    }else if(fabsl(tmin) > 1.9 && fabsl(tmin) < 2.1) {
      dofullcycle = 2;
    }
  }else {
    dofullcycle = 0;
  }

  // Start with plotting effects on the pulsar for individual pulsars
  pgplot_options.box.drawlabels = 1;
  pgplot_options.box.drawtitle = 1;
  pgplot_options.viewport.dontclose = 1;
  pgplot_options.viewport.aspectratio = 1.0;
  strcpy(pgplot_options.viewport.plotDevice, device);


  // Derive the masses to be used by the respective mass functions of the different companions
  inner_mass[0] = mpsr;
  for(companion = 0; companion < TEMPO3_MaxNrCompanions; companion++) {
    if(parfile->companion_defined[companion]) {
      if(hierarchical != 0 && companion > 0) {
	derive_orbital_parameters(parfile, paramset, companion-1, parfile->parameter[TEMPO3_PARAM_T0+companion-1], inclination, inner_mass[companion-1], 0, &xp, &yp, &xc, &yc, &companionmass, 0, dplus);
	inner_mass[companion] = inner_mass[companion-1] + companionmass;
      }else {
	inner_mass[companion] = mpsr;
      }
    }
  }

  // plottype == 0 = resdual plot
  // plottype == 1 = psr orbit in plane sky plot
  // plottype == 2 = psr orbit plot
  // plottype == 3 = companion orbit plot
  for(plottype = 0; plottype < 4; plottype++) {
    range_min2 = range_max2 = 0;
    // if companion == TEMPO3_MaxNrCompanions -> combined effect  - So companion is the current companion to consider in the plots
    for(companion = 0; companion <= TEMPO3_MaxNrCompanions; companion++) {
      // See if plot can be made for this companion
      makeplot = 0;
      if(companion == TEMPO3_MaxNrCompanions && nrcompanions > 1) {  // Only make "total" plots if there are multiple companions
	makeplot = 1;
      }else {
	if(parfile->companion_defined[companion])
	  makeplot = 1;
      }
      if(makeplot) {
	// Find range to be plotted if range was not set
	if(companion < TEMPO3_MaxNrCompanions) {
	  c1 = companion;
	  c2 = companion;
	}else {
	  c1 = 0;
	  c2 = TEMPO3_MaxNrCompanions-1;
	}
	for(curc = c1; curc <= c2; curc++) { // Loop over all companions + one extra, which is total effect (if multiple companions are defined)
	  if(parfile->companion_defined[curc]) {
	    if(dofullcycle == 1 || dofullcycle == 2) {   // If wants to cover full orbit if data-range is limited, adjust tmin and tmax
	      if(parfile->parameter[TEMPO3_PARAM_ECC+curc] >= 1.0) {
		dummy1 = parfile->parameter[TEMPO3_PARAM_T0+curc]-20.0*parfile->parameter[TEMPO3_PARAM_PB+curc];
		dummy2 = dummy1 + 40.0*parfile->parameter[TEMPO3_PARAM_PB+curc];
	      }else {
		dummy1 = parfile->parameter[TEMPO3_PARAM_T0+curc]-parfile->parameter[TEMPO3_PARAM_PB+curc];
		dummy2 = dummy1 + 2.0*parfile->parameter[TEMPO3_PARAM_PB+curc];
	      }
	      if(curc == c1 || (dummy1 < tmin)) {
		tmin = dummy1;
	      }
	      if(curc == c1 || (dummy2 > tmax)) {
		tmax = dummy2;
	      }
	    }
	  }
	}
	if(dofullcycle == 2) {  // Make sure more than one orbit is covered to fit in all data points
	  if(toas != NULL) {
	    for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
	      if(toas->deleted[i] == 0) {
		if(toas->mjdx[i] < tmin)
		  tmin = toas->mjdx[i];
		if(toas->mjdx[i] > tmax)
		  tmax = toas->mjdx[i];
	      }
	    }
	  }
	}
	// Check if orbits are resolved
	for(curc = c1; curc <= c2; curc++) {
	  if(parfile->companion_defined[curc]) {
	    if((double)nrpoints_line/(double)((tmax-tmin)/(double)parfile->parameter[TEMPO3_PARAM_PB+curc]) < 10) {
	      printwarning(0, "WARNING plotOrbit(): Orbit of companion %d is barely or not at all resolved with the resolution of at least some of the plots", curc+1);
	    }
	  }
	}
	// Make the timing effect plot
	if(plottype == 0) {
	  strcpy(pgplot_options.box.ylabel, "TOA delay [s]");
	  // Draw line of timing effect as function of time
	  for(i = 0; i < nrpoints_line; i++) {
	    xgraph[i] = (tmax-tmin)*(long double)i/(long double)(nrpoints_line-1);
	    if(companion == TEMPO3_MaxNrCompanions) {
	      ygraph[i] = -ephemeris_get_torb((long double)(tmin + xgraph[i])*TEMPO3_SecPerDay, parfile, paramset, 0, -1);
	    }else {
	      ygraph[i] = -ephemeris_get_torb((long double)(tmin + xgraph[i])*TEMPO3_SecPerDay, parfile, paramset, 0, companion);
	    }
	    // Centre w.r.t. T0 if the full cycle is plotted, otherwise it is centred w.r.t. tmin
	    if(dofullcycle && companion != TEMPO3_MaxNrCompanions) {
	      //	      printf("XXXXX %f %Lf %Lf\n", xgraph[i], tmin,  parfile->parameter[TEMPO3_PARAM_T0+companion]);
	      xgraph[i] += tmin -  parfile->parameter[TEMPO3_PARAM_T0+companion];
	      sprintf(pgplot_options.box.xlabel, "t-%Lf [days]", parfile->parameter[TEMPO3_PARAM_T0+companion]);
	    }else {
	      sprintf(pgplot_options.box.xlabel, "t-%Lf [days]", tmin);
	    }
	  }
	  if(companion == TEMPO3_MaxNrCompanions) {
	    sprintf(pgplot_options.box.title, "Total effect");
	  }else {
	    sprintf(pgplot_options.box.title, "Effect of companion %d", companion+1);
	  }
	  if(strlen(orbit_dump_prefix) > 0) {
	    char filename[MaxFilenameLength];
	    FILE *fout;
	    sprintf(filename, "%s.predicteddelay%d", orbit_dump_prefix, companion);
	    printf("Output graph to %s\n", filename);
	    fout = fopen(filename, "w");
	    if(fout == NULL) {
	      printwarning(0, "WARNING: Cannot open %s for writing.", filename);
	    }else {
	      fprintf(fout, "# %s\n", pgplot_options.box.title);
	      fprintf(fout, "# Predicted delay (%s) as function of time (%s)\n", pgplot_options.box.ylabel, pgplot_options.box.xlabel);
	      for(i = 0; i < nrpoints_line; i++) {
		fprintf(fout, "%e %e\n", xgraph[i], ygraph[i]);
	      }
	      fclose(fout);
	    }
	  }
	  //int pgplotGraph1(&pgplot_viewport_def viewport, float *data, float *datax, float *sigma, int nrx, float xmin, float xmax, int dontsetranges, float xmin_show, float xmax_show, float ymin_show, float ymax_show, int forceMinZero, pgplot_box_def pgplotbox, int hist, int noline, int pointtype, int color, int boxcolor, regions_definition *regions, verbose_definition verbose);
	  pgplot_options.viewport.noclear = 0;
	  if(pgplotGraph1(&pgplot_options, ygraph, xgraph, NULL, nrpoints_line, -1, -1, 0, -1, -1, -1, -1, 0, 0, 0, 1, 0, 1, 1, NULL, -1, verbose) == 0) {
	    printerror(0, "ERROR plotOrbit(): Cannot plot graph.");
	    exit(0);
	  }

	  // Draw specified time
	  ppgsch(1.5*markersize);
	  ppgsci(2);
	  if(companion == TEMPO3_MaxNrCompanions) {
	    if(dofullcycle && companion != TEMPO3_MaxNrCompanions) {
	      ppgpt1(t-parfile->parameter[TEMPO3_PARAM_T0+companion], -ephemeris_get_torb((long double)(t)*TEMPO3_SecPerDay, parfile, paramset, 0, -1), 17);
	    }else {
	      ppgpt1(t-tmin, -ephemeris_get_torb((long double)(t)*TEMPO3_SecPerDay, parfile, paramset, 0, -1), 17);
	    }
	  }else {
	    if(dofullcycle && companion != TEMPO3_MaxNrCompanions) {
	      ppgpt1(t-parfile->parameter[TEMPO3_PARAM_T0+companion], -ephemeris_get_torb((long double)(t)*TEMPO3_SecPerDay, parfile, paramset, 0, companion), 17);
	    }else {
	      ppgpt1(t-tmin, -ephemeris_get_torb((long double)(t)*TEMPO3_SecPerDay, parfile, paramset, 0, companion), 17);
	    }
	  }
	  ppgsci(1);
	  ppgsch(1.0);

	  // Draw points corresponding with TOA's
	  pgplot_options.viewport.noclear = 1;
	  pgplot_options.viewport.dontopen = 1;
	  if(toas != NULL) {
	    idx = 0;
	    for(i = 0; i < toas->nrtoas + toas->nrtoasdeleted; i++) {
	      if(toas->deleted[i] == 0) {
		xgraph[idx] = toas->mjdx[i] - tmin;
		if(companion == TEMPO3_MaxNrCompanions) {
		  ygraph[idx] = -ephemeris_get_torb((long double)(tmin + xgraph[i])*TEMPO3_SecPerDay, parfile, paramset, 0, -1);
		}else {
		  ygraph[idx] = -ephemeris_get_torb((long double)(tmin + xgraph[i])*TEMPO3_SecPerDay, parfile, paramset, 0, companion);
		}
		idx++;
	      }
	      if(dofullcycle && companion != TEMPO3_MaxNrCompanions) {
		xgraph[i] += tmin -  parfile->parameter[TEMPO3_PARAM_T0+companion];
	      }
	    }
	    if(pgplotGraph1(&pgplot_options, ygraph, xgraph, NULL, toas->nrtoas, -1, -1, 1, -1, -1, -1, -1, 0, 0, 1, 1, 17, 3, 1, NULL, -1, verbose) == 0) {
	      printerror(0, "ERROR plotOrbit(): Cannot plot graph.");
	      exit(0);
	    }
	    if(strlen(orbit_dump_prefix) > 0) {
	      char filename[MaxFilenameLength];
	      FILE *fout;
	      sprintf(filename, "%s.toas_predicteddelay%d", orbit_dump_prefix, companion);
	      printf("Output graph to %s\n", filename);
	      fout = fopen(filename, "w");
	      if(fout == NULL) {
		printwarning(0, "WARNING: Cannot open %s for writing.", filename);
	      }else {
		fprintf(fout, "# %s\n", pgplot_options.box.title);
		fprintf(fout, "# Predicted delay at TOA's (%s) as function of time (%s)\n", pgplot_options.box.ylabel, pgplot_options.box.xlabel);
		for(i = 0; i < toas->nrtoas; i++) {
		  fprintf(fout, "%e %e\n", xgraph[i], ygraph[i]);
		}
		fclose(fout);
	      }
	    }
	  }
	}else if(plottype == 1 || plottype == 2 || plottype == 3) {
	  if(companion == TEMPO3_MaxNrCompanions) {
	    if(plottype == 1)
	      sprintf(pgplot_options.box.title, "Pulsar orbit projected on l.o.s.");
	    else if(plottype == 2)
	      sprintf(pgplot_options.box.title, "Pulsar orbit in orbital plane");
	    else if(plottype == 3)
	      sprintf(pgplot_options.box.title, "Companion orbits (copy of prev. plot?)");
	  }else {
	    if(plottype == 1)
	      sprintf(pgplot_options.box.title, "Effect of companion %d on pulsar projected on l.o.s.", companion+1);
	    else if(plottype == 2)
	      sprintf(pgplot_options.box.title, "Effect of companion %d on pulsar in orbital plane", companion+1);
	    else if(plottype == 3)
	      sprintf(pgplot_options.box.title, "Orbit of companion %d (if inclination is %.2Lf deg)", companion+1, inclination*180.0/M_PI);
	  }
	  if(plottype == 1) {
	    sprintf(pgplot_options.box.ylabel, "Line of sight distance away from COM [lt-s]");
	    sprintf(pgplot_options.box.xlabel, "Position in plane sky projected on line of nodes [lt-s] \\(2235) sin(i)");
	  }else if(plottype == 2) {
	    sprintf(pgplot_options.box.ylabel, "Y [lt-s] \\(2235) sin(i)");
	    sprintf(pgplot_options.box.xlabel, "X (position projected on line of nodes) [lt-s] \\(2235) sin(i)");
	  }else if(plottype == 3) {
	    sprintf(pgplot_options.box.ylabel, "Y [lt-s]");
	    sprintf(pgplot_options.box.xlabel, "X (position projected on line of nodes) [lt-s]");
	  }

	  // Plot orbit
	  for(i = 0; i < nrpoints_line; i++) {
	    xgraph[i] = 0;
	    ygraph[i] = 0;
	  }
	  for(curc = c1; curc <= c2; curc++) {
	    range_min = range_max = 0;
	    if(parfile->companion_defined[curc]) {  // Make plots for individual companions
	      for(i = 0; i < nrpoints_line; i++) {
		dummy1 = (tmax-tmin)*(long double)i/(long double)(nrpoints_line-1);
		timearray[i] = dummy1+tmin;
		derive_orbital_parameters(parfile, paramset, curc, dummy1+tmin, inclination, inner_mass[curc], 0, &xp, &yp, &xc, &yc, &companionmass, finderrors, dplus);
		if(plottype == 3) {
		  xgraph[i] = xc;
		  ygraph[i] = yc;
		}else {
		  xgraph[i] += xp;
		  ygraph[i] += yp;
		}
	      }
	      for(i = 0; i < nrpoints_line; i++) {
		if(fabs(xgraph[i]) > range_max) {
		  range_min = -fabs(xgraph[i]);
		  range_max = fabs(xgraph[i]);
		}
	      }
	      if(range_min < range_min2)
		range_min2 = range_min;
	      if(range_max > range_max2)
		range_max2 = range_max;
	      pgplot_options.viewport.noclear = 0;
	      if(plottype == 3) {
		if(c1 != c2) {   // Combined plot - remember largest size
		  if(curc != c1)
		    pgplot_options.viewport.noclear = 1;  // overplot next orbit
		  range_min = range_min2;
		  range_max = range_max2;
		}
		if(pgplotGraph1(&pgplot_options, ygraph, xgraph, NULL, nrpoints_line, range_min*1.1, range_max*1.1, pgplot_options.viewport.noclear, range_min*1.1, range_max*1.1, range_min*1.1, range_max*1.1, 0, 0, 1, 1, -1, 1, 1, NULL, -1, verbose) == 0) {
		  printerror(0, "ERROR plotOrbit(): Cannot plot graph.");
		  exit(0);
		}
		// Mark equidistant times along orbit
		ppgsch(1.0*markersize);
		for(i = 0; i < nrpoints_line; i++) {
		  if(i%larger_points_in_orbit == 0) {
		    ppgpt1(xgraph[i], ygraph[i], 17);
		  }
		}
		// Mark periastron
		if(c1 == c2) {
		  derive_orbital_parameters(parfile, paramset, curc, -1, inclination, inner_mass[curc], 0, &xp, &yp, &xc, &yc, &companionmass, finderrors, dplus);
		  ppgsch(1.5*markersize);
		  ppgsci(4);
		  ppgpt1(xc, yc, 5);
		  ppgsch(1.0*markersize);
		  //		  float xch, ych;
		  //		  ppgqcs(4, &xch, &ych);
		  //		  ppgptxt(xp-0.25*ych, yp-0.25*ych, 0.0, 1.0, "P");
		}
		// Mark ascending node
		if(c1 == c2) {
		  derive_orbital_parameters(parfile, paramset, curc, -2, inclination, inner_mass[curc], 0, &xp, &yp, &xc, &yc, &companionmass, finderrors, dplus);
		  ppgsch(1.5*markersize);
		  ppgsci(4);
		  ppgpt1(xc, yc, 2);
		  ppgsch(1.0*markersize);
		  ppgsci(1);
		  //		  float xch, ych;
		  //		  ppgqcs(4, &xch, &ych);
		  //		  ppgptxt(xp-0.25*ych, yp-0.25*ych, 0.0, 1.0, "P");
		}
		// Mark specified time
		derive_orbital_parameters(parfile, paramset, curc, t, inclination, inner_mass[curc], 0, &xp, &yp, &xc, &yc, &companionmass, finderrors, dplus);
		ppgsch(1.5*markersize);
		ppgsci(2);
		ppgpt1(xc, yc, 17);
	       
	      }
	    }
	  }

	  if(strlen(orbit_dump_prefix) > 0  && (plottype == 2 || (plottype == 3 && companion != TEMPO3_MaxNrCompanions))) {
	    FILE *fout;
	    char filename[MaxFilenameLength];
	    if(plottype == 2) {
	      if(companion != TEMPO3_MaxNrCompanions)
		sprintf(filename, "%s.projPSR%d", orbit_dump_prefix, companion+1);
	      else
		sprintf(filename, "%s.projPSRall", orbit_dump_prefix);
	    }else {
	      if(companion != TEMPO3_MaxNrCompanions)
		sprintf(filename, "%s.companion%d", orbit_dump_prefix, companion+1);
	      else
		sprintf(filename, "%s.companionall", orbit_dump_prefix);
	    }
	    printf("Output graph to %s\n", filename);
	    fout = fopen(filename, "w");
	    if(fout == NULL) {
	      printwarning(0, "WARNING: Cannot open %s for writing.", filename);
	    }else {
	      fprintf(fout, "# %s\n", pgplot_options.box.title);
	      if(plottype == 2) {
		fprintf(fout, "# MJD X*sin(i) Y*sin(i)\n");
	      }else {
		fprintf(fout, "# MJD X Y\n");
	      }
	      for(i = 0; i < nrpoints_line; i++) {
		fprintf(fout, "%e %e %e\n", timearray[i], xgraph[i], ygraph[i]);
	      }
	      fclose(fout);
	    }
	  }

	  if(plottype != 3) {
	    if(c1 != c2) {   // Combined plot - remember largest size
	      range_min = range_min2;
	      range_max = range_max2;
	    }
	    if(pgplotGraph1(&pgplot_options, ygraph, xgraph, NULL, nrpoints_line, range_min*1.1, range_max*1.1, 0, range_min*1.1, range_max*1.1, range_min*1.1, range_max*1.1, 0, 0, 1, 1, -1, 1, 1, NULL, -1, verbose) == 0) {
	      printerror(0, "ERROR plotOrbit(): Cannot plot graph.");
	      exit(0);
	    }
	    // Mark equidistant times along orbit
	    ppgsch(1.0*markersize);
	    for(i = 0; i < nrpoints_line; i++) {
	      if(i%larger_points_in_orbit == 0) {
		ppgpt1(xgraph[i], ygraph[i], 17);
	      }
	    }
	    // Mark periastron and the ascending node
	    for(curc = c1; curc <= c2; curc++) {
	      if(parfile->companion_defined[curc]) {
		if(c1 == c2) {
		  derive_orbital_parameters(parfile, paramset, curc, -1, inclination, inner_mass[curc], 0, &xp, &yp, &xc, &yc, &companionmass, finderrors, dplus);
		  ppgsch(1.5*markersize);
		  ppgsci(4);
		  ppgpt1(xp, yp, 5);
		  derive_orbital_parameters(parfile, paramset, curc, -2, inclination, inner_mass[curc], 0, &xp, &yp, &xc, &yc, &companionmass, finderrors, dplus);
		  ppgpt1(xp, yp, 2);
		  ppgsch(1.0*markersize);
		  ppgsci(1);
		  //		  float xch, ych;
		  //		  ppgqcs(4, &xch, &ych);
		  //		  ppgptxt(xp-0.25*ych, yp-0.25*ych, 0.0, 1.0, "P");
		}
	      }
	    }
		// Mark 
		if(c1 == c2) {
		  //		  float xch, ych;
		  //		  ppgqcs(4, &xch, &ych);
		  //		  ppgptxt(xp-0.25*ych, yp-0.25*ych, 0.0, 1.0, "P");
		}

	    // Mark specified time
	    xgraph[0] = ygraph[0] = 0;
	    for(curc = c1; curc <= c2; curc++) {
	      if(parfile->companion_defined[curc]) {
		if(c1 == c2 && plottype == 1)
		  derive_orbital_parameters(parfile, paramset, curc, t, inclination, inner_mass[curc], 1, &xp, &yp, &xc, &yc, &companionmass, finderrors, dplus);
		else
		  derive_orbital_parameters(parfile, paramset, curc, t, inclination, inner_mass[curc], 0, &xp, &yp, &xc, &yc, &companionmass, finderrors, dplus);
		xgraph[0] += xp;
		ygraph[0] += yp;
	      }
	    }
	    ppgsch(1.5*markersize);
	    ppgsci(2);
	    ppgpt1(xgraph[0], ygraph[0], 17);
	  }

	  // Bark barycentre
	  ppgsci(1);
	  ppgpt1(0, 0, 2);
	  ppgsch(1.0*markersize);
	}
	i = pgplot_device_type(NULL, verbose); 
	if(i <= 2) {
	  printf("Press a key to continue (in terminal)\n");
	  fflush(stdout);
	  pgetch();
	}
      }  // End of if statement if makeplot is set
    }  // End loop over companion number
  }   // End loop over plottype

  free(xgraph);
  free(ygraph);
  free(timearray);
  ppgclos();
  ppgslct(deviceID);
  return;
}



// If mjd == -1, the periastron position is returned
// If mjd == -2, the ascending node position is returned
// innermass is either the mass of the pulsar, or the mass of all inner bodies (if -hierarchical is specified)
void derive_orbital_parameters(tempo3_parameters_def *parfile, int *paramset, int companion, long double mjd, long double inclination, long double innermass, int doprint, long double *xpsini, long double *ypsini, long double *xc, long double *yc, long double *companionmass, int showerrors, long double *dplus)
{
  long double pb, pb_err, t0, pbdot, ecc, ecc_err, asini_psr, asini_psr_err, omega_psr, omegaB, omegaB_err, M, dummy, phase, E, AT, rpsini, fn, fn_err, m2, m2_err, atot, atot_err, a2, a2_err, r2, epoch_an, specific_energy, vrel, mjd1, mjd2, dt, xc_old, yc_old, xpsini_old, ypsini_old, rtot, v2, vpsr, value, err;
  double GM          = 1.3271243999e20;      /* Gravitational constant * mass sun */
  if(parfile->companion_defined[companion]) {

    // Input parameters
    pb = 0;
    pb_err = 0;
    if(paramset[TEMPO3_PARAM_PB+companion]) {
      pb = parfile->parameter[TEMPO3_PARAM_PB+companion];       // days
      if(showerrors)
	pb_err = dplus[TEMPO3_PARAM_PB+companion];
    }
    if(paramset[TEMPO3_PARAM_T0+companion])
      t0 = parfile->parameter[TEMPO3_PARAM_T0+companion];       // days
    else
      t0 = 0;
    if(paramset[TEMPO3_PARAM_PBDOT+companion])
      pbdot = parfile->parameter[TEMPO3_PARAM_PBDOT+companion]; // s/s = days/day
    else
      pbdot = 0;
    if(paramset[TEMPO3_PARAM_ECC+companion]) {
      ecc = parfile->parameter[TEMPO3_PARAM_ECC+companion];     // s/s = days/day
      if(showerrors)
	ecc_err = dplus[TEMPO3_PARAM_ECC+companion];
    }else {
      ecc = 0;
      ecc_err = 0;
    }
    if(paramset[TEMPO3_PARAM_A1+companion]) {
      asini_psr = parfile->parameter[TEMPO3_PARAM_A1+companion];       // lt-s // Really p
      if(showerrors)
	asini_psr_err = dplus[TEMPO3_PARAM_A1+companion];
    }else {
      asini_psr = 0;
      asini_psr_err = 0;
    }
    if(paramset[TEMPO3_PARAM_OM+companion])
      omega_psr = parfile->parameter[TEMPO3_PARAM_OM+companion]*M_PI/180.0;    // rad
    else
      omega_psr = 0;

    omegaB = 2.0*M_PI/pb;             // Mean angular frequency [rad/day]
    if(showerrors)
      omegaB_err = omegaB*pb_err/pb;

    // Derive the epoch of the ascending node
    AT = -omega_psr; // Location of the ascending node, since periastron (AT=0) happens per definition omega_psr later than the ascending node
    if(ecc < 1.0) {
      epoch_an = pb*( 2.0*atanl(sqrtl((1.0-ecc)/(1.0+ecc))*tanl(0.5*AT))- ecc*sqrtl(1.0-ecc*ecc)*sinl(AT)/(1.0+ecc*cosl(AT)))/(2.0*M_PI) + t0;
    }else if(ecc > 1.0) {
      epoch_an = pb*(   (ecc*sqrtl(ecc*ecc-1.0)*sinl(AT))/(1.0+ecc*cosl(AT)) -logl((sqrtl(ecc+1.0)+sqrtl(ecc-1.0)*tanl(0.5*AT))/(sqrtl(ecc+1.0)-sqrtl(ecc-1.0)*tanl(0.5*AT))) )/(2.0*M_PI) + t0;
    }else {
      epoch_an = tanl(0.5*AT);
      epoch_an = 0.5*epoch_an+epoch_an*epoch_an*epoch_an/6.0;  // Barker's equation, which gives mean anomaly
      epoch_an = epoch_an/omegaB+t0;
    }
    if(mjd < 0) {
      if(fabsl(mjd) < 1.1 && fabsl(mjd) > 0.9) {
	mjd = t0; // Periastron per definition occurs at t0
      }else if(fabsl(mjd) < 2.1 && fabsl(mjd) > 1.9) {
	mjd = epoch_an;
      }
    }

    // Do calculation at two slightly different times if doprint is set, so orbital velocity can be measured to do some additional checks
    mjd1 = mjd;
    mjd2 = mjd;
    dt = 1e-8*pb;
    if(doprint) {
      mjd1 -= dt;
    }
    for(mjd = mjd1; mjd <= mjd2+0.01*dt; mjd += dt) {
      // Derived parameters
      M = omegaB*(mjd - t0);            // Mean anomaly [rad]
      if(ecc < 1.0) {
	M = derotate_rad_longdouble(M);   // Make the Mean anomaly a number between 0 and 2pi
	dummy = (mjd - t0)/pb;
	phase = M - M_PI*pbdot*dummy*dummy;
	phase = derotate_rad_longdouble(phase);
      }else {
	phase = M;  // For parabolic/hyperbolic orbits M is non-periodic
      }
      // Solve E-ecc*sin(E) = phase
      E = getEccentricAnomaly(phase, ecc);                      // Eccentric anomaly [rad]
      if(ecc < 1.0)
	AT = 2.0*atanl(sqrtl((1.0+ecc)/(1.0-ecc))*tanl(0.5*E));   // True anomaly [rad]
      else if(ecc == 1.0)
	AT = E;   // That is what function returns
      else
	AT = 2.0*atanl(sqrtl((ecc+1.0)/(ecc-1.0))*tanhl(0.5*E));   // True anomaly [rad]
      // Distance CM -> PSR
      if(ecc < 1.0) {
	//    rpsini = asini_psr*(1.0-ecc*cos(E));    // Should do the same as the next equation
	rpsini = asini_psr*(1.0-ecc*ecc)/(1.0+ecc*cos(AT));
      }else if(ecc == 1.0) {
	rpsini = asini_psr/(1.0+cos(AT));  // asini_psr really psini
      }else {
	rpsini = asini_psr*(ecc*ecc-1.0)/(1.0+ecc*cos(AT));
      }
      
      // Mass function
      fn  = omegaB*omegaB*asini_psr*asini_psr*asini_psr/GM;
      fn *= TEMPO3_SPEEDOFLIGHT*TEMPO3_SPEEDOFLIGHT*TEMPO3_SPEEDOFLIGHT/(TEMPO3_SecPerDay*TEMPO3_SecPerDay);
      fn_err = 0;
      // Solve mass function for the companion mass
      m2 = solve_mass_function_m2(fn, sin(inclination), innermass);
      m2_err = 0;
      if(showerrors) {
	fn_err = sqrtl(fn*fn*(9.0*asini_psr_err*asini_psr_err/(asini_psr*asini_psr) + 4.0*pb_err*pb_err/(pb*pb))); // Also used above
	m2_err = m2*(innermass+m2)*fn_err/(fn*(3.0*innermass+m2));  // Also used above
      }
      *companionmass = m2; // Pass information on to whoever called this function
	  
      // Sum of semi-major axes between pulsar and companion follows from 3rd law of Kepler
      // Note that this is the actual value (without a sin(i)), since i was already used to get m2.
      dummy = GM*(innermass+m2)/(omegaB*omegaB);
      dummy *= TEMPO3_SecPerDay*TEMPO3_SecPerDay;            // Convert in SI units
      atot = powl(dummy, 1.0/3.0)/TEMPO3_SPEEDOFLIGHT; // in lt-s, which is ptot for parabolic orbits
      //      if(m2 < 1e-4)
      //	printf("XXXXXX %Le %Le %Le\n", m2, omegaB, atot);

      if(showerrors)
	atot_err = atot*sqrtl(4.0*(omegaB_err/omegaB)*(omegaB_err/omegaB)+(m2_err/(innermass+m2))*(m2_err/(innermass+m2)))/3.0;

      // semi-major axis of companion in lt-s
      a2 = atot*innermass/(innermass+m2);   // Is p2 for parabolic orbits
      a2_err = 0;
      if(showerrors)
	a2_err = a2*sqrtl(atot_err*atot_err/(atot*atot) + m2_err*m2_err/((innermass+m2)*(innermass+m2)));
      // Current distance CM -> Companion in lt-s
      if(ecc < 1.0) {
	r2 = a2*(1.0-ecc*ecc)/(1.0+ecc*cos(AT)); // Here I assume that the true anomaly of the companion is the same as that of pulsar: i.e. they are in their periastron/apastron at the same time. I think that must be true to keep psr, comp and CM on a single line.
      }else if(ecc == 1.0) {
	r2 = a2/(1.0+cos(AT)); // Here I assume that the true anomaly of the companion is the same as that of pulsar: i.e. they are in their periastron/apastron at the same time. I think that must be true to keep psr, comp and CM on a single line.
      }else {
	r2 = a2*(ecc*ecc-1.0)/(1.0+ecc*cos(AT)); // Here I assume that the true anomaly of the companion is the same as that of pulsar: i.e. they are in their periastron/apastron at the same time. I think that must be true to keep psr, comp and CM on a single line.
      }

      // rpsini is really the projected rp*sin(i), since a1 really is simi-major axis * sin(i)
      // The x-axis is choosen in the plane of the sky, the y-axis is in the plane of the orbit, orthogonal to the x-axis.
      // rpsini should be replaced with rpsini/sin(i) to get actual shape of orbit.
      // To get orbit projected in plane sky, the x value remains unchanged (as it is already in the plane of the sky), but the y-axis gets multiplied again with sin(i) to do the projection.
      *xpsini = rpsini*cos(AT+omega_psr);
      *ypsini = rpsini*sin(AT+omega_psr);
      
      // In lt-s
      *xc = r2*cos(AT+omega_psr+M_PI);
      *yc = r2*sin(AT+omega_psr+M_PI);


      if(doprint == 0) {
	break;
      }

      if(mjd > mjd1) {
	printf("Companion %d\n", companion + 1);
	printf("  Input: Pb=%Lf +- %Lf days, T0=%Lf, Pb_dot=%Le, asini=%Le +- %Le lt-s, e=%Le +- %Le, om=%Lf deg\n", pb, pb_err, t0, pbdot, asini_psr, asini_psr_err, ecc, ecc_err, omega_psr*180.0/M_PI);
	printf("  Epoch of passage ascending node       = %Lf\n", epoch_an);
	printf("  At mjd %Lf:\n", mjd);
	printf("    Mean anomaly                        = %Lf deg\n", M*180.0/M_PI);
	if(ecc != 1.0)
	  printf("    Eccentric anomaly                   = %Lf deg\n", E*180.0/M_PI);
	printf("    True anomaly                        = %Lf deg\n", AT*180.0/M_PI);
	printf("  Results independent of the assumed mass of the pulsar\n");
	//      asini_psr = parfile->parameter[TEMPO3_PARAM_A1+companion];       // lt-s
	if(ecc < 1.0) {
	  value = asini_psr*(1.0-ecc);
	  err = sqrt((1.0-ecc)*(1.0-ecc)*asini_psr_err*asini_psr_err + asini_psr*asini_psr*ecc_err*ecc_err);
	  if(showerrors)
	    printf("    Periastron distance*sin(i) (PSR-CM) = %Lf +- %Lf lt-s = %Lf +- %Lf AU = %Lf +- %Lf km\n", value, err, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, err*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0, err*TEMPO3_SPEEDOFLIGHT/1000.0);
	  else
	    printf("    Periastron distance*sin(i) (PSR-CM) = %Lf lt-s = %Lf AU = %Lf km\n", value, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0);
	  value = asini_psr*(1.0+ecc);
	  err = sqrt((1.0+ecc)*(1.0+ecc)*asini_psr_err*asini_psr_err + asini_psr*asini_psr*ecc_err*ecc_err);
	  if(showerrors)
	    printf("    Apastron distance*sin(i) (PSR-CM)   = %Lf +- %Lf lt-s = %Lf +- %Lf AU = %Lf +- %Lf km\n", value, err, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, err*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0, err*TEMPO3_SPEEDOFLIGHT/1000.0);
	  else
	    printf("    Apastron distance*sin(i) (PSR-CM)   = %Lf lt-s = %Lf AU = %Lf km\n", value, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0);
	}else if(ecc == 1.0) {
	  value = 0.5*asini_psr;
	  err = 0.5*asini_psr_err;
	  if(showerrors)
	    printf("    Periastron distance*sin(i) (PSR-CM) = %Lf +- %Lf lt-s = %Lf +- %Lf AU = %Lf +- %Lf km\n", value, err, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, err*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0, err*TEMPO3_SPEEDOFLIGHT/1000.0);
	  else
	    printf("    Periastron distance*sin(i) (PSR-CM) = %Lf lt-s = %Lf AU = %Lf km\n", value, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0);
	  printf("    Apastron distance*sin(i) (PSR-CM)   = Inf (per definition)\n");
	  if(showerrors == 0) {
	    printf("    Maximum true anomaly                = 180.0 deg (per definition)\n");
	  }else {
	    printf("    Maximum true anomaly                = 180.0 +- 0.0 deg (per definition)\n");
	  }
	}else {
	  value = asini_psr*(ecc-1.0);
	  err = sqrt((1.0-ecc)*(1.0-ecc)*asini_psr_err*asini_psr_err + asini_psr*asini_psr*ecc_err*ecc_err);
	  if(showerrors)
	    printf("    Periastron distance*sin(i) (PSR-CM) = %Lf +- %Lf lt-s = %Lf +- %Lf AU = %Lf +- %Lf km\n", value, err, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, err*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0, err*TEMPO3_SPEEDOFLIGHT/1000.0);
	  else
	    printf("    Periastron distance*sin(i) (PSR-CM) = %Lf lt-s = %Lf AU = %Lf km\n", value, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0);
	  printf("    Apastron distance*sin(i) (PSR-CM)   = Inf (per definition)\n");
	  if(showerrors == 0) {
	    printf("    Maximum true anomaly                = %Lf deg\n", acosl(-1.0/ecc)*180.0/M_PI);
	  }else {
	    printf("    Maximum true anomaly                = %Lf +- %Lf deg\n", acosl(-1.0/ecc)*180.0/M_PI, ecc_err/(ecc*sqrtl(ecc*ecc-1.0))*180.0/M_PI);
	  }
	}
	printf("    Current distance*sin(i) (PSR-CM)    = %Lf lt-s = %Lf AU\n", rpsini, rpsini*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU);
	printf("  Assuming pulsar mass of %Lf Msun and inclination i=%Lf deg\n", innermass, inclination*180.0/M_PI);
	if(showerrors) {
	  if(m2 < 1e-3)
	    printf("    Companion mass                      = %Le +- %Le Msun = %Lf +- %Lf Mearth\n", m2, m2_err, m2*1.989e6/5.972, m2_err*1.989e6/5.972);
	  else
	    printf("    Companion mass                      = %Lf +- %Lf Msun = %Lf +- %Lf Mearth\n", m2, m2_err, m2*1.989e6/5.972, m2_err*1.989e6/5.972);
	}else {
	  printf("    Companion mass                      = %Lf Msun = %Lf Mearth\n", m2, m2*1.989e6/5.972);
	}
	if(ecc < 1.0) {
	  value = a2;
	  printf("    Semi-major axis (Comp-CM)           = %Lf lt-s = %Lf AU = %Lf km\n", value, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0);
	  value = a2*(1.0-ecc);
	  err = sqrtl((1.0-ecc)*(1.0-ecc)*a2_err*a2_err+a2*a2*ecc_err*ecc_err);
	  if(showerrors)
	    printf("    Periastron distance (Comp-CM)       = %Lf +- %Lf lt-s = %Lf +- %Lf AU = %Lf +- %Lf km\n", value, err, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, err*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0, err*TEMPO3_SPEEDOFLIGHT/1000.0);
	  else
	    printf("    Periastron distance (Comp-CM)       = %Lf lt-s = %Lf AU = %Lf km\n", value, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0);
	  value = a2*(1.0+ecc);
	  err = sqrtl((1.0+ecc)*(1.0+ecc)*a2_err*a2_err+a2*a2*ecc_err*ecc_err);
	  if(showerrors)
	    printf("    Apastron distance (Comp-CM)         = %Lf +- %Lf lt-s = %Lf +- %Lf AU = %Lf +- %Lf km\n", value, err, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, err*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0, err*TEMPO3_SPEEDOFLIGHT/1000.0);
	  else
	    printf("    Apastron distance (Comp-CM)         = %Lf lt-s = %Lf AU = %Lf km\n", value, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0);
	}else if(ecc == 1.0) {
	  value = 0.5*a2;
	  err = 0.5*a2_err;
	  if(showerrors)
	    printf("    Periastron distance (Comp-CM)       = %Lf +- %Lf lt-s = %Lf +- %Lf AU = %Lf +- %Lf km\n", value, err, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, err*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0, err*TEMPO3_SPEEDOFLIGHT/1000.0);
	  else
	    printf("    Periastron distance (Comp-CM)       = %Lf lt-s = %Lf AU = %Lf km\n", value, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0);
	  printf("    Apastron distance (Comp-CM)         = Inf (per definition)\n");
	}else {
	  value = a2*(ecc-1.0);
	  err = sqrtl((1.0-ecc)*(1.0-ecc)*a2_err*a2_err+a2*a2*ecc_err*ecc_err);
	  if(showerrors)
	    printf("    Periastron distance (Comp-CM)       = %Lf +- %Lf lt-s = %Lf +- %Lf AU = %Lf +- %Lf km\n", value, err, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, err*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0, err*TEMPO3_SPEEDOFLIGHT/1000.0);
	  else
	    printf("    Periastron distance (Comp-CM)       = %Lf lt-s = %Lf AU = %Lf km\n", value, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0);
	  printf("    Apastron distance (Comp-CM)         = Inf (per definition)\n");
	  value = a2*ecc*sqrtl(ecc*ecc-1.0);
	  err = value*sqrtl(a2_err*a2_err/(a2*a2)+((2.0*ecc*ecc-1.0)/(ecc*ecc-1.0))*(ecc_err*ecc_err/(ecc*ecc)));
	  if(showerrors)
	    printf("    Impact parameter (Comp-CM)          = %Lf +- %Lf lt-s = %Lf +- %Lf AU = %Lf +- %Lf km\n", value, err, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, err*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0, err*TEMPO3_SPEEDOFLIGHT/1000.0);
	  else
	    printf("    Impact parameter (Comp-CM)          = %Lf lt-s = %Lf AU = %Lf km\n", value, value*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, value*TEMPO3_SPEEDOFLIGHT/1000.0);
	}
	printf("    Current distance (Comp-CM)          = %Lf lt-s = %Lf AU\n", r2, r2*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU);
	rtot = rpsini/sin(inclination) + r2;
	//	v = sqrtl(GM*(innermass+m2)*(2.0/r2-1.0/a2)/TEMPO3_SPEEDOFLIGHT);
	//	v = sqrtl(GM*(innermass)*(2.0/r2-1.0/a2)/TEMPO3_SPEEDOFLIGHT);
	if(ecc < 1.0)
	  vrel = sqrtl(GM*(innermass+m2)*(2.0/rtot-1.0/atot)/TEMPO3_SPEEDOFLIGHT);
	else if(ecc == 1.0)	
	  vrel = sqrtl(2.0*GM*(innermass+m2)/(rtot*TEMPO3_SPEEDOFLIGHT));
	else
	  vrel = sqrtl(GM*(innermass+m2)*(2.0/rtot+1.0/atot)/TEMPO3_SPEEDOFLIGHT);
	printf("    Current relative orbital velocity   = %Le m/s\n", vrel);
	v2 = (innermass/(innermass+m2))*vrel;
	printf("    Current orbital velocity            = %Le m/s\n", v2);
	if(ecc < 1.0)
	  specific_energy = -0.5*GM*(innermass+m2)/(atot*TEMPO3_SPEEDOFLIGHT);
	else
	  specific_energy = 0.5*GM*(innermass+m2)/(atot*TEMPO3_SPEEDOFLIGHT);
	printf("    Specific orbital energy             = %Le J/kg\n", specific_energy);
	if(ecc > 1.0) {
	  long double vinf, vinf_err;
	  vinf = sqrtl(GM*(innermass+m2)/(atot*TEMPO3_SPEEDOFLIGHT));
	  vinf_err = vinf*sqrtl(atot_err*atot_err/(atot*atot)+m2_err*m2_err/((innermass+m2)*(innermass+m2)))/2.0;
	  if(showerrors)
	    printf("    Excess relative velocity at infinity  = %Le +- %Le m/s\n", vinf, vinf_err);
	  else
	    printf("    Excess relative velocity at infinity  = %Le m/s\n", vinf);
	  err = vinf*sqrtl(vinf_err*vinf_err/(vinf*vinf)+m2_err*m2_err/((innermass+m2)*(innermass+m2)));
	  if(showerrors)
	    printf("    Excess velocity companion at infinity = %Le +- %Le m/s\n", vinf*innermass/(innermass+m2), err);
	  else
	    printf("    Excess velocity companion at infinity = %Le m/s\n", vinf*innermass/(innermass+m2));
	  value = vinf*m2/(innermass+m2);
	  err = value*sqrtl(vinf_err*vinf_err/(vinf*vinf)+(innermass*m2_err/(m2*m2))*(innermass*m2_err/(m2*m2)));
	  if(showerrors)
	    printf("    Excess velocity pulsar at infinity    = %Le +- %Le m/s\n", value, err);
	  else
	    printf("    Excess velocity pulsar at infinity    = %Le m/s\n", value);
	}
	long double gamma;
	gamma = atan2l(ecc*sin(AT), 1.0+ecc*cos(AT));
	printf("    Current flight path angle           = %Le deg\n", gamma*180.0/M_PI);
	printf("    Current orbital position            = (%Le , %Le) lt-s = (%Le , %Le) AU\n", r2*cos(AT+M_PI), r2*sin(AT+M_PI), r2*cos(AT+M_PI)*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, r2*sin(AT+M_PI)*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU);
	printf("    Current orbital position of pulsar  = (%Le , %Le) lt-s = (%Le , %Le) AU\n", (rpsini/sin(inclination))*cos(AT), (rpsini/sin(inclination))*sin(AT), (rpsini/sin(inclination))*cos(AT)*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, (rpsini/sin(inclination))*sin(AT)*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU);
	printf("    Current orbital velocity            = (%Le , %Le) m/s\n", -v2*sin(AT-gamma+M_PI), v2*cos(AT-gamma+M_PI));
	vpsr = (m2/(innermass+m2))*vrel;
	printf("    Current orbital velocity of pulsar  = (%Le , %Le) m/s\n", -vpsr*sin(AT-gamma), vpsr*cos(AT-gamma));

	printf("  Some internal checks to see if things are consistent:\n");
	if(ecc < 1.0)
	  printf("    Check of Kepler's equation:         %Lf - %Lf rad  = %Le rad (should be zero)\n", E-ecc*sin(E), phase, E-ecc*sin(E)- phase);
	else if(ecc > 1.0)
	  printf("    Check of Kepler's equation:         %Lf - %Lf rad  = %Le rad (should be zero)\n", ecc*sinh(E)-E, phase, ecc*sinh(E) - E - phase);
	if(ecc < 1.0) {
	  printf("    Check consistency of anomalies (x): %Lf - %Lf lt-s = %Le lt-s (should be zero)\n", rpsini*cosl(AT), asini_psr*(cosl(E)-ecc), rpsini*cosl(AT) - asini_psr*(cosl(E)-ecc));
	  printf("    Check consistency of anomalies (y): %Lf - %Lf lt-s = %Le lt-s (should be zero)\n", rpsini*sinl(AT), asini_psr*sqrtl(1.0-ecc*ecc)*sinl(E), rpsini*sinl(AT) - asini_psr*sqrtl(1-ecc*ecc)*sinl(E));
	}else if(ecc > 1.0) {
	  //XXXXXXXXXXXXXXXXXXXXXXXXXXX
	  printf("    Check consistency of anomalies (x): %Lf - %Lf lt-s = %Le lt-s (should be zero)\n", rpsini*cosl(AT), asini_psr*(ecc-coshl(E)), rpsini*cosl(AT) - asini_psr*(ecc-coshl(E)));
	  printf("    Check consistency of anomalies (y): %Lf - %Lf lt-s = %Le lt-s (should be zero)\n", rpsini*sinl(AT), asini_psr*sqrtl(ecc*ecc-1.0)*sinhl(E), rpsini*sinl(AT) - asini_psr*sqrtl(ecc*ecc-1.0)*sinhl(E));
	}
	// MJD in seconds
	// whichcompanion = -1 = all effects combined, otherwise only that of one companion (counting from 0)
	long double torb;
	torb = -ephemeris_get_torb(epoch_an*TEMPO3_SecPerDay, parfile, paramset, 0, companion);
	printf("    Check epoch ascending node in timing model                       %Le s (should be zero)\n", torb);
	torb = -ephemeris_get_torb(mjd*TEMPO3_SecPerDay, parfile, paramset, 0, companion);
	printf("    Check of timing model at given time %Lf - %Lf s    = %Le s (should be zero)\n", torb, *ypsini, torb - *ypsini);
	printf("    Check momentum conservation:        %Lf - %Lf lt-s = %Le lt-s (should be zero)\n", asini_psr*innermass/sin(inclination), a2*m2, asini_psr*innermass/sin(inclination) - a2*m2);
	printf("    Check 2 momentum conservation:      %Lf - %Lf lt-s = %Le lt-s (should be zero)\n", rpsini*innermass/sin(inclination), r2*m2, rpsini*innermass/sin(inclination) - r2*m2);
	long double vmeas, xrel, yrel, xrel_old, yrel_old;
	vmeas = sqrtl((*xc - xc_old)*(*xc - xc_old)+(*yc - yc_old)*(*yc - yc_old))*TEMPO3_SPEEDOFLIGHT/(dt*TEMPO3_SecPerDay);
	printf("    Check current orbital velocity      %Le - %Le m/s  = %Le m/s (should be zero but follows from interpolation)\n", vmeas, v2, vmeas-v2);
	xrel = *xc - *xpsini/sin(inclination);
	yrel = *yc - *ypsini/sin(inclination);
	xrel_old = xc_old - xpsini_old/sin(inclination);
	yrel_old = yc_old - ypsini_old/sin(inclination);
	vmeas = sqrtl((xrel-xrel_old)*(xrel-xrel_old) + (yrel-yrel_old)*(yrel-yrel_old))*TEMPO3_SPEEDOFLIGHT/(dt*TEMPO3_SecPerDay);
	// Note: if relative orbital velocity is ok, it implies that the specific orbital energy will be concerved in the orbit
	printf("    Check relative orbital velocity     %Le - %Le m/s  = %Le m/s (should be zero but follows from interpolation)\n", vmeas, vrel, vmeas-vrel);
      }
      xc_old = *xc;
      yc_old = *yc;
      xpsini_old = *xpsini;
      ypsini_old = *ypsini;
    }
  }else {
    printerror(0, "derive_orbital_parameters(): Requested companion is not defined.");
    exit(0);
  }
}

//double *mrq_alpha, *mrq_alpha_inv -> just temporary space
//long double *mrq_deriv -> just temporary space
void finderrors_tempo3(long double *dplus, long double *dmin, tempo3_parameters_def *parfile, toa_data_def toas, residual_def residuals, int *fixed, int *paramset, int *paramlinked, long double *xstart, long double *dx, long double *xfit, double *mrq_alpha, double *mrq_alpha_inv, long max_nr_fitted_params, long double *mrq_deriv, int dumpcovar, verbose_definition verbose_state)
{
  int nrtoasincluded, i, j, k;
  int *fitparam_list;
  long double chi2w, rms;
  double sigma;
  tempo3_parameters_def parfile_tmp;
  allocate_ephemeris_par_only(&parfile_tmp, 0, verbose_state);
  
  copy_parfile(&parfile_tmp, parfile);
  calcChi2(&residuals, toas, &rms, parfile->parameter[TEMPO3_PARAM_F0], &chi2w);


  fitparam_list = malloc(TEMPO3_PARAM_TOTAL*sizeof(int));
  if(fitparam_list == NULL) {
    printerror(verbose_state.debug, "finderrors_tempo3: Memory allocation error");
    exit(0);
  }

  long nr_actual_fit_parameters;
  nr_actual_fit_parameters = 0;
  for(i = 0; i < TEMPO3_PARAM_TOTAL; i++) {
    if(fixed[i] == 0) {
      fitparam_list[nr_actual_fit_parameters] = i;  // Keep a list of the actually used parameters
      nr_actual_fit_parameters++;
    }
  }
  if(max_nr_fitted_params > 0) {
    if(nr_actual_fit_parameters > max_nr_fitted_params) {
      printerror(verbose_state.debug, "ERROR tempo3: Maximum number of parameters that can be fit simultaneously is exceeded (%ld > %ld).\n", nr_actual_fit_parameters, max_nr_fitted_params);
      exit(0);
    }
  }

  /* Initialise mrq_alpha matrix with zeros  */
  for(j = 0; j < nr_actual_fit_parameters*nr_actual_fit_parameters; j++) {
    mrq_alpha[j] = 0;
  }
  
  memcpy(xstart, xfit, TEMPO3_PARAM_TOTAL*sizeof(long double));
  ephemeris_setLinkedValues_par_version(&parfile_tmp, paramlinked);
  
  double remember_factor;
  // In not-weighted mode, sigma is always the same, and shouldn't be re-calculated every time in loop below
  if(parfile->mode == 0) {
    // Count nr of toas
    nrtoasincluded = 0;
    for(j = 0; j < toas.nrtoas+toas.nrtoasdeleted; j++) {
      if(toas.deleted[j] == 0) 
	if(toas.mjd[j] >= residuals.mjdmin && toas.mjd[j] <= residuals.mjdmax) 
	  nrtoasincluded += 1;
    }
    // all errors were taken to be 1
    sigma = 1.0; 
    // Scale all errorbars such that rms would be 1 (rms is really variance, so have to devide by the number of points)
    sigma *= sqrt(rms/nrtoasincluded); 
    sigma = sigma*sigma;
  }else {
    remember_factor = 1e-6*parfile->parameter[TEMPO3_PARAM_F0]*sqrt(chi2w/(long double)(toas.nrtoas-nr_actual_fit_parameters));
    printwarning(0, "WARNING: Be aware that all errors are scaled to make the reduced chi2 = 1 in the error calculation.");
  }
  printwarning(0, "WARNING: Be aware that the errorbars are calculated under the assumption that the timing solution is linear.");
  for(j = 0; j < toas.nrtoas+toas.nrtoasdeleted; j++) {
    /* Calculate derivatives of function y at position x */
    if(toas.deleted[j] == 0) {
      if(toas.mjd[j] >= residuals.mjdmin && toas.mjd[j] <= residuals.mjdmax) {
	for(i = 0; i < nr_actual_fit_parameters; i++) {
	  k = fitparam_list[i];
	  long double stepsize;
	  if(dx != NULL) {
	    stepsize = 0.001*dx[k];
	  }else {
	    stepsize = -1;
	  }
	  // PATRICK: CAN USE freqSite RATHER THAN freqSSB TO REVERT TO OLD BEHAVIOUR
	  mrq_deriv[i] =  calc_derivative(k, xstart, &parfile_tmp, stepsize, toas.mjd[j], toas.freqSSB[j], toas.ssb1[j], toas.ssb2[j], toas.ssb3[j], toas.site[j], toas.flags[j], paramlinked, paramset);
	}
	for(i = 0; i < nr_actual_fit_parameters; i++) {
	  for(k = 0; k < nr_actual_fit_parameters; k++) {
	    if(parfile->mode != 0) {
	      //		      if(j == 0)
	      //			fprintf(stderr, "XXXXX %f (%f) %Lf\n", toas.err[j]*1e-6, toas.mjdx[j], chi2w/(toas.nrtoas-nrparamsused));
	      //		      sigma = toas.err[j]*1e-6*parfile->parameter[TEMPO3_PARAM_F0]; // Convert errorbar in turns
	      //sigma = toas.err[j]*1e-6*sqrtl(chi2w/(long double)(toas.nrtoas-nrparamsused)); // Convert errorbar in turns
	      sigma = toas.err[j]*remember_factor; // Convert errorbar in turns
	      sigma = sigma*sigma;
	    }
	    if(sigma == 0)
	      sigma = 1e-10;
	    // When not squaring here, can avoid doing square root later on? Probably not because elements are summed
	    // So mrq_deriv is d(phase)/d(parameter) = basis function value at point of TOA = X.
	    // So mrq_alpha = alpha nrparams X nrparams matrix of NR.
	    mrq_alpha[nr_actual_fit_parameters*i+k] += mrq_deriv[i]*mrq_deriv[k]/(sigma);
	  }
	}
      }
    }
  }
  
  // Check if there are zero's on diagonal, which means that the matrix cannot be inverted.
  // If there are zero's on the diagonal, it probably means that the fit function is independent of the parameter. However, since alpha is a sum over the number of data points, couldn't an element of alpha being zero by chance without it being a problem? If you want to check if the matrix is invertable, check if the determinant is zero.
  for(i = 0; i < nr_actual_fit_parameters; i++) {
    if(mrq_alpha[nr_actual_fit_parameters*i+i] == 0) {
      printerror(0, "ERROR, CANNOT DETERMINE ERRORBARS!\n");
      // Cannot remember why I did this, just to get rediculus error-bars if co-variant parameters are defined, but avoiding a crash?
      for(k = 0; k < nr_actual_fit_parameters; k++) {
	for(j = 0; j < nr_actual_fit_parameters; j++) {
	  // mrq_alpha[1+k][1+j]
	  if(k == j) {
	    mrq_alpha[nr_actual_fit_parameters*k+j] = 1e-99;
	  }else {
	    mrq_alpha[nr_actual_fit_parameters*k+j] = 0;
	  }
	}
      }
      break;
    }
  }
  
  // Define the solution vectors to be a unit matrix
  /*
    for(i = 0; i < nr_actual_fit_parameters; i++) {
    for(j = 0; j < nr_actual_fit_parameters; j++) {
    //	    mrq_covar[1+i][1+j]
    if(i == j)
    mrq_covar[i*nr_actual_fit_parameters+j] = 1;
    else
    mrq_covar[i*nr_actual_fit_parameters+j] = 0;
    }
    }
  */
  
  /*	fprintf(stderr, "\nXXXXXXX %d\n", nrparamsused); */
  /*	for(i = 0; i < nrparamsused; i++) {
	for(j = 0; j < nrparamsused; j++) {
	fprintf(stderr, "mrq_alpha[%d][%d] = %Le\n", i+1, j+1, mrq_alpha[1+i][1+j]);
	}
	}
	for(i = 0; i < nrparamsused; i++) {
	for(j = 0; j < nrparamsused; j++) {
	fprintf(stderr, "mrq_covar[%d][%d] = %Le\n", i+1, j+1, mrq_covar[1+i][1+j]);
	  }
	}*/
  int errors_found;
  errors_found = 1;
  if(linalg_solve_matrix_eq(mrq_alpha, nr_actual_fit_parameters, 0, NULL, 0, NULL, 0, mrq_alpha_inv, 1, 1, verbose_state) != 0) {  // Failed
    printwarning(0, "WARNING tempo3: Error bars could not be determined. Two or more parameters might be 100%% covariant. Error-bars set to -1.");
    errors_found = 0;
  }
  /*	fprintf(stderr, "\nXXXXXXX\n"); */
  
  /*	for(i = 0; i < nrparamsused; i++) {
	for(j = 0; j < nrparamsused2; j++) {
	printf("mrq_alpha[%d][%d] = %Le\n", i+1, j+1, mrq_alpha[1+i][1+j]);
	}
	}*/
  for(i = 0; i < nr_actual_fit_parameters; i++) {
    k = fitparam_list[i];
    if(errors_found)
      dmin[k] = dplus[k] = sqrt(mrq_alpha_inv[nr_actual_fit_parameters*i+i]);
    else
      dmin[k] = dplus[k] = -1;
  }

  if(dumpcovar) {
    if(dumpcovar == 2) {
      for(i = 0; i < nr_actual_fit_parameters; i++) {
	double norm;
	norm = sqrt(mrq_alpha_inv[nr_actual_fit_parameters*i+i]);
	for(j = 0; j < nr_actual_fit_parameters; j++) {
	  mrq_alpha_inv[nr_actual_fit_parameters*i+j] /= norm;    // Divide row by normalisation
	  mrq_alpha_inv[nr_actual_fit_parameters*j+i] /= norm;    // Divide column by normalisation
	}
      }
    }
    printf("\nCovariance matrix:\n");
    printf("     ");
    for(i = 0; i < nr_actual_fit_parameters; i++) {
      k = fitparam_list[i];
      printf("%13s ", tempo3_parfile_descr.identifiers[k][0]);
    }
    printf("\n");
    for(i = 0; i < nr_actual_fit_parameters; i++) {
      k = fitparam_list[i];
      printf("%6s:", tempo3_parfile_descr.identifiers[k][0]);
      for(j = 0; j < nr_actual_fit_parameters; j++) {
	if(mrq_alpha_inv[nr_actual_fit_parameters*i+j] >= 0)
	  printf(" ");
	printf(" %e", mrq_alpha_inv[nr_actual_fit_parameters*i+j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  free_ephemeris_par_only(&parfile_tmp);
  free(fitparam_list);
}

