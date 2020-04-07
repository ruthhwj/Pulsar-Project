//START REGION RELEASE
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include "psrsalsa.h"

#define AmoebaAlgorithm    0   // 1=NR

//START REGION DEVELOP

//START REGION RELEASE
/* Read a text file with the vonMises functions. 

  Return 1: OK
         0: Cannot open file
*/
int readVonMisesModel(char *filename, vonMises_collection_definition *components, verbose_definition verbose)
{
  int i;
  FILE *fin;
  if(verbose.verbose) printf("Opening %s for reading\n", filename);
  fin = fopen(filename, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readVonMisesModel: Cannot open %s", filename);
    return 0;
  }
  for(components->nrcomponents = 0; components->nrcomponents < maxNrVonMisesComponents; (components->nrcomponents)++) {
    i = fscanf(fin, "%lf %lf %lf", &(components->centre[components->nrcomponents]), &(components->concentration[components->nrcomponents]), &(components->height[components->nrcomponents]));
    if(i != 3) {
      break;
    }else if(verbose.verbose) {
      printf("  component %d: phase=%e concentration=%e amplitude=%e\n", components->nrcomponents, components->centre[components->nrcomponents], components->concentration[components->nrcomponents], components->height[components->nrcomponents]);
    }
  }
  if(verbose.verbose) printf("Closing %s\n", filename);
  fclose(fin);
  if(components->nrcomponents == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readVonMisesModel: No components found in %s", filename);
    return 0;
  }
  return 1;
}
//START REGION DEVELOP

// When concentration is found to be negative, it can be transformed in a different function with a positive concentration.
void vonMises_simplify_parameters(vonMises_collection_definition *components)
{
  int i;
  for(i = 0; i < components->nrcomponents; i++) {
    // If concentration is negative, it can be written as a different von Mises function with a positive concentration
    if(components->concentration[i] < 0) {
      components->centre[i] += 0.5;
      if(components->centre[i] >= 1)
	components->centre[i] -= 1.0;
      components->concentration[i] *= -1.0;
      components->height[i] *= exp(2.0*components->concentration[i]);
    }
  }
}

/* Write a text file with the vonMises functions.

   If filename = NULL, the stdout is used as input instead.

  Return 1: OK
         0: Cannot open/write file
*/
int writeVonMisesModel(char *filename, vonMises_collection_definition *components, verbose_definition verbose)
{
  int n;
  FILE *fout;
  if(filename != NULL) {
    if(verbose.verbose) printf("Opening %s for writing\n", filename);
    fout = fopen(filename, "w");
    if(fout == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeVonMisesModel: Cannot open %s", filename);
      return 0;
    }
  }else {
    fout = stdout;
  }
  for(n = 0; n < components->nrcomponents; n++) {
    fprintf(fout, "%e %e %e\n", components->centre[n], components->concentration[n], components->height[n]);
    if(verbose.verbose) {
      if(n == 0) {
	printf("y = A*exp((cos(2.0*M_PI*(x-x0))-1.0)*concentration)\n");
	printf("  component: x0 concentration H\n");
      }
      printf("  component %d: %e %e %e\n", n, components->centre[n], components->concentration[n], components->height[n]);
    }
  }
  if(filename != NULL) {
    if(verbose.verbose) printf("Closing %s\n", filename);
    fclose(fout);
  }
  return 1;
}

//START REGION RELEASE
/* Calculate a vonMises function at a certain phase by applying a
   phase shift. */
double calcVonMisesFunction2(double centre, double concentration, double height, double phase, double shift)
{
  double y;
  y = exp((cos(2.0*M_PI*(phase-centre-shift))-1.0)*concentration) * height;
  return y;
}

//START REGION DEVELOP
/* Calculate the derivative of a vonMises function (with respect to
   phase) at a certain phase by applying a phase shift. */
double calcVonMisesFunction2_deriv(double centre, double concentration, double height, double phase, double shift)
{
  double y;
  y = calcVonMisesFunction2(centre, concentration, height, phase, shift);
  y *= -2.0*M_PI*concentration*sin(2.0*M_PI*(phase-centre-shift));
  return y;
}

//START REGION RELEASE


/* Integrate a vonMises function over its entire domain. 
Here the 2pi in the formula corresponds to a full rotation in radians, so the units are intensity*radians. 
Replace the 2pi with 360 or the number of bins across the profile to get it in different units.
*/
double integrateVonMisesFunction2(double concentration, double height)
{
  double area;
  //  fprintf(stderr, "CONCENTRATION: %lf\n", concentration);
  //  area = 2*M_PI*height*gsl_sf_bessel_I0(concentration)*exp(-concentration);
  // I0_scaled = I0*exp(-fabs(concentration)). For large values of the concentration this avoids overflows from happening.
  area = 2*M_PI*height*gsl_sf_bessel_I0_scaled(concentration);
  return area;
}


/* Returns the width of a single vonMises function given its
   concentration in radians. The width is measured at a level ampfrac,
   which is a number between 0 and 1. A fwhm corresponds to ampfrac =
   0.5 and a W10 to an ampfrac = 0.1. For example the fwhm is defined
   at half the peak amplitude. So the actual minimum of the function
   (which is only asymptotically approaching zero for narrow
   distributions) is not considered. Returns sqrt(-1) if the fwhm doesn't
   exists. This happens if the concentration is very low, which means
   that the function never goes to zero. */
double widthVonMisesFunction2(double concentration, double ampfrac)
{
  // exp((cos(2.0*M_PI*(phase-centre-shift))-1.0)*concentration) = 0.5
  // (cos(2.0*M_PI*(phase-centre-shift))-1.0)*concentration = ln(0.5)
  // cos(2.0*M_PI*(phase-centre-shift)) = ln(0.5)/concentration + 1.0
  double value = log(ampfrac)/concentration + 1.0;
  if(value < -1.0 || value >= 1.0) {
    return sqrt(-1);
  }
  return 2.0*acos(value);
}


/* Calculate the profile out of the model by at a certain phase by
   applying a phase shift. */
double calcVonMisesFunction(vonMises_collection_definition *components, double phase, double shift)
{
  int n;
  double y;
  y = 0;
  if(components->nrcomponents > 0) {
    for(n = 0; n < components->nrcomponents; n++) {
      y += calcVonMisesFunction2(components->centre[n], components->concentration[n], components->height[n], phase, shift);
    }
  }
  return y;
}
//START REGION DEVELOP

/* Calculate the derivative (with respect to phase) of the profile (the model) at a certain
   phase by applying a phase shift. */
double calcVonMisesFunction_deriv(vonMises_collection_definition *components, double phase, double shift)
{
  int n;
  double y;
  y = 0;
  if(components->nrcomponents > 0) {
    for(n = 0; n < components->nrcomponents; n++) {
      y += calcVonMisesFunction2_deriv(components->centre[n], components->concentration[n], components->height[n], phase, shift);
    }
  }
  return y;
}

/* Calculate the integral of the model over its entire domain.  Here
the 2pi in the formula corresponds to a full rotation in radians, so
the units are intensity*radians.  Replace the 2pi with 360 or the
number of bins across the profile to get it in different units.*/
double integrateVonMisesFunction(vonMises_collection_definition *components)
{
  int n;
  double area;
  area = 0;
  if(components->nrcomponents > 0) {
    for(n = 0; n < components->nrcomponents; n++) {
      area += integrateVonMisesFunction2(components->concentration[n], components->height[n]);
    }
  }
  return area;
}

//START REGION RELEASE
/* Calculate the profile out of the model by applying a phase
   shift. If the normalize flag is set the profile will be
   normalized. */
void calcVonMisesProfile(vonMises_collection_definition *components, int nrbins, float *profile, double shift, int normalize)
{
  int i;
  double x, Imax = -1;

  for(i = 0; i < nrbins; i++) {
    profile[i] = 0;
  }
  for(i = 0; i < nrbins; i++) {
    x = i/(double)nrbins;
    profile[i] = calcVonMisesFunction(components, x, shift);
  }
  if(normalize) {
    for(i = 0; i < nrbins; i++) {
      if(profile[i] > Imax)
	Imax = profile[i];
    }
    for(i = 0; i < nrbins; i++) {
      profile[i] /= Imax;
    }
  }
}
//START REGION DEVELOP

/* For a given profile with nrbins bins, calculate the rms of the
   profile - model after applying a phase shift to the model. If the
   normalize flag is set the profile will be normalized. */
void calcVonMisesProfile_resid_rms(vonMises_collection_definition *components, int nrbins, float *profile, double shift, double *rms)
{
  int i;
  double x, y;

  *rms = 0;

  for(i = 0; i < nrbins; i++) {
    x = i/(double)nrbins;
    y = profile[i] - calcVonMisesFunction(components, x, shift);
    *rms += y*y;
  }
  *rms /= (double)nrbins;
  *rms = sqrt(*rms);
}

//START REGION DEVELOP
//START REGION RELEASE


/* Find the best shift in phase to match the profile (the return
   value). verbose-1 is number of spaces before output. */
float correlateVonMisesFunction(vonMises_collection_definition *components, int nrbins, float *profile, verbose_definition verbose)
{
  float correl_max, *profile2;
  int i;
  /*, I;
    int i, j, istart, iend, di; */
  int ishift;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Correlating template with profile\n");
    verbose.indent += 2;
  }
  profile2 = (float *)malloc(nrbins*sizeof(float));
  if(profile2 == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR correlateVonMisesFunction: Memory allocation error.");
    return 0;
  }

  calcVonMisesProfile(components, nrbins, profile2, 0, 0);
  find_peak_correlation(profile, profile2, nrbins, 0, 1, 1, &ishift, &correl_max, verbose);

  /*
  ymax = 0;
  ymin = 0;
  ishift = 0;
  // If many bins, do a rough search first 
  if(nrbins > 250) {
    di = nrbins/40;
    for(i = 0; i < nrbins; i++) {
      calcVonMisesProfile(components, nrbins, profile2, i/(float)nrbins, 0);
      I = 0;
      for(j = 0; j < nrbins; j++) {
	I += (profile[j]*profile2[j]);
      }
      if(I > ymax || i == 0) {
	ymax = I;
	ishift = i;
      }
      if(fabs(I) < ymin || i == 0) {
	ymin = fabs(I);
      }
      i += di;
    }
    istart = ishift - di;
    iend = ishift + di;
    if(verbose.verbose) printf("  Found a shift of approx %d bins (%f phase).\n", ishift, ishift/(float)nrbins);
  }else {
    istart = 0;
    iend = nrbins-1;
  }

  if(verbose.verbose) printf("  Search between bins %d - %d.\n", istart, iend);
  // Do the correlation to shift the standard with respect to the profile 
  for(i = 0; i <= iend; i++) {
    calcVonMisesProfile(components, nrbins, profile2, i/(float)nrbins, 0);
    I = 0;
    for(j = 0; j < nrbins; j++) {
      I += (profile[j]*profile2[j]);
    }
    if(I > ymax || i == 0) {
      ymax = I;
      ishift = i;
    }
    if(fabs(I) < ymin || i == 0) {
      ymin = fabs(I);
    }
  }
  */
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Found a shift of %d bins (%f phase) and max/min correlation is %f.\n", ishift, ishift/(float)nrbins, correl_max);
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("so you could do a shift by %d bins (%f phase) to rotate profile to align with model.\n", -ishift, -ishift/(float)nrbins);
  }
  free(profile2);
  return ishift/(float)nrbins;
}

//START REGION DEVELOP

/* 
   Return:  1 Found rising slope
           -1 Found declining slope
	    0 Nothing found after istart
*/
int getboundry_single(float *profile, int nrbins, int istart, float y, int *bin)
{
  int i, slope, found;
  /* slope = 1 -> search for rising slope */
  if(profile[istart] < y)
    slope = 1;
  else 
    slope = -1;
  found = 0;
  for(i = istart; i < nrbins; i++) {
    if(profile[i] > y && slope == 1) {
      *bin = i-1;
      if(*bin < istart)
	*bin = istart;
      found = 1;
      break;
    }
    if(profile[i] < y && slope == -1) {
      *bin = i+1;
      if(*bin >= nrbins)
	*bin = nrbins;
      found = 1;
      break;
    }
  }
  return found*slope;
}


void find_boundaries(float *profile, int nrbins, float y, pulselongitude_regions_definition *regions)
{
  int ret, lastbin, bin, lastbin_set, firstpulse;
  lastbin = -1;
  firstpulse = 1;
  clearPulselongitudeRegion(regions);
  lastbin_set = 0;
  do {
    ret = getboundry_single(profile, nrbins, lastbin+1, y, &bin);
    /*    printf("ret=%d i=%d\n", ret, bin); */
    if(firstpulse == 1 && ret == -1) {
      regions->left_bin[regions->nrRegions] = 0;
      regions->right_bin[regions->nrRegions] = bin;
      regions->bins_defined[regions->nrRegions] = 1;
      (regions->nrRegions)++;
      lastbin = bin;
      firstpulse = 0;
    }else if(ret == 1) {
      lastbin = bin;
      lastbin_set = 1;
      firstpulse = 0;
    }else if(ret == -1) {
      regions->left_bin[regions->nrRegions] = lastbin;
      regions->right_bin[regions->nrRegions] = bin;
      regions->bins_defined[regions->nrRegions] = 1;
      (regions->nrRegions)++;
      lastbin = bin;
      lastbin_set = 0;
      firstpulse = 0;
    }
    if(regions->nrRegions == MAX_pulselongitude_regions)
      ret = 0;
  }while(ret != 0);
  if(lastbin_set && (regions->nrRegions < MAX_pulselongitude_regions)) {
    regions->left_bin[regions->nrRegions] = lastbin;
    regions->right_bin[regions->nrRegions] = nrbins-1;
    regions->bins_defined[regions->nrRegions] = 1;
    (regions->nrRegions)++;
  }
}


int internal_fitvonmises_nrcomponents;
int internal_fitvonmises_nrbins;
int internal_fitvonmises_plotallsteps;
int internal_fitvonmises_fitbaseline;
int internal_fitvonmises_fitdummy;

float *internal_fitvonmises_profile;
float *internal_fitvonmises_fitprofile;

int internal_fitvonmises_relative_phases = 0; // If set, the phases of the components after component 1 are relative to that of the first component.
int internal_fitvonmises_relative_amps = 0;   // If set, the amplitudes of the components after component 1 are relative to that of the first component.
int internal_fitvonmises_debug = 0;           // If set, debug messages are enabled in internal_fitvonmises_funk()
int internal_fitvonmises_showed_overflow_warning = 0; // Prevent many warnings to be produced.
int internal_fitvonmises_avoid_neg_components = 0; // If set, negative amplitudes get a big chi2 penalty during the fitting
pulselongitude_regions_definition *internal_fitvonmises_onpulse = NULL;  // The onpulse region, if fitting is to be limited to a certain onpulse region.
/*
params[0] = baseline (if internal_fitvonmises_fitbaseline is set)
params[1] = center_1
params[2] = concentration_1
params[3] = height_1
....
params[..] = dummy valiable (if internal_fitvonmises_fitdummy is set)
 */
float internal_fitvonmises_funk(float *params)
{
  static vonMises_collection_definition *components = NULL;  // Note that it saves some non-dynamic memory allocation. Data is never freed though....
  int n, i, got_neg_component, got_component_outside_onpulse_region;
  double chi2, chi, y;
  verbose_definition verbose;
  cleanVerboseState(&verbose);
  verbose.debug = internal_fitvonmises_debug;

  if(components == NULL) {
    components = malloc(sizeof(vonMises_collection_definition));
    if(components == NULL) {
      fprintf(stderr, "ERROR internal_fitvonmises_funk: Cannot allocate memory\n");
      exit(0);
    }
  }

  got_neg_component = 0;
  got_component_outside_onpulse_region = 0;
  components->nrcomponents = internal_fitvonmises_nrcomponents;
  for(n = 0; n < internal_fitvonmises_nrcomponents; n++) {

    double centre = params[3*n+internal_fitvonmises_fitbaseline];
    //    printf("comp %d: orig x=%f (relphases = %d)\n", n+1, centre, internal_fitvonmises_relative_phases);
    if(internal_fitvonmises_relative_phases && n > 0) {
      centre += components->centre[0];
    }
    components->centre[n] = centre;

    components->concentration[n] = params[3*n+1+internal_fitvonmises_fitbaseline];

    components->height[n] = params[3*n+2+internal_fitvonmises_fitbaseline];
    if(components->height[n] <= 0.0) {
      got_neg_component = 1;
    }
    if(internal_fitvonmises_onpulse != NULL) {
      if(internal_fitvonmises_onpulse->nrRegions > 0) {
	centre -= floor(centre);
	centre *= internal_fitvonmises_nrbins;
	if(centre >= internal_fitvonmises_nrbins) {
	  centre -= internal_fitvonmises_nrbins;
	}
	if(checkRegions(centre, internal_fitvonmises_onpulse, 0, verbose) == 0) {
	  got_component_outside_onpulse_region = 1;
	}
      }
    }

    if(internal_fitvonmises_relative_amps && n > 0) {
      components->height[n] *= params[2+internal_fitvonmises_fitbaseline];
    }

    if(internal_fitvonmises_plotallsteps) {
      ppgsci(4);
      y = calcVonMisesFunction2(components->centre[n], components->concentration[n], components->height[n], 0, 0);
      if(internal_fitvonmises_fitbaseline)
	y += params[0];
      ppgmove(0, y);
      for(i = 1; i < internal_fitvonmises_nrbins; i++) {
	y = calcVonMisesFunction2(components->centre[n], components->concentration[n], components->height[n], i/(double)internal_fitvonmises_nrbins, 0);
	if(internal_fitvonmises_fitbaseline)
	  y += params[0];
	ppgdraw(i*360.0/(double)internal_fitvonmises_nrbins, y);
      }
      //      printf("comp %d: x=%f c=%f h=%f\n", n+1, components->centre[n], components->concentration[n], components->height[n]);
    }
    //    printf("comp %d: x=%f c=%f h=%f\n", n+1, components->centre[n], components->concentration[n], components->height[n]);
  }

  /*
  printf("Concentrations are: ");
  for(n = 0; n < internal_fitvonmises_nrcomponents; n++) {
    printf(" %f", components->concentration[n]);
  }
  printf("\n");
  */

  calcVonMisesProfile(components, internal_fitvonmises_nrbins, internal_fitvonmises_fitprofile, 0, 0);

  chi2 =0;
  for(n = 0; n < internal_fitvonmises_nrbins; n++) {
    int accept_bin = 1;
    if(internal_fitvonmises_onpulse != NULL) {
      if(internal_fitvonmises_onpulse->nrRegions > 0) {
	if(checkRegions(n, internal_fitvonmises_onpulse, 0, verbose) == 0) {
	  accept_bin = 0;
	}
      }
    }
    if(accept_bin) {
      chi = (internal_fitvonmises_fitprofile[n] - internal_fitvonmises_profile[n]);
      if(internal_fitvonmises_fitbaseline)
	chi += params[0];
      chi2 += chi*chi;
    }
  }

  // Dummy variable adds to chi2 (as x^2), but doesn't affect shape of mode.
  if(internal_fitvonmises_fitdummy) {
    n = 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline;
    chi2 += params[n]*params[n];
    //    printf("XXXXXX dummy=%e\n", params[n]);
    if(isnan(params[n]) || isinf(params[n])) {
      fflush(stdout);
      printerror(verbose.debug, "internal_fitvonmises_funk: chi2=%f produced by dummy variable\n", chi2);
    }
  }

  if(isnan(chi2) || isinf(chi2)) {
    // Something wrong. Probably concentration did became too negative.
    int tooneg;
    tooneg = 0;
    //    printf("Concentrations are: ");
    for(n = 0; n < internal_fitvonmises_nrcomponents; n++) {
      // -40 will probably overflow a floating point calculation of the von Mises function.
      //      printf(" %f", components->concentration[n]);
      if(components->concentration[n] < -15 && isinf(components->concentration[n]) == 0) {
	tooneg = 1;
      }
    }
    //    printf("\n");
    if(tooneg) {
      if(internal_fitvonmises_showed_overflow_warning == 0) {
	fflush(stdout);
	printwarning(verbose.debug, "internal_fitvonmises_funk: At least one trial resulted in a too negative concentration causing a floating-point overflow. These solution will be ignored. This warning will only be shown once.");
	internal_fitvonmises_showed_overflow_warning = 1;
      }
    }

    /*
    if(tooneg == 0) {
      printf("Concentrations are: ");
      for(n = 0; n < internal_fitvonmises_nrcomponents; n++) {
	printf(" %f", components->concentration[n]);
      }
      printf("\n");
    }
    */

    // Ignore the showing the error if it is because of a too negative concentration. Instead, return a very high chi2.
    if(tooneg == 0 || verbose.verbose || verbose.debug) {
      fflush(stdout);
      printerror(verbose.debug, "internal_fitvonmises_funk: chi2=%f", chi2);
      if(verbose.debug) {
	writeVonMisesModel(NULL, components, verbose);
	if(internal_fitvonmises_fitbaseline) {
	  printf("Baseline=%lf\n", params[0]);
	}
	calcVonMisesProfile(components, internal_fitvonmises_nrbins, internal_fitvonmises_fitprofile, 0, 0);
	printf("Fitted function:");
	for(n = 0; n < internal_fitvonmises_nrbins; n++) {
	  printf(" %f", internal_fitvonmises_fitprofile[n]);
	}
	printf("\n");
	printf("Data:");
	for(n = 0; n < internal_fitvonmises_nrbins; n++) {
	  printf(" %f", internal_fitvonmises_profile[n]);
	}
	printf("\n");
	/*
	  float test;
	  printf("Test:\n");
	  test = calcVonMisesFunction(components, 3.0/(float)internal_fitvonmises_nrbins, 0);
	  printf("%f\n", test);
	  test = calcVonMisesFunction2(components->centre[0], components->concentration[0], components->height[0], 0, 0);
	  printf("%f\n", test);
	  test = exp((cos(2.0*M_PI*(-components->centre[0]))-1.0)*components->concentration[0]) * components->height[0];
	  printf("%f\n", test);
	  test = ((cos(2.0*M_PI*(-components->centre[0]))-1.0)*components->concentration[0]);
	  printf("%f\n", test);
	*/
      }
    }
    chi2 = 1e100;
    //        exit(0);
  }

  // Try to avoid negative concentrations. They can be written as positive concentrations, but can more easily explodes the exp() when becoming too negative.
  /*
  for(n = 0; n < internal_fitvonmises_nrcomponents; n++) {
    if(components->concentration[n] < 0) {
      return sqrt(chi2)*1e10;
    }
  }
  */

  if((got_neg_component && internal_fitvonmises_avoid_neg_components) || got_component_outside_onpulse_region) {
    //    printf("Penalty!!!\n");
    chi2 *= 1e10;
  }

  return sqrt(chi2);
}


/* Could improve things if educatedguess is not set by trying both
   signs of the maximum/minimum. Sometimes there is a narrow negative
   spike and a positive component which is then not fitted.


superverbose = 1: extra print statements
                  2: plot best intermediate steps
                  3: plot all intermediate steps
  return 0=error, 1=ok, 99=no new component found 
*/
int fitvonmises_addcomponent(vonMises_collection_definition *components, int refine, int educatedguess, float sigmalimit, int fitvonMises_Nr_X_trials, int fitvonMises_Nr_C_trials, float *baseline, pulselongitude_regions_definition *onpulse, verbose_definition verbose, int superverbose)
{
  float *xstart, *dx, *xfit, best_baseline;
  int *fixed;
  float ymax, ymin, av, chi2, chi2_best, rms, pwr, y, s2n;
  int ix, ic, i, nfunk, n, ret;
  int bin, pulsewidth;
  float snrbest, E_best;
  vonMises_collection_definition *best;

  xstart = malloc((3*maxNrVonMisesComponents+1)*sizeof(float));
  dx     = malloc((3*maxNrVonMisesComponents+1)*sizeof(float));
  xfit   = malloc((3*maxNrVonMisesComponents+1)*sizeof(float));
  fixed  = malloc((3*maxNrVonMisesComponents+1)*sizeof(int));
  best   = malloc(sizeof(vonMises_collection_definition));
  if(xstart == NULL || dx == NULL || xfit == NULL || fixed == NULL || best == NULL) {
    printerror(verbose.debug, "ERROR fitvonmises_addcomponent: Cannot allocate memory");
    return 0;
  }

  /* Make sure refine is 0 or 1 */
  if(refine) {
    refine = 1;
  }

  if(refine == 0) {
    /* Calculate residual profile */
    calcVonMisesProfile(components, internal_fitvonmises_nrbins, internal_fitvonmises_fitprofile, 0, 0);
    for(i = 0; i < internal_fitvonmises_nrbins; i++) {
      internal_fitvonmises_fitprofile[i] = internal_fitvonmises_profile[i]-internal_fitvonmises_fitprofile[i];
    }
    /* Determine maximum amplitude */
    ymax = internal_fitvonmises_fitprofile[0];
    ymin = internal_fitvonmises_fitprofile[0];
    av = 0;
    for(i = 0; i < internal_fitvonmises_nrbins; i++) {
      av += internal_fitvonmises_profile[i];
      if(internal_fitvonmises_fitprofile[i] > ymax)
	ymax = internal_fitvonmises_fitprofile[i];
      if(internal_fitvonmises_fitprofile[i] < ymin)
	ymin = internal_fitvonmises_fitprofile[i];
      /* Make all values positive so that boxcarFindpeak can find maximum if power is negative */
      /*      internal_fitvonmises_fitprofile[i] = fabs(internal_fitvonmises_fitprofile[i]); */
    }
    av /= (float)internal_fitvonmises_nrbins;
    if(fabs(ymin) > fabs(ymax))
      ymax = ymin;

    if(superverbose)
      printf("Residual: amplitude=%f, av=%f\n", ymax, av);
  }

  /*
  memset(xstart, 0, (3*maxNrVonMisesComponents+1)*sizeof(float));
  memset(dx, 0, (3*maxNrVonMisesComponents+1)*sizeof(float));
  memset(fixed, 0, (3*maxNrVonMisesComponents+1)*sizeof(int));
  */

  /* Set the parameters corresponding to the components previously found. */
  internal_fitvonmises_nrcomponents = components->nrcomponents + 1 - refine;
  for(i = 0; i < components->nrcomponents; i++) {
    xstart[3*i+internal_fitvonmises_fitbaseline]   = components->centre[i];
    xstart[3*i+1+internal_fitvonmises_fitbaseline] = components->concentration[i];
    xstart[3*i+2+internal_fitvonmises_fitbaseline] = components->height[i];
    dx[3*i+internal_fitvonmises_fitbaseline]   = 0.01;
    dx[3*i+1+internal_fitvonmises_fitbaseline] = 0.1*xstart[3*i+1+internal_fitvonmises_fitbaseline];
    dx[3*i+2+internal_fitvonmises_fitbaseline] = 0.1*xstart[3*i+2+internal_fitvonmises_fitbaseline];
    fixed[3*i+internal_fitvonmises_fitbaseline]   = 1-refine;
    fixed[3*i+1+internal_fitvonmises_fitbaseline] = 1-refine;
    fixed[3*i+2+internal_fitvonmises_fitbaseline] = 1-refine;
  }
  if(internal_fitvonmises_fitbaseline) {
    xstart[0] = 0;
    dx[0] = 0;
    fixed[0] = 0;
  }

  /* If only refining there are no trials */
  if(refine) {
    fitvonMises_Nr_X_trials = 1;
    fitvonMises_Nr_C_trials = 1;
  }

  /* If educateguess there are no position trials */
  if(educatedguess) {
    fitvonMises_Nr_X_trials = 1;
  }

  bin = 0;
  pulsewidth = 0;
  snrbest = 0;
  E_best = 0;

  int have_onpulse = 0;
  if(onpulse != NULL) {
    if(onpulse->nrRegions > 0) {
      have_onpulse = 1;
    }
  }
  int posOrNeg;
  if(internal_fitvonmises_avoid_neg_components) {
    posOrNeg = 0;
  }else {
    posOrNeg = 1;
  }
  for(ix = 0; ix < fitvonMises_Nr_X_trials; ix++) {
    for(ic = 0; ic < fitvonMises_Nr_C_trials; ic++) {
      /* Define parameters of the new component */
      if(refine == 0) {
	if(educatedguess) {
	  if(ic == 0) {
	    if(have_onpulse == 0) {
	      boxcarFindpeak(internal_fitvonmises_fitprofile, internal_fitvonmises_nrbins, NULL, &bin, &pulsewidth, &snrbest, &E_best, 0, posOrNeg, 0, 1, -1, have_onpulse, 0, verbose);
	    }else {
	      boxcarFindpeak(internal_fitvonmises_fitprofile, internal_fitvonmises_nrbins, onpulse, &bin, &pulsewidth, &snrbest, &E_best, 0, posOrNeg, 0, 1, -1, have_onpulse, 0, verbose);
	    }
	  }
	  if(superverbose) {
	    printf("boxcar (ix=%d ic=%d): bin=%d width=%d snr=%f E=%f\n", ix, ic, bin, pulsewidth, snrbest, E_best); 
	  }
	  xstart[3*(internal_fitvonmises_nrcomponents - 1)+internal_fitvonmises_fitbaseline]   = (bin+0.5*pulsewidth)/(float)internal_fitvonmises_nrbins;
	  dx[3*(internal_fitvonmises_nrcomponents - 1)+internal_fitvonmises_fitbaseline]   = 0.05;
	  xstart[3*(internal_fitvonmises_nrcomponents - 1)+2+internal_fitvonmises_fitbaseline] = E_best/(float)pulsewidth;
	}else {
	  xstart[3*(internal_fitvonmises_nrcomponents - 1)+internal_fitvonmises_fitbaseline]   = ix/(float)fitvonMises_Nr_X_trials;
	  dx[3*(internal_fitvonmises_nrcomponents - 1)+internal_fitvonmises_fitbaseline]   = 0.5/(float)fitvonMises_Nr_X_trials;
	  xstart[3*(internal_fitvonmises_nrcomponents - 1)+2+internal_fitvonmises_fitbaseline] = ymax;
	}
	xstart[3*(internal_fitvonmises_nrcomponents - 1)+1+internal_fitvonmises_fitbaseline] = 1.0+3.5*ic/(float)(fitvonMises_Nr_C_trials-1);
	xstart[3*(internal_fitvonmises_nrcomponents - 1)+1+internal_fitvonmises_fitbaseline] = pow(10,xstart[3*(internal_fitvonmises_nrcomponents - 1)+1+internal_fitvonmises_fitbaseline]);

	dx[3*(internal_fitvonmises_nrcomponents - 1)+1+internal_fitvonmises_fitbaseline] = xstart[3*(internal_fitvonmises_nrcomponents - 1)+1+internal_fitvonmises_fitbaseline]/2.0;
	dx[3*(internal_fitvonmises_nrcomponents - 1)+2+internal_fitvonmises_fitbaseline] = fabs(xstart[3*(internal_fitvonmises_nrcomponents - 1)+2+internal_fitvonmises_fitbaseline]/3.5);

	fixed[3*(internal_fitvonmises_nrcomponents - 1)+internal_fitvonmises_fitbaseline]   = 0;
	fixed[3*(internal_fitvonmises_nrcomponents - 1)+1+internal_fitvonmises_fitbaseline] = 0;
	fixed[3*(internal_fitvonmises_nrcomponents - 1)+2+internal_fitvonmises_fitbaseline] = 0;
      
	if(internal_fitvonmises_fitbaseline) {
	  xstart[0] = av;
	  dx[0] = av*0.1;
	  fixed[0] = 0;
	}
	if(superverbose) {
	  if((ix == 0 && ic == 0) || chi2 < chi2_best) {
	    printf("xstart = pos=%f con=%f amp=%f (ymax=%f)\n", xstart[3*(internal_fitvonmises_nrcomponents - 1)+internal_fitvonmises_fitbaseline], xstart[3*(internal_fitvonmises_nrcomponents - 1)+1+internal_fitvonmises_fitbaseline], xstart[3*(internal_fitvonmises_nrcomponents - 1)+2+internal_fitvonmises_fitbaseline], ymax);
	    printf("dx = %f %f %f\n", dx[3*(internal_fitvonmises_nrcomponents - 1)+internal_fitvonmises_fitbaseline], dx[3*(internal_fitvonmises_nrcomponents - 1)+1+internal_fitvonmises_fitbaseline], dx[3*(internal_fitvonmises_nrcomponents - 1)+2+internal_fitvonmises_fitbaseline]);
	    printf("fixed = %d %d %d\n", fixed[3*(internal_fitvonmises_nrcomponents - 1)+internal_fitvonmises_fitbaseline], fixed[3*(internal_fitvonmises_nrcomponents - 1)+1+internal_fitvonmises_fitbaseline], fixed[3*(internal_fitvonmises_nrcomponents - 1)+2+internal_fitvonmises_fitbaseline]);


	    /*	    ppgclos();
      pgplot_viewport_def viewport;
      pgplot_clear_viewport_def(&viewport);
      strcpy(viewport.plotDevice, "1242353/xs");
      viewport.dontclose = 1;
	    if(pgplotGraph1(viewport, internal_fitvonmises_fitprofile, NULL, NULL, internal_fitvonmises_nrbins, 0, 360, 0, 360, 0, "Pulse longitude", "Intensity", "", 0, 0, 1, 1, NULL, -1, debug) == 0) {
    fflush(stdout);
	      printerror(verbose.debug, "fitvonMises: ERROR cannot open plot device\n");
	      return 0;
	    }
	    int dummy;
	    scanf("%d", &dummy);
	    */
	  }
	}
      }

      if(superverbose > 2)
	internal_fitvonmises_plotallsteps = 1;
      else
	internal_fitvonmises_plotallsteps = 0;

      /* Do the fitting */
      ret = doAmoeba(AmoebaAlgorithm, xstart, dx, fixed, xfit, &chi2, 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline, &internal_fitvonmises_funk, 1e-3/(1.0+9.0*refine), &nfunk, 0*verbose.verbose, 0, 3.0, NULL, NULL);
      if(ret != 0) {
	fflush(stdout);
	printerror(verbose.debug, "fitvonMises: Amoeba failed (return code = %d).", ret);
	free(xstart);
	free(dx);
	free(xfit);
	free(fixed);
	free(best);
	return 0;      
      }

      if(superverbose) {
	printf("Trial x=%f y=%f: chi2=%f after %d steps.\n", xstart[3*(internal_fitvonmises_nrcomponents - 1)+internal_fitvonmises_fitbaseline], xfit[3*(internal_fitvonmises_nrcomponents - 1)+2+internal_fitvonmises_fitbaseline], chi2, nfunk);
	printf("    x=%f c=%f h=%f chi=%f\n", xfit[3*(internal_fitvonmises_nrcomponents-1)+internal_fitvonmises_fitbaseline], xfit[3*(internal_fitvonmises_nrcomponents-1)+1+internal_fitvonmises_fitbaseline], xfit[3*(internal_fitvonmises_nrcomponents-1)+2+internal_fitvonmises_fitbaseline], chi2);
      }

      if(superverbose > 1) {
	internal_fitvonmises_plotallsteps = 1;
	internal_fitvonmises_funk(xfit); 
	internal_fitvonmises_plotallsteps = 0;
      }

      if((ix == 0 && ic == 0) || chi2 < chi2_best) {
	chi2_best = chi2;
	best->nrcomponents = internal_fitvonmises_nrcomponents;
	/* 	if(verbose.verbose) printf("  Found\n"); */
	for(i = 0; i < internal_fitvonmises_nrcomponents; i++) {
	  best->centre[i] = xfit[3*i+internal_fitvonmises_fitbaseline];
	  best->concentration[i] = xfit[3*i+1+internal_fitvonmises_fitbaseline];
	  best->height[i] = xfit[3*i+2+internal_fitvonmises_fitbaseline];
	  if(internal_fitvonmises_fitbaseline) {
	    best_baseline = xfit[0];
	    if(superverbose) {
	      printf("Baseline = %f\n", best_baseline);
	    }
	  }
	  /*	  if(verbose.verbose) printf("    %f %f %f\n", xfit[3*i], xfit[3*i+1], xfit[3*i+2]); */
	}
      }
    }
  }

  if(refine == 0) {
    /* Calculate new residual profile */  
    calcVonMisesProfile(best, internal_fitvonmises_nrbins, internal_fitvonmises_fitprofile, 0, 0);
    for(i = 0; i < internal_fitvonmises_nrbins; i++) {
      internal_fitvonmises_fitprofile[i] = internal_fitvonmises_profile[i]-internal_fitvonmises_fitprofile[i];
    }
    rms = 0;
    pwr = 0;
    n = 0;
    for(i = 0; i < internal_fitvonmises_nrbins; i++) {
      rms += internal_fitvonmises_fitprofile[i]*internal_fitvonmises_fitprofile[i];
      y = calcVonMisesFunction2(best->centre[internal_fitvonmises_nrcomponents - 1], best->concentration[internal_fitvonmises_nrcomponents - 1], best->height[internal_fitvonmises_nrcomponents - 1], i/(float)internal_fitvonmises_nrbins, 0);
      y = fabs(y);
      pwr += y;
      if(y >= fabs(0.5*best->height[internal_fitvonmises_nrcomponents - 1]))
	n++;
    }
    rms = sqrt(rms/(float)internal_fitvonmises_nrbins);
    if(n == 0)
      n = 1;
    s2n = pwr/(rms*sqrt(n));
    
    if(s2n < sigmalimit) {
      if(verbose.verbose) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("Discarded component: pwr=%f fwhm=%d s2n=%f (< %f) rms=%e\n", pwr, n, s2n, sigmalimit, rms);
      }
      if(superverbose) {
	for(i = 0; i < internal_fitvonmises_nrcomponents; i++) {
	  printf("    x=%f c=%f h=%f\n", best->centre[i], best->concentration[i], best->height[i]);
	}
      }
      free(xstart);
      free(dx);
      free(xfit);
      free(fixed);
      free(best);
      return 99;
    }
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("New component:       pwr=%f fwhm=%d s2n=%f rms=%e\n", pwr, n, s2n, rms);
    }
  }

  if(internal_fitvonmises_fitbaseline)
    *baseline = best_baseline;
  memcpy(components, best, sizeof(vonMises_collection_definition));
  free(xstart);
  free(dx);
  free(xfit);
  free(fixed);
  free(best);
  return 1;
}

/* Fit a vonMises function through the pulse profile in
   psrdata. MaxNrComp defines the maximum nr of components which are
   allowed to be fitted and sigmalimit the limit up to where
   components are fitted (5 seems to work reasonably good). If
   educatedguess it will use the boxcarFindpeak algorithm to find
   initial guesses (this option seems to work best, so you should use
   it by default). If avoid_neg_components is non-zero, then
   components with negative amplitudes are avoided. If fitbaseline is
   set it is not assumed that the baseline is removed and its level is
   stored in baseline. The larger precision is set to the finer trials
   are tried. Something like 2 seems reasonable. If onpulse != NULL
   (and the number of regions is non-zero), the components must be
   centered within the specified range.  verbose-1 is number of spaces
   before output.

   superverbose = 1: extra print statements
                  2: plot best intermediate steps
                  3: plot all intermediate steps

Return 1 is success 
*/
int fitvonmises(datafile_definition psrdata, vonMises_collection_definition *components, int MaxNrComp, float sigmalimit, int educatedguess, int fitbaseline, int precision, int avoid_neg_components, float *baseline, pulselongitude_regions_definition *onpulse, verbose_definition verbose, int superverbose, char *plotdevice)
{
  int n, ret, fitvonMises_Nr_X_trials, fitvonMises_Nr_C_trials, moreaccurate;
  verbose_definition verbosedebug;
  pgplot_options_definition *pgplot_options;
  copyVerboseState(verbose, &verbosedebug);

  internal_fitvonmises_debug = verbose.debug;
  internal_fitvonmises_avoid_neg_components = avoid_neg_components;
  internal_fitvonmises_onpulse = onpulse;
  if(onpulse != NULL) {
    region_frac_to_int(onpulse, psrdata.NrBins, 0);
  }

  if(verbose.debug == 0) {
    verbosedebug.verbose = 0;
  }


  /*
  memset(components, 0, sizeof(vonMises_def));
  */

  fitvonMises_Nr_X_trials = 4;
  fitvonMises_Nr_C_trials = 2;
  moreaccurate = 0;
  if(fitbaseline)
    internal_fitvonmises_fitbaseline = 1;
  else
    internal_fitvonmises_fitbaseline = 0;
  internal_fitvonmises_fitdummy = 0;  // Only used in fitvonmises_refine_model()
  if(precision < 1)
    precision = 1;
  if(psrdata.NrFreqChan != 1 || psrdata.NrSubints != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitvonmises: This is not a pulse profile.");
    return 0;
  }
  if(psrdata.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitvonmises: only works if data is loaded into memory.");
    return 0;
  }

  internal_fitvonmises_fitprofile = malloc(psrdata.NrBins*sizeof(float));
  internal_fitvonmises_profile = psrdata.data;
  internal_fitvonmises_plotallsteps = 0;
  if(internal_fitvonmises_fitprofile == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "fitvonMises: Memory allocation error.\n");
    return 0;
  }

  components->nrcomponents = 0;
  internal_fitvonmises_nrbins = psrdata.NrBins;

  if(plotdevice != NULL) {
    pgplot_options = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
    if(pgplot_options == NULL) {
      printerror(verbose.debug, "ERROR fitvonMises: Memory allocation error");
      return 0;
    }
    pgplot_clear_options(pgplot_options);
    strcpy(pgplot_options->box.xlabel, "Pulse longitude");
    strcpy(pgplot_options->box.ylabel, "Intensity");
    strcpy(pgplot_options->box.title, "");
    strcpy(pgplot_options->viewport.plotDevice, plotdevice);
    pgplot_options->viewport.dontopen  = 0;
    pgplot_options->viewport.dontclose = 1;
    if(pgplotGraph1(pgplot_options, psrdata.data, NULL, NULL, psrdata.NrBins, 0, 360, 0, 0, 360, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, verbosedebug) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "fitvonMises: ERROR cannot open plot device");
      return 0;
    }
  }

  for(n = 0; n < MaxNrComp; n++) {
    if(plotdevice != NULL) {
      ppgpage();
      pgplot_options->viewport.dontopen  = 1;
      pgplot_options->viewport.dontclose = 1;
      if(pgplotGraph1(pgplot_options, psrdata.data, NULL, NULL, psrdata.NrBins, 0, 360, 0, 0, 360, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, verbosedebug) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "fitvonMises: ERROR cannot open plot device");
	return 0;
      }
    }
    ret = fitvonmises_addcomponent(components, 0, educatedguess, sigmalimit, fitvonMises_Nr_X_trials, fitvonMises_Nr_C_trials, baseline, onpulse, verbose, superverbose);
    if(ret == 0) {
      fflush(stdout);
      printerror(verbose.debug, "fitvonMises: ERROR while adding a new component");
      return 0;
    }
    if(ret == 99) {
  /* Makes finer trials when failing to find a component. The larger
     the number the more accurate */
      if(moreaccurate < precision) {
	moreaccurate++;
	fitvonMises_Nr_X_trials *= 3;
	fitvonMises_Nr_C_trials *= 2;
	n -= 1;
      }else {
	break;
      }
    }else {
      if(plotdevice != NULL) {
	ppgpage();
	pgplot_options->viewport.dontopen  = 1;
	pgplot_options->viewport.dontclose = 1;
	if(pgplotGraph1(pgplot_options, psrdata.data, NULL, NULL, psrdata.NrBins, 0, 360, 0, 0, 360, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, verbosedebug) == 0) {
	  fflush(stdout);
	  printerror(verbose.debug, "fitvonMises: ERROR cannot open plot device");
	  return 0;
	}
      }
      if(fitvonmises_addcomponent(components, 1, educatedguess, sigmalimit, fitvonMises_Nr_X_trials, fitvonMises_Nr_C_trials, baseline, onpulse, verbose, superverbose) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "fitvonMises: ERROR while adding a new component");
	return 0;
      }
    }
  }

  vonMises_simplify_parameters(components);

  if(verbose.verbose) {
    for(n = 0; n < components->nrcomponents; n++) {
      printf("    comp %d: x=%f c=%f h=%f\n", n+1, components->centre[n], components->concentration[n], components->height[n]);
    }
  }

  if(plotdevice != NULL) {
    free(pgplot_options);
    ppgclos();
  }
  free(internal_fitvonmises_fitprofile);
  return 1;
}



/* Fit a given vonMises function through the pulse profile in psrdata
   and refine its parameters. If fitbaseline is set it is not assumed
   that the baseline is removed and its level is stored in
   baseline. If avoid_neg_components is non-zero, then components with
   negative amplitudes are avoided.  If fixamp, fixwidth and fixphase
   are set, the amplitudes, widths and/or phases are held fixed. If
   fixrelphase is set, the components cannot shift relative to each
   other. With fixrelamp, the components cannot change height relative
   to each other. If there is only one fit parameter (for instance
   when fitting an overall scaling, without a baseline variation), a
   dummy fit parameter is introduced to avoid the downhill-simplex
   method to fail. Not elegant, but it works.

   return 0=error,    1=ok */
int fitvonmises_refine_model(datafile_definition psrdata, vonMises_collection_definition *components, int fitbaseline, int avoid_neg_components, float *baseline, int fixamp, int fixwidth, int fixphase, int fixrelamp, int fixrelphase, verbose_definition verbose)
{
  float *xstart, *dx, *xfit;
  float chi2, av; //, ymax, ymin, maxmodelheight;
  int *fixed, nfunk;
  int i, superverbose;


  internal_fitvonmises_debug = verbose.debug;
  internal_fitvonmises_avoid_neg_components = avoid_neg_components;

  internal_fitvonmises_relative_phases = fixrelphase;
  internal_fitvonmises_relative_amps   = fixrelamp;
  superverbose = 0;
  if(verbose.verbose) {
    printf("fitvonmises_refine_model: Start refining model with %d components.\n", components->nrcomponents);
  }

  if(fixamp && fixwidth && fixphase) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitvonmises_refine_model: Cannot simultaneously fix the phases, concentrations and amplitudes of the components.");
    return 0;
  }

  if(psrdata.NrFreqChan != 1 || psrdata.NrSubints != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitvonmises_refine_model: This is not a pulse profile.");
    return 0;
  }
  if(psrdata.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitvonmises_refine_model: only works if data is loaded into memory.");
    return 0;
  }

  internal_fitvonmises_nrbins = psrdata.NrBins;
  internal_fitvonmises_fitprofile = malloc(psrdata.NrBins*sizeof(float));
  internal_fitvonmises_profile = psrdata.data;
  if(internal_fitvonmises_fitprofile == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fitvonMises_refine_model: Memory allocation error.");
    return 0;
  }


  xstart = malloc((3*maxNrVonMisesComponents+1)*sizeof(float));
  dx     = malloc((3*maxNrVonMisesComponents+1)*sizeof(float));
  xfit   = malloc((3*maxNrVonMisesComponents+1)*sizeof(float));
  fixed  = malloc((3*maxNrVonMisesComponents+1)*sizeof(int));
  if(xstart == NULL || dx == NULL || xfit == NULL || fixed == NULL) {
    printerror(verbose.debug, "ERROR fitvonmises_refine_model: Cannot allocate memory");
    return 0;
  }

  if(fitbaseline)
    internal_fitvonmises_fitbaseline = 1;
  else
    internal_fitvonmises_fitbaseline = 0;

  // Not sure why scaling is based on what is in the residual profile? Is causing problems. Scaling now done if fitvonMises. Could move scaling from there to this function.
  av = 0;
  for(i = 0; i < internal_fitvonmises_nrbins; i++) {
    av += internal_fitvonmises_profile[i];
  }
  av /= (float)internal_fitvonmises_nrbins;

  // Scale model if required
  if(fixamp == 0) {  // Scale the model so it roughly fits the input data.
    float min, max, value, scale;
    min = max = psrdata.data[0];
    for(i = 0; i < psrdata.NrBins; i++) {
      if(psrdata.data[i] < min)
	min = psrdata.data[i];
      if(psrdata.data[i] > max)
	max = psrdata.data[i];
    }
    scale = max - min;  // Estimation of height of profile
    for(i = 0; i < psrdata.NrBins; i++) {
      value = calcVonMisesFunction(components, i/(float)psrdata.NrBins, 0);
      if(i == 0 || value > max)
	max = value;
    }
    scale /= max; // Express height profile in terms of units height template
    if(verbose.debug)
      printf("Scaling input model amplitudes with factor %e\n", scale);
    for(i = 0; i < components->nrcomponents; i++) {
      components->height[i] *= scale;
    }
  }

  if(fixphase == 0) {   // Shift the model so it roughly fits the input data.
    float value;
    value = correlateVonMisesFunction(components, psrdata.NrBins, psrdata.data, verbose);
    if(verbose.debug)
      printf("Shifting input model %e in phase\n", value);
    for(i = 0; i < components->nrcomponents; i++) {
      components->centre[i] += value;
    }
  }


  // Calculate residual profile
  /*
  calcVonMisesProfile(*components, internal_fitvonmises_nrbins, internal_fitvonmises_fitprofile, 0, 0);
  for(i = 0; i < internal_fitvonmises_nrbins; i++) 
    internal_fitvonmises_fitprofile[i] = internal_fitvonmises_profile[i]-internal_fitvonmises_fitprofile[i];
  // Determine maximum aplitude
  ymax = internal_fitvonmises_fitprofile[0];
  ymin = internal_fitvonmises_fitprofile[0];
  av = 0;
  for(i = 0; i < internal_fitvonmises_nrbins; i++) {
    av += internal_fitvonmises_profile[i];
    if(internal_fitvonmises_fitprofile[i] > ymax)
      ymax = internal_fitvonmises_fitprofile[i];
    if(internal_fitvonmises_fitprofile[i] < ymin)
      ymin = internal_fitvonmises_fitprofile[i];
  }
  av /= (float)internal_fitvonmises_nrbins;
  if(fabs(ymin) > fabs(ymax))
    ymax = ymin;

  maxmodelheight = fabs(components->height[0]);
  for(i = 0; i < components->nrcomponents; i++) {
    if(fabs(components->height[i]) > maxmodelheight)
      maxmodelheight = fabs(components->height[i]);
  }

  if(verbose.verbose) {
    printf("    Residual: amplitude=%f, av=%f\n", ymax, av);
    printf("    Scaling model by %f\n", ymax/maxmodelheight);
  }

  if(fixamp == 0) {
    for(i = 0; i < components->nrcomponents; i++) {
      components->height[i] *= ymax/maxmodelheight;
    }
  }
  */

  if(internal_fitvonmises_fitbaseline) {
    xstart[0] = av;
    dx[0] = av*0.1;
    fixed[0] = 0;
  }

  /* Set the parameters corresponding to the components previously found. */
  internal_fitvonmises_nrcomponents = components->nrcomponents;
  int nrfitparameters;
  nrfitparameters = 0;
  for(i = 0; i < components->nrcomponents; i++) {
    xstart[3*i+internal_fitvonmises_fitbaseline]   = components->centre[i];
    if(fixrelphase && i > 0)
      xstart[3*i+internal_fitvonmises_fitbaseline] -= components->centre[0];

    xstart[3*i+1+internal_fitvonmises_fitbaseline] = components->concentration[i];

    xstart[3*i+2+internal_fitvonmises_fitbaseline] = components->height[i];
    if(fixrelamp && i > 0)
      xstart[3*i+2+internal_fitvonmises_fitbaseline] /= components->height[0];

    if(verbose.debug) {
      printf("Refine component %d with phi=%f, con=%f, amp=%f\n", i+1, components->centre[i], components->concentration[i], components->height[i]);
    }
    if(fixphase || (fixrelphase && i > 0)) {
      fixed[3*i+internal_fitvonmises_fitbaseline]   = 1;
      dx[3*i+internal_fitvonmises_fitbaseline]   = 0;
    }else {
      fixed[3*i+internal_fitvonmises_fitbaseline]   = 0;
      dx[3*i+internal_fitvonmises_fitbaseline]   = 0.01;
      nrfitparameters++;
    }
    if(fixwidth) {
      dx[3*i+1+internal_fitvonmises_fitbaseline] = 0;
      fixed[3*i+1+internal_fitvonmises_fitbaseline] = 1;
    }else {
      dx[3*i+1+internal_fitvonmises_fitbaseline] = 0.1*xstart[3*i+1+internal_fitvonmises_fitbaseline];
      fixed[3*i+1+internal_fitvonmises_fitbaseline] = 0;
      nrfitparameters++;
    }
    if(fixamp || (fixrelamp && i > 0)) {
      fixed[3*i+2+internal_fitvonmises_fitbaseline] = 1;
      dx[3*i+2+internal_fitvonmises_fitbaseline] = 0;
    }else {
      fixed[3*i+2+internal_fitvonmises_fitbaseline] = 0;
      dx[3*i+2+internal_fitvonmises_fitbaseline] = 0.1*xstart[3*i+2+internal_fitvonmises_fitbaseline];
      nrfitparameters++;
    }
  }

  if(nrfitparameters == 1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING fitvonMises_refine_model: Adding a dummy variable to construct a multi-parameter fit as required by the downhill-simplex method.");
    i = 3*components->nrcomponents+internal_fitvonmises_fitbaseline;
    internal_fitvonmises_fitdummy = 1;
    xstart[i] = 0;
    dx[i] = 0;
    fixed[i] = 0;
  }else {
    internal_fitvonmises_fitdummy = 0;
  }


  if(superverbose > 2)
    internal_fitvonmises_plotallsteps = 1;
  else
    internal_fitvonmises_plotallsteps = 0;

  if(verbose.verbose) {
    printf("    Start the fitting with %d params, of which ", 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy);
    int n = 0;
    for(i = 0; i < 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy; i++) {
      if(fixed[i] == 0) {
	n++;
      }
    }
    printf("%d are free parameters.", n);
    if(internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy) {
      printf(" This includes ");
      if(internal_fitvonmises_fitbaseline) {
	printf("a baseline offset");
	if(internal_fitvonmises_fitdummy) {
	  printf(" and");
	}
      }
      if(internal_fitvonmises_fitdummy) {
	printf("a dummy variable");
      }
    }
    printf("\n");
    if(verbose.debug) {
      printf("Initial values are:     ");
      for(i = 0; i < 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy; i++) {
	if(fixed[i] == 0) {
	  printf("%e ", xstart[i]);
	}
      }
      printf("\n");
      printf("Initial step sizes are: ");
      for(i = 0; i < 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy; i++) {
	if(fixed[i] == 0) {
	  printf("%e ", dx[i]);
	}
      }
      printf("\n");
    }
  }
  /* Do the fitting */
  int ret;
  ret = doAmoeba(AmoebaAlgorithm, xstart, dx, fixed, xfit, &chi2, 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline+internal_fitvonmises_fitdummy, &internal_fitvonmises_funk, 1e-4, &nfunk, 0*verbose.verbose, 0, 3.0, NULL, NULL);
  if(ret != 0) {
    fflush(stdout);
    printerror(verbose.debug, "fitvonMises_refine_model: Amoeba failed.");
    if(ret == 1) {
      printerror(verbose.debug, "fitvonMises_refine_model: Maximum number of itterations exceeded.");
    }else if(ret == 2) {
      printerror(verbose.debug, "fitvonMises_refine_model: Memory allocation error.");
    }else if(ret == 2) {
      printerror(verbose.debug, "fitvonMises_refine_model: Not enough free parameters to perform fit.");
    }else if(ret == 2) {
      printerror(verbose.debug, "fitvonMises_refine_model: Algorithm is not available.");
    }else if(ret == 2) {
      printerror(verbose.debug, "fitvonMises_refine_model: Other fit fail error (such as singular matrix).");
    }else if(ret == 2) {
      printerror(verbose.debug, "fitvonMises_refine_model: Undefined error code.");
    }
    free(xstart);
    free(dx);
    free(xfit);
    free(fixed);
    return 0;      
  }

  if(verbose.verbose) {
    printf("    chi2=%f after %d steps.\n", chi2, nfunk);
  }

  double offset;
  offset = xfit[internal_fitvonmises_fitbaseline]-components->centre[0];
  for(i = 0; i < internal_fitvonmises_nrcomponents; i++) {
    //    printf("XXXXXX changing comp %d from %f to %f\n", i, components->centre[i], xfit[3*i+internal_fitvonmises_fitbaseline]);
    if(fixrelphase) {
      components->centre[i] += offset;
    }else {
      components->centre[i] = xfit[3*i+internal_fitvonmises_fitbaseline];
    }
    //    if(fixrelphase && i > 0) {
    //      printf("XXXXXX changing comp %d by adding %f to %f\n", i, components->centre[internal_fitvonmises_fitbaseline], components->centre[i]);
    //      components->centre[i] += components->centre[internal_fitvonmises_fitbaseline];
    //    }
    components->concentration[i] = xfit[3*i+1+internal_fitvonmises_fitbaseline];
    components->height[i] = xfit[3*i+2+internal_fitvonmises_fitbaseline];
    if(fixrelamp && i > 0)
      components->height[i] *= xfit[2+internal_fitvonmises_fitbaseline];
  }
  if(internal_fitvonmises_fitbaseline) {
    *baseline = xfit[0];
    if(verbose.verbose) {
      printf("    Baseline = %f\n", *baseline);
    }
  }
  if(internal_fitvonmises_fitdummy) {
    i = 3*internal_fitvonmises_nrcomponents+internal_fitvonmises_fitbaseline;
    if(verbose.verbose) {
      printf("    Dummy = %e (should be zero)\n", xfit[i]);
    }
  }

  free(internal_fitvonmises_fitprofile);
  if(verbose.verbose) {
    printf("fitvonmises_refine_model: done.\n");
  }
  internal_fitvonmises_relative_phases = 0;
  internal_fitvonmises_relative_amps = 0;

  vonMises_simplify_parameters(components);

  free(xstart);
  free(dx);
  free(xfit);
  free(fixed);
  return 1;
}



//START REGION DEVELOP

/* Calculate a shape parameter based on the profiled model after
   applying a phase shift. phaseprecision, a number << 1, indicates
   how precise the shape parameter needs to be. More precision = more
   computations. The possibilities for shapepar are defined in
   psrsalsa_defines.h, and start with SHAPEPAR_. The shape parameter
   measure might depend on one or mor additional input values (a range
   for example), which are provided as shapepar_aux. The measurement
   is returned as measurement.
*/
void calcVonMisesProfile_shape_parameter(vonMises_collection_definition *components, double shift, double phaseprecision, int shapepar, double *shapepar_aux, double *measurement, verbose_definition verbose)
{
  int cont;
  long i, nrbins, bin1, bin2;
  double x, y, xmax, ymax;
  
  if(verbose.debug) {
    printf("Find shape parameter %d\n", shapepar);
  }
  if(components->nrcomponents <= 0) {
    printerror(verbose.debug, "ERROR calcVonMisesProfile_shape_parameter: Analytic model is undefined, so no shape parameter could be measured.");
    *measurement = sqrt(-1);
    return;
  }

  // Do a rough search for overall peak location, as that knowledge is required for all shape parameters
  nrbins = 2048;
  bin1 = 0;
  bin2 = nrbins-1;
  do {
    for(i = bin1; i < bin2; i++) {
      x = i/(double)nrbins;
      y = calcVonMisesFunction(components, x, shift);
      if(i == 0 || y > ymax) {
	ymax = y;
	xmax = x;
      }
    }
    if(verbose.debug) {
      printf("  Found maximum at phase %lf (max=%e)\n", xmax, ymax);
    }
    if(1.0/phaseprecision > nrbins)
      cont = 1;
    else
      cont = 0;
    nrbins *= 10;
    bin1 = xmax*nrbins-20;
    bin2 = xmax*nrbins+20;
  }while(cont);
  if(shapepar == SHAPEPAR_PEAKPHASE) {
    *measurement = xmax;
    return;
  }else if(shapepar == SHAPEPAR_PEAKAMP) {
    *measurement = ymax;
    return;
  }

  // Search for widths at given intensity level
  double yrequest, xedge1, xedge2;
  int search_edges = 1;
  if(shapepar == SHAPEPAR_W10) {
    yrequest = 0.1*ymax;
  }else if(shapepar == SHAPEPAR_W25) {
    yrequest = 0.25*ymax;
  }else if(shapepar == SHAPEPAR_W50) {
    yrequest = 0.5*ymax;
  }else if(shapepar == SHAPEPAR_W75) {
    yrequest = 0.75*ymax;
  }else if(shapepar == SHAPEPAR_W90) {
    yrequest = 0.9*ymax;
  }else {
    search_edges = 0;
  }

  /*
    printerror(verbose.debug, "ERROR calcVonMisesProfile_shape_parameter: Requested shape parameter is not implemented");
    *measurement = sqrt(-1);
    return;
  */
  // I think the same code is repeated below, so make changes in both if required
  if(search_edges) {
    int direction;
    int found;
    for(direction = -1; direction <= 1; direction += 2) {  // Do a search in two directions to find both edges
      found = 0;
      nrbins = 2048;
      if(direction == -1) {
	bin1 = (xmax+0.5)*nrbins;  // Start half a turn after maximum, end half a turn earlier
	bin2 = (xmax-0.5)*nrbins;
      }else {
	bin1 = (xmax-0.5)*nrbins;  // Start half a turn before maximum, end half a turn later
	bin2 = (xmax+0.5)*nrbins;
      }
      do {
	// Need to find continue until intensity drops below yrequest (i.e. if there is an interpulse) -> found = 1 
	// Need then to continue until intensity exceeds yrequest -> found = 2
	i = bin1;
	do {
	  x = i/(double)nrbins;
	  y = calcVonMisesFunction(components, x, shift);
	  //	printf("XXXX bin=%ld (search range %ld .. %ld) x=%e y=%e found=%d\n", i, bin1, bin2, x, y, found);
	  if(found == 0 && y <= yrequest) {
	    found = 1;
	  }
	  if(found == 1 && y >= yrequest) {
	    found = 2;
	    break;
	  }
	  i += direction;
	  if(direction == 1) {
	    if(i <= bin2)
	      cont = 1;
	    else
	      cont = 0;
	  }else {
	    if(i >= bin2)
	      cont = 1;
	    else
	      cont = 0;
	  }
	}while(cont);
	if(found < 2) {
	  break;
	}
	if(verbose.debug) {
	  printf("  Found edge at phase %lf (y=%e is close to %e)\n", x, y, yrequest);
	}
	if(1.0/phaseprecision > nrbins) {
	  cont = 1;
	  found = 0;
	}else {
	  cont = 0;
	}
	nrbins *= 10;
	if(direction == 1) {
	  bin1 = x*nrbins-20;
	  bin2 = x*nrbins+20;
	}else {
	  bin1 = x*nrbins+20;
	  bin2 = x*nrbins-20;
	}
      }while(cont);
      if(verbose.debug && direction == -1) {
	printf("  Found one edge, try to find next\n");
      }
      if(found < 2)
	break;
      if(direction == -1) {
	xedge2 = x;
      }else if(direction == 1) {
	xedge1 = x;
      }
    }  // Other direction
    if(found == 0) {
      printwarning(verbose.debug, "WARNING calcVonMisesProfile_shape_parameter: Requested width at intensity level (%e) is 360 deg since the intensity level appears to be always higher.\n", yrequest);
      *measurement = 1;
      return;
    }else if(found == 1) {
      printwarning(verbose.debug, "WARNING calcVonMisesProfile_shape_parameter: Requested width at intensity level (%e) is 0 deg since the intensity level appears to be always lower.\n", yrequest);
      *measurement = 0;
      return;
    }else {
      *measurement = xedge2-xedge1;
      return;
    }
  }

  int search_other_peak = 0;
  if(shapepar == SHAPEPAR_PEAKSEPPHASE) {
    search_other_peak = 1;
  }else if(shapepar == SHAPEPAR_PEAKAMPRATIO || shapepar == SHAPEPAR_PEAKAMPRATIO_RECI) {
    search_other_peak = 1;
  }
  // I think the same code is repeated below, so make changes in both if required
  if(search_other_peak) {
    int found;
    found = 0;
    nrbins = 2048;
    if(shapepar_aux[1] > shapepar_aux[0]) {
      bin1 = (xmax+shapepar_aux[0])*nrbins;  // Start with given offset from main peak
      bin2 = (xmax+shapepar_aux[1])*nrbins;
    }else {
      bin1 = (xmax+shapepar_aux[1])*nrbins;  // Start with given offset from main peak
      bin2 = (xmax+shapepar_aux[0])*nrbins;
    }
    do {
      // Need to find a sign change in the derivative of the fit function.
      for(i = bin1; i <= bin2; i++) {
	double dydx;
	int oldsign;
	x = i/(double)nrbins;
	dydx = calcVonMisesFunction_deriv(components, x, shift);
	if(i == bin1) {
	  if(dydx >= 0)
	    oldsign = 1;
	  else
	    oldsign = -1;
	}
	if(dydx <= 0 && oldsign == 1) {
	  //	  printf("Found neg slope = %e, while it was positive before\n", dydx);
	  found = 1;
	  break;
	}
	if(oldsign == -1 && dydx >= 0)
	  oldsign = 1;
	//	printf("XXXX bin=%ld (search range %ld .. %ld) x=%e dydx=%e oldsign=%d found=%d\n", i, bin1, bin2, x, dydx, oldsign, found);
      }
      if(found  == 0) {
	break;
      }
      if(verbose.debug) {
	printf("  Found secondary maximum at phase %lf\n", x);
      }
      if(1.0/phaseprecision > nrbins) {
	cont = 1;
	found = 0;
      }else {
	cont = 0;
      }
      nrbins *= 10;
      bin1 = x*nrbins-20;
      bin2 = x*nrbins+20;
    }while(cont);
    if(found == 0) {
      printwarning(verbose.debug, "WARNING calcVonMisesProfile_shape_parameter: No peak found in phase offset range %lf .. %lf from phase of main peak at %lf. Shape parameter undefined.\n", shapepar_aux[0], shapepar_aux[1], xmax);
      *measurement = sqrt(-1);
      return;
    }else {
      if(shapepar == SHAPEPAR_PEAKSEPPHASE) {
	*measurement = x-xmax;
      }else if(shapepar == SHAPEPAR_PEAKAMPRATIO || shapepar == SHAPEPAR_PEAKAMPRATIO_RECI) {
	double y;
	y = calcVonMisesFunction(components, x, shift);
	if(shapepar == SHAPEPAR_PEAKAMPRATIO_RECI) {
	  *measurement = ymax/y;
	}else {
	  *measurement = y/ymax;
	}
      }
      return;
    }
  }

  if(shapepar == SHAPEPAR_W10_MAXAMP || shapepar == SHAPEPAR_W25_MAXAMP || shapepar == SHAPEPAR_W50_MAXAMP || shapepar == SHAPEPAR_W75_MAXAMP || shapepar == SHAPEPAR_W90_MAXAMP) {
    // Identify which component has the largest amplitude
    int nmax = 0;
    double amp_max = components->height[0];
    if(components->nrcomponents > 1) {
      int n;
      for(n = 1; n < components->nrcomponents; n++) {
	if(components->height[n] > amp_max) {
	  nmax = n;
	  amp_max = components->height[n];
	}
      }
    }
    if(shapepar == SHAPEPAR_W10_MAXAMP) {
      *measurement = widthVonMisesFunction2(components->concentration[nmax], 0.1)/(2.0*M_PI);
    }else if(shapepar == SHAPEPAR_W25_MAXAMP) {
      *measurement = widthVonMisesFunction2(components->concentration[nmax], 0.25)/(2.0*M_PI);
    }else if(shapepar == SHAPEPAR_W50_MAXAMP) {
      *measurement = widthVonMisesFunction2(components->concentration[nmax], 0.5)/(2.0*M_PI);
    }else if(shapepar == SHAPEPAR_W75_MAXAMP) {
      *measurement = widthVonMisesFunction2(components->concentration[nmax], 0.75)/(2.0*M_PI);
    }else if(shapepar == SHAPEPAR_W90_MAXAMP) {
      *measurement = widthVonMisesFunction2(components->concentration[nmax], 0.9)/(2.0*M_PI);
    }
    if(*measurement < 0) {
      *measurement = sqrt(-1);
    }
    return;
  }

  printerror(verbose.debug, "ERROR calcVonMisesProfile_shape_parameter: Requested shape parameter is not implemented");
  *measurement = sqrt(-1);
  return;
}

/* Print the shape parameter, plus a description if showdescr != 0. The error is quoted when >= 0. */
void print_shape_par(FILE *fout, int showdescr, int shapepar, double measurement, double error)
{
  if(showdescr) {
    switch(shapepar) {
    case SHAPEPAR_PEAKPHASE: fprintf(fout, "Pulse phase of peak: "); break;
    case SHAPEPAR_PEAKAMP: fprintf(fout, "Amplitude of peak: "); break;
    case SHAPEPAR_W10: fprintf(fout, "Width at 10%% of maximum amplitude (in phase): "); break;
    case SHAPEPAR_W25: fprintf(fout, "Width at 25%% of maximum amplitude (in phase): "); break;
    case SHAPEPAR_W50: fprintf(fout, "Width at 50%% of maximum amplitude (in phase): "); break;
    case SHAPEPAR_W75: fprintf(fout, "Width at 75%% of maximum amplitude (in phase): "); break;
    case SHAPEPAR_W90: fprintf(fout, "Width at 90%% of maximum amplitude (in phase): "); break;
    case SHAPEPAR_W10_MAXAMP: fprintf(fout, "Width at 10%% of amplitude vonMises tallest component (in phase): "); break;
    case SHAPEPAR_W25_MAXAMP: fprintf(fout, "Width at 25%% of amplitude vonMises tallest component (in phase): "); break;
    case SHAPEPAR_W50_MAXAMP: fprintf(fout, "Width at 50%% of amplitude vonMises tallest component (in phase): "); break;
    case SHAPEPAR_W75_MAXAMP: fprintf(fout, "Width at 75%% of amplitude vonMises tallest component (in phase): "); break;
    case SHAPEPAR_W90_MAXAMP: fprintf(fout, "Width at 90%% of amplitude vonMises tallest component (in phase): "); break;
    case SHAPEPAR_PEAKSEPPHASE: fprintf(fout, "Peak separation (in phase): "); break;
    case SHAPEPAR_PEAKAMPRATIO: fprintf(fout, "Peak amplitude ratio: "); break;
    case SHAPEPAR_PEAKAMPRATIO_RECI: fprintf(fout, "Reciprocal peak amplitude ratio: "); break;
    default: printerror(0, "ERROR print_shape_par: Requested shape parameter is not implemented"); break;
    }
  }
  if(shapepar == SHAPEPAR_PEAKPHASE || shapepar == SHAPEPAR_W10 || shapepar == SHAPEPAR_W25 || shapepar == SHAPEPAR_W50 || shapepar == SHAPEPAR_W75 || shapepar == SHAPEPAR_W90  || shapepar == SHAPEPAR_W10_MAXAMP  || shapepar == SHAPEPAR_W25_MAXAMP  || shapepar == SHAPEPAR_W50_MAXAMP  || shapepar == SHAPEPAR_W75_MAXAMP  || shapepar == SHAPEPAR_W90_MAXAMP || shapepar == SHAPEPAR_PEAKSEPPHASE) {
    if(error >= 0) {
      fprintf(fout, "%e +- %e (phase) = %e +- %e (deg)\n", measurement, error, measurement*360.0, error*360.0);
    }else {
      fprintf(fout, "%e (phase) = %e (deg)\n", measurement, measurement*360.0);
    }
  }else {
    if(error >= 0) {
      fprintf(fout, "%e +- %e\n", measurement, error);
    }else {
      fprintf(fout, "%e\n", measurement);
    }
  }
}

//START REGION DEVELOP
