//START REGION RELEASE
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_cdf.h>
#include "psrsalsa.h"
//START REGION DEVELOP
#include <stdio.h>
#include <stdlib.h>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_interp.h"
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>

//START REGION RELEASE
/* Generates a random unsigned int. There are 60,000,000 unique values
   that can be generated. Can be used for instance to generate unique
   temporary file names */
long randomUnsignedInt()
{
  time_t seconds;
  struct timeval precisetime;
  /*
    Get value from system clock and place in seconds variable
    (seconds since January 1, 1970).
  */
  time(&seconds);
  gettimeofday(&precisetime,0x0);

  return (long)seconds*(long)precisetime.tv_usec;
}
//START REGION DEVELOP



//START REGION RELEASE
/* Generates a random idnum to be used in Numerical Recipies random
   number generator. There are 60,000,000 unique values that can be
   genetated. */
void randomize_idnum(long *idnum)
{
  
  *idnum = -randomUnsignedInt();
}
//START REGION DEVELOP

/*
  Converts a number from a uniform distribution (between 0 and 1) to a
  number from a sinusoidal distribution (between 0 and 180, peaking at
  90).
 */
double uniform_to_sinusoidal_180(double uniform)
{
  return 180*asin(2.0*uniform-1.0)/M_PI + 90.0;
}

/*
  Converts a number from a uniform distribution (between 0 and 1) to a
  number from a sinusoidal distribution (between 0 and 90, peaking at
  90).
 */
double uniform_to_sinusoidal_90(double uniform)
{
  double value;
  value = uniform_to_sinusoidal_180(uniform);
  if(value > 90)
    value = 180 - value;
  return value;
}

//START REGION RELEASE
/* Add zeropad zero's before and after data and then add extra zero's
   to make the data length a power of two. If circularpad is set the
   padding is not done with zero's but with the edges of the data
   itself. When duplicate is set the profile is first duplicated,
   which can be usefull when correlating profiles which occur at the
   edge. Then it correlates the data and returns the lag number where
   the maximum was found. correl_max indices how many times the
   maximum of the correlation function was higher than the
   minimum. Returns 1 if succesful. verbose-1 is number of spaces
   before output.*/
int find_peak_correlation(float *data1, float *data2, int ndata, int zeropad, int circularpad, int duplicate, int *lag, float *correl_max, verbose_definition verbose)
{
  int i, lag_max;
  int npoints;
  float *paddata1, *paddata2, *ans, ans_max, ans_min;
  
  if(duplicate)
    duplicate = ndata;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("%d points in the data\n", ndata);
    if(zeropad != 0) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Padding at least %d points before and after data\n", zeropad);
    }
    if(duplicate) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Duplicating data to avoid wrap problems enabled.\n");
    }
  }
  i = (int) (log10(1.0 * (ndata+2*zeropad+duplicate))/0.301031);
  npoints = pow(2.0,(i+1));
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Going to zero-pad it to %d points\n", npoints);
  }

  paddata1 = (float *)malloc(npoints*sizeof(float));
  paddata2 = (float *)malloc(npoints*sizeof(float));
  ans = (float *)malloc(2*npoints*sizeof(float));
  if(paddata1 == NULL || paddata2 == NULL || ans == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR find_peak_correlation: Memory allocation error.");
    return 0;
  }

  zeropad = (npoints - duplicate - ndata)/2;
  for(i = 0; i < ndata; i++) {
    paddata1[i+zeropad] = data1[i];
    paddata2[i+zeropad] = data2[i];
    if(duplicate) {
      paddata1[i+ndata+zeropad] = data1[i];
      paddata2[i+ndata+zeropad] = data2[i];
    }
  }
  for(i = 0; i < zeropad; i++) {
    paddata1[i] = data1[i-zeropad+ndata];
    paddata2[i] = data2[i-zeropad+ndata];
  }
  for(i = ndata+duplicate+zeropad; i < npoints; i++) {
    paddata1[i] = data1[i-duplicate-zeropad-ndata];
    paddata2[i] = data2[i-duplicate-zeropad-ndata];
  }

  /*
  pgplot_viewport_def viewport;
  pgplot_clear_viewport_def(&viewport);
  sprintf(viewport.plotDevice, "2/xs");
  pgplot_box_def pgplotbox;
  clear_pgplot_box(&pgplotbox);
  sprintf(pgplotbox.title, "data1");
  pgplotbox.drawtitle = 1;
  pgplotGraph1(viewport, data1, NULL, NULL, ndata, 0, (float)ndata, 0, (float)ndata, 0, 0, 0, pgplotbox, 0, 1, 1, 1, NULL); 
  sprintf(viewport.plotDevice, "3/xs");
  sprintf(pgplotbox.title, "paddata1");
  pgplotGraph1(viewport, paddata1, NULL, NULL, ndata, 0, (float)ndata, 0, (float)ndata, 0, 0, 0, pgplotbox, 0, 1, 1, 1, NULL); 
  sprintf(viewport.plotDevice, "4/xs");
  sprintf(pgplotbox.title, "paddata2");
  pgplotGraph1(viewport, paddata2, NULL, NULL, ndata, 0, (float)ndata, 0, (float)ndata, 0, 0, 0, pgplotbox, 0, 1, 1, 1, NULL); 
  */
  
  /* NR code:   correl(paddata1-1, paddata2-1, npoints, ans-1); */
  if(crosscorrelation_fft(paddata1, paddata2, npoints, ans, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR find_peak_correlation: Cross correlation failed.");
    return 0;
  }

  lag_max = 0;
  ans_max = ans[0];
  ans_min = ans[0];
  for(i = 1; i < npoints; i++) {
    if(ans[i] > ans_max) {
      ans_max = ans[i];
      lag_max = i;
    }
    if(ans[i] < ans_min)
      ans_min = ans[i];
  }
  if(lag_max >= npoints/2)
    lag_max -= npoints;
  *correl_max = ans_max/ans_min;
  *lag = lag_max;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Found a lag of %d (correlation %f higher)\n", *lag, *correl_max);
  }
  /*
    for(i = ndata/2; i <= ndata-1; i++) 
      printf("%ld %f\n", i-ndata, ans[i]);
    for(i = 0; i < ndata/2; i++) 
      printf("%ld %f\n", i, ans[i]);
  */
  /*  pgplotGraph1("2/xs", 0, 0, ans, npoints, 0, npoints, 0, "bin", "I", "ans", NULL, 0, 0); */

  free(paddata1);
  free(paddata2);
  free(ans);
  return 1;
}
//START REGION DEVELOP



/*
  Interpolate the functions to get dataY(x). The value y is found via
  linear interpolation. dy is the estimated error. Order indicates
  what type of interpolation is done. 2 indicates a 2-point
  interpolation (linear interpolation).

  returns 1 on success, 0 if out of bounds.

 */
/*
int interpolate_NR(float *dataX, float *dataY, long nrdatapoints, float x, float *y, float *dy, long int order, int verbose)
{
  unsigned long idx;
  long idx_start;

  if(order < 2) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR interpolate: Cannot do less than a 2-point interpolation");
    return 0;
  }
  if(order > nrdatapoints) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR interpolate: Cannot do a %ld-point interpolation with %ld data points", order, nrdatapoints);
    return 0;
  }


  locate(dataX-1, nrdatapoints, x, &idx);
  if(idx == 0 || idx == nrdatapoints) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR interpolate: Out of bounds. Value %f appears to be less than %f or more than %f", x, dataX[0], dataX[nrdatapoints-1]);
    return 0;
  }
  idx--;
  if(verbose.verbose) 
    printf("interpolate: Value x=%f is in between %f and %f (index %ld and %ld)\n", x, dataX[idx], dataX[idx+1], idx, idx+1);

  idx_start = idx-(order-1)/2;
  if(idx_start < 0)
    idx_start = 0;
  if(idx_start + order > nrdatapoints)
    idx_start = nrdatapoints - order;
  if(verbose.verbose)
    printf("interpolate: take start index = %ld for the %ld-point interpolation function\n", idx_start, order);
  polint(&dataX[idx_start-1], &dataY[idx_start-1], order, x, y, dy);
  if(verbose.verbose)
    printf("interpolate: found y = %f with %f uncertainty\n", *y, *dy);

  return 1;
}
*/

/*
  Interpolate the functions to get dataY(x). The value y is found via
  linear interpolation. Order indicates what type of interpolation is
  done. 2 indicates a 2-point interpolation (linear
  interpolation). The data points are assumed to be stored in
  increasing order.

  returns 1 on success, 0 if out of bounds.

  interpolate() will call interpolate_double() 
*/
int interpolate_double(double *dataX, double *dataY, long nrdatapoints, double x, double *y, long int order, verbose_definition verbose)
{
  unsigned long idx;
  long idx_start;

  if(order < 2) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR interpolate: Cannot do less than a 2-point interpolation");
    return 0;
  }
  if(order > nrdatapoints) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR interpolate: Cannot do a %ld-point interpolation with %ld data points", order, nrdatapoints);
    return 0;
  }

  if(x < dataX[0] || x > dataX[nrdatapoints-1]) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR interpolate: Out of bounds. Value %lf appears to be less than %lf or more than %lf", x, dataX[0], dataX[nrdatapoints-1]);
    return 0;
  }

  idx = gsl_interp_bsearch (dataX, x, 0, nrdatapoints);

  //  locate(dataX-1, nrdatapoints, x, &idx);
  if(verbose.verbose) 
    printf("interpolate: Value x=%lf is in between %lf and %lf (index %ld and %ld)\n", x, dataX[idx], dataX[idx+1], idx, idx+1);

  idx_start = idx-(order-1)/2;
  if(idx_start < 0)
    idx_start = 0;
  if(idx_start + order > nrdatapoints)
    idx_start = nrdatapoints - order;
  if(verbose.verbose)
    printf("interpolate: take start index = %ld for the %ld-point interpolation function\n", idx_start, order);

  gsl_interp *interp;
  gsl_interp_accel *interp_accel;

  if(order == 2) {
    interp = gsl_interp_alloc(gsl_interp_linear, order);
  }else {
#if GSL_VERSION_NUMBER >= 101
    interp = gsl_interp_alloc(gsl_interp_polynomial, order);
#else
    printerror(verbose.debug, "ERROR interpolate: Interpolation with an order > 2 is not supported in GSL 1.0");
    exit(0);
#endif
  }

  if(gsl_interp_init(interp, &dataX[idx_start], &dataY[idx_start], order) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR interpolate: Error in gsl_interp_init");
    return 0;
  }


  interp_accel = gsl_interp_accel_alloc();

  if(gsl_interp_eval_e (interp, &dataX[idx_start], &dataY[idx_start], x, interp_accel, y) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR interpolate: Error in gsl_interp_eval_e");
    return 0;
  }

  if(verbose.verbose)
    printf("interpolate: found y = %e\n", *y);

  gsl_interp_accel_free(interp_accel);
  gsl_interp_free(interp);

  return 1;
}

/* This wrapper function will call interpolate_double. */
int interpolate(float *dataX, float *dataY, long nrdatapoints, float x, float *y, long int order, verbose_definition verbose)
{
  double *doubleX, *doubleY, ydouble;
  int ret;
  long i;

  doubleX = (double *)malloc(nrdatapoints*sizeof(double));
  doubleY = (double *)malloc(nrdatapoints*sizeof(double));
  if(doubleX == NULL || doubleY == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR interpolate: Cannot allocate memory");
    return 0;
  }
  for(i = 0; i < nrdatapoints; i++) {
    doubleX[i] = dataX[i];
    doubleY[i] = dataY[i];
  }
  ret = interpolate_double(doubleX, doubleY, nrdatapoints, x, &ydouble, order, verbose);
  *y = ydouble;
  free(doubleX);
  free(doubleY);
  return ret;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Calculates in which bin number the value x belongs, given the bin
   width dx and the minimum value of the distribution min_x. Set
   centered_at_zero to one to make the centre of one bin equal to
   zero. exta_phase will shift the bin_centre by a certain fraction of
   a bin. */
long calculate_bin_number(double x, double dx, double min_x, int centered_at_zero, double extra_phase)
{
  long bin, step, binzero;
  if(min_x < 0)
    step = -min_x/dx+10;
  else
    step = 0;
  bin     = (    x+(step+0.5*centered_at_zero-extra_phase)*dx)/dx;
  binzero = (min_x+(step+0.5*centered_at_zero-extra_phase)*dx)/dx;
  bin -= binzero;
  return bin;
}

//START REGION DEVELOP
//START REGION RELEASE

// Like calculate_bin_number(), but returns to central value of bin binnr
double calculate_bin_location(long binnr, double dx, double min_x, int centered_at_zero, double extra_phase)
{
  long step, binzero;
  double x;
  if(min_x < 0)
    step = -min_x/dx+10;
  else
    step = 0;
  binzero = (min_x+(step+0.5*centered_at_zero-extra_phase)*dx)/dx;

  binnr = binnr + binzero;
  x = binnr*dx - (step+0.5*centered_at_zero-extra_phase)*dx;
  x += 0.5*dx; // From left edge to centre
  return x;
}

//START REGION DEVELOP
//START REGION RELEASE

// Like calculate_bin_number(), but returns to bin width dx required to make x fall in the centre of bin binnr. Can be useful to determine required binwidth to end up with a given number of bins.
double calculate_required_bin_width(double x, long binnr, double min_x, int centered_at_zero, double extra_phase, verbose_definition verbose)
{
  long lastbin, ok, timesinloop;
  double dx, offset;
  /*
    long step;
  if(min_x < 0)
    step = -min_x/dx+10;
  else
    step = 0;
  */
  //  x = (binnr + floor((min_x+(step+0.5*centered_at_zero-extra_phase)*dx)/dx))*dx - (step+0.5*centered_at_zero-extra_phase)*dx + 0.5*dx;
  // Ignored the floor, which makes the guess not accurate enough.
  //  x = (binnr + (min_x+(step+0.5*centered_at_zero-extra_phase)*dx)/dx)*dx - (step+0.5*centered_at_zero-extra_phase)*dx + 0.5*dx;
  //  x = binnr*dx + (min_x/dx+(step+0.5*centered_at_zero-extra_phase))*dx - (step+0.5*centered_at_zero-extra_phase)*dx + 0.5*dx;
  //  x = binnr*dx + min_x+(step+0.5*centered_at_zero-extra_phase)*dx - (step+0.5*centered_at_zero-extra_phase)*dx + 0.5*dx;
  //  x = binnr*dx + min_x + 0.5*dx;
  //  dx = (x - min_x)/(double)(binnr+0.5)
  // + 0.5 is requisted phase, i.e. at middle of the bin. I think you can make it larger to make fall later in bin, but is failing anyways because of floor thing.
  dx = (x - min_x)/(double)(binnr+0.5);

  timesinloop = 0;
  offset = 0;
  do {
    ok = 1;
    dx = (x - min_x)/(double)(binnr+0.5) + offset;
    lastbin = calculate_bin_number(x, dx, min_x, centered_at_zero, extra_phase);
    //    fprintf(stderr, "dx = %e; binnr = %ld\n", dx, lastbin);
    if(lastbin < binnr) {
      offset += (x-min_x-0.5*dx)/(double)(binnr+0.5) - dx;
      ok = 0;
    }else if(lastbin > binnr) {
      offset += (x-min_x+0.5*dx)/(double)(binnr+0.5) -dx;
      ok = 0;
    }
    timesinloop++;
    if(timesinloop > 10) {
      printerror(verbose.debug, "ERROR calculate_required_bin_width: Cannot find suitable binsize.\n");
      return dx;
    }
  }while(ok == 0);
  return dx;
}

//START REGION DEVELOP
//START REGION RELEASE

/*
  Given data in the range [min_x_data, max_x_data], determine the
  required binning. If rangex_set is non-zero, the range is adjusted
  to [rangex_min, rangex_max]. nrbins is the requested number of bins,
  which is used when nrbins_specified is set. Otherwise *dx is
  assumed. If centered_at_zero is set, the mid-point of one bin
  (possibly outside range) will coincide with zero. Otherwise an edge
  will coincide with zero. extra_phase changes where zero falls with
  respect to the bin-centre (cyclic behaviour with period of 1).

  The return parameters are the actual range that will be covered
  [min_x, max_x], which might be different from what is specified if
  randex_set is set (to avoid rounding errors). Otherwise, it will be
  set to [min_x_data, max_x_data]. The bin-size is returned in dx.
  
  Return value:
  0: ok
  1: Determination of bin width failed.
  2: Fist bin not in correct location
 */
int set_binning_histogram(double min_x_data, double max_x_data, int rangex_set, double rangex_min, double rangex_max, int nrbins_specified, long nrbins, int centered_at_zero, double extra_phase, double *min_x, double *max_x, double *dx, verbose_definition verbose)
{
  int reset;
  do { // Loop over the determination of the proper range to be use, as the -rangex and -rangey options might need to be adjusted when using the -range options if they fall at border between bins to avoid rounding errors.
    reset = 0;
    // The range is initially set to the range as found in the data
    *min_x = min_x_data;
    *max_x = max_x_data;
      
    if(rangex_set) {
      *min_x = rangex_min;
      *max_x = rangex_max;
    }
      
    // Guess binsizes, will need to be refined depending on -zero option etc
    if(nrbins_specified) {
      *dx = calculate_required_bin_width(*max_x, nrbins-1, *min_x, centered_at_zero, extra_phase, verbose);
      //START REGION DEVELOP
      long lastbin;
      lastbin = calculate_bin_number(*max_x, *dx, *min_x, centered_at_zero, extra_phase);
      if(lastbin != nrbins - 1) {
	printerror(verbose.debug, "ERROR set_binning_histogram: determination of bin width failed.\n");
	return 1;
      }
      //START REGION RELEASE
      if(verbose.verbose)
	fprintf(stdout, "Going to use binsize %e.\n", *dx);
    }
	
    // Check if minimum value really falls in bin zero
    {
      long firstbin;
      firstbin = calculate_bin_number(*min_x, *dx, *min_x, centered_at_zero, extra_phase);
      if(firstbin != 0) {
	printerror(verbose.debug, "ERROR set_binning_histogram: Expected first bin to be number zero (it is %ld).\n", firstbin);
	return 2;
      }
    }
	
    // Check if -rangex borders fall at edge of bin
    long i;
    double diff;
    i = calculate_bin_number(*max_x, *dx, *min_x, centered_at_zero, extra_phase);
    // diff is the offset of the max value from the bin centre
    diff = *max_x - calculate_bin_location(i, *dx, *min_x, centered_at_zero, extra_phase);
    //Express as a phase
    diff = diff/(*dx);
    if(verbose.debug) {
      printf("Current set max value (%e) is falling at phase=%e w.r.t. centre of last generated bin.\n", *max_x, diff);
    }
    if(diff > -0.501 && diff < -0.499) { // Maximum falls between last and one-but-last bin, so no danger of overflow, essentially the last bin is (almost) not used, so can make max half a bin smaller.
      if(rangex_set) {
	rangex_max -= 0.5*(*dx);                   // To avoid rounding errors, adjust range, which might be especially useful with the -trunc option
	printwarning(verbose.debug, "WARNING set_binning_histogram: Adjusting maximum value of the specified range to %e to avoid rounding errors. Going to reset choosen binning.", rangex_max);
	reset = 1;
      }
    }
    if(diff > 0.499 && diff < 0.501) { // Maximum falls at outer edge of last bin, so danger of overflow. Increase the max by half a bin, so we're sure enough bins are generated to accomodate all points.
      if(rangex_set) {
	rangex_max += 0.5*(*dx);                   // To avoid rounding errors, adjust range, which might be especially useful with the -trunc option
	printwarning(verbose.debug, "WARNING set_binning_histogram: Adjusting maximum value of the specified range to %e to avoid rounding errors. Going to reset choosen binning.", rangex_max);
	reset = 1;
      }else {
	*max_x += 0.5*(*dx);
	if(nrbins_specified) {
	  printwarning(verbose.debug, "WARNING set_binning_histogram: Nr of bins might be different from what was requested to avoid rounding errors.");
	}
      }
    }
    
    i = calculate_bin_number(*min_x, *dx, *min_x, centered_at_zero, extra_phase);
    // diff is the offset of the min value from the bin centre
    diff = *min_x - calculate_bin_location(i, *dx, *min_x, centered_at_zero, extra_phase);
    //Express as a phase
    diff = diff/(*dx);
    if(verbose.debug) {
      printf("Current set min value (%e) is falling at phase=%e w.r.t. centre of first generated bin.\n", *min_x, diff);
    }
    if(diff > -0.501 && diff < -0.499) { // Minimum falls at outer edge for first bin, so danger of overflow
      if(rangex_set) {
	rangex_min -= 0.5*(*dx);                   // To avoid rounding errors, adjust range, which might be especially useful with the -trunc option
	printwarning(verbose.debug, "WARNING set_binning_histogram: Adjusting minimum value of the specified range to %e to avoid rounding errors. Going to reset choosen binning.", rangex_min);
	reset = 1;
      }else {
	*min_x -= 0.5*(*dx);
	if(nrbins_specified) {
	  printwarning(verbose.debug, "WARNING set_binning_histogram: Nr of bins might be different from what was requested to avoid rounding errors.");
	}
      }
    }
    if(diff > 0.499 && diff < 0.501) { // Maximum falls between first and second bin, so no danger of overflow
      if(rangex_set) {
	rangex_min += 0.5*(*dx);                   // To avoid rounding errors, adjust range, which might be especially useful with the -trunc option
	printwarning(verbose.debug, "WARNING set_binning_histogram: Adjusting minimum value of the specified range to %e to avoid rounding errors. Going to reset choosen binning.", rangex_min);
	reset = 1;
      }
    }
    if(reset) {
      printwarning(verbose.debug, "WARNING set_binning_histogram: Re-adjusting choosen binning.\n", rangex_max);
      //	  free(data_x);
      //	  if(twoDmode) {
      //	    free(data_y);
      //	  }
    }
  }while(reset == 1);
  return 0;
}

//START REGION DEVELOP

void swap_double_array_elements(double *array, long index1, long index2)
{
  double tmp_double;
  tmp_double = array[index1];
  array[index1] = array[index2];
  array[index2] = tmp_double;
}


/*
  Fits cubic B-spline functions to get dataY(x), and dataYerr defines
  the errorbar on the data points (set to NULL if not used). The value
  y is found fitting cubic B-spline basis functions with uniform
  breakpoints. The number of fit coefficients is set by
  nrcoefficients. Set y=NULL and/or yerr=NULL if you're not interested
  in the error on the y estimate. The fit is returned in splineFit
  (unless it is set to NULL, memory will be allocated) for future use,
  such as with cubicBspline_eval(). Use cubicBspline_free() to free up
  this memory. dataX is not required to be sorted. If nrtrials is set,
  instead of assuming uniform gridding of the nodes (which is first
  trial), random distributions of node positions are explored
  (nrtrials nr of times) to find a better fit. The variable
  frac_replace defines the fraction of grid-points which are replaced
  by random values in each itteration. verbose-1 is number of spaces
  before output. If debug is set, more information is printed.

  returns 1 on success, 0 on error.

  cubicBspline_fit() will call cubicBspline_fit_double() 
*/
int cubicBspline_fit_double(double *dataX, double *dataY, double *dataYerr, long nrdatapoints, int nrcoefficients, double x, double *y, double *yerr, splineFit_def *splineFit, long nrtrials, double frac_replace, verbose_definition verbose)
{
#if GSL_VERSION_NUMBER < 109
  printerror(verbose.debug, "ERROR cubicBspline_fit_double: Not supported for GSL < 1.9");  
  exit(0);
#else

  size_t i, j;
  int order, nbreak;
  gsl_bspline_workspace *spline_workspace;
  gsl_vector *spline;
  gsl_vector *xarray, *yarray, *warray;
  gsl_matrix *X, *cov;
  double xmin, xmax, xvalue;
  long trial, idnum;
  double chisq, chisq_best, dof, yerror, yvalue, *breakpoints_array, *breakpoints_array_best;
  gsl_vector *c, *breakpoints;
  gsl_multifit_linear_workspace *multifit_workspace;
  gsl_rng *rand_num_gen;
  const gsl_rng_type *rand_num_gen_type;
  //  double Rsq, tss;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Smooth %ld datapoints with cubic B-spline basis functions with %d coefficients\n", nrdatapoints, nrcoefficients);
  }

  order = 4;  // for cubic B-spline basis functions

  //Determine minimum and maximum value of dataX
  xmin = dataX[0];
  xmax = dataX[0];
  for(i = 0; i < nrdatapoints; i++) {
    if(dataX[i] < xmin)
      xmin = dataX[i];
    if(dataX[i] > xmin)
      xmax = dataX[i];
  }

  // If x-value is used, check if in range
  if(y != NULL || yerr != NULL) {
    if(x < xmin || x > xmax) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR cubicBspline_fit_double: Out of bounds. Value %lf appears to be less than %lf or more than %lf", x, xmin, xmax);
      return 0;
    }
  }

  // According to gsl manual
  // nrcoefficients = nbreak + order - 2
  // nbreak = nrcoefficients + 2 - order 
  nbreak = nrcoefficients + 2 - order;
  if(nbreak < 2) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cubicBspline_fit_double: Specify more coefficients in order fit a cubic spline.");
    return 0;
  }
  if(nbreak <= 2 && nrtrials) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cubicBspline_fit_double: Specify more coefficients in order to allow non-random gridding. The grid-points are currently fixed by range in data.");
    return 0;
  }
  breakpoints = gsl_vector_alloc(nbreak);
  breakpoints_array = (double *)malloc(nbreak*sizeof(double));
  breakpoints_array_best = (double *)malloc(nbreak*sizeof(double));
  if(breakpoints_array == NULL || breakpoints_array_best == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cubicBspline_fit_double: Memory allocation error");
    return 0;
  }

  // allocate workspace for spline
  spline_workspace = gsl_bspline_alloc(order, nbreak);

  spline = gsl_vector_alloc(nrcoefficients);
  xarray = gsl_vector_alloc(nrdatapoints);
  yarray = gsl_vector_alloc(nrdatapoints);
  warray = gsl_vector_alloc(nrdatapoints);
  X = gsl_matrix_alloc(nrdatapoints, nrcoefficients);

  // Prepare the data to be fitted
  for(i = 0; i < nrdatapoints; i++) {
    gsl_vector_set(xarray, i, dataX[i]);
    gsl_vector_set(yarray, i, dataY[i]);
    if(dataYerr != NULL)
      gsl_vector_set(warray, i, 1.0 / (dataYerr[i] * dataYerr[i]));
    else
      gsl_vector_set(warray, i, 1.0);
  }

  // Set-up random number generator if required
  if(nrtrials) {
    gsl_rng_env_setup();    /* Set the default generators, can be influenced by environment variables */
    rand_num_gen_type = gsl_rng_default;
    rand_num_gen = gsl_rng_alloc(rand_num_gen_type);
    randomize_idnum(&idnum);
    gsl_rng_set(rand_num_gen, idnum);
  }

  for(trial = 0; trial <= nrtrials; trial++) {
    if(trial == 0) { // First trial: use uniform breakpoints on the provided range
      // The following is identical to: gsl_bspline_knots_uniform(xmin, xmax, spline_workspace);
      for (i = 0; i < nbreak; i++) {
	xvalue = xmin + i*(xmax-xmin)/(double)(nbreak-1);
	breakpoints_array[i] = xvalue;
      }
    }else {
      // Start with current best grid
      for (i = 0; i < nbreak; i++) {
	breakpoints_array[i] = breakpoints_array_best[i];
      }
      // Want to keep end points fixed, put them in locations 0 and 1. At this point, the array should be sorted
      swap_double_array_elements(breakpoints_array, 1, nbreak-1);
      long nrswap, iswap;
      nrswap = (1.0-frac_replace)*(nbreak-2);
      for(iswap = 0; iswap < nrswap; iswap++) {
	i = gsl_rng_uniform_int(rand_num_gen, nbreak-2);
	swap_double_array_elements(breakpoints_array, 2+iswap, 2+i);
      }
      // Assign random grid points for others
      for (i = 2+nrswap; i < nbreak; i++) {
	breakpoints_array[i] = xmin + gsl_rng_uniform(rand_num_gen)*(xmax-xmin);
      }
      gsl_sort(breakpoints_array, 1, nbreak);   // Breakpoints need to be sorted
    }

    for (i = 0; i < nbreak; i++) {
      gsl_vector_set(breakpoints, i, breakpoints_array[i]);
    }
    gsl_bspline_knots(breakpoints, spline_workspace);

    if(verbose.debug) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Trying breakpoints:");
      for (i = 0; i < nbreak; i++) {
	printf(" %f", breakpoints_array[i]);
      }
      printf("\n");
    }

    // construct the fit matrix X 
    for (i = 0; i < nrdatapoints; ++i) {
      double xi = gsl_vector_get(xarray, i);
      // compute B_j(xi) for all j
      gsl_bspline_eval(xi, spline, spline_workspace);
      // fill in row i of X
      for (j = 0; j < nrcoefficients; j++) {
	double Bj = gsl_vector_get(spline, j);
	gsl_matrix_set(X, i, j, Bj);
      }
    }
        

    c = gsl_vector_alloc(nrcoefficients);
    cov = gsl_matrix_alloc(nrcoefficients, nrcoefficients);
    multifit_workspace = gsl_multifit_linear_alloc(nrdatapoints, nrcoefficients);
    // do the fit
    gsl_multifit_wlinear(X, warray, yarray, c, cov, &chisq, multifit_workspace);
    
    dof = nrdatapoints - nrcoefficients;
    //  tss = gsl_stats_wtss(warray->data, 1, yarray->data, 1, yarray->size);
    //  Rsq = 1.0 - chisq / tss;
    
    if((verbose.verbose && trial == 0) || verbose.debug) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Fit has a reduced chisquare of %lf\n", chisq / dof);
      if(nrtrials > 0) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  Starting to exploring %ld trial griddings\n", nrtrials);
      }
    }

    // Store the fit, if necessary, for future use.
    if(trial == 0 || chisq < chisq_best) {
      if(splineFit != NULL) {
	splineFit->order = order;
	splineFit->nbreak = nbreak;
	splineFit->x_start = xmin;
	splineFit->x_end = xmax;
	if(trial == 0) {
	  splineFit->coefficients = (double *)malloc(nrcoefficients*sizeof(double));
	  splineFit->breakpoints = (double *)malloc(nbreak*sizeof(double));
	  if(splineFit->coefficients == NULL || splineFit->breakpoints == NULL) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR cubicBspline_fit_double: Memory allocation error");
	    return 0;
	  }
	}
	for(i = 0; i < nrcoefficients; i++) {
	  splineFit->coefficients[i] = gsl_vector_get(c, i);
	}
	for(i = 0; i < nbreak; i++) {
	  splineFit->breakpoints[i] = breakpoints_array[i];
	}
      }
      // If x-value is used, obtain the fit value + error
      if(y != NULL || yerr != NULL) {
	gsl_bspline_eval(x, spline, spline_workspace);
	gsl_multifit_linear_est(spline, c, cov, &yvalue, &yerror);
	if(y != NULL)
	  *y = yvalue;
	if(yerr != NULL)
	  *yerr = yerror;
      }
      for(i = 0; i < nbreak; i++) {
	breakpoints_array_best[i] = breakpoints_array[i];
      }
      chisq_best = chisq;
    }
  }  // End of trial loop

  if(nrtrials && verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  After trying non-uniform gridding a reduced chisquare of %lf was found\n", chisq_best / dof);
  }

  gsl_bspline_free(spline_workspace);
  gsl_vector_free(spline);
  gsl_vector_free(xarray);
  gsl_vector_free(yarray);
  gsl_vector_free(warray);
  gsl_matrix_free(X);
  gsl_matrix_free(cov);
  gsl_vector_free(c);
  gsl_vector_free(breakpoints);
  gsl_multifit_linear_free(multifit_workspace);
  if(nrtrials)
    gsl_rng_free(rand_num_gen);

#endif
  return 1;
}

/* This wrapper function will call cubicBspline_fit_double. */
int cubicBspline_fit(float *dataX, float *dataY, float *dataYerr, long nrdatapoints, int nrcoefficients, float x, float *y, float *yerr, splineFit_def *splineFit, long nrtrials, float frac_replace, verbose_definition verbose)
{
  double *doubleX, *doubleY, *doubleYerr, ydouble, yerrdouble;
  int ret;
  long i;

  doubleX = (double *)malloc(nrdatapoints*sizeof(double));
  doubleY = (double *)malloc(nrdatapoints*sizeof(double));
  if(dataYerr != NULL) {
    doubleYerr = (double *)malloc(nrdatapoints*sizeof(double));
    if(doubleYerr == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR interpolate: Memory allocation error");
      return 0;
    }
  }else {
    doubleYerr = NULL;
  }
  if(doubleX == NULL || doubleY == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR interpolate: Memory allocation error");
    return 0;
  }
  for(i = 0; i < nrdatapoints; i++) {
    doubleX[i] = dataX[i];
    doubleY[i] = dataY[i];
    if(dataYerr != NULL)
      doubleYerr[i] = dataYerr[i];
  }
  ret = cubicBspline_fit_double(doubleX, doubleY, doubleYerr, nrdatapoints, nrcoefficients, x, &ydouble, &yerrdouble, splineFit, nrtrials, frac_replace, verbose);
  *y = ydouble;
  if(yerr != NULL)
    *yerr = yerrdouble;
  free(doubleX);
  free(doubleY);
  if(dataYerr != NULL)
    free(doubleYerr);

  return ret;
}

// Free up allocated space after spline was allocated with cubicBspline_fit_double
void cubicBspline_free(splineFit_def *splineFit)
{
  splineFit->order = 0;
  splineFit->nbreak = 0;
  splineFit->x_start = 0;
  splineFit->x_end = 0;
  free(splineFit->coefficients);
  splineFit->coefficients = NULL;
  free(splineFit->breakpoints);
  splineFit->breakpoints = NULL;
}

/*
  Using a fit made by cubicBspline_fit(), calculate the value of the
  cubic B-spline functions at x.

  returns 1 on success, 0 if out of bounds.

  cubicBspline_fit() will call cubicBspline_fit_double() 
*/
int cubicBspline_eval_double(splineFit_def splineFit, double x, double *y, verbose_definition verbose)
{
#if GSL_VERSION_NUMBER < 109
  printerror(verbose.debug, "ERROR cubicBspline_eval_double: Not supported for GSL < 1.9");  
  exit(0);
#else
  gsl_bspline_workspace *spline_workspace;
  gsl_vector *spline;
  int nrcoefficients, i;
  gsl_vector *breakpoints;

  // If x-value is used, check if in range
  if(x < splineFit.x_start || x > splineFit.x_end) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cubicBspline_eval_double: Out of bounds. Value %lf appears to be less than %lf or more than %lf", x, splineFit.x_start, splineFit.x_end);
    return 0;
  }

  // allocate workspace for spline
  spline_workspace = gsl_bspline_alloc(splineFit.order, splineFit.nbreak);

  nrcoefficients = splineFit.nbreak -2 + splineFit.order;
  spline = gsl_vector_alloc(nrcoefficients);

  // Assign breakpoints
  // Assume uniform breakpoints were used on the provided range */
  //  gsl_bspline_knots_uniform(splineFit.x_start, splineFit.x_end, spline_workspace);
  breakpoints = gsl_vector_alloc(splineFit.nbreak);
  for (i = 0; i < splineFit.nbreak; i++) {
    gsl_vector_set(breakpoints, i, splineFit.breakpoints[i]);
  }
  gsl_bspline_knots(breakpoints, spline_workspace);


  // Get the basis functions
  gsl_bspline_eval(x, spline, spline_workspace);

  // To obtain the actual value, multiply basis function value with fit coefficient
  *y = 0;
  for(i = 0; i < nrcoefficients; i++) {
    *y += gsl_vector_get(spline, i) * splineFit.coefficients[i];
  }

  gsl_bspline_free(spline_workspace);
  gsl_vector_free(spline);
  gsl_vector_free(breakpoints);
 
  return 1;
#endif
}

//  Wrapper for cubicBspline_fit_double() 
int cubicBspline_eval(splineFit_def splineFit, double x, float *y, verbose_definition verbose)
{
  int ret;
  double yvalue;
  ret = cubicBspline_eval_double(splineFit, x, &yvalue, verbose);
  *y = yvalue;
  return ret;
}

//START REGION DEVELOP
//START REGION RELEASE

double kstest_cdf_flat(double x, double min_x, double max_x)
{
  if(x <= min_x)
    return 0;
  if(x >= max_x)
    return 1;  
  return (x-min_x)/(max_x-min_x);
}

double kstest_cdf_sin(double x)
{
  if(x <= 0)
    return 0;
  if(x >= 90)
    return 1;
  return 1-cos(x*M_PI/180.0);
}

//START REGION DEVELOP
//START REGION RELEASE

/*
  If data2 != NULL, n2 > 0: Given two data-sets, calculate the maximum
  difference between their cummulative distribution, and the
  corresponding probability that the distributions are not sampled
  from the same underlying distribution. A small probability means the
  distributions are statistically shown to be different. The input
  arrays data1 and data2 will be sorted after the call to kstest.

  If data2 == NULL, n2 > 0:
    - If cdf_type == 0: a pointer to the function cdf should be
      provided, which is the cummulative distribution of the distribution
      which is assumed to be the parent distribution from which the points
      data1 are drawn.
    - If cdf_type == 1: Instead use a flat distribution starting/ending 
      at maximum value of data1.
    - If cdf_type == 2: Instead use a sin distribution from 0 to 90 deg.
    - If cdf_type == 3: Instead use a flat distribution starting/ending 
      at values input_value1 and input_value2.

  This implementation should be similar to that implemented in
  Numerical Recipes in C, second edition.

  */
void kstest(double *data1, long n1, double *data2, long n2, int cdf_type, double input_value1, double input_value2, double (*cdf)(double), double *max_diff, double *prob, verbose_definition verbose)
{
  long i1, i2;
  double effective_n, ks_statistic, sign, cur_term, last_term, coeff;
  int converged;

  gsl_sort(data1, 1, n1);   // Sort the input array to make the cummulative distribution
  if(n2 > 0 && data2 != NULL)
    gsl_sort(data2, 1, n2);   // Sort the other input array to make the cummulative distribution

  // Find the largest separation between the cdf's, which is the test statistic
  *max_diff = 0;
  if(n2 > 0 && data2 != NULL) {
    i1 = 0;
    i2 = 0;
    while(i1 < n1 && i2 < n2) { // Step through the cummulative distributions
      double diff;
      if(data1[i1] == data2[i2]) {   // Overlapping points, make step in both cdf's
	i1++;
	i2++;
      }else if(data1[i1] < data2[i2]) { // If bin cdf of first distribution starts before that of second distribution
	i1++;                           // Make step in the first distribution to catch up
      }else {                           // If bin cdf of second distribution starts before that of first distribution
	i2++;                           // Make step in the first distribution to catch up
      }
      diff = fabs(i1/(double)n1 - i2/(double)n2); // Difference in cdf of current points. The x-value of one distribution is always in between those of the current point of the other distribution and the next.
      if(diff > *max_diff)          // See if difference in cdf at current point is larger than previous largest difference
	*max_diff = diff;
    }
    effective_n=n1*n2/(double)(n1+n2);  // weighting of statistic depends on the effective number of points
  }else { // Simpler case: the cdf function to compare with is a continuous defined function
    double cdf_right, cdf_left, cdf_model;
    cdf_left = 0;               // cdf histogram starts at zero
    for(i1 = 0; i1 < n1; i1++) { // Just loop over input array = steps in cdf
      cdf_right = (i1+1)/(double)n1;  // The value of the cdf at right-hand side of bin (stepping up)
      if(cdf_type == 0) {
	cdf_model = (*cdf)(data1[i1]); // The model value at current point
      }else if(cdf_type == 1) {
	cdf_model =  kstest_cdf_flat(data1[i1], data1[0], data1[n1-1]);
	if(i1 == 0) {
	  printwarning(verbose.debug, "WARNING kstest: Probability will be overestimated, since the minimum/maximum of the uniform distribution is based on the input values.");    
	}
      }else if(cdf_type == 3) {
	cdf_model =  kstest_cdf_flat(data1[i1], input_value1, input_value2);
      }else if(cdf_type == 2) {
	cdf_model =  kstest_cdf_sin(data1[i1]);
      }else {
	fflush(stdout);
	printerror(verbose.debug, "ERROR kstest: Undefined type of cdf is specified.");
	exit(0);
      }
      if(fabs(cdf_right-cdf_model) > *max_diff)   // If model difference is after before step up
	*max_diff = fabs(cdf_right-cdf_model);
      if(fabs(cdf_left-cdf_model) > *max_diff)    // If model difference is larger before step up
	*max_diff = fabs(cdf_left-cdf_model);
      cdf_left = cdf_right;   // Value cdf at right-hand side bin is that of the left-hand side of next bin 
    }
    effective_n=n1;  // weighting of statistic depends on just the number of points in the distribution
  }

  if(effective_n < 4) {
    printwarning(verbose.debug, "WARNING kstest: Number of data-points is too low to make use of approximations used in this implementation of the KS-test.");    
  }

  effective_n=sqrt(effective_n);  // or actually the square-root of the effective number of points
  ks_statistic = (*max_diff)*(effective_n+0.12+0.11/effective_n); // Approximate interpolation to give right behaviour as function of effective_n. In principle need a look-up table.
  coeff = -2.0*ks_statistic*ks_statistic;  // The coefficients in the exponentials in the sum we're going to calculate below.

  *prob = 0;     // The probability is the sum of a series
  sign = 1;      // The sign of the terms in the sum will be alternating, first term is positive.
  last_term = 0; // Last term undefined, but set to zero to avoid exiting the loop prematurely.
  converged = 0; // Will be set if series is converged
  for(i1 = 1; i1 <= 100; i1++) {  // Add at most 100 terms together, otherwise the probability is very high, possibly even 1.
    cur_term = sign*2.0*exp(coeff*i1*i1); // The answer is an series with alternating sign
    *prob += cur_term;               // Add the next term to the series
    // Quit series if the current term becomes very small compared to previous term,
    // or if the current term is small in comparison to the everal sum. If that is the
    // case the series is practically converged. Need this condition, because for identical
    // distributions the series does not converge properly.
    if(fabs(cur_term) <= 1e-5*fabs(last_term) || fabs(cur_term) <= 1e-10*(*prob)) {
      converged = 1;
      break;
    }
    last_term = cur_term; // Remember last term to use in test convergence
    sign = -sign;         // Sign of terms is alternating
  }
  if(!converged)
    *prob = 1; // Didn't converge, which means the probability should be very large, possibly even 1.

  //  *prob = 1-0.997300203937;
  //  *prob = 1-0.682689492137;
  //  *prob = 0;
  //  *prob = 1;

  if(verbose.verbose) {
    printf("KS-test statistic max_diff: %lf = %e\n", *max_diff, *max_diff);
    printf("KS-test probability:        %lf = %e\nA small probability means the two sets of points are drawn from a different distribution.\n", *prob, *prob);
    printf("The null hypothesis ");
    if(n2 > 0 && data2 != NULL) {
      printf("(data sets are drawn from the same distribution) ");
    }else {
      printf("(data set is drawn from the specified distribution) ");
    }
#if GSL_VERSION_NUMBER >= 104
    printf("can be rejected at the %.2lf sigma level.\n", gsl_cdf_gaussian_Pinv(0.5*(1.0-*prob)+0.5, 1.0));
#else
    printf("can be rejected at the XXXX sigma level (need GSL >= 1.4 to get this number).\n");
#endif
  }
  //  printf("XXXX %lf\n", erf(1.0/sqrt(2)));
}

//START REGION DEVELOP
