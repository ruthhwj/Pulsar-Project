//START REGION RELEASE
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include "psrsalsa.h"

//START REGION DEVELOP
#include "gsl/gsl_roots.h"
#include "gsl/gsl_errno.h"

/* Based on Russells arbitrary_bin_shift function
   (bin_shift.cc). Returns 0 if failed, otherwise 1. */

//START REGION RELEASE

void print_fftw_version_used(FILE *stream)
{
  fprintf(stream, "%s (library)", fftwf_version);
}

int rotateSinglepulse(float *data, int npts, float epsilon, verbose_definition verbose)
{
  int i, npts2;
  float fac, dtheta;
  fftwf_complex *dataFFT;
  fftwf_plan plan1, plan2;

  npts2 = npts/2+1;
  dataFFT = (fftwf_complex *)fftwf_malloc(npts2*sizeof(fftwf_complex));
  if(dataFFT == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR rotateSinglepulse: fftwf_malloc failed.");
    return 0;
  }
  plan1 = fftwf_plan_dft_r2c_1d(npts, data, dataFFT, FFTW_ESTIMATE);
  plan2 = fftwf_plan_dft_c2r_1d(npts, dataFFT, data, FFTW_ESTIMATE);

  /* Do the forward FFT */
  fftwf_execute(plan1);

  /* multiply by exponential... The scale factor is to correct for the
     DFT normalization */
  fac = 1.0/(float)npts;
  dtheta = -2.0*M_PI*epsilon/(float)npts;
  for (i=0; i < npts2; i++) {
    dataFFT[i] *= fac*(cos(i*dtheta) + I*sin(i*dtheta));
  }

  /* Do the inverse FFT */
  fftwf_execute(plan2);

  fftwf_destroy_plan(plan1);
  fftwf_destroy_plan(plan2);
  fftwf_free(dataFFT); 
  return 1;
}

//START REGION DEVELOP


/* Smooth each pulse independent of each other using fft algorithm by
   zapping nrSmoothFreq frequency bins. If highpass is set, the low
   frequencies are zapped rathern than the low frequencies.  Returns 0
   if failed, otherwise 1. Verbose determines the number of spaces
   before output. */
int fftSmooth(float *data, long npts, int highpass, long nrSmoothFreq, verbose_definition verbose)
{
  int i;
  long npts2, freqlow;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
  }
  if(highpass) {
    if(verbose.verbose) {
      printf("Smooth data by zapping frequency channels 0 - %ld\n", nrSmoothFreq);
    }
    return fftZap(data, npts, 0, nrSmoothFreq, verbose);
  }else {
    npts2 = npts/2 + 1;
    freqlow = npts2-nrSmoothFreq;
    if(freqlow < 0)
      freqlow = 0;
    if(verbose.verbose) {
      printf("Smooth data by zapping frequency channels %ld - %ld\n", freqlow, npts2-1);
    }
    return fftZap(data, npts, freqlow, npts2-1, verbose);
  }
}

/* Zaps given range of frequency bins.  Returns 0 if failed, 
otherwise 1.
*/
int fftZap(float *data, long npts, long freq1, long freq2, verbose_definition verbose)
{
  long i, npts2;
  fftwf_complex *dataFFT;
  fftwf_plan plan1, plan2;

  if(sizeof(int) == 4) {
    if(npts > 2147483644) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR fftZap: requested fft too long.");
      return 0;
    }
  }else if(npts > 32764) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fftZap: requested fft too long.");
    return 0;
  }
  npts2 = npts/2+1;
  dataFFT = (fftwf_complex *)fftwf_malloc(npts2*sizeof(fftwf_complex));
  if(dataFFT == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fftZap: fftwf_malloc failed.");
    return 0;
  }
  plan1 = fftwf_plan_dft_r2c_1d(npts, data, dataFFT, FFTW_ESTIMATE);
  plan2 = fftwf_plan_dft_c2r_1d(npts, dataFFT, data, FFTW_ESTIMATE);


  /* Do the forward FFT */
  fftwf_execute(plan1);

  for (i = freq1; i <= freq2; i++) {
    //    printerror(debug, "XXXXX %d", i);
    dataFFT[i] = (0+0*I);
  }

  /* Do the inverse FFT */
  fftwf_execute(plan2);

  fftwf_destroy_plan(plan1);
  fftwf_destroy_plan(plan2);
  fftwf_free(dataFFT); 
  return 1;
}

/* Replace data by its (real) fft with length npts/2. Last sample is
   ignored if the npts is odd. Returns 0 if failed, otherwise 1. */
int replacefft(float *data, long npts, verbose_definition verbose)
{
  long i, npts2;
  fftwf_complex *dataFFT;
  fftwf_plan plan1;

  if(sizeof(int) == 4) {
    if(npts > 2147483644) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR replacefft: requested fft too long.");
      return 0;
    }
  }else if(npts > 32764) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR replacefft: requested fft too long.");
    return 0;
  }
  /* Make npts an even number */
  npts /= 2;
  npts *= 2;

  /* Array size needed by fftw3 */
  npts2 = npts/2+1;
  dataFFT = (fftwf_complex *)fftwf_malloc(npts2*sizeof(fftwf_complex));
  if(dataFFT == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR replacefft: fftwf_malloc failed.");
    return 0;
  }
  plan1 = fftwf_plan_dft_r2c_1d(npts, data, dataFFT, FFTW_ESTIMATE);

  /* Do the forward FFT */
  fftwf_execute(plan1);

  /* Ignore the npts2/2+1 bin. Not 100% sure what's in there, but it
     is at least partially padded. */
  for (i=0; i < npts/2; i++) {
    data[i] = cabs(dataFFT[i]);
  }

  fftwf_destroy_plan(plan1);
  fftwf_free(dataFFT); 
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE
/* Calculates the cross-correlation using fft's. The two data-sets
   should have the same length (ndata), not necessarily a power of
   two. You could use the function crosscorrelation_fft_padding
   instead, which allows you to zero-pad the data first. The
   cross-correlation function is returned in cc. Return value is 0 if
   there is an error. */
int crosscorrelation_fft(float *data1, float *data2, int ndata, float *cc, verbose_definition verbose)
{
  int i, npts2;
  float fac;
  fftwf_complex *dataFFT1, *dataFFT2;
  fftwf_plan plan1, plan2, plan3;

  npts2 = ndata/2+1;
  dataFFT1 = (fftwf_complex *)fftwf_malloc(npts2*sizeof(fftwf_complex));
  dataFFT2 = (fftwf_complex *)fftwf_malloc(npts2*sizeof(fftwf_complex));
  if(dataFFT1 == NULL || dataFFT2 == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR crosscorrelation_fft: fftwf_malloc failed.");
    return 0;
  }
  plan1 = fftwf_plan_dft_r2c_1d(ndata, data1, dataFFT1, FFTW_ESTIMATE);
  plan2 = fftwf_plan_dft_r2c_1d(ndata, data2, dataFFT2, FFTW_ESTIMATE);
  plan3 = fftwf_plan_dft_c2r_1d(ndata, dataFFT1, cc, FFTW_ESTIMATE);

  /* Do the forward FFTs */
  fftwf_execute(plan1);
  fftwf_execute(plan2);

  /* multiply fft by complex conjugate of other fft. The scale factor
     is to correct for the DFT normalization */
  fac = 1.0/(float)ndata;
  for (i=0; i < npts2; i++) {
    dataFFT1[i] *= fac*conj(dataFFT2[i]);
  }

  /* Do the backward FFT */
  fftwf_execute(plan3);

  //  for (i=0; i < npts2; i++) {
  //    cc[i] = creal(1.0/(dataFFT1[i]*conj(dataFFT2[i])));
  //  }


  fftwf_destroy_plan(plan1);
  fftwf_destroy_plan(plan2);
  fftwf_destroy_plan(plan3);
  fftwf_free(dataFFT1); 
  fftwf_free(dataFFT2); 
  return 1;
}
//START REGION DEVELOP
//START REGION RELEASE


/* Calculates the length of the output of the function
crosscorrelation_fft_padding (i.e. what cclength will be set to). This
includes the zero-padding in the calculation.  */
int crosscorrelation_fft_padding_cclength(int ndata, int extrazeropad)
{
  int i, ndata_padded;
  i = (int) (log10(1.0 * (ndata+extrazeropad))/log10(2.0));
  ndata_padded = pow(2.0,(i+1));
  // If ndata is exactly a power of 2, make sure ndata_padded doesn't come out as twice ndata
  if(ndata_padded/2 == ndata)
    ndata_padded = ndata;
  return ndata_padded;
}

/* Calculates the cross-correlation using fft's. The two data-sets
   should have the same length (ndata). If extrazeropad is set, at
   least this amount of zero's are pasted after the data set. If the
   length of the data (plus any extra zeropadding) is not a power of
   two, it will be zero-padded further to a length that is a power of
   two. The cross-correlation function is returned in cc (memory will
   be allocated), which is cclength points long. Return value is 0 if
   there is an error. Verbose determines the number of spaces before
   output. */
int crosscorrelation_fft_padding(float *data1, float *data2, int ndata, int extrazeropad, float **cc, int *cclength, verbose_definition verbose)
{
  float *padded1, *padded2;
  int i, ndata_padded;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Calculating cross correlation of data with length %d\n", ndata);
    if(extrazeropad != 0) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Padding at least %d zero's after data\n", extrazeropad);
    }
  }
  ndata_padded = crosscorrelation_fft_padding_cclength(ndata, extrazeropad);
  if(ndata_padded == ndata) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Data length is already a power of two\n");
    }
    padded1 = data1;
    padded2 = data2;
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Zero-pad it to %d points\n", ndata_padded);
    }
    padded1 = calloc(ndata_padded, sizeof(float));
    padded2 = calloc(ndata_padded, sizeof(float));
    if(padded1 == NULL || padded2 == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR crosscorrelation_fft_padding: Memory allocation error"); 
      return 0;
    }
    memcpy(padded1, data1, ndata*sizeof(float));
    memcpy(padded2, data2, ndata*sizeof(float));
  }
  *cc = (float *)malloc(ndata_padded*sizeof(float));
  *cclength = ndata_padded;
  if(*cc == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR crosscorrelation_fft_padding: Memory allocation error"); 
    return 0;
  }

  /* NR code:   correl(padded1-1, padded2-1, ndata_padded, (*cc)-1); */

  if(crosscorrelation_fft(padded1, padded2, ndata_padded, *cc, verbose) == 0)
    return 0;

  /*  if(lagoutput == 0) {
    for(i = 0; i < ndata; i++) 
      printf("%d %f\n", i, ans[i]);
  }else {
    for(i = ndata/2; i <= ndata-1; i++) 
      printf("%d %f\n", i-ndata, ans[i]);
    for(i = 0; i < ndata/2; i++) 
      printf("%d %f\n", i, ans[i]);
      }*/

  if(ndata_padded != ndata) {
    free(padded1);
    free(padded2);
  }
  return 1;
}

//START REGION DEVELOP

/*
  Determine delay between data1 and data2, which should each be nrbins
  long arrays. The delay is found by determining the location of the
  maximum in the cross-correlation function, which will have a
  precision of 1 bin. The delay is returned in phase. Verbose
  determines the number of spaces before output.

  Return 0 on error, 1 if successful.
 */
int crosscorrelation_find_delay(float *data1, float *data2, long nrbins, double *delay, verbose_definition verbose)
{
  float *cc_function, maxcc;
  int cclength, i;
  if(crosscorrelation_fft_padding(data1, data2, nrbins, 0, &cc_function, &cclength, verbose) == 0) {
    return 0;
  }
  maxcc = cc_function[0];
  *delay = 0;
  for(i = 0; i < cclength; i++) {
    if(cc_function[i] > maxcc) {
      maxcc = cc_function[i];
      *delay = -i/(double)nrbins;
    }
  }
  free(cc_function);
  return 1;
}

typedef struct {
  long nrbins;
  float *template_phases, *data_phases, *template_ampl, *data_ampl;
}internal_FPG_params;


double internal_FPG_root_function(double delay, void *params)
{
  internal_FPG_params *params2 = params;
  long k;
  double ret;
  ret = 0;
  for(k = 0; k < params2->nrbins; k++) {
    ret += k*params2->template_ampl[k]*params2->data_ampl[k]*sin(params2->data_phases[k]-params2->template_phases[k]-k*delay);
  }
  return ret;
}

/* Calculates the delay using the Fourier phase gradient method
   descibed in Taylor 1992. The two data-sets should have the same
   length (ndata), not necessarily a power of two. The delay is
   returned as a delay in phase. Verbose determines the nr of spaces
   before output.

   Return value:
      0 - Critical error occured
      1 - Converged
      2 - Not converged, not no critical error occured
*/
int delay_Fourier_phase_gradient(float *data1, float *data2, int ndata, double *delay, verbose_definition verbose)
{
  int i;
  fftwf_complex *dataFFT1, *dataFFT2;
  fftwf_plan plan1, plan2;
  internal_FPG_params params;
  verbose_definition verbose2;

  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Find delay between two %d long arrays using FFT phase gradient\n", ndata);
  }
  params.nrbins = ndata/2+1;
  dataFFT1 = (fftwf_complex *)fftwf_malloc(params.nrbins*sizeof(fftwf_complex));
  dataFFT2 = (fftwf_complex *)fftwf_malloc(params.nrbins*sizeof(fftwf_complex));
  if(dataFFT1 == NULL || dataFFT2 == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR delay_Fourier_phase_gradient: fftwf_malloc failed.");
    return 0;
  }
  plan1 = fftwf_plan_dft_r2c_1d(ndata, data1, dataFFT1, FFTW_ESTIMATE);
  plan2 = fftwf_plan_dft_r2c_1d(ndata, data2, dataFFT2, FFTW_ESTIMATE);

  params.template_phases = (float *)malloc(params.nrbins*sizeof(float));
  params.data_phases = (float *)malloc(params.nrbins*sizeof(float));
  params.template_ampl = (float *)malloc(params.nrbins*sizeof(float));
  params.data_ampl = (float *)malloc(params.nrbins*sizeof(float));
  if(params.template_phases == NULL || params.data_phases == NULL || params.template_ampl == NULL || params.data_ampl == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR delay_Fourier_phase_gradient: Cannot allocate memory.");
    return 0;
  }

  /* Do the forward FFTs */
  fftwf_execute(plan1);
  fftwf_execute(plan2);

  // Get the phases and amplitudes of the FFT, used in fuction passed to root finder
  for (i=0; i < params.nrbins; i++) {
    params.data_phases[i] = carg(dataFFT1[i]);
    params.template_phases[i] = carg(dataFFT2[i]);
    params.data_ampl[i] = cabs(dataFFT1[i]);
    params.template_ampl[i] = cabs(dataFFT2[i]);
  }

  // Get initial guess for the delay, so we can search for the root bracketed with a small range
  if(crosscorrelation_find_delay(data1, data2, ndata, delay, verbose2) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "Error delay_Fourier_phase_gradient: Finding delay using cross-correlation failed.");
    return 0;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  Initial guess for delay: %lf phase\n", *delay);
  }

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  long iter, max_iter;
  int status;
  double x_lo, x_hi, val1, val2;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  F.function = &internal_FPG_root_function;
  F.params = &params;


  x_lo = 2.0*M_PI*(*delay-1.0/ndata);
  x_hi = 2.0*M_PI*(*delay+1.0/ndata);

  val1 = internal_FPG_root_function(x_lo, &params);
  val2 = internal_FPG_root_function(x_hi, &params);
  if((val1 > 0 && val2 > 0) || (val1 < 0 && val2 < 0)) {
    printwarning(verbose.debug, "WARNING delay_Fourier_phase_gradient: Cannot find suitable start point. Maybe observation is a non-detection? Using maximum of cross-correlation instead.");
    //Avoid an error message, just use maximum of cross-correlation.
    status = GSL_SUCCESS;
  }else {
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);
    /*
      printf ("delay_Fourier_phase_gradient: using %s method\n", gsl_root_fsolver_name (s));
      printf ("delay_Fourier_phase_gradient: %5s [%9s, %9s] %9s %9s\n",
      "iter", "lower", "upper", "root",
      "err(est)");
    */
    iter = 0;
    max_iter = 100;
    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      *delay = gsl_root_fsolver_root(s);
      x_lo = gsl_root_fsolver_x_lower(s);
      x_hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(x_lo, x_hi, 1e-9, 0);
      /*
	if(status == GSL_SUCCESS)
	printf ("Converged:\n");
	printf ("%5ld [%.7f, %.7f] %.7f %+.7f\n",
	iter, x_lo, x_hi,
	*delay,
	x_hi - x_lo);
      */
    }while(status == GSL_CONTINUE && iter < max_iter);
    *delay /= 2.0*M_PI;
    /*
      for(i = 0; i < 100; i++) {
      printf("XXXXX root: %e %e\n", 0.01*2*M_PI*i, internal_FPG_root_function(0.01*2*M_PI*i, &params));
      }
    */

    gsl_root_fsolver_free(s);
  }

  free(params.template_phases);
  free(params.data_phases);
  free(params.template_ampl);
  free(params.data_ampl);
  fftwf_destroy_plan(plan1);
  fftwf_destroy_plan(plan2);
  fftwf_free(dataFFT1); 
  fftwf_free(dataFFT2); 
  if(status == GSL_SUCCESS) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Found a delay of: %lf phase\n", *delay);
    }
    return 1;
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Did not converge, no delay found\n");
    }
    return 2;
  }
}


/* This function takes a pulse of data, npts bins long and convolves
   it with an exponential tail (fast-rise-slow-decay) function. The
   timeconst is in bins. Note that the end of the scatter tail ends up
   at the start of the pulse.

   Returns 0 if failed, otherwise 1. */
int convolveScatterTail_singlepulse(float *data, int npts, float timeconst, verbose_definition verbose)
{
  int i, npts2;
  float fac, *scatter;
  double norm, value;
  fftwf_complex *dataFFT, *scatterFFT;
  fftwf_plan plan1, plan2, plan3;

  //Allocate memory
  npts2 = npts/2+1;
  scatter = (float *)malloc(npts*sizeof(float));
  dataFFT = (fftwf_complex *)fftwf_malloc(npts2*sizeof(fftwf_complex));
  scatterFFT = (fftwf_complex *)fftwf_malloc(npts2*sizeof(fftwf_complex));
  if(dataFFT == NULL || scatter == NULL || scatterFFT == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR convolveScatterTail: Memory allocation failed.");
    return 0;
  }

  //  printerror(debug, "XXXXX timeconst = %f", timeconst);

  // The convolve function is: exp(-n/timeconst)
  norm = 0;
  for(i = 0; i < npts; i++) {
    value = exp(-(float)i/timeconst);
    scatter[i] = value;
    norm += value;
  }
  norm = 1.0/norm;
  for(i = 0; i < npts; i++) {
    scatter[i] *= norm;   
  }

  //Create fft plans
  plan1 = fftwf_plan_dft_r2c_1d(npts, data, dataFFT, FFTW_ESTIMATE);
  plan2 = fftwf_plan_dft_r2c_1d(npts, scatter, scatterFFT, FFTW_ESTIMATE);
  plan3 = fftwf_plan_dft_c2r_1d(npts, dataFFT, data, FFTW_ESTIMATE);

  // Do the forward FFT of data and scatter tail
  fftwf_execute(plan1);
  fftwf_execute(plan2);

  /* Convolution is multiplication of fft's. The scale factor is to
     correct for the DFT normalization */
  fac = 1.0/(float)npts;
  for (i=0; i < npts2; i++) {
    dataFFT[i] *= fac*scatterFFT[i];
    //    dataFFT[i] = fac*scatterFFT[i];
  }

  // Do the inverse FFT
  fftwf_execute(plan3);

  /*
  for(i = 0; i < npts; i++) {
    data[i] = exp(-(float)i/timeconst);
  }
  */

  fftwf_destroy_plan(plan1);
  fftwf_destroy_plan(plan2);
  fftwf_destroy_plan(plan3);
  fftwf_free(dataFFT); 
  fftwf_free(scatterFFT); 
  free(scatter);
  return 1;
}



/* This function takes a pulse of data, npts bins long and convolves
   it with tophat function of full width timeconst (in bins). Verbose
   level detirmines nr of spaces before output.

   Returns 0 if failed, otherwise 1. */
int tophat_smooth(float *data, int npts, float timeconst, verbose_definition verbose)
{
  int i, npts2;
  float fac, *boxcarFFT;
  double norm, value;
  fftwf_complex *dataFFT;
  fftwf_plan plan1, plan2;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Smooting with a width of %f bins\n", timeconst);
  }

  //Allocate memory
  npts2 = npts/2+1;
  boxcarFFT = (float *)malloc(npts2*sizeof(float));
  dataFFT = (fftwf_complex *)fftwf_malloc(npts2*sizeof(fftwf_complex));
  if(dataFFT == NULL || boxcarFFT == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR tophat_smooth: Memory allocation failed.");
    return 0;
  }

  //  printerror(debug, "XXXXX timeconst = %f", timeconst);

  // The fft of convolve function is: sinc(-timeconst*n/2) times a constant, which is purely real
  norm = 0;
  for(i = 0; i < npts2; i++) {
    if(i == 0) {
      value = 1.0;
    }else {
      value = i*timeconst*M_PI/(double)npts;
      value = sin(value)/value;
    }
    boxcarFFT[i] = value;
    norm += value*value;
  }
  // Not entirely sure why, but result is best if not normalised (energy conserved)
  /*
  norm = 1;
  for(i = 0; i < npts2; i++) {
    boxcarFFT[i] *= norm;   
  }
  */

  //Create fft plans
  plan1 = fftwf_plan_dft_r2c_1d(npts, data, dataFFT, FFTW_ESTIMATE);
  plan2 = fftwf_plan_dft_c2r_1d(npts, dataFFT, data, FFTW_ESTIMATE);

  // Do the forward FFT of data 
  fftwf_execute(plan1);

  /* Convolution is multiplication of fft's. The scale factor is to
     correct for the DFT normalization */
  fac = 1.0/(float)npts;
  for (i=0; i < npts2; i++) {
    dataFFT[i] *= fac*(boxcarFFT[i]+0*I);
  }

  // Do the inverse FFT
  fftwf_execute(plan2);

  fftwf_destroy_plan(plan1);
  fftwf_destroy_plan(plan2);
  fftwf_free(dataFFT); 
  free(boxcarFFT);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done\n");
  }
  return 1;
}
