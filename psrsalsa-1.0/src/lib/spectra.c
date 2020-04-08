//START REGION RELEASE
#include <math.h>
#include <string.h>
#include "psrsalsa.h"
//START REGION DEVELOP
/* Comment this line out if you don't want to use fftw3 */
//START REGION RELEASE

#define USEFFTW3     1  

//START REGION DEVELOP
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

//START REGION RELEASE
#ifdef USEFFTW3 
  #include <complex.h>
  #include <fftw3.h>
#else 
  #include "nr.h"
  #include "nrutil.h"
#endif



/* 
   Calculate the 2DFS of data, which is an array of nrx by nry
   points. The on-pulse region should be defined and should have a
   length which is a power of two. The region number to calculate the
   2DFS for can be specified. Only multiples of fft_size are used in
   the y-direction). The 2DFS's of length fft_size are added together
   in 2dfs, which should be an array of nrx by 1+fft_size/2
   points. All the defined (on-pulse) regions are ignored for the
   off-pulse noise subtraction. The off-pulse subtraction will be done
   using an off-pulse region with a with equal to the selected
   region. This off-pulse region start at bin 0, if it fits. Otherwise
   the start bin of the off-pulse region is increased until a suitable
   off-pulse region is identified. If there is no suitable off-pulse
   region, no off-pulse noise subtraction will be applied.  If verbose
   is set it shows the number of blocks that is calculated.

   0 = nothing done 
   1 = ok
*/
int calc2DFS(float *data, long nry, long nrx, unsigned long fft_size, float *twodfs, pulselongitude_regions_definition *onpulse, int region, verbose_definition verbose)
{
  unsigned long nr_fftblocks;
  float junk_float;
  long junk_int, i, nf, nb, nb2, np, nrx2, bin_offpulse_left;
  int ok;

  #ifdef USEFFTW3 
    float *inputdata, pwr;
    fftwf_complex *fftdata; 
    fftwf_plan plan;
  #else
    float ***inputdata, **speq;
  #endif

  // Not strictly necessary, but avoids warning when optimizing code
  ok = 0;
  pwr = 0;

  junk_float = log(fft_size)/log(2);
  junk_int = junk_float;
  junk_float = pow(2, junk_int);
  if(fabs(junk_float-fft_size) > 0.1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: fft length is not a power of two!");
    return 0;
  }
  if(onpulse->nrRegions < region) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: selected region is not defined");
    return 0;
  }
  if(onpulse->bins_defined[region] == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: Region is not defined in bins");
    //    printRegions(onpulse);
    return 0;
  }

  nrx2 = onpulse->right_bin[region]-onpulse->left_bin[region]+1;

  /* FFTW3 doesn't appear to mind if horizontal range is a power of two, but it should be an even number in order to avoid confusion about the labelling of the x-axis. */
  #ifdef USEFFTW3 
    junk_int = nrx2 / 2;
    junk_int *= 2;
    junk_int -= nrx2;
    if(junk_int != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calc2DFS: Length onpulse region should be an even number.");
      return 0;
    }
  #else
    junk_float = log(nrx2)/log(2);
    junk_int = junk_float;
    junk_float = pow(2, junk_int);
    if(fabs(junk_float-nrx2) > 0.1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calc2DFS: Length onpulse region is not a power of two. Please recompile using fftw3 support.");
      return 0;
    }
  #endif

  /* Calculate nr of spectra we're going to calculate */
  nr_fftblocks = nry/fft_size;
  if(nr_fftblocks == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: Cannot calculate 2dfs for %ld pulses (fft size = %ld)", nry, fft_size);
    return 0;
  }
  if(verbose.verbose) printf("Calculating 2DFS (%ld blocks)\n", nr_fftblocks);
  if(onpulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: Onpulse region is undefined.");
    return 0;
  }
  if(onpulse->nrRegions == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calc2DFS: Onpulse region is undefined.");
    return 0;
  }

  #ifdef USEFFTW3 
    inputdata = (float *)malloc(nrx2*fft_size*sizeof(float));
    fftdata = (fftwf_complex *)fftwf_malloc(nrx2*(1+fft_size/2)*sizeof(fftwf_complex));
    if(fftdata == NULL || inputdata == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calc2DFS: fftwf_malloc failed.");
      return 0;
    }
    plan = fftwf_plan_dft_r2c_2d(nrx2, fft_size, inputdata, fftdata, FFTW_ESTIMATE); 
  #else
    inputdata = f3tensor(1,1,1,nrx2,1,fft_size);
    speq = matrix(1,1,1,2*nrx2); 
    if(inputdata == NULL || speq == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calc2DFS: Cannot allocate memory");
      return 0;
    }
  #endif

  bin_offpulse_left = -1;
  if(onpulse != NULL) {
    for(bin_offpulse_left = 0; bin_offpulse_left < nrx; bin_offpulse_left++) {
      ok = 1;
      for(i = 0; i < nrx2; i++) {
	if(i+bin_offpulse_left >= nrx) {
	  ok = 0;
	  break;
	}
	if(checkRegions(i+bin_offpulse_left, onpulse, 0, verbose) != 0) {
	  ok = 0;
	  break;
	}
      }
      if(ok) {
	break;
      }
    }
    if(ok == 0)
      bin_offpulse_left = -1;
  }
  if(bin_offpulse_left >= 0) {
    if(verbose.verbose) printf("  Found suitable offpulse region (%ld %ld).\n", bin_offpulse_left, bin_offpulse_left+nrx2-1);
  }else {
    if(verbose.verbose) printwarning(verbose.debug, "  WARNING calc2DFS: Didn't found suitable offpulse region.");
  }


  /* Clean memory, so we can add the spectra */
  for(i = 0; i < nrx2*(1+fft_size/2); i++)
    twodfs[i] = 0;

  /* Loop over the blocks of which we're calculating the spectra */
  for(nf = 0; nf < nr_fftblocks; nf++) {
    /* Loop over the bins */
    for(nb = 0; nb < nrx2; nb++) {
      /* Loop over pulses */
      for(np = 0; np < fft_size; np++) {
	/* Prepare input data */
#ifdef USEFFTW3 
	inputdata[nb*fft_size + np]  = data[(nf*fft_size+np)*nrx + nb + onpulse->left_bin[region]]; 
#else
	inputdata[1][nb+1][np+1] = data[(nf*fft_size+np)*nrx + nb + onpulse->left_bin[region]]; 
#endif
      }
    }
    

#ifdef USEFFTW3 
    fftwf_execute(plan);
    for(nb = 0; nb < nrx2; nb++) {
      nb2 = nb+nrx2/2;
      if(nb2 >= nrx2)
	nb2 -= nrx2;
      for(np = 0; np < fft_size/2+1; np++) {
	if(np != 0) {                                   
//START REGION DEVELOP
/* I assumed this is to ignore DC offset. It didn't had the curly brackets, making it add the square of an undefined number for the first point. Caused problems if it is a nan. */
//START REGION RELEASE
	  pwr = cabs(fftdata[nb*(fft_size/2+1)+np]);
	  twodfs[np*nrx2+nb2] += pwr*pwr;
	}
      }
    }
#else
    rlft3(inputdata, speq, 1, nrx2, fft_size, 1);
    for(nb = 0; nb < nrx2; nb++) {
      nb2 = nb+nrx2/2;
      if(nb2 >= nrx2)
	nb2 -= nrx2;
      for(np = 0; np < fft_size/2; np++) {
	if(np != 0) 
	  twodfs[np*nrx2+nb2] += inputdata[1][nb+1][2*np+1]*inputdata[1][nb+1][2*np+1]+inputdata[1][nb+1][2*np+2]*inputdata[1][nb+1][2*np+2]; 
      }
      twodfs[(fft_size/2)*nrx2+nb2] += speq[1][2*nb+1]*speq[1][2*nb+1]+speq[1][2*nb+2]*speq[1][2*nb+2];
    }
#endif


    /**********************************************
     *         Subtract offpulse region           *
     **********************************************/ 

    if(bin_offpulse_left >= 0) {
      for(nb = 0; nb < nrx2; nb++) {
	/* Loop over pulses */
	for(np = 0; np < fft_size; np++) {
	  /* Prepare input data */
#ifdef USEFFTW3 
	  inputdata[nb*fft_size + np]  = data[(nf*fft_size+np)*nrx + nb + bin_offpulse_left]; 
#else
	  inputdata[1][nb+1][np+1] = data[(nf*fft_size+np)*nrx + nb + bin_offpulse_left]; 
#endif
	}
      }
      
#ifdef USEFFTW3 
    fftwf_execute(plan);
    for(nb = 0; nb < nrx2; nb++) {
      nb2 = nb+nrx2/2;
      if(nb2 >= nrx2)
	nb2 -= nrx2;
      for(np = 0; np < fft_size/2+1; np++) {
	if(np != 0) 
	  pwr = cabs(fftdata[nb*(fft_size/2+1)+np]);
	  twodfs[np*nrx2+nb2] -= pwr*pwr;
      }
    }
#else
      rlft3(inputdata, speq, 1, nrx2, fft_size, 1);
      for(nb = 0; nb < nrx2; nb++) {
      	nb2 = nb+nrx2/2;
	if(nb2 >= nrx2)
	  nb2 -= nrx2;
	for(np = 0; np < fft_size/2; np++) {
	  if(np != 0) 
	    twodfs[np*nrx2+nb2] -= inputdata[1][nb+1][2*np+1]*inputdata[1][nb+1][2*np+1]+inputdata[1][nb+1][2*np+2]*inputdata[1][nb+1][2*np+2]; 
	}
	twodfs[(fft_size/2)*nrx2+nb2] -= speq[1][2*nb+1]*speq[1][2*nb+1]+speq[1][2*nb+2]*speq[1][2*nb+2];
      }
#endif 
    }

    if(nr_fftblocks > 1 && verbose.nocounters == 0) {
      printf("Block %ld of the %ld     \r", nf+1, nr_fftblocks);
      fflush(stdout);
    }
  }
 
  if(nr_fftblocks > 1 && verbose.nocounters == 0)   
    printf("Done                       \n");

#ifdef USEFFTW3 
  fftwf_destroy_plan(plan);
  fftwf_free(fftdata); 
  free(inputdata);
#else
  free_matrix(speq, 1,1,1,2*nrx2);
  free_f3tensor(inputdata,1,1,1,nrx2,1,fft_size);
#endif
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* 
   Note there is a wrapper function calcLRFS_padding() which allows
   you to do zero-padding. Calculate the LRFS of data, which is an
   array of nrx by nry points. The powers of the fft's of length
   fft_size are added together in lrfs, which should be an array of
   nrx by fft_size/2+1 points (which corresponds to the frequencies
   0,1/fft_size,2/fft_size, ..., 0.5). If subtractDC is specified, the
   DC component is subtracted from the LRFS. If
   avrg_offpulse_lrfs_power != NULL, it will be set to the average
   offpulse power in the obtained lrfs (so the average power per
   spectral/longitude bin in *lrfs. phase_track should be nrx points
   long (if calcPhaseTrack is specified), which contains the
   vertically collapsed subpulse phase for each bin in degrees. The
   range in fluctuation frequencies specified by freq_min and freq_max
   is used for the track calculation. If phase_track_phases is not
   NULL, the phase offsets between the phase_tracks in subsequent
   fft_blocks are stored (so it should be an array of nry/fft_size
   floats). If track_only_first_region is zero, it means that the full
   onpulse region is used to align the phases obtained from the
   individual blocks. If set to 1, only the first on-pulse region is
   used to do the correlations, while the rest of the regions are
   still used for the noise subtraction of the LRFS. The noise
   subtraction is based on all not-selected pulse-longitude bins. If
   calcsubpulseAmplitude is set (subpulseAmplitude should have a
   length nrx otherwise), the power of the drifting subpulses is
   stored for each longitude bin (as the stddev in selected frequency
   range divided by the maximum stddev in the profile). These are the
   averages of the different spectral channels for a given block, not
   weighted by anything. The final subpulse phase track is a weighted
   average however. Set the mask_freqs flag to mask out everything
   outside the frequency range. var_rms is the rms (after subtracting
   the mean, using the separate blocks of data rather than the final
   lrfs for which the different blocks are summed) of the samples in
   the lrfs in the off-pulse region.
//START REGION DEVELOP
   If inverseFFT is set it calculates the inverse transform. In
   that case the data is overwritten. I'm not sure if the
   normalization is done correctly. You might be better of to select a
   single bin.
//START REGION RELEASE 
   If regions is defined (on-pulse), than the off-pulse noise
   contribution is subtracted and the rms of the offpulse rms is
   calculated (used by calcModindex). If verbose is set it shows the
   number of blocks that is calculated. If argc > 0 and argv != NULL,
   it will search the command line for the -p3zap option and zap these
   regions from the lrfs (can be specified in cpp or in bins). The 
   debug option gives more output.

   0 = nothing done 
   1 = ok
*/
int calcLRFS(float *data, long nry, long nrx, unsigned long fft_size, float *lrfs, int subtractDC, float *avrg_offpulse_lrfs_power, float *phase_track, float *phase_track_phases, int calcPhaseTrack, float freq_min, float freq_max, int track_only_first_region, float *subpulseAmplitude, int calcsubpulseAmplitude, int mask_freqs, int inverseFFT, pulselongitude_regions_definition *regions, float *var_rms, int argc, char **argv, verbose_definition verbose)
{
  long fftblock, binnr, pulsenr;
  unsigned long nr_fftblocks;
  float *data1, *lrfs_tmp, pwr, pwrtot, freq, var_mean;
  double avr_power_subtracted_from_lrfs;
  long i, j, n, tot_var_rms_samples;
  float zapmin, zapmax, p3;
  int k;
  float *phase_track_complex_template, *phase_track_complex;
  //  int freq_bin_selected;
  long ok, itteration, nrphasetracks, nspecbins;
#ifdef USEFFTW3
  fftwf_plan plan1;
  fftwf_plan plan2;
#endif


  if(regions != NULL) {
    *var_rms = 0;
    var_mean = 0;
    tot_var_rms_samples = 0;
    avr_power_subtracted_from_lrfs = 0;
  }

  // Shouldn't be necessary, but gets rid of warnings when optimising code
  phase_track_complex_template = NULL;
  phase_track_complex = NULL;
  nspecbins = 0;

//START REGION DEVELOP
  /* Calculate nr of spectra we're going to calculate (per longitude bin). */
//START REGION RELEASE
  nr_fftblocks = nry/fft_size;
  if(verbose.verbose) printf("Calculating LRFS (%ld blocks)\n", nr_fftblocks);
  if(nr_fftblocks == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calcLRFS: Cannot calculate lrfs for %ld pulses (smaller than fft size = %ld)", nry, fft_size);
    return 0;
  }

//START REGION DEVELOP
  /* Allocate some temporary memory. data1 will contain a column from
     the pulse-stack (we have to do a transpose) which will go into
     the fft routine. lrfs_tmp will contain the lrfs of a single block
     of data. The extra 2 values is because there will be
     (fft_size/2+1) complex valued spectral bins. */
//START REGION RELEASE
  data1 = (float *)malloc((fft_size+2)*sizeof(float));
  lrfs_tmp = (float *)malloc(nrx*(fft_size/2+1)*sizeof(float));
  if(data1 == NULL || lrfs_tmp == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calcLRFS: Cannot allocate memory");
    return 0;
  }


#ifdef USEFFTW3
//START REGION DEVELOP
  /* Do an inplace transform, like realft of NR. */
//START REGION RELEASE
  plan1 = fftwf_plan_dft_r2c_1d(fft_size, data1, (fftwf_complex *)data1, FFTW_ESTIMATE);
  plan2 = fftwf_plan_dft_c2r_1d(fft_size, (fftwf_complex *)data1, data1, FFTW_ESTIMATE);
//START REGION DEVELOP
  /* FFTW routine takes integer as fft length */
//START REGION RELEASE
  if(fft_size > 2147483640) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calcLRFS: requested fft too long.");
    return 0;
  }
#endif

//START REGION DEVELOP
  /* Calculate the number of spectral bins selected and store the last
     frequency bin (under the assumption at this point that there is
     only one spectral bin selected, otherwise this is not very
     useful). */
//START REGION RELEASE
  if(calcPhaseTrack || inverseFFT || calcsubpulseAmplitude) {
    nspecbins = 0;
    for(i = 0; i <= fft_size/2; i++) {
      freq = i/(float)fft_size;
      if(freq >= freq_min && freq <= freq_max) {
	nspecbins++;
	//	freq_bin_selected = i;
      }
    }
    if(verbose.verbose) {
      printf("  Using %ld spectral bins for the calculation of the subpulse track.\n", nspecbins);
    }
    if(nspecbins == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calcLRFS: Cannot calculate phase track for %ld frequency bins", nspecbins);
      return 0;
    }else if(verbose.debug) {
      for(i = 0; i <= fft_size/2; i++) {
	freq = i/(float)fft_size;
	if(freq >= freq_min && freq <= freq_max) {
	  printf("    This is frequency bin with a frequency of %ld/%ld=%f cpp (P3 = %f P1).\n", i, fft_size, freq, 1.0/freq);
	}
      }	
    }
  }

//START REGION DEVELOP
  /* If we're going to calculate the phase track, we'll have to store
     all the complex phases (for each longitude bin, each fft block
     and each frequency bin. This is required so we can coherently add
     the phases later. */
//START REGION RELEASE
  if(calcPhaseTrack || calcsubpulseAmplitude) {
    phase_track_complex = (float *)malloc(2*nrx*nspecbins*nr_fftblocks*sizeof(float));
    phase_track_complex_template = (float *)malloc(2*nrx*sizeof(float));
    if(phase_track_complex == NULL || phase_track_complex_template == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calcLRFS: Cannot allocate memory");
      return 0;
    }
//START REGION DEVELOP
    /* To check for memory uninitialized I set these arrays to zero, it should not be necessary */
    /*    for(i = 0; i < 2*nrx*nspecbins*nr_fftblocks; i++)
      phase_track_complex[i] = 0;
    */
    /*    for(i = 0; i < 2*nrx; i++)
      phase_track_complex_template[i] = 0;
    */
//START REGION RELEASE
  }

//START REGION DEVELOP
  /* Empty the lrfs, so we can add the spectra of the individual
     blocks later. */
//START REGION RELEASE
  for(i = 0; i < nrx*(1+fft_size/2); i++)
    lrfs[i] = 0;

//START REGION DEVELOP
  /* This is a counter that loops over the number of phase tracks
     which have been calculated. There will be nspecbins*nr_fftblocks
     phase tracks in total. */
//START REGION RELEASE
  nrphasetracks = 0;
//START REGION DEVELOP
  /* Loop over the blocks of pulses for which we're calculating the spectra */
//START REGION RELEASE
  for(fftblock = 0; fftblock < nr_fftblocks; fftblock++) {

//START REGION DEVELOP
    /* Set the temporary array for the lrfs of this block to zero */
//START REGION RELEASE
    for(i = 0; i < nrx*(fft_size/2+1); i++)
      lrfs_tmp[i] = 0;
//START REGION DEVELOP
    /* Loop over the longitude bins */
//START REGION RELEASE
    for(binnr = 0; binnr < nrx; binnr++) {
//START REGION DEVELOP
      /* Get column out of pulse-stack block and calculate DC component */
//START REGION RELEASE
      pwrtot = 0;
      for(pulsenr = 0; pulsenr < fft_size; pulsenr++) {
	data1[pulsenr] = data[(fftblock*fft_size+pulsenr)*nrx + binnr];
	pwrtot += data1[pulsenr];
      }

//START REGION DEVELOP
      /* Calculate fft of the column of (real valued) data resulting in a complex spectrum */
//START REGION RELEASE
#ifdef USEFFTW3
      fftwf_execute(plan1);
#else
      realft(data1-1, fft_size, 1);
//START REGION DEVELOP
      /* Unpack the (real valued) Nyquest frequency, which is stored
	 as the imaginary part of the real valued DC bin. */
      /* Nyquist freq. was not calculated correctly before, this is now fixed. */
//START REGION RELEASE
      data1[2*(fft_size/2)]   = data1[1];
      data1[2*(fft_size/2)+1] = 0;
      data1[1] = 0;
#endif

//START REGION DEVELOP
      /* Set the phase track for this particular block and phase bin to zero */
//START REGION RELEASE
      if(calcPhaseTrack) {
	phase_track_complex[2*(nrx*nrphasetracks + binnr)]     = 0;
	phase_track_complex[2*(nrx*nrphasetracks + binnr)+1]   = 0;
      }
    
//START REGION DEVELOP
      /* Loop over spectral bins */
//START REGION RELEASE
      for(i = 0; i <= fft_size/2; i++) {
//START REGION DEVELOP
	/* Frequency in cpp of current spectral channel */
//START REGION RELEASE
	freq = i/(float)fft_size;
//START REGION DEVELOP
	/* Calculate the amplitude squared of the complex fft. */
//START REGION RELEASE
	pwr  = data1[2*i]*data1[2*i] + data1[2*i+1]*data1[2*i+1];
//START REGION DEVELOP
	/* Subtract DC component using what we measured from the pulse
	   stack. If everything is correct, the DC channel should be
	   zero per definition. */
//START REGION RELEASE
	if(i == 0 && subtractDC) {
	  pwr -= pwrtot*pwrtot;
	}
//START REGION DEVELOP
	/* Calculate power of lrfs if within specified frequency
	   range. If mask_freqs is not set the power is always added
	   to lrfs. */
//START REGION RELEASE
	if(mask_freqs == 0 || (freq >= freq_min && freq <= freq_max)) {
	  lrfs_tmp[i*nrx+binnr]   = pwr;
//START REGION DEVELOP
	  /* If in specified frequency range, store the phase
	     information to calculate the phase track later on. */
//START REGION RELEASE
	  if(calcPhaseTrack) {
	    if(freq >= freq_min && freq <= freq_max) {
	      phase_track_complex[2*(nrx*nrphasetracks + binnr)]     += data1[2*i];
	      phase_track_complex[2*(nrx*nrphasetracks + binnr)+1]   += data1[2*i+1];
	      nrphasetracks++;
	    }
	  }
//START REGION DEVELOP
	}
//	}else if(inverseFFT) { /* Set zapped freq. channels to zero if we want to calculate inverse fft of result. */
	if(inverseFFT) { // Set zapped freq. channels to zero if we want to calculate inverse fft of result.
	  // replace else if with if, make sure the else if is still being achieved and add -p3zap
	  if(mask_freqs && (freq < freq_min || freq > freq_max)) {
	    data1[2*i] = 0;
	    data1[2*i+1] = 0;
	  }
	  // Code related to -p3zap used later again
	  if(argc > 0 && argv != NULL) {
	    int index;
	    for(index = 0; index < argc-1; index++) {
	      if(strcmp(argv[index], "-p3zap") == 0) {
		j = sscanf(argv[index+1], "%f %f", &zapmin, &zapmax);
		if(j != 2) {
		  fflush(stdout);
		  printerror(verbose.debug, "Cannot parse -p3zap option.");
		  exit(0);
		}
		for(j = 0; j <= (fft_size/2); j++) {
		  if(zapmin > 0.9 || zapmax > 0.9) {
		    p3 = j;   // If the boundary exceeds 0.5 cpp, then the input parameters are assumed to be in bins
		  }else {
		    p3 = j/(float)(fft_size);   
		  }
		  if(p3 >= zapmin && p3 <= zapmax) {
		    data1[2*j] = 0;
		    data1[2*j+1] = 0;
		  }
		}
	      }
	    }
	  }

	  //START REGION RELEASE
	}
      }           /* End loop over spectral channels. */

//START REGION DEVELOP
      if(inverseFFT) {
	/*	if(inverseFFT == 2) {
	  if(nspecbins != 1) {
	  fflush(stdout);
	    printwarning(verbose.debug, "WARNING: frequency shifting of the LRFS unlikely to be sensible if more than one frequency channel is selected!!!!!!!!!!!!.");
	  }
	  for(i = 0; i <= fft_size/2; i++) {
	    j = i - freq_bin_selected;
	    if(j < 0)
	      j += fft_size/2;
	    if(i == freq_bin_selected) {
	      data1[2*j] = data1[2*i];
	      data1[2*j+1] = data1[2*i+1];
	    }else {
	      data1[2*j] = 0;
	      data1[2*j+1] = 0;
	    }
	  }
	  data1[0] = sqrt(data1[1]*data1[1] + data1[0]*data1[0]);
	  data1[1] = 0;
	  }else {*/
	  /*}*/

#ifdef USEFFTW3
	fftwf_execute(plan2);
#else
	/* Pack the real valued Nyquist freq in the imaginary part of
	   the real valued DC. */
	data1[1] = data1[2*(fft_size/2)];
	realft(data1-1, fft_size, -1);
#endif

	for(pulsenr = 0; pulsenr < fft_size; pulsenr++) {
	  /* The 2/fft_size is a normalisation factor for NR, 1/fft_size for fftw */
#ifdef USEFFTW3
	  data[(fftblock*fft_size+pulsenr)*nrx + binnr] = data1[pulsenr]/(float)fft_size;
#else
	  data[(fftblock*fft_size+pulsenr)*nrx + binnr] = data1[pulsenr]*2.0/(float)fft_size;
#endif
	}
      } /* End of inverseFFT if statement */
//START REGION DEVELOP
      /* This is necessary because most inner loop is not longitude bins. */
//START REGION RELEASE
      nrphasetracks -= nspecbins;
    } // Loop over bin number

//START REGION DEVELOP
    /* If regions is specified, subtract offpulse spectral power from
       the lrfs. */
//START REGION RELEASE
    if(regions != NULL) {
      if(regions->nrRegions > 0) {
	for(i = 0; i <= fft_size/2; i++) {
	  pwrtot = 0;
	  n = 0;
	  for(binnr = 0; binnr < nrx; binnr++) {
	    if(checkRegions(binnr, regions, 0, verbose) == 0) {
	      pwrtot += lrfs_tmp[i*nrx+binnr];
	      n++;
	      *var_rms += lrfs_tmp[i*nrx+binnr]*lrfs_tmp[i*nrx+binnr];
	      var_mean += lrfs_tmp[i*nrx+binnr];
	      tot_var_rms_samples++;
	    }
	  }
	  if(n > 0) {
	    pwrtot /= (float)n;
	  }else if(i == 0 && fftblock == 0) {
	    fflush(stdout);
	    printwarning(verbose.debug, "WARNING: onpulse region is defined, but no offpulse region is available!");
	  }
	  for(binnr = 0; binnr < nrx; binnr++) {
	    lrfs_tmp[i*nrx+binnr] -= pwrtot;
	  }
	  avr_power_subtracted_from_lrfs += pwrtot;
	}
      }
    }
  

//START REGION DEVELOP
    /* Add the power to lrfs if within specified frequency range. If
       mask_freqs is not set the power is always added to lrfs. The
       power is the square of the spectral coefficients, hence
       according to Parseval's theorem this is proportional to the
       variance of the original time samples. This suggest we can add
       the different blocks, as the total variance is given by the sum
       of the squares of the time samples. */
//START REGION RELEASE
    for(binnr = 0; binnr < nrx; binnr++) {
      for(i = 0; i <= fft_size/2; i++) {
	lrfs[i*nrx+binnr] += lrfs_tmp[i*nrx+binnr];
      }
    }
  
//START REGION DEVELOP
    /* After finishing nspecbins phase tracks, goto next block. */
//START REGION RELEASE
    nrphasetracks += nspecbins;
    if(verbose.verbose && verbose.nocounters == 0) {
      printf("  Block %ld of the %ld     \r", fftblock+1, nr_fftblocks);
      fflush(stdout);
    }
  }


  if(regions != NULL) {
    *var_rms /= (float)tot_var_rms_samples;
    var_mean /= (float)tot_var_rms_samples;
    *var_rms = sqrt((*var_rms)-var_mean*var_mean);
  }

//START REGION DEVELOP
  /* Now calculate the phase track. We have all the phases as
     calculated above. The problem now is that the phases will be
     different for different blocks and frequency channels. Therefore
     adding the complex numbers of different blocks will not work,
     because it will average out to zero. The phase offset between
     different blocks should be independent of longitude bin. So we
     can correlate a template with all the individual phase tracks to
     find this phase offset. Then by derotating the individual phase
     tracks we can coherently add the phase tracks to obtain a better
     template etc. As initial guess we simply incoherently add all the
     available phase tracks. I hope there should always be some signal
     left. */
//START REGION RELEASE

  if(calcPhaseTrack || calcsubpulseAmplitude) {

    int itteration_max = 100;
    double complex corr, phase, *phases;
    float total_phase_offset;
  
    if(verbose.verbose) printf("  Coherently add phase tracks\n");

//START REGION DEVELOP
    /* Reserve memory to store the phases of the individual tracks. */
//START REGION RELEASE
    phases = (double complex *)malloc(nrphasetracks*sizeof(double complex));
    if(phases == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calcLRFS: Cannot allocate memory");
      return 0;
    }
    for(i = 0; i < nrphasetracks; i++) {
      phases[i] = 1;
    }
  
//START REGION DEVELOP
    /* Do one extra loop, because last loop is only to make new
       template. It is not doing a new correlation anymore. */
//START REGION RELEASE
    for(itteration = 0; itteration < itteration_max+1; itteration++) {
//START REGION DEVELOP
      /* Keep track of the total phase offset so we know when to stop. */
//START REGION RELEASE
      total_phase_offset = 0;

//START REGION DEVELOP
      /* Add the individual phase tracks together by assuming they are
	 in phase. At first that will not be true, but we will
	 itteratively keep derotating the phase tracks to make this a
	 coherent addition. */
//START REGION RELEASE
      for(binnr = 0; binnr < 2*nrx; binnr++) {
	phase_track_complex_template[binnr] = 0;
	for(i = 0; i < nrphasetracks; i++) {
	  phase_track_complex_template[binnr] += phase_track_complex[binnr+i*2*nrx];
	}
      }

//START REGION DEVELOP
      /* As long as we're not in the last loop, we can calculate new
	 phases. Otherewise there is not much use as the addition is
	 done above. */
//START REGION RELEASE
      if(itteration < itteration_max) {
//START REGION DEVELOP
	/* Loop over all phase tracks. */
//START REGION RELEASE
	for(i = 0; i < nrphasetracks; i++) {

//START REGION DEVELOP
	  /* Now calculate the correlation between the phase track and
	     the template: corr =
	     phase_track_complex_template*conj(phase_track_complex) In
	     effect this is an inproduct giving the phase offset as a
	     complex number.*/
//START REGION RELEASE
	  corr = 0;
	  for(binnr = 0; binnr < nrx; binnr++) {
	    ok = 1;
	    if(regions != NULL) {
	      if(regions->nrRegions > 0) {
		if(track_only_first_region) {
		  if(checkRegions(binnr, regions, 1, verbose) == 0) {
		    ok = 0;
		  }
		}else {
		  if(checkRegions(binnr, regions, 0, verbose) == 0) {
		    ok = 0;
		  }
		}
	      }
	    }
	    if(ok) {
	      corr += (phase_track_complex_template[2*binnr]+I*phase_track_complex_template[2*binnr+1])*(phase_track_complex[2*binnr+i*2*nrx]-I*phase_track_complex[2*binnr+1+i*2*nrx]);
	    }
	  }
//START REGION DEVELOP
	  /* If zero, there was no signal, so there is nothing to do */
//START REGION RELEASE
	  if(cabs(corr) > 0) {
//START REGION DEVELOP
	    /* Make the correlation have a length 1, so it becomes a
	       pure rotation. */
//START REGION RELEASE
	    corr /= cabs(corr);
	    total_phase_offset += fabs(carg(corr));
//START REGION DEVELOP
	    /*  Now adjust the phases of the individual phase tracks so
		we can coherently add them at the start of the outer
		loop. This is simply done by multiplying the complex
		subpulse phase by the complex phase offset (corr).
	    */
//START REGION RELEASE
	    for(binnr = 0; binnr < nrx; binnr++) {
	      phase = phase_track_complex[2*binnr+i*2*nrx] + I*phase_track_complex[2*binnr+1+i*2*nrx];
	      phase *= corr;
	      phase_track_complex[2*binnr+i*2*nrx] = creal(phase);
	      phase_track_complex[2*binnr+1+i*2*nrx] = cimag(phase);
//START REGION DEVELOP
	      /* Also store the total phase rotation per phase track */
//START REGION RELEASE
	    }
	    phases[i] *= corr;
	  }
	}
	total_phase_offset *= 180.0/(float)(M_PI*nrphasetracks);
	if(verbose.debug) printf("  Phase offset: %f\n", total_phase_offset);

	/* If we reach an average offset of 0.00001 deg, quit the loop
	   (after coherently adding the phases one more time). */
	if(total_phase_offset <= 0.00001) {
	  itteration = itteration_max-1;
	}
      }
    }
    
    /* Final step: Store the phases of the template as real
       numbers. */
    for(binnr = 0; binnr < nrx; binnr++) {
      phase_track[binnr] = polar_angle_rad(phase_track_complex_template[2*binnr], phase_track_complex_template[2*binnr+1])*180.0/M_PI;
      phase_track[binnr] = derotate_deg(phase_track[binnr]);
    }

    if(phase_track_phases != NULL) {
      if(nspecbins > 1) {
	printwarning(verbose.debug, "  WARNING: More than one P3 modulation frequency channel is selected. The subpulse phase tracks for the individual frequency bins will be averaged, but not with an optimal weighting. You probably want to use shorter fft transforms to ensure that only one frequency bin is selected for the computation of the phase track.");
      }
      for(i = 0; i < nry/fft_size; i++) {
	phase_track_phases[i] = 0;
	for(j = 0; j < nspecbins; j++) {
	  phase_track_phases[i] += carg(phases[i*nspecbins+j])*180.0/(M_PI*(float)nspecbins); 
	}
//START REGION DEVELOP
	/*	printf("%f phases\n", phase_track_phases[i]);*/
//START REGION RELEASE
      }
    }
    free(phases);
  }  // End of phase track/amplidude if statement
  
  // Code related to -p3zap used earlier
  if(argc > 0 && argv != NULL) {
    for(i = 0; i < argc-1; i++) {
      if(strcmp(argv[i], "-p3zap") == 0) {
	j = sscanf(argv[i+1], "%f %f", &zapmin, &zapmax);
	if(j != 2) {
	  fflush(stdout);
	  printerror(verbose.debug, "Cannot parse -p3zap option.");
	  exit(0);
	}
	for(j = 0; j <= (fft_size/2); j++) {
	  if(zapmin > 0.9 || zapmax > 0.9) {
	    p3 = j;   // If the boundary exceeds 0.5 cpp, then the input parameters are assumed to be in bins
	  }else {
	    p3 = j/(float)(fft_size);   
	  }
	  if(p3 >= zapmin && p3 <= zapmax) {
	    for(k = 0; k < nrx; k++) {
	      lrfs[j*nrx+k] = 0;
	    }
	  }
	}
      }
    }
  }

  if(calcsubpulseAmplitude) {
    for(binnr = 0; binnr < nrx; binnr++)
      subpulseAmplitude[binnr] = 0;
    float peakpwr;
    peakpwr = 0;
    for(binnr = 0; binnr < nrx; binnr++) {
//START REGION DEVELOP
      /* I do not totally understand the normalisation when there are more than one phase tracks, but it is consistent with Russells software and it prevents getting more than 100% subpulse amplitude when power is offset from bin centre. */
//START REGION RELEASE
      subpulseAmplitude[binnr] = sqrt(phase_track_complex_template[2*binnr]*phase_track_complex_template[2*binnr]+phase_track_complex_template[2*binnr+1]*phase_track_complex_template[2*binnr+1])/sqrt(nspecbins);
      pwrtot = 0;
      for(i = 0; i < nr_fftblocks*fft_size; i++) {
	pwrtot += data[i*nrx+binnr];
      }
      if(binnr == 0 || pwrtot > peakpwr)
	peakpwr = pwrtot;
    }
    // Normalise with power of peak intensity pulse profile
    for(binnr = 0; binnr < nrx; binnr++) {
      subpulseAmplitude[binnr] /= peakpwr;
//START REGION DEVELOP
/* Russells magic factor, not entirely sure why. The subpulse amplitude in (if power is in one frequency bin) comes out to be 1, independent of bin length. */
//START REGION RELEASE
      subpulseAmplitude[binnr] *= 2;   
    }
  }

  if(verbose.verbose && verbose.nocounters == 0) 
    printf("Done                       \n");

  avr_power_subtracted_from_lrfs /= (double)(fft_size/2+1);
  if(avrg_offpulse_lrfs_power != NULL) {
    *avrg_offpulse_lrfs_power = avr_power_subtracted_from_lrfs;
  }
  if(verbose.debug) {
    if(regions != NULL) {
      if(regions->nrRegions > 0) {
	printf("  Average (off-pulse) power subtracted from LRFS per frequency/pulse longitude bin: %e\n", avr_power_subtracted_from_lrfs);
      }
    }
  }


  free(data1);
  free(lrfs_tmp);
  if(calcPhaseTrack || calcsubpulseAmplitude) {
    free(phase_track_complex);
    free(phase_track_complex_template);
  }
#ifdef USEFFTW3
  fftwf_destroy_plan(plan1);
  fftwf_destroy_plan(plan2);
#endif


  return 1;
}
//START REGION DEVELOP

// This function is the same as calcLRFS, but it takes as an extra argument the number of subintegrations to add to the end of the dataset filled with zeros. The only extra parameter is nr_zero_pad_subints_to_add. So this is a wrapper function that does zero-padding.
int calcLRFS_padding(float *data, long nry, long nrx, unsigned long fft_size, float *lrfs, int subtractDC, float *avrg_offpulse_lrfs_power, float *phase_track, float *phase_track_phases, int calcPhaseTrack, float freq_min, float freq_max, int track_only_first_region, float *subpulseAmplitude, int calcsubpulseAmplitude, int mask_freqs, int inverseFFT, pulselongitude_regions_definition *regions, float *var_rms, int argc, char **argv, verbose_definition verbose, int nr_zero_pad_subints_to_add)
{
  int ret;
  long sample;
  float *padded_data;
  padded_data = malloc(nrx*(nry+nr_zero_pad_subints_to_add)*sizeof(float));
  if(padded_data == NULL) {
    printerror(verbose.debug, "ERROR calcLRFS_padding: Cannot allocate memory");
    return 0;
  }
  memcpy(padded_data, data, nrx*nry*sizeof(float));
  for(sample = nry*nrx; sample < nrx*(nry+nr_zero_pad_subints_to_add); sample++) {
    padded_data[sample] = 0.0;
  }
  ret = calcLRFS(padded_data, nry+nr_zero_pad_subints_to_add, nrx, fft_size, lrfs, subtractDC, avrg_offpulse_lrfs_power, phase_track, phase_track_phases, calcPhaseTrack, freq_min, freq_max, track_only_first_region, subpulseAmplitude, calcsubpulseAmplitude, mask_freqs, inverseFFT, regions, var_rms, argc, argv, verbose);
  free(padded_data);
  return ret;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Calculate std deviation (sigma) and the modulation index. The
   profile will be normalized in the process. If regions is set it
   will calculate error-bars as well (an array for rms_sigma and
   rms_modind). Calculate the profile for a integer number of
   fft-blocks and set nrpulses to a multiple of fft_size. You can pass
   on var_rms from calcLRFS to obtain a higher precision (and smaller
   numbers) in the error-bar calculation if there are baseline
   variations. Not sure if it is better, but it is more consistent
   with Russell's software. Not sure if it is really working well. In
   fact, I don't believe the errorbars are correct in Russells
   software either (I simulated a 2D sinusoid convolved with a top-hat
   function as a test). You should do bootstrapping to obtain reliable
   errorbars. If *avrg_offpulse_lrfs_power is != NULL, then this
   value as provided by calcLRFS is used to report the offpulse
   normalised variance as spectrally determined.
*/
void calcModindex(float *lrfs, float *profile, long nrx, unsigned long fft_size, unsigned long nrpulses, float *sigma, float *rms_sigma, float *modind, float *rms_modind, pulselongitude_regions_definition *regions, float var_rms, float *avrg_offpulse_lrfs_power, verbose_definition verbose)
{
  long b, f, n, nrblocks;
  float var, max, rms, av_sigma, rms_var;

  /* Calculate variance profile */
  for(b = 0; b < nrx; b++) {
    var = 0;
    for(f = 0; f <= fft_size/2; f++) {
      var += lrfs[f*nrx+b];
    }
    sigma[b] = var;
  }
  nrblocks = nrpulses / (fft_size);
  /* Find maximum of profile. */
  max = profile[0];
  for(b = 0; b < nrx; b++) {
    if(max < profile[b])
      max = profile[b];
  }

  if(verbose.debug) {
    printf("Modulation index calculation: profile normalisation = %e\n", max);
  }
  if(verbose.verbose) {
    if(avrg_offpulse_lrfs_power != NULL) {
      float offpulse_stddev;
      int nrspecchannels = 1+fft_size/2;
      offpulse_stddev = (*avrg_offpulse_lrfs_power)*nrspecchannels;    // Offpulse lrfs power after summing freqs
      offpulse_stddev *= nrpulses*2.0/(float)fft_size;                 // FFT dependent normalisation factor to get variance
      offpulse_stddev = sqrt(offpulse_stddev);                         // Get std dev
      offpulse_stddev /= max;                                          // Normalisation
      printf("Spectrally determined normalised off-pulse standard deviation = %e\n", offpulse_stddev);      
    }
  }

  /* Normalize profile and variance profile */
  for(b = 0; b < nrx; b++) {
    profile[b] /= max;
    sigma[b] /= max*max;
//START REGION DEVELOP
    /* Not sure I understand this scaling. Formula is basically
       Parseval's theorem: integral fft^2 = integral f^2. So sigma^2 =
       Sum(y^2)/N-y_av = Sum(LRFS)-y_av. Because DC component is taken
       out in LRFS it shouldn't be done again. The 2/fft_size probably
       comes from the way the fft is calculated. NR says that if you
       do the inverse transform you have to multiply by 2/n,
       suggesting that something happens with the normalization. Note
       that by taking the square root of the variance profile we
       obtained the standard deviation profile. */
//START REGION RELEASE
    sigma[b] = sqrt(fabs(sigma[b])*nrpulses*2.0/(float)fft_size);
  }
  rms = 0;
  rms_var = 0;
  av_sigma = 0;
  n = 0;
  if(regions != NULL) {
    /* Calculate average sigma, which has to be subtracted to get a
       rms, because sigma is per definition a positive number */
    if(regions->nrRegions > 0) {
      for(b = 0; b < nrx; b++) {
	if(checkRegions(b, regions, 0, verbose) == 0) {
	  rms_var += (sigma[b])*(sigma[b])*(sigma[b])*(sigma[b]);
	  rms += profile[b]*profile[b];
	  av_sigma += sigma[b]*sigma[b];  /* Really a variance */
	  n++;
	}
      }
      if(n > 0) {
	av_sigma /= (float)n;
	rms_var /= (float)n;
	rms /= (float)n;
//START REGION DEVELOP
	/* Hmmm, average is non-zero, so maybe the non-zeroness of the
	   offpulse variance should add to the errorbar????? Also,
	   Russells code seem to produce different error-bars for
	   different bins. Not sure how that works. */
	/*	*rms_sigma = sqrt(*rms_sigma-av_sigma*av_sigma);  */
	/*	rms_var = sqrt(rms_var-av_sigma*av_sigma);  */

	/* Do the error propagation from sqrt(x), hence, assume that
	   the errors on the variance is independent of pulse
	   longitude, making the error of the standard deviation a
	   function of pulse longitude [sigma_sqrt_x/sqrt(x) = 0.5*sigma_x/x] */
//START REGION RELEASE
	for(b = 0; b < nrx; b++) {
//START REGION DEVELOP
	  /* This is more in line with what Russell does, which underestimates errorbars??? */
//START REGION RELEASE
	  rms_sigma[b] = 0.5*sqrt(rms_var-av_sigma*av_sigma)/sigma[b];          
//START REGION DEVELOP
	  /*	  rms_sigma[b] = sqrt(sqrt(rms_var-av_sigma*av_sigma));	  */   /* This is what I used to do, which overestimates errorbars???? */
//START REGION RELEASE
	}

	rms = sqrt(rms);
      }
    }

    if(var_rms > 0) {
      var_rms /= max*max;
      /* var_rms is per point in a single lrfs block, so calculate the overall rms */
      var_rms *= sqrt(n*(fft_size/2)*nrblocks);
      *rms_sigma = sqrt(fabs(var_rms*var_rms)*nrpulses*2.0/(float)fft_size);
    }
  }


  for(b = 0; b < nrx; b++) {
    if(profile[b] != 0)
      modind[b] = sigma[b]/profile[b];
    else
      modind[b] = -1e10;
    if(regions != NULL) {
      if(regions->nrRegions > 0) {
//START REGION DEVELOP
	/* While checking out error propagation for Sarah, I noticed
	   that error propagation was not done properly. If rms is
	   really what I think it is, you should add the terms is
	   quadrature? */
	/*		rms_modind[b] = fabs(modind[b]*(*rms_sigma/sigma[b] + rms/profile[b])); */
	/*	rms_modind[b] = fabs(modind[b]*modind[b]*(pow(*rms_sigma/sigma[b], 2.0) + pow(rms/profile[b], 2.0))); */
//START REGION RELEASE
	rms_modind[b] = sqrt(fabs(modind[b]*modind[b]*(pow(*rms_sigma/sigma[b], 2.0) + pow(rms/profile[b], 2.0)))); 
      }
    }
  }
}

//START REGION DEVELOP


/* Calculate std deviation (sigma) and the modulation index directly
   from the pulse stack. No errors are calculated and no baseline is
   subtracted. If regions is defined the variance of the off-pulse
   region is subtracted before calculating the modulation index. */
void calcModindexSimple(float *pulsestack, long nrx, unsigned long nrpulses, float *sigma, float *modind, pulselongitude_regions_definition *regions, verbose_definition verbose)
{
  long b, p, n;
  float var, intensity, av_sigma;

  for(b = 0; b < nrx; b++) {
    var = 0;
    intensity = 0;
    for(p = 0; p < nrpulses; p++)
      intensity += pulsestack[p*nrx+b];
    intensity /= (float)nrpulses;
    for(p = 0; p < nrpulses; p++) 
      var += (pulsestack[p*nrx+b]-intensity)*(pulsestack[p*nrx+b]-intensity);
    var /= (float)nrpulses;
    sigma[b] = var;
  }


  av_sigma = 0;
  if(regions != NULL) {
    /* Calculate average sigma, which has to be subtracted to get a
       rms, because sigma is per definition a positive number */
    if(regions->nrRegions > 0) {
      n = 0;
      av_sigma = 0;
      for(b = 0; b < nrx; b++) {
	if(checkRegions(b, regions, 0, verbose) == 0) {
	  av_sigma += sigma[b];  /* Really a variance */
	  n++;
	}
      }
      if(n > 0) {
	av_sigma /= (float)n;
      }
    }
  }

  /*  av_sigma = 0; */

  for(b = 0; b < nrx; b++) {
    intensity = 0;
    for(p = 0; p < nrpulses; p++)
      intensity += pulsestack[p*nrx+b];
    intensity /= (float)nrpulses;
    sigma[b] -= av_sigma;
    sigma[b] = sqrt(fabs(sigma[b]))*sigma[b]/fabs(sigma[b]);
    modind[b] = sigma[b]/intensity;
  }
}

//START REGION DEVELOP



/* 
   Calculate the longitude resolved cross-correlation map of the data
   for the specified lag. The map is written out to lrcc, which should
   be nrx*nrx floats in size. If noSubtractMean is specified, the mean
   profiles is not subtracted from the LRCC. 

   lrcc_ij = sqrt(sum_i sum_j sum_t stack_i(t)*stack_j(t+lag))

   0 = nothing done 
   1 = ok
*/
int calcLRCC(float *data, long nry, long nrx, int lag, float *lrcc, int noSubtractMean, verbose_definition verbose)
{
  float *dataMean;
  long i, j, n, nstart, nend;

  if(verbose.verbose) printf("Calculating LRCC with a lag of %d\n", lag);

  /* Calculate mean intensity of each bin if required. */
  if(noSubtractMean == 0) {
    dataMean = (float *)malloc(nrx*sizeof(float));
    if(dataMean == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calcLRCC: Cannot allocate memory");
      return 0;
    }
    if(verbose.verbose) printf("  Calculating mean\n");
    for(i = 0; i < nrx; i++)
      dataMean[i] = 0;
    for(n = 0; n < nry; n++) {
      for(i = 0; i < nrx; i++) {
	dataMean[i] += data[n*nrx+i];
      }
    }
    for(i = 0; i < nrx; i++) {
      dataMean[i] /= (float)nry;
    }
  }
  
  /* Empty the lrcc, so we can add the correlations later. */
  for(i = 0; i < nrx*nrx; i++)
    lrcc[i] = 0;

  if(lag >= 0) {
    nstart = 0;
    nend = nry-lag;
  }else {
    nstart = -lag;
    nend = nry;
  }
  for(n = nstart; n < nend; n++) {
    for(i = 0; i < nrx; i++) {
      for(j = 0; j < nrx; j++) {
	if(noSubtractMean) {
	  lrcc[i*nrx+j] += data[n*nrx+i]*data[(n+lag)*nrx+j];
	}else {
	  lrcc[i*nrx+j] += (data[n*nrx+i]-dataMean[i])*(data[(n+lag)*nrx+j]-dataMean[j]);
	}
      }
    }
    if(verbose.verbose && verbose.nocounters == 0) {
      printf("  Pulse %ld of the %ld     \r", (n-nstart)+1, nry);
      fflush(stdout);
    }
  }
  /* Take the square root to make the scale more like intensity, rather than intensity squared. */
  for(i = 0; i < nrx; i++) {
    for(j = 0; j < nrx; j++) {
      if(lrcc[i*nrx+j] >= 0) {
	lrcc[i*nrx+j] = sqrt(lrcc[i*nrx+j]);
      }else {
	lrcc[i*nrx+j] = -sqrt(-lrcc[i*nrx+j]);
      }
    }
  }



  if(verbose.verbose && verbose.nocounters == 0) 
    printf("Done                       \n");
  if(noSubtractMean == 0) 
    free(dataMean);
  return 1;
}

//START REGION DEVELOP
/* 
   Calculate the HRFS of data, which is an array of nrx by nry
   points. The data is analysed per block of pulses_per_block pulses and the
   resulting power spectra are added. The resulting spectrum is
   hrfs_length points in length (memory will be allocated). If regions
   is defined (on-pulse), than the off-pulse bins are set to zero
   before doing calculating the hrfs.

   0 = nothing done 
   1 = ok
*/
int calcHRFS(float *data, long nry, long nrx, long pulses_per_block, float **hrfs, long *hrfs_length, pulselongitude_regions_definition *regions, verbose_definition verbose)
{
  unsigned long nr_fftblocks;
  long i, fft_size, fftblock, binnr, pulsenr;
  int add_extra_zero;
  float *datablock, pwr;
#ifdef USEFFTW3
  fftwf_plan fftplan;
#endif

  // Calculate nr of spectra we're going to average.
  nr_fftblocks = nry/pulses_per_block;
  if(verbose.verbose) printf("Calculating HRFS (%ld blocks) using %ld/%ld pulses\n", nr_fftblocks, nr_fftblocks*pulses_per_block, nry);
  if(nr_fftblocks == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calcHRFS: Cannot calculate hrfs for %ld pulses (smaller than specified block size = %ld)", nry, pulses_per_block);
    return 0;
  }

  // Calculated size of final spectrum and make sure it is even.
  fft_size = nrx*pulses_per_block;
  i = fft_size / 2;
  i *= 2;
  if(i != fft_size) {
    if(verbose.verbose) {
      fflush(stdout);
      printf("  An extra zero is added to each block of %ld points to make it even\n", fft_size);
    }
    add_extra_zero = 1;
    fft_size++;
  }else {
    add_extra_zero = 0;
  }


  if(regions != NULL) {
    if(regions->bins_defined[0] == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calcHRFS: Region is not defined in bins");
      //    printRegions(onpulse);
      return 0;
    }
  }

  // Allocate some temporary memory. datablock will contain a block of pulses from
  // the pulse-stack which will go into
  // the fft routine. hrfs will contain the result (added power spectra) of each block
  // of data. The extra 2 values is because there will be
  // (fft_size/2+1) complex valued spectral bins. 
  datablock = (float *)malloc((fft_size+2)*sizeof(float));
  *hrfs = (float *)malloc((fft_size/2+1)*sizeof(float));
  *hrfs_length = (fft_size/2+1);
  if(datablock == NULL || *hrfs == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calcHRFS: Cannot allocate memory");
    return 0;
  }
  if(verbose.verbose) printf("  output is a fft of length %ld\n", *hrfs_length);

  // Empty the hrfs, so we can add the spectra of the individual
  // blocks later.
  for(i = 0; i < (1+fft_size/2); i++)
    (*hrfs)[i] = 0;

#ifdef USEFFTW3
  // Do an inplace transform, like realft of NR.
  fftplan = fftwf_plan_dft_r2c_1d(fft_size, datablock, (fftwf_complex *)datablock, FFTW_ESTIMATE);
  // FFTW routine takes integer as fft length 
  if(fft_size > 2147483640) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calcHRFS: requested fft too long.");
    return 0;
  }
#endif

  // Loop over the blocks of pulses for which we're calculating the spectra 
  for(fftblock = 0; fftblock < nr_fftblocks; fftblock++) {
    // Extract a block of data
    for(pulsenr = 0; pulsenr < pulses_per_block; pulsenr++) {
      for(binnr = 0; binnr < nrx; binnr++) {
	datablock[pulsenr*nrx+binnr] = data[(fftblock*pulses_per_block+pulsenr)*nrx + binnr];
	if(regions != NULL) {
	  if(checkRegions(binnr, regions, 0, verbose) == 0) {
	    datablock[pulsenr*nrx+binnr] = 0;
	  }
	}
      }
    }
    // Input array is not an even number, so add an extra zero to make it even
    if(add_extra_zero) {
      datablock[pulses_per_block*nrx] = 0;
    }

    // Calculate fft of the (real valued) data resulting in a complex spectrum 
#ifdef USEFFTW3
    fftwf_execute(fftplan);
#else
    realft(datablock-1, fft_size, 1);
    // Unpack the (real valued) Nyquest frequency, which is stored
    // as the imaginary part of the real valued DC bin. 
    datablock[2*(fft_size/2)]   = datablock[1];
    datablock[2*(fft_size/2)+1] = 0;
    datablock[1] = 0;
#endif

    // Loop over spectral bins 
    for(i = 0; i <= fft_size/2; i++) {
      // Calculate the amplitude squared of the complex fft. 
      pwr  = datablock[2*i]*datablock[2*i] + datablock[2*i+1]*datablock[2*i+1];
      // Add it to final spectrum
      (*hrfs)[i] += pwr;
    }

  } // End of fftblock loop


#ifdef USEFFTW3
  fftwf_destroy_plan(fftplan);
#endif

  free(datablock);

  return 1;
}
//START REGION DEVELOP



//START REGION DEVELOP
/* 
   Calculate the LRAC (longitude resolved auto-correlation) of
   data. If lrac_keepbothhalves is set, both positive and negative
   lags are stored, which is an array of nrx by cclength points (as
   returned by crosscorrelation_fft_padding_cclength() when adding at
   least nry extra zero-padding zero's). If lrac_keepbothhalves is
   zero, only zero lag + positive lags are stored (the return value of
   crosscorrelation_fft_padding_cclength()/2 points). This length is
   returned in cclength. The data will be zero-padded in the nry
   direction to make it a power of 2. If remove_wf is set, the
   triangular shaped effect because of the window function is removed
   from the AC. If remove_zerolag the zero lag coefficients (often
   respondsible for a large spike) is set to zero. The debug option
   gives more output.

   0 = error
   1 = ok
*/
int calcLRAC(float *data, long nry, long nrx, float *lrac, int remove_wf, int remove_spike, int lrac_keepbothhalves, long *cclength, verbose_definition verbose)
{
  long i, binnr, pulsenr, extrazeropad;
  int cclength2;
  float *cc_tmp, *data1;
  double av, central_spike;
  verbose_definition verbose2;

  copyVerboseState(verbose, &verbose2);
  verbose2.verbose = verbose.debug;

  if(verbose.debug) printf("LRAC: Pad with at least %ld zero's to avoid any wrap-effects\n", nry);
  extrazeropad = nry;

  *cclength = crosscorrelation_fft_padding_cclength(nry, extrazeropad);
  if(lrac_keepbothhalves == 0) {
    *cclength = (*cclength)/2;
  }
  if(verbose.verbose) printf("Calculating LRAC (%ld X %ld points)\n", nrx, *cclength);

  // Enough memory to hold the cross-correlation function of one bin, so the output can be transposed
  data1 = (float *)malloc(nry*sizeof(float));
  if(data1 == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR calcLRAC: Cannot allocate memory");
    return 0;
  }

  for(binnr = 0; binnr < nrx; binnr++) {
    /* Get column out of pulse-stack */
    central_spike = 0;
    av = 0;
    for(pulsenr = 0; pulsenr < nry; pulsenr++) {
      data1[pulsenr] = data[pulsenr*nrx + binnr];
      central_spike += data1[pulsenr]*data1[pulsenr];
      av += data1[pulsenr];
    }
    av /= (double)nry;
    // Calculate the auto-correlation
    if(crosscorrelation_fft_padding(data1, data1, nry, extrazeropad, &cc_tmp, &cclength2, verbose2) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR calcLRAC: Error computing the auto-correlation function.\n");
      return 0;
    }
    if( (lrac_keepbothhalves && (*cclength != cclength2)) || (lrac_keepbothhalves == 0 && (*cclength != cclength2/2)) ){
      fflush(stdout);
      printerror(verbose.debug, "ERROR calcLRAC: There is an unexpected mismatch between the expected and actual length of the auto-correlation output.\n");
      return 0;
    }

    // Unpack the cc function so it starts with the most negative lag, rather than lag 0.
    if(lrac_keepbothhalves) {
      for(i = cclength2/2; i <= cclength2-1; i++) {
	//Note: lag number = i-cclength2
	if(remove_wf) {
	  if(i - cclength2 > -(cclength2-nry)) {
	    if(i-cclength2 > -nry) {     // If lag is larger than length data, the cc should be zero without any zero's
	      cc_tmp[i] -= av*av*(i-cclength2); // Remove slope
	      cc_tmp[i] -= nry*av*av; // Remove offset of slope-bit of the triangle
	    }
	  }else {
	    cc_tmp[i] -= (2*nry-cclength2)*av*av;  // Remove flat bit after triangle
	  }
	}
	lrac[(i-cclength2/2)*nrx+binnr] = cc_tmp[i];
      }
    }
    for(i = 0; i < cclength2/2; i++) {
      //Note: lag number = i
      if(remove_wf) {
	if(i == 0) {
	  if(remove_spike) {
	    cc_tmp[i] -= central_spike;  // Remove central spike
	  }
	}else {
	  if(i < cclength2-nry) {
	    if(i < nry) {  // If lag is larger than length data, the cc should be zero without any zero's
	      cc_tmp[i] += av*av*i; // Remove slope
	      cc_tmp[i] -= nry*av*av;  // Remove offset of slope-bit of the triangle
	    }
	  }else {
	    cc_tmp[i] -= (2*nry-cclength2)*av*av; // Remove flat bit after triangle
	  }
	}
      }
      if(lrac_keepbothhalves) {
	lrac[(i+cclength2/2)*nrx+binnr] = cc_tmp[i];
      }else {
	lrac[i*nrx+binnr] = cc_tmp[i];
      }
    }

    free(cc_tmp);
  }

  //  if(verbose.verbose && verbose.nocounters == 0) 
  //    printf("Done                       \n");

  free(data1);

  return 1;
}

//START REGION DEVELOP

typedef struct {
  float centre, width, s2n;  // location are in cpp
}p3classify_peak_info;

// padding_up_to: if the nr of pulses analysed is less than this number, the data will be zero-padded to this length. When set to zero, no zero-padding will be applied
// produce_plot_when_debug
//   0 = make no plots
//   1 = make pulse stack + initial paps plot
//   2 = add final paps plot 
// Return 0 on error, 1 on success
int p3classify_core(float *data, long nstart, long nend, long nrx, long nry, long padding_up_to, float cpp1, float cpp2, float *lrfs, float *paps, pulselongitude_regions_definition *regions, int argc, char **argv, p3classify_peak_info *peak, verbose_definition verbose, int produce_plot_when_debug)
{
  long i, j, ret;
  int ok;
  long fft_size, spec_bin, zero_pulses_to_add;
  float var_rms;
  pgplot_options_definition *pgplot_options;
  pgplot_options = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
  if(pgplot_options == NULL) {
    printerror(verbose.debug, "ERROR p3classify_core: Memory allocation error");
    return 0;
  }
  pgplot_clear_options(pgplot_options);


  int oldverbose = verbose.verbose;
  int olddebug = verbose.debug;
  verbose.verbose = 0;

  // Do the lrfs, potentially with zero-padding enables
  zero_pulses_to_add = padding_up_to - (nend - nstart + 1);
  if(zero_pulses_to_add <= 0) {
    fft_size = nend - nstart + 1;  // Take the full pulse-range as fft-size, so the data is not split in separate blocks ensuring that the full range of pulses is used. Note that fftw3 allows dft's of any size.
    if(calcLRFS(&data[nstart*nrx], fft_size, nrx, fft_size, lrfs, 1, NULL, NULL, NULL, 0, 0.0, 0.0, 0, NULL, 0, 0, 0, regions, &var_rms, argc, argv, verbose) == 0) {
      printerror(verbose.debug, "ERROR p3classify_core: Cannot compute lrfs.");
      free(pgplot_options);
      return 0;
    }
  }else {
    //    printf("DO ZERO PADDING: %ld\n", zero_pulses_to_add);

    fft_size = padding_up_to;  // Take the full pulse-range as fft-size, including the added padding.
    if(calcLRFS_padding(&data[nstart*nrx], nend - nstart + 1, nrx, fft_size, lrfs, 1, NULL, NULL, NULL, 0, 0.0, 0.0, 0, NULL, 0, 0, 0, regions, &var_rms, argc, argv, verbose, zero_pulses_to_add) == 0) {
      printerror(verbose.debug, "ERROR p3classify_core: Cannot compute lrfs.");
      free(pgplot_options);
      return 0;
    }
  }
  verbose.verbose = oldverbose;

  // Integrate the lrfs within the onpulse region to get the fluctuation power as function of frequency
  for(spec_bin = 0; spec_bin < fft_size/2+1; spec_bin++) {
    paps[spec_bin] = 0;
  }
  for(i = 0; i < nrx; i++) {
    ok = 1;
    if(regions != NULL) {
      if(regions->nrRegions > 0) {
	if(regions->bins_defined[0]) {
	  if(checkRegions(i, regions, 0, verbose) == 0) {
	    ok = 0;
	  }
	}
      }
    }
    if(ok) {
      for(spec_bin = 0; spec_bin < fft_size/2+1; spec_bin++) {
	paps[spec_bin] += lrfs[spec_bin*nrx + i];
      }
    }
  }

  pulselongitude_regions_definition cpp_range;
  if(initPulselongitudeRegion(&cpp_range, verbose) == 0) {
    printerror(verbose.debug, "ERROR p3classify_core: Cannot initialise onpulse region.");
    free(pgplot_options);
    return 0;
  }
  cpp_range.nrRegions = 1;
  cpp_range.bins_defined[0] = 1;
  cpp_range.left_bin[0]  = cpp1*(fft_size/2+1)/0.5; // I'm not 100% about this conversion. I assumed this to be good enough here and in the following.
  cpp_range.right_bin[0] = cpp2*(fft_size/2+1)/0.5;
  

  // Find and subtract baseline
  double baseline;
  long nr_used_bins;
  float zapmin, zapmax;
  int loop;
  baseline = 0;
  nr_used_bins = 0;
  for(loop = 0; loop < 2; loop++) { // First loop: determine the baseline. Second loop: subtract it.
    for(i = 0; i < fft_size/2+1; i++) {
      ok = 1;
      if(loop == 0) {
	if(cpp_range.nrRegions > 0) {
	  if(cpp_range.bins_defined[0]) {
	    if(checkRegions(i, &cpp_range, 0, verbose) != 0) {
	      ok = 0;
	    }
	  }
	}
      }

      for(j = 1; j < argc-1; j++) {
	if(strcmp(argv[j], "-p3zap") == 0) {
	  ret = sscanf(argv[j+1], "%f %f", &zapmin, &zapmax);
	  if(ret != 2) {
	    printerror(verbose.debug, "ERROR p3classify_core: Cannot parse -p3zap option, specify two values.");
	    free(pgplot_options);
	    return 0;
	  }
	  if(zapmin > 0.9 || zapmax > 0.9) {  // If the boundary exceeds 0.5 cpp, then the input parameters are assumed to be in bins
	    printerror(verbose.debug, "ERROR p3classify_core: The -p3zap option should be specified in cpp.");
	    free(pgplot_options);
	    return 0;
	  }
	  // Convert zapmin/zapmax to bin numbers
	  zapmin = (fft_size/2+1)*zapmin/0.5;
	  zapmax = (fft_size/2+1)*zapmax/0.5;
	  if(i >= zapmin && i <= zapmax) {
	    ok = 0;
	  }
	}
      }
      if(loop == 0) {
	if(ok) {
	  baseline += paps[i];
	  nr_used_bins++;
	}
      }else {
	if(ok) {
	  paps[i] -= baseline;
	}
      }
    }
    if(nr_used_bins == 0) {
      printerror(verbose.debug, "ERROR p3classify_core: There are no spectral bins that can be used to determine the baseline.");
      free(pgplot_options);
      return 0;
    }
    if(loop == 0) {
      baseline /= (double)nr_used_bins;
    }
  }

  int peak_start, peak_width;
  float peak_snr, peak_integral;
  verbose.verbose = 0;
  if(boxcarFindpeak(paps, fft_size/2+1, &cpp_range, &peak_start, &peak_width, &peak_snr, &peak_integral, 0, 0, 0, 1, cpp_range.right_bin[0] - cpp_range.left_bin[0], 1, 1, verbose) == 0) {
    printerror(verbose.debug, "ERROR p3classify_core: Peak searching failed.");
    free(pgplot_options);
    return 0;
  }
  verbose.verbose = oldverbose;

  peak->centre = 0.5*(peak_start+0.5*peak_width)/(fft_size/2+1);
  peak->width = 0.5*peak_width/(fft_size/2+1);
  peak->s2n = peak_snr;

  if(verbose.debug && produce_plot_when_debug) {

    if(produce_plot_when_debug == 1) { // Make pulse stack plot
      long first_bin_to_plot, last_bin_to_plot;
      
      // Find out pulse longitude region in the onpulse region, if defined
      first_bin_to_plot = nrx;
      last_bin_to_plot = -1;
      if(regions != NULL) {
	if(regions->nrRegions > 0) {
	  if(regions->bins_defined[0]) {
	    
	    for(i = nrx-1; i >= 0; i--) {
	      if(checkRegions(i, regions, 0, verbose) != 0) {
		first_bin_to_plot = i;
	      }
	    }
	    for(i = 0; i < nrx; i++) {
	      if(checkRegions(i, regions, 0, verbose) != 0) {
		last_bin_to_plot = i;
	      }
	    }
	  }
	}
      }
      if(first_bin_to_plot == nrx) {
	first_bin_to_plot = 0;
      }
      if(last_bin_to_plot == -1) {
	last_bin_to_plot = nrx-1;
      }
      
      pgplot_options->viewport.xsize = 0.3;
      pgplot_options->viewport.ysize = 1.2;
      pgplot_options->viewport.dyplot = -0.1;
      pgplot_options->viewport.dontclose = 1;
      
      long nrpulsestoplot;

      nrpulsestoplot = 3*(nend - nstart + 1);
      if(nstart + nrpulsestoplot > nry)
	nrpulsestoplot = nry - nstart;

      verbose.debug = 0;
      pgplotMap(pgplot_options, &data[nstart*nrx], nrx, nrpulsestoplot, 0, nrx-1, first_bin_to_plot, last_bin_to_plot, nstart, nstart+nrpulsestoplot-1, nstart, nstart+nrpulsestoplot-1, PPGPLOT_HEAT, 0, 0, 0, NULL, 0, 0, 1.0, 0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, verbose);
      verbose.debug = olddebug;
    }

    if(produce_plot_when_debug == 1) { // Make first paps plot
      pgplot_options->viewport.dontopen = 1;
      pgplot_options->viewport.noclear = 1;
      pgplot_options->viewport.dxplot = 0.35;
      pgplot_options->viewport.ysize = 0.2;
      pgplot_options->viewport.xsize = 0.6;
      pgplot_options->viewport.dyplot = 0.72;
      pgplot_options->viewport.dontclose = 1;
      pgplot_options->box.box_labelsize = 0.7;
      pgplotGraph1(pgplot_options, paps, NULL, NULL, fft_size/2+1, 0.0, 0.5, 0, 0.0, 0.5, 0.0, 0.0, 0, 1, 0, 1, 0, 1, 1, NULL, -1, verbose);
    
      printf("  Peak location: %f cpp (%f cpp wide = %.1f bins) with a S/N=%f\n", peak->centre, peak->width, (fft_size/2+1)*peak->width/0.5, peak->s2n);
    }else {
      pgplot_options->viewport.dontopen = 1;
      pgplot_options->viewport.noclear = 1;
      pgplot_options->viewport.dxplot = 0.35;
      pgplot_options->viewport.ysize = 0.2;
      pgplot_options->viewport.xsize = 0.6;
      pgplot_options->viewport.dyplot = 0.50;
      pgplot_options->viewport.dontclose = 0;
      pgplot_options->box.box_labelsize = 0.7;
      pgplotGraph1(pgplot_options, paps, NULL, NULL, fft_size/2+1, 0.0, 0.5, 0, 0.0, 0.5, 0.0, 0.0, 0, 1, 0, 1, 0, 1, 1, NULL, -1, verbose);
    
      printf("  Improved peak location: %f cpp (%f cpp wide = %.1f bins) with a S/N=%f\n", peak->centre, peak->width, (fft_size/2+1)*peak->width/0.5, peak->s2n);
    }
  }
  freePulselongitudeRegion(&cpp_range);
  free(pgplot_options);
  return 1;
}

// plottype
//     0 = No plot generated
//     1 = Generate plot
//     2 = Add 2nd paps
// Return 0 = error
int p3classify_slide(float *data, long nrx, long nry, long padding_up_to, long nstart, long nend, int first_itteration, int plottype, float cpp1, float cpp2, pulselongitude_regions_definition *regions, int argc, char **argv, float *orig_s2n, float *cur_s2n, float *max_s2n_found, long *nstart_when_s2n_is_max, long *nend_when_s2n_is_max, float *lrfs, float *paps, verbose_definition verbose)
{
  int produce_plot_when_debug;
  p3classify_peak_info peak;

  if(first_itteration) {
    //    if(plottype)
    produce_plot_when_debug = plottype;
    *nstart_when_s2n_is_max = nstart;
    *nend_when_s2n_is_max = nend;
    *max_s2n_found = -1;
    printf("  Initial guess: nstart=%ld nend=%ld\n", nstart, nend);
  }else {
    produce_plot_when_debug = 0;
  }
  if(p3classify_core(data, nstart, nend, nrx, nry, padding_up_to, cpp1, cpp2, lrfs, paps, regions, argc, argv, &peak, verbose, produce_plot_when_debug) == 0) {
    printerror(verbose.debug, "ERROR p3classify_slide: Classifying failed.");
    return 0;
  }
  *cur_s2n = peak.s2n; // Keep the s/n of this particular guess of the start/end point
  if(peak.s2n >= *max_s2n_found) {
    *max_s2n_found = peak.s2n;
    *nstart_when_s2n_is_max = nstart;
    *nend_when_s2n_is_max = nend;
  }
  if(first_itteration) {
    *orig_s2n = peak.s2n;
  }
  return 1;
}

/*
  If argc > 0 and argv != NULL,
   it will search the command line for the -p3zap option and zap these
   regions from the lrfs (can be specified in cpp or in bins). The 
   debug option gives more output.

 Return 0 on error
*/
int p3classify(float *data, long nry, long nrx, long padding_up_to, long nblock, long min_length, float cpp1, float cpp2, float *output, pulselongitude_regions_definition *regions, int argc, char **argv, verbose_definition verbose)
{
  float max_length_in_blocks = 3;  // The user specifies nblock, the initial guess of the mode duration. The mode could be longer, but not longer than this number of blocks.

  if(verbose.verbose) printf("Starting classifying P3 as function of time (initial block size=%ld, Npad=%ld)\n", nblock, padding_up_to);
#ifndef USEFFTW3
  // We need fftw3 because that will work for any size of the fft, unlike NR which requires a power of 2.
  printerror(verbose.debug, "ERROR p3classify: This function relies on fftw3, but its use was not enabled during compilation.");
  return 0;
#endif

  float *lrfs, *paps, *s2ngraph;
  lrfs = malloc((nry/2+1)*nrx*sizeof(float));   // Probably way too much memory of what is going to be used, but this should always work.
  paps = malloc((nry/2+1)*sizeof(float));  // Memory for the phase-averaged power spectrum (side panel of lrfs integrated over the on-pulse region)
  s2ngraph = malloc(max_length_in_blocks*nblock*sizeof(float));
  if(lrfs == NULL || paps == NULL || s2ngraph == NULL) {
    printerror(verbose.debug, "ERROR p3classify: Memory allocation error.");
    return 0;
  }

  long nstart_cur, nend_cur;
  long current_block;

  // First attempt: take a block with fixed size nblock
  nstart_cur = 0;
  current_block = 0;
  while(nstart_cur < nry) {

    nend_cur   = nstart_cur + nblock - 1;  // Calculate the current guess for the end of the block given currenst start pulse number.
    if(nend_cur >= nry) {
      nend_cur = nry - 1;
    }

    printf("Determining start/end pulse numbers for stretch number %ld\n", current_block+1);
    

    int first_itteration, plottype;
    long nstart, nend;
    
    
    float max_s2n_found, cur_s2n, orig_s2n, prev_max_s2n; long nstart_when_s2n_is_max, nend_when_s2n_is_max;
    
    for(nstart = 0; nstart < max_length_in_blocks*nblock; nstart++) {
      s2ngraph[nstart] = 0;
    }

    // Step 2: try if making nend larger improve things
    first_itteration = 1;
    nstart = nstart_cur;

    //nend_cur = 10;

    for(nend = nend_cur; nend < nry; nend++) {  // nend can never exceed the number of available subints
      if(nend - nstart + 1 > max_length_in_blocks*nblock) // There is a maximum mode duration to consider
	break;
      plottype = 1;
      if(p3classify_slide(data, nrx, nry, padding_up_to, nstart, nend, first_itteration, plottype, cpp1, cpp2, regions, argc, argv, &orig_s2n, &cur_s2n, &max_s2n_found, &nstart_when_s2n_is_max, &nend_when_s2n_is_max, lrfs, paps, verbose) == 0) {
	printerror(verbose.debug, "ERROR p3classify: Classifying failed.");
	free(lrfs);
	free(paps);
	free(s2ngraph);
	return 0;
      }
      s2ngraph[nend-nstart] = cur_s2n;
      //      if(verbose.debug)
	//	printf("   nstart = %ld nend = %ld s2n = %e\n", nstart, nend, cur_s2n);
      first_itteration = 0;
    }
    printf("  Increasing nend to %ld improved the s/n from %f to %f\n", nend_when_s2n_is_max, orig_s2n, max_s2n_found);
    nend_cur = nend_when_s2n_is_max;
    prev_max_s2n = max_s2n_found; // Remember the best S/N found when increasing the length of the mode

    
    // Step 3: try instead if making nend smaller than what we initially had improve things
    first_itteration = 1;
    nstart = nstart_cur;
    for(nend = nstart_cur + nblock - 1; nend >= nstart+min_length-1; nend--) {  // nend can never result in a mode size which is less than min_length pulses
      plottype = 0;
      if(p3classify_slide(data, nrx, nry, padding_up_to, nstart, nend, first_itteration, plottype, cpp1, cpp2, regions, argc, argv, &orig_s2n, &cur_s2n, &max_s2n_found, &nstart_when_s2n_is_max, &nend_when_s2n_is_max, lrfs, paps, verbose) == 0) {
	printerror(verbose.debug, "ERROR p3classify: Classifying failed.");
	free(lrfs);
	free(paps);
	free(s2ngraph);
	return 0;
      }
      s2ngraph[nend-nstart] = cur_s2n;
      first_itteration = 0;
    }
    printf("  Decreasing nend to %ld improved the s/n from %f to %f\n", nend_when_s2n_is_max, orig_s2n, max_s2n_found);
    if(verbose.debug) { // Make s/n evolution plot
      pgplot_options_definition *pgplot_options;
      pgplot_options = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
      if(pgplot_options == NULL) {
	printerror(verbose.debug, "ERROR p3classify: Memory allocation error");
	return 0;
      }
      pgplot_clear_options(pgplot_options);
      
      pgplot_options->viewport.dontopen = 1;
      pgplot_options->viewport.noclear = 1;
      pgplot_options->viewport.dxplot = 0.35;
      pgplot_options->viewport.ysize = 0.2;
      pgplot_options->viewport.xsize = 0.6;
      pgplot_options->viewport.dyplot = 0.25;
      pgplot_options->viewport.dontclose = 1;
      sprintf(pgplot_options->box.xlabel, "nend");
      sprintf(pgplot_options->box.ylabel, "spectral S/N");
      pgplot_options->box.box_labelsize = 0.7;
      pgplotGraph1(pgplot_options, s2ngraph, NULL, NULL, max_length_in_blocks*nblock, nstart, nstart+floor(max_length_in_blocks*nblock)-1, 0, nstart+min_length-1, nstart+floor(max_length_in_blocks*nblock)-1, 0.0, 0.0, 0, 1, 0, 1, 0, 1, 1, NULL, -1, verbose);
      free(pgplot_options);
      //      printf("Press a key to continue\n");
      //      pgetch();
    }
    if(max_s2n_found > prev_max_s2n) {
      printf("    Changing nend to %ld\n", nend_when_s2n_is_max);
      nend_cur = nend_when_s2n_is_max;
    }else {
      printf("    Keep nend to be %ld\n", nend_cur);
    }
    
    // Step 4: try if making nstart larger than what we initially had improve things
    first_itteration = 1;
    nend = nend_cur;
    for(nstart = 0; nstart < max_length_in_blocks*nblock; nstart++) {
      s2ngraph[nstart] = 0;
    }
    for(nstart = nstart_cur; nstart <= nend-min_length+1; nstart++) {  // nstart can never result in a mode size which is less than min_length pulses
      plottype = 0;
      if(p3classify_slide(data, nrx, nry, padding_up_to, nstart, nend, first_itteration, plottype, cpp1, cpp2, regions, argc, argv, &orig_s2n, &cur_s2n, &max_s2n_found, &nstart_when_s2n_is_max, &nend_when_s2n_is_max, lrfs, paps, verbose) == 0) {
	printerror(verbose.debug, "ERROR p3classify: Classifying failed.");
	free(lrfs);
	free(paps);
	free(s2ngraph);
	return 0;
      }
      s2ngraph[nstart-nstart_cur] = cur_s2n;
      first_itteration = 0;
    }
    printf("  Increasing nstart to %ld improved the s/n from %f to %f\n", nstart_when_s2n_is_max, orig_s2n, max_s2n_found);
    if(verbose.debug) { // Make s/n evolution plot
      pgplot_options_definition *pgplot_options;
      pgplot_options = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
      if(pgplot_options == NULL) {
	printerror(verbose.debug, "ERROR p3classify: Memory allocation error");
	return 0;
      }
      pgplot_clear_options(pgplot_options);
      pgplot_options->viewport.dontopen = 1;
      pgplot_options->viewport.noclear = 1;
      pgplot_options->viewport.dxplot = 0.35;
      pgplot_options->viewport.ysize = 0.2;
      pgplot_options->viewport.xsize = 0.6;
      pgplot_options->viewport.dyplot = 0.0;
      pgplot_options->viewport.dontclose = 1;
      sprintf(pgplot_options->box.xlabel, "nstart");
      sprintf(pgplot_options->box.ylabel, "spectral S/N");
      pgplot_options->box.box_labelsize = 0.7;
      pgplotGraph1(pgplot_options, s2ngraph, NULL, NULL, nend-min_length+1-nstart_cur+1, nstart_cur, nend-min_length+1, 0, nstart_cur, nend-min_length+1, 0.0, 0.0, 0, 1, 0, 1, 0, 1, 1, NULL, -1, verbose);
      free(pgplot_options);
      //      printf("Press a key to continue\n");
      //      pgetch();
    }
    nstart_cur = nstart_when_s2n_is_max;
    
    
    if(verbose.debug) {
      printf("Showing spectrum of the final result\n");
      plottype = 2;
      first_itteration = 1;
      nstart = nstart_cur;
      nend = nend_cur;
      if(p3classify_slide(data, nrx, nry, padding_up_to, nstart, nend, first_itteration, plottype, cpp1, cpp2, regions, argc, argv, &orig_s2n, &cur_s2n, &max_s2n_found, &nstart_when_s2n_is_max, &nend_when_s2n_is_max, lrfs, paps, verbose) == 0) {
	printerror(verbose.debug, "ERROR p3classify: Classifying failed.");
	free(lrfs);
	free(paps);
	free(s2ngraph);
	return 0;
      }
    }

    if(verbose.verbose) printf("  Found a stretch between %ld and %ld with s/n %f\n", nstart_cur, nend_cur, max_s2n_found);


    nstart_cur = nend_cur+1;
    current_block++;
  }  // End of while loop

  free(lrfs);
  free(paps);
  free(s2ngraph);
  return 1;
}

