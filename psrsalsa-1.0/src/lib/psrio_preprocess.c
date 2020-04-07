//START REGION RELEASE
#include <math.h>
#include <string.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_fit.h"
#include "psrsalsa.h"

/* Create a rebinned clone. The clone should not exist yet. Return 1 =
   success. Verbose level determines nr of spaces before output. If
   nocounters is not set, the progresss is shown. */
int preprocess_rebin(datafile_definition original, datafile_definition *clone, long NrBins, verbose_definition verbose)
{
  long p, f, n;
  int i;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Rebinning data to %ld bins\n", NrBins);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rebin: Rebinning only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rebin: Cannot handle PA data.");
    return 0;
  }
  if(NrBins > original.NrBins) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rebin: Cannot rebin to a larger amount of bins.");
    return 0;
  }
  if(original.NrBins % NrBins != 0) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_rebin: Rebinning from %ld to %ld bins implies that separate bins are not entirely independent.", original.NrBins, NrBins);
  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->fixedtsamp *= original.NrBins/(double)NrBins;
  clone->NrBins = NrBins;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  if(clone->data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rebin: Memory allocation error.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < clone->NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(n = 0; n < clone->NrSubints; n++) {
	if(!rebinPulse(&(original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))]), original.NrBins, &(clone->data[clone->NrBins*(p+clone->NrPols*(f+n*clone->NrFreqChan))]), clone->NrBins, 1, verbose)) {
	  return 0;
	}
        if(verbose.verbose && verbose.nocounters == 0) {
	  long doprint;
	  doprint = 1;
	  if(original.NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	    doprint = 0;
	  if(doprint) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*((p*clone->NrFreqChan+f)*clone->NrSubints+n))/(float)(clone->NrSubints*clone->NrFreqChan*clone->NrPols));
	    fflush(stdout);
	  }
	}
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done                               \n");
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE


/* Returns 1 if there are NaN's in the data, or zero otherwise. If
   generate_warning is set, a warning is generated. Verbose level
   determines nr of spaces before output. */
int preprocess_checknan(datafile_definition original, int generate_warning, verbose_definition verbose)
{
  long p, f, n, b;
  int i;
  float I;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Checking data for NaN's\n");
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_checknan: only works if data is loaded into memory.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	for(b = 0; b < original.NrBins; b++) {
	  if(readPulsePSRData(&original, n, p, f, b, 1, &I, verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_checknan: Cannot read data.");
	    exit(-1);
	  }
	  if(isnan(I)) {
	    if(generate_warning) {
	      printwarning(verbose.debug, "WARNING Found a NaN in file %s. First occurance is at pulse number %ld, freq channel %ld, pol channel %ld and bin number %ld.", original.filename, n+1, f+1, p+1, b+1);
	    }
	    return 1;
	  }
	}
      }
    }
  }
  return 0;
}

//START REGION DEVELOP


/* Function sets all NaN's to zero. Returns 1 if there are NaN's in
   the data, or zero otherwise. Verbose level determines nr of spaces
   before output. */
int preprocess_removenan(datafile_definition original, verbose_definition verbose)
{
  long p, f, n, b;
  int i, ret;
  float I;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Removing any NaN's from data if they occur\n");
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_removenan: only works if data is loaded into memory.");
    return 0;
  }


  ret = 0;
  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	for(b = 0; b < original.NrBins; b++) {
	  if(readPulsePSRData(&original, n, p, f, b, 1, &I, verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_removenan: Cannot read data.");
	    exit(-1);
	  }
	  if(isnan(I)) {
	    I = 0;
	    if(writePulsePSRData(&original, n, p, f, b, 1, &I, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR preprocess_removenan: Cannot write data.");
	      exit(-1);
	    }
	    ret = 1;
	  }
	}
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done             \n");
  }
  return ret;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Returns 1 if there are INF's in the data, or zero otherwise. If
   generate_warning is set, a warning is generated. Verbose level
   determines nr of spaces before output. */
int preprocess_checkinf(datafile_definition original, int generate_warning, verbose_definition verbose)
{
  long p, f, n, b;
  int i;
  float I;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Checking data for INF's\n");
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_checkinf: only works if data is loaded into memory.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	for(b = 0; b < original.NrBins; b++) {
	  if(readPulsePSRData(&original, n, p, f, b, 1, &I, verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_checkinf: Cannot read data.");
	    exit(-1);
	  }
	  if(isinf(I) != 0) {
	    if(generate_warning) {
	      printwarning(verbose.debug, "WARNING Found a INF (or -INF) in file %s. First occurance is at pulse number %ld, freq channel %ld, pol channel %ld and bin number %ld.", original.filename, n+1, f+1, p+1, b+1);
	    }
	    return 1;
	  }
	}
      }
    }
  }
  return 0;
}

//START REGION DEVELOP


/* Function sets all INF's (or -INF's) to zero. Returns 1 if there are INF's in
   the data, or zero otherwise. Verbose level determines nr of spaces
   before output. */
int preprocess_removeinf(datafile_definition original, verbose_definition verbose)
{
  long p, f, n, b;
  int i, ret;
  float I;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Removing any INF's from data if they occur\n");
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_removeinf: only works if data is loaded into memory.");
    return 0;
  }


  ret = 0;
  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	for(b = 0; b < original.NrBins; b++) {
	  if(readPulsePSRData(&original, n, p, f, b, 1, &I, verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_removeinf: Cannot read data.");
	    exit(-1);
	  }
	  if(isinf(I) != 0) {
	    I = 0;
	    if(writePulsePSRData(&original, n, p, f, b, 1, &I, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR preprocess_removeinf: Cannot write data.");
	      exit(-1);
	    }
	    ret = 1;
	  }
	}
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done             \n");
  }
  return ret;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Rotates each pulse independent of each other using an fft
   algorithm. shiftPhase sets a constant phase offset to be applied to
   each subint. If addslope is set, an additional offset is applied
   which is zero for the first subint, but for following subints it is
   increasing by slope.

   Return 1 = success. Verbose level determines nr of spaces before
   output. If nocounters is set, no progress is shown. */
int preprocess_fftshift(datafile_definition original, float shiftPhase, int addslope, float slope, verbose_definition verbose)
{
  long p, f, n;
  int i;
  float offset;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(addslope == 0)
      printf("Rotating data by %f phase\n", shiftPhase);
    else
      printf("Rotating data by %f phase + %e/subint\n", shiftPhase, slope);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_fftshift: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_fftshift: Cannot handle PA data.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	if(addslope) {
	  offset = (shiftPhase+n*slope)*original.NrBins;
	}else {
	  offset = shiftPhase*original.NrBins;
	}
	if(rotateSinglepulse(&(original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))]), original.NrBins, offset, verbose) == 0)
	  return 0;
        if(verbose.verbose && verbose.nocounters == 0) {
	  long doprint;
	  doprint = 1;
	  if(original.NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	    doprint = 0;
	  if(doprint) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*((p*original.NrFreqChan+f)*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan*original.NrPols));
	    fflush(stdout);
	  }
	}
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(verbose.verbose) printf("  done                            \n");
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE


/* Create a clone with only the selected frequency channel. The clone
   should not exist yet. Verbose level determines nr of spaces before
   output. Return 1 = success */
int preprocess_channelselect(datafile_definition original, datafile_definition *clone, long chanelnr, verbose_definition verbose)
{
  long p, f, n;
  int i;
  float *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Selecting frequency channel %ld\n", chanelnr);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_chanelselect: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_chanelselect: Cannot handle PA data.");
    return 0;
  }
  if(chanelnr < 0 || chanelnr >= original.NrFreqChan) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_chanelselect: Invalid frequency chanel number.");
    return 0;
  }

  if(original.NrSubints > 1 && original.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_chanelselect: Selecting single frequency channels from a multi-subint dataset is only implemented when the frequency channels are uniformely separated.");
    return 0;
  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->NrFreqChan = 1;
  //  clone->freq_cent = original.freq_cent-0.5*original.bw+original.channelbw*chanelnr;
  clone->freqMode = FREQMODE_UNIFORM;
  if(clone->freqlabel_list != NULL) {
    free(clone->freqlabel_list);
    clone->freqlabel_list = NULL;
  }
  //  clone->freq_list = malloc(2*sizeof(double));
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_chanelselect: Memory allocation error.");
    return 0;
  }
  set_centre_frequency(clone, get_weighted_channel_freq(original, 0, chanelnr, verbose), verbose);
  double bw;
  if(get_channelbandwidth(original, &bw, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_chanelselect: Gatting bandwidth failed.");
    return 0;
  }
  if(set_bandwidth(clone, bw, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_chanelselect: Bandwidth changing failed.");
    return 0;
  }


  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < clone->NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < clone->NrSubints; n++) {
	if(f == chanelnr) {
	  if(readPulsePSRData(&original, n, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_chanelselect: Error reading pulse.");
	    return 0;
	  }
	  if(writePulsePSRData(clone, n, p, 0, 0, clone->NrBins, pulse, verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_chanelselect: Error writing pulse.");
	    return 0;
	  }
	}
      }
    }
  }
  free(pulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done\n");
  }
  return 1;
}


//START REGION DEVELOP
//START REGION RELEASE

/* Create a clone with only the selected polarization channel. The
   clone should not exist yet. Verbose level determines nr of spaces
   before output. 

   Return 1 = success */
int preprocess_polselect(datafile_definition original, datafile_definition *clone, long polnr, verbose_definition verbose)
{
  long p, f, n;
  int i;
  float *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Selecting polarization channel %ld\n", polnr);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_polselect: only works if data is loaded into memory.");
    return 0;
  }
  //  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "ERROR preprocess_polselect: Cannot handle PA data.");
  //    return 0;
  //  }
  if(polnr < 0 || polnr >= original.NrPols) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_polselect: Invalid polarization chanel number.");
    return 0;
  }

  /*
  if(original.NrPols == 1) {
    cleanPSRData(clone, verbose);
    copy_params_PSRData(original, clone, verbose);
    clone->data = original.data;
     In swap funcion the free of the original will distroy clone/
  }else {
  */
  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  //  if(clone->gentype == GENTYPE_ILVPA) {
  //    clone->gentype = GENTYPE_SUBINTEGRATIONS;
  //    fflush(stdout);
  //    printwarning(verbose.debug, "WARNING preprocess_polselect: Changing gentype from ILVPA to SUBINTEGRATIONS.");
  //  }

  clone->NrPols = 1;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_polselect: Memory allocation error.");
    return 0;
  }
    
    
    /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(n = 0; n < clone->NrSubints; n++) {
	if(p == polnr) {
	  if(readPulsePSRData(&original, n, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_polselect: Error reading pulse.");
	    return 0;
	  }
	  if(writePulsePSRData(clone, n, 0, f, 0, clone->NrBins, pulse, verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_polselect: Error writing pulse.");
	    return 0;
	  }
	}
      }
    }
  }
  
  free(pulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done              \n");
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Create a clone with only the selected pulses. The clone should not
   exist yet. Verbose level determines nr of spaces before
   output. Return 1 = success */
int preprocess_pulsesselect(datafile_definition original, datafile_definition *clone, long nskip, long nread, verbose_definition verbose)
{
  long p, f, n, i;
  float *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Selecting pulses %ld-%ld\n", nskip, nread+nskip-1);
  }
  /*
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_pulsesselect: only works if data is loaded into memory.");
    return 0;
  }
  */
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_pulsesselect: Cannot handle PA data.");
    return 0;
  }
  if(nskip < 0 || nskip >= original.NrSubints) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_pulsesselect: Invalid nskip.");
    return 0;
  }
  if(nread < 0 || nskip+nread > original.NrSubints) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_pulsesselect: Invalid nread.");
    return 0;
  }
  if(original.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_pulsesselect: Frequency channels are not necessarily uniformly separated.");
    return 0;
  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->NrSubints = nread;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  clone->format = MEMORY_format;
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_pulsesselect: Memory allocation error.");
    return 0;
  }


  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(n = nskip; n < nskip+nread; n++) {
	if(readPulsePSRData(&original, n, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_pulsesselect: Error reading pulse.");
	  return 0;
	}
	if(writePulsePSRData(clone, n-nskip, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_pulsesselect: Error writing pulse.");
	  return 0;
	}
	if(original.format != MEMORY_format) {
	  if(verbose.verbose && verbose.nocounters == 0) {
	    long doprint;
	    doprint = 1;
	    if(clone->NrFreqChan > 4 && n != nskip)  // Avoid excessive nr of prints
	      doprint = 0;
	    if(doprint) {
	      for(i = 0; i < verbose.indent; i++)      
		printf(" ");
	      printf("  %.1f%%     \r", (100.0*(((n-nskip)+nread*f+p*nread*clone->NrFreqChan))/(float)(nread*clone->NrFreqChan*clone->NrPols)));
	      fflush(stdout);
	    }
	  }
	}
      }
    }
  }
  free(pulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done            \n");
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Create a clone with only a multiple of blocksize
   subintegrations. The clone should not exist yet. Verbose level
   determines nr of spaces before output. Return 1 = success */
int preprocess_blocksize(datafile_definition original, datafile_definition *clone, int blocksize, verbose_definition verbose)
{
  int i;
  long nread;
  verbose_definition verbose2;

  nread = original.NrSubints/blocksize;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Selecting %ld blocks of %d subints\n", nread, blocksize);
  }

  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;

  return preprocess_pulsesselect(original, clone, 0, nread*blocksize, verbose2);
}

//START REGION DEVELOP

/* Create a clone with an inverted bin-axis and subint-axis. The clone
   should not exist yet. Verbose level determines nr of spaces before
   output. Return 1 = success */
int preprocess_invertXY(datafile_definition original, datafile_definition *clone, verbose_definition verbose)
{
  long p, f, n, b;
  int i;
  float *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Inverting bins/subints\n");
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertXY: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    if((original.gentype != GENTYPE_PADIST && original.gentype != GENTYPE_ELLDIST) || original.NrPols != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_invertXY: Cannot handle PA data.");
      return 0;
    }
  }
  if(original.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertXY: Frequency channels are not necessarily uniformly separated.");
    return 0;
  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->NrBins = original.NrSubints;
  clone->NrSubints = original.NrBins;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((original.NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertXY: Memory allocation error.");
    return 0;
  }


  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	/*	fprintf(stderr, "preprocess_invertXY: XXXXXXXXXXX\n"); */
	if(readPulsePSRData(&original, n, p, f, 0, original.NrBins, pulse, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_invertXY: Error reading pulse.");
	  return 0;
	}
	/*	   fprintf(stderr, "preprocess_invertXY: YYYYYYYYYYYY\n"); */
	for(b = 0; b < original.NrBins; b++) {
	  if(writePulsePSRData(clone, b, p, f, n, 1, &pulse[b], verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_invertXY: Error writing pulse.");
	    return 0;
	  }
	}
      }
    }
  }
  free(pulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done        \n");
  }
  return 1;
}

/* Create a clone with an inverted frequency-axis (the order of the
   channels is inverted). Verbose level determines nr of spaces before
   output. The clone should not exist yet. Return 1 = success */
int preprocess_invertF(datafile_definition original, datafile_definition *clone, verbose_definition verbose)
{
  long p, f, n;
  int i;
  float *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Inverting frequency channel order\n");
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertF: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertF: Cannot handle PA data.");
    return 0;
  }
  if(original.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertF: Cannot handle non-uniform frequency sampling.");
    return 0;
  }
  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  double bw = get_bandwidth(original, verbose);
  if(set_bandwidth(clone, -1.0*bw, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertF: Bandwidth changing failed.");
    return 0;
  }
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((original.NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertF: Memory allocation error.");
    return 0;
  }


  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	/*	fprintf(stderr, "preprocess_invertF: XXXXXXXXXXX\n"); */
	if(readPulsePSRData(&original, n, p, f, 0, original.NrBins, pulse, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_invertF: Error reading pulse.");
	  return 0;
	}
	/*	   fprintf(stderr, "preprocess_invertF: YYYYYYYYYYYY\n"); */
	if(writePulsePSRData(clone, n, p, original.NrFreqChan-f-1, 0, original.NrBins, pulse, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_invertF: Error writing pulse.");
	  return 0;
	}
      }
    }
  }
  free(pulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done               \n");
  }
  return 1;
}

/* Create a clone with an inverted frequency-axis and subint-axis. In
   fact this is just a header opperation. Verbose level determines nr
   of spaces before output. Return 1 = success */
int preprocess_invertFP(datafile_definition *original, verbose_definition verbose)
{
  int i;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Inverting frequency channels/subints\n");
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertfp: only works if data is loaded into memory.");
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertfp: Cannot handle PA data.");
    return 0;
  }
  if(original->NrSubints != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertfp: only works when there is one subint present.");
    return 0;
  }
  if(original->freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertfp: Frequency channels are not necessarily uniformly separated.");
    return 0;
  }

  original->NrSubints = original->NrFreqChan;
  original->NrFreqChan = 1;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done            \n");
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Create a clone with an inverted bin-axis and frequency channel
   axis. The clone should not exist yet. Verbose level determines nr
   of spaces before output. Return 1 = success */
int preprocess_invertFX(datafile_definition original, datafile_definition *clone, verbose_definition verbose)
{
  long p, f, n, b;
  int i;
  float *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Inverting frequency channels/bins\n");
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertFX: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    if((original.gentype != GENTYPE_PADIST && original.gentype != GENTYPE_ELLDIST) || original.NrPols != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_invertFX: Cannot handle PA data.");
      return 0;
    }
  }
  if(original.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertFX: Frequency channels are not necessarily uniformly separated.");
    return 0;
  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->NrBins = original.NrFreqChan;
  clone->NrFreqChan = original.NrBins;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((original.NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_invertFX: Memory allocation error.");
    return 0;
  }


  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	/*	fprintf(stderr, "preprocess_invertXY: XXXXXXXXXXX\n"); */
	if(readPulsePSRData(&original, n, p, f, 0, original.NrBins, pulse, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_invertFX: Error reading pulse.");
	  return 0;
	}
	/*	   fprintf(stderr, "preprocess_invertXY: YYYYYYYYYYYY\n"); */
	for(b = 0; b < original.NrBins; b++) {
	  if(writePulsePSRData(clone, n, p, b, f, 1, &pulse[b], verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_invertFX: Error writing pulse.");
	    return 0;
	  }
	}
      }
    }
  }
  free(pulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done        \n");
  }
  return 1;
}


//START REGION DEVELOP
//START REGION RELEASE

/* Create a clone with a different order of the axis that makes it
   possible to directly plot raw filterbank data. It will have subints
   on the horizontal axis (each nrbins long) and the frequency
   channels on the vertical axis. This block of data for the first
   polarization is followed by the next polarization.

   Return 1 = success. Verbose
   level of verbose determines nr of spaces before output. */
int preprocess_transposeRawFBdata(datafile_definition original, datafile_definition *clone, verbose_definition verbose)
{
  long p, f, n;
  float *pulse;
  int i;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Transposing data\n");
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_transposeRawFBdata: only works if data is loaded into memory.");
    return 0;
  }
  //  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "ERROR preprocess_transposeRawFBdata: Cannot handle PA data.");
  //    return 0;
  //  }
  if(original.isTransposed) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_transposeRawFBdata: Data already appears to be transposed.");
    return 0;
  }
  //  if(original.freqMode != FREQMODE_UNIFORM) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "ERROR preprocess_transposeRawFBdata: Frequency channels are not necessarily uniformly separated.");
  //    return 0;
  //  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_transposeRawFBdata: Memory allocation error.");
    return 0;
  }


  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(n = 0; n < clone->NrSubints; n++) {
	if(readPulsePSRData(&original, n, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_transposeRawFBdata: Error reading pulse.");
	  return 0;
	}
	memcpy(&(clone->data[original.NrBins*(n+original.NrSubints*(f+p*original.NrFreqChan))]), pulse, sizeof(float)*(clone->NrBins));
      }
    }
  }

  clone->isTransposed = 1;
  free(pulse);
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Create a clone with nrpulses succesive pulses added together
   (i.e. if set to the total nr of pulses you get a pulse profile). If
   nrpulses is negative each subint is written out -nrpulses times. If
   complete is set, each subint are thrown away to make each subint
   the sum of the same number of input subints. The clone should not
   exist yet. Return 1 = success. Verbose level determines nr of
   spaces before output. If nocounters is set the progress is not
   shown. */
int preprocess_addsuccessivepulses(datafile_definition original, datafile_definition *clone, long nrpulses, int complete, verbose_definition verbose)
{
  long p, f, n, n2, b;
  float *pulse, *addedpulse;
  int i, use_depar;
  datafile_definition clone_depar;
  verbose_definition verbose2;

  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;

  use_depar = 0;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(nrpulses >= 0)
      printf("Add each %ld succesive subints\n", nrpulses);
    else
      printf("Write out each subint %ld times\n", -nrpulses);
  }
  /*
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: only works if data is loaded into memory.");
    return 0;
  }
  */
  // ILVPA ends up with angles outside range etc, normalisation not easy, as could have many zero's, so nr of detections different for each bin. Adding errorbars even more difficult.
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Cannot handle position angle data when adding subints. Add subints using Stokes parameter data.");
    return 0;
  }
  if(nrpulses > original.NrSubints || nrpulses == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Invalid number of subints to add.");
    return 0;
  }
  if(original.freqMode != FREQMODE_UNIFORM) {
    if(nrpulses != 1) {  // If no time-scrunching is done, this isn't an issue
      /*
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Frequency channels are not necessarily uniformly separated.");
      return 0;
      */
      if(preprocess_dedisperse(&original, 0, 1, -1.0, verbose) != 2) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Dedispersion failed.");
	return 0;
      }
    }
  }

  // If four polarizations, do parallactic angle correction first if required.
  if(original.NrPols == 4 && original.isDePar == -1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_addsuccessivepulses: Parallactic angle correction state is unknown, no correction will be done.");
  }else if(original.NrPols == 4 && original.isDePar == 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Correcting for parallactic angle changes first.\n");
    }
    if(preprocess_corrParAng(&original, &clone_depar, 0, verbose2) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Parallactic angle correction failed.");
      return 0;
    }
    use_depar = 1;
  }

  cleanPSRData(clone, verbose);
  if(use_depar)
    copy_params_PSRData(clone_depar, clone, verbose);
  else
    copy_params_PSRData(original, clone, verbose);
  if(nrpulses > 0) {
    clone->NrSubints = original.NrSubints/nrpulses;
    if(complete == 0) {
      if(clone->NrSubints * nrpulses != original.NrSubints) {
	printwarning(verbose.debug, "WARNING preprocess_addsuccessivepulses: Last subint has a different duration.");
	//Things like observation duration will not be correct in the output.
	clone->NrSubints += 1;  // Since it was rounded down, make sure no data is lost
      }
    }
  }else {
    clone->NrSubints = original.NrSubints*(-nrpulses);
  }
  clone->tsubMode = TSUBMODE_TSUBLIST;
  if(clone->tsub_list != NULL)
    free(clone->tsub_list);
  clone->tsub_list = (double *)malloc(clone->NrSubints * sizeof(double));
  if(clone->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Memory allocation error");
    return 0;
  }
  //  if(nrpulses > 0) {
  //    clone->fixedtsub = nrpulses*get_tsub(original, 0, verbose);
  //  }else {
  //    clone->fixedtsub = -get_tsub(original, 0, verbose)/(double)nrpulses;
  //  }
  clone->format = MEMORY_format;
     
  if(clone->gentype == GENTYPE_PULSESTACK) {
    clone->gentype = GENTYPE_SUBINTEGRATIONS;
    //    if(fabs(get_tobs(*clone, verbose)) < 0.0001) { THIS SHOULD ALREADY BE DONE WHEN READING IN HEADER
    //      clone->tsubMode = TSUBMODE_FIXEDTSUB;
    //      if(nrpulses > 0)
    //	clone->fixedtsub = nrpulses*get_period(original, 0, verbose);
    //      else
    //	clone->fixedtsub = -get_period(original, 0, verbose)/(double)nrpulses;
    //      //      clone->tobs = original.NrSubints * get_period(original, 0, verbose);
    //    }
  }
  if(clone->gentype == GENTYPE_SUBINTEGRATIONS && clone->NrSubints == 1) {
    clone->gentype = GENTYPE_PROFILE;
  }
  if(clone->gentype != GENTYPE_PROFILE && clone->gentype != GENTYPE_SUBINTEGRATIONS && clone->gentype != GENTYPE_PULSESTACK && clone->gentype != GENTYPE_UNDEFINED && clone->gentype != GENTYPE_POLNCAL) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_addsuccessivepulses: Unsure about adding subints for a %s file. Setting gentype to undefined.", returnGenType_str(clone->gentype));
    clone->gentype = GENTYPE_UNDEFINED;
  }

  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  addedpulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL || addedpulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Memory allocation error.");
    return 0;
  }


  /* Loop over pulses and frequency channels etc. */
  double curtsub;
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      curtsub = 0;
      for(n = 0; n < clone->NrSubints; n++) {
	if(nrpulses > 0) {
	  for(b = 0; b < original.NrBins; b++)
	    addedpulse[b] = 0;
	  for(n2 = 0; n2 < nrpulses; n2++) {
	    if(n*nrpulses+n2 < original.NrSubints) {
	      if(use_depar) {
		if(readPulsePSRData(&clone_depar, n*nrpulses+n2, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Error reading pulse.");
		  return 0;
		}
		curtsub += get_tsub(clone_depar, n*nrpulses+n2, verbose);
		//		if(f == 0 && p == 0) {
		//		  fprintf(stderr, "+=%lf (%ld)\n", get_tsub(clone_depar, n*nrpulses+n2, verbose), n*nrpulses+n2);
		//		}
	      }else {
		if(readPulsePSRData(&original, n*nrpulses+n2, p, f, 0, clone->NrBins, pulse, verbose) != 1) {
		  fflush(stdout);
		  printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Error reading pulse.");
		  return 0;
		}
		curtsub += get_tsub(original, n*nrpulses+n2, verbose);
	      }
	      for(b = 0; b < original.NrBins; b++)
		addedpulse[b] += pulse[b];
	    }
	  }
	}else {
	  if(use_depar) {
	    if(readPulsePSRData(&clone_depar, n/(-nrpulses), p, f, 0, clone->NrBins, addedpulse, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Error reading pulse.");
	      return 0;
	    }
	    curtsub += get_tsub(clone_depar, n/(-nrpulses), verbose)/(double)(-nrpulses);
	  }else {
	    if(readPulsePSRData(&original, n/(-nrpulses), p, f, 0, clone->NrBins, addedpulse, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Error reading pulse.");
	      return 0;
	    }
	    curtsub += get_tsub(original, n/(-nrpulses), verbose)/(double)(-nrpulses);
	  }
	}
	clone->tsub_list[n] = curtsub;
	//	if(f == 0 && p == 0)
	  //	  fprintf(stderr, "XXXXX n=%ld %lf\n", n, clone->tsub_list[n]);
	if(writePulsePSRData(clone, n, p, f, 0, clone->NrBins, addedpulse, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_addsuccessivepulses: Error writing pulse.");
	  return 0;
	}
	curtsub = 0;
        if(verbose.verbose && verbose.nocounters == 0) {
	  long doprint;
	  doprint = 1;
	  if(clone->NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	    doprint = 0;
	  if(doprint) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*((p*clone->NrFreqChan+f)*clone->NrSubints+n))/(float)(clone->NrSubints*clone->NrFreqChan*clone->NrPols));
	    fflush(stdout);
	  }
	}
      }
    }
  }
  free(pulse);
  free(addedpulse);
  if(use_depar) {
    closePSRData(&clone_depar, 0, verbose);    
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(verbose.verbose) printf("  done                            \n");
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Create a clone with nrfreq succesive frequency channels added
   together (i.e. if set to the total nr of frequency channels you get
   a dedispersed time series). The data should already be dedispersed
   and de-Faraday rotated. If nrfreq is negative each frequency
   channel is written out -nrfreq times. The clone should not exist
   yet. Return 1 = success, 0 on error. Verbose level determines nr of
   spaces before output. This function is almost identical to
   preprocess_addsuccessivepulses. If fzapMask is != NULL, the
   frequency channels that are != 0 in the array are ignored. If
   nocounters is set the progress is not shown. */
int preprocess_addsuccessiveFreqChans(datafile_definition original, datafile_definition *clone, long nrfreq, int *fzapMask, verbose_definition verbose)
{
  long p, f, n, n2, b;
  float *pulse, *addedpulse;
  int i, ok;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(nrfreq >= 0)
      printf("Add %ld succesive frequency channels\n", nrfreq);
    else
      printf("Write out each frequency channel %ld times\n", -nrfreq);
  }
  /*
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: only works if data is loaded into memory.\n");
    return 0;
    }*/
  // ILVPA ends up with angles outside range etc, normalisation not easy, as could have many zero's, so nr of detections different for each bin. Adding errorbars even more difficult.
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Cannot handle position angle data when adding frequency channels. Add subints using Stokes parameter data.");
    return 0;
  }
  if(nrfreq > original.NrFreqChan || nrfreq == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Invalid number of frequency channels to add.");
    return 0;
  }
  //  if(original.freqMode != FREQMODE_UNIFORM) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Frequency channels are not necessarily uniformly separated.");
  //    return 0;
  //  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  if(nrfreq > 0) {
    clone->NrFreqChan = original.NrFreqChan/nrfreq;
    if(clone->NrFreqChan * nrfreq != original.NrFreqChan) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_addsuccessiveFreqChans: Last channel will be the sum of a different number of channels compared to the others. The frequency labeling will no longer be correct.");
    }
  }else {
    clone->NrFreqChan = original.NrFreqChan*(-nrfreq);
  }
  clone->format = MEMORY_format;
  //  clone->channelbw = clone->bw/(double)clone->NrFreqChan;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  //  printf("XXXXXXX alloc %p\n", clone->data);
  pulse = (float *)malloc((clone->NrBins)*sizeof(float));
  addedpulse = (float *)malloc((clone->NrBins)*sizeof(float));
  if(clone->data == NULL || pulse == NULL || addedpulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Memory allocation error.");
    return 0;
  }
  if(original.freqMode == FREQMODE_FREQTABLE) {
    if(clone->freqlabel_list != NULL)
      free(clone->freqlabel_list);
    clone->freqlabel_list = malloc((clone->NrFreqChan)*(clone->NrSubints)*sizeof(double));
    if(clone->freqlabel_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Memory allocation error.");
      return 0;
    }
  }
  double newfreq;  // The new frequency of the channel
  long newfreq_nradded; // newfreq is initially a cummulative sum, so keep track of how many entries
  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < clone->NrPols; p++) {
    for(n = 0; n < clone->NrSubints; n++) {
      for(f = 0; f < clone->NrFreqChan; f++) {
	newfreq = 0;
	newfreq_nradded = 0;
	if(nrfreq > 0) {  // Sum channels
	  for(b = 0; b < original.NrBins; b++)
	    addedpulse[b] = 0;
	  for(n2 = 0; n2 < nrfreq; n2++) {
	    ok = 1;
	    if(fzapMask != NULL) {
	      if(fzapMask[f*nrfreq+n2] != 0) {
		ok = 0;
	      }
	    }
	    if(ok) {
	      if(readPulsePSRData(&original, n, p, f*nrfreq+n2, 0, clone->NrBins, pulse, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Error reading pulse.");
		return 0;
	      }
	      if(original.freqMode == FREQMODE_FREQTABLE) {
		newfreq += get_weighted_channel_freq(original, n, f*nrfreq+n2, verbose);
		newfreq_nradded++;
	      }
	      for(b = 0; b < original.NrBins; b++)
		addedpulse[b] += pulse[b];
	    }
	  }
	  newfreq /= (double)newfreq_nradded;
	}else { // Duplicate channels
	  ok = 1;
	  if(fzapMask != NULL) {
	    if(fzapMask[f/(-nrfreq)] != 0) {
	      ok = 0;
	    }
	  }
	  if(ok) {
	    if(readPulsePSRData(&original, n, p, f/(-nrfreq), 0, clone->NrBins, addedpulse, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Error reading pulse.");
	      return 0;
	    }
	  }else {
	    for(b = 0; b < original.NrBins; b++)
	      addedpulse[b] = 0;
	  }
	  if(original.freqMode == FREQMODE_FREQTABLE) {
	    newfreq = get_weighted_channel_freq(original, n, f/(-nrfreq), verbose);
	  }
	}
	if(original.freqMode == FREQMODE_FREQTABLE) {
	  if(set_weighted_channel_freq(clone, n, f, newfreq, verbose) == 0) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Setting frequency labeling failed.");
	    return 0;
	  }
	}
	if(writePulsePSRData(clone, n, p, f, 0, clone->NrBins, addedpulse, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_addsuccessiveFreqChans: Error writing pulse.");
	  return 0;
	}
        if(verbose.verbose && verbose.nocounters == 0) {
	  long doprint;
	  doprint = 1;
	  if(clone->NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	    doprint = 0;
	  if(doprint) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*((p*clone->NrFreqChan+f)*clone->NrSubints+n))/(float)(clone->NrSubints*clone->NrFreqChan*clone->NrPols));
	    fflush(stdout);
	  }
	}
      }
    }
  }
  free(pulse);
  free(addedpulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done                            \n");
  }
  return 1;
}


//START REGION DEVELOP
//START REGION RELEASE

/* Subtracts baseline from data. If onpulse != NULL, the onpulse
   region will be excluded from the calculation of the averages. If
   baseline != NULL, the subtracted baseline values are stored to be
   used in preprocess_restore_debase(). If remove_slope = 0, a simple
   offset is subtracted. If 1, a linear gradient is removed from the
   offpulse region. Memory will be allocated to store this
   information.

   Returns 0 on error. */
int preprocess_debase(datafile_definition *original, pulselongitude_regions_definition *onpulse, float **baseline, int remove_shape, verbose_definition verbose)
{
  long p, f, n, j, nrOffpulseBins, offpulse_bin_nr;
  int i;
  float avrg, *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Subtracting baseline\n");
  }

  if(onpulse != NULL) {
    region_frac_to_int(onpulse, original->NrBins, 0);
  }

  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_debase: only works if data is loaded into memory.");
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_debase: Cannot handle PA data.");
    return 0;
  }


  pulse = (float *)malloc((original->NrBins)*sizeof(float));
  if(pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_debase: Memory allocation error.");
    return 0;
  }
  if(baseline != NULL) {
    *baseline = (float *)malloc((original->NrPols)*(original->NrFreqChan)*(original->NrSubints)*sizeof(float));
    if(*baseline == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_debase: Memory allocation error.");
      return 0;
    }
  }
  double *profilex_double, *profile_double; //, , *profileSigma_double;
  if(remove_shape != 0) {
    profilex_double = (double *)malloc((original->NrBins)*sizeof(double));
    profile_double = (double *)malloc((original->NrBins)*sizeof(double));
    if(profilex_double == NULL || profile_double == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_debase: Memory allocation error.");
      return 0;
    }
    offpulse_bin_nr = 0;
    for(j = 0; j < original->NrBins; j++) {   // Store which bin numbers have offpulse data, to be used in linear fit
      if(onpulse != NULL) {
	if(checkRegions(j, onpulse, 0, verbose) == 0) {
	  profilex_double[offpulse_bin_nr++] = j;
	}
      }else {
	profilex_double[offpulse_bin_nr++] = j;
      }
    }
  }
  /* Loop over pulses and frequency channels etc. */
  long baseline_index = 0;
  int ok;
  for(p = 0; p < original->NrPols; p++) {
    for(f = 0; f < original->NrFreqChan; f++) {
      for(n = 0; n < original->NrSubints; n++) {
	if(readPulsePSRData(original, n, p, f, 0, original->NrBins, pulse, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_debase: Error reading data.");
	  return 0;
	}
	nrOffpulseBins = 0;
	avrg = 0;
	offpulse_bin_nr = 0;
	for(j = 0; j < original->NrBins; j++) {
	  ok = 0;
	  if(onpulse != NULL) {
	    if(checkRegions(j, onpulse, 0, verbose) == 0 || onpulse->nrRegions == 0) {
	      ok = 1;
	    }
	  }else { // Same as above, but for all bins
	    ok = 1;
	  }
	  if(ok) {
	    avrg += pulse[j];
	    nrOffpulseBins++;
	    if(remove_shape != 0) {
	      profile_double[offpulse_bin_nr++] = pulse[j]; 
	    }
	  }
	}
	if(nrOffpulseBins > 0) {
	  if(remove_shape == 0) {
	    avrg /= (float)nrOffpulseBins;
	    /*float min, max;
	      min=max=pulse[0];*/
	    for(j = 0; j < original->NrBins; j++) {
	      pulse[j] -= avrg; 
	      /*	    if(pulse[j] > max)
			    max=pulse[j];
			    if(pulse[j] < min)
			    min=pulse[j];*/
	    }
	  }else if(remove_shape == 1) {
	    // Note: this is very similar to the code in pmod()
	    double a, b, cov00, cov01, cov11, sumsq;
	    gsl_fit_linear(profilex_double, 1, profile_double, 1, offpulse_bin_nr, &a, &b, &cov00, &cov01, &cov11, &sumsq);
	    //	    printf("XXXXX fit: %lf + %lf*x\n", a, b);
	    for(j = 0; j < original->NrBins; j++) {
	      /*
		if(l == 0 && k == 0) {
		printf("%f %f %f\n", (float)j, profileI[j], a+profilex_double[j]*b);
		}*/
	      pulse[j] -= a+(float)j*b;  
	    }
	    avrg = 0;
	  }
	  /*	  printf("XXXXXX %f %f\n", min, max); */
	  if(writePulsePSRData(original, n, p, f, 0, original->NrBins, pulse, verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_debase: Error writing data.");
	    return 0;
	  }
	  if(baseline != NULL) {
	    (*baseline)[baseline_index++] = avrg;
	  }
	}else {
	  if(baseline != NULL) {
	    (*baseline)[baseline_index++] = 0.0;
	  }
	}
        if(verbose.verbose && verbose.nocounters == 0) {
	  long doprint;
	  doprint = 1;
	  if(original->NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	    doprint = 0;
	  if(doprint) {
	    for(i = 0; i < verbose.indent; i++)
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*((p*original->NrFreqChan+f)*original->NrSubints+n))/(float)(original->NrSubints*original->NrFreqChan*original->NrPols));
	    fflush(stdout);
	  }
	}

      }
    }
  }
  original->isDebase = 1;
  free(pulse);
  if(remove_shape != 0) {
    free(profilex_double);
    free(profile_double);
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(verbose.verbose) printf("  done                            \n");
  }
  return 1;
}

//START REGION RELEASE
//START REGION DEVELOP

/* Restores the baseline as was previously subtracted with preprocess_debase(). 

   Returns 0 on error. */
int preprocess_restore_debase(datafile_definition *original, float *baseline, verbose_definition verbose)
{
  long p, f, n, j;
  int i;
  float *pulse;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Restoring baseline\n");
  }

  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_restore_debase: only works if data is loaded into memory.");
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_restore_debase: Cannot handle PA data.");
    return 0;
  }

  pulse = (float *)malloc((original->NrBins)*sizeof(float));
  if(pulse == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_restore_debase: Memory allocation error.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  long baseline_index = 0;
  for(p = 0; p < original->NrPols; p++) {
    for(f = 0; f < original->NrFreqChan; f++) {
      for(n = 0; n < original->NrSubints; n++) {
	if(baseline[baseline_index] != 0.0) {
	  if(readPulsePSRData(original, n, p, f, 0, original->NrBins, pulse, verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_restore_debase: Error reading data.");
	    return 0;
	  }
	  for(j = 0; j < original->NrBins; j++) {
	    pulse[j] += baseline[baseline_index]; 
	  }
	  if(writePulsePSRData(original, n, p, f, 0, original->NrBins, pulse, verbose) != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR preprocess_restore_debase: Error writing data.");
	    return 0;
	  }
	}
	baseline_index++;
        if(verbose.verbose && verbose.nocounters == 0) {
	  long doprint;
	  doprint = 1;
	  if(original->NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	    doprint = 0;
	  if(doprint) {
	    for(i = 0; i < verbose.indent; i++)
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*((p*original->NrFreqChan+f)*original->NrSubints+n))/(float)(original->NrSubints*original->NrFreqChan*original->NrPols));
	    fflush(stdout);
	  }
	}
      }
    }
  }
  original->isDebase = 0;
  free(pulse);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(verbose.verbose) printf("  done                            \n");
  }
  return 1;
}


//START REGION DEVELOP
//START REGION RELEASE

static int preprocess_addNoise_random_nr_generater_initialized = 0;
static gsl_rng *preprocess_addNoiserand_num_gen;
static const gsl_rng_type *rand_num_gen_type;

/* Makes a clone with Gaussian noise with the given standard deviation
   added to the data. The clone should not exist yet. Verbose level
   determines nr of spaces before output. If nocounters is not set,
   the progresss is shown. Return 1 = success. */
int preprocess_addNoise(datafile_definition original, datafile_definition *clone, float rms, verbose_definition verbose)
{
  long p, f, n, b;
  int i;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Adding noise to data with rms %f\n", rms);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addNoise: Adding noise only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    if(original.NrPols == 1) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_addNoise: The polarization state suggests this is not necessarily a Stokes parameter, hence the error distribution is not necessarily Gaussian, which is assumed to be the case. This might be a problem.");
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_addNoise: The polarization channels appear to include L and PA, hence the error distribution is not necessarily Gaussian.");
      return 0;
    }
  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  if(clone->data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_addNoise: Memory allocation error.");
    return 0;
  }

  if(preprocess_addNoise_random_nr_generater_initialized == 0) {
    gsl_rng_env_setup();    /* Set the default generators, can be influenced by environment variables */
    rand_num_gen_type = gsl_rng_default;
    preprocess_addNoiserand_num_gen = gsl_rng_alloc(rand_num_gen_type);
    preprocess_addNoise_random_nr_generater_initialized = 1;
    // If randomizing seed, make also a non-randomized option, such that in testing the results are reproducable
    // see preprocess_shuffle how to implement randomize
  }

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < clone->NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(n = 0; n < clone->NrSubints; n++) {
	if(writePulsePSRData(clone, n, p, f, 0, clone->NrBins, &(original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))]), verbose) != 1) {
	  return 0;
	}
	for(b = 0; b < clone->NrBins; b++) {
	  /* NR code:	  clone->data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] += rms*gasdev(&preprocess_addNoise_idnum); */
	  clone->data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] += gsl_ran_gaussian(preprocess_addNoiserand_num_gen, rms);
	}
	if(verbose.verbose && verbose.nocounters == 0) {
	  long doprint;
	  doprint = 1;
	  if(clone->NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	    doprint = 0;
	  if(doprint) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*((p*clone->NrFreqChan+f)*clone->NrSubints+n))/(float)(clone->NrSubints*clone->NrFreqChan*clone->NrPols));
	    fflush(stdout);
	  }
	}
      }
    }
  }
  /* Don't free the random number generator, as we want each call to
     be unique, as it is used for instance to bootstrap data etc.

     gsl_rng_free (rand_num_gen); 
  */

  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(verbose.verbose) printf("done                              \n");
  }
  return 1;
}

//START REGION RELEASE

static int preprocess_shuffle_random_nr_generater_initialized = 0;
static gsl_rng *preprocess_shufflerand_num_gen;
static const gsl_rng_type *rand_num_gen_type;

/* Makes a clone with the subints placed in a random order. The clone
   should not exist yet. If fixseed is set, the random number
   generator seed is initialised with a fixed seed to make results
   reproducable. Verbose level determines nr of spaces before
   output. If nocounters is not set, the progresss is shown.

   Return 1 = success, 0 = error. */
int preprocess_shuffle(datafile_definition original, datafile_definition *clone, int fixseed, verbose_definition verbose)
{
  long p, f, n;
  int i;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Shuffeling order of subints\n");
  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->format = MEMORY_format;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  if(clone->data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_shuffle: Memory allocation error.");
    return 0;
  }
  if(original.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_shuffle: Frequency channels are not necessarily uniformly separated.");
    return 0;
  }

  if(preprocess_shuffle_random_nr_generater_initialized == 0) {
    long idnum;
    gsl_rng_env_setup();    /* Set the default generators, can be influenced by environment variables */
    rand_num_gen_type = gsl_rng_default;
    preprocess_shufflerand_num_gen = gsl_rng_alloc(rand_num_gen_type);
    if(fixseed)
      idnum = 1;
    else
      randomize_idnum(&idnum);
    gsl_rng_set(preprocess_shufflerand_num_gen, idnum);
    preprocess_shuffle_random_nr_generater_initialized = 1;
  }

  long *subintlist;
  subintlist = (long *)malloc(clone->NrSubints*sizeof(long));
  if(subintlist == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_shuffle: Memory allocation error.");
    return 0;
  }
  for (n = 0; n < clone->NrSubints; n++) {
    subintlist[n] = n;
  }
  gsl_ran_shuffle (preprocess_shufflerand_num_gen, subintlist, clone->NrSubints, sizeof(long));  // Shuffle the order

  /* Loop over pulses and frequency channels etc. */
  for(n = 0; n < clone->NrSubints; n++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(p = 0; p < clone->NrPols; p++) {
	if(readPulsePSRData(&original, n, p, f, 0, clone->NrBins, &(clone->data[clone->NrBins*(p+clone->NrPols*(f+subintlist[n]*clone->NrFreqChan))]), verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_shuffle: Read error.");
	  free(subintlist);
	  return 0;
	}
	if(verbose.verbose && verbose.nocounters == 0) {
	  long doprint;
	  doprint = 1;
	  if(clone->NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	    doprint = 0;
	  if(doprint) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*((p*clone->NrFreqChan+f)*clone->NrSubints+n))/(float)(clone->NrSubints*clone->NrFreqChan*clone->NrPols));
	    fflush(stdout);
	  }
	}
      }
    }
  }
  /* Don't free the random number generator, as we want each call to
     be unique, as it is used for instance to bootstrap data etc.

     gsl_rng_free (rand_num_gen); 
  */

  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(verbose.verbose) printf("done                              \n");
  }
  free(subintlist);
  return 1;
}

//START REGION DEVELOP


/* Subtract a given RVM model (defined by alpha, beta, pa0 and l0 in
   degrees) from the data by rotating the Q & U parameters. The clone
   should not exist yet. If inplace is set, the variable clone is
   ignored, and original is modified instead. Verbose level determines
   nr of spaces before output. If nocounters is not set, the progresss
   is shown.

   Return 1 = success, 0 = error. */
int preprocess_subtractRVM(datafile_definition *original, datafile_definition *clone, int inplace, float alpha, float beta, float pa0, float l0, verbose_definition verbose)
{
  int b, n;
  float *angle_array, longitude;

  if(verbose.verbose) {
    for(n = 0; n < verbose.indent; n++)      
      printf(" ");
    printf("Subtracting RVM with alpha=%f deg, beta=%f deg, pa0=%f deg, l0=%f deg\n", alpha, beta, pa0, l0);
  }

  angle_array = malloc(original->NrBins*sizeof(float));
  if(angle_array == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_subtractRVM: Memory allocation failed.");
    return 0;
  }

  for(b = 0; b < original->NrBins; b++) {
    longitude = get_pulse_longitude(*original, 0, b, verbose);
    //    printf("Bin %d gives longitude=%f deg\n", b, longitude);
    angle_array[b] =  -2.0*paswing(alpha, beta, longitude, pa0, l0, 0, NULL, NULL, 1e10, 0)*M_PI/180.0;
  }

  int ret;
  ret = preprocess_rotateStokes(original, clone, inplace, -1, 0, angle_array, 1, 2, verbose);
  if(ret == 0) {
    return 0;
  }

  free(angle_array);
  return 1;
}


/* Subtract the average (over pulse number) PA from each sample of
   each bin from the data by rotating the Q & U parameters. The clone
   should not exist yet. If inplace is set, the variable clone is
   ignored, and original is modified instead. Verbose level determines
   nr of spaces before output. If nocounters is not set, the progresss
   is shown.

   Return 1 = success, 0 = error. */
int preprocess_subtractAveragePA(datafile_definition *original, datafile_definition *clone, int inplace, verbose_definition verbose)
{
  int b, n;
  float *angle_array;
  datafile_definition profile;

  if(verbose.verbose) {
    for(n = 0; n < verbose.indent; n++)      
      printf(" ");
    printf("Subtracting average PA-swing from data\n");
  }

  angle_array = malloc(original->NrBins*sizeof(float));
  if(angle_array == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_subtractAveragePA: Memory allocation failed.");
    return 0;
  }


  //  if(!preprocess_addsuccessivepulses(*original, &profile, original->NrSubints, 0, verbose)) {
  if(preprocess_make_profile(*original, &profile, 0, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_subtractAveragePA: Cannot form a pulse profile from the input data.");
    return 0;
  }

  // Make everything the onpulse region. Note that we don't care about the error bars. We just want the PA as function of pulse longitude. So the selection is irrelevant.
  pulselongitude_regions_definition onpulse;
  if(initPulselongitudeRegion(&onpulse, verbose) == 0) {
    printerror(verbose.debug, "ERROR preprocess_subtractAveragePA: Initialising onpulse region failed.");
    return 0;
  }
  onpulse.nrRegions = 1;
  onpulse.bins_defined[0] = 1;
  onpulse.left_bin[0] = 0;
  onpulse.right_bin[0] = original->NrBins-1;

  if(make_paswing_fromIQUV(&profile, 0, onpulse, 0, 0, 1.0, 1.0, 1, 0.0, 0.0, NULL, 1.0, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_subtractAveragePA: Cannot form the PA swing of the average profile.");
    return 0;
  }

  if(readPulsePSRData(&profile, 0, 3, 0, 0, original->NrBins, angle_array, verbose) != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_subtractAveragePA: Extracting PA swing of the average profile failed.");
    return 0;
  }

  for(b = 0; b < original->NrBins; b++) {
    angle_array[b] =  -2.0*angle_array[b]*M_PI/180.0;
  }

  int ret;
  ret = preprocess_rotateStokes(original, clone, inplace, -1, 0, angle_array, 1, 2, verbose);
  if(ret == 0) {
    return 0;
  }

  free(angle_array);
  closePSRData(&profile, 0, verbose);
  freePulselongitudeRegion(&onpulse);
  return 1;
}



//START REGION RELEASE

/* Rotate stokes parameter stokes1 into stokes parameter stokes2. So
   if stokes1 = 1 and stokes2 = 2 Stokes Q is rotated into U by a
   certain angle in degrees. If angle_array is set to NULL, the value
   angle is used for all pulse longitude bins. angle should be in
   degrees, angle_array in radians. Otherwise, angle_array should
   point to an array specifying the angle for each pulse longitude
   bin. The clone should not exist yet. If inplace is set, the
   variable clone is ignored, and original is modified instead. Set
   subint to the subint you want to modify (counting from zero), or to
   -1 if you want to correct all subints with the specified
   angle. Verbose level determines nr of spaces before output. If
   nocounters is not set, the progresss is shown.

   Return 1 = success, 0 = error. */
int preprocess_rotateStokes(datafile_definition *original, datafile_definition *clone, int inplace, int subint, float angle, float *angle_array, int stokes1, int stokes2, verbose_definition verbose)
{
  int i;
  long f, n, b, p;
  float x, y;
  if(verbose.verbose) {
    for(n = 0; n < verbose.indent; n++)      
      printf(" ");
    printf("Rotating ");
    switch(stokes1) {
    case 0: printf("I"); break;
    case 1: printf("Q"); break;
    case 2: printf("U"); break;
    case 3: printf("V"); break;
    default: 
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_rotateStokes: Non-existing stokes parameter was defined.");
      return 0;
    }
    printf(" and ");
    switch(stokes2) {
    case 0: printf("I"); break;
    case 1: printf("Q"); break;
    case 2: printf("U"); break;
    case 3: printf("V"); break;
    default: 
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_rotateStokes: Non-existing stokes parameter was defined.");
      return 0;
    }
    if(angle_array == NULL) 
      printf(" by %f degrees\n", angle);
    else
      printf(" in a pulse longitude dependent way\n");
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rotateStokes: Rotating Stokes parameters only works if data is loaded into memory.");
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rotateStokes: Cannot handle PA data.");
    return 0;
  }
  if(original->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rotateStokes: There should be 4 polarization channels in the data.");
    return 0;
  }
  if(original->poltype == POLTYPE_COHERENCY) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_rotateStokes: Date is written as coherency parameters. Please convert into Stokes parameters first.");
    return 0;
  }else if(original->poltype == POLTYPE_UNKNOWN) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_rotateStokes: It is assumed the data are Stokes parameters.");
  }else if(original->poltype != POLTYPE_STOKES) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_rotateStokes: Cannot process this polarization type of the data.");
    return 0;
  }
  if(verbose.debug) {
    for(n = 0; n < verbose.indent; n++)      
      printf(" ");
    if(original->isDebase != 1) {
      printf("WARNING preprocess_rotateStokes: Baseline variations in time/freq will be introduced if baseline is not removed\n");
    }
  }

  if(inplace == 0) {
    cleanPSRData(clone, verbose);
    copy_params_PSRData(*original, clone, verbose);
    clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
    if(clone->data == NULL) {
      fflush(stdout); 
      printerror(verbose.debug, "ERROR preprocess_rotateStokes: Memory allocation error.");
      return 0;
    }
  }

  if(angle_array == NULL) { // Otherwise the array is assumed to be in radians.
    angle *= M_PI/180.0; 
  }
  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < original->NrFreqChan; f++) {
    for(n = 0; n < original->NrSubints; n++) {
      /*	if(writePulsePSRData(clone, n, p, f, 0, clone->NrBins, &(original->data[original->NrBins*(p+original->NrPols*(f+n*original->NrFreqChan))]), verbose) != 1) {
	return 0;
	}*/
      // Continue processing if this is the selected subint, if all subints need to be processed, or if clone needs to be generated
      if(n == subint || subint < 0 || inplace == 0) {
	for(b = 0; b < original->NrBins; b++) {
	  if(angle_array != NULL) {
	    angle = angle_array[b];
	  }
	  for(p = 0; p < 4; p++) {
	    /* Copy the two unaffected polarization channels to clone if using clone
	       or copy all polarization channels if this is not the selected subint
	    */
	    if((p != stokes1 && p != stokes2 && inplace == 0) || (inplace == 0 && subint != n)) {
	      clone->data[original->NrBins*(p+original->NrPols*(f+n*original->NrFreqChan))+b] = original->data[original->NrBins*(p+original->NrPols*(f+n*original->NrFreqChan))+b];
	    }else if(p == stokes1) {  // Copy two channels which are rotated
	      x = original->data[original->NrBins*(stokes1+original->NrPols*(f+n*original->NrFreqChan))+b]*cos(angle)-original->data[original->NrBins*(stokes2+original->NrPols*(f+n*original->NrFreqChan))+b]*sin(angle);
	      y = original->data[original->NrBins*(stokes1+original->NrPols*(f+n*original->NrFreqChan))+b]*sin(angle)+original->data[original->NrBins*(stokes2+original->NrPols*(f+n*original->NrFreqChan))+b]*cos(angle);
	      if(inplace == 0) {
		clone->data[original->NrBins*(stokes2+original->NrPols*(f+n*original->NrFreqChan))+b] = y;
		clone->data[original->NrBins*(stokes1+original->NrPols*(f+n*original->NrFreqChan))+b] = x;
	      }else {
		original->data[original->NrBins*(stokes2+original->NrPols*(f+n*original->NrFreqChan))+b] = y;
		original->data[original->NrBins*(stokes1+original->NrPols*(f+n*original->NrFreqChan))+b] = x;
	      }
	    }
	  }
	}
      }
    }
    if(verbose.verbose && verbose.nocounters == 0) {
      long doprint;
      doprint = 1;
      if(original->NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	doprint = 0;
      if(doprint) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  %.1f%%     \r", (100.0*(f*original->NrSubints+n))/(float)(original->NrSubints*original->NrFreqChan));
	fflush(stdout);
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(verbose.verbose) printf("  done                              \n");
  }
  return 1;
}

//START REGION DEVELOP

/* Create a clone in which the data outside range defined by onpulse
is thrown away.  The clone should not exist yet. Verbose level
determines nr of spaces before output. If nocounters is not set, the
progresss is shown. Return 1 = success */
int preprocess_gate(datafile_definition original, datafile_definition *clone, pulselongitude_regions_definition onpulse, verbose_definition verbose)
{
  int i;
  long p, f, n;
  region_frac_to_int(&onpulse, original.NrBins, 0);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Gating data");
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_gate: Gating only works if data is loaded into memory.");
    return 0;
  }
  //  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "ERROR preprocess_gate: Cannot handle PA data.");
  //    return 0;
  //  }
  if(onpulse.nrRegions == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_gate: No region is defined.");
    return 0;
  }
  if(onpulse.bins_defined[0] == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_gate: Region is not defined in bins.");
    return 0;
  }
  if(onpulse.left_bin[0] < 0 || onpulse.right_bin[0] >= original.NrBins || onpulse.nrRegions == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_gate: Specified range is not valid.");
    return 0;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  using bins %d - %d\n", onpulse.left_bin[0], onpulse.right_bin[0]);
  }

  cleanPSRData(clone, verbose);
  copy_params_PSRData(original, clone, verbose);
  clone->NrBins = onpulse.right_bin[0]-onpulse.left_bin[0]+1;
  clone->data = (float *)malloc((clone->NrBins)*(clone->NrPols)*(clone->NrFreqChan)*(clone->NrSubints)*sizeof(float));
  if(clone->data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_gate: Memory allocation error.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < clone->NrPols; p++) {
    for(f = 0; f < clone->NrFreqChan; f++) {
      for(n = 0; n < clone->NrSubints; n++) {
	if(writePulsePSRData(clone, n, p, f, 0, clone->NrBins, &(original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+onpulse.left_bin[0]]), verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR preprocess_gate: write error.");
	  return 0;
	}
      }
    }
    if(verbose.verbose && verbose.nocounters == 0) {
      long doprint;
      doprint = 1;
      if(clone->NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	doprint = 0;
      if(doprint) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  %.1f%%     \r", (100.0*((p*clone->NrFreqChan+f)*clone->NrSubints+n))/(float)(clone->NrSubints*clone->NrFreqChan*clone->NrPols));
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done                              \n");
  }
  return 1;
}

/* Smooth each pulse independent of each other using the fft
   algorithm. If highpass is set, the low frequencies are zapped
   rathern than the low frequencies. Verbose level determines nr of
   spaces before output. If nocounters is not set, the progresss is
   shown. Return 1 = success, 0 on error. */
int preprocess_fftSmooth(datafile_definition original, long nrZapFreqs, int highpass, verbose_definition verbose)
{
  int i;
  long p, f, n;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Smoothing: Zapping %ld/%ld fft frequency channels\n", nrZapFreqs, original.NrBins/2);
    verbose.indent += 2;
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_fftSmooth: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_fftSmooth: Cannot handle PA data.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	if(fftSmooth(&(original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))]), original.NrBins, highpass, nrZapFreqs, verbose) == 0)
	  return 0;
      }
    }
    if(verbose.verbose && verbose.nocounters == 0) {
      long doprint;
      doprint = 1;
      if(original.NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	doprint = 0;
      if(doprint) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  %.1f%%     \r", (100.0*((p*original.NrFreqChan+f)*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan*original.NrPols));
	fflush(stdout);
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("done                              \n");
  }
  return 1;
}

/* Zap each pulse independent of each other using fft
   algorithm. Verbose level determines nr of spaces before output. If
   nocounters is not set, the progresss is shown. Return 1 =
   success */
int preprocess_fftZap(datafile_definition original, long freq1, long freq2, verbose_definition verbose)
{
  int i;
  long p, f, n;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Zapping (%ld -> %ld)/%ld fft frequency channels\n", freq1, freq2, original.NrBins/2);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_fftZap: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_fftZap: Cannot handle PA data.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	if(fftZap(&(original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))]), original.NrBins, freq1, freq2, verbose) == 0)
	  return 0;
      }
    }
    if(verbose.verbose && verbose.nocounters == 0) {
      long doprint;
      doprint = 1;
      if(original.NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	doprint = 0;
      if(doprint) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  %.1f%%     \r", (100.0*((p*original.NrFreqChan+f)*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan*original.NrPols));
	fflush(stdout);
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(verbose.verbose) printf("  done                              \n");
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* If global=0: Normalizes each subint/channel independent of each
   other (based on first polarization channel).

   If global is set: Normalizes each subint/channel with the same
   scale factor based on the peak value found in the first
   polarization channel.

   If onpulse == NULL the whole pulse longitude range is searched for the
   peak value, otherwise only the onpulse region. Return 1 =
   success */
int preprocess_norm(datafile_definition original, float normvalue, pulselongitude_regions_definition *onpulse, int global, verbose_definition verbose)
{
  int ok, first, itt;
  long p, f, n, b, i;
  float max, fac, globalmax;
  pulselongitude_regions_definition onpulse_converted;
  if(onpulse != NULL) {
    //    memcpy(&onpulse_converted, onpulse, sizeof(regions_definition));
    if(initPulselongitudeRegion(&onpulse_converted, verbose) == 0) {
      printerror(verbose.debug, "ERROR preprocess_norm: Initialising onpulse region failed.");
      return 0;
    }
    copyPulselongitudeRegion(*onpulse, &onpulse_converted);
    region_frac_to_int(&onpulse_converted, original.NrBins, 0);
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(global)
      printf("Normalize data to %f (same scaling for each channels/subint)\n", normvalue);
    else
      printf("Normalize data to %f (individual channels/subints)\n", normvalue);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_norm: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_norm: Cannot handle PA data.");
    return 0;
  }
  
  // Define 2 itterations (only used when global is set): first itteration determines scale factor, second applies it.
  for(itt = 0; itt < 2; itt++) {
    if(global == 0)
      itt = 1;
    // So itt=0 means find global max first
    // itt=1 means apply single chan/subint max, or apply global max

    if(global && verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      if(itt == 0)
	printf("  Find normalisation contant\n");
      else
	printf("  Applying normalisation contant %e                   \n", fac);
    }

    /* Loop over pulses and frequency channels etc. */
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	p = 0;
	max = original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))];
	if(f == 0 && n == 0 && itt == 0)
	  globalmax = max;
	first = 1;
	for(b = 0; b < original.NrBins; b++) {
	  ok = 0;
	  if(onpulse == NULL) {
	    ok = 1;
	  }else {
	    if(checkRegions(b, &onpulse_converted, 0, verbose) != 0 || onpulse_converted.nrRegions == 0) {
	      ok = 1;
	    }
	  }
	  if((original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] > max || first == 1) && ok == 1) {
	    max = original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b];
	    first = 0;
	  }
	}
	if(itt == 0) {
	  if(max > globalmax) {
	    globalmax = max;
	    if(globalmax != 0)
	      fac = normvalue/globalmax;
	    else
	      fac = 1;
	  }
	}else if(global == 0) {
	  if(max != 0)
	    fac = normvalue/max;
	  else
	    fac = 1;
	}
	if(itt == 1) {
	  for(p = 0; p < original.NrPols; p++) {
	    for(b = 0; b < original.NrBins; b++) {
	      original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] *= fac;
	    }
	  }
	}
	if(verbose.verbose && verbose.nocounters == 0) {
	  long doprint;
	  doprint = 1;
	  if(original.NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	    doprint = 0;
	  if(doprint) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*(f*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan));
	    fflush(stdout);
	  }
	}
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done                              \n");
  }
  if(onpulse != NULL) {
    freePulselongitudeRegion(&onpulse_converted);
  }
  return 1;
}


/* 
   All values > clipvalue are set to the clipvalue. 
   All values < -clipvalue are set to -clipvalue. 
   Return 1 = success */
int preprocess_clip(datafile_definition original, float clipvalue, verbose_definition verbose)
{
  long p, f, n, b, i;
  float value;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Clip all samples exceeding %e (will be set to threshold value)\n", clipvalue);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_clip: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_clip: Clipping of data that represents PA values and error-bars might not be what you want.");
  }
  
  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < original.NrFreqChan; f++) {
    for(n = 0; n < original.NrSubints; n++) {
      for(p = 0; p < original.NrPols; p++) {
	for(b = 0; b < original.NrBins; b++) {
	  value = original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b];
	  if(value > clipvalue) {
	    original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] = clipvalue;
	  }else if(value < -clipvalue) {
	    original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] = -clipvalue;
	  }
	}
      }
      if(verbose.verbose && verbose.nocounters == 0) {
	long doprint;
	doprint = 1;
	if(original.NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	  doprint = 0;
	if(doprint) {
	  for(i = 0; i < verbose.indent; i++)      
	    printf(" ");
	  printf("  %.1f%%     \r", (100.0*(f*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan));
	  fflush(stdout);
	}
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done                              \n");
  }
 
  return 1;
}

//START REGION DEVELOP

/* Normalizes each subint/channel independent of each other (based on
   first polarization channel) by making the off-pulse rms =
   normvalue. If onpulse == NULL, the whole pulse longitude range is
   used to determine the rms, otherwise only the not selected region
   is used. Return 1 = success */
int preprocess_normRMS(datafile_definition original, float normvalue, pulselongitude_regions_definition *onpulse, verbose_definition verbose)
{
  int nodebase, polchan;
  long p, f, n, b, i;
  float fac, baseline, rms;
  pulselongitude_regions_definition onpulse_converted;

  if(onpulse != NULL) {
    if(initPulselongitudeRegion(&onpulse_converted, verbose) == 0) {
      printerror(verbose.debug, "ERROR preprocess_normRMS: Initialising onpulse region failed.");
      return 0;
    }
    copyPulselongitudeRegion(*onpulse, &onpulse_converted);
    //    memcpy(&onpulse_converted, onpulse, sizeof(regions_definition));
    region_frac_to_int(&onpulse_converted, original.NrBins, 0);
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Normalize offpulse rms to %f (individual channels/subints)\n", normvalue);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_normRMS: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_normRMS: Cannot handle PA data.");
    return 0;
  }

  nodebase = 0; // Subtract baseline when determining rms
  polchan = 0;  // Use first polarization channel

  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < original.NrFreqChan; f++) {
    for(n = 0; n < original.NrSubints; n++) {
      if(onpulse != NULL)
	offpulseStats(&(original.data[original.NrBins*(polchan+original.NrPols*(f+n*original.NrFreqChan))]), original.NrBins, &baseline, &rms, &onpulse_converted, nodebase, verbose);
      else
	offpulseStats(&(original.data[original.NrBins*(polchan+original.NrPols*(f+n*original.NrFreqChan))]), original.NrBins, &baseline, &rms, NULL, nodebase, verbose);
      if(rms != 0)
	fac = normvalue/rms;
      else
	fac = 1;
      for(p = 0; p < original.NrPols; p++) {
	for(b = 0; b < original.NrBins; b++) {
	  original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] *= fac;
	}
      }
      if(verbose.verbose && verbose.nocounters == 0) {
	long doprint;
	doprint = 1;
	if(original.NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	  doprint = 0;
	if(doprint) {
	  for(i = 0; i < verbose.indent; i++)      
	    printf(" ");
	  printf("  %.1f%%     \r", (100.0*(f*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan));
	  fflush(stdout);
	}
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done                              \n");
  }
  if(onpulse != NULL) {
    freePulselongitudeRegion(&onpulse_converted);
  }

  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Convert coherency parameters into stokes parameters. Verbose level
   determines nr of spaces before output.  Return 1 = success, 0 =
   error. */
int preprocess_stokes(datafile_definition *original, verbose_definition verbose)
{
  long f, n, b, i;
  float I, Q, U, V, c1, c2, c3, c4;
  int basis;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Forming Stokes parameters\n");
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: only works if data is loaded into memory.");
    return 0;
  }
  if(original->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Expected 4 polarization channels, only %ld present.", original->NrPols);
    return 0;
  }
  if(original->poltype == POLTYPE_STOKES) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_stokes: Data already in Stokes parameters, nothing will be done.");
    return 1;
  }else {
    if(original->poltype != POLTYPE_COHERENCY) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_stokes: Expected poltype=%d (coherency parameters), got %d.", POLTYPE_COHERENCY, original->poltype);
      return 0;
    }
  }
  if(original->feedtype == FEEDTYPE_LINEAR || original->feedtype == FEEDTYPE_INV_LINEAR) {
    basis = 1;
  }else if(original->feedtype == FEEDTYPE_CIRCULAR || original->feedtype == FEEDTYPE_INV_CIRCULAR) {
    basis = 2;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Expected feedtype=%d, this basis is not implemented.", original->feedtype);
    return 0;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(basis == 1)
      printf("  Using a linear basis\n");
    else if(basis == 2)
      printf("  Using a circular basis\n");
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Cannot handle PA data.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < original->NrFreqChan; f++) {
    for(n = 0; n < original->NrSubints; n++) {
      for(b = 0; b < original->NrBins; b++) {
	c1 = original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b];
	c2 = original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b];
	c3 = original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b];
	c4 = original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b];
	if(basis == 1) {  // Linear
	  I = c1+c2;
	  Q = c1-c2;
	  U = 2.0*c3;
	  V = 2.0*c4;
	}else if(basis == 2) {  // Circular
	  I = c1+c2;
	  Q = 2.0*c3;
	  U = 2.0*c4;
	  V = c1-c2;
	}
	original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b] = I;
	original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b] = Q;
	original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b] = U;
	original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b] = V;
      }
    }
  }
  original->poltype = POLTYPE_STOKES;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done                    \n");
  }
  return 1;
}

/* Convert Stokes parameters into coherency parameters. Verbose level
   determines nr of spaces before output.  Return 1 = success */
int preprocess_coherency(datafile_definition *original, verbose_definition verbose)
{
  /* If you change this function, also change preprocess_stokes() */
  long f, n, b, i;
  float I, Q, U, V, c1, c2, c3, c4;
  int basis;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Forming coherency parameters\n");
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: only works if data is loaded into memory.");
    return 0;
  }
  if(original->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Expected 4 polarization channels, only %ld present.", original->NrPols);
    return 0;
  }
  if(original->poltype != POLTYPE_STOKES) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Expected poltype=%d (Stokes parameters), got %d.", POLTYPE_STOKES, original->poltype);
    return 0;
  }
  if(original->feedtype == FEEDTYPE_LINEAR || original->feedtype == FEEDTYPE_INV_LINEAR) {
    basis = 1;
  }else if(original->feedtype == FEEDTYPE_CIRCULAR || original->feedtype == FEEDTYPE_INV_CIRCULAR) {
    basis = 2;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Expected feedtype=%d, this basis is not implemented.", original->feedtype);
    return 0;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(basis == 1)
      printf("  Using a linear basis\n");
    else if(basis == 2)
      printf("  Using a circular basis\n");
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_stokes: Cannot handle PA data.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < original->NrFreqChan; f++) {
    for(n = 0; n < original->NrSubints; n++) {
      for(b = 0; b < original->NrBins; b++) {
	I = original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b];
	Q = original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b];
	U = original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b];
	V = original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b];
	if(basis == 1) {  // Linear
	  /*
	    I = c1+c2;
	    Q = c1-c2;
	    U = 2.0*c3;
	    V = 2.0*c4;
	  */
          c1 = 0.5*(I+Q);
	  c2 = 0.5*(I-Q);
	  c3 = 0.5*U;
	  c4 = 0.5*V;
	}else if(basis == 2) {  // Circular
	  /*
	    I = c1+c2;
	    Q = 2.0*c3;
	    U = 2.0*c4;
	    V = c1-c2;
	  */
          c1 = 0.5*(I+V);
	  c2 = 0.5*(I-V);
	  c3 = 0.5*Q;
	  c4 = 0.5*U;
	}
	original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b] = c1;
	original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b] = c2;
	original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b] = c3;
	original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b] = c4;
      }
    }
  }
  original->poltype = POLTYPE_COHERENCY;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done                      \n");
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Multiplies data values with a certain value after applying an
   offset. Verbose level determines nr of spaces before output. If
   nocounters is not set, the progresss is shown.  Return 1 =
   success */
int preprocess_scale(datafile_definition original, float factor, float offset, verbose_definition verbose)
{
  long p, f, n, b, i;
  float sample;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Scale data with %f and offset %f\n", factor, offset);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_scale: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_scale: Cannot handle PA data.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < original.NrFreqChan; f++) {
    for(n = 0; n < original.NrSubints; n++) {
      for(p = 0; p < original.NrPols; p++) {
	for(b = 0; b < original.NrBins; b++) {
	  sample = original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b];
	  original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] = (sample+offset)* factor;
	}
      }
      if(verbose.verbose && verbose.nocounters == 0) {
	long doprint;
	doprint = 1;
	if(original.NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	  doprint = 0;
	if(doprint) {
	  for(i = 0; i < verbose.indent; i++)      
	    printf(" ");
	  printf("  %.1f%%     \r", (100.0*(f*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan));
	  fflush(stdout);
	}
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done                              \n");
  }
  return 1;
}

//START REGION DEVELOP

/* Replaces data with fft. Verbose level determines nr of spaces
   before output. If nocounters is not set, the progresss is
   shown. Return 1 = success */
int preprocess_fft(datafile_definition *original, verbose_definition verbose)
{
  int i;
  long p, f, n;
  pulselongitude_regions_definition onpulse;
  datafile_definition clone;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Calculating fft of data\n");
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_fft: only works if data is loaded into memory.");
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_fft: Cannot handle PA data.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original->NrPols; p++) {
    for(f = 0; f < original->NrFreqChan; f++) {
      for(n = 0; n < original->NrSubints; n++) {
	if(replacefft(&(original->data[original->NrBins*(p+original->NrPols*(f+n*original->NrFreqChan))]), original->NrBins, verbose) == 0)
	  return 0;
	if(verbose.verbose && verbose.nocounters == 0) {
	  long doprint;
	  doprint = 1;
	  if(original->NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	    doprint = 0;
	  if(doprint) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*((p*original->NrFreqChan+f)*original->NrSubints+n))/(float)(original->NrSubints*original->NrFreqChan*original->NrPols));
	    fflush(stdout);
	  }
	}
      }
    }
  }

  if(initPulselongitudeRegion(&onpulse, verbose) == 0) {
    printerror(verbose.debug, "ERROR preprocess_fft: Initialising onpulse region failed");
    return 0;
  }
  onpulse.nrRegions = 1;
  onpulse.left_bin[0] = 0;
  onpulse.right_bin[0] = original->NrBins/2;
  onpulse.bins_defined[0] = 1;
  if(preprocess_gate(*original, &clone, onpulse, verbose) == 0) {
    printerror(verbose.debug, "ERROR preprocess_fft: gating failed");
    return 0;
  }
  swap_orig_clone(original, &clone, verbose);

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done                              \n");
  }
  freePulselongitudeRegion(&onpulse);
  return 1;
}

/* Replaces signal with a test signal. Return 1 = success. Verbose
   level determines nr of spaces before output. If nocounters is not
   set, the progresss is shown. */
int preprocess_testsignal(datafile_definition original, verbose_definition verbose)
{
  long p, f, n, b;
  int i;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Injecting test signal\n");
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_testsignal: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_testsignal: Cannot handle PA data.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	for(b = 0; b < original.NrBins; b++) {
	  original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))+b] = (f/(float)original.NrFreqChan+1)*(sin(M_PI*b/(float)(original.NrBins)))+0.1*n;
	}
	if(verbose.verbose && verbose.nocounters == 0) {
	  long doprint;
	  doprint = 1;
	  if(original.NrFreqChan > 4 && n != 0)  // Avoid excessive nr of prints
	    doprint = 0;
	  if(doprint) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*((p*original.NrFreqChan+f)*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan*original.NrPols));
	    fflush(stdout);
	  }
	}
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(verbose.verbose) printf("done                              \n");
  }

  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE


/* Rotates each pulse independent of each other in Fourier space based
   on the expected dispersive delay w.r.t. the reference frequency
   specified in datafile_definition. If undo is set, the dispersion
   delay is re-introduced into the data if the data was de-dispersed
   previously. Otherwise, if update is zero, the nothing is done if
   the data is already dedispersed. If update is nonzero, dedispersion
   will be applied again to change the reference frequency from
   freq_ref (-1 or 1e10 = inf freq) to the reference frequency
   specified in datafile_definition.

   Return 0 = error, 1 = command was ignored because of warnings (such
   as data already being dedispersed), 2 = dispersive delay was taken
   out. Verbose level determines nr of spaces before output. */
int preprocess_dedisperse(datafile_definition *original, int undo, int update, double freq_ref, verbose_definition verbose)
{
  long p, f, n;
  int i, inffreq, inffreq_old;
  long double dt, dt_samples;
  double freq;

  if(undo && update) {
    printerror(verbose.debug, "ERROR preprocess_dedisperse (%s): Cannot update the reference frequency and re-dedisperse the data simultaneously.", original->filename);    
    return 0;
  }

  if(original->freq_ref < -1.1) {
    if(undo == 0) {
      printwarning(verbose.debug, "WARNING preprocess_dedisperse (%s): Reference frequency is unknown. The reference frequency is set to infinite frequency.", original->filename);
      original->freq_ref = 1e10;
    }else {
      printerror(verbose.debug, "ERROR preprocess_dedisperse (%s): Cannot re-dedisperse the data if the reference frequency of the current dedispersion is unknown.", original->filename);    
      return 0;
    }
  }

  if((original->freq_ref > -1.1 && original->freq_ref < -0.9) || (original->freq_ref > 0.99e10 && original->freq_ref < 1.01e10))
    inffreq = 1;
  else
    inffreq = 0;

  if((freq_ref > -1.1 && freq_ref < -0.9) || (freq_ref > 0.99e10 && freq_ref < 1.01e10))
    inffreq_old = 1;
  else if(freq_ref < 0) {
    printerror(verbose.debug, "ERROR preprocess_dedisperse (%s): Requested frequency is invalid (%f).", original->filename, freq_ref);    
    return 0;
  }else
    inffreq_old = 0;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(undo == 0) {
      printf("De-dispersing %s with DM=%f with reference frequency=", original->filename, original->dm);
    }else {
      printf("Re-dedispersing %s with DM=%f with reference frequency=", original->filename, original->dm);
    }
    if(inffreq)
      printf("infinity\n");
    else
      printf("%f MHz\n", original->freq_ref);
  }

  if(original->isDeDisp == 1 && undo == 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Data in %s is already dedispersed\n", original->filename);
    }
  }else if(original->isDeDisp == 0 && undo) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Data in %s is already dispersed\n", original->filename);
    }
  }else if(original->isDeDisp == -1) {
    fflush(stdout);
    if(undo == 0) {
      printwarning(verbose.debug, "WARNING preprocess_dedisperse (%s): unknown dedispersion state. Data is assumed to be already dedispersed.", original->filename);
    }else {
      printwarning(verbose.debug, "WARNING preprocess_dedisperse (%s): unknown dedispersion state. Data is assumed to be already dispersed.", original->filename);
    }
  }
  if(undo == 0) {
    if(original->isDeDisp == 1 || original->isDeDisp == -1) {
      if(update) {
	if(verbose.verbose) {
	  for(i = 0; i < verbose.indent; i++)      
	    printf(" ");
	  printf("  Updating reference frequency from ");
	  if(inffreq_old)
	    printf("infinity");
	  else
	    printf("%f MHz", freq_ref);
	  if(inffreq)
	    printf(" to infinity\n");
	  else
	    printf(" to %f MHz\n", original->freq_ref);
	}
	if((inffreq == 1 && inffreq_old == 1) || fabs(original->freq_ref - freq_ref) < 0.001) {
	  if(verbose.verbose) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("  Reference frequency is the same, nothing done\n");
	  }
	  return 1;
	}
      }else {
	return 1;
      }
    }else {
      update = 0;  // Make sure update is only set when update reference frequency needs to be done
    }
  }else {   // if undo is set
    if(original->isDeDisp == 0 || original->isDeDisp == -1) {
      return 1;
    }
  }

  dt_samples = 0; // To avoid compiler warning
  // If updating reference frequency, there is an overall freq. independent shift to be applied.
  if(update) {
    freq = get_weighted_channel_freq(*original, 0, 0, verbose);   // Doesn't really matter, pick a frequency
    dt = -calcDMDelay(freq, freq_ref, inffreq_old, original->dm); // To undo old dedispersion
    dt += calcDMDelay(freq, original->freq_ref, inffreq, original->dm); // To apply new dedispersion
    dt_samples = dt/get_tsamp(*original, 0, verbose);
  }

  if(fabs(original->dm) < 1e-6) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_dedisperse (%s): DM appears to be not set", original->filename);
    return 1;
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    f = original->NrFreqChan/2;
    freq = get_weighted_channel_freq(*original, 0, f, verbose);
    if(update == 0) {
      if(undo == 0) {
	printf("  Rotating central frequency channel of first subint (%lf MHz) by %Lf phase\n", freq, calcDMDelay(freq, original->freq_ref, inffreq, original->dm)/(original->NrBins*get_tsamp(*original, 0, verbose)));
      }else {
	printf("  Rotating central frequency channel of first subint (%lf MHz) by %Lf phase\n", freq, -calcDMDelay(freq, original->freq_ref, inffreq, original->dm)/(original->NrBins*get_tsamp(*original, 0, verbose)));
      }
    }else {
      printf("  Rotating central frequency channel of first subint (%lf MHz) by %Lf phase\n", freq, dt_samples/(double)(original->NrBins));
    }
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_dedisperse (%s): only works if data is loaded into memory.", original->filename);
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_dedisperse (%s): Cannot handle PA data.", original->filename);
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original->NrPols; p++) {
    for(f = 0; f < original->NrFreqChan; f++) {
      for(n = 0; n < original->NrSubints; n++) {
	if(update == 0) {
	  dt = calcDMDelay(get_weighted_channel_freq(*original, n, f, verbose), original->freq_ref, inffreq, original->dm);
	  if(undo)
	    dt *= -1.0;
	  dt /= get_tsamp(*original, 0, verbose);   /* Convert shift from seconds in samples */
	}else {
	  dt = dt_samples;
	}
	if(rotateSinglepulse(&(original->data[original->NrBins*(p+original->NrPols*(f+n*original->NrFreqChan))]), original->NrBins, -dt, verbose) == 0)
	  return 0;
      }
    }
  }
  if(undo == 0) {
    original->isDeDisp = 1;
  }else {
    original->isDeDisp = 0;
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done              \n");
  }
  return 2;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Compensate for the effect of Faraday rotation with the reference
   frequency determined in datafile_definition. If undo is set, a
   Faraday rotation is applied rather than corrected for. If update is
   set, the reference frequency is changed from freq_ref (-1 = inf
   freq) to what is set in datafile_definition. Cannot do an undo and
   update simultaneously. Verbose determines nr of spaces before
   output.

   If rm_table is set to NULL, the RM from the header is
   used. Otherwise rm_table is interpretted as a list of RM's for each
   pulse longitude bin.

   Return:
     0 = error, 
     1 = ignored because of warning (i.e. already de-Faraday rotated),
     2 = applied de-Faraday rotation. 
*/
int preprocess_deFaraday(datafile_definition *original, int undo, int update, double freq_ref, double *rm_table, verbose_definition verbose)
{
  long f, n, b;
  int i, inffreq, inffreq_old;
  float dphi, phi, L, *pulseQ, *pulseU;
  verbose_definition verbose2;

  if(original->freq_ref < -1.1) {
    printwarning(verbose.debug, "WARNING preprocess_deFaraday (%s): Reference frequency is unknown. The reference frequency is set to infinite frequency.", original->filename);
    original->freq_ref = 1e10;
  }

  if((original->freq_ref > -1.1 && original->freq_ref < -0.9) || (original->freq_ref > 0.99e10 && original->freq_ref < 1.01e10))
    inffreq = 1;
  else
    inffreq = 0;

  if((freq_ref > -1.1 && freq_ref < -0.9) || (freq_ref > 0.99e10 && freq_ref < 1.01e10))
    inffreq_old = 1;
  else if(freq_ref < 0) {
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Requested original reference frequency is invalid (%f).", original->filename, freq_ref);
    return 0;
  }else
    inffreq_old = 0;

  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(undo == 0) {
      if(rm_table == NULL)
	printf("De-Faraday rotate frequency channels of %s with RM=%lf with reference frequency=", original->filename, original->rm);
      else
	printf("De-Faraday rotate frequency channels of %s with RM's from a table with reference frequency=", original->filename);
      if(inffreq)
	printf("infinity\n");
      else
	printf("%f MHz\n", original->freq_ref);
    }else {
      if(rm_table == NULL) 
	printf("Undo de-Faraday rotation of frequency channels of %s with RM=%lf\n", original->filename, original->rm);
      else
	printf("Undo de-Faraday rotation of frequency channels of %s with RM's from a table\n", original->filename);
    }
  }
  if(update && undo) {
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Cannot undo Faraday rotation and update reference frequency simultaneously.", original->filename);
    return 0;
  }

  if(original->isDeFarad == 1 && undo == 0) {
    if(update == 0) {
      if(verbose.verbose) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  Data in %s is already de-Faraday rotated, will not correct data.\n",  original->filename);
      }
      return 1;
    }else {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Updating reference frequency from ");
      if(inffreq_old)
	printf("infinity");
      else
	printf("%f MHz", freq_ref);
      if(inffreq)
	printf(" to infinity\n");
      else
	printf(" to %f MHz\n", original->freq_ref);
      if((inffreq == 1 && inffreq_old == 1) || fabs(original->freq_ref - freq_ref) < 0.001) {
	if(verbose.verbose) {
	  for(i = 0; i < verbose.indent; i++)      
	    printf(" ");
	  printf("  Reference frequency is the same, nothing done\n");
	}
	return 1;
      }
    }
  }
  if(original->isDeFarad == 0 && update) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Data in %s is not de-Faraday rotated, cannot update reference frequency.\n",  original->filename);
      return 1;
    }
  }
  if(original->isDeFarad == 0 && undo) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    printf("  Data in %s is not de-Faraday rotated, will not undo correction.\n",  original->filename);
    return 1;
  }
  if(original->isDeFarad == -1) {
    if(undo == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_deFaraday (%s): unknown de-Faraday rotation state. Data will not be de-Faraday rotated.", original->filename);
    }else {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_deFaraday (%s): unknown de-Faraday rotation state. De-Faraday rotation will not be undone.", original->filename);
    }
    return 1;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Remove Faraday rotation from data which still contains the four Stokes parameters.", original->filename);
    return 0;
  }
  if(original->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Expected four polarization channels.", original->filename);
    return 0;
  }
  if(original->poltype != POLTYPE_STOKES) {
    fflush(stdout);
    //    printwarning(verbose.debug, "WARNING preprocess_deFaraday (%s): Expected Stokes parameters, but poltype = %d != %d. Trying to convert data to Stokes parameters.", original->filename, original->poltype, POLTYPE_STOKES);

    if(preprocess_stokes(original, verbose2) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Conversion into Stokes parameters failed.", original->filename);
      return 0;
    }

    if(original->poltype != POLTYPE_STOKES) {
      printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Conversion into Stokes parameters failed.", original->filename);
      return 0;
    }
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): only works if data is loaded into memory.", original->filename);
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Cannot handle PA data.", original->filename);
    return 0;
  }
  if(fabs(original->rm) < 1e-3) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_deFaraday (%s): RM extremely small (%e). Is it really set? De-Faraday rotation is not applied.", original->filename, original->rm);
    return 1;
  }
  if(verbose.debug) {
    for(n = 0; n < verbose.indent; n++)      
      printf(" ");
    if(original->isDebase != 1) {
      printf("WARNING preprocess_deFaraday: Baseline variations in time/freq will be introduced if baseline is not removed\n");
    }
  }

  dphi = 0; // To avoid compiler warning
  if(update) {  // If updating reference frequency, there is an overall freq. independent shift to be applied. So here this overall shift is calculated.
    if(rm_table != NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Cannot update the reference frequency when the RM is specified speperately for each pulse longitude bin.", original->filename);
      return 0;
    }

    double freq;
    freq = get_weighted_channel_freq(*original, 0, 0, verbose); // Frequency doesn't matter, so just pick one.
    dphi = -calcRMAngle(freq, freq_ref, inffreq_old, original->rm);   // To undo previous correction
    dphi += calcRMAngle(freq, original->freq_ref, inffreq, original->rm);  // To apply new correction
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(rm_table == NULL) {
      if(update == 0)
	dphi = calcRMAngle(get_weighted_channel_freq(*original, 0, original->NrFreqChan/2, verbose), original->freq_ref, inffreq, original->rm);
      printf("  Rotating Q&U of central frequency channel of first subint (%f MHz) by %f rad = %f deg\n", get_weighted_channel_freq(*original, 0, original->NrFreqChan/2, verbose), 2.0*dphi, 2.0*dphi*180.0/M_PI);
    }else {
      printf("  Rotating Q&U\n");
    }
  }

  pulseQ = (float *)malloc(original->NrBins*sizeof(float));
  pulseU = (float *)malloc(original->NrBins*sizeof(float));
  if(pulseQ == NULL || pulseU == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Cannot allocate memory.", original->filename);
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < original->NrFreqChan; f++) {
    for(n = 0; n < original->NrSubints; n++) {
      if(update == 0) {
	if(rm_table == NULL) // Deal with the bin-depedent RM below
	  dphi = calcRMAngle(get_weighted_channel_freq(*original, n, f, verbose), original->freq_ref, inffreq, original->rm);
	if(undo)
	  dphi *= -1.0;
      }
      //    printf("XXXXX %ld: %f\n", f, dphi);
      if(readPulsePSRData(original, n, 1, f, 0, original->NrBins, pulseQ, verbose) != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Read error.", original->filename);
	return 0;
      }
      if(readPulsePSRData(original, n, 2, f, 0, original->NrBins, pulseU, verbose) != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Read error.", original->filename);
	return 0;
      }
      for(b = 0; b < original->NrBins; b++) {
	if(rm_table != NULL) {
	  if(update == 0) {
	    dphi = calcRMAngle(get_weighted_channel_freq(*original, n, f, verbose), original->freq_ref, inffreq, rm_table[b]);
	    if(undo)
	      dphi *= -1.0;
	  }
	}

	L = sqrt(pulseQ[b]*pulseQ[b]+pulseU[b]*pulseU[b]);
	phi = atan2(pulseU[b], pulseQ[b]);
	phi -= 2.0*dphi;
	pulseQ[b] = L*cos(phi);
	pulseU[b] = L*sin(phi);
      }
      if(writePulsePSRData(original, n, 1, f, 0, original->NrBins, pulseQ, verbose) != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Write error.", original->filename);
	return 0;
      }
      if(writePulsePSRData(original, n, 2, f, 0, original->NrBins, pulseU, verbose) != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR preprocess_deFaraday (%s): Write error.", original->filename);
	return 0;
      }
    }
  }
  if(undo == 0)
    original->isDeFarad = 1;
  else
    original->isDeFarad = 0;
  free(pulseQ);
  free(pulseU);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done              \n");
  }
  return 2;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Compensate for the effect of parallactic angle. If clone is set to
   NULL, the original is modified (if read into memory). Otherwise a
   clone is generated. Verbose level determines nr of spaces before
   output.  If nocounters is not set, the progresss is shown. If undo
   is set, the data is distorted rather than corrected.

   Return 0 = error, 
   1 = command is ignored because already parallactic angle corrected, 
   2 = ignored because of warnings (such as unknown location of telescope), 
   3 = parallactic angle got changed. 

   Note that a clone is generated or data is modified only when the
   return code is 3.
  */
int preprocess_corrParAng(datafile_definition *original, datafile_definition *clone, int undo, verbose_definition verbose)
{
  int i;
  long n;
  double parang, parang2;
  verbose_definition verbose2;

  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("De-rotate data to ");
    if(undo)
      printf("add the effect of parallactic angle\n");
    else
      printf("remove the effect of parallactic angle\n");
  }
  if(original->isDePar == 1 && undo == 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    printf("  Data in %s is already parallactic angle corrected, so cannot do correction.\n", original->filename);
    return 1;
  }
  if(original->isDePar == 0 && undo) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    printf("  Data in %s is not parallactic angle corrected, so cannot undo correction.\n", original->filename);
    return 1;
  }
  if(original->isDePar == -1 && undo == 0) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): unknown parallactic angle state. Parallactic angle correction is ignored.", original->filename);
    return 2;
  }
  if(original->isDePar == -1 && undo) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): unknown parallactic angle state. Will not undo Parallactic angle correction.", original->filename);
    return 2;
  }

  if(verbose.verbose) {
    if(data_parang(*original, -1, &parang, verbose) == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): Cannot calculate parallactic angle, ignoring parallactic angle correction.", original->filename);
      return 2;
    }
    parang *= 180.0/M_PI;
    if(undo)
      parang *= -1.0;
    for(i = 0; i < verbose.indent; i++)
      printf(" ");
    printf("  Correcting %s for a parallactic angle (which is %f deg for midpoint of observation)\n", original->filename, parang);
    if(original->NrSubints > 1) {
      if(data_parang(*original, 0, &parang, verbose) == 0) {
	fflush(stdout);
	printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): Cannot calculate parallactic angle, ignoring parallactic angle correction.", original->filename);
	return 2;
      }
      if(data_parang(*original, original->NrSubints-1, &parang2, verbose) == 0) {
	fflush(stdout);
	printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): Cannot calculate parallactic angle, ignoring parallactic angle correction.", original->filename);
	return 2;
      }
      parang = parang2 - parang;
      parang *= 180.0/M_PI;
      if(undo)
	parang *= -1.0;
      for(i = 0; i < verbose.indent; i++)
	printf(" ");
      printf("  Total change in parallactic angle throughout observation is %f deg.\n", parang);
    }
  }

  // If requested, make a clone which has all data in memory
  if(clone != NULL) {
    if(make_clone(*original, clone, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_corrParAng: Cannot make clone of data, cannot remove effect of parallactic angle.");
      return 0;
    }
  }


  // If data are coherency parameters, convert into Stokes parameters
  if(original->poltype != POLTYPE_STOKES) {
    if(clone != NULL) {
      if(preprocess_stokes(clone, verbose2) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR preprocess_corrParAng (%s): Cannot form the Stokes parameters.", original->filename);
	return 0;
      }
    }else {
      if(preprocess_stokes(original, verbose2) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR preprocess_corrParAng (%s): Cannot form the Stokes parameters.", original->filename);
	return 0;
      }
    }
  }

  for(n = 0; n < original->NrSubints; n++) {
    if(data_parang(*original, n, &parang, verbose) == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING preprocess_corrParAng (%s): Cannot calculate parallactic angle, ignoring parallactic angle correction.", original->filename);
      return 2;
    }
    parang *= 180.0/M_PI;
    if(undo)
      parang *= -1.0;
    if(verbose.verbose && original->NrSubints > 1 && n == 0) {
      for(i = 0; i < verbose.indent; i++)
	printf(" ");
      printf("For first subint:\n");
    }
    if(n != 0) { // Don't output angle for each subint
      verbose2.verbose = 0;
      verbose2.nocounters = 0;
    }
    if(clone != NULL) {
      if(preprocess_rotateStokes(clone, NULL, 1, n, 2.0*parang, NULL, 1, 2, verbose2) != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR preprocess_corrParAng (%s): Cannot rotate Stokes Q and U.", original->filename);
	return 0;
      }
    }else {
      if(preprocess_rotateStokes(original, NULL, 1, n, 2.0*parang, NULL, 1, 2, verbose2) != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR preprocess_corrParAng (%s): Cannot rotate Stokes Q and U.", original->filename);
	return 0;
      }
    }
  }
  if(undo) {
    if(clone != NULL)
      clone->isDePar = 0;
    else
      original->isDePar = 0;
  }else {
    if(clone != NULL)
      clone->isDePar = 1;
    else
      original->isDePar = 1;
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done          \n");
  }

  return 3;
}

//START REGION DEVELOP

/* Change the coherency parameters in such a way that a cable-swap is
   undone. Return 1 = success */
int preprocess_swapcables(datafile_definition *original, verbose_definition verbose)
{
  int swaptoStokes;
  long f, n, b, i;
  float c1, c2, c3, c4, c1new, c2new, c3new, c4new;
  swaptoStokes = 0;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Undo a cable swap\n");
  }
  if(original->format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_swapcables: only works if data is loaded into memory.");
    return 0;
  }
  if(original->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_swapcables: Expected 4 polarization channels, only %ld present.", original->NrPols);
    return 0;
  }
  if(original->poltype == POLTYPE_ILVPAdPA || original->poltype == POLTYPE_PAdPA || original->poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_swapcables: Cannot handle PA data.");
    return 0;
  }
  /* Got Stokes parameters, get it into coherency parameters */
  if(original->poltype == POLTYPE_STOKES) {
    if(preprocess_coherency(original, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_swapcables: Cannot convert Stokes into coherency parameters.");
      return 0;
    }
    swaptoStokes = 1;
  }
  if(original->poltype != POLTYPE_COHERENCY) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_swapcables: Expected poltype=%d (coherency parameters), got %d.", POLTYPE_COHERENCY, original->poltype);
    return 0;
  }

  if(original->cableSwap == -1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_swapcables: State of cables during observations is unknown.");
  }else if(original->cableSwap == 0) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_swapcables: Header suggests that cables were not swapped during the observation.");
  }
  if(original->cableSwapcor == 1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_swapcables: Header suggests that a cable swap was already applied.");
  }
  if(original->isDeFarad == 1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING preprocess_swapcables: Header suggests that Faraday rotation has already removed. This is probably unlikely to be a good idea when applying a cable swap.");    
  }
  original->cableSwapcor = 1;

  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < original->NrFreqChan; f++) {
    for(n = 0; n < original->NrSubints; n++) {
      for(b = 0; b < original->NrBins; b++) {
	c1 = original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b];
	c2 = original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b];
	c3 = original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b];
	c4 = original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b];
	/* I think the following is true:
	   c1 = XX*            -> YY*            = c2
	   c2 = YY*            -> XX*            = c1
	   c3 = (XY* + YX*)/2  -> (YX* + XY*)/2  = c3
	   c4 = (YX* - XY*)/2i -> (XY* - YX*)/2i = -c4;

Lin: Stokes V swaps sign and Q as well (PA-swing reverses)
	  I = c1+c2;
	  Q = c1-c2;
	  U = 2.0*c3;
	  V = 2.0*c4;

Circ: Stokes V swaps sign and U as well (PA-swing reverses)
	  I = c1+c2;
	  Q = 2.0*c3;
	  U = 2.0*c4;
	  V = c1-c2;
	*/
	c1new = c2;
	c2new = c1;
	c3new = c3;
	c4new = -c4;
	original->data[original->NrBins*(0+original->NrPols*(f+n*original->NrFreqChan))+b] = c1new;
	original->data[original->NrBins*(1+original->NrPols*(f+n*original->NrFreqChan))+b] = c2new;
	original->data[original->NrBins*(2+original->NrPols*(f+n*original->NrFreqChan))+b] = c3new;
	original->data[original->NrBins*(3+original->NrPols*(f+n*original->NrFreqChan))+b] = c4new;
      }
    }
  }

  if(swaptoStokes) {
    if(preprocess_stokes(original, verbose) != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_swapcables: Cannot convert coherency parameters into Stokes parameters.");
      return 0;
    }
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  done\n");
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Using the original data, construct a profile. Dedispersion,
   de-Faraday rotation etc is taken care off. If stokesI is set,
   Stokes parameters are formed and only StokesI is stored in the
   profile.  The profile should not exist yet. Return 1 = success, 0
   on error. Verbose level determines nr of spaces before output. If
   nocounters is not set, the progresss is shown. */
int preprocess_make_profile(datafile_definition original, datafile_definition *profile, int stokesI, verbose_definition verbose)
{
  int i;
  datafile_definition clone, clone2;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Constructing ");
    if(stokesI)
      printf("Stokes I ");
    printf("profile\n");
    verbose.indent += 2;
  }

  if(make_clone(original, &clone, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_make_profile: Cannot make clone of data, cannot construct profile");
    return 0;
  }

  if(original.NrPols == 4 && stokesI) {
    if(original.poltype != POLTYPE_STOKES) {
      if(preprocess_stokes(&clone, verbose) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR preprocess_make_profile: forming Stokes parameters failed, cannot construct profile");
	return 0;
      }
    }
    if(original.NrPols > 1) {
      //  fflush(stdout); fprintf(stderr, "XXXXXX\n");
      if(preprocess_polselect(clone, &clone2, 0, verbose) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR preprocess_make_profile: Selecting Stokes I failed, cannot construct profile");
	return 0;
      }
      //      fflush(stdout); fprintf(stderr, "XXXXXX %d %d\n", clone.dumpOnClose, clone2.dumpOnClose);
      swap_orig_clone(&clone, &clone2, verbose);
      //  fflush(stdout); fprintf(stderr, "XXXXXX\n");
    }
  }

  if(clone.NrFreqChan > 1) {
    if(preprocess_dedisperse(&clone, 0, 0, 0, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_make_profile: de-dispersing failed, cannot construct profile");
      return 0;
    }

    if(preprocess_addsuccessiveFreqChans(clone, &clone2, clone.NrFreqChan, NULL, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_make_profile: summing frequency channels failed, cannot construct profile");
      return 0;
    }
    swap_orig_clone(&clone, &clone2, verbose);
  }

  if(preprocess_addsuccessivepulses(clone, profile, clone.NrSubints, 0, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_make_profile: summing subints failed, cannot construct profile");
    return 0;
  }

  closePSRData(&clone, 0, verbose);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Constructing profile done\n");
  }
  return 1;
}

//START REGION DEVELOP

/* Smooth the data by convolving each subint/freq. channel in Fourier
   space with a tophat function with a full width width (in
   bins). Return 0 = error, 1 = success. Verbose level determines nr
   of spaces before output. */
int preprocess_smooth(datafile_definition original, float width, verbose_definition verbose)
{
  long p, f, n;
  int i;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Smoothing by convolving data with tophat function of width %f\n", width);
    verbose.indent += 2;
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_smooth (%s): only works if data is loaded into memory.", original.filename);
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_smooth (%s): Cannot handle PA data.", original.filename);
    return 0;
  }

  verbose_definition noverbose;
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;

  /* Loop over pulses and frequency channels etc. */
  for(p = 0; p < original.NrPols; p++) {
    for(f = 0; f < original.NrFreqChan; f++) {
      for(n = 0; n < original.NrSubints; n++) {
	if(tophat_smooth(&(original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))]), original.NrBins, width, noverbose) == 0)
	  return 0;
      }
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("done              \n");
  }
  return 1;
}

void get_mean_rms(float *values, int *bins_selected, int nrbins, float *mean, float *rms)
{
  int i, nrselected;
  double av, av2;
  av = 0;
  av2 = 0;
  nrselected = 0;
  for(i = 0; i < nrbins; i++) {
    if(bins_selected[i]) {
      av += (double)values[i];
      av2 += (double)values[i]*(double)values[i];
      nrselected++;
    }
  }
  if(nrselected != 0) {
    av /= (double)nrselected;
    av2 /= (double)nrselected;
    *mean = av;
    *rms = sqrt(av2-av*av);
    //    printf("XXXXX %lf %lf %lf %f\n", av2, av, av2-av*av, *rms);
  }else {
    *mean = 0;
    *rms = 1;
  }
}

int widest_notselected(int *bins_selected, int nrbins)
{
  int i, widest, width, status;
  widest = 0;
  width = 0;
  status = 0;
  for(i = 0; i < nrbins; i++) {
    if(bins_selected[i] == 0 && status == 0) {   // New region
      width = 1;
      status = 1;
    }else if(bins_selected[i] == 0 && status == 1) {   // continuing region
      width++;
    }else if(bins_selected[i] == 1 && status == 1) {   // stopping region
      status = 0;
    }else if(bins_selected[i] == 1 && status == 0) {
    }
    if(width > widest)
      widest = width;
  }
  if(status) {    // If ending while counting, continue from start until 
    for(i = 0; i < nrbins; i++) {
      if(bins_selected[i] == 0) {   // continuing region
	width++;
      }else {
	break;
      }
      if(width > widest)
	widest = width;
    }
  }
  return widest;
}

void remove_narrow_notselected(int *bins_selected, int nrbins, int widest, int minwidthbins)
{
  int i, j, width, status;

  if(minwidthbins > widest)
    minwidthbins = widest - 1;

  status = 0;
  for(i = 0; i < nrbins; i++) {
    if(bins_selected[i] == 0 && status == 0) {   // New region
      width = 1;
      status = 1;
    }else if(bins_selected[i] == 0 && status == 1) {   // continuing region
      width++;
    }else if(bins_selected[i] == 1 && status == 1) {   // stopping region
      status = 0;
      //      printf("XXXXX %d <= %d for i=%d?\n", width, minwidthbins, i);
      if(width <= minwidthbins) {   // Remove region
	for(j = i-1; j >= 0; j--) {
	  if(bins_selected[j] == 0)
	    bins_selected[j] = 1;
	  else
	    break;
	}
      }
    }else if(bins_selected[i] == 1 && status == 0) {
    }
  }
  if(status) {    // If ending while counting, continue from start until 
    for(i = 0; i < nrbins; i++) {
      if(bins_selected[i] == 0) {   // continuing region
	width++;
      }else {
	if(width <= minwidthbins) {   // Remove region
	  for(j = 0; j < nrbins; j++) {
	    if(bins_selected[j] == 0)
	      bins_selected[j] = 1;
	    else
	      break;
	  }
	  for(j = nrbins-1; j >= 0; j--) {
	    if(bins_selected[j] == 0)
	      bins_selected[j] = 1;
	    else
	      break;
	  }
	}
	break;
      }
    }
  }
}

// Return 0 on error, 1 on success
int expand_notselected(float *profile, int *bins_selected, int nrbins, float mean, float rms, float sigma2, pulselongitude_regions_definition *onpulse, verbose_definition verbose)
{
  int i, reg, quitloop;
  float dy;
  if(selection_to_region(bins_selected, nrbins, onpulse, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR expand_notselected: Cannot construct region (selection_to_region failed at first call).");
    return 0;
  }

  // Expand not-selected regions
  for(reg = 0; reg < onpulse->nrRegions; reg++) {
    i = onpulse->right_bin[reg];
    if(i != nrbins-1) {
      quitloop = 0;
      do {
	if(i < 0) {
	  if(quitloop == 0) {
	    i += nrbins;
	    quitloop = 1;
	  }else {
	    quitloop = 2;
	  }
	}
	dy = fabs(profile[i] - mean);
	//	printf("XXXXX %d %lf, %lf\n", i, dy, dy/(sigma2*rms));
	if(dy > sigma2*rms)
	  bins_selected[i] = 0;
	else
	  quitloop = 2;
	i--;
      }while(quitloop < 2);
    }
    
    i = onpulse->left_bin[reg];
    if(i != 0) {
      quitloop = 0;
      do {
	if(i >= nrbins) {
	  if(quitloop == 0) {
	    i -= nrbins;
	    quitloop = 1;
	  }else {
	    quitloop = 2;
	  }
	}
	dy = fabs(profile[i] - mean);
	//            printf("XXXXX %d %lf, %lf\n", i, dy, dy/(sigma2*rms));
	if(dy > sigma2*rms)
	  bins_selected[i] = 0;
	else
	  quitloop = 2;
	i++;
      }while(quitloop < 2);
    }
  }
  
  if(selection_to_region(bins_selected, nrbins, onpulse, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR expand_notselected: Cannot construct region (selection_to_region failed at call 2).");
    return 0;
  }
  return 1;
}

/* Estimate the onpulse region automatically of the DEDISPERSED
   profile. It actually finds the offpulse region, so the onpulse
   region is probably a bit too large, so it is best used to use the
   rest as baseline. Set smooth_width to the width (in rotational
   phase) of the smoothing function to determine the minimum in the
   profile, or to a negative value to take the default value (might
   change in future). sigma is the nr of times the profile intensity
   has to be deviating from the mean with the rms of the offpulse
   region to count as a detection (set to negative to take the
   default). minwidthbins removes not selected regions of a maximum of
   this width from the not-selected regions (set to negative to take
   the default). sigma2 is similar to sigma, but after applying the
   minwidthbins option, the not-selected regions are expanded using
   the sigma2 limit (which should be lower than sigma) to ensure there
   is no low-level emission included in not-selected regions. The same
   process is repeated after smoothing the profile with
   smooth_width2. Set nritts to the number of itterations to reach
   final result (set to negative to take the default). Return 1 =
   success, 0 on error. Verbose level determines nr of spaces before
   output. If nocounters is not set, the progresss is shown. */
int preprocess_autoOnpulse(datafile_definition original, float smooth_width, float smooth_width2, float sigma, float sigma2, int minwidthbins, int nritts, pulselongitude_regions_definition *onpulse, verbose_definition verbose)
{
  //Comment following line out to avoid debug mode, leave it in to enable debugging
  //  #define Enable_preprocess_autoOnpulse_debugflag 1

  int i, bin_ymin, *bins_selected, itt, widest;
  float ymin, mean, rms, dy;
  datafile_definition profile, profile_smooth;
  verbose_definition verbose2;

#ifdef Enable_preprocess_autoOnpulse_debugflag
  pgplot_viewport_def viewport;
  pgplot_box_def pgplotbox;
#endif

  if(smooth_width < 0)
    smooth_width = 0.1;
  if(smooth_width2 < 0)
    smooth_width2 = 0.05;
  if(sigma < 0)
    sigma = 3;
  if(sigma2 < 0)
    sigma2 = 1;
  if(minwidthbins < 0)
    minwidthbins = 2;
  if(nritts < 0)
    nritts = 3;

  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Estimating off-pulse region using smoothwidth=%f smoothwidth2=%f sigma=%f sigma2=%f minwidth=%d\n", smooth_width, smooth_width2, sigma, sigma2, minwidthbins);
    verbose.indent += 2;
    verbose2.indent += 2;
  }

  clearPulselongitudeRegion(onpulse);

  //  fflush(stdout); fprintf(stderr, "XXXXXX");
  if(preprocess_make_profile(original, &profile, 1, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Making profile failed.");
    return 0;
  }
  //  fflush(stdout); fprintf(stderr, "XXXXXX");


#ifdef Enable_preprocess_autoOnpulse_debugflag
  pgplot_clear_viewport_def(&viewport);
  clear_pgplot_box(&pgplotbox);
  strcpy(pgplotbox.ylabel, "Stokes I");
  strcpy(pgplotbox.title, "Original profile");
  pgplotGraph1(viewport, profile.data, NULL, NULL, profile.NrBins, 0, profile.NrBins-1, 0, 0, profile.NrBins-1, 0, 0, 0, pgplotbox, 0, 0, 0, 1, 1, NULL);
#endif

  if(make_clone(profile, &profile_smooth, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Making clone failed.");
    return 0;
  }

  if(tophat_smooth(profile_smooth.data, profile.NrBins, profile.NrBins*smooth_width, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Smoothing of profile failed.");
    return 0;
  }

#ifdef Enable_preprocess_autoOnpulse_debugflag
  pgplot_clear_viewport_def(&viewport);
  clear_pgplot_box(&pgplotbox);
  strcpy(pgplotbox.ylabel, "Stokes I");
  strcpy(pgplotbox.title, "Smooth version of profile");
  pgplotGraph1(viewport, profile_smooth.data, NULL, NULL, profile.NrBins, 0, profile.NrBins-1, 0, 0, profile.NrBins-1, 0, 0, 0, pgplotbox, 1, 0, 0, 1, 1, NULL);
#endif

  for(i = 0; i < profile_smooth.NrBins; i++) {
    if(i == 0 || profile_smooth.data[i] < ymin) {
      bin_ymin = i;
      ymin = profile_smooth.data[i];
    }
  }
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Found minimum in smoothed profile at bin %d\n", bin_ymin);
  }

  bins_selected = calloc(profile.NrBins, sizeof(int));
  if(bins_selected == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Cannot allocate memory.");
    return 0;
  }

  for(i = 0; i < profile.NrBins; i++) {
    if((double)abs(i - bin_ymin) < 0.5*profile.NrBins*smooth_width) {
      bins_selected[i] = 1;
    }
    if((double)labs(i - bin_ymin + profile.NrBins) < 0.5*profile.NrBins*smooth_width) {
      bins_selected[i] = 1;
    }
    if((double)labs(i - bin_ymin - profile.NrBins) < 0.5*profile.NrBins*smooth_width) {
      bins_selected[i] = 1;
    }
  }

  if(selection_to_region(bins_selected, profile.NrBins, onpulse, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Cannot construct region (selection_to_region failed at call 1).");
    return 0;
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Initial guess for offpulse region:\n");
    printRegions(onpulse, verbose);
  }


#ifdef Enable_preprocess_autoOnpulse_debugflag
  pgplot_clear_viewport_def(&viewport);
  clear_pgplot_box(&pgplotbox);
  strcpy(pgplotbox.ylabel, "Stokes I");
  strcpy(pgplotbox.title, "Initial guess offpulse region");
  pgplotGraph1(viewport, profile.data, NULL, NULL, profile.NrBins, 0, profile.NrBins-1, 0, 0, profile.NrBins-1, 0, 0, 0, pgplotbox, 1, 0, 0, 1, 1, onpulse);
#endif

  for(itt = 0; itt < nritts; itt++) {
    get_mean_rms(profile.data, bins_selected, profile.NrBins, &mean, &rms);
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  itt %d: mean=%f, rms=%f\n", itt, mean, rms);
    }

    // Select all < sigma points
    for(i = 0; i < profile.NrBins; i++) {
      dy = fabs(profile.data[i] - mean);
      if(dy < sigma*rms)
	bins_selected[i] = 1;
      else
	bins_selected[i] = 0;
      //      printf("XXXXX %d: %d %f (%lf, %lf)\n", i, bins_selected[i], profile.data[i], dy, dy/(sigma*rms));
    }

    widest = widest_notselected(bins_selected, profile.NrBins);
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  itt %d: widest not selected region=%d bins using new rms\n", itt, widest);
    }

#ifdef Enable_preprocess_autoOnpulse_debugflag
    pgplot_clear_viewport_def(&viewport);
    clear_pgplot_box(&pgplotbox);
    strcpy(pgplotbox.ylabel, "Stokes I");
    strcpy(pgplotbox.title, "Profile after de-selecting high rms points (step 1)");
    if(selection_to_region(bins_selected, profile.NrBins, onpulse, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Cannot construct region (selection_to_region failed at call 2).");
      return 0;
    }
    pgplotGraph1(viewport, profile.data, NULL, NULL, profile.NrBins, 0, profile.NrBins-1, 0, 0, profile.NrBins-1, 0, 0, 0, pgplotbox, 1, 0, 0, 1, 1, onpulse);
    printf("  itt %d: Step 1:\n", itt);
    printRegions(onpulse, verbose);
#endif



    // Expand non-selected regions using the non-smoothed profile    
    if(expand_notselected(profile.data, bins_selected, profile.NrBins, mean, rms, sigma2, onpulse, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Expanding regions failed.");
      return 0;
    }
    
    // And again using a smoothed version
    memcpy(profile_smooth.data, profile.data, profile.NrBins);
    if(tophat_smooth(profile_smooth.data, profile.NrBins, profile.NrBins*smooth_width2, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Smoothing of profile failed.");
      return 0;
    }
    get_mean_rms(profile_smooth.data, bins_selected, profile.NrBins, &mean, &rms);
    if(expand_notselected(profile_smooth.data, bins_selected, profile.NrBins, mean, rms, sigma2, onpulse, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Expanding regions failed.");
      return 0;
    }

#ifdef Enable_preprocess_autoOnpulse_debugflag
    pgplot_clear_viewport_def(&viewport);
    clear_pgplot_box(&pgplotbox);
    strcpy(pgplotbox.ylabel, "Stokes I");
    strcpy(pgplotbox.title, "Profile after trying to expand on-pulse region (step 2)");
    if(selection_to_region(bins_selected, profile.NrBins, onpulse, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Cannot construct region (selection_to_region failed at call 3).");
      return 0;
    }
    pgplotGraph1(viewport, profile.data, NULL, NULL, profile.NrBins, 0, profile.NrBins-1, 0, 0, profile.NrBins-1, 0, 0, 0, pgplotbox, 1, 0, 0, 1, 1, onpulse);
    printf("  itt %d: step 2\n", itt);
    printRegions(onpulse, verbose);
#endif




    remove_narrow_notselected(bins_selected, profile.NrBins, widest, minwidthbins);

    get_mean_rms(profile.data, bins_selected, profile.NrBins, &mean, &rms);
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  itt %d: mean=%f, rms=%f after removing narrow spikes\n", itt, mean, rms);
    }

#ifdef Enable_preprocess_autoOnpulse_debugflag
    pgplot_clear_viewport_def(&viewport);
    clear_pgplot_box(&pgplotbox);
    strcpy(pgplotbox.ylabel, "Stokes I");
    strcpy(pgplotbox.title, "Profile after removed very narrow onpulse regions (step 3)");
    if(selection_to_region(bins_selected, profile.NrBins, onpulse, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Cannot construct region (selection_to_region failed at call 4).");
      return 0;
    }
    pgplotGraph1(viewport, profile.data, NULL, NULL, profile.NrBins, 0, profile.NrBins-1, 0, 0, profile.NrBins-1, 0, 0, 0, pgplotbox, 1, 0, 0, 1, 1, onpulse);
    printf("  itt %d: Step 3\n", itt);
    printRegions(onpulse, verbose);
#endif

    //    for(i = 0; i < profile.NrBins; i++) {
    //      printf("XXXXX %d %d\n", i, bins_selected[i]);
    //    }



  }

  // Invert selection
  for(i = 0; i < profile.NrBins; i++) {
    if(bins_selected[i])
      bins_selected[i] = 0;
    else
      bins_selected[i] = 1;
  }

  if(selection_to_region(bins_selected, profile.NrBins, onpulse, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_autoOnpulse: Cannot construct region (selection_to_region failed at call 5).");
    return 0;
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Final estimate for offpulse region:\n");
    verbose_definition verbose2;
    copyVerboseState(verbose, &verbose2);
    verbose2.indent = verbose.indent + 2;
    printRegions(onpulse, verbose2);
  }

#ifdef Enable_preprocess_autoOnpulse_debugflag
  pgplot_clear_viewport_def(&viewport);
  clear_pgplot_box(&pgplotbox);
  strcpy(pgplotbox.ylabel, "Stokes I");
  strcpy(pgplotbox.title, "Profile");
  pgplotGraph1(viewport, profile.data, NULL, NULL, profile.NrBins, 0, profile.NrBins-1, 0, 0, profile.NrBins-1, 0, 0, 0, pgplotbox, 1, 0, 0, 1, 1, onpulse);
#endif

  closePSRData(&profile, 0, verbose);
  closePSRData(&profile_smooth, 0, verbose);
  free(bins_selected);
  return 1;
}


//START REGION RELEASE

/*
  Replace the reference frequency with the new frequency new_freq (set
  to -1 or 1e10 if infinite frequency). If the data is not
  dedispersed/de-Faraday rotated, nothing is done except changing the
  header parameter. If the data is dedispersed/de-Faraday rotated, it
  will be adjusted to make it consistent with the new reference
  frequency.

  Return
    0 = error.
    1 = success.
*/
int preprocess_changeRefFreq(datafile_definition *original, double freq_ref_new, verbose_definition verbose)
{
  double freq_ref_cur;
  int i, inffreq_cur, inffreq_new;
  verbose_definition verbose2;

  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 4;

  freq_ref_cur = original->freq_ref;
  inffreq_cur = 0;
  if((original->freq_ref > -1.1 && original->freq_ref < -0.9) || (original->freq_ref > 0.99e10 && original->freq_ref < 1.01e10)) {
    inffreq_cur = 1;
  }else if(original->freq_ref < 0) {  // Unknown reference frequency.
    if(original->isDeDisp || original->isDeFarad) {
      printwarning(verbose.debug, "WARNING preprocess_changeRefFreq: Current reference frequency is unknown, it is assumed it was infinite frequency");
    }
    inffreq_cur = 1;
    freq_ref_cur = 1e10; // This is the code used to indicate infinite frequency.
  }

  if((freq_ref_new > -1.1 && freq_ref_new < -0.9) || (freq_ref_new > 0.99e10 && freq_ref_new < 1.01e10)) {
    inffreq_new = 1;
  }else {
    inffreq_new = 0;
  }


  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Updating reference frequency of %s from ", original->filename);
    if(inffreq_cur)
      printf("infinity");
    else
      printf("%f MHz", freq_ref_cur);
    if(inffreq_new)
      printf(" to infinity\n");
    else
      printf(" to %f MHz\n", freq_ref_new);
  }

  original->freq_ref = freq_ref_new;

  if(original->isDeDisp == 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Data not yet dispersed, so no dedispersion done.\n");
    }
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Re-dedispersing data.\n");
    }
    if(preprocess_dedisperse(original, 0, 1, freq_ref_cur, verbose2) == 0) {
      printerror(verbose.debug, "WARNING preprocess_changeRefFreq: Re-dedispersion failed");
      return 0;
    }
  }

  if(original->isDeFarad == 0) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Data not yet de-Faraday rotated, so no de-Faraday rotation done.\n");
    }
  }else {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  Re-de-Farady rotating data.\n");
    }
    //    printf("XXXXX %lf\n", freq_ref_cur);
    if(preprocess_deFaraday(original, 0, 1, freq_ref_cur, NULL, verbose2) == 0) {
      printerror(verbose.debug, "WARNING preprocess_changeRefFreq: Re-de-Faraday rotating data failed");
      return 0;
    }
  }

  return 1;
}

//START REGION DEVELOP

/* The 0-DM profile is subtracted from the individual non-zapped
   frequency channels. Any dedispersion is removed from the data
   first. This default behaviour is obtained if both
   posdips_stddev_thresh and negdips_stddev_thresh are set to a
   negative value.

   If float negdips_stddev_thresh is set to a positive number, then
   only pulse phases where there is a dip in the intensity exceeding
   negdips_stddev_thresh times the standard deviation will be
   affected.

   posdips_stddev_thresh does something similar to positive spikes.

   if either or both posdips_stddev_thresh and negdips_stddev_thresh
   are positive (or zero), then mode refines the behaviour:

   mode=1: Subtract freq-avrg from each channel
   mode=2: Set each channel to its avrg (zero per def)

   if either or both posdips_stddev_thresh and negdips_stddev_thresh
   are positive (or zero), then maxw refines the behaviour:

   When searching for peaks, if maxw is a positive integer, it will do
   a boxcar search for the most significant peak of various widths in
   the 0-DM profile. Only if the most significant peak is more narrow
   than maxw bins, it will be removed. At either side of the found
   signal extrabins extra bins are removed. In addition, if
   checkdedispersed is non-zero, the width of a signal will be
   compared with the width after dedispersion. Actual pulsar signals
   should become narrower, undispersed signals wider. So only if the
   width after dedispersion becomes wider by minwidening bins (allowed
   to be negative), the signal is regarded to be RFI signal.

   Return 1 = success, 0 on error. 
*/
int preprocess_zero_dming(datafile_definition *original, float posdips_stddev_thresh, float negdips_stddev_thresh, int mode, int maxw, int extrabins, int checkdedispersed, int minwidening, verbose_definition verbose)
{
  int i, process_allbins, posOrNeg;
  long p, f, n, b;
  datafile_definition subints_averages, subints_averages_dedispersed;
  verbose_definition verbose_indent;
  verbose_definition verbose_noverbose; // Disable verbose, unless debug is set

  if(negdips_stddev_thresh < 0.0 && posdips_stddev_thresh < 0.0) {
    process_allbins = 1;  // Means that all bins are processed
    posOrNeg = 0;  // In fact not used in this mode
  }else {
    process_allbins = 0;  // Means that bins are processed dependent of the detection of spikes
    if(mode != 1 && mode != 2) {
      printerror(verbose.debug, "ERROR zero_dming: mode is expected to be 1 or 2.");
      return 0;
    }
    if(posdips_stddev_thresh >= 0.0 && negdips_stddev_thresh < 0.0) {
      posOrNeg = 0;
    }else if(posdips_stddev_thresh >= 0.0 && negdips_stddev_thresh >= 0.0) {
      posOrNeg = 1;
    }else if(posdips_stddev_thresh < 0.0 && negdips_stddev_thresh >= 0.0) {
      posOrNeg = 2;
    }
  }
  if(process_allbins != 0 && maxw > 0) {
    printerror(verbose.debug, "ERROR zero_dming: Cannot set a maximum width if not defining a significance threshold.");
    return 0;
  }
  if(maxw > 0 && posOrNeg == 1 && posdips_stddev_thresh != negdips_stddev_thresh) {
    if(posdips_stddev_thresh > negdips_stddev_thresh) {
      negdips_stddev_thresh = posdips_stddev_thresh;
    }else {
      posdips_stddev_thresh = negdips_stddev_thresh;
    }
    printwarning(verbose.debug, "WARNING zero_dming: Cannot have separate thresholds for positive/negative dips when maxw pararameter is positive.");
  }
  if(maxw <= 0 && checkdedispersed) {
    printerror(verbose.debug, "ERROR zero_dming: Cannot set the dedispersed width if no maximum width maxw is specified.");
    return 0;
  }
  if(checkdedispersed) {
    printwarning(verbose.debug, "WARNING zero_dming: Cleaning might not work properly at pulse edges because of wrap issues?");
  }


  copyVerboseState(verbose, &verbose_indent);
  copyVerboseState(verbose, &verbose_noverbose);
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Zero-DMing data");
    if(negdips_stddev_thresh >= 0.0) {
      printf(" at negative dips exceeding %f sigma", negdips_stddev_thresh);
      if(posdips_stddev_thresh >= 0.0) {
	printf(" and");
      }
    }
    if(posdips_stddev_thresh >= 0.0) {
      if(negdips_stddev_thresh >= 0.0) {
	printf(" and");
      }
      printf(" for positive spikes exceeding %f sigma", posdips_stddev_thresh);
    }
    if(maxw > 0) {
      printf(" with a maximum width of %d bins (undedispersed)", maxw);
    }
    if(checkdedispersed) {
      printf(" while ensuring that the width increases by >= %d bins after dedispersion", minwidening);
    }
    if(extrabins) {
      printf(" and removing %d extra bins at either side", extrabins);
    }
    printf("\n");
    // Increase indentation level
    verbose_indent.indent += 4;
    verbose_noverbose.indent += 4;
    if(verbose.debug == 0) {  // Disable verbose, unless debug is set
      verbose_noverbose.verbose = 0;
    }
  }

  if(original->format != MEMORY_format) {
    printerror(verbose.debug, "ERROR zero_dming: only works if data is loaded into memory.");
    return 0;
  }

  if(original->NrFreqChan <= 1) {
    printerror(verbose.debug, "ERROR zero_dming: Expected data to have frequency resolution.");
    return 0;
  }
  if(original->NrPols > 1) {
    if(original->poltype != POLTYPE_STOKES) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING zero_dming: Need to think about how to deal with non-Stokes data.");
      if(negdips_stddev_thresh >= 0.0 || posdips_stddev_thresh >= 0.0) {
	printwarning(verbose.debug, "WARNING zero_dming: First polarization channel is used to determine which pulse phases to correct.");
      }
    }
  }


  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  - Make sure the data is undispersed\n");
  }
  if(preprocess_dedisperse(original, 1, 0, 0, verbose_indent) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR zero_dming: remove dedispersion failed, cannot work on undedispersed data.");
    return 0;
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  - Make sure the average signal in each channel is zero.\n");
    printwarning(verbose.debug, "WARNING zero_dming: Baseline subtraction should be done EXCLUDING spikes and pulsar signal ideally.....");
  }
  float *baseline;
  if(preprocess_debase(original, NULL, &baseline, 0, verbose_indent) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR zero_dming: removing baseline failed, cannot work on undedispersed data.");
    return 0;
  }

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  - Make not dedispersed frequency-scrunched template(s)\n");
  }
  if(preprocess_addsuccessiveFreqChans(*original, &subints_averages, original->NrFreqChan, NULL, verbose_indent) == 0) {
    printerror(verbose.debug, "ERROR zero_dming: Doing a frequency scrunch failed.");
    return 0;
  }

  if(checkdedispersed) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("  - Make dedispersed frequency-scrunched template(s)\n");
    }
    datafile_definition clone;
    if(make_clone(*original, &clone, verbose_indent) == 0) {
      printerror(verbose.debug, "ERROR zero_dming: Cannot make clone of data, cannot make a dedispersed equivalent dataset.");
      return 0;
    }
    // Not dedispersed, so this should be change in header parameter only
    if(preprocess_changeRefFreq(&clone, get_centre_frequency(clone, verbose_indent), verbose_indent) == 0) {
      printerror(verbose.debug, "ERROR zero_dming: Cannot change reference frequency for the dedispersed equivalent dataset.");
      return 0;
    }
    // Set reference freq to mid-point in the band, such that the peak location should remain the same
    if(preprocess_dedisperse(&clone, 0, 0, 0, verbose_indent) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR zero_dming: remove dedispersion failed, cannot work on undedispersed data.");
      return 0;
    }
    if(preprocess_addsuccessiveFreqChans(clone, &subints_averages_dedispersed, original->NrFreqChan, NULL, verbose_indent) == 0) {
      printerror(verbose.debug, "ERROR zero_dming: Doing a frequency scrunch failed.");
      return 0;
    }
    closePSRData(&clone, 0, verbose);
  }

  /*
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  - Scaling frequency-scrunched profile(s)\n");
  }
  if(preprocess_scale(subints_averages, 1.0/(original->NrFreqChan), 0.0, verbose_indent) != 1) {
    printerror(verbose.debug, "ERROR zero_dming: Scaling data failed.");
    return 0;
  }
  */

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  - Subtracting template(s) from data\n");
  }
  float templatesample, origsample;
  long nrzapped;
  int *iszapped;
  int *processbin;
  float rms, avrg;

  iszapped = malloc(original->NrFreqChan*sizeof(int));
  processbin = malloc(original->NrBins*sizeof(int));
  if(iszapped == NULL || processbin == NULL) {
    printerror(verbose.debug, "ERROR zero_dming: Memory allocation error.");
    return 0;
  }


  for(n = 0; n < original->NrSubints; n++) {
    int pulsarsignal = 0;
    if(nrchannels_iszapped(*original, n, &nrzapped, verbose_indent) == 0) {
      printerror(verbose.debug, "ERROR zero_dming: Determining the number of zapped channels failed.");
      free(processbin);
      free(iszapped);
      return 0;
    }
    // Determine average and rms of first polarization channel of the template
    if(process_allbins == 0) {
      if(rmsPSRData(subints_averages, &rms, &avrg, NULL, 0, n, 0, 0, verbose_indent) == 0) {
	printerror(verbose.debug, "ERROR zero_dming: Determining rms of channel failed.");
	free(processbin);
	free(iszapped);
	return 0;
      }
    }
    // Make an lookup table if frequency channel is zapped
    for(f = 0; f < original->NrFreqChan; f++) {
      if(pulse_iszapped(*original, n, f, &(iszapped[f]), verbose_indent) == 0) {
	printerror(verbose.debug, "ERROR zero_dming: Determining the zap state of channel failed.");
	free(processbin);
	free(iszapped);
	return 0;
      }
    }
 
    if(process_allbins) {  // If not looking for excess sample values, then correct bin regardless of anything
      for(b = 0; b < original->NrBins; b++) {
	processbin[b] = 1;
      }
    }else if(process_allbins == 0 && maxw <= 0) {  // If looking for excess sample values, and no maxw is defined, then base the decision to zap the bin on the sample value
      for(b = 0; b < original->NrBins; b++) {
	processbin[b] = 0;
      }
      for(b = 0; b < original->NrBins; b++) {
	if(readPulsePSRData(&subints_averages, n, 0, 0, b, 1, &templatesample, verbose_indent) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR zero_dming: Cannot read data.");
	  free(processbin);	
	  free(iszapped);
	  return 0;
	}
	if(negdips_stddev_thresh >= 0.0) {
	  if((avrg-templatesample)/rms >= negdips_stddev_thresh) { // Significant negative deflection
	    long b2;
	    for(b2 = b - extrabins; b2 <= b + extrabins; b2++) {
	      if(b2 >= 0 && b2 < original->NrBins) {
		processbin[b2] = 1;
	      }
	    }
	  }
	}
	if(posdips_stddev_thresh >= 0.0) {
	  if((templatesample-avrg)/rms >= posdips_stddev_thresh) { // Significant positive deflection
	    long b2;
	    for(b2 = b - extrabins; b2 <= b + extrabins; b2++) {
	      if(b2 >= 0 && b2 < original->NrBins) {
		processbin[b2] = 1;
	      }
	    }
	  }
	}
      }
   }else if(process_allbins == 0 && maxw > 0) {  // If looking for excess sample values, and maxw is defined, then base the decision to zap he bin on a boxcar search and the width of the most significant peak found

      int peak_bin_start, peakwidth, peak_bin_start_dispersed, peakwidth_dispersed;
      float peak_snr, peak_snr_dispersed;

      // Obtains a pointer current pulse data in 0DM template
      float *pulse_ptr;
      if(get_pointer_PulsePSRData(&subints_averages, n, 0, 0, 0, &pulse_ptr, verbose_indent) == 0) {
	printerror(verbose.debug, "ERROR zero_dming: Data access error.");
	free(processbin);
	free(iszapped);
	return 0;
      }
      if(boxcarFindpeak(pulse_ptr, original->NrBins, NULL, &peak_bin_start, &peakwidth, &peak_snr, NULL, 0, posOrNeg, 0, 1, -1, 0, 0, verbose_noverbose) == 0) {
	printerror(verbose.debug, "ERROR zero_dming: Boxcar peak search failed.");
	free(processbin);
	free(iszapped);
	return 0;
      }
      if(verbose.debug) {
	for(i = 0; i < verbose.indent; i++)      
	  printf(" ");
	printf("  - boxcar search found something at bin=%d with width %d bins and S/N=%lf.\n", peak_bin_start, peakwidth, peak_snr);
	printwarning(verbose.debug, "WARNING zero_dming: Baseline subtraction should be done EXCLUDING spikes and pulsar signal ideally.....");	
      }

      if(checkdedispersed) {  // Check if the found pulse could be a real pulsar signal
	  if(get_pointer_PulsePSRData(&subints_averages_dedispersed, n, 0, 0, 0, &pulse_ptr, verbose_indent) == 0) {
	  printerror(verbose.debug, "ERROR zero_dming: Data access error.");
	  free(processbin);
	  free(iszapped);
	  return 0;
	}
	if(boxcarFindpeak(pulse_ptr, original->NrBins, NULL, &peak_bin_start_dispersed, &peakwidth_dispersed, &peak_snr_dispersed, NULL, 0, posOrNeg, 0, 1, -1, 0, 0, verbose_noverbose) == 0) {
	  printerror(verbose.debug, "ERROR zero_dming: Boxcar peak search failed.");
	  free(processbin);
	  free(iszapped);
	  return 0;
	}
	if(verbose.debug) {
	  for(i = 0; i < verbose.indent; i++)      
	    printf(" ");
	  printf("  - boxcar search found something at bin=%d with width %d bins and S/N=%lf in dedispersed data.\n", peak_bin_start_dispersed, peakwidth_dispersed, peak_snr_dispersed);
	  printwarning(verbose.debug, "WARNING zero_dming: Baseline subtraction should be done EXCLUDING spikes and pulsar signal ideally.....");	
	}
	if(peakwidth_dispersed - peakwidth < minwidening) {
	  pulsarsignal = 1;
	  if(verbose.debug) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("    Signal did narrow enough after dedispersion to be considered pulsar signal.\n");
	  }
	}else {
	  pulsarsignal = 0;
	  if(verbose.debug) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	    printf("    Signal didn't narrow enough after dedispersion to be considered pulsar signal.\n");
	  }
	}
	// Compare if centres got moved by too much
	if(fabs(peak_bin_start_dispersed+0.5*peakwidth_dispersed-(peak_bin_start+0.5*peakwidth)) > 0.5*peakwidth) {
	  pulsarsignal = 1;
	  //	    if(verbose.debug) {
	  //	      for(i = 0; i < verbose.indent; i++)      
	  //		printf(" ");
	  printwarning(verbose.debug, "Signal did move too much in time between dedispersed and non-dedispersed datasets, suggestive of different peaks being detected in profile: no correction is done for subint %ld (counting from zero).\n", n);
	  //	    }
	}
      }

      for(b = 0; b < original->NrBins; b++) {
	processbin[b] = 0;
	if(pulsarsignal == 0) { // Not widened enough after dedispersion, not considered real RFI
	  if(peakwidth <= maxw) {
	    if(b >= peak_bin_start-extrabins) {
	      if(b < peak_bin_start+peakwidth+extrabins) {
		if(posOrNeg == 2) {  // Looking for negative dips
		  if(peak_snr >= negdips_stddev_thresh) {
		    processbin[b] = 1;
		  }
		}else {
		  if(peak_snr >= posdips_stddev_thresh) {  // Note posdips_stddev_thresh=negdips_stddev_thresh if looking for both +&-
		    processbin[b] = 1;
		  }
		}
	      }
	    }
	  }
	}
      }
    }  // End of if(process_allbins == 0 && maxw > 0)
    if(verbose.debug && checkdedispersed) {
      printf("Dispersed and dedispersed templates are:");
      for(b = 0; b < original->NrBins; b++) {
	float *ptr1, *ptr2;
	if(get_pointer_PulsePSRData(&subints_averages, n, 0, 0, 0, &ptr1, verbose_indent) == 0) {
	  printerror(verbose.debug, "ERROR zero_dming: Data access error.");
	  free(processbin);
	  free(iszapped);
	  return 0;
	}
	if(get_pointer_PulsePSRData(&subints_averages_dedispersed, n, 0, 0, 0, &ptr2, verbose_indent) == 0) {
	  printerror(verbose.debug, "ERROR zero_dming: Data access error.");
	  free(processbin);
	  free(iszapped);
	  return 0;
	}
	printf("%ld %e %e\n", b, ptr1[b], ptr2[b]);
      }
    }

    for(b = 0; b < original->NrBins; b++) {
      if(processbin[b]) {
	for(p = 0; p < original->NrPols; p++) {
	  if(p != 0) {   // If p == 0, then templatesample is already determined
	    if(readPulsePSRData(&subints_averages, n, p, 0, b, 1, &templatesample, verbose_indent) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR zero_dming: Cannot read data.");
	      free(processbin);
	      free(iszapped);
	      return 0;
	    }
	  }
	  for(f = 0; f < original->NrFreqChan; f++) {
	    //	    iszapped[f] = 0;  // Useful to see which bins are corrected (since zapped channels are now affected as well)
	    if(iszapped[f] == 0) {
	      if(readPulsePSRData(original, n, p, f, b, 1, &origsample, verbose_indent) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR zero_dming: Cannot read data.");
		free(processbin);
		free(iszapped);
		return 0;
	      }
	      if(process_allbins == 0 && mode == 2) {
		origsample = 0; // Set signal to zero
	      }else {
		origsample -= templatesample/(float)(original->NrFreqChan-nrzapped); // Subtract out average signal
	      }
	      //	      origsample = 0;
	      if(writePulsePSRData(original, n, p, f, b, 1, &origsample, verbose_indent) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR zero_dming: Cannot write data.");
		free(processbin);
		free(iszapped);
		return 0;
	      }
	    }
	  }
	  if(verbose.verbose && verbose.nocounters == 0) {
	    for(i = 0; i < verbose_indent.indent; i++)      
	      printf(" ");
	    printf("  %.1f%%     \r", (100.0*((n*original->NrBins+b)*original->NrPols+p))/(float)(original->NrSubints*original->NrBins*original->NrPols));
	    fflush(stdout);
	  }
	}
      }
    }
  }
  
  closePSRData(&subints_averages, 0, verbose);
  if(checkdedispersed) {
    closePSRData(&subints_averages_dedispersed, 0, verbose);
  }


  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  - Restoring baseline in data\n");
  }
  if(preprocess_restore_debase(original, baseline, verbose_indent) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR zero_dming: restoring baseline failed.");
    return 0;
  }
  free(baseline);

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("  - done             \n");
  }
  free(processbin);
  free(iszapped);
  return 1;
}

//START REGION DEVELOP
