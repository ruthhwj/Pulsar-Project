//START REGION RELEASE
#include <stdio.h>
#include <math.h>
#include "psrsalsa.h"

/* Input Ipulse (NrBins long), output Ipulse2 (NrBins2 long). 
   If noDependencyWarning is set, there is no warning generated 
   when the rebinning is not done with an integer amount (this
   is used when warning is already generated in wrapper function,
   and no warning should be generated for each subint for instance).
   The resulting intensities are weighted such that the mean remains the same.
   Returns 1 if succesfull. */
int rebinPulse(float *Ipulse, long NrBins, float *Ipulse2, long NrBins2, int noDependencyWarning, verbose_definition verbose)
{
  long j, i1, i2;
  float x, x2;
  if(noDependencyWarning == 0) {
    if(NrBins % NrBins2 != 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING rebinPulse: Rebinning from %ld to %ld bins implies that separate bins are not entirely independent.", NrBins, NrBins2);
    }
  }
  /* Clear pulse */
  for(j = 0; j < NrBins2; j++) 
    Ipulse2[j] = 0;
  for(j = 0; j < NrBins; j++) {
    x = (j)/(float)NrBins;
    x *= NrBins2;
    x2 = (j+1)/(float)NrBins;
    x2 *= NrBins2;
    i1 = x;
    i2 = x2;
//START REGION DEVELOP
    /*   fprintf(stderr, "%f %f %ld %ld (nrbins=%ld)\n", x, x2, i1, i2, NrBins2); */
    /* For each original bin j, work out the bin range i1,i2 (x,
       x2) in which it falls in the rebinned data. Because the
       nr of bins becomes smaller, the bin falls either entirely
       in one rebinned bin, or it falls in two output
       bins. However, due to rounding errors it can fall in
       three bins as well if the in/output nr of bins are
       (exactly) the same. */
//START REGION RELEASE
    if(i1 == i2) {
      Ipulse2[i1] += Ipulse[j]*(x2-x);
    }else if(i2-i1 == 1) {
      Ipulse2[i1] += Ipulse[j]*(i2-x);
      if(i2 < NrBins2)
	Ipulse2[i2] += Ipulse[j]*(x2-i2);
    }else if(i2-i1 == 2) {
      Ipulse2[i1]   += Ipulse[j]*(i1+1-x);
      Ipulse2[i1+1] += Ipulse[j];
      if(i2 < NrBins2)
	Ipulse2[i2]   += Ipulse[j]*(x2-i2);
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR rebinPulse: Error in rebinning function.");
      return 0;
    }
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE


/* Takes data from fin and write it to fout (which will be opened and
   the header should be written out, as well as a history if argc !=
   0). The data is shifted with shift bins and the shift is carried
   over to the next pulse. If circularShift is set, no pulses are lost
   because the last pulse is used to fill in the first pulse. oformat
   can be set to MEMORY_format in which case the data will be written
   to memory. fout should only be declared, but it should be closed
   afterwards.  If verbose2 is set it lets you know what parts of the
   pulse are written out. Return 1 = success */
int continuous_shift(datafile_definition fin, datafile_definition *fout, int shift, int circularShift, char *output_name, int oformat, int argc, char **argv, verbose_definition verbose, int verbose2)
{
  int i;
  long nout, p, f, n;
  float *Ipulse, *Ifirst, *Ilast;
  verbose_definition verbose_counters_verbose2;

  if(fin.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR continuous_shift: Frequency channels need to be uniformely separated.");
    return 0;
  }

  copyVerboseState(verbose, &verbose_counters_verbose2);
  if(verbose2)
    verbose_counters_verbose2.nocounters = 0;
  else
    verbose_counters_verbose2.nocounters = 1;

  /* Copy header parameters to output header */
  cleanPSRData(fout, verbose);
  copy_params_PSRData(fin, fout, verbose);

  fout->NrSubints = fin.NrSubints-1;         /* Leave out one pulse to shift profile if no circular shift */
  if(circularShift != 0)
    (fout->NrSubints)++;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Continuous shift: shift=%d bins circularShift=%d.\n", shift, circularShift);
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Output data will contain %ld pulses.\n", fout->NrSubints);
  }
  if(fout->NrFreqChan > 1 && fout->isDeDisp == 0) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING continuous_shift: You might want to dedisperse the data first.");
  }
  if(fout->NrSubints == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR continuous_shift: If there is only one pulse circular shifting must be used.");
    return 0;
  }

  /* Open output file and write data */
  if(oformat == MEMORY_format) {
    fout->format = oformat;
    fout->data = (float *)malloc(fout->NrPols*fout->NrBins*fout->NrSubints*fout->NrFreqChan*sizeof(float));
    if(fout->data == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR continuous_shift: Cannot allocate memory.");
      return 0;
    }
  }else {
    if(!openPSRData(fout, output_name, oformat, 1, 0, 0, verbose_counters_verbose2))
      return 0;
    int cmdOnly = 0;
    if(!writeHeaderPSRData(fout, argc, argv, cmdOnly, verbose))
      return 0;
    //    writeHistory(fout, argc, argv, verbose);
  }

  /* Allocate memory */
  Ipulse = (float *)malloc(fout->NrPols*fout->NrBins*sizeof(float));
  Ifirst = (float *)malloc(fout->NrPols*fout->NrBins*sizeof(float));
  Ilast = (float *)malloc(fout->NrPols*fout->NrBins*sizeof(float));
  if(Ipulse == NULL || Ifirst == NULL || Ilast == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR continuous_shift: Cannot allocate memory.");
    return 0;
  }

  if(shift < 0)
    shift += fout->NrBins;

  if(shift >= fout->NrBins || shift < 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR continuous_shift: Shift is not valid.");
    return 0;
  }
    

//START REGION DEVELOP
  /* I think shift is per defenition positive in the code above, so
     part of the code is actually never used!!! */
//START REGION RELEASE

  for(p = 0; p < fin.NrPols; p++) {
    /* Write out shifted stack */
    for(n = 0; n < fin.NrSubints; n++) {
      /* nout is the pulsenr in output of the first bin of pulse n of the input.*/
      if(shift >= 0)
	nout = n;
      else
	nout = n-1;
      for(f = 0; f < fin.NrFreqChan; f++) {
	if(readPulsePSRData(&fin, n, p, f, 0, fin.NrBins, Ipulse, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR continuous_shift: Read error.");
	  return 0;
	}
	if(n == 0) {                                       /* Write last part of first pulse */
	  if(shift == 0 && circularShift != 0) {            /* If no shift, skip first pulse if not circular */
	    if(verbose2) printf("\ncontinuous_shift: Write out whole pulse %ld\n", n);
	    if(shift >= 0) {
	      if(writePulsePSRData(fout, nout, p, f, shift, fin.NrBins-shift, Ipulse, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		return 0;
	      }
	      if(writePulsePSRData(fout, nout+1, p, f, 0, shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		return 0;
	      }
	    }else {
	      if(writePulsePSRData(fout, nout, p, f, fin.NrBins-shift, shift, Ipulse, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		return 0;
	      }
	      if(writePulsePSRData(fout, nout+1, p, f, 0, fin.NrBins-shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		return 0;
	      }
	    }
	  }
	  if(shift > 0) {
	    if(circularShift != 0) {             /* If circular shift, write out last part of last pulse */
	      if(verbose2) printf("\ncontinuous_shift: Write %d bins of last pulse\n", shift);
	      /*		pumawrite(Ilast+fout->NrBins-shift+p*fout->NrBins, sizeof(float), shift, fout); */
	      if(readPulsePSRData(&fin, fin.NrSubints-1, p, f, 0, fin.NrBins, Ilast, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR continuous_shift: Read error.");
		return 0;
	      }
	      if(writePulsePSRData(fout, nout, p, f, 0, shift, Ilast+fout->NrBins-shift, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		return 0;
	      }
	      /* FIX */ /* Only write out whole pulse if there are more than one pulses */
	      if(n != fin.NrSubints-1) {
		if(verbose2) printf("\ncontinuous_shift: Write out whole pulse %ld\n", n);
		if(shift >= 0) {
		  if(writePulsePSRData(fout, nout, p, f, shift, fin.NrBins-shift, Ipulse, verbose) != 1) {
		    fflush(stdout);
		    printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		    return 0;
		  }
		  if(writePulsePSRData(fout, nout+1, p, f, 0, shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
		    fflush(stdout);
		    printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		    return 0;
		  }
		}else {
		  if(writePulsePSRData(fout, nout, p, f, fin.NrBins-shift, shift, Ipulse, verbose) != 1) {
		    fflush(stdout);
		    printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		    return 0;
		  }
		  if(writePulsePSRData(fout, nout+1, p, f, 0, fin.NrBins-shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
		    fflush(stdout);
		    printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		    return 0;
		  }
		}
	      }
	    }else {
	      /* no circular, shift > 0, first pulse */
	      if(verbose2) printf("\ncontinuous_shift: Write %d bins of pulse %ld (freq = %ld)\n", shift, n, f);
	      /*		pumawrite(Ipulse+fout->NrBins-shift, sizeof(float), shift, fout); */
	      if(writePulsePSRData(fout, nout, p, f, 0, shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		return 0;
	      }
	    }
	  }else if(shift < 0) {
	    if(verbose2) printf("\ncontinuous_shift: Write %ld bins of pulse %ld (freq = %ld)\n", fout->NrBins+shift, n, f);
	    /*	      pumawrite(Ipulse-shift, sizeof(float), fout->NrBins+shift, fout); */
	    if(writePulsePSRData(fout, nout, p, f, 0, fout->NrBins+shift, Ipulse-shift, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR continuous_shift: Write error.");
	      return 0;
	    }
	  }
	}
	if(n == fin.NrSubints-1) {           /* Only write out first bit of one last pulse */
	  if(shift >= 0) {
	    if(circularShift == 0 && f == 0)
	      nout -= 1;
	    if(verbose2) printf("\ncontinuous_shift: Write %ld bins of pulse %ld\n", fout->NrBins-shift, n);
	    /*	      pumawrite(Ipulse, sizeof(float), fout->NrBins-shift, fout); */
	    if(writePulsePSRData(fout, nout, p, f, shift, fout->NrBins-shift, Ipulse, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR continuous_shift: Write error.");
	      return 0;
	    }
	  }else if(shift < 0) {
	    if(circularShift != 0) {             /* If circular shift, write out first part of first pulse */
	      /* FIX: Only write out whole pulse if more than one pulse*/
	      if(n != 0) {
		if(verbose2) printf("\ncontinuous_shift: continuous_shift: Write out whole pulse %ld (freq = %ld)\n", n, f);
		/*		  pumawrite(Ipulse, sizeof(float), fout->NrBins, fout); */
		if(shift >= 0) {
		  if(writePulsePSRData(fout, nout, p, f, shift, fin.NrBins-shift, Ipulse, verbose) != 1) {
		    fflush(stdout);
		    printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		    return 0;
		  }
		  if(writePulsePSRData(fout, nout+1, p, f, 0, shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
		    fflush(stdout);
		    printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		    return 0;
		  }
		}else {
		  if(writePulsePSRData(fout, nout, p, f, fin.NrBins-shift, shift, Ipulse, verbose) != 1) {
		    fflush(stdout);
		    printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		    return 0;
		  }
		  if(writePulsePSRData(fout, nout+1, p, f, 0, fin.NrBins-shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
		    fflush(stdout);
		    printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		    return 0;
		  }
		}
	      }
	      if(verbose2) printf("\ncontinuous_shift: Write %d bins of first pulse\n", -shift);
	      if(readPulsePSRData(&fin, 0, p, f, 0, fin.NrBins, Ifirst, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR continuous_shift: Read error.");
		return 0;
	      }
	      /*		pumawrite(Ifirst+p*fout->NrBins, sizeof(float), -shift, fout); */
	      if(writePulsePSRData(fout, nout+1, p, f, fin.NrBins-shift, -shift, Ipulse, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		return 0;
	      }
	    }else {	      
	      /* last pulse, no circular, shift < 0 */
	      if(verbose2) printf("\ncontinuous_shift: Write %d bins of pulse %ld\n", -shift, n);
	      /*		pumawrite(Ipulse, sizeof(float), -shift, fout); */
	      if(writePulsePSRData(fout, nout, p, f, fin.NrBins-shift, -shift, Ipulse, verbose) != 1) {
		fflush(stdout);
		printerror(verbose.debug, "ERROR continuous_shift: Write error.");
		return 0;
	      }
	    }
	  }
	}
	if(n != fin.NrSubints-1 && n != 0) {
	  if(n == 1)
	    if(verbose2) printf("\ncontinuous_shift: Write %ld pulses\n", fout->NrSubints-2);
	  /*	    pumawrite(Ipulse, sizeof(float), fout->NrBins, fout); */
	  if(shift >= 0) {
	    if(circularShift == 0 && f == 0)
	      nout -= 1;
	      /*	      printf("Write %ld %d %ld\n", nout, shift, fin.NrBins-shift); */
	    if(writePulsePSRData(fout, nout, p, f, shift, fin.NrBins-shift, Ipulse, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR continuous_shift: Write error (pulse=%ld pol=%ld freq=%ld, start=%d, length=%ld).", nout, p, f, shift, fin.NrBins-shift);
	      return 0;
	    }
	      /*	      printf("Write %ld %d %d\n", nout+1, 0, shift); */
	    if(writePulsePSRData(fout, nout+1, p, f, 0, shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR continuous_shift: Write error (pulse=%ld pol=%ld freq=%ld, start=%d, length=%d).", nout, p, f, 0, shift);
	      return 0;
	    }
	  }else {
	    if(writePulsePSRData(fout, nout, p, f, fin.NrBins-shift, shift, Ipulse, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR continuous_shift: Write error (pulse=%ld pol=%ld freq=%ld, start=%ld, length=%d).", nout, p, f, fin.NrBins-shift, shift);
	      return 0;
	    }
	    if(writePulsePSRData(fout, nout+1, p, f, 0, fin.NrBins-shift, Ipulse+fout->NrBins-shift, verbose) != 1) {
	      fflush(stdout);
	      printerror(verbose.debug, "ERROR continuous_shift: Write error (pulse=%ld pol=%ld freq=%ld, start=%d, length=%ld).", nout, p, f, 0, fin.NrBins-shift);
	      return 0;
	    }
	  }
	}
	if(verbose2 && verbose.nocounters == 0) {
	  printf("continuous_shift: Writing file: %.1f%%     \r", (100.0*(n+1))/(float)(fout->NrSubints-2));
	  /*	  fprintf(stderr, "\nn=%ld, f=%ld\n", n, f); */
	  fflush(stdout);
	}
      }
    }
  }
  free(Ipulse);
  free(Ifirst);
  free(Ilast);
  if(verbose2 && verbose.nocounters == 0)
    printf("\n");
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

// Returns the fold period for subint (counting from zero) in datafile.
// At the moment, only a fixed period folding model is supported, 
// so at the moment, subint is not used.
//
// Returns: 0 = ok
//          1 = not folded, which also results in a warning. (period = 0)
//          2 = error  (period = 0)
int get_period(datafile_definition datafile, long subint, double *period, verbose_definition verbose)
{
  *period = 0;
  if(datafile.isFolded != 1) {
    printwarning(verbose.debug, "WARNING get_period: Data does not appear to be folded");
    return 1;
  }
  if(datafile.foldMode != FOLDMODE_FIXEDPERIOD) {
    printerror(verbose.debug, "ERROR get_period: Unknown folding mode");
    return 2;
  }
  *period = datafile.fixedPeriod;
  return 0;
}

// Returns the fold period for subint (counting from zero) in datafile.
// At the moment, only a fixed period folding model is supported, 
// so at the moment, subint is not used.
double get_tsamp(datafile_definition datafile, long subint, verbose_definition verbose)
{
  if(datafile.tsampMode == TSAMPMODE_LONGITUDELIST) {
    printerror(verbose.debug, "ERROR get_tsamp: Data is not defined to have a regular sampling interval.");
    exit(0);
  }
  if(datafile.tsampMode != TSAMPMODE_FIXEDTSAMP) {
    printerror(verbose.debug, "ERROR get_tsamp: Unknown sampling mode");
    exit(0);
  }
  return datafile.fixedtsamp;
}

// Get the pulse longitude (in degrees)
double get_pulse_longitude(datafile_definition datafile, long subint, long binnr, verbose_definition verbose)
{
  if(datafile.tsampMode == TSAMPMODE_LONGITUDELIST) {
    if(datafile.tsamp_list == NULL) {
      printerror(verbose.debug, "ERROR get_pulse_longitude: Longitudes appear to be undifined in the data");
      exit(0);
    }
    return datafile.tsamp_list[binnr];
  }else if(datafile.tsampMode == TSAMPMODE_FIXEDTSAMP) {
    if(datafile.fixedtsamp == 0) {
      printerror(verbose.debug, "ERROR get_pulse_longitude (%s): Sampling time appears to be set to zero.", datafile.filename);
      exit(0);
    }

    double longitude;
    double period;
    int ret;
    ret = get_period(datafile, subint, &period, verbose);
    if(ret != 0) {
      printerror(verbose.debug, "ERROR get_pulse_longitude (%s): Cannot obtain period", datafile.filename);
      exit(0);
    }
    longitude = binnr;
    longitude *= datafile.fixedtsamp;
    longitude /= period;
    return 360.0*longitude;
  }else {
    printerror(verbose.debug, "ERROR get_pulse_longitude: Unknown sampling mode");
    exit(0);
  }
  return datafile.fixedtsamp;
}

// Attempt to convert for instance a pulse longitude list in a fixed sampling time
// Return: 1 = successful, 0 = No changes could be made
int convert_to_fixed_tsamp(datafile_definition *datafile, verbose_definition verbose)
{
  int ret;
  long binnr;
  double period, delta, delta_exp, diff;
  if(verbose.debug)
    printf("Entering convert_to_fixed_tsamp()\n");
  if(datafile->tsampMode == TSAMPMODE_FIXEDTSAMP) {
    printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: Data already has a fixed sampling time.");
    return 0;
  }
  ret = get_period(*datafile, 0, &period, verbose);
  if(ret == 2) {
    printerror(verbose.debug, "ERROR convert_to_fixed_tsamp (%s): Cannot obtain period", datafile->filename);
    return 0;
  }
  if(period <= 0.0) {
    printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: Period is not set, command is ignored.");
    return 0;
  }
  if(datafile->tsamp_list == NULL) {
    printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: The tsamp_list does not appear to be initialised, command is ignored.");
    return 0;
  }
  if(datafile->NrBins <= 1) {
    printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: More than 1 bin is required, command is ignored.");
    return 0;
  }

  delta_exp = get_pulse_longitude(*datafile, 0, 1, verbose);
  delta_exp -= get_pulse_longitude(*datafile, 0, 0, verbose);
  if(delta_exp <= 0.0) {
    printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: Sampling time seems to be negative or zero, command is ignored.");
    return 0;
  }
  if(datafile->NrBins > 2) {
    for(binnr = 1; binnr < datafile->NrBins-1; binnr++) {
      delta = get_pulse_longitude(*datafile, 0, binnr+1, verbose);
      delta -= get_pulse_longitude(*datafile, 0, binnr, verbose);
      diff = fabs(delta-delta_exp)/delta_exp;
      if(diff > 1.0001 && diff < 0.9999) {
	printwarning(verbose.debug, "WARNING convert_to_fixed_tsamp: Sampling time does not appear to be regular.");
	return 0;
      }
    }
  }

  free(datafile->tsamp_list);
  datafile->tsamp_list = NULL;
  datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
  datafile->fixedtsamp = period*delta_exp/360.0;
  if(verbose.debug) {
    printf("  convert_to_fixed_tsamp: new sampling time is %lf sec\n", datafile->fixedtsamp);
    printf("Exiting convert_to_fixed_tsamp()\n");
  }
  return 1;
}

// Returns the length (in seconds) of the subint (counting from zero)
// in datafile. At the moment, only a fixed subint durations are
// supported, subint is not used.
double get_tsub(datafile_definition datafile, long subint, verbose_definition verbose)
{
  if(datafile.tsubMode == TSUBMODE_FIXEDTSUB) {
    return datafile.tsub_list[0];
  }else if(datafile.tsubMode == TSUBMODE_TSUBLIST) {
    //    fprintf(stderr, "XXXX %lf\n", datafile.tsub_list[subint]);
    return datafile.tsub_list[subint];
  }else {
    printerror(verbose.debug, "ERROR get_tsub: Unknown subint duration mode", datafile.tsubMode);
    if(verbose.debug) {
      printerror(verbose.debug, "Subint duration mode is set to %d", datafile.tsubMode);
    }
    exit(0);
  }
}

// Returns the length of the observation in seconds.
double get_tobs(datafile_definition datafile, verbose_definition verbose)
{
  long i;
  double tobs;

  // The following data products have subints defined, which are not really subints
  if(datafile.gentype == GENTYPE_LRFS || datafile.gentype == GENTYPE_2DFS || datafile.gentype == GENTYPE_S2DFSP3 || datafile.gentype == GENTYPE_S2DFSP2 || datafile.gentype == GENTYPE_P3FOLD || datafile.gentype == GENTYPE_LRCC || datafile.gentype == GENTYPE_RMMAP || datafile.gentype == GENTYPE_PADIST) {
    return get_tsub(datafile, 0, verbose);
  }else if(datafile.gentype == GENTYPE_RECEIVERMODEL || datafile.gentype == GENTYPE_RECEIVERMODEL2) {  // These have no associated observing length defined
    return 0;
  }else {
    tobs = 0;
    for(i = 0; i < datafile.NrSubints; i++)
      tobs += get_tsub(datafile, i, verbose);
    return tobs;
  }
}

// Returns the mid-mjd of the subint (counting from zero)
// in datafile.
long double get_mjd_subint(datafile_definition datafile, long subint, verbose_definition verbose)
{
  int i;
  long double mjd = datafile.mjd_start;
  long double offset;
  offset = 0;
  if(subint > 0) {
    for(i = 0; i < subint; i++)
      offset += get_tsub(datafile, i, verbose);
  }
  offset += 0.5*get_tsub(datafile, subint, verbose);
  mjd += offset/(3600.0*24.0);
  return mjd;
}

// Returns the channel bandwidth (can be negative) for datafile
//
// Return 1 = ok, 0 = failed
int get_channelbandwidth(datafile_definition datafile, double *channelbw, verbose_definition verbose)
{
  //  if(datafile.freqMode != FREQMODE_UNIFORM) {
  //    printerror(verbose.debug, "ERROR get_channelbandwidth: Unknown frequency channel bandwidth (freqMode = %d)", datafile.freqMode);
  //    return 0;
  //  }
  //  if(subint < 0 || channel < 0 || subint >= datafile.NrSubints || channel >= datafile.NrFreqChan) {
    // Sometimes frequency in header is set before number of subints is known.
  //    if(subint != 0 || channel != 0) {
  //      printerror(verbose.debug, "ERROR get_channelbandwidth: subint/channel out of range.");
  //      return 0;
  //    }
  //  }
  if(datafile.isTransposed == 0)
    *channelbw = get_bandwidth(datafile, verbose)/(double)datafile.NrFreqChan;
  else
    *channelbw = get_bandwidth(datafile, verbose)/(double)datafile.NrSubints;
  return 1;
}

// Sets the channel bandwidth (can be negative) for datafile
void set_channelbandwidth(datafile_definition *datafile, double channelbw, verbose_definition verbose)
{
  //  if(datafile->freqMode != FREQMODE_UNIFORM) {
  //    printerror(verbose.debug, "ERROR set_channelbandwidth: Not implemented for frequency mode %d", datafile->freqMode);
  //    exit(0);
  //  }
  //  if(subint < -1 || channel < -1 || subint >= datafile->NrSubints || channel >= datafile->NrFreqChan) {
    // Sometimes frequency in header is set before number of subints is known.
  //    if(subint != 0 || channel != 0) {
  //      printerror(verbose.debug, "ERROR set_channelbandwidth: subint/channel out of range.");
  //      exit(0);
  //    }
  //  }
  //  if(channel != -1 || subint != -1) {
  //    printerror(verbose.debug, "ERROR set_channelbandwidth: Not implemented for individual frequency channels/subints.");
  //    exit(0);
  //  }
  datafile->bandwidth = channelbw*datafile->NrFreqChan;
}

// Returns the overall bandwidth (can be negative) for datafile, which is the difference between the first and last channel, including their widths. The channels are therefore assumed to be in order.
double get_bandwidth(datafile_definition datafile, verbose_definition verbose)
{
  //  if(datafile.freqMode != FREQMODE_UNIFORM) {
  //    printerror(verbose.debug, "ERROR get_bandwidth: Unknown observing frequency");
  //    exit(0);
  //  }
  //  if(subint < 0 || subint >= datafile.NrSubints) {
    // Sometimes frequency in header is set before number of subints is known.
  //    if(subint != 0) {
  //      printerror(verbose.debug, "ERROR set_bandwidth: subint out of range.");
  //      exit(0);
  //    }
  //  }
  //  if(verbose.debug) {
  //    printf("Obtaining bandwidth from array starting at %p\n", datafile.freq_list);
  //  }
  return datafile.bandwidth;
}

// Sets the channel bandwidth (can be negative) for datafile
//
// Return 1 = ok, 0 = failed
int set_bandwidth(datafile_definition *datafile, double bw, verbose_definition verbose)
{
  //  if(datafile->freqMode != FREQMODE_UNIFORM) {
  //    printerror(verbose.debug, "ERROR set_bw: Not implemented for frequency mode %d", datafile->freqMode);
  //    return 0;
  //  }
  //  if(subint < -1 || subint >= datafile->NrSubints) {
    // Sometimes frequency in header is set before number of subints is known.
  //    if(subint != 0) {
  //      printerror(verbose.debug, "ERROR set_bw: subint out of range.");
  //      return 0;
  //    }
  //  }
  //  if(subint != -1) {
  //    printerror(verbose.debug, "ERROR set_bw: Not implemented for individual subints");
  //    return 0;
  //  }
  datafile->bandwidth = bw;
  return 1;
}

// Returns the centre frequency, i.e. (nonweighted) middle of the overall band.
double get_centre_frequency(datafile_definition datafile, verbose_definition verbose)
{
  //  if(datafile.freqMode != FREQMODE_UNIFORM) {
  //    printerror(verbose.debug, "ERROR get_centre_frequency: Unknown observing frequency");
  //    exit(0);
  //  }
  //  if(subint < 0 || subint >= datafile.NrSubints) {
    // Sometimes frequency in header is set before number of subints is known.
  //    if(subint != 0) {
  //      printerror(verbose.debug, "ERROR get_centre_frequency: subint out of range (%ld while number of subints is %ld).", subint, datafile.NrSubints);
  //      exit(0);
  //    }
  //  }
  return datafile.centrefreq;
}

// Set the centre frequency, i.e. (nonweighted) middle of the overall band.
void set_centre_frequency(datafile_definition *datafile, double freq, verbose_definition verbose)
{
  //  if(datafile->freqMode != FREQMODE_UNIFORM) {
  //    printerror(verbose.debug, "ERROR set_centre_frequency: Not implemented for frequency mode %d", datafile->freqMode);
  //    exit(0);
  //  }
  //  if(subint != -1) {
  //    printerror(verbose.debug, "ERROR set_centre_frequency: Not implemented for individual subints.");
  //    exit(0);
  //  }
  //  if(subint < -1 || subint >= datafile->NrSubints) {
    // Sometimes frequency in header is set before number of subints is known.
  //    if(subint != 0) {
  //      printerror(verbose.debug, "ERROR set_centre_frequency: subint out of range (%ld while number of subints is %ld).", subint, datafile->NrSubints);
  //      exit(0);
  //    }
  //  }
  datafile->centrefreq = freq;
}

/* Returns the midpoint of the frequency  channel number in MHz. */
double get_nonweighted_channel_freq(datafile_definition psrdata, long channel, verbose_definition verbose)
{
  double freq_bottom, cfreq;
  //  if(psrdata.freqMode != FREQMODE_UNIFORM) {
  //    printerror(verbose.debug, "ERROR get_nonweighted_channel_freq: Unknown observing frequency");
  //    if(verbose.debug) {
  //      printerror(0, "freqMode is set to %d, expected %d.\n", psrdata.freqMode, FREQMODE_UNIFORM);
  //    }
  //    exit(0);
  //  }
  if(psrdata.isTransposed == 0) {
    if(channel < 0 || channel >= psrdata.NrFreqChan) {
      // Sometimes frequency in header is set before number of subints is known.
      //    if(channel != 0) {
      printerror(verbose.debug, "ERROR get_nonweighted_channel_freq: channel out of range (channel=%ld/%ld).", channel, psrdata.NrFreqChan);
      exit(0);
      //    }
    }
  }else { 
    if(channel < 0 || channel >= psrdata.NrSubints) {
      // Sometimes frequency in header is set before number of subints is known.
      //    if(channel != 0) {
      printerror(verbose.debug, "ERROR get_nonweighted_channel_freq: channel out of range (channel=%ld/%ld).", channel, psrdata.NrFreqChan);
      exit(0);
      //    }
    }
  }

  double chanbw, cfreq_total;
  if(get_channelbandwidth(psrdata, &chanbw, verbose) == 0) {
    printerror(verbose.debug, "ERROR get_nonweighted_channel_freq (%s): Cannot obtain channel bandwidth.", psrdata.filename);
    exit(0);
  }
  cfreq_total = get_centre_frequency(psrdata, verbose);
  if(psrdata.isTransposed == 0) {
    freq_bottom = cfreq_total - 0.5*(double)psrdata.NrFreqChan*chanbw;
    cfreq = freq_bottom + chanbw*(channel+0.5);
    //  printf("XXXXXX cfreq=%lf channelbw=%lf channel=%ld bottom_freq=%lf channel centre freq = %lf\n", psrdata.freq_cent, get_channelbw(psrdata), channel, freq_bottom, cfreq);
  }else {
    freq_bottom = cfreq_total - 0.5*(double)psrdata.NrSubints*chanbw;
    cfreq = freq_bottom + chanbw*(channel+0.5);
  }
  return cfreq;
}

/* Returns the frequency label (if FREQMODE_FREQTABLE) subint/channel number in MHz. Otherwise the same as get_nonweighted_channel_freq(). */
double get_weighted_channel_freq(datafile_definition psrdata, long subint, long channel, verbose_definition verbose)
{
  if(psrdata.freqMode != FREQMODE_FREQTABLE) {
    return get_nonweighted_channel_freq(psrdata, channel, verbose);
  }

  if(subint < 0 || channel < 0 || subint >= psrdata.NrSubints || channel >= psrdata.NrFreqChan) {
    // Sometimes frequency in header is set before number of subints is known.
    //    if(subint != 0 || channel != 0) {
      printerror(verbose.debug, "ERROR get_weighted_channel_freq: subint/channel out of range (subint=%ld/%ld, channel=%ld/%ld).", subint, psrdata.NrSubints, channel, psrdata.NrFreqChan);
      exit(0);
      //    }
  }
  return psrdata.freqlabel_list[subint*psrdata.NrFreqChan+channel];
}

/* Sets the frequency label (if FREQMODE_FREQTABLE) subint/channel number in MHz. Return 0 on error. */
int set_weighted_channel_freq(datafile_definition *psrdata, long subint, long channel, double freq, verbose_definition verbose)
{
  if(psrdata->freqMode != FREQMODE_FREQTABLE) {
    printerror(verbose.debug, "ERROR set_weighted_channel_freq: freqMode is set to %d, expected %d.\n", psrdata->freqMode, FREQMODE_FREQTABLE);
    return 0;
  }
  if(subint < 0 || channel < 0 || subint >= psrdata->NrSubints || channel >= psrdata->NrFreqChan) {
    // Sometimes frequency in header is set before number of subints is known.
    //    if(subint != 0 || channel != 0) {
      printerror(verbose.debug, "ERROR set_weighted_channel_freq: subint/channel out of range (subint=%ld/%ld, channel=%ld/%ld).", subint, psrdata->NrSubints, channel, psrdata->NrFreqChan);
      return 0;
      //    }
  }
  psrdata->freqlabel_list[subint*psrdata->NrFreqChan+channel] = freq;
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Get the parallactic angle in radians for the midpoint of the
   specified subintnr (counting from zero and assumed to be of equal
   length). Set subintnr to -1 to obtain the parallactic angle at the
   midpoint of the observation. The ra and dec are assumed to be not
   yet precessed and to be in J2000. If successful (returns 1,
   otherwise 0).*/
int data_parang(datafile_definition data, long subintnr, double *parang, verbose_definition verbose)
{
  if(fabs(data.telescope_X) < 1e-6 && fabs(data.telescope_Y) < 1e-6 && fabs(data.telescope_Z) < 1e-6) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR data_parang: Position of telescope appears to be undefined.");   
    return 0;
  }
  if(fabs(data.ra) < 1e-6 && fabs(data.dec) < 1e-6) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR data_parang: Position of source appears to be undefined.");   
    return 0;
  }
  //  double telescope_long = observatory_long_geocentric(data);
  //  double telescope_lat = observatory_lat_geocentric(data);
  double telescope_long = observatory_long_geodetic(data);
  double telescope_lat = observatory_lat_geodetic(data);
  long double mjd;
  if(subintnr < 0) {
    mjd = data.mjd_start + 0.5*get_tobs(data, verbose)/(double)(3600.0*24.0);  // Go to middle of observation
  }else {
    mjd = get_mjd_subint(data, subintnr, verbose);
    //    mjd = data.mjd_start;
    //    double subint_length = (data.tobs/(double)(3600.0*24.0))/(double)data.NrSubints;
    //    mjd += subintnr*subint_length; // Skip requested nr of subints
    //    mjd += 0.5*subint_length;      // Go to middle of subint
  }
  *parang = calc_parang(telescope_long, telescope_lat, data.ra, data.dec, mjd, 1);
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* By looking at Stokes I to be either just positive or negative, it
   is determined if the baseline is subtracted or not. Returns 0 if
   the baseline appears not to be subtracted, and 1 if the baseline is
   (potentially) subtracted. */
int check_baseline_subtracted(datafile_definition data, verbose_definition verbose)
{
  long f, n, b;
  float sample, miny, maxy;
  miny = maxy = 0;
  if(data.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR check_baseline_subtracted: Data should be loaded into memory.");
    return 0;
  }
  for(f = 0; f < data.NrFreqChan; f++) {
    for(n = 0; n < data.NrSubints; n++) {
      for(b = 0; b < data.NrBins; b++) {
	if(readPulsePSRData(&data, n, 0, f, b, 1, &sample, verbose) != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR check_baseline_subtracted: Cannot read data.");
	  exit(-1);
	}
	if(b == 0) {
	  miny = sample;
	  maxy = sample;
	}else {
	  if(sample < miny)
	    miny = sample;
	  if(sample > maxy)
	    maxy = sample;
	}
      }
      if(maxy < 0 || miny > 0) {
	return 0;
      }
    }
  }
  return 1;
}

//START REGION DEVELOP

/* This function converts the timeconst (e-fold time) in seconds, as
   defined at frequency reffreq (MHz), varying with frequency with a
   powerlaw-index pwrlawindex (typically -4 for scattering in ISM)
   into a scatter time in bins for frequency channel freqchan. */
double scatter_efoldtime_bins(datafile_definition original, long subint, long freqchan, double timescale, double reffreq, double pwrlawindex, verbose_definition verbose)
{
  double freq, efoldtime;
  freq = get_weighted_channel_freq(original, subint, freqchan, verbose);
  //Get efoldtime in seconds
  efoldtime = timescale*pow(freq/reffreq, pwrlawindex);
  //and in bins
  efoldtime /= get_tsamp(original, 0, verbose);
  if(verbose.debug) {
    printf("e-Fold time at freq=%f MHz = %f bins\n", freq, efoldtime);
  }
  return efoldtime;
}

/* This function takes a pulse of data, npts bins long and convolves
   it with an exponential tail (fast-rise-slow-decay) function. The
   timeconst (e-fold time) is in seconds, as defined at frequency
   reffreq (MHz). The timescale varies with frequency with a
   powerlaw-index pwrlawindex (typically -4 for scattering in
   ISM). Each pulse is treated seperately, so the end of scatter tails
   end up at the start of the same pulse/subint. Return 1 =
   success. Verbose level determines nr of spaces before output. If
   nocounters is set, no progress is shown. */
int convolveScatterTail(datafile_definition original, float timescale, float reffreq, float pwrlawindex, verbose_definition verbose)
{
  long p, f, n;
  int i;
  float efoldtime;
  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Adding scatter tail to dataset (%f sec at %f MHz)\n", timescale, reffreq);
  }
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_convolveScatterTail: only works if data is loaded into memory.");
    return 0;
  }
  if(original.poltype == POLTYPE_ILVPAdPA || original.poltype == POLTYPE_PAdPA || original.poltype == POLTYPE_ILVPAdPATEldEl) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR preprocess_convolveScatterTail: Cannot handle PA data.");
    return 0;
  }

  /* Loop over pulses and frequency channels etc. */
  for(f = 0; f < original.NrFreqChan; f++) {
    for(p = 0; p < original.NrPols; p++) {
      for(n = 0; n < original.NrSubints; n++) {
	efoldtime = scatter_efoldtime_bins(original, n, f, timescale, reffreq, pwrlawindex, verbose);
	if(convolveScatterTail_singlepulse(&(original.data[original.NrBins*(p+original.NrPols*(f+n*original.NrFreqChan))]), original.NrBins, efoldtime, verbose) == 0) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR convolveScatterTail: convolve failed");
	  return 0;
	}
        if(verbose.verbose && verbose.nocounters == 0) {
	  for(i = 0; i < verbose.indent; i++)      
	    printf(" ");
	  printf("processed %.1f%%     \r", (100.0*((f*original.NrPols+p)*original.NrSubints+n))/(float)(original.NrSubints*original.NrFreqChan*original.NrPols));
	  fflush(stdout);
	}
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    if(verbose.verbose) printf(" done                            \n");
  }
  return 1;
}
