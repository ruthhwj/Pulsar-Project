//START REGION RELEASE
#include <math.h>
#include <string.h>
#include "psrsalsa.h"

// The boxcarwidths to use
#define NrBoxCarWidths 29
int BoxCars[NrBoxCarWidths] = {1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,40,50,60,75,100,250,500,1000,1500,2000,2500,3000};

//START REGION DEVELOP
//START REGION RELEASE


/*
This function adds the bins between (and including) bin1 and bin2 by
taking into acount the specified baseline. If squared is set the
signal is squared first.
 */
float integratePulseEnergy(float *pulse, int bin1, int bin2, float baseline, int squared)
{
  int i;
  float E;
  E = 0;
  for(i = bin1; i <= bin2; i++) {
    if(squared != 0)
      E += (pulse[i]-baseline)*(pulse[i]-baseline);
    else
      E += (pulse[i]-baseline);
  }
  return E;
}

//START REGION DEVELOP
//START REGION RELEASE

/*
This function determines the baseline, rms of the offpulse region
OUTSIDE the region defined by onpulse (Set to NULL if not used). If
nodebase is set, the baseline will be assumed to be already
subtracted, which can be useful when a running baseline has been used.
 */
void offpulseStats(float *pulse, int nrBins, float *baseline, float *rms, pulselongitude_regions_definition *onpulse, int nodebase, verbose_definition verbose)
{
  int i, NrOffpulseBins;
  float E;
  E = 0;
  NrOffpulseBins = 0;
  for(i = 0; i < nrBins; i++) {
    if(checkRegions(i, onpulse, 0, verbose) == 0) {
      NrOffpulseBins++;
      E += pulse[i];
    }
  }
  *baseline = E/(float)(NrOffpulseBins);
  if(nodebase)
    *baseline = 0;

  *rms = 0;
  for(i = 0; i < nrBins; i++) {
    if(checkRegions(i, onpulse, 0, verbose) == 0) {
      *rms += (pulse[i]-(*baseline))*(pulse[i]-(*baseline));
    }
  }
  *rms = sqrt(*rms);
  *rms /= sqrt(NrOffpulseBins);
}


//START REGION DEVELOP
//START REGION RELEASE

void boxcarFindpeak_core(int width, float *pulse, int nrBins, int *bin, int *pulsewidth, float *snrbest, float *E_best, int posOrNeg, int squared, float baseline, float rms, int *firsttime, int *allowedWidths, verbose_definition verbose)
{					       
  float E, snr;
  int b, ok;
  long nrtrials;

  if(width < nrBins) {
    if(verbose.verbose)
      printf("boxcarFindpeak: Try out width %d: ", width);
    nrtrials = 0;
    for(b = 0; b < nrBins-width; b++) {
      /*
      ok = 1;
      for(b2 = b; b2 < b+width; b2++) {
	if(checkRegions(b2, &onpulse_search, 0) == 0) {
	  ok = 0;
	}
      }
      */
      ok = 1;
      if(allowedWidths[b] < b+width-1)
	ok = 0;
      if(ok == 1) {
	nrtrials++;
	//	fprintf(stderr, "XXXXX before: %f\n", *snrbest);
	E = integratePulseEnergy(pulse, b, b+(width-1), baseline, squared);
	snr = E/(rms*sqrt(width));
	ok = 0;
	if(posOrNeg == 0) {
	  if(snr > *snrbest || *firsttime == 1)
	    ok = 1;
	}else if(posOrNeg == 1) {
	  if(fabs(snr) > *snrbest || *firsttime == 1) {
	    ok = 1;
	    snr = fabs(snr);
	  }
	}else if(posOrNeg == 2) {
	  if((snr < 0 && fabs(snr) > *snrbest) || *firsttime == 1) {
	    ok = 1;
	    snr = -snr;
	  }
	}else {
	  printerror(verbose.debug, "ERROR boxcarFindpeak: posOrNeg variable is set to an unrecognized value.");
	  exit(0);
	}
	if(ok) {
	  *firsttime = 0;
	  *bin = b;
	  *pulsewidth = width; 
	  *snrbest = snr;
	  if(E_best != NULL) {
	    *E_best = E;
	  }
	  //	  fprintf(stderr, "XXXXX after: %f %f\n", *E_best, *snrbest);
	}
      }
    }
    //Note: If width > 1, I'm not sure how to calculate truely independent trials when sliding by 1 bin at the time
    if(verbose.verbose)
      printf("%ld trials (not necessarily independent).\n", nrtrials);
  }
}

//START REGION DEVELOP
//START REGION RELEASE

/* Does a boxcar search for the highest s/n peak in pulse. onpulse
   (set to NULL if not used) is used to determine off-pulse statistics
   (baseline and rms). Returns the first bin of the box and its
   pulsewidth, the s/n and the integral over the box (if E_best !=
   NULL). If squared is set the signal is squared. posOrNeg defines if
   the algorithm searches for both positive or negative peaks.

   posOrNeg: 0 = positive peaks only
   posOrNeg: 1 = positive & negative peaks
   posOrNeg: 2 = negative peaks only

   If allwidths is set all bin widths are tried instead of a fixed nr
   of widths (extremely slow). If maxwidth >= 0, it indicates that
   only widths smaller than that number should be tried. If refine is
   set, try all bin-widths around the optimum of the finite nr of
   boxcars search if s2n > 2. The only_onpulse parameter is used to
   limit the search range. If it is set to zero the full pulse phase
   range is searched. If it is set to 1 only the onpulse region is
   searched. If it is set to 2 only the first onpulse region is
   considered, while the other onpulse regions are define the region
   excluded from the off-pulse region. If nodebase is set, the
   baseline will be assumed to be already subtracted, which can be
   useful when a running baseline has been used. Returns 1 if
   successful.  */
int boxcarFindpeak(float *pulse, int nrBins, pulselongitude_regions_definition *onpulse, int *bin, int *pulsewidth, float *snrbest, float *E_best, int squared, int posOrNeg, int allwidths, int refine, int maxwidth, int only_onpulse, int nodebase, verbose_definition verbose)
{
  float baseline, rms;
  int w, w1, w2, dw, width, NrWidths, firsttime, *allowedWidths, b, b2;
  pulselongitude_regions_definition onpulse_search;

  if(initPulselongitudeRegion(&onpulse_search, verbose) == 0) {
    printerror(verbose.debug, "ERROR boxcarFindpeak: Initialising onpulse region failed.");
    return 0;
  }
  offpulseStats(pulse, nrBins, &baseline, &rms, onpulse, nodebase, verbose);

  *snrbest = 0;
  firsttime = 1;
  if(allwidths) {
    NrWidths = nrBins;
  }else {
    NrWidths = NrBoxCarWidths;
  }
  /* Default is full range */
  onpulse_search.nrRegions = 1;
  onpulse_search.left_bin[0] = 0;
  onpulse_search.right_bin[0] = nrBins-1;
  onpulse_search.bins_defined[0] = 1;
  if(onpulse != NULL) {
    copyPulselongitudeRegion(*onpulse, &onpulse_search);
    //    memcpy(&onpulse_search, onpulse, sizeof(regions_definition));
    if(only_onpulse == 0) {
      onpulse_search.nrRegions = 1;
      onpulse_search.left_bin[0] = 0;
      onpulse_search.right_bin[0] = nrBins-1;
      onpulse_search.bins_defined[0] = 1;
    }
    if(only_onpulse == 2)
      onpulse_search.nrRegions = 1;
  }else {
    if(only_onpulse != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR boxcarFindpeak: should define onpulse region when limiting search range.");
      return 0;
    }
  }
  allowedWidths = (int *)malloc(nrBins*sizeof(int));
  if(allowedWidths == 0) {
    fflush(stdout);
    printerror(verbose.debug, "boxcarFindpeak: Cannot allocate memory.");
    return 0;
  }
  /* Find out what boxcar widths will work for which start bin
     number. This will greatly speed up the process compared to
     finding this out in the for each trial in boxcarFindpeak_core. */
  for(b = 0; b < nrBins; b++) {
    allowedWidths[b] = -1;
    for(b2 = b; b2 < nrBins; b2++) {
      if(checkRegions(b2, &onpulse_search, 0, verbose) == 0) {
	break;
      }else {
	allowedWidths[b] = b2;
      }
    }
  }

  *pulsewidth = 0;
  for(w = 0; w < NrWidths; w++) {
    if(allwidths)
      width = w;
    else
      width = BoxCars[w];
    if(width > 0) {
      if(maxwidth <= 0 || width <= maxwidth)
	boxcarFindpeak_core(width, pulse, nrBins, bin, pulsewidth, snrbest, E_best, posOrNeg, squared, baseline, rms, &firsttime, allowedWidths, verbose);
    }
  }					       
  /*
  if(refine == 1 && allwidths == 0) {
    for(w = 0; w < NrWidths; w++) {
      if(BoxCars[w] == *pulsewidth)
	break;
    }
    if(w != 0 && w != NrBoxCarWidths-1) {
      printf("Refining %d %d\n", BoxCars[w-1], BoxCars[w+1]);
      for(width = BoxCars[w-1]; width < BoxCars[w+1]; width++) {
	boxcarFindpeak_core(width, pulse, nrBins, bin, pulsewidth, snrbest, E_best, posOrNeg, squared, baseline, rms, &firsttime);	
      }
    }
    }*/
  if(refine == 1 && allwidths == 0) {
    for(w = 0; w < NrWidths; w++) {
      if(BoxCars[w] == *pulsewidth)
	break;
    }
    if(w > 0 && w < NrBoxCarWidths-1) {
      w1 = BoxCars[w-1];
      w2 = BoxCars[w+1];
      if(w2 > nrBins)
	w2 = nrBins;
      if(w1 < 1)
	w1 = 1;
      if(w2-w1 > 1) {
	/*	printf("Refining %d %d\n", w1, w2); */
	do {
	  dw = (w2-w1)/7;
	  if(dw < 1)
	    dw = 1;
	  for(w = w1+dw; w <= w2-dw; ) {
	    if(maxwidth <= 0 || w <= maxwidth)
	      boxcarFindpeak_core(w, pulse, nrBins, bin, pulsewidth, snrbest, E_best, posOrNeg, squared, baseline, rms, &firsttime, allowedWidths, verbose);
	    w += dw;
	  }
	  w1 = *pulsewidth - dw;
	  w2 = *pulsewidth + dw;
	  if(w2 > nrBins)
	    w2 = nrBins;
	  if(w1 < 1)
	    w1 = 1;
	  /*	  printf("  re-Refining %d %d\n", w1, w2); */
	}while(dw != 1);
      }
    }
  }
  freePulselongitudeRegion(&onpulse_search);
  free(allowedWidths);
  return 1;
}

//START REGION DEVELOP

/*
  Determines the single pulse rms when excluded the onpulse region (if
  not set to NULL). You can select the polarization channel as
  well. If nodebase is set, the baseline will be assumed to be already
  subtracted, which can be useful when a running baseline has been
  used. This is the average determined from analysing each
  channel/subint individually.
 */
float determine_singlepulse_rms(datafile_definition psrdata, int polchan, pulselongitude_regions_definition *onpulse, int nodebase, verbose_definition verbose)
{
  long f, n;
  float baseline, rms;
  double rms_tot;
  rms_tot = 0;
  for(f = 0; f < psrdata.NrFreqChan; f++) {
    for(n = 0; n < psrdata.NrSubints; n++) {
      offpulseStats(&(psrdata.data[psrdata.NrBins*(polchan+psrdata.NrPols*(f+n*psrdata.NrFreqChan))]), psrdata.NrBins, &baseline, &rms, onpulse, nodebase, verbose);
      rms_tot += rms*rms;
    }
  }
  rms_tot = sqrt(rms_tot/((double)((psrdata.NrFreqChan)*(psrdata.NrSubints))));
  return rms_tot;
}
//START REGION DEVELOP
