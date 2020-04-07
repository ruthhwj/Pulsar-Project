//START REGION RELEASE
#include <math.h>
#include <string.h>
#include "psrsalsa.h"

#define carousel_interpolation_limit 0.005           // Define the limit on the Gaussian intensity when performing the interpolation. The lower it is, the further the 2D gaussian extends.
#define carousel_from_p3fold_smoothing_multiplier 10 // Number of times to overlay P3 stacks to smooth the output

/* Folds data onto nr_p3_bins with a period foldp3. Offset is a
   positive start phase for the folding in pulse numbers. nrcounts
   will be filled in with a weight that represents the amount of data
   that went into each p3 phase bin. If noNormalise is set the
   resulting fold (map) is not normalised (flattened) by the
   weights. By default each pulse is split over two p3 phase bins,
   depending how accurate it folds on the integer number of p3 phase
   bins. However, if smoothwidth is positive, then the for a given
   pulse the weight for each p3 phase is set by a gaussian with this
   width. As another alternative, if noSmooth is set, all power ends
   up in a single bin. If slope is set (in degrees per bin) a subpulse
   phase slope is subtracted. Offset2 (in degrees) is an additional
   offset2 in subpulse phase to be applied. */
void foldP3_simple(float *data, long nry, long starty, long nrx, float *map, float *nrcounts, int nr_p3_bins, float foldp3, float offset, int noNormalise, int noSmooth, float smoothWidth, float slope, float offset2, 
//START REGION DEVELOP
		  FILE *pulse_list_file, 
//START REGION RELEASE
		   int debug)
{
  float p3, j_frac, weight, weightnext, dp3;
  int i, j, b, jnext;

  /* Clear the map */
  for(i = 0; i < nr_p3_bins; i++) {
    for(b = 0; b < nrx; b++) {
      nrcounts[i*nrx+b] = 0;
      map[i*nrx+b] = 0;
    }
  }

  for(i = starty; i < nry; i++) {
    if(smoothWidth <= 0) {
      for(b = 0; b < nrx; b++) {
	p3 = 360.0*(float)(i+offset)/foldp3-slope*b + offset2;
	p3 = derotate_deg(p3);
	if(p3 == 360.0)  // Sometimes can be exactly 360 deg, resulting in a buffer overflow
	  p3 = 0;
	if(b == 0) {
	  if(debug) {
	    printf("DEBUG foldP3_simple: pulse=%d (block=%ld ... %ld) folded at P3=%f P, with an offset=%f P and additional offset2=%f deg: Subpulse phase of pulse longitude bin 0 = %f deg\n", i, starty, nry-1, foldp3, offset, offset2, p3);
	  }
//START REGION DEVELOP
	  if(pulse_list_file != NULL) {
	    fprintf(pulse_list_file, "%d %f\n", i, p3);
	  }
//START REGION RELEASE
	}
	j = (p3/360.0)*nr_p3_bins;
	if(noSmooth) {
	  j_frac = 0;
	}else {
	  j_frac = (p3/360.0)*nr_p3_bins - j;
	}
	weight = 1.0-j_frac;
	if(j >= nr_p3_bins || j < 0)
	  fprintf(stderr, "ERROR foldP3_simple: Buffer overflow for i=%d, j=%d (nr_p3_bins=%d) p3=%f\n", i, j, nr_p3_bins, p3);
	map[j*nrx+b] += weight*data[i*nrx+b];
	nrcounts[j*nrx+b] += weight;
	
	if(noSmooth == 0) {
	  jnext = j+1;
	  if(jnext == nr_p3_bins)
	    jnext = 0;
	  weightnext = j_frac;
	  map[jnext*nrx+b] += weightnext*data[i*nrx+b];
	  nrcounts[jnext*nrx+b] += weightnext;
	}
      }
      /*    printf("XXXX %f = %d %f, weight = %f, %f, j=%d,%d\n", p3*nr_p3_bins, j, j_frac, weight, weightnext, j, jnext); */
    }else {
      for(j = 0; j < nr_p3_bins; j++) {
	for(b = 0; b < nrx; b++) {
	  p3 = 360.0*(float)(i+offset)/foldp3-slope*b + offset2;
	  p3 = derotate_deg(p3);
	  if(p3 == 360.0)  // Sometimes can be exactly 360 deg, resulting in a buffer overflow
	    p3 = 0;

	  if(b == 0 && j == 0) {
	    if(debug) {
	      printf("DEBUG foldP3_simple: pulse=%d (block=%ld ... %ld) folded at P3=%f P, with an offset=%f P and additional offset2=%f deg: Subpulse phase of pulse longitude bin 0 = %f deg\n", i, starty, nry-1, foldp3, offset, offset2, p3);
	    }
//START REGION DEVELOP
	    if(pulse_list_file != NULL) {
	      fprintf(pulse_list_file, "%d %f\n", i, p3);
	    }
//START REGION RELEASE
	  }
	  p3 *= nr_p3_bins/360.0;
	  dp3 = fabs(p3-j);
	  if(fabs(p3-j+nr_p3_bins) < dp3) {
	    dp3 = fabs(p3-j+nr_p3_bins);
	  }
	  if(fabs(p3-j-nr_p3_bins) < dp3) {
	    dp3 = fabs(p3-j-nr_p3_bins);
	  }
	  weight = exp(-(dp3*dp3/(smoothWidth*smoothWidth)));
	  map[j*nrx+b] += weight*data[i*nrx+b];
	  nrcounts[j*nrx+b] += weight;
	}
      }
    }
  }

  if(noNormalise == 0) {
    for(i = 0; i < nr_p3_bins; i++) {
      for(b = 0; b < nrx; b++) {
	if(nrcounts[i*nrx+b] > 0) {
	  map[i*nrx+b] /= nrcounts[i*nrx+b];
	}
      }
    }
  }
}

//START REGION DEVELOP

//START REGION RELEASE
/* Folds the data (an array of nrx phase bins by nry pulses) on a
   given p3 value (expressed in pulse periods). The resulting map
   (which should already be allocated) has is nrx by nr_p3_bins points
   in size. The resulting map is normalised by the number of times a
   specific bin in the P3 cycle is filled in. If refine is set the
   routine will align the blocks first by doing a cross-correlation in
   order to try to correct for P3 changes etc. If refine is larger
   than 1, it will first produce a template, which then will be used
   to calculate more accurate offsets etc. The larger the value of
   refine, the better the result should be. Set cyclesperblock to set
   the number of P3 cycles that are folded first before attempting to
   align the blocks (number of pulses folded first is
   foldp3*cyclesperblock rounded down). If onpulse is set, it will
   only use the onpulse region to align the blocks. By default each
   pulse is split over two p3 phase bins, depending how accurate it
   folds on the integer number of p3 phase bins. However, if
   smoothWidth is positive, then the for a given pulse the weight for
   each p3 phase is set by a gaussian with this width (in p3 phase
   bins). As another alternative, if noSmooth is set, all power ends
   up in a single bin. If slope is nonzero (in degrees per bin) a
   subpulse phase slope is subtracted. Offset (in degrees) is
   added. If offsetFileName is set, the found offsets are written out
   to a file (readwriteOffsets = 1) or read in and used
   (readwriteOffsets = 2). Note that verbose.debug is only passed on
   to foldP3_simple once the optimum offset has been determined (for
   the given itteration).

   0 = nothing done 
   1 = ok
*/
int foldP3(float *data, long nry, long nrx, float *map, int nr_p3_bins, float foldp3, int refine, int cyclesperblock, int noSmooth, float smoothWidth, float slope, float subpulse_offset, pulselongitude_regions_definition *onpulse
//START REGION DEVELOP
	   , char *offsetFileName, int readwriteOffsets, char *pulse_list_filename, int writeP3foldPulseList
//START REGION RELEASE
, verbose_definition verbose)
{
  float *blockmap, correl, maxcorrel, *template, *nrcounts, *nrcounts_block;
  int i, b, offset, *bestoffset, itt, ok;
  long startpulse, pulsesleft, dN, blockcounter;
//START REGION DEVELOP
  long blockcounter_tmp;
  FILE *offsetFile, *pulse_list_file;
//START REGION RELEASE



  if(cyclesperblock < 1) {
    fflush(stdout);
    printerror(verbose.debug, "foldP3: cyclesperblock (%d) makes no sense", cyclesperblock);
    return 0;
  }

  dN = foldp3*cyclesperblock;   /* Block size that will be folded per step */


//START REGION DEVELOP
  pulse_list_file = NULL;
  if(writeP3foldPulseList) {
    if(verbose.verbose) 
      printf("Writing list of subpulse phases for each pulse to %s\n", pulse_list_filename);
    if(readwriteOffsets != 2) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR foldP3: Cannot open write out the subpulse for each pulse unless the previously determined offsets for each block are being read in from a file rather than being determined from the data.");
      return 0;      
    }
    pulse_list_file = fopen(pulse_list_filename, "w");
    if(pulse_list_file == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR foldP3: Cannot open %s for reading", pulse_list_filename);
      return 0;
    }
  }
  if(offsetFileName != NULL && readwriteOffsets == 2) {
    if(verbose.verbose) 
      printf("Folding data onto P3 using phases from %s\n", offsetFileName);
    offsetFile = fopen(offsetFileName, "r");
    if(offsetFile == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR foldP3: Cannot open %s for reading", offsetFileName);
      return 0;
    }
    if(refine <= 0)
      refine = 1;   /* Set refine, as we are not doing a fixed folding. */
    if(refine > 1) {
      printwarning(verbose.debug, "  WARNING: MULTIPLE ITTERATIONS MAKES NO SENSE WHEN LOADING PHASE OFFSETS FROM A FILE.");
      refine = 1;
    }
  }else 
//START REGION RELEASE
    {
      if(verbose.verbose) 
	printf("Folding data onto P3=%f using %ld subint blocks of data (refine=%d cyclesperblock=%d smoothWidth=%f slope=%f deg/bin offset=%f)\n", foldp3, dN, refine, cyclesperblock, smoothWidth, slope, subpulse_offset);
    }

  nrcounts = (float *)malloc(nrx*nr_p3_bins*sizeof(float));
  if(nrcounts == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "foldP3: cannot allocate memory");
    return 0;
  }
  bestoffset = (int *)malloc((1+nry/dN)*sizeof(int));
  if(bestoffset == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "foldP3: cannot allocate memory");
    return 0;
  }
  

 


  if(refine <= 0) {  /* Just simply fold the data. */
    foldP3_simple(data, nry, 0, nrx, map, nrcounts, nr_p3_bins, foldp3, 0, 0, noSmooth, smoothWidth, slope, subpulse_offset, 
//START REGION DEVELOP
		  pulse_list_file, 
//START REGION RELEASE
		  0*verbose.debug);
  }else {
    /* This will be the map of the current block, which will be added
       to the result after determining a phase offset */
    blockmap = (float *)malloc(nr_p3_bins*nrx*sizeof(float));
    if(blockmap == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "foldP3: cannot allocate memory");
      return 0;
    }
    nrcounts_block = (float *)malloc(nrx*nr_p3_bins*sizeof(float));
    if(nrcounts_block == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "foldP3: cannot allocate memory");
      return 0;
    }
    if(refine > 1) {
      template = (float *)malloc(nr_p3_bins*nrx*sizeof(float));
      if(template == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "foldP3: cannot allocate memory");
	return 0;
      }
    }

    for(itt = 0; itt < refine; itt++) {
      /* Clear the final map */
      for(i = 0; i < nr_p3_bins; i++) {
	for(b = 0; b < nrx; b++) {
	  nrcounts[i*nrx+b] = 0;
	  map[i*nrx+b] = 0;
	}
      }
      
      startpulse = 0;
      pulsesleft = nry;
      blockcounter = 0;
      while(pulsesleft > foldp3*cyclesperblock) {
	
//START REGION DEVELOP
	if(offsetFileName != NULL && readwriteOffsets == 2) {
	  ok = fscanf(offsetFile, "%ld %d\n", &blockcounter_tmp, &bestoffset[blockcounter]);
	  if(ok != 2) {
	    fflush(stdout);
	    printerror(verbose.debug, "\nERROR foldP3: Reading from %s failed", offsetFileName);
	    return 0;
	  }
	  if(blockcounter_tmp != blockcounter) {
	    fflush(stdout);
	    printerror(verbose.debug, "\nERROR foldP3: Unexpected block number: %ld != %ld", blockcounter_tmp, blockcounter);
	    return 0;
	  }
	}else 
//START REGION RELEASE
	{
	  maxcorrel = 0;
	  for(offset = 0; offset < nr_p3_bins; offset++) {
	    foldP3_simple(data, startpulse+dN, startpulse, nrx, blockmap, nrcounts_block, nr_p3_bins, foldp3, offset*foldp3/(float)nr_p3_bins, 1, noSmooth, smoothWidth, slope, subpulse_offset, 
//START REGION DEVELOP
		  pulse_list_file, 
//START REGION RELEASE
0*verbose.debug);
	    correl = 0;
	    for(i = 0; i < nr_p3_bins; i++) {
	      for(b = 0; b < nrx; b++) {
		ok = 0;
		if(onpulse == NULL) {
		  ok = 1;
		}else if(checkRegions(b, onpulse, 0, verbose)) {
		  ok = 1;
		}
		if(ok) {
		  if(itt > 0)
		    correl +=  template[i*nrx+b]*template[i*nrx+b]*blockmap[i*nrx+b]*blockmap[i*nrx+b];
		  else
		    correl +=  map[i*nrx+b]*map[i*nrx+b]*blockmap[i*nrx+b]*blockmap[i*nrx+b];
		}
	      }
	    }
	    /*	  printf("Itt = %d Offset = %d Correl = %f\n", itt, offset, correl); */
	    if(correl > maxcorrel || offset == 0) {
	      maxcorrel = correl;
	      bestoffset[blockcounter] = offset;
	    }
	  
	  }
	  /*	printf("Correl = %f (offset %d)\n", maxcorrel, bestoffset[blockcounter]); */
	}
	
	foldP3_simple(data, startpulse+dN, startpulse, nrx, blockmap, nrcounts_block, nr_p3_bins, foldp3, bestoffset[blockcounter]*foldp3/(float)nr_p3_bins, 0, noSmooth, smoothWidth, slope, subpulse_offset, 
//START REGION DEVELOP
		  pulse_list_file, 
//START REGION RELEASE
verbose.debug);
	for(i = 0; i < nr_p3_bins; i++) {
	  for(b = 0; b < nrx; b++) {
	    nrcounts[i*nrx+b] += nrcounts_block[i*nrx+b];
	    map[i*nrx+b] += (float)blockmap[i*nrx+b];
	  }
	}
	startpulse += dN;
	pulsesleft -= dN;
	if(verbose.verbose && verbose.nocounters == 0) {
	  printf("  Itteration %d/%d, pulse %ld/%ld     \r", itt+1, refine, startpulse, nry);
	  fflush(stdout);
	}
	blockcounter++;
      }
      if(refine > 1) {
	for(i = 0; i < nr_p3_bins; i++) {
	  for(b = 0; b < nrx; b++) {
	    template[i*nrx+b] = (float)map[i*nrx+b];
	  }
	}    
      }  
    }

    if(verbose.verbose && verbose.nocounters == 0) {
      printf("\n");
    }

//START REGION DEVELOP
    if(offsetFileName != NULL && readwriteOffsets == 2) {
      fclose(offsetFile);
    }

    if(offsetFileName != NULL && readwriteOffsets == 1) {
      printf("  Writing out phases to %s\n", offsetFileName);
      offsetFile = fopen(offsetFileName, "w");
      if(offsetFile == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR foldP3: Cannot open %s for writing", offsetFileName);
	return 0;
      }
      startpulse = 0;
      pulsesleft = nry;
      blockcounter = 0;
      while(pulsesleft > dN) {
	startpulse += dN;
	pulsesleft -= dN;
	fprintf(offsetFile, "%ld %d\n", blockcounter, bestoffset[blockcounter]);
	blockcounter++;
      }      
      fclose(offsetFile);
    }
//START REGION RELEASE

    /* Refold the data with found offsets, but now subtract slope */
    /*
    if(slope != 0) {
      if(verbose.verbose && verbose.nocounters == 0) {
	printf("\n  Refold the data to subtract slope\n");
      }
      
      for(i = 0; i < nr_p3_bins; i++) {
	for(b = 0; b < nrx; b++) {
	  nrcounts[i*nrx+b] = 0;
	  map[i*nrx+b] = 0;
	}
      }
      
      startpulse = 0;
      pulsesleft = nry;
      blockcounter = 0;
      while(pulsesleft > dN) {
	foldP3_simple(data, startpulse+dN, startpulse, nrx, blockmap, nrcounts_block, nr_p3_bins, foldp3, bestoffset[blockcounter]*foldp3/(float)nr_p3_bins, 0, noSmooth, smoothWidth, slope, subpulse_offset, 
//START REGION DEVELOP
		  pulse_list_file, 
//START REGION RELEASE
verbose.debug);
	for(i = 0; i < nr_p3_bins; i++) {
	  for(b = 0; b < nrx; b++) {
	    nrcounts[i*nrx+b] += nrcounts_block[i*nrx+b];
	    map[i*nrx+b] += (float)blockmap[i*nrx+b];
	  }
	}
	startpulse += dN;
	pulsesleft -= dN;
	if(verbose.verbose && verbose.nocounters == 0) {
	  printf("   pulse %ld/%ld     \r", startpulse, nry);
	  fflush(stdout);
	}
	blockcounter++;
      }
    }*/

    free(blockmap);
    free(nrcounts_block);
    if(refine > 1) {
      free(template);
    }
  }

//START REGION DEVELOP
  if(pulse_list_file != NULL) {
    fclose(pulse_list_file);
  }
//START REGION RELEASE
  free(nrcounts);
  free(bestoffset);
  if(verbose.verbose)
    printf("  Done\n");
  return 1;
}

//START REGION DEVELOP

float deg2rad(float degrees) { // Function to convert from degrees to radians
  float radians = degrees * (M_PI / 180.0);
  return radians;
}


float rad2deg(float radians) { // Function to convert from radians to degrees
  float degrees = radians * (180.0 / M_PI);
  return degrees;
}

float RCalc(float phi, carousel_parameter_list plist) { // Function to calculate R
  float tempR = asin( sqrt((sin((phi-plist.phi0)/2.0)*sin((phi-plist.phi0)/2.0)*sin(plist.alpha)*sin(plist.alpha+plist.beta))+(sin(plist.beta/2.0)*sin(plist.beta/2.0))) );
  float R = 2.0*tempR;
  return R;
}


float ThetaCalc(int k, float phi, float R, carousel_parameter_list plist) { // Function to calculate Theta
  float theta_rot, theta_trans;
  theta_rot = (2 * M_PI / plist.P3_hat) * ( k - plist.k0 + ((phi - plist.phi0)/(2 * M_PI)));
  theta_trans = asin(  sin(plist.alpha + plist.beta) * sin(phi - plist.phi0) / sin(R) );		
  float Theta = plist.A*theta_rot + plist.B*theta_trans;
  return Theta;
}

void stack2car(int subint, int bin, float *px, float *py, carousel_parameter_list plist) {
  // This function takes bin and subint numbers and
  // maps them to R and Theta on the carousel, before
  // mapping those in turn to a position in the 
  // output image.

  float phi, R, Theta, x, y;

  // first calculate the pulse phase phi, and convert it to radians:
  phi = ((float)bin - (((float) plist.NrBins - 1.0) / 2.0)) * (plist.subint_width / (float) plist.NrBins);
  phi = deg2rad(phi);

  // calculate R and Theta
  R = RCalc(phi, plist);
  Theta = ThetaCalc(subint, phi, R, plist);

  // Now map R and Theta onto a point in Cartesian coordinates.
  // Let us define increasing Theta in an anticlockwise direction from the negative y-axis.

  x = R * cos(Theta - (M_PI/2.0));
  y = R * sin(Theta - (M_PI/2.0));

  // convert x and y to points in relation to the extent (R_max) of the image
  x = x / plist.R_max;
  y = y / plist.R_max;

  // Now convert x and y to a point on the image of (SIZE)x(SIZE) pixels
  // Shift x and y so they lie between 0 and 1 
  x = (x + 1.0)/2.0;
  y = (y + 1.0)/2.0;
  // Map that to a pixel number (indexing starts at 0)
  px[0] = x * (float)(plist.SIZE - 1);
  py[0] = y * (float)(plist.SIZE - 1);
}

int roundNumber(float number) {
  int integerPart = (int) number;
  float fractionalPart = number - (float) integerPart;
  if (fractionalPart >= 0.5) {
    return (int) ceil(number);
  } else {
    return (int) floor(number);
  }
}


int indexconv4D(int subint, int freq, int pol, int bin, int NrFreqs, int NrPols, int NrBins) { // Function to convert a 3D coordinate [subint, freq, pol, bin] to an index in a 1D array
  int INDEX = bin + NrBins*(pol + NrPols*(freq + NrFreqs * subint));
  return INDEX;
}

float interp_multiplier(int i, int j, float x, float y, float sigma) {
  float multiplier = (1.0 / (sigma * sqrt( 2.0 * M_PI))) * exp( -((i + 0.5 - x)*(i + 0.5 - x) + (j + 0.5 - y)*(j + 0.5 - y))/(2*sigma*sigma));
  return multiplier;
}

void placeGaussian(float x, float y, int bin, int subint, float *output_image, float *artefact_map, datafile_definition datain, carousel_parameter_list plist) {  
  int gaussx, gaussy, pol, freq;
  float sigma = 1.0;
  int delta = roundNumber(sqrt(-2.0 * sigma * sigma * log(carousel_interpolation_limit)));
  int xx = roundNumber(x); // Round x and y to the nearest integer
  int yy = roundNumber(y);
  
  for (gaussx = (xx - delta); gaussx < (xx + delta+1); gaussx++) {
    for (gaussy = (yy - delta); gaussy < (yy + delta+1); gaussy++) {
      if (gaussx >=0 && gaussx< plist.SIZE && gaussy >=0 && gaussy < plist.SIZE) {
	// Do all polarisation and frequency channels whilst we're at it!
	for(pol = 0; pol < plist.NrPols; pol++){
	  for(freq = 0; freq < plist.NrFreqs; freq++){
	    int outindex = indexconv4D(gaussy, freq, pol, gaussx, plist.NrFreqs, plist.NrPols, plist.SIZE);
	    int inindex = indexconv4D(subint, freq, pol, bin, plist.NrFreqs, plist.NrPols, plist.NrBins);
	    //	    if(plist.debug) {
	    //	      printf("XXXXXX interp_multiplier=%f data=%f -> index=%d\n", interp_multiplier(gaussx, gaussy, x, y, sigma), datain.data[inindex], outindex);
	    //	    }
	    output_image[outindex] += datain.data[inindex]*interp_multiplier(gaussx, gaussy, x, y, sigma);
	    artefact_map[outindex] += 1.0*interp_multiplier(gaussx, gaussy,x, y, sigma);
	  }
	}
      }
    }
  }
}

float max_element(float *array, int ElementNr) {
  int i;
  float largest = array[0];
  for (i = 1; i < ElementNr; i++) {
    if (array[i] > largest) {
      largest = array[i];
    }
  }
  return largest;
}

float plotCarousel(carousel_parameter_list plist, char *filename_ptr, float *output_image, datafile_definition datain) {
  int i;
  // Convert all angles to radians
  plist.alpha = deg2rad(plist.alpha);
  plist.beta = deg2rad(plist.beta);
  plist.phi0 = deg2rad(plist.phi0);
  plist.R_max = deg2rad(plist.R_max);  

  // Calculate the angular width of each subint
  float subint_period = plist.NrBins * plist.sampling_period;
  plist.subint_width = 360.0 * (subint_period / plist.P1);

  // Generate an  array for the artefact map
  float *artefact_map;
  artefact_map = (float *)malloc(sizeof(float)*plist.SIZE*plist.SIZE*plist.NrPols*plist.NrFreqs);
  if (artefact_map == NULL) {
    printerror(plist.debug, "ERROR: Could not assign memory for the artefact map.");
  }
  for(i = 0; i < plist.SIZE*plist.SIZE*plist.NrPols*plist.NrFreqs; i++) { // Fill it with zeroes
    artefact_map[i] = 0.0;
    output_image[i] = 0.0;
  }

  // Calculate the (minimum) number of pulses required to plot the entire carousel
  int subints2plot;
  if(plist.p3car == 1){
    subints2plot = ceil(plist.N * plist.P3)*carousel_from_p3fold_smoothing_multiplier;
    printf("Plotting %d pulses to generate the full carousel.\n", subints2plot);
  } else {
    subints2plot = plist.num_pulses;
  }
  // START MAPPING BLOCK ------------------------------------------------------------------------
  if(plist.debug) printf("    Mapping the raw data\n");
  int bin, subint;
  float px;
  float py;

  for(subint = 0; subint < subints2plot; subint++) { // Loop through subints
    if(plist.app_pointer->verbose_state.nocounters == 0) {
      printf("%.2f%% complete...          \r", ((float)subint*100.0/(float)subints2plot));
      fflush(stdout);
    }

    if(plist.p3car == 1){
      plist.whichsubint = (int) plist.NrSubints * ((subint/plist.P3) - (int)(subint/plist.P3));
    }else{
      plist.whichsubint = plist.first_pulse + subint;
    }

    for(bin = 0; bin < plist.NrBins; bin++) { // Loop through bins
	stack2car(subint, bin, &px, &py, plist);
	//	if(plist.debug) printf("    Carousel mapping: subint=%d, bin=%d, px=%f, py=%f\n", subint, bin, px, py);
	placeGaussian(px, py, bin, plist.whichsubint, output_image, artefact_map, datain, plist);
    }
  }

  // Clean the image and return it
  for(i = 0; i < plist.SIZE*plist.SIZE*plist.NrPols*plist.NrFreqs; i++) {
    if(artefact_map[i] !=0 ){
      output_image[i] /=  artefact_map[i];
    }
  }

  free(artefact_map);
  // Plot the data
  //pgplotCAROUSEL(output_image, plist, filename_ptr, "?"); // He's dead Jim.
  

  // Free up memory and return the maximum element
  float mxelmnt = max_element(output_image, plist.SIZE*plist.SIZE);
  return mxelmnt;
}
