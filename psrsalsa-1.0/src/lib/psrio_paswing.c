#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_sort_float.h>
#include <gsl/gsl_statistics_float.h>
#include <gsl/gsl_integration.h>
#include "psrsalsa.h"

/* Filters out all non-significant PA-points, i.e. those with an PA
   error which is negative. Returns the total number of significant
   points. */
int filterPApoints(datafile_definition *datafile, verbose_definition verbose)
{
  int dPa_polnr;
  long i, j, nrpoints;
  float *olddata;
  if(datafile->poltype != POLTYPE_ILVPAdPA && datafile->poltype != POLTYPE_PAdPA && datafile->poltype != POLTYPE_ILVPAdPATEldEl) {
    printerror(verbose.debug, "ERROR filterPApoints: Data doesn't appear to have poltype ILVPAdPA, PAdPA or ILVPAdPATEldEl."); 
    return 0;
  }
  if(datafile->poltype == POLTYPE_ILVPAdPA && datafile->NrPols != 5) {
    printerror(verbose.debug, "ERROR filterPApoints: 5 polarization channels were expected, but there are only %ld.", datafile->NrPols); 
    return 0;
  }else if(datafile->poltype == POLTYPE_ILVPAdPATEldEl && datafile->NrPols != 8) {
    printerror(verbose.debug, "ERROR filterPApoints: 8 polarization channels were expected, but there are only %ld.", datafile->NrPols); 
    return 0;
  }else if(datafile->poltype == POLTYPE_PAdPA && datafile->NrPols != 2) {
    printerror(verbose.debug, "ERROR filterPApoints: 2 polarization channels were expected, but there are only %ld.", datafile->NrPols); 
    return 0;
  }
  if(datafile->NrSubints > 1 || datafile->NrFreqChan > 1) { // As otherwise different bins are removed from different subints.
    printerror(verbose.debug, "ERROR filterPApoints: Can only do this opperation if there is one subint and one frequency channel."); 
    return 0;
  }
  if(datafile->tsampMode != TSAMPMODE_LONGITUDELIST) {
    printerror(verbose.debug, "ERROR filterPApoints: Expected pulse longitudes to be defined."); 
    return 0;
  }

  if(datafile->poltype == POLTYPE_ILVPAdPA || datafile->poltype == POLTYPE_ILVPAdPATEldEl) {
    dPa_polnr = 4;
  }else if(datafile->poltype == POLTYPE_PAdPA) {
    dPa_polnr = 1;
  }


  // Determine the number of points to keep
  nrpoints = 0;
  for(i = 0; i < datafile->NrBins; i++) {
    if(datafile->data[i+dPa_polnr*datafile->NrBins] > 0) {
      nrpoints++;
    }
  }
  if(verbose.verbose)
    printf("Keeping %ld significant PA points\n", nrpoints);

  olddata = datafile->data;
  datafile->data = (float *)malloc(nrpoints*datafile->NrPols*sizeof(float));
  if(datafile->data == NULL) {
    printerror(verbose.debug, "ERROR filterPApoints: Memory allocation error."); 
    return 0;
  }

  // alloc new data, but tsamp_list not necessary, as list is equal or smaller

  j = 0;
  for(i = 0; i < datafile->NrBins; i++) {
    if(olddata[i+dPa_polnr*datafile->NrBins] > 0) {
      if(datafile->poltype == POLTYPE_ILVPAdPA) {
	datafile->data[j+0*nrpoints] = olddata[i+0*datafile->NrBins];
	datafile->data[j+1*nrpoints] = olddata[i+1*datafile->NrBins];
	datafile->data[j+2*nrpoints] = olddata[i+2*datafile->NrBins];
	datafile->data[j+3*nrpoints] = olddata[i+3*datafile->NrBins];
	datafile->data[j+4*nrpoints] = olddata[i+4*datafile->NrBins];
      }else if(datafile->poltype == POLTYPE_ILVPAdPATEldEl) {
	datafile->data[j+0*nrpoints] = olddata[i+0*datafile->NrBins];
	datafile->data[j+1*nrpoints] = olddata[i+1*datafile->NrBins];
	datafile->data[j+2*nrpoints] = olddata[i+2*datafile->NrBins];
	datafile->data[j+3*nrpoints] = olddata[i+3*datafile->NrBins];
	datafile->data[j+4*nrpoints] = olddata[i+4*datafile->NrBins];
	datafile->data[j+5*nrpoints] = olddata[i+5*datafile->NrBins];
	datafile->data[j+6*nrpoints] = olddata[i+6*datafile->NrBins];
	datafile->data[j+7*nrpoints] = olddata[i+7*datafile->NrBins];
      }else {
	datafile->data[j+0*nrpoints] = olddata[i+0*datafile->NrBins];
	datafile->data[j+1*nrpoints] = olddata[i+1*datafile->NrBins];
      }
      datafile->tsamp_list[j] = datafile->tsamp_list[i];
      j++;
    }
  }
  free(olddata);
  datafile->NrBins = nrpoints;
  return datafile->NrBins;
}

//START REGION RELEASE

/* Provide Q and U data and the onpulse region. If extended is
   non-zero, the total polarization sqrt(Q^2+U^2+V^2) and ellipticity
   + error will be computed as well. Each subint and frequency channel
   is calculated independent of each other. Set normalize if you want
   to normalize intensities to the peak value. Set correctQV to -1 if
   you want to swap sign of Q and V (or to +1 if you don't want to
   make a change). Stokes V is divided by correctV. The pulse
   longitudes are calculated (i.e. data will be in
   TSAMPMODE_LONGITUDELIST, unless nolongitudes in nonzero). loffset
   can be used to shift the data (has no effect in nolongitudes is
   used). Whatever extended is set to, it will always calculate total
   degree of linear polarization and calculates pulse longitudes. By
   default (correctLbias = 0) the median value of the offpulse L is
   subtracted from the data to do the "L bias" correction. If
   correctLbias is set to 1 the much better Wardle & Kronberg debias
   method is used. If set to -1 no bias is subtracted at all. paoffset
   (in degrees) is an offset that will be applied to all PA values. If
   rms_file is not the NULL pointer, the off-pulse rmses are
   determined from this dataset rather than from datafile. In both
   cases onpulse is used to identify the offpulse region in the
   relevant dataset. When applying the rmses on the actual dataset to
   be processed, a scaling can be applied via the variable
   rebin_factor, which indicates how many bins should be summed to get
   a sampling time of datafile. Returns 0 if
   not successful. */
int make_paswing_fromIQUV(datafile_definition *datafile, int extended, pulselongitude_regions_definition onpulse, int normalize, int correctLbias, float correctQV, float correctV, int nolongitudes, float loffset, float paoffset, datafile_definition *rms_file, float rebin_factor, verbose_definition verbose)
{
  int indent, rms_file_specified;
  long i, j, NrOffpulseBins, pulsenr, freqnr, output_nr_pols;
  float ymax, baseline_intensity, RMSQ, RMSU, *Loffpulse, *Poffpulse, medianL, medianP, *newdata, *newdata_rms;

  rms_file_specified = 1;
  if(rms_file == NULL) {  // Make rms_file point to the actual datafile when no seperate offpulse file is provided
    rms_file = datafile;
    rms_file_specified = 0;
  }

  if(verbose.verbose) {
    for(indent = 0; indent < verbose.indent; indent++) printf(" ");
    printf("Constructing PA and degree of linear polarization");
    if(extended)
      printf(", total polarization and ellipticity");
    if(rms_file_specified)
      printf(" (using a seperate file to determine the off-pulse rms)");
    printf("\n");
    for(indent = 0; indent < verbose.indent; indent++) printf(" ");
    printf("  Reference frequency for PA is ");
    if(datafile->isDeFarad) {
      if((datafile->freq_ref > -1.1 && datafile->freq_ref < -0.9) || (datafile->freq_ref > 0.99e10 && datafile->freq_ref < 1.01e10))
	printf("infinity\n");
      else if(datafile->freq_ref < 0)
	printf("unknown\n");
      else
	printf("%f MHz\n", datafile->freq_ref);
    }else {
      if(datafile->NrFreqChan == 1)
	printf("%lf MHz\n", get_centre_frequency(*datafile, verbose));
      else
	printf("observing frequencies of individual frequency channels\n");
    }
    for(indent = 0; indent < verbose.indent; indent++) printf(" ");
    printf("  ");
    switch(correctLbias) {
    case -1: printf("No L de-bias applied"); break;
    case 0: printf("De-bias L using median noise subtraction"); break;
    case 1: printf("De-bias L using Wardle & Kronberg correction"); break;
    default: printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Undefined L de-bias method specified."); return 0;
    }
    if(correctQV != 1 || correctV != 1)
      printf(", Q correction factor %f, V correction factor %f", 1.0/correctQV, 1.0/(correctQV*correctV));
    if(normalize)
      printf(", output is normalised");
    if(loffset != 0)
      printf(", pulse longitude shifted by %f deg\n", loffset);
    if(paoffset != 0)
      printf(", PA shifted by %f deg\n", paoffset);
    printf("\n");
    if(extended) {
      printwarning(verbose.debug, "WARNING make_paswing_fromIQUV: Total polarization is computed, and the median off-pulse value is subtracted (probably not a very good idea)."); 
    }
  }

  if(datafile->NrPols != 4 || rms_file->NrPols != 4) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Expected 4 input polarizations.");
    return 0;
  }

  if(datafile->poltype != POLTYPE_STOKES) {
    if(datafile->poltype == POLTYPE_UNKNOWN) {
      printwarning(verbose.debug, "WARNING make_paswing_fromIQUV: Polarization state unknown, it is assumed the data are Stokes parameters.");
    }else {
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Convert data into Stokes parameters first.");
      return 0;
    }
  }
  if(rms_file_specified) {
    if(rms_file->poltype != POLTYPE_STOKES) {
      if(rms_file->poltype == POLTYPE_UNKNOWN) {
	printwarning(verbose.debug, "WARNING make_paswing_fromIQUV: Polarization state of the data to be used to determine the off-pulse rms is unknown, it is assumed the data are Stokes parameters.");
      }else {
	printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Convert data to be used to determine the off-pulse rms into Stokes parameters first.");
	return 0;
      }
    }
  }

  //  if((datafile->NrFreqChan > 1 || datafile->NrSubints > 1) && datafile->freqMode != FREQMODE_UNIFORM) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "WARNING make_paswing_fromIQUV: The frequency channels are not necessarily uniformly separated, or the same from subint to subint.");
  //  }

  if(datafile->tsampMode != TSAMPMODE_FIXEDTSAMP) {   // Otherwise tsamp_list is already in use, so cannot generate it from scratch. Could swap arrays, like data array.
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Expects input to have a regular sampling.");
    return 0;
  }


  if(correctQV == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: correctQV is set to zero, you probably want this to be 1.");
    return 0;    
  }

  if(correctV == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: correctV is set to zero, you probably want this to be 1.");
    return 0;    
  }

  if(datafile->isDebase == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Please remove baseline first, i.e. use pmod -debase.");
    return 0;
  }else if(datafile->isDebase != 1) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING make_paswing_fromIQUV: Unknown baseline state. It is assumed the baseline has already removed from the data.");
  }

  if(rms_file_specified) {
    if(rms_file->isDebase == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Please remove baseline first, i.e. use pmod -debase.");
      return 0;
    }else if(rms_file->isDebase != 1) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING make_paswing_fromIQUV: Unknown baseline state. It is assumed the baseline has already removed from the data.");
    }
  }

  if(rms_file_specified) {
    if(datafile->NrSubints != rms_file->NrSubints) {
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Number of subintegrations is different in the data to be used to determine the off-pulse rms.");
	return 0;
    }
    if(datafile->NrFreqChan != rms_file->NrFreqChan) {
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Number of frequency channels is different in the data to be used to determine the off-pulse rms (%ld != %ld).", rms_file->NrFreqChan, datafile->NrFreqChan);
	return 0;
    }
    if(correctLbias == 0) {   // The calculation of the median of L is not necessarily trivial
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Subtracting the median of L is not supported when a separate file is used for the offpulse statistics.");
      return 0;
    }
    if(extended) {   // The calculation of the median of P is not necessarily trivial
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Subtracting the median of P is not supported when a separate file is used for the offpulse statistics.");
      return 0;
    }
  }

  /* Allocate memory */
  if(datafile->offpulse_rms != NULL) {
    free(datafile->offpulse_rms);
  }
  Loffpulse = (float *)malloc((rms_file->NrBins)*sizeof(float));
  Poffpulse = (float *)malloc((rms_file->NrBins)*sizeof(float));
  if(extended) {
    output_nr_pols = 8;
  }else {
    output_nr_pols = 5;
  }
  newdata = (float *)malloc(datafile->NrBins*datafile->NrSubints*datafile->NrFreqChan*output_nr_pols*sizeof(float));  // Alloc memory to hold I,L,V,Pa,dPa,tot pol, Ell, dEll
  if(rms_file_specified) {
    newdata_rms = (float *)malloc(rms_file->NrBins*rms_file->NrSubints*rms_file->NrFreqChan*output_nr_pols*sizeof(float));  // Alloc memory to hold I,L,V,Pa,dPa,tot pol, Ell, dEll
  }else {
    newdata_rms = newdata;  // Make newdata_rms point to wherever the offpulse L and T are stored
  }
  datafile->offpulse_rms = (float *)malloc(datafile->NrSubints*datafile->NrFreqChan*output_nr_pols*sizeof(float));

  if(Loffpulse == NULL  || Poffpulse  == NULL || newdata == NULL || datafile->offpulse_rms == NULL || newdata_rms == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Memory allocation error.");
    return 0;
  }
  if(nolongitudes == 0) {
    datafile->tsamp_list = (double *)malloc(datafile->NrBins*sizeof(double));
    if(datafile->tsamp_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_paswing_fromIQUV: Memory allocation error.");
      return 0;
    }
    // Generate pulse longitudes
    for(j = 0; j < (datafile->NrBins); j++) {
      datafile->tsamp_list[j] = get_pulse_longitude(*datafile, 0, j, verbose);
      datafile->tsamp_list[j] += loffset;
      //    fprintf(stderr, "XXXXX generating long %lf (j=%ld)\n", datafile->tsamp_list[j], j);
    }
  }



  if(normalize && (datafile->NrSubints > 1 || datafile->NrFreqChan > 1)) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING make_paswing_fromIQUV: Normalization will cause all subintegrations/frequency channels to be normalised individually. This may not be desired.");
  }

  long sindex_I, sindex_Q, sindex_U, sindex_V, sindex_I_rms, sindex_Q_rms, sindex_U_rms, sindex_V_rms, newindex_I, newindex_L, newindex_V, newindex_Pa, newindex_dPa, newindex_T, newindex_Ell, newindex_dEll, newindex_L_rms, newindex_T_rms;

  /* Now loop over the subintegrations and calculate PA etc. */
  for(pulsenr = 0; pulsenr < datafile->NrSubints; pulsenr++) {
    for(freqnr = 0; freqnr < datafile->NrFreqChan; freqnr++) {

      // Calculate start index of pulse in Stokes I-V in the data to make code more readable and potentially a bit faster
      sindex_I = datafile->NrBins*(0+datafile->NrPols*(freqnr+pulsenr*datafile->NrFreqChan));
      sindex_Q = datafile->NrBins*(1+datafile->NrPols*(freqnr+pulsenr*datafile->NrFreqChan));
      sindex_U = datafile->NrBins*(2+datafile->NrPols*(freqnr+pulsenr*datafile->NrFreqChan));
      sindex_V = datafile->NrBins*(3+datafile->NrPols*(freqnr+pulsenr*datafile->NrFreqChan));
      newindex_I   = datafile->NrBins*(0+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr));
      newindex_L   = datafile->NrBins*(1+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr));
      newindex_V   = datafile->NrBins*(2+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr));
      newindex_Pa  = datafile->NrBins*(3+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr));
      newindex_dPa = datafile->NrBins*(4+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr));
      if(extended) {
	newindex_T    = datafile->NrBins*(5+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr));
	newindex_Ell  = datafile->NrBins*(6+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr));
	newindex_dEll = datafile->NrBins*(7+output_nr_pols*(freqnr+datafile->NrFreqChan*pulsenr));
      }
      if(rms_file_specified) {
	sindex_I_rms = rms_file->NrBins*(0+rms_file->NrPols*(freqnr+pulsenr*rms_file->NrFreqChan));
	sindex_Q_rms = rms_file->NrBins*(1+rms_file->NrPols*(freqnr+pulsenr*rms_file->NrFreqChan));
	sindex_U_rms = rms_file->NrBins*(2+rms_file->NrPols*(freqnr+pulsenr*rms_file->NrFreqChan));
	sindex_V_rms = rms_file->NrBins*(3+rms_file->NrPols*(freqnr+pulsenr*rms_file->NrFreqChan));
	newindex_L_rms   = rms_file->NrBins*(1+output_nr_pols*(freqnr+rms_file->NrFreqChan*pulsenr));
	if(extended) {
	  newindex_T_rms    = rms_file->NrBins*(5+output_nr_pols*(freqnr+rms_file->NrFreqChan*pulsenr));
	}
      }else {
	sindex_I_rms = sindex_I;
	sindex_Q_rms = sindex_Q;
	sindex_U_rms = sindex_U;
	sindex_V_rms = sindex_V;
	newindex_L_rms   = newindex_L;
	if(extended) {
	  newindex_T_rms    = newindex_T;
	}
      }
      
      /* Determine normalisation constant if required */
      if(normalize == 0) {
	ymax = 1;
      }else {
	ymax = datafile->data[sindex_I];
	for(j = 1; j < (datafile->NrBins); j++) {
	  if(datafile->data[sindex_I+j] > ymax)
	    ymax = datafile->data[sindex_I+j];
	}
      }
      if(ymax == 0)  // Disable normalisation of zapped channels/subints
	ymax = 1;
      
      for(j = 0; j < (datafile->NrBins); j++) {
	datafile->data[sindex_I + j] /= ymax;
	datafile->data[sindex_Q + j] /= correctQV*ymax;
	datafile->data[sindex_U + j] /= ymax;
	datafile->data[sindex_V + j] /= correctV*correctQV*ymax;
	// L = sqrt(Q^2+U^2)
	newdata[j+newindex_L] = sqrt((datafile->data[sindex_Q+j])*(datafile->data[sindex_Q+j])+(datafile->data[sindex_U+j])*(datafile->data[sindex_U+j]));
	if(extended) {
	  // totpol = sqrt(Q^2+U^2+V^2)
	  newdata[j+newindex_T] = sqrt(newdata[j+newindex_L]*newdata[j+newindex_L]+datafile->data[sindex_V+j]*datafile->data[sindex_V+j]);
	}
      }
      // Scale the data to be used to determine the offpulse rms in the same way
      if(rms_file_specified) {
	for(j = 0; j < (rms_file->NrBins); j++) {
	  rms_file->data[sindex_I_rms + j] /= ymax;
	  rms_file->data[sindex_Q_rms + j] /= correctQV*ymax;
	  rms_file->data[sindex_U_rms + j] /= ymax;
	  rms_file->data[sindex_V_rms + j] /= correctV*correctQV*ymax;
	  // L = sqrt(Q^2+U^2)
	  newdata_rms[j+newindex_L_rms] = sqrt((rms_file->data[sindex_Q_rms+j])*(rms_file->data[sindex_Q_rms+j])+(rms_file->data[sindex_U_rms+j])*(rms_file->data[sindex_U_rms+j]));
	  if(extended) {
	  // totpol = sqrt(Q^2+U^2+V^2)
	    newdata_rms[j+newindex_T_rms] = sqrt(newdata[j+newindex_L_rms]*newdata[j+newindex_L_rms]+rms_file->data[sindex_V_rms+j]*rms_file->data[sindex_V_rms+j]);
	  }
	}
      }

      // 5 polarizations are defined, only three are used
      datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = 0;  /* I */
      datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = 0;  /* L */
      datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = 0;  /* V */
      datafile->offpulse_rms[3+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = -1;  /* Pa */
      datafile->offpulse_rms[4+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = -1;  /* dPa */
      if(extended) {
	datafile->offpulse_rms[5+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = 0;  /* totpol */
	datafile->offpulse_rms[6+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = -1;  /* Ell */
	datafile->offpulse_rms[7+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = -1;  /* dEll */
      }

      baseline_intensity = 0;
      RMSQ = 0;
      RMSU = 0;
      NrOffpulseBins = 0;

      for(i = 0; i < (rms_file->NrBins); i++) {
	Loffpulse[i] = 0;  // This array will be L as function of pulse longitude after cutting out the on-pulse region, which allows the median offpulse value to be calculated
	if(checkRegions(i, &onpulse, 0, verbose) == 0) {
	  NrOffpulseBins++;
	  baseline_intensity += rms_file->data[sindex_I_rms+i];
	  // RMS in Q & U
	  RMSQ += (rms_file->data[sindex_Q_rms+i])*(rms_file->data[sindex_Q_rms+i]);
	  RMSU += (rms_file->data[sindex_U_rms+i])*(rms_file->data[sindex_U_rms+i]);
	  // RMS in I, L, V & totpol
	  datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] += (rms_file->data[sindex_I_rms+i])*(rms_file->data[sindex_I_rms+i]);
	  datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] += (newdata_rms[i+newindex_L_rms])*(newdata_rms[i+newindex_L_rms]);
	  datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] += (rms_file->data[sindex_V_rms+i])*(rms_file->data[sindex_V_rms+i]);
	  if(extended) {
	    datafile->offpulse_rms[5+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] += (newdata_rms[i+newindex_L_rms])*(newdata_rms[i+newindex_L_rms]) + (rms_file->data[sindex_V_rms+i])*(rms_file->data[sindex_V_rms+i]);
	  }
	  // Keep a copy of only the L in the off-pulse region
	  Loffpulse[NrOffpulseBins-1] = newdata_rms[i+newindex_L_rms];
	  if(extended)
	    Poffpulse[NrOffpulseBins-1] = newdata_rms[i+newindex_T_rms]; 
	}
      }
      baseline_intensity /= (float)NrOffpulseBins;
      RMSQ = sqrt(RMSQ/(float)NrOffpulseBins);
      RMSU = sqrt(RMSU/(float)NrOffpulseBins);
      datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = sqrt(datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/(float)NrOffpulseBins);
      datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = sqrt(datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/(float)NrOffpulseBins);
      datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = sqrt(datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/(float)NrOffpulseBins);
      if(extended) {
	datafile->offpulse_rms[5+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] = sqrt(datafile->offpulse_rms[5+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/(float)NrOffpulseBins);
      }
      if(rms_file_specified) { // Scale all rms'es according to the specified rebinning
	// This assumes that the rebinning is done in such a way that the mean is preserved. If it would be simply summed, the rms would become bigger rather than smaller.
	float scale = 1.0/sqrt(rebin_factor);
	RMSQ *= scale;
	RMSU *= scale;
	datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] *= scale;
	datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] *= scale;
	datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] *= scale;
	if(extended) {
	  datafile->offpulse_rms[5+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)] *= scale;
	}
      }

      if(verbose.verbose) {
	if((freqnr == 0 && pulsenr == 0) || verbose.debug) {
	  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	  fprintf(stdout, "  PA conversion output for subint %ld frequency channel %ld:\n", pulsenr, freqnr);
	  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	  fprintf(stdout, "    Average baseline Stokes I: %f\n", baseline_intensity);
	  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	  fprintf(stdout, "    RMS I:                  %f\n", datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]);
	  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	  fprintf(stdout, "    RMS Q:                  %f\n", RMSQ);
	  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	  fprintf(stdout, "    RMS U:                  %f\n", RMSU);
	  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	  fprintf(stdout, "    RMS V:                  %f\n", datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]);
	  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	  fprintf(stdout, "    RMS L (before de-bias): %f\n", datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]);
	  if(extended) {
	    for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	    fprintf(stdout, "    RMS sqrt(Q^2+U^2+V^2):  %f\n", datafile->offpulse_rms[5+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]);
	  }
	}
      }
      /* NR code:
	 sort(NrOffpulseBins, &Loffpulse[-1]); 
	 medianL = Loffpulse[NrOffpulseBins/2];   
      */
      
      //      for(i = 0; i < NrOffpulseBins; i++) {
      //	printf("XXXX %ld (%d): %f\n", i, extended, Loffpulse[i]);
      //      }
      gsl_sort_float (Loffpulse, 1, NrOffpulseBins);
      medianL = gsl_stats_float_median_from_sorted_data(Loffpulse, 1, NrOffpulseBins);

      fflush(stdout);
      if((verbose.verbose && freqnr == 0 && pulsenr == 0) || verbose.debug) {
	for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	fprintf(stdout, "    Median L: %f\n", medianL);
      }
      if(extended) {
	/* NR code:
	   sort(NrOffpulseBins, &Poffpulse[-1]); 
	   medianP = Poffpulse[NrOffpulseBins/2];   
	*/
	//	for(i = 0; i < NrOffpulseBins; i++) {
	//	  printf("XXXX %ld (%d): %f\n", i, extended, Poffpulse[i]);
	//	}
	gsl_sort_float (Poffpulse, 1, NrOffpulseBins);
	medianP = gsl_stats_float_median_from_sorted_data(Poffpulse, 1, NrOffpulseBins);
	fflush(stdout);
	if((verbose.verbose && freqnr == 0 && pulsenr == 0) || verbose.debug) {
	  for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	  fprintf(stdout, "    Median sqrt(Q^2+U^2+V^2): %f\n", medianP);
	}
      }

      for(j = 0; j < (datafile->NrBins); j++) {
	if(correctLbias == 1) {
	  float junk = (0.5*(RMSQ+RMSU)/newdata[j+newindex_L]);
	  if(junk < 1)
	    newdata[j+newindex_L] *= sqrt(1.0-junk*junk);
	  else
	    newdata[j+newindex_L] = 0.0;
	}else if(correctLbias == 0) {
	  newdata[j+newindex_L] -= medianL;
	}
	if(extended)
	  newdata[j+newindex_T] -= medianP;
      }
      
      //      if(datafile->NrPApoints > 0) {
      //      if(datafile->NrPApoints != datafile->NrBins*datafile->NrSubints*datafile->NrFreqChan) {
      //	fflush(stdout);
      //	printerror(verbose.debug, "ERROR make_paswing_fromIQUV: NrPApoints != NrBins*NrSubints*NrFreqChan (%ld != %ld), no correct memory allocation done?", datafile->NrPApoints, datafile->NrBins*datafile->NrSubints*datafile->NrFreqChan);
      //	return 0;
      //      }
      for(i = 0; i < (datafile->NrBins); i++) {
	/* PA = 1/2 atan(U/Q) */
	/* PA = 90.0*atan2(U,Q)/M_PI; */
	newdata[i+newindex_Pa] = 90.0*atan2(datafile->data[sindex_U+i],datafile->data[sindex_Q+i])/M_PI;
	if(paoffset) {
	  newdata[i+newindex_Pa] += paoffset;
	  newdata[i+newindex_Pa] = derotate_180_small_double(newdata[i+newindex_Pa]);
	}

	//	  if(datafile->data_pa[i+datafile->NrBins*(freqnr+datafile->NrFreqChan*pulsenr)] > 90)
	//	    printf("XXXXX %f\n", datafile->data_pa[i+datafile->NrBins*(freqnr+datafile->NrFreqChan*pulsenr)]);
	if(datafile->data[i+sindex_Q] != 0 || datafile->data[i+sindex_U] != 0) {
	  /* dPA = sqrt(Q^2s_u^2+U^2s_q^2)/(2*(Q^2+U^2)). */
	  /* dPA = sqrt((Q*RMSU)*(Q*RMSU) + (U*RMSQ)*(U*RMSQ));*/
	  newdata[i+newindex_dPa] = sqrt((datafile->data[sindex_Q+i]*RMSU)*(datafile->data[sindex_Q+i]*RMSU) + (datafile->data[sindex_U+i]*RMSQ)*(datafile->data[sindex_U+i]*RMSQ));
	  /* dPA /= 2.0*(Q*Q + U*U); */
	  newdata[i+newindex_dPa] /= 2.0*(datafile->data[i+sindex_Q]*datafile->data[i+sindex_Q] + datafile->data[i+sindex_U]*datafile->data[i+sindex_U]); 
	  newdata[i+newindex_dPa] *= 180.0/M_PI;
	}else {  // Zapped data, set error to zero
	  newdata[i+newindex_dPa] = 0;
	}

	if(extended) {
//START REGION RELEASE
	  newdata[i+newindex_Ell] = 0;
	  //	  newdata[i+newindex_Ell] = 90.0*atan2(datafile->data[sindex_V+i], sqrt(datafile->data[sindex_Q+i]*datafile->data[sindex_Q+i]+datafile->data[sindex_U+i]*datafile->data[sindex_U+i]))/M_PI;
//START REGION DEVELOP
	  newdata[i+newindex_Ell] = 90.0*asin(datafile->data[sindex_V+i]/sqrt(datafile->data[sindex_Q+i]*datafile->data[sindex_Q+i]+datafile->data[sindex_U+i]*datafile->data[sindex_U+i]+datafile->data[sindex_V+i]*datafile->data[sindex_V+i]))/M_PI;
	  if(datafile->data[sindex_Q+i] == 0 && datafile->data[sindex_U+i] == 0 && datafile->data[sindex_V+i] == 0) {   // Zapped data, set error to zero
//START REGION RELEASE
	    newdata[i+newindex_dEll] = 0;
//START REGION DEVELOP
	  }else {
	    newdata[i+newindex_dEll]  = datafile->data[sindex_V+i]*datafile->data[sindex_V+i]*(datafile->data[sindex_Q+i]*datafile->data[sindex_Q+i]*RMSQ*RMSQ+datafile->data[sindex_U+i]*datafile->data[sindex_U+i]*RMSU*RMSU);
	    newdata[i+newindex_dEll] += (datafile->data[sindex_Q+i]*datafile->data[sindex_Q+i]+datafile->data[sindex_U+i]*datafile->data[sindex_U+i])*(datafile->data[sindex_Q+i]*datafile->data[sindex_Q+i]+datafile->data[sindex_U+i]*datafile->data[sindex_U+i])*datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]*datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)];
	    newdata[i+newindex_dEll] /= 4.0*(datafile->data[sindex_Q+i]*datafile->data[sindex_Q+i]+datafile->data[sindex_U+i]*datafile->data[sindex_U+i]);
	    newdata[i+newindex_dEll] = sqrt(newdata[i+newindex_dEll]);
	    newdata[i+newindex_dEll] /= (datafile->data[sindex_Q+i]*datafile->data[sindex_Q+i]+datafile->data[sindex_U+i]*datafile->data[sindex_U+i]+datafile->data[sindex_V+i]*datafile->data[sindex_V+i]);
	    newdata[i+newindex_dEll] *= 180.0/M_PI;
	  }
//START REGION RELEASE
	}

      }
	//      }else {
	//	fflush(stdout);
	//	printerror(verbose.debug, "ERROR make_paswing_fromIQUV: PA array not defined, should allocate memory first.");
	//	return 0;
	//      }

      // Set I, V as well in new dataset
      for(j = 0; j < (datafile->NrBins); j++) {
	newdata[j+newindex_I] = datafile->data[j+sindex_I];
	newdata[j+newindex_V] = datafile->data[j+sindex_V];
      }

//START REGION DEVELOP
      double sumI = 0;
      long nrIbins = 0;
      double sumL = 0;
      long nrLbins = 0;
      double sumV = 0;
      for(i = 0; i < (datafile->NrBins); i++) {
	if(checkRegions(i, &onpulse, 0, verbose) == 1) {
	  sumI += newdata[i+newindex_I];
	  sumV += newdata[i+newindex_V];
	  //	  printf("XXXX %f\t%lf\n", newdata[i+newindex_V], sumV);
	  nrIbins++;
	  if(newdata[i+newindex_L] != 0.0) {
	    sumL += newdata[i+newindex_L];
	    nrLbins++;
	  }
	}
      }
      if(verbose.debug) {
	for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	fprintf(stdout, "    Considering first selected onpulse region only\n");
	for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	fprintf(stdout, "    Sum of I (%ld onpulse bins):              %f +- %f\n", nrIbins, sumI, sqrt(nrIbins)*datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]);
	for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	fprintf(stdout, "    Sum of L (%ld onpulse bins where L != 0): %f +- %f\n", nrLbins, sumL, sqrt(nrIbins)*datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]);
	double LoverI = sumL/sumI;
	// (L/I)*sqrt(((RMS_L/L)^2+(RMS_I/I)^2)*nrbins)
	double LoverIerr = LoverI*sqrt(((datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/sumL)*(datafile->offpulse_rms[1+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/sumL) + (datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/sumI)*(datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/sumI))*nrIbins);
	for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	fprintf(stdout, "    %%L: %f +- %f (positively biased and error propagating only valid for high significance points)\n", 100.0*LoverI, 100.0*LoverIerr);
	for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	fprintf(stdout, "    Sum of V (%ld onpulse bins):              %f +- %f\n", nrIbins, sumV, sqrt(nrIbins)*datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]);
	double VoverI = sumV/sumI;
	// (V/I)*sqrt(((RMS_V/V)^2+(RMS_I/I)^2)*nrbins)
	//	printf("XXXX %f\n", datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]);
	//	printf("XXXX %f\n", datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]);
	double VoverIerr = fabs(VoverI)*sqrt(((datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/sumV)*(datafile->offpulse_rms[2+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/sumV) + (datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/sumI)*(datafile->offpulse_rms[0+output_nr_pols*(freqnr + datafile->NrFreqChan*pulsenr)]/sumI))*nrIbins);
	for(indent = 0; indent < verbose.indent; indent++) printf(" ");
	fprintf(stdout, "    %%V: %f +- %f\n", 100.0*VoverI, 100.0*VoverIerr);
      }
//START REGION RELEASE

    } // End of frequency channel loop
  }   // End of subint loop

  // Free Stokes data and replace it with ILVPAdPA data
  free(datafile->data);
  datafile->data = newdata;
  if(rms_file_specified) {
    free(newdata_rms);
  }
  if(nolongitudes == 0) {
    datafile->tsampMode = TSAMPMODE_LONGITUDELIST;
  }
  datafile->NrPols = output_nr_pols;
  if(extended) {
    datafile->poltype = POLTYPE_ILVPAdPATEldEl;
  }else {
    datafile->poltype = POLTYPE_ILVPAdPA;
  }

  free(Loffpulse);
  free(Poffpulse);
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

// Returns 0 on error, 1 if ok.
int writePPOLHeader(datafile_definition datafile, int argc, char **argv, verbose_definition verbose)
{
  char *txt;
  txt = malloc(10000);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePPOLHeader: Memory allocation error.");
    return 0;
  }
  constructCommandLineString(txt, 10000, argc, argv, verbose);
  fprintf(datafile.fptr_hdr, "#ppol file: %s\n", txt);
  free(txt);
  return 1;
}


//START REGION DEVELOP
//START REGION RELEASE

int readPPOLHeader(datafile_definition *datafile, int extended, verbose_definition verbose)
{
  float dummy_float;
  int ret, maxlinelength, nrwords;
  char *txt, *ret_ptr, *word_ptr;

  maxlinelength = 2000;
  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLHeader: Memory allocation error.");
    return 0;
  }

  // Assume data is folded, some information not known from the header
  datafile->isFolded = 1;
  datafile->foldMode = FOLDMODE_FIXEDPERIOD;
  datafile->fixedPeriod = 0;
  datafile->tsampMode = TSAMPMODE_LONGITUDELIST;
  datafile->fixedtsamp = 0;
  datafile->tsubMode = TSUBMODE_FIXEDTSUB;
  if(datafile->tsub_list != NULL)
    free(datafile->tsub_list);
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLHeader: Memory allocation error");
    return 0;
  }
  datafile->tsub_list[0] = 0;

  datafile->NrSubints = 1;
  datafile->NrFreqChan = 1;
  datafile->datastart = 0;
  //  datafile->longitudes_defined = 1;

  rewind(datafile->fptr);
  ret = fread(txt, 1, 3, datafile->fptr);
  txt[3] = 0;
  if(ret != 3) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLHeader: cannot read from file.");
    free(txt);
    return 0;
  }
  if(strcmp(txt, "#pp") != 0
//START REGION DEVELOP
&& strcmp(txt, "#ma") != 0
//START REGION RELEASE
) {
    fflush(stdout);
    printwarning(verbose.debug, "WARNING readPPOLHeader: File does not appear to be in PPOL or PPOLSHORT format. I will try to load file, but this will probably fail. Did you run ppol first?");
  }

  skipallhashedlines(datafile);
  datafile->NrBins = 0;
  dummy_float = 0;
  do {
    ret_ptr = fgets(txt, maxlinelength, datafile->fptr);
    if(ret_ptr != NULL) {
      if(txt[0] != '#') {
	if(extended) {
	  word_ptr = pickWordFromString(txt, 2, &nrwords, 1, ' ', verbose);
	  if(nrwords != 10 && nrwords != 14) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR readPPOLHeader: Line should have 10 or 14 words, got %d", nrwords);
	    if(nrwords == 3)
	      printerror(verbose.debug, "                             Maybe file is in format %s?", returnFileFormat_str(PPOL_SHORT_format));
	    printerror(verbose.debug, "                             Line: '%s'.", txt);
	    free(txt);
	    return 0;
	  }
	  if(nrwords == 10) {
	    datafile->poltype = POLTYPE_ILVPAdPA;
	    datafile->NrPols = 5;  /* I,L,V,Pa,dPa */
	  }else {
	    datafile->poltype = POLTYPE_ILVPAdPATEldEl;
	    datafile->NrPols = 8;  /* I,L,V,Pa,dPa, totpol, Ell, dEll */
	  }
	}else {
	  word_ptr = pickWordFromString(txt, 1, &nrwords, 1, ' ', verbose);
	  if(nrwords != 3) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR readPPOLHeader: Line should have 3 words, got %d", nrwords);
	    if(nrwords == 10)
	      printerror(verbose.debug, "                             Maybe file is in format %s?",  returnFileFormat_str(PPOL_format));
	    printerror(verbose.debug, "                             Line: '%s'.", txt);
	    free(txt);
	    return 0;
	  }
	}
	ret = sscanf(word_ptr, "%f", &dummy_float);
	if(ret != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readPPOLHeader: Cannot interpret as a float: '%s'.", txt);
	  free(txt);
	  return 0;
	}
	if(dummy_float >= 360) {
	  fflush(stdout);
	  printwarning(verbose.debug, "WARNING: IGNORING POINTS AT PULSE LONGITUDES > 360 deg.");
	}else {
	  (datafile->NrBins)++;
	}
      }
    }
  }while(ret_ptr != NULL && dummy_float < 360);


  if(extended == 0) {
    datafile->poltype = POLTYPE_PAdPA;
    datafile->NrPols = 2;  /* Just PA,dPa */
  }


  //  datafile->NrPApoints = datafile->NrBins;
  fflush(stdout);
  if(verbose.verbose) fprintf(stdout, "Going to load %ld points from %s\n", datafile->NrBins, datafile->filename); 
  if(datafile->NrBins == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLHeader: No data in %s", datafile->filename); 
    free(txt);
    return 0;
  }
  fseek(datafile->fptr, datafile->datastart, SEEK_SET);  
  free(txt);

  // Allocate rms list & pulse longitude list
  if(datafile->offpulse_rms != NULL) {
    free(datafile->offpulse_rms);
    datafile->offpulse_rms = NULL;
  }
  if(extended) {
    datafile->offpulse_rms = (float *)malloc(datafile->NrSubints*datafile->NrFreqChan*datafile->NrPols*sizeof(float));
    if(datafile->offpulse_rms == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPPOLHeader: Memory allocation error"); 
      return 0;
    }
  }
  datafile->tsamp_list = (double *)malloc(datafile->NrBins*sizeof(double));
  if(datafile->tsamp_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLHeader: Memory allocation error"); 
    return 0;
  }

  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/*
It will ignore all points after 360 degrees longitude. This is because
it is important for the pa-fitting that there is only one profile
loaded, and not two rotation periods. 
*/
int readPPOLfile(datafile_definition *datafile, float *data, int extended, float add_longitude_shift, verbose_definition verbose)
{
  int maxlinelength;
  long i, k, dummy_long;
  char *txt, *ret_ptr;

  if(datafile->NrBins == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLfile: No data in %s", datafile->filename); 
    return 0;
  }
  maxlinelength = 2000;
  txt = malloc(maxlinelength);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPPOLfile: Memory allocation error.");
    return 0;
  }

  fseek(datafile->fptr, datafile->datastart, SEEK_SET);  
  k = 0;
  if(extended) {
    datafile->offpulse_rms[3] = -1;  // RMS on PA and dPA polarization are undefined
    datafile->offpulse_rms[4] = -1;
    if(datafile->NrPols == 8) {
      datafile->offpulse_rms[6] = -1;
      datafile->offpulse_rms[7] = -1;
    }
  }
  for(i = 0; i < datafile->NrBins; i++) {
    ret_ptr = fgets(txt, maxlinelength, datafile->fptr);
    if(ret_ptr == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPPOLfile: Cannot read next line, should not happen after successfully reading in header");
      free(txt);
      return 0;
    }
    if(txt[0] != '#') {
      if(extended == 0) {
	sscanf(txt, "%lf %f %f", &(datafile->tsamp_list[k]), &(data[k]), &(data[k+datafile->NrBins]));
      }else {
	if(datafile->NrPols == 8) {
	  sscanf(txt, "%ld %lf %f %f %f %f %f %f %f %f %f %f %f %f", &dummy_long, &(datafile->tsamp_list[k]), &(data[k]), &(datafile->offpulse_rms[0]), &(data[k+datafile->NrBins]), &(datafile->offpulse_rms[1]), &(data[k+2*datafile->NrBins]), &(datafile->offpulse_rms[2]), &(data[k+3*datafile->NrBins]), &(data[k+4*datafile->NrBins]), &(data[k+5*datafile->NrBins]), &(datafile->offpulse_rms[5]), &(data[k+6*datafile->NrBins]), &(data[k+7*datafile->NrBins]));
	}else {
	  sscanf(txt, "%ld %lf %f %f %f %f %f %f %f %f", &dummy_long, &(datafile->tsamp_list[k]), &(data[k]), &(datafile->offpulse_rms[0]), &(data[k+datafile->NrBins]), &(datafile->offpulse_rms[1]), &(data[k+2*datafile->NrBins]), &(datafile->offpulse_rms[2]), &(data[k+3*datafile->NrBins]), &(data[k+4*datafile->NrBins]));
	}
      }
      datafile->tsamp_list[k] += add_longitude_shift;
      if(datafile->tsamp_list[k] >= 0 && datafile->tsamp_list[k] < 360) {
	//	if(data_dpa[k] > 0 || onlysignificantPA == 0)
	  k++;
      }else {
	fflush(stdout);
	printwarning(verbose.debug, "WARNING readPPOLfile: IGNORING POINTS AT PULSE LONGITUDES > 360 deg.");
      }
    }
  }
  if(k != datafile->NrBins) {
    fflush(stdout);
    printerror(verbose.debug, "WARNING readPPOLfile: The nr of bins read in is different as determined from header. Something is wrong.");
    return 0;
  }
  //  datafile->NrBins = k;
  fflush(stdout);
  if(verbose.verbose) fprintf(stdout, "readPPOLfile: Accepted %ld points\n", datafile->NrBins); 
  /*
  if(verbose.verbose) {
    for(i = 0; i < datafile->NrBins; i++) 
      fprintf(stderr, "%f %f %f\n", data_long[i], data_pa[i], data_dpa[i]);
      }*/
  /*  fclose(datafile->fptr);  */
	    free(txt);
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/*
Set fptr in datafile to to stdout is possible. If onlysignificantPA is
set only the significant PA's are written out. Can decide to write out
two profiles and or only the significant PA values. If extended == 0,
only the PA and its errorbar are written out. Return value of 1 is OK.
*/
int writePPOLfile(datafile_definition datafile, float *data, int extended, int onlysignificantPA, int twoprofiles, float PAoffset, verbose_definition verbose)
{
  long j;

  if(datafile.poltype != POLTYPE_ILVPAdPA && datafile.poltype != POLTYPE_PAdPA && datafile.poltype != POLTYPE_ILVPAdPATEldEl) {
    printerror(verbose.debug, "ERROR writePPOLfile: Data doesn't appear to have poltype ILVPAdPA, PAdPA or ILVPAdPATEldEl (it is %d).", datafile.poltype); 
    return 0;
  }
  if(datafile.poltype == POLTYPE_ILVPAdPA && datafile.NrPols != 5) {
    printerror(verbose.debug, "ERROR writePPOLfile: 5 polarization channels were expected, but there are %ld.", datafile.NrPols); 
    return 0;
  }else if(datafile.poltype == POLTYPE_PAdPA && datafile.NrPols != 2) {
    printerror(verbose.debug, "ERROR writePPOLfile: 2 polarization channels were expected, but there are %ld.", datafile.NrPols); 
    return 0;
  }else if(datafile.poltype == POLTYPE_ILVPAdPATEldEl && datafile.NrPols != 8) {
    printerror(verbose.debug, "ERROR writePPOLfile: 8 polarization channels were expected, but there are %ld.", datafile.NrPols); 
    return 0;
  }
  if(datafile.NrSubints > 1 || datafile.NrFreqChan > 1) { // As otherwise different bins are removed from different subints.
    printerror(verbose.debug, "ERROR writePPOLfile: Can only do this opperation if there is one subint and one frequency channel."); 
    return 0;
  }
  if(datafile.tsampMode != TSAMPMODE_LONGITUDELIST) {
    printerror(verbose.debug, "ERROR writePPOLfile: Expected pulse longitudes to be defined."); 
    return 0;
  }

  int pa_offset, dpa_offset;
  if(datafile.poltype == POLTYPE_ILVPAdPA || datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
    pa_offset = 3;
    dpa_offset = 4;
  }else if(datafile.poltype == POLTYPE_PAdPA) {
    pa_offset = 0;
    dpa_offset = 1;
  }

  for(j = 0; j < datafile.NrBins; j++) {
    if(data[j+dpa_offset*datafile.NrBins] > 0 || onlysignificantPA == 0) {
      if(extended) {
	/*
	fprintf(stderr, "XXXXX j=%ld NrBins=%ld\n", j, datafile.NrBins);
	fprintf(datafile.fptr, "%ld\n", j);
	fprintf(datafile.fptr, "%e\n", datafile.tsamp_list[j]);
	fprintf(datafile.fptr, "%e\n", data[j]);
	fprintf(datafile.fptr, "%e\n", datafile.offpulse_rms[0]);
	fprintf(datafile.fptr, "%e\n", data[j+datafile.NrBins]);
	fprintf(datafile.fptr, "%e\n", datafile.offpulse_rms[1]);
	fprintf(datafile.fptr, "%e\n", data[j+2*datafile.NrBins]);
	fprintf(datafile.fptr, "%e\n", datafile.offpulse_rms[2]);
	fprintf(datafile.fptr, "%e\n", data[j+3*datafile.NrBins]+PAoffset);
	fprintf(datafile.fptr, "%e\n", data[j+dpa_offset*datafile.NrBins]);
	*/
	
	fprintf(datafile.fptr, "%ld %e %e %e %e %e %e %e %e %e", j, datafile.tsamp_list[j], data[j], datafile.offpulse_rms[0], data[j+datafile.NrBins], datafile.offpulse_rms[1], data[j+2*datafile.NrBins], datafile.offpulse_rms[2], data[j+pa_offset*datafile.NrBins]+PAoffset, data[j+dpa_offset*datafile.NrBins]);
	if(datafile.poltype == POLTYPE_ILVPAdPATEldEl)
	  fprintf(datafile.fptr, " %e %e %e %e", data[j+5*datafile.NrBins], datafile.offpulse_rms[5], data[j+6*datafile.NrBins], data[j+7*datafile.NrBins]);
	fprintf(datafile.fptr, "\n");
      }else {
	fprintf(datafile.fptr, "%e %e %e\n", datafile.tsamp_list[j], data[j+pa_offset*datafile.NrBins]+PAoffset, data[j+dpa_offset*datafile.NrBins]);
      }
    }
  }
  if(twoprofiles) {
    for(j = 0; j < datafile.NrBins; j++) {
      if(data[j+dpa_offset*datafile.NrBins] > 0 || onlysignificantPA == 0) {
	if(extended) {
	  fprintf(datafile.fptr, "%ld %e %e %e %e %e %e %e %e %e", j, datafile.tsamp_list[j]+360, data[j], datafile.offpulse_rms[0], data[j+datafile.NrBins], datafile.offpulse_rms[1], data[j+2*datafile.NrBins], datafile.offpulse_rms[2], data[j+pa_offset*datafile.NrBins]+PAoffset, data[j+dpa_offset*datafile.NrBins]);
	  if(datafile.poltype == POLTYPE_ILVPAdPATEldEl)
	    fprintf(datafile.fptr, " %e %e %e %e", data[j+5*datafile.NrBins], datafile.offpulse_rms[5], data[j+6*datafile.NrBins], data[j+7*datafile.NrBins]);
	  fprintf(datafile.fptr, "\n");
	}else {
	  fprintf(datafile.fptr, "%e %e %e\n", datafile.tsamp_list[j]+360, data[j]+PAoffset, data[j+dpa_offset*datafile.NrBins]);
	}
      }
    }
  }
  return 1;
}

//START REGION DEVELOP

/*
  Formula from Everett & Weisberg 2001, although they didn't invent the formula.
  
  P0 = L/sigmaI
  eta0 = (P0/sqrt(2))*cos(2*dPA)
  G(dPA) = (1/sqrt(M_PI)+eta0*exp(eta0*eta0)*(1+erff(eta0)))/sqrt(M_PI)
*/
double internal_pa_probability_distribution_P0;

double pa_probability_distribution(double dpa)
{
  double P0, eta0, G;
  P0 = internal_pa_probability_distribution_P0;
  eta0 = (P0/sqrt(2.0))*cos(2.0*dpa);
  G = (1.0/sqrt(M_PI)+eta0*exp(eta0*eta0)*(1.0+erff(eta0)))*exp(-0.5*P0*P0)/sqrt(M_PI);
  return G;
}

double pa_probability_distribution_gsl(double dpa, void *params)
{
  return pa_probability_distribution(dpa);
}

/* P0 = L/sigma (the std. deviation of any of the Stokes parameters,
   they should be equal). L is sqrt(Q^2+U^2) for the pa point, after
   subtracting the bias in L. The 1 sigma errorbar of the pa-point is
   returned in radians. Formula is from Naghizadeh-Khouei & Clarke
   1993.*/
double realisticPAerrorbar(double P0, verbose_definition verbose)
{
  double dpa, ddpa, I, abserr;
  int direction;
  size_t neval;
  gsl_function F;

  internal_pa_probability_distribution_P0 = P0;
  dpa = 0.01;
  ddpa = 0.01;
  direction = 1;
  F.function = &pa_probability_distribution_gsl;
  F.params = NULL;
  do {
    /*    I = qromb(pa_probability_distribution, -dpa, dpa); */
    /* NR code:    I = qsimp(pa_probability_distribution, -dpa, dpa); */
    if(gsl_integration_qng (&F, -dpa, dpa, 0, 1e-6, &I, &abserr, &neval) != 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR realisticPAerrorbar: gsl_integration_qng Failed");
      exit(0);
    }
    /*   fprintf(stderr, "XXXX I=%f dpa=%f ddpa=%e factor %f higher\n", I, dpa, ddpa, ddpa/(0.00001*M_PI/180.0)); */
    if(I < 0.682689492137086 && direction == 1) {  /* Keep increasing dpa */
      dpa += ddpa;
    }else if(I < 0.682689492137086 && direction == -1) {  /* Reverse direction ddpa */
      ddpa *= 0.5;
      direction *= -1;
    }else if(I > 0.682689492137086 && direction == 1) {   /* Reverse direction ddpa */
      ddpa *= 0.5;
      direction *= -1;
    }else if(I > 0.682689492137086 && direction == -1) {   /* Continue decreasing ddpa */
      dpa -= ddpa;
    }else {                                           /* Must be very close to answer if gets here */
      ddpa *= 0.5;
    }
  }while(ddpa > 0.00001*M_PI/180.0);
  /*
  I = qsimp(pa_probability_distribution, -M_PI*0.5, M_PI*0.5);
  fprintf(stderr, "XXXX Itot=%f\n", I);
  */
  return dpa;
}

/*
  Makes a PA (or ellipticity) distribution with nrbins pa bins. If
  ellipticity != 0, the ellipticities are used rather than the pa. The
  input data should already be in the form of I, L, V, PA and its
  error bar. The memory for the output data will be allocated and
  dataout can be completely uninitialized. If normalise is set, the
  distribution will be normalised such that the pulse longitude bin
  column with most counts will have a sum of 1.

  In addition, a "PA-mask" can be supplied. This is an array of data,
  with the dimensions as the pa-distribution to be produced, i.e. the
  number of "subint" should be equal to the nrbins variable and the
  number of pulse longitude bins should match those of datain. If this
  pointer to pamask is set to NULL, it will be ignored. Otherwise, the
  input data will be masked in the following way:

  For each PA bin (i.e. subint)/pulse-longitude bin combination,
  pamask will have a recorded value in the first polarization
  channel. If the PA value of a given sample in the input data
  (datain, the pulse stack) matches the value given by pamask_value
  (+- 0.01), the input data will not be altered. Otherwise, that
  particular sample is set to zero. If pamask_value is NaN, the output
  is either the corresponding value as specified in the mask file, or
  zero when the PA is not significant.

  Returns 0 on error
 */
int make_pa_distribution(datafile_definition datain, datafile_definition *dataout, int nrbins, int normalise, datafile_definition *pamask, float pamask_value, int ellipticity, verbose_definition verbose)
{
  long i, j, f, nrpointsadded, nrpointsadded_max, binnr;
  float dpa;

  if(datain.NrSubints <= 1 && datain.NrFreqChan <= 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_pa_distribution: Need more than a single subints and frequency channel to make a PA distribution");
    return 0;
  }

  if(datain.poltype != POLTYPE_ILVPAdPA && datain.poltype != POLTYPE_PAdPA && datain.poltype != POLTYPE_ILVPAdPATEldEl) {
    printerror(verbose.debug, "ERROR make_pa_distribution: Data doesn't appear to have poltype ILVPAdPA, ILVPAdPATEldEl or PAdPA."); 
    return 0;
  }
  if(ellipticity && datain.poltype != POLTYPE_ILVPAdPATEldEl) {
    printerror(verbose.debug, "ERROR make_pa_distribution: The data isn't of polarization type ILVPAdPATEldEl, while the ellipticity distribution was requested."); 
    return 0;
  }

  if(datain.poltype == POLTYPE_ILVPAdPA && datain.NrPols != 5) {
    printerror(verbose.debug, "ERROR make_pa_distribution: 5 polarization channels were expected, but there are %ld.", datain.NrPols); 
    return 0;
  }else if(datain.poltype == POLTYPE_ILVPAdPATEldEl && datain.NrPols != 8) {
    printerror(verbose.debug, "ERROR make_pa_distribution: 8 polarization channels were expected, but there are %ld.", datain.NrPols); 
    return 0;
  }else if(datain.poltype == POLTYPE_PAdPA && datain.NrPols != 2) {
    printerror(verbose.debug, "ERROR make_pa_distribution: 2 polarization channels were expected, but there are %ld.", datain.NrPols); 
    return 0;
  }

  if(pamask != NULL) {
    if(datain.NrBins != pamask->NrBins) {
      printerror(verbose.debug, "ERROR make_pa_distribution: Applying a PA mask only works when the input data has the same number of pulse longitude bins compared to that of the provided mask. (the input data has %ld pulse longitude bins, while the mask has %ld).", datain.NrBins, pamask->NrBins); 
      return 0;    
    }
    if(nrbins != pamask->NrSubints) {
      printerror(verbose.debug, "ERROR make_pa_distribution: Applying a PA mask only works when generating a PA-distribution with an equal number of PA bins compared to that of the provided mask. (now %ld pa bins are requested, while the mask has %ld pa-bins defined).", nrbins, pamask->NrSubints); 
      return 0;    
    }
    if(pamask->NrFreqChan > 1) {
      printerror(verbose.debug, "ERROR make_pa_distribution: Applying a PA mask only works when the mask has one frequency channel defined. There are currently %ld channels defined).", pamask->NrFreqChan); 
      return 0;    
    }
  }

  cleanPSRData(dataout, verbose);
  copy_params_PSRData(datain, dataout, verbose);
  dataout->format = MEMORY_format;
  dataout->NrSubints = nrbins;
  dataout->NrPols = 1;
  //  dataout->NrPApoints = 0;
  dataout->NrFreqChan = 1;
  //  dataout->longitudes_defined = 0;
  //  dataout->bins_defined = 0;
  if(ellipticity == 0)
    dataout->gentype = GENTYPE_PADIST;
  else
    dataout->gentype = GENTYPE_ELLDIST;
  dataout->tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
  if(dataout->tsub_list != NULL)
    free(dataout->tsub_list);
  dataout->tsub_list = (double *)malloc(sizeof(double));
  if(dataout->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_pa_distribution: Memory allocation error");
    return 0;
  }
  dataout->tsub_list[0] = get_tobs(datain, verbose);

  dataout->yrangeset = 1;
  if(ellipticity == 0) {
    dataout->yrange[0] = -90.0+0.5*180.0/(float)(nrbins);
    dataout->yrange[1] = 90.0-0.5*180.0/(float)(nrbins);
  }else {
    dataout->yrange[0] = -45.0+0.5*90.0/(float)(nrbins);
    dataout->yrange[1] = 45.0-0.5*90.0/(float)(nrbins);
  }

  dataout->data = (float *)calloc(dataout->NrSubints*dataout->NrBins*dataout->NrPols*dataout->NrFreqChan, sizeof(float));
  if(dataout->data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_pa_distribution: Cannot allocate memory for data.");
    return 0;
  }

  if(ellipticity == 0) {
    dpa = 180.0/(float)nrbins;
  }else {
    dpa = 90.0/(float)nrbins;
  }


  int pa_chan, dpa_chan;
  if(datain.poltype == POLTYPE_ILVPAdPA || datain.poltype == POLTYPE_ILVPAdPATEldEl) {
    if(ellipticity == 0) {
      pa_chan = 3;
      dpa_chan = 4;
    }else {
      pa_chan = 6;
      dpa_chan = 7;
    }
  }else {
    pa_chan = 0;
    dpa_chan = 1;
  }

  nrpointsadded_max = 0;
  for(j = 0; j < datain.NrBins; j++) {
    nrpointsadded = 0;
    for(i = 0; i < datain.NrSubints; i++) {
      for(f = 0; f < datain.NrFreqChan; f++) {
	float paerr;
	paerr = datain.data[j+datain.NrBins*(dpa_chan+datain.NrPols*(f+datain.NrFreqChan*i))];
	if(paerr > 0) {
	  float pa = datain.data[j+datain.NrBins*(pa_chan+datain.NrPols*(f+datain.NrFreqChan*i))];
	  if(ellipticity == 0) {
	    if(pa == 90.0)  // In the unlikely event the pa is exactly 90, put it in first bin to avoid an overflow of the array
	      binnr = 0;
	    else
	      binnr = (pa + 90.0)/dpa;
	  }else {
	    if(pa == 45.0)  // In the unlikely event the pa is exactly 45, put it in first bin to avoid an overflow of the array
	      binnr = 0;
	    else
	      binnr = (pa + 45.0)/dpa;
	  }
	  //	  fprintf(stderr, "XXXXXX %f %ld\n", pa, binnr);
	  if(binnr < 0 || binnr >= nrbins) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR make_pa_distribution: %ld %f %f BUG!!!!!!!!!!!!!!!", binnr, datain.data[j+datain.NrBins*(pa_chan+datain.NrPols*(f+datain.NrFreqChan*i))], dpa);
	    return 0;
	  }
	  dataout->data[j+dataout->NrBins*binnr] += 1.0;
	  nrpointsadded++;
	}
	if(pamask != NULL) {
	  float value;
	  if(paerr <= 0) {
	    value = 0; // The value is irrelevant in condition below, but if paerr < 0, binnr is not necessarily initialised above.
	  }else {
	    value = pamask->data[j+pamask->NrBins*(0+pamask->NrPols*(0+pamask->NrFreqChan*binnr))];
	  }
	  if(isnan(pamask_value)) {
	    int curpol;
	    for(curpol = 0; curpol < datain.NrPols; curpol++) {
	      datain.data[j+datain.NrBins*(curpol+datain.NrPols*(f+datain.NrFreqChan*i))] = value;
	    }
	  }else {
	    if(value < pamask_value-0.01 || value > pamask_value+0.01 || paerr <= 0) {
	      int curpol;
	      for(curpol = 0; curpol < datain.NrPols; curpol++) {
		datain.data[j+datain.NrBins*(curpol+datain.NrPols*(f+datain.NrFreqChan*i))] = 0;
	      }
	    }
	  }
	}
      }
    }
    if(nrpointsadded > nrpointsadded_max)
      nrpointsadded_max = nrpointsadded; 
  }
  // Normalise the histogram, if there are data points.
  if(normalise) {
    if(nrpointsadded_max > 0) {
      for(i = 0; i < nrbins; i++) {
	for(j = 0; j < datain.NrBins; j++) {
	  dataout->data[j+datain.NrBins*i] /= (float)nrpointsadded_max;
	}
      }
    }
  }
  return 1;
}

//START REGION DEVELOP

/*
  Fit a PA distribution with one or more von Mises functions. The
  memory will be allocated and fit (containing the fit) and resid
  (containing the residual) can be completely uninitialized.

  maxNrComp = The maximum allowed number of functions to fit.
  sigmalimit = the limit up to where components are fitted (5 seems to work reasonably good). 
  precision = The larger precision is set to the finer trials are tried. Something like 2 seems reasonable.

  Returns 0 on error
 */
int fit_pa_distribution(datafile_definition datain, datafile_definition *fit, datafile_definition *resid, int maxNrComp, float sigmalimit, int precision, verbose_definition verbose)
{
  long i, j;
  datafile_definition padist;
  // Define as pointer to avoid to allow dynamic memory allocation
  vonMises_collection_definition *vonMises;
  float baseline;
  verbose_definition noverbose;

  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;
  if(verbose.debug) {
    noverbose.debug = 1;
    noverbose.verbose = 1;
  }

  if(datain.gentype != GENTYPE_PADIST) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fit_pa_distribution: Expected the input of this function to to have a gentype set to PA distribution.");
    return 0;
  }
  if(datain.NrPols != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fit_pa_distribution: Expected the input of this function to have a single polarization channel");
    return 0;
  }
  if(datain.NrFreqChan != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fit_pa_distribution: Expected the input of this function to have a single frequency channel");
    return 0;
  }

  cleanPSRData(fit, verbose);
  cleanPSRData(resid, verbose);
  copy_params_PSRData(datain, fit, verbose);
  copy_params_PSRData(datain, resid, verbose);
  fit->format = MEMORY_format;
  fit->tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
  resid->format = MEMORY_format;
  resid->tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
  if(fit->tsub_list != NULL)
    free(fit->tsub_list);
  if(resid->tsub_list != NULL)
    free(resid->tsub_list);
  fit->tsub_list = (double *)malloc(sizeof(double));
  resid->tsub_list = (double *)malloc(sizeof(double));
  vonMises = (vonMises_collection_definition *)malloc(sizeof(vonMises_collection_definition));
  if(fit->tsub_list == NULL || resid->tsub_list == NULL || vonMises == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fit_pa_distribution: Memory allocation error");
    return 0;
  }
  fit->tsub_list[0] = get_tobs(datain, verbose);
  resid->tsub_list[0] = get_tobs(datain, verbose);

  fit->data = (float *)calloc(fit->NrSubints*fit->NrBins*fit->NrPols*fit->NrFreqChan, sizeof(float));
  resid->data = (float *)calloc(fit->NrSubints*fit->NrBins*fit->NrPols*fit->NrFreqChan, sizeof(float));
  if(fit->data == NULL || resid->data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fit_pa_distribution: Cannot allocate memory for data.");
    return 0;
  }

  cleanPSRData(&padist, verbose);
  copy_params_PSRData(datain, &padist, verbose);
  padist.NrSubints = 1;
  padist.NrBins = datain.NrSubints;

  padist.data =   (float *)malloc(datain.NrSubints*sizeof(float));
  if(padist.data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR fit_pa_distribution: Cannot allocate memory for data.");
    return 0;
  }

  for(j = 0; j < datain.NrBins; j++) {
    float sum;
    sum = 0;
    for(i = 0; i < datain.NrSubints; i++) {
      // Extract a column out of the pa distribution for the current pulse longitude
      padist.data[i] = datain.data[j+datain.NrBins*i];
      sum += padist.data[i];
      //      printf("XXXXX %f\n", padist.data[i]);
    }
    // If sum == 0, then there is no data and fit will fail, so no point even in trying.
    if(sum == 0) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING fit_pa_distribution: Ignoring empty bin %d.", j);
      for(i = 0; i < datain.NrSubints; i++) {
	resid->data[j+resid->NrBins*i] = datain.data[j+datain.NrBins*i];
      }
    }else {
      //      int maxNrComp = 10; // Allow at most 2 components to be defined.
      //      float sigmalimit = 3.0; // (5 seems to work reasonably good
      int educatedguess = 1;
      //      int precision = 3;
      int ret;
      ret = fitvonmises(padist, vonMises, maxNrComp, sigmalimit, educatedguess, 1, precision, 0, &baseline, NULL, noverbose, 0, NULL);
      if(ret != 1) {
	fflush(stdout);
	printwarning(verbose.debug, "WARNING fit_pa_distribution: Fitting PA distribution failed for bin %d.", j);
	for(i = 0; i < datain.NrSubints; i++) {
	  resid->data[j+resid->NrBins*i] = datain.data[j+datain.NrBins*i];
	}
      }else {
	/* Calculate the profile out of the model by applying a phase
   shift. If the normalize flag is set the profile will be
   normalized. */
	calcVonMisesProfile(vonMises, datain.NrSubints, padist.data, 0, 0);
	for(i = 0; i < datain.NrSubints; i++) {
	  resid->data[j+resid->NrBins*i] = datain.data[j+datain.NrBins*i]-padist.data[i]-baseline;
	  fit->data[j+fit->NrBins*i] = padist.data[i]+baseline;
	}
      }
    }

  }

  closePSRData(&padist, 0, verbose);
  return 1;
}


/* 
   Makes a projection of the the Poincare sphere on a
   map of nrx by nry points (idealy with a ration 2:1). The area
   outside the spherical surface will be filled with value
   background. 

   projection = 1 Hammer-Aitoff projection. 
   projection = 2 3D sphere
   projection = 3 Simple linear longitude-latitude map

   For projection 1 and 2 the horizontal axis runs from -2.25 to 2.25
   and the vertical axis from -1.125 to 1.125. For projection 1 the
   poles are at (0, +1) and (0, -1) and the equator extends from (-2,
   0) to (+2, 0). For projection 2 the sphere has a radius of 1. If a
   point lies on the front-side, the sphere is centred at x=-1.125,
   y=0. If the point lies on the far side the point ends on a sphere
   with an x-offset in the opposite direction. For projection 3 the
   horizontal axis runs from -180 to 180 degrees and the vertical axis
   from -90 to 90 degrees.


   The baseline of the data is expected to be subtracted. If weighting
   = 0, the map is basically a count map. If set to 1, weighting
   w.r.t. the polarized power is applied.

   If binnr is non-negative, only the selected pulse longitude bin is
   considered. If binnr == -1, the onpulse region is
   considered. Otherwise all bins are considered.

   Normally longitude 0 and latitude 0 defines the centre of the
   projection, but this grid can be rotated by using rot_long and
   rot_lat (both in degrees), which are added to the longitude and
   latitude lines which are drawn.

   If subtract_pa_data != NULL, the PA-swing from this datafile from
   the data is subtracted. So this should be 1 subit, 1 freq channel,
   converted in PA etc, and the number of bins should match.

   Return 0 on error, 1 on success

*/
int make_projection_map_fromIQUV(datafile_definition datafile, float *map, int nrx, int nry, float background, int binnr, pulselongitude_regions_definition onpulse, int weighting, int projection, float rot_long, float rot_lat, datafile_definition *subtract_pa_data, verbose_definition verbose)
{
  int ok;
  long i, xi, yi, pulsenr;
  float longitude, L, P, latitude, x, y, weight;

  if(datafile.NrFreqChan != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR make_HammerAitoff_map_fromIQUV: Expected 1 frequency channel.");
    return 0;
  }

  if(datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
    if(datafile.NrPols != 8) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_HammerAitoff_map_fromIQUV: Expected 8 input polarizations when reading in data containing PA's and ellipticities.");
      return 0;
    }
  }else {
    if(datafile.NrPols != 4) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_HammerAitoff_map_fromIQUV: Expected 4 input polarizations.");
      return 0;
    }
  }

  if(subtract_pa_data != NULL) {
    if(subtract_pa_data->NrBins != datafile.NrBins) {
      printerror(verbose.debug, "ERROR make_HammerAitoff_map_fromIQUV: The reference PA-swing has a different number of pulse phase bins compared to the input data.");
      return 0;
    }
    if(subtract_pa_data->NrFreqChan != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_HammerAitoff_map_fromIQUV: The reference PA-swing should have a single frequency channel.");
      return 0;
    }
    if(subtract_pa_data->NrSubints != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR make_HammerAitoff_map_fromIQUV: The reference PA-swing should have a single frequency channel.");
      return 0;
    }
  }

  rot_long *= M_PI/180.0;
  rot_lat *= M_PI/180.0;

  /* Clear map, so we can add polarization points later. */
  for(xi = 0; xi < nrx; xi++) {
    for(yi = 0; yi < nry; yi++) {
      map[xi+nrx*yi] = 0;
    }
  }

  /* Now loop over the subintegrations and calculate PA etc. */
  for(pulsenr = 0; pulsenr < datafile.NrSubints; pulsenr++) {
    for(i = 0; i < (datafile.NrBins); i++) {
      ok = 1;
      if(binnr >= 0) {
	if(i != binnr)
	  ok = 0;
      }else if(binnr == -1) {
	if(checkRegions(i, &onpulse, 0, verbose) == 0) {
	  ok = 0;
	  /*
	    if(pulsenr == 0)
	    printf("XXX rejected bin %ld\n", i);
	  */
	}
      }
      if(ok) {
	/* tan longitude = U/Q -> so a number between -pi and pi. */
	if(datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
	  longitude = 2.0*datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+3*(datafile.NrBins)]*M_PI/180.0;
	  L = datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+1*(datafile.NrBins)];
	  /* tan latitude = V/L */
	  latitude = 2.0*datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+6*(datafile.NrBins)]*M_PI/180.0;
	  if(datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+4*(datafile.NrBins)] < 0 || datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+7*(datafile.NrBins)] < 0) {
	    longitude = latitude = sqrt(-1.0);
	  }
	}else {
	  longitude = atan2(datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+2*(datafile.NrBins)],datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+(datafile.NrBins)]);
	  L = sqrt(datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+2*(datafile.NrBins)]*datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+2*(datafile.NrBins)] + datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+(datafile.NrBins)]*datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+(datafile.NrBins)]);
	  /* tan latitude = V/L */
	  latitude = atan(datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+3*(datafile.NrBins)]/L);
	}
	if(subtract_pa_data != NULL) {
	  int pachannel_subtract_fin;
	  if(subtract_pa_data->poltype == POLTYPE_ILVPAdPATEldEl) {
	    pachannel_subtract_fin = 3;
	  }else {
	    pachannel_subtract_fin = subtract_pa_data->NrPols-2;
	  }
	  longitude -= 2.0*subtract_pa_data->data[i+subtract_pa_data->NrBins*(pachannel_subtract_fin)]*M_PI/180.0;  // Times to since PA is only over a domain of 180 deg
	}
	if(!isnan(latitude)) {
	  if(weighting) {
	    if(datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
	      P = datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+5*(datafile.NrBins)];
	    }else {
	      P = sqrt(L*L+datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+3*(datafile.NrBins)]*datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+3*(datafile.NrBins)]);
	    }
	  }
	  if(projection == 1) {
	    projectionHammerAitoff_xy(longitude, latitude, rot_long, rot_lat, &x, &y);
	    weight = 1;
	  }else if(projection == 2) {
	    projection_sphere_xy(longitude, latitude, rot_long, rot_lat, &x, &y, &weight);
	  }else if(projection == 3) {
	    projection_longlat_xy(longitude, latitude, rot_long, rot_lat, &x, &y);
	    /* Rescale coordinates such that they line up with other projections */
	    x /= 80.0;
	    y /= 80.0;
	  }
	  xi = 0.5*nrx + x*nrx/4.5;
	  yi = 0.5*nry + y*nry/2.25;
	  //	  if(map[xi+nrx*yi] < 1000)
	  if(weighting == 0)
	    map[xi+nrx*yi] += 1.0*weight;
	  else
	    map[xi+nrx*yi] += P*weight;
	  //	  else 
	  //	    printf("XXXX (%f %f) = (%ld %ld)  L=%f V=%f\n", longitude, latitude, xi, yi, L, datafile.data[i+datafile.NrBins*datafile.NrPols*pulsenr+3*(datafile.NrBins)]);
	}
      }
    }
  }

  return 1;
}

//START REGION DEVELOP

/*
  Given the input datafile (POLTYPE_ILVPAdPA), determine two datasets
  (mode1 and mode2) which contain the two polarization modes which
  combined reproduces datafile. The data is stored as IQUV data and
  mode1 and mode2 are assumed to be uninitialized.

  The variable method indicates how the splitting is done: 

  1 = Depolarization is because of two orthogonal linear modes. mode1
  is the the stronger mode at a given sample in the data. Stokes V is
  divided equally over the two modes.

  Returns 1 on success, 0 on error.
 */
int decompose_polarization_modes(datafile_definition datafile, datafile_definition *mode1, datafile_definition *mode2, int method, verbose_definition verbose)
{
  if(datafile.format != MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR decompose_polarization_modes: only works if data is loaded into memory.");
    return 0;
  }
  if(datafile.poltype != POLTYPE_ILVPAdPA) {
    printerror(verbose.debug, "ERROR decompose_polarization_modes: Data doesn't appear to have poltype ILVPAdPA."); 
    return 0;
  }
  int pa_chan, L_chan, V_chan, I_chan; // , dpa_chan
  I_chan = 0;
  L_chan = 1;
  V_chan = 2;
  pa_chan = 3;
  //  dpa_chan = 4;

  if(datafile.NrPols != 5) {
    printerror(verbose.debug, "ERROR decompose_polarization_modes: 5 polarization channels were expected, but there are %ld.", datafile.NrPols); 
    return 0;
  }
  

  cleanPSRData(mode1, verbose);
  cleanPSRData(mode2, verbose);
  copy_params_PSRData(datafile, mode1, verbose);
  copy_params_PSRData(datafile, mode2, verbose);
  mode1->format = MEMORY_format;
  mode2->format = MEMORY_format;
  mode1->poltype = POLTYPE_STOKES;
  mode2->poltype = POLTYPE_STOKES;
  mode1->NrPols = 4;
  mode2->NrPols = 4;
  mode1->data = (float *)malloc(datafile.NrSubints*datafile.NrBins*datafile.NrPols*datafile.NrFreqChan*sizeof(float));
  mode2->data = (float *)malloc(datafile.NrSubints*datafile.NrBins*datafile.NrPols*datafile.NrFreqChan*sizeof(float));
  if(mode1->data == NULL || mode2->data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR decompose_polarization_modes: Cannot allocate memory for data.");
    return 0;
  }


  long binnr, subint, freqchan;
  for(subint = 0; subint < datafile.NrSubints; subint++) {
    for(freqchan = 0; freqchan < datafile.NrFreqChan; freqchan++) {
      for(binnr = 0; binnr < datafile.NrBins; binnr++) {
	float I, L, V, pa;
	I  = datafile.data[binnr+datafile.NrBins*(I_chan +datafile.NrPols*(freqchan+datafile.NrFreqChan*subint))];
	L  = datafile.data[binnr+datafile.NrBins*(L_chan +datafile.NrPols*(freqchan+datafile.NrFreqChan*subint))];
	V  = datafile.data[binnr+datafile.NrBins*(V_chan +datafile.NrPols*(freqchan+datafile.NrFreqChan*subint))];
	pa = datafile.data[binnr+datafile.NrBins*(pa_chan+datafile.NrPols*(freqchan+datafile.NrFreqChan*subint))];
	//	printf("XXXXX: %f\n", pa);
	float i1, i2, q1, q2, u1, u2, v1, v2;
	if(method == 1) {
	  float L1, L2;
	  L1 = 0.5*(I-fabs(V)+L);
	  L2 = I-fabs(V)-L1;
	  q1 = L1*cos(2.0*pa*M_PI/180.0);
	  q2 = L2*cos(2.0*(pa+90.0)*M_PI/180.0);
	  u1 = L1*sin(2.0*pa*M_PI/180.0);
	  u2 = L2*sin(2.0*(pa+90.0)*M_PI/180.0);
	  v1 = 0.5*V;
	  v2 = 0.5*V;
	  //	  i1 = sqrt(L1*L1+v1*v1);
	  //	  i2 = sqrt(L2*L2+v2*v2);
	  i1 = L1;
	  i2 = L2;
	}else {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR decompose_polarization_modes: Invalid method is specified.");
	  return 0;
	}
	mode1->data[binnr+mode1->NrBins*(0+mode1->NrPols*(freqchan+mode1->NrFreqChan*subint))] = i1;
	mode2->data[binnr+mode2->NrBins*(0+mode2->NrPols*(freqchan+mode2->NrFreqChan*subint))] = i2;
	mode1->data[binnr+mode1->NrBins*(1+mode1->NrPols*(freqchan+mode1->NrFreqChan*subint))] = q1;
	mode2->data[binnr+mode2->NrBins*(1+mode2->NrPols*(freqchan+mode2->NrFreqChan*subint))] = q2;
	mode1->data[binnr+mode1->NrBins*(2+mode1->NrPols*(freqchan+mode1->NrFreqChan*subint))] = u1;
	mode2->data[binnr+mode2->NrBins*(2+mode2->NrPols*(freqchan+mode2->NrFreqChan*subint))] = u2;
	mode1->data[binnr+mode1->NrBins*(3+mode1->NrPols*(freqchan+mode1->NrFreqChan*subint))] = v1;
	mode2->data[binnr+mode2->NrBins*(3+mode2->NrPols*(freqchan+mode2->NrFreqChan*subint))] = v2;
	//	printf("IQUV mode2: %f %f %f %f\n", i2, q2, u2, v2);
      }
    }
  }

  return 1;
}

