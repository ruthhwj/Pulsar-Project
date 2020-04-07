#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <math.h>
#include <string.h>
#include "psrsalsa.h"

// Type can be:
//   1: Maybe this is the calibrated version?
//   2: This is a more simple format: one line header per pulse
int writeAOAsciiHeader(datafile_definition datafile, int type, verbose_definition verbose)
{
  float degpersamp;
  double period;

  if(type != 1) {
    printerror(verbose.debug, "ERROR writeAOAsciiHeader: Writing of type %d data is not implemented", type);
    return 0;
  }


  if(get_period(datafile, 0, &period, verbose) == 2) {
    printerror(verbose.debug, "ERROR writeAOAsciiHeader (%s): Cannot obtain period", datafile.filename);
    return 0;
  }

  if(datafile.freqMode != FREQMODE_UNIFORM) {
    printerror(verbose.debug, "ERROR writeAOAsciiHeader (%s): Frequency channels do not appear to be uniformly distributed", datafile.filename);
    return 0;
  }

  degpersamp = (360.0 * get_tsamp(datafile, 0, verbose) / period);
  if(fabs(degpersamp) < 0.001) {
    printwarning(verbose.debug, "The degrees per sample was set to: 360*%f/%f = %f", get_tsamp(datafile, 0, verbose), period, degpersamp);
    degpersamp = 360.0/(float)datafile.NrBins;
    printwarning(verbose.debug, "I set it to: 360/%ld = %f", datafile.NrBins, degpersamp);
  }
  /*  fprintf(datafile.fptr_hdr, "%ld  %d  %f  %d\n", datafile.NrSubints, datafile.NrBins, degpersamp, datafile.NrPols); */
  /* Always wite out four polarizations, although it could be crap */
  fprintf(datafile.fptr_hdr, " %ld  %ld  %f  %d\n", datafile.NrSubints, datafile.NrBins, degpersamp, 4);
  if(datafile.psrname == NULL) {
    fprintf(datafile.fptr_hdr, " B????+??\n");
  }else if(strlen(datafile.psrname) == 0) {
    fprintf(datafile.fptr_hdr, " B????+??\n");
  }else {
    fprintf(datafile.fptr_hdr, " %s\n", datafile.psrname);
  }
  if(datafile.observatory == NULL) {
    fprintf(datafile.fptr_hdr, " ?/");
  }else if(strlen(datafile.observatory) == 0) {
    fprintf(datafile.fptr_hdr, " ?/");
  }else {
    fprintf(datafile.fptr_hdr, " %s/", datafile.observatory);
  }
  if(datafile.institute == NULL) {
    fprintf(datafile.fptr_hdr, "?\n");
  }else if(strlen(datafile.institute) == 0) {
    fprintf(datafile.fptr_hdr, "?\n");
  }else {
    fprintf(datafile.fptr_hdr, "%s\n", datafile.institute);
  }

  /* No idea why 0 should be there, but it must (at least for the recent file format) */
  fprintf(datafile.fptr_hdr, " %d 0  %f  %f  %12.10f  %16.10Lf\n", (int)datafile.mjd_start, get_centre_frequency(datafile, verbose), get_bandwidth(datafile, verbose), period, datafile.mjd_start);
  return 1;
}

int writeAOAsciifile(datafile_definition datafile, float *data, int type, verbose_definition verbose)
{
  long n, f, p, i;
  float sample, SigI, AverageI, ScaleVal, maxI;

  if(type != 1) {
    printerror(verbose.debug, "ERROR writeAOAsciifile: Writing of type %d data is not implemented", type);
    return 0;
  }

  if(datafile.NrFreqChan != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeAOAsciifile: Only one freq channel is supported in this data format");
    return 0;
  }
  SigI = 0;
  AverageI = 0;
  ScaleVal = 1;
  f = 0;
  maxI = data[0];
  for(n = 0; n < datafile.NrSubints; n++) {
    for(i = 0; i < datafile.NrBins; i++) {
      for(p = 0; p < datafile.NrPols; p++) {
	if(fabs(data[i+datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))]) > maxI) 
	  maxI = fabs(data[i+datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))]);
      }
    }
  }
  printwarning(verbose.debug, "WARNING writeAOAsciifile: Data scaled by 1/%e", maxI); 

  for(n = 0; n < datafile.NrSubints; n++) {
    if(verbose.verbose && verbose.nocounters == 0) printf("writeAOAsciifile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
    for(i = 0; i < datafile.NrBins; i++) {
      for(p = 0; p < datafile.NrPols; p++) {
	sample = data[i+datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))];
	if(datafile.NrPols == 1) {
	  fprintf(datafile.fptr, "%f  %f  %f  %f ", sample/maxI, 0.0,0.0,0.0);
	}else {
	  if(p == 3)
	    fprintf(datafile.fptr, "%f ", sample);
	  else
	    fprintf(datafile.fptr, "%f  ", sample);
	}
      }
      if(((i+1)%2) == 0) fprintf(datafile.fptr, "\n");
    }
    fprintf(datafile.fptr, "%f  %f  %ld\n", SigI, AverageI/ScaleVal, n+1);
  }
  if(verbose.verbose && verbose.nocounters == 0) printf("writeAOAsciifile: done                 \n");  
  return 1;
}

int readAOAsciiHeader(datafile_definition *datafile, int type, verbose_definition verbose)
{
  double mjdfrac, deg_per_sample;
  long c;
  int ret, i;
  char *txt, *txt2;
  long filepos;

  txt = malloc(1001);
  if(txt == NULL) {
    printerror(verbose.debug, "ERROR readAOAsciiHeader: Memory allocation error");
    return 0;
  }

  if(type < 1 || type > 2) {
    printerror(verbose.debug, "ERROR readAOAsciiHeader: Reading of type %d data is not implemented", type);
    free(txt);
    return 0;
  }

  if(type == 1) {
    ret = fscanf(datafile->fptr,"%ld %ld %lf %ld %s", &(datafile->NrSubints), &(datafile->NrBins), &deg_per_sample, &(datafile->NrPols), txt);
    if(ret != 5) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot read header! (nrarguments=%d)", ret);
      free(txt);
      return 0;
    }
    if(set_psrname_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Setting pulsar name failed.");
      free(txt);
      return 0;
    }
    
    do {
      i = fgetc(datafile->fptr);
    }while(i != EOF && i != '\n');
    //  fgets(txt, 900, datafile->fptr);
    fscanf(datafile->fptr, "%s", txt);
    if(set_scanID_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Setting scan ID failed.");
      free(txt);
      return 0;
    }
    do {
      i = fgetc(datafile->fptr);
    }while(i != EOF && i != '\n');  // Skip rest of line
    /*  printf("%d '%s'\n", i, datafile->scanid); */
    
    datafile->isFolded = 1;
    datafile->foldMode = FOLDMODE_FIXEDPERIOD;
    datafile->tsubMode = TSUBMODE_FIXEDTSUB;
    if(datafile->tsub_list != NULL)
      free(datafile->tsub_list);
    datafile->tsub_list = (double *)malloc(sizeof(double));
    if(datafile->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Memory allocation error");
      free(txt);
      return 0;
    }
    datafile->tsub_list[0] = 0;
    
    datafile->freqMode = FREQMODE_UNIFORM;
    double freq, bw;
    ret = fscanf(datafile->fptr,"%Lf %ld %lf %lf %lf %lf", &(datafile->mjd_start), &c, &freq, &bw, &(datafile->fixedPeriod), &mjdfrac);
    if(ret != 6) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot read header! (nrarguments=%d)", ret);
      free(txt);
      return 0;
    }
    // Skip rest of line as there is a mysterious other integer in at least some cases, but not always....
    do {
      i = fgetc(datafile->fptr);
    }while(i != EOF && i != '\n'); // Rest of line
    set_centre_frequency(datafile, freq, verbose);
    if(set_bandwidth(datafile, bw, verbose) == 0) {
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Bandwidth changing failed.");
      return 0;
    }
    datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
    double period;
    if(get_period(*datafile, 0, &period, verbose) == 2) {
      printerror(verbose.debug, "ERROR readAOAsciiHeader (%s): Cannot obtain period", datafile->filename);
      free(txt);
      return 0;
    }
    datafile->fixedtsamp = period*deg_per_sample/360.0;
    /*
      if (verbose) printf("npulses=%ld\nnbin=%ld\ndeg per sample=%f\nnpol=%ld\n", *npulses, *nbin, *deg_per_sample,*npol); 
      if (verb) printf("Pulsar=%s\n", pulsar_name);
      if (verb) printf("Comment=%s\n", comment);
      if (verb) printf("MJD=%ld\n%ld\nFreq=%f\nBandwidth=%f\nPulse period=%f\nMJD=%f\n",*mjd, c, *freq,*bw,*pulse_period,*mjdfrac);
    */
    datafile->NrFreqChan = 1;
    datafile->dm = -1;
    datafile->rm = 0;

  }else if(type == 2) {
    int nrwords;
    double value;
    // Read first line
    //# 53156.0 68955.9999996 0.6451611911 1 339.500 52.667 1024 3 1 J0815+0939
    if(fgets(txt, 1000, datafile->fptr) == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot read first line of header");
      free(txt);
      return 0;
    }
    filepos = ftell(datafile->fptr);
    if(txt[0] != '#' || txt[1] != ' ') {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Header doesn't start with a # and a space");
      free(txt);
      return 0;
    }
    txt2 = pickWordFromString(txt, 1, &nrwords, 0, ' ', verbose);
    if(nrwords != 11) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Expected 11 words in header line, got %d", nrwords);
      free(txt);
      return 0;
    }
    txt2 = pickWordFromString(txt, 2, &nrwords, 0, ' ', verbose);
    i = sscanf(txt2, "%Lf", &(datafile->mjd_start));
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 2nd word");
      free(txt);
      return 0;
    }
    txt2 = pickWordFromString(txt, 3, &nrwords, 0, ' ', verbose);
    i = sscanf(txt2, "%lf", &mjdfrac);
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 3rd word");
      free(txt);
      return 0;
    }
    datafile->isFolded = 1;
    datafile->foldMode = FOLDMODE_FIXEDPERIOD;
    txt2 = pickWordFromString(txt, 4, &nrwords, 0, ' ', verbose);
    i = sscanf(txt2, "%lf", &(datafile->fixedPeriod));
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 4th word");
      free(txt);
      return 0;
    }
    txt2 = pickWordFromString(txt, 5, &nrwords, 0, ' ', verbose); // Not sure what this is supposed to be. It was 1 for a pulse stack
    i = sscanf(txt2, "%lf", &value);
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 5th word");
      free(txt);
      return 0;
    }
    if(value < 1-1e-3 || value > 1+1e-3) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Expected 5th word to be 1");
      free(txt);
      return 0;
    }
    datafile->freqMode = FREQMODE_UNIFORM;
    txt2 = pickWordFromString(txt, 6, &nrwords, 0, ' ', verbose);
    double freq;
    i = sscanf(txt2, "%lf", &freq);
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 6th word");
      free(txt);
      return 0;
    }
    set_centre_frequency(datafile, freq, verbose);
    txt2 = pickWordFromString(txt, 7, &nrwords, 0, ' ', verbose);
    i = sscanf(txt2, "%lf", &(datafile->dm));
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 7th word");
      free(txt);
      return 0;
    }
    txt2 = pickWordFromString(txt, 8, &nrwords, 0, ' ', verbose);
    i = sscanf(txt2, "%ld", &(datafile->NrBins));
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 8th word");
      free(txt);
      return 0;
    }
    txt2 = pickWordFromString(txt, 9, &nrwords, 0, ' ', verbose); // Not sure what this is supposed to be. It was 3 for a pulse stack
    i = sscanf(txt2, "%lf", &value);
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 9th word");
      free(txt);
      return 0;
    }
    if(value < 3-1e-3 || value > 3+1e-3) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Expected 9th word to be 3");
      free(txt);
      return 0;
    }
    txt2 = pickWordFromString(txt, 10, &nrwords, 0, ' ', verbose); // Not sure what this is supposed to be. It was 1 for a pulse stack
    i = sscanf(txt2, "%lf", &value);
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 10th word");
      free(txt);
      return 0;
    }
    if(value < 1-1e-3 || value > 1+1e-3) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Expected 10th word to be 1");
      free(txt);
      return 0;
    }
    txt2 = pickWordFromString(txt, 11, &nrwords, 0, ' ', verbose); // Not sure what this is supposed to be. It was 1 for a pulse stack
    if(set_psrname_PSRData(datafile, txt2, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Setting pulsar name failed.");
      free(txt);
      return 0;
    }

    datafile->NrFreqChan = 1;

    datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
    double period;
    if(get_period(*datafile, 0, &period, verbose) == 2) {
      printerror(verbose.debug, "ERROR readAOAsciiHeader (%s): Cannot obtain period", datafile->filename);
      free(txt);
      return 0;
    }
    datafile->fixedtsamp = period/(double)datafile->NrBins;
    // Read in 1st line of data
    if(fgets(txt, 1000, datafile->fptr) == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot read first line of header");
      free(txt);
      return 0;
    }

    txt2 = pickWordFromString(txt, 1, &nrwords, 0, ' ', verbose);
    datafile->NrPols = nrwords - 1;
    if(datafile->NrPols <= 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Expected at least two words on 2nd line");
      free(txt);
      return 0;
    }
    i = sscanf(txt2, "%lf", &value);
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 1st word from 2nd line");
      free(txt);
      return 0;
    }
    if(value < 1-1e-3 || value > 1+1e-3) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Expected 1st word of 2nd line to be 1");
      free(txt);
      return 0;
    }

    datafile->tsubMode = TSUBMODE_FIXEDTSUB;
    if(datafile->tsub_list != NULL)
      free(datafile->tsub_list);
    datafile->tsub_list = (double *)malloc(sizeof(double));
    if(datafile->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readAOAsciiHeader: Memory allocation error");
      free(txt);
      return 0;
    }
    datafile->tsub_list[0] = 0;
    // Goto next header (of 2nd subint)
    if(skipLinesInFile(datafile->fptr, datafile->NrBins-1, verbose) == 0) {
      datafile->NrSubints = 1;
    }else {
      if(fgets(txt, 1000, datafile->fptr) == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot read 2nd header line");
	free(txt);
	return 0;
      }

      txt2 = pickWordFromString(txt, 3, &nrwords, 0, ' ', verbose);
      if(nrwords != 11) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readAOAsciiHeader: 2nd header line does not appear to have 11 words");
	free(txt);
	return 0;
      }
      double mjdfrac2;
      i = sscanf(txt2, "%lf", &mjdfrac2);
      if(i != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 3rd word on 2nd header line");
	free(txt);
	return 0;
      }
      datafile->tsub_list[0] = mjdfrac2 - mjdfrac;

      // Determine nr lines in file
      long nrlines;
      rewind(datafile->fptr);
      if(ascii_file_stats(datafile->fptr, 0, &nrlines, 10000, 0, NULL, verbose) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot obtain file statistics");
	free(txt);
	return 0;
      }
      datafile->NrSubints = nrlines/(datafile->NrBins+1);
      if(nrlines - datafile->NrSubints*(datafile->NrBins+1) != 0) {
	printwarning(verbose.debug, "ERROR readAOAsciiHeader: The number of lines in the file is suspicious, so the data might be corrupt");
      }
    }


/* Assuming the lines in the text file (fin is the file pointer and
   the file is already opened) are less than maxlinelength characters
   long, determine the nr of lines in the file. If the line starts
   with character skipChar (maybe set to '#'), the line is ignored
   (set to zero to ignore this feature). If autoNrColumns is set, the
   number of columns (nrColumns) is determined as well. An error is
   generated in the nr of columns appear to change. If autoNrColumns
   is set to zero, the nr of columns is compared with that in
   *nrColumns, unless that is set to the NULL pointer. If nrColumns is
   set to a negative number, each line is expected to have at least
   -nrColumns of columns. The function returns 1 on success, 0 on
   error. */
//int ascii_file_stats(FILE *fin, char skipChar, long *nrlines, int maxlinelength, int autoNrColumns, int *nrColumns, verbose_definition verbose);

    //    


    //datafile->NrSubints

    // Go back to position of first line of data
    fseek(datafile->fptr, filepos, SEEK_SET);
  }else {
    printerror(verbose.debug, "ERROR readAOAsciiHeader: Reading of type %d data is not implemented", type);
    return 0;
  }
  datafile->mjd_start += mjdfrac/86400.0;
  free(txt);
  return 1;
}



int readAOAsciifile(datafile_definition datafile, float *data, int type, verbose_definition verbose)
{
  int ret, intvalue, nrwords;
  long p, n, i, j, f;
  float sample, c1, c2, c3;
  char *txt, *txt2;
  if(type < 1 || type > 2) {
    printerror(verbose.debug, "ERROR readAOAsciiHeader: Reading of type %d data is not implemented", type);
    return 0;
  }
  if(datafile.NrFreqChan != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readAOAsciifile: Only one frequency channel is supported.");
    return 0;
  }

  txt = malloc(1001);
  if(txt == NULL) {
    printerror(verbose.debug, "ERROR readAOAsciiHeader: Memory allocation error");
    return 0;
  }
  if(type == 1) {
    if(verbose.verbose) printf("Reading file in AO ascii format 1.");
    f = 0;
    for(p = 0; p < datafile.NrPols; p++) {
      fseek(datafile.fptr, datafile.datastart, SEEK_SET);
      if(verbose.verbose) printf("\n  Polarization %ld/%ld\n", p+1, datafile.NrPols); 
      for(n = 0; n < datafile.NrSubints; n++) {
	if(verbose.verbose && verbose.nocounters == 0) printf("  readAOAsciifile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
	for(i = 0; i < datafile.NrBins; i++) {
	  /*	    if(verb) printf("bin %d/%d\n", i+1, nbin); */
	  for(j = 0; j < datafile.NrPols; j++) {
	    fscanf(datafile.fptr, "%f", &sample);
	    if(j == p) {
	      data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] = sample;
	    }
	  }
	}
	fscanf(datafile.fptr, "%f %f %f", &c1, &c2, &c3); /* Ignore RMS values */
	/*      if(p == 0) fprintf(tsys_file, "%f %f %f\n", c1, c2, c3); */
      }
    }
  }else if(type == 2) {
    rewind(datafile.fptr);
    f = 0;
    for(n = 0; n < datafile.NrSubints; n++) {
      // Read header line
      if(fgets(txt, 1000, datafile.fptr) == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot read header line of subint %ld", n);
	free(txt);
	return 0;
      }
      for(i = 0; i < datafile.NrBins; i++) {
	// Read data line
	if(fgets(txt, 1000, datafile.fptr) == NULL) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot read data line of subint %ld and bin %ld", n, i);
	  free(txt);
	  return 0;
	}
	txt2 = pickWordFromString(txt, 1, &nrwords, 0, ' ', verbose);
	if(nrwords != datafile.NrPols+1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readAOAsciiHeader: Expected %d words in data line, got %d", datafile.NrPols+1, nrwords);
	  free(txt);
	  return 0;
	}
	ret = sscanf(txt2, "%d", &intvalue);
	if(ret != 1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret 1st words in data line");
	  free(txt);
	  return 0;
	}
	if(intvalue != i+1) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readAOAsciiHeader: Bin numbering appears to be out of order");
	  free(txt);
	  return 0;
	}
	for(p = 0; p < datafile.NrPols; p++) {
	  txt2 = pickWordFromString(txt, 2+p, &nrwords, 0, ' ', verbose);
	  ret = sscanf(txt2, "%f", &(data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i]));
	  if(ret != 1) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR readAOAsciiHeader: Cannot interpret dample in data line");
	    free(txt);
	    return 0;
	  }
	}
      }
    }
  }

  if(verbose.verbose && verbose.nocounters == 0) printf("  Reading is done.                  \n");
  free(txt);
  return 1;
}
