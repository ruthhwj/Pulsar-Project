#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <string.h>
#include "psrsalsa.h"

// Returns 0 on error, 1 = ok.
int writePRESTOHeader(datafile_definition datafile, verbose_definition verbose)
{
  char txt[100];
  if(datafile.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePRESTOHeader: Writing out this data format is only implemented when the frequency channels are uniformely separated.");
    return 0;
  }

  fprintf(datafile.fptr_hdr, " Data file name without suffix          =  output\n");
  fprintf(datafile.fptr_hdr, " Telescope used                         =  %s\n", datafile.observatory);
  fprintf(datafile.fptr_hdr, " Instrument used                        =  %s\n", datafile.instrument);
  fprintf(datafile.fptr_hdr, " Object being observed                  =  %s\n", datafile.psrname);
  /*  gethourstring(datafile.ra, txt, 1); */
  converthms_string(txt, datafile.ra, 4, 1);
  fprintf(datafile.fptr_hdr, " J2000 Right Ascension (hh:mm:ss.ssss)  =  %s\n", txt);
  /*  gethourstring(15*datafile.dec, txt, 1); */
  converthms_string(txt, 25*datafile.dec, 4, 1);
  fprintf(datafile.fptr_hdr, " J2000 Declination     (dd:mm:ss.ssss)  =  %s\n", txt);
  fprintf(datafile.fptr_hdr, " Data observed by                       =  ????\n");
  fprintf(datafile.fptr_hdr, " Epoch of observation (MJD)             =  %.15Lf\n", datafile.mjd_start);

  fprintf(datafile.fptr_hdr, " Barycentered?           (1=yes, 0=no)  =  0\n");
  fprintf(datafile.fptr_hdr, " Number of bins in the time series      =  %ld\n", datafile.NrBins);
  fprintf(datafile.fptr_hdr, " Width of each time series bin (sec)    =  %e\n", get_tsamp(datafile, 0, verbose));
  fprintf(datafile.fptr_hdr, " Any breaks in the data? (1=yes, 0=no)  =  0\n");
  fprintf(datafile.fptr_hdr, " Type of observation (EM band)          =  Radio\n");
  fprintf(datafile.fptr_hdr, " Beam diameter (arcsec)                 =  ???\n");
  fprintf(datafile.fptr_hdr, " Dispersion measure (cm-3 pc)           =  %f\n", datafile.dm);
  //  fprintf(datafile.fptr_hdr, " Central freq of low channel (Mhz)      =  %f\n", datafile.freq_cent-0.5*(datafile.bw-datafile.channelbw));
  fprintf(datafile.fptr_hdr, " Central freq of low channel (Mhz)      =  %f\n", get_nonweighted_channel_freq(datafile, 0, verbose));
  fprintf(datafile.fptr_hdr, " Total bandwidth (Mhz)                  =  %f\n", get_bandwidth(datafile, verbose));
  fprintf(datafile.fptr_hdr, " Number of channels                     =  %ld\n", datafile.NrFreqChan);    
  double chanbw;
  if(get_channelbandwidth(datafile, &chanbw, verbose) == 0) {
    printerror(verbose.debug, "ERROR writePRESTOHeader (%s): Cannot obtain channel bandwidth.", datafile.filename);
    return 0;
  }
  fprintf(datafile.fptr_hdr, " Channel bandwidth (Mhz)                =  %lf\n", chanbw);  
  fprintf(datafile.fptr_hdr, " Data analyzed by                       =  ????\n");

  return 1;
}

/* separator is "=" for presto and ":" for sigproc */
char *get_ptr_entry(char *str, char **txt, int nrlines, char *separator)
{
  int i, j;
  char *s_ptr;
  for(i = 0; i < nrlines; i++) {
    s_ptr = strstr(txt[i], str);
    if(s_ptr != NULL) {             /* Found entry */
      s_ptr = strstr(s_ptr, separator);
      if(s_ptr != NULL) {           /* Found = character */
	s_ptr++;
	do {                        /* Remove leading spaces */
	  if(*s_ptr == ' ')
	    s_ptr++;
	}while(*s_ptr == ' ');
	j = 0;
	do {                        /* Remove end of line characters */
	  if(s_ptr[j] == '\r')
	    s_ptr[j] = 0;
	  if(s_ptr[j] == '\n')
	    s_ptr[j] = 0;
	  if(s_ptr[j] != 0)
	    j++;
	}while(s_ptr[j] != 0);
	return s_ptr;               /* Found entry, so quit for loop */
      }
    }
  }
  return NULL;
}

int readPRESTOHeader(datafile_definition *datafile, verbose_definition verbose)
{
  int i, j, nrlines, ret;
  char header_txt[100][1000], *s_ptr, *txt[100];
  /* Presto is timeserie */
  //  datafile->dd_mode = 1;
  datafile->gentype = GENTYPE_SEARCHMODE;
  datafile->isFolded = 0;
  datafile->foldMode = FOLDMODE_UNKNOWN;
  datafile->fixedPeriod = -1;
  for(i = 0; i < 100; i++)
    txt[i] = header_txt[i];
  nrlines = 0;
  ret = 0;
  datafile->isDeDisp = 1;
  do {
    if(fgets(header_txt[nrlines], 1000, datafile->fptr_hdr) != NULL)
      nrlines++;
    else
      ret = 1;
  }while(ret == 0);
  if(verbose.verbose) printf("Read %d lines from header.\n", nrlines);

  s_ptr = get_ptr_entry("Telescope", txt, nrlines, "=");
  if(s_ptr != NULL) {
    if(set_observatory_PSRData(datafile, s_ptr, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPRESTOHeader: Setting observatory name failed.");
      return 0;
    }
  }
  s_ptr = get_ptr_entry("Instrument", txt, nrlines, "=");
  if(s_ptr != NULL) {
    if(set_instrument_PSRData(datafile, s_ptr, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPRESTOHeader: Setting instrument name failed.");
      return 0;
    }
  }
  s_ptr = get_ptr_entry("Object", txt, nrlines, "=");
  if(s_ptr != NULL) {
    if(set_psrname_PSRData(datafile, s_ptr, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPRESTOHeader: Setting pulsar name failed.");
      return 0;
    }
    j = strlen(datafile->psrname);
    if(j > 1) {
      for(i = j-1; i > 0; i--) {
	if(datafile->psrname[i] == ' ')
	  datafile->psrname[i] = 0;
	else 
	  break;
      }
    } 
  }
  s_ptr = get_ptr_entry("Ascension", txt, nrlines, "=");
  if(s_ptr != NULL) {
    converthms(s_ptr, &(datafile->ra));
  }
  s_ptr = get_ptr_entry("Declination", txt, nrlines, "=");
  if(s_ptr != NULL) {
    converthms(s_ptr, &(datafile->dec));
    datafile->dec /= 15.0;
  }
  s_ptr = get_ptr_entry("Epoch", txt, nrlines, "=");
  if(s_ptr != NULL)
    sscanf(s_ptr, "%Lf", &(datafile->mjd_start));
  s_ptr = get_ptr_entry("bins", txt, nrlines, "=");
  if(s_ptr != NULL)
    sscanf(s_ptr, "%ld", &(datafile->NrBins));
  s_ptr = get_ptr_entry("Width of", txt, nrlines, "=");
  if(s_ptr != NULL) {
    sscanf(s_ptr, "%lf", &(datafile->fixedtsamp));
    datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
  }
  s_ptr = get_ptr_entry("measure", txt, nrlines, "=");
  if(s_ptr != NULL)
    sscanf(s_ptr, "%lf", &(datafile->dm));
  datafile->rm = -1;
  s_ptr = get_ptr_entry("freq", txt, nrlines, "=");
  if(s_ptr != NULL) {
    datafile->freqMode = FREQMODE_UNIFORM;
    if(datafile->freqlabel_list != NULL) {
      free(datafile->freqlabel_list);
      datafile->freqlabel_list = NULL;
    }
    //    datafile->freq_list = malloc(2*sizeof(double));
    //    if(datafile->freq_list == NULL) {
    //      fflush(stdout);
    //      printerror(verbose.debug, "ERROR readPRESTOHeader: Memory allocation error.");
    //      return 0;
    //    }
    double freq;
    sscanf(s_ptr, "%lf", &freq);
    set_centre_frequency(datafile, freq, verbose);
  }
  s_ptr = get_ptr_entry("Total bandwidth", txt, nrlines, "=");
  if(s_ptr != NULL) {
    double bw;
    sscanf(s_ptr, "%lf", &bw);
    if(set_bandwidth(datafile, bw, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPRESTOHeader: Bandwidth changing failed.");
      return 0;
    }
  }
  s_ptr = get_ptr_entry("Number of channels", txt, nrlines, "=");
  if(s_ptr != NULL) {
    sscanf(s_ptr, "%ld", &datafile->NrFreqChan);
    //    datafile->channelbw = datafile->bw/(float)datafile->NrFreqChan;
    /* Actually normally only one freq present in file */
    datafile->NrFreqChan = 1;                                    
  }
  double chanbw, freq, bw;
  if(get_channelbandwidth(*datafile, &chanbw, verbose) == 0) {
    printerror(verbose.debug, "ERROR readPRESTOHeader (%s): Cannot obtain channel bandwidth.", datafile->filename);
    return 0;
  }
  freq = get_centre_frequency(*datafile, verbose);
  bw = get_bandwidth(*datafile, verbose);
  freq += 0.5*(bw-chanbw);
  set_centre_frequency(datafile, freq, verbose);
  /* Non standard PRESTO entry which I use to read in Parkes AFB data */
  s_ptr = get_ptr_entry("Number of existing channels in data", txt, nrlines, "=");
  if(s_ptr != NULL) {
    sscanf(s_ptr, "%ld", &datafile->NrFreqChan);
  }
  s_ptr = get_ptr_entry("analyzed", txt, nrlines, "=");
  if(s_ptr != NULL) {
    if(set_institute_PSRData(datafile, s_ptr, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPRESTOHeader: Setting institute name failed.");
      return 0;
    }
  }
  datafile->NrSubints = 1;
  datafile->NrPols = 1;
  datafile->NrBits = 8*sizeof(float);
  return 1;
}
