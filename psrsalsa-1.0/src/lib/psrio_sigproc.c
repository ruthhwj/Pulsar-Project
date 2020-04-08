//START REGION RELEASE
#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "psrsalsa.h"

//START REGION DEVELOP
static int ingridflavour = 0;
static char *ingridflavour_parfile = NULL;

//START REGION RELEASE
// A sigproc header file is basically an int (length of string),
// followed by the string. This defines what is written after this
// which is read in by readSigprocHeader(). id should not yet be
// allocated. If != NULL, it is freed first.
//
// Return 0 = error, 1 = success
int readSigprocHeader_readParamID(FILE *fin, char **id, verbose_definition verbose)
{
  int idlength;
  if(fread(&idlength, sizeof(int), 1, fin) != 1) {
    printerror(verbose.debug,"ERROR readSigprocHeader_readParamID: Read error");
    return 0;
  }
  if(idlength < 1 || idlength > 80) {
    printerror(verbose.debug,"ERROR readSigprocHeader_readParamID: ID length (%d) makes no sense", idlength);
    return 0;
  }
  if(*id != NULL) {
    free(*id);
  }
  *id = malloc(idlength+1);
  if(*id == NULL) {
    printerror(verbose.debug,"ERROR readSigprocHeader_readParamID: Memory allocation error");
    return 0;
  }
  if(fread(*id, 1, idlength, fin) != idlength) {
    printerror(verbose.debug,"ERROR readSigprocHeader_readParamID: Read error");
    return 0;
  }
  (*id)[idlength] = 0;
  return 1;
}

//START REGION DEVELOP
/* Add the Ingrid flavour to the sigproc format. This adds the pulse
   phase of the start mjd to the header. Set parfile to the location
   of the par file. The name of the parfile should be stay in the same
   memory location, even after calling this funtion. */
void add_IngridFlavour_SigprocAsciiFormat(int val, char *parfile)
{
  ingridflavour = val;
  ingridflavour_parfile = parfile;
}

/* Should be defined in psrio_presto.c */
//char *get_ptr_entry(char *str, char **txt, int nrlines, char *separator);

//START REGION RELEASE
int readSigprocHeader(datafile_definition *datafile, verbose_definition verbose)
{
  char *id;
  int dummy_int, nifs;
  double freq_chan1, freq_chanbw, dummy_double;
  long int dummy_long;

  // Default sigproc format is not folded
  datafile->gentype = GENTYPE_SEARCHMODE;
  datafile->isFolded = 0;
  datafile->foldMode = FOLDMODE_UNKNOWN;
  datafile->fixedPeriod = -1;
  datafile->tsubMode = TSUBMODE_FIXEDTSUB;
  if(datafile->tsub_list != NULL)
    free(datafile->tsub_list);
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    printerror(verbose.debug,"ERROR readSigprocHeader: Memory allocation error");
    exit(-1);
  }
  datafile->tsub_list[0] = 0;

  id = NULL;
  freq_chan1 = -1;
  freq_chanbw = 0;
  nifs = 1;

  if(readSigprocHeader_readParamID(datafile->fptr, &id, verbose) == 0) {
    printerror(verbose.debug, "ERROR readSigprocHeader: Cannot read ID");
    return 0;
  }
  if(verbose.debug) {
    printf("DEBUG: ID=%s\n", id);
  }
  if(strcmp(id, "HEADER_START") != 0) {
    printerror(verbose.debug, "ERROR readSigprocHeader: First ID is not the expected ID (%s != HEADER_START).", id);
    return 0;
  }

  do {
    if(readSigprocHeader_readParamID(datafile->fptr, &id, verbose) == 0) {
      printerror(verbose.debug, "ERROR readSigprocHeader: Cannot read ID");
      return 0;
    }
    if(verbose.debug) {
      printf("DEBUG: ID=%s\n", id);
    }
    if(strcmp(id, "machine_id") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      char txt[100];
      sprintf(txt, "%d", dummy_int);
      if(set_instrument_PSRData(datafile, txt, verbose) == 0) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Setting instrument failed");
	return 0;
      }
    }else if(strcmp(id, "telescope_id") == 0) { // Note - sigproc has its own codes. See aliases.c
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      char txt[100];
      sprintf(txt, "%d", dummy_int);
      if(set_observatory_PSRData(datafile, txt, verbose) == 0) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Setting observatory failed");
	return 0;
      }
    }else if(strcmp(id, "data_type") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      if(dummy_int == 3) {
	datafile->gentype = GENTYPE_SUBINTEGRATIONS;
	datafile->isFolded = 1;
      }
    }else if(strcmp(id, "fch1") == 0) {
      if(fread(&freq_chan1, sizeof(double), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%lf\n", id, freq_chan1);
      }
    }else if(strcmp(id, "foff") == 0) {
      if(fread(&freq_chanbw, sizeof(double), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%lf\n", id, freq_chanbw);
      }
    }else if(strcmp(id, "nchans") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      datafile->NrFreqChan = dummy_int;
    }else if(strcmp(id, "source_name") == 0) {
      char *txt;
      txt = NULL;
      if(readSigprocHeader_readParamID(datafile->fptr, &txt, verbose) == 0) {
	printerror(verbose.debug, "ERROR readSigprocHeader: Cannot read string");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%s\n", id, txt);
      }
      if(set_psrname_PSRData(datafile, txt, verbose) == 0) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Setting observatory failed");
	return 0;
      }
      free(txt);
    }else if(strcmp(id, "src_raj") == 0) {
      int hour, min;
      double sec;
      if(fread(&(datafile->ra), sizeof(double), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%lf\n", id, datafile->ra);
      }
      hour = datafile->ra/10000.0;
      datafile->ra -= hour*10000.0;
      min = datafile->ra/100.0;
      datafile->ra -= min*100.0;
      sec = datafile->ra;
      datafile->ra = 2.0*M_PI*((double)(hour+(double)(min+(double)sec/60.0)/60.0)/24.0);
    }else if(strcmp(id, "src_dej") == 0) {
      int deg, min, sign;
      double sec;
      if(fread(&(datafile->dec), sizeof(double), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%lf\n", id, datafile->dec);
      }
      if(datafile->dec < 0) {
	sign = -1;
	datafile->dec = -datafile->dec;
      }else {
	sign = 1;
      }
      deg = datafile->dec/10000.0;
      datafile->dec -= deg*10000.0;
      min = datafile->dec/100.0;
      datafile->dec -= min*100.0;
      sec = datafile->dec;
      datafile->dec = (double)sign*M_PI*((double)(deg+(double)(min+(double)sec/60.0)/60.0))/180.0;
    }else if(strcmp(id, "refdm") == 0) {
      if(fread(&dummy_double, sizeof(double), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%lf\n", id, dummy_double);
      }
      datafile-> dm = dummy_double;
    }else if(strcmp(id, "nbits") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      datafile->NrBits = dummy_int;
      if(datafile->NrBits != 32 && datafile->NrBits != 8) {
	printerror(verbose.debug, "ERROR readSigprocHeader: Can only handle 32-bit or 8-bit data. Got %d bit data.", datafile->NrBits);
	return 0;
      }
    }else if(strcmp(id, "nifs") == 0) {
      if(fread(&nifs, sizeof(int), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%d\n", id, nifs);
      }
    }else if(strcmp(id, "tstart") == 0) {
      if(fread(&dummy_double, sizeof(double), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%lf\n", id, dummy_double);
      }
      datafile->mjd_start = dummy_double;
    }else if(strcmp(id, "tsamp") == 0) {
      if(fread(&dummy_double, sizeof(double), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%lf\n", id, dummy_double);
      }
      datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
      datafile->fixedtsamp = dummy_double;
    }else if(strcmp(id, "npuls") == 0) {
      if(fread(&(dummy_long), sizeof(long int), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%ld\n", id, dummy_long);
      }
      // Not sure about this, maybe it is the number of pulses in subint?
      // Does not appear to be the number of subints
      //      datafile->NrSubints = dummy_long;
    }else if(strcmp(id, "nbins") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%d\n", id, dummy_int);
      }
      datafile->NrBins = dummy_int;
    }else if(strcmp(id, "period") == 0) {
      if(fread(&dummy_double, sizeof(double), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%lf\n", id, dummy_double);
      }
      datafile->isFolded = 1;
      datafile->foldMode = FOLDMODE_FIXEDPERIOD;
      datafile->fixedPeriod = dummy_double;
    }else if(strcmp(id, "barycentric") == 0 || strcmp(id, "pulsarcentric") == 0 || strcmp(id, "nsamples") == 0) {
      if(fread(&dummy_int, sizeof(int), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%d\n", id, dummy_int);
      }
    }else if(strcmp(id, "az_start") == 0 || strcmp(id, "za_start") == 0) {
      if(fread(&dummy_double, sizeof(double), 1, datafile->fptr) != 1) {
	printerror(verbose.debug,"ERROR readSigprocHeader: Read error");
	return 0;
      }
      if(verbose.debug) {
	printf("DEBUG: %s=%lf\n", id, dummy_double);
      }
    }else if(strcmp(id, "HEADER_END") == 0) {
    }else {
      printerror(verbose.debug,"ERROR readSigprocHeader: ID=%s is not recognized", id);
      return 0;
    }
  }while(strcmp(id, "HEADER_END") != 0);
  free(id);

  if(datafile->NrFreqChan <= 0) {
    // Nr of freq channels not set in sigproc when only one channel.
    printwarning(verbose.debug, "WARNING readSigprocHeader: The number of frequency channels do not appear to be defined. It is assumed only one frequency channel is present.");
    datafile->NrFreqChan = 1;
  }

  if(datafile->isFolded) {
    printwarning(verbose.debug, "WARNING readSigprocHeader: For folded data the full baseline is assumed to be stored.");
    datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
    if(datafile->fixedPeriod <= 0) {
      printwarning(verbose.debug, "WARNING readSigprocHeader: Period is not set in the header, assume it is 1 sec.");
      datafile->foldMode = FOLDMODE_FIXEDPERIOD;
      datafile->fixedPeriod = 1;
    }
    datafile->fixedtsamp = datafile->fixedPeriod/(double)datafile->NrBins;
  }

  if(freq_chan1 > 0) {
    datafile->freqMode = FREQMODE_UNIFORM;
    if(datafile->freqlabel_list != NULL) {
      free(datafile->freqlabel_list);
      datafile->freqlabel_list = NULL;
    }
    //    datafile->freq_list = malloc(2*sizeof(double));
    //    if(datafile->freq_list == NULL) {
    //      fflush(stdout);
    //      printerror(verbose.debug, "ERROR readSigprocHeader: Memory allocation error.");
    //      return 0;
    //    }
    if(set_bandwidth(datafile, datafile->NrFreqChan*freq_chanbw, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readSigprocHeader: Bandwidth changing failed.");
      return 0;
    }
    set_centre_frequency(datafile, freq_chan1+0.5*(datafile->NrFreqChan*freq_chanbw), verbose);
  }else {
    printwarning(verbose.debug, "WARNING readSigprocHeader: The frequency of the observation does not appear to be defined.");
    if(datafile->freq_ref > 0) {
      // Freq of obs is not set in sigproc when only one channel. Reference freq is not the same as reference freq though..... But better as nothing, and it is a warning, so ok for now.
      printwarning(verbose.debug, "WARNING readSigprocHeader: It is assumed that the frequency of the observationdefines the reference frequency.");
      datafile->freqMode = FREQMODE_UNIFORM;
      if(datafile->freqlabel_list != NULL) {
	free(datafile->freqlabel_list);
	datafile->freqlabel_list = NULL;
      }
      //      datafile->freq_list = malloc(2*sizeof(double));
      //      if(datafile->freq_list == NULL) {
      //	fflush(stdout);
      //	printerror(verbose.debug, "ERROR readSigprocHeader: Memory allocation error.");
      //	return 0;
      //      }
      if(set_bandwidth(datafile, 0.0, verbose) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readSigprocHeader: Bandwidth changing failed.");
	return 0;
      }
      set_centre_frequency(datafile, datafile->freq_ref, verbose);
    }
  }

  // Find out filesize
  off_t curpos, datasize;
  curpos = ftello(datafile->fptr);
  fseeko(datafile->fptr, 0, SEEK_END);
  datasize = ftello(datafile->fptr);
  datasize -= curpos;
  fseeko(datafile->fptr, curpos, SEEK_SET);
  if(verbose.debug) {
    printf("DEBUG: data size=%ld byte\n", (long)datasize);
  }
  
  if(datafile->isFolded && datafile->NrBins <= 0) {
    printwarning(verbose.debug, "WARNING readSigprocHeader: For folded data the number of bins should be defined in header. Assume the data is not folded");
    datafile->isFolded = 0;
  }

  if(datafile->isFolded == 0) { // Find out number of bins
    long double tmp;
    tmp = datasize*8.0;                     // Size in bits
    tmp /= (long double)(datafile->NrBits); // Nr of samples
    tmp /= (long double)nifs;               // Per IF
    tmp /= (long double)datafile->NrFreqChan;   // Actual number of samples
    datafile->NrBins = tmp;
    if(verbose.debug) {
      printf("DEBUG: Implied number of bins = %ld\n", datafile->NrBins);
    }    
    // Assume there is only one subint
    datafile->NrSubints = 1;
  }else {
    //    datasize + headersize = nrsubints*(nfreq*nbin*nrbits/8+headersize)
    long double tmp;
    long double subintsize;
    subintsize = ((long double)(datafile->NrFreqChan)*(long double)(datafile->NrBins)*(long double)(datafile->NrBits)/8.0 + (long double)curpos);
    tmp = datasize + curpos;
    tmp /= subintsize;
    datafile->NrSubints = roundl(tmp);
  }

  // Maybe need to multiply with npuls, although that might be different for different subints
  datafile->tsub_list[0] = datafile->NrBins * datafile->fixedtsamp;

  /*
  datafile->rm = 0;
  s_ptr = get_ptr_entry("Reference frequency", txt, nrlines, ":");
  if(s_ptr != NULL) {
    sscanf(s_ptr, "%lf", &(datafile->freq_ref));
  }
  */

  printwarning(verbose.debug, "WARNING readSigprocHeader: Assuming there is only one polarization channel in the data");
  datafile->NrPols = 1;

  return 1;
}
//START REGION DEVELOP


//START REGION RELEASE
int readSigprocfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  //  char txt[10000];
  long n, f, i, p;
  float *sample_f;
  unsigned char *sample_b;
  int ret;

  if(datafile.NrPols > 1) {
    printerror(verbose.debug, "readSigprocfile: Data should have just one polarization");
    return 0;
  }
  //  if(datafile.NrSubints > 1) {
  //    printerror(verbose.debug, "readSigprocfile: Data should have just one subint");
  //    return 0;
  //  }
  //  if(datafile.NrFreqChan > 1) {
  //    printerror(verbose.debug, "readSigprocfile: Data should have just one frequency channel");
  //    return 0;
  //  }
  if(datafile.NrBits != 32 && datafile.NrBits != 8) {
    printerror(verbose.debug, "ERROR readSigprocfile: Can only handle 32-bit or 8-bit data. Got %d bit data.", datafile.NrBits);
    return 0;
  }
  
  if(datafile.NrBits == 32) {
    sample_f = malloc(datafile.NrFreqChan*sizeof(float));
    if(sample_f == NULL) {
      printerror(verbose.debug, "readSigprocfile: Memory allocation error");
      return 0;
    }
  }else if(datafile.NrBits == 8) {
    sample_b = malloc(datafile.NrFreqChan);
    if(sample_b == NULL) {
      printerror(verbose.debug, "readSigprocfile: Memory allocation error");
      return 0;
    }
  }
  for(n = 0; n < datafile.NrSubints; n++) {
    //    if(verbose.verbose) printf("readSigprocfile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
    for(i = 0; i < datafile.NrBins; i++) {
      if(datafile.NrBits == 32) {
	ret = fread(sample_f, sizeof(float), datafile.NrFreqChan, datafile.fptr);
      }else if(datafile.NrBits == 8) {
	ret = fread(sample_b, 1, datafile.NrFreqChan, datafile.fptr);
      }
      if(ret != datafile.NrFreqChan) {
	printerror(verbose.debug, "ERROR readSigprocfile: Cannot read data (sample %ld of subint %ld).", i, n);
	return 0;
      }
      p = 0;
      for(f = 0; f < datafile.NrFreqChan; f++) {
	if(datafile.NrBits == 32) {
	  data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] = sample_f[f];
	}else if(datafile.NrBits == 8) {
	  data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] = sample_b[f];
	}
      }
    }
    if(n != datafile.NrSubints - 1) { // Skip next header, as each subint has a header
      fseeko(datafile.fptr, datafile.datastart, SEEK_CUR);
    }
  }
  if(datafile.NrBits == 32) {
    free(sample_f);
  }else if(datafile.NrBits == 8) {
    free(sample_b);
  }
  if(verbose.verbose) printf("Reading is done.                           \n");
  return 1;
}

// Return 1 on success
int readPulseSigprocData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  long i;
  float *sample_f;
  unsigned char *sample_b;
  int ret;
  off_t offset;

  if(datafile.NrPols > 1) {
    printerror(verbose.debug, "readPulseSigprocData: Data should have just one polarization");
    return 0;
  }
  if(datafile.NrBits != 32 && datafile.NrBits != 8) {
    printerror(verbose.debug, "ERROR readPulseSigprocData: Can only handle 32-bit or 8-bit data. Got %d bit data.", datafile.NrBits);
    return 0;
  }
  
  if(datafile.NrBits == 32) {
    sample_f = malloc(datafile.NrFreqChan*sizeof(float));
    if(sample_f == NULL) {
      printerror(verbose.debug, "readPulseSigprocData: Memory allocation error");
      return 0;
    }
  }else if(datafile.NrBits == 8) {
    sample_b = malloc(datafile.NrFreqChan);
    if(sample_b == NULL) {
      printerror(verbose.debug, "readPulseSigprocData: Memory allocation error");
      return 0;
    }
  }

  offset = (pulsenr*datafile.NrBins+binnr)*datafile.NrFreqChan+freq; // Nr of samples to skip
  if(datafile.NrBits == 32) {           
    offset *= 4;
  }
  offset += (pulsenr+1)*datafile.datastart;  // Each subint has a header
  fseeko(datafile.fptr, offset, SEEK_SET);

  for(i = 0; i < nrSamples; i++) {
    if(datafile.NrBits == 32) {
      ret = fread(sample_f, sizeof(float), 1, datafile.fptr);
      pulse[i] = *sample_f;
    }else if(datafile.NrBits == 8) {
      ret = fread(sample_b, 1, 1, datafile.fptr);
      pulse[i] = *sample_b;
    }
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readPulseSigprocData: Cannot read data.");
      return 0;
    }
    if(datafile.NrFreqChan > 1) {
      if(i < nrSamples-1) {
	if(datafile.NrBits == 32) {
	  fseeko(datafile.fptr, (datafile.NrFreqChan-1)*sizeof(float), SEEK_CUR);
	}else {
	  fseeko(datafile.fptr, (datafile.NrFreqChan-1), SEEK_CUR);
	}
      }
    }
  }
  if(datafile.NrBits == 32) {
    free(sample_f);
  }else if(datafile.NrBits == 8) {
    free(sample_b);
  }
  return 1;
}

//START REGION DEVELOP

//START REGION RELEASE
int readSigprocASCIIHeader(datafile_definition *datafile, verbose_definition verbose)
{

  /* # MJD_INT SEC_AFTER_MIDNIGHT FOLDPERIOD 1 BOTTOMFREQ DM NRBINS SITE???? 1 SRCNAME */

  int j, c;
  double sec, freq;
  long double mjd;
  char tmp1[10000];
  char tmp2[1000];
  char tmp3[1000];
  char tmp4[1000];
  char tmp5[1000];
  char *txtptr;
  long pos;
  j = fscanf(datafile->fptr_hdr, "%s %Lf %lf %lf %s %lf %lf %ld %s %s %s", tmp1, &mjd, &sec, &(datafile->fixedPeriod), tmp2, &freq,  &(datafile->dm), &(datafile->NrBins), tmp4, tmp3, tmp5);
  if(set_psrname_PSRData(datafile, tmp5, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readSigprocASCIIHeader: Setting pulsar name failed.");
    return 0;
  }

  if(j != 11) {
    printerror(verbose.debug,"ERROR readSigprocASCIIHeader: Error reading first line (%d != 11).", j);
    return 0;
  }
  datafile->isFolded = 1;
  datafile->foldMode = FOLDMODE_FIXEDPERIOD;
  if(set_observatory_PSRData(datafile, tmp4, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readSigprocASCIIHeader: Setting observatory name failed.");
    return 0;
  }
  if(strcmp(tmp1, "#") != 0) {
    printerror(verbose.debug, "readSigprocASCIIHeader: I do not understand header (%s != #).", tmp1);
    return 0;
  }
  if(strcmp(tmp2, "1") != 0) {
    printerror(verbose.debug, "readSigprocASCIIHeader: I do not understand header (%s != 1).", tmp2);
    return 0;
  }
  if(strcmp(tmp3, "1") != 0) {
    printerror(verbose.debug, "readSigprocASCIIHeader: I do not understand header (%s != 1).", tmp3);
    return 0;
  }
//START REGION DEVELOP
  if(ingridflavour) {   /* Read in and ignore the reference phase */
    j = fscanf(datafile->fptr_hdr, "%s", tmp1);
    if(j != 1) {
      printerror(verbose.debug,"ERROR readSigprocASCIIHeader: Error reading first line (%d != 1).", j);
      exit(-1);
    }
  }
//START REGION RELEASE
  datafile->mjd_start = mjd;
  datafile->mjd_start += sec/(double)(24.0*60.0*60.0);
  datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
  datafile->fixedtsamp = 0;
  datafile->tsubMode = TSUBMODE_FIXEDTSUB;
  if(datafile->tsub_list != NULL)
    free(datafile->tsub_list);
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    printerror(verbose.debug,"ERROR readSigprocASCIIHeader: Memory allocation error");
    exit(-1);
  }
  datafile->tsub_list[0] = 0;

  /* Not sure this is really supposed to be the centre frequency. I
     thought the header program of sigproc claims it is the frequency
     of the first channel. Don't know what the bandwidth is, so take
     central frequency. */
  datafile->freqMode = FREQMODE_UNIFORM;
  if(datafile->freqlabel_list != NULL) {
    free(datafile->freqlabel_list);
    datafile->freqlabel_list = NULL;
  }
  //  datafile->freq_list = malloc(2*sizeof(double));
  //  if(datafile->freq_list == NULL) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "ERROR readSigprocASCIIHeader: Memory allocation error.");
  //    return 0;
  //  }
  if(set_bandwidth(datafile, 0.0, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readSigprocASCIIHeader: Bandwidth changing failed.");
    return 0;
  }
  set_centre_frequency(datafile, freq, verbose);

  printwarning(verbose.debug, "WARNING readSigprocASCIIHeader: Assuming there is only one polarization channel in the data");
  datafile->NrPols = 1;

  pos = ftell(datafile->fptr_hdr);
  datafile->NrSubints = 0; 
  do {
    c = fgetc(datafile->fptr_hdr);
    if(c == '\n')
      datafile->NrSubints++;
  }while(c != EOF);
  datafile->NrSubints /= (datafile->NrBins + 1);  /* There is 1 header line per pulse */
  fseek(datafile->fptr_hdr, pos, SEEK_SET);

  fgets(tmp1, 9999, datafile->fptr_hdr);
  if(strlen(tmp1) < 3)
    fgets(tmp1, 9999, datafile->fptr_hdr);
  datafile->NrFreqChan = 0;
  do {
    if(datafile->NrFreqChan == 0) {
      txtptr = strtok(tmp1, " ");
      if(
//START REGION DEVELOP
	 ingridflavour == 0 && 
//START REGION RELEASE
	 strcmp(txtptr, "1")) {
	printerror(verbose.debug, "ERROR readSigprocASCIIHeader: Expected different bin number ('%s' != %d).", txtptr, 1);
	return 0;
//START REGION DEVELOP
      }else if(ingridflavour != 0 && strcmp(txtptr, "0")) {
	printerror(verbose.debug, "ERROR readSigprocASCIIfile: Expected different bin number ('%s' != %d).", txtptr, 0);
	return 0;
//START REGION RELEASE
      }
      txtptr = strtok(NULL, " ");
    }else {
      txtptr = strtok(NULL, " ");
    }
    if(txtptr != NULL) {
      datafile->NrFreqChan ++;
    }
  }while(txtptr != NULL);
  fseek(datafile->fptr_hdr, pos, SEEK_SET);

  return 1;
}
//START REGION DEVELOP

void mpolyco_t1(char  *fname,char  *psrname,int imjd,double fmjd,char  *nsite,int *nspan,int *ncoeff, char *par_directory);
void ppolyco(char *unfname,int imjd,double frmjd,double *pobs,double *phobs);

//START REGION RELEASE
int writeSigprocASCIIHeader(datafile_definition datafile, verbose_definition verbose)
{
  /* # MJD_INT SEC_AFTER_MIDNIGHT FOLDPERIOD 1 BOTTOMFREQ DM NRBINS SITE???? 1 SRCNAME */
  int int_mjd;
  double sec, freq;
//START REGION DEVELOP
  int nspan, ncoeff;
  double pobs, refph;
  char site[1];
  long i;
//START REGION RELEASE

  if(datafile.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeSigprocASCIIHeader: Writing out this data format is only implemented when the frequency channels are uniformely separated.");
    return 0;
  }

  int_mjd = datafile.mjd_start;
  sec = datafile.mjd_start;
  sec -= int_mjd;
  sec *= (24.0*60.0*60.0);
  /* Not sure this is really supposed to be the centre frequency. I
     thought the header program of sigproc claims it is the frequency
     of the first channel, but BW is not written out to header, so I'm
     not sure. I take central frequency. */
  freq = get_centre_frequency(datafile, verbose);
  
  fprintf(datafile.fptr_hdr, "# %d %lf %f %d %lf %f %ld %s %d %s", int_mjd, sec, datafile.fixedPeriod, 1, freq,  datafile.dm, datafile.NrBins, datafile.observatory, 1, datafile.psrname);
//START REGION DEVELOP
  if(ingridflavour) {
    /* Must first generate the polyco at the midpoint */
    sprintf(site,"i");  
    extern int verb2;   /* Enable verbose mode in the mpolyco function */
    verb2 = 1;
    mpolyco_t1(datafile.filename,datafile.psrname,int_mjd,datafile.mjd_start-int_mjd,site,&nspan,&ncoeff, ingridflavour_parfile);
    ppolyco(datafile.filename,int_mjd,datafile.mjd_start-int_mjd,&pobs,&refph);
    i = refph;
    refph -= i;
    if(refph < -1.0)   /* In case something goes wrong with rounding negative numbers. */
      refph += 1.0;
    printwarning(verbose.debug, "WARNING: ASSUMING SITE IS WSRT AND THAT DATA IS REDUCED WITH fpuma.");
    
    fprintf(datafile.fptr_hdr, " %lf\n", refph);
  }else {
//START REGION RELEASE
    fprintf(datafile.fptr_hdr, "\n");
//START REGION DEVELOP
  }
//START REGION RELEASE
  return 1;
}
//START REGION DEVELOP

//START REGION RELEASE
int readSigprocASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  char txt[10000];
  long n, f, i, p;
  float sample;
  int ret;

  if(datafile.NrPols > 1) {
    printerror(verbose.debug, "readSigprocASCIIfile: Data should have just one polarization");
    return 0;
  }
  for(n = 0; n < datafile.NrSubints; n++) {
    if(verbose.verbose && verbose.nocounters == 0) printf("readSigprocASCIIfile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
    for(i = 0; i < datafile.NrBins; i++) {
      for(f = 0; f < datafile.NrFreqChan; f++) {
	ret = fscanf(datafile.fptr, "%s", txt);
	if(ret != 1) {
	  printerror(verbose.debug, "ERROR readSigprocASCIIfile: Cannot read data.");
	  return 0;
	}
	if(f == 0) {
	  if(
//START REGION DEVELOP
	     ingridflavour == 0 && 
//START REGION RELEASE
	     atoi(txt) != i+1) {
	    printerror(verbose.debug, "ERROR readSigprocASCIIfile: Expected different bin number ('%s' != %ld).", txt, i+1);
	    return 0;
//START REGION DEVELOP
	  }else if(ingridflavour != 0 && atoi(txt) != i) {
	    printerror(verbose.debug, "ERROR readSigprocASCIIfile: Expected different bin number ('%s' != %ld).", txt, i);
	    return 0;
//START REGION RELEASE
	  }
	  ret = fscanf(datafile.fptr, "%s", txt);
	  if(ret != 1) {
	    printerror(verbose.debug, "ERROR readSigprocASCIIfile: Cannot read data.");
	    return 0;
	  }
	}
	sscanf(txt, "%f", &sample);
	p = 0;
	data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] = sample;
      }
    }
    /* Skip header */
    fgets(txt, 9999, datafile.fptr_hdr);
    if(strlen(txt) < 2)
      fgets(txt, 9999, datafile.fptr_hdr);
  }
  if(verbose.verbose && verbose.nocounters == 0) printf("Reading is done.                           \n");
  return 1;
}
//START REGION DEVELOP

//START REGION RELEASE
int writeSigprocASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long n, f, i, p;
  float sample;

  if(datafile.NrPols > 1) {
    printerror(verbose.debug, "writeSigprocASCIIfile: Data should have just one polarization");
    return 0;
  }
  p = 0;
  for(n = 0; n < datafile.NrSubints; n++) {
    if(verbose.verbose && verbose.nocounters == 0) printf("writeSigprocASCIIfile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
    for(i = 0; i < datafile.NrBins; i++) {
      for(f = 0; f < datafile.NrFreqChan; f++) {
	if(f == 0) {
//START REGION DEVELOP
	  if(ingridflavour == 0)
//START REGION RELEASE
	    fprintf(datafile.fptr, "%5ld", i+1);
//START REGION DEVELOP
	  else
	    fprintf(datafile.fptr, "%5ld", i);
//START REGION RELEASE
	}
	sample = data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i]; 
	fprintf(datafile.fptr, " %f", sample);
      }
      fprintf(datafile.fptr, "\n");
    }
    if(n != datafile.NrSubints - 1) {
      if(writeSigprocASCIIHeader(datafile, verbose) == 0) {
	printerror(verbose.debug, "writeSigprocASCIIfile: Writing header line failed.");
	return 0;
      }
    }
  }
  if(verbose.verbose && verbose.nocounters == 0) printf("Writing is done.                           \n");
  return 1;
}
//START REGION DEVELOP
