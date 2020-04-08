#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include "psrsalsa.h"

int readParkesHeader(datafile_definition *datafile, verbose_definition verbose)
{
  long i, h, m;
  double s, channelbw;
  char txt[1000];

  /* Parkes data is a timeserie */
  //  datafile->dd_mode = 1;
  datafile->gentype = GENTYPE_SEARCHMODE;
  datafile->isFolded = 0;
  datafile->foldMode = FOLDMODE_UNKNOWN;
  datafile->fixedPeriod = -1;

  fscanf(datafile->fptr_hdr, "%s", txt);
  fscanf(datafile->fptr_hdr, "%s", txt);
  fscanf(datafile->fptr_hdr, "%s", txt);
  fscanf(datafile->fptr_hdr, "%ld", &i);
  fscanf(datafile->fptr_hdr, "%ld:%ld:%lf", &h, &m, &s);
  datafile->mjd_start += i+ (double)h/24.0+(double)m/1440.0 + s/86400.0;
  fscanf(datafile->fptr_hdr, "%s", txt);
  fscanf(datafile->fptr_hdr, "%s", txt);
  fscanf(datafile->fptr_hdr, "%s", txt);
  fscanf(datafile->fptr_hdr, "%s", txt);
  fscanf(datafile->fptr_hdr, "%s", txt);
  fscanf(datafile->fptr_hdr, "%s", txt);
  fscanf(datafile->fptr_hdr, "%s", txt);
  fseek(datafile->fptr_hdr, 208, SEEK_SET);
  fscanf(datafile->fptr_hdr, "%lf", &channelbw);
  fscanf(datafile->fptr_hdr, "%s", txt);
  fscanf(datafile->fptr_hdr, "%ld", &(datafile->NrFreqChan));
  fscanf(datafile->fptr_hdr, "%s", txt);
  fscanf(datafile->fptr_hdr, "%s", txt);
  datafile->freqMode = FREQMODE_UNIFORM;
  //  if(datafile->freq_list != NULL) {
  //    free(datafile->freq_list);
  //  }
  //  datafile->freq_list = malloc(2*sizeof(double));
  //  if(datafile->freq_list == NULL) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "ERROR readParkesHeader: Memory allocation error.");
  //    return 0;
  //  }
  double freq;
  fscanf(datafile->fptr_hdr, "%lf", &freq);
  freq = (freq-0.5*channelbw)+(0.5*channelbw*(datafile->NrFreqChan));
  set_centre_frequency(datafile, freq, verbose);
  if(set_bandwidth(datafile, channelbw*(datafile->NrFreqChan), verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readParkesHeader: Bandwidth changing failed.");
    return 0;
  }
  datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
  fscanf(datafile->fptr_hdr, "%lf", &(datafile->fixedtsamp));
  datafile->fixedtsamp *= 0.001;
  fseek(datafile->fptr_hdr, 490, SEEK_SET);
  fscanf(datafile->fptr_hdr, "%ld", &i);
  if(verbose.verbose) printf("There are %ld blocks of 8192 samples in the data.\n", i);
  /* i = Nr blocks */
  datafile->NrBins = i*8192;
  datafile->NrBits = 16;
  /*
  i = 87;
  printf("i=%ld\nbytes=%ld\n", i, i*8192*2+i*8);
  */
  fseek(datafile->fptr_hdr, 528, SEEK_SET);
  fscanf(datafile->fptr_hdr, "%s", txt);
  if(set_psrname_PSRData(datafile, txt, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readParkesHeader: Setting pulsar name failed.");
    return 0;
  }

  datafile->NrPols = 1;
  
  /* Simons dedispersed timeserie is frequency scrunched */
  /* Patrick: Disable for now to try to read in sc_td output. */
  /*  datafile->NrFreqChan = 1;  */
  
  /*
  fscanf(datafile->fptr_hdr, "%s", txt);
  printf("txt=%s\n", txt);
  */

  return 1;
}

