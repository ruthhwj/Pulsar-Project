#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <math.h>
#include <string.h>
#include "psrsalsa.h"


int readGMRTasciiHeader(datafile_definition *datafile, verbose_definition verbose)
{
  double dummy_d;
  int ret;

  ret = fscanf(datafile->fptr,"%ld %ld %lf %lf %lf %lf", &(datafile->NrSubints), &(datafile->NrBins), &(datafile->fixedtsamp), &(datafile->fixedPeriod), &dummy_d, &dummy_d);
  if(ret != 6) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readGMRTasciiHeader: Cannot read header! (nrarguments=%d)", ret);
    return 0;
  }
  datafile->isFolded = 1;
  datafile->foldMode = FOLDMODE_FIXEDPERIOD;
  datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
  datafile->NrFreqChan = 1;
  datafile->NrPols = 1;
  datafile->fixedtsamp *= 0.001; // convert to seconds
  datafile->fixedPeriod *= 0.001; // convert to seconds
  datafile->tsubMode = TSUBMODE_FIXEDTSUB;
  if(datafile->tsub_list != NULL)
    free(datafile->tsub_list);
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readGMRTasciiHeader: Memory allocation error");
    return 0;
  }
  datafile->tsub_list[0] = 0;
  datafile->freqMode = FREQMODE_UNIFORM;
  if(datafile->freqlabel_list != NULL) {
    free(datafile->freqlabel_list);
    datafile->freqlabel_list = NULL;
  }
  //  datafile->freq_list = malloc(2*sizeof(double));
  //  if(datafile->freq_list == NULL) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "ERROR readGMRTasciiHeader: Memory allocation error");
  //    return 0;
  //  }
  set_centre_frequency(datafile, 0.0, verbose);
  if(set_bandwidth(datafile, 0.0, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readGMRTasciiHeader: Bandwidth changing failed.");
    return 0;
  }
  if(set_psrname_PSRData(datafile, "Unknown", verbose) == 0) { // Set with something, allowing DESH ascii format to work
    fflush(stdout);
    printerror(verbose.debug, "ERROR readGMRTasciiHeader: Setting pulsar name failed.");
    return 0;
  }
  if(set_institute_PSRData(datafile, "Unknown", verbose) == 0) { // Set with something, allowing DESH ascii format to work
    fflush(stdout);
    printerror(verbose.debug, "ERROR readGMRTasciiHeader: Setting pulsar name failed.");
    return 0;
  }
  return 1;
}



int readGMRTasciifile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long p, n, i, f, ret, pulsenr, binnr;
  float sample;
  double dummy_d;

  if(datafile.NrFreqChan != 1) {
    printerror(verbose.debug, "ERROR readGMRTasciifile: Only one frequency channel is supported.");
    return 0;
  }
  if(datafile.NrPols != 1) {
    printerror(verbose.debug, "ERROR readGMRTasciifile: Only one polarization channel is supported.");
    return 0;
  }

  if(verbose.verbose) printf("Reading file in GMRT folded ascii format.\n");
  f = 0;
  p = 0;
  fseek(datafile.fptr, datafile.datastart, SEEK_SET);
  for(n = 0; n < datafile.NrSubints; n++) {
    if(verbose.verbose) printf("  readGMRTasciifile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
    ret = fscanf(datafile.fptr,"%ld %lf %lf %lf %lf %lf", &(pulsenr), &dummy_d, &dummy_d, &dummy_d, &dummy_d, &dummy_d);
    if(ret != 6) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readGMRTasciifile: Cannot read subint header! (nrarguments=%ld)", ret);
      return 0;
    }
    if(pulsenr != n+1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readGMRTasciifile: Subint header for subint %ld does start with a different number (%ld)", n+1, pulsenr);
      return 0;      
    }

    for(i = 0; i < datafile.NrBins; i++) {
      ret = fscanf(datafile.fptr, "%ld %ld", &pulsenr, &binnr);
      if(ret != 2) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readGMRTasciifile: Cannot read bin header! (nrarguments=%ld, subint=%ld, bin=%ld)", ret, n+1, i+1);
	return 0;
      }
      if(pulsenr != n+1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readGMRTasciifile: Bin header for subint %ld does start with a different number (%ld)", n+1, pulsenr);
	return 0;      
      }
      if(binnr != i+1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readGMRTasciifile: Bin header for bin %ld does start with a different number (%ld) for subint %ld", i+1, binnr, n+1);
	return 0;      
      }
      ret = fscanf(datafile.fptr, "%f", &sample);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readGMRTasciifile: Cannot read intensity sample! (nrarguments=%ld)", ret);
	return 0;
      }
      data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] = sample;
    }
  }
  if(verbose.verbose) printf("  Reading is done.                  \n");
  return 1;
}
