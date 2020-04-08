//START REGION RELEASE
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include "psrsalsa.h"

//START REGION DEVELOP
#include "gsl/gsl_rng.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_randist.h"

//START REGION RELEASE
int readWSRTHeader(datafile_definition *datafile, verbose_definition verbose);
int readPulseWSRTData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse);
int writePulseWSRTData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse);
int writeWSRTHeader(datafile_definition datafile, verbose_definition verbose);
int writePuMafile(datafile_definition datafile, float *data, verbose_definition verbose);
int readPuMafile(datafile_definition datafile, float *data, verbose_definition verbose);
int readPSRFITSHeader(datafile_definition *datafile, int readnoscales, int nowarnings, verbose_definition verbose);
int readFITSpulse(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
int readFITSfile(datafile_definition *datafile, float *data, verbose_definition verbose);
int writePSRFITSHeader(datafile_definition *datafile, verbose_definition verbose);
int writeFITSpulse(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
int writeFITSfile(datafile_definition datafile, float *data, verbose_definition verbose);
int readPSRCHIVE_ASCIIHeader(datafile_definition *datafile, verbose_definition verbose);
int readPSRCHIVE_ASCIIfilepulse(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
int readPSRCHIVE_ASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose);
int writePSRCHIVE_ASCIIHeader(datafile_definition datafile, verbose_definition verbose);
int writePSRCHIVE_ASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose);
//START REGION DEVELOP
int writeAOAsciiHeader(datafile_definition datafile, int type, verbose_definition verbose);
int readAOAsciiHeader(datafile_definition *datafile, int type, verbose_definition verbose);
int writeAOAsciifile(datafile_definition datafile, float *data, int type, verbose_definition verbose);
int readAOAsciifile(datafile_definition datafile, float *data, int type, verbose_definition verbose);
//START REGION RELEASE
int readEPNHeader(datafile_definition *datafile, int what, verbose_definition verbose);
int readPulseEPNData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
int writeEPNfile(datafile_definition datafile, float *data, verbose_definition verbose);
int readEPNfile(datafile_definition *datafile, float *data, verbose_definition verbose, long request_only_one_pulse);
int readEPNsubHeader(datafile_definition *datafile, float *scale, float *offset, verbose_definition verbose);
//START REGION DEVELOP
int readParkesHeader(datafile_definition *datafile, verbose_definition verbose);
int writePRESTOHeader(datafile_definition datafile, verbose_definition verbose);
int readPRESTOHeader(datafile_definition *datafile, verbose_definition verbose);
//START REGION RELEASE
int readSigprocHeader(datafile_definition *datafile, verbose_definition verbose);
// Return 1 on success
int readPulseSigprocData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose);
//START REGION DEVELOP

//START REGION RELEASE
int readPPOLHeader(datafile_definition *datafile, int extended, verbose_definition verbose);
int writePPOLHeader(datafile_definition datafile, int argc, char **argv, verbose_definition verbose);
int readHistoryFITS(datafile_definition *datafile, verbose_definition verbose);
int writeHistoryFITS(datafile_definition datafile, verbose_definition verbose);
int readHistoryPSRData(datafile_definition *datafile, verbose_definition verbose);
int writeHistoryPSRData(datafile_definition *datafile, int argc, char **argv, int cmdOnly, verbose_definition verbose);
int writeHistoryPuma(datafile_definition datafile, verbose_definition verbose);
int readHistoryPuma(datafile_definition *datafile, verbose_definition verbose);
int readSigprocfile(datafile_definition datafile, float *data, verbose_definition verbose);
int readSigprocASCIIHeader(datafile_definition *datafile, verbose_definition verbose);
int writeSigprocASCIIHeader(datafile_definition datafile, verbose_definition verbose);
int writeSigprocASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose);
int readSigprocASCIIfile(datafile_definition datafile, float *data, verbose_definition verbose);
//START REGION DEVELOP
int readGMRTasciiHeader(datafile_definition *datafile, verbose_definition verbose);
int readGMRTasciifile(datafile_definition datafile, float *data, verbose_definition verbose);

//START REGION RELEASE

/*
  Return 1 if data format is recognized, or zero otherwise.
*/
int isValidPSRDATA_format(int format) 
{
  if(format == PUMA_format)
    return 1;
//START REGION DEVELOP
  if(format == AO_ASCII_1_format)
    return 1;
  if(format == AO_ASCII_2_format)
    return 1;
  if(format == PRESTO_format)
    return 1;
  if(format == PARKESFB_format)
    return 1;
//START REGION RELEASE
  if(format == PSRCHIVE_ASCII_format)
    return 1;
  if(format == EPN_format)
    return 1;
  if(format == FITS_format)
    return 1;
  if(format == SIGPROC_format)
    return 1;
//START REGION RELEASE
  if(format == PPOL_format)
    return 1;
  if(format == PPOL_SHORT_format)
    return 1;
  if(format == SIGPROC_ASCII_format)
    return 1;
//START REGION DEVELOP
  if(format == GMRT_ASCII_format)
    return 1;
//START REGION RELEASE
  if(format == PSRSALSA_BINARY_format)
    return 1;
  if(format == MEMORY_format)
    return 1;
  printerror(0, "ERROR isValidPSRDATA_format: specified data format is not recognized.");
  return 0;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Print the available data formats to device (could be for instance
   stdio). nrspaces defines the number of spaces before each line.*/
void printPSRDataFormats(FILE *printdevice, int nrspaces)
{
  int i, nrspaces2;
  nrspaces2 = nrspaces + 17; 
//START REGION DEVELOP
  nrspaces2 += 3;
//START REGION RELEASE
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
//START REGION DEVELOP
  fprintf(printdevice, "%2d ", PUMA_format);
//START REGION RELEASE
  fprintf(printdevice, "(PUMA)         - WSRT PuMa format\n");
//START REGION DEVELOP
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "%2d (AOASCII1)     - Folded Arecibo ASCII data\n", AO_ASCII_1_format);
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "%2d (AOASCII2)     - Another type of folded Arecibo ASCII data\n", AO_ASCII_2_format);
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "%2d (PRESTO)       - PRESTO timeseries\n", PRESTO_format);
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "%2d (PARKESFB)     - Parkes analogue fb timeseries (Simon's format)\n", PARKESFB_format);
//START REGION RELEASE
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
//START REGION DEVELOP
  fprintf(printdevice, "%2d ", PSRCHIVE_ASCII_format);
//START REGION RELEASE
  fprintf(printdevice, "(ASCII)        - PSRCHIVE ascii dump file (generated by for instance\n");
  for(i = 0; i < nrspaces2; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "\"pdv -t\"). File has limited header information.\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
//START REGION DEVELOP
  fprintf(printdevice, "%2d ", EPN_format);
//START REGION RELEASE
  fprintf(printdevice, "(EPN)          - EPN format\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
//START REGION DEVELOP
  fprintf(printdevice, "%2d ", FITS_format);
//START REGION RELEASE
  fprintf(printdevice, "(PSRFITS)      - PSRFITS format (generated by for instance\n");
  for(i = 0; i < nrspaces2; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "\"pam -a PSRFITS\"). Note that the data files written out\n");
  for(i = 0; i < nrspaces2; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "do not 100%% conform with the PSRFITS definition, so\n");
  for(i = 0; i < nrspaces2; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "proper behavior in other software cannot be guaranteed.\n");
  for(i = 0; i < nrspaces2; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "Especially timing experiments are not recommended.\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
//START REGION DEVELOP
  fprintf(printdevice, "%2d ", SIGPROC_format);
//START REGION RELEASE
  fprintf(printdevice, "(SIGPROC)      - SIGPROC binary format\n");
//START REGION RELEASE
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
//START REGION DEVELOP
  fprintf(printdevice, "%2d ", PPOL_format);
//START REGION RELEASE
  fprintf(printdevice, "(PPOL)         - PPOL format\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
//START REGION DEVELOP
  fprintf(printdevice, "%2d ", PPOL_SHORT_format);
//START REGION RELEASE
  fprintf(printdevice, "(PPOLSHORT)    - PPOL SHORT format (longitude, pa, pa error)\n");
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
//START REGION DEVELOP
  fprintf(printdevice, "%2d ", SIGPROC_ASCII_format);
//START REGION RELEASE
  fprintf(printdevice, "(SIGPROCASCII) - Sigproc ascii format\n");
//START REGION DEVELOP
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
  fprintf(printdevice, "%2d (GMRT)         - GMRT ASCII folded data\n", GMRT_ASCII_format);
//START REGION RELEASE
  for(i = 0; i < nrspaces; i++) fprintf(printdevice, " ");
//START REGION DEVELOP
  fprintf(printdevice, "%2d ", PSRSALSA_BINARY_format);
//START REGION RELEASE
  fprintf(printdevice, "(PSRSALSA)     - PSRSALSA binary format\n");
}
//START REGION DEVELOP

//START REGION RELEASE
/* Parse a string to find the data format for the -oformat/-iformat
   command line option. Returns 0 on error, otherwise the data
   format. */
int parsePSRDataFormats(char *cmd)
{
  if(strcasecmp(cmd, "ASCII") == 0 || atoi(cmd) == PSRCHIVE_ASCII_format)
    return PSRCHIVE_ASCII_format;
  else if(strcasecmp(cmd, "PSRFITS") == 0 || strcasecmp(cmd, "FITS") == 0 || atoi(cmd) == FITS_format)
    return FITS_format;
  else if(strcasecmp(cmd, "PSRSALSA") == 0 || strcasecmp(cmd, "SALSA") == 0 || atoi(cmd) == PSRSALSA_BINARY_format)
    return PSRSALSA_BINARY_format;
  else if(strcasecmp(cmd, "PUMA") == 0 || atoi(cmd) == PUMA_format)
    return PUMA_format;
//START REGION DEVELOP
  else if(strcasecmp(cmd, "DESH") == 0 || strcasecmp(cmd, "AOASCII") == 0 || strcasecmp(cmd, "ASCIIAO") == 0 || strcasecmp(cmd, "AOASCII1") == 0 || strcasecmp(cmd, "ASCIIAO1") == 0 || atoi(cmd) == AO_ASCII_1_format)
    return AO_ASCII_1_format;
  else if(strcasecmp(cmd, "AOASCII2") == 0 || strcasecmp(cmd, "ASCIIAO2") == 0 || atoi(cmd) == AO_ASCII_2_format)
    return AO_ASCII_2_format;
  else if(strcasecmp(cmd, "PRESTO") == 0 || atoi(cmd) == PRESTO_format)
    return PRESTO_format;
  else if(strcasecmp(cmd, "PARKESFB") == 0 || atoi(cmd) == PARKESFB_format)
    return PARKESFB_format;
//START REGION RELEASE
  else if(strcasecmp(cmd, "EPN") == 0 || atoi(cmd) == EPN_format)
    return EPN_format;
  else if(strcasecmp(cmd, "SIGPROC") == 0 || atoi(cmd) == SIGPROC_format)
    return SIGPROC_format;
  else if(strcasecmp(cmd, "PPOL") == 0 || strcasecmp(cmd, "PASWING") == 0 || atoi(cmd) == PPOL_format)
    return PPOL_format;
  else if(strcasecmp(cmd, "PPOLSHORT") == 0 || strcasecmp(cmd, "PPOL_SHORT") == 0 || strcasecmp(cmd, "PASWINGSHORT") == 0 || atoi(cmd) == PPOL_SHORT_format)
    return PPOL_SHORT_format;
  else if(strcasecmp(cmd, "SIGPROCASCII") == 0 || strcasecmp(cmd, "SIGPROC_ASCII") == 0 || atoi(cmd) == SIGPROC_ASCII_format)
    return SIGPROC_ASCII_format;
//START REGION DEVELOP
  else if(strcasecmp(cmd, "GMRT") == 0 || atoi(cmd) == GMRT_ASCII_format)
    return GMRT_ASCII_format;
//START REGION RELEASE
  else {
    fflush(stdout);
    printerror(0, "parsePSRDataFormats: Cannot parse '%s' as a valid data format", cmd);
    return 0;
  }
  return 0;
}

//START REGION DEVELOP
/*
void addquotesargv(char *argv, char *txt)
{
  int i;
  char txt1[1005], txt2[1005];
  i = sscanf(argv, "%s %s", txt1, txt2);
  if(i == 2) 
    sprintf(txt, "'%s'", argv);
  else
    sprintf(txt, "%s", argv);
}
*/

/* 
Takes the command line and stores it in the scanid field of the datafile struct.
This function is now obsolete.
*/
/*
void putCmdlineInScanid(datafile_definition *datafile, int argc, char **argv, verbose_definition verbose)
{
  int i;
  char txt[1002];
  for(i = 0; i < argc; i++) {
    if(strlen(argv[i]) > 1000) {
      fflush(stdout);
      printwarning(0, "WARNING putCmdlineInScanid: command line option too long.");
    }else {
      addquotesargv(argv[i], txt);
      if(strlen(txt) + strlen(datafile->scanID) + 2 > 1000) {
	fflush(stdout);
	printwarning(0, "WARNING putCmdlineInScanid: command line options too long.");
	break;
      }else {
	char txt3[1002];
	strcpy(txt3, datafile->scanID);
	strcat(txt3, txt);
	strcat(txt3, " ");
	if(set_scanID_PSRData(datafile, txt3, verbose) == 0) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR putCmdlineInScanid: Setting scanID name failed.");
	  exit(0);
	}
      }
    }
  }
}
*/

//START REGION RELEASE

/* Clean the struct which is probably filled with random junk at the
   time of declaration. */
void cleanPSRData(datafile_definition *datafile, verbose_definition verbose) 
{
  memset(datafile, 0, sizeof(datafile_definition));
  datafile->fptr = NULL;
  datafile->fptr_hdr = NULL;
  datafile->fits_fptr = NULL;
  datafile->scales = NULL;
  datafile->offsets = NULL;
  datafile->weights = NULL;
  datafile->weight_stats_set = 0;
  datafile->data = NULL;
  datafile->format = 0;
  datafile->version = 0;
  datafile->opened_flag = 0;
  datafile->enable_write_flag = 0;
  datafile->dumpOnClose = 0;
  //  datafile->dd_mode = 0;
  datafile->filename = malloc(1);
  //  printf("XXXXXX allocated %p\n", datafile->filename);
  if(datafile->filename == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cleanPSRData: Memory allocation error.");
    exit(0);
  }
  datafile->filename[0] = 0;
  datafile->NrSubints = 0;
  datafile->NrBins = 0;
  datafile->NrBits = 0;
  datafile->NrPols = 0;
  datafile->NrFreqChan = 0;
  datafile->isFolded = -1;
  datafile->foldMode = FOLDMODE_UNKNOWN;
  datafile->fixedPeriod = 0;
  datafile->tsampMode = TSAMPMODE_UNKNOWN;
  datafile->fixedtsamp = 0;
  datafile->tsamp_list = NULL;
  datafile->tsubMode = TSUBMODE_UNKNOWN;
  datafile->tsub_list = (double *)malloc(sizeof(double));
  if(datafile->tsub_list == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cleanPSRData: Memory allocation error.");
    exit(0);
  }
  //  printf("XXXXX alloc %p\n", datafile->tsub_list);
  datafile->tsub_list[0] = 0;
  datafile->freq_ref = -2;  // Undefined
  datafile->freqMode = FREQMODE_UNKNOWN;
  datafile->freqlabel_list = NULL; //(double *)malloc(sizeof(double));
  //  if(datafile->freq_list == NULL) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "ERROR cleanPSRData: Memory allocation error.");
  //    exit(0);
  //  }
  //  datafile->freq_list[0] = 0;
  //  datafile->uniform_freq_cent = 0;
  //  datafile->uniform_bw = 0;
  datafile->bandwidth = 0;
  datafile->centrefreq = 0;
  datafile->ra = 0;
  datafile->dec = 0;
  datafile->dm = 0;
  datafile->rm = 0;
  datafile->mjd_start = 0;
  datafile->psrname = malloc(1);
  datafile->observatory = malloc(1);
  datafile->institute = malloc(1);
  datafile->instrument = malloc(1);
  datafile->scanID = malloc(1);
  datafile->observer = malloc(1);
  datafile->projectID = malloc(1);
  if(datafile->psrname == NULL || datafile->observatory == NULL || datafile->institute == NULL || datafile->instrument == NULL || datafile->scanID == NULL || datafile->observer == NULL || datafile->projectID == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR cleanPSRData: Memory allocation error.");
    exit(0);
  }
  datafile->psrname[0] = 0;
  datafile->observatory[0] = 0;
  datafile->institute[0] = 0;
  datafile->instrument[0] = 0;
  datafile->scanID[0] = 0;
  datafile->observer[0] = 0;
  datafile->projectID[0] = 0;
  //  datafile->NrPApoints = 0;
  //  datafile->rmspoints_defined = 0;
  //  datafile->data_pa = NULL;
  //  datafile->data_dpa = NULL;
  //  datafile->longitudes_defined = 0;
  //  datafile->data_long = NULL;
  //  datafile->bins_defined = 0;
  //  datafile->data_bin = NULL;
  //  datafile->rmsvalues = NULL;
  datafile->offpulse_rms = NULL;
  datafile->feedtype = FEEDTYPE_UNKNOWN;
  //  datafile->fd_sang = 0;
  //  datafile->fd_xyph = 0;
  datafile->poltype = POLTYPE_UNKNOWN;
  datafile->datastart = 0;
  datafile->isTransposed = 0;
  datafile->gentype = GENTYPE_UNDEFINED;
  datafile->xrangeset = 0;
  datafile->xrange[0] = 0;
  datafile->xrange[1] = 0;
  datafile->yrangeset = 0;
  datafile->yrange[0] = 0;
  datafile->yrange[1] = 0;
  datafile->isDeDisp  = -1;
  datafile->isDeFarad = -1;
  datafile->isDePar = -1;
  datafile->isDebase = -1;
  //  datafile->telescope_long = 0;
  //  datafile->telescope_lat = 0;
  datafile->telescope_X = 0;
  datafile->telescope_Y = 0;
  datafile->telescope_Z = 0;
  datafile->cableSwap = -1;
  datafile->cableSwapcor = -1;
  datafile->history.timestamp = NULL;
  datafile->history.cmd = NULL;
  datafile->history.user = NULL;
  datafile->history.hostname = NULL;
  datafile->history.nextEntry = NULL;
  /*  
      datafile-> = 0;
  */
}

//START REGION DEVELOP
//START REGION RELEASE

/* Copy the struct to another, except the file pointers. No memory is
   allocated to hold data and the data pointer is not set. The
   destination should have been initialised with cleanPSRData() at
   some point before calling this function.

If return value = 0, then there is an error. Returns 1 on success.
*/
int copy_params_PSRData(datafile_definition datafile_source, datafile_definition *datafile_dest, verbose_definition verbose) 
{
  datafile_dest->fptr = NULL;
  datafile_dest->fptr_hdr = NULL;
  datafile_dest->fits_fptr = NULL;
  datafile_dest->scales = NULL;
  datafile_dest->offsets = NULL;
  datafile_dest->weights = NULL;
  datafile_dest->weight_stats_set = datafile_source.weight_stats_set;
  datafile_dest->weight_stats_zeroweightfound = datafile_source.weight_stats_zeroweightfound;
  datafile_dest->weight_stats_differentweights = datafile_source.weight_stats_differentweights;
  datafile_dest->weight_stats_negativeweights = datafile_source.weight_stats_negativeweights;
  datafile_dest->weight_stats_weightvalue = datafile_source.weight_stats_weightvalue;
  datafile_dest->offpulse_rms = NULL;
  datafile_dest->format = datafile_source.format;
  datafile_dest->version = datafile_source.version;
  /* Don't copy open status and write status */
  /*  datafile_dest->opened_flag = datafile_source.opened_flag; 
      datafile_dest->enable_write_flag = datafile_source.enable_write_flag; 
      datafile_dest->dumpOnClose = datafile_source.dumpOnClose; 
  */
  //  datafile_dest->dd_mode = datafile_source.dd_mode;
  datafile_dest->NrSubints = datafile_source.NrSubints;
  datafile_dest->NrBins = datafile_source.NrBins;
  datafile_dest->NrBits = datafile_source.NrBits;
  datafile_dest->NrPols = datafile_source.NrPols;
  datafile_dest->NrFreqChan = datafile_source.NrFreqChan;
  datafile_dest->isFolded = datafile_source.isFolded;
  datafile_dest->foldMode = datafile_source.foldMode;
  datafile_dest->fixedPeriod = datafile_source.fixedPeriod;
  datafile_dest->tsampMode = datafile_source.tsampMode;
  if(datafile_source.tsampMode == TSAMPMODE_LONGITUDELIST && datafile_source.tsamp_list != NULL) {
    if(datafile_dest->tsamp_list != NULL)
      free(datafile_dest->tsamp_list);
    datafile_dest->tsamp_list = (double *)malloc(datafile_source.NrBins*sizeof(double));
    if(datafile_dest->tsamp_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
      return 0;    
    }
    memcpy(datafile_dest->tsamp_list, datafile_source.tsamp_list, datafile_source.NrBins*sizeof(double));
  }
  datafile_dest->fixedtsamp = datafile_source.fixedtsamp;
  datafile_dest->tsubMode = datafile_source.tsubMode;
  if(datafile_dest->tsub_list != NULL) {
    //    printf("XXXXX free %p\n", datafile_dest->tsub_list);
    free(datafile_dest->tsub_list);
  }
  if(datafile_source.tsubMode == TSUBMODE_TSUBLIST) {
    datafile_dest->tsub_list = (double *)malloc(datafile_source.NrSubints*sizeof(double));
    if(datafile_dest->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
      return 0;    
    }
    //    printf("XXXXX alloc %p\n", datafile_dest->tsub_list);
    memcpy(datafile_dest->tsub_list, datafile_source.tsub_list, datafile_source.NrSubints*sizeof(double));
  }else {
    datafile_dest->tsub_list = (double *)malloc(sizeof(double));
    if(datafile_dest->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
      return 0;    
    }
    //    printf("XXXXX alloc %p\n", datafile_dest->tsub_list);
    memcpy(datafile_dest->tsub_list, datafile_source.tsub_list, sizeof(double));
  }
  datafile_dest->freq_ref = datafile_source.freq_ref;
  datafile_dest->freqMode = datafile_source.freqMode;
  if(datafile_dest->freqlabel_list != NULL) {
    //    printf("XXXXX free %p\n", datafile_dest->freq_list);
    free(datafile_dest->freqlabel_list);
    datafile_dest->freqlabel_list = NULL;
  }
  datafile_dest->bandwidth = datafile_source.bandwidth;
  datafile_dest->centrefreq = datafile_source.centrefreq;
  if(datafile_source.freqMode == FREQMODE_FREQTABLE) {
    datafile_dest->freqlabel_list = (double *)malloc(datafile_source.NrFreqChan*datafile_source.NrSubints*sizeof(double));
    if(datafile_dest->freqlabel_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
      return 0;    
    }
    //    printf("XXXXX alloc %p\n", datafile_dest->tsub_list);
    memcpy(datafile_dest->freqlabel_list, datafile_source.freqlabel_list, datafile_source.NrFreqChan*datafile_source.NrSubints*sizeof(double));
  }else {
    datafile_dest->freqlabel_list = NULL;
  }
  //  datafile_dest->uniform_freq_cent = datafile_source.uniform_freq_cent;
  //  datafile_dest->uniform_bw = datafile_source.uniform_bw;
  datafile_dest->ra = datafile_source.ra;
  datafile_dest->dec = datafile_source.dec;
  datafile_dest->dm = datafile_source.dm;
  datafile_dest->rm = datafile_source.rm;
  datafile_dest->mjd_start = datafile_source.mjd_start;
  //  datafile_dest->filename = NULL;
  if(set_filename_PSRData(datafile_dest, datafile_source.filename, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting file name failed.");
    return 0;
  }
  //  datafile_dest->psrname = NULL;
  if(set_psrname_PSRData(datafile_dest, datafile_source.psrname, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting pulsar name failed.");
    return 0;
  }
  //  datafile_dest->observatory = NULL;
  if(set_observatory_PSRData(datafile_dest, datafile_source.observatory, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting observatory name failed.");
    return 0;
  }
  if(set_institute_PSRData(datafile_dest, datafile_source.institute, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting institute name failed.");
    return 0;
  }
  if(set_instrument_PSRData(datafile_dest, datafile_source.instrument, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting instrument name failed.");
    return 0;
  }
  if(set_scanID_PSRData(datafile_dest, datafile_source.scanID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting scan ID name failed.");
    return 0;
  }
  if(set_observer_PSRData(datafile_dest, datafile_source.observer, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting name of observer failed.");
    return 0;
  }
  if(set_projectID_PSRData(datafile_dest, datafile_source.projectID, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR copy_paramsPSRData: Setting project ID failed.");
    return 0;
  }
  //  datafile_dest->telescope_long = datafile_source.telescope_long;
  //  datafile_dest->telescope_lat = datafile_source.telescope_lat;
  datafile_dest->telescope_X = datafile_source.telescope_X;
  datafile_dest->telescope_Y = datafile_source.telescope_Y;
  datafile_dest->telescope_Z = datafile_source.telescope_Z;

  //  datafile_dest->NrPApoints = datafile_source.NrPApoints;
  //  datafile_dest->rmspoints_defined = datafile_source.rmspoints_defined;

  //  datafile_dest->longitudes_defined = datafile_source.longitudes_defined;
  //  datafile_dest->bins_defined = datafile_source.bins_defined;

  datafile_dest->feedtype = datafile_source.feedtype;
  //  datafile_dest->fd_sang = datafile_source.fd_sang;
  //  datafile_dest->fd_xyph = datafile_source.fd_xyph;
  datafile_dest->poltype = datafile_source.poltype;

  datafile_dest->isTransposed = datafile_source.isTransposed;
  datafile_dest->gentype    = datafile_source.gentype;
  datafile_dest->xrangeset  = datafile_source.xrangeset;
  datafile_dest->xrange[0]  = datafile_source.xrange[0];
  datafile_dest->xrange[1]  = datafile_source.xrange[1];
  datafile_dest->yrangeset  = datafile_source.yrangeset;
  datafile_dest->yrange[0]  = datafile_source.yrange[0];
  datafile_dest->yrange[1]  = datafile_source.yrange[1];
  datafile_dest->isDeDisp = datafile_source.isDeDisp;
  datafile_dest->isDeFarad = datafile_source.isDeFarad;
  datafile_dest->isDePar = datafile_source.isDePar;
  datafile_dest->isDebase = datafile_source.isDebase;
  datafile_dest->cableSwap = datafile_source.cableSwap;
  datafile_dest->cableSwapcor = datafile_source.cableSwapcor;

  datafile_history_entry_definition *dest_hist, *source_hist;
  dest_hist = &(datafile_dest->history);
  source_hist = &(datafile_source.history);
  int ok = 1;
  do {
    dest_hist->timestamp = NULL;
    dest_hist->cmd = NULL;
    dest_hist->user = NULL;
    dest_hist->hostname = NULL;
    dest_hist->nextEntry = NULL;
    if(source_hist->timestamp != NULL) {
      dest_hist->timestamp = malloc(strlen(source_hist->timestamp)+1);
      if(dest_hist->timestamp == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
	return 0;
      }
      strcpy(dest_hist->timestamp, source_hist->timestamp);
    }
    if(source_hist->cmd != NULL) {
      dest_hist->cmd = malloc(strlen(source_hist->cmd)+1);
      if(dest_hist->cmd == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
	return 0;
      }
      strcpy(dest_hist->cmd, source_hist->cmd);
    }
    if(source_hist->user != NULL) {
      dest_hist->user = malloc(strlen(source_hist->user)+1);
      if(dest_hist->user == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
	return 0;
      }
      strcpy(dest_hist->user, source_hist->user);
    }
    if(source_hist->hostname != NULL) {
      dest_hist->hostname = malloc(strlen(source_hist->hostname)+1);
      if(dest_hist->hostname == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
	return 0;
      }
      strcpy(dest_hist->hostname, source_hist->hostname);
    }
    if(source_hist->nextEntry != NULL) {
      dest_hist->nextEntry = malloc(sizeof(datafile_history_entry_definition));
      if(dest_hist->nextEntry == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR copy_paramsPSRData: Memory allocation error.");
	return 0;
      }
      dest_hist = dest_hist->nextEntry;
      source_hist = source_hist->nextEntry;
    }else {
      ok = 0;
    }
  }while(ok);
  /*
    datafile_dest-> = datafile_source.;
  */
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* No seeking is done (important to allow for instance multiple files
   are dumped in a single dump file), header is written out
   immediately to the header file pointer defined in the datafile
   struct, so file should already be opened. The datastart variable in
   datafile is set to the position where the data can be found.

   Return 1 on success, 0 on error.
 */
int writePSRSALSAHeader(datafile_definition *datafile, verbose_definition verbose)
{
  int ret, dummyi;
  char identifier[] = "PSRSALSAdump";
  int version = 2;


  // Write identifier and version number
  ret = fwrite(identifier, 12, 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&version, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  
  // Write strings. These can have variable length. Their size is written out first, followed by the string itself
  dummyi = strlen(datafile->psrname);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->psrname, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }

  dummyi = strlen(datafile->scanID);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->scanID, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }

  dummyi = strlen(datafile->observer);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->observer, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }

  dummyi = strlen(datafile->projectID);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->projectID, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }

  dummyi = strlen(datafile->observatory);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->observatory, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }

  dummyi = strlen(datafile->instrument);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->instrument, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }

  dummyi = strlen(datafile->institute);
  ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fwrite(datafile->institute, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }

  // Write out some other general information
  ret = fwrite(&(datafile->telescope_X), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->telescope_Y), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->telescope_Z), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->NrBits), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->isDeDisp), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->isDeFarad), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->isDePar), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->isDebase), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->cableSwap), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->cableSwapcor), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->dm), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->rm), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->freq_ref), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->feedtype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->poltype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->mjd_start), sizeof(long double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }

  // The dimensions of the data
  ret = fwrite(&(datafile->NrSubints), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->NrBins), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->NrPols), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->NrFreqChan), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }

  // Information about folding/time resolution
  ret = fwrite(&(datafile->isFolded), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->foldMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->fixedPeriod), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->tsampMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->fixedtsamp), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(datafile->tsampMode == TSAMPMODE_LONGITUDELIST) {
    ret = fwrite(datafile->tsamp_list, sizeof(double), datafile->NrBins, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  ret = fwrite(&(datafile->tsubMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(datafile->tsubMode == TSUBMODE_TSUBLIST) {
    ret = fwrite(datafile->tsub_list, sizeof(double), datafile->NrSubints, datafile->fptr_hdr);
    if(ret != datafile->NrSubints) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s (return value is %ld)", datafile->filename, ret);
      return 0;
    }
  }else {
    ret = fwrite(datafile->tsub_list, sizeof(double), 1, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }

  // Information about pointing
  ret = fwrite(&(datafile->ra), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->dec), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }


  // Information about the frequency coverage
  ret = fwrite(&(datafile->freqMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->centrefreq), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->bandwidth), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }

  if(datafile->freqMode == FREQMODE_FREQTABLE) {
    ret = fwrite(datafile->freqlabel_list, sizeof(double), datafile->NrSubints*datafile->NrFreqChan, datafile->fptr_hdr);
    if(ret != datafile->NrSubints*datafile->NrFreqChan) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }else if(datafile->freqMode != FREQMODE_UNKNOWN && datafile->freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Frequency mode is not implemented");
    return 0;
  }

  // Data reduction information.
  ret = fwrite(&(datafile->gentype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->isTransposed), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->xrangeset), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(datafile->xrangeset) {
    ret = fwrite(datafile->xrange, sizeof(float), 2, datafile->fptr_hdr);
    if(ret != 2) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  ret = fwrite(&(datafile->yrangeset), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  if(datafile->yrangeset) {
    ret = fwrite(datafile->yrange, sizeof(float), 2, datafile->fptr_hdr);
    if(ret != 2) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }


  // Other information, which is a bit messy at the moment


  /*
  ret = fwrite(&(datafile->fd_sang), sizeof(float), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  ret = fwrite(&(datafile->fd_xyph), sizeof(float), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  */

  datafile->datastart = ftell(datafile->fptr_hdr);

  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Writes out the history table with command lines etc. 
   Returns 1 on success, 0 on error */
int writeHistoryPSRSALSA(datafile_definition *datafile, verbose_definition verbose)
{
  // Write history, start by counting nr of history lines
  datafile_history_entry_definition *curHistoryEntry;
  int nrHistoryLines, dummyi;
  size_t ret;
  char *history_id = "HISTORY";
  curHistoryEntry = &(datafile->history);
  nrHistoryLines = 0;
  do {
    if(curHistoryEntry->timestamp != NULL || curHistoryEntry->cmd != NULL || curHistoryEntry->user != NULL || curHistoryEntry->hostname != NULL || curHistoryEntry->nextEntry != NULL) {
      nrHistoryLines++;
      curHistoryEntry = curHistoryEntry->nextEntry; // Goto next entry
    }
  }while(curHistoryEntry != NULL);
  if(verbose.debug) {
    printf("Start writing history in PSRSALSA binary format (%d lines)\n", nrHistoryLines);
  }
  ret = fwrite(history_id, 1, strlen(history_id), datafile->fptr_hdr);  // Write out history identifier, allowing us to determine if there is a history table in the file
  if(ret != strlen(history_id)) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }  
  ret = fwrite(&(nrHistoryLines), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
    return 0;
  }
  //Now write out the actual history information
  curHistoryEntry = &(datafile->history);
  do {
    if(curHistoryEntry->timestamp != NULL || curHistoryEntry->cmd != NULL || curHistoryEntry->user != NULL || curHistoryEntry->hostname != NULL || curHistoryEntry->nextEntry != NULL) {
      if(curHistoryEntry->timestamp == NULL) {
	dummyi = 0;
      }else {
	dummyi = strlen(curHistoryEntry->timestamp);
      }
      ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
	return 0;
      }
      if(dummyi > 0) {
	ret = fwrite(curHistoryEntry->timestamp, sizeof(char), dummyi, datafile->fptr_hdr);
	if(ret != dummyi) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
	  return 0;
	}
      }

      if(curHistoryEntry->user == NULL) {
	dummyi = 0;
      }else {
	dummyi = strlen(curHistoryEntry->user);
      }
      ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
	return 0;
      }
      if(dummyi > 0) {
	ret = fwrite(curHistoryEntry->user, sizeof(char), dummyi, datafile->fptr_hdr);
	if(ret != dummyi) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
	  return 0;
	}
      }

      if(curHistoryEntry->hostname == NULL) {
	dummyi = 0;
      }else {
	dummyi = strlen(curHistoryEntry->hostname);
      }
      ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
	return 0;
      }
      if(dummyi > 0) {
	ret = fwrite(curHistoryEntry->hostname, sizeof(char), dummyi, datafile->fptr_hdr);
	if(ret != dummyi) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
	  return 0;
	}
      }

      if(curHistoryEntry->cmd == NULL) {
	dummyi = 0;
      }else {
	dummyi = strlen(curHistoryEntry->cmd);
      }
      ret = fwrite(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
	return 0;
      }
      if(dummyi > 0) {
	ret = fwrite(curHistoryEntry->cmd, sizeof(char), dummyi, datafile->fptr_hdr);
	if(ret != dummyi) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
	  return 0;
	}
      }

      curHistoryEntry = curHistoryEntry->nextEntry; // Goto next entry
    }
  }while(curHistoryEntry != NULL);
  datafile->datastart = ftell(datafile->fptr_hdr);
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* No seeking is done (important to allow for instance multiple files
   are dumped in a single dump file), header is read immediately from
   the header file pointer defined in the datafile struct, so file
   should already be opened. If nohistory_expected is nonzero, no
   warning is generated if there is no history (but a warning is
   generated if there is).

   Return 1 on success, 0 on error.
 */
int readPSRSALSAHeader(datafile_definition *datafile, int nohistory_expected, verbose_definition verbose)
{
  int ret, dummyi;
  char identifier[13], *txt;
  int version, maxversion_supported = 2;

  txt = malloc(10000);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error", datafile->filename);
    return 0;
  }

  // Read identifier and version number
  ret = fread(identifier, 1, 12, datafile->fptr_hdr);
  if(ret != 12) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  identifier[12] = 0;
  if(strcmp(identifier, "PSRSALSAdump") != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: File does not appear to be in the expected format", datafile->filename);
    return 0;
  }

  ret = fread(&version, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(version < 1 || version > maxversion_supported) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: File %s is in an unsupported version number (%d). The maximum supported version in your installation is %d. An update might be required to read this file.", datafile->filename, version, maxversion_supported);
    return 0;
  }
  if(verbose.debug) {
    printf("  PSRSALSA file is version %d\n", version);
  }

  // Read strings. These can have variable length. Their size is written out first, followed by the string itself
  ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(dummyi > 9999) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    txt[dummyi] = 0;
    if(set_psrname_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting pulsar name failed for %s", datafile->filename);
      return 0;
    }
  }

  ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(dummyi > 9999) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    txt[dummyi] = 0;
    if(set_scanID_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting scan ID failed for %s", datafile->filename);
      return 0;
    }
  }

  if(version > 1) {
    ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    if(dummyi > 9999) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
      return 0;
    }
    if(dummyi > 0) {
      ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
      if(ret != dummyi) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
	return 0;
      }
      txt[dummyi] = 0;
      if(set_observer_PSRData(datafile, txt, verbose) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting scan ID failed for %s", datafile->filename);
	return 0;
      }
    }

    ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    if(dummyi > 9999) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
      return 0;
    }
    if(dummyi > 0) {
      ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
      if(ret != dummyi) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
	return 0;
      }
      txt[dummyi] = 0;
      if(set_projectID_PSRData(datafile, txt, verbose) == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting scan ID failed for %s", datafile->filename);
	return 0;
      }
    }
  }

  ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(dummyi > 9999) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    txt[dummyi] = 0;
    if(set_observatory_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting observatory name failed for %s", datafile->filename);
      return 0;
    }
  }


  ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(dummyi > 9999) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    txt[dummyi] = 0;
    if(set_instrument_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting instrument name failed for %s", datafile->filename);
      return 0;
    }
  }

  ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(dummyi > 9999) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Maximum string length exceeded for %s", datafile->filename);
    return 0;
  }
  if(dummyi > 0) {
    ret = fread(txt, sizeof(char), dummyi, datafile->fptr_hdr);
    if(ret != dummyi) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    txt[dummyi] = 0;
    if(set_institute_PSRData(datafile, txt, verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Setting institute name failed for %s", datafile->filename);
      return 0;
    }
  }


  // Read in some other general information
  ret = fread(&(datafile->telescope_X), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->telescope_Y), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->telescope_Z), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->NrBits), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->isDeDisp), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->isDeFarad), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->isDePar), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->isDebase), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->cableSwap), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->cableSwapcor), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->dm), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->rm), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->freq_ref), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->feedtype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->poltype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->mjd_start), sizeof(long double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }

  // The dimensions of the data
  ret = fread(&(datafile->NrSubints), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->NrBins), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->NrPols), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->NrFreqChan), sizeof(long), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }

  // Information about folding
  ret = fread(&(datafile->isFolded), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->foldMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->fixedPeriod), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->tsampMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->fixedtsamp), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  datafile->tsamp_list = NULL;
  if(datafile->tsampMode == TSAMPMODE_LONGITUDELIST) {
    datafile->tsamp_list = (double *)malloc(sizeof(double)*datafile->NrBins);
    if(datafile->tsamp_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading %s", datafile->filename);
      return 0;
    }
    ret = fread(datafile->tsamp_list, sizeof(double), datafile->NrBins, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writePSRSALSAHeader: Write error to %s", datafile->filename);
      return 0;
    }
  }
  ret = fread(&(datafile->tsubMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(datafile->tsub_list != NULL)
    free(datafile->tsub_list);
  if(datafile->tsubMode == TSUBMODE_TSUBLIST) {
    datafile->tsub_list = (double *)malloc(sizeof(double)*datafile->NrSubints);
    if(datafile->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading %s", datafile->filename);
      return 0;
    }
    ret = fread(datafile->tsub_list, sizeof(double), datafile->NrSubints, datafile->fptr_hdr);
    if(ret != datafile->NrSubints) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
  }else {
    datafile->tsub_list = (double *)malloc(sizeof(double));
    if(datafile->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading %s", datafile->filename);
      return 0;
    }
    ret = fread(datafile->tsub_list, sizeof(double), 1, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
  }

  // Information about pointing
  ret = fread(&(datafile->ra), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->dec), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }


  // Information about the frequency coverage
  ret = fread(&(datafile->freqMode), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->centrefreq), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->bandwidth), sizeof(double), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(datafile->freqMode == FREQMODE_FREQTABLE) {
    if(datafile->freqlabel_list != NULL) {
      free(datafile->freqlabel_list);
    }
    datafile->freqlabel_list = malloc(datafile->NrSubints*datafile->NrFreqChan*sizeof(double));
    if(datafile->freqlabel_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error");
      return 0;
    }
    ret = fread(datafile->freqlabel_list, sizeof(double), datafile->NrSubints*datafile->NrFreqChan, datafile->fptr_hdr);
    if(ret != datafile->NrSubints*datafile->NrFreqChan) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
  }

  // Data reduction information.
  ret = fread(&(datafile->gentype), sizeof(int), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->isTransposed), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->xrangeset), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(datafile->xrangeset) {
    ret = fread(datafile->xrange, sizeof(float), 2, datafile->fptr_hdr);
    if(ret != 2) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
  }
  ret = fread(&(datafile->yrangeset), sizeof(char), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  if(datafile->yrangeset) {
    ret = fread(datafile->yrange, sizeof(float), 2, datafile->fptr_hdr);
    if(ret != 2) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
  }


  // Other information, which is a bit messy at the moment


  /*
  ret = fread(&(datafile->fd_sang), sizeof(float), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  ret = fread(&(datafile->fd_xyph), sizeof(float), 1, datafile->fptr_hdr);
  if(ret != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
    return 0;
  }
  */

  int readhistory;
  off_t filepos;
  filepos = ftello(datafile->fptr_hdr); // Remember position before attempting to read the history, in case it is not there
  readhistory = 1;
  ret = fread(identifier, 1, 7, datafile->fptr_hdr);
  if(ret != 7) {
    fflush(stdout);
    if(nohistory_expected == 0) {
      printwarning(verbose.debug, "WARNING readPSRSALSAHeader: No history table found in %s", datafile->filename);
    }
    readhistory = 0;
    fseeko(datafile->fptr_hdr, filepos, SEEK_SET);
  }else {
    identifier[7] = 0;
    if(strcmp(identifier, "HISTORY") != 0) {
      fflush(stdout);
      if(nohistory_expected == 0) {
	printwarning(verbose.debug, "WARNING readPSRSALSAHeader: No history table found in %s", datafile->filename);
      }
      if(verbose.debug) {
	printf("First bytes are '%c' '%c' and '%c'.\n", identifier[0], identifier[1], identifier[2]);
      }
      readhistory = 0;
      fseeko(datafile->fptr_hdr, filepos, SEEK_SET);
    }
  }
  if(nohistory_expected && readhistory) {
    printwarning(verbose.debug, "WARNING readPSRSALSAHeader: History table found in %s, while this was not expected.", datafile->filename);
  }

  if(readhistory) {
    // Read history, start by reading the nr of history lines
    datafile_history_entry_definition *curHistoryEntry;
    int curline, nrHistoryLines;
    ret = fread(&nrHistoryLines, sizeof(int), 1, datafile->fptr_hdr);
    if(ret != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
      return 0;
    }
    
    //Now read in the actual history information
    curHistoryEntry = &(datafile->history);
    for(curline = 0; curline < nrHistoryLines; curline++) {
      
      ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
	return 0;
      }
      if(dummyi > 0) {
	curHistoryEntry->timestamp = malloc(dummyi+1);
	if(curHistoryEntry->timestamp == NULL) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading in %s", datafile->filename);
	  return 0;
	}
	ret = fread(curHistoryEntry->timestamp, sizeof(char), dummyi, datafile->fptr_hdr);
	if(ret != dummyi) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
	  return 0;
	}
	curHistoryEntry->timestamp[dummyi] = 0;
      }
      
      ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
	return 0;
      }
      if(dummyi > 0) {
	curHistoryEntry->user = malloc(dummyi+1);
	if(curHistoryEntry->user == NULL) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading in %s", datafile->filename);
	  return 0;
	}
	ret = fread(curHistoryEntry->user, sizeof(char), dummyi, datafile->fptr_hdr);
	if(ret != dummyi) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
	  return 0;
	}
	curHistoryEntry->user[dummyi] = 0;
      }
      
      ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
	return 0;
      }
      if(dummyi > 0) {
	curHistoryEntry->hostname = malloc(dummyi+1);
	if(curHistoryEntry->hostname == NULL) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading in %s", datafile->filename);
	  return 0;
	}
	ret = fread(curHistoryEntry->hostname, sizeof(char), dummyi, datafile->fptr_hdr);
	if(ret != dummyi) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
	  return 0;
	}
	curHistoryEntry->hostname[dummyi] = 0;
      }
      
      ret = fread(&dummyi, sizeof(int), 1, datafile->fptr_hdr);
      if(ret != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
	return 0;
      }
      if(dummyi > 0) {
	curHistoryEntry->cmd = malloc(dummyi+1);
	if(curHistoryEntry->cmd == NULL) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading in %s", datafile->filename);
	  return 0;
	}
	ret = fread(curHistoryEntry->cmd, sizeof(char), dummyi, datafile->fptr_hdr);
	if(ret != dummyi) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readPSRSALSAHeader: Read error from %s", datafile->filename);
	  return 0;
	}
	curHistoryEntry->cmd[dummyi] = 0;
      }
      
      // Add another entry if required
      if(curline < nrHistoryLines - 1) {
	curHistoryEntry->nextEntry = malloc(sizeof(datafile_history_entry_definition));
	if(curHistoryEntry->nextEntry == NULL) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readPSRSALSAHeader: Memory allocation error while reading in %s", datafile->filename);
	  return 0;
	}
	curHistoryEntry = curHistoryEntry->nextEntry; // Goto next entry
	curHistoryEntry->timestamp = NULL;  // Will be filled above in next loop, if they are defined.
	curHistoryEntry->cmd = NULL;
	curHistoryEntry->user = NULL;
	curHistoryEntry->hostname = NULL;
	curHistoryEntry->nextEntry = NULL;
      }
    }
  } // End of reading in history table

  // Obtain header length by reading current position
  datafile->datastart = ftell(datafile->fptr_hdr);

  //  if(verbose.debug) {
  //    printf("  XXXX frequency information stored in %p\n", datafile->freq_list);
  //  }

  free(txt);
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

int readPulsePSRSALSAData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  // Data is organized as subint, freq, pol, bin, like in memory format
  long long filepos;
  size_t ret;
  filepos = datafile.NrBins*(polarization+datafile.NrPols*(freq+pulsenr*datafile.NrFreqChan))+binnr;
  filepos *= sizeof(float);
  filepos += datafile.datastart;
  fseeko(datafile.fptr, filepos, SEEK_SET);
  ret = fread(pulse, sizeof(float), nrSamples, datafile.fptr);
  if(ret != nrSamples) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPulsePSRSALSAData: File read failed.");
    return 0;
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

int writePulsePSRSALSAData(datafile_definition datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  // Data is organized as subint, freq, pol, bin, like in memory format
  long long filepos;
  size_t ret;
  filepos = datafile.NrBins*(polarization+datafile.NrPols*(freq+pulsenr*datafile.NrFreqChan))+binnr;
  filepos *= sizeof(float);
  filepos += datafile.datastart;
  fseeko(datafile.fptr, filepos, SEEK_SET);
  ret = fwrite(pulse, sizeof(float), nrSamples, datafile.fptr);
  if(ret != nrSamples) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePulsePSRSALSAData: File write failed.");
    return 0;
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

int readPSRSALSAfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  // Data is organized as subint, freq, pol, bin, like in memory format
  long n, f, p;
  size_t ret;
  if(verbose.verbose) {
    printf("Start reading PSRSALSA binary file\n");
  }
  fseeko(datafile.fptr, datafile.datastart, SEEK_SET);
  for(n = 0; n < datafile.NrSubints; n++) {
    for(f = 0; f < datafile.NrFreqChan; f++) {
      for(p = 0; p < datafile.NrPols; p++) {
	if(verbose.verbose && verbose.nocounters == 0) 
	  printf("  Progress reading PSRSALSA binary file (%.1f%%)\r", 100.0*(p+(f+n*datafile.NrFreqChan)*datafile.NrPols)/(float)(datafile.NrSubints*datafile.NrFreqChan*datafile.NrPols));
	ret = fread(&data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))], sizeof(float), datafile.NrBins, datafile.fptr);
	if(ret != datafile.NrBins) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR readPSRSALSAfile: File read failed.");
	  return 0;
	}
      }
    }
  }
  if(verbose.verbose) printf("  Reading is done.                                \n");
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

int writePSRSALSAfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  long n, f, p;
  size_t ret;
  fseeko(datafile.fptr, datafile.datastart, SEEK_SET);
  for(n = 0; n < datafile.NrSubints; n++) {
    for(f = 0; f < datafile.NrFreqChan; f++) {
      for(p = 0; p < datafile.NrPols; p++) {
	if(verbose.verbose && verbose.nocounters == 0) 
	  printf("  Progress writing PSRSALSA binary file (%.1f%%)\r", 100.0*(p+(f+n*datafile.NrFreqChan)*datafile.NrPols)/(float)(datafile.NrSubints*datafile.NrFreqChan*datafile.NrPols));
	ret = fwrite(&data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))], sizeof(float), datafile.NrBins, datafile.fptr);
	if(ret != datafile.NrBins) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR writePSRSALSAfile: File write failed.");
	  return 0;
	}
      }
    }
  }
  if(verbose.verbose) printf("  Writing is done.                                   \n");
  return 1;
}


//START REGION DEVELOP
//START REGION RELEASE

/* Copy the string source in dest. If no memory is allocated (!= NULL)
to hold the string it is allocated. If already momory is allocated, it
will be freed first. If source is the NULL pointer, the pointer in the
data struct will point to an empty string. A \n or \r will be
stripped. This functions are is used to handle the filename string in
the datafile_definition struct.

If return value = 0, then there is a memory allocation
error. Returns 1 on success.  */
int set_string_PSRData(char **dest, char *source, verbose_definition verbose)
{
  int i;
  if(*dest != NULL) {
    //    printf("XXXXX free %p\n", *dest);
    free(*dest);
    *dest = NULL;
  }
  if(source != NULL) {
    *dest = malloc(strlen(source)+1);
    if(*dest == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR set_string_PSRData: Memory allocation error.");
      return 0;
    }
    strcpy(*dest, source);
    for(i = 0; i < 2; i++) { // Strip any trailing return characters
      if(strlen(*dest) > 0) {
	if((*dest)[strlen(*dest)-1] == '\n' || (*dest)[strlen(*dest)-1] == '\r') {
	  (*dest)[strlen(*dest)-1] = 0;
	}
      }
    }
  }else {
    *dest = malloc(1);
    if(*dest == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR set_string_PSRData: Memory allocation error.");
      return 0;
    }
    (*dest)[0] = 0;
  }
  return 1;
}

/* Copy the string filename in the data structure. If no memory is
allocated to hold the filename it is allocated. If already momory is
allocated, it will be freed first. If filename is the NULL pointer,
the pointer in the data struct will point to an empty string.

If return value = 0, then there is a memory allocation
error. Returns 1 on success.  */
int set_filename_PSRData(datafile_definition *datafile_dest, char *filename, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->filename), filename, verbose);
}

int set_psrname_PSRData(datafile_definition *datafile_dest, char *psrname, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->psrname), psrname, verbose);
}

int set_observatory_PSRData(datafile_definition *datafile_dest, char *observatory, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->observatory), observatory, verbose);
}

int set_institute_PSRData(datafile_definition *datafile_dest, char *institute, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->institute), institute, verbose);
}

int set_instrument_PSRData(datafile_definition *datafile_dest, char *instrument, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->instrument), instrument, verbose);
}

int set_scanID_PSRData(datafile_definition *datafile_dest, char *scanID, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->scanID), scanID, verbose);
}

int set_observer_PSRData(datafile_definition *datafile_dest, char *observer, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->observer), observer, verbose);
}

int set_projectID_PSRData(datafile_definition *datafile_dest, char *projectID, verbose_definition verbose)
{
  return set_string_PSRData(&(datafile_dest->projectID), projectID, verbose);
}
      
//START REGION DEVELOP
//START REGION RELEASE

/* Tries to guess what format the file is. Verbose level determines 
   nr of spaces before output. If debug is set, more information is
   printed. If noerror is set, no error is generated if data format
   cannot be determined. Returns format:
     -3 is an other error
     -2 is a known type which cannot be read
     -1 is file opening/read error
      0 is unknown type

//START REGION DEVELOP
   Formats that should be recognized are: 
         PUMA_format, 
         PSRCHIVE_ASCII_format, 
         EPN_format,
         FITS_format, 
	 PSRSALSA_BINARY_format
         AO_ASCII_1_format
         AO_ASCII_2_format
	 PPOL_format
	 PPOL_SHORT_format
	 SIGPROC_format
//START REGION RELEASE
*/
int guessPSRData_format(char *filename, int noerror, verbose_definition verbose)
{
  FILE *fin;
  char *txt;
  int i, indent, ok;
  float fdummy;
//START REGION DEVELOP
  int n, j;
//START REGION RELEASE

  txt = malloc(10001);
  if(txt == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR guessPSRData_format: Cannot allocate memory");
    return -3;    
  }

  if(verbose.verbose) {
    for(indent = 0; indent < verbose.indent; indent++)      
      printf(" ");
    printf("Trying to determine data format of file '%s'\n", filename);
    if(verbose.debug) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      printf("  Opening '%s'\n", filename);
    }
  }
  fin = fopen(filename, "r");
  if(fin == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR guessPSRData_format: opening file '%s' failed.", filename);
    free(txt);
    return -1;
  }
  i = fread(txt, 1, 3, fin);
  txt[3] = 0;
  if(i != 3) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR guessPSRData_format: reading file '%s' failed.", filename);
    free(txt);
    return -1;
  }
  if(verbose.debug) {
    printf("  First three bytes have values %d, %d and %d.\n", txt[0], txt[1], txt[2]);
    printf("  This corresponds to characters '%c', '%c' and '%c'.\n", txt[0], txt[1], txt[2]);
  }
  if(strcmp(txt, "PSR") == 0) {  // Check psrsalsa binary format
    if(verbose.debug) {
      if(verbose.verbose) {
	for(indent = 0; indent < verbose.indent; indent++)      
	  printf(" ");
      }
      printf("  file might be in PSRSALSA binary format.\n");
    }
    // Read a bit more
    rewind(fin);
    i = fread(txt, 1, 12, fin);
    txt[12] = 0;
    if(i != 12) {
      fflush(stdout);
      if(verbose.debug) {
	if(verbose.verbose) {
	  for(indent = 0; indent < verbose.indent; indent++)      
	    printf(" ");
	  printf("  file is not in PSRSALSA binary format.\n");
	}
      }
      rewind(fin);
      i = fread(txt, 1, 3, fin);
      txt[3] = 0;
    }else {
      if(strcmp(txt, "PSRSALSAdump") == 0) {
	if(verbose.verbose) {
	  for(indent = 0; indent < verbose.indent; indent++)      
	    printf(" ");
	  printf("  file is in PSRSALSA binary format.\n");
	}
	fclose(fin);
	free(txt);
	return PSRSALSA_BINARY_format;
      }else {
	fflush(stdout);
	if(verbose.debug) {
	  for(indent = 0; indent < verbose.indent; indent++)      
	    printf(" ");
	  printf("  file is not in PSRSALSA binary format.\n");
	}
	rewind(fin);
	i = fread(txt, 1, 3, fin);
	txt[3] = 0;
      }
    }
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    printf("  file is not in PSRSALSA binary format.\n");
  }
  
  if(strcmp(txt, "DPC") == 0) {  // Check for puma format
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      printf("  file is in WSRT PuMa 1 format.\n");
    }
    fclose(fin);
    free(txt);
    return PUMA_format;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    printf("  file is not in WSRT PuMa 1 format.\n");
  }
  if(strcmp(txt, "Fil") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    if(verbose.verbose) 
      printf("  file is in PSRCHIVE ASCII format.\n");
    fclose(fin);
    free(txt);
    return PSRCHIVE_ASCII_format;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    printf("  file is not PSRCHIVE ASCII format.\n");
  }
  if(strcmp(txt, "EPN") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    if(verbose.verbose) 
      printf("  file is in EPN format.\n");
    fclose(fin);
    free(txt);
    return EPN_format;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    printf("  file is not EPN format.\n");
  }
  if(txt[0] == 12 && txt[1] == 0 && txt[2] == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    if(verbose.verbose) 
      printf("  file is in SIGPROC format.\n");
    fclose(fin);
    free(txt);
    return SIGPROC_format;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    printf("  file is not SIGPROC format.\n");
  }
  if(strcmp(txt, "SIM") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    if(verbose.verbose) 
      printf("  file is in PSRFITS format.\n");
    fclose(fin);
    free(txt);
    return FITS_format;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    printf("  file is not PSRFITS format.\n");
  }
//START REGION RELEASE
  if(strcmp(txt, "#pp") == 0
//START REGION DEVELOP
     || strcmp(txt, "#ma") == 0
//START REGION RELEASE
     ) {
    if(verbose.verbose && verbose.debug) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    if(verbose.debug) 
      printf("  File might be in PPOL or PPOLSHORT format.\n");
    if(fgets(txt, 10000, fin) == NULL) {
      printf("ERROR guessPSRData_format: Cannot read first line\n");
      free(txt);
      return 0;
    }
    if(fgets(txt, 10000, fin) == NULL) {
      printf("ERROR guessPSRData_format: Cannot read second line\n");
      free(txt);
      return 0;
    }
    ok = 0;
    i = sscanf(txt, "%f %f %f", &fdummy, &fdummy, &fdummy);
    if(i == 3) {
      ok = 1;
      i = sscanf(txt, "%f %f %f %f %f %f %f %f %f", &fdummy, &fdummy, &fdummy, &fdummy, &fdummy, &fdummy, &fdummy, &fdummy, &fdummy);
      if(i == 9)
	ok = 2;
    }
    fclose(fin);
    if(ok == 2) {
      if(verbose.verbose) {
	for(indent = 0; indent < verbose.indent; indent++)      
	  printf(" ");
      }
      if(verbose.verbose) printf("  file is in PPOL format.\n");
      free(txt);
      return PPOL_format;
    }else if(ok == 1) {
      if(verbose.verbose) {
	for(indent = 0; indent < verbose.indent; indent++)      
	  printf(" ");
      }
      if(verbose.verbose) printf("  file is in PPOLSHORT format.\n");
      free(txt);
      return PPOL_SHORT_format;
    }else if(verbose.debug) {
      if(verbose.verbose) {
	for(indent = 0; indent < verbose.indent; indent++)      
	  printf(" ");
      }
      printf("  file is not in PPOL or PPOLSHORT format.\n");
    }
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    printf("  file is not in PPOL or PPOLSHORT format.\n");
  }
//START REGION DEVELOP
  if(txt[0] == ' ') {  /* Could be AO ASCII format 1 */
    if(verbose.verbose && verbose.debug) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    if(verbose.debug) 
      printf("  file could be AO ASCII format 1 (starting with space).\n");
    i = fread(txt, 1, 10000, fin);
    ok = 1;
    n = 0;
    for(j = 0; j < i; j++) {
      if(txt[j] < 0) {
	ok = 0;
	break;
      }
      if(txt[j] == ' ')
	n++;
    }
    if((float)i/(float)n > 10)          /* If not enough spaces, not AO ASCII format 1 */
      ok = 0;
    if(ok) {
      if(verbose.verbose && verbose.debug) {
	for(indent = 0; indent < verbose.indent; indent++)      
	  printf(" ");
      }
      if(verbose.debug) 
	printf("  file could still be AO ASCII format 1.\n");
      rewind(fin);
      txt[0] = 0;
      fscanf(fin, "%s", txt);
      fscanf(fin, "%s", txt);
      n = 0;
      sscanf(txt, "%d", &n);
      if(n > 10 && n < 10000) {    /* If outside this range, it is unlikely to be a valid number of bins */
	if(verbose.debug) {
	  printf("  file could still be AO ASCII format 1 since the number of bins appears to be %d, which is plausible.\n", n);
	}
	fscanf(fin, "%s", txt);
	fscanf(fin, "%s", txt);
	i = 0;
	sscanf(txt, "%d", &i);
	n *= i;
	if(i >= 1 && i <= 4) {    /* Only proceed when nr polarizations make sense */
	  if(verbose.debug) {
	    printf("  file could still be AO ASCII format 1 since the number of polarizations appears to be %d, which is plausible.\n", i);
	  }
	  fscanf(fin, "%s", txt); //Pulsar name
	  do {
	    i = fgetc(fin);
	  }while(i != EOF && i != '\n'); // Rest of line
	  do {
	    i = fgetc(fin);
	  }while(i != EOF && i != '\n'); // Skip next line as well

	  for(j = 0; j < 6; j++) {
	    txt[0] = 0;
	    fscanf(fin, "%s", txt);
	    //	    if(j < 3)
	    //	      fprintf(stderr, "XXXXX first entry = %s\n", txt);
	  }
	  do {
	    i = fgetc(fin);
	  }while(i != EOF && i != '\n'); // Skip next line as well, since in some cases there appears to be a 7th value (integer)
	  //	  fgets(txt, 100, fin);   
	  /*	  fscanf(fin, "%s", txt);
	  fscanf(fin, "%s", txt);
	  fscanf(fin, "%s", txt);
	  fscanf(fin, "%s", txt);
	  fscanf(fin, "%s", txt);
	  fscanf(fin, "%s", txt); */
	  //	  for(j = 0; j < n+9; j++) {
	  for(j = 0; j < n+3; j++) {
	    txt[0] = 0;
	    fscanf(fin, "%s", txt);
	    //	    if(j < 3)
	    //	      fprintf(stderr, "XXXXX first entry = %s\n", txt);
	  }
	  /*	  fscanf(fin, "%s", txt);
	  fscanf(fin, "%s", txt);
	  fscanf(fin, "%s", txt); */
	  if(strcmp(txt, "1") == 0 || strcmp(txt, "1.00000000") == 0) {
	    if(verbose.verbose) {
	      for(indent = 0; indent < verbose.indent; indent++)      
		printf(" ");
	    }
	    if(verbose.verbose) printf("  file is in AO ASCII format 1.\n");
	    fclose(fin);
	    free(txt);
	    return AO_ASCII_1_format;
	  }else {
	    if(verbose.debug) {
	      if(verbose.verbose) {
		for(indent = 0; indent < verbose.indent; indent++)      
		  printf(" ");
	      }
	      printf("  file is not in AO ASCII format 1 (expected '1', got %s).\n", txt);
	    }
	  }
	}else {
	  if(verbose.debug) {
	    if(verbose.verbose) {
	      for(indent = 0; indent < verbose.indent; indent++)      
		printf(" ");
	    }
	    if(verbose.debug) printf("  file is not in AO ASCII format 1 (nr polarizations=%d).\n", i);
	  }
	}
      }else {
	if(verbose.debug) {
	  if(verbose.verbose) {
	    for(indent = 0; indent < verbose.indent; indent++)      
	      printf(" ");
	  }
	  if(verbose.verbose) printf("  file is not in AO ASCII format 1 (nr bins=%d).\n", n);
	}
      }
    }
    if(verbose.verbose && verbose.debug) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    if(verbose.debug) 
      printf("  file is not AO ASCII format 1.\n");
    rewind(fin);
    i = fread(txt, 1, 3, fin);   /* Leave things they were before the test */
    txt[3] = 0;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    printf("  file is not AO ASCII format 1.\n");
  }

  if(txt[0] == '#' && txt[1] == ' ') {  /* Could be AO ASCII format 2 */
    if(verbose.verbose && verbose.debug) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    if(verbose.debug) 
      printf("  file could be AO ASCII format 2 (starting with # and space).\n");
    rewind(fin);
    if(fgets(txt, 10000, fin) == NULL) {
      printf("ERROR guessPSRData_format: Cannot read first line\n");
      free(txt);
      return 0;
    }
    if(strlen(txt) < 200) {
      char *txt2;
      int nrwords;
      txt2 = pickWordFromString(txt, 1, &nrwords, 0, ' ', verbose);
      if(nrwords == 11) {
	float value;
	txt2 = pickWordFromString(txt, 2, &nrwords, 0, ' ', verbose);
	i = sscanf(txt2, "%f", &value);
	if(i == 1) {
	  if(value >= 39126 || value < 416787) {
	    txt2 = pickWordFromString(txt, 8, &nrwords, 0, ' ', verbose);
	    i = sscanf(txt2, "%f", &value);
	    if(i == 1) {
	      i = value;
	      float diff = fabs(value-i);
	      if(diff < 1e-3) {
		if(verbose.verbose) printf("  file is in AO ASCII format 2.\n");
		fclose(fin);
		free(txt);
		return AO_ASCII_2_format;
	      }else {
		if(verbose.verbose) {
		  for(indent = 0; indent < verbose.indent; indent++)
		    printf(" ");
		}
		printf("  file is not AO ASCII format 2 (expected 8th word to be an integer: diff with integer is %f).\n", diff);
	      }
	    }else {
	      if(verbose.verbose) {
		for(indent = 0; indent < verbose.indent; indent++)
		  printf(" ");
	      }
	      printf("  file is not AO ASCII format 2 (expected 8th word to be a value).\n");
	    }
	  }else {
	    if(verbose.verbose) {
	      for(indent = 0; indent < verbose.indent; indent++)
		printf(" ");
	    }
	    printf("  file is not AO ASCII format 2 (expected %f to be a MJD).\n", value);
	  }
	}else {
	  if(verbose.verbose) {
	    for(indent = 0; indent < verbose.indent; indent++)
	      printf(" ");
	  }
	  printf("  file is not AO ASCII format 2 (expected 2nd word to be a value).\n");
	}
      }else { // Nr words test failed
	if(verbose.verbose) {
	  for(indent = 0; indent < verbose.indent; indent++)
	    printf(" ");
	}
	printf("  file is not AO ASCII format 2 (expected 11 words on first line, got %d).\n", nrwords);
      }
    }else {  // if strlen failed
      if(verbose.verbose) {
	for(indent = 0; indent < verbose.indent; indent++)
	  printf(" ");
      }
      printf("  file is not AO ASCII format 2 (abnormal long first line).\n");
    }
    rewind(fin);
    i = fread(txt, 1, 3, fin);   /* Leave things they were before the test */
    txt[3] = 0;
  }else if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    printf("  file is not AO ASCII format 2.\n");
  }

//START REGION RELEASE

  // Check if in timer format (warning only)
  int istimer = 1;
  if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    printf("  check if in TIMER format.\n");
  }
  if(fseek(fin, 32+32+8, SEEK_SET) != 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR guessPSRData_format: Seek in file '%s' failed.", filename);
    free(txt);
    return -1;
  }
  i = fread(&fdummy, sizeof(float), 1, fin);
  if(i != 1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR guessPSRData_format: reading file '%s' failed.", filename);
    free(txt);
    return -1;
  }
  if(verbose.debug) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
    }
    printf("    got version = %f\n", fdummy);
  }
  if(fabs(fdummy) < 1e-10) {
    i = fread(&fdummy, sizeof(float), 1, fin);
    if(i != 1) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR guessPSRData_format: reading file '%s' failed.", filename);
      free(txt);
      return -1;
    }
    if(verbose.debug) {
      if(verbose.verbose) {
	for(indent = 0; indent < verbose.indent; indent++)      
	  printf(" ");
      }
      printf("    got minor version = %f\n", fdummy);
    }
    if(fabs(fdummy) < 1e-10) {
      if(fseek(fin, 2*sizeof(int), SEEK_CUR) != 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR guessPSRData_format: Seek in file '%s' failed.", filename);
	free(txt);
	return -1;
      }
      i = fread(txt, 1, 16, fin);
      if(i != 16) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR guessPSRData_format: reading file '%s' failed.", filename);
	free(txt);
	return -1;
      }
      txt[16] = 0;
      if(verbose.debug) {
	if(verbose.verbose) {
	  for(indent = 0; indent < verbose.indent; indent++)      
	    printf(" ");
	}
	printf("    got date = '%s'\n", txt);
      }
      if(txt[2] != '-' || txt[5] != '-') {
	istimer = 0;
      }
    }else {
      istimer = 0;
    }
  }else {
    istimer = 0;
  }



  fflush(stdout);
  fclose(fin);
  free(txt);
  if(istimer) {
    printerror(verbose.debug, "ERROR guessPSRData_format: file '%s' appears to be in TIMER format.", filename);
    printerror(verbose.debug, "ERROR guessPSRData_format: If you have psrchive, you can use 'pam -a psrfits -e fits %s' to convert to PSRFITS.", filename);
    return -2;
  }else {
    if(noerror == 0) {
      printerror(verbose.debug, "ERROR guessPSRData_format: determining file type of '%s' failed.", filename);
    }else {
      if(verbose.verbose) {
	for(indent = 0; indent < verbose.indent; indent++)      
	  printf(" ");
	printf("  file is in an undetermined format.\n");
      }
    }
  }
  return 0;
}

//START REGION DEVELOP

//START REGION RELEASE
/* Open filename which is in the specified format. Set the
   enable_write flag to open the file with write permissions. If
   enable_write is not set, cleanPSRData() is called, meaning that if
   this is done in the main progamme, there will be a memory leak
   since the struct will be re-initialised.  read_in_memory is set the
   file is read into memory and the file is closed. You can acces the
   data as if the file is still open. If the format is set to zero it
   will try to guess the format. If enable_write is not set, the
   datafile structure is always cleared before any reading/writing is
   done. If nowarnings == 0, warnings are shown. If 1, no warnings are
   generated for incomplete headers (the checks done AFTER the
   file-format specific reading has been done), which can be useful
   for applications which do not rely on the header information. If 2,
   all warnings are ignored while reading in the header. If nocounters
   is set it will not write out any counters (useful for log
   files). Verbose level determines nr of spaces before output. If
   debug is set, more information is printed. Returns 0 on error. */
int openPSRData(datafile_definition *datafile, char *filename, int format, int enable_write, int read_in_memory, int nowarnings, verbose_definition verbose)
{
  int status = 0, iomode, i;   /* CFITSIO status value MUST be initialized to zero! */
  char open_mode[100], filename2[1000];
  verbose_definition verbose2;
//START REGION DEVELOP
  char hdrname[1000];
//START REGION RELEASE


//  fprintf(stderr, "XXXXX open=%d format=%d\n", datafile->opened_flag, format);

  //Forget everything (header params etc) which might be set by a previous file which has been read in
  if(enable_write == 0)
    cleanPSRData(datafile, verbose);

  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;

  if(verbose.verbose) {
    for(i = 0; i < verbose.indent; i++)      
      printf(" ");
    printf("Opening file '%s'\n", filename);
  }
  if(set_filename_PSRData(datafile, filename, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR openPSRData: Making copy of data parameters failed.");
    return 0;
  }
  datafile->opened_flag = 0;
  /*  printf("Opening format %s\n", returnFileFormat_str(format)); */
  datafile->format = format;
  if(enable_write) {
    datafile->enable_write_flag = 1;
  }

  if(format <= 0) {
    if(access(filename, F_OK) != 0) {
      fflush(stdout);
      if(enable_write == 0) {
	printerror(verbose.debug, "ERROR openPSRData: Cannot open file %s for reading.", filename);
      }else {
	printerror(verbose.debug, "ERROR openPSRData: File %s requested to be opened in write mode without file format being specified. File doesn't exist yet, so file format cannot be determined from existing file.", filename);
      }
      return 0;	
    }
    format = guessPSRData_format(filename, 0, verbose2);
  }
  if(format == -1) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR openPSRData:\n  Error reading file '%s' - file format couldn't be determined.", filename);
    return 0;
  }
  if(format == -2) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR openPSRData:\n  The file type of '%s' isn't supported for reading in.", filename);
    return 0;
  }
  if(format == 0 || format == -3) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR openPSRData:\n  Error determining file type of '%s'. It is either not supported, or you might need to specify the format on the command line.", filename);
    return 0;
  }
  datafile->format = format;



  open_mode[0] = 0;
  if(enable_write) {
    /* Also allow reading which is required when writing out the history (at least in puma format) */
    strcat(open_mode, "w+");
    if(verbose.debug) {
      if(verbose.verbose) {
	for(i = 0; i <= verbose2.indent; i++)      
	  printf(" ");
      }
      printf("Opening file '%s' for writing\n", filename);
    }
    if(format == PSRCHIVE_ASCII_format) {
      if(verbose.debug) {
	if(verbose.verbose) {
	  for(i = 0; i <= verbose2.indent; i++)      
	    printf(" ");
	}
	printf("  The write commands will be buffered and executed once the file is being closed.\n");
      }
      datafile->dumpOnClose = 1; // The writing itself should be done in writePSRData() on closing, but not if that function was used by program rather than writing individual pulses, although, that is calling write individual pulses, so should be ok. However, make sure that data does not exist in memory twice!
    }else {
      datafile->dumpOnClose = 0;
    }
  }else {
    strcat(open_mode, "r");
    if(verbose.debug) {
      if(verbose.verbose) {
	for(i = 0; i <= verbose2.indent; i++)      
	  printf(" ");
      }
      printf("Opening file '%s' for reading\n", filename);
    }
  }

  /* WSRT, Presto and Parkes data are binary formats */
  if(format == PUMA_format || format == PSRSALSA_BINARY_format ||
//START REGION DEVELOP
     format == PRESTO_format || format == PARKESFB_format || 
//START REGION RELEASE
     format == SIGPROC_format) {
    strcat(open_mode, "b");
    datafile->fptr = fopen(filename, open_mode);
    if(datafile->fptr == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR openPSRData: Cannot open file %s: %s", filename, strerror(errno));
      datafile->opened_flag = 0;
    }else {
      datafile->opened_flag = 1;
    }
  }else if(format == PSRCHIVE_ASCII_format || format == PPOL_format || format == PPOL_SHORT_format || format == SIGPROC_ASCII_format || format == EPN_format
//START REGION DEVELOP
 || format ==  AO_ASCII_1_format || format ==  AO_ASCII_2_format || format == GMRT_ASCII_format
//START REGION RELEASE
       ) {
    datafile->fptr = fopen(filename, open_mode);
    if(datafile->fptr == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR openPSRData: Cannot open file %s: %s", filename, strerror(errno));
      datafile->opened_flag = 0;
    }else {
      datafile->opened_flag = 1;
    }
  }else if(format == FITS_format) {
      /* Table 6 should be the subintegrations */
    if(enable_write) {
      iomode = READWRITE;
      /* exlamation mark caused file to be overwritten. */
      sprintf(filename2, "!%s", filename);
      if (!fits_create_file(&(datafile->fits_fptr), filename2, &status)) {
	if(verbose.debug) {
	  if(verbose.verbose) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	  }
	  printf("  '%s' opened\n", filename);
	}
	datafile->opened_flag = 1;
      }else {
	datafile->opened_flag = 0;
      }
    }else {
      strcpy(filename2, filename);
      iomode = READONLY;
      if (!fits_open_file(&(datafile->fits_fptr), filename2, iomode, &status)) {
	if(verbose.debug) {
	  if(verbose.verbose) {
	    for(i = 0; i < verbose.indent; i++)      
	      printf(" ");
	  }
	  printf("  '%s' opened\n", filename);
	}
	datafile->opened_flag = 1;
      }else {
	datafile->opened_flag = 0;
      }
    }
    if(status) {
      fflush(stdout);
      fits_report_error(stderr, status); /* print any error message */
    }
  }

//START REGION DEVELOP
  /* PRESTO and Parkes uses separate header file */
  if(format == PRESTO_format || format == PARKESFB_format) {   
    strcpy(hdrname, filename);
    hdrname[strlen(hdrname)-4] = 0;
    if(format == PRESTO_format)
      strcat(hdrname, ".inf");
    else if(format == PARKESFB_format)
      strcat(hdrname, ".hdr");
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    if(verbose.verbose) printf("  Opening header file %s\n", hdrname);
    if(enable_write)
      datafile->fptr_hdr = fopen(hdrname, "w+");
    else
      datafile->fptr_hdr = fopen(hdrname, "r");
    if(datafile->fptr_hdr == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR openPSRData: Cannot open file %s: %s", hdrname, strerror(errno));
      return 0;
    }
  }else {
//START REGION RELEASE
    datafile->fptr_hdr = datafile->fptr;
//START REGION DEVELOP
  }
//START REGION RELEASE
  switch(format) {
    //  case PUMA_format: datafile->datastart = sizeof(Header_type); break;
  case PUMA_format: datafile->datastart = 4504; break;
  default: 
    datafile->datastart = 0; 
    break;
  }

//START REGION RELEASE
  if(read_in_memory && datafile->opened_flag) {
//START REGION DEVELOP
      /* EPN format no longer treated differently. Used to be that:
	 EPN format is handled differently, because the header should
	 not be read separately and the memory should not be allocated
	 yet and the file should be close. */
    /*    if(datafile->format == EPN_format) {
      fclose(datafile->fptr);
      if(readEPNfile(datafile, verbose, -1)) {
	datafile->format = MEMORY_format;
	datafile->opened_flag = 1;
      }else {
    fflush(stdout);
	printerror(verbose.debug, "ERROR openPSRData: Cannot read data.");
	return 0;
      }
    }else 
    */
//START REGION RELEASE
    if(readHeaderPSRData(datafile, 0, nowarnings, verbose2)) {
      if(datafile->NrPols != 0) {
	long datasize = datafile->NrSubints*datafile->NrBins*datafile->NrPols*datafile->NrFreqChan*sizeof(float);
	datafile->data = (float *)malloc(datasize);
	//	fprintf(stderr, "\nXXXXX Allocating %ld bytes for data\n", datafile->NrSubints*datafile->NrBins*datafile->NrPols*datafile->NrFreqChan*sizeof(float));
	if(datafile->data == NULL) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR openPSRData: Cannot allocate memory (data=%ld bytes=%.3fGB).", datasize, datasize/1073741824.0);
	  closePSRData(datafile, 0, verbose2);
	  return 0;
	}
      }
//START REGION DEVELOP
/*
      if(format == PPOL_format || format == PPOL_SHORT_format) {
	datafile->data_pa = malloc(datafile->NrSubints*datafile->NrBins*datafile->NrFreqChan*sizeof(float)); 
	datafile->data_dpa = malloc(datafile->NrSubints*datafile->NrBins*datafile->NrFreqChan*sizeof(float));
	datafile->data_long = malloc(datafile->NrBins*sizeof(float));
	//	fprintf(stderr, "\nXXXXX Allocating %ld bytes for data\n", datafile->NrSubints*datafile->NrBins*datafile->NrFreqChan*sizeof(float));
	//	fprintf(stderr, "\nXXXXX Allocating %ld bytes for data\n", datafile->NrSubints*datafile->NrBins*datafile->NrFreqChan*sizeof(float));
	//	fprintf(stderr, "\nXXXXX Allocating %ld bytes for data\n", datafile->NrBins*sizeof(float));
	if(datafile->rmspoints_defined) {
	  datafile->rmsvalues = malloc(datafile->NrSubints*datafile->rmspoints_defined*datafile->NrFreqChan*sizeof(float)); 
	  //	  fprintf(stderr, "\nXXXXX Allocating %ld bytes for data\n", datafile->NrSubints*datafile->rmspoints_defined*datafile->NrFreqChan*sizeof(float));
	  if(datafile->rmsvalues == NULL) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR openPSRData: Cannot allocate memory (rmsvalues).");
	    closePSRData(datafile, verbose2);
	    return 0;
	  }
	}
	if(datafile->data_pa == NULL || datafile->data_dpa == NULL || datafile->data_long == NULL) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR openPSRData: Cannot allocate memory (pa/dpa/longitudes).");
	  closePSRData(datafile, verbose2);
	  return 0;
	}
	if(format == PPOL_format) {
	  datafile->data_bin = malloc(datafile->NrBins*sizeof(long int));
	  //	  fprintf(stderr, "\nXXXXX Allocating %ld bytes for data\n", datafile->NrBins*sizeof(long int));
	  if(datafile->data_bin == NULL) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR openPSRData: Cannot allocate memory (bins).");
	    closePSRData(datafile, verbose2);
	    return 0;
	  }
	}
      }
*/
//START REGION RELEASE
      if(readPSRData(datafile, datafile->data, verbose2)) {
	closePSRData(datafile, 2, verbose2);  // Perserve header information+data, as the data will be available in MEMORY_format
	datafile->format = MEMORY_format;
	datafile->opened_flag = 1;

      }else {
	fflush(stdout);
	printerror(verbose.debug, "ERROR openPSRData: Cannot read data.");
	closePSRData(datafile, 0, verbose2);
	return 0;
      }
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR openPSRData: Cannot read header.");
      closePSRData(datafile, 0, verbose2);
      return 0;
    }
  }

  return datafile->opened_flag;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Close the files and free the memory if required. Data is written to
   file if writing is buffered (useful for ascii formats). To be
   specific:

   - If requested by the dumpOnClose flag, the data is written to the file and memory is released.
   - If the file is opened it is closed.
   - If the format is in MEMORY_format the allocated memory (data) is released
   - If perserve_info is 1, memory associated with things like the pulsar name will not be freed. If set to 2, also the data memory is not released.

Returns 0 on success. 
*/
int closePSRData(datafile_definition *datafile, int perserve_info, verbose_definition verbose)
{
  // If freeing memory for history & strings, it will cause a problem for swap_orig_clone.
  int indent;
  int status = 0;

  if(verbose.debug) {
    printf("Closing file '%s'\n", datafile->filename);
  }

  
  if(datafile->opened_flag) {
    if(datafile->dumpOnClose) {
      if(verbose.debug) {
	printf("  - Dumping data to file before closing\n");
      }
      if(verbose.verbose) {
	for(indent = 0; indent < verbose.indent; indent++)      
	  printf(" ");
	printf("Writing buffer to file %s before it is closed\n", datafile->filename);
      }
      if(datafile->data == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR closePSRData: Although output is buffered, no memory is allocated?");
	return 1;
      }
      if(writePSRData(datafile, datafile->data, verbose) != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR closePSRData: Writing of buffered data failed.");
	return 1;
      }
    }
//START REGION DEVELOP
    if(datafile->format == PRESTO_format || datafile->format == PARKESFB_format) {
      if(verbose.debug) {
	printf("  - Releasing file pointer to separate header file\n");
      }
      fclose(datafile->fptr_hdr);
    }
//START REGION RELEASE
    if(datafile->format == FITS_format) {
      if(verbose.debug) {
	printf("  - Releasing FITS file pointer\n");
      }
      fits_close_file(datafile->fits_fptr, &status);
      if(status) {
	fflush(stdout);
	fits_report_error(stderr, status); /* print any error message */
      }
      if(verbose.debug) {
	printf("  - Releasing memory related to scales/offsets/weights\n");
      }
      free(datafile->scales);
      free(datafile->offsets);
      free(datafile->weights);
      datafile->scales = NULL;
      datafile->offsets = NULL;
      datafile->weights = NULL;
    }else if(datafile->format != MEMORY_format){
      if(verbose.debug) {
	printf("  - Releasing file pointer\n");
      }
      fclose(datafile->fptr);
    }
    datafile->opened_flag = 0;
  }   // If file is not open
   
  if(perserve_info != 2) {
    if(datafile->data != NULL) {
      if(verbose.debug) {
	printf("  - Releasing memory containing data\n");
      }
      //      printf("XXXXXXX releasing %p\n", datafile->data);
      free(datafile->data);
      datafile->data = NULL;  // Important if datafile struct is re-used. Needs to forget that memory was allocated in the past.
    }
  }

  if(perserve_info == 0) {
    if(verbose.debug) {
      printf("  - Releasing header related memory\n");
    }
    //    printf("XXXXX free %p\n", datafile->filename);
    free(datafile->filename); // Is always allocated, at least when cleanPSRData() was used
    datafile->filename = NULL;

    free(datafile->psrname);
    datafile->psrname = NULL;
    free(datafile->observatory);
    datafile->observatory = NULL;
    free(datafile->institute);
    datafile->institute = NULL;
    free(datafile->instrument);
    datafile->instrument = NULL;
    free(datafile->scanID);
    datafile->scanID = NULL;
    free(datafile->observer);
    datafile->observer = NULL;
    free(datafile->projectID);
    datafile->projectID = NULL;

    if(datafile->tsamp_list != NULL) {
      free(datafile->tsamp_list); 
      datafile->tsamp_list = NULL;
    }
    if(datafile->tsub_list != NULL) {
      //      printf("XXXXX free %p\n", datafile->tsub_list);
      free(datafile->tsub_list); 
      datafile->tsub_list = NULL;
    }
    if(datafile->freqlabel_list != NULL) {
      //      printf("XXXXX free %p\n", datafile->freq_list);
      free(datafile->freqlabel_list); 
      datafile->freqlabel_list = NULL;
    }
    if(datafile->offpulse_rms != NULL) {
      free(datafile->offpulse_rms);
      datafile->offpulse_rms = NULL;
    }

    datafile_history_entry_definition *history_ptr, *history_ptr_next;
    int firsthistoryline;
    history_ptr = &(datafile->history);
    firsthistoryline = 1;
    do {
      if(history_ptr->timestamp != NULL) {
	free(history_ptr->timestamp);
	history_ptr->timestamp = NULL;
      }
      if(history_ptr->cmd != NULL) {
	free(history_ptr->cmd);
	history_ptr->cmd = NULL;
      }
      if(history_ptr->user != NULL) {
	free(history_ptr->user);
	history_ptr->user = NULL;
      }
      if(history_ptr->hostname != NULL) {
	free(history_ptr->hostname);
	history_ptr->hostname = NULL;
      }
      history_ptr_next = history_ptr->nextEntry;
      history_ptr->nextEntry = NULL;
      if(firsthistoryline == 1) {
	firsthistoryline = 0;   // First line is not dynamically allocated
      }else {
	free(history_ptr);
      }
      history_ptr = history_ptr_next;
    }while(history_ptr_next != NULL);

  }
 
  return status;
}

//START REGION DEVELOP
//START REGION RELEASE

static char * internal_gentype_string_undefined = "Not set";
static char * internal_gentype_string_profile = "Profile";
static char * internal_gentype_string_pulsestack = "Pulsestack";
static char * internal_gentype_string_subints = "Subints";
//static char * internal_gentype_string_freqphase = "Freq. vs. phase";
static char * internal_gentype_string_searchmode = "Search mode";
static char * internal_gentype_string_bandpass = "Bandpass";
static char * internal_gentype_string_dynspec = "Dynamic spectrum";
static char * internal_gentype_string_polcalib = "Pol. calibrator";
static char * internal_gentype_string_lrfs = "LRFS";
static char * internal_gentype_string_2dfs = "2DFS";
static char * internal_gentype_string_s2dfsp3 = "S2DFS P3 map";
static char * internal_gentype_string_s2dfsp2 = "S2DFS P2 map";
static char * internal_gentype_string_p3fold = "P3 fold";
static char * internal_gentype_string_polarmap = "Polar map";
static char * internal_gentype_string_lrcc = "LRCC";
static char * internal_gentype_string_lrac = "LRAC";
static char * internal_gentype_string_padist = "PA distr.";
static char * internal_gentype_string_elldist = "Ell distr.";
//static char * internal_gentype_string_ILVPA = "I,L,V,PA data";
static char * internal_gentype_string_recmodel = "Receiver model";
static char * internal_gentype_string_recmodel2 = "Receiver model with chi^2";
static char * internal_gentype_string_rmmap = "RM map";
static char * internal_gentype_string_penergy = "penergy output";
static char * internal_gentype_string_hrfs_unfolded = "HRFS (unfolded)";
static char * internal_gentype_string_hrfs = "HRFS";
static char * internal_gentype_string_bug = "BUG, undefined????";

/* return a string with the description of the gentype */
char *returnGenType_str(int gentype)
{
  switch(gentype) {
  case GENTYPE_UNDEFINED: return internal_gentype_string_undefined;  break;
  case GENTYPE_PROFILE: return internal_gentype_string_profile; break;
  case GENTYPE_PULSESTACK: return internal_gentype_string_pulsestack; break;
  case GENTYPE_SUBINTEGRATIONS: return internal_gentype_string_subints; break;
    //  case GENTYPE_FREQPHASE: return internal_gentype_string_freqphase; break;
  case GENTYPE_SEARCHMODE: return internal_gentype_string_searchmode; break;
  case GENTYPE_BANDPASS: return internal_gentype_string_bandpass; break;
  case GENTYPE_DYNAMICSPECTRUM: return internal_gentype_string_dynspec; break;
  case GENTYPE_POLNCAL: return internal_gentype_string_polcalib; break;
  case GENTYPE_LRFS: return internal_gentype_string_lrfs; break;
  case GENTYPE_2DFS: return internal_gentype_string_2dfs; break;
  case GENTYPE_S2DFSP3: return internal_gentype_string_s2dfsp3; break;
  case GENTYPE_S2DFSP2: return internal_gentype_string_s2dfsp2; break;
  case GENTYPE_P3FOLD: return internal_gentype_string_p3fold; break;
  case GENTYPE_POLARMAP: return internal_gentype_string_polarmap; break;
  case GENTYPE_LRCC: return internal_gentype_string_lrcc; break;
  case GENTYPE_LRAC: return internal_gentype_string_lrac; break;
  case GENTYPE_PADIST: return internal_gentype_string_padist; break;
  case GENTYPE_ELLDIST: return internal_gentype_string_elldist; break;
    //  case GENTYPE_ILVPA: return internal_gentype_string_ILVPA; break;
  case GENTYPE_RECEIVERMODEL: return internal_gentype_string_recmodel; break;
  case GENTYPE_RECEIVERMODEL2: return internal_gentype_string_recmodel2; break;
  case GENTYPE_RMMAP: return internal_gentype_string_rmmap; break;
  case GENTYPE_PENERGY: return internal_gentype_string_penergy; break;
  case GENTYPE_HRFS_UNFOLDED: return internal_gentype_string_hrfs_unfolded; break;
  case GENTYPE_HRFS: return internal_gentype_string_hrfs; break;
  default: return internal_gentype_string_bug; break;
  }
}

/* Print the description of the gentype to the destination (e.g. stdout) */
void printGenType(int gentype, FILE *destination)
{
  fprintf(destination, "%s", returnGenType_str(gentype));
}


//START REGION DEVELOP
static char * internal_format_string_AOascii1      = "AO ascii type 1";
static char * internal_format_string_AOascii2      = "AO ascii type 2";
static char * internal_format_string_Presto        = "Presto";
static char * internal_format_string_ParkesFB      = "ParkesAFB";
static char * internal_format_string_GMRTAscii     = "GMRT ascii";
//START REGION RELEASE
static char * internal_format_string_PuMa          = "PuMa";
static char * internal_format_string_EPN           = "EPN";
static char * internal_format_string_PPOL          = "PAswing";
static char * internal_format_string_PPOL_short    = "PAswing (short)";
static char * internal_format_string_PSRCHIVEAscii = "PSRCHIVE ascii";
static char * internal_format_string_PSRFITS       = "PSRfits";
static char * internal_format_string_SIGPROC       = "Sigproc";
static char * internal_format_string_SIGPROCAscii  = "Sigproc (ascii)";
static char * internal_format_string_PSRSALSA      = "PSRSALSA binary";
static char * internal_format_string_Memory        = "Loaded in RAM";
static char * internal_format_string_bug           = "BUG, undefined????";

/* return a string with the description of the file format */
char *returnFileFormat_str(int format)
{
  switch(format) {
  case PUMA_format: return internal_format_string_PuMa;  break;
//START REGION DEVELOP
  case AO_ASCII_1_format: return internal_format_string_AOascii1;  break;
  case AO_ASCII_2_format: return internal_format_string_AOascii2;  break;
  case PRESTO_format: return internal_format_string_Presto;  break;
  case PARKESFB_format: return internal_format_string_ParkesFB;  break;
//START REGION RELEASE
  case EPN_format: return internal_format_string_EPN;  break;
  case SIGPROC_format: return internal_format_string_SIGPROC;  break;
  case PPOL_format: return internal_format_string_PPOL;  break;
  case PPOL_SHORT_format: return internal_format_string_PPOL_short;  break;
  case SIGPROC_ASCII_format: return internal_format_string_SIGPROCAscii;  break;
//START REGION DEVELOP
  case GMRT_ASCII_format: return internal_format_string_GMRTAscii;  break;
//START REGION RELEASE
  case PSRCHIVE_ASCII_format: return internal_format_string_PSRCHIVEAscii;  break;
  case FITS_format: return internal_format_string_PSRFITS;  break;
  case PSRSALSA_BINARY_format: return internal_format_string_PSRSALSA;  break;
  case MEMORY_format: return internal_format_string_Memory;  break;
  default: return internal_format_string_bug; break;
  }
}

/* Automatically done in readHeaderPSRData if verbose is set. If
   update is set, the first line indicates these are updated header
   commands (used after using -header command line). Indent indicates
   the nr of spaces before output. If debug is set, the output is
   more precise, but less concise. */
void printHeaderPSRData(datafile_definition datafile, int update, verbose_definition verbose)
{
  int i;
  char txt[1000];
  for(i = 0; i < verbose.indent; i++) printf(" ");
  if(update == 0)
    printf("========================= Dump of header parameters =========================\n");
  else
    printf("===================== Dump of updated header parameters =====================\n");
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("Pulsar=%s  ", datafile.psrname);
  printf("Ra=");
  converthms_string(txt, datafile.ra*12.0/M_PI, 2, 2);
  printf("%s ", txt);
  if(verbose.debug) {
    printf("= %f rad = %f deg = %f hours \n", datafile.ra, datafile.ra*180.0/M_PI, datafile.ra*12.0/M_PI);
    for(i = 0; i < verbose.indent; i++) printf(" ");
  }
  printf("Dec=");
  converthms_string(txt, datafile.dec*180.0/M_PI, 2, 3);
  printf("%s  ", txt);
  if(verbose.debug) {
    printf(" = %f rad = %f deg\n", datafile.dec, datafile.dec*180.0/M_PI);
    for(i = 0; i < verbose.indent; i++) printf(" ");
  }
  printf("\n");
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("tobs=%.1f sec", get_tobs(datafile, verbose));
  if(datafile.tsubMode == TSUBMODE_FIXEDTSUB) {
    printf(" (fixed subint length)\n");
  }else if(datafile.tsubMode == TSUBMODE_TSUBLIST) {
    printf(" (variable subint length)\n");
  }else {
    printf(" (unknown subint length)\n");
  }
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("DM=%f  RM=%f", datafile.dm, datafile.rm);
  printf(" RefFreq=");
  if((datafile.freq_ref < -0.9 && datafile.freq_ref >= -1.1) || (datafile.freq_ref > 0.99e10 && datafile.freq_ref < 1.01e10)){
    printf("infinity");
  }else if(datafile.freq_ref >= 0) {
    if(verbose.debug == 0)
      printf("%f MHz", datafile.freq_ref);
    else
      printf("%.9e MHz", datafile.freq_ref);
  }else {
    printf("Unknown");
  }
  printf("\n");
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("de-dispersed=");
  if(datafile.isDeDisp == 1)
    printf("yes");
  else if(datafile.isDeDisp == 0)
    printf("no");
  else
    printf("unknown");    
  printf(" de-Faraday=");
  if(datafile.isDeFarad == 1)
    printf("yes");
  else if(datafile.isDeFarad == 0)
    printf("no");
  else
    printf("unknown");    
  printf(" de-par. angle=");
  if(datafile.isDePar == 1)
    printf("yes");
  else if(datafile.isDePar == 0)
    printf("no");
  else
    printf("unknown");    
  printf(" de-baselined=");
  if(datafile.isDebase == 1)
    printf("yes\n");
  else if(datafile.isDebase == 0)
    printf("no\n");
  else
    printf("unknown\n");
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("NrSubints=%ld NrBins=%ld NrPols=%ld NrFreqChan=%ld NrBits=%d\n",datafile.NrSubints,datafile.NrBins,datafile.NrPols, datafile.NrFreqChan,datafile.NrBits);
  for(i = 0; i < verbose.indent; i++) printf(" ");
  mjd2dateString(datafile.mjd_start, txt, 0, 1, " ");

  if(datafile.isFolded == 1) {
    double period;
    if(get_period(datafile, 0, &period, verbose) == 2) {
      printerror(verbose.debug, "ERROR printHeaderPSRData: Cannot obtain period");
      exit(0);
    }
    if(verbose.debug == 0)
      printf("period=%lf sec", period);
    else
      printf("period=%.9e sec", period);
  }else {
    printf("Not folded");
  }
  if(datafile.tsampMode == TSAMPMODE_LONGITUDELIST) {
    printf(" SampTime=longitude array");
  }else {
    if(verbose.debug == 0)
      printf(" SampTime=%f sec", get_tsamp(datafile, 0, verbose));
    else
      printf(" SampTime=%.9e sec", get_tsamp(datafile, 0, verbose));
  }
  if(verbose.debug == 0)
    printf(" mjd=%.9Lf (%s)\n", datafile.mjd_start, txt);
  else
    printf(" mjd=%.9Lf (%s)\n", datafile.mjd_start, txt);
  for(i = 0; i < verbose.indent; i++) printf(" ");
  if(datafile.freqMode == FREQMODE_UNKNOWN) {
    printwarning(verbose.debug, "OBSERVING FREQUENCY IS NOT SET");
  }else {
    double bw, chanbw, freq_cent;
    bw = get_bandwidth(datafile, verbose);
    if(get_channelbandwidth(datafile, &chanbw, verbose) == 0) {
      printerror(verbose.debug, "ERROR printHeaderPSRData: Cannot obtain channel bandwidth");
      exit(0);
    }
    freq_cent = get_centre_frequency(datafile, verbose);
    if(verbose.debug == 0)
      printf("freq_cent=%f MHz bw=%f MHz channelbw=%f MHz\n", freq_cent, bw, chanbw);
    else
      printf("freq_cent=%.9e MHz bw=%.9e MHz channelbw=%.9e MHz\n", freq_cent, bw, chanbw);
  }
  for(i = 0; i < verbose.indent; i++) printf(" ");
  if(verbose.debug) {
    printf("observatory=%s  ITRF (X,Y,Z)=(%lf,%lf,%lf) m\n", datafile.observatory, datafile.telescope_X, datafile.telescope_Y, datafile.telescope_Z);
    for(i = 0; i < verbose.indent; i++) printf(" ");
    printf("ITRF derived geocentric long=%lf deg lat=%lf deg\n", observatory_long_geocentric(datafile)*180.0/M_PI, observatory_lat_geocentric(datafile)*180.0/M_PI);
    for(i = 0; i < verbose.indent; i++) printf(" ");
    double longitude, latitude, height;
    tempo2_ITRF_to_GRS80(datafile.telescope_X, datafile.telescope_Y, datafile.telescope_Z, &longitude, &latitude, &height);
    printf("GRS80 derived geodetic  long=%lf deg lat=%lf deg height=%lf m\n", longitude*180.0/M_PI, latitude*180.0/M_PI, height);
  }else {
    printf("observatory=%s  long=%lf deg lat=%lf deg (geodetic derived)\n", datafile.observatory, observatory_long_geodetic(datafile)*180.0/M_PI, observatory_lat_geodetic(datafile)*180.0/M_PI);
  }
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("projectID=%s  observer=%s  institute=%s  instrument=%s\n", datafile.projectID, datafile.observer, datafile.institute, datafile.instrument);
//START REGION DEVELOP
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("Cable swap during obs: ");
  if(datafile.cableSwap == 0) {
    printf("No ");
  }else if(datafile.cableSwap == 1) {
    printf("Yes ");
  }else {
    printf("Unknown ");
  }
  printf("  correction applied: ");
  if(datafile.cableSwapcor == 0) {
    printf("No\n");
  }else if(datafile.cableSwapcor == 1) {
    printf("Yes\n");
  }else {
    printf("Unknown\n");
  }
//START REGION RELEASE
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("feedtype=%d ", datafile.feedtype);
  if(datafile.feedtype == FEEDTYPE_LINEAR) {
    printf("(positive linear)");
  }else if(datafile.feedtype == FEEDTYPE_INV_LINEAR) {
    printf("(negative linear)");
  }else if(datafile.feedtype == FEEDTYPE_CIRCULAR) {
    printf("(positive circular)"); 
  }else if(datafile.feedtype == FEEDTYPE_INV_CIRCULAR) {
    printf("(negative circular)");
  }else {
    printf("(Unknown)");
  }
  printf("  poltype=%d ", datafile.poltype);
  if(datafile.poltype == POLTYPE_UNKNOWN) {
    printf("(undefined)\n");
  }else if(datafile.poltype == POLTYPE_STOKES) {
    printf("(Stokes)\n");
  }else if(datafile.poltype == POLTYPE_COHERENCY) {
    printf("(coherency)\n");
  }else if(datafile.poltype == POLTYPE_ILVPAdPA) {
    printf("(I,L,V,Pa and error)\n");
  }else if(datafile.poltype == POLTYPE_ILVPAdPATEldEl) {
    printf("(I,L,V,Pa+error,tot pol,ell+error)\n");
  }else if(datafile.poltype == POLTYPE_PAdPA) {
    printf("(Pa and error)\n");
  }else {
    printf("BUG!!!!\n");
  }
//START REGION DEVELOP
//  for(i = 0; i < verbose.indent; i++) printf(" ");
//  printf("fd_sang=%f  fd_xyph=%f\n", datafile.fd_sang, datafile.fd_xyph);
//START REGION RELEASE
//  for(i = 0; i < verbose.indent; i++) printf(" ");
  /*  printf("paswing=");
  if(!needstoberemoved)
    printf("YES");
  else
    printf("NO");
  printf("  longitudes=");
  if(!needstoberemoved)
    printf("YES");
  else
    printf("NO");
  printf("  bins=");
  if(!needstoberemoved)
    printf("YES  ");
  else
    printf("NO  ");
  */
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("FileName=");
  if(datafile.filename == NULL)
    printf(" ");
  else
    printf("%s", datafile.filename);
  printf("  ScanID=%s  Header length = %ld bytes\n", datafile.scanID, (long)datafile.datastart);
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("GenType=");
  printGenType(datafile.gentype, stdout);
  if(datafile.xrangeset)
    printf("  xrange=%f %f", datafile.xrange[0], datafile.xrange[1]);
  if(datafile.yrangeset)
    printf("  yrange=%f %f", datafile.yrange[0], datafile.yrange[1]);
  printf(" rms values=");
  if(datafile.offpulse_rms != NULL)
    printf("YES\n");
  else
    printf("NO\n");
  for(i = 0; i < verbose.indent; i++) printf(" ");
  printf("====================== End of dump of header parameters ======================\n");
}

//START REGION DEVELOP
//START REGION RELEASE

// Set the ITRF telescope location XYZ by name of the observatory.
// If verbose is set, the location is printed, with verbose-1 spaces in front of output.
// Return 1 if observatory was recognized, or 0 if not.
int setITRFlocation_by_name(datafile_definition *datafile, char *observatory, verbose_definition verbose)
{
  int indent;
  if(strcasecmp(observatory, "PARKES") == 0 || strcasecmp(observatory, "PKS") == 0 || strcasecmp(observatory, "7") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this Parkes telescope data");
    }
    //      datafile->telescope_lat = M_PI*(-32.99994444444444444444)/180.0;
    //      datafile->telescope_long = M_PI*(148+15/60.0+48.636/3600.0)/180.0;
    // These coordinates are from a fits file
    //      datafile->telescope_X = -4554231.6;
    //      datafile->telescope_Y = 2816759.1;
    //      datafile->telescope_Z = -3454036.1;
    // These coordinates are from tempo2
    datafile->telescope_X = -4554231.5;
    datafile->telescope_Y = 2816759.1;
    datafile->telescope_Z = -3454036.3;
  }else if(strcasecmp(observatory, "JODRELL") == 0 || strcasecmp(observatory, "JB") == 0 || strcasecmp(observatory, "JBO") == 0 || strcasecmp(observatory, "LOVELL") == 0 || strcasecmp(observatory, "JBDFB") == 0 || strcasecmp(observatory, "JBDFB") == 0 || strcasecmp(observatory, "JBODFB") == 0 || strcasecmp(observatory, "8") == 0 || strcasecmp(observatory, "q") == 0 || strcasecmp(observatory, "JB_MKII") == 0 || strcasecmp(observatory, "JBMK2") == 0 || strcasecmp(observatory, "h") == 0 || strcasecmp(observatory, "JB42") == 0 || strcasecmp(observatory, "JB_42ft") == 0) {
    //      datafile->telescope_lat = M_PI*(53.233901)/180.0;
    //      datafile->telescope_long = -M_PI*(2+18/60.0 + 25.74/3600.0)/180.0;
    if(strcasecmp(observatory, "JODRELL") == 0 || strcasecmp(observatory, "JB") == 0 || strcasecmp(observatory, "JBO") == 0 || strcasecmp(observatory, "LOVELL") == 0 || strcasecmp(observatory, "JBDFB") == 0 || strcasecmp(observatory, "JBDFB") == 0 || strcasecmp(observatory, "JBODFB") == 0 || strcasecmp(observatory, "8") == 0 || strcasecmp(observatory, "q") == 0) {
      if(verbose.verbose) {
	for(indent = 0; indent < verbose.indent; indent++)      
	  printf(" ");
	fflush(stdout);
	printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Jodrell Bank (Lovell) data");
      }
      // Taken from tempo2 (with revised Jodrell coordinates)
      datafile->telescope_X = 3822252.643;
      datafile->telescope_Y = -153995.683;
      datafile->telescope_Z = 5086051.443;
    }else if(strcasecmp(observatory, "JB_MKII") == 0 || strcasecmp(observatory, "JBMK2") == 0 || strcasecmp(observatory, "h") == 0) {
      if(verbose.verbose) {
	for(indent = 0; indent < verbose.indent; indent++)      
	  printf(" ");
	fflush(stdout);
	printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Jodrell Bank (MKII) data");
      }
      // Taken from tempo2 (with revised Jodrell coordinates)
      datafile->telescope_X = 3822473.365;
      datafile->telescope_Y = -153692.318;
      datafile->telescope_Z = 5085851.303;
    }else if(strcasecmp(observatory, "JB42") == 0 || strcasecmp(observatory, "JB_42ft") == 0) {
      if(verbose.verbose) {
	for(indent = 0; indent < verbose.indent; indent++)      
	  printf(" ");
	fflush(stdout);
	printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Jodrell Bank (42 ft telescope) data");
      }
      // Taken from tempo2 (with revised Jodrell coordinates)
      datafile->telescope_X = 3822294.825;
      datafile->telescope_Y = -153862.275;
      datafile->telescope_Z = 5085987.071;
    }else {
      fflush(stdout);
      printerror(verbose.debug, "ERROR setITRFlocation_by_name: BUG!!!!!!");
      return 0;
    }
  }else if(strcasecmp(observatory, "WSRT") == 0 || strcasecmp(observatory, "i") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is the WSRT data");
    }
    //      datafile->telescope_lat = M_PI*(52.91691528690951429326)/180.0;
    //      datafile->telescope_long = M_PI*(6.60417111186232741460)/180.0;
    datafile->telescope_X = 3828445.659;
    datafile->telescope_Y = 445223.600000;
    datafile->telescope_Z = 5064921.5677;
  }else if(strcasecmp(observatory, "GBT") == 0 || strcasecmp(observatory, "1") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is GBT data");
    }
    datafile->telescope_X = 882589.65;
    datafile->telescope_Y = -4924872.32;
    datafile->telescope_Z = 3943729.348;
  }else if(strcasecmp(observatory, "ARECIBO") == 0 || strcasecmp(observatory, "AO") == 0 || strcasecmp(observatory, "3") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Arecibo data");
    }
    datafile->telescope_X = 2390490.0;
    datafile->telescope_Y = -5564764.0;
    datafile->telescope_Z = 1994727.0;
  }else if(strcasecmp(observatory, "EFFELSBERG") == 0 || strcasecmp(observatory, "eff") == 0 || strcasecmp(observatory, "g") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Effelsberg data");
    }
    datafile->telescope_X = 4033949.5;
    datafile->telescope_Y = 486989.4;
    datafile->telescope_Z = 4900430.8;
  }else if(strcasecmp(observatory, "GMRT") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is GMRT data");
    }
    datafile->telescope_X = 1656342.30;
    datafile->telescope_Y = 5797947.77;
    datafile->telescope_Z = 2073243.16;
  }else if(strcasecmp(observatory, "NANCAY") == 0 || strcasecmp(observatory, "ncy") == 0 || strcasecmp(observatory, "NUPPI") == 0 || strcasecmp(observatory, "ncyobs") == 0 || strcasecmp(observatory, "OP") == 0 || strcasecmp(observatory, "obspm") == 0 || strcasecmp(observatory, "f") == 0 || strcasecmp(observatory, "w") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Nancay data");
    }
    datafile->telescope_X = 4324165.81;
    datafile->telescope_Y = 165927.11;
    datafile->telescope_Z = 4670132.83;
  }else if(strcasecmp(observatory, "HARTRAO") == 0 || strcasecmp(observatory, "Hartebeesthoek") == 0 || strcasecmp(observatory, "hart") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is HARTRAO data");
    }
    datafile->telescope_X = 5085442.780;
    datafile->telescope_Y = 2668263.483;
    datafile->telescope_Z = -2768697.034;
  }else if(strcasecmp(observatory, "NANSHAN") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Nanshan data");
    }
    datafile->telescope_X = -228310.702;
    datafile->telescope_Y = 4631922.905;
    datafile->telescope_Z = 4367064.059;
  }else if(strcasecmp(observatory, "LOFAR") == 0 || strcasecmp(observatory, "t") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is central LOFAR data");
    }
    datafile->telescope_X = 3826577.462;
    datafile->telescope_Y = 461022.624;
    datafile->telescope_Z = 5064892.526;
  }else if(strcasecmp(observatory, "DE601LBA") == 0 || strcasecmp(observatory, "DE601LBH") == 0 || strcasecmp(observatory, "EFlfrlba") == 0 || strcasecmp(observatory, "EFlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Effelsberg international LOFAR station data (low band antenna's)");
    }
    datafile->telescope_X = 4034038.635;
    datafile->telescope_Y = 487026.223;
    datafile->telescope_Z = 4900280.057;
  }else if(strcasecmp(observatory, "DE601HBA") == 0 || strcasecmp(observatory, "DE601") == 0 || strcasecmp(observatory, "EFlfrhba") == 0 || strcasecmp(observatory, "EFlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Effelsberg international LOFAR station data (high band antenna's)");
    }
    datafile->telescope_X = 4034101.901;
    datafile->telescope_Y = 487012.401;
    datafile->telescope_Z = 4900230.210;
  }else if(strcasecmp(observatory, "DE602LBA") == 0 || strcasecmp(observatory, "DE602LBH") == 0 || strcasecmp(observatory, "UWlfrlba") == 0 || strcasecmp(observatory, "UWlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Garching/Unterweilenbach international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 4152561.068;
    datafile->telescope_Y = 828868.725;
    datafile->telescope_Z = 4754356.878;
  }else if(strcasecmp(observatory, "DE602HBA") == 0 || strcasecmp(observatory, "DE602") == 0 || strcasecmp(observatory, "UWlfrhba") == 0 || strcasecmp(observatory, "UWlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Garching/Unterweilenbach international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 4152568.416;
    datafile->telescope_Y = 828788.802;
    datafile->telescope_Z = 4754361.926;
  }else if(strcasecmp(observatory, "DE603LBA") == 0 || strcasecmp(observatory, "DE603LBH") == 0 || strcasecmp(observatory, "TBlfrlba") == 0 || strcasecmp(observatory, "TBlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Tautenburg international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 3940285.328;
    datafile->telescope_Y = 816802.001;
    datafile->telescope_Z = 4932392.757;
  }else if(strcasecmp(observatory, "DE603HBA") == 0 || strcasecmp(observatory, "DE603") == 0 || strcasecmp(observatory, "TBlfrhba") == 0 || strcasecmp(observatory, "TBlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Tautenburg international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 3940296.126;
    datafile->telescope_Y = 816722.532;
    datafile->telescope_Z = 4932394.152;
  }else if(strcasecmp(observatory, "DE604LBA") == 0 || strcasecmp(observatory, "DE604LBH") == 0 || strcasecmp(observatory, "POlfrlba") == 0 || strcasecmp(observatory, "POlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Potsdam international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 3796327.609;
    datafile->telescope_Y = 877591.315;
    datafile->telescope_Z = 5032757.252;
  }else if(strcasecmp(observatory, "DE604HBA") == 0 || strcasecmp(observatory, "DE604") == 0 || strcasecmp(observatory, "POlfrhba") == 0 || strcasecmp(observatory, "POlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Potsdam international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 3796380.254;
    datafile->telescope_Y = 877613.809;
    datafile->telescope_Z = 5032712.272;
  }else if(strcasecmp(observatory, "DE605LBA") == 0 || strcasecmp(observatory, "DE605LBH") == 0 || strcasecmp(observatory, "JUlfrlba") == 0 || strcasecmp(observatory, "JUlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Julich international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 4005681.742;
    datafile->telescope_Y = 450968.282;
    datafile->telescope_Z = 4926457.670;
  }else if(strcasecmp(observatory, "DE605HBA") == 0 || strcasecmp(observatory, "DE605") == 0 || strcasecmp(observatory, "JUlfrhba") == 0 || strcasecmp(observatory, "JUlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Julich international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 4005681.407;
    datafile->telescope_Y = 450968.304;
    datafile->telescope_Z = 4926457.940;
  }else if(strcasecmp(observatory, "FR606LBA") == 0 || strcasecmp(observatory, "FR606LBH") == 0 || strcasecmp(observatory, "FRlfrlba") == 0 || strcasecmp(observatory, "FRlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Nancey international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 4323980.155;
    datafile->telescope_Y = 165608.408;
    datafile->telescope_Z = 4670302.803;
  }else if(strcasecmp(observatory, "FR606HBA") == 0 || strcasecmp(observatory, "FR606") == 0 || strcasecmp(observatory, "FRlfrhba") == 0 || strcasecmp(observatory, "FRlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Nancey international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 4324017.054;
    datafile->telescope_Y = 165545.160;
    datafile->telescope_Z = 4670271.072;
  }else if(strcasecmp(observatory, "SE607LBA") == 0 || strcasecmp(observatory, "SE607LBH") == 0 || strcasecmp(observatory, "ONlfrlba") == 0 || strcasecmp(observatory, "ONlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Onsala international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 3370287.366;
    datafile->telescope_Y = 712053.586;
    datafile->telescope_Z = 5349991.228;
  }else if(strcasecmp(observatory, "SE607HBA") == 0 || strcasecmp(observatory, "SE607") == 0 || strcasecmp(observatory, "ONlfrhba") == 0 || strcasecmp(observatory, "ONlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Onsala international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 3370272.092;
    datafile->telescope_Y = 712125.596;
    datafile->telescope_Z = 5349990.934;
  }else if(strcasecmp(observatory, "UK608LBA") == 0 || strcasecmp(observatory, "UK608LBH") == 0 || strcasecmp(observatory, "UKlfrlba") == 0 || strcasecmp(observatory, "UKlfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Chilbolton international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 4008438.796;
    datafile->telescope_Y = -100310.064;
    datafile->telescope_Z = 4943735.554;
  }else if(strcasecmp(observatory, "UK608HBA") == 0 || strcasecmp(observatory, "UK608") == 0 || strcasecmp(observatory, "UKlfrhba") == 0 || strcasecmp(observatory, "UKlfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Chilbolton international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 4008462.280;
    datafile->telescope_Y = -100376.948;
    datafile->telescope_Z = 4943716.600;
  }else if(strcasecmp(observatory, "FI609LBA") == 0 || strcasecmp(observatory, "FI609LBH") == 0 || strcasecmp(observatory, "Filfrlba") == 0 || strcasecmp(observatory, "Filfrlbh") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Kilpisjarvi international LOFAR station data (low band antenna's) data");
    }
    datafile->telescope_X = 2136833.225;
    datafile->telescope_Y = 810088.740;
    datafile->telescope_Z = 5935285.279;
  }else if(strcasecmp(observatory, "FI609HBA") == 0 || strcasecmp(observatory, "FI609") == 0 || strcasecmp(observatory, "Filfrhba") == 0 || strcasecmp(observatory, "Filfr") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is Kilpisjarvi international LOFAR station data (high band antenna's) data");
    }
    datafile->telescope_X = 2136819.1940;
    datafile->telescope_Y = 810039.5757;
    datafile->telescope_Z = 5935299.0536;
  }else if(strcasecmp(observatory, "KAT7") == 0 || strcasecmp(observatory, "k7") == 0 || strcasecmp(observatory, "KAT-7") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is KAT-7 data");
    }
    datafile->telescope_X = 5109943.1050;
    datafile->telescope_Y = 2003650.7359;
    datafile->telescope_Z = -3239908.3195;
  }else if(strcasecmp(observatory, "MEERKAT") == 0) {  // Thanks to Maciej this telescope has been added
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is MeerKAT data");
    }
    // Coordinates set by Maciej
    //    datafile->telescope_X = 5109318.8410;
    //    datafile->telescope_Y = 2006836.3673;
    //    datafile->telescope_Z = -3238921.7749;
    // Changed to coordinates from OzStar on 28/3/2019
    datafile->telescope_X = 5109360.133;
    datafile->telescope_Y = 2006852.586;
    datafile->telescope_Z = -3238948.127;
  }else if(strcasecmp(observatory, "LWA") == 0 || strcasecmp(observatory, "LWA1") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is LWA data");
    }
    datafile->telescope_X = -1602196.60;
    datafile->telescope_Y = -5042313.47;
    datafile->telescope_Z = 3553971.51;
  }else if(strcasecmp(observatory, "FAST") == 0) {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "setITRFlocation_by_name: Guessing this is FAST data");
    }
    datafile->telescope_X = -1668557.0;
    datafile->telescope_Y = 5506838.0;
    datafile->telescope_Z = 2744934.0;
  }else {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fflush(stdout);
      printwarning(verbose.debug, "WARNING setITRFlocation_by_name: Cannot guess location of telescope '%s'", observatory);
    }
    return 0;
  }
  return 1;
}

/*
  Finds out if there are weights read in that are zero
  (zeroweightfound != 0), if there are weights that are different from
  each other which means 0 + one other value is allowed
  (differentweights != 0), if there are negative weights found
  (negativeweights != 0). The non-zero weight value is returned as
  weightvalue. If weights are undefined, the return values are as if
  all weights were set to 1.
*/
void determineWeightsStat(datafile_definition *datafile)
{
  int firstnonzeroweight = 1;
  long n;
  datafile->weight_stats_zeroweightfound = 0;
  datafile->weight_stats_differentweights = 0;
  datafile->weight_stats_negativeweights = 0;
  if(datafile->weights == NULL) {
    datafile->weight_stats_weightvalue = 1.0;
    return;
  }
  for(n = 0; n < datafile->NrSubints*datafile->NrFreqChan; n++) {
    //    printf("XXXXXX %ld: %f\n", n, datafile->weights[n]);
    if(datafile->weights[n] < 0) {
      datafile->weight_stats_negativeweights = 1;
    }
    if(datafile->weights[n] == 0.0) {
      datafile->weight_stats_zeroweightfound = 1;
    }else {
      if(firstnonzeroweight) {
	datafile->weight_stats_weightvalue = datafile->weights[n];
	firstnonzeroweight = 0;
      }else if(datafile->weights[n] != datafile->weight_stats_weightvalue) {
	datafile->weight_stats_differentweights = 1;
      }
    }
  }
  datafile->weight_stats_set = 1;
}

/* Read the header. If readnoscales is set, the scales/offsets/weights
  of the subints are not read in (if supported/used by file
  format). Useful when only interested in header information.  Verbose
  level determines nr of spaces before output. If debug is set, more
  information is printed. If nowarnings == 0, warnings are shown. If
  1, no warnings are generated for incomplete headers (the checks done
  AFTER the file-format specific reading has been done), which can be
  useful for applications which do not rely on the header
  information. If 2, all warnings are ignored while reading in the
  header. Returns 0 on error, 1 on success. */
int readHeaderPSRData(datafile_definition *datafile, int readnoscales, int nowarnings, verbose_definition verbose)
{
  int i, nowarnings2, ret = 0;
  verbose_definition verbose2;

  if(nowarnings == 2) {
    nowarnings2 = 1;
  }else {
    nowarnings2 = 0;
  }

  if(verbose.debug) {
    printf("Entering readHeaderPSRData()\n");
  }

  copyVerboseState(verbose, &verbose2);
  verbose2.indent = verbose.indent + 2;

  if(datafile->format == PSRSALSA_BINARY_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading PSRSALSA binary header\n");
    }
    ret = readPSRSALSAHeader(datafile, 0, verbose);
  }else if(datafile->format == PUMA_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading WSRT header\n");
    }
    ret = readWSRTHeader(datafile, verbose);
//START REGION DEVELOP
  }else if(datafile->format == PRESTO_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading PRESTO header\n");
    }
    ret = readPRESTOHeader(datafile, verbose);
  }else if(datafile->format == PARKESFB_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading Parkes header\n");
    }
    ret = readParkesHeader(datafile, verbose);
//START REGION RELEASE
  }else if(datafile->format == PSRCHIVE_ASCII_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading PSRCHIVE ascii header\n");
    }
    ret = readPSRCHIVE_ASCIIHeader(datafile, verbose);
  }else if(datafile->format == EPN_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      //      printf("Reading EPN header (will be done again automatically when reading the whole file)\n");
      printf("Reading EPN header\n");
    }
    ret = readEPNHeader(datafile, 1, verbose);
    if(ret == 1) {
      float scale, offset;
      ret = 0;
      ret = readEPNsubHeader(datafile, &scale, &offset, verbose);
    }
    /*    if(verbose.verbose) 
	  printf("readHeaderPSRData (%s): Reading EPN file will be done automatically when reading the whole file.\n", datafile->filename);
	  return 0;*/
//START REGION DEVELOP
  }else if(datafile->format == AO_ASCII_1_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading AO ascii format 1\n");
    }
    ret = readAOAsciiHeader(datafile, 1, verbose);
  }else if(datafile->format == AO_ASCII_2_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading AO ascii format 2\n");
    }
    ret = readAOAsciiHeader(datafile, 2, verbose);
//START REGION RELEASE
  }else if(datafile->format == FITS_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading PSRFITS header\n");
    }
    ret = readPSRFITSHeader(datafile, readnoscales, nowarnings2, verbose);
    if(verbose.debug) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading PSRFITS header done\n");
    }
  }else if(datafile->format == SIGPROC_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading sigproc header\n");
    }
    ret = readSigprocHeader(datafile, verbose);
 }else if(datafile->format == PPOL_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading ppol file header\n");
    }
    ret = readPPOLHeader(datafile, 1, verbose);
  }else if(datafile->format == PPOL_SHORT_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading short ppol file header\n");
    }
    ret = readPPOLHeader(datafile, 0, verbose);
 }else if(datafile->format == SIGPROC_ASCII_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading sigproc ascii header\n");
    }
    ret = readSigprocASCIIHeader(datafile, verbose);
//START REGION DEVELOP
  }else if(datafile->format == GMRT_ASCII_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      printf("Reading GMRT folded ascii format\n");
    }
    ret = readGMRTasciiHeader(datafile, verbose);
//START REGION RELEASE
  }else if(datafile->format == MEMORY_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHeaderPSRData (%s): The file is already read in, cannot read header.", datafile->filename);
    return 0;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Data type (%s) not implemented in readHeaderPSRData", datafile->filename, returnFileFormat_str(datafile->format));
    return 0;
  }
  if(ret == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Error reading header.", datafile->filename);
    return 0;
  }
  if(datafile->format != FITS_format) {
    datafile->datastart = ftell(datafile->fptr);
  }else {
    datafile->datastart = 0;
  }

  if(convert_if_uniform_frequency_spacing(datafile, nowarnings2, verbose) == 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Attempt to convert dataset in an equally spaced frequency channel format failed.", datafile->filename);
    return 0;    
  }


  /* If longitude and latitude of the site appear to be unset, try to guess it */
  //  if(fabs(datafile->telescope_long) < 1e-6 && fabs(datafile->telescope_lat) < 1e-6) {
  if(fabs(datafile->telescope_X) < 1e-6 && fabs(datafile->telescope_Y) < 1e-6 && fabs(datafile->telescope_Z) < 1e-6) {
    fflush(stdout);
    if(strlen(datafile->observatory) == 0) {
      if(nowarnings == 0) {
	printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): Telescope location not set", datafile->filename);
      }
    }else {
      if(nowarnings == 0) {
	printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): Telescope location not set, try to determine location from name '%s'", datafile->filename, datafile->observatory);
      }
      // Always output what telescope is selected, as this should result in a warning as telescope location ought to be defined in the file.
      if(setITRFlocation_by_name(datafile, datafile->observatory, verbose2) == 0) {
	fflush(stdout);
	if(nowarnings == 0) {
	  printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): Telescope location not recognized by name", datafile->filename);
	}
      }
    }
  }

  double period;
  if(datafile->isFolded == 1) {
    if(get_period(*datafile, 0, &period, verbose) == 2) {
      printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Cannot obtain period", datafile->filename);
      return 0;
    }
  }else {
    period = -1;
  }

  if(period < 0.001 && datafile->isFolded != 0) {
    fflush(stdout);
    if(nowarnings == 0) {
      printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The period does not appear to be set in the header. Consider using the -header option.", datafile->filename);
      if(verbose.debug) {
	printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): isFolded=%d\n", datafile->filename, datafile->isFolded);
      }
    }
  }
  if(datafile->tsampMode != TSAMPMODE_LONGITUDELIST) {  // If longitudes are defined, no sampling time is required (PA-swing files for instance)
    if(datafile->tsampMode == TSAMPMODE_UNKNOWN) {
      fflush(stdout);
      if(nowarnings == 0) {
	printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The sampling time is unknown.", datafile->filename);      
      }
    }else {
      if(get_tsamp(*datafile, 0, verbose) < 0.0000001 || get_tsamp(*datafile, 0, verbose) >= 100) {
	fflush(stdout);
	if(nowarnings == 0) {
	  printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The sampling time does not appear to be set correctly in the header. Consider using the -header option.", datafile->filename);
	}
      }
    }
  }
  if(datafile->isFolded == 0 && datafile->gentype == GENTYPE_UNDEFINED)
    datafile->gentype = GENTYPE_SEARCHMODE;

  // Check if full period is stored, if regularly sampled data.
  if(datafile->tsampMode != TSAMPMODE_LONGITUDELIST && datafile->tsampMode != TSAMPMODE_UNKNOWN) {
    double diff_nbin;
    double period;
    if(datafile->isFolded == 1) {
      if(get_period(*datafile, 0, &period, verbose) == 2) {
	printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Cannot obtain period", datafile->filename);
	return 0;
      }
    }else {
      period = -1;
    }
    if(get_tsamp(*datafile, 0, verbose) != 0)
      diff_nbin = period/get_tsamp(*datafile, 0, verbose) - datafile->NrBins;
    else
      diff_nbin = 12345;
    
    if(datafile->isFolded && datafile->gentype != GENTYPE_2DFS && datafile->gentype != GENTYPE_RECEIVERMODEL  && datafile->gentype != GENTYPE_RECEIVERMODEL2) {
      if(verbose.debug) {
	fflush(stdout);
	fprintf(stderr, "readHeaderPSRData (%s): Check if full period is stored - Tsamp=%lf, NBIN=%ld: Period=%lf suggests nr of bins is off by %lf.\n", datafile->filename, get_tsamp(*datafile, 0, verbose), datafile->NrBins, period, diff_nbin);
      }
      
      if(diff_nbin < -0.5 || diff_nbin > 0.5) {
	fflush(stdout);
	if(nowarnings == 0) {
	  printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The sampling and period suggest that not the whole rotational phase range is stored. If not correct, consider using the -header option.", datafile->filename);
	}
      }
    }
  }

  if(datafile->gentype == GENTYPE_PULSESTACK) {
    double tobs_expected;
    double period;
    if(get_period(*datafile, 0, &period, verbose) == 2) {
      printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Cannot obtain period", datafile->filename);
      return 0;
    }
    tobs_expected = datafile->NrSubints * period;
    if(fabs(get_tsub(*datafile, 0, verbose)) < 0.0001) {  // If tobs not set, set it.
      datafile->tsubMode = TSUBMODE_FIXEDTSUB;
      if(datafile->tsub_list != NULL)
	free(datafile->tsub_list);
      datafile->tsub_list = (double *)malloc(sizeof(double));
      if(datafile->tsub_list == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR readHeaderPSRData (%s): Memory allocation error.", datafile->filename);
	return 0;
      }
      datafile->tsub_list[0] = period;
    }else if(get_tobs(*datafile, verbose)/tobs_expected > 1.02 || get_tobs(*datafile, verbose)/tobs_expected < 0.98) { // See if close enough
      fflush(stdout);
      if(nowarnings == 0) {
	printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The period and nr of pulses appear to be incompatable with tobs given that these are thought to be single pulses. Gentype is set to undefined. Consider using the -header option to fix problem.", datafile->filename);
      }
      datafile->gentype = GENTYPE_UNDEFINED;
    }
  }

  // If reference frequency is not set
  if(datafile->freq_ref < -1.1) {
    if(datafile->isDeDisp || datafile->isDeFarad) {
      if(nowarnings == 0) {
	printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): The reference frequency used in dedispersion/de-Faraday rotation is unknown.", datafile->filename);
      }
    }else {  // At this point the reference frequency is unimportant, so set to infinity.
      datafile->freq_ref = -1;
    }
  }
  if(datafile->freq_ref < -0.9 && datafile->freq_ref >= -1.1) {  // Corresponds to infinite frequency
    datafile->freq_ref = 1e10; // Now the prefered encoding of an infinite frequency
  }

  if(readnoscales == 0) {
    determineWeightsStat(datafile);
    if(datafile->weight_stats_negativeweights) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): Found negative weights in input file, so probably the file is corrupted. You might want to use the -absweights option.", datafile->filename);
    }
    if(datafile->weight_stats_zeroweightfound && verbose.verbose) {
      printf("FITS file contains zero-weighted data. By default this is applied, unless the -noweights option is used.\n");
    }
    if(datafile->weight_stats_differentweights == 0 && (verbose.debug || datafile->weight_stats_weightvalue != 1.0)) {
      printf("FITS file contains uniform weights with a value %f", datafile->weight_stats_weightvalue);
      if(datafile->weight_stats_zeroweightfound)
	printf(" (appart from zero weights).\n");
      else
	printf("\n");
    }
    if(datafile->weight_stats_differentweights) {
      fflush(stdout);
      printwarning(verbose.debug, "WARNING readHeaderPSRData (%s): FITS file contains data with different weights. Make sure it is used as desired.", datafile->filename);
    }
  }else {
    datafile->weight_stats_set = 1;
    datafile->weight_stats_zeroweightfound = 0;
    datafile->weight_stats_differentweights = 0;
    datafile->weight_stats_negativeweights = 0;
    datafile->weight_stats_weightvalue = 1.0;
  }

  if(readHistoryPSRData(datafile, verbose2) == 0) {
    printwarning(verbose.debug, "WARNING: Reading history failed.");
  }
  //  showHistory(*datafile, verbose2);




  if(verbose.verbose) {
    printHeaderPSRData(*datafile, 0, verbose2);
  }

  if(verbose.debug) {
    printf("Exiting readHeaderPSRData()\n");
  }
  return 1;
}


//START REGION DEVELOP
//START REGION RELEASE

/* Write the header, including the present command line (if argc >
   0). If cmdOnly is set, only the command is writen out. This makes
   the history table without timestamps and therefore results are
   exactly reproducable.  Verbose level determines nr of spaces before
   output. 
   Returns 1 if successful or 0 on error. */
int writeHeaderPSRData(datafile_definition *datafile, int argc, char **argv, int cmdOnly, verbose_definition verbose)
{
  int i, ret;  // , debug
  //  debug = 0;
  if(datafile->format == PSRSALSA_BINARY_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      if(verbose.verbose) printf("Write PSRSALSA header.\n");
    }
    ret = writePSRSALSAHeader(datafile, verbose);
    if(verbose.debug) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      if(verbose.verbose) printf("Writing PSRSALSA header done.\n");
    }
  }else if(datafile->format == PUMA_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      if(verbose.verbose) printf("Write PuMa header.\n");
    }
    ret = writeWSRTHeader(*datafile, verbose);
    if(verbose.debug) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      if(verbose.verbose) printf("Writing PuMa header done.\n");
    }
//START REGION DEVELOP
  }else if(datafile->format == AO_ASCII_1_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    if(verbose.verbose) printf("Write AO ascii type 1 header.\n");
    ret = writeAOAsciiHeader(*datafile, 1, verbose);
  }else if(datafile->format == AO_ASCII_2_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    if(verbose.verbose) printf("Write AO ascii type 2 header.\n");
    ret = writeAOAsciiHeader(*datafile, 2, verbose);
  }else if(datafile->format == PRESTO_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    if(verbose.verbose) printf("Write Presto header.\n");
    ret = writePRESTOHeader(*datafile, verbose);
//START REGION RELEASE
  }else if(datafile->format == PSRCHIVE_ASCII_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    if(verbose.verbose) printf("Write PSRCHIVE_ASCII header.\n");
    ret = writePSRCHIVE_ASCIIHeader(*datafile, verbose);
  }else if(datafile->format == FITS_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    if(verbose.verbose) printf("Write PSRFITS header.\n");
    ret = writePSRFITSHeader(datafile, verbose);
    if(verbose.debug) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
      if(verbose.verbose) printf("Writing PSRFITS header done.\n");
    }
  }else if(datafile->format == EPN_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    if(verbose.verbose) {
      printwarning(verbose.debug, "WARNING writeHeaderPSRData: Write EPN header not necessary, will when writing data.");
    }
    ret = 1;
  }else if(datafile->format == SIGPROC_ASCII_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    if(verbose.verbose) printf("Write Sigproc ascii header.\n");
    ret = writeSigprocASCIIHeader(*datafile, verbose);
  }else if(datafile->format == PPOL_format || datafile->format == PPOL_SHORT_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)      
	printf(" ");
    }
    if(verbose.verbose) printf("Write ppol file header.\n");
    ret = writePPOLHeader(*datafile, argc, argv, verbose);
//START REGION DEVELOP
  }else if(datafile->format == GMRT_ASCII_format) {
    if(verbose.verbose) {
      for(i = 0; i < verbose.indent; i++)
	printf(" ");
    }
    if(verbose.verbose) printf("Write GMRT folded ascii header.\n");
    printerror(verbose.debug, "ERROR writeHeaderPSRData: Writing of this type of header is not implemented");
    ret = 0;
//START REGION RELEASE
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writeHeaderPSRData: Writing of this type of header is not implemented");
    ret = 0;
  }
  if(ret != 0) {
    //    showHistory(*datafile, verbose);
    if(datafile->format == PSRSALSA_BINARY_format || datafile->format == FITS_format || datafile->format == PUMA_format) {
      writeHistoryPSRData(datafile, argc, argv, cmdOnly, verbose);
    }
  }
  if(ret == 0) {
    printerror(verbose.debug, "ERROR writeHeaderPSRData: Writing header failed.");
  }
  return ret;
}

//START REGION DEVELOP

/* Obtains a pointer pulse_ptr which points to a single pulse (or
   subint) starting at binnr etc. The behaviour is unpredictable if
   you read in beyond a subint/polarization/frequency channel. An
   error will be generated if the data was not yet read into memory.

  Returns 1 if successful, 0 on error. */
int get_pointer_PulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, float **pulse_ptr, verbose_definition verbose)
{
  if(datafile->format != MEMORY_format) {
    printerror(verbose.debug, "ERROR get_pointer_PulsePSRData: Data does not appear to exist in memory yet.");
    return 0;
  }

  *pulse_ptr = &(datafile->data[datafile->NrBins*(polarization+datafile->NrPols*(freq+pulsenr*datafile->NrFreqChan))+binnr]);
  return 1;
}

//START REGION RELEASE

/* Read a single pulse (or subint) starting at binnr and with length
   nrSamples. The behaviour of this function is unpredictable if you
   read in beyond a subint/polarization/frequency channel. The
   reading should not nessesarily be from start to end. If the data is
   psrfits, and psrfits_set_use_weighted_freq function is used to set
   the weighted frequency flag, the centre frequency is set to the
   wieghted frequency of the first pulse read in. Returns 1 if
   successful. */
int readPulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  if(datafile->format == PSRSALSA_BINARY_format)
    return readPulsePSRSALSAData(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  else if(datafile->format == PUMA_format)
    return readPulseWSRTData(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse);
  else if(datafile->format == FITS_format) 
    return readFITSpulse(datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  else if(datafile->format == PSRCHIVE_ASCII_format) 
    return readPSRCHIVE_ASCIIfilepulse(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  else if(datafile->format == EPN_format) 
    return readPulseEPNData(datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  else if(datafile->format == SIGPROC_format) 
    return readPulseSigprocData(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  else if(datafile->format == MEMORY_format) {
    memcpy(pulse, &datafile->data[datafile->NrBins*(polarization+datafile->NrPols*(freq+pulsenr*datafile->NrFreqChan))+binnr], sizeof(float)*nrSamples);
    return 1;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPulsePSRData: Reading of this format is not implemented. Maybe converting the data in a different format will solve this issue.");
  }
  return 0;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Write a single pulse/subint starting at binnr and with length
   nrSamples. The behaviour of this function is unpredictable if you
   read in beyond a subint/polarization/frequency channel. The
   writing should not nessesarily be from start to end. Return 1 if
   successful. */
int writePulsePSRData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  if(pulsenr < 0 || binnr < 0 || freq < 0) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePulsePSRData: Parameters outside boundaries.");
    return 0;
  }
  if(datafile->format == MEMORY_format || datafile->dumpOnClose) {
    // If not allocated memory previously, allocate it not allowing the buffering. Instead of writing data to file, write it to this buffer
    if(datafile->dumpOnClose && datafile->data == NULL) {
      long datasize = datafile->NrSubints*datafile->NrBins*datafile->NrPols*datafile->NrFreqChan*sizeof(float);
      datafile->data = (float *)malloc(datasize);
      //	fprintf(stderr, "\nXXXXX Allocating %ld bytes for data\n", datafile->NrSubints*datafile->NrBins*datafile->NrPols*datafile->NrFreqChan*sizeof(float));
      if(datafile->data == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR writePulsePSRData: Cannot allocate memory (data=%ld bytes=%.3fGB).", datasize, datasize/1073741824.0);
	datafile->dumpOnClose = 0;  // To avoid more writing when closing file
	verbose_definition verbose2;
	copyVerboseState(verbose, &verbose2);
	if(verbose.debug == 0)
	  verbose2.verbose = 0;
	verbose2.nocounters = 1;
	closePSRData(datafile, 0, verbose2);
	return 0;
      }else if(verbose.debug) {
	printf("DEBUG: Allocated %ld bytes of memory for memory buffering.\n", datasize);
      }
    }
    //    if(verbose.debug) {
    //      printf("DEBUG: writing %ld bytes to memory for buffering at position %ld.\n", sizeof(float)*nrSamples, datafile->NrBins*(polarization+datafile->NrPols*(freq+pulsenr*datafile->NrFreqChan))+binnr);
    //    }
    memcpy(&datafile->data[datafile->NrBins*(polarization+datafile->NrPols*(freq+pulsenr*datafile->NrFreqChan))+binnr], pulse, sizeof(float)*nrSamples);
  }else if(datafile->format == PSRSALSA_BINARY_format) {
    return writePulsePSRSALSAData(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  }else if(datafile->format == PUMA_format) {
    return writePulseWSRTData(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse);
  }else if(datafile->format == FITS_format) {
    return writeFITSpulse(*datafile, pulsenr, polarization, freq, binnr, nrSamples, pulse, verbose);
  }else if(datafile->format == PSRCHIVE_ASCII_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePulsePSRData: Writing out individual subintegrations is not implemented for ASCII formats.");
    return 0;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePulsePSRData: Writing in this format (%s) is not implemented", returnFileFormat_str(datafile->format));
    return 0;
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Read the whole observation at once in memory. The data loops over
   pulse number, frequency channel, polarization channel and bin nr.
   If the data is psrfits, and psrfits_set_use_weighted_freq function
   is used to set the weighted frequency flag, the centre frequency is
   set to the wieghted frequency of the first pulse read in. If
   nocounters is set, no counters are shown. Return 0 if there is an
   error. */
int readPSRData(datafile_definition *datafile, float *data, verbose_definition verbose)
{
  if(datafile->format == PSRSALSA_BINARY_format)
    return readPSRSALSAfile(*datafile, data, verbose);
  else if(datafile->format == PUMA_format)
    return readPuMafile(*datafile, data, verbose);
  else if(datafile->format == PSRCHIVE_ASCII_format)
    return readPSRCHIVE_ASCIIfile(*datafile, data, verbose);
//START REGION DEVELOP
  else if(datafile->format == AO_ASCII_1_format)
    return readAOAsciifile(*datafile, data, 1, verbose);
  else if(datafile->format == AO_ASCII_2_format)
    return readAOAsciifile(*datafile, data, 2, verbose);
//START REGION RELEASE
  else if(datafile->format == EPN_format)
    return readEPNfile(datafile, data, verbose, -1);   // No debug information is printed
  else if(datafile->format == FITS_format) {
    int ret;
    ret = readFITSfile(datafile, data, verbose);
    return ret;
  }else if(datafile->format == PPOL_format)
    return readPPOLfile(datafile, data, 1, 0, verbose);
  else if(datafile->format == PPOL_SHORT_format)
    return readPPOLfile(datafile, data, 0, 0, verbose);
  else if(datafile->format == SIGPROC_ASCII_format)
    return readSigprocASCIIfile(*datafile, data, verbose);
  else if(datafile->format == SIGPROC_format)
    return readSigprocfile(*datafile, data, verbose);
//START REGION DEVELOP
  else if(datafile->format == GMRT_ASCII_format)
    return readGMRTasciifile(*datafile, data, verbose);
//START REGION RELEASE
  else 
    fflush(stdout);
    printerror(verbose.debug, "ERROR readPSRData: Reading whole dataset is not supported for this type of data.");
  return 0;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Write the whole dataset at once from memory. The data loops over
   pulse number, frequency channel, polarization channel and bin
   nr. Set pa data pointers to something random if not used. Returns 1
   if successful. If nocounters is set, no counters are shown. */
int writePSRData(datafile_definition *datafile, float *data, verbose_definition verbose)
{
  if(verbose.verbose) printf("Writing %ld x %ld x %ld x %ld samples\n", datafile->NrSubints, datafile->NrFreqChan, datafile->NrBins, datafile->NrPols);
  datafile->dumpOnClose = 0;  // To avoid more writing when closing file as data is already written. Data is written out in this function.
  if(datafile->format == PSRSALSA_BINARY_format) {
    return writePSRSALSAfile(*datafile, data, verbose);
  }else if(datafile->format == PUMA_format) {
    return writePuMafile(*datafile, data, verbose);
  }else if(datafile->format == EPN_format) {
    return writeEPNfile(*datafile, data, verbose);
//START REGION DEVELOP
  }else if(datafile->format == AO_ASCII_1_format) {
    return writeAOAsciifile(*datafile, data, 1, verbose);
  }else if(datafile->format == AO_ASCII_2_format) {
    return writeAOAsciifile(*datafile, data, 2, verbose);
//START REGION RELEASE
  }else if(datafile->format == PSRCHIVE_ASCII_format) {
    return writePSRCHIVE_ASCIIfile(*datafile, data, verbose);
  }else if(datafile->format == FITS_format)  {
    return writeFITSfile(*datafile, data, verbose);
  }else if(datafile->format == PPOL_format) {
    return writePPOLfile(*datafile, data, 1, 0, 0, 0, verbose);
  }else if(datafile->format == PPOL_SHORT_format) {
    return writePPOLfile(*datafile, data, 0, 0, 0, 0, verbose);
  }else if(datafile->format == SIGPROC_ASCII_format) {
    return writeSigprocASCIIfile(*datafile, data, verbose);
//START REGION DEVELOP
  }else if(datafile->format == GMRT_ASCII_format) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRData: Writing whole dataset is not supported for this type of data (%s).", returnFileFormat_str(datafile->format));
    return 0;
//START REGION RELEASE
  }else  {
    fflush(stdout);
    printerror(verbose.debug, "ERROR writePSRData: Writing whole dataset is not supported for this type of data (%s).", returnFileFormat_str(datafile->format));
  }
  return 0;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Read file and produce a pulse profile. Set zapMask to NULL is you
   don't want to zap pulses. Set polchan to select the polarization
   channel. Returns 1 if successful. */
int read_profilePSRData(datafile_definition datafile, float *profileI, int *zapMask, int polchan, verbose_definition verbose)
{
  return read_partprofilePSRData(datafile, profileI, zapMask, polchan, 0, datafile.NrSubints, verbose);
}

//START REGION DEVELOP
//START REGION RELEASE

/* Identical to read_profilePSRData, but now you can set nread and nskip. Returns 1 if successful. */
int read_partprofilePSRData(datafile_definition datafile, float *profileI, int *zapMask, int polchan, long nskip, long nread, verbose_definition verbose)
{
  int zap;
  long i, j, k; //, l;
  float *data;
  if(verbose.verbose) printf("Generating average pulse profile (polarization channel %d)\n", polchan);

  //  if(datafile.freqMode != FREQMODE_UNIFORM) {
  //    fflush(stdout);
  //    printerror(verbose.debug, "ERROR read_partprofilePSRData: Frequency channels need to be uniformely separated.");
  //    return 0;
  //  }


  data = (float *)malloc(datafile.NrBins*sizeof(float));
  if(data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR read_partprofilePSRData: Cannot allocate memory.");
    return 0;
  }
  for(j = 0; j < datafile.NrBins; j++) {
    profileI[j] = 0;
  }
  for(k = 0; k < datafile.NrFreqChan; k++) {
    for(i = nskip; i < nskip+nread; i++) {
      if(readPulsePSRData(&datafile, i, polchan, k, 0, datafile.NrBins, data, verbose) == 0)
	return 0;
      zap = 0;
      if(zapMask != NULL) {
	if(zapMask[i] != 0)
	  zap = 1;
      }
      if(zap == 0) {
	for(j = 0; j < datafile.NrBins; j++) {
	  profileI[j] += data[j];
	}
      }
      /*
      l = (i-nskip)/10;
      if(verbose && (10*l == i) && (l > 0)) {
	printf("read_partprofilePSRData: processed %.1f%%     \r",(100.0*(i-nskip))/(float)(nread));
	fflush(stdout);
	}*/
    }
  }
  free(data);
  /*  if(verbose.verbose) printf("  Done                                      \n"); */
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Read data file and produce a list of rms'es and averages. Set
   zapMask to NULL is you don't want to zap pulses. The defined
   regions are NOT included, unless invert is set. If onlyI is set
   only the total intensity is generated. Set polchan to select the
   polarization channel. Set freqchan to -1 if you want to get average
   of all channels instead of a specific channel. Set rms to NULL if
   you're not interested in the rms'es, likewise for avrg. Return 0 if
   there is an error. This function is very similar to rmsPSRData(). */
int read_rmsPSRData(datafile_definition datafile, float *rms, float *avrg, int *zapMask, pulselongitude_regions_definition *regions, int invert, int polchan, int freqchan, verbose_definition verbose)
{
  long i, j, k, freq0, freq1;
  float *data;
  double *rms_double, *avrg_double;
  int nrOffpulseBins, zap, zap2;

  /*  if(verbose.verbose) printf("Generating list of rms'es\n"); */
  data = (float *)malloc(datafile.NrBins*sizeof(float));
  /* Must store rms'es and avrg'es in double precission, or else the
     baseline subtraction of the rms values will not be precise enough
     if there is a large DC components. */
  rms_double = (double *)malloc(datafile.NrSubints*sizeof(double));
  avrg_double = (double *)malloc(datafile.NrSubints*sizeof(double));
  if(data == NULL || rms_double == NULL || avrg_double == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR read_rmsPSRData: Cannot allocate memory.");
    return 0;
  }
  for(i = 0; i < datafile.NrSubints; i++) {
    rms_double[i] = 0;
    avrg_double[i] = 0;
  }
  nrOffpulseBins = datafile.NrBins;
  if(regions != NULL) {
    for(j = 0; j < datafile.NrBins; j++) {
      if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0))
	nrOffpulseBins--;
    }
  }

  if(nrOffpulseBins == 0) {
    printerror(verbose.debug, "ERROR read_rmsPSRData: An off-pulse rms was requested, but everything is defined as onpulse.");
    return 0;
  }

  if(freqchan < 0) {
    freq0 = 0;
    freq1 = datafile.NrFreqChan;
    freqchan = 0;
  }else {
    freq0 = freqchan;
    freq1 = freqchan+1;
    freqchan = 1;
  }
  for(k = freq0; k < freq1; k++) {
    for(i = 0; i < datafile.NrSubints; i++) {
      if(readPulsePSRData(&datafile, i, polchan, k, 0, datafile.NrBins, data, verbose) == 0)
	return 0;
      zap = 0;
      if(zapMask != NULL) {
	if(zapMask[i] != 0)
	  zap = 1;
      }
      for(j = 0; j < datafile.NrBins; j++) {
	zap2 = zap;
	if(regions != NULL) {
	  if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0))
	    zap2 = 1;
	}
	if(zap2 == 0) {
	  /*	    printf("I=%f\n", data[j]);  */
	  rms_double[i] += data[j]*data[j];
	  avrg_double[i] += data[j];
	}
      }
       /*
      l = i/10;
     if((10*l == i) && (l > 0)) {
	printf("read_rmsPSRData: processed %.1f%%     \r",(100.0*i)/(float)(datafile.NrSubints));
	fflush(stdout);
	}*/
    }
  }
  if(freqchan == 0) {
    freqchan = datafile.NrFreqChan;
  }else {
    freqchan = 1;
  }
  for(i = 0; i < datafile.NrSubints; i++) {
    /*    printf("rms=%f av=%f\n", rms_double[i], avrg_double[i]); */
    avrg_double[i] /= (double)(nrOffpulseBins*freqchan);
    /* sum (I-a)^2 = sum (I*I -2aI + a*a) = sum (I*I) -2a*N*a + a*a*N = sum(I*I) - N*a*a */
    rms_double[i] -= nrOffpulseBins*freqchan*avrg_double[i]*avrg_double[i];
    /*    printf("rms2=%f %d\n", rms_double[i], nrOffpulseBins); */
    rms_double[i] /= (double)(nrOffpulseBins*freqchan);
    if(rms != NULL) 
      rms[i] = sqrt(rms_double[i]);
    if(avrg != NULL)
      avrg[i] = avrg_double[i];
  }
  free(data);
  free(rms_double);
  free(avrg_double);
  /*  printf("  Done                              \n"); */
  return 1;
}

//START REGION DEVELOP

/* Read data file and obtain an rms and average for a given subint,
   polchan and freqchan. The defined regions are NOT included, unless
   invert is set. Set rms to NULL if you're not interested in the rms,
   likewise for avrg. Return 0 if there is an error. This function is
   very similar to read_rmsPSRData().*/
int rmsPSRData(datafile_definition datafile, float *rms, float *avrg, pulselongitude_regions_definition *regions, int invert, long subint, long polchan, long freqchan, verbose_definition verbose)
{
  long j;
  float *data;
  double rms_double, avrg_double;
  int nrOffpulseBins, zap;

  /*  if(verbose.verbose) printf("Generating list of rms'es\n"); */
  data = (float *)malloc(datafile.NrBins*sizeof(float));
  if(data == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR rmsPSRData: Cannot allocate memory.");
    return 0;
  }
  rms_double = 0;
  avrg_double = 0;
  nrOffpulseBins = datafile.NrBins;
  if(regions != NULL) {
    for(j = 0; j < datafile.NrBins; j++) {
      if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0))
	nrOffpulseBins--;
    }
  }

  if(readPulsePSRData(&datafile, subint, polchan, freqchan, 0, datafile.NrBins, data, verbose) == 0)
    return 0;
  for(j = 0; j < datafile.NrBins; j++) {
    zap = 0;
    if(regions != NULL) {
      if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0))
	zap = 1;
    }
    if(zap == 0) {
      /*	    printf("I=%f\n", data[j]);  */
      rms_double += data[j]*data[j];
      avrg_double += data[j];
    }
  }

  avrg_double /= (double)(nrOffpulseBins);
  /* sum (I-a)^2 = sum (I*I -2aI + a*a) = sum (I*I) -2a*N*a + a*a*N = sum(I*I) - N*a*a */
  rms_double -= nrOffpulseBins*avrg_double*avrg_double;
  /*    printf("rms2=%f %d\n", rms_double[i], nrOffpulseBins); */
  rms_double /= (double)(nrOffpulseBins);
  if(rms != NULL) 
    *rms = sqrt(rms_double);
  if(avrg != NULL)
    *avrg = avrg_double;
  free(data);
  /*  printf("  Done                              \n"); */
  return 1;
}

/* Read data file and produce a list of correlations and averages. The
   file pointer should be set to the start of the data. Set zapMask to
   NULL is you don't want to zap pulses. The defined regions are NOT
   included, unless invert is set. Set polchan to select the
   polarization channel.  */
int read_correlPSRData(datafile_definition datafile, float *correl, float *avrg, int *zapMask, pulselongitude_regions_definition *regions, int invert, int polchan, verbose_definition verbose)
{
  long i, j, k;
  float Ilast, *data;
  double *correl_double, *avrg_double;
  int nrOffpulseBins, zap, zap2;

  /*
// Not sure why this check would be relevant. Only put in so I didn't had to think of consequences?
  if(datafile.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR read_correlPSRData: Frequency channels need to be uniformely separated.");
    return 0;
  }
  */

  /*  if(verbose.verbose) printf("Generating list of correlations\n"); */

  data = (float *)malloc(datafile.NrBins*sizeof(float));
  /* Must store rms'es and avrg'es in double precission, or else the
     baseline subtraction of the rms values will not be precise enough
     if there is a large DC components. */
  correl_double = (double *)malloc(datafile.NrSubints*sizeof(double));
  avrg_double = (double *)malloc(datafile.NrSubints*sizeof(double));
  if(data == NULL || correl_double == NULL || avrg_double == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR read_correlPSRData: Cannot allocate memory.");
    return 0;
  }
  
  for(i = 0; i < datafile.NrSubints; i++) {
    correl_double[i] = 0;
    avrg_double[i] = 0;
  }
  nrOffpulseBins = datafile.NrBins;
  if(regions != NULL) {
    for(j = 0; j < datafile.NrBins; j++) {
      if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0))
	nrOffpulseBins--;
    }
  }
  for(k = 0; k < datafile.NrFreqChan; k++) {
    for(i = 0; i < datafile.NrSubints; i++) {
      if(readPulsePSRData(&datafile, i, polchan, k, 0, datafile.NrBins, data, verbose) == 0)
	return 0;
      zap = 0;
      if(zapMask != NULL) {
	if(zapMask[i] != 0) {
	  /* printf("Z %ld\n", i); */
	  zap = 1;
	}
      }
      for(j = 0; j < datafile.NrBins; j++) {
	if(j == 0)
	  Ilast = 0;
	zap2 = zap;
	if(regions != NULL) {
	  if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0)) {
	    zap2 = 1;
	    Ilast = 0;
	  }
	}
	if(zap2 == 0) {
	  /*	  printf("I=%f\n", I); */
	  Ilast = data[j];
	  correl_double[i] += Ilast*data[j];
	  avrg_double[i] += data[j];
	  Ilast = data[j];
	}
      }
      /*
      l = i/10;
      if((10*l == i) && (l > 0)) {
	printf("read_correlPSRData: processed %.1f%%     \r",(100.0*i)/(float)(datafile.NrSubints));
	fflush(stdout);
	}*/
    }
  }
  for(i = 0; i < datafile.NrSubints; i++) {
    avrg_double[i] /= (float)(nrOffpulseBins*datafile.NrFreqChan);
    /* sum (I-a)^2 = sum (I*I -2aI + a*a) = sum (I*I) -2a*N*a + a*a*N = sum(I*I) - N*a*a */
    correl_double[i] -= nrOffpulseBins*datafile.NrFreqChan*avrg_double[i]*avrg_double[i];
    correl_double[i] /= (float)(nrOffpulseBins*datafile.NrFreqChan);
    correl[i] = sqrt(fabs(correl_double[i]));
  }
  /*  printf("  Done                                     \n");*/
  free(correl_double);
  free(avrg_double);
  free(data);
  return 1;
}

/* Read data file and produce a list of peak rms'es and averages
   (after adding addbins together). Set zapMask to NULL is you don't
   want to zap pulses. The defined regions are NOT included, unless
   invert is set. Set polchan to select the polarization channel. The
   defined regions are NOT included, unless invert is set.

   If return value = 0, then there is a memory allocation error.
*/
int read_peakrmsPSRData(datafile_definition datafile, float *rms, float *avrg, int *zapMask, pulselongitude_regions_definition *regions, int invert, int addbins, int polchan, verbose_definition verbose)
{
  long i, j, k, m, n;
  float *pulse, *pulse2, *data;
  double *avrg_double, *rms_double, *rms2_double, I, avrg_tmp;
  int nrOffpulseBins, zap, zap2;

  /*
// Not sure why this check would be relevant. Only put in so I didn't had to think of consequences?
  if(datafile.freqMode != FREQMODE_UNIFORM) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR read_peakrmsPSRData: Frequency channels need to be uniformely separated.");
    return 0;
  }
  */

  data = (float *)malloc(datafile.NrBins*sizeof(float));
  /* Must store rms'es and avrg'es in double precission, or else the
     baseline subtraction of the rms values will not be precise enough
     if there is a large DC components. */
  rms_double = (double *)malloc(datafile.NrSubints*sizeof(double));
  avrg_double = (double *)malloc(datafile.NrSubints*sizeof(double));
  if(data == NULL || rms_double == NULL || avrg_double == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR read_peakrmsPSRData: Cannot allocate memory.");
    return 0;
  }

  nrOffpulseBins = datafile.NrBins;
  if(regions != NULL) {
    for(j = 0; j < datafile.NrBins; j++) {
      if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0))
	nrOffpulseBins--;
    }
  }

  pulse  = (float *)malloc(datafile.NrBins*sizeof(float));
  pulse2 = (float *)malloc(datafile.NrBins*sizeof(float));
  rms2_double = (double *)malloc(datafile.NrSubints*sizeof(double));
  if(pulse == NULL || pulse2 == NULL || rms2_double == NULL) {
    fflush(stdout);
    printerror(verbose.debug, "ERROR read_peakrmsPSRData: Cannot allocate memory.");
    return 0;
  }

  /*  if(verbose.verbose) printf("Generating list of peak rms'es\n"); */
  for(i = 0; i < datafile.NrSubints; i++) {
    rms2_double[i] = 0;
    rms_double[i] = 0;
    avrg_double[i] = 0;
  }

  /*  int doplot=0;
  printf("plot?\n");
  scanf("%d", &doplot);
  */
  for(k = 0; k < datafile.NrFreqChan; k++) {
    for(i = 0; i < datafile.NrSubints; i++) {
      if(readPulsePSRData(&datafile, i, polchan, k, 0, datafile.NrBins, data, verbose) == 0)
	return 0;
      zap = 0;
      if(zapMask != NULL) {
	if(zapMask[i] != 0) 
	  zap = 1;
      }
      avrg_tmp = 0;
      for(j = 0; j < datafile.NrBins; j++) {
	zap2 = zap;
	if(regions != NULL) {
	  if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0)) {
	    zap2 = 1;
	  }
	}
	if(zap2 == 0) {
	  avrg_tmp += data[j];
	  pulse[j] = data[j];
	}else
	  pulse[j] = 0;
      }
      /* Replace zapped bins by average */
      avrg_tmp /= (float)(nrOffpulseBins);
      for(j = 0; j < datafile.NrBins; j++) {
	zap2 = zap;
	if(regions != NULL) {
	  if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0)) {
	    zap2 = 1;
	  }
	}
	if(zap2) 
	  pulse[j] = avrg_tmp;
      }

      if(addbins < 1) 
	addbins = 1;
      /*      if(doplot && i == 2578) {
	pgplotGraph1("2/xs", 0, 0, pulse, datafile.NrBins, 0, datafile.NrBins, 0, "Bins", "I", "pulse", NULL, 0, 0);
	printf("Adding %d\n", addbins);
	}*/
      for(j = 0; j < datafile.NrBins; j++) {
	pulse2[j] = 0;
	for(m = 0; m < addbins; m++) {
	  n = j+m;
	  if(n >= datafile.NrBins)
	    n -= datafile.NrBins;
	  pulse2[j] += pulse[n];
	}
      }
      /*      if(doplot && i == 2578)
	pgplotGraph1("3/xs", 0, 0, pulse2, datafile.NrBins, 0, datafile.NrBins, 0, "Bins", "I", "pulse2", NULL, 0, 0);
      */
      for(j = 0; j < datafile.NrBins; j++) {
	zap2 = 0;
	if(regions != NULL) {
	  if((checkRegions(j, regions, 0, verbose) != 0 && invert == 0) || (checkRegions(j, regions, 0, verbose) == 0 && invert != 0)) {
	    zap2 = 1;
	  }
	}
	if(zap2 == 0) {
	  if(pulse2[j]*pulse2[j] > rms_double[i])
	    rms_double[i] = pulse2[j]*pulse2[j];
	  rms2_double[i] += pulse2[j]*pulse2[j];
	  avrg_double[i] += pulse2[j];
	}
      }
      /*
      l = i/10;
      if((10*l == i) && (l > 0)) {
	printf("read_peakrmsPSRData: processed %.1f%%     \r",(100.0*i)/(float)(datafile.NrSubints));
	fflush(stdout);
	}*/
    }
  }
  for(i = 0; i < datafile.NrSubints; i++) {
    if(i == 2578) {
      /*      printf("%ld: %lf %lf %lf\n", i, avrg_double[i], rms_double[i], rms2_double[i]); */
    }
    /*   rms2 is the sigma^2 of the noise and rms is the sigma^2 of the peak */
    avrg_double[i] /= (float)(nrOffpulseBins*datafile.NrFreqChan);
    /* sum (I-a)^2 = sum (I*I -2aI + a*a) = sum (I*I) -2a*N*a + a*a*N = sum(I*I) - N*a*a */
    rms2_double[i] -= nrOffpulseBins*datafile.NrFreqChan*avrg_double[i]*avrg_double[i];
    rms2_double[i] /= (float)(nrOffpulseBins*datafile.NrFreqChan);
    /*
    if(i == 2578) {
      printf("%ld: %lf %lf %lf\n", i, avrg_double[i], rms_double[i], rms2_double[i]);
      } */
    /* (I-a)^2 = (I*I -2aI + a*a) */
    I = rms_double[i];
    rms_double[i] = I -2.0*avrg_double[i]*sqrt(I) + avrg_double[i]*avrg_double[i];
    rms_double[i] /= (float)(datafile.NrFreqChan);
    avrg[i] = avrg_double[i];
    if(zapMask != NULL) {
      if(zapMask[i] != 0) {
	rms[i] = 0;
      }else {
	if(rms2_double[i] != 0)
	  rms[i] = sqrt(rms_double[i]/rms2_double[i]);
	else
	  rms[i] = 0;
	if(rms_double[i]/rms2_double[i] < 0) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR read_peakrmsPSRData: Error (pulse %ld)", i);
	}
      }
    }else {
      if(rms2_double[i] != 0)
	rms[i] = sqrt(rms_double[i]/rms2_double[i]);
      else
	rms[i] = 0;
    }
  }
  free(pulse);
  free(pulse2);
  free(rms2_double);
  free(rms_double);
  free(avrg_double);
  free(data);
  /*  printf("  Done                                     \n"); */
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Parse the command line for -header or -headerUFL options and applies the changes
   to psrdata. Returns 0 if parse error */
int PSRDataHeader_parse_commandline(datafile_definition *psrdata, int argc, char **argv, verbose_definition verbose)
{
  int i, j, ok;
  char identifier[100], value[100], txt[100];

  ok = 0;
  for(i = 1; i < argc - 1; i++) {
    if(strcmp(argv[i], "-header") == 0 || strcmp(argv[i], "-headerUFL") == 0) {
      ok = 1;
    }
  }
  if(ok && verbose.verbose) {
    printf("Changing header parameters:\n");
  }
  for(i = 1; i < argc - 1; i++) {
    if(strcmp(argv[i], "-header") == 0) {
      i++;
      j = sscanf(argv[i], "%s %s", identifier, value);
      if(j != 2) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline:  Cannot parse -header option.");
	printHeaderCommandlineOptions(stderr);
	return 0;
      }else {
	if(strcasecmp(identifier,"name") == 0 || strcasecmp(identifier, "psrname") == 0 || strcasecmp(identifier, "pulsar") == 0 || strcasecmp(identifier, "psr") == 0) {
	  if(set_psrname_PSRData(psrdata, value, verbose) == 0) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Setting pulsar name failed.");
	    return 0;
	  }
	  if(verbose.verbose) printf("  hdr.psrname = %s\n", psrdata->psrname);
	}else if(strcasecmp(identifier, "freq") == 0 || strcasecmp(identifier, "cfreq") == 0 || strcasecmp(identifier, "freq_cent") == 0) {
	  if(psrdata->freqMode != FREQMODE_UNIFORM) {
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Centre frequency can only be changed in the header if the data has uniformly distributed channels.");
	      return 0;
	  }else {
	    double freq;
	    sscanf(value, "%lf", &freq);
	    set_centre_frequency(psrdata, freq, verbose);
	    if(verbose.verbose) printf("  hdr.freq_cent = %lf MHz\n", freq);
	  }
	}else if(strcasecmp(identifier, "bw") == 0) {
	  if(psrdata->freqMode != FREQMODE_UNIFORM) {
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Bandwidth can only be changed in the header if the data has uniformly distributed channels.");
	      return 0;
	  }else {
	    double bw;
	    sscanf(value, "%lf", &bw);
	    if(set_bandwidth(psrdata, bw, verbose) == 0) {
	      printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Bandwidth changing failed.");
	      return 0;
	    }
	    if(verbose.verbose) printf("  hdr.bw = %f MHz\n", bw);
	    double chanbw;
	    chanbw = bw/(double)psrdata->NrFreqChan;
	    if(verbose.verbose) printf("  hdr.channelbw = %f MHz\n", chanbw);
	  }
	}else if(strcasecmp(identifier,"chbw") == 0 || strcasecmp(identifier,"chanbw") == 0 || strcasecmp(identifier,"channelbw") == 0) {
	  if(psrdata->freqMode != FREQMODE_UNIFORM) {
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Bandwidth can only be changed in the header if the data has uniformly distributed samples.");
	      return 0;
	  }else {
	    double chanbw;
	    sscanf(value, "%lf", &chanbw);
	    double bw;
	    bw = chanbw*psrdata->NrFreqChan;
	    if(set_bandwidth(psrdata, bw, verbose) == 0) {
	      printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Bandwidth changing failed.");
	      return 0;
	    }
	    if(verbose.verbose) printf("  hdr.bw = %f MHz\n", bw);
	    if(verbose.verbose) printf("  hdr.channelbw = %f MHz\n", chanbw);
	  }
	}else if(strcasecmp(identifier, "freqref") == 0 || strcasecmp(identifier, "freq_ref") == 0 || strcasecmp(identifier, "ref_freq") == 0 || strcasecmp(identifier, "reffreq") == 0) {
	  psrdata->freq_ref = atof(value);
	  if(verbose.verbose) printf("  hdr.freq_ref = %f MHz\n", psrdata->freq_ref);
	  if(psrdata->freq_ref > -1.01 && psrdata->freq_ref < -0.99) {
	    psrdata->freq_ref = 1e10; // THe new preferred way to encode infinite frequency
	  }
	}else if(strcasecmp(identifier, "p0") == 0 || strcasecmp(identifier, "P0") == 0 || strcasecmp(identifier, "period") == 0) {
	  if(psrdata->isFolded != 1 || psrdata->foldMode != FOLDMODE_FIXEDPERIOD) {
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Period can only be changed in the header if there is a fixed period throughout the whole dataset.");
	      return 0;
	  }else {
	    psrdata->fixedPeriod = atof(value);
	    double period;
	    if(get_period(*psrdata, 0, &period, verbose) == 2) {
	      printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline (%s): Cannot obtain period", psrdata->filename);
	      return 0;
	    }
	    if(verbose.verbose) printf("  hdr.period = %lf s\n", period);
	  }
	}else if(strcasecmp(identifier, "dt") == 0 || strcasecmp(identifier, "tsamp") == 0 || strcasecmp(identifier, "samptime") == 0) {
	  if(psrdata->tsampMode != TSAMPMODE_FIXEDTSAMP) {
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Sampling time can only be changed in the header if there is a fixed sampling time throughout the whole dataset.");
	      return 0;
	  }else {
	    psrdata->fixedtsamp = atof(value);
	    if(verbose.verbose) printf("  hdr.SampTime = %lf s\n", get_tsamp(*psrdata, 0, verbose));
	  }
	}else if(strcasecmp(identifier, "tsub") == 0 || strcasecmp(identifier, "tsubint") == 0 || strcasecmp(identifier, "t_sub") == 0) {
	  char *substring;
	  int nrwords, ret, nsub;
	  substring = pickWordFromString(argv[i], 2, &nrwords, 1, ' ', verbose);
	  if(nrwords == 2) {   // Single value given
	    ret = sscanf(substring, "%lf", &(psrdata->tsub_list[0]));
	    if(ret != 1) {
	      printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Cannot parse %s as a value.", substring);
	    }
	    psrdata->tsubMode = TSUBMODE_FIXEDTSUB;
	    if(verbose.verbose) printf("  hdr.tsub = %lf s (fixed for each subint)\n", get_tsub(*psrdata, 0, verbose));
	  }else {
	    if(nrwords != 1+psrdata->NrSubints) {
	      printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: When setting individual subint durations, provide a list of space separated durations which defines a duration for each subint (got %d durations, need %d values). Example: -header 'tsub 30.0 30.0 30.0' when there are three subints present.", nrwords-1, psrdata->NrSubints);
	      return 0;
	    }
	    if(psrdata->tsub_list != NULL)
	      free(psrdata->tsub_list);
	    psrdata->tsub_list = (double *)malloc(psrdata->NrSubints*sizeof(double));
	    if(psrdata->tsub_list == NULL) {
	      printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Memory allocation error.");
	      return 0;
	    }
	    psrdata->tsubMode = TSUBMODE_TSUBLIST;
	    if(verbose.verbose) printf("  hdr.tsub = ");
	    for(nsub = 0; nsub < psrdata->NrSubints; nsub++) {
	      substring = pickWordFromString(argv[i], 2+nsub, &nrwords, 1, ' ', verbose);
	      ret = sscanf(substring, "%lf", &(psrdata->tsub_list[nsub]));
	      if(ret != 1) {
		printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Cannot parse %s as a value.", substring);
	      }
	      if(verbose.verbose) {
		if(nsub != 0)
		  printf(",");
		printf("%lf", get_tsub(*psrdata, nsub, verbose));
	      }
	    }
	    if(verbose.verbose) {
	      printf(" sec\n");
	    }
	  }
	}else if(strcasecmp(identifier,"mjd") == 0) {
	  psrdata->mjd_start = atof(value);
	  if(verbose.verbose) printf("  hdr.mjd = %Lf\n", psrdata->mjd_start);
	}else if(strcasecmp(identifier,"length") == 0 || strcasecmp(identifier,"tobs") == 0 || strcasecmp(identifier,"dur") == 0) {
	  psrdata->tsubMode = TSUBMODE_FIXEDTSUB;
	  if(psrdata->tsub_list != NULL)
	    free(psrdata->tsub_list);
	  psrdata->tsub_list = (double *)malloc(sizeof(double));
	  if(psrdata->tsub_list == NULL) {
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Memory allocation error.");
	    return 0;
	  }
	  psrdata->tsub_list[0] = atof(value)/(double)psrdata->NrSubints;
	  if(verbose.verbose) printf("  hdr.tobs = %lf\n", get_tobs(*psrdata, verbose));
	  printwarning(verbose.debug, "WARNING PSRDataHeader_parse_commandline: Assuming equal subint lengths.");
	}else if(strcasecmp(identifier,"loc") == 0 || strcasecmp(identifier,"location") == 0) {
	  //	  printf("XXXXX location = %s\n", value);
	  int ret;
	  ret = sscanf(value, "%lf,%lf,%lf", &(psrdata->telescope_X), &(psrdata->telescope_Y), &(psrdata->telescope_Z));
	  //	  if(ret != 3) {
	  //	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: In option %s of -header, '%s' expected to be of the form XVALUE,YVALUE,ZVALUE (without spaces).", identifier, value);
	  //	    return 0;	    
	  //	  }
	  if(ret != 3) {
	    if(verbose.verbose) {
	      printf("  Looking up ITRF coordinates for site '%s'\n", value);
	    }
	    if(setITRFlocation_by_name(psrdata, value, verbose) == 0) {
	      printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: '%s' not recognized as a telescope location. Consider entering location as XVALUE,YVALUE,ZVALUE (without spaces)", value);
	      return 0;
	    }
	  }
	  if(verbose.verbose) printf("  hdr.telescope_X = %lf m\n", psrdata->telescope_X);
	  if(verbose.verbose) printf("  hdr.telescope_Y = %lf m\n", psrdata->telescope_Y);
	  if(verbose.verbose) printf("  hdr.telescope_Z = %lf m\n", psrdata->telescope_Z);
	}else if(strcasecmp(identifier,"locationGEO") == 0 || strcasecmp(identifier,"locGEO") == 0) {
	  //	  printf("XXXXX location = %s\n", value);
	  int ret;
	  double longitude, latitude, height;
	  ret = sscanf(value, "%lf,%lf,%lf", &longitude, &latitude, &height);
	  if(ret != 3) {
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: In option %s of -header, '%s' expected to be of the form LONGITUDE,LATITUDE,HEIGHT (without spaces).", identifier, value);
	    return 0;	    
	  }
	  if(verbose.verbose) {
	    printf("  Derive ITRF X,Y,Z coordinates from geodetic GRS80 longitude=%f deg, latitude=%f deg, height=%f m\n", longitude, latitude, height);
	  }
	  longitude *= M_PI/180.0;
	  latitude *= M_PI/180.0;
	  tempo2_GRS80_to_ITRF(longitude, latitude, height, &(psrdata->telescope_X), &(psrdata->telescope_Y), &(psrdata->telescope_Z));
	  if(verbose.verbose) printf("  hdr.telescope_X = %lf m\n", psrdata->telescope_X);
	  if(verbose.verbose) printf("  hdr.telescope_Y = %lf m\n", psrdata->telescope_Y);
	  if(verbose.verbose) printf("  hdr.telescope_Z = %lf m\n", psrdata->telescope_Z);
	}else if(strcasecmp(identifier,"scan") == 0 || strcasecmp(identifier,"scanid") == 0) {
	  if(set_scanID_PSRData(psrdata, value, verbose) == 0) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Setting scan ID failed.");
	    return 0;
	  }
	  if(verbose.verbose) printf("  hdr.ScanID = %s\n", psrdata->scanID);
	}else if(strcasecmp(identifier,"project") == 0 || strcasecmp(identifier,"projectid") == 0 || strcasecmp(identifier,"projid") == 0) {
	  if(set_projectID_PSRData(psrdata, value, verbose) == 0) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Setting project ID failed.");
	    return 0;
	  }
	  if(verbose.verbose) printf("  hdr.projectID = %s\n", psrdata->projectID);
	}else if(strcasecmp(identifier, "observer") == 0 || strcasecmp(identifier, "observers") == 0) {
	  if(set_observer_PSRData(psrdata, value, verbose) == 0) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Setting name of observer(s) failed.");
	    return 0;
	  }
	  if(verbose.verbose) printf("  hdr.observer = %s\n", psrdata->observer);
	}else if(strcasecmp(identifier,"observatory") == 0 || strcasecmp(identifier,"telescope") == 0) {
	  if(set_observatory_PSRData(psrdata, value, verbose) == 0) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: Setting observatory name failed.");
	    return 0;
	  }
	  if(verbose.verbose) printf("  hdr.observatory = %s\n", psrdata->observatory);
	  int found;
	  found = 0;
	  for(j = 1; j < argc - 1; j++) {
	    sscanf(argv[j], "%s", txt);
	    if(strcasecmp(txt,"loc") == 0 || strcasecmp(txt,"location") == 0 || strcasecmp(txt,"locationGEO") == 0 || strcasecmp(txt,"locGEO") == 0) {
	      found = 1;
	      break;
	    }
	  }
	  if(found == 0) {
	    printwarning(verbose.debug, "WARNING PSRDataHeader_parse_commandline: When changing the name of the telescope, you might want to consider to change the location as well.");
	  }
	}else if(strcasecmp(identifier,"nrpulses") == 0 || strcasecmp(identifier,"npulses") == 0 || strcasecmp(identifier,"pulses") == 0 || strcasecmp(identifier,"nrsub") == 0  || strcasecmp(identifier,"nsub") == 0 || strcasecmp(identifier, "nsubint") == 0 || strcasecmp(identifier, "subints") == 0) {
	  psrdata->NrSubints = atol(value);
	  if(verbose.verbose) printf("  hdr.NrSubints = %ld\n", psrdata->NrSubints);
	}else if(strcasecmp(identifier,"nrbin") == 0 || strcasecmp(identifier,"nbin") == 0) {
	  psrdata->NrBins = atol(value);
	  if(verbose.verbose) printf("  hdr.NrBins = %ld\n", psrdata->NrBins);
	}else if(strcasecmp(identifier,"nrbits") == 0 || strcasecmp(identifier,"nbits") == 0) {
	  psrdata->NrBits = atoi(value);
	  if(verbose.verbose) printf("  hdr.NrBits = %d\n", psrdata->NrBits);
	}else if(strcasecmp(identifier,"nrchan") == 0 || strcasecmp(identifier,"nchan") ==  0 || strcasecmp(identifier,"nrfreq") == 0 || strcasecmp(identifier,"nfreq") == 0 || strcasecmp(identifier,"nrfreqchan") == 0 || strcasecmp(identifier,"nfreqchan") == 0) {
	  psrdata->NrFreqChan = atol(value);
	  if(verbose.verbose) printf("  hdr.NrFreqChan = %ld\n", psrdata->NrFreqChan);
	}else if(strcasecmp(identifier,"nrpol") == 0 || strcasecmp(identifier,"npol") == 0 || strcasecmp(identifier,"nrpols") == 0 || strcasecmp(identifier,"npols") == 0) {
	  psrdata->NrPols = atol(value);
	  if(verbose.verbose) printf("  hdr.NrPols = %ld\n", psrdata->NrPols);
	}else if(strcasecmp(identifier,"fdtype") == 0 || strcasecmp(identifier,"fd_type") == 0 || strcasecmp(identifier,"feedtype") == 0) {
	  psrdata->feedtype = atoi(value);
	  if(verbose.verbose) printf("  hdr.feedtype = %d\n", psrdata->feedtype);
	  //	}else if(strcasecmp(identifier,"fd_sang") == 0 || strcasecmp(identifier,"fdsang") == 0 || strcasecmp(identifier,"sang") == 0) {
	  //	  psrdata->fd_sang = atof(value);
	  //	  if(verbose.verbose) printf("  hdr.fd_sang = %f\n", psrdata->fd_sang);
	  //	}else if(strcasecmp(identifier,"xyph") == 0 || strcasecmp(identifier,"fd_xyph") == 0 || strcasecmp(identifier,"fdxyph") == 0) {
	  //	  psrdata->fd_xyph = atof(value);
	  //	  if(verbose.verbose) printf("  hdr.fd_xyph = %f\n", psrdata->fd_xyph);
	}else if(strcasecmp(identifier,"ra") == 0) {
	  if(strstr(value, ":") != NULL) {
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: In option %s of -header, a single value in degrees is expected.", identifier);
	    return 0;
	  }
	  psrdata->ra = atof(value)*M_PI/180.0;
	  if(verbose.verbose) printf("  hdr.ra = %f rad = %f deg = %f hours = ", psrdata->ra, psrdata->ra*180.0/M_PI, psrdata->ra*12.0/M_PI);
	  converthms_string(txt, psrdata->ra*12.0/M_PI, 2, 2);
	  printf("%s\n", txt);
	}else if(strcasecmp(identifier, "dec") == 0) {
	  if(strstr(value, ":") != NULL) {
	    printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: In option %s of -header, a single value in degrees is expected.", identifier);
	    return 0;
	  }
	  psrdata->dec = atof(value)*M_PI/180.0;
	  if(verbose.verbose) printf("  hdr.dec = %f rad = %f deg = ", psrdata->dec, psrdata->dec*180.0/M_PI);
	  converthms_string(txt, psrdata->dec*180.0/M_PI, 2, 3);
	  printf("%s\n", txt);
	}else if(strcasecmp(identifier,"poltype") == 0 || strcasecmp(identifier,"pol_type") == 0) {
	  psrdata->poltype = atoi(value);
	  if(verbose.verbose) printf("  hdr.poltype = %d\n", psrdata->poltype);
	}else if(strcasecmp(identifier,"dedisp") == 0 || strcasecmp(identifier,"dedispersed") == 0 || strcasecmp(identifier,"isdedisp") == 0 || strcasecmp(identifier,"isdedispersed") == 0) {
	  psrdata->isDeDisp = atoi(value);
	  if(verbose.verbose) printf("  hdr.isDeDisp = %d\n", psrdata->isDeDisp);
	}else if(strcasecmp(identifier,"defarad") == 0 || strcasecmp(identifier,"isdefarad") == 0) {
	  psrdata->isDeFarad = atoi(value);
	  if(verbose.verbose) printf("  hdr.isDeFarad = %d\n", psrdata->isDeFarad);
	  if(fabs(psrdata->rm) < 1e-6) {
	    fflush(stdout);
	    printwarning(verbose.debug, "WARNING PSRDataHeader_parse_commandline: Be cautious with changing the de-Faraday rotation state while the RM is not defined.");
	  }
	}else if(strcasecmp(identifier,"depar") == 0 || strcasecmp(identifier,"isdepar") == 0) {
	  psrdata->isDePar = atoi(value);
	  if(verbose.verbose) printf("  hdr.isDePar = %d\n", psrdata->isDePar);
	}else if(strcasecmp(identifier,"debase") == 0 || strcasecmp(identifier,"isdebase") == 0) {
	  psrdata->isDebase = atoi(value);
	  if(verbose.verbose) printf("  hdr.isDebase = %d\n", psrdata->isDebase);
	}else if(strcasecmp(identifier,"dm") == 0) {
	  psrdata->dm = atof(value);
	  if(verbose.verbose) printf("  hdr.dm = %f\n", psrdata->dm);
	}else if(strcasecmp(identifier,"rm") == 0) {
	  psrdata->rm = atof(value);
	  if(verbose.verbose) printf("  hdr.rm = %f\n", psrdata->rm);
	}else if(strcasecmp(identifier,"cableswap") == 0) {
	  psrdata->cableSwap = atoi(value);
	  if(verbose.verbose) printf("  hdr.cableswap = %d\n", psrdata->cableSwap);
	}else if(strcasecmp(identifier,"cableswapcor") == 0) {
	  psrdata->cableSwapcor = atoi(value);
	  if(verbose.verbose) printf("  hdr.cableswapcor = %d\n", psrdata->cableSwapcor);
	}else if(strcasecmp(identifier,"gentype") == 0 || strcasecmp(identifier,"type") == 0) {
	  int orig_gentype = psrdata->gentype;
	  sscanf(value, "%d", &(psrdata->gentype));
	  if(verbose.verbose) printf("  hdr.gentype = %d\n", psrdata->gentype);
	  if(orig_gentype == GENTYPE_SEARCHMODE) {
	    psrdata->isFolded = 1;
	    psrdata->foldMode = FOLDMODE_FIXEDPERIOD;
	    if(psrdata->fixedPeriod <= 0)
	      psrdata->fixedPeriod = 1.0;
	    fflush(stdout);
	    printwarning(verbose.debug, "WARNING PSRDataHeader_parse_commandline:  Period is not set. Use -header \"period value\" to set it to something appropriate");
	  }
	  if(psrdata->gentype == GENTYPE_SEARCHMODE) {
	    psrdata->fixedPeriod = 0;
	    psrdata->foldMode = FOLDMODE_UNKNOWN;
	    psrdata->isFolded = 0;
	  }
	}else if(strcasecmp(identifier,"yrange") == 0) {
	  //	  printf("XXXXX location = %s\n", value);
	  int ret;
	  double value1, value2;
	  ret = sscanf(value, "%lf,%lf", &value1, &value2);
	  if(ret == 2) {
	    psrdata->yrangeset = 1;
	    psrdata->yrange[0] = value1;
	    psrdata->yrange[1] = value2;
	    if(verbose.verbose) printf("  hdr.yrange1 = %lf\n", psrdata->yrange[0]);
	    if(verbose.verbose) printf("  hdr.yrange2 = %lf\n", psrdata->yrange[1]);
	  }else {
	    if(strcasecmp(value, "x") == 0 || strcasecmp(value, "undefined") == 0 || strcasecmp(value, "empty") == 0 || strcasecmp(value, "nothing") == 0) {
	      psrdata->yrangeset = 0;
	      if(verbose.verbose) printf("  hdr.yrange = undefined\n");
	    }else {
	      printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline: In option %s of -header, '%s' expected to be of the form VALUE1,VALUE2 (without spaces) or the word UNDEFINED.", identifier, value);
	      return 0;	    	      
	    }
	  }
	}else {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR PSRDataHeader_parse_commandline:  '%s' not recognized as a header parameter.", identifier);
	  printHeaderCommandlineOptions(stderr);
	  return 0;	  
	}
      }
    }else if(strcmp(argv[i], "-headerUFL") == 0) {
      force_uniform_frequency_spacing(psrdata, verbose);
    }
  }
  if(ok && verbose.verbose) {
    verbose_definition verbose2;
    copyVerboseState(verbose, &verbose2);
    verbose2.verbose = 0;
    printHeaderPSRData(*psrdata, 1, verbose2);
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

void printHeaderCommandlineOptions(FILE *printdevice)
{
  fprintf(printdevice, "Valid options for the -header option are:\n");
  fprintf(printdevice, "  bw           band width.\n");
//START REGION DEVELOP
  fprintf(printdevice, "  cableswap    Cables swapped during obs (0 or 1).\n");
  fprintf(printdevice, "  cableswapcor Cable swap corrected (0 or 1).\n");
//START REGION RELEASE
  fprintf(printdevice, "  chbw         Channel band width.\n");
  fprintf(printdevice, "  dec          Declination (single number, in degrees).\n");
  fprintf(printdevice, "  dm           Dispersion Measure.\n");
  fprintf(printdevice, "  dt           Sampling time.\n");
  fprintf(printdevice, "  fdtype       1=lin 2=circ, sign is handiness.\n");
  fprintf(printdevice, "  freq         Centre frequency (i.e. defines frequency labeling).\n");
  fprintf(printdevice, "  debase       Data is baseline-subtracted (0 or 1).\n");
  fprintf(printdevice, "  dedisp       Data is de-dispersed (0 or 1).\n");
  fprintf(printdevice, "  defarad      Data is de-faraday rotated (0 or 1).\n");
  fprintf(printdevice, "  depar        Data is parallactic angle corrected (0 or 1).\n");
  fprintf(printdevice, "  gentype      Identifies the type of data, use -gentypelist to see options.\n");
  fprintf(printdevice, "  length       The duration of the observation.\n");
  fprintf(printdevice, "  location     The ITRF location of telescope in meters, e.g.\n");
  fprintf(printdevice, "               -header \"loc 3822252.643,-153995.683,5086051.443\" for the Lovell.\n");  
  fprintf(printdevice, "               or -header \"loc Lovell\" for recognized site names.\n");  
  fprintf(printdevice, "  locationGEO  Derive ITRF location of telescope via geodetic GRS80 coordinates:\n");
  fprintf(printdevice, "               -header \"locGEO long,lat,height\" where longitude/latitude are.\n");  
  fprintf(printdevice, "               in degrees, and height in meters.\n");  
  fprintf(printdevice, "  mjd          Start MJD of the observation.\n");
  fprintf(printdevice, "  name         Name of the pulsar.\n");
  fprintf(printdevice, "  nbin         Nr of time-bins (use at own risk).\n");
  fprintf(printdevice, "  nbits        Nr of bits (use at own risk).\n");
  fprintf(printdevice, "  nchan        Nr of frequency channels (use at own risk).\n");
  fprintf(printdevice, "  npol         Nr of polarization channels (use at own risk).\n");
  fprintf(printdevice, "  nsub         Nr of sub integrations in observation (use at own risk).\n");
  fprintf(printdevice, "  observatory  Change the telescope (name only, not location).\n");
  fprintf(printdevice, "  p0           Period.\n");
  fprintf(printdevice, "  poltype      Polarization type (%d = undefined, %d=Stokes, %d=Coherency,\n", POLTYPE_UNKNOWN, POLTYPE_STOKES, POLTYPE_COHERENCY);
  fprintf(printdevice, "               %d=I,L,V,Pa+error, %d=Pa+error, %d=I,L,V,Pa+error,tot pol,ell+error).\n", POLTYPE_ILVPAdPA, POLTYPE_PAdPA, POLTYPE_ILVPAdPATEldEl);
  fprintf(printdevice, "  ra           Right ascension (single number, in degrees).\n");
  fprintf(printdevice, "  reffreq      Reference frequency (for dedispersion/de Faraday rotation).\n");
  fprintf(printdevice, "               -1=Infinite freq -2=Unknown\n");
  fprintf(printdevice, "  rm           Rotation Measure.\n");
//START REGION DEVELOP
//  fprintf(printdevice, "  sang         Angle of receiver.\n");
//START REGION RELEASE
  fprintf(printdevice, "  scan         Scan ID of the observation.\n");
  fprintf(printdevice, "  tsub         Followed by single number: set subint durations to this number\n");
  fprintf(printdevice, "               in seconds.\n");
  fprintf(printdevice, "               If followed by space separated numbers: set individual subint\n");
  fprintf(printdevice, "               durations\n");
  fprintf(printdevice, "  yrange       Set the yrange to the range provided. Specify this as\n");
  fprintf(printdevice, "               -header \"yrange VALUE1,VALUE2\". To unset the use of\n");
  fprintf(printdevice, "               this parameter use -header \"yrange undefined\".\n");

//START REGION DEVELOP
//  fprintf(printdevice, "  xyph         Other angle of receiver.\n");
//START REGION RELEASE
}

//START REGION RELEASE

// Shows the options of the -header gentype option
void printHeaderGentypeOptions(FILE *printdevice)
{
  fprintf(printdevice, "Valid options for the -header \"gentype number\" option are:\n");
  fprintf(printdevice, "  %3d           Undefined.\n", GENTYPE_UNDEFINED);
  fprintf(printdevice, "  %3d           Profile.\n", GENTYPE_PROFILE);
  fprintf(printdevice, "  %3d           Pulse stack.\n", GENTYPE_PULSESTACK);
  fprintf(printdevice, "  %3d           Subintegration data.\n", GENTYPE_SUBINTEGRATIONS);
  fprintf(printdevice, "  %3d           Seach mode data (not folded).\n", GENTYPE_SEARCHMODE);
  fprintf(printdevice, "  %3d           Bandpass.\n", GENTYPE_BANDPASS);
  fprintf(printdevice, "  %3d           Dynamic spectrum.\n", GENTYPE_DYNAMICSPECTRUM);
  fprintf(printdevice, "  %3d           Polarization calibration signal (noise diode).\n", GENTYPE_POLNCAL);
  fprintf(printdevice, "  %3d           LRFS.\n", GENTYPE_LRFS);
  fprintf(printdevice, "  %3d           2DFS.\n", GENTYPE_2DFS);
  fprintf(printdevice, "  %3d           Sliding 2DFS (P3).\n", GENTYPE_S2DFSP3);
  fprintf(printdevice, "  %3d           Sliding 2DFS (P2).\n", GENTYPE_S2DFSP2);
  fprintf(printdevice, "  %3d           P3 fold.\n", GENTYPE_P3FOLD);
//START REGION DEVELOP
  fprintf(printdevice, "  %3d           Polar map.\n", GENTYPE_POLARMAP);
  fprintf(printdevice, "  %3d           LRCC.\n", GENTYPE_LRCC);
  fprintf(printdevice, "  %3d           RM map.\n", GENTYPE_RMMAP);
  fprintf(printdevice, "  %3d           Position angle distribution.\n", GENTYPE_PADIST);
  fprintf(printdevice, "  %3d           Ellipticity angle distribution.\n", GENTYPE_ELLDIST);
//START REGION RELEASE
  fprintf(printdevice, "  %3d           Receiver model.\n", GENTYPE_RECEIVERMODEL);
  fprintf(printdevice, "  %3d           Receiver model with chi^2 and nfree.\n", GENTYPE_RECEIVERMODEL2);
}

//START REGION DEVELOP

/* Allocates memory for paswing, pulse longitudes and bins and RMS
   values. Returns 0 if not successful */
/*
int alloc_paswing_PSRdata(datafile_definition *psrdata, int paswing, int longitudes, int bins, int rmsValues, verbose_definition verbose)
{
  int debug;
  debug = 0;
  //  if(psrdata->NrSubints > 1) {
  //      fflush(stdout);
  //    printerror(verbose.debug, "ERROR alloc_paswing_PSRdata: Requires only one pulse!");
  //    return 0;
  //  }
  //  if(psrdata->NrFreqChan != 1) {
  //      fflush(stdout);
  //    printerror(verbose.debug, "ERROR alloc_paswing_PSRdata: Requires only one frequency channel!");
  //    return 0;
  //  }

  if(paswing) {
    psrdata->data_pa = (float *)malloc((psrdata->NrSubints*psrdata->NrFreqChan*psrdata->NrBins)*sizeof(float));
    psrdata->data_dpa = (float *)malloc((psrdata->NrSubints*psrdata->NrFreqChan*psrdata->NrBins)*sizeof(float));  
    if(psrdata->data_pa == NULL || psrdata->data_dpa == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR alloc_paswing: Cannot allocate memory"); 
      return 0;    
    } 
    psrdata->NrPApoints = psrdata->NrSubints*psrdata->NrFreqChan*psrdata->NrBins;
  }
  if(longitudes) {
    psrdata->data_long = (float *)malloc((psrdata->NrBins)*sizeof(float));
    if(psrdata->data_long == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR alloc_paswing: Cannot allocate memory"); 
      return 0;    
    } 
    psrdata->longitudes_defined = 1;
  }
  if(bins) {
    psrdata->data_bin = (long int *)malloc((psrdata->NrBins)*sizeof(long int));
    if(psrdata->data_bin == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR alloc_paswing: Cannot allocate memory"); 
      return 0;    
    } 
    psrdata->bins_defined = 1;
  }
  if(rmsValues) {
    psrdata->rmspoints_defined = 3;
    psrdata->rmsvalues = (float *)malloc(psrdata->rmspoints_defined*(psrdata->NrSubints*psrdata->NrFreqChan)*sizeof(float));
    if(psrdata->rmsvalues == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR alloc_paswing: Cannot allocate memory"); 
      return 0;    
    } 
  }
  return 1;
}
*/

//START REGION DEVELOP
//START REGION RELEASE

/* This function makes a clone of the data and puts it in a (as yet
   non-existent) clone. So the clone will contain have all data
   available in memory. Returns 1 on success, 0 on error. The data
   must have been read into memory first. */
int make_clone(datafile_definition original, datafile_definition *clone, verbose_definition verbose)
{
  int debug;
  debug = 0;
  cleanPSRData(clone, verbose);
  if(original.format != MEMORY_format) {
    fflush(stdout);
    printerror(debug, "ERROR make_clone: Only works if data is already loaded into memory.");
    return 0;
  }
  //  if(original.NrPApoints > 0) {
  //    fflush(stdout);
  //    printerror(debug, "ERROR make_clone: Cannot handle PA data.");
  //    return 0;
  //  }

  copy_params_PSRData(original, clone, verbose);

  clone->data = (float *)malloc((original.NrBins)*(original.NrPols)*(original.NrFreqChan)*(original.NrSubints)*sizeof(float));
  if(clone->data == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR make_clone: Memory allocation error.");
    return 0;
  }
  memcpy(clone->data, original.data, (original.NrBins)*(original.NrPols)*(original.NrFreqChan)*(original.NrSubints)*sizeof(float));

  if(original.offpulse_rms != NULL) {
    clone->offpulse_rms = (float *)malloc(original.NrPols*original.NrFreqChan*original.NrSubints*sizeof(float));
    if(clone->offpulse_rms == NULL) {
      fflush(stdout);
      printerror(debug, "ERROR make_clone: Memory allocation error.");
      return 0;
    }
    memcpy(clone->offpulse_rms, original.offpulse_rms, original.NrPols*original.NrFreqChan*original.NrSubints*sizeof(float));
  }

  clone->opened_flag = 1;
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Moves the clone to original and free up the original. */
void swap_orig_clone(datafile_definition *original, datafile_definition *clone, verbose_definition verbose)
{
  verbose_definition verbose2;
  copyVerboseState(verbose, &verbose2);
  verbose2.verbose = 0;
  verbose2.nocounters = 1;
  closePSRData(original, 0, verbose);  // Release all memory in original, and replace all header information with that in the clone
// USED TO DO 1 AND SAY: Perserves the header information (pointers), as that will be passed on to the clone, but release memory. Not entirely sure where the close is happening though and why before the move.
  memmove(original, clone, sizeof(datafile_definition));
}

//START REGION DEVELOP
//START REGION RELEASE

/* Writes the history to the file, including the present command line
   (if argc > 0). If cmdOnly is set, only the command is writen
   out. This makes the history table without timestamps and therefore
   results are exactly reproducable. The file should be opened with
   write permission AND THE HEADER SHOULD ALREADY BE WRITTEN OUT AND
   NOT BE MODIFIED ANYMORE.  Returns 1 on success, 0 on error */
int writeHistoryPSRData(datafile_definition *datafile, int argc, char **argv, int cmdOnly, verbose_definition verbose)
{
  int ret;
  char txt[10000], txt2[1000], *username_ptr, hostname[1000]; //, username[1000]
  time_t curtime;
  datafile_history_entry_definition *curHistoryEntry;

  if(verbose.verbose)
    fprintf(stdout, "Writing history\n");

  if(argc > 0) {
    // Add a new entry if required
    curHistoryEntry = &(datafile->history);
    // Check if the first entry is not empty. If so, we need to find last entry and add a new one to it.
    if(curHistoryEntry->timestamp != NULL || curHistoryEntry->cmd != NULL || curHistoryEntry->user != NULL || curHistoryEntry->hostname != NULL || curHistoryEntry->nextEntry != NULL) {
      // Find last entry
      while(curHistoryEntry->nextEntry != NULL) {
	curHistoryEntry = curHistoryEntry->nextEntry;
      }
      // Add entry
      curHistoryEntry->nextEntry = malloc(sizeof(datafile_history_entry_definition));
      if(curHistoryEntry->nextEntry == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error"); 
	return 0;
      }
      curHistoryEntry = curHistoryEntry->nextEntry; // Goto next entry and make it empty
      curHistoryEntry->timestamp = NULL;
      curHistoryEntry->cmd = NULL;
      curHistoryEntry->user = NULL;
      curHistoryEntry->hostname = NULL;
      curHistoryEntry->nextEntry = NULL;
    }

    constructCommandLineString(txt, 10000, argc, argv, verbose);
    curHistoryEntry->cmd = malloc(strlen(txt)+1);
    if(curHistoryEntry->cmd == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error"); 
      return 0;
    }
    strcpy(curHistoryEntry->cmd, txt);


    /*
      txt[0] = 0;
      for(i = 0; i < argc; i++) {
      if(strlen(txt) + strlen(argv[i]) + 4 > 9999) {
      printwarning(verbose.debug, "WARNING writeHistoryPSRData: Truncating command line which is too long.");
      break;
      }
      if(strchr(argv[i], ' ') == NULL) {
      strcat(txt, argv[i]);
      }else {
      strcat(txt, "\"");
      strcat(txt, argv[i]);
      strcat(txt, "\"");
      }
      if(i != argc-1)
      strcat(txt, " ");
      }
    */
    /*  strcpy(txt2, asctime(gmtime(time(NULL)))); */
    if(cmdOnly == 0) {
      curtime = time(NULL);
      strcpy(txt2, asctime(gmtime(&curtime))); 
      if(txt2[strlen(txt2)-1] == '\n')
	txt2[strlen(txt2)-1] = 0;
      if(txt2[strlen(txt2)-1] == '\r')
	txt2[strlen(txt2)-1] = 0;
      if(txt2[strlen(txt2)-1] == '\n')
	txt2[strlen(txt2)-1] = 0;
      curHistoryEntry->timestamp = malloc(strlen(txt2)+1);
      if(curHistoryEntry->timestamp == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error"); 
	return 0;
      }
      strcpy(curHistoryEntry->timestamp, txt2);
      
      /*  printf("Time stamp: '%s'\n", txt2); */
      
      //      username_ptr = who_am_i();
      if(getUsername(&username_ptr, verbose) == 0) {
	//strcpy(username, username_ptr);
	//      }else {
	fflush(stdout); 
	printwarning(verbose.debug, "writeHistoryPSRData: Cannot identify user.");
	username_ptr = malloc(8);
	if(username_ptr == NULL) {
	  printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error");
	  return 0;
	}
	sprintf(username_ptr, "Unknown");
      }
      curHistoryEntry->user = malloc(strlen(username_ptr)+1);
      if(curHistoryEntry->user == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error"); 
	return 0;
      }
      strcpy(curHistoryEntry->user, username_ptr);
      free(username_ptr);
      
      getMachinename(hostname, 1000, verbose);
      curHistoryEntry->hostname = malloc(strlen(hostname)+1);
      if(curHistoryEntry->hostname == NULL) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR writeHistoryPSRData: Memory allocation error"); 
	return 0;
      }
      strcpy(curHistoryEntry->hostname, hostname);
      
      /*  printf("Machine: %s\n", hostname); */
    
      //  printf("XXXXX '%s'\nXXXXXXX '%s'\n", txt, txt2);
    }
  }

  if(datafile->format == PSRSALSA_BINARY_format) {
    ret = writeHistoryPSRSALSA(datafile, verbose);
  }else if(datafile->format == FITS_format) {
    ret = writeHistoryFITS(*datafile, verbose);
  }else if(datafile->format == PUMA_format) {
    ret = writeHistoryPuma(*datafile, verbose);
  }else {
    fflush(stdout);
    printwarning(verbose.debug, "writeHistoryPSRData: Writing a history is not supported in this file format.");
    ret = 0;
  }
  if(verbose.verbose)
    fprintf(stdout, "  done\n");
  return ret;
}

//START REGION RELEASE

/* Reads the content of the history that is stored in the file.
   Returns:
     1  on success
     0  on error when history should be supported for format.
     -1 no history was read, since it isn't supported for this format. It should therefore not result in termination of program. 
*/
int readHistoryPSRData(datafile_definition *datafile, verbose_definition verbose)
{
  int ret, indent, doread;
  if(verbose.verbose) {
    for(indent = 0; indent < verbose.indent; indent++)      
      printf(" ");
    fprintf(stdout, "Reading history\n");
  }
  //    printf("History line %d: ", rownr+1);
  // date    printf("%s ", txt2);
    // user   printf("%s ", txt2);
    // hostname   printf("%s ", txt2);
    // cmd   printf("%s\n", txt2);


  doread = 0;
  if(datafile->format == FITS_format) {
    ret = readHistoryFITS(datafile, verbose);
    doread = 1;
  }else if(datafile->format == PUMA_format) {
    ret = readHistoryPuma(datafile, verbose);
    doread = 1;
  }else {
    //    fflush(stdout);
    //    printerror(verbose.debug, "ERROR readHistoryPSRData: Reading the history is not supported in this file format.");
    ret = -1;
  }

  if(doread) {
    if(ret == 0) {
      printwarning(verbose.debug, "WARNING: Reading history failed.");
    }
    
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fprintf(stdout, "  done\n");
    }
  }else {
    if(verbose.verbose) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      fprintf(stdout, "  skipped\n");
    }
  }
  return ret;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Show the content of the history that is stored in the file.
   Returns 1 on success, 0 on error, for instance no history supported. */
int showHistory(datafile_definition datafile, verbose_definition verbose)
{
  datafile_history_entry_definition *curHistoryEntry;
  long rownr;
  int indent;

  curHistoryEntry = &(datafile.history);
  rownr = 0;
  do {
    if(curHistoryEntry->timestamp != NULL || curHistoryEntry->cmd != NULL || curHistoryEntry->user != NULL || curHistoryEntry->hostname != NULL || curHistoryEntry->nextEntry != NULL) {
      for(indent = 0; indent < verbose.indent; indent++)      
	printf(" ");
      printf("History line %ld: ", rownr+1);
      if(curHistoryEntry->timestamp != NULL)
	printf("%s ", curHistoryEntry->timestamp);
      if(curHistoryEntry->user != NULL)
	printf("%s ", curHistoryEntry->user);
      if(curHistoryEntry->hostname != NULL)
	printf("%s ", curHistoryEntry->hostname);
      if(curHistoryEntry->cmd != NULL)
	printf("%s\n", curHistoryEntry->cmd);
      curHistoryEntry = curHistoryEntry->nextEntry; // Goto next entry
      rownr++;
    }
  }while(curHistoryEntry != NULL && rownr != 0);  // Keep on going as long as the next history line is defined. Quit when the first entry does not contain any information, as otherwise this would be an infinite loop.
  if(rownr > 0)
    return 1;
  else 
    return 0;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Rewinds the file and skil lines starting with a #. The start
   datastart field in the datafile struct is adjusted. */
int skipallhashedlines(datafile_definition *datafile)
{
  long pos;
  char c;
  int ret, debug;
  debug = 0;
  rewind(datafile->fptr);
  do {
    pos = ftell(datafile->fptr);
    ret = fscanf(datafile->fptr, "%c", &c);
    if(ret != 1) {
      fflush(stdout);
      printerror(debug, "ERROR skipallhashedlines: Cannot read data.");
      return 0;
    }
    if(c == '#') {
      do {
	fscanf(datafile->fptr, "%c", &c);
      }while(c != '\n');
      c = '#';
    }
  }while(c == '#');
  datafile->datastart = (long long)pos;
  fseek(datafile->fptr, datafile->datastart, SEEK_SET);  
  return 1;
}


//START REGION DEVELOP
//START REGION RELEASE

/* The string text will be parsed and certain keywords will replaced
   with header parameters. The return value is a pointer of the new
   string (memory is allocated, so it should be freed after use), or
   NULL if an error occured. Use the function str_list_replace_keys()
   to list all possible keywords. */
char *str_replace_header_params(datafile_definition data, char *text, verbose_definition verbose)
{
  int i, debug;
  char *newtext, *newtext2, headerparam[1000];

  if(verbose.debug) {
    fflush(stdout);
    printf("Entering str_replace_header_params()\n");
  }

  debug = 0;

  sprintf(headerparam, "%Lf", data.mjd_start);
  newtext = str_replace(text, "%MJD", headerparam, verbose);
  if(newtext == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }

  sprintf(headerparam, "%ld", data.NrSubints);
  newtext2 = str_replace(newtext, "%NRSUBINTS", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;

  sprintf(headerparam, "%ld", data.NrFreqChan);
  newtext2 = str_replace(newtext, "%NRFREQCHAN", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;

  double period;
  int ret;
  if(data.isFolded) {
    ret = get_period(data, 0, &period, verbose);
    if(ret == 2) {
      printerror(verbose.debug, "ERROR str_replace_header_params (%s): Cannot obtain period", data.filename);
      return 0;
    }
  }else {
    ret = 1;
    period = 0;
  }
  sprintf(headerparam, "%lf", period);
  newtext2 = str_replace(newtext, "%PERIOD", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;

  if(data.freqMode != FREQMODE_UNKNOWN) {
    sprintf(headerparam, "%lf", get_centre_frequency(data, verbose));
  }
  newtext2 = str_replace(newtext, "%FREQ", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;

  if((data.freq_ref > -1.1 && data.freq_ref < -0.9) || (data.freq_ref > 0.99e10 && data.freq_ref < 1.01e10)) {
    sprintf(headerparam, "Inf");
  }else if(data.freq_ref >= 0) {
    sprintf(headerparam, "%lf", data.freq_ref);
  }else {
    sprintf(headerparam, "Unknown");
  }
  newtext2 = str_replace(newtext, "%REFFREQ", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;

  if(data.freqMode != FREQMODE_UNKNOWN) {
    sprintf(headerparam, "%lf", get_bandwidth(data, verbose));
  }
  newtext2 = str_replace(newtext, "%BW", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;

  sprintf(headerparam, "%lf", data.dm);
  newtext2 = str_replace(newtext, "%DM", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;

  sprintf(headerparam, "%lf", data.rm);
  newtext2 = str_replace(newtext, "%RM", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;

//START REGION DEVELOP
  if(data.cableSwapcor == 1)
    sprintf(headerparam, "YES");
  else if(data.cableSwapcor == 0)
    sprintf(headerparam, "NO");
  else
    sprintf(headerparam, "UNKNOWN");
  newtext2 = str_replace(newtext, "%CABLESWAPCOR", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
//START REGION RELEASE

  strncpy(headerparam, data.filename, 999);
  for(i = 0; i < strlen(headerparam); i++) {
    if(headerparam[i] == '.') {
      headerparam[i] = 0;
      break;
    }
  }
  newtext2 = str_replace(newtext, "%FILEBASE", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;

  strncpy(headerparam, data.filename, 999);
  newtext2 = str_replace(newtext, "%FILE", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;

  strncpy(headerparam, data.psrname, 999);
  newtext2 = str_replace(newtext, "%PSR", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;

//START REGION DEVELOP
  if(data.cableSwap == 1)
    sprintf(headerparam, "YES");
  else if(data.cableSwap == 0)
    sprintf(headerparam, "NO");
  else
    sprintf(headerparam, "UNKNOWN");
  newtext2 = str_replace(newtext, "%CABLESWAP", headerparam, verbose);
  if(newtext2 == NULL) {
    fflush(stdout);
    printerror(debug, "ERROR str_replace_header_params: Cannot replace text");
    return NULL;
  }
  free(newtext);
  newtext = newtext2;
//START REGION RELEASE

  if(verbose.debug) {
    fflush(stdout);
    printf("Exiting str_replace_header_params()\n");
  }

  return newtext;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Lists the available keywords for the str_replace_header_params()
   function. nrspaces spaces are printed in front of the output. */
void str_list_replace_keys(int nrspaces)
{
  int i;
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%BW             - Bandwidth in MHz\n");
//START REGION DEVELOP
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%CABLESWAP      - Original observation affected by cable swap?\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%CABLESWAPCOR   - Data corrected for cable swap?\n");
//START REGION RELEASE
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%DM             - DM\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%FILE           - File name\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%FILEBASE       - File name (everything before first .)\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%FREQ           - Centre frequency in MHz\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%MJD            - MJD\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%NRFREQCHAN     - Nr of frequency channels\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%NRSUBINTS      - Nr of subintegrations\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%PERIOD         - Period in seconds\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%PSR            - Name of the pulsar\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%RM             - RM\n");
  if(nrspaces > 0)
    for(i = 0; i < nrspaces; i++)
      printf(" ");
  printf("%%REFFREQ        - Reference frequency for dedispersion/de-Faraday rotation in MHz\n");
}


// Set the struct to the default state (no verbose)
void cleanVerboseState(verbose_definition *verbose_state)
{
  verbose_state->verbose = 0;
  verbose_state->debug = 0;
  verbose_state->nocounters = 0;
  verbose_state->indent = 0;
}

// Copy verbose state
void copyVerboseState(verbose_definition verbose_state_src, verbose_definition *verbose_state_dst)
{
  verbose_state_dst->verbose = verbose_state_src.verbose;
  verbose_state_dst->debug = verbose_state_src.debug;
  verbose_state_dst->nocounters = verbose_state_src.nocounters;
  verbose_state_dst->indent = verbose_state_src.indent;
}


/* If the data is in FREQMODE_FREQTABLE mode, check if the data really
   makes use of this feature. If all the channels are equally spaced
   in line with the given bandwidth and centre-frequency, then change
   the data to FREQMODE_UNIFORM, which is compatible with more
   functionality. 

   Return -
   0: Error
   1: Not converted, because input data not in FREQMODE_FREQTABLE mode
   2: Not converted, because frequency channels are not uniformly separated
   3: Converted
*/

int convert_if_uniform_frequency_spacing(datafile_definition *datafile, int nowarnings, verbose_definition verbose) 
{
  long subint, freqchan;
  double freq, freq_first, freq_prev, df_expected, actual_centre_freq;
  if(verbose.debug) {
    printf("Check if data has uniformly distributed channels\n");
  }
  if(datafile->freqMode != FREQMODE_FREQTABLE) {
    if(verbose.debug) {
      printf("  Data does not have frequency defined for separate channels/subintegrations, so nothing is done\n");
    }
    return 1;
  }
  
  // Alright, find out if uniformly separated. Done by checking that all channels are equally spaced.
  freq_first = df_expected = actual_centre_freq = 0.0;
  for(subint = 0; subint < datafile->NrSubints; subint++) {
    freq_prev = 0.0;
    for(freqchan = 0; freqchan < datafile->NrFreqChan; freqchan++) {
      freq = get_weighted_channel_freq(*datafile, subint, freqchan, verbose);
      if(subint == 0)
	actual_centre_freq += freq;   // Keep track of cumulative frequency to be used later.
      if(subint == 0 && freqchan == 0) {
	freq_first = freq;   // Frequency of the first channel
	df_expected = -freq;
      }else if(subint == 0 && freqchan == 1) {
	df_expected += freq;          // Define channel separation of first two channels of first two subints as the expected separation
	freq_prev = freq;
      }else if(freqchan == 0 && subint != 0) { // Not possible to determine actual frequency separation
	if(fabs(freq - freq_first) > 1e-6) {
	  if(nowarnings == 0) {
	    printwarning(verbose.debug, "WARNING (%s): Frequency channels do not appear to be equally spaced (%lf != %lf for first channels of subint 0 and %ld). Some operations might might not work with this data. Warnings for other channels are suppressed.", datafile->filename, freq, freq_first, subint);
	  }
	  return 2;
	}
	freq_prev = freq;
      }else {
	if(fabs(freq-freq_prev - df_expected) > 1e-6) {  // Check for equally spacing
	  if(nowarnings == 0) {
	    printwarning(verbose.debug, "WARNING (%s): Frequency channels do not appear to be equally spaced (%lf - %lf != %lf for channel %ld and subint %ld). Some operations might might not work with this data. Warnings for other channels are suppressed.", datafile->filename, freq, freq_prev, df_expected, freqchan, subint);
	  }
	  return 2;
	}
	freq_prev = freq;
      }
    }
  }

  // All channels are equally spaced and aligned between different subints. So can change to FREQMODE_UNIFORM.
  if(verbose.debug) {
    printf("  The frequency channels appear to be uniformly separated, so convert in data-set described by only a bandwidth and a centre frequency.\n");
  }
  actual_centre_freq /= (double)datafile->NrFreqChan;
  double hdr_freq_cent;
  hdr_freq_cent = get_centre_frequency(*datafile, verbose);

  // Check if centre frequency needs to be adjusted
  if(fabs(actual_centre_freq-hdr_freq_cent) > 1e-6) {  // Apparently the centre frequency is different
    if(hdr_freq_cent > 0 && (hdr_freq_cent < 0.99e10 || hdr_freq_cent > 1.01e10)) {     // There is already a centre frequency defined there should be warnings, otherwise just do the change

      // zero'th channel is assumed to be DC and therefore should dropped in calculation of OBSFREQ
      double chanbw;
      if(get_channelbandwidth(*datafile, &chanbw, verbose) == 0) {
	printerror(verbose.debug, "ERROR (%s): Cannot obtain channel bandwidth.", datafile->filename);
	return 0;
      }
      if(fabs(hdr_freq_cent-actual_centre_freq-0.5*fabs(chanbw)) < 1e-6) {
	if(verbose.debug) {
	  printf("  (%s): Updating centre frequency from %lf to %lf MHz as suggested by the frequency list in the fits table, corresponding to an offset expected from a dropped DC channel.\n", datafile->filename, hdr_freq_cent, actual_centre_freq);
	}
      }else {
	if(nowarnings == 0) {
	  printwarning(verbose.debug, "WARNING (%s): Updating centre frequency from %lf to %lf MHz as suggested by the frequency list in the data.", datafile->filename, hdr_freq_cent, actual_centre_freq);
	}
      }
    }else {
      if(verbose.debug)
	printf("  %s: Use subint frequency table frequency as centre frequency = %f MHz.\n", datafile->filename, actual_centre_freq);
    }
    set_centre_frequency(datafile, actual_centre_freq, verbose);
  }

  datafile->freqMode = FREQMODE_UNIFORM;
  free(datafile->freqlabel_list);
  datafile->freqlabel_list = NULL;


  if(datafile->NrFreqChan > 1) {
    df_expected = (freq_prev-freq_first)/(double)(datafile->NrFreqChan-1);  // Possibly a marginal more precise calculation
    double chanbw;
    if(get_channelbandwidth(*datafile, &chanbw, verbose) == 0) {
      printerror(verbose.debug, "ERROR (%s): Cannot obtain channel bandwidth.", datafile->filename);
      return 0;
    }
    if(fabs(df_expected-chanbw) > 1e-6) {  //Need to update bandwidth
      fflush(stdout);
      if(nowarnings == 0) {
	printwarning(verbose.debug, "WARNING (%s): Updating channel bandwidth from %lf to %lf MHz suggested by the frequency list in the subint table.", datafile->filename, chanbw, df_expected);
      }
      if(set_bandwidth(datafile, df_expected*datafile->NrFreqChan, verbose) == 0) {
	printerror(verbose.debug, "ERROR (%s): Bandwidth changing failed.", datafile->filename);
	return 0;
      }
    }
  }

  return 3;
}

/* If the data is in FREQMODE_FREQTABLE mode, check if the data really
   makes use of this feature. If all the channels are equally spaced
   in line with the given bandwidth and centre-frequency, then change
   the data to FREQMODE_UNIFORM, which is compatible with more
   functionality. 

   Return -
   0: Error
   1: Not converted, because input data not in FREQMODE_FREQTABLE mode
   2: Not converted, because frequency channels are not uniformly separated
   3: Converted
*/

/*
  If the data is in FREQMODE_FREQTABLE mode, force all channels to
   become equally spaced with a single bandwidth and
   centre-frequency. Then change the data to FREQMODE_UNIFORM, which
   is compatible with more functionality. NO CHECK IS DONE IF THIS IS
   REASONALBLE!

   Return -
   1: Not converted, because input data not in FREQMODE_FREQTABLE mode
   3: Converted
 */
int force_uniform_frequency_spacing(datafile_definition *datafile, verbose_definition verbose) 
{
  if(verbose.debug) {
    printf("Force channels to be uniformly distributed in frequency\n");
  }
  if(datafile->freqMode != FREQMODE_FREQTABLE) {
    if(verbose.debug) {
      printf("  Data does not have frequency defined for separate channels/subintegrations, so nothing is done\n");
    }
    return 1;
  }

  free(datafile->freqlabel_list);
  datafile->freqlabel_list = NULL;
  datafile->freqMode = FREQMODE_UNIFORM;

  printwarning(verbose.debug, "WARNING (%s): Re-labelling frequency channels to be uniformly distributed and forcing them to be the same for each subintegration. Expect all frequency dependent effects to be wrong.", datafile->filename);

  return 3;
}

/*
    double cfreq, dfreq, df_expected, freqlast;
    double freq;
    int i, convert, suppresswarning;
    suppresswarning = 0;
    cfreq = freqlast = 0;
 	

    if(convert) {
      if(fabs(cfreq-hdr_freq_cent) > 1e-6) {
	if(hdr_freq_cent > 0 && (hdr_freq_cent < 0.99e10 || hdr_freq_cent > 1.01e10)) {
	  // zero'th channel is assumed to be DC and therefore should dropped in calculation of OBSFREQ
	  double chanbw;
	  if(get_channelbandwidth(*datafile, &chanbw, verbose) == 0) {
	    printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot obtain channel bandwidth.", datafile->filename);
	    exit(0);
	  }
	  if(fabs(hdr_freq_cent-cfreq-0.5*fabs(chanbw)) < 1e-6) {
	    if(verbose.debug) {
	      printf("  readPSRFITSHeader (%s): Updating centre frequency from %lf to %lf MHz as suggested by the frequency list in the fits table, corresponding to an offset expected from a dropped DC channel.\n", datafile->filename, hdr_freq_cent, cfreq);
	    }
	  }else {
	    printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Updating centre frequency from %lf to %lf MHz as suggested by the frequency list in the fits table.", datafile->filename, hdr_freq_cent, cfreq);
	  }
	}else {
	  //	    printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Set centre frequency to %lf MHz as suggested by the frequency list in the subint table.", datafile->filename, cfreq);
	  if(verbose.debug)
	    printf("  readPSRFITSHeader (%s): Use subint frequency table frequency as centre frequency = %f MHz.\n", datafile->filename, cfreq);
	}
	datafile->freqMode = FREQMODE_UNIFORM;
	set_centre_frequency(datafile, cfreq, verbose);
      }



      if(datafile->NrFreqChan > 1) {
	dfreq = (freq_prev-freq_first)/(double)(datafile->NrFreqChan-1);
	double chanbw;
	if(get_channelbandwidth(*datafile, &chanbw, verbose) == 0) {
	  printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Cannot obtain channel bandwidth.", datafile->filename);
	  exit(0);
	}
	if(fabs(dfreq-chanbw) > 1e-6) {
	  fflush(stdout);
	  printwarning(verbose.debug, "WARNING readPSRFITSHeader (%s): Updating channel bandwidth from %lf to %lf MHz suggested by the frequency list in the subint table.", datafile->filename, chanbw, dfreq);
	  if(set_bandwidth(datafile, dfreq*datafile->NrFreqChan, verbose) == 0) {
	    printerror(verbose.debug, "ERROR readPSRFITSHeader (%s): Bandwidth changing failed.", datafile->filename);
	    exit(0);
	  }
	}
      }
    }
  }
}
*/

//START REGION DEVELOP

/* If the data is constant and zero for all polarizations in a given
   single pulse (or subint) in a given frequency channel, then
   iszapped is set to 1, otherwise it is set to 0.

   Returns 1 if successful, or zero if an error occured.
 */
int pulse_iszapped(datafile_definition datafile, long subintnr, int freq, int *iszapped, verbose_definition verbose)
{
  long p, b;
  float sample, min, max;
  for(p = 0; p < datafile.NrPols; p++) {
    for(b = 0; b < datafile.NrBins; b++) {
      if(readPulsePSRData(&datafile, subintnr, p, freq, b, 1, &sample, verbose) != 1) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR pulse_iszapped: Cannot read data.");
	return 0;
      }
      if(b == 0 && p == 0) {
	min = max = sample;
      }else {
	if(sample < min)
	  min = sample;
	if(sample > max)
	  max = sample;
      }
    }
  }
  if(min == 0.0 && max == 0.0) {
    *iszapped = 1;
  }else {
    *iszapped = 0;
  }
  return 1;
}

/* Using pulse_iszapped(), find the number of zapped frequency
   channels for a given subint.

   Returns 1 if successful, or zero if an error occured.
 */
int nrchannels_iszapped(datafile_definition datafile, long subintnr, long *nrzapped, verbose_definition verbose)
{
  long f;
  int iszapped;
  *nrzapped = 0;
  for(f = 0; f < datafile.NrFreqChan; f++) {
    if(pulse_iszapped(datafile, subintnr, f, &iszapped, verbose) == 0) {
      printerror(verbose.debug, "ERROR nrchannels_iszapped: Cannot determine if frequency channel is zapped.");
      return 0;
    }
    if(iszapped) {
      (*nrzapped) += 1;
    }
  }
  return 1;
}

