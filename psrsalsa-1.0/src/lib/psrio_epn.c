//START REGION RELEASE

#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <string.h>
#include <math.h>
#include "psrsalsa.h"

//Each pulse starts with a long header (written by write_epn_longheader())
//Each channel starts with a short header (written by write_epn_shortheader())
//The data is written as 4 hex bytes = 2^16 levels. Each channel is written as a multiple of 80 bytes, possibly ending with zero's


// Removes any white spaces at end of string
void stripspaces(char *txt)
{
  while(strlen(txt) > 0) {
    if(txt[strlen(txt)-1] == ' ')
      txt[strlen(txt)-1] = 0;
    else
      break;
  }
}

void fillspaces(char *txt, int n)
{
  memset(txt, ' ', n);
  txt[n] = 0;
}

double parse_unit_freq(char *txt, verbose_definition verbose)
{
  int i;
  stripspaces(txt);
  for(i = 0; i < 8; i++) {
    if(txt[i] != ' ')
      break;
  }
  if(i >= 7)
    i = 0;
  if(strcasecmp(&txt[i], "MHz") == 0) { // Nothing needs to be done, already in MHz
    return 1;
  }else if(strcasecmp(&txt[i], "GHz") == 0) { // Nothing needs to be done, already in MHz
    return 1000;
  }else {
    printwarning(verbose.debug, "WARNING readEPNsubHeader: Unknown frequency unit '%s', assume it is in GHz.", txt);
    return 1000;
  }
}

// version can be 6.0 or 6.3 (*100)
void write_epn_longheader(datafile_definition datafile, int version, verbose_definition verbose)
{
  int counter, nrlines, i;
  double dbl;
  char txt[100];

  nrlines = datafile.NrBins/20;     // Nr of 80 byte lines of data per channel
  if(datafile.NrBins > nrlines*20)
    nrlines++;

  counter=6+datafile.NrPols*datafile.NrFreqChan*(nrlines+2);  // Nr of 80 byte lines per pulse
  // First line
  if(version == 600) {
    fillspaces(txt, 68);   // History
    sprintf(txt, "  Written using libpsrsalsa");
    txt[strlen(txt)] = ' ';
    txt[68] = 0;
    fprintf(datafile.fptr, "EPN 6.00%4d%s", counter, txt);
  }else {
    fillspaces(txt, 66);   // History
    sprintf(txt, "  Written using libpsrsalsa");
    txt[strlen(txt)] = ' ';
    txt[66] = 0;
    fprintf(datafile.fptr, "EPN 6.30%6d%s", counter, txt);
  }

  // second line
  // get the name without the PSR prefix
  fillspaces(txt, 12);   // Jname
  if(strncmp(datafile.psrname, "PSR ", 4) == 0) {
    strncpy(txt, &datafile.psrname[4], 10);
  }else if(strncmp(datafile.psrname, "PSR_", 4) == 0 ) {
    strncpy(txt, &datafile.psrname[4], 10);
  }else if(strncmp(datafile.psrname, "PSR", 3) == 0 ) {
    strncpy(txt, &datafile.psrname[3], 10);
  }else {
    strncpy(txt, &datafile.psrname[0], 10);
  }
  txt[12] = 0;
  fprintf(datafile.fptr, "%12s", txt); // Jname
  fprintf(datafile.fptr, "%12s", txt); // Common nname

  i = get_period(datafile, 0, &dbl, verbose);
  if(i == 2) {
    printerror(verbose.debug, "ERROR write_epn_longheader (%s): Cannot obtain period", datafile.filename);
    exit(0);
  }else if (i == 1) {
    printwarning(verbose.debug, "WARNING write_epn_longheader: Error obtaining the fold period");
  }
  if(dbl < 0) {
    printwarning(verbose.debug, "WARNING write_epn_longheader: Error obtaining the fold period");
    dbl = 0;
  }
  sprintf(txt, "%16.12lf", dbl);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%8.3f", datafile.dm);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%10.3f", datafile.rm);
  fprintf(datafile.fptr, "%s", txt);
  fillspaces(txt, 6);                  // catref
  fprintf(datafile.fptr, "%s", txt);
  fillspaces(txt, 8);                  // bibref
  fprintf(datafile.fptr, "%s", txt);
  fillspaces(txt, 8);                  // Expansion
  fprintf(datafile.fptr, "%s", txt);

  // Third line
  converthms_string(txt, datafile.ra*12.0/M_PI, 3, 4);
  fprintf(datafile.fptr, "%s", txt);
  if(datafile.ra >= 0)
    fprintf(datafile.fptr, "+"); // No + sign by default of converthms_string()
  converthms_string(txt, datafile.dec*180.0/M_PI, 3, 4);
  fprintf(datafile.fptr, "%s", txt);

  fillspaces(txt, 8);
  if(datafile.observatory != NULL) {
    strncpy(txt, datafile.observatory, 8);
    //    datafile.observatory[strlen(datafile.observatory)] = ' ';
    //    txt[8] = 0;
  }
  fprintf(datafile.fptr, "%8s", txt);

  sprintf(txt, "%10.3Lf", datafile.mjd_start);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%8.3f", 0.0);  // opos = position angle
  fprintf(datafile.fptr, "%s", txt);
  fprintf(datafile.fptr, " ");   // PA flag (absolute PA)
  fprintf(datafile.fptr, " ");   // TIM flag (absolute time stamp UTC)
  fillspaces(txt, 31);                // Expansion
  fprintf(datafile.fptr, "%s", txt);

  // Fourth line
  sprintf(txt, "%17.5lf", datafile.telescope_X);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%17.5lf", datafile.telescope_Y);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%17.5lf", datafile.telescope_Z);
  fprintf(datafile.fptr, "%s", txt);
  fillspaces(txt, 29);                // Expansion
  fprintf(datafile.fptr, "%s", txt);

  // Fifth line
  sprintf(txt, "01121978");            // Creation date
  fprintf(datafile.fptr, "%s", txt);

  i = 0;
  sscanf(datafile.scanID, "%d", &i);
  sprintf(txt, "%04d", i);
  fprintf(datafile.fptr, "%s", txt);
  if(version == 600) {
    sprintf(txt, "%04d", 0);
    fprintf(datafile.fptr, "%s", txt);
  }else {
    sprintf(txt, "%05d", 0);
    fprintf(datafile.fptr, "%s", txt);
  }
  sprintf(txt, "%02ld", datafile.NrPols);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%04ld", datafile.NrFreqChan);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%04ld", datafile.NrBins);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%12.6lf", 1e6*get_tsamp(datafile, 0, verbose));
  fprintf(datafile.fptr, "%s", txt);
  fprintf(datafile.fptr, "%s", txt);  // Set tres = tbin
  sprintf(txt, "%06d", 1);            // Set the number of pulses per integration = 1 (as each pulse is a separate block)
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%04d", 0);            // Bin number start cal signal
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%04d", 0);            // Bin number length cal signal
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, " ");                  // Not flux calibrated
  fprintf(datafile.fptr, "%s", txt);
  if(version == 600) {
    fillspaces(txt, 15);                // Expansion
    fprintf(datafile.fptr, "%s", txt);
  }else {
    fillspaces(txt, 14);                // Expansion
    fprintf(datafile.fptr, "%s", txt);
  }
  fprintf(datafile.fptr, "--------------------------------------------------------------------------------");
}


//idfield, 'I' for example, will be truncated at 8 chars
//nband = Ordinal number of current stream
void write_epn_shortheader(datafile_definition datafile, int version, char *idfield, int nband, verbose_definition verbose)
{
  char txt[100];

  if(datafile.freqMode != FREQMODE_UNIFORM) {
    printerror(verbose.debug, "ERROR write_epn_shortheader (%s): Frequency channels do not appear to be uniformly distributed.", datafile.filename);
    exit(0);
  }

  fillspaces(txt, 8);
  strncpy(txt, idfield, 8);
  txt[8] = 0;
  fprintf(datafile.fptr, "%s", txt);

  sprintf(txt, "%04d", nband);            // Bin number length cal signal
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%04d", 1);                // Nr of streams averaged
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%12.8f", get_centre_frequency(datafile, verbose)*0.001);                
  fprintf(datafile.fptr, "%s", txt);
  if(version >= 610) {
    sprintf(txt, "%8s", "GHz");                
    fprintf(datafile.fptr, "%s", txt);
  }
  sprintf(txt, "%12.6f", get_bandwidth(datafile, verbose));                
  fprintf(datafile.fptr, "%s", txt);
  if(version >= 610) {
    sprintf(txt, "%8s", "MHz");                
    fprintf(datafile.fptr, "%s", txt);
  }
  sprintf(txt, "%17.5f", 0.0); // Tstart in us after epoch
  fprintf(datafile.fptr, "%s", txt);
  if(version >= 620) {
    fprintf(datafile.fptr, " ");  // plflag
    sprintf(txt, "%6.1f", 0.0);       // plval = %6.1f
    fprintf(datafile.fptr, "%s", txt);
  }else if(version == 610) {
    fillspaces(txt, 7);
    fprintf(datafile.fptr, "%s", txt);
  }else if(version == 600) {
    fillspaces(txt, 23);
    fprintf(datafile.fptr, "%s", txt);
  }
}

void write_epn_data(datafile_definition datafile, float scale, float offset, float rms, unsigned int *iprofile, verbose_definition verbose)
{
  int i;
  char txt[100];
  double dbl;

  sprintf(txt, "%#12G", scale);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%#12G", offset);
  fprintf(datafile.fptr, "%s", txt);
  sprintf(txt, "%#12G", rms);
  fprintf(datafile.fptr, "%s", txt);

  i = get_period(datafile, 0, &dbl, verbose);
  if(i == 2) {
    printerror(verbose.debug, "ERROR write_epn_data (%s): Cannot obtain period", datafile.filename);
    exit(0);
  }else if(i == 1) {
    printwarning(verbose.debug, "WARNING write_epn_data: Error obtaining the fold period");
  }
  if(dbl < 0) {
    printwarning(verbose.debug, "WARNING write_epn_data: Error obtaining the fold period");
    dbl = 0;
  }
  sprintf(txt, "%16.12lf", dbl);
  fprintf(datafile.fptr, "%s", txt);

  fillspaces(txt, 28);
  fprintf(datafile.fptr, "%s", txt);

  for(i = 0; i < datafile.NrBins; i++)
    fprintf(datafile.fptr, "%04X", iprofile[i]);
  for(i = 0; i < 20-(datafile.NrBins-(datafile.NrBins/20)*20); i++)
    fprintf(datafile.fptr, "0000");
}

/* Based on inverted write_epn_data.
ReadLargeHeader_psrsalsa() the same as what=1, except stored in datafile.
This is the main header
ReadSmallHeader_psrsalsa() the same as what=2, except stored in datafile. version must be set om datafile.
subheader
Returns 0 on error, 1 if sucess.
 */
int readEPNHeader(datafile_definition *datafile, verbose_definition verbose)
{
  int nlines, ret, i;
  double dbl_version;
  char txt[1000];

  //    verbose.debug = 1;

  /* first line */
    ret = fread(txt, 1, 3, datafile->fptr);
    txt[3] = 0;
    if(ret != 3) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    if(strcmp(txt, "EPN") != 0) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 5, datafile->fptr);
    txt[5] = 0;
    if(ret != 5) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%lf", &dbl_version);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    datafile->version = round(dbl_version*10);
    if(verbose.debug) {
      printf("DEBUG readEPNHeader: EPN version = %lf (%d)\n", dbl_version, datafile->version);
    }
  if((datafile->version < 60) || (datafile->version > 63)) {
    printerror(verbose.debug, "ERROR readEPNHeader: EPN versions 6.0 - 6.3 are supported. This is version %lf.", dbl_version);
    return 0;
  }

    int counter;
    if(datafile->version == 63) {
      ret = fread(txt, 1, 6, datafile->fptr);
      txt[6] = 0;
      if(ret != 6) {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
      ret = sscanf(txt, "%d", &counter);
      if(ret != 1) {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
      /*  printf("readEPNHeader: Counter = %d\n", counter); */
      // This is a history line
      ret = fread(txt, 1, 66, datafile->fptr);
      txt[66] = 0;
      if(ret != 66) {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
    }else {
      ret = fread(txt, 1, 4, datafile->fptr);
      txt[4] = 0;
      if(ret != 4) {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
      ret = sscanf(txt, "%d", &counter);
      if(ret != 1) {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
      /*  printf("readEPNHeader: Counter = %d\n", counter); */
      // This is a history line
      ret = fread(txt, 1, 68, datafile->fptr);
      txt[68] = 0;
      if(ret != 68) {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
    }
    
    // Second line
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    stripspaces(txt);
    for(i = 0; i < 12; i++) {
      if(txt[i] != ' ')
	break;
    }
    if(i >= 11)
      i = 0;
    if(set_psrname_PSRData(datafile, &txt[i], verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readEPNHeader: Setting pulsar name failed.");
      return 0;
    }
    // The common name, rather than the J-name above
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 16, datafile->fptr);
    txt[16] = 0;
    if(ret != 16) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    datafile->isFolded = 1;
    datafile->foldMode = FOLDMODE_FIXEDPERIOD;
    ret = sscanf(txt, "%lf", &datafile->fixedPeriod);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 8, datafile->fptr);
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%lf", &datafile->dm);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 10, datafile->fptr);
    txt[10] = 0;
    if(ret != 10) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%lf", &datafile->rm);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 6, datafile->fptr); //catref
    txt[6] = 0;
    if(ret != 6) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 8, datafile->fptr); // bibref
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 8, datafile->fptr); // Expansion
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    int hour, min;
    double sec, sign;
    ret = fread(txt, 1, 10, datafile->fptr); // Position ra=(2h, 2m, 6s) and dec=(1sign,2d,2m,6s)
    txt[10] = 0;
    if(ret != 10) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    sscanf(txt, "%2d%2d%lf", &hour, &min, &sec);
    datafile->ra = (hour+(min+sec/60.0)/60.0)*M_PI/12.0;

    ret = fread(txt, 1, 11, datafile->fptr); // Position ra=(2h, 2m, 6s) and dec=(1sign,2d,2m,6s)
    txt[11] = 0;
    if(ret != 11) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    sscanf(txt, "%3d%2d%lf", &hour, &min, &sec);
    if(hour < 0) {
      sign = -1.0;
      hour = -hour;
    }else {
      sign = +1.0;
    }
    datafile->dec = sign*(hour+(min+sec/60.0)/60.0)*M_PI/180.0;
    
    ret = fread(txt, 1, 8, datafile->fptr);
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    for(i = 0; i < 8; i++) {
      if(txt[i] != ' ')
	break;
    }
    if(i >= 7)
      i = 0;
    if(set_observatory_PSRData(datafile, &txt[i], verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readEPNHeader: Setting observatory name failed.");
      return 0;
    }
    
    ret = fread(txt, 1, 10, datafile->fptr);
    txt[10] = 0;
    if(ret != 10) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%Lf", &datafile->mjd_start);
    if(ret != 1) {
      datafile->mjd_start = 0;
    }
    ret = fread(txt, 1, 8, datafile->fptr); // opos
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 1, datafile->fptr); // paflag
    txt[1] = 0;
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 1, datafile->fptr); // timflag
    txt[1] = 0;
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 31, datafile->fptr); // Expansion
    txt[31] = 0;
    if(ret != 31) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 17, datafile->fptr); // Xtel
    txt[17] = 0;
    if(ret != 17) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 17, datafile->fptr);  // Ytel
    txt[17] = 0;
    if(ret != 17) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 17, datafile->fptr);  // Ztel
    txt[17] = 0;
    if(ret != 17) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 29, datafile->fptr);  // Expansion
    txt[29] = 0;
    if(ret != 29) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    
    // Line 4:
    
    /* This is the day/month/year/scanno/subscan */
    if(datafile->version > 60) {
      ret = fread(txt, 1, 17, datafile->fptr);
      txt[17] = 0;
      if(ret != 17) {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
    }else {
      ret = fread(txt, 1, 16, datafile->fptr);
      txt[16] = 0;
      if(ret != 16) {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
    }
    if(set_scanID_PSRData(datafile, &txt[8], verbose) == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readEPNHeader: Setting scan ID failed.");
      return 0;
    }
    ret = fread(txt, 1, 2, datafile->fptr);
    txt[2] = 0;
    if(ret != 2) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%ld", &datafile->NrPols);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    if(datafile->NrPols <= 0) {
      printwarning(verbose.debug, "WARNING readEPNHeader: The number of polarizations (%ld) does not make sense. I assume it should be one???", datafile->NrPols);
      datafile->NrPols = 1;
    }
    datafile->NrBits = 16;  // True for all EPN files (4 bytes, but 2^16 levels as is hex)
    ret = fread(txt, 1, 4, datafile->fptr);
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%ld", &datafile->NrFreqChan);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 4, datafile->fptr);
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%ld", &datafile->NrBins);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 12, datafile->fptr);
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    datafile->tsampMode = TSAMPMODE_FIXEDTSAMP;
    datafile->tsubMode = TSUBMODE_FIXEDTSUB;
    if(datafile->tsub_list != NULL)
      free(datafile->tsub_list);
    datafile->tsub_list = (double *)malloc(sizeof(double));
    if(datafile->tsub_list == NULL) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR readEPNHeader: Memory allocation error");
      return 0;
    }
    //    ret = sscanf(txt, "%lf", datafile->tsub_list);
    ret = sscanf(txt, "%lf", &(datafile->fixedtsamp));
    datafile->tsub_list[0] = datafile->fixedPeriod;
    datafile->fixedtsamp *= 1e-6;
    //  printf("xxxxx %s\n", txt);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 12, datafile->fptr);  // tres
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 6, datafile->fptr);
    txt[6] = 0;
    if(ret != 6) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%ld", &datafile->NrSubints);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 4, datafile->fptr); // ncal
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 4, datafile->fptr);  // lcal
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 1, datafile->fptr);  // fluxflag
    txt[1] = 0;
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
      return 0;
    }
    
    // Expansion
    if(datafile->version > 60) {
      ret = fread(txt, 1, 14, datafile->fptr);
      txt[14] = 0;
      if(ret != 14) {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
    }else {
      ret = fread(txt, 1, 15, datafile->fptr);
      txt[15] = 0;
      if(ret != 15) {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
    }
    
    for(i = 0; i < 80; i++) {
      ret = fread(txt, 1, 1, datafile->fptr);
      if(ret != 1) {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
      if(txt[0] != '-') {
	printerror(verbose.debug, "ERROR readEPNHeader: This file is not an EPN file.");
	return 0;
      }
    }
    
    // Is this ok? what if nrbins = multiple 20
    nlines=datafile->NrBins/20+1;
    if(counter != 6+datafile->NrPols*datafile->NrFreqChan*(nlines+2)) {
      printerror(verbose.debug, "ERROR readEPNHeader: Something wrong with EPN file, counter doesn't make sense (%d != %ld). npol=%ld nfreq=%ld nbin=%ld", counter, 6+datafile->NrPols*datafile->NrFreqChan*(nlines+2), datafile->NrPols, datafile->NrFreqChan, datafile->NrBins);
      return 0;
    }


    // Determine nr of subints
    int length_of_block;
    long size;
    length_of_block = counter*80;
    datafile->datastart = ftell(datafile->fptr);
    fseek(datafile->fptr, 0, SEEK_END);
    size = ftell(datafile->fptr);
    fseek(datafile->fptr, datafile->datastart, SEEK_SET);
    if(size < 0) {
      printerror(verbose.debug, "ERROR readEPNHeader: Determining file size failed.");
      return 0;
    }
    datafile->NrSubints = size/length_of_block;
    if(size != (long)datafile->NrSubints*(long)length_of_block) {
      printerror(verbose.debug, "ERROR readEPNHeader: File size does not match an integer number of subints.");
      return 0;
    }
  return 1;
}



/* Based on inverted write_epn_data.
Returns 0 on error, 1 if sucess.
 */
int readEPNsubHeader(datafile_definition *datafile, float *scale, float *offset, verbose_definition verbose)
{
  int ret, version;
  double dbl_version;
  char txt[1000];
  double freq, bw;

  //  verbose.debug = 1;
  version = datafile->version;
  version = 63;  // Overwrite to make compatible with files I generate, which are a mixture of 6.0 and 6.3 to make compatible with psrchive.

  /* first line */
  if((version < 60) || (version > 63)) {
    dbl_version = (version)*0.1;
    printerror(verbose.debug, "ERROR readEPNsubHeader: EPN versions 6.0 - 6.3 are supported. This is version %lf.", dbl_version);
    return 0;
  }

    ret = fread(txt, 1, 8, datafile->fptr);  //idfield
    txt[8] = 0;
    if(ret != 8) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 4, datafile->fptr);  //nband
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 4, datafile->fptr);  //navg
    txt[4] = 0;
    if(ret != 4) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 12, datafile->fptr);  //f0
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    if(datafile->freqMode == FREQMODE_UNKNOWN) {
      ret = sscanf(txt, "%lf", &freq);   // Here it is assumed there is only one freq present
      if(ret != 1) {
	printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
	return 0;
      }
    }
    if(version > 60) {   // uf
      ret = fread(txt, 1, 8, datafile->fptr);
      txt[8] = 0;
      if(ret != 8) {
	printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
	return 0;
      }
      stripspaces(txt);
    }else {
      txt[0] = 0;
    }
    if(datafile->freqMode == FREQMODE_UNKNOWN) {
      freq *= parse_unit_freq(txt, verbose);
      //      datafile->freqMode = FREQMODE_UNIFORM;   // Do not set this, otherwise BW will not be set
    }


    ret = fread(txt, 1, 12, datafile->fptr);  //df
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    if(datafile->freqMode == FREQMODE_UNKNOWN) {
      ret = sscanf(txt, "%lf", &bw);   // Here it is assumed there is only one freq present
      if(ret != 1) {
	printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
	if(verbose.debug) {
	  printerror(verbose.debug, "ERROR readEPNsubHeader: Cannot parse '%s' as a number.", txt);	  
	}
	return 0;
      }
    }
    if(version > 60) {   // ud
      ret = fread(txt, 1, 8, datafile->fptr);
      txt[8] = 0;
      if(ret != 8) {
	printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
	return 0;
      }
    }else {
      txt[0] = 0;
    }
    if(datafile->freqMode == FREQMODE_UNKNOWN) {
      bw *= parse_unit_freq(txt, verbose);
      datafile->freqMode = FREQMODE_UNIFORM;
      if(datafile->freqlabel_list != NULL) {
      	free(datafile->freqlabel_list);
	datafile->freqlabel_list = NULL;
      }
      //      datafile->freq_list = malloc(2*sizeof(double));
      //      if(datafile->freq_list == NULL) {
      //	printerror(verbose.debug, "ERROR readEPNsubHeader: Memory allocation error.");
      //	return 0;
      //      }
      set_centre_frequency(datafile, freq, verbose);
      if(set_bandwidth(datafile, bw, verbose) == 0) {
	printerror(verbose.debug, "ERROR readEPNsubHeader: Bandwidth changing failed.");
	return 0;
      }
    }

    ret = fread(txt, 1, 17, datafile->fptr);  //tstart
    txt[17] = 0;
    if(ret != 17) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    if(version == 60) {   // Extension
      ret = fread(txt, 1, 23, datafile->fptr);
      txt[23] = 0;
      if(ret != 23) {
	printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
	return 0;
      }
    }else if(version == 61) {   // Extension
      ret = fread(txt, 1, 7, datafile->fptr);
      txt[7] = 0;
      if(ret != 7) {
	printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
	return 0;
      }
    }else if(version >= 62) {
      ret = fread(txt, 1, 1, datafile->fptr); // plflag
      txt[1] = 0;
      if(ret != 1) {
	printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
	return 0;
      }
      ret = fread(txt, 1, 6, datafile->fptr); //plval
      txt[6] = 0;
      if(ret != 6) {
	printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
	return 0;
      }
    }

    ret = fread(txt, 1, 12, datafile->fptr);  // scale
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%f", scale);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 12, datafile->fptr);  // offset
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = sscanf(txt, "%f", offset);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 12, datafile->fptr);  // rms
    txt[12] = 0;
    if(ret != 12) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 16, datafile->fptr);  // papp
    txt[16] = 0;
    if(ret != 16) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
    ret = fread(txt, 1, 28, datafile->fptr);  // Expansion
    txt[28] = 0;
    if(ret != 28) {
      printerror(verbose.debug, "ERROR readEPNsubHeader: This file is not an EPN file.");
      return 0;
    }
  return 1;
}

// Data written in format long header format 6.0, but short header format 6.3 to make it compatible with PSRCHIVE. I didn't find an official description of the formats beyond 5.0
int writeEPNfile(datafile_definition datafile, float *data, verbose_definition verbose)
{
  float dmin, dmax, scale, offset;
  long n, p, i, f;
  int maxint;
  unsigned int *iprofile;

  maxint  = (1 << 16);

  if(datafile.NrFreqChan != 1) {
    printerror(verbose.debug, "ERROR writeEPNfile: Only one freq channel is supported in this data format");
    return 0;
  }
  iprofile  = (unsigned int *)calloc(datafile.NrBins, sizeof(unsigned int));
  if(iprofile == NULL) {
    printerror(verbose.debug, "ERROR writeEPNfile: Cannot allocate memory.");
    return 0;
  }
 
  f = 0;  // Only one freq channel present
  for(n = 0; n < datafile.NrSubints; n++) {
    if(verbose.verbose && verbose.nocounters == 0) {
      printf("writeEPNfile: pulse %ld/%ld\r", n+1, datafile.NrSubints);
    }
    //    for(p = 0; p < datafile.NrPols; p++) {
    //      memcpy(&datablock[p*datafile.NrBins], &data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))], sizeof(float)*datafile.NrBins);
    //    }
    //    maxmin1(datablock, datafile.NrPols*datafile.NrBins, &dmax, &dmin);
    for(p = 0; p < datafile.NrPols; p++) {
      for(i = 0; i < datafile.NrPols*datafile.NrBins; i++) {
	if((i == 0 && p == 0) || dmax < data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i]) {
	  dmax = data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i];
	}
	if((i == 0 && p == 0) || dmin > data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i]) {
	  dmin = data[datafile.NrBins*(p+datafile.NrPols*(f+n*datafile.NrFreqChan))+i];
	}
      }
    }
    if(dmax == dmin)
      dmax = dmin+1;
    offset = dmin;
    scale = (float)(maxint - 1)/(dmax - dmin);
    //   
    for(i = 0; i < datafile.NrBins; i++)
      iprofile[i] = (int)((data[datafile.NrBins*(0+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] - offset) * scale);
    write_epn_longheader(datafile, 600, verbose);
    write_epn_shortheader(datafile, 630, "I       ", 1, verbose);
    write_epn_data(datafile, scale, offset, 0.0, iprofile, verbose);
    if(datafile.NrPols == 4) {
      for(i = 0; i < datafile.NrBins; i++) 
	iprofile[i] = (int)((data[datafile.NrBins*(1+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] - offset) * scale);
      write_epn_shortheader(datafile, 630, "Q       ", 2, verbose);
      write_epn_data(datafile, scale, offset, 0.0, iprofile, verbose);

      for(i = 0; i < datafile.NrBins; i++) 
	iprofile[i] = (int)((data[datafile.NrBins*(2+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] - offset) * scale);
      write_epn_shortheader(datafile, 630, "U       ", 3, verbose);
      write_epn_data(datafile, scale, offset, 0.0, iprofile, verbose);

      for(i = 0; i < datafile.NrBins; i++) 
	iprofile[i] = (int)((data[datafile.NrBins*(3+datafile.NrPols*(f+n*datafile.NrFreqChan))+i] - offset) * scale);
      write_epn_shortheader(datafile, 630, "V       ", 4, verbose);
      write_epn_data(datafile, scale, offset, 0.0, iprofile, verbose);
    }
    //    epn.tstart += (1.0E+06 * get_period(datafile, 0, verbose));
  }
  free(iprofile);
  printf("  Done                       \n");
  return 1;
}



  //Returns 0 on error, 1 if sucess.
int readPulseEPNData(datafile_definition *datafile, long pulsenr, int polarization, int freq, int binnr, long nrSamples, float *pulse, verbose_definition verbose)
{
  //  long i, j, k, l, f;
  float scale, offset;
  int n;
  long pos, ret, nrrecords;

  // The data is written as 80 byte lines.
  nrrecords = datafile->NrBins/20;  // nr of full lines
  if(datafile->NrBins > nrrecords*20) // If incomplete line, it is completed with junk
    nrrecords += 1;

  if(datafile->NrFreqChan > 1) {
    printerror(verbose.debug, "ERROR readPulseEPNData: Only one freq channel is supported in this data format");
    return 0;
  }

  //  Assume the individual pulses are in order
  pos = 480+(160+nrrecords*80)*(datafile->NrPols*datafile->NrFreqChan); // Size of one subint = header size + subheadersizes + data sizes
  pos *= pulsenr; // Start of block = subint
  pos += 480; // Skip header of subint
  pos += (160+nrrecords*80)*polarization; // Skip earlier polarizations. Should have something with NrFreqChan/freq here
  fseek(datafile->fptr, pos, SEEK_SET);

  //Returns 0 on error, 1 if sucess.
  if(readEPNsubHeader(datafile, &scale, &offset, verbose) == 0) {
    printerror(verbose.debug, "ERROR readPulseEPNData: Reading subheader failed.");
    return 0;
  }

  // Skip the requested number of bins
  fseek(datafile->fptr, 4*binnr, SEEK_CUR);
  unsigned int iprofile; 
  //  iprofile = malloc(4*nrSamples);
  //  if(iprofile == NULL) {
  //    printerror(verbose.debug, "ERROR readPulseEPNData: Memory allocation error.");
  //    return 0;
  //  }
  //  ret = fread(iprofile, 4, nrSamples, datafile->fptr);
  //  if(ret != nrSamples) {
  //    printerror(verbose.debug, "ERROR readPulseEPNData: Read error.");
  //    return 0;
  //  }
  for(n = 0; n < nrSamples; n++) {
    ret = fscanf(datafile->fptr,"%04X", &iprofile);
    if(ret != 1) {
      printerror(verbose.debug, "ERROR readPulseEPNData: Reading data failed.");
      return 0;
    }
    //    pulse[n] = (float)iprofile[n] / scale + offset;
    pulse[n] = (float)iprofile / scale + offset;
    //    printf("XXXX %f\n", pulse[n]);
  }

  //  printf("XXXX scale=%f offset=%f\n", scale, offset);
  
  //  free(iprofile);
  return 1;
}


int readEPNfile(datafile_definition *datafile, float *data, verbose_definition verbose, long request_only_one_pulse)
{
  long n, f, p;
  if(verbose.verbose) {
    printf("Start reading EPN file\n");
  }
  for(n = 0; n < datafile->NrSubints; n++) {
    for(p = 0; p < datafile->NrPols; p++) {
      for(f = 0; f < datafile->NrFreqChan; f++) {
        if(verbose.verbose && verbose.nocounters == 0) 
          printf("  Progress reading EPN file (%.1f%%)\r", 100.0*(n+(f+p*datafile->NrFreqChan)*datafile->NrSubints)/(float)(datafile->NrSubints*datafile->NrFreqChan*datafile->NrPols));
	if(readPulseEPNData(datafile, n, p, f, 0, datafile->NrBins, &data[datafile->NrBins*(p+datafile->NrPols*(f+n*datafile->NrFreqChan))], verbose) == 0) {
	  printerror(verbose.debug, "ERROR readEPNfile: Read error.");
	  return 0;
	}
      }
    }
  }
  if(verbose.verbose) printf("  Reading is done.                           \n");
  return 1;
}


//START REGION DEVELOP
