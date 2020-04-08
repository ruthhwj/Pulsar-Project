#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "psrsalsa.h"

void SHOWREVISIONINFO_prog() {
#include "prebin.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);

psrsalsaApplication application;

int main(int argc, char **argv)
{
  int read_wholefile, index, j, dummyi;
  long rebin, n, f, p;
  float *Ipulse, *Ipulse2;
  char output_name[1000];
  datafile_definition fin, fout, clone;

  initApplication(&application, "prebin", "[options] inputfile");
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.oformat = FITS_format;
  application.switch_formatlist = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_showrevision = 1;
  application.switch_nocounters = 1;

  rebin = 0;
  read_wholefile = 1;
  strcpy(output_name, "rebin.gg");

  if(argc < 2) {
    printApplicationHelp(&application);
    fprintf(stdout, "\nOther options:\n");
    fprintf(stdout, "-rebin       Rebin to this number of bins.\n");
    fprintf(stdout, "-output      Specify output filename, default is output.dat\n");
    fprintf(stdout, "-memsave     Try to use less memory. Not all conversions will work with\n             or without this switch.\n");
    return 0;
  }else {
    for(dummyi = 1; dummyi < argc; dummyi++) {
      index = dummyi;
      if(processCommandLine(&application, argc, argv, &index)) {
	dummyi = index;
      }else if(strcmp(argv[dummyi], "-rebin") == 0) {
	j = sscanf(argv[dummyi+1], "%ld", &rebin);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "prebin: Error parsing -rebin option\n");
	  return 0;
	}
        dummyi++;
      }else if(strcmp(argv[dummyi], "-output") == 0) {
	strcpy(output_name, argv[++dummyi]);
      }else if(strcmp(argv[dummyi], "-memsave") == 0) {
	read_wholefile = 0;
      }else if(dummyi < argc-1) {
	printerror(application.verbose_state.debug, "prebin: Unknown option (%s).\n", argv[dummyi]);
	return 0;
      }
    }
  }

  /* Clear struct */
  cleanPSRData(&fin, application.verbose_state);
  cleanPSRData(&fout, application.verbose_state);

  if(rebin == 0) {
    printerror(application.verbose_state.debug, "prebin: You should specify the -rebin option.\n\n");
    return 0;
  }

  if(application.iformat <= 0) {
    application.iformat = guessPSRData_format(argv[argc-1], 0, application.verbose_state);
    if(application.iformat == -2 || application.iformat == -3)
      return 0;
  }
  if(isValidPSRDATA_format(application.iformat) == 0) {
    printerror(application.verbose_state.debug, "ERROR prebin: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", argv[argc-1]);
    return 0;
  }

  /* Open input file and read header */
  if(!openPSRData(&fin, argv[argc-1], application.iformat, 0, read_wholefile, 0, application.verbose_state))
    return 0;
  if(read_wholefile == 0)
    if(!readHeaderPSRData(&fin, 0, 0, application.verbose_state))
      return 0;

  /* Search commandline for header parameters to overwrite header
     parameters. */
  if(PSRDataHeader_parse_commandline(&fin, argc, argv, application.verbose_state) == 0)
    return 0;

  /* Copy header parameters to output header */
  copy_params_PSRData(fin, &fout, application.verbose_state);
  fout.NrBins = rebin;
  fout.fixedtsamp = get_tsamp(fin, 0, application.verbose_state)*fin.NrBins/(double)rebin;
  fout.tsampMode = TSAMPMODE_FIXEDTSAMP;


  if(fout.NrBins > fin.NrBins) {
    printerror(application.verbose_state.debug, "prebin: Cannot rebin to a larger number of bins (%ld vs %ld).\n", fout.NrBins, fin.NrBins);
    return 0;
  }

  /* Open output file and write data */
  if(!openPSRData(&fout, output_name, application.oformat, 1, 0, 0, application.verbose_state))
    return 0;
  if(!writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state))
    return 0;
  //  appendHistoryLine(fout, argc, argv, application.verbose_state);

  /* Allocate memory */

  Ipulse = (float *)malloc(fin.NrBins*sizeof(float));
  Ipulse2 = (float *)malloc(fout.NrBins*sizeof(float));
  if(Ipulse == NULL || Ipulse2 == NULL) {
    printerror(application.verbose_state.debug, "prebin: Cannot allocate memory.\n");
    return 0;
  }

  if(read_wholefile == 1) {
    if(!preprocess_rebin(fin, &clone, rebin, application.verbose_state))
      return 0;
    swap_orig_clone(&fin, &clone, application.verbose_state); 
    /* Write out all the data, it's as simple as that. */
    if(writePSRData(&fout, fin.data, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "prebin: Cannot write data\n");
      return 0;
    }
  }else {
    if(fin.NrBins % fout.NrBins != 0) {
      fflush(stdout);
      printwarning(application.verbose_state.debug, "WARNING prebin: Rebinning from %ld to %ld bins implies that separate bins are not entirely independent.", fin.NrBins, fout.NrBins);
    }
    /* Loop over pulses and frequency channels */
    for(p = 0; p < fin.NrPols; p++) {
      for(f = 0; f < fin.NrFreqChan; f++) {
	for(n = 0; n < fin.NrSubints; n++) {
	  long i;
	  i = 0;  // For some reason in code below there was an i (which was uninitialised). Assume it should be zero in code below.
	  if(readPulsePSRData(&fin, n, p, f, i, fin.NrBins, Ipulse, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "prebin: Cannot read data\n");
	    return 0;
	  }
	  if(rebinPulse(Ipulse, fin.NrBins, Ipulse2, fout.NrBins, 1, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "prebin: Cannot rebin data.\n");
	    return 0;
	  }
	  if(writePulsePSRData(&fout, n, p, f, i, fout.NrBins, Ipulse2, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "prebin: Cannot write data\n");
	    return 0;
	  }
	}
	printf("Processed %.1f%%     \r", (100.0*(n+1))/(float)(fin.NrSubints));
	fflush(stdout);
      }
    }
    printf("\n");
  }

  free(Ipulse);
  free(Ipulse2);
  /* Free resources and close files */
  closePSRData(&fin, 0, application.verbose_state);
  closePSRData(&fout, 0, application.verbose_state);
  terminateApplication(&application);
  return 0;
}
