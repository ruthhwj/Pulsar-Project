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
#include "pcenter.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);


void PlotProfile(int NrBins, float *Ipulse, char *xlabel, char *ylabel, char *title, int color, int clearPage);
void ShiftProfile(int shift, int NrBins, float *Iprofile, float *outputProfile);
void circular_smoothing(float *serie, float *new_serie, long nrbins, int smooth);

int main(int argc, char **argv)
{
  int j, index, read_wholefile;
  datafile_definition fin, fout;
  long p, n, f;

  char output_name[1000], PlotDevice[100], c;
  int b, circularShift, noinput; //, fouriershift, deviceID; 
  int shift, bin1, bin2, autocentre, halfphase, smooth, shiftPhaseSet;  //, smooth2
  float *Ismooth, *Ipulse, *Iprofile, *shiftedProfile, x, y, I, shiftPhase;
  psrsalsaApplication application;

  initApplication(&application, "pcenter", "[options] inputfile");
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.oformat = 1;
  application.switch_formatlist = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_showrevision = 1;
  application.switch_nocounters = 1;

  strcpy(PlotDevice, "?");
  strcpy(output_name, "centred.gg");
  circularShift = 0;
  noinput = 0;
  shift = 0;
  autocentre = 0;
  halfphase = 0;
  smooth = 0;
  //  smooth2 = 0;
  shiftPhaseSet = 0;
  shiftPhase = 0;
  //  fouriershift = 0;
  read_wholefile = 1;
  bin1 = bin2 = 0;

  if(argc < 2) {
    printf("Program to add center the pulse in pulse phase (without using a timing solution). Usage:\n\n");
    printApplicationHelp(&application);
    printf("\n\nOptional options:\n");
    printf("-output      Output name. Default is \"%s\"\n", output_name);
    printf("-c           Turn on circular shifting (so that no pulse is lost).\n");
    printf("-peak        Put peak at phase zero.\n");
    printf("-peakhalf    Put peak at phase 0.5.\n");
    printf("-smoothprof  Smooth profile with this number of bins.\n");
    printf("-smooth      Smooth output with this number of bins.\n");
    printf("-nrbins      No graphical input, shifting by this number of bins.\n");
    printf("-phase       No graphical input, shifting by this phase.\n");
    printf("-memsave     Try to use less memory. Not all conversions will work with\n               or without this switch.\n");
    return 0;
  }else {
    for(j = 1; j < argc; j++) {
      index = j;
      if(processCommandLine(&application, argc, argv, &index)) {
	j = index;
      }else if(strcmp(argv[j], "-output") == 0) {
	strcpy(output_name,argv[j+1]);
        j++;
      }else if(strcmp(argv[j], "-memsave") == 0) {
	read_wholefile = 0;
      }else if(strcmp(argv[j], "-peak") == 0) {
	autocentre = 1;
	noinput = 1;
      }else if(strcmp(argv[j], "-peakhalf") == 0) {
	autocentre = 1;
	noinput = 1;
	halfphase = 1;
      }else if(strcmp(argv[j], "-c") == 0) {
	circularShift = 1;
      }else if(strcmp(argv[j], "-nrbins") == 0) {
	noinput = 1;
	shift = atoi(argv[j+1]);
	j++;
      }else if(strcmp(argv[j], "-phase") == 0) {
	noinput = 1;
	shiftPhase = atof(argv[j+1]);
	shiftPhaseSet = 1;
	j++;
      }else if(strcmp(argv[j], "-smoothprof") == 0) {
	smooth = atoi(argv[j+1]);
	j++;
      }else if(strcmp(argv[j], "-smooth") == 0) {
        printerror(application.verbose_state.debug, "Option (no longer) implemented?\n");
	return 0;
	//	smooth2 = atoi(argv[j+1]);
	j++;
      }else if(j < argc-1){
        printerror(application.verbose_state.debug, "Unknown option: %s\n", argv[j]);
	return 0;
      }
    }
  }

  printf("Circular shifting is ");
  if(circularShift)
    printf("on\n");
  else
    printf("off\n");

  /* Clear struct */
  cleanPSRData(&fin, application.verbose_state);

  if(application.iformat <= 0) {
    application.iformat = guessPSRData_format(argv[argc-1], 0, application.verbose_state);
    if(application.iformat == -2 || application.iformat == -3)
      return 0;
  }
  if(isValidPSRDATA_format(application.iformat) == 0) {
      printerror(application.verbose_state.debug, "ERROR pconv: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", argv[argc-1]);
    return 0;
  }

  /* Open input file and read header */
  if(!openPSRData(&fin, argv[argc-1], application.iformat, 0, read_wholefile, 0, application.verbose_state))
    return 0;
  if(read_wholefile == 0) {
    if(!readHeaderPSRData(&fin, 0, 0, application.verbose_state))
      return 0;
  }

  /* Search commandline for header parameters to overwrite header
     parameters. */
  if(PSRDataHeader_parse_commandline(&fin, argc, argv, application.verbose_state) == 0)
    return 0;

  if(noinput == 0) {
    //    deviceID = 
    ppgopen(PlotDevice);
    ppgask(0);
    ppgslw(1);
  }

  if(shiftPhaseSet) {
    shift = shiftPhase*fin.NrBins;
  }
  if(shift >= fin.NrBins)
    shift -= fin.NrBins;
  if(shift < 0)
    shift += fin.NrBins;

  /* Allocate memory */
  Ipulse = (float *)malloc(fin.NrPols*fin.NrBins*sizeof(float));
  Iprofile = (float *)malloc(fin.NrPols*fin.NrBins*sizeof(float));
  Ismooth = (float *)malloc(fin.NrPols*fin.NrBins*sizeof(float));
  shiftedProfile = (float *)malloc(fin.NrPols*fin.NrBins*sizeof(float));
  if(Ipulse == NULL || Iprofile == NULL || Ismooth == NULL || shiftedProfile  == NULL) {
    printerror(application.verbose_state.debug, "Cannot allocate memory.\n");
    return 0;
  }

  memset(Iprofile, 0, fin.NrBins*sizeof(float));
  for(p = 0; p < fin.NrPols; p++) {
    for(n = 0; n < fin.NrSubints; n++) {
      for(f = 0; f < fin.NrFreqChan; f++) {
	if(readPulsePSRData(&fin, n, p, f, 0, fin.NrBins, Ipulse, application.verbose_state) != 1) {
	  printerror(application.verbose_state.debug, "pcenter: Read error.\n");
	  return 0;
	}
	/* Store first and last pulse for circular shifting */
	/*	if(n == 0)
	  memcpy(&Ifirst[p*fin.NrBins], Ipulse, fin.NrBins*sizeof(float));
	if(n == hdr.redn.NTimeInts-1)
	memcpy(&Ilast[p*fin.NrBins], Ipulse, fin.NrBins*sizeof(float));*/
	for(b = 0; b < fin.NrBins; b++) {
	  Iprofile[b+p*fin.NrBins] += Ipulse[b];
	}
      }
    }
    if(application.verbose_state.verbose && application.verbose_state.nocounters == 0) {
      printf("pcenter: Reading file: %.1f%%     \r", (100.0*(n+1))/(float)(fin.NrSubints-2));
      fflush(stdout);
    }
  }
  if(application.verbose_state.verbose && application.verbose_state.nocounters == 0) {
    printf("\nReading is done.                        \n");
  }

  if(smooth) {
    circular_smoothing(Iprofile, Ismooth, fin.NrBins, smooth);
    memcpy(Iprofile, Ismooth, fin.NrBins*sizeof(float));
  }
  /*  if(smooth2) {
    for(j = 0; j < fin.NrPols; j++) {
      circular_smoothing(Ifirst+j*fin.NrBins, Ismooth, fin.NrBins, smooth2);
      memcpy(Ifirst+j*fin.NrBins, Ismooth, fin.NrBins*sizeof(float));
      circular_smoothing(Ilast+j*fin.NrBins, Ismooth, fin.NrBins, smooth2);
      memcpy(Ilast+j*fin.NrBins, Ismooth, fin.NrBins*sizeof(float));
    }
    }*/

  if(autocentre) {
    I = Iprofile[0];
    shift = 0;
    for(b = 0; b < fin.NrBins; b++) {
      if(Iprofile[b] > I) {
	I = Iprofile[b];
	shift = b;
      }
    }
    shift = -shift;
    if(halfphase)
      shift += fin.NrBins/2;
    if(shift >= fin.NrBins)
      shift -= fin.NrBins;
    if(shift < 0)
      shift += fin.NrBins;
  }

  if(application.verbose_state.verbose) printf("Shift is %d\n", shift);

  /* Plot profile of first file */
  if(noinput == 0)
    PlotProfile(fin.NrBins, Iprofile, "Bin number", "I", "Pulse profile from first pulse stack", 1, 1);      

  /* Shift profile graphically*/
  do {                              
    if(noinput == 0) {
      /* Plot profile of second file */
      ShiftProfile(shift, fin.NrBins, Iprofile, shiftedProfile);
      PlotProfile(fin.NrBins, Iprofile, "Bin number", "Intensity", "Click to shift profile, press s to stop", 1, 1);      
      PlotProfile(fin.NrBins, shiftedProfile, "", "", "", 2, 0);      

      /* Ask user what shift should be, or take it from command line */
      j = 0;
      do {
	if(j == 0) 
	  ppgband(0, 0, 0.0, 0.0, &x, &y, &c);
	else
	  ppgband(4, 0, bin1, 0.0, &x, &y, &c);  
	if(c == 65) {    /* left mouse button */
	  if(j == 0)
	    bin1 = x;
	  else
	    bin2 = x;
	  j++;
	}else if(c == 115 || c ==83 ) {    /* s or S */
	  j = 10;
	}
      }while(j < 2);
      if(j < 10) {
	shift += bin2 - bin1;
	if(shift >= fin.NrBins)
	  shift -= fin.NrBins;
	if(shift < 0)
	  shift += fin.NrBins;
      }
    }else {
      j = 10;
    }
  }while(j < 10);


  if(continuous_shift(fin, &fout, shift, circularShift, output_name, application.oformat, argc, argv, application.verbose_state, application.verbose_state.verbose && (application.verbose_state.nocounters == 0)) == 0)
    return 0;
    
  closePSRData(&fin, 0, application.verbose_state);
  closePSRData(&fout, 0, application.verbose_state);

  free(Ipulse);
  free(Ismooth);
  free(Iprofile);
  free(shiftedProfile);

  if(noinput == 0)
    ppgend();
  terminateApplication(&application);
  return 0;
}
 

void PlotProfile(int NrBins, float *Ipulse, char *xlabel, char *ylabel, char *title, int color, int clearPage)
{
  long j;
  float ymin, ymax;
  
  ymin = ymax = Ipulse[0];
  for(j = 1; j < NrBins; j++) {
    if(Ipulse[j] > ymax)
      ymax = Ipulse[j];
    if(Ipulse[j] < ymin)
      ymin = Ipulse[j];
  }

  if(clearPage) {
    ppgpage();
    ppgsci(1);
    ppgsvp(0.1, 0.9, 0.1, 0.9);
    ppgswin(0,NrBins-1,-0.1,1.1);
    ppgbox("bcnsti",0.0,0,"bcnti",0.0,0);
    ppglab(xlabel, ylabel, title);
  }
  ppgsci(color);
  ppgmove(0, Ipulse[0]/ymax);
  for(j = 1; j < NrBins; j++) {
    ppgsci(color);
    ppgdraw(j, Ipulse[j]/ymax);
  }
  ppgsci(1);
}


void ShiftProfile(int shift, int NrBins, float *Iprofile, float *outputProfile)
{
  int b, b2;
  for(b = 0; b < NrBins; b++) {
    b2 = b+shift;
    if(b2 < 0)
      b2 += NrBins;
    if(b2 >= NrBins)
      b2 -= NrBins;
    outputProfile[b2] = Iprofile[b];
  }
}

void circular_smoothing(float *serie, float *new_serie, long nrbins, int smooth)
{
  long i, j;
  for(i = 0; i < nrbins; i++) {
    new_serie[i] = 0;
    for(j = i-smooth; j <= i+smooth; j++) {
      if(j < 0)
	new_serie[i] += serie[j+nrbins];
      else if(j >= nrbins)
	new_serie[i] += serie[j-nrbins];
      else
	new_serie[i] += serie[j];
    }
  }
}




