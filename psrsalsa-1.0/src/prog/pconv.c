//START REGION RELEASE
/* Make sure that programs can handle large files */
#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <stdio.h>
#include <string.h>
#include <psrsalsa.h>
//START REGION DEVELOP
//#include "libpuma.h"
//#include "fpuma.h"
//START REGION RELEASE

int pumawrite(void *prptr, int size, int nelem,FILE *out);

//START REGION DEVELOP
void SHOWREVISIONINFO_prog() {
#include "pconv.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);
//START REGION RELEASE

//START REGION DEVELOP
/* 
   This is an internal epn function, which is used for low level
   writing when the program runs in single pulse per time mode. Using
   these functions is not the standard way to handle pulsar data using
   psrsalsalib.
*/
//void write_epn_psrsalsa(FILE *fptr, EPNpsrsalsa epn, int IOSwitch, int version);
//int fill_headerEPNstruct(datafile_definition datafile, headerEPN *epn, verbose_definition verbose);
//START REGION RELEASE

void print_help();

int main(int argc, char **argv)
{
  datafile_definition fin, fout;
  int read_wholefile, indx;
  long i, n, p, f, n1, n2;
  float *pulseData, sample;
  char outputname[MaxFilenameLength], *dummy_ptr;
  psrsalsaApplication application;

//START REGION DEVELOP
  int sctd;
  long k;
  int searchmode;
//START REGION RELEASE

  initApplication(&application, "pconv", "[options] inputfile(s)");
//START REGION DEVELOP
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
//START REGION RELEASE
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.switch_ext = 1;
  application.switch_output = 1;
  application.switch_fchan = 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_formatlist = 1;
  application.switch_tscr = 1;
  application.switch_tscr_complete = 1;
  application.switch_TSCR = 1;
  application.switch_FSCR = 1;
  application.switch_fscr = 1;
  application.switch_rebin = 1;
  application.switch_scale = 1;
  application.switch_polselect = 1;
  application.switch_nskip = 1;
  application.switch_nread = 1;
  application.switch_noweights = 1;
  application.switch_useweights = 1;
  application.switch_uniformweights = 1;
  application.switch_norm = 1;
  application.switch_debase = 1;
  application.switch_onpulse = 1;
  application.switch_onpulsef = 1;
  application.switch_nocounters = 1;
  application.switch_stokes = 1;
  application.switch_deparang = 1;
  application.switch_history_cmd_only = 1;
//START REGION DEVELOP
  application.switch_testsignal = 1;
  application.switch_gate = 1;
  application.switch_useweightedfreq = 1;
  application.switch_absweights = 1;
  application.switch_normRMS = 1;
  application.switch_showrevision = 1;
//START REGION RELEASE

 /*
  application.switch_ = 1;
  */
  read_wholefile = 1;
//START REGION DEVELOP
  sctd = 0;
  searchmode = 0;
//START REGION RELEASE



  if(argc <= 1) {
    printApplicationHelp(&application);
    print_help();
    terminateApplication(&application);
    return 0;
  }
  for(indx = 1; indx < argc; indx++) {
    if(processCommandLine(&application, argc, argv, &indx)) {
    }else if(strcmp(argv[indx], "-memsave") == 0) {
      read_wholefile = 0;
//START REGION DEVELOP
    }else if(strcmp(argv[indx], "-sctd") == 0) {
      sctd = 1;
    }else if(strcmp(argv[indx], "-ingrid") == 0) {
      add_IngridFlavour_SigprocAsciiFormat(1, argv[++indx]);
    }else if(strcmp(argv[indx], "-searchmode") == 0) {
      searchmode = 1;
//START REGION RELEASE
    }else {
      /* If the option is not recognized, assume it is a filename */
      if(argv[indx][0] == '-') {
	printerror(application.verbose_state.debug, "ERROR pconv: Unknown option (%s).\n\nRun pconv without command line arguments to show help", argv[indx]);
	terminateApplication(&application);
	return 0;
      }else {
	if(applicationAddFilename(indx, application.verbose_state) == 0)
	  return 0;
      }
      /*
    }else if(indx < argc-1) {
      printerror(application.verbose_state.debug, "pconv: Unknown option (%s).", argv[indx]);
      return 0;*/
    }
  }

  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }

  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pconv: No files specified");
    return 0;
  }

  if(isValidPSRDATA_format(application.oformat) == 0) {
    printerror(application.verbose_state.debug, "ERROR pconv: Please specify a valid output format with the -oformat option.");
    terminateApplication(&application);
    return 0;
  }

  /* Clear struct */
  //  cleanPSRData(&fin, application.verbose_state);
  cleanPSRData(&fout, application.verbose_state);

  /* If there are files specified with the -i option or if there are
     files at the end of the command line, then process them. */
  while((dummy_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {


    if(application.iformat <= 0) {
      application.iformat = guessPSRData_format(dummy_ptr, 0, application.verbose_state);
      if(application.iformat == -2 || application.iformat == -3) {
	closePSRData(&fout, 0, application.verbose_state);
	terminateApplication(&application);
	return 0;
      }
    }
    if(isValidPSRDATA_format(application.iformat) == 0) {
      printerror(application.verbose_state.debug, "ERROR pconv: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", dummy_ptr);
      closePSRData(&fout, 0, application.verbose_state);
      terminateApplication(&application);
      return 0;
    }

  /* Make sure that some conversions will work for which -memsave is or isn't working */
  /*
  if((application.iformat == EPN_format) && read_wholefile == 0) {
    printwarning(application.verbose_state.debug, "WARNING pconv: Ignoring -memsave for this format.");
    read_wholefile = 1;
  }
  */
//START REGION DEVELOP
  if((application.iformat == AO_ASCII_1_format) && read_wholefile == 0) {
    printwarning(application.verbose_state.debug, "WARNING pconv: Ignoring -memsave for this format.");
    read_wholefile = 1;
  }
  if(application.iformat == PARKESFB_format && read_wholefile == 1) {
    printwarning(application.verbose_state.debug, "WARNING pconv: Enabling -memsave for this format.");
    read_wholefile = 0;
  }
  if((application.iformat == PRESTO_format) && read_wholefile == 1) {
    printwarning(application.verbose_state.debug, "WARNING pconv: Enabling -memsave for this format.");
    read_wholefile = 0;
  }
//START REGION RELEASE
  //  if((application.iformat == SIGPROC_format) && read_wholefile == 1) {
  //    printwarning(application.verbose_state.debug, "WARNING pconv: Enabling -memsave for this format.");
  //    read_wholefile = 0;
  //  }

  /* Open input file and read header */
//START REGION DEVELOP
  if(application.useweightedfreq)
    psrfits_set_use_weighted_freq(1);
//START REGION RELEASE
  if(!openPSRData(&fin, dummy_ptr, application.iformat, 0, read_wholefile, 0, application.verbose_state))
    return 0;


  if(read_wholefile == 0) {
    if(application.fchan_select != -1) {
      printerror(application.verbose_state.debug, "ERROR pconv: -fchan option doesn't work with -memsave option.");
      return 0;
    }
    if(!readHeaderPSRData(&fin, 0, 0, application.verbose_state))
      return 0;
  }

  /* Search commandline for header parameters to overwrite header
     parameters. */
  if(PSRDataHeader_parse_commandline(&fin, argc, argv, application.verbose_state) == 0)
    return 0;

  /* Convert the specified onpulse regions which were defined as
     fractions on the commandline. This is now possible because we
     know the number of bins. */
  region_frac_to_int(&(application.onpulse), fin.NrBins, 0);
  
  if(read_wholefile != 0) {
    for(i = 1; i < argc; i++) {
      if(strcmp(argv[i], "-header") == 0) {
	printwarning(application.verbose_state.debug, "WARNING pconv: If using the -header option, be aware it applied BEFORE the preprocessing.");
	break;
      }
    }
    if(preprocessApplication(&application, &fin) == 0) {
      return 0;
    }
  }

  /* Copy header parameters to output header */
  copy_params_PSRData(fin, &fout, application.verbose_state);
  /* Make sure the period is set to something */
  double period;
  int ret;
  if(fin.isFolded) {
    ret = get_period(fin, 0, &period, application.verbose_state);
    if(ret == 2) {
      printerror(application.verbose_state.debug, "ERROR pconv (%s): Cannot obtain period", fin.filename);
      return 0;
    }
  }else {
    ret = 1;
    period = -1;
  }
  if(ret == 0 && period >= 0 && period < 0.001) {
    fout.fixedPeriod = 1;
    fout.foldMode = FOLDMODE_FIXEDPERIOD;
  }

  if(getOutputName(&application, dummy_ptr, outputname, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pconv: Changing filename failed");
    return 0;
  }

  n1 = application.nskip;
  if(application.nread > 0)
    n2 = application.nskip + application.nread;
  else
    n2 = fin.NrSubints;

  fout.NrSubints = n2 - n1;


//START REGION DEVELOP
  if(searchmode) {
    fout.isFolded = 0;
    fout.foldMode = FOLDMODE_UNKNOWN;
    fout.gentype = GENTYPE_SEARCHMODE;
  }
//START REGION RELEASE

  /* Open output file and write data */
  if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
    return 0;
  if(!writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state))
    return 0;
  //  appendHistoryLine(fout, argc, argv, application.verbose_state);

  /* If this switch is set the data is read/written on a pulse to
     pulse basis. Not all conversions are implemented, mostly because
     the order of reading and writing is really important. For most
     formats it is much more easy to read in everything and then write
     it out. */
  if(read_wholefile == 1) {
    /* Write out all the data, it's as simple as that. */
    if(writePSRData(&fout, fin.data, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR pconv: Cannot write data");
      return 0;
    }
  }else {
//START REGION DEVELOP
    if(application.oformat == AO_ASCII_1_format) {
      float SigI, AverageI, ScaleVal;
      SigI = 0;
      AverageI = 0;
      ScaleVal = 1;
      for(n = n1; n < n2; n++) {
	if(application.verbose_state.verbose) printf("pulse %ld/%ld\r", n+1-n1, fout.NrSubints);
	for(i = 0; i < fin.NrBins; i++) {
	  for(p = 0; p < fin.NrPols; p++) {
	    if(readPulsePSRData(&fin, n, p, 0, i, 1, &sample, application.verbose_state) == 0) {
	      printerror(application.verbose_state.debug, "ERROR pconv: Cannot read data");
	      return 0;
	    }
	    if(fin.NrPols == 1) {
	      fprintf(fout.fptr, "%f  %f  %f  %f ", sample/500000.0, 0.0,0.0,0.0);
	    }else {
	      if(p == 3)
		fprintf(fout.fptr, "%f ", sample);
	      else
		fprintf(fout.fptr, "%f  ", sample);
	    }
	  }
	  if(((i+1)%2) == 0) fprintf(fout.fptr, "\n");
	}
	fprintf(fout.fptr, "%f  %f  %ld\n", SigI, AverageI/ScaleVal, n+1);
      }
    }else if(application.oformat == PRESTO_format) {
      if(fin.NrPols > 1) {
	printwarning(application.verbose_state.debug, "WARNING pconv: only stokes I is converted.");
      }
      if(fin.NrSubints > 1) {
	printwarning(application.verbose_state.debug, "WARNING pconv: only first pulse is converted.");
      }
      printwarning(application.verbose_state.debug, "WARNING pconv: Assumed that data is barycentered.");
      for(i = 0; i < fin.NrBins; i++) {
	if(readPulsePSRData(&fin, 0, 0, 0, i, 1, &sample, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pconv: Cannot read data");
	  return 0;
	}
	fwrite(&sample, sizeof(float), 1, fout.fptr);
      }
      /*    }else if(application.oformat == EPN_format) {
      float dmin, dmax;
      int maxint, IOSwitch;
      EPNpsrsalsa epn;

      //      writeEPNHeader(fout, &epn, application.verbose_state);
      fill_headerEPNstruct(fout, &epn, application.verbose_state);

      maxint  = (1 << 16);
      pulseData = (float *)malloc(fin.NrBins*sizeof(float)*fin.NrPols);
      epn.iprofile  = (unsigned int *)   calloc(4*fin.NrBins, sizeof(unsigned int));
      if(pulseData == NULL || epn.iprofile == NULL) {
	printerror(application.verbose_state.debug, "ERROR pconv: Cannot allocate memory.");
	return 0;
      }
      for(n = n1; n < n2; n++) {
	if(application.verbose_state.verbose) printf("pulse %ld/%ld\r", n+1-n1, fout.NrSubints);
	for(p = 0; p < fin.NrPols; p++) {
	  if(readPulsePSRData(&fin, n, p, 0, 0, fin.NrBins, &pulseData[p*fin.NrBins], application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pconv: Cannot read data");
	    return 0;
	  }
	}
	maxmin1(pulseData, fin.NrPols*fin.NrBins, &dmax, &dmin);
	//	epn.rms    = SigI[i]; 
      // Let's see if we get away with this 
	epn.rms    = 0;
	epn.offset = dmin;
	if(dmax == dmin)
	  dmax = dmin+1;
	epn.scale  = (float) (maxint - 1) / (dmax - dmin);
	strcpy(epn.idfield, "I       ");             
	epn.nband   = p+1;
	epn.navg  = 1;
	IOSwitch = 0;
	for(i = 0; i < fin.NrBins; i++) 
	  epn.iprofile[i] = (int)((pulseData[i] - epn.offset) * epn.scale);
	write_epn_psrsalsa(fout.fptr, epn, IOSwitch, 600);
	if(fin.NrPols == 4) {
	  IOSwitch = 1;
	  strcpy(epn.idfield, "Q       ");             epn.nband   = i+1;
	  for(i = 0; i < fin.NrBins; i++) 
	    epn.iprofile[i] = (int)((pulseData[fin.NrBins+i] - epn.offset) * epn.scale);
	  write_epn_psrsalsa(fout.fptr, epn, IOSwitch, 600);
	  strcpy(epn.idfield, "U       ");             epn.nband   = i+1;
	  for(i = 0; i < fin.NrBins; i++) 
	    epn.iprofile[i] = (int)((pulseData[2*fin.NrBins+i] - epn.offset) * epn.scale);
	  write_epn_psrsalsa(fout.fptr, epn, IOSwitch, 600);
	  strcpy(epn.idfield, "V       ");             epn.nband   = i+1;
	  for(i = 0; i < fin.NrBins; i++) 
	    epn.iprofile[i] = (int)((pulseData[3*fin.NrBins+i] - epn.offset) * epn.scale);
	  write_epn_psrsalsa(fout.fptr, epn, IOSwitch, 600);
	}
	//      for (i=0;i<epn.nbins;i++) printf("%04X",epn.iprofile[i]); 
	epn.tstart += (1.0E+06 * get_period(fin, 0, application.verbose_state));
	}*/
    }else if(application.iformat == PRESTO_format && application.oformat == PUMA_format) {
      /* If there are frequency channels present, I made the hacked presto format having frequency channel in the inner loop, while puma format has it as an outer loop. */
      if(fin.NrFreqChan == 1) {
	for(i = 0; i < fin.NrBins; i++) {
	  fread(&sample, sizeof(float), 1, fin.fptr);
	  pumawrite(&sample, sizeof(float), 1, fout.fptr);
	}      
      }else {
	float *pulseData2;
#define PrestoPumaNrtimesamplesPerBlock   100
	pulseData = (float *)malloc(fin.NrFreqChan*PrestoPumaNrtimesamplesPerBlock*sizeof(float));
	pulseData2 = (float *)malloc(fin.NrFreqChan*PrestoPumaNrtimesamplesPerBlock*sizeof(float));
	if(pulseData == NULL || pulseData2 == NULL) {
	  printerror(application.verbose_state.debug, "ERROR pconv: Cannot allocate memory.");
	  return 0;
	}
	long j, k, binsavailable, curposinreadblock;
	long long filepos;
	printf("\n");
	k = 0;
	binsavailable = 0;
	curposinreadblock = 0;
	for(i = 0; i < fin.NrBins; i++) {
	  if(binsavailable == 0) {
	    binsavailable = fin.NrBins-i;
	    if(binsavailable > PrestoPumaNrtimesamplesPerBlock)
	      binsavailable = PrestoPumaNrtimesamplesPerBlock;
	    fread(pulseData, sizeof(float), fin.NrFreqChan*binsavailable, fin.fptr);
	    curposinreadblock = 0;
	  }
	  for(j = 0; j < fin.NrFreqChan; j++) {
	    pulseData2[j*PrestoPumaNrtimesamplesPerBlock+k] = pulseData[j+curposinreadblock*fin.NrFreqChan];
	  }
	  k++;
	  curposinreadblock++;
	  binsavailable--;
	  if(k == PrestoPumaNrtimesamplesPerBlock) {
	    filepos = (i-(k-1))*sizeof(float) + 4504;
	    for(j = 0; j < fin.NrFreqChan; j++) {
	      fseeko(fout.fptr, filepos, SEEK_SET);
	      pumawrite(&pulseData2[j*PrestoPumaNrtimesamplesPerBlock], sizeof(float), k, fout.fptr);
	      filepos += fin.NrBins*sizeof(float);                /* Begin next frequency channel */
	    }
	    k = 0;
	  }
	  if(application.verbose_state.nocounters == 0 && i%PrestoPumaNrtimesamplesPerBlock == 0) {
	    printf("\r%.2lf%%       ", 100.0*i/(double)(fin.NrBins));
	  }
	}
	filepos = (i-(k))*sizeof(float) + 4504;
	for(j = 0; j < fin.NrFreqChan; j++) {
	  fseeko(fout.fptr, filepos, SEEK_SET);
	  pumawrite(&pulseData2[j*PrestoPumaNrtimesamplesPerBlock], sizeof(float), k, fout.fptr);
	  filepos += fin.NrBins*sizeof(float);                /* Begin next frequency channel */
	}
	k = 0;
      }
      printf("\n");
   }else
 //START REGION RELEASE
      if(application.iformat == SIGPROC_format && application.oformat == PUMA_format) {
      long j;
      unsigned char *sample_b;
      float *sample_f;
      void *data_ptr;

      if(fin.NrSubints > 1) {
	printerror(application.verbose_state.debug, "ERROR pconv: Non-folded sigproc data can be converted in PuMa format with the -memsave option. This data has %ld subints.", fin.NrSubints);
	return 0;
      }

      if(fin.NrBits == 32) {
	data_ptr = malloc(4*fin.NrFreqChan);
	sample_f = data_ptr;
      }else if(fin.NrBits == 8) {
	data_ptr = malloc(fin.NrFreqChan);
	sample_b = data_ptr;
      }else {
	printerror(application.verbose_state.debug, "ERROR pconv: Only 32-bit or 8-bit sigproc data can be converted in PuMa format. This is %ld bit data.", fin.NrBits);
	return 0;
      }
      if(data_ptr == NULL) {
	printerror(application.verbose_state.debug, "ERROR pconv: Cannot allocate memory.");
	return 0;
      }

      //  for(p = 0; p < datafile.NrPols; p++) {
      //    for(f = 0; f < datafile.NrFreqChan; f++) {
      //      for(n = 0; n < datafile.NrSubints; n++) {

      long long filepos;
      long polarization = 0;
      long pulsenr = 0;
      
      for(i = 0; i < fin.NrBins; i++) {
	for(j = 0; j < fin.NrFreqChan; j++) {
	  if(fin.NrBits == 32) {
	    fread(&sample_f[j], sizeof(float), 1, fin.fptr);
	  }else if(fin.NrBits == 8) {
	    fread(&sample_b[j], sizeof(unsigned char), 1, fin.fptr);
	  }
	}
	filepos = polarization*fin.NrFreqChan*fin.NrSubints*fin.NrBins;  /* Begin polarization block */
	filepos += pulsenr*fin.NrBins;                      /* Begin first pulse */
	filepos += i;
	filepos *= sizeof(float);
	filepos += fout.datastart;
	for(j = 0; j < fin.NrFreqChan; j++) {
	  fseeko(fout.fptr, filepos, SEEK_SET);
	  if(fin.NrBits == 32) {
	    sample = sample_f[j];
	  }else if(fin.NrBits == 8) {
	    sample = sample_b[j];
	  }
	  pumawrite(&sample, sizeof(float), 1, fout.fptr);
	  filepos += fin.NrSubints*fin.NrBins*sizeof(float);               
	}  
      }
      free(data_ptr);
//START REGION DEVELOP
    }else if(application.iformat == PARKESFB_format && application.oformat == PUMA_format) {
      short int sample2, l, l2;
      long nrblocks, nsamps_written, nsamps_read;
      long skipsamples = 1878*1635;
      long readsamples = 100*1635;
      long freqsadd = 1;
      float sample_avrg = 0;
      char sample3;

      skipsamples = 0;
      readsamples = 10000;


      l2 = 0;
      if(sizeof(sample2) != 2) {
	printerror(application.verbose_state.debug, "ERROR pconv: sample2 is %d bytes, and it should be 2 bytes.", (int)sizeof(sample2));
	return 0;
      }
      nrblocks = fin.NrBins*fin.NrFreqChan/8192;
      if(application.verbose_state.verbose) printf("Going to process %ld samples (%ld blocks).\n", fin.NrBins, nrblocks);
      nsamps_written = 0;
      nsamps_read = 0;
      /*    for(j = 0; j < nrblocks; j++) { */
      do {
	if(sctd == 0) {
	  fread(&sample2, 2, 1, fin.fptr);
	  fread(&sample2, 2, 1, fin.fptr);
	}else {
	  fread(&sample3, 1, 1, fin.fptr);
	  fread(&sample3, 1, 1, fin.fptr);
	}
	for(i = 0; i < 8192; i++) {
	  if(sctd == 0) {
	    k = fread(&sample2, 2, 1, fin.fptr);
	  }else {
	    k = fread(&sample3, 1, 1, fin.fptr);
	  }
	  if(k != 1) {
	    printerror(application.verbose_state.debug, "ERROR pconv: Unexpected EOF.");
	    break;
	  }
	  if(sctd == 0) {
	    sample = sample2;
	    pumawrite(&sample, sizeof(float), 1, fout.fptr);
	    nsamps_written++;
	  }else {
	    for(l = 0; l < 8; l++) {
	      if((sample3 & (1 << l)) != 0) {
		sample = 1;
		/*		printf("HOI 1\n");*/
	      }else {
		sample = 0;
		/*		printf("HOI 0\n");*/
	      }
	      sample_avrg += sample;
	      l2++;
	      if(l2 >= freqsadd) {
		sample_avrg /= (float)freqsadd;
		/* Probably wrong way around. binnr should be inner loop. */
		if(nsamps_read/fout.NrFreqChan >= skipsamples && nsamps_read/fout.NrFreqChan < skipsamples+readsamples) {
		  /*		pumawrite(&sample, sizeof(float), 1, fout.fptr); */
		  fwrite(&sample_avrg, sizeof(float), 1, fout.fptr);
		  nsamps_written++;
		}
		l2 = 0;
		sample_avrg = 0;
	      }
	      nsamps_read++;
	    }
	  }
	}
	if(k != 1) 
	  break;
	if(sctd == 0) {
	  fread(&sample2, 2, 1, fin.fptr);
	  fread(&sample2, 2, 1, fin.fptr);
	}else {
	  fread(&sample3, 1, 1, fin.fptr);
	  fread(&sample3, 1, 1, fin.fptr);
	}
      }while(k == 1);
      if(sctd == 0) {
	fout.NrFreqChan = 1;
	fout.NrBins = nsamps_written;
      }else {
	fout.NrFreqChan /= freqsadd;
	fout.NrBins = nsamps_written/fout.NrFreqChan;
      }
      fseeko(fout.fptr, 0, SEEK_SET);
      writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state);
      //      appendHistoryLine(fout, argc, argv, application.verbose_state);
      if(application.verbose_state.verbose) printf("Processed %ld samples.\n", nsamps_written);
 //START REGION RELEASE
    }else if(application.iformat == PSRCHIVE_ASCII_format && application.oformat == PUMA_format) {
      char txt[100];
      for(n = 0; n < fin.NrSubints; n++) {
	if(application.verbose_state.verbose) printf("pulse %ld/%ld\r", n+1-n1, fout.NrSubints);
	for(f = 0; f < fin.NrFreqChan; f++) {
	  for(i = 0; i < fin.NrBins; i++) {
	    fscanf(fin.fptr, "%s", txt);
	    if(i == 0 && f == 0 && n == 0) {
	      if(strcmp(txt, "0") != 0) {
		printerror(application.verbose_state.debug, "ERROR pconv: PSRCHIVE ascii header format has probably been changed.");
		return 0;
	      }
	    }
	    fscanf(fin.fptr, "%s", txt);
	    if(i == 0 && f == 0 && n == 0) {
	      if(strcmp(txt, "0") != 0) {
		printerror(application.verbose_state.debug, "ERROR pconv: PSRCHIVE ascii header format has probably been changed.");
		return 0;
	      }
	    }
	    fscanf(fin.fptr, "%s", txt);
	    if(i == 0 && f == 0 && n == 0) {
	      if(strcmp(txt, "0") != 0) {
		printerror(application.verbose_state.debug, "ERROR pconv: PSRCHIVE ascii header format has probably been changed.");
		return 0;
	      }
	    }
	    for(p = 0; p < fin.NrPols; p++) {
	      fscanf(fin.fptr, "%f", &sample);
	      if(n >= n1 && n < n2) {
		if(writePulsePSRData(&fout, n-n1, p, f, i, 1, &sample, application.verbose_state) == 0) {
		  printerror(application.verbose_state.debug, "ERROR pconv: Cannot write data");
		  return 0;
		}
	      }
	    }
	  }
	}
      }
    }else {    /*if(application.iformat == 7 && application.oformat == 1) {*/
      pulseData = (float *)malloc(fin.NrBins*sizeof(float));
      if(pulseData == NULL) {
	printerror(application.verbose_state.debug, "ERROR pconv: Cannot allocate memory.");
	return 0;
      }
      /* For reading the EPN format it is much quicker if pulse number is in the outer loop */
      for(n = n1; n < n2; n++) {
	if(application.verbose_state.verbose) {
	  printf("pulse %ld/%ld\r", n+1-n1, fout.NrSubints);
	  fflush(stdout);
	}
	/*	for(i = 0; i < fin.NrBins; i++) { */
	for(p = 0; p < fin.NrPols; p++) {
	  for(f = 0; f < fin.NrFreqChan; f++) {
	    if(readPulsePSRData(&fin, n, p, f, 0, fin.NrBins, pulseData, application.verbose_state) == 0) {
	      printerror(application.verbose_state.debug, "ERROR pconv: Cannot read individual pulses from data.");
	      return 0;
	    }
	    if(writePulsePSRData(&fout, n-n1, p, f, 0, fin.NrBins, pulseData, application.verbose_state) == 0) {
	      printerror(application.verbose_state.debug, "ERROR pconv: Cannot write individual pulses.");
	      return 0;
	    }
	  }
	}
      }
      free(pulseData);
    }
  }

  /* Free resources and close files */
  closePSRData(&fin, 0, application.verbose_state);
  closePSRData(&fout, 0, application.verbose_state);
  
  }
  terminateApplication(&application);
  return 0;
}

//START REGION DEVELOP
//START REGION RELEASE

void print_help()
{
  fprintf(stdout, "Other options:\n");
  fprintf(stdout, "  -memsave     Try to use less memory. Not all conversions will work with\n               or without this switch.\n");
//START REGION DEVELOP
  fprintf(stdout, "  -sctd        In combination with -iformat 4: input is sc_td data instead of \n");
  fprintf(stdout, "                 pdmd data (experimental and not working).\n");
  fprintf(stdout, "  -ingrid      Adds the Ingrid flavour to the sigproc ascii format using this \n");
  fprintf(stdout, "                 directory to find the parfile.\n");
  fprintf(stdout, "  -searchmode  Change the input data to be classified as in searchmode.\n");  
//START REGION RELEASE
  fprintf(stdout, "\n");
//START REGION DEVELOP
  fprintf(stdout, "Supported formats for reading are:\n");
  fprintf(stdout, "  FITS file\n");
  fprintf(stdout, "  PSRCHIVE ascii dump file\n");
  fprintf(stdout, "  EPN format\n");
  fprintf(stdout, "  PuMa file\n");
  fprintf(stdout, "  Desh arecibo ascii gated file\n");
  fprintf(stdout, "\nSupported formats for writing are:\n");
  fprintf(stdout, "  PSRCHIVE FITS file\n");
  fprintf(stdout, "  PSRCHIVE ascii dump file\n");
  fprintf(stdout, "  EPN format\n");
  fprintf(stdout, "  PuMa gated file\n");
  fprintf(stdout, "  Desh arecibo ascii gated file\n");
  fprintf(stdout, "\nOther more limited conversions:\n");
  fprintf(stdout, "  PuMa timeserie     -> PRESTO timeserie\n");
  fprintf(stdout, "  Sigproc binary     -> PuMa\n");
  fprintf(stdout, "  Simon's Parkes FB  -> PuMa\n");
  fprintf(stdout, "  PRESTO timeserie   -> PuMa timeserie\n");
  fprintf(stdout, "  Parkes timeserie   -> PuMa timeserie\n");
//START REGION RELEASE
  printf("\n");
  printCitationInfo();
}

//START REGION DEVELOP
