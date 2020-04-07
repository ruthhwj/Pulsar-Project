#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "psrsalsa.h"

void SHOWREVISIONINFO_prog() {
#include "psplit.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);

int main(int argc, char **argv)
{
  //  FILE *fin, *fout;
  char output_name[1000], txt[1000], txt2[1000], *inputname;
  //  int NrBins, NrPol, NrFreqChan;
  // long NrPulses, 
  long BlockSize, BlockNumber, NrBlocksToOutput;
  float *Ipulse, *Ipulse2;
  long n, p, f, i, pulse_in;
  //  long long puma_filepos;
  int IndividualChannelFlag, polsplit;
  int NoPolarizationFlag;
  //  int , Verbose;
  int RunningBlock, BlockSplitting;
  datafile_definition fin, fout;

  psrsalsaApplication application;
  initApplication(&application, "psplit", "[options] inputfiles");
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.switch_formatlist = 1;
  application.switch_nocounters = 1;
  application.switch_headerlist = 1;
  application.switch_header = 1;
  
  BlockSize = -1;
  NrBlocksToOutput = -1;
  IndividualChannelFlag = 0;
  NoPolarizationFlag = 0;
  //  Verbose = 0;
  RunningBlock = 0;
  BlockSplitting = 0;
  polsplit = 0;
	
  if(argc < 2) {
    printf("Program to split data files in separate files. Usage:\n\n");
    printApplicationHelp(&application);
    printf("Optional options:\n");
    printf("-f Write out frequency channels seperately.\n");
    printf("-polsplit  Write out polarization channels seperately.\n");
    printf("-I Only writes out first polarization channel.\n");
    printf("-n BlockSize (default is all pulses).\n");
    printf("-r Use in combination with -n option. Each output block of -n subints is offset by one subint, rather than -n subints.\n");
    printf("-b Number of blocks to output (default is all).\n");
    //    printf("-1 Write out blocks to single file.\n");   // Obsolete -> -tscr
    //    printf("-s Shift pulse stack by this number of bins.\n"); // Can be done with pmod before psplit
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      int dummy_int;
      dummy_int = i;
      if(processCommandLine(&application, argc, argv, &dummy_int)) {
	i = dummy_int;
      }else if(strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "-N") == 0) {
	BlockSize = atol(argv[i+1]);
	i++;
	BlockSplitting = 1;
      }else if(strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "-B") == 0) {
	NrBlocksToOutput = atol(argv[i+1]);
	i++;
	//      }else if(strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "-S") == 0) {
	//	ShiftFlag = 1;
	//	NrShiftBins = atoi(argv[i+1]);
	//	i++;
      }else if(strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "-F") == 0) {
	IndividualChannelFlag = 1;
      }else if(strcmp(argv[i], "-polsplit") == 0) {
	polsplit = 1;
      }else if(strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "-I") == 0) {
	NoPolarizationFlag = 1;
	//      }else if(strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "-V") == 0) {
	//	Verbose = 1;
      }else if(( strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "-R") == 0 ) &&  BlockSplitting == 1) {
	RunningBlock = 1;
      }else {
	/* If the option is not recognized, assume it is a filename */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "ERROR psplit: Unknown option %s. Run psplit without command-line options to get help.", argv[i]);
	  terminateApplication(&application);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }

  if(IndividualChannelFlag != 0 && polsplit == 1) {
    printerror(application.verbose_state.debug, "ERROR psplit: Cannot use -f and -polsplit simultaneously.\n");
    return 0;      
  }

  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  int nrinputfiles;
  nrinputfiles = numberInApplicationFilenameList(&application, argv, application.verbose_state);
  if(nrinputfiles < 1) {
    printerror(application.verbose_state.debug, "ERROR psplit: Need at least one input files");
    return 0;
  }

  //  cleanPSRData(&fin, application.verbose_state);
  while((inputname = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
    if(openPSRData(&fin, inputname, application.iformat, 0, 0, 0, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR psplit: Cannot open %s\n", inputname);
      return 0;
    }
    if(readHeaderPSRData(&fin, 0, 0, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR psplit: Cannot read header of file %s\n", inputname);
      return 0;
    }
    /* Search commandline for header parameters to overwrite header
     parameters. */
    if(PSRDataHeader_parse_commandline(&fin, argc, argv, application.verbose_state) == 0)
      return 0;
    printf("Data contains %ld bins, %ld pulses, %ld polarizations and %ld frequencies.\n", fin.NrBins, fin.NrSubints, fin.NrPols, fin.NrFreqChan);
    if(BlockSize <= 0) {
      BlockSize = fin.NrSubints;
      //      if(ShiftFlag)
      //	BlockSize--;
    }

    // I suspect this was only checked for to avoid having to think about it while making other changes
    if(IndividualChannelFlag != 0 && fin.freqMode != FREQMODE_UNIFORM) {
      if(fin.NrSubints > 1) {
	printerror(application.verbose_state.debug, "ERROR psplit: Frequency channels are not uniformly separated. The data cannot be split in individual frequency channels when there are multiple subintegrations, unless the frequency channels have the same uniform labelling from subint to subint.");
	return 0;
      }
      printwarning(application.verbose_state.debug, "WARNING psplit: Frequency channels are not uniformly separated. Be aware that the centre frequencies of the splitted data might not be set differently from what you expect.\n");
    }

    cleanPSRData(&fout, application.verbose_state);
    copy_params_PSRData(fin, &fout, application.verbose_state);

    if(application.oformat != -1) {
      fout.format = application.oformat;
      //      fout.data = NULL; // Important to do, if it is an ascii file. It is not done in copy_params_PSRData or cleanPSRData, which means that otherwise memory is assumed to be allocated.
    }

    if(NoPolarizationFlag != 0 || polsplit != 0) { /*Set header if only stokes I is outputted */
      fout.NrPols = 1;
    }
    if(IndividualChannelFlag != 0) { /* Set header if freq. channels are outputted seperately */
      fout.NrFreqChan = 1;
    }

    Ipulse = (float *)malloc(fin.NrBins*sizeof(float));
    Ipulse2 = (float *)malloc(fin.NrBins*sizeof(float));
    if(Ipulse == NULL || Ipulse2 == NULL) {
      printerror(application.verbose_state.debug, "ERROR psplit: Memory allocation error\n");
      return 0;
    }

    BlockNumber = 0;
    pulse_in = 0;
    if(NrBlocksToOutput < 0) /* If not set, make sure it outputs all pulses */
      NrBlocksToOutput = fin.NrSubints;

    do {
      long outern, outer1, outer2;
      long innern, inner1, inner2;
      if(polsplit != 0) {
	outer1 = 0;
	outer2 = fin.NrPols;
	if(NoPolarizationFlag != 0)
	  outer2 = 1;

	inner1 = 0;
	inner2 = fin.NrFreqChan;
      }else {
	inner1 = 0;
	inner2 = fin.NrPols;
	if(NoPolarizationFlag != 0)
	  inner2 = 1;

	outer1 = 0;
	outer2 = fin.NrFreqChan;
      }

      for(outern = outer1; outern < outer2; outern++) {
	for(innern = inner1; innern < inner2; innern++) {
	  if(polsplit != 0) {
	    p = outern;
	    f = innern;
	  }else {
	    f = outern;
	    p = innern;
	  }


	  sprintf(txt, "block%05ld.gg", BlockNumber);
	  if(change_filename_extension(inputname, output_name, txt, 999, application.verbose_state) != 1) {
	    printerror(application.verbose_state.debug, "ERROR psplit: Changing extension failed\n");
	    return 0;
	  }
	  if(IndividualChannelFlag != 0) { /* Make filename for frequency channel */
	    sprintf(txt, "freq%05ld.gg", f);
	    strcpy(txt2, output_name);
	    if(change_filename_extension(txt2, output_name, txt, 999, application.verbose_state) != 1) {
	      printerror(application.verbose_state.debug, "ERROR psplit: Changing extension failed\n");
	      return 0;
	    }
	    fout.freqMode = FREQMODE_UNIFORM;
	    if(fout.freqlabel_list != NULL) {
	      free(fout.freqlabel_list);
	      fout.freqlabel_list = NULL;
	    }
	    //	  fout.freq_list = malloc(2*sizeof(double));
	    //	  if(fout.freq_list == NULL) {
	    //	    fflush(stdout);
	    //	    printerror(application.verbose_state.debug, "ERROR psplit: Memory allocation error.");
	    //	    return 0;
	    //	  }
	    set_centre_frequency(&fout, get_nonweighted_channel_freq(fin, f, application.verbose_state), application.verbose_state);
	    double chanbw;
	    if(get_channelbandwidth(fin, &chanbw, application.verbose_state) == 0) {
	      printerror(application.verbose_state.debug, "ERROR psplit (%s): Cannot obtain channel bandwidth.", fin.filename);
	      return 0;
	    }
	    if(set_bandwidth(&fout, chanbw, application.verbose_state) == 0) {
	      printerror(application.verbose_state.debug, "ERROR psplit (%s): Bandwidth changing failed.", fin.filename);
	      return 0;
	    }
	  }
	  if(polsplit) {
	    sprintf(txt, "pol%05ld.gg", p);
	    strcpy(txt2, output_name);
	    if(change_filename_extension(txt2, output_name, txt, 999, application.verbose_state) != 1) {
	      printerror(application.verbose_state.debug, "ERROR psplit: Changing extension failed\n");
	      return 0;
	    }
	  }

	  if((IndividualChannelFlag != 0 && p == 0) || (polsplit != 0 && f == 0) || (f == 0 && p == 0)) {
	    fout.NrSubints = BlockSize; /* Write out header for (part) of block */
	    if(pulse_in + BlockSize > fin.NrSubints) {
	      fout.NrSubints = fin.NrSubints-pulse_in;
	    }
	    if(fout.tsub_list != NULL)
	      free(fout.tsub_list);
	    fout.tsub_list = (double *)malloc(fout.NrSubints*sizeof(double));
	    fout.tsubMode = TSUBMODE_TSUBLIST;
	    for(n = 0; n < fout.NrSubints; n++) {
	      fout.tsub_list[n] = get_tsub(fin, pulse_in+n, application.verbose_state);
	    }
	    // Open the output file
	    if(openPSRData(&fout, output_name, fout.format, 1, 0, 0, application.verbose_state) == 0) {
	      printerror(application.verbose_state.debug, "ERROR psplit: Cannot open %s", output_name);
	      return 0;
	    }
	    if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) != 1) {
	      printerror(application.verbose_state.debug, "ERROR psplit: Cannot write header to %s", output_name);
	      return 0;
	    }
	    //	    if(application.verbose_state.debug) fprintf(stderr, "Opening file %s with %ld frequency channels\n", output_name, fout.NrFreqChan);
	    //	  if(f == 1) {
	    //	    closePSRData(&fout, application.verbose_state);
	    //	    return 0;
	    //	  }
	  }
	  if(IndividualChannelFlag != 0 && application.verbose_state.debug)
	    printf("Freq channel: %ld/%ld\n", f+1, fin.NrFreqChan);
	  if(polsplit != 0 && application.verbose_state.debug)
	    printf("Pol channel: %ld/%ld\n", p+1, fin.NrPols);
	  for(n = 0; n < fout.NrSubints; n++) {
	    if(application.verbose_state.verbose && application.verbose_state.nocounters == 0) {
	      if(application.verbose_state.debug) {
		printf("(pol=%ld freq=%ld) %f%%        \r", p, f, 100.0*(p*fout.NrSubints+n+ 1)/(float)(fout.NrSubints*fin.NrPols));
	      }else {
		printf("%f%%        \r", 100.0*(p*fout.NrSubints+n+ 1)/(float)(fout.NrSubints*fin.NrPols));
	      }
	    }
	    //	      puma_filepos = PuMaFilepos(pulse_in+n, f, p, fin.NrSubints, fin.NrBins, fin.NrFreqChan, data_start);
	    //	    ret = fseeko(fin, puma_filepos, SEEK_SET); /* Seek for pulse in input */
	    //						if(ret != 0) {
	    //							fprintf(stderr, "Cannot seek in puma file (byte %lld).\n", puma_filepos);
	    //							return -1;
	    //						}
	    //						if(IndividualChannelFlag == 0)
	    //							puma_filepos = PuMaFilepos(n, f, p, hdr.redn.NTimeInts, fin.NrBins, fin.NrFreqChan, data_start);
	    //						else
	    //							puma_filepos = PuMaFilepos(n, 0, p, hdr.redn.NTimeInts, fin.NrBins, 1, data_start);
	    //						ret = fseeko(fout, puma_filepos, SEEK_SET); /* Seek for pulse in output */
	    //						if(ret != 0) {
	    //							fprintf(stderr, "Cannot seek in puma file (byte %lld).\n", puma_filepos);
	    //							return -1;
	    //						}
	    if(readPulsePSRData(&fin, pulse_in+n, p, f, 0, fin.NrBins, Ipulse, application.verbose_state) != 1) {
	      printerror(application.verbose_state.debug, "ERROR psplit: Read error");
	      return 0;
	    }
	    //	    if(ShiftFlag) {
	      //	      pumawrite(Ipulse+NrShiftBins, sizeof(float), fin.NrBins-NrShiftBins, fout);
	      //	      puma_filepos = PuMaFilepos(pulse_in+n+1, f, p, fin.NrSubints, fin.NrBins, fin.NrFreqChan, data_start);
	      //	      ret = fseeko(fin, puma_filepos, SEEK_SET); /* Seek for pulse in input */
	      //	      if(ret != 0) {
		//		fprintf(stderr, "Cannot seek in puma file (byte %lld).\n", puma_filepos);
		//		return -1;
		//	      }
	      //	      pumaread(Ipulse, sizeof(float), fin.NrBins, fin);
	      /* int junk;
		 for(junk=0; junk < fin.NrBins; junk++)
		 Ipulse[junk] = 0; */
	      //	      pumawrite(Ipulse, sizeof(float), NrShiftBins, fout);
	      //	    }else {
/* Write a single pulse/subint starting at binnr and with length
   nrSamples. The behaviour of this function is unpredictable if you
   read in beyond a subint/polarization/frequency channel. The
   writing should not nessesarily be from start to end. Return 1 if
   successful. */
	    if(IndividualChannelFlag == 0 && polsplit == 0) {
	      if(writePulsePSRData(&fout, n, p, f, 0, fin.NrBins, Ipulse, application.verbose_state) != 1) {
		printerror(application.verbose_state.debug, "ERROR psplit: Write error");
		return 0;
	      }
	    }else if(IndividualChannelFlag != 0) {
	      if(writePulsePSRData(&fout, n, p, 0, 0, fin.NrBins, Ipulse, application.verbose_state) != 1) {
		printerror(application.verbose_state.debug, "ERROR psplit: Write error");
		return 0;
	      }
	    }else if(polsplit != 0) {
	      if(application.verbose_state.debug) fprintf(stderr, "Writing subint %ld, freq chan %ld, pol chan %ld\n", n, f, p);
	      if(writePulsePSRData(&fout, n, 0, f, 0, fin.NrBins, Ipulse, application.verbose_state) != 1) {
		printerror(application.verbose_state.debug, "ERROR psplit: Write error");
		return 0;
	      }
	    }
	      //	    }
	  } // End of subint loop
	} // End of inner loop (default is polarization) loop
	if(polsplit != 0) { /* If polarizations are split, outer loop are the polarizations */
	  closePSRData(&fout, 1, application.verbose_state);  // Header information is probably required for next block
	  if(application.verbose_state.verbose)
	    printf("\n");
	}
	if(polsplit == 0 && (IndividualChannelFlag != 0 || (f == fin.NrFreqChan - 1))) { /* If freq separation or writen out last channel */
	  closePSRData(&fout, 1, application.verbose_state);  // Header information is probably required for next block
	  if(application.verbose_state.verbose)
	    printf("\n");
	}
      } // End of outer loop (default is frequency) loop
      if (RunningBlock != 0) {
	pulse_in += 1;
      }else {
	pulse_in += BlockSize;
      }
      BlockNumber++;
      if(BlockNumber >= NrBlocksToOutput) /* If done, make sure we quit the loop */
	pulse_in = fin.NrSubints+100;
    }while( (RunningBlock != 0 && pulse_in <= (fin.NrSubints-BlockSize))  || (RunningBlock == 0 && pulse_in < fin.NrSubints) );
    closePSRData(&fin, 0, application.verbose_state);
    free(Ipulse);
    free(Ipulse2);
  }
  closePSRData(&fout, 0, application.verbose_state);
  terminateApplication(&application);
  return 0;
}
