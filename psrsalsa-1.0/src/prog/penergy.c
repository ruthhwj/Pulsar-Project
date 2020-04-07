//START REGION RELEASE
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "psrsalsa.h"

int writePSRCHIVE_ASCIIHeader(datafile_definition datafile, verbose_definition verbose);

//START REGION DEVELOP
void SHOWREVISIONINFO_prog() {
#include "penergy.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);
//START REGION RELEASE


int main(int argc, char **argv)
{ 
  FILE    *ofile;
  float   energy, on_totenergy, off_totenergy,on_peakenergy;
  float   off_peakenergy, on_rmsenergy, off_rmsenergy, s2n;
  float   *Ipulse, *output_values; 
  float   snrTresh, maxbin, maxbin2;
  long    i, j, polnr, freqnr, subintnr, okflag, freq1, freq2;
  int     binnr, filename, init_nrRegions, polmode, burstmode, refine, allwidths, posOrNeg, only_onpulse, nodebase;
//START REGION DEVELOP
  int force_binrange, selectMoreOnpulseRegions;
  int find_amplitudes, oldstyle;
//START REGION RELEASE
  int     used_onpulse_bins, individual_bin_mode;
  int     used_offpulse_bins;
  int     output2file, suppressNotFiniteS2Nwarnings;
  char    *oname, *filename_ptr, output_suffix[100];
  datafile_definition datain, opfile, pulse_profile;
  psrsalsaApplication application;
  pgplot_options_definition pgplot_options;
  //  long i;
  int index;
  int guessing_format;

  initApplication(&application, "penergy", "[options] inputfile(s)");
//START REGION DEVELOP
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
//START REGION RELEASE
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.oformat = FITS_format;
  application.switch_formatlist = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
//START REGION DEVELOP
  application.switch_showrevision = 1;
//START REGION RELEASE
  application.switch_forceUniformFreqLabelling = 1;
  application.switch_onpulse = 1; 
  application.switch_onpulsef = 1;
  application.switch_nocounters = 1;
  application.switch_device = 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_tscr = 1;
  application.switch_tscr_complete = 1;
  application.switch_TSCR = 1;
  application.switch_fscr = 1;
  application.switch_FSCR = 1;
  application.switch_debase = 1;
  application.switch_history_cmd_only = 1;
  application.switch_rebin = 1;
  application.switch_nskip = 1;
  application.switch_nread = 1;
  application.switch_scale = 1;
  application.switch_noweights = 1;
  application.switch_useweights = 1;
  application.switch_uniformweights = 1;

  snrTresh = 1;
  output2file = 1;
  individual_bin_mode = 0;
  filename = 0;
  strcpy(output_suffix, "en");
  suppressNotFiniteS2Nwarnings = 0;
  polmode = 0;
  burstmode = 0;
  freq1 = freq2 = -1;
  maxbin = 1e30;
  maxbin2 = 0;
  refine = 0;
  allwidths = 0;
  posOrNeg = 0;
  only_onpulse = 0;
  nodebase = 1;  // Set to zero if this command line option is released
  pgplot_clear_options(&pgplot_options);  
//START REGION DEVELOP
  nodebase = 0;  // In release version default is to not remove baseline, while in develop version baseline is removed by default.
  selectMoreOnpulseRegions = 0;
  force_binrange = 0;
  find_amplitudes = 0;
  oldstyle = 0;
//START REGION RELEASE

  if(argc < 2) {
    printf("Program for calculating pulse energy statistics. Use pdist to make distributions. The program can be run in two modes:\n\n1) Without the -burst option: Calculate burst statistics using the first selected on-pulse region. The not selected regions are used as the off-pulse regions only used for S/N determinations. When -debase is used, the baseline is determined from the off-pulse region as well. Unless specified otherwise, an ascii file will be generated.\n\n2) With the burst option: By default this program tries a fixed number of widths as a matched filter to search for bursts in the signal. Unless -only_onpulse or -only_onpulse1 is specified, the off-pulse region is only used to get the rms of the noise.\n\n");
    printApplicationHelp(&application);
    printf("Action options mode 1:\n\n");
    printf("-dynspec           Output a PSRFITS file (change format with -oformat) with over\n");
    printf("                   onpulse region integrated energies rather than an ascii file.\n");
    printf("-dynspec2          Output a PSRFITS file (change format with -oformat) with S/N\n");
    printf("                   rather than an ascii file.\n");
    printf("-b                 Output pulse energy stats for individual bins rather than\n");
    printf("                   considering the (first) onpulse region as a whole.\n");
//START REGION DEVELOP
    printf("-force             Force the number of used onpulse bins in the selected onpulse\n");
    printf("                   range to be this fixed amount, independent to actual number of\n");
    printf("                   selected bins.\n");
    printf("-onpulsegr         Enables selecting more on-pulse regions graphically than\n");
    printf("                   defined by -onpulse\n");
//START REGION RELEASE
    printf("-pol               Analyse all polarization channels (only first otherwise)\n");
    printf("\nAction options mode 2 (use -burst):\n\n");
    printf("-burst             Enable this mode\n");
    printf("-burst_all         Slower than -burst_refine, as all pulse widths are tried.\n");
//START REGION DEVELOP
    printf("-burst_amp         Prints out the value/location of the bin with the heighest\n");
    printf("                   intensity for each subintegration.\n");
//START REGION RELEASE
    printf("-burst_maxbin      Specify maximum width in bins allowed for a detection. If the\n");
    printf("                   brightest detection is wider, the subint will be ignored.\n");
    printf("-burst_maxbin2     Similar to -burst_maxbin, except that if there is a detection\n");
    printf("                   with a smaller width, it is accepted.\n");
    printf("-burst_neg         Allow the detection of dips as well as spikes.\n");
    printf("-burst_onpulse     Only look for bursts in the onpulse region(s).\n");
    printf("-burst_onpulse1    Only look for bursts in the first specified onpulse region.\n");
    printf("-burst_refine      Refine detection trying other widths close to the found\n");
    printf("                   detection. This is slower. Run penergy without this option\n");
    printf("                   (with -v) to see which pulse widths are considered otherwise.\n");
    printf("-burst_snr         Specify treshhold S/N [default is %.1f].\n", snrTresh);
//START REGION DEVELOP
    printf("-burst_nodebase    Assume the baseline is already removed, otherwise the\n");  // If releasing, it should also work for mode 1)
    printf("                   off-pulse region is used.\n");
//START REGION RELEASE
    printf("\nOutput file options:\n\n");
    printf("-ext              Specify suffix, default is '%s'\n", output_suffix);
    printf("-freq \"l h\"       Specify freq. channels to use (from l to h, default is all).\n");
    printf("-output filename  Write output to this file [def=.%s extension].\n", output_suffix);
//START REGION DEVELOP
    printf("-oldstyle         Use the old style of ascii output.\n");
//START REGION RELEASE
    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about using pulse energy distributions can be found in:\n");
    printf(" - Weltevrede et al. 2006, A&A, 458, 269\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-output") == 0) {
	filename = i+1;
	i++;
      }else if(strcmp(argv[i], "-dynspec2") == 0) {
	//	output2puma = 1;
	output2file = 3;
      }else if(strcmp(argv[i], "-dynspec") == 0) {
	//	output2puma = 2;
	output2file = 2;
      }else if(strcmp(argv[i], "-b") == 0) {
	individual_bin_mode = 1;
      }else if(strcmp(argv[i], "-ext") == 0) {
	if(strlen(argv[i+1]) > 99) {
	  printerror(application.verbose_state.debug, "ERROR penergy: Error parsing option '%s' - string too long", argv[i]);
	  return 0;
	}
	strcpy(output_suffix, argv[i+1]);
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-force") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &force_binrange, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR penergy: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-onpulsegr") == 0) {
	selectMoreOnpulseRegions = 1; 
      }else if(strcmp(argv[i], "-oldstyle") == 0) {
	oldstyle = 1; 
//START REGION RELEASE
      }else if(strcmp(argv[i], "-pol") == 0) {
	polmode = 1; 
      }else if(strcmp(argv[i], "-burst") == 0) {
	burstmode = 1; 
      }else if(strcmp(argv[i], "-burst_refine") == 0) {
	refine = 1;
      }else if(strcmp(argv[i], "-burst_all") == 0) {
	allwidths = 1;
      }else if(strcmp(argv[i], "-burst_neg") == 0) {
	posOrNeg = 1;
      }else if(strcmp(argv[i], "-burst_onpulse") == 0) {
	only_onpulse = 1;
      }else if(strcmp(argv[i], "-burst_onpulse1") == 0) {
	only_onpulse = 2;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-burst_amp") == 0) {
	find_amplitudes = 1;
      }else if(strcmp(argv[i], "-burst_nodebase") == 0) {
	nodebase = 1;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-burst_snr") == 0 || strcmp(argv[i], "-burst_s2n") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &snrTresh, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR penergy: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-burst_maxbin") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &maxbin, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR penergy: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-burst_maxbin2") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &maxbin2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR penergy: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-freq") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld %ld", &freq1, &freq2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR penergy: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else {
	/* If the option is not recognized, assume it is a filename */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "ERROR penergy: Unknown option: %s\n\nRun penergy without command line arguments to show help", argv[i]);
	  terminateApplication(&application);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }

  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR penergy: No input file(s) specified");
    return 0;
  }

  init_nrRegions = application.onpulse.nrRegions;

  guessing_format = 0;

  while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {

    if(application.verbose_state.verbose)
      printf("Reading from %s\n", filename_ptr);

    oname = (char *) calloc(strlen(filename_ptr)+strlen(output_suffix)+8, 1);
    if(oname == NULL) {
      printerror(application.verbose_state.debug, "ERROR penergy: Memory allocation error.\n");
      return 0;
    }
    //    cleanPSRData(&datain, application.verbose_state);
    /* Opening the file and checking all is well */
    if(guessing_format) {
      application.iformat = -1;
    }
    if(application.iformat <= 0) {
      application.iformat = guessPSRData_format(filename_ptr, 0, application.verbose_state);
      if(application.iformat == -2 || application.iformat == -3)
	return 0;
      guessing_format = 1;
    }
    if(isValidPSRDATA_format(application.iformat) == 0) {
      printerror(application.verbose_state.debug, "ERROR penergy: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", filename_ptr);
      return 0;
    }
    if(openPSRData(&datain, filename_ptr, application.iformat, 0, 1, 0, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR penergy: Error opening data");
      return 0;
    }

    /* Search commandline for header parameters to overwrite header
     parameters. */
    if(PSRDataHeader_parse_commandline(&datain, argc, argv, application.verbose_state) == 0)
      return 0;

    /* Convert the specified onpulse regions which were defined as
       fractions on the commandline. This is now possible because we
       know the number of bins. */
    region_frac_to_int(&(application.onpulse), datain.NrBins, 0);
    
    for(j = 1; j < argc; j++) {
      if(strcmp(argv[j], "-header") == 0) {
	printwarning(application.verbose_state.debug, "WARNING: If using the -header option, be aware it applied BEFORE the preprocessing.");
      }
    }

    //    fprintf(stderr, "XXXXXX %s\n", filename_ptr);
    if(preprocessApplication(&application, &datain) == 0) {
      return 0;
    }
    
    regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);

    // Default: process all frequency channels
    if(freq1 < 0 || freq2 < 0) {
      freq1 = 0;
      freq2 = datain.NrFreqChan-1;
    }

    if(application.onpulse.nrRegions == 0
//START REGION DEVELOP
       || selectMoreOnpulseRegions || force_binrange
//START REGION RELEASE
       ) {
      if(preprocess_make_profile(datain, &pulse_profile, 1, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR penergy: Cannot construct pulse profile");
	return 0;
      }
    }

    if(application.onpulse.nrRegions == 0
//START REGION DEVELOP
       || selectMoreOnpulseRegions
//START REGION RELEASE
       ) {
      strcpy(pgplot_options.viewport.plotDevice, application.pgplotdevice);
      strcpy(pgplot_options.box.xlabel, "Bin");
      strcpy(pgplot_options.box.ylabel, "Intensity");
      sprintf(pgplot_options.box.title, "Select on-pulse region of %s (%s)", datain.psrname, datain.filename);
      if(burstmode == 0) {
	printwarning(application.verbose_state.debug, "\nYou can select multiple onpulse regions, but only the first is used to calculate onpulse statistics. Non-selected regions are used for the off-pulse statistics.\n\n");
      }else if (only_onpulse == 0) {
	printwarning(application.verbose_state.debug, "\nThe non-selected regions are only used for the off-pulse statistics, unless -burst_onpulse or -burst_onpulse1 is used.\n\n");
      }
      selectRegions(pulse_profile.data, datain.NrBins, &pgplot_options, 0, 0, 0, &(application.onpulse), application.verbose_state);
//START REGION DEVELOP
      if(force_binrange) {
	application.onpulse.right_bin[0] = application.onpulse.left_bin[0]+force_binrange-1;
	if(application.onpulse.right_bin[0] >= datain.NrBins) {
	  printerror(application.verbose_state.debug, "penergy: Selected range is out of bounds");
	  return 0;
	}
      }
//START REGION RELEASE
      regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);
    }

    if(application.onpulse.nrRegions == 0) {
      printerror(application.verbose_state.debug, "ERROR penergy: At least one on-pulse region needs to be defined. Use the -onpulse option or select at least one region in pgplot window.");
      return 0;
    }
    if(application.onpulse.bins_defined[0] == 0) {
      printerror(application.verbose_state.debug, "ERROR penergy: The onpulse region is not defined in bins.");
      return 0;
    }

    Ipulse = malloc(datain.NrBins*sizeof(float));            // Input: read NrBins of data at a time
    if(Ipulse == NULL) { // Iprofile == NULL || 
      printerror(application.verbose_state.debug, "ERROR penergy: Cannot allocate memory");
      return 0;
    }

    if(burstmode == 0) {  // The original penergy code. The pfindbursts code is the other option.
      output_values = malloc(datain.NrSubints*sizeof(float));  // Need construct full NrBins worth of output data to write to a dynamic spectrum file
      if(output_values == NULL) { // Iprofile == NULL || 
	printerror(application.verbose_state.debug, "ERROR penergy: Cannot allocate memory");
	return 0;
      }
      
//START REGION DEVELOP
      if(force_binrange) {
	if(application.onpulse.nrRegions != 1) {
	  printerror(application.verbose_state.debug, "ERROR penergy: Only one region is expected when the -force option is selected.");
	  return 0;
	}
	strcpy(pgplot_options.viewport.plotDevice, application.pgplotdevice);
	strcpy(pgplot_options.box.xlabel, "Bin");
	strcpy(pgplot_options.box.ylabel, "Intensity");
	sprintf(pgplot_options.box.title, "The selected region after applying -force option, press S to continue");
	selectRegions(pulse_profile.data, datain.NrBins, &pgplot_options, 1, 0, 0, &(application.onpulse), application.verbose_state);
	printf("Only showing selected range. Do not change the selected range.\n");
	fflush(stdout);
	if(application.onpulse.right_bin[0] - application.onpulse.left_bin[0] != force_binrange-1) {
	  printerror(application.verbose_state.debug, "ERROR penergy: Selected range is not the amount specified with -force option.");
	  return 0;
	}
      }
//START REGION RELEASE

      if(application.onpulse.nrRegions == 0) {
	printerror(application.verbose_state.debug, "ERROR penergy: At least one on-pulse region needs to be defined. Use the -onpulse option or select at least one region in pgplot window.");
	return 0;
      }
      if(application.onpulse.bins_defined[0] == 0) {
	printerror(application.verbose_state.debug, "ERROR penergy: The onpulse region is not defined in bins.");
	return 0;
      }
      if(application.onpulse.nrRegions > 1) {
	printwarning(application.verbose_state.debug, "WARNING penergy: When selecting multiple onpulse regions, only the first is used to calculate onpulse statistics. Non-selected regions are used for the off-pulse statistics.");
      }
      
      for(binnr = application.onpulse.left_bin[0]; binnr <= application.onpulse.right_bin[0]; binnr++) {   /* Bin loop is here, in case each bin is a separate output file */

      // Construct output filename
	strcpy(oname, filename_ptr);
	if(filename != 0) {
	  strcpy(oname, argv[filename]);
	}else {
	  if(individual_bin_mode != 0) {
	    sprintf(oname,"%s.%05d.%s", filename_ptr, binnr, output_suffix);
	  }else {
	    strcat(oname,".");
	    strcat(oname, output_suffix);
	  }
	}
      
	// Open output file 
	if(output2file == 1) {
	  printf("Output Ascii to %s\n",oname);
	  ofile = fopen(oname,"w+");
	  if(ofile == NULL) {
	    printerror(application.verbose_state.debug, "Error opening '%s'", oname);
	    perror("");
	    return 0;
	  }
//START REGION DEVELOP
	  if(oldstyle) {
	    fprintf(ofile,"#Used on-pulse region: bin %d,%d\n", application.onpulse.left_bin[0], application.onpulse.right_bin[0]);
	    fprintf(ofile,"#Used off-pulse region: everyting except");
	    for(j = 0; j < application.onpulse.nrRegions; j++) {
	      fprintf(ofile," %d,%d", application.onpulse.left_bin[j], application.onpulse.right_bin[j]);
	    }
	    fprintf(ofile,"\n");
	    fprintf(ofile,"#SUBINT POLN FREQ: Peak intensity (on & off-pulse) Integrated energy (on & off-pulse) RMS (on & off-pulse) S/N\n");
	  }else {
//START REGION RELEASE
	    cleanPSRData(&opfile, application.verbose_state);
	    copy_params_PSRData(datain, &opfile, application.verbose_state);
	    //	    fprintf(stderr, "XXXX copy_params %d %d\n", output2file, oldstyle);
	    opfile.fptr_hdr = ofile;
	    opfile.NrBins = datain.NrSubints;          // Each bin is a subintegration
	    opfile.NrSubints = datain.NrPols;
	    if(polmode == 0)
	      opfile.NrSubints = 1;               // If only writing Stokes I
	    opfile.NrFreqChan = freq2-freq1+1;
	    opfile.NrPols = 7;    // Nr of columns of data

	    if(set_filename_PSRData(&opfile, "penergy", application.verbose_state) == 0) {
	      printerror(application.verbose_state.debug, "ERROR penergy: Setting filename failed");
	      return 0;
	    }
	    if(writePSRCHIVE_ASCIIHeader(opfile, application.verbose_state) == 0) {
	      printerror(application.verbose_state.debug, "ERROR penergy: Writing header line failed");
	      return 0;
	    }
//START REGION DEVELOP
	  }
//START REGION RELEASE
	}else if(output2file > 1) {  // Write out data as pulsar data
	  printf("Output to %s file %s\n", returnFileFormat_str(application.oformat), oname);
	  cleanPSRData(&opfile, application.verbose_state);
	  copy_params_PSRData(datain, &opfile, application.verbose_state);
	  
	  if(openPSRData(&opfile, oname, application.oformat, 1, 0, 0, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR penergy: Cannot open %s", oname);
	    return 0;
	  }
	  opfile.NrPols = 1;   /* Only write Stokes I */
	  opfile.gentype = GENTYPE_DYNAMICSPECTRUM;  // Output will be a dynamic spectrum
	  opfile.NrBins = opfile.NrSubints;          // Each bin is a subintegration
	  opfile.NrSubints = 1;
	  opfile.tsub_list[0] = get_tobs(datain, application.verbose_state);
	  opfile.tsubMode = TSUBMODE_FIXEDTSUB;
	  opfile.tsampMode = TSAMPMODE_FIXEDTSAMP;
	  opfile.fixedtsamp = get_tobs(datain, application.verbose_state)/(double)opfile.NrBins;
	  if(writeHeaderPSRData(&opfile, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR penergy: Cannot write header to %s", oname);
	    return 0;
	  }
	}
      
	on_totenergy = 0;
	off_totenergy = 0;
	on_peakenergy = 0;
	off_peakenergy = 0;
	on_rmsenergy = 0;
	off_rmsenergy = 0;
	s2n = 0;
	long nrpolstoconsider;
	nrpolstoconsider = datain.NrPols;
	if(polmode == 0)
	  nrpolstoconsider = 1;
	for(polnr = 0; polnr < nrpolstoconsider; polnr++) {
	  for(freqnr = freq1; freqnr <= freq2; freqnr++) {
	    for(subintnr = 0; subintnr < datain.NrSubints; subintnr++) {
	      if(readPulsePSRData(&datain, subintnr, polnr, freqnr, 0, datain.NrBins, Ipulse, application.verbose_state) != 1) {
		printerror(application.verbose_state.debug, "ERROR penergy: read error, technically not possible");
		return 0;
	      }
	      used_onpulse_bins = used_offpulse_bins = 0;
	      for(j = 0; j < datain.NrBins; j++) {
		energy = Ipulse[j];	
		okflag = 0;
		if(individual_bin_mode != 0 && j == binnr) {
		  okflag = 1;
		}else if(individual_bin_mode == 0 && checkRegions(j, &(application.onpulse), 1, application.verbose_state) != 0)
		  okflag = 1;
		if(okflag != 0) {
		  on_totenergy += energy;
		  on_rmsenergy += energy*energy;
		  if(energy > on_peakenergy)
		    on_peakenergy=energy;
		  used_onpulse_bins++;
		}
		okflag = 0;
		if(checkRegions(j, &(application.onpulse), 0, application.verbose_state) == 0)
		  okflag = 1;
		if(okflag != 0) {
		  off_totenergy += energy;
		  off_rmsenergy += energy*energy;
		  if(energy > off_peakenergy)
		    off_peakenergy=energy;
		  used_offpulse_bins++;
		}
	      }
	      if(freqnr == 0 && subintnr == 0) 
		printf("Using %d on-pulse and %d off-pulse bins\n", used_onpulse_bins, used_offpulse_bins);
	      on_rmsenergy = sqrt(on_rmsenergy)/sqrt(used_onpulse_bins);
	      if(used_offpulse_bins > 0) {
		off_rmsenergy = sqrt(off_rmsenergy)/sqrt(used_offpulse_bins);
		s2n = on_totenergy/(off_rmsenergy*sqrt(used_onpulse_bins));
	      }else {
		off_rmsenergy = -1;
		s2n = -1;
	      }

	      if(output2file == 1) {
//START REGION DEVELOP
		if(oldstyle) {
		  fprintf(ofile, "%ld %ld %ld: %e %e %e %e %e %e %e\n",subintnr,polnr,freqnr,on_peakenergy, off_peakenergy,on_totenergy, off_totenergy, on_rmsenergy, off_rmsenergy, s2n);
		}else {
//START REGION RELEASE
		  fprintf(ofile, "%ld %ld %ld %e %e %e %e %e %e %e\n", polnr, freqnr-freq1, subintnr, on_peakenergy, off_peakenergy,on_totenergy, off_totenergy, on_rmsenergy, off_rmsenergy, s2n);
//START REGION DEVELOP
		}
//START REGION RELEASE
	      }else if(output2file > 1) { // If constructing psrdata file, need to store NrBins of output data first (for FITS)
		if(output2file == 3) {
		  if(isnan(s2n) || isinf(s2n)) {
		    if(suppressNotFiniteS2Nwarnings == 0) {
		      fflush(stdout);
		      printwarning(application.verbose_state.debug, "WARNING penergy: s2n = %f will be replaced with 0.", s2n);
		      suppressNotFiniteS2Nwarnings = 1;
		    }
		    s2n = 0;
		  }
		  output_values[subintnr] = s2n;
		}else {
		  output_values[subintnr] = on_totenergy;
		}
	      }
	      
	      on_totenergy = 0;
	      off_totenergy = 0;
	      on_peakenergy = 0;
	      off_peakenergy = 0;
	      on_rmsenergy = 0;
	      off_rmsenergy = 0;
	      s2n = 0;
	    } // End of subint loop
	    
	    // Write out complete output "subint"
	    if(output2file > 1) {
	      if(writePulsePSRData(&opfile, 0, polnr, freqnr, 0, opfile.NrBins, output_values, application.verbose_state) != 1) {
		printerror(application.verbose_state.debug, "ERROR penergy: Write error");
		return 0;
	      }
	    }
	    
	  } // End of frequency loop
	} // End of polarization loop
	
	if(output2file == 1) {
//START REGION DEVELOP
	  if(oldstyle == 0) {
//START REGION RELEASE
	    char *txt;
	    txt = malloc(10000);
	    if(txt == NULL) {
	      fflush(stdout);
	      printerror(application.verbose_state.debug, "ERROR penergy: Memory allocation error.");
	      return 0;
	    }
	    constructCommandLineString(txt, 10000, argc, argv, application.verbose_state);
	    fprintf(ofile, "#%s\n", txt);
	    free(txt);
	    fprintf(ofile,"#Used on-pulse region: bin %d,%d\n", application.onpulse.left_bin[0], application.onpulse.right_bin[0]);
	    fprintf(ofile,"#Used off-pulse region: everyting except");
	    for(j = 0; j < application.onpulse.nrRegions; j++) {
	      fprintf(ofile," %d,%d", application.onpulse.left_bin[j], application.onpulse.right_bin[j]);
	    }
	    fprintf(ofile,"\n");
	    fprintf(ofile,"#POLN FREQ SUBINT Peak intensity (on & off-pulse) Integrated energy (on & off-pulse) RMS (on & off-pulse) S/N\n");
	    if(application.verbose_state.verbose) {
	      printf("Writing done, see end of file for explanation of the columns.\n");
	    }
//START REGION DEVELOP
	  }
//START REGION RELEASE
	  fclose(ofile);
	}
	if(output2file > 1 || (output2file == 1
//START REGION DEVELOP
			       && oldstyle == 0
//START REGION RELEASE
			       )) {
	  //	  fprintf(stderr, "XXXX closing opfile\n");
	  closePSRData(&opfile, 0, application.verbose_state);
	}
	/* Makes sure we quit the loop over the bins if we only want to output one file */
	if(individual_bin_mode == 0)      
	  break;
      } // End of bin loop
    }else {    // If burstmode is set
      if(polmode) {
	printerror(application.verbose_state.debug, "ERROR penergy: The -burst mode only opperates on the first polarization channel");
	return 0;
      }


      for(freqnr = freq1; freqnr <= freq2; freqnr++) {

	// Construct output filename
	strcpy(oname, filename_ptr);
	if(filename != 0) {
	  strcpy(oname, argv[filename]);
	}else {
	  if(datain.NrFreqChan > 1) {
	    sprintf(oname,"%s.%05ld.%s", filename_ptr, freqnr, output_suffix);
	  }else {
	    strcat(oname,".");
	    strcat(oname, output_suffix);
	  }
	}
	if(application.verbose_state.verbose)
	  printf("Output to: %s\n", oname);
	ofile = fopen(oname, "w");
	if(ofile == NULL) {
	  printf("Cannot open %s\n\n", oname);
	  return 0;
	}
//START REGION DEVELOP
	if(find_amplitudes) {
	  fprintf(ofile, "#pulsenr, bin, intensity\n");
	}else {
//START REGION RELEASE
	  fprintf(ofile, "#pulsenr, bin (left edge), pulse phase (center), width (bins), snr, integrated energy\n");
//START REGION DEVELOP
	}
//START REGION RELEASE

	for(subintnr = 0; subintnr <= (datain.NrSubints-1); subintnr++) {
	  if(application.verbose_state.verbose && (subintnr % 10 == 0) && subintnr > 0 && application.verbose_state.nocounters == 0) {
	    printf("\r%f%%        ", 100.0*(subintnr+1)/(float)((datain.NrSubints-1)+1.0));
	    fflush(stdout);
	  }
	  //	  fprintf(stderr, "XXXX: subintnr=%ld freqnr=%ld NrBins=%ld\n", subintnr, freqnr, datain.NrBins);
	  if(readPulsePSRData(&datain, subintnr, 0, freqnr, 0, datain.NrBins, Ipulse, application.verbose_state) != 1) {
	    printerror(application.verbose_state.debug, "ERROR penergy: Reading failed.\n");
	    return 0;
	  }
//START REGION DEVELOP
	  if(find_amplitudes) {
	    int firsttime, ok;
	    long binnr_max;
	    float max;
	    firsttime = 1;
	    for(binnr = 0; binnr < datain.NrBins; binnr++) {
	      ok = 0;
	      if(only_onpulse == 0)
		ok = 1;
	      if(checkRegions(binnr, &(application.onpulse), 1, application.verbose_state) != 0)
		ok = 1;
	      if(checkRegions(binnr, &(application.onpulse), 0, application.verbose_state) != 0 && only_onpulse == 1)
		ok = 1;
	      if(ok) {
		if(Ipulse[binnr] > max || firsttime == 1) {
		  max = Ipulse[binnr];
		  firsttime = 0;
		  binnr_max = binnr;
		}
	      }
	    }
	    fprintf(ofile, "%ld %ld %e\n", subintnr, binnr_max, max);
	  }else {
//START REGION RELEASE
	    float E, snr;
	    int  binnr_max, pulsewidth;
	    verbose_definition verbose2;
	    copyVerboseState(application.verbose_state, &verbose2);
	    verbose2.verbose = application.verbose_state.verbose;
	    if(freqnr != 0 || subintnr != 0)
	      verbose2.verbose = 0;
	    boxcarFindpeak(Ipulse, datain.NrBins, &(application.onpulse), &binnr_max, &pulsewidth, &snr, &E, 0, posOrNeg, allwidths, refine, maxbin2, only_onpulse, nodebase, verbose2);
	    //	fprintf(stderr, "snr=%f  pulsewidth=%d\n", snr, pulsewidth);
	    if(snr > snrTresh && pulsewidth <= maxbin)
	      fprintf(ofile, "%ld %d %e %d %e %e\n", subintnr, binnr_max, ((float)binnr_max+0.5*(float)pulsewidth)/(float)datain.NrBins, pulsewidth, snr, E);
      }
//START REGION DEVELOP
    }
//START REGION RELEASE
    fflush(stdout);
    if(application.verbose_state.verbose) printf("Done\n");
    fclose(ofile);
  }
    }
    closePSRData(&datain, 0, application.verbose_state);
    free(Ipulse);
    if(burstmode == 0)
      free(output_values);
    free(oname);

    /* If there are more files to be opened, reset the number of regions to the initial ones that were specified */
    application.onpulse.nrRegions = init_nrRegions;
  } // Loop over input files

  terminateApplication(&application);
  return 0;
}
