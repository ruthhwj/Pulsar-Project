//START REGION RELEASE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "psrsalsa.h"

#define MaxNrPolarizations 5  // These are the max nr of polarizations that can be handled if the data-set needs to be separated in separate polarizations first (i.e. p3 folding)


//START REGION DEVELOP
void SHOWREVISIONINFO_prog() {
#include "pfold.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);


void display_parameters(carousel_parameter_list plist, verbose_definition verbose);

//START REGION RELEASE
int main(int argc, char **argv)
{
  int index, originalNrPols, originalNrPolsP3;
  long i, p;
  int p3_fold_flag, p3_fold_refine, p3_fold_cpb, p3_fold_nbin, p3_fold_onpulse_flag;
  int write_flag, zoom_flag, zoom_flag1, selectMoreOnpulseRegions;
  float p3_fold, p3_fold_smoothWidth, p3fold_dphase, p3fold_nosmooth, slope;
  float xmin, xmax, xmin_zoom, xmax_zoom, *profileI, *p3foldmap, *p3foldmap2;
  char onpulseselectdevice[1000], p3fold_device[1000], outputname[1000];
  psrsalsaApplication application;
  pgplot_options_definition pgplot_options;
  verbose_definition noverbose;
  datafile_definition fin[MaxNrPolarizations], clone, fout;
//START REGION DEVELOP
  int p3_fold_inputfile, p3_fold_pulse_list_file, readwriteP3foldOffsets, writeP3foldPulseList;
  long j;
  carousel_parameter_list plist; // plist should be globally accessible?

//START REGION RELEASE
  initApplication(&application, "pfold", "[options] inputfile");
//START REGION DEVELOP
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
//START REGION RELEASE
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_formatlist = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.switch_nocounters = 1;
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_onpulse = 1;
  application.switch_onpulsef = 1;
  //  application.switch_itf = 1;
  application.switch_rebin = 1;
  application.switch_nskip = 1;
  application.switch_nread = 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_conshift= 1;
  application.switch_circshift= 1;
  //  application.switch_shuffle = 1;
  application.switch_libversions = 1;
  application.switch_history_cmd_only = 1;
//START REGION DEVELOP
  application.switch_polselect = 1;
  application.switch_ppgplot = 1;
  application.switch_showrevision = 1;


//START REGION RELEASE
  write_flag = 0;
  zoom_flag = 0;
  zoom_flag1 = 0;
  selectMoreOnpulseRegions = 0;
  slope = 0;
  p3fold_dphase = 0;
  p3fold_nosmooth = 0;
  p3_fold_flag = 0;
  p3_fold_cpb = 1;
  p3_fold_refine = 1;
  p3_fold_smoothWidth = -1;
  p3_fold_onpulse_flag = 1;
  sprintf(onpulseselectdevice, "?");
  sprintf(p3fold_device, "?");

  pgplot_clear_options(&pgplot_options);
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;
//START REGION DEVELOP
//START REGION RELEASE
  application.oformat = FITS_format;
//START REGION DEVELOP
  p3_fold_inputfile = 0;
  p3_fold_pulse_list_file = 0;
  readwriteP3foldOffsets = 0;
  writeP3foldPulseList = 0;
  // Assign default flag values
  plist.maxbrightness = -1.0;
  plist.first_pulse = -1;
  plist.last_pulse = -1;
  plist.A = plist.B = 1;
  plist.car = 0;
  plist.p3car = 0;
  plist.k0 = 0;
  plist.phi0 = 0;

//START REGION RELEASE
  if(argv[argc-1][0] == '-' &&  strcmp(argv[argc-1], "-formatlist") != 0 &&  strcmp(argv[argc-1], "-headerlist") != 0) {
    printerror(application.verbose_state.debug, "pfold: Last command line option is expected to be a file name.\nRun pfold without command line arguments to show help");
    terminateApplication(&application);
    return 0;
  }

  if(argc < 2) {
    printf("Program to fold (folded) single pulse data thereby visualising the subpulse modulation cycle.\n\n");
    printApplicationHelp(&application);
    printf("General options:\n");
    printf("  -w                  Write out the results to files.\n");
    printf("\nOutput options related to p3 folding:\n");
    printf("  -p3fold             \"P3 n\": Fold the data using this P3 value in pulse\n"); 
    printf("                      periods and the P3 cycle is divided in n bins in the final\n");
    printf("                      result. If n > P3, the different bins are not independent\n");
    printf("                      and you might want to use -p3fold_smooth option to make\n");
    printf("                      the effective resolution equal to P3.\n");
    printf("  -p3fold_dphase      Add this subpulse phase in degrees (can also use -slope).\n");
    printf("  -p3fold_norefine    Do not attemt to align subsequent blocks, i.e. fixed\n");
    printf("                      period folding\n");
    printf("  -p3fold_nritt       Set the number of itterations, which produces a template\n");
    printf("                      first thereby producing better results. Default is %d.\n", p3_fold_refine);
    printf("  -p3fold_cpb         Set the number of cycles per block used in the cross\n");
    printf("                      correlation used to compensate for P3 variations. More\n");
    printf("                      means more signal to correlate (more precise alignment of\n");
    printf("                      the blocks, less means less smearing in each block because\n");
    printf("                      of P3 variation within the block. Default is %d.\n", p3_fold_cpb);
    printf("  -p3fold_smooth      Replace the tophat weight used to assign power to the P3\n");
    printf("                      bins with a Gausian weight with this width in pulse\n");
    printf("                      periods. This could make oversampling look nicer and\n");
    printf("                      reduce the effective resolution. Example: if P3=10P, you\n");
    printf("                      could set n in the -p3fold option to 20, resulting in\n");
    printf("                      oversampling with a factor 2. By setting -p3fold_smooth\n");
    printf("                      to 2, the effective resolution is reduced by a factor 2\n");
    printf("                      because each input pulse is smeared out with this width\n");
    printf("                      in pulse periods.\n");
//START REGION DEVELOP
    printf("  -p3fold_nosmooth    Disable any weighting, all power ends up in a single bin.\n");
    printf("                      By default a given input pulse is split over the two\n");
    printf("                      nearest otput rows of the P3 fold.\n");
//START REGION RELEASE
    printf("  -p3fold_noonpulse   Ignore selected pulse longitude range, but use the full\n"); 
    printf("                      range when doing the cross correlations\n");
//START REGION DEVELOP
    printf("  -p3fold_writephases Write the found phase offsets to a file.\n");
    printf("  -p3fold_readphases  Read offsets from specified file (made with\n");
    printf("                      -p3fold_writephases) rather than finding them via\n");
    printf("                      correlations.\n");
    printf("  -p3fold_writepulselist filename   - Write out a file with the subpulse phases\n");
    printf("                      for each pulse. Requires the use of -p3fold_readphases.\n");
//START REGION RELEASE
    printf("  -slope              Subtract slope from subpulse phases (in degrees subpulse\n");
    printf("                      phase per degree pulse longitude).\n");
//START REGION DEVELOP
    printf("\nOutput options related to carousel mapping:\n");
    printf("  -carousel \"alpha beta size[pixels] radial_extent first_pulse last_pulse P3_hat\"\n");
    printf("                      Generate a carousel map with these parameters.\n");  // NEEDS A DESCRIPTION OF EACH PARAMETER
    printf("  -p3carousel \"alpha beta size=pixels radial_extent #_of_sparks\"\n");
    printf("                     Plot a carousel with these parameters from P3-folded data.\n");
    printf("  -abrot       Abnormal (\"retrograde\")carousel rotation.\n");
    printf("  -phi0        Introduce a shift in pulse phase (degrees). By default pulse longitue 180 is assumed to be in the meridional plane. This specifies an additional offset.\n");
//START REGION RELEASE
    printf("\nGraphics options:\n");
    printf("  -onpulsed           Set pgplot device for the selection of the onpulse region.\n");
    printf("  -p3foldd            Set pgplot device for the P3 fold map.\n");
//START REGION DEVELOP
    printf("  -zoom               Use selected region to zoom in.\n");
    printf("  -zoomw              Use selected region to zoom in, but show some offpulse as\n");
    printf("                      well.\n");
    printf("  -zoom1              Use only first selected region to zoom in.\n");
//START REGION RELEASE
    printf("  -onpulsegr          Enables graphical selection of additional on-pulse regions\n");
    printf("                      to those defined with the -onpulse option.\n");

    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else if(argc > 2 || strcmp(argv[argc-1], "-formatlist") == 0 || strcmp(argv[argc-1], "-headerlist") == 0) {
    int lastindex;
    lastindex = argc-1;
    if(strcmp(argv[argc-1], "-formatlist") == 0 || strcmp(argv[argc-1], "-headerlist") == 0)
      lastindex++;
    for(i = 1; i <= lastindex; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-onpulsed") == 0) {
	strcpy(onpulseselectdevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-p3foldd") == 0) {
	strcpy(p3fold_device, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-w") == 0) {
	write_flag = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-zoom") == 0) {
	zoom_flag = 1;
      }else if(strcmp(argv[i], "-zoomw") == 0) {
	zoom_flag = 2;
      }else if(strcmp(argv[i], "-zoom1") == 0) {
	if(zoom_flag == 0)
	  zoom_flag = 1;
	zoom_flag1 = 1;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-onpulsegr") == 0) {
	selectMoreOnpulseRegions = 1; 
      }else if(strcmp(argv[i], "-p3fold") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %d", &p3_fold, &p3_fold_nbin, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	p3_fold_flag = 1;
	i++;
      }else if(strcmp(argv[i], "-p3fold_nritt") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &p3_fold_refine, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-p3fold_cpb") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &p3_fold_cpb, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-p3fold_smooth") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &p3_fold_smoothWidth, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-p3fold_nosmooth") == 0) {
	p3fold_nosmooth = 1;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-p3fold_norefine") == 0) {
	p3_fold_refine = 0;
      }else if(strcmp(argv[i], "-p3fold_noonpulse") == 0) {
	p3_fold_onpulse_flag = 0;
      }else if(strcmp(argv[i], "-p3fold_dphase") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &p3fold_dphase, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-p3fold_writephases") == 0) {
	readwriteP3foldOffsets = 1;
      }else if(strcmp(argv[i], "-p3fold_readphases") == 0) {
	readwriteP3foldOffsets = 2;
	p3_fold_inputfile = i+1;
	i++;
      }else if(strcmp(argv[i], "-p3fold_writepulselist") == 0) {
	writeP3foldPulseList = 1;
	p3_fold_pulse_list_file = i+1;
	i++;
      } else if(strcmp(argv[i], "-carousel") == 0) {			// plot a carousel			
	j = sscanf(argv[i+1], "%f %f %d %f %d %d %lf", &plist.alpha, &plist.beta, &plist.SIZE, &plist.R_max, &plist.first_pulse, &plist.last_pulse, &plist.P3_hat);
	if(j != 7) {
	  printerror(application.verbose_state.debug, "ERROR: Error parsing option '%s'", argv[i]);
	  return 0;
	}
	plist.car = 1;
	i++;							
      } else if(strcmp(argv[i], "-p3carousel") == 0) {			// plot a carousel for P3-folded data			
	j = sscanf(argv[i+1], "%f %f %d %f %d", &plist.alpha, &plist.beta, &plist.SIZE, &plist.R_max, &plist.N);
	if(j != 5) {
	  printerror(application.verbose_state.debug, "ERROR: Error parsing option '%s'", argv[i]);
	  return 0;
	}
	plist.p3car = 1;
	i++;							
      } else if(strcmp(argv[i], "-abrot") == 0) { // Flag for "abnormal" carousel rotation (A = +1, B = -1). Default is "normal" (A = B = +1)
	plist.B = -1;
      } else if(strcmp(argv[i], "-phi0") == 0) {			// Input for phi0
	j = sscanf(argv[i+1], "%f", &plist.phi0);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "ERROR: Error parsing option '%s'", argv[i]);
	  return 0;
	}
	i++;
	

	
//START REGION RELEASE
      }else if(strcmp(argv[i], "-slope") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &slope, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pfold: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else {
	/* If the option is not recognized, assume it is a filename */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "ERROR pfold: Unknown option: %s\nRun 'pfold' without command line arguments to show help", argv[i]);
	  return 0;
	}
	else {
	  if(applicationAddFilename(i, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pfold: applicationAddFilename() failed", argv[i]);
	    terminateApplication(&application);
	    return 0;
	  }
	}
      }
    }
  }
//START REGION DEVELOP

  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
 
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) != 1) {
    printerror(application.verbose_state.debug, "ERROR: Please specify one and only one input file.\nCurrently there are %d input files recognized on the command line", numberInApplicationFilenameList(&application, argv, application.verbose_state));
    return 0; // Kill the program if there are more than one filenames
  }


  if(application.ppgplot_name != NULL) {
    if(pgopenoutputfile(application.ppgplot_name) == 0) {
      printerror(application.verbose_state.debug, "ERROR pfold: Cannot open %s", application.ppgplot_name);
      return 0;
    }
  }

//START REGION RELEASE
 
  /* Open data and read in profile */
  char *filename_ptr;
  filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state);
  if(filename_ptr == NULL) {
    printerror(application.verbose_state.debug, "ERROR pfold: Cannot obtain next file name to process.");
    return 0;
  }
  for(i = 0; i < MaxNrPolarizations; i++)
    cleanPSRData(&fin[i], application.verbose_state);
  if(application.iformat <= 0) {
    application.iformat = guessPSRData_format(argv[argc-1], 0, application.verbose_state);
    if(application.iformat == -2 || application.iformat == -3)
      return 0;
  }
  if(isValidPSRDATA_format(application.iformat) == 0) {
    printerror(application.verbose_state.debug, "ERROR pfold: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", filename_ptr);
    return 0;
  }
  /* Load data file into fin[0], the data will be split later into the separate polarizations */
  closePSRData(&fin[0], 0, application.verbose_state);  // When reading in only a data-set, make sure it is not allocated
  if(!openPSRData(&fin[0], argv[argc-1], application.iformat, 0, 1, 0, application.verbose_state))
    return 0;
  /* Search commandline for header parameters to overwrite header
     parameters. */
  if(PSRDataHeader_parse_commandline(&fin[0], argc, argv, application.verbose_state) == 0)
    return 0;
  for(i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-header") == 0) {
      printwarning(application.verbose_state.debug, "WARNING pfold: If using the -header option, be aware it applied BEFORE the preprocessing.");
    }
  }
  if(preprocessApplication(&application, &fin[0]) == 0) {
    return 0;
  }

  // Can check for period here if required

  double period;
  int ret;
  ret = get_period(fin[0], 0, &period, application.verbose_state);
  if(ret == 2) {
    printerror(application.verbose_state.debug, "ERROR pfold (%s): Cannot obtain period", fin[0].filename);
    return 0;
  }
  if(period < 0.001) {
    printerror(application.verbose_state.debug, "ERROR pfold: The period does not appear to be set in the header. Consider using the -header option.");
    return 0;
  }
  if(get_tsamp(fin[0], 0, application.verbose_state) < 0.0000001 || get_tsamp(fin[0], 0, application.verbose_state) >= 10) {
    printerror(application.verbose_state.debug, "ERROR pfold: The sampling time does not appear to be set correctly in the header. Consider using the -header option.");
    return 0;
  }

  if(fin[0].isDebase == 0) {
    printerror(application.verbose_state.debug, "ERROR pfold: Baseline is not subtracted. Use pmod -debase first.");
    return 0;
  }else if(fin[0].isDebase != 1) {
    printwarning(application.verbose_state.debug, "WARNING pfold:  It is not known if baseline is already subtracted. Use pmod -debase first.");
  }
  if(check_baseline_subtracted(fin[0], application.verbose_state) == 0) {
    printwarning(application.verbose_state.debug, "WARNING pfold: Baseline does not appear to be subtracted. Use pmod -debase first.");
  }
  /* Convert the specified onpulse regions which were defined as
     fractions on the commandline. This is now possible because we
     know the number of bins. */
  region_frac_to_int(&(application.onpulse), fin[0].NrBins, 0);

  /* Split the data in individual polarization channels, as required for p3 folding */
  if(p3_fold_flag) {
    if(fin[0].NrPols > MaxNrPolarizations) {
      printerror(application.verbose_state.debug, "ERROR pfold: Maximum supported input parameters is exceeded.\n");
      return 0;
    }
    originalNrPols = fin[0].NrPols;
    //  printf("XXXXXX Before selecting polarization 0: %p\n", fin[4].tsub_list);
    if(fin[0].NrPols > 0) {
      for(i = fin[0].NrPols-1; i >= 0; i--) {
	if(preprocess_polselect(fin[0], &clone, i, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pfold: Error selecting Stokes parameter %ld.", i);
	  return 0;
	}
	swap_orig_clone(&fin[i], &clone, application.verbose_state); 
      }
    }
    //  printf("XXXXXX After selecting polarization 0: %p\n", fin[4].tsub_list);
  }

  /* Copy header parameters, which is usefull when we write out the data. */
  cleanPSRData(&fout, application.verbose_state);
  copy_params_PSRData(fin[0], &fout, application.verbose_state);


  /* Do selection of the onpulse/offpulse regions */
  if(p3_fold_flag) {
    profileI = (float *)malloc(fin[0].NrBins*sizeof(float));
    if(profileI == NULL) {
      printerror(application.verbose_state.debug, "ERROR pfold: Cannot allocate memory");
      return 0;
    }
    read_profilePSRData(fin[0], profileI, NULL, 0, noverbose);
    
    xmin = 0;
    ret = get_period(fin[0], 0, &period, application.verbose_state);
    if(ret == 2) {
      printerror(application.verbose_state.debug, "ERROR pfold (%s): Cannot obtain period", fin[0].filename);
      return 0;
    }
    xmax = 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period;
    if(application.onpulse.nrRegions == 0 || selectMoreOnpulseRegions) {
      strcpy(pgplot_options.viewport.plotDevice, onpulseselectdevice);
      strcpy(pgplot_options.box.xlabel, "Bin");
      strcpy(pgplot_options.box.ylabel, "Intensity");
      strcpy(pgplot_options.box.title, "Select on-pulse region");
      selectRegions(profileI, fin[0].NrBins, &pgplot_options, 0, 0, 1, &application.onpulse, application.verbose_state);
    }else {
      if(strcmp(onpulseselectdevice, "?") == 0)
	printf("Specify plotting device to show the profile showing the selected regions: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, onpulseselectdevice);
      strcpy(pgplot_options.box.xlabel, "Phase[deg]");
      strcpy(pgplot_options.box.ylabel, "Intensity");
      strcpy(pgplot_options.box.title, fin[0].psrname);
      if(pgplotGraph1(&pgplot_options, profileI, NULL, NULL, fin[0].NrBins, xmin, xmax, 0, xmin, xmax, 0, 0, 0, 1, 0, 1, 0, 1, 1, &application.onpulse, -1, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pfold: Unable to open plotdevice.\n");
	return 0;
      }
    }

    /* Convert the specified onpulse regions which were defined as
       fractions on the commandline. This is now possible because we
       know the number of bins. */
    region_int_to_frac(&(application.onpulse), 1.0/(float)fin[0].NrBins, 0);
    regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);
    xmin_zoom = xmin;
    xmax_zoom = xmax;
    if(zoom_flag) {
      if(application.onpulse.nrRegions > 0) {
	if(zoom_flag1) {
	  if(application.onpulse.bins_defined[0] == 0) {
	    printerror(application.verbose_state.debug, "ERROR pfold: region not defined in bins");
	    return 0;
	  }
	  xmin_zoom = application.onpulse.left_bin[0];
	  xmax_zoom = application.onpulse.right_bin[0];
	}else {
	  xmin_zoom = 0;
	  for(i = 0; i < fin[0].NrBins; i++) {
	    if(checkRegions(i, &application.onpulse, 0, application.verbose_state) != 0) {
	      xmin_zoom = i;
	      break;
	    }
	  }
	  xmax_zoom = fin[0].NrBins-1;
	  for(i = fin[0].NrBins-1; i >= 0; i--) {
	    if(checkRegions(i, &application.onpulse, 0, application.verbose_state) != 0) {
	      xmax_zoom = i;
	      break;
	    }
	  }
	}
	if(zoom_flag == 2) {
	  i = xmax_zoom - xmin_zoom;
	  xmin_zoom -= i;
	  xmax_zoom += i;
	  if(xmin_zoom < 0)
	    xmin_zoom = 0;
	  if(xmax_zoom >= fin[0].NrBins)
	    xmax_zoom = fin[0].NrBins-1;
	}
	xmin_zoom = (xmax-xmin)*xmin_zoom/(float)fin[0].NrBins+xmin;
	xmax_zoom = (xmax-xmin)*xmax_zoom/(float)fin[0].NrBins+xmin;
      }
    }
  }


//START REGION RELEASE
  if(p3_fold_flag) {
    p3foldmap = malloc(originalNrPols*fin[0].NrBins * p3_fold_nbin * sizeof(float));
    if(p3foldmap == NULL) {
      printerror(application.verbose_state.debug, "ERROR pfold: Cannot allocate memory");
      return 0;
    }

//START REGION DEVELOP
    if(writeP3foldPulseList) {
      if(readwriteP3foldOffsets != 2) {
	printerror(application.verbose_state.debug, "ERROR pfold: The -p3fold_writepulselist option requires the use of -p3fold_readphases.");
	return 0;	
      }
    }
    if(readwriteP3foldOffsets == 2) {
      strncpy(outputname, argv[p3_fold_inputfile], 999);
    }else {
      if(change_filename_extension(argv[argc-1], outputname, "p3foldphases", 1000, application.verbose_state) == 0) {
	return 0;
      }
    }
//START REGION RELEASE
    originalNrPolsP3 = originalNrPols;
//START REGION DEVELOP
    if(readwriteP3foldOffsets != 2 && originalNrPols > 1) {
      printwarning(application.verbose_state.debug, "WARNING pfold: Only first polarization chanel is folded. You could use -p3fold_readphases in combination with the -polselect option.");
      originalNrPolsP3 = 1;
    }else if(1 == 0) {
//START REGION RELEASE
      if(originalNrPols > 1) {
	printwarning(application.verbose_state.debug, "WARNING pfold: Only first polarization chanel is folded.");
	originalNrPolsP3 = 1;
      }
//START REGION DEVELOP
    }
//START REGION RELEASE
//    printf("XXXXXX Before P3 fold: %p\n", fin[4].tsub_list);
    for(i = 0; i < originalNrPolsP3; i++) {
      if(p3_fold_onpulse_flag) {
	if(foldP3(fin[i].data, fin[i].NrSubints, fin[i].NrBins, &p3foldmap[i*fin[0].NrBins * p3_fold_nbin], p3_fold_nbin,  p3_fold, p3_fold_refine, p3_fold_cpb, p3fold_nosmooth, p3_fold_smoothWidth, slope*360.0/(float)fin[i].NrBins, p3fold_dphase, &application.onpulse
//START REGION DEVELOP
		  , outputname, readwriteP3foldOffsets, argv[p3_fold_pulse_list_file], writeP3foldPulseList
//START REGION RELEASE
, application.verbose_state) == 0) {
	  return 0;
	}
      }else {
	if(foldP3(fin[i].data, fin[i].NrSubints, fin[i].NrBins, &p3foldmap[i*fin[0].NrBins * p3_fold_nbin], p3_fold_nbin,  p3_fold, p3_fold_refine, p3_fold_cpb, p3fold_nosmooth, p3_fold_smoothWidth, slope*360.0/(float)fin[i].NrBins, p3fold_dphase, NULL
//START REGION DEVELOP
		  , outputname, readwriteP3foldOffsets, argv[p3_fold_pulse_list_file], writeP3foldPulseList
//START REGION RELEASE
, application.verbose_state) == 0) {
	  return 0;
	}
      }
    }
    if(strcmp(p3fold_device, "?") == 0)
      printf("Specify plotting device to show the P3 fold map: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, p3fold_device);
      strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.ylabel, "P3 [pulse periods]");
      strcpy(pgplot_options.box.title, "P3 fold");
      ret = get_period(fin[0], 0, &period, application.verbose_state);
      if(ret == 2) {
	printerror(application.verbose_state.debug, "ERROR pfold (%s): Cannot obtain period", fin[0].filename);
	return 0;
      }
      pgplotMap(&pgplot_options, p3foldmap, fin[0].NrBins, p3_fold_nbin, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period, xmin_zoom, xmax_zoom, 0.5*p3_fold/(float)p3_fold_nbin, 0.5*p3_fold/(float)p3_fold_nbin + p3_fold*(p3_fold_nbin-1)/(float)p3_fold_nbin, 0, p3_fold, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, application.verbose_state);
    if(write_flag) {
      fout.NrSubints = p3_fold_nbin;
      fout.NrBins = fin[0].NrBins;
      fout.NrPols = originalNrPolsP3;
      fout.gentype = GENTYPE_P3FOLD;
      fout.tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
      if(fout.tsub_list != NULL)
	free(fout.tsub_list);
      fout.tsub_list = (double *)malloc(sizeof(double));
      if(fout.tsub_list == NULL) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR pfold: Memory allocation error");
	return 0;
      }
      fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
      fout.yrangeset = 1;
      fout.yrange[0] = 0.5*p3_fold/(float)p3_fold_nbin;
      fout.yrange[1] = 0.5*p3_fold/(float)p3_fold_nbin + p3_fold*(p3_fold_nbin-1)/(float)p3_fold_nbin;

      if(change_filename_extension(argv[argc-1], outputname, "p3fold", 1000, application.verbose_state) == 0) {
	return 0;
      }
      if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state)) {
	printerror(application.verbose_state.debug, "ERROR pfold: Unable to open file for writing.\n");
	return 0;
      }
      if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pfold: Unable to write header.\n");
	return 0;
      }
      //      appendHistoryLine(fout, argc, argv, application.verbose_state);

      /* Transpose the data */
      p3foldmap2 = malloc(originalNrPols*fin[0].NrBins * p3_fold_nbin * sizeof(float));
      if(p3foldmap2 == NULL) {
	printerror(application.verbose_state.debug, "ERROR pfold: Cannot allocate memory");
	return 0;
      }
      for(i = 0; i < p3_fold_nbin; i++) {
	for(p = 0; p < originalNrPolsP3; p++) {
	  memcpy(&p3foldmap2[(originalNrPolsP3*i+p)*fin[0].NrBins], &p3foldmap[(p*p3_fold_nbin+i)*fin[0].NrBins], fin[0].NrBins*sizeof(float));
	}
      }
      if(writePSRData(&fout, p3foldmap2, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pfold: Unable to write data.\n");
	return 0;
      }
      free(p3foldmap2);
      closePSRData(&fout, 1, application.verbose_state); // Keep history, as header can be written out again
      fout.gentype = GENTYPE_UNDEFINED;
      fout.yrangeset = 0;
      fout.xrangeset = 0;
      fout.NrPols = 1;
    }
    
    free(p3foldmap);
  }  // End of if(p3_fold_flag)

//START REGION DEVELOP
  if(plist.car == 1 || plist.p3car == 1) {
    // Set debugging flags appropriately
    plist.verbose = application.verbose_state.verbose;
    plist.debug = application.verbose_state.debug;
    plist.app_pointer = &application;

    // read in various parameters from the input data file 
    // Maybe stick this in a function at some point?
    plist.NrBins = fin[0].NrBins;
    plist.NrSubints = fin[0].NrSubints;
    plist.NrPols = fin[0].NrPols;
    plist.NrFreqs = fin[0].NrFreqChan;

    // Only automatically assign these if there is no user input
    ret = get_period(fin[0], 0, &(plist.P1), application.verbose_state);
    if(ret == 2) {
      printerror(application.verbose_state.debug, "ERROR pfold (%s): Cannot obtain period", fin[0].filename);
      return 0;
    }
    if(plist.verbose) printf("  Reading P1 from input file.\n");

    plist.sampling_period = get_tsamp(fin[0], 0, application.verbose_state);


    if(plist.p3car == 1) {
      if(fin[0].gentype == GENTYPE_UNDEFINED) {
	printwarning(application.verbose_state.debug, "WARNING pfold: The type of data is undefined. It will be assumed this is a P3 fold.");
      }else if(fin[0].gentype != GENTYPE_P3FOLD) {
	printerror(application.verbose_state.debug, "WARNING pfold: The input file does not appear to be a P3 fold.");
	return 0;
      }
      if(fin[0].yrangeset == 0) {
	printerror(application.verbose_state.debug, "WARNING pfold: The yrange variable does not appear to be set in the input data, which should not happen for a P3 fold.");
	return 0;
      }
      plist.P3 = plist.NrSubints * fin[0].yrange[0] * 2.0; // Calculate P3 from the input data
      plist.P3_hat = plist.P3 * (float)plist.N;
      printf("P3 = %.3f\n", plist.P3);
    }else if(plist.car == 1) {
      if(fin[0].gentype == GENTYPE_UNDEFINED) {
	printwarning(application.verbose_state.debug, "WARNING pfold: The type of data is undefined. It will be assumed this is a pulse stack.");
      }else if(fin[0].gentype == GENTYPE_SUBINTEGRATIONS) {
	printwarning(application.verbose_state.debug, "WARNING pfold: It will be assumed that every subintegration is an individual pulse.");
      }else if(fin[0].gentype != GENTYPE_PULSESTACK) {
	printerror(application.verbose_state.debug, "WARNING pfold: The input file does not appear to be a P3 fold.");
	return 0;
      }
    }
 
    if ((plist.first_pulse != -1) && (plist.last_pulse != -1)) {
      plist.num_pulses = plist.last_pulse - plist.first_pulse + 1;
      plist.NrSubints = plist.num_pulses;
    } else { 
      plist.first_pulse = 0;
      plist.last_pulse = plist.NrSubints - 1;
      plist.num_pulses = plist.NrSubints;
    }
   
    // Function to display parameters if the flag -p was included
    display_parameters(plist, application.verbose_state);

    // Generate an array of size (SIZE)x(SIZE)x(NrPols)x(NrFreqs) to store the output data
    float *output_image;
    output_image = (float *)malloc(sizeof(float)*plist.SIZE*plist.SIZE*plist.NrPols*plist.NrFreqs);
    if (output_image == NULL) {
      printerror(plist.debug, "ERROR: Could not assign memory for the output image.");
      return 0;
    }

    if (plist.car == 1) { // Carousel Plotting(float)plist.SIZE)
      printf("\nRendering carousel for %s.\n", filename_ptr);   
      if(plist.debug) printf("    [%d] Called plot_carousel\n", __LINE__+1);
      plotCarousel(plist, filename_ptr, output_image, fin[0]);
      
    } else if (plist.p3car == 1) { // Carousel plotting of P3-folded data
      printf("\nRendering carousel from P3-folded data (input = %s).\n", filename_ptr);
      ret = plotCarousel(plist, filename_ptr, output_image, fin[0]);
      if(plist.debug) {
	printf("PlotCarousel returned %d\n", ret);
	//	for(i = 0; i < plist.SIZE*plist.SIZE; i++) {
	//	  printf("XXXXX %f\n", output_image[i]);
	//	}
      }
    } else if (plist.stack == 1) { // Stack Plotting
      // Write a new function to plot a pulse stack, call it plot_stack or something
    } else if (plist.movie == 1) {	// Movie Plotting
      // Removed this functionality whilst the rest of the code is being streamlined -rewrite later!
    } else {
      printerror(plist.debug, "ERROR pfold: Something seriously wrong with the code.");
      return 0;
    }

    if(write_flag) {
      printf("Writing carousel to file.\n");

      fout.NrBins = plist.SIZE;
      fout.NrSubints = plist.SIZE;
      fout.NrFreqChan = plist.NrFreqs;
      fout.NrPols = plist.NrPols;
      fout.xrangeset = 1;
      fout.yrangeset = 1;
      fout.xrange[0] = -plist.R_max;
      fout.xrange[1] = plist.R_max;
      fout.yrange[0] = -plist.R_max;
      fout.yrange[1] = plist.R_max;
      fout.gentype = GENTYPE_POLARMAP;
      fout.tsubMode = TSUBMODE_FIXEDTSUB;

      if(fout.tsub_list != NULL) {
	free(fout.tsub_list);
      }
      fout.tsub_list = malloc(sizeof(double));
      if (fout.tsub_list == NULL) {
	printerror(plist.debug, "ERROR pfold: Could not assign memory for fout.tsub_list.");
      }
      fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state)/plist.SIZE; // Note: given the number of "subintegrations" the reported observation duration should still be correctly reported.


      // Change the input filename and replaces the extension with "carousel
      char outputname[1000];
      if(change_filename_extension(filename_ptr, outputname, "carousel", 1000, application.verbose_state) == 0) {
	return 0;
      }

      // Open an actual file, enabling writing, on the harddisk with the given name and the given format
      if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state)) {return 0;}
      // Write out the header information
      if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pfold: Unable to write header.\n");
	return 0;
      }
      // Writes out all data which was previously already was put in memory to disk
      if(writePSRData(&fout, output_image, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pfold: Unable to write data.\n");
	return 0;
      }


    }
    free(output_image);
  } // End of carousel stuff
//START REGION RELEASE


  closePSRData(&fout, 0, application.verbose_state);  // Release header 
  for(i = 0; i < MaxNrPolarizations; i++)
    closePSRData(&fin[i], 0, application.verbose_state);
  if(p3_fold_flag) {
    free(profileI);
  }

  ppgend();
//START REGION DEVELOP
  if(application.ppgplot_name != NULL) {
    pgcloseoutputfile();
  }

//START REGION RELEASE

  terminateApplication(&application);
  return 0;
}


//START REGION DEVELOP
void display_parameters(carousel_parameter_list plist, verbose_definition verbose) {
  if (verbose.verbose) {
    // This needs updating a fair amount...
    printf("\nCurrent parameters:\n==============================\n");
    printf("alpha	= %f degrees\n", plist.alpha);
    printf("beta	= %f degrees\n", plist.beta);
    printf("A	= %d\n", plist.A);
    printf("B	= %d\n", plist.B);
    printf("k0	= %d\n", plist.k0);
    printf("phi0	= %f degrees\n", plist.phi0);
    printf("P1	= %f seconds\n", plist.P1);
    printf("P3_hat	= %f x P1\n", plist.P3_hat);
    printf("R_max	= %f degrees\n", plist.R_max);
    printf("SIZE	= %d pixels\n", plist.SIZE);
    printf("first_pulse = %d\n", plist.first_pulse);
    printf("last_pulse = %d\n", plist.last_pulse);
    printf("sampling_period = %f seconds\n\n", plist.sampling_period);
    //    printf("pulsepf = %d\n", plist.pulsepf);
    //    printf("frame_limiter = %d\n", plist.frame_limiter);
  }
}

//START REGION RELEASE
