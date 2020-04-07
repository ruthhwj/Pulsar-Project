//START REGION DEVELOP
//START REGION RELEASE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_randist.h"
#include "psrsalsa.h"

// If isnormal is not defined on your system, uncomment the following
// line to disable this test of the results.
// #define isnormal(x) ( 1 )

//START REGION DEVELOP
void SHOWREVISIONINFO_prog() {
#include "rmsynth.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);
//START REGION RELEASE

double get_FWHM_RMsynthesis(datafile_definition datain, verbose_definition verbose);

int main(int argc, char **argv)
{
  int i, j, index, nrrmsteps, status, parabola, pgplot_mapdevice, pgplot_rmdevice, pgplot_profiledevice;
  int collapse, write_ascii, bootstrap, bootstrap_itt, nrOnpulseBins, not_normal_result_warning;
//START REGION DEVELOP
  int write_map, showPAprofile;
  float offpulsedata_rebinfactor;
  char offpulsedata_name[MaxFilenameLength];
  datafile_definition data_offpulse;
//START REGION RELEASE
  long binnr, freqchannelnr, polnr, idnum;
  float rmmin, rmmax, sample;
  float *rmsynth_array, *cmap, *rms_channels;  
  float *singlespectrum;
  double *singlespectrum_double, *rmgrid, *sigmagrid, *rmestimate, *rmcurestimate;
  double *rm_av, *rm_square, rmsigma, expectedRMerror, ftol;
  char outputname[MaxFilenameLength];
  FILE *ofile;
  psrsalsaApplication application;
  datafile_definition datain, clone;
  pgplot_options_definition pgplot_options;
  fitfunc_collection_type function;
  datafile_definition profiledata;
  gsl_rng *rand_num_gen;
  const gsl_rng_type *rand_num_gen_type;

  initApplication(&application, "rmsynth", "[options] inputfile");
//START REGION DEVELOP
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
//START REGION RELEASE

  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_rebin = 1;
  application.switch_rot = 1;
  application.switch_onpulse = 1;
  application.switch_onpulsef = 1;
  application.switch_onpulse2 = 1;
  application.switch_onpulsef2 = 1;
  application.switch_nocounters = 1;
  //  application.switch_debase = 1;    -- Not allowed because of complication that -onpulse2 should be used, which is defined by user after the pre-process options are applied. But the preprocess options should be applied first in case the profile is affected.
  application.switch_stokes = 1;
  application.oformat = FITS_format;
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_dedisperse = 1;
  application.switch_fixseed = 1;
//START REGION DEVELOP
  application.switch_history_cmd_only = 1;
//START REGION RELEASE

  nrrmsteps = 0;
  parabola = 0;
  pgplot_mapdevice = 0;
  pgplot_rmdevice = 0;
  pgplot_profiledevice = 0;
  collapse = 1;
  write_ascii = 0;
  bootstrap = 0;
  not_normal_result_warning = 0;
  ftol = -1;  // Means take the default values, which depends on the fitting algorithm used
//START REGION DEVELOP
  write_map = 0;
  showPAprofile = 0;
  offpulsedata_name[0] = 0;
//START REGION RELEASE

  if(argc < 2) {
    printApplicationHelp(&application);
    fprintf(stdout, "Other options:\n");
    fprintf(stdout, "  -device       pgplot device for the Faraday depth plot.\n");
    fprintf(stdout, "  -device2      pgplot device for the longitude resolved map.\n");
    fprintf(stdout, "  -device3      pgplot device for the profile (used with -bootstrap).\n");
    fprintf(stdout, "  -rm           \"min max steps\" Set the minimum and maximum rm to search over and the nr of rm steps to use.\n");
    fprintf(stdout, "  -parabola     Fit data with parabola rather than estimated reponds.\n");
    fprintf(stdout, "  -long         Get a RM for each bin independently.\n");
    fprintf(stdout, "  -ascii        An ascii file the found RM values is generated.\n");
//START REGION DEVELOP
    fprintf(stdout, "  -asciispec    Adds the RM spectrum to the output of -ascii.\n");
    fprintf(stdout, "  -writemap     The RM map is written out to a file.\n");
//START REGION RELEASE
    fprintf(stdout, "  -bootstrap    Bootstrap the errorbars on the RM this number of times.\n");
//START REGION DEVELOP
    fprintf(stdout, "  -showprofile  Show de-Faraday rotated profile.\n");
    fprintf(stdout, "  -offpulsedata \"file rebinfactor\"\n");
    fprintf(stdout, "                The specified file will be used to determine the noise\n");
    fprintf(stdout, "                characteristics. This can be useful if the actual data to be\n");
    fprintf(stdout, "                used to measure the RM is re-binned, thereby reducing the number\n");
    fprintf(stdout, "                of offpulse bins, hence reducing the precision of the rms.\n");
    fprintf(stdout, "                Set rebinfactor to the ratio of the sampling times of the actual\n");
    fprintf(stdout, "                data and the offpulse data file (>1 if actual data is rebinned).\n");
    fprintf(stdout, "                Any selection of the onpulse (-onpulse2) region should be related\n");
    fprintf(stdout, "                to the file to be used for the noise calculation. It is assumed\n");
    fprintf(stdout, "                that the rebinning has been done in such a way that the mean\n");
    fprintf(stdout, "                remains the same, as is the case with pmod -rebin.\n");
//START REGION RELEASE
    fprintf(stdout, "  -ftol         Set the fractional tolerance of the fitting algorithm used.\n");
//START REGION DEVELOP
    fprintf(stdout, "\n\nNote: specifying an onpulse region on the command line can greatly improve the efficiency.\n");
//START REGION RELEASE
    printf("\nPlease use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about RM syntesis method as implemented here can be found in: Ilie et al. 2018, accepted for publication in MNRAS, astro-ph/1811.12831\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-rm") == 0) {
	j = sscanf(argv[i+1], "%f %f %d", &rmmin, &rmmax, &nrrmsteps);
	if(j != 3) {
	  printerror(application.verbose_state.debug, "Error parsing %s option", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-bootstrap") == 0) {
	j = sscanf(argv[i+1], "%d", &bootstrap);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "Error parsing %s option", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-ftol") == 0) {
	j = sscanf(argv[i+1], "%lf", &ftol);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "Error parsing %s option", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-parabola") == 0) {
	parabola = 1;
      }else if(strcmp(argv[i], "-long") == 0) {
	collapse = 0;
      }else if(strcmp(argv[i], "-ascii") == 0) {
	write_ascii = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-showprofile") == 0) {
	showPAprofile = 1;
      }else if(strcmp(argv[i], "-asciispec") == 0) {
	write_ascii = 2;
      }else if(strcmp(argv[i], "-writemap") == 0) {
	write_map = 1;
      }else if(strcmp(argv[i], "-offpulsedata") == 0) {
	int ret;
	ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%1000s %f", offpulsedata_name, &offpulsedata_rebinfactor, NULL);
	if(ret != 2) {
	  printerror(application.verbose_state.debug, "ERROR rmsynth: Error parsing '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-device") == 0) {
	pgplot_rmdevice = i+1;
	i++;
      }else if(strcmp(argv[i], "-device2") == 0) {
	pgplot_mapdevice = i+1;
	i++;
      }else if(strcmp(argv[i], "-device3") == 0) {
	pgplot_profiledevice = i+1;
	i++;
      }else {
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "rmsynth: Unknown option: %s\n\nRun rmsynth without command line arguments to show help", argv[i]);
	  terminateApplication(&application);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }

  if(application.dodebase) {
    application.dodebase = 2;  // To ensure -onpulse2 is used rather than -onpulse
  }

  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "rmsynth: No files specified");
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) > 1) {
    printerror(application.verbose_state.debug, "rmsynth: Does not support more than one input file");
    return 0;
  }

  if(nrrmsteps == 0) {
    printerror(application.verbose_state.debug, "Use the -rm option to specify range of rm's to explore");
    return 0;
  }
  if(nrrmsteps <= 1) {
    printerror(application.verbose_state.debug, "Use the -rm option to specify more than one rm to explore");
    return 0;
  }

  gsl_rng_env_setup();    /* Set the default generators, can be influenced by environment variables */
  rand_num_gen_type = gsl_rng_default;
  rand_num_gen = gsl_rng_alloc(rand_num_gen_type);
  if(application.fixseed)
    idnum = 1;
  else
    randomize_idnum(&idnum);
  gsl_rng_set(rand_num_gen, idnum);

//START REGION DEVELOP
  // Read in separate data to be used for the offpulse statistics, if requested.
  if(offpulsedata_name[0] != 0) {
    if(!openPSRData(&data_offpulse, offpulsedata_name, application.iformat, 0, 1, 0, application.verbose_state))
      return 0;
    if(application.verbose_state.verbose) {
      fflush(stdout);
      printwarning(application.verbose_state.debug, "Offpulse input data contains %ld bins, %ld pulses, %ld polarizations and %ld frequencies.", (data_offpulse.NrBins), data_offpulse.NrSubints, (data_offpulse.NrPols), data_offpulse.NrFreqChan);
    }

    if(data_offpulse.isDeFarad == -1) {
      printwarning(application.verbose_state.debug, "WARNING rmsynth: De-Faraday rotation state is unknown. Assume the data is not yet de-Faraday rotated.\n");
      datain.isDeFarad = 0;
    }

    if(data_offpulse.isDebase == 0) {
      printerror(application.verbose_state.debug, "ERROR rmsynth: Baseline is not subtracted. Use pmod -debase first.\n");
      return 0;
    }else if(data_offpulse.isDebase != 1) {
      printwarning(application.verbose_state.debug, "WARNING rmsynth: It is not known if baseline is already subtracted, which might affect the fitting.\n");
    }
    if(check_baseline_subtracted(data_offpulse, application.verbose_state) == 0) {
      printwarning(application.verbose_state.debug, "WARNING rmsynth: Baseline does not appear to be subtracted, which might affect the fitting.\n");
    }

  }
//START REGION RELEASE


  /* Read in data */
  //  cleanPSRData(&datain, application.verbose_state);
  if(!openPSRData(&datain, argv[argc-1], 0, 0, 1, 0, application.verbose_state))
    return 0;

  /* Search commandline for header parameters to overwrite header
     parameters. */
  if(PSRDataHeader_parse_commandline(&datain, argc, argv, application.verbose_state) == 0)
    return 0;

  if(datain.isDeFarad == -1) {
    printwarning(application.verbose_state.debug, "WARNING rmsynth: De-Faraday rotation state is unknown. Assume the data is not yet de-Faraday rotated.\n");
    datain.isDeFarad = 0;
  }

  if(datain.isDebase == 0) {
    printerror(application.verbose_state.debug, "ERROR rmsynth: Baseline is not subtracted. Use pmod -debase first.\n");
    return 0;
  }else if(datain.isDebase != 1) {
    printwarning(application.verbose_state.debug, "WARNING rmsynth: It is not known if baseline is already subtracted, which might affect the fitting.\n");
  }
  if(check_baseline_subtracted(datain, application.verbose_state) == 0) {
    printwarning(application.verbose_state.debug, "WARNING rmsynth: Baseline does not appear to be subtracted, which might affect the fitting.\n");
  }

     

  /* Convert the specified onpulse regions which were defined as
     fractions on the commandline. This is now possible because we
     know the number of bins. */
  region_frac_to_int(&(application.onpulse), datain.NrBins, 0);

//START REGION DEVELOP
  if(offpulsedata_name[0] != 0) {
    region_frac_to_int(&(application.onpulse2), data_offpulse.NrBins, 0);
  }else {
//START REGION RELEASE
    region_frac_to_int(&(application.onpulse2), datain.NrBins, 0);
//START REGION DEVELOP
  }
//START REGION RELEASE


  for(i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-header") == 0) {
      printwarning(application.verbose_state.debug, "WARNING: If using the -header option, be aware it applied BEFORE the preprocessing.");
    }
  }
  // See if there are any preprocess commands on command line
  if(preprocessApplication(&application, &datain) == 0) {
    return 0;
  }

  // Take out dispersive delay.
  if(preprocess_dedisperse(&datain, 0, 0, 0, application.verbose_state) == 0)
    return 0;
//START REGION DEVELOP
  if(offpulsedata_name[0] != 0) {
    if(preprocess_dedisperse(&data_offpulse, 0, 0, 0, application.verbose_state) == 0)
      return 0;
  }
//START REGION RELEASE

//START REGION DEVELOP
  // Make a frequency scrunched clone to show as a pulse profile
  if(offpulsedata_name[0] != 0) {
    if(!preprocess_addsuccessiveFreqChans(data_offpulse, &profiledata, datain.NrFreqChan, NULL, application.verbose_state))
      return 0;
  }else {
//START REGION RELEASE
    if(!preprocess_addsuccessiveFreqChans(datain, &profiledata, datain.NrFreqChan, NULL, application.verbose_state))
      return 0;
//START REGION DEVELOP
  }
//START REGION RELEASE

  pgplot_clear_options(&pgplot_options);
  if(pgplot_profiledevice)
    strcpy(pgplot_options.viewport.plotDevice, argv[pgplot_profiledevice]);
  else
    strcpy(pgplot_options.viewport.plotDevice, "?");
  strcpy(pgplot_options.box.ylabel, "Stokes I");
  strcpy(pgplot_options.box.title, "Select on-pulse region (for noise calculation only)");

  // Ask the user to select an onpulse region, if it is not specified on command line
  if(application.onpulse2.nrRegions == 0) {
    selectRegions(profiledata.data, profiledata.NrBins, &pgplot_options, 0, 0, 0, &(application.onpulse2), application.verbose_state);
  }else {
    pgplotGraph1(&pgplot_options, profiledata.data, NULL, NULL, profiledata.NrBins, 0, profiledata.NrBins-1, 0, 0, profiledata.NrBins-1, 0, 0, 0, 1, 0, 1, 0, 1, 1, &(application.onpulse2), -1, application.verbose_state);
  }
//START REGION DEVELOP
  if(offpulsedata_name[0] != 0) {
    region_int_to_frac(&(application.onpulse2), 1.0/(float)profiledata.NrBins, 0);
  }else {
//START REGION RELEASE
    region_int_to_frac(&(application.onpulse2), 1.0/(float)datain.NrBins, 0);
//START REGION DEVELOP
  }
//START REGION RELEASE
  regionShowNextTimeUse(application.onpulse2, "-onpulse2", "-onpulsef2", stdout);

//START REGION DEVELOP
  // After the offpulse selection has been made, make a new profile of the actual input data to be shown in the graph showing RM vs pulse longitude.
  if(offpulsedata_name[0] != 0) {
    closePSRData(&profiledata, 0, application.verbose_state);
    if(!preprocess_addsuccessiveFreqChans(datain, &profiledata, datain.NrFreqChan, NULL, application.verbose_state))
      return 0;
  }
//START REGION RELEASE

  if(bootstrap) {
    rms_channels = (float *)calloc(datain.NrFreqChan, sizeof(float));
    rm_av     = (double *)calloc(datain.NrBins, sizeof(double));
    rm_square = (double *)calloc(datain.NrBins, sizeof(double));
    if(rm_av  == NULL || rm_square  == NULL || rms_channels == NULL) {
      printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot allocate memory");
      return 0;
    }
    
    for(i = 0; i < datain.NrFreqChan; i++) {
//START REGION DEVELOP
      if(offpulsedata_name[0] != 0) {
	// Note we will assume that the rms'es are the same for all Stokes parameters, so only keep rms of Stokes I.
	if(read_rmsPSRData(data_offpulse, &rms_channels[i], NULL, NULL, &(application.onpulse2), 0, 0, i, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot determine RMS");
	  return 0;	
	}
	// Note that the rebinning is keeping the mean the same (as happens with pmod -rebin). This implies that the rms should
	// go down after rebinning.
	rms_channels[i] /= sqrt(offpulsedata_rebinfactor);
      }else {
//START REGION RELEASE
	if(read_rmsPSRData(datain, &rms_channels[i], NULL, NULL, &(application.onpulse2), 0, 0, i, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot determine RMS");
	  return 0;	
	}
//START REGION DEVELOP
      }
//START REGION RELEASE
      //      printf("RMS channel %d = %e\n", i, rms_channels[i]);
    }
  }

  // Make a clone to be used in bootstrapping
  cleanPSRData(&clone, application.verbose_state);
  if(copy_params_PSRData(datain, &clone, application.verbose_state) == 0)
    return 0;
  //  memmove(&clone, &datain, sizeof(datafile_definition));
  clone.data = (float *)malloc(datain.NrBins*datain.NrFreqChan*datain.NrPols*sizeof(float));
  if(clone.data == NULL) {
    printerror(application.verbose_state.debug, "ERROR rmsynth: Memory allocation error");
    return 0;	
  }

  // Makes sure the memory gets allocated at first itteration
  rmsynth_array = NULL;

  // Map to show in pgplot of RM versus pulse longitude
  cmap = (float *)malloc(nrrmsteps*datain.NrBins*sizeof(float));

  // The best found RMs as function of pulse longitude
  rmestimate = (double *)calloc(datain.NrBins, sizeof(double));
  rmcurestimate = (double *)calloc(datain.NrBins, sizeof(double));

  // The following array will contain the RM synthesis power as function of RM (in float and double format)
  singlespectrum = malloc(nrrmsteps*sizeof(float));
  singlespectrum_double = malloc(nrrmsteps*sizeof(double));
  
  // The following array will contain the error bars on the RM synthesis points. In fact, will all be set to 1.
  sigmagrid = malloc(nrrmsteps*sizeof(double));
  
  // The following array will contain the RM values corresponding to the singlespectrum array
  rmgrid = malloc(nrrmsteps*sizeof(double));
    
  // Check if allocation worked
  if(cmap == NULL || rmestimate == NULL || rmcurestimate == NULL || rmgrid == NULL || singlespectrum == NULL || singlespectrum_double == NULL || sigmagrid == NULL) {
    printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot allocate memory");
    return 0;
  }


  //Loop over the nr of bootstrap itterations to determine an errorbar
  for(bootstrap_itt = 0; bootstrap_itt <= bootstrap; bootstrap_itt++) {

    if(bootstrap_itt == 1 && application.verbose_state.verbose)
      printf("Starting bootstrap\n");

    if(application.verbose_state.debug)
      printf("  itteration %d: Making clone of data with added noise\n", bootstrap_itt+1);
    for(freqchannelnr = 0; freqchannelnr < datain.NrFreqChan; freqchannelnr++) {
      for(binnr = 0; binnr < datain.NrBins; binnr++) {
	for(polnr = 0; polnr < datain.NrPols; polnr++) {
	  if(readPulsePSRData(&datain, 0, polnr, freqchannelnr, binnr, 1, &sample, application.verbose_state) != 1) {
	    printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot read data.");
	    return 0;
	  }
	  //	  printf("XXXX sample = %e\n", sample);
	  //	  exit(0);
	  if(bootstrap_itt > 0)
	    sample += gsl_ran_gaussian(rand_num_gen, rms_channels[freqchannelnr]);
	  if(writePulsePSRData(&clone, 0, polnr, freqchannelnr, binnr, 1, &sample, application.verbose_state) != 1) {
	    printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot write data.");
	    return 0;
	  }
	}
      }
    }


    // Do the rm synthesis. Results will be stored in rmsynth_array,
    // which will be allocated in first itteration of bootstrap_itt
    verbose_definition verbose2;
    copyVerboseState(application.verbose_state, &verbose2);
    verbose2.verbose = application.verbose_state.verbose;
    if(bootstrap_itt != 0)
      verbose2.verbose = 0;
    if(application.verbose_state.debug)
      verbose2.verbose = 1;
    if(application.verbose_state.debug)
      printf("  itteration %d: Applying RM synthesis\n", bootstrap_itt+1);
    if(rmSynthesis(clone, rmmin, rmmax, &rmsynth_array, nrrmsteps, &(application.onpulse), verbose2) == 0)
      return 0;
    if((application.verbose_state.verbose && bootstrap_itt == 0) || application.verbose_state.debug)
      printf("RM synthesis calculation done\n");
    
  
    // Show plot of RM map in first itteration
    if(bootstrap_itt == 0) {
      // Reorder rmsynth_array in order to show a map with pgplot
      for(i = 0; i < datain.NrBins; i++) {
	for(j = 0; j < nrrmsteps; j++) {
	  cmap[j*datain.NrBins + i] = rmsynth_array[2*(j*datain.NrBins+i)];
	}
      }
      
//START REGION DEVELOP
      // Write the RM map if requested
      if(write_map) {
	datafile_definition fout;
	cleanPSRData(&fout, application.verbose_state);
	copy_params_PSRData(datain, &fout, application.verbose_state);
	fout.NrSubints = nrrmsteps;
	fout.NrFreqChan = 1;
	fout.NrPols = 1;
	fout.gentype = GENTYPE_RMMAP;
	fout.tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
	free(fout.tsub_list);
	fout.tsub_list = (double *)malloc(sizeof(double));
	if(fout.tsub_list == NULL) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR rmsynth: Memory allocation error");
	  return 0;
	}
	fout.tsub_list[0] = get_tobs(datain, application.verbose_state);
	fout.yrangeset = 1;
	fout.yrange[0] = rmmin;
	fout.yrange[1] = rmmax;
	if(change_filename_extension(argv[argc-1], outputname, "RMmap", 1000, application.verbose_state) == 0) {
	  return 0;
	}
	if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	  return 0;
	if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR rmsynth: Unable to write header.\n");
	  return 0;
	}
	//	appendHistoryLine(fout, argc, argv, application.verbose_state);
	if(writePSRData(&fout, cmap, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR rmsynth: Unable to write data.\n");
	  return 0;
	}
	closePSRData(&fout, 0, application.verbose_state);
      }
//START REGION RELEASE

      // Draw the rm-pulse longitude map
      pgplot_clear_options(&pgplot_options);
      if(pgplot_mapdevice)
	strcpy(pgplot_options.viewport.plotDevice, argv[pgplot_mapdevice]);
      else
	strcpy(pgplot_options.viewport.plotDevice, "?");
      pgplot_options.viewport.ysize = 0.7;
      pgplot_options.viewport.dontclose = 1;
      strcpy(pgplot_options.box.xlabel, "Pulse longitude (bin)");
      strcpy(pgplot_options.box.ylabel, "RM (rad/m\\u2\\d)");
      pgplotMap(&pgplot_options, cmap, datain.NrBins, nrrmsteps, 0, datain.NrBins-1, 0-0.5, datain.NrBins-1 + 0.5, rmmin, rmmax, rmmin, rmmax, PPGPLOT_HEAT, 0, 0, 0, NULL, 0, 0, 1.0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, application.verbose_state);
      
      // Draw the profile above the map
      pgplot_clear_options(&pgplot_options);
      if(pgplot_mapdevice)
	strcpy(pgplot_options.viewport.plotDevice, argv[pgplot_mapdevice]);
      else
	strcpy(pgplot_options.viewport.plotDevice, "?");
      pgplot_options.viewport.dontopen = 1;
      pgplot_options.viewport.noclear = 1;
      pgplot_options.viewport.ysize = 0.25;
      pgplot_options.viewport.dxplot = 0.05;
      pgplot_options.viewport.xsize = 0.75;
      pgplot_options.viewport.dyplot = 0.62;
      strcpy(pgplot_options.box.ylabel, "Pulse profile");
      strcpy(pgplot_options.box.title, "RM synthesis map");
      strcpy(pgplot_options.box.box_xopt, "bcsti");
      
      // Ask the user to select an onpulse region, if it is not specified on command line
      if(application.onpulse.nrRegions == 0) {
	selectRegions(profiledata.data, profiledata.NrBins, &pgplot_options, 0, 0, 0, &(application.onpulse), application.verbose_state);
      }else {
	pgplotGraph1(&pgplot_options, profiledata.data, NULL, NULL, profiledata.NrBins, 0, profiledata.NrBins-1, 0, 0, profiledata.NrBins-1, 0, 0, 0, 1, 0, 1, 0, 1, 1, &(application.onpulse), -1, application.verbose_state);
      }
      ppgend();
      region_int_to_frac(&(application.onpulse), 1.0/(float)datain.NrBins, 0);
      regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);
    }  // End of plotting RM map
    
    
    
    // Draw the RM synthesis power as function of RM, set device etc only once 
    if(bootstrap_itt == 0) {
      pgplot_clear_options(&pgplot_options);
      if(pgplot_rmdevice)
	strcpy(pgplot_options.viewport.plotDevice, argv[pgplot_rmdevice]);
      else
	strcpy(pgplot_options.viewport.plotDevice, "?");
    }
    
    //Loop over pulse longitude bin number (in fact, if collapse is set loop will be terminated after one itteration
    for(binnr = 0; binnr < datain.NrBins; binnr++) {
      
      // Fill the double arrays
      if(collapse) {
	if(application.verbose_state.debug)
	  printf("  itteration %d: Collapsing RM synthesis spectrum\n", bootstrap_itt+1);
    
	// Add pulse longitude bins together in a single spectrum
	collapseRMSynthesisArray(rmsynth_array, nrrmsteps, datain.NrBins, application.onpulse, singlespectrum, application.verbose_state);
	for(i = 0; i < nrrmsteps; i++) {
	  singlespectrum_double[i] = singlespectrum[i];
	}
      }else {
	if(application.verbose_state.debug)
	  printf("  itteration %d: Extracting RM synthesis spectrum for bin %ld\n", bootstrap_itt+1, binnr);
	for(i = 0; i < nrrmsteps; i++) {
	  singlespectrum[i] = rmsynth_array[2*(i*datain.NrBins+binnr)];
	  singlespectrum_double[i] = rmsynth_array[2*(i*datain.NrBins+binnr)];
	}
      }
      
      if(checkRegions(binnr, &application.onpulse, 0, application.verbose_state) != 0 || application.onpulse.nrRegions == 0 || collapse) {
	if(application.verbose_state.debug)
	  printf("  itteration %d: Plotting spectrum for bin %ld\n", bootstrap_itt+1, binnr);
	for(i = 0; i < nrrmsteps; i++) {
	  rmgrid[i] = rmmin + (rmmax-rmmin)*i/(float)(nrrmsteps-1);
	  sigmagrid[i] = 1.0;
	}
	
	pgplot_options.viewport.dontclose = 1;
	strcpy(pgplot_options.box.ylabel, "RM synthesis power");
	strcpy(pgplot_options.box.xlabel, "RM (rad/m\\u2\\d)");
	pgplotGraph1(&pgplot_options, singlespectrum, NULL, NULL, nrrmsteps, rmmin, rmmax, 0, rmmin, rmmax, 0, 0, 0, 0, 0, 1, -10, 1, 1, NULL, -1, application.verbose_state);
	pgplot_options.viewport.dontopen = 1;   // If looping over binnr, we don't have to open device again
	

	/* See if want to fit with parabola, or sync-like function */
	if(parabola) {
	  if(application.verbose_state.debug)
	    printf("  itteration %d: Start fitting of parabola for bin %ld\n", bootstrap_itt+1, binnr);
	  /*
	    ampl*(x-x1)*(x-x2)+base;
	    ampl*(x*x + x*(-x1-x2)+x1*x2)+base;
	    (ampl)*x*x + x*ampl*(-x1-x2)+(ampl*x1*x2+base);
	  */
	  // Define a parabola to be fitted
	  function.nrfuncs = 3;
	  function.func[0].type = FUNC_POLYNOMAL;
	  function.func[0].param[0] = 0;
	  function.func[0].start[0] = singlespectrum[0];
	  function.func[0].value[0] = singlespectrum[0];
	  function.func[0].fit_flag[0] = 1;
	  function.func[1].type = FUNC_POLYNOMAL;
	  function.func[1].param[0] = 1;
	  function.func[1].start[0] = 0.0;
	  function.func[1].value[0] = 0.0;
	  function.func[1].fit_flag[0] = 1;
	  function.func[2].type = FUNC_POLYNOMAL;
	  function.func[2].param[0] = 2;
	  function.func[2].start[0] = 0.0;
	  function.func[2].value[0] = 0.0;
	  function.func[2].fit_flag[0] = 1;

	  // Do the fit of a parabola
	  if(application.verbose_state.debug) {
	    //	    print_fitfunctions(function, 1, 0, -1);
	    //	    printf("Initial guess of fit function: %e + %e*x + %e*x^2\n", function.func[0].start[0], function.func[0].start[1], function.func[0].start[2]);
	  }
	  if(ftol < 0) {
	    ftol = 1e-5;
	  }
	  fit_levmar(1, &function, rmgrid, singlespectrum_double, sigmagrid, nrrmsteps, 0, 1, 0, ftol, 1000, &status, 0, 0, application.verbose_state);
	  if(status != 0) {
	    printerror(application.verbose_state.debug, "Fitting failed: error %d", status);
	    switch(status) {
	    case 1: printerror(application.verbose_state.debug, "Memory allocation error"); return 0;
	    case 2: printerror(application.verbose_state.debug, "Max number of itterations reached"); return 0;
	    case 3: printwarning(application.verbose_state.debug, "WARNING: Machine precision limit reached"); break;
	    case 4: printerror(application.verbose_state.debug, "Cannot determine suitable trial step."); return 0;
	    case 10: printerror(application.verbose_state.debug, "Shouldn't happen, not documented error code of gsl"); return 0;
	    }
	  }
	  if(application.verbose_state.debug) {
	    printf("  itteration %d: Start fitting of parabola for bin %ld done\n", bootstrap_itt+1, binnr);
	    print_fitfunctions(&function, 0, 0, -1, application.verbose_state);
	  }
	  
	  ppgbbuf();
	  ppgsci(2);
	  for(i = 0; i < nrrmsteps; i++) {
	    if(i == 0)
	      ppgmove(rmgrid[i], evaluate_fitfunc_collection(&function, rmgrid[i], application.verbose_state));
	    else
	      ppgdraw(rmgrid[i], evaluate_fitfunc_collection(&function, rmgrid[i], application.verbose_state));
	  }
	  ppgebuf();
	  
	  /*
	    a*x*x+b*x+c
	    2ax + b = 0
	    x = -b/2a
	  */
	  rmcurestimate[binnr] = -0.5*function.func[1].value[0]/function.func[2].value[0];
	}else {
	  float best_rm, best_offset, best_scale, *rmsynth_responds;
	  if(datain.freqMode != FREQMODE_UNIFORM) {
	    fflush(stdout);
	    printwarning(application.verbose_state.debug, "WARNING rmsynth: Frequency range is non-uniform. The instrumental responds is determined using the non-weighted frequencies.");
	  }
	  
	  if(application.verbose_state.debug) 
	    printf("  itteration %d: Fitting expected instrumental shape for bin %ld\n", bootstrap_itt+1, binnr);
	  double chanbw;
	  if(get_channelbandwidth(datain, &chanbw, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR rmsynth (%s): Cannot obtain channel bandwidth.", datain.filename);
	    return 0;
	  }
	  if(ftol < 0) {
	    ftol = 0.001;
	  }
	  if(rmSynthesis_fitInstrumentalResponds(singlespectrum, rmmin, rmmax, nrrmsteps, &best_rm, &best_offset, &best_scale, datain.NrFreqChan, chanbw, get_nonweighted_channel_freq(datain, 0, application.verbose_state), ftol, application.verbose_state) != 0) {
	    printerror(application.verbose_state.debug, "ERROR rmsynth: Fitting of instrumental responds shape failed");
	    return 0;
	  }
	  if(application.verbose_state.debug) 
	    printf("  itteration %d: Generating instrumental shape for bin %ld\n", bootstrap_itt+1, binnr);
	  rmsynth_responds = NULL;
	  if(rmSynthesis_instrument_responds(datain.NrFreqChan, chanbw, get_nonweighted_channel_freq(datain, 0, application.verbose_state), rmmin, rmmax, &rmsynth_responds, NULL, 0, nrrmsteps, best_rm, 0, application.verbose_state) == 0) {
	    return 0;
	  }
	  rmcurestimate[binnr] = best_rm;

	  if(application.verbose_state.debug) 
	    printf("  itteration %d: Plotting instrumental shape for bin %ld\n", bootstrap_itt+1, binnr);
	  ppgbbuf();
	  ppgsci(2);
	  for(i = 0; i < nrrmsteps; i++) {
	    rmsynth_responds[i] *=  best_scale;
	    rmsynth_responds[i] += best_offset;
	    if(i == 0)
	      ppgmove(rmgrid[i], rmsynth_responds[i]);
	    else
	      ppgdraw(rmgrid[i], rmsynth_responds[i]);
	  }
	  ppgebuf();
	  free(rmsynth_responds);
	} // End of if parabola or instrument responds
	
	if(!isnormal(rmcurestimate[binnr])) {
	  printwarning(application.verbose_state.debug, "Not a normal number for RM: %f", rmcurestimate[binnr]);
	  not_normal_result_warning = 1;
	  //	  pgetch();
	}
	/*
	printf("\n\nBest estimate of RM ");
	if(collapse == 0)
	  printf("for bin %d ", binnr);
	printf("is %f\n", rmestimate);
	*/
	
      }else {
	if(application.verbose_state.debug)
	  printf("  itteration %d: Ignoring bin %ld\n", bootstrap_itt+1, binnr);
      } // End of if in onpulse region statement
      if(bootstrap_itt > 0) {
	rm_av[binnr] += rmcurestimate[binnr];
	rm_square[binnr] += rmcurestimate[binnr]*rmcurestimate[binnr];
      }else {
	rmestimate[binnr] = rmcurestimate[binnr];
      }
      if(collapse)  // If fitting collapsed faraday depth, there is only one fit that needs to be done
	break;
     
      if(application.verbose_state.debug) 
	printf("\n");
      if(application.verbose_state.nocounters == 0 && bootstrap && application.verbose_state.verbose && bootstrap_itt > 0) {
	printf("\r%.2f%%       ", 100.0*((bootstrap_itt-1)*datain.NrBins+binnr)/(float)(datain.NrBins*bootstrap));
	fflush(stdout);
      }
      if(application.verbose_state.debug) 
	printf("\n");
    } // End of bin number loop
  
  }  // End of bootstrap loop
  if(application.verbose_state.debug) 
    printf("\n");
  if(bootstrap && application.verbose_state.verbose) {
    printf("\rDone                    \n");
  }

  //  if(bootstrap == 0 || showPAprofile) {
  if(application.verbose_state.verbose)
    printf("Creating de-Faraday rotated profile\n");

  if(application.verbose_state.verbose)
    printf("de-Faraday rotate data\n");
  if(collapse) {
    datain.rm = rmestimate[0];
  }else {
    nrOnpulseBins = 0;
    datain.rm = 0;
    for(binnr = 0; binnr < datain.NrBins; binnr++) {
      if(checkRegions(binnr, &application.onpulse, 0, application.verbose_state) != 0 || application.onpulse.nrRegions == 0) {
	datain.rm += rmestimate[binnr];
	//	printf("Adding RM=%f for bin %ld (total = %f)\n", rmestimate[binnr], binnr, datain.rm);
	nrOnpulseBins++;
      }
    }
    datain.rm /= (double)nrOnpulseBins;
    printwarning(application.verbose_state.debug, "WARNING: To get the analytic errorbar, the average of the longitude-resolved RM was used = %f", datain.rm);
  }
  int ret = preprocess_deFaraday(&datain, 0, 0, 0, NULL, application.verbose_state);
  if(ret == 1) {
    if(fabs(datain.rm) < 1e-3) {
      printwarning(application.verbose_state.debug, "WARNING rmsynth: de-Faraday rotation wasn't applied because of the very low RM. Something appears to be strange.");
      ret = 2;
    }else {
      printerror(application.verbose_state.debug, "de-Faraday rotation failed, cannot determine analytic errorbar");
      return 0;
    }
  }
  if(ret != 2) {
    printerror(application.verbose_state.debug, "de-Faraday rotation failed, cannot determine analytic errorbar");
    return 0;
  }
  if(application.verbose_state.verbose)
    printf("Adding frequency channels\n");
  closePSRData(&clone, 0, application.verbose_state);
  //  free(clone.data); // The rest will be freed with datain
  //  cleanPSRData(&clone, application.verbose_state);
  if(!preprocess_addsuccessiveFreqChans(datain, &clone, datain.NrFreqChan, NULL, application.verbose_state))
    return 0;

  if(application.verbose_state.verbose)
    printf("Calculating L from Q and U\n");
  //  if(alloc_paswing_PSRdata(&clone, 1, 1, 1, 1, application.verbose_state) == 0) {
  //    printerror(application.verbose_state.debug, "Memory allocation error");
  //    return 0;
  //  }
//START REGION DEVELOP
  if(offpulsedata_name[0] != 0) {
    datafile_definition clone2;
    if(!preprocess_addsuccessiveFreqChans(data_offpulse, &clone2, datain.NrFreqChan, NULL, application.verbose_state))
      return 0;
    if(make_paswing_fromIQUV(&clone, 0, application.onpulse2, 0, 1, 1, 1, 0, 0.0, 0.0, &clone2, offpulsedata_rebinfactor, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "Creation of L profile failed, cannot determine analytic errorbar");
      return 0;
    }
    closePSRData(&clone2, 0, application.verbose_state);
  }else {
//START REGION RELEASE
    if(make_paswing_fromIQUV(&clone, 0, application.onpulse2, 0, 1, 1, 1, 0, 0.0, 0.0, NULL, 1.0, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "Creation of L profile failed, cannot determine analytic errorbar");
      return 0;
    }
//START REGION DEVELOP
  }
//START REGION RELEASE

//START REGION DEVELOP
  if(showPAprofile) {
    pgplot_clear_options(&pgplot_options);
    strcpy(pgplot_options.box.title, "RM corrected profile");
    pgplotPAplot(clone, 0, 0, 0, &pgplot_options, "Pulse longitude [deg]", "I,Linear,V", "PA [deg]", "", 0, 360, 0, 0, 0, -90, 90, 0, 3, 1, 1, 0, 0, "-text", "-herrorbar", "-herrorbar2", "-verrorbar", "-verrorbar2", argc, argv, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NULL, NULL, NULL, -90, 90, 1.0, NULL, 1.0, 0, application.verbose_state);
  }
  //START REGION RELEASE

//START REGION DEVELOP
  if(offpulsedata_name[0] == 0) {
//START REGION RELEASE
    double sumL, fwhm, s2nL;
    sumL = 0;
    nrOnpulseBins = 0;
    for(binnr = 0; binnr < datain.NrBins; binnr++) {
      if(checkRegions(binnr, &application.onpulse2, 0, application.verbose_state) != 0 || application.onpulse2.nrRegions == 0) {
	sumL += clone.data[datain.NrBins+binnr];
	nrOnpulseBins++;
      }
    }
    printwarning(application.verbose_state.debug, "WARNING: The analytic errorbar is based on the total L (within -onpulse2) as determined using the overall RM (using -onpulse)."); 
    
    if(collapse) {
      s2nL = sumL/(sqrt(nrOnpulseBins)*clone.offpulse_rms[1]);
    }else {
      s2nL = (sumL/(double)nrOnpulseBins)/(sqrt(nrOnpulseBins)*clone.offpulse_rms[1]);
      printwarning(application.verbose_state.debug, "WARNING: The analytic errorbar as computed is per definition the same for each pulse longitude bin.");
    }
    if(application.verbose_state.verbose) {
      printf("Total L = %f in %d bins\n", sumL, nrOnpulseBins);
      printf("RMS L   = %f\n", clone.offpulse_rms[1]);
      printf("S/N L   = %f\n", s2nL);
    }
    
    fwhm = get_FWHM_RMsynthesis(datain, application.verbose_state);
    expectedRMerror = 0.5*fwhm/s2nL;
    if(application.verbose_state.verbose) {
      printf("Expected FWHM RM synthesis peak = %f\n", fwhm);
      
      printf("Expected errorbar on RM = %f\n", expectedRMerror);
    }

//START REGION DEVELOP
  }else {
    expectedRMerror = sqrt(-1);
    printwarning(application.verbose_state.debug, "WARNING rmsynth: no analytic error is computed when the -offpulsedata option is used.");
  }
//START REGION RELEASE

  closePSRData(&clone, 0, application.verbose_state);


  if(collapse || write_ascii == 0) {
    for(binnr = 0; binnr < datain.NrBins; binnr++) {
      if(checkRegions(binnr, &application.onpulse, 0, application.verbose_state) != 0 || application.onpulse.nrRegions == 0 || collapse) {
	printf("RM ");
	if(collapse == 0)
	  printf("for bin %ld ", binnr);
	printf("is %lf", rmestimate[binnr]);
	if(bootstrap) {
	  rmsigma = sqrt(rm_square[binnr]/(double)bootstrap - rm_av[binnr]*rm_av[binnr]/(double)(bootstrap*bootstrap));
	  printf(" +- %f (from bootstrapping)", rmsigma);
	  if(!isnormal(rmsigma)) {
	    printwarning(application.verbose_state.debug, "Not a normal number for rmsigma: %f (%f %f)", rmsigma, rm_square[binnr], rm_av[binnr]);
	    not_normal_result_warning = 2;
	  }
	}
	printf(" +- %f (analytic prediction for uniform frequency channel weights)", expectedRMerror);

	
	printf("\n");
	if(collapse)
	  break;
      }
    }
  }

  if(write_ascii) {
//START REGION DEVELOP
    if(write_ascii == 2 && collapse == 0) {
      printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot use -asciispec with the -long option enabled.\n");
      return 0;
    }
    if(write_ascii == 2 && bootstrap != 0) {
      printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot use -asciispec with the -bootstrap option enabled.\n");
      return 0;
    }
//START REGION RELEASE
    if(change_filename_extension(argv[argc-1], outputname, "RMtable", 1000, application.verbose_state) == 0) 
      return 0;
    ofile = fopen(outputname, "w+");
    if(ofile == NULL) {
      printerror(application.verbose_state.debug, "ERROR rmsynth: Cannot open file %s for writing", outputname);
      return 0;
    }
    fprintf(ofile, "#");
    if(collapse == 0)
      fprintf(ofile, "binnr Profile ");
    fprintf(ofile, "RM RMerror(analytic)");
    if(bootstrap)
      fprintf(ofile, " RMerror(bootstrap)");
    fprintf(ofile, "\n");
    for(binnr = 0; binnr < datain.NrBins; binnr++) {
      if(checkRegions(binnr, &application.onpulse, 0, application.verbose_state) != 0 || application.onpulse.nrRegions == 0 || collapse) {
	if(collapse == 0)
	  fprintf(ofile, "%ld %f ", binnr, profiledata.data[binnr]);
	fprintf(ofile, "%lf ", rmestimate[binnr]);
	fprintf(ofile, "%lf ", expectedRMerror);
	if(bootstrap) {
	  rmsigma = sqrt(rm_square[binnr]/(float)bootstrap - rm_av[binnr]*rm_av[binnr]/(float)(bootstrap*bootstrap));
	  fprintf(ofile, " %f", rmsigma);
	}
	fprintf(ofile, "\n");
	if(collapse)
	  break;
      }
    }
//START REGION DEVELOP
    if(write_ascii == 2) {
      fprintf(ofile, "RM spectrum:\n");
      for(i = 0; i < nrrmsteps; i++) {
	fprintf(ofile, "%e %e\n", rmmin+(rmmax-rmmin)*i/(double)(nrrmsteps-1), singlespectrum[i]);
      }
    }
//START REGION RELEASE

    fclose(ofile);
  }
  ppgend();

  free(singlespectrum);
  free(singlespectrum_double);
  free(rmgrid);
  free(sigmagrid);
  free(rmsynth_array);
  free(cmap);
  free(rmestimate);
  free(rmcurestimate);
  closePSRData(&datain, 0, application.verbose_state);
  closePSRData(&profiledata, 0, application.verbose_state);
//START REGION DEVELOP
  if(offpulsedata_name[0] != 0) {
    closePSRData(&data_offpulse, 0, application.verbose_state);
  }
//START REGION RELEASE
  // Already closed above.
  //  if(bootstrap)
  //    closePSRData(&clone, 0, application.verbose_state);
  gsl_rng_free (rand_num_gen); 
 
  if(bootstrap) {
    free(rms_channels);
    free(rm_av);
    free(rm_square);
  }


  if(not_normal_result_warning == 1) {
    printwarning(application.verbose_state.debug, "\nNote that warnings were generated for one or more suspicious values found for the RM, suggesting something is wrong.");
  }else if(not_normal_result_warning == 2) {
    printwarning(application.verbose_state.debug, "\nNote that warnings were generated for one or more suspicious values found for the error on the RM, suggesting something is wrong.");
  }
   
  terminateApplication(&application);
  return 0;
}



double get_FWHM_RMsynthesis(datafile_definition datain, verbose_definition verbose)
{
  int freqnr;
  double lambda0, lambda, c, f, norm, sigmal2, fwhm, channelbw; //, dlambda2

  c = 299792458.0;


  // Get the weighted centre wavelength
  lambda0 = 0;
  norm = 0;
  if(get_channelbandwidth(datain, &channelbw, verbose) == 0) {
    printerror(verbose.debug, "ERROR get_FWHM_RMsynthesis (%s): Cannot obtain channel bandwidth.", datain.filename);
    exit(0);
  }
  for(freqnr = 0; freqnr < datain.NrFreqChan; freqnr++) {
    /*
      l = c/f
      dl = -c/f^2 df
      dl^2 = c^2/f^4 df^2 = 
    */
    f = 1e6*get_nonweighted_channel_freq(datain, freqnr, verbose);
    lambda = c/f;
    //    dlambda2 = c*c*channelbw*channelbw*1e12/(f*f*f*f);
    lambda0 += lambda*lambda; //*dlambda2;
    //    norm += dlambda2;
    norm += 1;
  }
  lambda0 /= norm;
  lambda0 = sqrt(lambda0);

  //  printf("lambda0 = %f\n", lambda0);

  /*
    In Brentjens & de Bruyn 2005 the error on RM is defined to be:

    sigmaPA/(sigmal2*sqrt(nrchannels-2))   (Eq. 52)


    This is the error propagation formula for fitting a straight line
    to PA(lambda^2), which should give the RM as the
    gradient. Therefore sigmaPA is the sigmaPA per channel.

    sigmaPA is defined to be

    sigmaPA = sigma/(2*L)     (Eq. 55)

    This is what you get from standard error propagation from the PA
    formula.

    Note that Eq. 52 can be written as

    0.5*[1/sigmal2]/(sqrt(nrchannels-2)*L/sigma) 
    
    The term between square brackets is apparently the FWHM of the RM
    synthesis peak. The L/sigma (per channel) is multiplied with the
    square root of the ~nr of channels, which is effectively the
    overall S/N of the linear polarization.

    sigmaRM = 0.5*FWHM/[S/N of L]

    where we used that in error propagation of the formula of L it
    follows that simga_L = sigma_I.
  */

  sigmal2 = 0;
  for(freqnr = 0; freqnr < datain.NrFreqChan; freqnr++) {
    f = 1e6*get_nonweighted_channel_freq(datain, freqnr, verbose);
    lambda = c/f;
    sigmal2 += pow(lambda, 4) - pow(lambda0, 4);
    //    printf("YYYYYY %f\n", pow(lambda, 4) - pow(lambda0, 4));
    //    sigmal2 += pow(lambda, 4);
  }
  //  printf("XXXX %f\n", sigmal2);
  sigmal2 /= (double)(datain.NrFreqChan-1);
  //  printf("XXXX %f\n", sigmal2);
  //  sigmal2 -= pow(lambda0, 4);
  //  printf("XXXX %f\n", sigmal2);
  //  sigmal2 = fabs(sigmal2);
  sigmal2 = sqrt(sigmal2);

  
  //  printf("XXXX %f\n", sigmal2);

  //  sigmal2 = sqrt(sigmal2/(double)(datain.NrFreqChan-1));

  fwhm = 1.0/(sigmal2);

  return fwhm;
}
