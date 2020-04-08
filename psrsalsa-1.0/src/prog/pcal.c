#include <stdio.h>
#include <time.h>
#include <string.h>
#include "psrsalsa.h"

void SHOWREVISIONINFO_prog() {
#include "pcal.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);

psrsalsaApplication application;

int main(int argc, char **argv)
{
  long i, j, filenumber, phasemodel, diffgainmodel;
  int index, apply_file, mtm_file, mtm_guess_file, show_file, neat, show_fixscale, mtm_leakage, mtmheaderfix, fixed[7], bruteforce, extremes, dowait, combine_file, first, last, analytic_gain, mtm_inverse, apply_ignoreGain, apply_normalise, ignore_dodgy_channels, noparang, apply_expected, refmodel_file, refmodel_file_open, dumpparam, dumpfile, dumpext, dumpresid, showmodel, apply_reconstruct, doMueller, mueller_fdtype, doMueller_row, show_file_overlay;
  char outputname[MaxFilenameLength], outputname_dump[MaxFilenameLength], mtm_dumpfilename[MaxFilenameLength], *filename_ptr, *show_file_device, *profileresid_device, *title;
  int mtm_sequencenr, mtm_totsequencenr, already_aligned, fitphase, fitdiffgain, fitgain, fitell, fitor, writefit, ret;
  datafile_definition receiversolution, *receiversolution_collection, receiversolution_combined, datafile, datafileout, template;
  int receiversolution_used;
  clock_t tstart, tend;
  double rmsSigmaLimit;
  FILE *refmodel_file_fptr;
  //  regions_defenition onpulseregion_fit;  // onpulseregion, 

  // Record start time of program
  tstart = clock();
  initApplication(&application, "pcal", "[options] observation(s)");
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_nocounters = 1;
  application.switch_showrevision = 1;
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.switch_ext = 1;
  application.switch_output = 1;
  application.switch_onpulse = 1; 
  application.switch_onpulsef = 1;
  application.switch_onpulse2 = 1; 
  application.switch_onpulsef2 = 1;
  application.switch_filelist = 1;
  application.switch_rot = 1;     // This is not affecting the template, so not very useful
  application.switch_rotdeg = 1;
  application.switch_history_cmd_only = 1;
  application.oformat = FITS_format;

  apply_file = 0;
  mtm_file = 0;
  mtm_guess_file = 0;
  show_file = 0;
  mtm_leakage = 0;
  mtmheaderfix = 0;
  bruteforce = 0;
  extremes = 0;
  for(i = 0; i < 7; i++)
    fixed[i] = 0;
  sprintf(mtm_dumpfilename, "dump.dat");
  mtm_sequencenr = 0;
  mtm_totsequencenr = 0;
  already_aligned = 0;
  combine_file = 0;
  dowait = 0;
  show_fixscale = 0;
  analytic_gain = 0;
  mtm_inverse = 0;
  apply_ignoreGain = 0;
  show_file_device = NULL;
  profileresid_device = NULL;
  ignore_dodgy_channels = 1;
  noparang = 0;
  apply_expected = 0;
  fitphase = 0;
  fitdiffgain = 0;
  fitgain = 0;
  fitell = 0;
  fitor = 0;
  writefit = 0;
  rmsSigmaLimit = -1;
  title = NULL;
  phasemodel = 0;
  diffgainmodel = 0;
  refmodel_file = 0;
  refmodel_file_open = 0;
  dumpparam = -1;
  dumpfile = 0;
  dumpext = 0;
  dumpresid = 0;
  showmodel = 0;
  doMueller = 0;
  mueller_fdtype = 1; // linear is default for output of Mueller matrix
  doMueller_row = 0; // Default is matrix-like output
  apply_reconstruct = 1;  // Default is to reconstruct data
  apply_normalise = 1; // Default is to normalise the Mueller matrix when applying it to data
  show_file_overlay = 0;
  receiversolution_collection = NULL;
  receiversolution_used = 0;
  //  clearRegion(&onpulseregion);
  //  clearRegion(&onpulseregion_fit);

  if(argc < 2) {
    printf("Program to do polarization calibration.\n\n");
    printApplicationHelp(&application);
    printf("MTM fit options:\n");
    printf("  -mtm             Use the specified template to fit the observation, resulting\n");
    printf("                   in receiver parameters (excluding leakage terms).\n");
    printf("  -mtmfull         Same as -mtm, but now obtain the full Mueller matrix.\n");
    printf("  -mtmguess        Use the speciefied receiver model as the initial guess. The\n");
    printf("                   gain will be ignored from the initial guess.\n");
    printf("  -mtmheaderfix    Fix some common encountered problems with the header of the\n");
    printf("                   template: i.e. assume the template are Stokes parameters.\n");
    printf("  -mtmfix          Fix some of the receiver parameters, so \"0 1 1 0 0 0 0\"\n");
    printf("                   fixes the differential gain and phase.\n");
    printf("  -mtmanalyticgain Experimental option. Instead of treating the gain as a free\n");
    printf("                   parameter, it is derived from the comparison between Stokes\n");
    printf("                   I of the template and the measurement. This might work better\n");
    printf("                   if the S/N is low. Use only together with -mtminverse.\n");
    printf("  -mtminverse      Experimental option. Do fit by correcting data rather than\n");
    printf("                   distorting template.\n");
    printf("  -mtmbrute        Try (1+2*bruteforce)**(nr fit parameters) initial guesses,\n");
    printf("                   which is in general a very large number!\n");
    printf("  -mtmextremes     Try 3**5=243 initial guesses with -mtmfull mode (3 otherwise)\n");
    printf("                   which covers the extreme posibilities for the angles.\n");
    printf("  -mtmprep         Prepare to run mtm fitting over multiple instances. Give\n");
    printf("                   \"dumpfilename nrinstances\"\n");
    printf("  -mtmrun          Run mtm fitting for one of the instances prepared with\n");
    printf("                   -mtmprep. Give \"dumpfilename instancenumber\", starting\n");
    printf("                   counting from 1\n");
    printf("  -noalign         Assume template and data are already aligned, so do not a\n");
    printf("                   cross-correlation.\n");
    //    printf("  -debugdev        Specify device for the before/after plots (use with -debug).\n");
    printf("  -mtmprofdev      Specify device for the before/after profile plots.\n");
    printf("  -wait            Do not immediately quit. Can be useful with -mtmrun.\n");
    printf("  -noparang        Disable paralactic angle corrections to be made.\n");

    printf("\nCalibrate data options:\n");
    printf("  -apply           Specify receiver solution file (normalised by default) and\n");
    printf("                   apply it to the observation to correct it.\n");
    printf("  -apply_nonorm    By default, the solution is normalised such that the rms of\n");
    printf("                   Stokes I is unchanged (if the rms is identical for the four\n");
    printf("                   Stokes parameters. This option will disable this.\n");
    printf("  -apply_nogain    Use with -apply. When specified, the gain is ignored, which\n");
    printf("                   implies -apply_nonorm\n");
    printf("  -apply_dodgy     Use with -apply. When specified, the dodgy looking leakage\n");
    printf("                   solutions are not ignored.\n");
    printf("  -apply_inverse   Use with -apply. When specified, the inverse of the Mueller\n");
    printf("                   matrix is applied (i.e. to distort observation).\n");

    printf("\nReceiver solution options:\n");
    printf("  -combine fname   The specified output file name will contain the best\n");
    printf("                   solution for each channel as found in the input file(s).\n");
    printf("                   If used in combination with -fitgain etc., fname is the\n");
    printf("                   first (possibly only) input file with the receiver\n");
    printf("                   parameters replaced by the fit.\n");
    printf("  -dump            specify a parameter to dump:\n");
    printf("                      0 - Gain\n");
    printf("                      1 - Differential gain [%%]\n");
    printf("                      2 - Differential phase [deg]\n");
    printf("                      3 - Elipticity 0 [deg]\n");
    printf("                      4 - Orientation 0 [deg]\n");
    printf("                      5 - Elipticity 1 [deg]\n");
    printf("                      6 - Orientation 1 [deg]\n");
    printf("                      7 - Chi-square\n");
    printf("                      8 - Nfree\n");
    printf("  -dumpall         Dump all parameters\n");
    printf("  -dumpresid       Same as -dump, except that if a fit is made (or reference\n");
    printf("                   model is specified), it is subtracted\n");
    printf("  -dumpresidall    Same as -dumpresid, but all parameters\n");
    printf("  -dumpresid2      Same as -dumpresid, but reference model is not subtracted.\n");
    printf("  -dumpresid2all   Same as -dumpresidall, but reference model is not subtracted.\n");
    printf("  -dumpfile        Redirect output to this file rather than to the stdout\n");
    printf("  -dumpext         Similar to -dumpfile, but specify output file extension\n");
    printf("  -mueller         Show the Mueller matrix (i.e. obs = Mueller * intrinsic).\n");
    printf("  -imueller        Show the inverse of the Mueller matrix (i.e. intrinsic = Mueller * obs).\n");
    printf("  -mueller_unity   Show the Mueller matrix multiplied by its inverse (should be unity matrix).\n");
    printf("  -mueller_row     Use in combination with -mueller or -imueller: formats output in single rows.\n");
    printf("  -mueller_cir     Use in combination with -mueller or -imueller: a circular basis rather than a linear basis is assumed.\n");
    printf("  -show            View receiver solutions (can be used in combination with -mtm\n");
    printf("                   options).\n");
    printf("  -showneat        View receiver solutions, but have multiple params per plot.\n");
    printf("  -showneat2       View receiver solutions on a single page.\n");
    printf("  -showdevice      Plot the receiver solution on this device.\n");
    printf("  -showoverlay     Don't produce separate graphs for separate receiver solutions.\n");
    printf("                   This option implies -showneat2.\n");
    printf("  -title           Set the title for the receiver solution plot. See\n");
    printf("                   -titlekeywords for possible keywords to use.\n");
    printf("  -titlekeywords   Lists to keywords you can use in the -title option.\n");
    printf("  -showmodel       Adds text to the panels identifying type of model being fit.\n");
    printf("  -showfix         Fix the y-range of some plots. This option implies -showneat.\n");
    printf("  -fitgain n       Fit gain with nth order polynomial.\n");
    printf("  -fitdiffgain n   Fit differential gain with nth order polynomial.\n");
    printf("  -fitphase        Fit differential phase with nth order polynomial.\n");
    printf("  -fitel           Fit ellipticities with nth order polynomial.\n");
    printf("  -fitor           Fit orientations with nth order polynomial.\n");
    printf("  -fitsigma        Do not include points which are more than the specified\n");
    printf("                   factor away from the solution.\n");
    printf("  -writefit        Write out fit to file (can use together with -ext).\n");
    printf("  -phasemodel      Before fitting, subtract the specified reference model\n");
    printf("                   number for the differential phase. Specify -1 to search\n");
    printf("                   for optimum available reference model.\n");
    printf("  -diffgainmodel   Before fitting, subtract the specified reference model\n");
    printf("                   number for the differential gain. Specify -1 to search\n");
    printf("                   for optimum available reference model.\n");
    printf("  -refmodelfile    Best to be used together with -filelist: Read the reference\n");
    printf("                   model numbers from a file, which is just a list of numbers.\n");
    printf("                   Each line contains a -phasemodel and -diffgainmodel number.\n");
    printf("  -refmodellist    Show which reference models are available\n");
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-mtmanalyticgain") == 0) {
	analytic_gain = 1;
      }else if(strcmp(argv[i], "-apply") == 0) {
	apply_file = i+1;
	i++;
      }else if(strcmp(argv[i], "-noparang") == 0) {
	noparang = 1;
      }else if(strcmp(argv[i], "-apply_nogain") == 0) {
	apply_ignoreGain = 1;
	apply_normalise = 0;
	apply_expected = 1;
      }else if(strcmp(argv[i], "-apply_nonorm") == 0) {
	apply_normalise = 0;
	apply_expected = 1;
      }else if(strcmp(argv[i], "-apply_dodgy") == 0) {
	ignore_dodgy_channels = 0;
	apply_expected = 1;
      }else if(strcmp(argv[i], "-apply_inverse") == 0 || strcmp(argv[i], "-apply_invert") == 0) {
	apply_expected = 1;
	apply_reconstruct = 0;
      }else if(strcmp(argv[i], "-show") == 0) {
	show_file = 1;
      }else if(strcmp(argv[i], "-showoverlay") == 0) {
	show_file = 3;
	show_file_overlay = 1;
      }else if(strcmp(argv[i], "-showdevice") == 0) {
	show_file_device = argv[i+1];
	i++;
      }else if(strcmp(argv[i], "-mueller") == 0) {
	if(doMueller == 2) {
	  printerror(application.verbose_state.debug, "pcal: -mueller and -imueller cannot be specified together");
	  return 0;
	}
	if(doMueller == 3) {
	  printerror(application.verbose_state.debug, "pcal: -mueller and -mueller_unity cannot be specified together");
	  return 0;
	}
	doMueller = 1;
      }else if(strcmp(argv[i], "-imueller") == 0) {
	if(doMueller == 1) {
	  printerror(application.verbose_state.debug, "pcal: -mueller and -imueller cannot be specified together");
	  return 0;
	}
	if(doMueller == 3) {
	  printerror(application.verbose_state.debug, "pcal: -imueller and -mueller_unity cannot be specified together");
	  return 0;
	}
	doMueller = 2;
      }else if(strcmp(argv[i], "-mueller_unity") == 0) {
	if(doMueller == 1) {
	  printerror(application.verbose_state.debug, "pcal: -mueller_unity and -mueller cannot be specified together");
	  return 0;
	}
	if(doMueller == 2) {
	  printerror(application.verbose_state.debug, "pcal: -mueller_unity and -imueller cannot be specified together");
	  return 0;
	}
	doMueller = 3;
      }else if(strcmp(argv[i], "-mueller_row") == 0) {
	doMueller_row = 1;
      }else if(strcmp(argv[i], "-mueller_cir") == 0) {
	mueller_fdtype = 2;
      }else if(strcmp(argv[i], "-title") == 0) {
	title = argv[i+1];
	i++;
      }else if(strcmp(argv[i], "-titlekeywords") == 0) {
	str_list_replace_keys(0);
	return 0;
	//      }else if(strcmp(argv[i], "-debugdev") == 0) {
      }else if(strcasecmp(argv[i], "-mtmprofdev") == 0 || strcasecmp(argv[i], "-mtmprofiledev") == 0) {
	profileresid_device = argv[i+1];
	i++;
      }else if(strcmp(argv[i], "-showneat") == 0) {
	show_file = 2;
      }else if(strcmp(argv[i], "-showneat2") == 0) {
	show_file = 3;
      }else if(strcmp(argv[i], "-showfix") == 0 || strcmp(argv[i], "-showfixed") == 0) {
	if(show_file < 2) 
	  show_file = 2;
	show_fixscale = 1;
      }else if(strcmp(argv[i], "-refmodelfile") == 0) {
	refmodel_file = i+1;
	i++;
      }else if(strcmp(argv[i], "-refmodellist") == 0) {
	refmodellist();
	return 0;
      }else if(strcmp(argv[i], "-showmodel") == 0) {
	showmodel = 1;
      }else if(strcmp(argv[i], "-mtmheaderfix") == 0) {
	mtmheaderfix = 1;
      }else if(strcmp(argv[i], "-combine") == 0) {
	combine_file = i+1;
	i++;
      }else if(strcmp(argv[i], "-mtmguess") == 0) {
	mtm_guess_file = i+1;
	i++;
      }else if(strcmp(argv[i], "-mtm") == 0) {
	mtm_file = i+1;
	i++;
      }else if(strcmp(argv[i], "-mtmfix") == 0) {
	j = sscanf(argv[i+1], "%d %d %d %d %d %d %d", &fixed[0], &fixed[1], &fixed[2], &fixed[3], &fixed[4], &fixed[5], &fixed[6]);
        i++;
      }else if(strcmp(argv[i], "-writefit") == 0) {
	writefit = 1;
      }else if(strcmp(argv[i], "-fitphase") == 0) {
	j = sscanf(argv[i+1], "%d", &fitphase);
	fitphase += 1;
        i++;
      }else if(strcmp(argv[i], "-phasemodel") == 0) {
	j = sscanf(argv[i+1], "%ld", &phasemodel);
        i++;
      }else if(strcmp(argv[i], "-diffgainmodel") == 0) {
	j = sscanf(argv[i+1], "%ld", &diffgainmodel);
        i++;
      }else if(strcmp(argv[i], "-fitdiffgain") == 0) {
	j = sscanf(argv[i+1], "%d", &fitdiffgain);
	fitdiffgain += 1;
        i++;
      }else if(strcmp(argv[i], "-fitgain") == 0) {
	j = sscanf(argv[i+1], "%d", &fitgain);
	fitgain += 1;
        i++;
      }else if(strcmp(argv[i], "-fitel") == 0 || strcmp(argv[i], "-fitell") == 0) {
	j = sscanf(argv[i+1], "%d", &fitell);
	fitell += 1;
        i++;
      }else if(strcmp(argv[i], "-fitor") == 0) {
	j = sscanf(argv[i+1], "%d", &fitor);
	fitor += 1;
        i++;
      }else if(strcmp(argv[i], "-fitsigma") == 0) {
	j = sscanf(argv[i+1], "%lf", &rmsSigmaLimit);
        i++;
      }else if(strcmp(argv[i], "-mtmbrute") == 0) {
	j = sscanf(argv[i+1], "%d", &bruteforce);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-mtmprep") == 0) {
	j = sscanf(argv[i+1], "%s %d", mtm_dumpfilename, &mtm_totsequencenr);
	if(j != 2) {
	  printerror(application.verbose_state.debug, "Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	mtm_sequencenr = -1;
        i++;
      }else if(strcmp(argv[i], "-mtmrun") == 0) {
	j = sscanf(argv[i+1], "%s %d", mtm_dumpfilename, &mtm_sequencenr);
	if(j != 2) {
	  printerror(application.verbose_state.debug, "Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	mtm_totsequencenr = 10;  /* Will be read in from dumpfile */
        i++;
      }else if(strcmp(argv[i], "-mtmbrute") == 0) {
	j = sscanf(argv[i+1], "%d", &bruteforce);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-dump") == 0) {
	j = sscanf(argv[i+1], "%d", &dumpparam);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-dumpall") == 0) {
	dumpparam = 99;
      }else if(strcmp(argv[i], "-dumpresidall") == 0) {
	dumpparam = 99;
	dumpresid = 1;
      }else if(strcmp(argv[i], "-dumpresid2all") == 0) {
	dumpparam = 99;
	dumpresid = 2;
      }else if(strcmp(argv[i], "-dumpresid") == 0) {
	j = sscanf(argv[i+1], "%d", &dumpparam);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	dumpresid = 1;
        i++;
      }else if(strcmp(argv[i], "-dumpresid2") == 0) {
	j = sscanf(argv[i+1], "%d", &dumpparam);
	if(j != 1) {
	  printerror(application.verbose_state.debug, "Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	dumpresid = 2;
        i++;
      }else if(strcmp(argv[i], "-dumpfile") == 0) {
	dumpfile = i+1;
        i++;
      }else if(strcmp(argv[i], "-dumpext") == 0) {
	dumpext = i+1;
        i++;
      }else if(strcmp(argv[i], "-mtmextremes") == 0) {
	extremes = 1;
      }else if(strcmp(argv[i], "-mtminverse") == 0) {
	mtm_inverse = 1;
      }else if(strcmp(argv[i], "-mtmfull") == 0) {
	mtm_file = i+1;
	i++;
	mtm_leakage = 1;
      }else if(strcmp(argv[i], "-noalign") == 0) {
	already_aligned = 1;
      }else if(strcmp(argv[i], "-wait") == 0) {
	dowait = 1;
      }else {
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "pcal: Unknown option: %s\n\nRun pcal without command line arguments to show help", argv[i]);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }

  if(apply_expected && apply_file == 0) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "pcal: Please also use the -apply option (in addition to -apply_gain and/or -apply_dodgy and/or -apply_inverse) if you want to apply a receiver model.");
    return 0;
  }

  if(dumpext && dumpfile) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "You cannot use -dumpfile and -dumpext together");
    return 0;
  }

  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "pcal: No files specified");
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) > 1) {
    fflush(stdout);
    printwarning(application.verbose_state.debug, "WARNING: Multiple files are being processed. If binning happens to be different in different observations I assume onpulse selection will not work correctly.");
  }

  if(show_file && (mtm_sequencenr || mtm_totsequencenr)) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR pcal: You cannot show receiver solution when running multiple instances.");
    return 0;
  }
  if(doMueller && (mtm_sequencenr || mtm_totsequencenr)) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR pcal: You cannot output Mueller matrices when running multiple instances.");
    return 0;
  }
  if((fitphase || fitdiffgain || fitgain || fitell || fitor) && (mtm_sequencenr || mtm_totsequencenr)) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR pcal: You cannot fit model to receiver solution when running multiple instances.");
    return 0;
  }

  //  cleanPSRData(&datafile, application.verbose_state);


  int profileresid_device_id;   // The device to be used to make before/after calibration profile plots

  first = 1;
  last = 0;
  filenumber = 0;
  while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
    /* Check if this is last file, which we need to know when plotting multiple receiver models. */
    if(filenumber == numberInApplicationFilenameList(&application, argv, application.verbose_state) - 1)
      last = 1;
    /* Load data file from provided file list (not a receiver solution) when we are applying reveicer solution or mtm file. */

    if(refmodel_file) {
      if(refmodel_file_open == 0) {
	refmodel_file_fptr = fopen(argv[refmodel_file], "r");
	if(refmodel_file_fptr == NULL) {
	  printerror(application.verbose_state.debug, "ERROR pcal: Cannot open file %s", argv[refmodel_file]);
	  return 0;
	}
	refmodel_file_open = 1;
      }
      ret = fscanf(refmodel_file_fptr, "%ld %ld", &phasemodel, &diffgainmodel);
      if(ret != 2) {
	printerror(application.verbose_state.debug, "ERROR pcal: Cannot read reference model number from file %s", argv[refmodel_file]);
	return 0;
      }
      printf("File %s: reference models %ld and %ld is used\n", argv[refmodel_file], phasemodel, diffgainmodel);
    }

    if((apply_file || mtm_file) && mtm_sequencenr <= 0) {  // Do not open file if going to read dump file
      if(!openPSRData(&datafile, filename_ptr, application.iformat, 0, 1, 0, application.verbose_state))
	return 0;
      if(datafile.poltype == POLTYPE_COHERENCY) {
	if(preprocess_stokes(&datafile, application.verbose_state) == 0)
	  return 0;
      }
      /* Search commandline for header parameters to overwrite header
	 parameters. */
      if(PSRDataHeader_parse_commandline(&datafile, argc, argv, application.verbose_state) == 0)
	return 0;
      for(i = 1; i < argc; i++) {
	if(strcmp(argv[i], "-header") == 0) {
	  fflush(stdout);
	  printwarning(application.verbose_state.debug, "WARNING: If using the -header option, be aware it applied BEFORE the preprocessing.");
	}
      }
      if(preprocessApplication(&application, &datafile) == 0) {
	return 0;
      }
    }
    /* Generate output filename */
    if(getOutputName(&application, filename_ptr, outputname, application.verbose_state) == 0) {   // Even do this while doing -mtmrun
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR pcal: Changing filename failed");
      return 0;
    }

    if(mtm_file || (show_file && show_file_overlay == 0) || doMueller || (combine_file == 0 && (fitphase || fitdiffgain || fitgain || fitell || fitor)) || dumpparam != -1) {
      if(mtm_file) {
	if(mtm_sequencenr <= 0) {  // Do not open template if going to read dump file
	  //cleanPSRData(&template, application.verbose_state);
	  if(!openPSRData(&template, argv[mtm_file], application.iformat, 0, 1, 0, application.verbose_state))
	    return 0;
	  if(mtmheaderfix) {
	    template.poltype = POLTYPE_STOKES;
	  }
	  if(template.poltype == POLTYPE_COHERENCY) {
	    if(preprocess_stokes(&template, application.verbose_state) == 0)
	      return 0;
	  }
	  if(preprocessApplication(&application, &template) == 0) {
	    return 0;
	  }
	}
	/* Load a receiver model as an initial guess if provided */
	datafile_definition *guess_ptr;
	if(mtm_guess_file) {
	  //  cleanPSRData(&receiversolution, application.verbose_state);
	  receiversolution_used = 1;
	  if(!openPSRData(&receiversolution, argv[mtm_guess_file], application.iformat, 0, 1, 0, application.verbose_state))
	    return 0;
	  guess_ptr = &receiversolution;
	  receiversolution_used = 1;
	  //	  printf("OPEN RECEIVERSOLUTION 1\n");
	}else {
	  guess_ptr = 0;
	}
	if(filenumber == 0) {
	  if(profileresid_device != NULL) {   // Open device ourselves, so the cal_mtm() can append plots for each file.
	    profileresid_device_id = ppgopen(profileresid_device);
	  }else {
	    profileresid_device_id = 0;
	  }
	}
	if(cal_mtm(&datafile, &template, &application.onpulse, &application.onpulse2, already_aligned, guess_ptr, mtm_leakage, fixed, analytic_gain, mtm_inverse, bruteforce, extremes, outputname, argc, argv, mtm_dumpfilename, mtm_sequencenr, mtm_totsequencenr, noparang, application.verbose_state, NULL, profileresid_device_id) == 0) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR pcal: MTM method failed");
	  return 0;      
	}
	if(last && profileresid_device_id) {
	  ppgslct(profileresid_device_id);
	  ppgend();
	}

	/* Set filename_ptr to output file which was written out, so it is loaded when showing the receiver solution */
	filename_ptr = outputname;
	closePSRData(&template, 0, application.verbose_state);
      }
      if((show_file && show_file_overlay == 0) || doMueller || (combine_file == 0 && show_file_overlay == 0 && (fitphase || fitdiffgain || fitgain || fitell || fitor)) || dumpparam != -1) {
	//  cleanPSRData(&receiversolution, application.verbose_state);
	receiversolution_used = 1;
	//	printf("OPEN RECEIVERSOLUTION 2\n");
	if(!openPSRData(&receiversolution, filename_ptr, application.iformat, 0, 1, 0, application.verbose_state))
	  return 0;
	analyticReceiverSolution_def analyticsolution;
	int havemodel;
	havemodel = 0;
	//	phasemodel = 0;  // THIS LINE DISABLES SUBTRACTING OF PHASE MODEL
	if(combine_file == 0 && (fitphase || fitdiffgain || fitgain || fitell || fitor)) {
	  if(fitReceiverModel(receiversolution, fitgain, fitdiffgain, fitphase, phasemodel, diffgainmodel, fitell, fitell, fitor, fitor, rmsSigmaLimit, &analyticsolution, application.verbose_state) == 0) {
	    fflush(stdout);
	    printerror(application.verbose_state.debug, "ERROR pcal: Fitting of an analytic solution to receiver model failed.");
	    return 0;
	  }
	  printf("The following analytic solution was fitted:\n");
	  printAnalyticReceiverSolution(analyticsolution, 2);
	  if(writefit) {
	    if(getOutputName(&application, filename_ptr, outputname, application.verbose_state) == 0) {
	      fflush(stdout);
	      printerror(application.verbose_state.debug, "ERROR pcal: Changing filename failed");
	      return 0;
	    }
	    if(writeAnalyticReceiverSolution(outputname, analyticsolution) == 0) {
	      fflush(stdout);
	      printerror(application.verbose_state.debug, "ERROR pcal: Writing analytic receiver solution failed");
	      return 0;
	    }
	  }
	  havemodel = 1;
	}

	if(show_file && show_file_overlay == 0) {   // If overlaying solution, plot them all at once after loop
	  // Show receiver solution
	  if(show_file == 2)
	    neat = 1;
	  else if(show_file == 3)
	    neat = 2;
	  else
	    neat = 0;
	  if(havemodel) {
	    if(showReceiverModel(&receiversolution, 1, neat, show_fixscale, &analyticsolution, show_file_device, first, last, show_file_overlay, title, showmodel, application.verbose_state) == 0) {
	      fflush(stdout);
	      printerror(application.verbose_state.debug, "ERROR pcal: Showing receiver model failed.");
	      return 0;
	    }
	  }else {
	    if(showReceiverModel(&receiversolution, 1, neat, show_fixscale, NULL, show_file_device, first, last, show_file_overlay, title, showmodel, application.verbose_state) == 0) {
	      fflush(stdout);
	      printerror(application.verbose_state.debug, "ERROR pcal: Showing receiver model failed.");
	      return 0;
	    }
	  }
	}

	if(doMueller) {
	  // Print Mueller matrices of solution
	  if(printMuellerFromReceiverModel(receiversolution, doMueller-1, 0, 0, mueller_fdtype, doMueller_row, 1, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pcal: Printing Mueller matrix elements of receiver model failed.");
	    return 0;
	  }
	}

	if(dumpparam >= 0) {
	  analyticReceiverSolution_def *analyticsolution_ptr, analyticsolution2;
	  float *xdata[9], *ydata[9], *yerrorbar[9];
	  long nrpoints;
	  int paramnr, test_value;
	  FILE *fout;
	  analyticsolution_ptr = NULL;
	  if(dumpresid && havemodel) {
	    memcpy(&analyticsolution2, &analyticsolution, sizeof(analyticReceiverSolution_def));
	    if(dumpresid == 2) {   // Ignore reference model if required
	      analyticsolution2.phaseRefmodel = 0;
	      analyticsolution2.diffgainRefmodel = 0;
	    }
	    analyticsolution_ptr = &analyticsolution2;
	  }
	  for(paramnr = 0; paramnr <= 8; paramnr++) {
	    test_value = get_polnr_in_receiversolution(receiversolution, paramnr);
	    if(dumpparam == paramnr || (dumpparam == 99 && test_value >= 0)) {
	      if(extractReceiverParameter(paramnr, receiversolution, 0, analyticsolution_ptr, &xdata[paramnr], &ydata[paramnr], &yerrorbar[paramnr], &nrpoints, application.verbose_state) == 0) {
		fflush(stdout);
		printerror(application.verbose_state.debug, "ERROR pcal: Dumping of receiver parameter failed.");
		return 0;
	      }
	    }
	  }
	  if(dumpfile > 0 || dumpext > 0) {
	    if(dumpfile > 0) {
	      strcpy(outputname_dump, argv[dumpfile]);
	    }else if(dumpext > 0) {
	      if(change_filename_extension(filename_ptr, outputname_dump, argv[dumpext], MaxFilenameLength, application.verbose_state) == 0) {
		fflush(stdout);
		printerror(application.verbose_state.debug, "ERROR pcal: Changing filename extension failed.");
		return 0;		
	      }
	    }
	    fout = fopen(outputname_dump, "w");
	    if(fout == NULL) {
	      fflush(stdout);
	      printerror(application.verbose_state.debug, "ERROR pcal: Cannot open file %s.", outputname_dump);
	      return 0;
	    }
	  }else {
	    fout = stdout;
	  }
	  fprintf(fout, "#Dump of parameter: ");
	  if(dumpparam == 0 || (dumpparam == 99 && get_polnr_in_receiversolution(receiversolution, 0) >= 0)) {
	    fprintf(fout, "Gain, ");
	  }
	  if(dumpparam == 1 || (dumpparam == 99 && get_polnr_in_receiversolution(receiversolution, 1) >= 0)) {
	    fprintf(fout, "Differential gain (gamma) [%%], ");
	  }
	  if(dumpparam == 2 || (dumpparam == 99 && get_polnr_in_receiversolution(receiversolution, 2) >= 0)) {
	    fprintf(fout, "Differential phase [deg], ");
	  }
	  if(dumpparam == 3 || (dumpparam == 99 && get_polnr_in_receiversolution(receiversolution, 3) >= 0)) {
	    fprintf(fout, "Elipticity 0 [deg], ");
	  }
	  if(dumpparam == 4 || (dumpparam == 99 && get_polnr_in_receiversolution(receiversolution, 4) >= 0)) {
	    fprintf(fout, "Orientation 0 [deg], ");
	  }
	  if(dumpparam == 5 || (dumpparam == 99 && get_polnr_in_receiversolution(receiversolution, 5) >= 0)) {
	    fprintf(fout, "Elipticity 1 [deg], ");
	  }
	  if(dumpparam == 6 || (dumpparam == 99 && get_polnr_in_receiversolution(receiversolution, 6) >= 0)) {
	    fprintf(fout, "Orientation 1 [deg], "); 
	  }
	  if(dumpparam == 7 || (dumpparam == 99 && get_polnr_in_receiversolution(receiversolution, 7) >= 0)) {
	    fprintf(fout, "Chi-square, ");
	  }
	  if(dumpparam == 8 || (dumpparam == 99 && get_polnr_in_receiversolution(receiversolution, 8) >= 0)) {
	    fprintf(fout, "Nfree, ");
	  }
	  fprintf(fout, " (freq [MHz] value error)\n");
	  for(i = 0; i < nrpoints; i++) {
	    if(dumpparam == 99)
	      test_value = 0;
	    else
	      test_value = dumpparam;
	    fprintf(fout, "%f", (xdata[test_value])[i]);
	    for(paramnr = 0; paramnr <= 8; paramnr++) {
	      test_value = get_polnr_in_receiversolution(receiversolution, paramnr);
	      if(dumpparam == paramnr || (dumpparam == 99 && test_value >= 0)) {
		fprintf(fout, " %f %f", (ydata[paramnr])[i], (yerrorbar[paramnr])[i]);
	      }
	    }
	    fprintf(fout, "\n");
	  }

	  if(dumpfile > 0) {
	    fclose(fout);
	  }

	  for(paramnr = 0; paramnr <= 8; paramnr++) {
	    free(xdata[paramnr]);
	    free(ydata[paramnr]);
	    free(yerrorbar[paramnr]);
	  }
	  
	}

      }
    }

    if(combine_file || show_file_overlay) {  // Read in a collection of receiver solutions
      // Sometimes (always) already opened? For now, just close it to avoid memory leak.
      if(receiversolution_used) {
	closePSRData(&receiversolution, 0, application.verbose_state);
	receiversolution_used = 0;
      }
      //  cleanPSRData(&receiversolution, application.verbose_state);
      //      printf("OPEN RECEIVERSOLUTION 3\n");
      if(!openPSRData(&receiversolution, filename_ptr, application.iformat, 0, 1, 0, application.verbose_state))
	return 0;
      receiversolution_used = 1;
      if(receiversolution.gentype != GENTYPE_RECEIVERMODEL && receiversolution.gentype != GENTYPE_RECEIVERMODEL2) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR pcal: Input data in %s does not appears to be a receiver solution.", filename_ptr);
	return 0;	
      }
      // Allocate memory to store all individual solutions
      if(first) {
	long nrfiles = numberInApplicationFilenameList(&application, argv, application.verbose_state);
	receiversolution_collection = malloc(nrfiles*sizeof(datafile_definition));
	if(receiversolution_collection == NULL) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR pcal: Memory allocation error.");
	  return 0;	
	}
	//	printf("Allocated %ld bytes for receiver solutions\n", nrfiles*sizeof(datafile_definition));
      }
      // Make a copy for later use
      cleanPSRData(&(receiversolution_collection[filenumber]), application.verbose_state);
      copy_params_PSRData(receiversolution, &(receiversolution_collection[filenumber]), application.verbose_state);
      (receiversolution_collection[filenumber]).data = malloc((receiversolution_collection[filenumber]).NrSubints*(receiversolution_collection[filenumber]).NrBins*(receiversolution_collection[filenumber]).NrPols*(receiversolution_collection[filenumber]).NrFreqChan*sizeof(float));
      if((receiversolution_collection[filenumber]).data == NULL) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR pcal: Memory allocation error.");
	return 0;
      }
      memcpy((receiversolution_collection[filenumber]).data, receiversolution.data, (receiversolution_collection[filenumber]).NrSubints*(receiversolution_collection[filenumber]).NrBins*(receiversolution_collection[filenumber]).NrPols*(receiversolution_collection[filenumber]).NrFreqChan*sizeof(float));

      closePSRData(&receiversolution, 0, application.verbose_state);
      receiversolution_used = 0;
    }

    if(apply_file) {
      if(!openPSRData(&receiversolution, argv[apply_file], application.iformat, 0, 1, 0, application.verbose_state))
	return 0;
      
      if(datafile.feedtype == FEEDTYPE_UNKNOWN) {
	printerror(application.verbose_state.debug, "ERROR pcal (%s): The data file has an unspecified type of feed (circular? linear?). This information is needed before a receiver solution can be applied.", filename_ptr);
	return 0;
      }
      
      if(receiversolution.feedtype == FEEDTYPE_UNKNOWN) {
	printwarning(application.verbose_state.debug, "WARNING pcal (%s): The receiver parameter file has an unspecified type of feed (circular? linear?). It is assumed this is the same as that of the data file.", filename_ptr);
	receiversolution.feedtype = datafile.feedtype;
      }

      if(receiversolution.feedtype != datafile.feedtype) {
	printerror(application.verbose_state.debug, "ERROR pcal (%s): The data file and receiver file have a different type of feed (circular? linear?) specified. This suggest that this is not what was intended.", filename_ptr);
	return 0;
      }

      if(datafile.isDeFarad != 0 && datafile.isDeFarad != 1) {
	printerror(application.verbose_state.debug, "ERROR pcal (%s): The data file has an unspecified de-Faraday rotation state, which must be known before receiver model can be applied.", filename_ptr);
	return 0;
      }
      if(datafile.isDePar != 0 && datafile.isDePar != 1) {
	printerror(application.verbose_state.debug, "ERROR pcal (%s): The data file has an unspecified parallactic angle removal state, which must be known before receiver model can be applied.", filename_ptr);
	return 0;
      }

      //      if(apply_reconstruct) {  // Default: correct the data
	if(datafile.isDeFarad) {
	  printwarning(application.verbose_state.debug, "WARNING pcal (%s): Faraday rotation has been removed from the data file, which will be undone before receiver model can be applied.", filename_ptr);
	  if(preprocess_deFaraday(&datafile, 1, 0, 0, NULL, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pcal (%s): Removal of de-Faraday rotation failed.", filename_ptr);
	    return 0;
	  }
	}
	if(datafile.isDePar) {
	  printwarning(application.verbose_state.debug, "WARNING pcal (%s): Parallactic angle effect has been removed from the data file, which will be undone before receiver model can be applied.", filename_ptr);
	  if(preprocess_corrParAng(&datafile, NULL, 1, application.verbose_state) != 3) {
	    printerror(application.verbose_state.debug, "ERROR pcal (%s): Removal of parallactic angle correction failed.", filename_ptr);
	    return 0;
	  }
	}
	//      }else { 
	/*
	// Ref freq is important. In my mtm function it is always set to infinity first, so do that here as well. For solutions from for instance PSRCHIVE this is I assume not always true.
	if(receiversolution.freq_ref >= 0) {
	  printwarning(application.verbose_state.debug, "WARNING pcal (%s): The reference frequency of the solution is not infinite. This probably means the correction is not valid.", filename_ptr);
	}
	if(datafile.freq_ref >= 0) {
	  printwarning(application.verbose_state.debug, "WARNING pcal (%s): The reference frequency will be set to infinity before the receiver model can be applied.", filename_ptr);
	  if(preprocess_changeRefFreq(&datafile, -1, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pcal (%s): Changing reference frequency failed.", filename_ptr);
	    return 0;
	  }
	}
	if(datafile.isDeFarad == 0) {
	  printwarning(application.verbose_state.debug, "WARNING pcal (%s): The effect of Faraday rotation will be removed from the data file before the receiver model can be applied.", filename_ptr);
	  if(preprocess_deFaraday(&datafile, 0, 0, 0, NULL, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pcal (%s): de-Faraday rotation failed.", filename_ptr);
	    return 0;
	  }
	}
	if(datafile.isDePar == 0) {
	  printwarning(application.verbose_state.debug, "WARNING pcal (%s): The effect of parallactic angle rotation will be removed from the data file first before the receiver model can be applied.", filename_ptr);
	  if(preprocess_corrParAng(&datafile, NULL, 0, application.verbose_state) != 3) {
	    printerror(application.verbose_state.debug, "ERROR pcal (%s): Removal of parallactic angle correction failed.", filename_ptr);
	    return 0;
	  }
	}
	*/
	//      }




      /* Apply calibration */
      if(applyReceiverModel(&datafile, receiversolution, apply_reconstruct, apply_ignoreGain, ignore_dodgy_channels, apply_normalise, application.verbose_state) == 0) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR pcal: Applying receiver model failed.");
	return 0;
      }

      cleanPSRData(&datafileout, application.verbose_state);
      copy_params_PSRData(datafile, &datafileout, application.verbose_state);
      if(!openPSRData(&datafileout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	return 0;
      if(!writeHeaderPSRData(&datafileout, argc, argv, application.history_cmd_only, application.verbose_state))
	return 0;
      //      appendHistoryLine(datafileout, argc, argv, application.verbose_state);
      if(writePSRData(&datafileout, datafile.data, application.verbose_state) == 0) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR pcal: Cannot write data");
	return 0;
      }
      
      
      closePSRData(&receiversolution, 0, application.verbose_state);
      receiversolution_used = 0;
      closePSRData(&datafileout, 0, application.verbose_state);
      printf("Calibrated data written to %s\n", outputname);
    }


    if(apply_file || mtm_file) {
      //      if(!show_file) 
      closePSRData(&datafile, 0, application.verbose_state);
    }
    if(receiversolution_used) {
      closePSRData(&receiversolution, 0, application.verbose_state);
      receiversolution_used = 0;
      //      printf("CLOSE RECEIVERSOLUTION\n");
    }
    if(mtm_sequencenr == -1) {
      printf("\nSuggestion for launching jobs:\n\n");
      for(i = 0; i < mtm_totsequencenr; i++) {
	if(i != 0)
	  printf("sleep 2s\n");
	printf("xterm -geometry 160x20+0+0 -e 'nice +19 pcal -mtmrun \"%s %ld\" -wait", mtm_dumpfilename, i+1);
	for(index = 1; index < argc; index++) {
	  if(strcmp(argv[index], "-mtmprep") == 0)
	    index++;
	  else
	    printf(" \"%s\"", argv[index]);
	}
	printf("' &\n");
      }
      printf("\n");
    }
    first = 0;
    filenumber++;
  }  // End of loop over input files


  datafile_definition *solution_to_write, solution_analytic_version;
  if(combine_file) {
    long nrfiles = numberInApplicationFilenameList(&application, argv, application.verbose_state);
    int ret;
    if(fitphase == 0 && fitdiffgain == 0 && fitgain == 0 && fitell == 0 && fitor == 0) {
      ret = combineReceiverModels(receiversolution_collection, nrfiles, COMBINERECEIVERMODELS_CHI2, &receiversolution_combined, application.verbose_state);
      solution_to_write = &receiversolution_combined;   // Normally the combined dataset is written
    }else { // If fitting, then append the solutions for fitting
      ret = combineReceiverModels(receiversolution_collection, nrfiles, COMBINERECEIVERMODELS_APPEND, &receiversolution_combined, application.verbose_state);
      // Add 1 receiver solution, so is effectively just a copy
      ret *= combineReceiverModels(receiversolution_collection, 1, COMBINERECEIVERMODELS_APPEND, &solution_analytic_version, application.verbose_state);
      solution_to_write = &solution_analytic_version;   // Except when data is replaced with fit, then the data in first solution is replaced with fit
    }
    if(ret == 0) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR pcal: Cannot combine receiver solutions.");
      return 0;
    }


    if(fitphase || fitdiffgain || fitgain || fitell || fitor) {
      analyticReceiverSolution_def analyticsolution;
      if(fitReceiverModel(receiversolution_combined, fitgain, fitdiffgain, fitphase, phasemodel, diffgainmodel, fitell, fitell, fitor, fitor, rmsSigmaLimit, &analyticsolution, application.verbose_state) == 0) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR pcal: Fitting of an analytic solution to receiver model failed.");
	return 0;
      }
      printf("The following analytic solution was fitted:\n");
      printAnalyticReceiverSolution(analyticsolution, 2);

      if(makeReceiverModelAnalytic(&solution_analytic_version, &analyticsolution, application.verbose_state) == 0) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR pcal: Replacing receiver solution with fit failed.");
	return 0;
      }

    }

  }

  if(show_file && show_file_overlay) {   // If overlaying solution, plot them all at once after loop
    // Show receiver solution
    if(show_file == 2)
      neat = 1;
    else if(show_file == 3)
      neat = 2;
    else
      neat = 0;
    /*
    if(havemodel) {
      if(showReceiverModel(&receiversolution, 1, neat, show_fixscale, &analyticsolution, show_file_device, first, last, show_file_overlay, title, showmodel, application.verbose_state) == 0) {
    	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR pcal: Showing receiver model failed.");
	return 0;
      }
    }else {
    */
    long nrfiles = numberInApplicationFilenameList(&application, argv, application.verbose_state);
    if(showReceiverModel(receiversolution_collection, nrfiles, neat, show_fixscale, NULL, show_file_device, 0, 0, show_file_overlay, title, showmodel, application.verbose_state) == 0) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR pcal: Showing receiver model failed.");
      return 0;
    }
      //    }
  }



  if(combine_file) {
    if(!openPSRData(solution_to_write, argv[combine_file], FITS_format, 1, 0, 0, application.verbose_state))
      return 0;
    if(!writeHeaderPSRData(solution_to_write, argc, argv, application.history_cmd_only, application.verbose_state))
      return 0;
    //    appendHistoryLine(receiversolution_combined, argc, argv, application.verbose_state);
    if(writePSRData(solution_to_write, solution_to_write->data, application.verbose_state) == 0) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR pcal: Cannot write data");
      return 0;
    }
    closePSRData(solution_to_write, 0, application.verbose_state);
    if(fitphase || fitdiffgain || fitgain || fitell || fitor) {
      closePSRData(&solution_analytic_version, 0, application.verbose_state);
    }
  }

  if(receiversolution_collection != NULL) {
    long solnr, nrfiles = numberInApplicationFilenameList(&application, argv, application.verbose_state);
    for(solnr = 0; solnr < nrfiles; solnr++) {
      closePSRData(&(receiversolution_collection[solnr]), 0, application.verbose_state);
    }
    free(receiversolution_collection);
  }
    
  double ttotal;
  tend = clock();
  ttotal = (tend-tstart)/(double)CLOCKS_PER_SEC;
  if(application.verbose_state.verbose)
    printf("Program took %.2lf seconds = %.2lf minutes = %.2lf hours to run\n", ttotal, ttotal/60.0, ttotal/3600.0);

  if(dowait) {
    fflush(stdout);
    fprintf(stderr, "Press a key to quit\n");
    pgetch();
  }
  terminateApplication(&application);
  return 0;
}
