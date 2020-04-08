#include <math.h>
#include <stdio.h>
#include <string.h>
#include "psrsalsa.h"

void SHOWREVISIONINFO_prog() {
#include "fitvonMises.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);

void prescale(vonMises_collection_definition *components, datafile_definition fin, verbose_definition verbose);
int plotstuff(datafile_definition fin, vonMises_collection_definition components, float baseline, float *profile, int xleft, int xright, char *title, char *plotDevice, pulselongitude_regions_definition *regions, int limitonpulse, verbose_definition verbose);

void print_shapepar_help()
{
  printf("Possible shape parameters:\n");
  printf("w10                Width at 10%% of the maximum\n");
  printf("w25                Width at 25%% of the maximum\n");
  printf("w50                Width at 50%% of the maximum\n");
  printf("w75                Width at 75%% of the maximum\n");
  printf("w90                Width at 90%% of the maximum\n");
  printf("w10_maxamp         Width at 10%% of the maximum for that vonMises component with the maximum amplitude\n");
  printf("w25_maxamp         Width at 25%% of the maximum for that vonMises component with the maximum amplitude\n");
  printf("w50_maxamp         Width at 50%% of the maximum for that vonMises component with the maximum amplitude\n");
  printf("w75_maxamp         Width at 75%% of the maximum for that vonMises component with the maximum amplitude\n");
  printf("w90_maxamp         Width at 90%% of the maximum for that vonMises component with the maximum amplitude\n");
  printf("peakphase          Phase of the main peak\n");
  printf("peakamp            Amplitude of the main peak\n");
  printf("peaksep            Seperation in phase between the main peak and a smaller peak. The offset interval from the main peak where a sub-peak is searched for needs to be provided, e.g. \"peaksep 0.1 0.2\" will look for a sub-peak between 0.1 and 0.2 phase after the main peak.\n");
  printf("peakampratio       As peaksep, but calculate ration of amplitudes.\n");
  printf("peakampratio_reci  Reciprocal (1 over) of peakampratio.\n");
}


int main(int argc, char **argv)
{
  datafile_definition fin, fout;
  int i, j, index, educatedguess, precision, writemodel, maxnrcomponents, refine, fitbaseline, show, shape_parameter, shape_parameter_nofit;
  int refine_fixamp, refine_fixwidth, refine_fixphase, refine_fixrelamp, refine_fixrelphase, output_ascii, output_stats, avoid_neg_components, doprescale, limitonpulse;
  long shape_parameter_error_nr_itt;
  float *profile, ymin, ymax, baseline, sigmalimit, shape_parameter_error_maxdev;
  char *filename_ptr;
  vonMises_collection_definition components;
  psrsalsaApplication application;
  char outputname[MaxFilenameLength];
  double shape_parameter_precision, shapepar_aux[2];

  initApplication(&application, "fitvonMises", "[options] inputfile(s) (note: removing baseline should give better results)");
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
  application.switch_rebin = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_iformat = 1;
  application.switch_formatlist = 1;
  application.switch_rot= 1;
  application.switch_rotdeg = 1;
  application.switch_ext = 1;
  application.switch_output = 1;
  application.switch_showrevision = 1;
  application.switch_nocounters = 1;
  application.switch_FSCR = 1;
  application.switch_TSCR = 1;
  application.switch_stokes = 1; 
  application.switch_onpulse = 1; 
  application.switch_onpulsegr = 1; 
  application.switch_debase = 1; 
  application.switch_polselect = 1; 
  application.switch_autoonpulse = 1; 
  application.switch_device = 1;
  application.switch_scale = 1;
#ifdef HAVEPGPLOT2 
  strcpy(application.pgplotdevice, "/pw");
#else
  strcpy(application.pgplotdevice, "/xs");
#endif
  strcpy(application.outputname, "model.txt");

  writemodel = 0;
  educatedguess = 1;
  precision = 2;
  sigmalimit = 5;
  maxnrcomponents = 10;
  refine = 0;
  refine_fixamp = 0;
  refine_fixwidth = 0;
  refine_fixphase = 0;
  refine_fixrelphase = 0;
  refine_fixrelamp = 0;
  fitbaseline = 1;
  show = 0;
  output_ascii = 0;
  output_stats = 0;
  avoid_neg_components = 0;
  shape_parameter = 0;
  shape_parameter_error_nr_itt = 0;
  shape_parameter_precision = 1e-6;
  shape_parameter_nofit = 0;
  doprescale = 0;
  limitonpulse = 0;
  shape_parameter_error_maxdev = -1;

  if(argc < 2) {
    printApplicationHelp(&application);
    printf("Options related to precision of fitting:\n");
    printf("  -max                 Define the maximum number of components in the model.\n");
    printf("                       Default is %d.\n", maxnrcomponents);
    printf("  -notsmart            Is a lot slower, not necessarily more precise\n");
    printf("  -precision           The higher the number, the more finer the trials are.\n");
    printf("                       Default is %d.\n", precision);
    printf("  -sigmalimit          The lower the number, the more components will be fitted.\n");
    printf("                       Default is %.1f.\n", sigmalimit);
    printf("\nOther options related to fitting:\n");
    printf("  -limitonpulse        Limit the fitting to the onpulse region (chi2\n");
    printf("                       optimisation is limited to onpulse region, and the\n");
    printf("                       component centres are forced inside the onpulse region).\n");
    printf("                       Use this option together with -onpulse and/or -onpulsegr.\n");
    printf("  -nofitbaseline       Do not fit for an overall baseline.\n");
    printf("  -noneg               Avoid components with negative amplitudes.\n");
    printf("  -refine filename     Instead of making a model from scratch, refine this model\n");
    printf("  -refine_fixamp       With -refine, do not allow components to change amplitude\n");
    printf("  -refine_fixphase     With -refine, do not allow components to move in phase\n");
    printf("  -refine_fixrelamp    With -refine, do not allow components to change amplitude\n");
    printf("                       relative to each other\n");
    printf("  -refine_fixrelphase  With -refine, do not allow components to move relative\n");
    printf("                       to each other\n");
    printf("  -refine_fixwidth     With -refine, do not allow components to move in width\n");
    printf("  -refine_prescale     Before refining the model, scale the model first to fit\n");
    printf("                       the data (baseline should be subtracted from the data).\n");
    printf("  -shapepar X          Derive shape parameter X. Run without X to get a list.\n");
    printf("                       The shape parameter is derived from a fit to data (maybe\n");
    printf("                       using -refine or fitting an analytic template from\n");
    printf("                       scratch). Use -shapepar_nofit if fitting is undesired.\n");
    printf("                       For some shape parameters extra input is required.\n");
    printf("                       Example: -shapepar \"peaksep 0.1 0.2\"\n");
    printf("  -shapeparerr N       Bootstrap the errorbars on the shape parameter with this\n");
    printf("                       number of itterations.\n");
    printf("  -shapeparerr_maxdev M  Error itterations which deviate more than M with respect\n");
    printf("                       to the initially determined shape parameter are excluded.\n");
    printf("  -shapeparprec N      Set the precision of the shape parameter determination.\n");
    printf("                       Default=%.2e.\n", shape_parameter_precision);
    //    printf("  -w                   Write out model (use -output or -ext)\n");
    printf("\nOther actions not related to fitting:\n");
    printf("  -ascii nbin          The input files are von-Mises models. They are written\n");
    printf("                       out as ascii files with nbins bins. Their name is\n");
    printf("                       defined with -output or -ext.\n");
    printf("  -shapepar_nofit X    As -shapepar, but the input are models rather than data.\n");
    printf("  -show filename       Show this model, together with a data file provided at\n");
    printf("                       the end of the command line.\n");
    printf("  -show_prescale       Similar to -refine_prescale.\n");
    printf("  -stats               The input files are von-Mises models. Some statistics\n");
    printf("                       of the model are shown.\n");
    terminateApplication(&application); 
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-precision") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &precision, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-nofitbaseline") == 0 || strcmp(argv[i], "-fitnobaseline") == 0) {
	fitbaseline = 0;
      }else if(strcmp(argv[i], "-noneg") == 0) {
	avoid_neg_components = 1;
      }else if(strcmp(argv[i], "-limitonpulse") == 0) {
	limitonpulse = 1;
      }else if(strcmp(argv[i], "-show") == 0) {
	show = i+1;
        i++;
      }else if(strcmp(argv[i], "-ascii") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &output_ascii, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	if(output_ascii <= 0) {
	  printerror(application.verbose_state.debug, "fitvonMises: Error parsing option '%s %s' - a positive integer needs to be provided.", argv[i], argv[i+1]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-stats") == 0) {
	output_stats = 1;
      }else if(strcmp(argv[i], "-refine") == 0) {
	refine = i+1;
        i++;
      }else if(strcmp(argv[i], "-sigmalimit") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &sigmalimit, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-max") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &maxnrcomponents, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcasecmp(argv[i], "-shapeparerr") == 0 || strcasecmp(argv[i], "-shapeparerror") == 0 || strcasecmp(argv[i], "-shapepar_err") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%ld", &shape_parameter_error_nr_itt, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcasecmp(argv[i], "-shapeparerr_maxdev") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%f", &shape_parameter_error_maxdev, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcasecmp(argv[i], "-shapeparprec") == 0 || strcasecmp(argv[i], "-shapeparprecision") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%lf", &shape_parameter_precision, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-shapepar") == 0 || strcmp(argv[i], "-shapepar_nofit") == 0) {
	char shapepar_string[101];
	int ret;
	if(strcmp(argv[i], "-shapepar_nofit") == 0) {
	  shape_parameter_nofit = 1;
	}
	ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%100s %lf %lf", shapepar_string, &(shapepar_aux[0]), &(shapepar_aux[1]), NULL);
	if(ret == 0) {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot parse '%s' option.", argv[i]);
	  print_shapepar_help();
	  terminateApplication(&application); 
	  return 0;
	}
	if(strcasecmp(shapepar_string, "peakphase") == 0) {
	  shape_parameter = SHAPEPAR_PEAKPHASE;
	}else if(strcasecmp(shapepar_string, "peakamp") == 0) {
	  shape_parameter = SHAPEPAR_PEAKAMP;
	}else if(strcasecmp(shapepar_string, "w10") == 0) {
	  shape_parameter = SHAPEPAR_W10;
	}else if(strcasecmp(shapepar_string, "w25") == 0) {
	  shape_parameter = SHAPEPAR_W25;
	}else if(strcasecmp(shapepar_string, "w50") == 0) {
	  shape_parameter = SHAPEPAR_W50;
	}else if(strcasecmp(shapepar_string, "w75") == 0) {
	  shape_parameter = SHAPEPAR_W75;
	}else if(strcasecmp(shapepar_string, "w90") == 0) {
	  shape_parameter = SHAPEPAR_W90;
	}else if(strcasecmp(shapepar_string, "w10_maxamp") == 0) {
	  shape_parameter = SHAPEPAR_W10_MAXAMP;
	}else if(strcasecmp(shapepar_string, "w25_maxamp") == 0) {
	  shape_parameter = SHAPEPAR_W25_MAXAMP;
	}else if(strcasecmp(shapepar_string, "w50_maxamp") == 0) {
	  shape_parameter = SHAPEPAR_W50_MAXAMP;
	}else if(strcasecmp(shapepar_string, "w75_maxamp") == 0) {
	  shape_parameter = SHAPEPAR_W75_MAXAMP;
	}else if(strcasecmp(shapepar_string, "w90_maxamp") == 0) {
	  shape_parameter = SHAPEPAR_W90_MAXAMP;
	}else if(strcasecmp(shapepar_string, "peaksep") == 0) {
	  shape_parameter = SHAPEPAR_PEAKSEPPHASE;
	}else if(strcasecmp(shapepar_string, "peakampratio") == 0) {
	  shape_parameter = SHAPEPAR_PEAKAMPRATIO;
	}else if(strcasecmp(shapepar_string, "peakampratio_reci") == 0) {
	  shape_parameter = SHAPEPAR_PEAKAMPRATIO_RECI;
	}else {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises: '%s' not recognized as a shape parameter.", shapepar_string);
	  print_shapepar_help();
	  return 0;
	}
	if(shape_parameter == SHAPEPAR_PEAKSEPPHASE || shape_parameter == SHAPEPAR_PEAKAMPRATIO || shape_parameter == SHAPEPAR_PEAKAMPRATIO_RECI) {
	  // Check if expected auxilary parameters are provided
	  if(ret != 3) {
	    printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot parse '%s' option. Expected two additional values to be specified for this shape parameter.", argv[i]);
	    print_shapepar_help();
	    return 0;
	  }
	}else {
	  if(ret != 1) {
	    printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot parse '%s' option. Expected no additional values to be specified for this shape parameter.", argv[i]);
	    print_shapepar_help();
	    return 0;
	  }
	}
	i++;
	//      }else if(strcmp(argv[i], "-w") == 0) {
	//	writemodel = 1;
      }else if(strcmp(argv[i], "-notsmart") == 0) {
	educatedguess = 0;
      }else if(strcmp(argv[i], "-refine_fixamp") == 0) {
	refine_fixamp = 1;
      }else if(strcmp(argv[i], "-refine_fixwidth") == 0) {
	refine_fixwidth = 1;
      }else if(strcmp(argv[i], "-refine_fixphase") == 0) {
	refine_fixphase = 1;
      }else if(strcmp(argv[i], "-refine_fixrelphase") == 0) {
	refine_fixrelphase = 1;
      }else if(strcmp(argv[i], "-refine_fixrelamp") == 0) {
	refine_fixrelamp = 1;
      }else if(strcmp(argv[i], "-refine_prescale") == 0 || strcmp(argv[i], "-refineprescale") == 0 || strcmp(argv[i], "-show_prescale") == 0 || strcmp(argv[i], "-showprescale") == 0) {
	doprescale = 1;
      }else {
	/* If the option is not recognized, assume it is a filename */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises: Unknown option '%s', run command without command-line options to get help.", argv[i]);
	  terminateApplication(&application); 
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }

  if(refine_fixrelphase && refine_fixphase) {
    printerror(application.verbose_state.debug, "ERROR fitvonMises: You cannot use -refine_fixphase and -refine_fixrelphase simultaneously.");
    return 0;
  }
  if(refine_fixamp && refine_fixrelamp) {
    printerror(application.verbose_state.debug, "ERROR fitvonMises: You cannot use -refine_fixrelamp and -refine_fixamp simultaneously.");
    return 0;
  }

  if(shape_parameter == 0 && shape_parameter_nofit) {
    printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot use the -shapepar_nofit without specifying a shape parameter with the -shapepar option.");
    return 0;
  }
  if(shape_parameter_error_nr_itt && shape_parameter_nofit) {
    printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot use the -shapeparerr when the -shapepar_nofit is used.");
    return 0;
  }

  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR fitvonMises: No input files specified");
    return 0;
  }

  //  if(writemodel == 0) {
  for(i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-output") == 0 || strcmp(argv[i], "-ext") == 0) {
      //	printwarning(application.verbose_state.debug, "WARNING fitvonMises: Command line suggests you want to write out model. Please use the -w option");
      writemodel = 1;
    }
  }
//}
  
  while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {

    if(application.verbose_state.verbose) {
      printf("fitvonMises: processing %s\n", filename_ptr);      
    }
    baseline = 0;
    /* Clear struct */
    if(output_ascii == 0 && output_stats == 0 && shape_parameter_nofit == 0) { // Input files are data
      //      cleanPSRData(&fin, application.verbose_state);
      /* Open input file and read header */
      if(!openPSRData(&fin, filename_ptr, application.iformat, 0, 1, 0, application.verbose_state)) {
	printerror(application.verbose_state.debug, "ERROR fitvonMises (%s): Reading pulsar data failed. Note that this data is not expected to be an analytic template.", filename_ptr);
	return 0;
      }
      profile = calloc(fin.NrBins, sizeof(float));
    }else {  // Input files are von-Mises models
      if(readVonMisesModel(filename_ptr, &components, application.verbose_state) == 0) {
	return 0;
      }
      profile = calloc(output_ascii, sizeof(float));
    }
    if(profile == NULL) {
      printerror(application.verbose_state.debug, "fitvonMises: Memory allocation error.");
      return 0;
    }

    if(output_ascii) {
      cleanPSRData(&fout, application.verbose_state);
      fout.NrSubints = 1;
      fout.NrBins = output_ascii;
      fout.NrPols = 1;
      fout.NrFreqChan = 1;
      fout.format = PSRCHIVE_ASCII_format;
      calcVonMisesProfile(&components, fout.NrBins, profile, 0, 0);
      if(getOutputName(&application, filename_ptr, outputname, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR fitvonMises (%s): Changing filename failed", filename_ptr);
	return 0;
      }
      if(openPSRData(&fout, outputname, fout.format, 1, 0, 0, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR fitvonMises: Error opening file.");
	return 0;
      }
      if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR fitvonMises: Error writing header.");
	return 0;
      }
      if(application.verbose_state.verbose) {
	printf("Writing out data\n");    
      }
      if(writePSRData(&fout, profile, application.verbose_state) != 1) {
	printerror(application.verbose_state.debug, "ERROR fitvonMises: Error writing data.");
	return 0;
      }
      closePSRData(&fout, 0, application.verbose_state);
      printf("Done                       \n");
    }else if(output_stats) {
      long n;
      double area, w50, w10;
      for(n = 0; n < components.nrcomponents; n++) {
	printf("Component %ld\n", n+1);
	//   y = exp((cos(2.0*M_PI*(phase-centre-shift))-1.0)*concentration) * height;
	printf("  Phase comp. %ld:     %f phase = %f degrees = %f radians\n", n+1, components.centre[n], 360.0*components.centre[n], 2.0*M_PI*components.centre[n]);
	printf("  Max ampl. comp. %ld: %e intensity\n", n+1, components.height[n]);
	printf("  Min ampl. comp. %ld: %e intensity\n", n+1, components.height[n]*exp(-2.0*components.concentration[n]));
	area = integrateVonMisesFunction2(components.concentration[n], components.height[n]);
	printf("  Area comp. %ld:      %lf intensity*radians or %lf intensity*degrees\n", n+1, area, area*180.0/M_PI);
	w50 = widthVonMisesFunction2(components.concentration[n], 0.5);
	printf("  FWHM comp. %ld:      %lf degrees = %lf phase = %lf radians\n", n+1, 180.0*w50/M_PI, w50/(2.0*M_PI), w50);
	w10 = widthVonMisesFunction2(components.concentration[n], 0.1);
	printf("  W10 comp. %ld:       %lf degrees = %lf phase = %lf radians\n", n+1, 180.0*w10/M_PI, w10/(2.0*M_PI), w10);
      }
      printf("Note: the FWHM/W10 calculation ignores that the minimum is not necessarily zero for wide components.\n");
    }else if(shape_parameter && shape_parameter_nofit) {
      double measurement;
      calcVonMisesProfile_shape_parameter(&components, 0.0, shape_parameter_precision, shape_parameter, shapepar_aux, &measurement, application.verbose_state);
      print_shape_par(stdout, 1, shape_parameter, measurement, -1);
    }else {               // Input files are data
      application.dotscr = fin.NrSubints;  // Make sure the data is fully scrunched
      application.dofscr = fin.NrFreqChan;
      for(i = 1; i < argc; i++) {
	if(strcmp(argv[i], "-header") == 0) {
	  printwarning(application.verbose_state.debug, "WARNING: If using the -header option, be aware it applied BEFORE the preprocessing.");
	}
      }
      if(preprocessApplication(&application, &fin) == 0)
	return 0;

      if(show) {
	if(readVonMisesModel(argv[show], &components, application.verbose_state) == 0) {
	  return 0;
	}
	if(doprescale) {
	  prescale(&components, fin, application.verbose_state);
	}
      }else {
	if(refine) {
	  if(readVonMisesModel(argv[refine], &components, application.verbose_state) == 0) {
	    return 0;
	  }
	  if(doprescale) {
	     prescale(&components, fin, application.verbose_state);
	  }
	  if(fitvonmises_refine_model(fin, &components, fitbaseline, avoid_neg_components, &baseline, refine_fixamp, refine_fixwidth, refine_fixphase, refine_fixrelamp, refine_fixrelphase, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR fitvonmises (%s): fitvonmises_refine_model failed", filename_ptr);
	    return 0;
	  }
	}else {
	  /* Return 1 is success */
	  char *plot_device;
#ifdef HAVEPGPLOT2 
	  plot_device = "/pw";
#else
	  plot_device = "/xs";
#endif

	  pulselongitude_regions_definition *onpulse;
	  if(limitonpulse) {
	    onpulse = &(application.onpulse);
	  }else {
	    onpulse = NULL;
	  }
	  if(fitvonmises(fin, &components, maxnrcomponents, sigmalimit, educatedguess, fitbaseline, precision, avoid_neg_components, &baseline, onpulse, application.verbose_state, 2*application.verbose_state.verbose, plot_device) == 0) {
	    printerror(application.verbose_state.debug, "ERROR fitvonMises (%s): fitvonmises failed", filename_ptr);
	    return 0;
	  }
	}
      }
      
      ymax = fin.data[0];
      ymin = fin.data[0];
      for(i = 0; i < fin.NrBins; i++) {
	if(fin.data[i] > ymax)
	  ymax = fin.data[i];
	if(fin.data[i] < ymin)
	  ymin = fin.data[i];
      }
      if(ymax < 0 || ymin > 0 || fabs(ymax) < fabs(ymin)) {
	printwarning(application.verbose_state.debug, "\nWARNING fitvonMises (%s): Baseline does not appear to be subtracted, which might affect the fitting.\n", filename_ptr);
      }
      if(application.verbose_state.verbose || show) {
	if(plotstuff(fin, components, baseline, profile, 0, fin.NrBins-1, filename_ptr, application.pgplotdevice, &(application.onpulse), limitonpulse, application.verbose_state) == 0)
	  return 0;
	calcVonMisesProfile(&components, fin.NrBins, profile, 0, 0);
	pulselongitude_regions_definition boundaries;
	if(initPulselongitudeRegion(&boundaries, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises: Cannot find boundaries of profile.");
	  return 0;
	}

	find_boundaries(profile, fin.NrBins, 0.01*ymax, &boundaries);
	
	for(i = 0; i < boundaries.nrRegions; i++) {
	  j = pgplot_device_type(application.pgplotdevice, application.verbose_state); 
	  if((j >= 1 && j <= 2) || j == 100) {
	    printf("Press a key to continue %d/%d\n", i+1, boundaries.nrRegions);
	    pgetch();
	  }
	  if(boundaries.bins_defined[i] == 0) {
	    printerror(application.verbose_state.debug, "ERROR fitvonMises: region is not defined in bins");
	    return 0;
	  }
	  //	  if(plotstuff(fin, components, baseline, profile, application.onpulse.left_bin[i]-fin.NrBins/30, application.onpulse.right_bin[i]+fin.NrBins/30, filename_ptr, application.pgplotdevice, &(application.onpulse), limitonpulse, application.verbose_state) == 0)
	  if(plotstuff(fin, components, baseline, profile, boundaries.left_bin[i]-fin.NrBins/30, boundaries.right_bin[i]+fin.NrBins/30, filename_ptr, application.pgplotdevice, &(application.onpulse), limitonpulse, application.verbose_state) == 0) {
	    return 0;
	  }
	}
      }
      if(shape_parameter) {
	double rms;
	double measurement;
	calcVonMisesProfile_resid_rms(&components, fin.NrBins, fin.data, 0, &rms);
	if(application.verbose_state.verbose) {
	  printf("Found rms of the residual = %e\n", rms);
	}
	calcVonMisesProfile_shape_parameter(&components, 0.0, shape_parameter_precision, shape_parameter, shapepar_aux, &measurement, application.verbose_state);

	double stddev;
	stddev = -1;
	if(shape_parameter_error_nr_itt > 0) {
	  double av, sq;
	  long itt, nrsuccessfulitt;
	  av = sq = 0;
	  nrsuccessfulitt = 0;
	  for(itt = 0; itt < shape_parameter_error_nr_itt; itt++) {
	    vonMises_collection_definition components_itt;
	    datafile_definition clone;
	    memcpy(&components_itt, &components, sizeof(vonMises_collection_definition));
	    verbose_definition verbosedebug;
	    copyVerboseState(application.verbose_state, &verbosedebug);
	    if(application.verbose_state.debug == 0) {
	      verbosedebug.verbose = 0;
	    }
	    if(preprocess_addNoise(fin, &clone, rms, verbosedebug) != 1) {
	      printerror(application.verbose_state.debug, "ERROR fitvonmises (%s): Adding noise to observation failed", filename_ptr);
	      return 0;
	    }
	    if(fitvonmises_refine_model(clone, &components_itt, 1, 0, &baseline, 0, 0, 0, 0, 0, verbosedebug) == 0) {
	      printerror(application.verbose_state.debug, "ERROR fitvonmises (%s): fitvonmises_refine_model failed", filename_ptr);
	      return 0;
	    }
	    /*
	      if(application.verbose_state.verbose) {
	      if(plotstuff(fin, components_itt, baseline, profile, 0, fin.NrBins-1, filename_ptr, application.pgplotdevice, &(application.onpulse), limitonpulse, application.verbose_state) == 0)
	      return 0;
	      }
	    */

	    double measurement2;
	    calcVonMisesProfile_shape_parameter(&components_itt, 0.0, shape_parameter_precision, shape_parameter, shapepar_aux, &measurement2, application.verbose_state);
	    if(application.verbose_state.verbose) {
	      print_shape_par(stdout, 1, shape_parameter, measurement2, -1);
	    }
	    closePSRData(&clone, 0, application.verbose_state);

	    if(isnan(measurement2) == 0) {
	      //	      printf("dev = %f (maxdev = %f)\n", fabs(measurement2-measurement), shape_parameter_error_maxdev);
	      if(shape_parameter_error_maxdev <= 0 || fabs(measurement2-measurement) < shape_parameter_error_maxdev) {
		av += measurement2;
		sq += measurement2*measurement2;
		nrsuccessfulitt++;
	      }
	    }
	  }
	  sq /= (double)(nrsuccessfulitt);
	  av /= (double)(nrsuccessfulitt);
	  stddev = sqrt(sq-av*av);
	  //	  printf("XXXXX %lf (av=%lf sq=%lf)\n", stddev, av, sq);
	  if(nrsuccessfulitt != shape_parameter_error_nr_itt) {
	    printwarning(application.verbose_state.debug, "WARNING fitvonMises (%s): Not all itterations resulted in a successful measurement of the shape parameter (nr fails is %lf%%), so error might not be reliable.", filename_ptr, 100.0*(shape_parameter_error_nr_itt-nrsuccessfulitt)/(double)shape_parameter_error_nr_itt);
	  }
	}
	print_shape_par(stdout, 1, shape_parameter, measurement, stddev);
      }


      /* Free resources and close files */
      closePSRData(&fin, 0, application.verbose_state);
      if(writemodel) {
	if(getOutputName(&application, filename_ptr, outputname, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fitvonMises (%s): Changing filename failed", filename_ptr);
	  return 0;
	}
	if(writeVonMisesModel(outputname, &components, application.verbose_state) == 0)
	  return 0;
      }else {
	printf("Specify the -ext or -output option to write out the model to a file.\n");
      }
    }   // End of not ascii output if statement

    free(profile);
  }

  terminateApplication(&application);
  return 0;
}


void prescale(vonMises_collection_definition *components, datafile_definition fin, verbose_definition verbose)
{
  double area, area_data;
  long i;
  if(verbose.verbose) {
    printf("Pre-scaling model before refinement");
  }
  area = integrateVonMisesFunction(components);
  if(verbose.verbose) {
    printf("  Integral of model: %e intensity*radians\n", area);
  }
  area_data = 0;
  for(i = 0; i < fin.NrBins; i++) {
    area_data += fin.data[i];
  }
  area_data *= 2*M_PI/(double)fin.NrBins;
  double scale_factor;
  scale_factor = area_data/area;
  if(verbose.verbose) {
    printf("  Integral of data:  %e intensity*radians  (assuming baseline is subtracted)\n", area_data);
    printf("  Scale factor: %e\n", scale_factor);
  }
  for(i = 0; i < components->nrcomponents; i++) {
    components->height[i] *= scale_factor;
  }
}

int plotstuff(datafile_definition fin, vonMises_collection_definition components, float baseline, float *profile, int xleft, int xright, char *title, char *plotDevice, pulselongitude_regions_definition *regions, int limitonpulse, verbose_definition verbose)
{
  int i, n, color;
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);


  strcpy(pgplot_options.viewport.plotDevice, plotDevice);
  pgplot_options.viewport.dontclose = 1;
  strcpy(pgplot_options.box.xlabel, "Pulse longitude bin");
  strcpy(pgplot_options.box.ylabel, "Intensity");
  strcpy(pgplot_options.box.title, title);
  if(pgplotGraph1(&pgplot_options, fin.data, NULL, NULL, fin.NrBins, 0, fin.NrBins-1, 0, xleft, xright, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, verbose) == 0) {
    printerror(verbose.debug, "Cannot open plot device\n");
    return 0;
  }

  // Calculate residuals
  calcVonMisesProfile(&components, fin.NrBins, profile, 0, 0);
  for(i = 0; i < fin.NrBins; i++) {
    profile[i] = fin.data[i]-profile[i] + 2*baseline;
  }
  pgplot_options.viewport.dontopen = 1;
  pgplot_options.viewport.noclear = 1;
  pgplot_options.box.drawbox = 0;
  pgplot_options.box.drawtitle = 0;
  pgplot_options.box.drawlabels = 0;
  if(limitonpulse == 0) {
    regions = NULL;
    color = 3;
  }else {
    color = 15;
  }
  if(pgplotGraph1(&pgplot_options, profile, NULL, NULL, fin.NrBins, 0, fin.NrBins-1, 1, xleft, xright, 0, 0, 0, 1, 0, 1, 0, color, 1, regions, 3, verbose) == 0) {
    printerror(verbose.debug, "Cannot open plot device\n");
    return 0;
  }

  //int pgplotGraph1(pgplot_options_definition *pgplot, float *data, float *datax, float *sigma, int nrx, float xmin, float xmax, int dontsetranges, float xmin_show, float xmax_show, float ymin_show, float ymax_show, int forceMinZero, int hist, int noline, int linewidth, int pointtype, int color, int boxcolor, pulselongitude_regions_definition *regions, int onpulsecolor, verbose_definition verbose)


  // Draw individual components
  color = 4;
  for(n = 0; n < components.nrcomponents; n++) {
    for(i = 0; i < fin.NrBins; i++) {
      profile[i] = calcVonMisesFunction2(components.centre[n], components.concentration[n], components.height[n], i/(float)fin.NrBins, 0) + baseline;
    }
    if(pgplotGraph1(&pgplot_options, profile, NULL, NULL, fin.NrBins, 0, fin.NrBins-1, 1, xleft, xright, 0, 0, 0, 0, 0, 5, 0, color, 1, NULL, -1, verbose) == 0) {
      printerror(verbose.debug, "Cannot open plot device\n");
      return 0;
    }
  }

  // Draw overall fit
  calcVonMisesProfile(&components, fin.NrBins, profile, 0, 0);
  for(i = 0; i < fin.NrBins; i++) {
    profile[i] += baseline;
  }
  color = 2;
  if(pgplotGraph1(&pgplot_options, profile, NULL, NULL, fin.NrBins, 0, fin.NrBins-1, 1, xleft, xright, 0, 0, 0, 0, 0, 1, 0, color, 1, NULL, -1, verbose) == 0) {
    printerror(verbose.debug, "Cannot open plot device\n");
    return 0;
  }

  ppgclos();

  return 1;
}
