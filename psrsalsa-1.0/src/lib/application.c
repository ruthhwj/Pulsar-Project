//START REGION RELEASE
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include "psrsalsa.h"

int   internal_application_cmdline_FilenameList[MaxNrApplicationFilenames];
int   internal_application_cmdline_Nrfilenames;
int   internal_application_cmdline_CurFilename;
int   internal_application_filelist_CurFilename;
long  internal_application_filelist_Nrfilenames;
int   internal_application_filelist_fileopen;   // 0 = not yet opened, 1 = open
FILE *internal_application_filelist_fptr;
char *internal_application_filelist_filename;

//START REGION DEVELOP

void SHOWREVISIONINFO_DEFAULT()
{
  printf("  Revision information is udefined.\n");
}
void (*SHOWREVISIONINFO)(void);

//START REGION RELEASE
void initApplication(psrsalsaApplication *application, char *name, char *genusage)
{
  verbose_definition verbose;
  cleanVerboseState(&verbose);

  internal_application_cmdline_Nrfilenames = 0;
  internal_application_cmdline_CurFilename = 0;
  internal_application_filelist_CurFilename = 0;
  internal_application_filelist_Nrfilenames = 0;
  internal_application_filelist_fileopen = 0;
  strcpy(application->progname, name);
  application->genusage = malloc(strlen(genusage)+1);
  if(application->genusage == NULL) {
    printerror(0, "ERROR initApplication: Cannot allocate memory.");
    exit(0);
  }
  strcpy(application->genusage, genusage);
  application->switch_verbose = 0;
  application->switch_debug = 0;
  application->switch_nocounters = 0;
  cleanVerboseState(&(application->verbose_state));
  application->switch_formatlist = 0;
  application->switch_iformat = 0;
  application->iformat = -1;
  application->switch_oformat = 0;
  application->oformat = -1;
  application->switch_headerlist = 0;
  application->switch_header = 0;
  application->switch_onpulse = 0;
  application->switch_onpulsef = 0;
  if(initPulselongitudeRegion(&(application->onpulse), verbose) == 0) {
    printerror(0, "ERROR initApplication: Cannot initialise onpulse region.");
    exit(0);
  }
  application->switch_polselect = 0;
  application->polselectnr = -1;
  application->switch_itf = 0;
  application->itf = 0;
  application->switch_rebin = 0;
  application->dorebin = 0;
  application->switch_nread = 0;
  application->nread = -1;
  application->switch_nskip = 0;
  application->nskip = 0;
  application->switch_conshift = 0;
  application->doconshift = 0;
  application->switch_circshift = 0;
  application->docircshift = 0;
  application->switch_rot = 0;
  application->switch_rotdeg = 0;
  application->doshiftphase = 0;
  application->shiftPhase_cmdline = 0;
  application->shiftPhase = 0;
  application->switch_filelist = 0;
  application->filelist = 0;
  application->switch_device = 0;
  strcpy(application->pgplotdevice, "?");
  application->switch_tscr = 0;
  application->dotscr = 0;
  application->switch_tscr_complete = 0;
  application->tscr_complete = 0;
  application->switch_TSCR = 0;
  application->doTSCR = 0;
  application->switch_fscr = 0;
  application->dofscr = 0;
  application->switch_FSCR = 0;
  application->doFSCR = 0;
  application->switch_dedisperse = 0;
  application->do_dedisperse = 0;
  application->switch_deFaraday = 0;
  application->do_deFaraday = 0;
  application->switch_changeRefFreq = 0;
  application->newRefFreq = -100;
  application->switch_stokes = 0;
  application->dostokes = 0;
  application->switch_coherence = 0;
  application->docoherence = 0;
  application->switch_noweights = 0;
  application->noweights = 0;
  application->switch_useweights = 0;
  application->useweights = 0;
  application->switch_uniformweights = 0;
  application->uniformweights = 0;
  application->switch_scale = 0;
  application->doscale = 0;
  application->switch_debase = 0;
  application->dodebase = 0;
  application->switch_onpulsegr = 0;
  application->doonpulsegr = 0;
  application->switch_size = 0;
  application->windowwidth = -1;
  application->windowheight = -1;
  application->switch_macro = 0;
  application->macro_ptr = NULL;
  application->switch_noplotsubset = 0;
  application->do_noplotsubset = 0;
  application->switch_cmaplist = 0;
  application->switch_cmap = 0;
  application->cmap = PPGPLOT_HEAT;
  application->switch_deparang = 0;
  application->switch_insertparang = 0;
  application->do_parang_corr = 0;
  application->switch_norm = 0;
  application->normvalue = 1;
  application->do_norm = 0;
  application->switch_normglobal = 0;
  application->do_normglobal = 0;
  application->switch_clip = 0;
  application->do_clip = 0;
  application->clipvalue = 0;
  application->switch_fchan = 0;
  application->fchan_select = -1;
  application->switch_history_cmd_only = 0;
  application->history_cmd_only = 0;
  application->switch_fixseed = 0;
  application->fixseed = 0;
  application->switch_template = 0;
  application->template_specified = 0;
  application->switch_templatedata = 0;
  application->template_data_index = 0;
  cleanPSRData(&(application->template_file), application->verbose_state);
  application->switch_align = 0;
  application->doalign = 0;
  application->switch_blocksize = 0;
  application->blocksize = 0;
  application->switch_ext = 0;
  application->extension = NULL;
  application->switch_output = 0;
  sprintf(application->outputname, "output.dat");
  application->switch_shuffle = 0;
  application->doshuffle = 0;
  application->switch_rotateStokes = 0;
  application->nr_rotateStokes = 0;
  application->switch_libversions = 0;
  application->switch_forceUniformFreqLabelling = 0;
  application->dodebase_slope = 0;

//START REGION DEVELOP
  application->switch_rotateQU = 0;
  application->dorotateQU = 0;
  application->switch_rotateUV = 0;
  application->dorotateUV = 0;
  application->switch_subtractRVM = 0;
  application->dosubtractRVM = 0;
  application->switch_subtractAvPA = 0;
  application->dosubtractAvPA = 0;
  application->switch_deFaradayTable = 0;
  application->do_deFaradayTable = 0;
  application->deFaradayTable_filename = NULL;

//The following parameter is used by pmod
//START REGION RELEASE
  application->fzapMask = NULL;
  //START REGION DEVELOP
// The following is used by getNextFilenameFromList(), so needs to be defined
//START REGION RELEASE
  application->doautot = 0;
  application->switch_onpulse2 = 0;
  application->switch_onpulsef2 = 0;
  if(initPulselongitudeRegion(&(application->onpulse2), verbose) == 0) {
    printerror(0, "ERROR initApplication: Cannot initialise onpulse region.");
    exit(0);
  }
//START REGION DEVELOP

  application->switch_invertxy = 0;
  application->invertxy = 0;
  application->switch_invertfp = 0;
  application->switch_invertf = 0;
  application->invertfp = 0;
  application->invertf = 0;
  application->switch_invertfx = 0;
  application->invertfx = 0;
  application->switch_ppgplot = 0;
  application->ppgplot_name = NULL;
  application->switch_autoonpulse = 0;
  application->do_autoonpulse = 0;
  application->switch_gate = 0;
  application->doGate = 0;
  application->switch_fftsmooth = 0;
  application->do_fftsmooth = -1;
  application->switch_smooth = 0;
  application->do_smooth = -1;
  application->switch_fftzap = 0;
  application->do_fftzap = 0;
  application->switch_fft = 0;
  application->dofft = 0;
  application->switch_template_onpulse = 0;
  application->template_onpulse = 0;
  application->switch_centert = 0;
  application->docentert = 0;
  application->switch_autot = 0;
  application->switch_normRMS = 0;
  application->donormRMS = 0;
  application->switch_normRMSi = 0;
  application->donormRMSi = 0;
  application->switch_testsignal = 0;
  application->do_testsignal = 0;
  application->switch_useweightedfreq = 0;
  application->useweightedfreq = 0;
  application->switch_absweights = 0;
  application->absweights = 0;
  application->switch_checknan = 0;
  application->dochecknan = 0;
  application->switch_removenan = 0;
  application->doremovenan = 0;
  application->switch_checkinf = 0;
  application->docheckinf = 0;
  application->switch_removeinf = 0;
  application->doremoveinf = 0;
  application->switch_showrevision = 0;
  application->switch_swapcables = 0;
  application->doswapcables = 0;
  application->switch_dedeFaraday = 0;
  application->switch_dededisperse = 0;
  application->switch_rotslope = 0;
  application->dorotslope = 0;
  application->switch_zerodming = 0;
  application->dozerodming = 0;
  application->switch_zerodming_adv = 0;
  application->dozerodming_adv = 0;
  application->switch_debase_slope = 0;
//START REGION RELEASE
//START REGION DEVELOP

  /*
application->switch_ = 0;
  */
  SHOWREVISIONINFO = SHOWREVISIONINFO_DEFAULT;
  //START REGION RELEASE
}

// Releases some memory
void terminateApplication(psrsalsaApplication *application)
{
  free(application->genusage);
  freePulselongitudeRegion(&(application->onpulse));
  freePulselongitudeRegion(&(application->onpulse2));
//START REGION DEVELOP
  if(application->deFaradayTable_filename != NULL)
    free(application->deFaradayTable_filename);
//START REGION RELEASE
  closePSRData(&(application->template_file), 0, application->verbose_state);
 }

void printCitationInfo()
{
  printf("If you make use of PSRSALSA, please cite \"Weltevrede 2016, A&A, 590, A109\" and refer to the following website: https://github.com/weltevrede/psrsalsa\n");
}

// argv[argv_index] = string to be parsed, so argv[argv_index-1] is the name of the option
// If check_only != 0, no errors messages are generated, so you can "try out" if a given format would work or not
// When parsing a string, the max length (plus termination character) should be provides: e.g. %100s
// minrequestedparameters = minimum number of parameters that should be successfully parsed. This is useful if the user can either provide for example 1 or 2 parameters. If <= 0, all parameters should be read in.
// Return value: nr of parameters successfully parsed if no errors are generated
//               0 is returned when an error occured
int parse_command_string(verbose_definition verbose, int argc, char **argv, int argv_index, int check_only, int minrequestedparameters, char *format, ...)
{
  va_list args;
  long i, curargnr, nrarguments;
  void *ptr;

  va_start(args, format); // Retrieve a pointer to the first variable in ..., starting after the defined fixed variable format

  if(argv_index >= argc) {
    fflush(stdout);
    if(check_only == 0) {
      printerror(verbose.debug, "ERROR parse_command_string: the command line option to be parsed does not appear to be followed by input values/text.", format);
    }
    return 0;
  }

  if(verbose.debug) {
    printf("parse_command_string: start parsing the \"%s\" option (check_only=%d min_req_params=%d)\n", argv[argv_index-1], check_only, minrequestedparameters);
    printf("parse_command_string: parsing \"%s\" as format \"%s\"\n", argv[argv_index], format);
  }
  nrarguments = 0;
  for(i = 0; i < strlen(format); i++) {
    if(format[i] == '%') {
      nrarguments++;
    }
  }
  if(verbose.debug) {
    printf("parse_command_string: expected number of arguments = %ld\n", nrarguments);
  }
  if(minrequestedparameters <= 0) {
    minrequestedparameters = nrarguments;
  }
  if(nrarguments == 0) {
    if(check_only == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR parse_command_string: format string '%s' suggest no variables need to be parsed.", format);
    }
    return 0;
  }
  for(i = 0; i < nrarguments; i++) {
    ptr = va_arg(args, void *);
    //    printf("%p\n", ptr);
    if(ptr == NULL) {
      if(check_only == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR parse_command_string: %ld pointers are expected to be passed to this function, but at least one of them appears to be the NULL pointer, suggestive of a bug in the programme.", nrarguments);
      }
      return 0;
    }
  }
  ptr = va_arg(args, void *);  // The pointers are expected to be followed by the null pointer to terminate the list of variables
  //  printf("%p\n", ptr);
  if(ptr != NULL) {
    if(check_only == 0) {
      fflush(stdout);
      printerror(verbose.debug, "ERROR parse_command_string: %ld pointers are expected to be passed to this function followed by the NULL pointer. However, the NULL pointer is not detected, suggestive of a bug in the programme.", nrarguments);
    }
    return 0;
  }
  va_end(args); // In order to rewind, destroy structure and start again

  int nrwords;
  ptr = pickWordFromString(argv[argv_index], 1, &nrwords, 1, ' ', verbose);
  if(nrwords != nrarguments) {
    if(nrwords < minrequestedparameters) {
      if(check_only == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR parse_command_string: string \"%s\" appears to contain %d words, while the program expects at least %d/%ld options to be provided.", argv[argv_index], nrwords, minrequestedparameters, nrarguments);
	if(nrarguments > 1) {
	  printerror(verbose.debug, "ERROR parse_command_string: When multiple arguments are expected, provide these within quotes, e.g. -option \"value1 value2\".");
	}
      }
      return 0;
    }
    if(nrwords > nrarguments) {
      if(check_only == 0) {
	fflush(stdout);
	printerror(verbose.debug, "ERROR parse_command_string: string \"%s\" appears to contain %d words, while the program expects at most %ld options to be provided.", argv[argv_index], nrwords, nrarguments);
      }
      return 0;
    }
  }
  minrequestedparameters = nrwords; // So minrequestedparameters is the number of words that are now expected to be read in

  va_start(args, format); // Retrieve a pointer to the first variable in ..., starting after the defined fixed variable format

  // Parse format string and obtain all types
  curargnr = 0;
  int expecttype, expectnumber, handlednumber;  //, expectstring
  long maxsize;
  expecttype = 0;
  expectnumber = 0;
  maxsize = 0;
  //  expectstring = 0;
  for(i = 0; i < strlen(format); i++) {
    handlednumber = 0;
    if(expecttype) { // Check if it is really a number (e.g. %100s) rather than a type (e.g. ld).
      if(format[i] >= '0' && format[i] <= '9') {
	handlednumber = 1;
	//	expecttype = 0; // Is a number, so it is not a type (yet)
	//	expectstring = 1; // The next type should be a string since a size is defined
	if(expectnumber == 0) {
	  maxsize = format[i] - '0';
	  expectnumber = 1;
	}else {
	  maxsize *= 10;
	  maxsize += format[i] - '0';
	}
	//	printf("XXXXX %ld\n", maxsize);
      }
    }else {
      expectnumber = 0;
    }
    if(expecttype && handlednumber == 0) {
      char *word_ptr, *word;
      if(curargnr > minrequestedparameters) { // Reached end of provided arguments
	return minrequestedparameters;
      }
      word_ptr = pickWordFromString(argv[argv_index], curargnr, &nrwords, 1, ' ', verbose);
      word = malloc(strlen(word_ptr)+1);
      if(word == NULL) {
	if(check_only == 0) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR parse_command_string: Memory allocation error");
	}
	return 0;      
      }
      sscanf(word_ptr, "%s", word); // read word, since word_ptr points to start of word we're interested in followed by the others.
      if(format[i] == 'c') {  // Character
	if(verbose.debug) {
	  printf("parse_command_string: Parsing \"%s\" as a character\n", word);
	}
	if(strlen(word) != 1) {
	  if(check_only == 0) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR parse_command_string: Cannot parse '%s' as a character, since it is not 1 character in length", word);
	  }
	  free(word);
	  return 0;      
	}
	char *char_ptr;
	char_ptr = va_arg(args, char *);
	char_ptr[0] = word[0];
	if(verbose.debug) {
	  printf("parse_command_string: Parsing \"%s\" as character '%c'\n", word, char_ptr[0]);	
	}
      }else if(format[i] == 's') {  // String
	  // XXXXXXXXXXXXXXXXXXX
	if(verbose.debug) {
	  printf("parse_command_string: Parsing \"%s\" as a string with maximum size %ld\n", word, maxsize);
	}

	if(strlen(word) + 1 >= maxsize) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR parse_command_string: Command line option '%s' exceeds the maximum variable length %ld", word, maxsize);
	  if(maxsize == 0) {
	    printerror(verbose.debug, "ERROR parse_command_string: This might be a programming error where %%s is passed on to parse_command_string() rathern than for example %%100s.");
	  }
	  free(word);
	  return 0;
	}
	char *string_ptr;
	string_ptr = va_arg(args, char *);
	strncpy(string_ptr, word, maxsize);
	if(verbose.debug) {
	  printf("parse_command_string: Parsing \"%s\" as %s\n", word, string_ptr);	
	}
      }else if((format[i] == 'l' && format[i+1] == 'f') || format[i] == 'f') {  // Floating point (float or double)
	if(verbose.debug) {
	  if(format[i] == 'l')
	    printf("parse_command_string: Parsing \"%s\" as a double\n", word);
	  else
	    printf("parse_command_string: Parsing \"%s\" as a floating point\n", word);
	}

	float *float_ptr;
	double *double_ptr;
	char *endptr;
	if(format[i] == 'l') {
	  double_ptr = va_arg(args, double *);
	  *double_ptr = strtod(word, &endptr);
	}else {
	  float_ptr = va_arg(args, float *);
	  *float_ptr = strtod(word, &endptr);
	}
	if(endptr == word || endptr == NULL || endptr != word+strlen(word)) {
	  if(check_only == 0) {
	    fflush(stdout);
	    printerror(verbose.debug, "ERROR parse_command_string: Cannot parse '%s' as a floating point", word);
	  }
	  free(word);
	  return 0;
	}
	if(verbose.debug) {
	  if(format[i] == 'l')
	    printf("parse_command_string: Parsing \"%s\" as %lf\n", word, *double_ptr);	
	  else
	    printf("parse_command_string: Parsing \"%s\" as %f\n", word, *float_ptr);	
	}
	if(format[i] == 'l')
	  i++;
      }else if((format[i] == 'l' && format[i+1] == 'd') || format[i] == 'd') {  // long int (or int)
	if(verbose.debug) {
	  if(format[i] == 'l')
	    printf("parse_command_string: Parsing \"%s\" as a long int\n", word);
	  else
	    printf("parse_command_string: Parsing \"%s\" as a int\n", word);
	}
	long int *long_ptr;
	int *int_ptr;
	char *endptr;
	if(format[i] == 'l') {
	  long_ptr = va_arg(args, long int *);
	  *long_ptr = strtol(word, &endptr, 0);  // use default base depending on how number is formatted, e.g. 0x123.
	}else {
	  int_ptr = va_arg(args, int *);
	  *int_ptr = strtol(word, &endptr, 0);  // use default base depending on how number is formatted, e.g. 0x123.
	}
	if(endptr == word || endptr == NULL || endptr != word+strlen(word)) {
	  if(check_only == 0) {
	    fflush(stdout);
	    if(format[i] == 'l') {
	      printerror(verbose.debug, "ERROR parse_command_string: Cannot parse '%s' as a long int", word);
	    }else {
	      printerror(verbose.debug, "ERROR parse_command_string: Cannot parse '%s' as a long", word);
	    }
	  }
	  free(word);
	  return 0;
	}
	if(verbose.debug) {
	  if(format[i] == 'l')
	    printf("parse_command_string: Parsing \"%s\" as %ld\n", word, *long_ptr);	
	  else
	    printf("parse_command_string: Parsing \"%s\" as %d\n", word, *int_ptr);	
	}
	if(format[i] == 'l')
	  i++;
      }else {
	if(check_only == 0) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR parse_command_string: format string \"%s\" contains a unrecognized data type '%c', suggestive of a bug in the software.", format, format[i]);
	}
	free(word);
	return 0;      
      }
      expecttype = 0;
      free(word);
    }else if(handlednumber == 0) {
      if(format[i] == '%') {
	curargnr++;
	expecttype = 1;  // Expect next characted to be a type of variable
      }else if(format[i] == ' ' || format[i] == '\t') {  // space or tab - just skip
      }else {
	if(check_only == 0) {
	  fflush(stdout);
	  printerror(verbose.debug, "ERROR parse_command_string: format string \"%s\" does not appear of the expected format, indicating a bug.", format);
	}
	return 0;
      }
    }
  }

  va_end(args);
  return curargnr;
}


//START REGION RELEASE
void printApplicationHelp(psrsalsaApplication *application)
{
  fprintf(stdout, "%s %s\n", application->progname, application->genusage);
  if(application->switch_iformat || application->switch_oformat || application->switch_formatlist || application->switch_headerlist || application->switch_header || application->switch_filelist || application->switch_noweights || application->switch_useweights || application->switch_uniformweights || application->switch_history_cmd_only || application->switch_ext || application->switch_output || application->switch_forceUniformFreqLabelling
  //START REGION DEVELOP
     || application->switch_useweightedfreq || application->switch_absweights
//START REGION RELEASE
) {
    fprintf(stdout, "\nGeneral Input/Output options:\n");
    if(application->switch_filelist) {
      fprintf(stdout, "  -filelist file    Specify file with input file names: only one file name\n");
      fprintf(stdout, "                    per line, and they can be commented out with a #.\n");
    }
    if(application->switch_ext)
      fprintf(stdout, "  -ext ext          Specify output extension ext\n");
    if(application->switch_output) {
      fprintf(stdout, "  -output file      ");
      if(application->switch_ext)
	fprintf(stdout, "instead ");
      fprintf(stdout, "specify output filename, default is %s\n", application->outputname);
    }
    if(application->switch_formatlist)
      fprintf(stdout, "  -formatlist       Show supported file formats\n");
    if(application->switch_iformat)
      fprintf(stdout, "  -iformat id       Specify input format (e.g.  -iformat PSRFITS)\n");
    if(application->switch_oformat)
      fprintf(stdout, "  -oformat id       Specify output format (e.g. -oformat PSRFITS)\n");
    if(application->switch_headerlist)
      fprintf(stdout, "  -headerlist       Show available parameters for the -header option\n");
    if(application->switch_header) 
      fprintf(stdout, "  -header           Change header parameter, e.g. -header 'name J0123+4567'\n");
    if(application->switch_forceUniformFreqLabelling) {
      fprintf(stdout, "  -headerUFL        Force uniform frequency labelling (from channel to\n");
      fprintf(stdout, "                    channel and from subint to subint). This will allow more\n");
      fprintf(stdout, "                    processes to work. However, this is not done in any\n");
      fprintf(stdout, "                    intelegent way, so expect all frequency dependent\n");
      fprintf(stdout, "                    effects to go wrong!\n");
    }
  //START REGION DEVELOP
    if(application->switch_useweightedfreq) {
      fprintf(stdout, "  -useweightedfreq  Use weighted centre freq rather than the centre freq of\n");
      fprintf(stdout, "                    psrfits file (assuming it is constant throughout the data)\n");
    }
//START REGION RELEASE
    if(application->switch_noweights) 
      fprintf(stdout, "  -noweights        Ignore the weights in the PSRFITS input file\n");
    if(application->switch_noweights) 
      fprintf(stdout, "  -useweights       Force usage of the weights in the PSRFITS input file\n");
    if(application->switch_noweights) 
      fprintf(stdout, "  -uniformweights   Ignore weights in the PSRFITS input file, except if zero\n");
//START REGION DEVELOP
    if(application->switch_absweights) 
      fprintf(stdout, "  -absweights       Take absolute value of weights in the PSRFITS input file\n");
//START REGION RELEASE
    if(application->switch_history_cmd_only) {
      fprintf(stdout, "  -history_cmd_only Write the history without timestamp, hence re-running the\n");
      fprintf(stdout, "                    same command might result in identical files.\n");
    }
//START REGION RELEASE
  }
  if(application->switch_templatedata || application->switch_align || application->switch_template 
//START REGION DEVELOP
     || application->switch_centert || application->switch_autot
//START REGION RELEASE
     ) {
    fprintf(stdout, "\nGeneral template options:\n");
    if(application->switch_align) {
      fprintf(stdout, "  -align             Rotate all subints with the same amount to align data\n");
      fprintf(stdout, "                     with template\n");
    }
//START REGION DEVELOP
    if(application->switch_autot) 
      fprintf(stdout, "  -autot             Fit a template automatically\n");
    if(application->switch_centert) 
      fprintf(stdout, "  -centert           Make peak template phase 0 (use -align to rotate data)\n");
//START REGION RELEASE
    if(application->switch_template) 
      fprintf(stdout, "  -template file     Use mathematical template (von-Mises ascii file)\n");
    if(application->switch_templatedata) 
      fprintf(stdout, "  -templatedata file Use this data file as a template profile\n");
  }
//START REGION RELEASE
  if(application->switch_polselect || application->switch_rebin || application->switch_nread || application->switch_nskip || application->switch_conshift || application->switch_circshift || application->switch_rot || application->switch_rotdeg || application->switch_tscr || application->switch_TSCR || application->switch_tscr_complete || application->switch_fscr || application->switch_FSCR || application->switch_dedisperse || application->switch_deFaraday || application->switch_stokes || application->switch_coherence || application->switch_changeRefFreq || application->switch_scale || application->switch_debase || application->switch_deparang || application->switch_insertparang || application->switch_norm || application->switch_normglobal || application->switch_fchan || application->switch_blocksize || application->switch_shuffle  || application->switch_clip || application->switch_rotateStokes
  //START REGION DEVELOP
      || application->switch_rotateQU || application->switch_rotateUV || application->switch_invertxy || application->switch_invertfp || application->switch_invertfx || application->switch_gate || application->switch_fftsmooth || application->switch_smooth|| application->switch_fftzap || application->switch_fft || application->switch_normRMS || application->switch_normRMSi || application->switch_testsignal || application->switch_checknan || application->switch_removenan || application->switch_checkinf || application->switch_removeinf || application->switch_dedeFaraday || application->switch_dededisperse || application->switch_invertf || application->switch_swapcables || application->switch_rotslope || application->switch_subtractRVM  || application->switch_subtractAvPA  || application->switch_deFaradayTable || application->switch_zerodming || application->switch_zerodming_adv  || application->switch_debase_slope
//START REGION RELEASE
) {
    fprintf(stdout, "\nGeneral preprocess options:\n");
    if(application->switch_blocksize) 
      fprintf(stdout, "  -blocksize b    Only read in a multiple of b subintegrations.\n");
    if(application->switch_clip) {
      fprintf(stdout, "  -clip           Specify threshhold value above which the value clipped\n");
      fprintf(stdout, "                  Clipping happens if sample > threshold or < -threshold.\n");
    }
    if(application->switch_dedisperse) 
      fprintf(stdout, "  -dedisp         Dedisperse the data.\n");
    if(application->switch_deFaraday) {
      fprintf(stdout, "  -defarad        De-Faraday rotate the data.\n");
    }
//START REGION DEVELOP
    if(application->switch_deFaradayTable) {
      fprintf(stdout, "  -defarad_table  De-Faraday rotate the data for each pulse longitude bin\n");
      fprintf(stdout, "                  separately using the specified ascii table as generated\n");
      fprintf(stdout, "                  by rmsyth.\n");
    }
//START REGION RELEASE
    if(application->switch_deparang)
      fprintf(stdout, "  -deparang       Remove parallactic angle effect from data\n");
//START REGION DEVELOP
    if(application->switch_dededisperse) {
      fprintf(stdout, "  -disp           Undo dedispersion of the data if currently corrected.\n");
    }
    if(application->switch_dedeFaraday) {
      fprintf(stdout, "  -farad          Undo de-Faraday rotation of the data if currently corrected.\n");
    }
    if(application->switch_checkinf) 
      fprintf(stdout, "  -checkinf       Checks data for INF's (or -INF's)\n");
    if(application->switch_removeinf) 
      fprintf(stdout, "  -removeinf      Replaces INF's with zero's in the data\n");
    if(application->switch_checknan) 
      fprintf(stdout, "  -checknan       Checks data for NaN's\n");
    if(application->switch_removenan) 
      fprintf(stdout, "  -removenan      Replaces NaN's with zero's in the data\n");
//START REGION RELEASE
    if(application->switch_coherence)
      fprintf(stdout, "  -coherence      Convert Stokes to coherency parameters\n");
    if(application->switch_debase) 
      fprintf(stdout, "  -debase         Subtract baseline from data (use with -onpulse)\n");
//START REGION DEVELOP
    if(application->switch_debase_slope) 
      fprintf(stdout, "  -debase_slope   As -debase, but takes out a gradient\n");
//START REGION RELEASE
    if(application->switch_fchan) 
      fprintf(stdout, "  -fchan f        Use frequency channel f only\n");
//START REGION DEVELOP
    if(application->switch_fft) 
      fprintf(stdout, "  -fft            Replace data with its fft\n");
    if(application->switch_fftsmooth) {
      fprintf(stdout, "  -fftsmooth f    Smooth each individual pulse up to frequency\n");
      fprintf(stdout, "                  f in Hz (high-pass)\n");
    }
    if(application->switch_fftzap)
      fprintf(stdout, "  -fftzap 'low high'  Zap this range of frequencies [Hz]\n");
    if(application->switch_gate) 
      fprintf(stdout, "  -gate           Throw away data outside range defined by -onpulse\n");
    if(application->switch_invertxy) 
      fprintf(stdout, "  -invertXY       Swap bin-axis and the subint-axis of the data\n");
    if(application->switch_invertfp)
      fprintf(stdout, "  -invertFP       Swap frequency-axis and the subint-axis of the data\n");
    if(application->switch_invertfx)
      fprintf(stdout, "  -invertFX       Swap frequency-axis and the pulse longitude axis of the data\n");
    if(application->switch_invertf)
      fprintf(stdout, "  -invertF        Invert the order of the frequency channels\n");
//START REGION RELEASE
    if(application->switch_norm) {
      fprintf(stdout, "  -norm           Normalize peak value of each subint/channel independently. If\n");
      fprintf(stdout, "                  an onpulse region is set, the peak of that region is used. The\n");
      fprintf(stdout, "                  normalization is determined by the first polarization channel.\n");
    }
    if(application->switch_normglobal) {
      fprintf(stdout, "  -norm_global    As -norm, but use global factor for all subints/channels.\n");
    }
//START REGION DEVELOP
    if(application->switch_normRMS) {
      fprintf(stdout, "  -normRMS        Normalize data using the offpulse rms (use -onpulse). A single\n");
      fprintf(stdout, "                  scale factor is used for the whole dataset, which makes the\n");
      fprintf(stdout, "                  average rms = 1 in the individual channels of the individual\n");
      fprintf(stdout, "                  subints for the first polarization channel.\n");
    }
    if(application->switch_normRMSi) {
      fprintf(stdout, "  -normRMSi       Normalize individual channels/subints using the offpulse rms.\n");
      fprintf(stdout, "                  (use -onpulse) based on the first polarization channel.\n");
    }
//START REGION RELEASE
    if(application->switch_nskip)
      fprintf(stdout, "  -nskip n        Skip n subintegrations in input data (default is 0)\n");
    if(application->switch_nread)
      fprintf(stdout, "  -nread n        Use n subintegrationss in input data (default is all)\n");
    if(application->switch_insertparang)
      fprintf(stdout, "  -parang         Insert parallactic angle effect (-deparang corrects data)\n");
//START REGION RELEASE
    if(application->switch_polselect)
      fprintf(stdout, "  -polselect p    Use polarization channel p (start counting from 0)\n");
    if(application->switch_rebin) {
      fprintf(stdout, "  -rebin n        Rebin to n bins (not necessarily a power of two).\n");
      fprintf(stdout, "                  The mean remains the same after rebinning.\n");
    }
    if(application->switch_conshift) {
      fprintf(stdout, "  -conshift       When applying rotations (-rot or -rotdeg), instead of rotating\n");
      fprintf(stdout, "                  subints independent of each other, rotate the end of one\n");
      fprintf(stdout, "                  subint to the start of another. One subint will be lost.\n");
    }
    if(application->switch_circshift) {
      fprintf(stdout, "  -circshift      idem, but makes the last subint spill over in the first.\n");
      fprintf(stdout, "                  No subints will be lost.\n");
    }
    if(application->switch_changeRefFreq) {
      fprintf(stdout, "  -reffreq f      Change the reference frequency for dedispersion/de-Faraday\n");      
      fprintf(stdout, "                  rotation to specified value f (in MHz, -1=inf frequency).\n");
      fprintf(stdout, "                  Data will be re-dedispersed/de-Farady rotated if required.\n");
    }      
    if(application->switch_rot)
      fprintf(stdout, "  -rot ph         Rotate each individual subint by ph pulse phase\n");
    if(application->switch_rotdeg)
      fprintf(stdout, "  -rotdeg ph      ditto, but in degrees\n");
    if(application->switch_rotateStokes) {
      fprintf(stdout, "  -rotateStokes   \"S1 S2 ph\" Rotate the polarization vectors in Stokes space\n");
      fprintf(stdout, "                  (Q,U,V) by ph degrees. S1 and S2 specify the direction of the\n");
      fprintf(stdout, "                  rotation, such that for instance \"Q U 10\" would be a\n");
      fprintf(stdout, "                  rotation about the Stokes V axis such that a vector in the\n");
      fprintf(stdout, "                  Q direction gets rotated towards the U axis. This option can\n");
      fprintf(stdout, "                  be used multiple times, and they are executed in the order\n");
      fprintf(stdout, "                  as specified on the command line.\n");
    }
//START REGION DEVELOP
    if(application->switch_rotateQU)
      fprintf(stdout, "  -rotateQU ph    Use -rotateStokes instead\n");
    if(application->switch_rotateUV)
      fprintf(stdout, "  -rotateUV ph    Use -rotateStokes instead\n");
    if(application->switch_rotslope) {
      fprintf(stdout, "  -rotslope slope Rotate each individual subint by slope pulse phase\n");
      fprintf(stdout, "                  per subint. So this introduces a slope in the pulse stack.\n");
    }
//START REGION RELEASE
    if(application->switch_scale)
      fprintf(stdout, "  -scale          \"scale offset\". output = scale*(input+offset)\n");
    if(application->switch_shuffle)
      fprintf(stdout, "  -shuffle        Shuffle the subints in a random order\n");
//START REGION DEVELOP
    if(application->switch_smooth) {
      fprintf(stdout, "  -smooth w       Smooth each individual pulse with a top-hat function\n");
      fprintf(stdout, "                  the specified width w in bins. This option acts as a\n");
      fprintf(stdout, "                  low-pass filter.\n");
    }
//START REGION RELEASE
    if(application->switch_stokes)
      fprintf(stdout, "  -stokes         Convert to Stokes parameters\n");
//START REGION DEVELOP
    if(application->switch_subtractAvPA) {
      fprintf(stdout, "  -subtractPA     Subtract the observed time-averaged PA-swing determined from\n");
      fprintf(stdout, "                  input data by rotating Stokes Q and U.\n");
    }
    if(application->switch_subtractRVM) {
      fprintf(stdout, "  -subtractRVM    Subtract the given RVM from the data by rotating Stokes\n");
      fprintf(stdout, "                  Q and U. Provide \"alpha beta pa0 l0\", all in degrees.\n");
    }
    if(application->switch_swapcables) {
      fprintf(stdout, "  -swapcables     Change the polarizations corresponding to a swap of the\n");
      fprintf(stdout, "                  two input polarizations channels in the back-end.\n");
    }
    if(application->switch_testsignal)
      fprintf(stdout, "  -testsignal     Substitute the data read in with a test signal\n");
//START REGION RELEASE
    if(application->switch_tscr) {
      fprintf(stdout, "  -tscr t         Add t successive subints together.\n");
      fprintf(stdout, "                  A negative number duplicates subints.\n");
    }
    if(application->switch_TSCR)
      fprintf(stdout, "  -TSCR           Add all successive subints together\n");
    if(application->switch_tscr_complete) {
      fprintf(stdout, "  -tscr_complete  Use in combination with -tscr. Throw away last subints to\n");
      fprintf(stdout, "                  ensure each subint is the sum of the same number of\n");
      fprintf(stdout, "                  input subints.\n");
    }
    if(application->switch_fscr) {
      fprintf(stdout, "  -fscr f         Add f successive frequency channels together.\n");
      fprintf(stdout, "                  Data will be dedispersed and de-Faraday rotated if required.\n");
      fprintf(stdout, "                  A negative number duplicates frequency channels.\n");
    }
    if(application->switch_FSCR) {
      fprintf(stdout, "  -FSCR           Add all successive frequency channels together.\n");
      fprintf(stdout, "                  Data will be dedispersed and de-Faraday rotated if required.\n");
    }
  //START REGION DEVELOP
    if(application->switch_zerodming) {
      fprintf(stdout, "  -zerodming      Remove a weighted frequency-scrunched non-dedispersed template\n");
      fprintf(stdout, "                  from the data to surpress non-dedispersed RFI.\n");
      fprintf(stdout, "  -zerodming_adv  \"mode poss negs maxw extrab dedisp minw\"\n");
      fprintf(stdout, "                  Similar to -zerodming, except that a significance poss and negs\n");
      fprintf(stdout, "                  can be defined to select time intervals which are corrected. If\n");
      fprintf(stdout, "                  the (freq-averaged) signal is exceeding poss times the rms the\n");
      fprintf(stdout, "                  time sample is corrected, or if the signal is a dip exceeding\n");
      fprintf(stdout, "                  negs times the rms. Set poss and/or negs to a negative nr to\n");
      fprintf(stdout, "                  disable this functionality.\n");
      fprintf(stdout, "                  mode=1: Subtract freq-avrg from each channel\n");
      fprintf(stdout, "                  mode=2: Set each channel to its avrg (zero per def)\n");
      fprintf(stdout, "                  When searching for peaks, if maxw is a positive integer, it\n");
      fprintf(stdout, "                  will do a boxcar search for the most significant peak of\n");
      fprintf(stdout, "                  various widths in the 0-DM profile. Only if the most significant\n");
      fprintf(stdout, "                  peak is more narrow than maxw bins, it will be removed.\n");
      fprintf(stdout, "                  At either side of the found signal extrab extra bins are removed.\n");
      fprintf(stdout, "                  If dedisp is non-zero, the width of a signal will be compared with\n");
      fprintf(stdout, "                  the width after dedispersion. Actual pulsar signals should\n");
      fprintf(stdout, "                  narrow. Only if the width becomes wider by minw\n");
      fprintf(stdout, "                  bins (allowed to be negative), the signal is regarded to be RFI\n");
      fprintf(stdout, "                  signal.\n");
      /*
	fprintf(stdout, "                  \n");
      */
    }
  //START REGION RELEASE
  }
  if(application->switch_itf || application->switch_device || application->switch_size || application->switch_noplotsubset || application->switch_cmaplist || application->switch_cmap 
  //START REGION DEVELOP
     || application->switch_ppgplot
  //START REGION RELEASE
     ) {
    fprintf(stdout, "\nGeneral plotting options:\n");
    if(application->switch_cmap) 
      fprintf(stdout, "  -cmap type    Select color map type\n");
    if(application->switch_cmaplist) 
      fprintf(stdout, "  -cmaplist     List available color map types\n");
  //START REGION RELEASE
    if(application->switch_device) 
      fprintf(stdout, "  -device dev   (or -dev) Specify PGPLOT plotting device dev\n");
    if(application->switch_noplotsubset) 
      fprintf(stdout, "  -noplotsubset Disable plotting subsets of data (this could increase file size)\n");
  //START REGION RELEASE
    if(application->switch_itf) {
      fprintf(stdout, "  -itf n        Set image transfer function for colour map plots\n");
      fprintf(stdout, "                (0=linear (default), 1=logarithmic, 2=square-root)\n");
    }
  //START REGION DEVELOP
    if(application->switch_ppgplot) 
      fprintf(stdout, "  -ppgplot      Output pgplot commands to this pgplot script file\n");
//START REGION RELEASE
    if(application->switch_size) 
      fprintf(stdout, "  -size         \"width height\". Specify resolution of plot device (in pixels).\n");
  }
  if(application->switch_onpulse || application->switch_onpulsef || application->switch_onpulsegr || application->switch_onpulse2 || application->switch_onpulsef2 
//START REGION DEVELOP
|| application->switch_template_onpulse || application->switch_autoonpulse
//START REGION RELEASE
) {
    fprintf(stdout, "\nGeneral data selection options:\n");
//START REGION DEVELOP
   if(application->switch_autoonpulse) 
     fprintf(stdout, "  -ao           automatically try to estimate a good on-pulse regions\n");
//START REGION RELEASE
   if(application->switch_onpulse) 
     fprintf(stdout, "  -onpulse      \"left right\" manually select on-pulse regions (in bins)\n");
   if(application->switch_onpulsef) 
     fprintf(stdout, "  -onpulsef     \"left right\" manually select on-pulse regions (in phase)\n");
   if(application->switch_onpulse2) {
     fprintf(stdout, "  -onpulse2     \"left right\" manually select second type of\n");
     fprintf(stdout, "                on-pulse regions (in bins)\n");
   }
   if(application->switch_onpulsef2) {
     fprintf(stdout, "  -onpulsef2    \"left right\" manually select second type of\n");
     fprintf(stdout, "                on-pulse regions (in phase)\n");
   }
//START REGION DEVELOP
   if(application->switch_template_onpulse) 
      fprintf(stdout, "  -onpulset     Use template to autodetect onpulse region\n");
//START REGION RELEASE
   if(application->switch_onpulsegr) 
      fprintf(stdout, "  -onpulsegr    Graphically select (additional) onpulse regions\n");
  }
  if(application->switch_verbose || application->switch_debug || application->switch_nocounters || application->switch_macro || application->switch_fixseed || application->switch_libversions
//START REGION DEVELOP
     || application->switch_showrevision
//START REGION RELEASE
     ) {
    fprintf(stdout, "\nOther general options:\n");
    if(application->switch_verbose) 
      fprintf(stdout, "  -v            Verbose mode (to get a better idea what is happening)\n");
    if(application->switch_debug) 
      fprintf(stdout, "  -debug        Enable more output (where implemented)\n");
    if(application->switch_fixseed) {
      fprintf(stdout, "  -fixseed      Do not randomy initialise seed of random number generator\n");
      fprintf(stdout, "                thereby making results reproducable.\n");
    }
    if(application->switch_nocounters)
      fprintf(stdout, "  -nocounters   Don't show counters etc (useful when generating log files)\n");
//START REGION DEVELOP
    if(application->switch_showrevision) 
      fprintf(stdout, "  -rev          Show revision information\n");
//START REGION RELEASE
    if(application->switch_libversions) {
      fprintf(stdout, "  -libversions  Show version information about libraries used by psrsalsa\n");      
    }
    if(application->switch_macro) {
      fprintf(stdout, "  -macro        Instead of taking commands from keyboard, read them from\n");
      fprintf(stdout, "                this macro file (put a ^ in front of symbol for the ctrl key)\n");
    }
  }
  fprintf(stdout, "\n");
}

//START REGION DEVELOP
//START REGION RELEASE

int processCommandLine(psrsalsaApplication *application, int argc, char **argv, int *index)
{
  if(strcmp(argv[*index], "-v") == 0 && application->switch_verbose) {
    application->verbose_state.verbose = 1;
    return 1;
  }else if(strcmp(argv[*index], "-debug") == 0 && application->switch_debug) {
    application->verbose_state.debug = 1;
    if(application->switch_verbose)  // Debug implies verbose
      application->verbose_state.verbose = 1;
    return 1;
  }else if(strcmp(argv[*index], "-iformat") == 0 && application->switch_iformat) {
    application->iformat = parsePSRDataFormats(argv[++(*index)]);
    if(application->iformat == 0)
      exit(0);
    return 1;
  }else if(strcmp(argv[*index], "-oformat") == 0 && application->switch_oformat) {
    application->oformat = parsePSRDataFormats(argv[++(*index)]);
    if(application->oformat == 0)
      exit(0);
    return 1;
  }else if(strcmp(argv[*index], "-formatlist") == 0 && application->switch_formatlist) {
    fprintf(stdout, "Supported file formats are:\n");
    printPSRDataFormats(stdout, 2);
    exit(0);
  }else if(strcmp(argv[*index], "-headerlist") == 0 && application->switch_headerlist) {
    printHeaderCommandlineOptions(stdout);
    exit(0);
  }else if(strcmp(argv[*index], "-gentypelist") == 0 && application->switch_headerlist) {
    printHeaderGentypeOptions(stdout);
    exit(0);
  }else if(strcmp(argv[*index], "-header") == 0 && application->switch_header) {
    (*index)++;
    return 1;
  }else if(strcmp(argv[*index], "-headerUFL") == 0 && application->switch_forceUniformFreqLabelling) {
    return 1;
//START REGION RELEASE
  }else if((strcmp(argv[*index], "-libversions") == 0 || strcmp(argv[*index], "-libversion") == 0 || strcmp(argv[*index], "-libraryversions") == 0) && application->switch_libversions) {
    showlibraryversioninformation(stdout);
    exit(0);
  }else if(strcmp(argv[*index], "-macro") == 0 && application->switch_macro) {
    printf("Opening macro '%s'\n", argv[++(*index)]);
    application->macro_ptr = fopen(argv[*index], "r");
    if(application->macro_ptr == NULL) {
      printerror(application->verbose_state.debug, "  Opening macro '%s' failed", argv[*index]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-cmaplist") == 0 && application->switch_cmaplist) {
    printCMAPCommandlineOptions(stdout);
    exit(0);
//START REGION DEVELOP
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-filelist") == 0 && application->switch_filelist) {
    (*index) += 1;
    application->filelist = *index;
    return 1;
  }else if(strcmp(argv[*index], "-nread") == 0 && application->switch_nread) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%ld", &(application->nread), NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-polselect") == 0 && application->switch_polselect) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%d", &(application->polselectnr), NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-nskip") == 0 && application->switch_nskip) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%ld", &(application->nskip), NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-template") == 0 && application->switch_template) {
    application->template_specified = ++(*index);
    if(readVonMisesModel(argv[application->template_specified], &(application->vonMises_components), application->verbose_state) == 0) {
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-templatedata") == 0 && application->switch_templatedata) {
    application->template_data_index = ++(*index);
    
    //    cleanPSRData(&(application->template_file), application->verbose_state);
    closePSRData(&(application->template_file), 0, application->verbose_state);
    if(openPSRData(&(application->template_file), argv[application->template_data_index], 0, 0, 1, 0, application->verbose_state) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot open template file.");
      exit(0);
    }
    datafile_definition clone;
    if(application->template_file.NrSubints > 1) {
      if(!preprocess_addsuccessivepulses(application->template_file, &clone, application->template_file.NrSubints, 1, application->verbose_state))
	return 0;
      swap_orig_clone(&(application->template_file), &clone, application->verbose_state);
    }
    if(application->template_file.NrFreqChan > 1) {
      if(!preprocess_dedisperse(&(application->template_file), 0, 0, 0, application->verbose_state))
	return 0;
      if(!preprocess_addsuccessiveFreqChans(application->template_file, &clone, application->template_file.NrFreqChan, NULL, application->verbose_state))
	return 0;
      swap_orig_clone(&(application->template_file), &clone, application->verbose_state); 
    }
    if(application->template_file.NrPols > 1) {
      printwarning(application->verbose_state.debug, "WARNING processCommandLine: Polarization channel 0 is used as a template.");
      if(preprocess_polselect(application->template_file, &clone, 0, application->verbose_state) == 0)
	return 0;
      swap_orig_clone(&(application->template_file), &clone, application->verbose_state);
    }
    return 1;
  }else if(strcmp(argv[*index], "-onpulse") == 0 && application->switch_onpulse) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%d %d", &(application->onpulse.left_bin[application->onpulse.nrRegions]), &(application->onpulse.right_bin[application->onpulse.nrRegions]), NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->onpulse.bins_defined[application->onpulse.nrRegions] = 1;
    (application->onpulse.nrRegions)++;
    if(application->onpulse.nrRegions == MAX_pulselongitude_regions) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "processCommandLine: To many regions selected.");
      exit(-1);
    }
    return 1;
  }else if(strcmp(argv[*index], "-onpulse2") == 0 && application->switch_onpulse2) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%d %d", &(application->onpulse2.left_bin[application->onpulse2.nrRegions]), &(application->onpulse2.right_bin[application->onpulse2.nrRegions]), NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->onpulse2.bins_defined[application->onpulse2.nrRegions] = 1;
    (application->onpulse2.nrRegions)++;
    if(application->onpulse2.nrRegions == MAX_pulselongitude_regions) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "processCommandLine: To many regions selected.");
      exit(-1);
    }
    return 1;
  }else if(strcmp(argv[*index], "-onpulsef") == 0 && application->switch_onpulsef) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f %f", &(application->onpulse.left_frac[application->onpulse.nrRegions]), &(application->onpulse.right_frac[application->onpulse.nrRegions]), NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->onpulse.frac_defined[application->onpulse.nrRegions] = 1;
    (application->onpulse.nrRegions)++;
    return 1;
  }else if(strcmp(argv[*index], "-onpulsef2") == 0 && application->switch_onpulsef2) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f %f", &(application->onpulse2.left_frac[application->onpulse2.nrRegions]), &(application->onpulse2.right_frac[application->onpulse2.nrRegions]), NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->onpulse2.frac_defined[application->onpulse2.nrRegions] = 1;
    (application->onpulse2.nrRegions)++;
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-ao") == 0 && application->switch_autoonpulse) {
    application->do_autoonpulse = 1;
    return 1;
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-ext") == 0 && application->switch_ext) {
    application->extension = argv[++(*index)];
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-ppgplot") == 0 && application->switch_ppgplot) {
    application->ppgplot_name = argv[++(*index)];
    return 1;
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-output") == 0 && application->switch_output) {
    strcpy(application->outputname, argv[++(*index)]);
    return 1;
  }else if((strcmp(argv[*index], "-device") == 0 || strcmp(argv[*index], "-dev") == 0) && application->switch_device) {
    strcpy(application->pgplotdevice, argv[++(*index)]);
    return 1;
  }else if(strcmp(argv[*index], "-fchan") == 0 && application->switch_fchan) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%d", &application->fchan_select, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-tscr_complete") == 0 && application->switch_tscr_complete) {
    application->tscr_complete = 1;
    return 1;
  }else if(strcmp(argv[*index], "-tscr") == 0 && application->switch_tscr) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%ld", &application->dotscr, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-fscr") == 0 && application->switch_fscr) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%ld", &application->dofscr, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-cmap") == 0 && application->cmap) {
    application->cmap = cmap_parse_commandline(argc, argv, application->verbose_state.debug);
    (*index)++;
    return 1;
  }else if(strcmp(argv[*index], "-rebin") == 0 && application->switch_rebin) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%ld", &application->rebin, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->dorebin = 1;
    return 1;
  }else if(strcmp(argv[*index], "-rot") == 0 && application->switch_rot) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f", &application->shiftPhase, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->shiftPhase_cmdline = application->shiftPhase;
    application->doshiftphase = 1;
    return 1;
  }else if(strcmp(argv[*index], "-rotdeg") == 0 && application->switch_rotdeg) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f", &application->shiftPhase, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->doshiftphase = 1;
    application->shiftPhase /= 360.0;
    application->shiftPhase_cmdline = application->shiftPhase;
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-fftsmooth") == 0 && application->switch_fftsmooth) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f", &application->do_fftsmooth, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-smooth") == 0 && application->switch_smooth) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f", &application->do_smooth, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-fftzap") == 0 && application->switch_fftzap) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f %f", &application->fftzap_low, &application->fftzap_high, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->do_fftzap = 1;
    return 1;
  }else if(strcmp(argv[*index], "-rotslope") == 0 && application->switch_rotslope) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f", &application->rotslope, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->dorotslope = 1;
    return 1;
//START REGION DEVELOP
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-itf") == 0 && application->switch_itf) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%d", &application->itf, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-size") == 0 && application->switch_size) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%d %d", &application->windowwidth, &application->windowheight, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-noplotsubset") == 0 && application->switch_noplotsubset) {
    application->do_noplotsubset = 1;
    return 1;
  }else if(strcmp(argv[*index], "-nocounters") == 0 && application->switch_nocounters) {
    application->verbose_state.nocounters = 1;
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-rev") == 0 && application->switch_showrevision) {
    printf("Revision info:\n");
    SHOWREVISIONINFO();
    printf("\n");
    return 1;
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-debase") == 0 && application->switch_debase) {
    application->dodebase = 1;
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-debase_slope") == 0 && application->switch_debase_slope) {
    application->dodebase_slope = 1;
    return 1;
  }else if(strcmp(argv[*index], "-testsignal") == 0 && application->switch_testsignal) {
    application->do_testsignal = 1;
    return 1;
  }else if(strcmp(argv[*index], "-gate") == 0 && application->switch_gate) {
    application->doGate = 1;
    return 1;
  }else if(strcmp(argv[*index], "-swapcables") == 0 && application->switch_swapcables) {
    application->doswapcables = 1;
    return 1;
  }else if(strcmp(argv[*index], "-zerodming") == 0 && application->switch_zerodming) {
    application->dozerodming = 1;
    return 1;
  }else if(strcmp(argv[*index], "-zerodming_adv") == 0 && application->switch_zerodming_adv) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%d %f %f %d %d %d %d", &application->zerodming_adv_mode, &application->zerodming_adv_possigma, &application->zerodming_adv_negsigma, &application->zerodming_adv_maxwidth, &application->zerodming_adv_extrabins, &application->zerodming_adv_checkdedispersed, &application->zerodming_adv_minwidening, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->dozerodming_adv = 1;
    return 1;
  }else if(strcmp(argv[*index], "-invertXY") == 0 && application->switch_invertxy) {
    application->invertxy = 1;
    return 1;
  }else if(strcmp(argv[*index], "-invertFP") == 0 && application->switch_invertfp) {
    application->invertfp = 1;
    return 1;
  }else if(strcmp(argv[*index], "-invertF") == 0 && application->switch_invertf) {
    application->invertf = 1;
    return 1;
  }else if(strcmp(argv[*index], "-invertFX") == 0 && application->switch_invertfp) {
    application->invertfx = 1;
    return 1;
  }else if(strcmp(argv[*index], "-fft") == 0 && application->switch_fft) {
    application->dofft = 1;
    return 1;
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-norm") == 0 && application->switch_norm) {
    application->do_norm = 1;
    return 1;
  }else if(strcmp(argv[*index], "-norm_global") == 0 && application->switch_normglobal) {
    application->do_normglobal = 1;
    return 1;
  }else if(strcmp(argv[*index], "-TSCR") == 0 && application->switch_TSCR) {
    application->doTSCR = 1;
    return 1;
  }else if(strcmp(argv[*index], "-FSCR") == 0 && application->switch_FSCR) {
    application->doFSCR = 1;
    return 1;
  }else if(strcmp(argv[*index], "-clip") == 0 && application->switch_clip) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f", &application->clipvalue, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->do_clip = 1;
    return 1;
  }else if(strcmp(argv[*index], "-reffreq") == 0 && application->switch_changeRefFreq) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%lf", &application->newRefFreq, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-dedisp") == 0 && application->switch_dedisperse) {
    application->do_dedisperse = 1;
    return 1;
  }else if(strcmp(argv[*index], "-defarad") == 0 && application->switch_deFaraday) {
    application->do_deFaraday = 1;
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-defarad_table") == 0 && application->switch_deFaradayTable) {
    application->do_deFaradayTable = 1;
    application->deFaradayTable_filename = malloc(strlen(argv[(*index)+1])+1);
    if(application->deFaradayTable_filename == NULL) {
      printerror(application->verbose_state.debug, "Cannot allocate memory.");
      exit(0);
    }
    strcpy(application->deFaradayTable_filename, argv[++(*index)]);
    return 1;
  }else if(strcmp(argv[*index], "-disp") == 0 && application->switch_dededisperse) {
    application->do_dedisperse = 2;
    return 1;
  }else if(strcmp(argv[*index], "-farad") == 0 && application->switch_dedeFaraday) {
    application->do_deFaraday = 2;
    return 1;
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-blocksize") == 0 && application->switch_blocksize) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%d", &application->blocksize, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-onpulset") == 0 && application->switch_template_onpulse) {
    application->template_onpulse = 1;
    return 1;
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-align") == 0 && application->switch_align) {
    application->doalign = 1;
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-centert") == 0 && application->switch_centert) {
    application->docentert = 1;
    return 1;
//START REGION DEVELOP
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-conshift") == 0 && application->switch_conshift) {
    application->doconshift = 1;
    return 1;
  }else if(strcmp(argv[*index], "-circshift") == 0 && application->switch_circshift) {
    application->doconshift = 1;
    application->docircshift = 1;
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-autot") == 0 && application->switch_autot) {
    application->doautot = 1;
    return 1;
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-onpulsegr") == 0 && application->switch_onpulsegr) {
    application->doonpulsegr = 1;
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-normRMS") == 0 && application->switch_normRMS) {
    application->donormRMS = 1;
    return 1;
  }else if(strcmp(argv[*index], "-normRMSi") == 0 && application->switch_normRMSi) {
    application->donormRMSi = 1;
    return 1;
  }else if(strcmp(argv[*index], "-checknan") == 0 && application->switch_checknan) {
    application->dochecknan = 1;
    return 1;
  }else if(strcmp(argv[*index], "-removenan") == 0 && application->switch_removenan) {
    application->doremovenan = 1;
    return 1;
  }else if(strcmp(argv[*index], "-checkinf") == 0 && application->switch_checkinf) {
    application->docheckinf = 1;
    return 1;
  }else if(strcmp(argv[*index], "-removeinf") == 0 && application->switch_removeinf) {
    application->doremoveinf = 1;
    return 1;
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-stokes") == 0 && application->switch_stokes) {
    application->dostokes = 1;
    return 1;
  }else if(strcmp(argv[*index], "-coherence") == 0 && application->switch_coherence) {
    application->docoherence = 1;
    return 1;
  }else if(strcmp(argv[*index], "-shuffle") == 0 && application->switch_shuffle) {
    application->doshuffle = 1;
    return 1;
  }else if(strcmp(argv[*index], "-fixseed") == 0 && application->switch_fixseed) {
    application->fixseed = 1;
    return 1;
  }else if(strcmp(argv[*index], "-scale") == 0 && application->switch_scale) {
    application->doscale = 1;
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f %f", &application->scale_scale, &application->scale_offset, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcasecmp(argv[*index], "-rotateStokes") == 0 && application->switch_rotateStokes) {
    if(application->nr_rotateStokes == maxNrRotateStokes) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Maximum number of uses of the '%s' option is exceeded.", argv[(*index)-1]);
      exit(0);
    }
    char s1, s2;
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%c %c %f", &s1, &s2, &(application->rotateStokesAngle[application->nr_rotateStokes]), NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    //    printf("Returned: %c %c %f\n", s1, s2, application->rotateStokesAngle[application->nr_rotateStokes]);
    if(s1 == 'i' || s1 == 'I') {
      application->rotateStokes1[application->nr_rotateStokes] = 0;
    }else if(s1 == 'q' || s1 == 'Q') {
      application->rotateStokes1[application->nr_rotateStokes] = 1;
    }else if(s1 == 'u' || s1 == 'U') {
      application->rotateStokes1[application->nr_rotateStokes] = 2;
    }else if(s1 == 'v' || s1 == 'V') {
      application->rotateStokes1[application->nr_rotateStokes] = 3;
    }else {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option: cannot interpret '%c' as a Stokes parameter.", argv[(*index)-1], s1);
      exit(0);
    }
    if(s2 == 'i' || s2 == 'I') {
      application->rotateStokes2[application->nr_rotateStokes] = 0;
    }else if(s2 == 'q' || s2 == 'Q') {
      application->rotateStokes2[application->nr_rotateStokes] = 1;
    }else if(s2 == 'u' || s2 == 'U') {
      application->rotateStokes2[application->nr_rotateStokes] = 2;
    }else if(s2 == 'v' || s2 == 'V') {
      application->rotateStokes2[application->nr_rotateStokes] = 3;
    }else {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option: cannot interpret '%c' as a Stokes parameter.", argv[(*index)-1], s2);
      exit(0);
    }
    if(application->rotateStokes1[application->nr_rotateStokes] == application->rotateStokes2[application->nr_rotateStokes]) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option: Cannot rotate a Stokes parameter towards the same Stokes parameter.", argv[(*index)-1]);
      exit(0);
    }
    (application->nr_rotateStokes)++;
    return 1;
//START REGION DEVELOP
  }else if(strcasecmp(argv[*index], "-subtractRVM") == 0 && application->switch_subtractRVM) {
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f %f %f %f", &(application->subtractRVM_alpha), &(application->subtractRVM_beta), &(application->subtractRVM_pa0), &(application->subtractRVM_l0), NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    application->dosubtractRVM = 1;
    return 1;
  }else if(strcasecmp(argv[*index], "-subtractPA") == 0 && application->switch_subtractAvPA) {
    application->dosubtractAvPA = 1;
    return 1;
  }else if(strcmp(argv[*index], "-rotateQU") == 0 && application->switch_rotateQU) {
    application->dorotateQU = 1;
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f", &application->rotateQUangle, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
  }else if(strcmp(argv[*index], "-rotateUV") == 0 && application->switch_rotateUV) {
    application->dorotateUV = 1;
    if(parse_command_string(application->verbose_state, argc, argv, ++(*index), 0, -1, "%f", &application->rotateUVangle, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "Cannot parse '%s' option.", argv[(*index)-1]);
      exit(0);
    }
    return 1;
//START REGION RELEASE
  }else if(application->switch_deparang && strcmp(argv[*index], "-deparang") == 0) {
    application->do_parang_corr = 1;
    return 1;
  }else if(application->switch_insertparang && strcmp(argv[*index], "-parang") == 0) {
    application->do_parang_corr = 2;
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-useweightedfreq") == 0 && application->switch_useweightedfreq) {
    application->useweightedfreq = 1;
    return 1;
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-noweights") == 0 && application->switch_noweights) {
    application->noweights = 1;
    psrfits_set_noweights(1);
    return 1;
  }else if(strcmp(argv[*index], "-uniformweights") == 0 && application->switch_uniformweights) {
    application->uniformweights = 1;
    psrfits_set_noweights(2);
    return 1;
  }else if(strcmp(argv[*index], "-useweights") == 0 && application->switch_useweights) {
    application->useweights = 1;
    psrfits_set_noweights(3);
    return 1;
//START REGION DEVELOP
  }else if(strcmp(argv[*index], "-absweights") == 0 && application->switch_absweights) {
    application->absweights = 1;
    psrfits_set_absweights(1);
    return 1;
//START REGION RELEASE
  }else if(strcmp(argv[*index], "-history_cmd_only") == 0 && application->switch_history_cmd_only) {
    application->history_cmd_only = 1;
    return 1;
  }else {
    return 0;
  }
  fflush(stdout);
  printwarning(application->verbose_state.debug, "WARNING processCommandLine: This line shouldn't be executed!");
  return 0;
}

//START REGION DEVELOP
//START REGION RELEASE
/*
It is assume outputname has a length MaxFilenameLength.

The output name is either from the outputname set with -output or the
extension set with -ext

Returns 1 on success, 0 on error.
 */
int getOutputName(psrsalsaApplication *application, char *filename, char *outputname, verbose_definition verbose)
{
  if(application->extension != NULL) {
    if(change_filename_extension(filename, outputname, application->extension, MaxFilenameLength, verbose) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "getOutputName: Changing filename failed");
      return 0;
    }
  }else {
    strcpy(outputname, application->outputname);
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/*
Applies any preprocess options specified on the command line. Returns 0 on error, 1 on success.
 */
int preprocessApplication(psrsalsaApplication *application, datafile_definition *psrdata)
{
  datafile_definition clone;
  int device, original_gentype, original_poltype, original_isDeDisp, original_isDeFarad, original_isDePar, original_isDebase;
  double original_freq_ref;
  float x;
  verbose_definition verbose1, verbose2;

  //de-dispersed=yes de-Faraday=no de-par. angle=no

  long i;
//START REGION DEVELOP
  float y, x1, x2, ymax, shift, baseline, rms;

//START REGION RELEASE

  original_gentype   = psrdata->gentype;
  original_poltype   = psrdata->poltype;
  original_isDeDisp  = psrdata->isDeDisp;
  original_isDeFarad = psrdata->isDeFarad;
  original_isDePar   = psrdata->isDePar;
  original_isDebase  = psrdata->isDebase;
  original_freq_ref  = psrdata->freq_ref;

  if(application->verbose_state.verbose) {
    printf("\nApplying preprocess options\n");
  }

  /* Verbose level determines nr spaces before output */
  copyVerboseState(application->verbose_state, &verbose1);
  copyVerboseState(application->verbose_state, &verbose2);
  verbose1.indent = application->verbose_state.indent + 2;
  verbose2.indent = application->verbose_state.indent + 4;

//START REGION DEVELOP

//START REGION DEVELOP

  if(application->do_testsignal) {
    if(preprocess_testsignal(*psrdata, verbose1) != 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Injecting test signal failed.");
      return 0;
    }
  }
  if(application->dochecknan) {
    if(preprocess_checknan(*psrdata, 1, verbose1) == 0)
      return 0;
    else {
      fflush(stdout);
      printerror(application->verbose_state.debug, "The data contains NaN's. You could use -removenan to remove these from the data, although you might want to reprocess the data.");
      exit(-1);
    }
  }
  if(application->doremovenan) {
    preprocess_removenan(*psrdata, verbose1);
  }
  if(application->docheckinf) {
    if(preprocess_checkinf(*psrdata, 1, verbose1) == 0)
      return 0;
    else {
      fflush(stdout);
      printerror(application->verbose_state.debug, "The data contains INF's. You could use -removeinf to remove these from the data, although you might want to reprocess the data.");
      exit(-1);
    }
  }
  if(application->doremoveinf) {
    preprocess_removeinf(*psrdata, verbose1);
  }
//START REGION DEVELOP
//START REGION RELEASE
  if(application->nskip != 0 || application->nread > 0) {
    if(application->nread <= 0)
      application->nread = psrdata->NrSubints-application->nskip;
    if(preprocess_pulsesselect(*psrdata, &clone, application->nskip, application->nread, verbose1) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Error selecting pulses.");
      return 0;
    }
    swap_orig_clone(psrdata, &clone, application->verbose_state); 
  }
//START REGION DEVELOP
  if(application->doswapcables) {
    if(preprocess_swapcables(psrdata, verbose1) == 0)
      return 0;
  }
//START REGION RELEASE
  if(application->dostokes) {
    if(preprocess_stokes(psrdata, verbose1) == 0)
      return 0;
  }
  if(application->docoherence) {
    if(preprocess_coherency(psrdata, verbose1) == 0)
      return 0;
  }
  if(application->nr_rotateStokes > 0) {
    for(i = 0; i < application->nr_rotateStokes; i++) {
      // Do inplace rotation
      if(preprocess_rotateStokes(psrdata, &clone, 1, -1, application->rotateStokesAngle[i], NULL, application->rotateStokes1[i], application->rotateStokes2[i], verbose1) == 0)
	return 0;
    }
  }
//START REGION DEVELOP
  if(application->dozerodming_adv) {
    if(preprocess_zero_dming(psrdata, application->zerodming_adv_possigma, application->zerodming_adv_negsigma, application->zerodming_adv_mode, application->zerodming_adv_maxwidth, application->zerodming_adv_extrabins, application->zerodming_adv_checkdedispersed, application->zerodming_adv_minwidening, verbose1) == 0)
      return 0;
  }else if(application->dozerodming) {  // If -zerodming and -zerodming_adv are both specified, -zerodming_adv was probably implied
    if(preprocess_zero_dming(psrdata, -1, -1, -1, -1, 0, 0, -1, verbose1) == 0)
      return 0;
  }
  if(application->dosubtractRVM) {
    if(preprocess_subtractRVM(psrdata, &clone, 1, application->subtractRVM_alpha, application->subtractRVM_beta, application->subtractRVM_pa0, application->subtractRVM_l0, verbose1) == 0)
      return 0;
  }
  if(application->dosubtractAvPA) {
    if(preprocess_subtractAveragePA(psrdata, &clone, 1, verbose1) == 0)
      return 0;
  }
  if(application->dorotateQU) {
    if(preprocess_rotateStokes(psrdata, &clone, 1, -1, application->rotateQUangle, NULL, 1, 2, verbose1) == 0)
      return 0;
    //    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(application->dorotateUV) {
    if(preprocess_rotateStokes(psrdata, &clone, 1, -1, application->rotateUVangle, NULL, 2, 3, verbose1) == 0)
      return 0;
    //    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
//START REGION RELEASE
  if(application->do_parang_corr > 0) {
    if(application->do_parang_corr == 2) { // Insert parallactic angle
      if(preprocess_corrParAng(psrdata, NULL, 1, verbose1) == 0)
	return 0;
    }else {  // Remove parallactic angle
      if(preprocess_corrParAng(psrdata, NULL, 0, verbose1) == 0)
	return 0;
    }
  }
  if(application->blocksize > 0) {
    if(preprocess_blocksize(*psrdata, &clone, application->blocksize, verbose1) == 0)
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(application->fchan_select != -1) {
    if(preprocess_channelselect(*psrdata, &clone, application->fchan_select, verbose1) == 0)
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
//START REGION DEVELOP
  if(application->invertxy) {
    if(preprocess_invertXY(*psrdata, &clone, verbose1) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Error inverting X and Y axis.");
      return 0;
    }
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
//START REGION DEVELOP
//START REGION RELEASE
  if(application->polselectnr >= 0) {
    if(preprocess_polselect(*psrdata, &clone, application->polselectnr, verbose1) == 0)
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
//START REGION DEVELOP
  if(application->invertfp) {
    if(preprocess_invertFP(psrdata, verbose1) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Error inverting frequency and subint axis.");
      return 0;
    }
  }
  if(application->invertf) {
    if(preprocess_invertF(*psrdata, &clone, verbose1) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Error inverting order of frequency channels.");
      return 0;
    }
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(application->invertfx) {
    if(preprocess_invertFX(*psrdata, &clone, verbose1) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Error inverting frequency/pulse longidude axes.");
      return 0;
    }
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
//START REGION RELEASE
  if(application->newRefFreq > -2) {  // If a new reference frequency is specified
    if(preprocess_changeRefFreq(psrdata, application->newRefFreq, verbose1) == 0) {
      printerror(application->verbose_state.debug, "preprocessApplication: Error changing reference frequency.");
      return 0;
    }    
  }
  if(application->doFSCR) {
    application->dofscr = psrdata->NrFreqChan;
    if(psrdata->NrFreqChan <= 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: -FSCR expect a positive number of channels. Maybe header parameters not loaded by main program?");
      return 0;
    }
  }
  if(application->do_dedisperse || application->dofscr) {
    if(application->do_dedisperse == 2) { // Undo dedispersion
      if(!preprocess_dedisperse(psrdata, 1, 0, 0, verbose1))
	return 0;
    }else {
      if(!preprocess_dedisperse(psrdata, 0, 0, 0, verbose1))
	return 0;
    }
  }
//START REGION DEVELOP
  if(application->do_deFaradayTable) {
    long i, j, nrdatapoints;
    double *rmtable, *rmtable2;
    int *bintable, showwarning;
    if(read_ascii_column_double(application->deFaradayTable_filename, 0, '#', 5, 0, &nrdatapoints, 3, 1.0, 0, &rmtable, NULL, NULL, NULL, application->verbose_state, 0) == 0) {
      printerror(application->verbose_state.debug, "  Reading from file '%s' failed", application->do_deFaradayTable);
      return 0;
    }
    if(read_ascii_column_int(application->deFaradayTable_filename, 0, '#', 5, 0, &nrdatapoints, 1, &bintable, NULL, NULL, NULL, application->verbose_state, 0) == 0) {
      printerror(application->verbose_state.debug, "  Reading from file '%s' failed", application->do_deFaradayTable);
      return 0;
    }
    rmtable2 = malloc(psrdata->NrBins*sizeof(double));
    if(rmtable2 == NULL) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Memory allocation error");
      return 0;
    }
    showwarning = 1;
    for(i = 0; i < psrdata->NrBins; i++) {
      for(j = 0; j < nrdatapoints; j++) {
	if(bintable[j] == i)
	  break;
      }
      if(j == nrdatapoints) {
	rmtable2[i] = psrdata->rm;
	if(showwarning) {
	  printwarning(application->verbose_state.debug, "preprocessApplication: RM table does not provide RM values for all pulse longitude bins. For the missing bins the RM of the header will be applied.");	  
	  showwarning = 0;
	}
      }else {
	rmtable2[i] = rmtable[j];
      }
    }
    //    for(i = 0; i < psrdata->NrBins; i++) {
    //      printf("%ld %lf\n", i, rmtable2[i]);
    //    }
    if(preprocess_deFaraday(psrdata, 0, 0, 0, rmtable2, verbose1) == 0)
      return 0;
    free(rmtable); 
    free(rmtable2);
    free(bintable);
  }
//START REGION RELEASE
  if(application->do_deFaraday || application->dofscr) {
    int skip;
    skip = 0;
    if(psrdata->NrPols != 4) {
      if(application->do_deFaraday == 0) {
	skip = 1;
      }
    }
    if(skip == 0) {
      if(application->do_deFaraday == 2 && application->dofscr == 0) { // Undo de-Faraday rotation
	if(!preprocess_deFaraday(psrdata, 1, 0, 0, NULL, verbose1))
	  return 0;
      }else {
	if(application->do_deFaraday == 2) {
	  fflush(stdout);
	  printerror(application->verbose_state.debug, "preprocessApplication: You cannot frequency scrunch and use the -farad option at the same time.");
	  return 0;
 	}
	if(!preprocess_deFaraday(psrdata, 0, 0, 0, NULL, verbose1))
	  return 0;
      }
    }
  }
  if(application->dofscr) {
    if(!preprocess_addsuccessiveFreqChans(*psrdata, &clone, application->dofscr, application->fzapMask, verbose1))
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state); 
  }
//START REGION DEVELOP
//START REGION RELEASE
  if(application->doTSCR) {
    application->dotscr = psrdata->NrSubints;
    if(psrdata->NrSubints <= 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: -TSCR expect a positive number of subints. Maybe header parameters not loaded by main program?");
      return 0;
    }
  }
  if(application->dotscr) {
    if(!preprocess_addsuccessivepulses(*psrdata, &clone, application->dotscr, application->tscr_complete, verbose1))
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state); 
  }
//START REGION DEVELOP
  if(application->doautot) {
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("Making template\n");
    }
    if(application->template_specified != 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Can't use -autot and -template option together.");
      return 0;
    }
    /* Make a time-scrunched clone. */
    if(!preprocess_addsuccessivepulses(*psrdata, &clone, psrdata->NrSubints, application->tscr_complete, verbose2)) {
      return 0;
    }
    if(fitvonmises(clone, &(application->vonMises_components), 10, 5, 1, 1, 2, 0, &baseline, NULL, verbose2, 0, NULL) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Generating template failed.");
      return 0;
    }
    application->template_specified = 1;
    /* Destroy the clone */
    closePSRData(&clone, 0, verbose2);
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("  done\n");
    }
  }
  if(application->docentert) {
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("Centring template\n");
    }
    if(application->template_specified == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Can only use -centert option together with the -template option.");
      return 0;
    }
    ymax = 0;
    for(i = 0; i < 1000*psrdata->NrBins; i++) {
      x = i/(float)(1000*psrdata->NrBins);
      y = calcVonMisesFunction(&(application->vonMises_components), x, 0);
      if(y > ymax) {
	ymax = y;
	shift = x;
      }
    }
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("  Applying %f phase shift to template\n", shift);
    }
    for(i = 0; i < application->vonMises_components.nrcomponents; i++)
      application->vonMises_components.centre[i] -= shift;
  }
//START REGION RELEASE
  if(application->doalign) {
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("Aligning data using template\n");
    }
    if(application->template_specified == 0 && application->template_data_index == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Can only use -align option together with a specified template on the command line.");
      return 0;
    }
    /*
    if(psrdata->NrFreqChan > 1) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Can only use -align option on frequency scrunched data.");
      return 0;
    }
    */
    /* Make a time-scrunched clone. */
    if(!preprocess_addsuccessivepulses(*psrdata, &clone, psrdata->NrSubints, application->tscr_complete, verbose2))
      return 0;
    if(!preprocess_dedisperse(&clone, 0, 0, 0, verbose2))
      return 0;
    if(psrdata->NrPols == 4) {
      if(!preprocess_deFaraday(&clone, 0, 0, 0, NULL, verbose2))
	return 0;
    }
    datafile_definition clone2;
    if(!preprocess_addsuccessiveFreqChans(clone, &clone2, clone.NrFreqChan, NULL, verbose2))
      return 0;
    if(application->template_specified) {  // If using a mathematical template
      /* Correlate profile with template to find offset */
      x = correlateVonMisesFunction(&(application->vonMises_components), clone2.NrBins, clone2.data, verbose2);
      /* Destroy the clone */
    }else {
      if(clone2.NrBins != application->template_file.NrBins) {
	fflush(stdout);
	printerror(application->verbose_state.debug, "preprocessApplication: The template and the data file have a different amount of bins (%ld != %ld).", clone2.NrBins, application->template_file.NrBins);
	return 0;
      }
      int lag;
      float correl_max;
      if(find_peak_correlation(clone2.data, application->template_file.data, clone2.NrBins, 0, 0, 1, &lag, &correl_max, verbose2) == 0) {
	return 0;
      }
      x = lag/(double)clone2.NrBins;
    }
    closePSRData(&clone2, 0, verbose2);
    closePSRData(&clone, 0, verbose2);
    /* Rotate the data */
    if(application->doshiftphase) {
      application->shiftPhase -= x;
    }else {
      application->doshiftphase = 1;
      application->shiftPhase = -x;
    }
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("  done       \n");
    }
  }
//START REGION DEVELOP
//START REGION RELEASE
  if(application->doshiftphase) {
    if(application->doconshift) {
      //      printf("XXXXX %f %ld\n",application->shiftPhase, psrdata->NrBins);
      i = application->shiftPhase*psrdata->NrBins;
     /* -rot in combination with -doalign can cause invalid shifts */
     if(i >= psrdata->NrBins)
       i -= psrdata->NrBins;
     if(i < 0)
       i += psrdata->NrBins;
     if(verbose1.verbose) {
       int j;
       for(j = 0; j < verbose1.indent; j++)      
	 printf(" ");
       printf("Rotating data by %ld bins\n", i);
     }
     if(continuous_shift(*psrdata, &clone, i, application->docircshift, "preprocessApplication", MEMORY_format, 0, NULL, verbose2, 0) == 0)
	return 0;
      swap_orig_clone(psrdata, &clone, application->verbose_state); 
    }else {
      if(preprocess_fftshift(*psrdata, application->shiftPhase, 0, 0, verbose2) == 0)
	return 0;
    }
  }
//START REGION DEVELOP
  if(application->dorotslope) {
    if(preprocess_fftshift(*psrdata, 0, 1, application->rotslope, verbose2) == 0)
      return 0;    
  }
  if(application->do_fftsmooth > 0) {
    x = 1.0/(get_tsamp(*psrdata, 0, verbose2)*(psrdata->NrBins));
    x = (application->do_fftsmooth)/x;
    /*    printf("XXXXXX %f\n", x); */
    /*0.5 * (application->do_fftsmooth) * (float)(psrdata->NrBins)*/
    if(preprocess_fftSmooth(*psrdata, x, 1, verbose1) == 0)
      return 0;
  }
  if(application->do_fftzap > 0) {
    x = 1.0/(get_tsamp(*psrdata, 0, verbose2)*(psrdata->NrBins));
    x1 = (application->fftzap_low)/x;
    x2 = (application->fftzap_high)/x;
    if(preprocess_fftZap(*psrdata, x1, x2, verbose1) == 0)
      return 0;
  }
  if(application->do_smooth > 0) {
    if(preprocess_smooth(*psrdata, application->do_smooth, verbose1) == 0)
      return 0;
  }
//START REGION RELEASE
  if(application->dorebin) {
    if(!preprocess_rebin(*psrdata, &clone, application->rebin, verbose1))
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state); 
  }
//START REGION DEVELOP
  if(application->do_autoonpulse) {
    if(preprocess_autoOnpulse(*psrdata, -1, -1, -1, -1, -1, -1, &(application->onpulse), verbose1) == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: auto onpulse selection failed.");
      return 0;      
    }
    regionShowNextTimeUse(application->onpulse, "-onpulse", "-onpulsef", stdout);
  }
  if(application->template_onpulse) {
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("Finding onpulse region via template\n");
    }
    if(application->template_specified == 0) {
      fflush(stdout);
      printerror(application->verbose_state.debug, "preprocessApplication: Can only use -onpulset option together with the -template option.");
      return 0;
    }
    /* Make a time-scrunched clone. */
    if(!preprocess_addsuccessivepulses(*psrdata, &clone, psrdata->NrSubints, application->tscr_complete, verbose2))
      return 0;
    /* Correlate profile with template to find offset */
    x = correlateVonMisesFunction(&(application->vonMises_components), clone.NrBins, clone.data, verbose2);
    /* Replace data in clone with the aligned template */
    calcVonMisesProfile(&(application->vonMises_components), clone.NrBins, clone.data, x, 1);
    /* Measure on-pulse region using the aligned template */
    find_boundaries(clone.data, clone.NrBins, 0.05, &(application->onpulse));
    /* Destroy the clone */
    closePSRData(&clone, 0, verbose2);
    printf("    ");
    regionShowNextTimeUse(application->onpulse, "-onpulse", "-onpulsef", stdout);
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("  done\n");
    }
  }
//START REGION RELEASE
  if(application->doonpulsegr) {
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("Select more onpulse regions\n");
    }
    /* Make a time-scrunched clone. */
    if(!preprocess_addsuccessivepulses(*psrdata, &clone, psrdata->NrSubints, application->tscr_complete, verbose2))
      return 0;
    ppgqid(&device);
    printf("Device = %d\n", device);
    pgplot_options_definition *pgplot_options;
    pgplot_options = (pgplot_options_definition *)malloc(sizeof(pgplot_options_definition));
    if(pgplot_options == NULL) {
      printerror(application->verbose_state.debug, "ERROR preprocessApplication: Memory allocation error");
      return 0;
    }
    pgplot_clear_options(pgplot_options);
    strcpy(pgplot_options->box.xlabel, "Bin");
    strcpy(pgplot_options->box.ylabel, "Intensity");
    strcpy(pgplot_options->box.title, "Select on-pulse region of ");
    strcat(pgplot_options->box.title, psrdata->psrname);
    selectRegions(clone.data, clone.NrBins, pgplot_options, 0, 0, 0, &(application->onpulse), verbose2); 
    if(device)
      ppgslct(device);
    /* Destroy the clone */
    closePSRData(&clone, 0, verbose2);
    regionShowNextTimeUse(application->onpulse, "-onpulse", "-onpulsef", stdout);
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("  done\n");
    }
    free(pgplot_options);
  }
//START REGION DEVELOP
//START REGION RELEASE
  if(application->dodebase || application->dodebase_slope) {
    if(application->dodebase == 2) {
//START REGION DEVELOP
      if(application->dodebase_slope) {
	if(!preprocess_debase(psrdata, &(application->onpulse2), NULL, 1, verbose1))
	  return 0;
      }else {
//START REGION RELEASE
	if(!preprocess_debase(psrdata, &(application->onpulse2), NULL, 0, verbose1))
	  return 0;
//START REGION DEVELOP
      }
//START REGION RELEASE
    }else {
//START REGION DEVELOP
      if(application->dodebase_slope) {
	if(!preprocess_debase(psrdata, &(application->onpulse), NULL, 1, verbose1))
	  return 0;
      }else {
//START REGION RELEASE
	if(!preprocess_debase(psrdata, &(application->onpulse), NULL, 0, verbose1))
	  return 0;
//START REGION DEVELOP
      }
//START REGION RELEASE
    }
  }
//START REGION DEVELOP
  if(application->doGate) {
    if(preprocess_gate(*psrdata, &clone, application->onpulse, verbose1) == 0)
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
  if(application->dofft) {
    /* Last sample is ignored if nrbins is odd */
    i = (psrdata->NrBins)/2;
    i *= 2;
    /* Change sampling time to sampling frequency in Hz */
    psrdata->fixedtsamp = 1.0/((float)i*get_tsamp(*psrdata, 0, verbose2)); 
    if(preprocess_fft(psrdata, verbose1) == 0)
      return 0;
    psrdata->tsampMode = TSAMPMODE_FIXEDTSAMP;
    if(psrdata->tsamp_list != NULL) {
      free(psrdata->tsamp_list);
      psrdata->tsamp_list = NULL;
    }
  }
//START REGION RELEASE
  if(application->do_norm) {
    if(preprocess_norm(*psrdata, application->normvalue, &(application->onpulse), 0, verbose1) == 0)
      return 0;
  }
  if(application->do_normglobal) {
    if(preprocess_norm(*psrdata, application->normvalue, &(application->onpulse), 1, verbose1) == 0)
      return 0;
  }
  if(application->do_clip) {
    if(preprocess_clip(*psrdata, application->clipvalue, verbose1) == 0)
      return 0;
  }
//START REGION DEVELOP
  if(application->donormRMS) {
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("Normalizing data using the off-pulse RMS using a single scale factor\n");
    }
    rms = determine_singlepulse_rms(*psrdata, 0, &(application->onpulse), 0, verbose2);
    if(preprocess_scale(*psrdata, 1.0/rms, 0, verbose2) == 0)
      return 0;
    if(verbose1.verbose) {
      for(i = 0; i < verbose1.indent; i++)      
	printf(" ");
      printf("  done\n");
    }
  }
  if(application->donormRMSi) {
    if(preprocess_normRMS(*psrdata,1.0,  &(application->onpulse), verbose2) == 0)
      return 0;
  }
//START REGION RELEASE
  if(application->doscale) {
    if(preprocess_scale(*psrdata, application->scale_scale, application->scale_offset, verbose1) == 0)
      return 0;
  }
  if(application->doshuffle) {
    if(preprocess_shuffle(*psrdata, &clone, application->fixseed, verbose1) == 0)
      return 0;
    swap_orig_clone(psrdata, &clone, application->verbose_state);
  }
// Processing done

  if(original_gentype != psrdata->gentype && application->verbose_state.debug) {
    for(i = 0; i < verbose1.indent; i++)      
      printf(" ");
    printf("Gentype of data changed from %s into %s\n", returnGenType_str(original_gentype), returnGenType_str(psrdata->gentype));
  }

  // Took out warning if only debase state has changed, as this should be obvious from the command line used.
  // || original_isDebase != psrdata->isDebase
  if(original_poltype != psrdata->poltype || original_isDeDisp != psrdata->isDeDisp || original_isDeFarad != psrdata->isDeFarad || original_isDePar != psrdata->isDePar || original_freq_ref != psrdata->freq_ref) {
    fflush(stdout);
    char *txt, *txt2;
    txt = malloc(10000);
    txt2 = malloc(10000);
    if(txt == NULL || txt2 == NULL) {
      printerror(application->verbose_state.debug, "ERROR preprocessApplication: Memory allocation error");
      return 0;
    }
    sprintf(txt, "WARNING: Note that after preprocessing the data are now:");
    if(original_poltype != psrdata->poltype) {
      if(psrdata->poltype == POLTYPE_STOKES) {
	strcat(txt, " Stokes parameters");
      }else if(psrdata->poltype == POLTYPE_COHERENCY) {
	strcat(txt, " Coherency parameters");
      }else if(psrdata->poltype == POLTYPE_ILVPAdPA) {
	strcat(txt, " I, L, V, PA and its error");
      }else if(psrdata->poltype == POLTYPE_ILVPAdPATEldEl) {
	strcat(txt, " I, L, V, PA and its error, total polarization, ellipticity and its error");
      }else if(psrdata->poltype == POLTYPE_PAdPA) {
	strcat(txt, " PA and its error");
      }else {
	strcat(txt, " unknown polarisation state");
      }
    }
    if(original_isDeDisp != psrdata->isDeDisp) {
      if(psrdata->isDeDisp) {
	strcat(txt, " dedispersed");
      }else {
	strcat(txt, " not dedispersed");
      }
    }
    if(original_isDeFarad != psrdata->isDeFarad) {
      if(psrdata->isDeFarad) {
	strcat(txt, " de-Faraday rotated");
      }else {
	strcat(txt, " not de-Faraday rotated");
      }
    }
    if(original_isDePar != psrdata->isDePar) {
      if(psrdata->isDePar) {
	strcat(txt, " de-parallactic angle rotated");
      }else {
	strcat(txt, " not de-parallactic angle rotated");
      }
    }
    if(original_isDebase != psrdata->isDebase) {
      if(psrdata->isDebase) {
	strcat(txt, " with baseline removed");
      }else {
	strcat(txt, " without removed baseline");
      }
    }

    if(original_freq_ref != psrdata->freq_ref) {
      if(psrdata->freq_ref > 0)
	sprintf(txt2, " with reference freq=%lf MHz", psrdata->freq_ref);
      else if((psrdata->freq_ref > -1.01 && psrdata->freq_ref < -0.99) || (psrdata->freq_ref > 0.99e10 && psrdata->freq_ref < 1.01e10))
	sprintf(txt2, " with reference freq=infinity");
      else
	sprintf(txt2, " with unknown reference freq");
      strcat(txt, txt2);
    }

    printwarning(application->verbose_state.debug, "%s", txt);
    free(txt);
    free(txt2);
  }

  if(verbose1.verbose) {
    printf("Preprocessing done\n\n");
  }

  if(application->verbose_state.debug) {
    printHeaderPSRData(*psrdata, 1, application->verbose_state);
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Call when you want to add string argi from argv to the list of
filenames.  returns 1 on success, 0 on error */
int applicationAddFilename(int argi, verbose_definition verbose)
{
  if(internal_application_cmdline_Nrfilenames < MaxNrApplicationFilenames) {
    internal_application_cmdline_FilenameList[internal_application_cmdline_Nrfilenames++] = argi;
  }else {
    fflush(stdout);
    printerror(verbose.debug, "applicationAddFilename: Too many file names on command-line (hard-coded limit is currently %d). Consider using -filelist.", MaxNrApplicationFilenames);
    return 0;
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* returns 1 if the filenames appear consecutive, or else returns
   0. If the list is empty returns 1. If argv != NULL, an error
   message is produced stating which options failed. */
int applicationFilenameList_checkConsecutive(char **argv, verbose_definition verbose)
{
  int i, j, indx_last;
  char txt[1000], txt2[1000];
  if(internal_application_cmdline_Nrfilenames < 2)
    return 1;
  indx_last = internal_application_cmdline_FilenameList[0];
  for(i = 1; i < internal_application_cmdline_Nrfilenames; i++) {
    if(internal_application_cmdline_FilenameList[i] != indx_last+1) {
      if(argv != NULL) {
	fflush(stdout);
	sprintf(txt, "Cannot parse command line. The following list of arguments (those not recognized as command line options) was interpreted as a list of file names:\n\n");
	for(j = 0; j < internal_application_cmdline_Nrfilenames; j++) {
	  fflush(stdout);
	  sprintf(txt2, "'%s' ", argv[internal_application_cmdline_FilenameList[j]]);
	  strcat(txt, txt2);
	}
	fflush(stdout);
	strcat(txt, "\n");
	strcat(txt, "\nThese should however be consecutive, so something appears to be wrong.");
	printerror(verbose.debug, "%s", txt);
      }
      return 0;
    }
    indx_last++;
  }
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE


int internal_open_filelist(psrsalsaApplication *application, char **argv, verbose_definition verbose)
{
  int nrColumns;
  internal_application_filelist_fptr = fopen(argv[application->filelist], "r");
  if(internal_application_filelist_fptr == NULL) {
    fflush(stdout);
    printerror(application->verbose_state.debug, "internal_open_filelist: Cannot open %s", argv[application->filelist]);
    return 0;
  }
  internal_application_filelist_fileopen = 1;

  //Determine the nr of lines and ensure there is at least one column (only first column is used)
  nrColumns = -1;
  if(ascii_file_stats(internal_application_filelist_fptr, '#', &internal_application_filelist_Nrfilenames, 2048, 0, &nrColumns, verbose) == 0) {
    fflush(stdout);
    printerror(application->verbose_state.debug, "internal_open_filelist: Determining the number of lines in %s failed", argv[application->filelist]);
    return 0;
  }
  rewind(internal_application_filelist_fptr);

  if(application->verbose_state.verbose) {
    printf("Opened input file name file '%s' with %ld file names\n", argv[application->filelist], internal_application_filelist_Nrfilenames);
  }

  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

/* Returns the number of filenames in the commandline + those
   specified in the file specified with -filelist. */
int numberInApplicationFilenameList(psrsalsaApplication *application, char **argv, verbose_definition verbose)
{
  long total;
  total = internal_application_cmdline_Nrfilenames;  // Those specified at the end of the filelist
  // Check if a list with input files is specified. 
  if(application->filelist) {
    // If so, open the file if not open yet
    if(internal_application_filelist_fileopen == 0) {
      if(internal_open_filelist(application, argv, verbose) == 0) {
	fflush(stdout);
	printerror(application->verbose_state.debug, "numberInApplicationFilenameList: Opening list with file names failed");
	return 0;
      }
    }
    total += internal_application_filelist_Nrfilenames;
  }
  //  printf("XXXXX There are %ld files specified\n", total);
  return total;
}

//START REGION DEVELOP
//START REGION RELEASE

/*
  Resets current position in the input filename list (either on the
  command-line and/or in the filelist specified with the -filelist
  option).
 */
void rewindFilenameList(psrsalsaApplication *application) {
  internal_application_cmdline_CurFilename = 0;
  internal_application_filelist_CurFilename = 0;
  if(application->filelist)
    rewind(internal_application_filelist_fptr);
}

//START REGION DEVELOP
//START REGION RELEASE

/* Returns NULL if the last file has been processed. */
char *getNextFilenameFromList(psrsalsaApplication *application, char **argv, verbose_definition verbose)
{
  int i;
  if(application->doautot != 0 && internal_application_cmdline_CurFilename != 0) {
    /* If auto template is specified, unset template so it can be
       re-generated for the next file. */
    application->template_specified = 0;
  }
  /* Reset the phase shift to that specified at the command line */
  application->shiftPhase = application->shiftPhase_cmdline;

  // First loop over files specified on the command line
  if(internal_application_cmdline_CurFilename < internal_application_cmdline_Nrfilenames) {
    return argv[internal_application_cmdline_FilenameList[internal_application_cmdline_CurFilename++]];
  }

  // Then loop over the ones in the specified file list
  if(application->filelist != 0) {
    if(internal_application_filelist_CurFilename <  internal_application_filelist_Nrfilenames) {

      if(internal_application_filelist_CurFilename == 0)
	internal_application_filelist_filename = malloc(2048);
      if(ascii_file_get_next_line(internal_application_filelist_fptr, internal_application_filelist_filename, 2048, '#', verbose) == 0) {
	fflush(stdout);
	printerror(application->verbose_state.debug, "ERROR getNextFilenameFromList: reading from %s failed", argv[application->filelist]);
	return 0;
      }
      // Take the first word as the filename
      for(i = 0; i < strlen(internal_application_filelist_filename); i++) {
	if(internal_application_filelist_filename[i] == '\n') {
	  internal_application_filelist_filename[i] = 0;
	  break;
	}else if(internal_application_filelist_filename[i] == '\r') {
	  internal_application_filelist_filename[i] = 0;
	  break;
	}else if(internal_application_filelist_filename[i] == ' ') {
	  internal_application_filelist_filename[i] = 0;
	  break;
	}
      }
      //      printf("XXXXXX next file is %s\n", internal_application_filelist_filename);
      internal_application_filelist_CurFilename++;
      return internal_application_filelist_filename;
    }
  }

  return NULL;
}

void showlibraryversioninformation(FILE *stream)
{
  fprintf(stream, "cfitsio version: ");  
  print_fitsio_version_used(stream);
  fprintf(stream, "\n");  

  fprintf(stream, "fftw version: ");  
  print_fftw_version_used(stream);
  fprintf(stream, "\n");

  fprintf(stream, "GSL version: ");  
  print_gsl_version_used(stream);
  fprintf(stream, "\n");

  fprintf(stream, "pgplot version: ");  
  print_pgplot_version_used(stream);
  fprintf(stream, "\n");
}

//START REGION DEVELOP
