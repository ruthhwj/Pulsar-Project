/* Make sure that programs can handle large files */
#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "psrsalsa.h"

void SHOWREVISIONINFO_prog() {
#include "pcalc.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);

psrsalsaApplication application;

void print_help();

int main(int argc, char **argv)
{
  datafile_definition fin1, fin2, fout;
  int read_wholefile, iformat1, iformat2, oformat, mask_set, ok_flag, index, normal0;
  long i, n, b, p, f, pulse1, pulse2;
  char outputname[MaxFilenameLength], *extension, *file1, *file2, *expression;
  char variables[2][100];
  float samples[2], answer, mask_value, mask_limit, *input_array1, *input_array2;

  initApplication(&application, "pcalc", "[options] -i inputfile1 -i inputfile2");
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
  application.switch_stokes = 1;
  application.switch_coherence = 1;
  application.switch_formatlist = 1;
  application.switch_verbose = 1;
  application.switch_nocounters = 1;
  application.switch_debug = 1;
  application.switch_showrevision = 1;
  application.switch_polselect = 1;
  application.switch_header = 1;
  application.switch_headerlist = 1;
  application.switch_history_cmd_only = 1;

  read_wholefile = 1;
  iformat1 = iformat2 = oformat = -1;
  sprintf(outputname, "output.dat");
  extension = NULL;
  file1 = file2 = NULL;
  expression = NULL;
  pulse1 = pulse2 = -1;
  mask_set = 0;
  normal0 = 0;

  /* Clear struct */
  //  cleanPSRData(&fin1, application.verbose_state);
  //  cleanPSRData(&fin2, application.verbose_state);
  cleanPSRData(&fout, application.verbose_state);

  if(argc <= 1) {
    print_help();
    closePSRData(&fout, 0, application.verbose_state);
    terminateApplication(&application);
    return 0;
  }

  for(i = 1; i < argc; i++) {
    index = i;
    if(processCommandLine(&application, argc, argv, &index)) {
      i = index;
    }else if(strcmp(argv[i], "-iformat1") == 0) {
      iformat1 = parsePSRDataFormats(argv[++(i)]);
      if(iformat1 == 0)
	return 0;
    }else if(strcmp(argv[i], "-iformat2") == 0) {
      iformat2 = parsePSRDataFormats(argv[++(i)]);
      if(iformat2 == 0)
	return 0;
    }else if(strcmp(argv[i], "-oformat") == 0) {
      oformat = parsePSRDataFormats(argv[++(i)]);
      if(oformat == 0)
	return 0;
    }else if(strcmp(argv[i], "-output") == 0) {
      strcpy(outputname, argv[++i]);
    }else if(strcmp(argv[i], "-i") == 0) {
      if(file1 == NULL)
	file1 = argv[++i];
      else if(file2 == NULL)
	file2 = argv[++i];
      else {
	printerror(application.verbose_state.debug, "pcalc: You can only specify two input files at most.\n");
	return 0;
      }
    }else if(strcmp(argv[i], "-ext") == 0) {
      extension = argv[++i];
    }else if(strcmp(argv[i], "-expr") == 0) {
      expression = argv[++i];
    }else if(strcmp(argv[i], "-memsave") == 0) {
      read_wholefile = 0;
    }else if(strcmp(argv[i], "-normal0") == 0) {
      normal0 = 1;
    }else if(strcmp(argv[i], "-prange") == 0) {
      n = sscanf(argv[i+1], "%ld %ld", &pulse1, &pulse2);
      if(n != 2) {
	printerror(application.verbose_state.debug, "pcalc: Error parsing -prange option\n");
	return 0;
      }
      i++;
    }else if(strcmp(argv[i], "-mask") == 0) {
      n = sscanf(argv[i+1], "%f %f", &mask_limit, &mask_value);
      if(n != 2) {
	printerror(application.verbose_state.debug, "pcalc: Error parsing -mask option\n");
	return 0;
      }
      mask_set = 1;
      i++;
    }else {
      printerror(application.verbose_state.debug, "pcalc: Unknown option (%s). Use -i twice to specify file names. Run pcalc without command line options to display the help.\n", argv[i]);
      closePSRData(&fout, 0, application.verbose_state);
      terminateApplication(&application);
      return 0;
    }
  }

  if(expression == NULL) {
    printerror(application.verbose_state.debug, "pcalc: Please specify an expression to calculate with the -expr option.\n\n");
    return 0;
  }
  if(file1 == NULL) {
    printerror(application.verbose_state.debug, "pcalc: Please specify at least one input file with the -i option.\n\n");
    return 0;
  }
  if(strstr(expression, "y") != NULL && file2 == NULL) {
    printerror(application.verbose_state.debug, "pcalc: Expression specified with -exp suggests that there should be two input files, but only one is specified with the -i option.\n\n");
    return 0;
  }


  if(iformat1 <= 0 && file1 != NULL)
    iformat1 = guessPSRData_format(file1, 0, application.verbose_state);
  if(iformat2 <= 0 && file2 != NULL) 
    iformat2 = guessPSRData_format(file2, 0, application.verbose_state);
  if(file1 != NULL) {
    if(iformat1 == -2 || iformat1 == -3)
      return 0;
    if(isValidPSRDATA_format(iformat1) == 0) {
      printerror(application.verbose_state.debug, "ERROR pcalc: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat1 option if the format is supported, but not automatically recognized.\n\n", file1);
      return 0;
    }
  }
  if(file2 != NULL) {
    if(iformat2 == -2 || iformat1 == -3)
      return 0;
    if(isValidPSRDATA_format(iformat2) == 0) {
      printerror(application.verbose_state.debug, "ERROR pcalc: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat2 option if the format is supported, but not automatically recognized.\n\n", file2);
      return 0;
    }
  }

  /* Open input file and read header */
  if(!openPSRData(&fin1, file1, iformat1, 0, read_wholefile, 0, application.verbose_state))
    return 0;
  if(oformat <= 0) {
    oformat = iformat1;
  }
  if(isValidPSRDATA_format(oformat) == 0) {
    printerror(application.verbose_state.debug, "ERROR pcalc: Please specify a valid output format with the -oformat option.\n\n");
    return 0;
  }
  if(oformat == FITS_format) {
    printwarning(application.verbose_state.debug, "WARNING pcalc: Be aware that when the data is written in psrfits format, there can be significant rounding errors when the the range of values generated is very large (for example when you have division by noise which is on average zero). You might consider using -oformat psrsalsa.");
  }
  if(read_wholefile == 0) {
    if(!readHeaderPSRData(&fin1, 0, 0, application.verbose_state))
      return 0;
    /* Search commandline for header parameters to overwrite header
     parameters. */
    if(PSRDataHeader_parse_commandline(&fin1, argc, argv, application.verbose_state) == 0)
      return 0;
    printwarning(application.verbose_state.debug, "WARNING pcalc: Pre-process options will be ignored when not whole file is read into memory.");
  }else {
    /* Search commandline for header parameters to overwrite header
     parameters. */
    if(PSRDataHeader_parse_commandline(&fin1, argc, argv, application.verbose_state) == 0)
      return 0;
    for(i = 1; i < argc; i++) {
      if(strcmp(argv[i], "-header") == 0) {
	printwarning(application.verbose_state.debug, "WARNING: If using the -header option, be aware it applied BEFORE the preprocessing.");
      }
    }
    if(preprocessApplication(&application, &fin1) == 0) {
      return 0;
    }
  }
  if(file2 != NULL) {
    if(!openPSRData(&fin2, file2, iformat2, 0, read_wholefile, 0, application.verbose_state))
      return 0;
    if(read_wholefile == 0) {
      /* Search commandline for header parameters to overwrite header
	 parameters. */
      if(PSRDataHeader_parse_commandline(&fin2, argc, argv, application.verbose_state) == 0)
      return 0;
      if(!readHeaderPSRData(&fin2, 0, 0, application.verbose_state))
	return 0;
    }else {
      /* Search commandline for header parameters to overwrite header
	 parameters. */
      if(PSRDataHeader_parse_commandline(&fin2, argc, argv, application.verbose_state) == 0)
	return 0;
      if(preprocessApplication(&application, &fin2) == 0) {
	return 0;
      }
    }
    if(fin1.NrSubints != fin2.NrSubints) {
      printerror(application.verbose_state.debug, "pcalc: Number of pulses are different in the input files.\n");
      closePSRData(&fin1, 0, application.verbose_state);
      closePSRData(&fin2, 0, application.verbose_state);
      closePSRData(&fout, 0, application.verbose_state);
      terminateApplication(&application);
      return 0;      
    }
    if(fin1.NrBins != fin2.NrBins) {
      printerror(application.verbose_state.debug, "pcalc: Number of bins are different in the input files.\n");
      closePSRData(&fin1, 0, application.verbose_state);
      closePSRData(&fin2, 0, application.verbose_state);
      closePSRData(&fout, 0, application.verbose_state);
      terminateApplication(&application);
      return 0;      
    }
    if(fin1.NrPols != fin2.NrPols) {
      printerror(application.verbose_state.debug, "pcalc: Number of polarizations are different in the input files.\n");
      closePSRData(&fin1, 0, application.verbose_state);
      closePSRData(&fin2, 0, application.verbose_state);
      closePSRData(&fout, 0, application.verbose_state);
      terminateApplication(&application);
      return 0;      
    }
    if(fin1.NrFreqChan != fin2.NrFreqChan) {
      printerror(application.verbose_state.debug, "pcalc: Number of frequency channels are different in the input files.\n");
      closePSRData(&fin1, 0, application.verbose_state);
      closePSRData(&fin2, 0, application.verbose_state);
      closePSRData(&fout, 0, application.verbose_state);
      terminateApplication(&application);
      return 0;      
    }
  }

  /* Copy header parameters to output header */
  if(pulse1 < 0) {
    pulse1 = 0;
    pulse2 = fin1.NrSubints -1;
  }
  copy_params_PSRData(fin1, &fout, application.verbose_state);
  fout.NrSubints = pulse2 - pulse1 + 1;

  if(extension != NULL) {
    if(change_filename_extension(file1, outputname, extension, 10000, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "pcalc: Changing filename failed\n");
      return 0;
    }
  }

  /* Open output file and write data */
  if(!openPSRData(&fout, outputname, oformat, 1, 0, 0, application.verbose_state))
    return 0;
  if(!writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state))
    return 0;
  //  appendHistoryLine(fout, argc, argv, application.verbose_state);

  sprintf(variables[0], "x");
  sprintf(variables[1], "y");


  input_array1 = (float *)malloc(fin1.NrBins*sizeof(float));
  if(input_array1 == NULL) {
    printerror(application.verbose_state.debug, "pcalc: Memory allocation error (%ld bytes)\n", fin1.NrBins*sizeof(float));
    return 0;
  }
  if(file2 != NULL) {
    input_array2 = (float *)malloc(fin2.NrBins*sizeof(float));
    if(input_array2 == NULL) {
      printerror(application.verbose_state.debug, "pcalc: Memory allocation error (%ld bytes)\n", fin2.NrBins*sizeof(float));
      return 0;
    }
  }

  /*
  samples[0] = 1.1;
  samples[1] = 2.2;
  calc_expressionf(expression, 2, variables, samples, &answer, 2);
  return 0;
  */
  for(n = pulse1; n <= pulse2; n++) {
    for(p = 0; p < fin1.NrPols; p++) {
      for(f = 0; f < fin1.NrFreqChan; f++) {
	if(readPulsePSRData(&fin1, n, p, f, 0, fin1.NrBins, input_array1, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "pcalc: Cannot read data\n");
	  return 0;
	}
	if(file2 != NULL) {
	  if(readPulsePSRData(&fin2, n, p, f, 0, fin2.NrBins, input_array2, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "pcalc: Cannot read data\n");
	    return 0;
	  }
	}
	for(b = 0; b < fin1.NrBins; b++) {
	  samples[0] = input_array1[b];
	  if(file2 != NULL) {
	    samples[1] = input_array2[b];
	  }else {
	    samples[1] = 0;
	  }

	  ok_flag = 1;
	  if(mask_set) {
	    if(samples[0] <= mask_limit)
	      ok_flag = 0;
	    if(samples[1] <= mask_limit)
	      ok_flag = 0;
	  }
	  if(ok_flag) {
	    if(calc_expressionf(expression, 2, variables, samples, &answer, 0) == 0) {
	      printerror(application.verbose_state.debug, "pcalc: Cannot calculate expression\n");
	      return 0;
	    }
	    //	    if(b == 500) {
	    //	      printf("%s: vars=%s and %s (%f and %f) gives %f\n", expression, variables[0], variables[1], samples[0], samples[1], answer);
	    //	    }
	  }else {
	    answer = mask_value;
	  }
	  if(normal0) {
	    if(isnan(answer))
	      answer = 0;
	    else if(isinf(answer))
	      answer = 0;
	  }
	  input_array1[b] = answer;
	}
	//	printf("subint %ld pol %ld freq %ld: wrote %f\n", n, p-pulse1, f, input_array1[500]);
	if(writePulsePSRData(&fout, n, p-pulse1, f, 0, fin1.NrBins, input_array1, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "pcalc: Cannot write data\n");
	  return 0;
	}
      }
    }
    b = (n)/10;
    if(application.verbose_state.verbose && 10*b == n) {
      if(application.verbose_state.nocounters == 0) {
	printf("processed %.1f%%     \r",(100.0*(n))/(float)(fin1.NrSubints));
	fflush(stdout);
      }
    }
  }
  if(application.verbose_state.verbose) printf("  Done                                      \n");


  /* Free resources and close files */
  closePSRData(&fin1, 0, application.verbose_state);
  if(file2 != NULL)
    closePSRData(&fin2, 0, application.verbose_state);
  closePSRData(&fout, 0, application.verbose_state);


  free(input_array1);
  if(file2 != NULL) {
    free(input_array2);
  }

  terminateApplication(&application);
  return 0;
}

void print_help()
{
  fprintf(stdout, "Program to combine two sets of pulsar data in various mathematical ways. Usage:\n\n");
  printApplicationHelp(&application);
  fprintf(stdout, "\nOther options:\n");
  fprintf(stdout, "-expr        The calculation. For instance \'x+y\' to get the sum of the two inputs.\n");
  fprintf(stdout, "-normal0     Filter the output such that nan's and inf's are set to zero.\n");
  fprintf(stdout, "-mask        \"limit value\" Ignore all input values below limit, and put value in the output.\n\n");
  fprintf(stdout, "-prange      Select pulse number range.\n");
  fprintf(stdout, "-iformat1    Specify format first input file\n");
  fprintf(stdout, "-iformat2    Specify format second input file\n");
  fprintf(stdout, "-oformat     Specify input format\n");
  fprintf(stdout, "-memsave     Try to use less memory.\n");
  fprintf(stdout, "-ext         Specify output extension\n");
  fprintf(stdout, "-output      instead you can specify output filename, default is output.dat\n");
  //  fprintf(stdout, "\nFormats are:\n");
  //  printPSRDataFormats(stdout, 2);
  fprintf(stdout, "\n");
  fprintf(stdout, "Supported functions are:\n");
  printCalcFunctions(stdout, 2);
  fprintf(stdout, "\n");
}
