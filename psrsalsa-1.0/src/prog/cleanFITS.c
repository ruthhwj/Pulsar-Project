#include <stdio.h>
#include <string.h>
#include <math.h>
#include "psrsalsa.h"

int main(int argc, char **argv)
{ 
  char *filename_ptr;
  int i, cleanspikes, cleanspikes_mode, cleanspikes_extrabins;
  int index, keepHeader;
  float cleanspikes_possig, cleanspikes_negsig;
  psrsalsaApplication application;
  datafile_definition datain;


  initApplication(&application, "cleanFITS", "[options] inputfile");

  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_nocounters = 1;
  //  application.switch_history_cmd_only = 1;
  application.switch_ext = 1;
  application.switch_output = 1;
  application.switch_debase_slope = 1;

  // Set the default output format
  application.oformat = FITS_format;

  cleanspikes = 0;
  keepHeader = 0;

  if(argc < 2) {
    printf("Program to RFI zap psrFITS files in more advanced ways than in possible with PSRCHIVE.\n\n");
    printApplicationHelp(&application);
    printf("Other options:\n");
    printf("  -cleanspikes    \"mode possig negsig extrabins\"\n");
    printf("                  This is an advanced version for zero DM-ing, where a \n");
    printf("                  non-dedispersed weighted frequency-scrunched template is\n");
    printf("                  removed from the data to surpress non-dedispersed RFI. A\n");
    printf("                  significance possig and negsig can be defined to select time\n");
    printf("                  intervals which are corrected. If the (freq-averaged) signal\n");
    printf("                  is exceeding possig times the rms the time sample is\n");
    printf("                  corrected, or if the signal is a dip exceeding negsig times\n");
    printf("                  the rms. Set possig and/or negsig to a negative value to\n");
    printf("                  disable removing spikes/dips from the data.\n");
    printf("                  mode=1: Subtract freq-avrg from each channel\n");
    printf("                  mode=2: Set each channel to its avrg\n");
    printf("                  At either side of the found signal extrabins extra bins are\n");
    printf("                  removed.\n");
    printf("  -keepHeader     By default the file generated will not be identical in\n");
    printf("                  terms of the header information, and shouldn't be used for\n"); 
    printf("                  timing purposes. By using this option a system call to\n");
    printf("                  fiddleFits will be done, which result in a PSRFITS file which\n");
    printf("                  is identical to the original, with only the data being\n");
    printf("                  replaced.\n");
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-cleanspikes") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %f %f %d", &cleanspikes_mode, &cleanspikes_possig, &cleanspikes_negsig, &cleanspikes_extrabins, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR: Error parsing option '%s'", argv[i]);
	  return 0;
	}
	cleanspikes = 1;
        i++;
      }else if(strcasecmp(argv[i], "-keepheader") == 0) {
	keepHeader = 1;
      }else {
	/* If the option is not recognized, assume it is a filename. It will be added to a list. It gives and error if the "file" starts with -, in which case it probably is a mistyped command line option rather than a file name */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "ERROR: Unknown option: %s\n\nRun cleanFITS without command line arguments to show help", argv[i]);
	  terminateApplication(&application);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0) {
	    return 0;
	  }
	}
      }
    }
  }

  // Checks if all the unrecognised options (should be files) appear as a list at the end of the command line
  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }

  // Check the number of files specified. In this example we want only one file.
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) < 1) {
    printerror(application.verbose_state.debug, "ERROR: Please specify at least one input file on the command line.");
    return 0;
  }


  if(cleanspikes == 0 && application.dodebase_slope == 0) {
    printerror(application.verbose_state.debug, "ERROR: No action was specified to be applied on the input data.");
    return 0;
  }

  // Loop over all the files.
  while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {

    // Open a file in a non-specified format (will be hopefully figured out automatically), writing is not enabled, and the file is read into memory.
    int nowarnings = 2;
    if(application.verbose_state.verbose) {
      nowarnings = 0;
    }
    //    printf("XXXX %d\n", nowarnings);
    // Do not read in memory, as we want to use determineWeightsStat(). If read in memory by openPSRData the memory will be released.
    if(openPSRData(&datain, filename_ptr, FITS_format, 0, 0, nowarnings, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR: Cannot open datafile '%s'.", filename_ptr);
      return 0;
    }
    //    int zeroweightfound, differentweights, negativeweights;
    //    float weightvalue;
    if(readHeaderPSRData(&datain, 0, nowarnings, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR: Cannot read header of '%s'.", filename_ptr);
      closePSRData(&datain, 0, application.verbose_state);
      return 0;
    }
    //    determineWeightsStat(&datain, &zeroweightfound, &differentweights, &negativeweights, &weightvalue);

    long datasize = datain.NrSubints*datain.NrBins*datain.NrPols*datain.NrFreqChan*sizeof(float);
    datain.data = (float *)malloc(datasize);
    if(datain.data == NULL) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR: Cannot allocate memory (data=%ld bytes=%.3fGB).", datasize, datasize/1073741824.0);
      closePSRData(&datain, 0, application.verbose_state);
      return 0;
    }

    if(readPSRData(&datain, datain.data, application.verbose_state)) {
      closePSRData(&datain, 2, application.verbose_state);  // Perserve header information+data, as the data will be available in MEMORY_format
      datain.format = MEMORY_format;
      datain.opened_flag = 1;
    }


    if(preprocessApplication(&application, &datain) == 0) {
      return 0;
    }
    
    if(cleanspikes) {
      if(datain.isDeDisp != 0) {
	printerror(application.verbose_state.debug, "ERROR: Unexpected dedispersions state for '%s'. Data should not have been dedispersed yet.", filename_ptr);
	return 0;
      }
      if(datain.weight_stats_differentweights) {
	printerror(application.verbose_state.debug, "ERROR: Weights are non-uniform in '%s'. Execution is stopped to prefent this information to be lost.", filename_ptr);
	return 0;
      }
      if(datain.weight_stats_zeroweightfound) {
	printwarning(application.verbose_state.debug, "WARNING: Zero-weighted data will be set to zero in '%s' after this opperation.", filename_ptr);
      }
      int orig_poltype = datain.poltype;
      if(orig_poltype == POLTYPE_COHERENCY) {
	if(preprocess_stokes(&datain, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR: Forming Stokes parameters failed for '%s'.", filename_ptr);
	  return 0;
	}
      }else if(orig_poltype == POLTYPE_STOKES) {
	if(application.verbose_state.verbose) {
	  printf("Data is already in Stokes, so no conversion is needed.\n");
	}
      }else {
	printerror(application.verbose_state.debug, "ERROR: Unexpected polarization state for '%s'.", filename_ptr);
	return 0;
      }
      // The data is modified
      if(preprocess_zero_dming(&datain, cleanspikes_possig, cleanspikes_negsig, cleanspikes_mode, -1, cleanspikes_extrabins, 0, 0, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR: Error opening data");
	return 0;
      }
      if(orig_poltype == POLTYPE_COHERENCY) {
	if(preprocess_coherency(&datain, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR: Forming coherency parameters failed for '%s'.", filename_ptr);
	  return 0;
	}
      }
    }

    
    char outputname[MaxFilenameLength], outputnametmp[MaxFilenameLength];
    if(getOutputName(&application, filename_ptr, outputname, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR: Changing filename failed");
      return 0;
    }
    if(keepHeader) {
      strcpy(outputnametmp, outputname);
      strcat(outputnametmp, ".cleanFITStmp54162");
    }else {
      strcpy(outputnametmp, outputname);
    }

    // Open an actual file, enabling writing, on the harddisk with the given name and the given format
    datafile_definition dataout;

    // When writing out data, start with initialisng some variables in the output data file struct
    cleanPSRData(&dataout, application.verbose_state);

    // Copy all header parameters to the new file
    // Note that the data will not be copied.
    // Note also that because the data is already read into memory, the output data will automatically be something that lives in memory as well.
    // After we constructed all required information, the file will be written to disk
    copy_params_PSRData(datain, &dataout, application.verbose_state);

    if(!openPSRData(&dataout, outputnametmp, application.oformat, 1, 0, 0, application.verbose_state))
      return 0;

    // Write out the header information
    if(writeHeaderPSRData(&dataout, argc, argv, 0, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR: Unable to write header.\n");
      return 0;
    }

    // Writes out all data which was previously already was put in memory to disk
    if(writePSRData(&dataout, datain.data, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR: Unable to write data.\n");
      return 0;
    }

    // Close the output file
    // Note that for some data formats (ascii I think), writing actually only happens when closing the file.
    closePSRData(&dataout, 0, application.verbose_state);

    // Finally, close the file, which will also free all associated memory.
    closePSRData(&datain, 0, application.verbose_state);

    if(keepHeader) {
      char cmd[MaxFilenameLength];
      sprintf(cmd, "fiddleFITS ");
      if(application.verbose_state.verbose) {
	strcat(cmd, "-v ");
      }
      strcat(cmd, "-replacesubint ");
      strcat(cmd, outputnametmp);
      strcat(cmd, " -overwrite ");
      strcat(cmd, outputname);

      if(application.verbose_state.verbose) {
	printf("Making copy of input data to: %s\n", outputname);
      }
      if(cp(filename_ptr, outputname, application.verbose_state) != 0) {
	printerror(application.verbose_state.debug, "ERROR cleanFITS: Copy failed (%s to %s).", filename_ptr, outputname);
	return 0;
      }
      if(application.verbose_state.verbose) {
	printf("Running: %s\n", cmd);
      }
      system(cmd);
      if(application.verbose_state.verbose) {
	printf("Removing: %s\n", outputnametmp);
      }
      remove(outputnametmp);
    }

  }

  terminateApplication(&application);
  return 0;
}
