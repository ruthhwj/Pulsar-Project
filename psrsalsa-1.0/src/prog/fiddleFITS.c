#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "psrsalsa.h"


extern void internalFITSscalePulse(float *pulse, long nrSamples, float *offset, float *scale, float maxvalue);


int main(int argc, char **argv)
{
  int status, iomode, anynul, whattodo, colnum, colnum_err, fullscales, fullscales2;
  int ncpar, nchan, npol, nbin, overwrite_flag;
  int highpass_nr;
  int obsmode, fitsfile2_name, reset_col_nr, chisq_defined, errors_defined;
  long i, j, k, n, nrows, datalength;
  char card[FLEN_CARD], comment[FLEN_COMMENT], value[FLEN_VALUE], txt[1000];
  float *data, *data2, *data3, *data4, reset_col_value, atten1, atten2;
  int *data_int;
  char output_suffix[MaxFilenameLength], output_name[MaxFilenameLength], input_name[MaxFilenameLength];
  FILE *outputfptr;
  fitsfile *fptr, *fptr2;
  psrsalsaApplication application;

  status = 0;       /* CFITSIO status value MUST be initialized to zero! */
  whattodo = 0;
  outputfptr = stdout;
  fitsfile2_name = -1;
  overwrite_flag = 0;
  strcpy(output_suffix, "fiddle");
  initApplication(&application, "fiddleFITS", "[options] inputfile");

  if(argc <= 2) {
    printf("Usage: fiddleFITS [options] fitsfile\n\n");
    printf("-obsmode CAL/PSR               Change obs mode in the header to CAL or PSR\n");
    printf("-pacv                          Extract calibration table from a pacv/pcm\n");
    printf("                               psrfits file.\n");
    printf("-resetgain                     Overwrite the gain to be 1 for a pacv/pcm style\n");
    printf("                               calibration table\n");
    printf("-resetcaltable \"Colnr Value\"   Reset column (counting from 1) in a a pacv/pcm\n");
    printf("                               style calibration table\n");
    printf("-output                        Output is directed to this text file (use with\n");
    printf("                               -pacv)\n");
    printf("-replacesubint FFread          Replaces subint data in fitsfile with that\n");
    printf("                               of FFread\n");
    printf("-atten \"dB_A dB_B\"             Correct the data for a given set of attenuations,\n");
    printf("                               the input file will be overwritten\n");
    printf("-highpass number               Highpass filter which removes the first number of\n");
    printf("                               frequency channels from the data. Each subint,\n");
    printf("                               frequency, and polarization channel is processed\n");
    printf("                               independently. If number is set to a negative\n");
    printf("                               number only the first harmonic is removed without\n");
    printf("                               calculating the full Fourier spectrum, which\n");
    printf("                               might be quicker.\n");
    printf("-v                             Enable verbose\n");
    printf("\nOutput options:\n\n");
    printf("-ext              Specify suffix, default is '%s'\n", output_suffix);
    printf("-overwrite        Overwrite the input file. Be very cautious!\n");
    return 0;
  }
  for(j = 1; j < argc-1; j++) {
    if(strcmp(argv[j], "-v") == 0) {
      application.verbose_state.verbose = 1;
    }else if(strcmp(argv[j], "-pacv") == 0) {
      whattodo = 1;
    }else if(strcmp(argv[j], "-output") == 0) {
      outputfptr = fopen(argv[j+1], "w");
      if(outputfptr == NULL) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot open '%s' for writing.\n", argv[j+1]);
	return 0;
      }
      j++;
    }else if(strcmp(argv[j], "-obsmode") == 0) {
      whattodo = 2;
      if((strcasecmp(argv[j+1],  "CAL") != 0) && (strcasecmp(argv[j+1],  "PSR") != 0)) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: obsmode options are only allowed to be CAL or PSR\n");
	 return 0;
      }
      obsmode = j+1;
      j++;
    }else if(strcmp(argv[j], "-replacesubint") == 0) {
      whattodo = 3;
      fitsfile2_name = j+1;
      j++;
    }else if(strcmp(argv[j], "-resetgain") == 0) {
      whattodo = 4;
    }else if(strcmp(argv[j], "-resetcaltable") == 0) {
      whattodo = 5;
      n = sscanf(argv[j+1], "%d %f", &reset_col_nr, &reset_col_value);
      if(n != 2) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot parse %s, expect two values.\n", argv[j]);
	return 0;
      }
      j++;
    }else if(strcmp(argv[j], "-atten") == 0) {
      whattodo = 6;
      n = sscanf(argv[j+1], "%f %f", &atten1, &atten2);
      if(n != 2) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot parse %s, expect two values.\n", argv[j]);
	return 0;
      }
      j++;
    }else if(strcmp(argv[j], "-andrew") == 0) {
      whattodo = 7;
      highpass_nr = -1;
    }else if(strcmp(argv[j], "-highpass") == 0) {
      whattodo = 7;
      n = sscanf(argv[j+1], "%d", &highpass_nr);
      if(n != 1) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot parse %s, expect one value.\n", argv[j]);
	return 0;
      }
      j++;
    }else if(strcmp(argv[j], "-overwrite") == 0) {
      overwrite_flag = 1;
    }else if(strcmp(argv[j], "-ext") == 0) {
      strcpy(output_suffix, argv[j+1]);
      j++;
    }else {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Unknown option '%s'.\n", argv[j]);
      return 0;
    }
  }

  if(whattodo == 0) {
    printerror(application.verbose_state.debug, "ERROR fiddleFITS: Specify on command line what you want to do with the fits file\n");
    return 0;
  }

  if(overwrite_flag == 0) {
    if(change_filename_extension(argv[argc-1], output_name, output_suffix, MaxFilenameLength, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot change extension in output name.");
      return 0;
    }
    if(application.verbose_state.verbose) printf("Creating copy of %s in %s\n", argv[argc-1], output_name);
    if(cp(argv[argc-1], output_name, application.verbose_state) != 0) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Copy failed (%s to %s).", argv[argc-1], output_name);
      return 0;
    }
    // Use copy as input below
    strcpy(input_name, output_name);
  }else {
    // Use input file directly
    strcpy(input_name, argv[argc-1]);
    if(application.verbose_state.verbose) {
      printwarning(application.verbose_state.debug, "Overwriting input file");
    }
  }

  /* Open files */
  iomode = READWRITE;
  if (!fits_open_file(&fptr, input_name, iomode, &status)) {
    if(application.verbose_state.verbose) printf("'%s' opened\n", input_name);
  }else {
    printf("ERROR fiddleFITS: Cannot open %s\n", input_name);
  }
  if (status) {
    fits_report_error(stderr, status); /* print any error message */
    return 0;
  }
  if(fitsfile2_name > 0) {
    iomode = READONLY;
    if (!fits_open_file(&fptr2, argv[fitsfile2_name], iomode, &status)) {
      if(application.verbose_state.verbose) printf("'%s' opened\n", argv[fitsfile2_name]);
    }else {
      printf("ERROR fiddleFITS: Cannot open %s\n", argv[fitsfile2_name]);
    }
    if (status) {
      fits_report_error(stderr, status); /* print any error message */
      return 0;
    }
  }

  if(whattodo == 1 || whattodo == 4 || whattodo == 5) {                                  /* Extract pacv information or reset gain or another column */
    if(fits_movnam_hdu(fptr, BINARY_TBL, "FEEDPAR", 0, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot move to FEEDPAR HDU.\n");
      return 0;
    } 

    if (fits_read_card(fptr,"NCPAR", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NCPAR keyword does not exist\n");
      return 0;
    }
    if(status) 
      fits_report_error(stderr, status); 
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &ncpar);

    if (fits_read_card(fptr,"NCHAN", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NCHAN keyword does not exist\n");
      return 0;
    }
    if(status) 
      fits_report_error(stderr, status); 
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &nchan);

    if(application.verbose_state.verbose) {
      if(whattodo == 1)
	printf("Extracting %d channels and %d parameters\n", nchan, ncpar);
      else if(whattodo == 4)
	printf("There are %d channels and %d parameters, going to reset the first parameter\n", nchan, ncpar);
      else if(whattodo == 5)
	printf("There are %d channels and %d parameters, going to reset parameter %d\n", nchan, ncpar, reset_col_nr);
    }

    if(whattodo == 5) {
      if(reset_col_nr < 1 || reset_col_nr > ncpar) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: specified column nr is invalid.\n");
	return 0;
      }
    }

    data = calloc(nchan*ncpar, sizeof(float));
    data2 = calloc(nchan, sizeof(float));
    data3 = calloc(nchan, sizeof(float));
    data4 = calloc(nchan*ncpar, sizeof(float));
    data_int = calloc(nchan, sizeof(int));
    if(data == NULL || data2 == NULL || data3 == NULL || data4 == NULL || data_int == NULL) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: cannot allocate memory\n");
      return 0;
    }

    if(fits_get_colnum (fptr, CASEINSEN, "DAT_FREQ", &colnum, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No frequency data in fits file?\n");
      return 0;
    }

    if(fits_read_col(fptr, TFLOAT, colnum, 1, 1, nchan, NULL, data2, &anynul, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read data\n");
      fits_report_error(stderr, status); 
      return 0;
    }

    chisq_defined = 1;
    if(fits_get_colnum (fptr, CASEINSEN, "CHISQ", &colnum, &status)) { 
      printwarning(application.verbose_state.debug, "WARNING fiddleFITS: No chi^2 values defined in fits file?\n");
      chisq_defined = 0;
      status = 0;
    }else {
      if(fits_read_col(fptr, TFLOAT, colnum, 1, 1, nchan, NULL, data3, &anynul, &status)) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read data\n");
	fits_report_error(stderr, status); 
	return 0;
      }

      if(fits_get_colnum (fptr, CASEINSEN, "NFREE", &colnum, &status)) { 
	printwarning(application.verbose_state.debug, "WARNING fiddleFITS: No degrees of freedom values in fits file?\n");
	chisq_defined = 0;
	status = 0;
      }else {
	if(fits_read_col(fptr, TINT, colnum, 1, 1, nchan, NULL, data_int, &anynul, &status)) {
	  printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read data\n");
	  fits_report_error(stderr, status); 
	  return 0;
	}
      }
    }
    errors_defined = 1;
    if(fits_get_colnum (fptr, CASEINSEN, "DATAERR", &colnum_err, &status)) { 
      printwarning(application.verbose_state.debug, "WARNING fiddleFITS: No error values defined in fits file?\n");
      errors_defined = 0;
      status = 0;
    }


    if(fits_get_colnum (fptr, CASEINSEN, "DATA", &colnum, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No data in fits file?\n");
      return 0;
    }

    fits_get_num_rows(fptr, &nrows, &status);
    /*
      int ncols;
    fits_get_num_cols(fptr, &ncols, &status);
    */
    for(i = 0; i < nrows; i++) {
      if(fits_read_col(fptr, TFLOAT, colnum, 1+i, 1, nchan*ncpar, NULL, data, &anynul, &status)) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read data\n");
	fits_report_error(stderr, status); 
	return 0;
      }
      if(errors_defined) {
	if(fits_read_col(fptr, TFLOAT, colnum_err, 1, 1, nchan*ncpar, NULL, data4, &anynul, &status)) {
	  printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read error data\n");
	  fits_report_error(stderr, status); 
	  return 0;
	}
      }
      
      if(whattodo == 1) {
	for(j = 0; j < nchan; j++) {
	  fprintf(outputfptr, "%ld %ld %f ", i, j, data2[j]);
	  if(chisq_defined)
	    fprintf(outputfptr, "%f ", data3[j]/(float)data_int[j]);
	  for(k = 0; k < ncpar; k++) {
	    fprintf(outputfptr, "%f ", data[ncpar*j+k]);
	    if(errors_defined)
	      fprintf(outputfptr, "%f ", data4[ncpar*j+k]);
	  }
	  fprintf(outputfptr, "\n");
	}
      }else {
	for(j = 0; j < nchan; j++) {   /* Assume gain is in the first column */
	  if(whattodo == 4)
	    data[ncpar*j+0] = 1;
	  else
	    data[ncpar*j+reset_col_nr-1] = reset_col_value;
	}	
	if(fits_write_col(fptr, TFLOAT, colnum, 1+i, 1, nchan*ncpar, data, &status) != 0) {
	  printerror(application.verbose_state.debug, "ERROR fiddleFITS: Error writing data.\n");
	  fits_report_error(stderr, status); 
	  return 0;
	}
      }
    }
    if(whattodo == 1) {
      printf("Output is subint number, channel number, frequency, ");
      if(chisq_defined)
	printf("reduced chi^2, ");
      printf("polarization calibration output");
      if(errors_defined)
	printf(" (value/error pairs)");
      printf("\n");
      for(i = 0; i < ncpar; i++) {
	sprintf(txt, "PAR_%04ld", i);
	if (fits_read_card(fptr,txt, card, &status)) {
	  printwarning(application.verbose_state.debug, "WARNING fiddleFITS: %s keyword does not exist\n", txt);
	  status = 0;
	}else {
	  fits_parse_value(card, value, comment, &status);
	  if(strcmp(value, "'gamma   '") == 0) {
	    printf("%s %s   Note: gamma=g0/g1-1=exp(2beta)-1\n", value, comment);
	  }else {
	    printf("%s %s\n", value, comment);
	  }
	}
      }
    }
    free(data);
    free(data2);
    free(data3);
    free(data4);
    free(data_int);
  }else if(whattodo == 2){ /* change obs_mode  */
    if (fits_read_card(fptr,"OBS_MODE", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: OBS_MODE keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    if(status) {
      fits_report_error(stderr, status); 
      return 0;
    }
    char obsmodestr[100];
    if(strcasecmp(argv[obsmode],  "CAL") == 0) {
      sprintf(obsmodestr, "CAL");
    }else {
      sprintf(obsmodestr, "PSR");
    }
    printf("current OBS_MODE is %s\n", value);
    printf("changing OBS_MODE to %s\n", obsmodestr);
    if(fits_update_key(fptr, TSTRING, "OBS_MODE", obsmodestr, "(PSR,CAL, SEARCH)", &status)){
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot change OBS_MODE\n");
      fits_report_error(stderr, status); 
      return 0;
    }
  }else if(whattodo == 3) {  /* Replace subint data */
    int colnum_s, colnum_o, colnum_w, colnum_d, colnum_s2, colnum_o2, colnum_w2, colnum_d2;
    if(outputfptr != stdout) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: You cannot use -output together with -replacesubint\n");
      return 0;
    }
    /* Read in header parameters from output file */
    if(fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: SUBINT table does not exist!\n");
      return 0;
    } 
    if (fits_read_card(fptr,"NPOL", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NPOL keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &npol);
    if (fits_read_card(fptr,"NBIN", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NBIN keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &nbin);
    if (fits_read_card(fptr,"NCHAN", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NCHAN keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &nchan);
    if (fits_read_card(fptr,"NAXIS2", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NAXIS2 keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%ld", &nrows);
    if(fits_get_colnum (fptr, CASEINSEN, "DAT_SCL", &colnum_s, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No scales in output fits file?\n");
      return 0;
    }
    if(fits_get_colnum (fptr, CASEINSEN, "DAT_OFFS", &colnum_o, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No offsets in output fits file?\n");
      return 0;
    }
    if(fits_get_colnum (fptr, CASEINSEN, "DAT_WTS", &colnum_w, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No weights in output fits file?\n");
      return 0;
    }
    if(fits_get_colnum (fptr, CASEINSEN, "DATA", &colnum_d, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No data in output fits file?\n");
      return 0;
    }
    sprintf(txt, "TFORM%d", colnum_s);
    if (fits_read_card(fptr,txt, card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: keyword %s does not exist\n", txt);
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    if(value[1] == 'E') {   /* Simply 1 float instead of an array which is written as '(number)E' */
      datalength = 1;
    }else {
      sscanf(value, "'%ld'", &datalength);
    }
    if(application.verbose_state.verbose)
      printf("Output datalength = %ld\n", datalength);
    if(datalength == npol*nchan) {
      fullscales = 1;
    }else if(datalength == nchan) {
      fullscales = 0;
    }else {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: %ld is a weird number for the number of scales.\n", datalength);
      return 0;
    }

    if(application.verbose_state.verbose) printf("Output data (%s) contains %ld subints, %d freq channels, %d pol channels and %d phase bins (fullscales=%d)\n", input_name, nrows, nchan, npol, nbin, fullscales);
  
    /* Check the header parameters of the input file */
    if(fits_movnam_hdu(fptr2, BINARY_TBL, "SUBINT", 0, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: SUBINT table does not exist!\n");
      return 0;
    } 
    if (fits_read_card(fptr2,"NPOL", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NPOL keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%ld", &i);
    if(npol != i) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Nr of polarizations doesn't match (%ld != %d)\n", i, npol);
      return 0;
    }
    if (fits_read_card(fptr2,"NBIN", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NBIN keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%ld", &i);
    if(nbin != i) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Nr of phase bins doesn't match (%ld != %d)\n", i, nbin);
      return 0;
    }
    if (fits_read_card(fptr2,"NCHAN", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NCHAN keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%ld", &i);
    if(nchan != i) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Nr of frequency channels doesn't match (%ld != %d)\n", i, nchan);
      return 0;
    }
    if (fits_read_card(fptr2,"NAXIS2", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NAXIS2 keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%ld", &i);
    if(nrows != i) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Nr of subints doesn't match (%ld != %ld)\n", i, nrows);
      return 0;
    }
    if(fits_get_colnum (fptr2, CASEINSEN, "DAT_SCL", &colnum_s2, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No scales in input fits file?\n");
      return 0;
    }
    if(fits_get_colnum (fptr2, CASEINSEN, "DAT_OFFS", &colnum_o2, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No offsets in input fits file?\n");
      return 0;
    }
    if(fits_get_colnum (fptr2, CASEINSEN, "DAT_WTS", &colnum_w2, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No weights in input fits file?\n");
      return 0;
    }
    if(fits_get_colnum (fptr2, CASEINSEN, "DATA", &colnum_d2, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No data in input fits file?\n");
      return 0;
    }
    sprintf(txt, "TFORM%d", colnum_s2);
    if (fits_read_card(fptr2,txt, card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: keyword %s does not exist\n", txt);
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    if(value[1] == 'E') {   /* Simply 1 float instead of an array which is written as '(number)E' */
      datalength = 1;
    }else {
      sscanf(value, "'%ld'", &datalength);
    }
    if(application.verbose_state.verbose)
      printf("Input datalength = %ld\n", datalength);
    if(datalength == npol*nchan) {
      fullscales2 = 1;
    }else if(datalength == nchan) {
      fullscales2 = 0;
    }else {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: %ld is a weird number for the number of scales in %s.\n", datalength, txt);
      return 0;
    }
    if(fullscales != fullscales2) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Nr of scales does not match\n");
      return 0;
    }
    if(application.verbose_state.verbose) printf("Header parameters of input data matches the output data\n");

    i = nchan*npol*nbin*sizeof(int);
    j = nrows*npol*nchan*sizeof(float);
    if(j > i)
      i = j;
    data = malloc(i);
    if(data == NULL) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Memory allocation failed\n");
      return 0;
    }

    if(application.verbose_state.verbose) printf("Start transferring data, scales, offsets and weights\n");
    for(n = 0; n < nrows; n++) {
      if(application.verbose_state.verbose) printf("Transferring subint %ld\n", n+1);
      if(fits_read_col(fptr2, TFLOAT, colnum_s2, 1+n, 1, nchan*(1+fullscales*(npol-1)), NULL, data, &anynul, &status)) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read data\n");
	fits_report_error(stderr, status); 
	return 0;
      }
      if(fits_write_col(fptr, TFLOAT, colnum_s, 1+n, 1, nchan*(1+fullscales*(npol-1)), data, &status) != 0) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Error writing data.\n");
	fits_report_error(stderr, status); 
	return 0;
      }
      if(fits_read_col(fptr2, TFLOAT, colnum_o2, 1+n, 1, nchan*(1+fullscales*(npol-1)), NULL, data, &anynul, &status)) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read data\n");
	fits_report_error(stderr, status); 
	return 0;
      }
      if(fits_write_col(fptr, TFLOAT, colnum_o, 1+n, 1, nchan*(1+fullscales*(npol-1)), data, &status) != 0) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Error writing data.\n");
	fits_report_error(stderr, status); 
	return 0;
      }
      if(fits_read_col(fptr2, TFLOAT, colnum_w2, 1+n, 1, nchan, NULL, data, &anynul, &status)) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read data\n");
	fits_report_error(stderr, status); 
	return 0;
      }
      if(fits_write_col(fptr, TFLOAT, colnum_w, 1+n, 1, nchan, data, &status) != 0) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Error writing data.\n");
	fits_report_error(stderr, status); 
	return 0;
      }

      if(fits_read_col(fptr2, TINT, colnum_d2, 1+n, 1, nchan*npol*nbin, NULL, data, &anynul, &status)) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read data\n");
	fits_report_error(stderr, status); 
	return 0;
      }
      if(fits_write_col(fptr, TINT, colnum_d, 1+n, 1, nchan*npol*nbin, data, &status)) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read data\n");
	fits_report_error(stderr, status); 
	return 0;
      }

      if (status) {
	fits_report_error(stderr, status); 
	return 0;
      }
    }
    if(application.verbose_state.verbose) printf("Transferring done\n");


    free(data);
  }else if(whattodo == 6) {  /* Change attenuation levels */
    int colnum_s;
    if(outputfptr != stdout) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: You cannot use -output together with -atten\n");
      return 0;
    }
    /* Read in header parameters from output file */
    if(fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: SUBINT table does not exist!\n");
      return 0;
    } 
    if (fits_read_card(fptr,"NPOL", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NPOL keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &npol);
    if(npol != 4) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Expect 4 polarization when changing the attenuation levels.\n");
      return 0;
    }

    if (fits_read_card(fptr,"POL_TYPE", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddlefits: POL_TYPE keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "'%s'", txt);
    for(i = strlen(txt)-1; i >= 0; i--) {
      if(i >= 0) {
	if(txt[i] == '\'')
	  txt[i] = 0;
      }
    }
    if(strcmp(txt, "AABBCRCI") != 0) {
      printerror(application.verbose_state.debug, "ERROR fiddlefits: Coherency parameters were expected, got %s\n", txt);
      return 0;
    }

    if (fits_read_card(fptr,"NCHAN", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NCHAN keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &nchan);
    if (fits_read_card(fptr,"NAXIS2", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NAXIS2 keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%ld", &nrows);
    if(fits_get_colnum (fptr, CASEINSEN, "DAT_SCL", &colnum_s, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No scales in output fits file?\n");
      return 0;
    }
    sprintf(txt, "TFORM%d", colnum_s);
    if (fits_read_card(fptr,txt, card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: keyword %s does not exist\n", txt);
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    if(value[1] == 'E') {   /* Simply 1 float instead of an array which is written as '(number)E' */
      datalength = 1;
    }else {
      sscanf(value, "'%ld'", &datalength);
    }
    if(application.verbose_state.verbose)
      printf("Output datalength = %ld\n", datalength);
    if(datalength == npol*nchan) {
      fullscales = 1;
    }else if(datalength == nchan) {
      fullscales = 0;
    }else {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: %ld is a weird number for the number of scales.\n", datalength);
      return 0;
    }
    if(fullscales == 0) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Changing attenuation settings is only working when each polarization has its own scale factor.\n");
      return 0;
    }

    if(application.verbose_state.verbose) printf("Output data (%s) contains %ld subints, %d freq channels and %d pol channels (fullscales=%d)\n", input_name, nrows, nchan, npol, fullscales);
  
    i = nrows*npol*nchan*sizeof(float);
    data = malloc(i);
    if(data == NULL) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Memory allocation failed\n");
      return 0;
    }

    if(application.verbose_state.verbose) printf("Start adjusting scales\n");
    float *scales, fac1, fac2;
    long f;
    scales = data;
    fac1 = pow(10.0, atten1/20.0);
    fac2 = pow(10.0, atten2/20.0);
    for(n = 0; n < nrows; n++) {
      if(application.verbose_state.verbose) printf("Adjusting subint %ld\n", n+1);
      if(fits_read_col(fptr, TFLOAT, colnum_s, 1+n, 1, nchan*(1+fullscales*(npol-1)), NULL, data, &anynul, &status)) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read data\n");
	fits_report_error(stderr, status); 
	return 0;
      }
      for(f = 0; f < nchan; f++) {
	scales[0*nchan+f] *= fac1*fac1;
	scales[1*nchan+f] *= fac2*fac2;
	scales[2*nchan+f] *= fac1*fac2;
	scales[3*nchan+f] *= fac1*fac2;
      }
      if(fits_write_col(fptr, TFLOAT, colnum_s, 1+n, 1, nchan*(1+fullscales*(npol-1)), data, &status) != 0) {
	printerror(application.verbose_state.debug, "ERROR fiddleFITS: Error writing data.\n");
	fits_report_error(stderr, status); 
	return 0;
      }
      if (status) {
	fits_report_error(stderr, status); 
	return 0;
      }
    }
    if(application.verbose_state.verbose) printf("Adjusting done\n");


    free(data);
  }else if(whattodo == 7) {  /* high-pass filter data */
    int colnum_s, colnum_o, colnum_w, colnum_d;
    if(outputfptr != stdout) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: You cannot use -output together with -highpass\n");
      return 0;
    }
    /* Read in header parameters from output file */
    if(fits_movnam_hdu(fptr, BINARY_TBL, "SUBINT", 0, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: SUBINT table does not exist!\n");
      return 0;
    } 
    if (fits_read_card(fptr,"NPOL", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NPOL keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &npol);
    if (fits_read_card(fptr,"NBIN", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NBIN keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &nbin);
    if (fits_read_card(fptr,"NCHAN", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NCHAN keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%d", &nchan);
    if (fits_read_card(fptr,"NAXIS2", card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: NAXIS2 keyword does not exist\n");
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    sscanf(value, "%ld", &nrows);
    if(fits_get_colnum (fptr, CASEINSEN, "DAT_SCL", &colnum_s, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No scales in output fits file?\n");
      return 0;
    }
    if(fits_get_colnum (fptr, CASEINSEN, "DAT_OFFS", &colnum_o, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No offsets in output fits file?\n");
      return 0;
    }
    if(fits_get_colnum (fptr, CASEINSEN, "DAT_WTS", &colnum_w, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No weights in output fits file?\n");
      return 0;
    }
    if(fits_get_colnum (fptr, CASEINSEN, "DATA", &colnum_d, &status)) { 
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: No data in output fits file?\n");
      return 0;
    }
    sprintf(txt, "TFORM%d", colnum_s);
    if (fits_read_card(fptr,txt, card, &status)) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: keyword %s does not exist\n", txt);
      return 0;
    }
    fits_parse_value(card, value, comment, &status);
    if(value[1] == 'E') {   /* Simply 1 float instead of an array which is written as '(number)E' */
      datalength = 1;
    }else {
      sscanf(value, "'%ld'", &datalength);
    }
    if(application.verbose_state.verbose)
      printf("Output datalength = %ld\n", datalength);
    if(datalength == npol*nchan) {
      fullscales = 1;
    }else if(datalength == nchan) {
      fullscales = 0;
    }else {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: %ld is a weird number for the number of scales.\n", datalength);
      return 0;
    }

    if(application.verbose_state.verbose) printf("Output data (%s) contains %ld subints, %d freq channels, %d pol channels and %d phase bins (fullscales=%d)\n", input_name, nrows, nchan, npol, nbin, fullscales);
  
    int *data_int;
    float *data_float, scales[4], offsets[4];
    long freqchan, polarization, istart;
    data_int   = malloc(npol*nbin*sizeof(int));
    data_float = malloc(npol*nbin*sizeof(float));

    if(data_int == NULL || data_float == NULL) {
      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Memory allocation failed\n");
      return 0;
    }

    if(application.verbose_state.verbose) printf("Start applying high-pass filter\n");
    for(n = 0; n < nrows; n++) {
      if(application.verbose_state.verbose) printf("Fiddling subint %ld\n", n+1);
      for(freqchan = 0; freqchan < nchan; freqchan++) {
	for(polarization = 0; polarization < npol; polarization++) {
	  /* First sample to be read in data */
	  istart = polarization*nbin*nchan+freqchan*nbin;
	  if(fits_read_col(fptr, TINT, colnum_d, 1+n, 1+istart, nbin, NULL, &data_int[polarization*nbin], &anynul, &status)) {
	    printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read data\n");
	    fits_report_error(stderr, status); 
	    return 0;
	  }
	  // Scales only written for each polarization if fullscales is set
	  if(polarization == 0 || fullscales) {
	    if(fits_read_col(fptr, TFLOAT, colnum_s, 1+n, 1+fullscales*polarization*nchan+freqchan, 1, NULL, &scales[polarization], &anynul, &status)) {
	      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read scales\n");
	      fits_report_error(stderr, status); 
	      return 0;
	    }else {
	      scales[polarization] = scales[0];
	    }
	  }
	  // Offsets only written for each polarization if fullscales is set
	  if(polarization == 0 || fullscales) {
	    if(fits_read_col(fptr, TFLOAT, colnum_o, 1+n, 1+fullscales*polarization*nchan+freqchan, 1, NULL, &offsets[polarization], &anynul, &status)) {
	      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot read scales\n");
	      fits_report_error(stderr, status); 
	      return 0;
	    }else {
	      offsets[polarization] = offsets[0];
	    }
	  }
	  // Construct floatint point numbers
	  for(i = 0; i < nbin; i++) {
	    //      int data[i] = (float pulse[i]-offset)/scale;
	    //      int data[i]*scale + offset = float pulse[i];
	    data_float[i+polarization*nbin] = data_int[i+polarization*nbin]*scales[polarization] + offsets[polarization];
	  }
	}

	// Do processing here for each polarization
	for(polarization = 0; polarization < npol; polarization++) {
	  if(highpass_nr < 0) {  // Use Andrews version of the high-pass filter
	    float as, ac, sum, dphi;
	    float as2, ac2, as3, ac3;
	    int nrignoredbins;
	    as = ac = sum = 0;
	    as2 = ac2 = as3 = ac3 = 0;
	    dphi=2.0*M_PI/(float)nbin;
	    nrignoredbins = 0;
	    for(i = 0; i < nbin; i++) {
	      //	      if(i >= 120 && i <= 140) {
	      //		nrignoredbins++;
	      //	      }else {
		as += data_float[i+polarization*nbin]*sin((i+0.5)*dphi);
		ac += data_float[i+polarization*nbin]*cos((i+0.5)*dphi);
		as2 += data_float[i+polarization*nbin]*sin(2.0*(i+0.5)*dphi);
		ac2 += data_float[i+polarization*nbin]*cos(2.0*(i+0.5)*dphi);
		as3 += data_float[i+polarization*nbin]*sin(3.0*(i+0.5)*dphi);
		ac3 += data_float[i+polarization*nbin]*cos(3.0*(i+0.5)*dphi);
		sum += data_float[i+polarization*nbin];
		//	      }
	    }
	    for(i = 0; i < nbin; i++) {
	      /*
	      data_float[i+polarization*nbin] -= (2.0/(float)nbin)*(as*sin((i+0.5)*dphi) + ac*cos((i+0.5)*dphi)) + sum/(float)nbin;
	      data_float[i+polarization*nbin] -= (2.0/(float)nbin)*(as2*sin(2.0*(i+0.5)*dphi) + ac2*cos(2.0*(i+0.5)*dphi));
	      data_float[i+polarization*nbin] -= (2.0/(float)nbin)*(as3*sin(3.0*(i+0.5)*dphi) + ac3*cos(3.0*(i+0.5)*dphi));
	      */
	      nrignoredbins = 0;
	      data_float[i+polarization*nbin] -= (2.0/(float)(nbin-nrignoredbins))*(as*sin((i+0.5)*dphi) + ac*cos((i+0.5)*dphi)) + sum/(float)(nbin-nrignoredbins);
	      //	      data_float[i+polarization*nbin] -= (2.0/(float)(nbin-nrignoredbins))*(as2*sin(2.0*(i+0.5)*dphi) + ac2*cos(2.0*(i+0.5)*dphi));
	      //	      data_float[i+polarization*nbin] -= (2.0/(float)(nbin-nrignoredbins))*(as3*sin(3.0*(i+0.5)*dphi) + ac3*cos(3.0*(i+0.5)*dphi));
	    }
	   }else {
	    if(fftSmooth(&data_float[polarization*nbin], nbin, 1, highpass_nr, application.verbose_state) != 1) {
	      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Applying high-pass filter failed\n");
	      return 0;
	    }
	    // Normalise output
	    float fac;
	    fac = 1.0/(float)nbin;
	    for(i = 0; i < nbin; i++) {
	      data_float[i+polarization*nbin] *= fac;
	    }
	  }
	}

	// Determine scaling
	if(fullscales == 0) {  // One scale for all polarizations
	  internalFITSscalePulse(data_float, npol*nbin, &offsets[0], &scales[0], 32767);
	  if(fits_write_col(fptr, TFLOAT, colnum_s, 1+n, 1+fullscales*polarization*nchan+freqchan, 1, &scales[0], &status) != 0) {
	    printerror(application.verbose_state.debug, "ERROR fiddleFITS: Error writing scales.\n");
	    fits_report_error(stderr, status); 
	    return 0;
	  }
	  if(fits_write_col(fptr, TFLOAT, colnum_o, 1+n, 1+fullscales*polarization*nchan+freqchan, 1, &offsets[0], &status) != 0) {
	    printerror(application.verbose_state.debug, "ERROR fiddleFITS: Error writing offsets.\n");
	    fits_report_error(stderr, status); 
	    return 0;
	  }
	  for(polarization = 0; polarization < npol; polarization++) {
	    for(i = 0; i < nbin; i++) {
	      //      int data[i] = (float pulse[i]-offset)/scale;
	      //      int data[i]*scale + offset = float pulse[i];
	      data_int[i+polarization*nbin] = (data_float[i+polarization*nbin]-offsets[0])/scales[0];
	    }
	  }
	}else {
	  for(polarization = 0; polarization < npol; polarization++) {
	    internalFITSscalePulse(&data_float[nbin*polarization], nbin, &offsets[polarization], &scales[polarization], 32767);
	    if(fits_write_col(fptr, TFLOAT, colnum_s, 1+n, 1+fullscales*polarization*nchan+freqchan, 1, &scales[polarization], &status) != 0) {
	      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Error writing scales.\n");
	      fits_report_error(stderr, status); 
	      return 0;
	    }
	    if(fits_write_col(fptr, TFLOAT, colnum_o, 1+n, 1+fullscales*polarization*nchan+freqchan, 1, &offsets[polarization], &status) != 0) {
	      printerror(application.verbose_state.debug, "ERROR fiddleFITS: Error writing scales.\n");
	      fits_report_error(stderr, status); 
	      return 0;
	    }
	    for(i = 0; i < nbin; i++) {
	      //      int data[i] = (float pulse[i]-offset)/scale;
	      //      int data[i]*scale + offset = float pulse[i];
	      data_int[i+polarization*nbin] = (data_float[i+polarization*nbin]-offsets[polarization])/scales[polarization];
	    }
	  }
	}
	// Now write out data
	for(polarization = 0; polarization < npol; polarization++) {
	  /* First sample to be read in data */
	  istart = polarization*nbin*nchan+freqchan*nbin;
	  if(fits_write_col(fptr, TINT, colnum_d, 1+n, 1+istart, nbin, &data_int[polarization*nbin], &status)) {
	    printerror(application.verbose_state.debug, "ERROR fiddleFITS: Cannot write data\n");
	    fits_report_error(stderr, status); 
	    return 0;
	  }
	}

	if (status) {
	  fits_report_error(stderr, status); 
	  return 0;
	}
      }  // End frequency channel loop
    } // End subint loop
    if(application.verbose_state.verbose) printf("Fiddling done\n");


    free(data_int);  
    free(data_float);  
  }else {
    printerror(application.verbose_state.debug, "ERROR fiddleFITS: Internal logic error\n");
    return 0;
  }

  /* Close files */
  fits_close_file(fptr, &status);
  if(status) 
    fits_report_error(stderr, status); 
  if(fitsfile2_name > 0) {
    fits_close_file(fptr2, &status);
    if(status) 
      fits_report_error(stderr, status); 
  }
  if(outputfptr != stdout) {
    fclose(outputfptr);
  }
  terminateApplication(&application);
  return 0;
}

