#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "sla_wrap.h"
#include "psrsalsa.h"
#include "psrsalsa_tempo3.h"

//Secret options: 
//
//-labdem_verbose  Enable verbose options to diagnose problems etc
//
//-labdem_test timfile barytable   Test out barycentring in the way the students do this

#define tmp_parfile "/tmp/42ft_tmp_parfile.par"
#define tmp_timfile "/tmp/42ft_tmp_timfile.tim"

int get_obsinfo_from_file(char *filename_ptr, double *period_observation, long double *startmjd_observation, psrsalsaApplication application);
int get_ephemeris_from_file(char *filename_ptr, char *filename_ptr_eph, ephemeris_def *parfile, verbose_definition verbose);
int make_bary_table(int year, ephemeris_def parfile, verbose_definition verbose);

// Not used function, but needs to be defined to use library
void print_ephemeris_gtk()
{
}


int main(int argc, char **argv)
{
  int i, labdem_test, index, guessperiod, barytable;
  char labdem_timfile[MaxFilenameLength], labdem_baryfile[MaxFilenameLength];
  psrsalsaApplication application;

  initApplication(&application, "42ft", "[options] inputfile(s)");
  initialise_tempo3_lib(application.verbose_state);


  // application.switch_verbose = 1;
  //  application.switch_debug = 1;
  guessperiod = 0;
  barytable = 0;
  labdem_test = 0;

  if(argc < 2) {
    printf("Program to aid the 42-ft lab experiment in Jodrell Bank. Run it with your data file from your observation as the last command-line option. Use one of the following additional command-line options to create input to help you with doing pulsar timing:\n\n");
    printApplicationHelp(&application);
    printf("For [options] choose one of the following:\n\n");
    printf("-guessperiod              For a given filename, produce an estimation for the pulsar period.\n");
    printf("-barytable                For a given filename, produce a table with positions of the barycentre.\n");
  }else {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-guessperiod") == 0) {
	guessperiod = 1;
      }else if(strcmp(argv[i], "-barytable") == 0) {
	barytable = 1;
      }else if(strcasecmp(argv[i], "-labdem_verbose") == 0) {
	application.verbose_state.verbose = 1;
	application.verbose_state.debug = 1;
      }else if(strcasecmp(argv[i], "-labdem_test") == 0 || strcasecmp(argv[i], "-labdem_check") == 0 || strcasecmp(argv[i], "-labdemtest") == 0 || strcasecmp(argv[i], "-labdemcheck") == 0) {
	labdem_test = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%1000s %1000s", &labdem_timfile, &labdem_baryfile, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR 42ft: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else {
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "42ft: Unknown option: %s\n\nRun 42ft without command line arguments to show help", argv[i]);
	  terminateApplication(&application);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }

  if(labdem_test == 0) {
    char *filename_ptr;
    while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
      double period_observation;
      long double startmjd_observation;
      if(get_obsinfo_from_file(filename_ptr, &period_observation, &startmjd_observation, application) == 0) {
	printerror(application.verbose_state.debug, "ERROR 42ft: Getting memory from file failed.");
	terminateApplication(&application);
	return 0;
      }
      if(application.verbose_state.verbose) {
	printf("\nThe data file reports a pulse period of %.10lf seconds\nStart MJD is %Lf\n\n", period_observation, startmjd_observation);
      }
      ephemeris_def parfile;
      if(get_ephemeris_from_file(filename_ptr, tmp_parfile, &parfile, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR 42ft: Getting ephemeris from file failed.");
	terminateApplication(&application);
	return 0;
      }
      
      
      // Make a single TOA to be sent to tempo2 to get SSB info etc
      toa_data_def jbo_toa;
      if(initialise_toas(&jbo_toa, 1, 0, 0, 0, 0, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR 42ft: Cannot initialise TOAs");
	return 0;
      }
      if(add_toa(&jbo_toa, 0, startmjd_observation, 123.0, 1400, 1400, 0, "FILENAME", NULL, "jb42", NULL, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR 42ft: Adding TOAs failed");
	return 0;
      }
      if(writetimfile(tmp_timfile, jbo_toa, NULL, NULL, 0, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR 42ft: Writing timfile failed");
	return 0;
      }
      free_toas(&jbo_toa);
      unsigned long *indx;
      if(loadtimfile(tmp_timfile, 0, 0, tmp_parfile, parfile, &jbo_toa, 1, &indx, 0, NULL, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR 42ft: Loading timfile failed");
	return 0;
      }
      free(indx); // Only 1 TOA, so sorting info not relevant
      if(application.verbose_state.verbose) {
	printf("\nBarycentric start MJD of observation is: %Lf\n", jbo_toa.mjd[0]);
	printf("The barycentric position is: %Lf %Lf %Lf\n", jbo_toa.ssb1[0]*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, jbo_toa.ssb2[0]*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, jbo_toa.ssb3[0]*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU);
      }
      

      long double spinfreq;
      if(evaluate_ephemeris_spinfrequency(parfile, startmjd_observation, &spinfreq, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR 42ft: Deriving spin-frequency from ephemeris failed.");
	return 0;
      }
      long double spinfreq_correct;
      if(evaluate_ephemeris_spinfrequency(parfile, jbo_toa.mjd[0], &spinfreq_correct, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR 42ft: Deriving spin-frequency from ephemeris failed.");
	return 0;
      }
      if(guessperiod) {
	if(application.verbose_state.verbose) {
	  printf("\nPeriod in datafile is                       %.10lf seconds\n", period_observation);
	  printf(  "Predicted period at start of observation is %.10Lf seconds (not barycentre corrected)\n", 1.0/spinfreq);
	  printf(  "Predicted period at start of observation is %.10Lf seconds (barycentre corrected)\n", 1.0/spinfreq_correct);
	}
	long double timeperturn = 2.0*3600.0;
	printf(  "Guessed initial period %.10Lf seconds\n", 1.0/(spinfreq_correct+1.0/timeperturn));
      }
      free_toas(&jbo_toa);

      if(barytable) {
	int year, month, day, hour, minute;
	float seconds;
	mjd2date(startmjd_observation, &year, &month, &day, &hour, &minute, &seconds);
	if(application.verbose_state.verbose) {
	  printf("MJD %Lf corresponds to a year %d\n", startmjd_observation, year);
	}
	//	year = 2007;
	if(make_bary_table(year, parfile, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR 42ft: Generating barycentre table failed");
	  return 0;
	}
      }

      if(free_ephemeris(&parfile, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR 42ft: Releasing ephemeris memory failed");
	return 0;
      }
      remove(tmp_parfile);
      remove(tmp_timfile);
    }
  }

  if(labdem_test == 1) {
    long nrTOAs, nrDaysTable;
    double *toas, *barytableX, *barytableY, *barytableZ;
    if(read_ascii_column_double(labdem_timfile, 0, 'C', 5, 0, &nrTOAs, 3, 1.0, 0, &toas, NULL, NULL, NULL, application.verbose_state, 1) == 0) {
      printerror(application.verbose_state.debug, "ERROR 42ft: Reading tim file failed.");
      return 0;
    }
    fprintf(stderr, "Read in %ld toas\n", nrTOAs);
    if(read_ascii_column_double(labdem_baryfile, 0, 0, 6, 0, &nrDaysTable, 4, 1.0, 0, &barytableX, NULL, NULL, NULL, application.verbose_state, 1) == 0) {
      printerror(application.verbose_state.debug, "ERROR 42ft: Reading tim file failed.");
      return 0;
    }
    if(read_ascii_column_double(labdem_baryfile, 0, 0, 6, 0, &nrDaysTable, 5, 1.0, 0, &barytableY, NULL, NULL, NULL, application.verbose_state, 1) == 0) {
      printerror(application.verbose_state.debug, "ERROR 42ft: Reading tim file failed.");
      return 0;
    }
    if(read_ascii_column_double(labdem_baryfile, 0, 0, 6, 0, &nrDaysTable, 6, 1.0, 0, &barytableZ, NULL, NULL, NULL, application.verbose_state, 1) == 0) {
      printerror(application.verbose_state.debug, "ERROR 42ft: Reading tim file failed.");
      return 0;
    }
    fprintf(stderr, "Read in %ld barycentre positions\n", nrDaysTable);
    free(toas);
  }

  cleanup_tempo3_lib(application.verbose_state);
  terminateApplication(&application);
  return 0;
}



// Return 0: Error
// Return 1: Ok
int get_obsinfo_from_file(char *filename_ptr, double *period_observation, long double *startmjd_observation, psrsalsaApplication application)
{
  datafile_definition datain;
  int nowarnings;
  if(application.verbose_state.verbose == 0) {
    nowarnings = 2;
  }else {
    nowarnings = 0;
  }
  if(openPSRData(&datain, filename_ptr, application.iformat, 0, 0, nowarnings, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR 42ft: Error opening data");
    return 0;
  }
  if(readHeaderPSRData(&datain, 1, nowarnings, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR 42ft: Error reading header");
    return 0;
  }
  if(get_period(datain, 0, period_observation, application.verbose_state) != 0) {
    printerror(application.verbose_state.debug, "ERROR 42ft: Fold period undefined?");
    return 0;
  }
  *startmjd_observation = get_mjd_subint(datain, 0, application.verbose_state);

  if(closePSRData(&datain, 0, application.verbose_state) != 0) {
    printerror(application.verbose_state.debug, "ERROR 42ft: Problem with releasing memory related to the pulsar data");
    return 0;
  }
  return 1;
}

int get_ephemeris_from_file(char *filename_ptr, char *filename_ptr_eph, ephemeris_def *parfile, verbose_definition verbose)
{
  char cmd[10000];
  // Get ephemeris from file
  remove(filename_ptr_eph);
  //  sprintf(cmd, "vap -E %s > %s", filename_ptr, filename_ptr_eph);
  sprintf(cmd, "vap -E %s | sed '{s/DE200/DE405/g}' > %s", filename_ptr, filename_ptr_eph);
  if(verbose.verbose) {
    printf("Running command: %s\n", cmd);
  }
  system(cmd);
  if(initialise_ephemeris(parfile, 0, verbose) == 0) {
    printerror(verbose.debug, "ERROR 42ft: Initialising parfile failed");
    //    remove(filename_ptr_eph);
    return 0;
    }
  if(read_ephemeris(filename_ptr_eph, parfile, 0, 1, verbose) == 0) {
    printerror(verbose.debug, "ERROR 42ft: Reading parfile failed");
    //    remove(filename_ptr_eph);
    return 0;
  }
  //  remove(filename_ptr_eph);
  if(verbose.verbose) {
    verbose.indent = 2;
    if(print_ephemeris(stdout, parfile, NULL, NULL, 1, 0, 1, 1, 0, 0, verbose) == 0) {
      printerror(verbose.debug, "ERROR 42ft: Showing parfile failed");
      return 0;
    }
  }
  verbose.indent = 0;
  return 1;
}


// Return 0 on error
// make_bary_table(year, application.verbose_state)
int make_bary_table(int year, ephemeris_def parfile, verbose_definition verbose)
{
  int ret, i;
  double startmjd_dbl, endmjd_dbl;
  int month, day, hour, minute;
  float seconds;
  unsigned long *indx;
  csla_cldj(year, 1, 1, &startmjd_dbl, &ret);
  if(ret != 0) {
    fprintf(stderr, "ERROR make_bary_table: Cannot calculate mjd.\n");
    return 0;
  }
  csla_cldj(year+1, 1, 1, &endmjd_dbl, &ret);
  if(ret != 0) {
    fprintf(stderr, "ERROR make_bary_table: Cannot calculate mjd.\n");
    return 0;
  }
  endmjd_dbl -= 1.0;
  if(verbose.verbose) {
    printf("Going to create table between MJD %lf and %lf\n", startmjd_dbl, endmjd_dbl);
  }
  long startmjd_int, endmjd_int;
  startmjd_int = round(startmjd_dbl);
  endmjd_int = round(endmjd_dbl);
	
      
  toa_data_def bary_toa;
  if(initialise_toas(&bary_toa, endmjd_int-startmjd_int+1, 0, 0, 0, 0, verbose) == 0) {
    printerror(verbose.debug, "ERROR make_bary_table: Cannot initialise TOAs");
    return 0;
  }
  for(i = startmjd_int; i <= endmjd_int; i++) {
    if(add_toa(&bary_toa, i-startmjd_int, i, 123.0, 1400, 1400, 0, "FILENAME", NULL, "@", NULL, verbose) == 0) {
      printerror(verbose.debug, "ERROR make_bary_table: Adding TOAs failed");
      return 0;
    }
  }
  if(writetimfile(tmp_timfile, bary_toa, NULL, NULL, 1, verbose) == 0) {
    printerror(verbose.debug, "ERROR make_bary_table: Writing timfile failed");
    return 0;
  }
  free_toas(&bary_toa);
  if(loadtimfile(tmp_timfile, 0, 0, tmp_parfile, parfile, &bary_toa, 1, &indx, 0, NULL, verbose) == 0) {
    printerror(verbose.debug, "ERROR make_bary_table: Loading timfile failed");
    return 0;
  }
  free(indx); // TOAs, already sorted, so sorting info not relevant
  
  for(i = startmjd_int; i <= endmjd_int; i++) {
    mjd2date(i, &year, &month, &day, &hour, &minute, &seconds);
    printf("%d %d %d %.9Lf %.9Lf %.9Lf\n", year, month, day, bary_toa.ssb1[i-startmjd_int]*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, bary_toa.ssb2[i-startmjd_int]*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU, bary_toa.ssb3[i-startmjd_int]*TEMPO3_SPEEDOFLIGHT/TEMPO3_AU);
  }
	
  free_toas(&bary_toa);
  return 1;
}
