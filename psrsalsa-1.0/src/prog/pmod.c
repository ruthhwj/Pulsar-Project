//START REGION RELEASE
#define _FILE_OFFSET_BITS 64
#define _USE_LARGEFILE 1
#define _USE_LARGEFILE 1
#define _LARGEFILE_SOURCE 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "psrsalsa.h"

#include "gsl/gsl_randist.h"

//START REGION DEVELOP
#include "gsl/gsl_rng.h"
#include <gsl/gsl_sort_float.h>
#include <gsl/gsl_statistics_float.h>
#include <gsl/gsl_fit.h>
/*#include <curses.h>*/

void SHOWREVISIONINFO_prog() {
#include "pmod.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);

//START REGION RELEASE
int readZapFile(char *zaplistname, int *zapMask, int zapSkipLines, int nrZapCols, int zapColumn, int zapColumn2, int inverseZap, verbose_definition verbose);
void make_blocks(long baseline_length, long blockSize, long nrPulses, long *nrOutputBlocks, int *zapMask, verbose_definition verbose);
//START REGION DEVELOP
int zapRMS(long nrPulses, float *rms, unsigned long *indx, int *zapMask, float rms_thresh, int subtractMinimum, int subtractMedian, float fracData, char *plotDevice, int *deviceOpened, verbose_definition verbose);


//START REGION RELEASE

int main(int argc, char **argv)
{
  int debase_flag, debase_offset_flag, index, deviceOpened, read_whole_file;
  int zapoption, inverseZap, fzapoption, finverseZap, zapColumn, zapColumn2, nrZapCols, zapSkipLines;
  int blockMode, remove_pulses_flag, prange_set;
  int nrPol, nrBins, NrFreqChan, addnoise_flag, removeOnPulse_flag;
  int selectMoreOnpulseRegions;
  int outputlist, filename;
  int *zapMask, *fzapMask;
  long i, j, k, l, m, tmp_used_pulses;
  long nrOutputBlocks, nrZapped, baseline_length, blockSize, firstPulseToKeep, lastPulseToKeep, nrPulses;
  long dataout2_pulse, pulse_nr_in_output, idnum;
  float *profileI, debase_offset_value, *baseline, *runningBaseline, *runningRMS, *rms, noiseRMS;
  char *filename_ptr, zaplistname[MaxFilenameLength], fzaplistname[MaxFilenameLength];
  char output_suffix[MaxFilenameLength], output_suffix2[MaxFilenameLength], output_name[MaxFilenameLength], output_name2[MaxFilenameLength];
  char txt[1000];
  datafile_definition datain, dataout, dataout2;
  psrsalsaApplication application;
  gsl_rng *rand_num_gen;
  const gsl_rng_type *rand_num_gen_type;
  pgplot_options_definition pgplot_options;
//START REGION DEVELOP
  int debase_slope_flag, debase_highpass_flag, debase_highpass_nrharmonics, clean_flag, poly, poly_replace;
  int correctbaseline, correctbaseline2;
  int cutoption, cutpower2, clippos_flag, clipneg_flag, threshold_flag, threshold_excess_flag;
  int status;
  long j2;
  float correctbaseline_bins, correctbaseline_amplitude, correctbaseline2_bins, correctbaseline2_amplitude, scatter_timescale, scatter_reffreq, scatter_pwrlaw, removeOffPulse_rms; //, clean_sigma;
  float rms_thresh1, rms_thresh2, rms_thresh3, I, clippos_val, clipneg_val, threshold_val, debase_sigma;
  double *profileI2_double, *profilex_double, *profileSigma_double;
  pulselongitude_regions_definition cutSelection;
  unsigned long *indx;



//START REGION RELEASE
  initApplication(&application, "pmod", "[options] inputfile(s)");
//START REGION DEVELOP
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
//START REGION RELEASE
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.oformat = FITS_format;
//START REGION DEVELOP
//  application.oformat = PUMA_format;
//START REGION RELEASE
  application.switch_formatlist = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_onpulse = 1; 
  application.switch_onpulsef = 1;
  application.switch_filelist = 1;
  application.switch_device = 1;
  application.switch_conshift= 1;
  application.switch_circshift= 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_rebin = 1;
  application.switch_tscr = 1;
  application.switch_TSCR = 1;
  application.switch_tscr_complete = 1;
  application.switch_fscr = 1;
  application.switch_FSCR = 1;
  application.switch_nocounters = 1;
  application.switch_dedisperse = 1;
  application.switch_deFaraday = 1;
  application.switch_polselect = 1;
  application.switch_stokes = 1;
  application.switch_coherence = 1;
  application.switch_noweights = 1;
  application.switch_useweights = 1;
  application.switch_uniformweights = 1;
  application.switch_scale = 1;
  application.switch_insertparang = 1;
  application.switch_deparang = 1;
  application.switch_changeRefFreq = 1;
  application.switch_history_cmd_only = 1;
  application.switch_fixseed = 1;
  application.switch_templatedata = 1;
  application.switch_template = 1;
  application.switch_align = 1;
  application.switch_shuffle = 1;
  application.switch_rotateStokes = 1;
  application.switch_libversions = 1;
//START REGION DEVELOP
  application.switch_iformat = 1;    /* Default output format is PuMa */  
  application.switch_oformat = 1;
  application.switch_invertxy = 1;
  application.switch_fftsmooth = 1;
  application.switch_fftzap = 1;
  application.switch_fft = 1;
  application.switch_template_onpulse = 1;
  application.switch_centert = 1;
  application.switch_autot= 1;
  application.switch_norm = 1;
  application.switch_normglobal = 1;
  application.switch_normRMS = 1;
  application.switch_normRMSi = 1;
  application.switch_invertfp = 1;
  application.switch_invertfx = 1;
  // Doesn't work with multiple files, probably because onpulse get changed. However, this is identical to the -l option which was already implemented, which is now renamed -gate.
  //  application.switch_gate = 1;
  application.switch_absweights = 1;
  application.switch_testsignal = 1;
  application.switch_checknan = 1;
  application.switch_removenan = 1;
  application.switch_checkinf = 1;
  application.switch_removeinf = 1;
  application.switch_rotateQU = 1;     // Use rotateStokes instead in release
  application.switch_rotateUV = 1;     // Use rotateStokes instead in release
  application.switch_showrevision = 1;
  application.switch_dedeFaraday = 1;
  application.switch_dededisperse= 1;
  application.switch_swapcables = 1;
  application.switch_smooth = 1;
  application.switch_autoonpulse = 1;
  application.switch_fchan = 1;
  application.switch_rotslope = 1;
  application.switch_subtractRVM = 1;
  application.switch_subtractAvPA = 1;
  application.switch_zerodming = 1;
  application.switch_zerodming_adv = 1;
  application.switch_nskip = 1;   // Don't put in release version, as it is essentially the same as the -prange option
  application.switch_nread = 1;   // Don't put in release version, as it is essentially the same as the -prange option
  application.switch_deFaradayTable = 1;

//START REGION DEVELOP
//START REGION RELEASE
  application.switch_forceUniformFreqLabelling = 1;

  /* Initialize standard values */
  debase_flag = 0;
  debase_offset_flag = 0;
  read_whole_file = 1;
  zapoption = 0;
  inverseZap = -1;  // -1 is undefined, 1=yes, 0=no
  finverseZap = -1;
  fzapoption = 0;
  nrZapCols = 1;
  zapColumn = 1;
  zapColumn2 = 0;  // Default - only read in one file
  zapSkipLines = 0;
  blockMode = 0;
  remove_pulses_flag = 0;
  prange_set = 0;
  debase_offset_value = 0;
  baseline_length = 0;
  nrOutputBlocks = 0;
  outputlist = 0;
  addnoise_flag = 0;
  noiseRMS = 0;
  selectMoreOnpulseRegions = 0;
  removeOnPulse_flag = 0;
  strcpy(output_suffix, "debase.gg");
  strcpy(output_suffix2, "zapped.gg");
  filename = 0;
  pgplot_clear_options(&pgplot_options);

//START REGION DEVELOP
  debase_sigma = -1;
  correctbaseline = 0;
  correctbaseline2 = 0;
  poly = 0;
  poly_replace = 0;
  clean_flag = 0;
  if(initPulselongitudeRegion(&cutSelection, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pmod: Cannot initialise onpulse region.");
    return 0;
  }
  cutSelection.nrRegions = 0;
  cutoption = 0;
  cutpower2 = 0;
  threshold_flag = 0;
  clippos_flag = 0;
  clipneg_flag = 0;
  threshold_excess_flag = 0;
  rms_thresh1 = 10;
  rms_thresh2 = 4;
  rms_thresh3 = 4;
  debase_slope_flag = 0;
  debase_highpass_flag = 0;
  scatter_timescale = -1;
  removeOffPulse_rms = -1;

//START REGION RELEASE
  if(argc < 2) {
    printf("Program to modify pulsar data in various ways. Usage:\n\n");
    printApplicationHelp(&application);
    printf("Action options:\n\n");
    printf("-debase              Subtract baseline from observation\n");
    printf("-debase_length       Specify number of pulses before and after the pulse to\n");
    printf("                     determine running baseline (default is -debase_length 0).\n");
    printf("                     Note that the number of pulses used is 2N+1\n");
//START REGION DEVELOP
    printf("-debase_sigma        Specify that pulses in excess of this number of sigma above\n");
    printf("                     the noise will be ignored in the running mean calculation.\n");
    printf("                     the rms of the noise is calculated over the debase_length.\n");
    printf("-debase_slope        Subtract linear baseline from each freqency channel/subint\n");
    printf("                     separately using the off-pulse region. This is done with\n");
    printf("                     linear regression, unless -debase_highpass is used, which\n");
    printf("                     implies Levenberg-Marquardt fitting is used.\n");
    printf("-debase_highpass n   Subtract n lower harmonics from the baseline of each freq\n");
    printf("                     channel/subint separately (so up to and including a\n");
    printf("                     wavelength nbin/n, where nbin is the number of longitude\n");
    printf("                     bins per subint. Only the off-pulse region is fitted\n");
    printf("                     (with Levenberg-Marquardt algorithm).\n");
//START REGION RELEASE
    printf("-debase_value        Subtract the specified value baseline (fixed value) rather\n");
    printf("                     than the determined baseline, i.e. -debase_value 1.234\n");
//START REGION DEVELOP
    printf("-clean               Try to clean observation automaticly\n");
    //    printf("-cleanradar s        Try to clean observation automaticly for undispersed negative dips\n");
    //    printf("                     exceeding s times the rms\n");
    printf("-correctbaseline     Try to correct negative dips after strong pulses by\n");
    printf("                     convolving the signal with an exponential recovery using\n");
    printf("                     this timescale (in bins) and fractional amplitude. The\n");
    printf("                     data should be debased first.\n");
    printf("-correctbaseline2    idem, but now correct dips BEFORE a strong pulse\n");
//START REGION RELEASE
    printf("-list                List the zap list to a file\n");
//START REGION DEVELOP
    printf("-gate                Cut longitude range to the onpulse region (reduce nr of bins).\n");
    printf("-gate2               Cut longitude range with 2**n samples\n");
    printf("-clip                Specify threshold value above which the value clipped.\n");
    printf("                     Clipping happens if sample > threshold or < -threshold.\n");
    printf("-clippos             Specify threshold value above which the value clipped.\n");
    printf("                     Clipping happens if sample > threshold only.\n");
    printf("-clipneg             Specify threshold value below which the value clipped.\n");
    printf("                     Clipping happens if sample < threshold only.\n");
    printf("-excess              Use together with -clip. It will write out the remainder after\n");
    printf("                     applying the clipping.\n");
    printf("-mask                Specify threshold value to make a mask\n");
    printf("                     If sample > threshold, the value will be threshold.\n");
    printf("                     If sample < -threshold, the value will be -threshold.\n");
    printf("                     If sample == threshold, the value will be 0.\n");
//START REGION RELEASE
    printf("-addnoise rms        Add white noise with RMS rms to the data.\n"); 
    printf("-onpulse_subst_noise Substitute the onpulse region with white noise based on the\n");
    printf("                     (running) off-pulse rms.\n");
//START REGION DEVELOP
    printf("-offpulse_subst_noise  RMS       Substitute the offpulse region with white noise\n");
    printf("                     with the given rms.\n");
    printf("-scatter             \"time reffreq index\" Convolve the data with a scatter tail\n");
    printf("                     with e-fold timescale in sec at the reference frequency in\n");
    printf("                     MHz, with a powerlaw dependence on frequency with index\n");
    printf("                     index (i.e. -4). Each pulse is treated seperately, so the\n");
    printf("                     end of scatter tail ends up at the start of the same\n"); 
    printf("                     pulse/subint.\n"); 
    printf("-poly N              Subtract a polynomial of order N from the off-pulse region\n");
    printf("-polyrep N           Same as -poly, but store the polynomial\n");

//START REGION RELEASE
    printf("\nData selection options:\n\n");
    printf("-zapfile file     Specify filename with pulse numbers to zap (first pulse is 0).\n");
    printf("                  Expected format: see -format.\n");
    printf("-zapfile_i        or alternatively, specify file with pulse numbers NOT to zap\n");
    printf("-zap              or alternatively, specify first and last pulse to zap, \n");
    printf("                  e.g. -zap \"0 10\" (can specify -zap multiple times)\n");
    printf("-prange           or alternatively, specify first and last pulse to keep\n");
    printf("-format           Specify column number, total number columns and lines to skip\n");
    printf("                  in zap file (default is \"%d %d %d\")\n", zapColumn, nrZapCols, zapSkipLines);
    printf("-fzapfile         Specify file with freq channels to zap (first channel is 0)\n");
    printf("-fzapfile_i       or alternatively, specify file with channels NOT to zap\n");
    printf("-zapfile2         If this flag is specified, the -zapfile or -zapfile_i or\n");
    printf("                  equivalent frequency channel zap options should specify a file\n");
    printf("                  with two columns, being the start and end values of ranges\n");
    printf("                  rather than individual values.\n");
    printf("-fzap             Specify first and last frequency channel to zap\n");
    printf("                  (can specify -fzap multiple times)\n");
    printf("-blocksize        Nr of succesive nonzapped pulses that should be writen out\n");
    printf("-remove           Remove zapped pulses instead of making them zero\n");

    printf("\nOn-pulse definition:\n\n");
//START REGION DEVELOP
    /*    printf("-C                Select on-pulse region (can use multiple times)\n"); */
//START REGION RELEASE
    printf("-onpulsegr        Enables selecting more on-pulse regions graphically than\n");
    printf("                  defined by -onpulse\n");
//START REGION RELEASE
    printf("\nOutput options:\n\n");
    printf("-ext              Specify suffix, default is '%s'\n", output_suffix);
    printf("-output filename  Write output to filename rather than changing the extension\n");
    printf("                  to '%s'.\n", output_suffix);
    printf("-memsave          Don't read the file in as a whole at the start of the program\n");
//START REGION DEVELOP
    /*
    printf("-w Output suffix. Default is \"debase.gg\"\n");
    printf("-a Specify number of pulses to read to generate average pulse (default is all)\n");
    printf("-p Output profile and baseline to file\n");
    printf("-q Do not output baseline to screen\n"); */
//START REGION RELEASE
    printf("\n");
    printCitationInfo();
//START REGION DEVELOP
    freePulselongitudeRegion(&cutSelection);
//START REGION RELEASE
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
      }else if(strcmp(argv[i], "-debase") == 0) {
	debase_flag = 1;      
      }else if(strcmp(argv[i], "-list") == 0) {
	outputlist = 1;
//START REGION DEVELOP  
// Something strange seems to happen when zapping pulses on command line, removing them and using this option. The pulses don't appear to be removed anymore.
      }else if(strcmp(argv[i], "-debase_sigma") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &debase_sigma, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-debase_value") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &debase_offset_value, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	debase_offset_flag = 1;
	debase_flag = 1;      
	i++;
      }else if(strcmp(argv[i], "-onpulse_subst_noise") == 0) {
	removeOnPulse_flag = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-offpulse_subst_noise") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &removeOffPulse_rms, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-scatter") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f", &scatter_timescale, &scatter_reffreq, &scatter_pwrlaw, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-debase_slope") == 0) {
	debase_slope_flag = 1;
      }else if(strcmp(argv[i], "-debase_highpass") == 0) {
	debase_highpass_flag = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &debase_highpass_nrharmonics, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-correctbaseline") == 0) {
	correctbaseline = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &correctbaseline_bins, &correctbaseline_amplitude, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-correctbaseline2") == 0) {
	correctbaseline2 = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &correctbaseline2_bins, &correctbaseline2_amplitude, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-debase_length") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &baseline_length, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	/* Make sure pulses at the start and end of the observation are removed when subtracting running baseline */
	if(baseline_length != 0)
	  remove_pulses_flag = 1;
        i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-poly") == 0 || strcmp(argv[i], "-polyrep") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &poly, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	if(strcmp(argv[i], "-polyrep") == 0) {
	  poly_replace = 1;
	}
        i++;
      }else if(strcmp(argv[i], "-clean") == 0) {
	clean_flag = 1;      
	//      }else if(strcmp(argv[i], "-cleanradar") == 0) {
	//	clean_flag = 2;      
	//	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &clean_sigma, NULL) == 0) {
	//	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	//	  return 0;
	//	}
	//	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-zapfile") == 0) {
	strcpy(zaplistname, argv[i+1]);
	if(zapoption) {
	  printerror(application.verbose_state.debug, "pmod: Can only either use -zapfile or -zapfile_i, and they cannot be used more than once.");
	  return 0;
	}
	if(inverseZap == 1) {
	  printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
	  return 0;
	}
	zapoption = 1;
	inverseZap = 0;
        i++;
      }else if(strcmp(argv[i], "-zapfile_i") == 0) {
	strcpy(zaplistname, argv[i+1]);
	zapoption = 1;
	inverseZap = 1;
	if(inverseZap == 0) {
	  printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-zapfile2") == 0) {
	if(zapColumn2 == 0) {
	  zapColumn2 = zapColumn+1;
	}
      }else if(strcmp(argv[i], "-fzapfile") == 0) {
	strcpy(fzaplistname, argv[i+1]);
	fzapoption = 1;
	if(finverseZap == 1) {
	  printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
	  return 0;
	}
	finverseZap = 0;
        i++;
      }else if(strcmp(argv[i], "-fzapfile_i") == 0) {
	strcpy(fzaplistname, argv[i+1]);
	fzapoption = 1;
	if(finverseZap == 0) {
	  printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
	  return 0;
	}
	finverseZap = 1;
        i++;
      }else if(strcmp(argv[i], "-blocksize") == 0) {
	blockMode = 1;
	remove_pulses_flag = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &blockSize, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-format") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d %d", &zapColumn, &nrZapCols, &zapSkipLines, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-prange") == 0) {
	prange_set = 1;
	if(inverseZap == 0) {
	  printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
	  return 0;
	}
	inverseZap = 1;
        i++;
      }else if(strcmp(argv[i], "-zap") == 0) {
	if(inverseZap == 1) {
	  printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
	  return 0;
	}
	/* Ignore for now, will deal with it later. */
	i++;
      }else if(strcmp(argv[i], "-fzap") == 0) {
	if(finverseZap == 1) {
	  printerror(application.verbose_state.debug, "pmod: Mixing of zap and inverse-zap options is not allowed.");
	  return 0;
	}
	/* Ignore for now, will deal with it later. */
	i++;
      }else if(strcmp(argv[i], "-addnoise") == 0) {
	addnoise_flag = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &noiseRMS, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-clip") == 0) {
	threshold_flag = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &threshold_val, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-clippos") == 0) {
	clippos_flag = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &clippos_val, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-clipneg") == 0) {
	clipneg_flag = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &clipneg_val, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-excess") == 0) {
	threshold_excess_flag = 1;
	/*	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &threshold_val, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	*/
      }else if(strcmp(argv[i], "-mask") == 0) {
	threshold_flag = 4;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &threshold_val, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
	/*
      }else if(strcmp(argv[i], "-C") == 0) {
	j = sscanf(argv[i+1], "%d %d", &(onpulseRegions.left[onpulseRegions.nrRegions]), &(onpulseRegions.right[onpulseRegions.nrRegions]));
	if(j != 2) {
	  printerror(application.verbose_state.debug, "pmod: Error parsing -C option");
	  return 0;
	}
	onpulseRegions.nrRegions += 1;
	if(onpulseRegions.nrRegions == maxNrRegions) {
	  printerror(application.verbose_state.debug, "pmod: To many regions selected.");
	  return 0;
	}
        i++;*/
//START REGION RELEASE
      }else if(strcmp(argv[i], "-ext") == 0) {
	strcpy(output_suffix, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-remove") == 0) {
	remove_pulses_flag = 1;
      }else if(strcmp(argv[i], "-onpulsegr") == 0) {
	selectMoreOnpulseRegions = 1; 
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-gate") == 0 || strcmp(argv[i], "-l") == 0) {   // -l was old option name, but make it the same as preprocess option which doesn't work correctly with multiple files as probably onpulse gets modified along the way
	cutoption = 1;
	cutpower2 = 0;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-memsave") == 0) {
	read_whole_file = 0;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-gate2") == 0 || strcmp(argv[i], "-L") == 0) {
	cutoption = 1;
	cutpower2 = 1;
	/*      }else if(i < argc-1) {    //Not sure why this was here, given that multiple input files are supported????
	printerror(application.verbose_state.debug, "pmod: Unknown option (%s).", argv[i]);
	return 0; */
//START REGION RELEASE
      }else {
	/* If the option is not recognized, assume it is a filename */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "pmod: Unknown option: %s\n\nRun pmod without command line arguments to show help", argv[i]);
//START REGION DEVELOP
	  freePulselongitudeRegion(&cutSelection);
//START REGION RELEASE
	  terminateApplication(&application);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }

  if(inverseZap == -1)
    inverseZap = 0;
  if(finverseZap == -1)
    finverseZap = 0;


  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pmod: No files specified");
    return 0;
  }

  strcpy(pgplot_options.viewport.plotDevice, application.pgplotdevice);
  pgplot_options.viewport.dontclose = 1;

  gsl_rng_env_setup();    /* Set the default generators, can be influenced by environment variables */
  rand_num_gen_type = gsl_rng_default;
  rand_num_gen = gsl_rng_alloc(rand_num_gen_type);
  if(application.fixseed)
    idnum = 1;
  else
    randomize_idnum(&idnum);
  //  idnum = 0;
  gsl_rng_set(rand_num_gen, idnum);
  
//START REGION RELEASE
  while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
    deviceOpened = 0;
//START REGION DEVELOP
    /*
      if(debase_flag == 0 && clean_flag == 0 && threshold_flag == 0) {
      printerror(application.verbose_state.debug, "Nothing to do.");
      return 0;
      } */
//START REGION RELEASE
    if(filename == 0) {
      if(change_filename_extension(filename_ptr, output_name, output_suffix, MaxFilenameLength, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pmod: Cannot change extension in output name. The input name is expected to have an extension (even in -output is used).");
	return 0;
      }
    }else {
      strcpy(output_name, argv[filename]);
    }
    if(change_filename_extension(filename_ptr, output_name2, output_suffix2, MaxFilenameLength, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR pmod: Cannot change extension in output name. The input name is expected to have an extension (even in -output is used).");
      return 0;
    }

    //    cleanPSRData(&datain, application.verbose_state);
    cleanPSRData(&dataout, application.verbose_state);
    cleanPSRData(&dataout2, application.verbose_state);
    
    if(application.iformat <= 0) {
      application.iformat = guessPSRData_format(filename_ptr, 0, application.verbose_state);
      if(application.iformat == -2 || application.iformat == -3)
	return 0;
    }
    if(isValidPSRDATA_format(application.iformat) == 0) {
      printerror(application.verbose_state.debug, "ERROR pmod: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", filename_ptr);
//START REGION DEVELOP
      freePulselongitudeRegion(&cutSelection);
//START REGION RELEASE
      closePSRData(&dataout, 0, application.verbose_state);
      closePSRData(&dataout2, 0, application.verbose_state);
      gsl_rng_free(rand_num_gen); 
      terminateApplication(&application);
      return 0;
    }
    i = openPSRData(&datain, filename_ptr, application.iformat, 0, read_whole_file, 0, application.verbose_state);
    if(i == 0) {
      printerror(application.verbose_state.debug, "ERROR pmod: Error opening data");
      return 0;
    }
    
 //START REGION DEVELOP
   /*
      if(invertXY && read_whole_file) {
      if(preprocess_invertXY(datain, &clone, application.verbose) == 0)
      return 0;
      swap_orig_clone(&datain, &clone);
      }else if(invertXY) {
      printerror(application.verbose_state.debug, "pmod: Cannot use the -invertXY option together with the -memsave option,");
      return 0;
      }*/
    
//START REGION RELEASE
    /* Read header */
    if(!read_whole_file) {
      if(readHeaderPSRData(&datain, 0, 0, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "pmod: Error reading header");
	return 0;
      }
    }

    /* Search commandline for header parameters to overwrite header
     parameters. */
    if(PSRDataHeader_parse_commandline(&datain, argc, argv, application.verbose_state) == 0)
      return 0;

    /* Convert the specified onpulse regions which were defined as
       fractions on the commandline. This is now possible because we
       know the number of bins. */
    region_frac_to_int(&(application.onpulse), datain.NrBins, 0);
    
    for(i = 1; i < argc; i++) {
      if(strcmp(argv[i], "-header") == 0) {
	fflush(stdout);
	printwarning(application.verbose_state.debug, "WARNING: If using the -header option, be aware it applied BEFORE the preprocessing.");
	break;
      }
    }
    // This has to be before the mallocs of zapMask and fzapMask
    if(preprocessApplication(&application, &datain) == 0) {
      return 0;
    }

//START REGION DEVELOP

    if(scatter_timescale > 0) {
      if(convolveScatterTail(datain, scatter_timescale, scatter_reffreq, scatter_pwrlaw, application.verbose_state) != 1) {
	printerror(application.verbose_state.debug, "ERROR: Adding scatter tail failed.");
	return 0;
      }
    }

//START REGION RELEASE

    /* Read in zap file(s) */
    zapMask = (int *)malloc(datain.NrSubints*sizeof(int));
    fzapMask = (int *)malloc(datain.NrFreqChan*sizeof(int));
    if(zapMask == NULL || fzapMask == NULL) {
      printerror(application.verbose_state.debug, "pmod: Memory allocation error.");
      return 0;
    }
    for(i = 0; i < datain.NrSubints; i++) 
      zapMask[i] = inverseZap;
    for(i = 0; i < datain.NrFreqChan; i++) 
      fzapMask[i] = finverseZap;
    if(argc > 2) {
      for(i = 1; i < argc-1; i++) {
	if(strcmp(argv[i], "-zap") == 0) {
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld %ld", &k, &l, NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	  for(j = k; j <= l; j++) {
	    if(j >= 0 && j < datain.NrSubints)
	      zapMask[j] = 1;
	  }
	  i++;
	}
      }
      for(i = 1; i < argc-1; i++) {
	if(strcmp(argv[i], "-fzap") == 0) {
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld %ld", &k, &l, NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	  for(j = k; j <= l; j++) {
	    if(j >= 0 && j < datain.NrFreqChan) {
	      fzapMask[j] = 1;
	    }
	  }
	  i++;
	}
      }
    }
    if(zapoption != 0) {
      if(readZapFile(zaplistname, zapMask, zapSkipLines, nrZapCols, zapColumn, zapColumn2, inverseZap, application.verbose_state) == 0)
	return 0;
      if(application.verbose_state.verbose) {
	j = 0;
	for(i = 0; i < datain.NrSubints; i++) {
	  if(zapMask[i] != 0)
	    j++;
	}
	printf("Zapfile zapped %ld pulses.\n", j);
      }
    }
    if(fzapoption != 0) {
      if(readZapFile(fzaplistname, fzapMask, zapSkipLines, nrZapCols, zapColumn, zapColumn2, finverseZap, application.verbose_state) == 0)
	return 0;
      if(application.verbose_state.verbose) {
	j = 0;
	for(i = 0; i < datain.NrFreqChan; i++) {
	  if(fzapMask[i] != 0)
	    j++;
	}
	printf("Zapfile zapped %ld/%ld frequency channels.\n", j, datain.NrFreqChan);
	if(j == 0) /* Make sure the mask is applied whith the freq. scrunch preprocess option */
	  application.fzapMask = NULL;
	else
	  application.fzapMask = fzapMask;
      }
    }
    /* Apply prange options */
    if(prange_set) {
      for(i = 1; i < argc-1; i++) {
	if(strcmp(argv[i], "-prange") == 0) {
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld %ld", &firstPulseToKeep, &lastPulseToKeep, NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pmod: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	  if(lastPulseToKeep >= datain.NrSubints) {
	    printerror(application.verbose_state.debug, "ERROR pmod: Invalid -prange option. Last subint to keep %ld >= %ld", lastPulseToKeep, datain.NrSubints);
	    return 0;
	  }
	  if(firstPulseToKeep < 0) {
	    printerror(application.verbose_state.debug, "ERROR pmod: Invalid -prange option. First subint to keep %ld < 0", firstPulseToKeep);
	    return 0;
	  }
	  for(j = firstPulseToKeep; j <= lastPulseToKeep; j++)
	    zapMask[j] = 0;
	  i++;
	}
      }
    }

//START REGION RELEASE
   

    /*  prheader(&hdr,datain.fptr); */
    /*  change_cmd_line(&hdr, argc, argv); */
    nrPol = datain.NrPols;
    nrBins = datain.NrBins;
    nrPulses= datain.NrSubints;
    NrFreqChan = datain.NrFreqChan;
    
    if(application.verbose_state.verbose) printf("Input data contains %d bins, %ld pulses, %d polarizations and %d frequencies.\n", nrBins, nrPulses, nrPol, NrFreqChan);
    profileI = (float *)malloc(nrBins*sizeof(float));
    rms      = (float *)malloc(nrPulses*nrPol*sizeof(float));
    runningBaseline = (float *)malloc(nrPulses*nrPol*sizeof(float));
    baseline = (float *)malloc(nrPulses*nrPol*sizeof(float));
    runningRMS = (float *)malloc(nrPulses*nrPol*sizeof(float));
//START REGION DEVELOP
    profileI2_double = (double *)malloc(nrBins*sizeof(double));
    profilex_double = (double *)malloc(nrBins*sizeof(double));
    profileSigma_double = (double *)malloc(nrBins*sizeof(double));
    indx     = (unsigned long *)malloc(nrPulses*sizeof(unsigned long));
//START REGION RELEASE
    if(profileI == NULL || rms == NULL || runningBaseline == NULL || baseline == NULL || runningRMS == NULL
//START REGION DEVELOP
       || profileI2_double == NULL || profilex_double == NULL || profileSigma_double == NULL || indx == NULL
//START REGION RELEASE
       ) {
      printerror(application.verbose_state.debug, "ERROR pmod: Memory allocation error.");
      return 0;
    }
    
//START REGION DEVELOP
    
    /*
    for(i = 0; i < nrPulses; i++) 
      printf("Zap mask pulse %ld: %d\n", i, zapMask[i]);
    */

    /* Apply user defined selected pulse range if selected */
    /*
      if(firstPulseToKeep >= 0) {
      for(i = 0; i < nrPulses; i++)
      zapMask[i] = 1;
      for(i = firstPulseToKeep; i <= lastPulseToKeep; i++)
      zapMask[i] = 0;
      blockMode = 1;
      blockSize = lastPulseToKeep-firstPulseToKeep+1;
      }
    */
    
    /* Remove pulses which have an extreme RMS including all bins, this
       in order to get a reasonable profile if there are a few extreme
       large pulses. */
    if(clean_flag == 1) {
      if(application.verbose_state.verbose) {
	printf("Find outliers with extreme large rms'es:\n");
      }
      if(read_rmsPSRData(datain, rms, runningBaseline, zapMask, NULL, 0, 0, -1, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "Error pmod: Cannot determine offpulse rms.");
	return 0;
      }
      zapRMS(nrPulses, rms, indx, zapMask, rms_thresh1, 1, 1, 0.5, application.pgplotdevice, &deviceOpened, application.verbose_state);
    }
    //    if(clean_flag == 2) {
    //      if(preprocess_zero_dming(&datain, -1, clean_sigma, 1, -1, 0, 0, -1, application.verbose_state) == 0) {
    //	printerror(application.verbose_state.debug, "Error pmod: Applying the -cleanradar option failed.");
    //	return 0;
    //      }
    //    }
    
//START REGION DEVELOP
//START REGION RELEASE

    /* Read in profile  */
    if(debase_flag || removeOnPulse_flag ||  selectMoreOnpulseRegions
//START REGION DEVELOP
       || removeOffPulse_rms >= 0 || clean_flag == 1 || poly > 0 || cutoption
//START REGION RELEASE
       ) {
      if(read_profilePSRData(datain, profileI, zapMask, 0, application.verbose_state) != 1) {
	printerror(application.verbose_state.debug, "Error pmod: Cannot form pulse profile");
	return 0;
      }
    }
    
//START REGION DEVELOP

    /* Find on-pulse region with template */
    /*
      if(application.template_specified) {
      blaat
      }
    */
    
    /* let user define pulse profile, if required */
    if(cutoption == 1) {
      if(application.onpulse.nrRegions == 0) {
	pgplot_options.viewport.dontopen = deviceOpened;
	strcpy(pgplot_options.box.xlabel, "Bin");
	strcpy(pgplot_options.box.ylabel, "Intensity");
	strcpy(pgplot_options.box.title, "Select region to cut of ");
	strcat(pgplot_options.box.title, datain.psrname);
	selectRegions(profileI, nrBins, &pgplot_options, 1, cutpower2, 0, &(application.onpulse), application.verbose_state);
	deviceOpened = 1;
      }else {
	if(application.onpulse.right_bin[0] >= nrBins) {
	  application.onpulse.right_bin[0] = nrBins - 1;
	}
	if(application.onpulse.left_bin[0] < 0) {
	  application.onpulse.left_bin[0] = 0;
	}
	if(cutpower2) {
	  if(application.onpulse.bins_defined[0] == 0) {
	    printerror(application.verbose_state.debug, "pmod ERROR: selected region is not defined in bins.");
	    return 0;
	  }
	  I = log(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1.0)/log(2.0);
	  j = I;
	  k = pow(2.0, j+1.0);
	  application.onpulse.left_bin[0] -= 0.5*(k - (application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1.0));
	  if(application.onpulse.left_bin[0] < 0)
	    application.onpulse.left_bin[0] = 0;
	  application.onpulse.right_bin[0] = application.onpulse.right_bin[0]+k-1;
	}
      }
    }
    if(correctbaseline || correctbaseline2) {
      if(debase_flag)
	printerror(application.verbose_state.debug, "pmod ERROR: Baseline will not be subtracted in baseline correction mode.");
      if(clean_flag == 1)
	printerror(application.verbose_state.debug, "pmod ERROR: cleaning will not be done in baseline correction mode.");
      /*    debase_flag = 1; */
    }else {
//START REGION RELEASE
      if((debase_flag && debase_offset_flag == 0) || removeOnPulse_flag || selectMoreOnpulseRegions
//START REGION DEVELOP
	 || removeOffPulse_rms >= 0 || clean_flag == 1 || poly > 0
//START REGION RELEASE
	 ) {
	if(selectMoreOnpulseRegions || application.onpulse.nrRegions == 0) {
	  pgplot_options.viewport.dontopen = deviceOpened;
	  strcpy(pgplot_options.box.xlabel, "Bin");
	  strcpy(pgplot_options.box.ylabel, "Intensity");
	  strcpy(pgplot_options.box.title, "Select onpulse region(s) (which will be ignored from baseline calculation) ");
	  strcat(pgplot_options.box.title, datain.psrname);
	  selectRegions(profileI, nrBins, &pgplot_options, 0, 0, 0, &(application.onpulse), application.verbose_state);
	  deviceOpened = 1;
	}else {
	  if(application.verbose_state.debug) {
	    pgplot_options.viewport.dontopen = deviceOpened;
	    strcpy(pgplot_options.box.xlabel, "Bin");
	    strcpy(pgplot_options.box.ylabel, "I");
	    strcpy(pgplot_options.box.title, "Profile");
	    pgplotGraph1(&pgplot_options, profileI, NULL, NULL, nrBins, 0, nrBins, 0, 0, nrBins, 0, 0, 0, 1, 0, 1, 0, 1, 1, &(application.onpulse), -1, application.verbose_state);
	    deviceOpened = 1;
	    printf("Press return to continue\n");
	    scanf("%c", txt);
	  }
	}
      }
//START REGION DEVELOP
    }
//START REGION RELEASE
    /* Convert the specified onpulse regions which were defined as
       fractions on the commandline. This is now possible because we
       know the number of bins. */
    region_int_to_frac(&(application.onpulse), 1.0/(float)datain.NrBins, 0);
    regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);
    
//START REGION DEVELOP
    if(clean_flag == 1) {
      if(application.verbose_state.verbose) {
	printf("Find highly correlated off-pulse regions\n");
      }
      /* Correl is now just RMS, because I set Ilast to I. */
      read_correlPSRData(datain, rms, runningBaseline, zapMask, &(application.onpulse), 0, 0, application.verbose_state);
      do {
	j = zapRMS(nrPulses, rms, indx, zapMask, rms_thresh2, 0, 1, 0.8, application.pgplotdevice, &deviceOpened, application.verbose_state);
	for(i = 0; i < nrPulses; i++)
	  if(zapMask[i])
	    rms[i] = 0;
      }while(j != 0);
      pgplot_options.viewport.dontopen = deviceOpened;
      strcpy(pgplot_options.box.xlabel, "Pulse number");
      strcpy(pgplot_options.box.ylabel, "Correl");
      strcpy(pgplot_options.box.title, "");
      pgplotGraph1(&pgplot_options, rms, NULL, NULL, nrPulses, 0, 0, nrPulses, 0, nrPulses, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
      deviceOpened = 1;
      
      if(application.verbose_state.verbose) printf("Find high off-pulse peaks\n");
      /* For zapped bins I is set to zero, what is not working for
	 smoothing. Moreover, the zapping is not checked when looking
	 for maxima, which is also not working in the case of a DC
	 component. */
      read_peakrmsPSRData(datain, rms, runningBaseline, zapMask, &(application.onpulse), 0, 1, 0, application.verbose_state);
      do {
	j = zapRMS(nrPulses, rms, indx, zapMask, rms_thresh3, 1, 1, 0.8, application.pgplotdevice, &deviceOpened, application.verbose_state);
	for(i = 0; i < nrPulses; i++)
	  if(zapMask[i])
	    rms[i] = 0;
      }while(j != 0);
      pgplot_options.viewport.dontopen = deviceOpened;
      strcpy(pgplot_options.box.xlabel, "Pulse number");
      strcpy(pgplot_options.box.ylabel, "Correl");
      strcpy(pgplot_options.box.title, "");
      pgplotGraph1(&pgplot_options, rms, NULL, NULL, nrPulses, 0, nrPulses, 0, 0, nrPulses, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
      deviceOpened = 1;
      
      if(application.verbose_state.verbose) printf("Find high off-pulse peaks of smoothed data\n");
      read_peakrmsPSRData(datain, rms, runningBaseline, zapMask, &(application.onpulse), 0, nrBins/50, 0, application.verbose_state);
      do {
	j = zapRMS(nrPulses, rms, indx, zapMask, rms_thresh3, 1, 1, 0.8, application.pgplotdevice, &deviceOpened, application.verbose_state);
	for(i = 0; i < nrPulses; i++) {
	  if(zapMask[i]) {
	    rms[i] = 0;
	  }
	}
      }while(j != 0);
      printf("-clean removed the following subints:");
      for(i = 0; i < nrPulses; i++) {
	if(zapMask[i]) {
	  printf(" %ld", i);
	}
      }
      printf("\n");
      pgplot_options.viewport.dontopen = deviceOpened;
      strcpy(pgplot_options.box.xlabel, "Pulse number");
      strcpy(pgplot_options.box.ylabel, "Correl");
      strcpy(pgplot_options.box.title, "");
      pgplotGraph1(&pgplot_options, rms, NULL, NULL, nrPulses, 0, nrPulses, 0, 0, nrPulses, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
      deviceOpened = 1;
    }   /* End of if(clean_flag == 1) */

//START REGION RELEASE
    /* Make sure that the pulses for which the running baseline is
       undefined are not zapped. This prevents that these pulses are not
       subtracted from the total number of pulses twice. */
    if(baseline_length) {
      for(i = 0; i < baseline_length; i++)
	zapMask[i] = 0;
      for(i = nrPulses-baseline_length; i < nrPulses; i++)
	zapMask[i] = 0;
    }
    /* Zap pulse ranges which are too short */
    if(blockMode != 0) 
      make_blocks(baseline_length, blockSize, nrPulses, &nrOutputBlocks, zapMask, application.verbose_state);
    /* Count the number of zapped pulses */
    if(application.verbose_state.debug) {
      if(baseline_length)
	printf("DEBUG: Zapping %ld subints because using a running baseline\n", 2*baseline_length);
    }
    nrZapped = 0;
    int prev_zap_state = -1;
    long prev_pulse_nr = 0;
    for(i = baseline_length; i < nrPulses-baseline_length; i++) {
      if(zapMask[i]) {
	nrZapped++;
      }
      if(application.verbose_state.debug) {
	if(prev_zap_state != zapMask[i]) {
	  if(zapMask[i]) {
	    if(i != baseline_length)
	      printf("%ld = %ld subints\n", i-1, i-prev_pulse_nr);
	    printf("DEBUG: Zapping subint %ld - ", i);
	    prev_pulse_nr = i;
	  }else {
	    if(i != baseline_length)
	      printf("%ld = %ld subints\n", i-1, i-prev_pulse_nr);
	    printf("DEBUG: Keeping subint %ld - ", i);
	    prev_pulse_nr = i;
	  }
	}
	prev_zap_state = zapMask[i];
      }
    }
    if(application.verbose_state.debug) {
      printf("%ld = %ld subints\n", nrPulses-baseline_length-1, nrPulses-baseline_length-prev_pulse_nr); 
    }
    if(application.verbose_state.verbose) printf("Total number of zapped pulses: %ld\n", nrZapped);
    
    if(outputlist) {    /* Write the list of zapped subintegrations to a file */
      FILE *fout_zap;
      fout_zap = fopen("zaplist.txt", "w");
      if(fout_zap == NULL) {
	printerror(application.verbose_state.debug, "pmod: Cannot open zaplist.txt");
	return 0;
      }
      for(i = 0; i < nrPulses; i++) {
	if(zapMask[i] != 0) {
	  fprintf(fout_zap, "%ld\n", i);
	}
      }
      fclose(fout_zap);
      printf("Zap list is written to zaplist.txt\n");
    }

    /* Copy pulsar parameters */
    copy_params_PSRData(datain, &dataout, application.verbose_state);
    copy_params_PSRData(datain, &dataout2, application.verbose_state);
    if(debase_flag) {
      dataout.isDebase = 1;
      dataout2.isDebase = 1;
    }
    dataout.NrSubints = datain.NrSubints - 2*baseline_length;
    if(remove_pulses_flag)
      dataout.NrSubints -= nrZapped;
 //START REGION DEVELOP
    if(cutoption == 1)
      dataout.NrBins = application.onpulse.right_bin[0] - application.onpulse.left_bin[0] + 1;
//START REGION RELEASE
    if(application.verbose_state.verbose) printf("Output data contains %ld bins, %ld pulses, %ld polarizations and %ld frequencies.\n", dataout.NrBins, dataout.NrSubints, dataout.NrPols, dataout.NrFreqChan);
    if(openPSRData(&dataout, output_name, application.oformat, 1, 0, 0, application.verbose_state) == 0) {
      printf("Cannot open %s\n\n", output_name);
      return 0;
    }
    if(writeHeaderPSRData(&dataout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
      printf("Cannot write header to %s\n\n", output_name);
      return 0;
    }
    //    appendHistoryLine(dataout, application.verbose_state);
    if(zapoption
 //START REGION DEVELOP
       || clean_flag == 1
//START REGION RELEASE
       ) {   /* Open second output file which will contain the zapped subints */
      dataout2.NrBins = nrBins;
      dataout2.NrSubints = nrZapped;
      if(dataout2.NrSubints > 0) {
	if(openPSRData(&dataout2, output_name2, application.oformat, 1, 0, 0, application.verbose_state) == 0) {
	  printf("Cannot open %s\n\n", output_name);
	  return 0;
	}
	if(writeHeaderPSRData(&dataout2, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	  printf("Cannot write header to %s\n\n", output_name);
	  return 0;
	}
	//	appendHistoryLine(dataout2, argc, argv, application.verbose_state);
      }
    }
 //START REGION DEVELOP

    if(debase_slope_flag == 0 && debase_highpass_flag == 0) {
//START REGION RELEASE
      for(k=0; k < nrPol; k++) {
	for(l = 0; l < NrFreqChan; l++) {
	  if(l == 0) {
	    if(application.verbose_state.verbose) {
	      printf("Processing polarization channel %ld (of the %d)\n", k+1, nrPol); 
	    }
	  }
	  /* First determine baseline of each pulse (average of off-pulse region) */
	  if(debase_flag || removeOnPulse_flag
//START REGION DEVELOP
|| removeOffPulse_rms >= 0
//START REGION RELEASE
) {
//START REGION DEVELOP
	    /* Very stupid ofcourse to read all the data to determine baseline
	       and than again to do the subtracting, but it works and it is
	       easy to implement. It also will not use huge amount of memory
	       (and swapspace), which could actually be quicker for very large
	       files. 

	       Read in table of rms and averages for the present
	       polarization channel and present frequency channel of
	       ALL subintegrations. Store each polarization
	       seperately, probably because I want to plot the
	       baseline of Stokes I afterwards.
	    */
//START REGION RELEASE
	    
//	    printf("XXXX onpulse = %d %d\n", application.onpulse.left_bin[0], application.onpulse.right_bin[0]);
	    if(read_rmsPSRData(datain, &rms[datain.NrSubints*k], &runningBaseline[datain.NrSubints*k], zapMask, &(application.onpulse), 0, k, l, application.verbose_state) == 0) {
	      printerror(application.verbose_state.debug, "Error pmod: Cannot determine offpulse rms.");
	      return 0;
	    }

	    /* Secondly, make the running baseline which will be subtracted */
	    if(l == 0 && application.verbose_state.verbose) {
	      if(datain.NrFreqChan > 1)
		printf("  Preparing baseline for first frequency channel\n");
	      else
		printf("  Preparing baseline\n");
	    }
	    for(i = 0; i < baseline_length; i++) {
	      baseline[i+datain.NrSubints*k] = runningBaseline[i+datain.NrSubints*k];
	      runningRMS[i+datain.NrSubints*k] = rms[i+datain.NrSubints*k];
	    }
	    for(i = nrPulses-baseline_length; i < nrPulses; i++) {
	      baseline[i+datain.NrSubints*k] = runningBaseline[i+datain.NrSubints*k];
	      runningRMS[i+datain.NrSubints*k] = rms[i+datain.NrSubints*k];
	    }
	    int itteration, ok;
	    float rms_old;
//START REGION DEVELOP
	    float baseline_old = 0;
//START REGION RELEASE
	    rms_old  = 0;
	    for(i = baseline_length; i < nrPulses-baseline_length; i++) {

	      for(itteration = 0; itteration < 2; itteration++) {
		baseline[i+datain.NrSubints*k] = 0;
		runningRMS[i+datain.NrSubints*k] = 0;
		tmp_used_pulses = 0;
		for(j = -baseline_length; j <= baseline_length; j++) {
		  /* Ignore zapped pulses */
		  if(runningBaseline[i+j+datain.NrSubints*k] != 0.0) {
//START REGION DEVELOP
		    if(debase_sigma < 0 || itteration == 0) {	
//START REGION RELEASE
		      ok = 1;
//START REGION DEVELOP
		    }else {
		      ok = 0;
		      if(fabs(runningBaseline[i+j+datain.NrSubints*k] - baseline_old)/rms_old < debase_sigma)
			ok = 1;
		    }
//START REGION RELEASE
		    if(ok) {
		      baseline[i+datain.NrSubints*k] += runningBaseline[i+j+datain.NrSubints*k];
		      runningRMS[i+datain.NrSubints*k] += rms[i+j+datain.NrSubints*k];
		      rms_old += runningBaseline[i+j+datain.NrSubints*k]*runningBaseline[i+j+datain.NrSubints*k];
		      tmp_used_pulses++;
		    }
		    //		}else {
		    //		  printf("Found zapped pulse\n");
		  }
		}
		if(tmp_used_pulses > 0) {
		  baseline[i+datain.NrSubints*k] /= (float)(tmp_used_pulses);
		  runningRMS[i+datain.NrSubints*k] /= (float)(tmp_used_pulses);
		  rms_old /= (float)(tmp_used_pulses);
		  rms_old -= baseline[i+datain.NrSubints*k]*baseline[i+datain.NrSubints*k];
		}
//START REGION DEVELOP
		baseline_old = baseline[i+datain.NrSubints*k];
		if(debase_sigma <= 0)
//START REGION RELEASE
		  break;
	      }
	      m = i/10;
	      if((10*m == i) && (m > 0) && application.verbose_state.nocounters == 0 && l == 0 && application.verbose_state.verbose) {
		printf("    %.1f%%     \r",(100.0*(i-baseline_length))/(float)(nrPulses-2.0*baseline_length));
		fflush(stdout);
	      }
	    }
	    if(l == 0 && application.verbose_state.verbose)
	      printf("    Done                      \n");
	  }

	  /* Do the subtracting of the baseline etc and writing out of the data */
	
	  if(debase_flag && l == 0)  {
	    if(application.verbose_state.verbose)
	      printf("  Subtracting baseline\n");
	  }
	  dataout2_pulse = 0;
	  pulse_nr_in_output = 0;
	  for(i = 0; i < nrPulses; i++) {
	    readPulsePSRData(&datain, i, k, l, 0, nrBins, profileI, application.verbose_state);
	    /*	    printf("XXX reading %ld %f\n", i, profileI[0]); */
	    if(debase_flag && debase_offset_flag == 0) {
	      for(j = 0; j < nrBins; j++) {
		// Do not remove baseline from pulses which appear to be zapped, as running baseline makes them non-zero
		if(runningBaseline[i+datain.NrSubints*k] != 0.0) {
		  profileI[j] -= baseline[i+datain.NrSubints*k];
		}
	      }
	    }
	    if(debase_flag && debase_offset_flag) {
	      for(j = 0; j < nrBins; j++) {
		if(runningBaseline[i+datain.NrSubints*k] != 0.0) {
		  profileI[j] -= debase_offset_value;
		}
	      }
	    }
	    if(removeOnPulse_flag) {
	      for(j = 0; j < nrBins; j++) {
		if(checkRegions(j, &(application.onpulse), 0, application.verbose_state) != 0) {
		  /* Old NR code:
		     profileI[j] = runningRMS[i+datain.NrSubints*k]*gasdev(&idnum);
		  */
		  profileI[j] = gsl_ran_gaussian(rand_num_gen, runningRMS[i+datain.NrSubints*k]);
		}
	      }
	    }
//START REGION DEVELOP
	    if(removeOffPulse_rms >= 0.0) {
	      for(j = 0; j < nrBins; j++) {
		if(checkRegions(j, &(application.onpulse), 0, application.verbose_state) == 0) {
		  profileI[j] = gsl_ran_gaussian(rand_num_gen, removeOffPulse_rms);
		}
	      }
	    }
	    if(correctbaseline) {
	      for(j = 0; j < nrBins; j++) {
		I = correctbaseline_amplitude*profileI[j];
		for(j2 = 1; j2 < nrBins; j2++) {
		  if(j+j2 < nrBins) {
		    /*		    profileI[j+j2] += 4e-3*profileI[j]*exp(-(j2)/50.0); */
		    profileI[j+j2] += I*exp(-((float)j2)/(float)correctbaseline_bins); 
		  }
		}
		/*	    profileI[j] -= 1e-3*correctbaseline_sum;
		  profileI[j] = -exp(-(j-1000)/50.0);
		  if(profileI[j] < -1)
		  profileI[j] = -1;*/
		/*		correctbaseline_sum += profileI[j]; */
	      }
	    }
	    if(correctbaseline2) {
	      for(j = nrBins-1; j >= 0; j--) {
		I = correctbaseline2_amplitude*profileI[j];
		for(j2 = 1; j2 < nrBins; j2++) {
		  if(j-j2 >= 0) {
		    profileI[j-j2] += I*exp(-((float)j2)/(float)correctbaseline2_bins); 
		  }
		}
	      }
	    }
	    if(poly > 0) {
	      /*
	      fitfunc_collection_type function;
	      function.nrfuncs = poly+1;
	      for(j = 0; j <= poly; j++) {
		function.func[j].type = FUNC_POLYNOMAL;
		function.func[j].param[0] = j;
		function.func[j].start[0] = 0;
		function.func[j].fit_flag[0] = 1;
	      }
	      function.nrfuncs = poly+1;
	      j2 = 0;
	      for(j = 0; j < nrBins; j++) {
		if(checkRegions(j, &(application.onpulse), 0, application.verbose_state) == 0) {
		  profileI2[j2] = profileI[j];
		  profilex[j2] = j;
		  profileSigma[j2] = 1;
		  j2++;
		}
	      }
	      print_fitfunctions(function, 0, 1, -1);
	      fit_levmar(1, &function, profilex, profileI2, profileSigma, j2, 0, 0, application.verbose, 0, application.verbose_state);
	      print_fitfunctions(function, 0, 1, -1);
	      for(j = 0; j < nrBins; j++) {
		profileI[j] -= evaluate_fitfunction_collection(function, j);
		if(application.verbose_state.verbose)
		  profileI[j] = evaluate_fitfunction_collection(function, j);
	      }
	      */
	      fitfunc_collection_type function;
	      function.nrfuncs = poly+1;
	      for(j = 0; j <= poly; j++) {
		function.func[j].type = FUNC_POLYNOMAL;
		function.func[j].param[0] = j;
		function.func[j].start[0] = 0;
		function.func[j].fit_flag[0] = 1;
	      }
	      function.nrfuncs = poly+1;
	      j2 = 0;
	      for(j = 0; j < nrBins; j++) {
		if(checkRegions(j, &(application.onpulse), 0, application.verbose_state) == 0) {
		  profileI2_double[j2] = profileI[j];
		  profilex_double[j2] = j;
		  profileSigma_double[j2] = 1;
		  j2++;
		}
	      }
	      print_fitfunctions(&function, 0, 1, -1, application.verbose_state);
	      fit_levmar(1, &function, profilex_double, profileI2_double, profileSigma_double, j2, 0, 0, 1e-4, 1e-4, 500, &status, 0, 0, application.verbose_state);
	      print_fitfunctions(&function, 0, 1, -1, application.verbose_state);
	      for(j = 0; j < nrBins; j++) {
		/*		printf("XXXX %Lf\n", evaluate_fitfunction_collection(function, j)); */
		if(poly_replace)
		  profileI[j] = evaluate_fitfunc_collection(&function, j, application.verbose_state);
		else
		  profileI[j] -= evaluate_fitfunc_collection(&function, j, application.verbose_state);
	      }

	    } // End of if(poly > 0)
//START REGION RELEASE
	    if(addnoise_flag) {
	      for(j = 0; j < nrBins; j++) {
		/* Old NR code:
		   profileI[j] += noiseRMS*gasdev(&idnum);
		*/
		profileI[j] += gsl_ran_gaussian(rand_num_gen, noiseRMS);
	      }
	    }
//START REGION DEVELOP
	    if(threshold_flag == 1) {
	      for(j = 0; j < nrBins; j++) {
		if(threshold_excess_flag == 0) {
		  if(profileI[j] > threshold_val)
		    profileI[j] = threshold_val;
		  else if(profileI[j] < -threshold_val)
		    profileI[j] = -threshold_val;
		}else {
		  if(profileI[j] > threshold_val)
		    profileI[j] -= threshold_val;
		  else if(profileI[j] < -threshold_val)
		    profileI[j] += -threshold_val;
		  else 
		    profileI[j] = 0;
		}
	      }
	    }
	    if(clippos_flag) {
	      for(j = 0; j < nrBins; j++) {
		if(threshold_excess_flag == 0) {
		  if(profileI[j] > clippos_val)
		    profileI[j] = clippos_val;
		}else {
		  if(profileI[j] > clippos_val)
		    profileI[j] -= clippos_val;
		  else 
		    profileI[j] = 0;
		}
	      }
	    }
	    if(clipneg_flag) {
	      for(j = 0; j < nrBins; j++) {
		if(threshold_excess_flag == 0) {
		  if(profileI[j] < clipneg_val)
		    profileI[j] = clipneg_val;
		}else {
		  if(profileI[j] < clipneg_val)
		    profileI[j] -= clipneg_val;
		  else 
		    profileI[j] = 0;
		}
	      }
	    }
	    if(threshold_flag == 4) {
	      for(j = 0; j < nrBins; j++) {
		if(profileI[j] < -threshold_val)
		  profileI[j] = -threshold_val;
		else if(profileI[j] > threshold_val)
		  profileI[j] = threshold_val;
		else
		  profileI[j] = 0;
	      }
	    }
//START REGION RELEASE	    
	    /* Write out zapped pulses to separate pulse-stack */
	    if(zapoption
//START REGION DEVELOP
	       || clean_flag == 1
//START REGION RELEASE	    
	       ) {
	      if(zapMask[i] != 0) {
		if(writePulsePSRData(&dataout2, dataout2_pulse, k, l, 0, nrBins, profileI, application.verbose_state) != 1) {
		  printerror(application.verbose_state.debug, "pmod: Writing data failed.");
		  return 0;
		}
		dataout2_pulse += 1;
	      }
	    }
	    if(zapMask[i] != 0) {
	      /*	      printf("XXX zapping %ld\n", i); */
	      for(j = 0; j < nrBins; j++)
		profileI[j] = 0;
	      baseline[i] = 0;
	    }
	    if(fzapMask[l] != 0) {
	      /*	      printf("YYY zapping %ld\n", i); */
	      for(j = 0; j < nrBins; j++)
		profileI[j] = 0;
	      baseline[i] = 0;
	    }
	    /* skip first and last pulses */
	    if(i >= baseline_length && i < nrPulses-baseline_length) {    
	      /* In block mode, do not write emty pulses */          
	      if(remove_pulses_flag == 0 || zapMask[i] == 0) {                     
		if(remove_pulses_flag == 0) {
		  pulse_nr_in_output = i;
		}
//START REGION DEVELOP
		if(cutoption == 1) {
		  if(writePulsePSRData(&dataout, pulse_nr_in_output, k, l, 0, dataout.NrBins, profileI+application.onpulse.left_bin[0], application.verbose_state) != 1) {
		    printerror(application.verbose_state.debug, "pmod: Writing data failed.");
		    return 0;
		  }
		}else {
//START REGION RELEASE
		  if(writePulsePSRData(&dataout, pulse_nr_in_output, k, l, 0, nrBins, profileI, application.verbose_state) != 1) {
		    printerror(application.verbose_state.debug, "pmod: Writing data failed.");
		    return 0;
		  }
//START REGION DEVELOP
		  /*		  printf("XXX writing %ld %f\n", i, profileI[0]); */
		}
//START REGION RELEASE
		pulse_nr_in_output++;
	      }
	    }
	    m = i/10;
	    if((10*m == i) && (m > 0) && application.verbose_state.nocounters == 0 && l != 0 && application.verbose_state.verbose) {
	      printf("   %.1f%%     \r",(100.0*(i+nrPulses*l))/(float)(nrPulses*NrFreqChan));
	      fflush(stdout);
	    }
	  }
	  if(NrFreqChan > 1 && l == 0 && application.verbose_state.verbose)
	    printf("  Processing other frequency channels\n");
	}
	if(application.verbose_state.verbose)
	  printf("    Done                             \n");
	
	if(k == 0 && (debase_flag && debase_offset_flag==0) && nrPulses-2*baseline_length > 1) {
	  strcpy(txt, "Baseline ");
	  strcat(txt, datain.psrname);
	  pgplot_options.viewport.dontopen = deviceOpened;
	  strcpy(pgplot_options.box.xlabel, "Pulse number");
	  strcpy(pgplot_options.box.ylabel, "Baseline");
	  strcpy(pgplot_options.box.title, txt);
	  pgplotGraph1(&pgplot_options, baseline+baseline_length, NULL, NULL, nrPulses-2*baseline_length, baseline_length, nrPulses-baseline_length, 0, baseline_length, nrPulses-baseline_length, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
	  deviceOpened = 1;
	}
      }
 //START REGION DEVELOP
    } // End of if(debase_slope_flag == 0)
  
    if(debase_slope_flag || debase_highpass_flag) {
      // Note: this is very similar to the code in preprocess_debase()
      // If only fitting for a linear slope, we can use linear regression
      if(debase_slope_flag && debase_highpass_flag == 0) {
	double a, b, cov00, cov01, cov11, sumsq;
	float value_float;
	long n;
	n = 0;
	for(j = 0; j < nrBins; j++) {
	  if(checkRegions(j, &(application.onpulse), 0, application.verbose_state) == 0)
	    profilex_double[n++] = j;
	}
	for(k=0; k < nrPol; k++) {
	  for(l = 0; l < NrFreqChan; l++) {
	    for(i = 0; i < nrPulses; i++) {
	      n = 0;
	      for(j = 0; j < nrBins; j++) {
		if(checkRegions(j, &(application.onpulse), 0, application.verbose_state) == 0) {
		  readPulsePSRData(&datain, i, k, l, j, 1, &value_float, application.verbose_state);
		  profileI2_double[n] = value_float; 
		  n++;
		}
	      }
	      
	      // For plotting purposes, subtract average 
	      //	  float avrg = 0;
	//		  for(j = 0; j < nrBins; j++)
	//		  avrg += profileI2_double[j];
	//		  avrg /= (float)nrBins;
	//		  for(j = 0; j < nrBins; j++)
	//		  profileI2_double[j] -= avrg;
	      
//	    for(j = 0; j < n; j++) {
//			    printf("%f %f FIT DATA\n", profilex_double[j], profileI2_double[j]);
//			    }
	      // Old NR code:
//		 fit(profilex_double-1, profileI2_double-1, n, NULL, 0, &a, &b, &siga, &sigb, &chi2, &q); 
	      gsl_fit_linear(profilex_double, 1, profileI2_double, 1, n, &a, &b, &cov00, &cov01, &cov11, &sumsq);
	      
	      readPulsePSRData(&datain, i, k, l, 0, nrBins, profileI, application.verbose_state); 
	      for(j = 0; j < nrBins; j++) {
		profileI[j] -= a+(float)j*b;  
	      }
	      if(writePulsePSRData(&dataout, i, k, l, 0, nrBins, profileI, application.verbose_state) != 1) {
		printerror(application.verbose_state.debug, "pmod: Writing data failed.");
		return 0;
	      }
	    }
	  }
	}
      }else { // End of if(debase_slope_flag && debase_highpass_flag == 0)
      // Fitting sinusoids: use MRQ fitting. This can be combined with fitting for a linear fit as well.
	float value_float;
	long n;
	int ret;
	n = 0;
	for(j = 0; j < nrBins; j++) {
	  if(checkRegions(j, &(application.onpulse), 0, application.verbose_state) == 0)
	    profilex_double[n++] = j;
	}
	for(k=0; k < nrPol; k++) {
	  for(l = 0; l < NrFreqChan; l++) {
	    for(i = 0; i < nrPulses; i++) {

	      // Construct the fit function
	      fitfunc_collection_type fitfunction;
	      n = 0;
	      fitfunction.func[n].type = FUNC_POLYNOMAL;  // Always fit for a constant offset
	      fitfunction.func[n].param[0] = 0;
	      fitfunction.func[n].start[0] = 0;
	      fitfunction.func[n].value[0] = 0;
	      fitfunction.func[n].fit_flag[0] = 1;
	      n++;
	      if(debase_slope_flag) {  // Do we want to fit for a slope?
		fitfunction.func[n].type = FUNC_POLYNOMAL;  // Always fit for a constant offset
		fitfunction.func[n].param[0] = 1;
		fitfunction.func[n].start[0] = 0;
		fitfunction.func[n].value[0] = 0;
		fitfunction.func[n].fit_flag[0] = 1;
		n++;
	      }
	      if(debase_highpass_nrharmonics > 0) {
		for(j = 0; j < debase_highpass_nrharmonics; j++) {
		  fitfunction.func[n].type = FUNC_SIN;
		  fitfunction.func[n].start[0] = 0;
		  fitfunction.func[n].start[1] = 2.0*M_PI*(double)(j+1)/(double)nrBins;
		  fitfunction.func[n].start[2] = 0;
		  fitfunction.func[n].value[0] = fitfunction.func[n].start[0];
		  fitfunction.func[n].value[1] = fitfunction.func[n].start[1];
		  fitfunction.func[n].value[2] = fitfunction.func[n].start[2];
		  fitfunction.func[n].fit_flag[0] = 1;
		  fitfunction.func[n].fit_flag[1] = 0;
		  fitfunction.func[n].fit_flag[2] = 0;
		  n++;		  

		  fitfunction.func[n].type = FUNC_COS;
		  fitfunction.func[n].start[0] = 0;
		  fitfunction.func[n].start[1] = 2.0*M_PI*(double)(j+1)/(double)nrBins;
		  fitfunction.func[n].start[2] = 0;
		  fitfunction.func[n].value[0] = fitfunction.func[n].start[0];
		  fitfunction.func[n].value[1] = fitfunction.func[n].start[1];
		  fitfunction.func[n].value[2] = fitfunction.func[n].start[2];
		  fitfunction.func[n].fit_flag[0] = 1;
		  fitfunction.func[n].fit_flag[1] = 0;
		  fitfunction.func[n].fit_flag[2] = 0;
		  n++;		  
		}
	      }
	      fitfunction.nrfuncs = n;

	      // Construct the x values (used bins)
	      n = 0;
	      for(j = 0; j < nrBins; j++) {
		if(checkRegions(j, &(application.onpulse), 0, application.verbose_state) == 0) {
		  readPulsePSRData(&datain, i, k, l, j, 1, &value_float, application.verbose_state);
		  profileI2_double[n] = value_float; 
		  n++;
		}
	      }

	      // Do fitting here
	      //	      gsl_fit_linear(profilex_double, 1, profileI2_double, 1, n, &a, &b, &cov00, &cov01, &cov11, &sumsq);


	      fit_levmar(1, &fitfunction, profilex_double, profileI2_double, NULL, n, 0, 1, 0.0, 0.0, 1000, &ret, 0, 0, application.verbose_state);
	      if(ret != 0 && ret != 3) {
		printerror(application.verbose_state.debug, "pmod: Fitting failed with code %d", ret);
		return 0;
	      }

	      // Subtract fit from data
	      readPulsePSRData(&datain, i, k, l, 0, nrBins, profileI, application.verbose_state); 
	      for(j = 0; j < nrBins; j++) {
		/*
		  if(l == 0 && k == 0) {
		  printf("%f %f %f\n", (float)j, profileI[j], a+profilex_double[j]*b);
		  }*/
		profileI[j] -= evaluate_fitfunc_collection(&fitfunction, j, application.verbose_state);

	      }
	      if(writePulsePSRData(&dataout, i, k, l, 0, nrBins, profileI, application.verbose_state) != 1) {
		printerror(application.verbose_state.debug, "pmod: Writing data failed.");
		return 0;
	      }
	    }
	  }
	}
      } // End of if(debase_slope_flag && debase_highpass_flag == 0)
      /* THE FOLLOWING DOESN'T QUITE WORK. REPLACING THE ONPULSE REGION WITH NOISE + ESTIMATED BASELINE DOES NOT QUITE REMOVE THE EFFECT OF THE ONPULSE REGION IN AN ITTERATIVE WAY AND SIGNIFICANT RESIDUALS ARE LEFT. INSTEAD DO ACTUAL FITTING, AS IS IMPLEMENTED ABOVE
      if(debase_highpass_flag) {
	float *profileCopy;       // Copy of original
	float *profileSmooth;     // Smoothed baseline
	float *profileCorrected;  // Pulse - smoothed baseline
	float factor;
	int itt;
	factor = 1.0/nrBins; // The normalisation factor
	profileCopy = (float *)malloc(nrBins*sizeof(float));
	profileSmooth = (float *)malloc(nrBins*sizeof(float));
	profileCorrected = (float *)malloc(nrBins*sizeof(float));
	if(profileCorrected == NULL || profileCopy == NULL || profileSmooth == NULL) {
	  printerror(application.verbose_state.debug, "ERROR pmod: Memory allocation error.");
	  return 0;
	}

	for(k=0; k < nrPol; k++) {
	  for(l = 0; l < NrFreqChan; l++) {
	    for(i = 0; i < nrPulses; i++) {
	      readPulsePSRData(&datain, i, k, l, 0, nrBins, profileI, application.verbose_state);
	      memcpy(profileCopy, profileI, nrBins*sizeof(float));
	      memcpy(profileCorrected, profileI, nrBins*sizeof(float));  // Make corrected equal to the original, to make itterative code below work for the first itteration
	      memset(profileSmooth, 0, nrBins*sizeof(float));  // Make smoothed equal to zero, to make itterative code below work for the first itteration
	      for(itt = 0; itt < 4; itt++) {
		float baseline_value, rms_value;
		offpulseStats(profileCorrected, nrBins, &baseline_value, &rms_value, &(application.onpulse), 0, application.verbose_state);
		if(application.verbose_state.debug) {
		  printf("nsub=%ld, nfreq=%ld, pol=%ld: avrg=%f rms=%f\n", i, l, k, baseline_value, rms_value);
		}
		// Make pulse equal to the copy, but replace pulse with smooth baseline + noise
		for(j = 0; j < nrBins; j++) {
		  if(checkRegions(j, &(application.onpulse), 0, application.verbose_state) != 0) {
		    profileI[j] = profileSmooth[j] + baseline_value + gsl_ran_gaussian(rand_num_gen, rms_value);
		  }else {
		    profileI[j] = profileCopy[j];
		  }
		}
		//		writePulsePSRData(&dataout, i, k, l, 0, nrBins, profileI, application.verbose_state); 
		if(fftSmooth(profileI, nrBins, 0, nrBins/2-debase_highpass_nrharmonics, application.verbose_state) == 0) {
		  printerror(application.verbose_state.debug, "ERROR pmod: Applying low-pass filter failed.");
		  return 0;
		}
		// Normalise result
		for(j = 0; j < nrBins; j++) {
		  profileI[j] *= factor;
		}
		memcpy(profileSmooth, profileI, nrBins*sizeof(float));		
		//writePulsePSRData(&dataout, i, k, l, 0, nrBins, profileSmooth, application.verbose_state); 
		for(j = 0; j < nrBins; j++) {
		  profileCorrected[j] = profileCopy[j] - profileSmooth[j];
		}
	      }
	      writePulsePSRData(&dataout, i, k, l, 0, nrBins, profileCorrected, application.verbose_state); 
	    }
	  }
	}
	free(profileCorrected);
	free(profileCopy);
      } // End of if(debase_highpass_flag) 
      */
    } // End of if(debase_slope_flag || debase_highpass_flag)
  
    /*
      int deviceID;
      deviceID = cpgopen(application.pgplotdevice);
      cpgask(0);
      cpgslw(1);
      
      if(output_postscript) {
      printf("Writing output to profile.ps\n");
      cpgopen("profile.ps/ps");
      cpgask(0);
      cpgslw(1);
      PlotProfile(nrBins, profileI, "bin number", "Intensity", txt, 1);
      cpgclos();
      cpgslct(deviceID);
      }
      PlotProfile(nrBins, profileI, "bin number", "Intensity", txt, 1);
      
      if(output_postscript) {
      printf("Writing output to baseline.ps\n");
      cpgopen("baseline.ps/ps");
      cpgask(0);
      cpgslw(1);
      PlotProfile(nrPulses-2*baseline_length, baseline+baseline_length, "pulse number", "Baseline", txt, 0);
	cpgclos();
	cpgslct(deviceID);
	}
  
*/
//START REGION RELEASE
    free(profileI);
    free(baseline); 
    free(rms);
    free(runningBaseline);
    free(runningRMS);
    free(zapMask);
    free(fzapMask);
//START REGION DEVELOP
    free(profileI2_double);
    free(profilex_double);
    free(profileSigma_double);
    free(indx);
//START REGION RELEASE
    closePSRData(&datain, 0, application.verbose_state);
    closePSRData(&dataout, 0, application.verbose_state);
    closePSRData(&dataout2, 0, application.verbose_state);
    if(deviceOpened) {
      ppgclos();
      ppgend();
    }
    printf("Writing of %s done\n", output_name);
  }  // End of while loop over input file names

//START REGION DEVELOP
  freePulselongitudeRegion(&cutSelection);
//START REGION RELEASE
  gsl_rng_free(rand_num_gen); 
  terminateApplication(&application);
  return 0;
}   // End of main()
 
//START REGION DEVELOP
//START REGION RELEASE

// If zapColumn2 > 0, two columns are read in being interpretted as start/end numbers
int readZapFile(char *zaplistname, int *zapMask, int zapSkipLines, int nrZapCols, int zapColumn, int zapColumn2, int inverseZap, verbose_definition verbose)
{
  FILE *fin_zap;
  long i, i2, j;
  char txt[100];

  if(zapColumn2 > 0) {
    if(zapColumn2 <= zapColumn) {
      printerror(verbose.debug, "pmod: column %d was expected to be larger than %d.", zapColumn2, zapColumn);
      return 0;
    }
  }

  fin_zap = fopen(zaplistname, "r");
  if(fin_zap == NULL) {
    printerror(verbose.debug, "pmod: Cannot open %s", zaplistname);
    return 0;
  }else {
    fprintf(stdout, "Opened zapfile '%s'\n", zaplistname);
  }
  if(zapSkipLines > 0) {
    printf("Skipping %d lines\n", zapSkipLines);
    for(i = 0; i < zapSkipLines; i++) {
      do {
	j = fgetc(fin_zap);
      }while(j !='\n');
    }
  }
  do {
    if(zapColumn > 1) {           /* Have to skip some columns */
      for(i = 0; i < zapColumn-1; i++) {
	j = fscanf(fin_zap, "%s", txt);
	if(j != 1)
	  break;
	/*	printf("SKIPPED %s\n", txt); */
      }
    }
    j = fscanf(fin_zap, "%ld", &i);              /* Read in pulse number */
    if(zapColumn2 > 0 && j == 1) {           
      if(zapColumn2 - zapColumn > 1) {           /* Have to skip some columns */
	for(i = zapColumn+1; i < zapColumn2; i++) {
	  j = fscanf(fin_zap, "%s", txt);
	  if(j != 1)
	    break;
	  /*	printf("SKIPPED %s\n", txt); */
	}
      }
      j = fscanf(fin_zap, "%ld", &i2);              /* Read in pulse number */
    }
    if(j == 1) {
      if(verbose.debug) {
	printf("DEBUG: read in value %ld", i); 
	if(zapColumn2 > 0)
	  printf(" and %ld", i2);
	printf(" (inverseZap = %d)\n", inverseZap);
      }
      if(zapColumn2 <= 0) {
	if(inverseZap == 0) {
	  zapMask[i] = 1;
	}else {
	  zapMask[i] = 0;
	}
      }else {
	long loopnr;
	for(loopnr = i; loopnr <= i2; loopnr++) {
	  if(inverseZap == 0) {
	    zapMask[loopnr] = 1;
	  }else {
	    zapMask[loopnr] = 0;
	  }
	}
      }
    }
    long nskip;
    if(zapColumn2 > 0)
      nskip = nrZapCols - zapColumn2;
    else
      nskip = nrZapCols - zapColumn;
    if(nskip > 0) {         /* Have to skip some columns */
      for(i = 0; i < nrZapCols - zapColumn; i++) {
	j = fscanf(fin_zap, "%s", txt);
	if(j != 1)
	  break;
      }
    }
  }while(j == 1);
  fclose(fin_zap);
  return 1;
}

//START REGION DEVELOP
//START REGION RELEASE

void make_blocks(long baseline_length, long blockSize, long nrPulses, long *nrOutputBlocks, int *zapMask, verbose_definition verbose)
{
  long i, j, k, nrZapped;
  nrZapped = 0;
  i = baseline_length;
  do {
    k = 1;
    for(j = i; j < i+blockSize; j++) {      /* Find out of next blockSize pulses exist */
      //      fprintf(stderr, "XXXX j=%ld i=%ld  blocksize=%ld  nrPulses=%ld baseline_length=%ld\n", j, i, blockSize, nrPulses, baseline_length);
      if(j == nrPulses-baseline_length) {
	k = -j;
	break;
      }
      if(zapMask[j]) {
	k = -j;
	break;
      }
    }
    if(k <= 0) {                       /* There are no blockSize pulses available */
      /*	printf("Zap pulse %ld t/m %ld\n", i, -k); */
      for(j = i; j <= -k; j++) {
	if(j < nrPulses) {
	  zapMask[j] = 1;
	  nrZapped++;
	}
      }
      i = -k+1;
    }else {                            /* There are is a complete block availavle */
      /*	printf("Found block starting from %ld\n", i); */
      i += blockSize;
      nrOutputBlocks++;
    }
  }while(i < nrPulses-baseline_length);
  if(verbose.verbose) printf("%ld pulses zapped in total to make integer number of %ld pulse blocks.\n", nrZapped, blockSize);
}

//START REGION DEVELOP

/*
Do not zap more than a fraction fracData of the data.

 */
int zapRMS(long nrPulses, float *rms, unsigned long *indx, int *zapMask, float rms_thresh, int subtractMinimum, int subtractMedian, float fracData, char *plotDevice, int *deviceOpened, verbose_definition verbose)
{
  long i, j, nrZapped, imedian;
  float rms2, *newrms, median;
  char txt[100];
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);

  strcpy(pgplot_options.viewport.plotDevice, plotDevice);
  pgplot_options.viewport.dontclose = 1; 
  pgplot_options.viewport.dontopen = *deviceOpened;

  newrms = (float *)malloc(nrPulses*sizeof(float));
  if(newrms == NULL) {
    printerror(verbose.debug, "pmod: Memory allocation error");
    return 0;
  }

  if(verbose.verbose) printf("  Sorting\n");
  /* Old NR code:
     indexx(nrPulses, rms-1, indx-1);
  */
  gsl_sort_float_index(indx, rms, 1, nrPulses);

  nrZapped = 0;
  for(i = 0; i < nrPulses; i++)
    if(zapMask[i] != 0)
      nrZapped++;

  imedian = indx[nrZapped+((nrPulses-nrZapped)/2)];
  median = rms[imedian];
  if(verbose.verbose) printf("  Median  RMS is pulse %ld (rms=%f)\n", imedian+1, median);
  if(verbose.verbose) printf("  Largest RMS is pulse %ld (rms=%f)\n", indx[nrPulses-1]+1, rms[indx[nrPulses-1]]);

  /*
  if(application.verbose_state.verbose) {
    pgplotGraph1(viewport, rms, NULL, NULL, nrPulses, 0, nrPulses, 0, "Pulse number", "RMS", "", 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
    *deviceOpened = 1;
    printf("Press return to continue to see result\n");
    scanf("%c", txt);
  }
  */

  /*
  median = 0;
  j = 0;
  for(i = 0; i < nrPulses; i++)
    if(zapMask[i] == 0) {
      median += rms[i];
      j++;
    }
  median /= (float)j;
  printf("Gemiddelde = %f\n", median);
  */

  /* Subtract minimum value and calculate rms of rms'es */
  if(subtractMinimum || subtractMedian) {
    rms2 = rms[indx[nrZapped]];
    if(subtractMedian)
      rms2 = median;
    printf("  Subtracting %f (%ld zapped)\n", rms2, nrZapped);
  }else {
    rms2 = 0;
  }
  for(i = 0; i < nrPulses; i++) {
    if(zapMask[i] == 0)
      newrms[i] = rms[i]-rms2;
    else
      newrms[i] = 0;
  }

  /* Calculate rms of rms'es */
  rms2 = 0;
  for(i = 0; i < nrPulses; i++) {
    if(zapMask[i] == 0) {
      rms2 += newrms[i]*newrms[i];
      /* if(verbose)      printf("%f\n", newrms[i]); */ 
    }
  }
  rms2 = sqrt(rms2/(float)(nrPulses-nrZapped));
  if(verbose.verbose) printf("  RMS of rms'es = %f (%ld)\n", rms2, nrPulses-nrZapped);

  /* Zap at most a fraction of the data at this initial stage */
  j = 0;
  for(i = 0; i < nrPulses; i++) {
    if((newrms[i] > rms_thresh*rms2 && zapMask[i] == 0) || (newrms[i] < -rms_thresh*rms2 && zapMask[i] == 0)) {
      zapMask[i] = 1;
      /* It is important to set these rms'es to zero for the next
	 round, or else the median is not calculated correctly. */
      rms[i] = 0;
      /* printf("Zapped %ld\n", indx[i]-1); */
      j++;
      nrZapped++;
    }
    if(nrZapped >= fracData*nrPulses)
      break;
  }
  if(verbose.verbose) printf("  Zapped %ld pulses\n", j);
  if(verbose.debug) {
    strcpy(pgplot_options.box.xlabel, "Pulse number");
    strcpy(pgplot_options.box.ylabel, "RMS");
    strcpy(pgplot_options.box.title, "");
    pgplotGraph1(&pgplot_options, newrms, NULL, NULL, nrPulses, 0, nrPulses, 0, 0, nrPulses, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, verbose);
    *deviceOpened = 1;
    ppgmove(0, rms2);
    ppgdraw(nrPulses, rms2);
    ppgsci(2);
    ppgmove(0, rms_thresh*rms2);
    ppgdraw(nrPulses, rms_thresh*rms2);
    ppgmove(0, -rms_thresh*rms2);
    ppgdraw(nrPulses, -rms_thresh*rms2);
    ppgsci(1);
    printf("Press return to continue\n");
    scanf("%c", txt);
  }
  free(newrms);
  return j;
}
