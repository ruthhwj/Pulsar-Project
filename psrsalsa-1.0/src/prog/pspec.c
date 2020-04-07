//START REGION RELEASE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "psrsalsa.h"

#define MaxNrPolarizations 5

//START REGION DEVELOP
void SHOWREVISIONINFO_prog() {
#include "pspec.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);

//START REGION RELEASE
int pgplotProfile(char *plotDevice, int windowwidth, int windowheight, float *profile, float *stddev, float *rms_stddev, float *modindex, float *rms_modindex, int nrx, float xmin, float xmax, char *xlabel, char *ylabel, char *title, int stddev_flag, int mod_flag, int zoom_flag, float xmin_zoom, float xmax_zoom, verbose_definition verbose);
//START REGION DEVELOP

//START REGION RELEASE
int main(int argc, char **argv)
{
  int fft_size, index, originalNrPols, selectMoreOnpulseRegions, powertwo, track_only_first_region;
  int profile_flag, lrfs_flag, stddev_flag, mod_flag, twodfs_flag, bootstrap, subtractDC, track_flag, amplitude_flag, ftrack_mask, inverseFFT, write_flag, modSimple_flag, zoom_flag, zoom_flag1, p2range_set, regionnr, s2dfs_p3_flag, s2dfs_p2_flag;
  long fft_blocks, junk_int;
  long i, j, k, l, p, nrpointsrms;
  float xmin, xmax, xmin_zoom, xmax_zoom, mod_sigma, stddev_sigma, sampleI, freq_min, freq_max, var_rms;
  float *profileI, *lrfs, *lrfs2, *stddev, *modindex, *rms_sigma, *rms_modindex, *twodfs, *clone_profileI, *phase_track, *phase_track_phases, *amplitude_profile, slope, track_dphase;
  float zapmin, zapmax, p2min, p2max, p2, p3, junk_float;
  double *stddev_av, *modindex_av, *stddev_square, *modindex_square, rms, avrg;
  char lrfsdevice[1000], onpulseselectdevice[1000], profiledevice[1000], trackdevice[1000], amplitudedevice[1000], twodfsdevice[1000], outputname[1000], txt[1000], s2dfs_p3_device[1000], s2dfs_p2_device[1000];
  FILE *fout_ascii;
  psrsalsaApplication application;
  datafile_definition fin[MaxNrPolarizations], clone, fout;
  pgplot_options_definition pgplot_options;
  verbose_definition noverbose;
//START REGION DEVELOP
//  int clip_lrfs_flag;
  int s_track_flag, s_modmap_flag, strack_mask, strack_2ndfile_flag, p3_subtract_flag;
  int lrcc_flag, lrcc_lag, lrcc_lagstop, lrcc_nosubmean, hrfs_flag, lrac_flag;
  int lrac_zap0, lrac_keepbothhalves;
  int trackphases_flag;
  int p3class, p3class_nblock, p3class_nmin, p3class_npad; float p3class_cpp1, p3class_cpp2;
  char s_track_device[1000], s_modmap_device[1000];
  char trackphasesdevice[1000], lrccdevice[1000], hrfsdevice[1000], lracdevice[1000];
  float *s_phase_track, *s_mod_map, p3_subtract;
  float strack_dphase, strack_dphase2;
  float *lrcc;

//START REGION RELEASE
  initApplication(&application, "pspec", "[options] inputfile");
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
  application.switch_itf = 1;
  application.switch_rebin = 1;
  application.switch_nskip = 1;
  application.switch_nread = 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_conshift= 1;
  application.switch_circshift= 1;
  application.switch_shuffle = 1;
  application.switch_libversions = 1;
  application.switch_history_cmd_only = 1;
//START REGION DEVELOP
  application.switch_polselect = 1;
  application.switch_ppgplot = 1;
  application.switch_showrevision = 1;


//START REGION RELEASE
  fft_size = 512;
  powertwo = 0;
  lrfs_flag = 0;
  profile_flag = 0;
  stddev_flag = 0;
  mod_flag = 0;
  mod_sigma = -1;
  stddev_sigma = -1;
  bootstrap = 0;
  subtractDC = 1;
  track_flag = 0;
  amplitude_flag = 0;
  freq_min = 0;
  freq_max = 0.5;
  ftrack_mask = 0;
  inverseFFT = 0;
  write_flag = 0;
  modSimple_flag = 0;
  zoom_flag = 0;
  zoom_flag1 = 0;
  slope = 0;
  track_dphase = 0;
  twodfs_flag = 0;
  p2range_set = 0;
  selectMoreOnpulseRegions = 0;
  pgplot_clear_options(&pgplot_options);
  sprintf(lrfsdevice, "?");
  sprintf(onpulseselectdevice, "?");
  sprintf(profiledevice, "?");
  sprintf(trackdevice, "?");
  sprintf(amplitudedevice, "?");
  sprintf(twodfsdevice, "?");
  sprintf(s2dfs_p3_device, "?");
  sprintf(s2dfs_p2_device, "?");
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;
  track_only_first_region = 0;
  s2dfs_p3_flag = 0;
  s2dfs_p2_flag = 0;

//START REGION DEVELOP
//  clip_lrfs_flag = 0;
  s_track_flag = 0;
  s_modmap_flag = 0;
  strack_mask = 0;
  strack_dphase = 0;
  strack_dphase2 = 0;
  strack_2ndfile_flag = 0;
  p3_subtract = 0;
  p3_subtract_flag = 0;
  trackphases_flag = 0;
  lrac_flag = 0;
  lrac_zap0 = 0;
  lrac_keepbothhalves = 0;
  p3class = 0;
  //  application.oformat = PUMA_format;
//START REGION RELEASE
  application.oformat = FITS_format;
//START REGION DEVELOP
  lrcc_flag = 0;
  lrcc_lag = 0;
  lrcc_lagstop = 0;
  lrcc_nosubmean = 0;
  hrfs_flag = 0;
  sprintf(s_track_device, "?");
  sprintf(s_modmap_device, "?");
  sprintf(trackphasesdevice, "?");
  sprintf(lrccdevice, "?");
  sprintf(hrfsdevice, "?");
  sprintf(lracdevice, "?");

//START REGION RELEASE
  if(argv[argc-1][0] == '-' &&  strcmp(argv[argc-1], "-formatlist") != 0 &&  strcmp(argv[argc-1], "-headerlist") != 0) {
    printerror(application.verbose_state.debug, "pspec: Last command line option is expected to be a file name (got %s).\nRun pspec without command line arguments to show help", argv[argc-1]);
    terminateApplication(&application);
    return 0;
  }

  if(argc < 2) {
    printf("Program to analyse (folded) single pulse data using mostly Fourier techniques.\nIt is assumed the data contains a single pulse in each subint and that the\nbaseline is subtracted (use pmod -debase).\n\n");
    printApplicationHelp(&application);
    printf("General options:\n");
    printf("  -nfft               Set size of fft's [default=%d].\n", fft_size);
//START REGION DEVELOP
    printf("  -C                  Select onpulse region (indentical to -onpulse).\n");
//START REGION RELEASE
    printf("  -powertwo           When manually selecting onpulse regions, they are forced\n");
    printf("                      to be a power of two bins wide.\n");
    printf("  -w                  Write out the results to files.\n");
    printf("  -bootstrap          Find error bars on the standard deviation, modulation\n");
    printf("                      index and subpulse phase by random adding noise to the\n");
    printf("                      data. This will be done for the specified number of times\n");
    printf("                      (larger value will be more precise, but takes longer). The\n");
    printf("                      error bars (although somewhat overestimated) are more\n");
    printf("                      accurate than the analytic approximation used by default.\n");
//START REGION RELEASE
    printf("\nOutput options:\n");
    printf("  -prof               Compute pulse profile.\n");
    printf("  -lrfs               Compute LRFS.\n");
    printf("  -DC                 Leave the DC channel in the LRFS.\n");
    printf("  -stddev             Compute standard deviation profile, normalised by the average\n");
    printf("                      intensity in the pulse longitude bin where the profile peaks.\n");
    printf("  -stddev_sigma       Specify sigma threshold for the stddev output to file.\n");
    printf("                      The plot (shown with -prof) only shows 3 sigma detections.\n");
    printf("  -mod                Compute modulation index profile.\n");
    printf("  -mod_sigma          Specify sigma threshold for the mod. index output to file.\n");
    printf("                      The plot (shown with -prof) only shows 3 sigma detections.\n");
//START REGION DEVELOP
    printf("  -modSimple          Calulate modulation index directly from pulse stack,\n");
    printf("                      rather than via the LRFS.\n");
    printf("  -p3class            Use spectral analysis as explained in Smits et al. 2005,\n");
    printf("                      A&A, 440, 683. to do the mode classification. The output\n");
    printf("                      is P3 as function of pulse number. The -p3class option\n");
    printf("                      should be followed by \"Nblock Nmin Npad cpp1 cpp2\".\n");
    printf("                      Nblock is the initial size of the sequences being\n");
    printf("                      analysed, before the start/end point are optimised. Nmin\n");
    printf("                      is the minimum nr of pulses allowed in a single stretch of\n");
    printf("                      a mode. Only spectral peaks between cpp1 and cpp2 cycles\n");
    printf("                      per period are considered to be a drift-mode to be\n");
    printf("                      identified. Npad defines the minimum data length in pulses\n");
    printf("                      to be used while calculating spectra. If the stretch of\n");
    printf("                      data is shorter than Npad the data will be zero-padded to\n");
    printf("                      Npad pulses.\n");
//START REGION RELEASE
    printf("  -track              Compute subpulse phase track (use with -freq).\n");
    printf("  -track_dphase       Add specified offset (in deg) to the subpulse phase track.\n");
    printf("  -track_firstregion  Only use the first selected onpulse region to find the\n");
    printf("                      alignments of the phases of the different fft blocks.\n");
    printf("                      The other onpulse regions are still used to subtract from\n");
    printf("                      the LRFS from which the phases are derived.\n");
    printf("  -slope              Subtract slope from subpulse phases (in degrees subpulse\n");
    printf("                      phase per degree pulse longitude).\n");
    printf("  -amplitude          Compute modulation amplitude (use with -freq).\n");
//START REGION DEVELOP
    printf("  -trackphases        Compute the relative phases of the subpulse phase tracks.\n");
//START REGION RELEASE
    printf("  -2dfs               Compute 2DFS.\n");
    printf("  -s2dfs_p3           Compute S2DFS (sliding 2DFS P3 map).\n");
    printf("  -s2dfs_p2           Compute S2DFS (sliding 2DFS P2 map)\n");
    printf("                      (for first selected region only).\n");
//START REGION DEVELOP
    printf("  -hrfs               Compute the (folded) HRFS.\n");
    printf("  -hrfs_notfolded     Compute not folded HRFS.\n");
    printf("  -lrac               Compute LRAC map (auto-correlations of the intensities\n");
    printf("                      for each pulse longitude bin)\n");
    printf("  -lrac_leavewf       No attempt is made to remove the window function effect\n");
    printf("                      because of zero padding of the data.\n");
    printf("  -lrac_zap0          Set coefficients for lag zero to zero.\n");
    printf("  -lrac_keepneg       Store the negative lags (no extra information).\n");
    printf("  -lrcc               Compute LRCC map (with this specified lag, or this range\n");
    printf("                      in lags).\n");
    printf("  -lrcc_nosubmean     Don't subtract mean becore calculating the LRCC map.\n");
    printf("  -smod               Compute sliding modulation index map.\n");
    printf("  -strack             Compute STRACK map.\n");
    printf("  -strack_mask        Mask all bins where the modulation index is not\n");
    printf("                      significant (works only if -smod option is specified).\n");
    printf("  -strack_mask2       Like -strack_mask option, but set undefined values to -1 \n");
    printf("                      instead of 0\n");
//START REGION DEVELOP
    printf("  -strack_dphase      Add this subpulse phase in degrees.\n");
    printf("  -dphase_2ndfile     Add this subpulse phase in degrees to second file which\n");
    printf("                      is then also written out.\n");
//START REGION RELEASE
    printf("  -freq               Define which fluctuation frequencies (in cpp) are used for\n");
    printf("                      the subpulse phase track/amplitude calculation\n");
//START REGION DEVELOP
    printf("                      (or the inverse fft).\n");
//START REGION RELEASE
    printf("                      Can only be used once on the command-line.\n");
//START REGION DEVELOP
    printf("  -freq_mask          Use together with -freq. With -freq_mask specified all\n");
    printf("                      other frequencies are removed in the lrfs. Note this only\n");
    printf("                      affects 2dfs when -inversefft is specified. Can be useful\n");
    printf("                      to check if the -freq option is selecting the correct.\n");
    printf("                      frequencies. Without this option specified the lrfs is\n");
    printf("                      unaffected.\n");
    printf("  -inversefft         Inverse fft the data using lrfs. Use with -freq and\n");
    printf("                      -freq_mask to select a single frequency range. By using\n");
    printf("                      -p3zap more complicated things can be done.\n");
//START REGION RELEASE
    printf("  -p2zap              \"P2min P2max\" Zap fluctuations in this P2 range in cpp.\n");
    printf("  -p3zap              \"P3min P3max\" Zap fluctuations in this P3 range.\n");
    printf("                      P3min and P3max can be specified as bins or in cpp.\n");
    printf("\nGraphics options:\n");
    printf("  -onpulsed           Set pgplot device for the selection of the onpulse region.\n");
    printf("  -profd              Set pgplot device for the pulse profile.\n");
    printf("  -lrfsd              Set pgplot device for the LRFS.\n");
//START REGION DEVELOP
    printf("  -hrfsd              Set pgplot device for the LRCC.\n");
    printf("  -lracd              Set pgplot device for the LRAC.\n");
    printf("  -lrccd              Set pgplot device for the LRCC.\n");
//START REGION RELEASE
    printf("  -trackd             Set pgplot device for the subpulse phase.\n");
    printf("  -amplituded         Set pgplot device for the modulation amplitude.\n");
//START REGION DEVELOP
    printf("  -trackphasesd       Set pgplot device for the relative subpulse phases.\n");
//START REGION RELEASE
    printf("  -2dfsd              Set pgplot device for the 2DFS.\n");
    printf("  -s2dfs_p3d          Set pgplot device for the S2DFS (P3 map).\n");
    printf("  -s2dfs_p2d          Set pgplot device for the S2DFS (P2 map).\n");
//START REGION DEVELOP
    printf("  -smodd              Set pgplot device for the sliding modulation index map.\n");
    printf("  -strackd            Set pgplot device for the STRACK map.\n");
//START REGION DEVELOP
    printf("  -zoom               Use selected region to zoom in.\n");
    printf("  -zoomw              Use selected region to zoom in, but show some offpulse as\n");
    printf("                      well.\n");
    printf("  -zoom1              Use only first selected region to zoom in.\n");
    printf("  -p2zoom             Zoom in to this region in P2 space (in cpp).\n");
//START REGION DEVELOP
    printf("  -p3                 Subtract this p3 value from STRACK map.\n");
//START REGION RELEASE
    printf("  -onpulsegr          Enables graphical selection of additional on-pulse regions\n");
    printf("                      to those defined with the -onpulse option.\n");

    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about the lrfs/2dfs/modulation index can be found in:\n");
    printf(" - Weltevrede et al. 2006, A&A, 445, 243\n");
    printf(" - Weltevrede et al. 2007, A&A, 469, 607\n");
    printf("More information about bootstrap/subpulse phase track & amplitude can be found in:\n");
    printf(" - Weltevrede et al. 2012, MNRAS, 424, 843\n");
    printf("More information about the sliding 2dfs can be found in:\n");
    printf(" - Serylak et al. 2009, A&A, 506, 865\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else if(argc > 2 || strcmp(argv[argc-1], "-formatlist") == 0 || strcmp(argv[argc-1], "-headerlist") == 0) {
    int lastindex;
    lastindex = argc-1;
    if(strcmp(argv[argc-1], "-formatlist") == 0 || strcmp(argv[argc-1], "-headerlist") == 0)
      lastindex++;
    for(i = 1; i < lastindex; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-nfft") == 0 || strcmp(argv[i], "-fft_size") == 0 || strcmp(argv[i], "-fftsize") == 0 || strcmp(argv[i], "-fft_length") == 0 || strcmp(argv[i], "-fftlength") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &fft_size, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-bootstrap") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &bootstrap, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-onpulsed") == 0) {
	strcpy(onpulseselectdevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-profd") == 0) {
	strcpy(profiledevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-lrfsd") == 0) {
	strcpy(lrfsdevice, argv[i+1]);
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-lracd") == 0) {
	strcpy(lracdevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-hrfsd") == 0) {
	strcpy(hrfsdevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-lrccd") == 0) {
	strcpy(lrccdevice, argv[i+1]);
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-trackd") == 0) {
	strcpy(trackdevice, argv[i+1]);
	i++; 
      }else if(strcmp(argv[i], "-amplituded") == 0) {
	strcpy(amplitudedevice, argv[i+1]);
	i++;
//START REGION DEVELOP
     }else if(strcmp(argv[i], "-smodd") == 0) {
	strcpy(s_modmap_device, argv[i+1]);
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-2dfsd") == 0) {
	strcpy(twodfsdevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-s2dfs_p3d") == 0) {
	strcpy(s2dfs_p3_device, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-s2dfs_p2d") == 0) {
	strcpy(s2dfs_p2_device, argv[i+1]);
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-strackd") == 0) {
	strcpy(s_track_device, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-trackphasesd") == 0) {
	strcpy(trackphasesdevice, argv[i+1]);
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-w") == 0) {
	write_flag = 1;
      }else if(strcmp(argv[i], "-prof") == 0 || strcmp(argv[i], "-profile") == 0) {
	profile_flag = 1;
      }else if(strcmp(argv[i], "-stddev") == 0) {
	stddev_flag = 1;
      }else if(strcmp(argv[i], "-mod") == 0) {
	mod_flag = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-modSimple") == 0) {
	modSimple_flag = 1;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-lrfs") == 0) {
	lrfs_flag = 1;
      }else if(strcmp(argv[i], "-track") == 0) {
	track_flag = 1;
      }else if(strcmp(argv[i], "-track_firstregion") == 0) {
	track_only_first_region = 1;
      }else if(strcmp(argv[i], "-amplitude") == 0) {
	amplitude_flag = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-freq_mask") == 0) {
	ftrack_mask = 1;
      }else if(strcmp(argv[i], "-strack") == 0) {
	s_track_flag = 1;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-2dfs") == 0) {
	twodfs_flag = 1;
      }else if(strcmp(argv[i], "-s2dfs_p3") == 0) {
	s2dfs_p3_flag = 1;
      }else if(strcmp(argv[i], "-s2dfs_p2") == 0) {
	s2dfs_p2_flag = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-smod") == 0) {
	s_modmap_flag = 1;
      }else if(strcmp(argv[i], "-zoom") == 0) {
	zoom_flag = 1;
      }else if(strcmp(argv[i], "-zoomw") == 0) {
	zoom_flag = 2;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-p2zap") == 0) {
	i++;
      }else if(strcmp(argv[i], "-p3zap") == 0) {
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-zoom1") == 0) {
	if(zoom_flag == 0)
	  zoom_flag = 1;
	zoom_flag1 = 1;
      }else if(strcmp(argv[i], "-trackphases") == 0) {
	trackphases_flag = 1;
      }else if(strcmp(argv[i], "-inversefft") == 0) {
	inverseFFT = 1;
      }else if(strcmp(argv[i], "-strack_mask") == 0) {
	strack_mask = 1;
      }else if(strcmp(argv[i], "-strack_mask2") == 0) {
	strack_mask = 2;
//START REGION DEVELOP
//START REGION RELEASE
      }else if(strcmp(argv[i], "-powertwo") == 0) {
	powertwo = 1;
      }else if(strcmp(argv[i], "-DC") == 0) {
	subtractDC = 0;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-lrcc_nosubmean") == 0) {
	lrcc_nosubmean = 1;
//START REGION DEVELOP
//START REGION RELEASE
      }else if(strcmp(argv[i], "-onpulsegr") == 0) {
	selectMoreOnpulseRegions = 1; 
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-hrfs_notfolded") == 0) {
	hrfs_flag = 1;
      }else if(strcmp(argv[i], "-hrfs") == 0) {
	hrfs_flag = 2;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-lrcc") == 0) {
	int ret;
	ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%d %d", &lrcc_lag, &lrcc_lagstop, NULL);
	if(ret == 1) {
	  lrcc_lagstop = lrcc_lag;
	}else if(ret != 2) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option. Expected one or two values to be provided.", argv[i]);
	  return 0;
	}
	/*
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 1, -1, "%d %d", &lrcc_lag, &lrcc_lagstop, NULL) == 0) { // Reading two options failed
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 1, -1, "%d", &lrcc_lag, NULL) == 0) { // Reading one options failed as well
	    printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option. Expected one or two values to be provided.", argv[i]);
	    return 0;
	  }else {
	    if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &lrcc_lag, NULL) == 0) {
	      printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	      return 0;
	    }
	    lrcc_lagstop = lrcc_lag;
	  }
	}else {
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &lrcc_lag, &lrcc_lagstop, NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	}
	//	j = sscanf(argv[i+1], "%d %d", &lrcc_lag, &lrcc_lagstop);
	//	if(j != 2) {
	//	  j = sscanf(argv[i+1], "%d", &lrcc_lag);
	//	  if(j != 1) {
	//	    printerror(application.verbose_state.debug, "Error parsing %s option, need one or two values", argv[i]);
	//	    return 0;
	//	  }
	//	  lrcc_lagstop = lrcc_lag;
	//	}
	*/
	lrcc_flag = 1;
	i++;
      }else if(strcmp(argv[i], "-lrac") == 0) {
	lrac_flag = 1;
      }else if(strcmp(argv[i], "-lrac_leavewf") == 0) {
	lrac_flag = 2;
      }else if(strcmp(argv[i], "-lrac_zap0") == 0) {
	lrac_zap0 = 1;
      }else if(strcmp(argv[i], "-lrac_keepneg") == 0) {
	lrac_keepbothhalves = 1;
      }else if(strcmp(argv[i], "-p3") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &p3_subtract, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	p3_subtract_flag = 1;
	i++;
      }else if(strcmp(argv[i], "-p3class") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d %d %f %f", &p3class_nblock, &p3class_nmin, &p3class_npad, &p3class_cpp1, &p3class_cpp2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	p3class = 1;
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-track_dphase") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &track_dphase, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-strack_dphase") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &strack_dphase, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-dphase_2ndfile") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &strack_dphase2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	strack_2ndfile_flag = 1;
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-slope") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &slope, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-freq") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &freq_min, &freq_max, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-p2zoom") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &p2min, &p2max, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	p2range_set = 1;
        i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-mod_sigma") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &mod_sigma, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-stddev_sigma") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &stddev_sigma, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-C") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &(application.onpulse.left_bin[application.onpulse.nrRegions]), &(application.onpulse.right_bin[application.onpulse.nrRegions]), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	application.onpulse.bins_defined[application.onpulse.nrRegions] = 1;
	application.onpulse.nrRegions += 1;
	if(application.onpulse.nrRegions == MAX_pulselongitude_regions) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Too many regions selected.");
	  return 0;
	}
        i++;
//START REGION RELEASE
      }else {
        printerror(application.verbose_state.debug, "Unknown option: %s", argv[i]);
	terminateApplication(&application);
	return 0;
      }
    }
  }
//START REGION DEVELOP

  if(application.ppgplot_name != NULL) {
    if(pgopenoutputfile(application.ppgplot_name) == 0) {
      printerror(application.verbose_state.debug, "ERROR pspec: Cannot open %s", application.ppgplot_name);
      return 0;
    }
  }

//START REGION RELEASE
  junk_float = log(fft_size)/log(2);
  junk_int = junk_float;
  junk_float = pow(2, junk_int);
  if(fabs(junk_float-fft_size) > 0.1) {
//START REGION DEVELOP
    if(hrfs_flag) {
      printwarning(application.verbose_state.debug, "WARNING pspec: fft length could be a power of two for the hrfs, but this exception is not implemented.");
    }
//START REGION RELEASE
    printerror(application.verbose_state.debug, "ERROR pspec: fft length is not a power of two.");
    return 0;
  }

//START REGION RELEASE

  /* Open data and read in profile */
  //  printf("XXXXX cleaning fin[i]\n");
  for(i = 0; i < MaxNrPolarizations; i++)
    cleanPSRData(&fin[i], application.verbose_state);
  //  cleanPSRData(&fin[0], application.verbose_state);
  //    cleanPSRData(&fin[1], application.verbose_state);
  //    cleanPSRData(&fin[2], application.verbose_state);
  //    cleanPSRData(&fin[3], application.verbose_state);
  //    cleanPSRData(&fin[4], application.verbose_state);
    //  printf("XXXXX cleaning fin[i] done\n");
  if(application.iformat <= 0) {
    application.iformat = guessPSRData_format(argv[argc-1], 0, application.verbose_state);
    if(application.iformat == -2 || application.iformat == -3)
      return 0;
  }
  if(isValidPSRDATA_format(application.iformat) == 0) {
    printerror(application.verbose_state.debug, "ERROR pspec: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", argv[argc-1]);
    return 0;
  }
  /* Load data file into fin[0], the data will be split later into the separate polarizations */
  //  printf("XXXXXX Closing fin[0]\n");
  closePSRData(&fin[0], 0, application.verbose_state);  // When reading in only a data-set, make sure it is not allocated
  //  printf("XXXXXX Closing fin[0] done\n");
  //  printf("XXXXXX Opening data\n");
  if(!openPSRData(&fin[0], argv[argc-1], application.iformat, 0, 1, 0, application.verbose_state))
    return 0;
  //  printf("XXXXXX Opening data done: %p\n", fin[4].tsub_list);
  /* Search commandline for header parameters to overwrite header
     parameters. */
  if(PSRDataHeader_parse_commandline(&fin[0], argc, argv, application.verbose_state) == 0)
    return 0;
  for(i = 1; i < argc; i++) {
    if(strcmp(argv[i], "-header") == 0) {
      printwarning(application.verbose_state.debug, "WARNING pspec: If using the -header option, be aware it applied BEFORE the preprocessing.");
    }
  }
  if(preprocessApplication(&application, &fin[0]) == 0) {
    return 0;
  }
  //  printf("XXXXXX After preprocess: %p\n", fin[4].tsub_list);
  double period;
  int ret;
  ret = get_period(fin[0], 0, &period, application.verbose_state);
  if(ret == 2) {
    printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
    return 0;
  }
  if(period < 0.001) {
    printerror(application.verbose_state.debug, "ERROR pspec: The period does not appear to be set in the header. Consider using the -header option.");
    return 0;
  }
  if(get_tsamp(fin[0], 0, application.verbose_state) < 0.0000001 || get_tsamp(fin[0], 0, application.verbose_state) >= 10) {
    printerror(application.verbose_state.debug, "ERROR pspec: The sampling time does not appear to be set correctly in the header. Consider using the -header option.");
    return 0;
  }
  if(fin[0].isDebase == 0) {
    printerror(application.verbose_state.debug, "ERROR pspec: Baseline is not subtracted. Use pmod -debase first.\n");
    return 0;
  }else if(fin[0].isDebase != 1) {
    printwarning(application.verbose_state.debug, "WARNING pspec:  It is not known if baseline is already subtracted. Use pmod -debase first.\n");
  }
  if(check_baseline_subtracted(fin[0], application.verbose_state) == 0) {
    printwarning(application.verbose_state.debug, "WARNING pspec: Baseline does not appear to be subtracted. Use pmod -debase first.\n");
  }
  /* Convert the specified onpulse regions which were defined as
     fractions on the commandline. This is now possible because we
     know the number of bins. */
  region_frac_to_int(&(application.onpulse), fin[0].NrBins, 0);
  if(fin[0].NrPols > MaxNrPolarizations) {
    printerror(application.verbose_state.debug, "ERROR pspec: Maximum supported input parameters is exceeded.\n");
    return 0;
  }

  /* Split the data in individual polarization channels */
  originalNrPols = fin[0].NrPols;
  //  printf("XXXXXX Before selecting polarization 0: %p\n", fin[4].tsub_list);
  if(fin[0].NrPols > 0) {
    for(i = fin[0].NrPols-1; i >= 0; i--) {
      if(preprocess_polselect(fin[0], &clone, i, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Error selecting Stokes parameter %ld.", i);
	return 0;
      }
      swap_orig_clone(&fin[i], &clone, application.verbose_state); 
    }
  }
  //  printf("XXXXXX After selecting polarization 0: %p\n", fin[4].tsub_list);

  /* Copy header parameters, which is usefull when we write out the data. */
  cleanPSRData(&fout, application.verbose_state);
  copy_params_PSRData(fin[0], &fout, application.verbose_state);


  /* For test purposes, see if scaling modindexes etc works correctly. */
  /*
  fin.NrSubints = 512;
  fin.NrSubints = 512*10;
  for(i = 1; i < 10; i++) {
    for(j = 0; j < 512*fin.NrBins; j++) {
      fin.data[i*512*fin.NrBins+j] = fin.data[j];
    }
    }
  */

  /* Do selection of the onpulse/offpulse regions */
  profileI = (float *)malloc(fin[0].NrBins*sizeof(float));
  if(profileI == NULL) {
    printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
    return 0;
  }
  fft_blocks = fin[0].NrSubints/fft_size;
  junk_int = fft_blocks*fft_size;
  if(application.verbose_state.verbose && (lrfs_flag || stddev_flag || mod_flag || track_flag  || amplitude_flag || modSimple_flag || profile_flag || twodfs_flag || s2dfs_p3_flag || s2dfs_p2_flag
//START REGION DEVELOP
					   || s_track_flag || s_modmap_flag || hrfs_flag
//START REGION RELEASE
					   )) { // For instance for P3 folding this message is irrelevant
    printf("Only using %ld of the %ld pulses for the spectra being generated (%ld blocks with fft size %d).\n", junk_int, fin[0].NrSubints, fft_blocks, fft_size);
  }
  if(junk_int == 0) {
    printerror(application.verbose_state.debug, "ERROR pspec: Not enough pulses, try a shorter fft length.");
    return 0;
  }
  read_partprofilePSRData(fin[0], profileI, NULL, 0, 0, fft_blocks*fft_size, noverbose);

  xmin = 0;
  ret = get_period(fin[0], 0, &period, application.verbose_state);
  if(ret == 2) {
    printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
    return 0;
  }
  xmax = 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period;
  if(application.onpulse.nrRegions == 0 || selectMoreOnpulseRegions) {
    strcpy(pgplot_options.viewport.plotDevice, onpulseselectdevice);
    strcpy(pgplot_options.box.xlabel, "Bin");
    strcpy(pgplot_options.box.ylabel, "Intensity");
    strcpy(pgplot_options.box.title, "Select on-pulse region");
    selectRegions(profileI, fin[0].NrBins, &pgplot_options, 0, powertwo, 1, &application.onpulse, application.verbose_state);
  }else {
    if(strcmp(profiledevice, "?") == 0)
      printf("Specify plotting device to show the profile showing the selected regions: \n  ");
    strcpy(pgplot_options.viewport.plotDevice, profiledevice);
    strcpy(pgplot_options.box.xlabel, "Phase[deg]");
    strcpy(pgplot_options.box.ylabel, "Intensity");
    strcpy(pgplot_options.box.title, fin[0].psrname);
    if(pgplotGraph1(&pgplot_options, profileI, NULL, NULL, fin[0].NrBins, xmin, xmax, 0, xmin, xmax, 0, 0, 0, 1, 0, 1, 0, 1, 1, &application.onpulse, -1, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR pspec: Unable to open plotdevice.\n");
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
	  printerror(application.verbose_state.debug, "ERROR pspec: region not defined in bins");
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

  stddev_av = modindex_av = stddev_square = modindex_square = NULL;  // To avoid compiler warnings

  if(lrfs_flag || stddev_flag || mod_flag || track_flag  || amplitude_flag || modSimple_flag || profile_flag) {
    lrfs = (float *)malloc(originalNrPols*(fft_size/2+1)*fin[0].NrBins*sizeof(float));
    stddev = (float *)calloc(fin[0].NrBins, sizeof(float));
    rms_sigma = (float *)calloc(fin[0].NrBins, sizeof(float));
    modindex = (float *)calloc(fin[0].NrBins, sizeof(float));
    rms_modindex = (float *)malloc(fin[0].NrBins*sizeof(float));
    phase_track = (float *)malloc(fin[0].NrBins*(bootstrap+2)*sizeof(float));   /* One for the actual phase track, the other for the stddev. */
    phase_track_phases = (float *)malloc((fin[0].NrSubints/fft_size)*sizeof(float));
    amplitude_profile = (float *)malloc(fin[0].NrBins*sizeof(float));
    if(lrfs == NULL  || stddev == NULL || rms_sigma == NULL || modindex == NULL || rms_modindex == NULL || phase_track == NULL || phase_track_phases == NULL || amplitude_profile == NULL) {
      printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
      return 0;
    }
    if(bootstrap > 0) {
      stddev_av = (double *)calloc(fin[0].NrBins, sizeof(double));
      modindex_av = (double *)calloc(fin[0].NrBins, sizeof(double));
      stddev_square = (double *)calloc(fin[0].NrBins, sizeof(double));
      modindex_square = (double *)calloc(fin[0].NrBins, sizeof(double));
      clone_profileI = (float *)malloc(fin[0].NrBins*sizeof(float));
      if(stddev_av  == NULL || modindex_av  == NULL || stddev_square  == NULL || modindex_square == NULL || clone_profileI == NULL) {
	printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
	return 0;
      }
      printf("Starting bootstrap\n");

      if(application.onpulse.nrRegions == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Cannot bootstrap without a selected region");
	return 0;
      }
      nrpointsrms = 0;
      rms = 0;
      avrg = 0;
      for(i = 0; i < fin[0].NrSubints; i++) {
	for(j = 0; j < fin[0].NrBins; j++) {
	  if(checkRegions(j, &application.onpulse, 0, application.verbose_state) == 0) {
	    nrpointsrms++;
	    if(readPulsePSRData(&fin[0], i, 0, 0, j, 1, &sampleI, application.verbose_state) != 1) {
	      printerror(application.verbose_state.debug, "ERROR pspec: Error reading data");
	      return 0;
	    }
	    rms += sampleI*sampleI;
	    avrg += sampleI;
	  }
	}
      }
      rms /= (double)nrpointsrms;
      avrg /= (double)nrpointsrms;
      rms -= avrg*avrg; 
      rms = sqrt(rms);
      if(application.verbose_state.verbose)
	printf("  Average off-pulse intensity = %e, rms = %e based on %ld points\n", avrg, rms, nrpointsrms);
      for(j = 0; j < fin[0].NrBins; j++) {
	modindex_av[j]  = 0;
	modindex_square[j]  = 0;
	stddev_av[j]  = 0;
	stddev_square[j] = 0;
      }
      for(i = 0; i < bootstrap; i++) {
	if(application.verbose_state.verbose && application.verbose_state.nocounters == 0) {
	  printf("\r  bootstrap step %ld/%d         ", i+1, bootstrap);
	  fflush(stdout);
	}
	//	cleanPSRData(&clone, application.verbose_state);
	if(preprocess_addNoise(fin[0], &clone, rms, noverbose) != 1) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Adding noise failed");
	  return 0;
	}
	/* Store the individual subpulse phase tracks, leave first tracks empty to fill with the calculation from the real data later and the stddev. */
	verbose_definition verbose2;
	copyVerboseState(application.verbose_state, &verbose2);
	if(application.verbose_state.verbose && (!application.verbose_state.nocounters))
	  verbose2.verbose = 1;
	else
	  verbose2.verbose = 0;
	if(calcLRFS(clone.data, fin[0].NrSubints, fin[0].NrBins, fft_size, lrfs, subtractDC, NULL, &phase_track[(i+2)*fin[0].NrBins], NULL, track_flag, freq_min, freq_max, track_only_first_region, NULL, 0, 0, 0, &application.onpulse, &var_rms, argc, argv, verbose2) == 0) { 
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate LRFS");
	  return 0;
	}

	read_partprofilePSRData(clone, clone_profileI, NULL, 0, 0, fft_blocks*fft_size, noverbose);
	// Disablse verbose in bootstrap itterations, unless debug is set
	int oldverbose = application.verbose_state.verbose;
	if(application.verbose_state.debug == 0) {
	  application.verbose_state.verbose = 0;
	}
	calcModindex(lrfs, clone_profileI, fin[0].NrBins, fft_size, fft_blocks*fft_size, stddev, rms_sigma, modindex, rms_modindex, &application.onpulse, var_rms, NULL, application.verbose_state);
	application.verbose_state.verbose = oldverbose;
//START REGION DEVELOP
	if(modSimple_flag)
	  calcModindexSimple(clone.data, fin[0].NrBins, fin[0].NrSubints, stddev, modindex, &application.onpulse, application.verbose_state);
//START REGION RELEASE
	for(j = 0; j < fin[0].NrBins; j++) {
	  modindex_av[j] += modindex[j];
	  modindex_square[j] += modindex[j]*modindex[j];
	  stddev_av[j] += stddev[j];
	  stddev_square[j] += stddev[j]*stddev[j];
	}
	clone.opened_flag = 1;   /* Required to free memory */
	closePSRData(&clone, 0, application.verbose_state);
      }
      printf("Bootstrap finished\n");
    }

    float avrg_offpulse_lrfs_power;
    for(i = 0; i < originalNrPols; i++) {
      int track_flag_pol, amplitude_flag_pol;
      track_flag_pol = track_flag;
      amplitude_flag_pol = amplitude_flag;
      if(i != 0) {
	track_flag_pol = 0;  // To make sure only first polarization channel is used for the subpulse phase
	amplitude_flag_pol = 0;
      }
      float *float_ptr;
      if(i == 0) {
	float_ptr = &avrg_offpulse_lrfs_power;
      }else {
	float_ptr = NULL;
      }
      if(calcLRFS(fin[i].data, fin[i].NrSubints, fin[i].NrBins, fft_size, &lrfs[i*(fft_size/2+1)*fin[0].NrBins], subtractDC, float_ptr, phase_track, phase_track_phases, track_flag_pol, freq_min, freq_max, track_only_first_region, amplitude_profile, amplitude_flag_pol, ftrack_mask, inverseFFT, &application.onpulse, &var_rms, argc, argv, application.verbose_state) == 0) {
	/* Plot only plot half the spectrum */
	printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate LRFS");
	return 0;
      }
    }
  
    /* Always calculate standard deviation and modulation index, even if only
       the profile is requested, because it always should be written
       out to the file to make output consistent with driftFigure */
    if(profile_flag || stddev_flag || mod_flag || modSimple_flag) {
//START REGION DEVELOP
      /* I tied -1 instead of var_rms. Also change later!!!! */
//START REGION RELEASE
//	avrg_offpulse_lrfs_power NEEDS TO BE PASSED ON TO calcModindex to report offpulse sigma

      calcModindex(lrfs, profileI, fin[0].NrBins, fft_size, fft_blocks*fft_size, stddev, rms_sigma, modindex, rms_modindex, &application.onpulse, var_rms, &avrg_offpulse_lrfs_power, application.verbose_state);

//START REGION DEVELOP
      if(modSimple_flag) {
	calcModindexSimple(fin[0].data, fin[0].NrBins, fin[0].NrSubints, stddev, modindex, &application.onpulse, application.verbose_state);
	for(j = 0; j < fin[0].NrBins; j++) {
	  rms_sigma[j] = 0;
	  rms_modindex[j] = 0;
	}
      }
      /*      printf("RMS standard deviation = %f\n", rms_sigma); Doesn't work anymore, is array. */
//START REGION RELEASE
      if(bootstrap > 0) {
	for(j = 0; j < fin[0].NrBins; j++) {
	  stddev_av[j] /= (double)bootstrap;
	  stddev_square[j] /= (double)bootstrap;
	  modindex_av[j] /= (double)bootstrap;
	  modindex_square[j] /= (double)bootstrap;
	  rms_sigma[j] = sqrt(stddev_square[j] - stddev_av[j]*stddev_av[j]);
	  rms_modindex[j] = sqrt(modindex_square[j] - modindex_av[j]*modindex_av[j]);
	  /* Better to take the average to be the best value. It turns
	     out that the modindex of the actual pulse stack is very
	     high, but that it is normal in the majority of the monte
	     carlo realizations (in noise with basically zero
	     intensity). The result can be a very large modulation
	     index with a relatively small error bar, making some
	     noise modulation indices significant. */
	  modindex[j] = modindex_av[j];
	  stddev[j] = stddev_av[j];
	}
      }
    }

    if(profile_flag) {
      if(strcmp(profiledevice, "?") == 0)
	printf("Specify plotting device to show the profile: \n  ");
      ret = get_period(fin[0], 0, &period, application.verbose_state);
      if(ret == 2) {
	printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
	return 0;
      }
      if(pgplotProfile(profiledevice, -1, -1, profileI, stddev, rms_sigma, modindex, rms_modindex, fin[0].NrBins, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period, "Phase[deg]", "Intensity", fin[0].psrname, stddev_flag, mod_flag|modSimple_flag, zoom_flag, xmin_zoom, xmax_zoom, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Unable to open plotdevice.\n");
	return 0;
      }
      if(write_flag) {
	 if(change_filename_extension(argv[argc-1], outputname, "profile", 1000, application.verbose_state) == 0) 
	   return 0;
	 fout_ascii = fopen(outputname, "w");
	 if(fout_ascii == NULL) {
	   printerror(application.verbose_state.debug, "ERROR pspec: Unable to open %s.", outputname);
	   return 0;
	 }
	 if(application.verbose_state.verbose) {
	   printf("Writing bin nr, intensity, stddev, rms stddev, modindex, rms modindex to %s\n", outputname);
	 }
	 for(i = 0; i < fin[0].NrBins; i++) {
	   fprintf(fout_ascii, "%ld %e", i, profileI[i]);
	   /*	   if(stddev_flag) */
	   if(stddev_sigma <= 0 || stddev[i]/rms_sigma[i] >= stddev_sigma)
	     fprintf(fout_ascii, " %e %e", stddev[i], rms_sigma[i]);
	   else
	     fprintf(fout_ascii, " %e %e", -1.0, 0.0);
	   /*	   if(mod_flag) */
	   if(mod_sigma <= 0 || modindex[i]/rms_modindex[i] >= mod_sigma)
	     fprintf(fout_ascii, " %e %e", modindex[i], rms_modindex[i]);
	   else
	     fprintf(fout_ascii, " %e %e", -1.0, 0.0);
	   fprintf(fout_ascii, "\n");
	 }
	 fclose(fout_ascii);
       }
    }
    if(lrfs_flag) {
      if(strcmp(lrfsdevice, "?") == 0)
	printf("Specify plotting device to show the LRFS: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, lrfsdevice);
      strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.ylabel, "P3 [cpp]");
      strcpy(pgplot_options.box.title, "LRFS");
  /*      pgplotMap(&pgplot_options, lrfs, fin[0].NrBins, 1+fft_size/2, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/get_period(fin[0], 0, application.verbose_state), xmin_zoom, xmax_zoom, 0, 0.5, 0, 0.5, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, "Pulse phase [degrees]", "P3 [cpp]", "LRFS", 1, 1, 1, -0.5, -0.5, 1.0, 1.0, 1, 1, 1, 0, 0, 0, 0, 0, application.verbose_state); */
      ret = get_period(fin[0], 0, &period, application.verbose_state);
      if(ret == 2) {
	printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
	return 0;
      }
      pgplotMap(&pgplot_options, lrfs, fin[0].NrBins, 1+fft_size/2, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period, xmin_zoom, xmax_zoom, 0, 0.5, 0, 0.5, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, application.verbose_state); 
      if(write_flag) {
	fout.NrSubints = 1+fft_size/2;
	fout.NrBins = fin[0].NrBins;
	fout.NrPols = originalNrPols;
	fout.gentype = GENTYPE_LRFS;
	fout.yrangeset = 1;
	fout.yrange[0] = 0;
	fout.yrange[1] = 0.5;
	fout.tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
	if(fout.tsub_list != NULL) {
	  free(fout.tsub_list);
	}
	fout.tsub_list = (double *)malloc(sizeof(double));
	if(fout.tsub_list == NULL) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
	  return 0;
	}
	fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
	//	fprintf(stderr, "XXXXX %lf\n", fout.fixedtsub);

	if(change_filename_extension(argv[argc-1], outputname, "lrfs", 1000, application.verbose_state) == 0) {
	  return 0;
	}
	if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	  return 0;
	if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
	  return 0;
	}
	//	appendHistoryLine(fout, argc, argv, application.verbose_state);

	/* Transpose the data for writing */
	lrfs2 = (float *)malloc(originalNrPols*(fft_size/2+1)*fin[0].NrBins*sizeof(float));
	if(lrfs2 == NULL) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
	  return 0;
	}

	for(i = 0; i < fft_size/2+1; i++) {
	  for(p = 0; p < originalNrPols; p++) {
	    memcpy(&lrfs2[(originalNrPols*i+p)*fin[0].NrBins], &lrfs[(p*(fft_size/2+1)+i)*fin[0].NrBins], fin[0].NrBins*sizeof(float));
	  }
	}
	if(writePSRData(&fout, lrfs2, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
	  return 0;
	}
	free(lrfs2);
	closePSRData(&fout, 1, application.verbose_state); // Keep history, as header can be written out again
	fout.gentype = GENTYPE_UNDEFINED;
	fout.yrangeset = 0;
	fout.xrangeset = 0;
	fout.NrPols = 1;
      }
    }

    if(track_flag) {
      for(j = 0; j < fin[0].NrBins; j++) {
	phase_track[fin[0].NrBins+j] = -1;   /* This will become the stddev, -1 if not calculated */
      }
      if(bootstrap > 0) {
	for(k = 0; k < 100; k++) {
	  for(j = 0; j < fin[0].NrBins; j++) {
	    phase_track[fin[0].NrBins+j] = 0;   /* This will become the stddev, has to be set to zero for each itteration */
	  }
	  for(i = 0; i < bootstrap; i++) {
	    for(j = 0; j < fin[0].NrBins; j++) {
	      float xdeg;
	      xdeg = phase_track[j] -  phase_track[(i+2)*fin[0].NrBins+j];
	      if(fabs(xdeg+360) < fabs(xdeg))
		xdeg += 360;
	      if(fabs(xdeg-360) < fabs(xdeg))
		xdeg -= 360;
	      phase_track[fin[0].NrBins+j] += xdeg*xdeg;
	    }
	  }
	  for(j = 0; j < fin[0].NrBins; j++) {
	    phase_track[fin[0].NrBins+j] = sqrt(phase_track[fin[0].NrBins+j]/(float)bootstrap);
	  }
	  /* Now we estimated the error bar, we can try to take out the
	     overall offset between the different bootstrap
	     itterations. They can have different subpulse phase offsets
	     (same offset for all bins) with respect of each other,
	     resulting in too big estimated error bars. So calculate and
	     correct the weighted offset for each itteration. */
	  float offset, weight;
	  for(i = 0; i < bootstrap; i++) {
	    offset = 0;
	    weight = 0;
	    for(j = 0; j < fin[0].NrBins; j++) {
	      if(checkRegions(j, &application.onpulse, 0, application.verbose_state) != 0) {
		float xdeg;
		xdeg = phase_track[(i+2)*fin[0].NrBins+j] - phase_track[j];
		if(fabs(xdeg+360) < fabs(xdeg))
		  xdeg += 360;
		if(fabs(xdeg-360) < fabs(xdeg))
		  xdeg -= 360;
		offset += (xdeg)/phase_track[fin[0].NrBins+j];
		weight += 1.0/phase_track[fin[0].NrBins+j];
	      }
	    }
	    offset /= weight;
	    if(application.verbose_state.verbose && application.verbose_state.nocounters == 0) {
	      printf("Itteration %9ld: Offset track %9ld: %20e\n", k+1, i, offset);
	    }
	    for(j = 0; j < fin[0].NrBins; j++) {
	      phase_track[(i+2)*fin[0].NrBins+j] -= offset; 
	      if(phase_track[(i+2)*fin[0].NrBins+j] < 0)
		phase_track[(i+2)*fin[0].NrBins+j] += 360;
	      if(phase_track[(i+2)*fin[0].NrBins+j] >= 360)
		phase_track[(i+2)*fin[0].NrBins+j] -= 360;
	    }
	  }	
	}
      }

      for(i = 0; i < fin[0].NrBins; i++) {
	float xdeg;
	xdeg = (xmax-xmin)*i/(float)(fin[0].NrBins-1);
	phase_track[i] -= slope*xdeg - track_dphase;
	phase_track[i] = derotate_deg(phase_track[i]);
      }
      if(bootstrap > 0) {
	for(j = 0; j < bootstrap; j++) {
	  for(i = 0; i < fin[0].NrBins; i++) {
	    float xdeg;
	    xdeg = (xmax-xmin)*i/(float)(fin[0].NrBins-1);
	    phase_track[(j+2)*fin[0].NrBins+i] -= slope*xdeg - track_dphase;
	    phase_track[(j+2)*fin[0].NrBins+i] = derotate_deg(phase_track[(j+2)*fin[0].NrBins+i]);
	  }
	}
      }

      if(strcmp(trackdevice, "?") == 0)
	printf("Specify plotting device to show the subpulse phase track: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, trackdevice);
      strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.ylabel, "Subpulse phase");
      strcpy(pgplot_options.box.title, "Subpulse phase track");
      if(bootstrap > 0) {
	pgplotGraph1(&pgplot_options, phase_track, NULL, &phase_track[fin[0].NrBins], fin[0].NrBins, xmin, xmax, 0, xmin_zoom, xmax_zoom, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
      }else {
	pgplotGraph1(&pgplot_options, phase_track, NULL, NULL, fin[0].NrBins, xmin, xmax, 0, xmin_zoom, xmax_zoom, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
      }
      if(write_flag) {
	 if(change_filename_extension(argv[argc-1], outputname, "track", 1000, application.verbose_state) == 0) 
	   return 0;
	 fout_ascii = fopen(outputname, "w");
	 if(fout_ascii == NULL) {
	   printerror(application.verbose_state.debug, "ERROR pspec: Unable to open %s.", outputname);
	   return 0;
	 }
	 for(i = 0; i < fin[0].NrBins; i++) {
	   /* fprintf(fout_ascii, "%ld %f %f %f %f %f\n", i, phase_track[i], phase_track[fin[0].NrBins+i], phase_track[2*fin[0].NrBins+i], phase_track[3*fin[0].NrBins+i], phase_track[4*fin[0].NrBins+i]); */
	   fprintf(fout_ascii, "%ld %f %e\n", i, phase_track[i], phase_track[fin[0].NrBins+i]); 
	 }
	 fclose(fout_ascii);
      }
//START REGION DEVELOP
      if(trackphases_flag) {
	if(strcmp(trackphasesdevice, "?") == 0)
	  printf("Specify plotting device to show the subpulse phase track phases: \n  ");
	strcpy(pgplot_options.viewport.plotDevice, trackphasesdevice);
	strcpy(pgplot_options.box.xlabel, "FFT block number");
	strcpy(pgplot_options.box.ylabel, "Subpulse phase track phase offset");
	strcpy(pgplot_options.box.title, "Subpulse phase track phase offsets");
	pgplotGraph1(&pgplot_options, phase_track_phases, NULL, NULL, fin[0].NrSubints/fft_size, 0, fin[0].NrSubints/fft_size - 1, 0, 0, fin[0].NrSubints/fft_size - 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
	if(write_flag) {
	  if(change_filename_extension(argv[argc-1], outputname, "trackphases", 1000, application.verbose_state) == 0) 
	    return 0;
	  fout_ascii = fopen(outputname, "w");
	  if(fout_ascii == NULL) {
	    printerror(application.verbose_state.debug, "ERROR pspec: Unable to open %s.", outputname);
	    return 0;
	  }
	  for(i = 0; i < fin[0].NrSubints/fft_size; i++) {
	    fprintf(fout_ascii, "%ld %f\n", i*fft_size, phase_track_phases[i]);
	  }
	  fclose(fout_ascii);
	}
      }
//START REGION RELEASE
    }  // End of if(track_flag)

    if(amplitude_flag) {
      if(strcmp(amplitudedevice, "?") == 0)
	printf("Specify plotting device to show the subpulse amplitude profile: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, amplitudedevice);
      pgplot_options.viewport.dontclose = 1;
      strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.ylabel, "Subpulse amplitude");
      strcpy(pgplot_options.box.title, "Subpulse amplitude profile");
      if(!(profile_flag || stddev_flag || mod_flag || modSimple_flag)) {  // If didn't calculate normalised profile before, normalise it now
	float imax;
	imax = 1;
	for(i = 0; i < fin[0].NrBins; i++) {
	  if(i == 0 || profileI[i] > imax)
	    imax = profileI[i];
	}
	for(i = 0; i < fin[0].NrBins; i++) {
	  profileI[i] /= imax;
	}
      }
      pgplotGraph1(&pgplot_options, profileI, NULL, NULL, fin[0].NrBins, xmin, xmax, 0, xmin_zoom, xmax_zoom, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
      /*pgplot_options.viewport.noclear=1; pgplot_options.viewport.dontclose = 1; pgplot_options.viewport.dontopen = 1; pgplotGraph1(&amplitudedevice, viewport, amplitude_profile, NULL, NULL, fin[0].NrBins, xmin, xmax, xmin_zoom, xmax_zoom, 0, "", "", "", 0, 1, 0, 2, 1, NULL, -1, application.verbose_state); */
      ppgsci(2);
      for(i = 0; i < fin[0].NrBins; i++) {
	float x;
	x = (xmax-xmin)*i/(float)(fin[0].NrBins-1);
	if(i == 0) {
	  ppgmove(x, amplitude_profile[i]);
	}else {
	  ppgdraw(x, amplitude_profile[i]);
	}
      }
      ppgsci(1);
      ppgclos();
      if(write_flag) {
	if(change_filename_extension(argv[argc-1], outputname, "amplitude", 1000, application.verbose_state) == 0) 
	  return 0;
	fout_ascii = fopen(outputname, "w");
	if(fout_ascii == NULL) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to open %s.", outputname);
	  return 0;
	}
	for(i = 0; i < fin[0].NrBins; i++) {
	  fprintf(fout_ascii, "%ld %f\n", i, amplitude_profile[i]);
	}
	fclose(fout_ascii);
      }
    }  // End of if(amplitude_flag)

//START REGION DEVELOP

    if(inverseFFT && write_flag) {
      fout.NrSubints = fin[0].NrSubints/(fft_size);
      fout.NrSubints *= fft_size;
      fout.NrBins = fin[0].NrBins;
      if(change_filename_extension(argv[argc-1], outputname, "inversefft", 1000, application.verbose_state) == 0) {
	return 0;
      }
      if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	return 0;
      if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
	return 0;
      }
      //      appendHistoryLine(fout, argc, argv, application.verbose_state);
      if(writePSRData(&fout, fin[0].data, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
	return 0;
      }
      closePSRData(&fout, 1, application.verbose_state); // Keep history, as header can be written out again
    }

//START REGION RELEASE
    free(lrfs);
    free(stddev);
    free(modindex);
    free(rms_sigma);
    free(rms_modindex);
    free(phase_track);
    free(phase_track_phases);
    free(amplitude_profile);
//START REGION DEVELOP
  }else if (ftrack_mask || inverseFFT) {
    printerror(application.verbose_state.debug, "To make the -freq_mask and -inversefft option do something, please specify the -lrfs, -stddev, -mod or -track option as well.");
    return 0;
//START REGION RELEASE
  }

  if(twodfs_flag) {
    for(regionnr = 0; regionnr < application.onpulse.nrRegions; regionnr++) {
      if(application.onpulse.bins_defined[regionnr] == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: region not defined in bins");
	return 0;
      }
      twodfs = (float *)malloc((1+fft_size/2)*(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1)*sizeof(float));
      if(twodfs == NULL) {
	printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
	return 0;
      }
      if(calc2DFS(fin[0].data, fin[0].NrSubints, fin[0].NrBins, fft_size, twodfs, &application.onpulse, regionnr, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate 2DFS");
	return 0;
      }
      for(i = 1; i < argc-1; i++) {
	if(strcmp(argv[i], "-p3zap") == 0) {
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &zapmin, &zapmax, NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	  //	  j = sscanf(argv[i+1], "%f %f", &zapmin, &zapmax);
	  //	  if(j != 2) {
	  //	    printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse -p3zap option, specify two values.");
	  //	    return 0;
	  //	  }
	  for(j = 0; j < (1+fft_size/2); j++) {
	    if(zapmin > 0.9 || zapmax > 0.9) {  // If the boundary exceeds 0.5 cpp, then the input parameters are assumed to be in bins
	      p3 = j; 
	    }else {
	      p3 = j/(float)fft_size;   /* not sure if exactly right, but anyway */
	    }
	    if(p3 >= zapmin && p3 <= zapmax) {
	      for(k = application.onpulse.left_bin[regionnr]; k <= application.onpulse.right_bin[regionnr]; k++) {
		twodfs[j*(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1)+(k-application.onpulse.left_bin[regionnr])] = 0;
	      }
	    }
	  }
	}

	if(strcmp(argv[i], "-p2zap") == 0) {
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &zapmin, &zapmax, NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	  //	  j = sscanf(argv[i+1], "%f %f", &zapmin, &zapmax);
	  //	  if(j != 2) {
	  //	    printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse -p2zap option, specify two values.");
	  //	    return 0;
	  //	  }
	  for(j = 0; j < (application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1); j++) {
	    p2 = j*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1)-fin[0].NrBins/2-0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1);
	    if(p2 >= zapmin && p2 <= zapmax) {
	      /*	      printf("XXXXXX: %f %ld: %e %e\n", p2, j, zapmin, zapmax); */
	      for(k = 0; k < (1+fft_size/2); k++) {
		twodfs[k*(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1)+(j)] = 0;
	      }
	    }
	  }
	}

      }
      if(p2range_set == 0) {
	p2min = -fin[0].NrBins/2-0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1);
	p2max = +fin[0].NrBins/2-0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1);
      }
      if(strcmp(twodfsdevice, "?") == 0)
	printf("Specify plotting device to show the 2DFS: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, twodfsdevice);
      strcpy(pgplot_options.box.xlabel, "P2 [cpp]");
      strcpy(pgplot_options.box.ylabel, "P3 [cpp]");
      strcpy(pgplot_options.box.title, "2DFS");
      /*      pgplotMap(&pgplot_options, twodfs, application.onpulse.right[regionnr]-application.onpulse.left[regionnr]+1, fft_size/2+1, -fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right[regionnr]-application.onpulse.left[regionnr]+1), fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right[regionnr]-application.onpulse.left[regionnr]+1), p2min, p2max, 0, 0.5, 0, 0.5, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, "P2 [cpp]", "P3 [cpp]", "2DFS", 1, 1, 1, -0.5, -0.5, 1.0, 1.0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state); */
      pgplotMap(&pgplot_options, twodfs, application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1, fft_size/2+1, -fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1), fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1), p2min, p2max, 0, 0.5, 0, 0.5, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state); 
      if(write_flag) {
	fout.NrSubints = fft_size/2+1;
	fout.NrBins = application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1;
	fout.gentype = GENTYPE_2DFS;
	fout.tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
	if(fout.tsub_list != NULL)
	  free(fout.tsub_list);
	fout.tsub_list = (double *)malloc(sizeof(double));
	if(fout.tsub_list == NULL) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
	  return 0;
	}
	fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
	fout.yrangeset = 1;
	fout.yrange[0] = 0;
	fout.yrange[1] = 0.5;
	fout.xrangeset = 1;
	fout.xrange[0] = -fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1);
	fout.xrange[1] = fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[regionnr]-application.onpulse.left_bin[regionnr]+1);

	sprintf(txt, "%d.2dfs", regionnr+1);
	if(change_filename_extension(argv[argc-1], outputname, txt, 1000, application.verbose_state) == 0) {
	  return 0;
	}
	if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	  return 0;
	if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
	  return 0;
	}
	//	appendHistoryLine(fout, argc, argv, application.verbose_state);
	if(writePSRData(&fout, twodfs, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
	  return 0;
	}
	closePSRData(&fout, 1, application.verbose_state); // Keep history, as header can be written out again
	fout.gentype = GENTYPE_UNDEFINED;
	fout.yrangeset = 0;
	fout.xrangeset = 0;
      }
      free(twodfs);
    }
  }   // End of (twodfs_flag)

//START REGION RELEASE

  if(s2dfs_p3_flag || s2dfs_p2_flag) {
    float *s2dfs_p3, *s2dfs_p2;
    printf("Calculating S2DFS\n");
    if(application.onpulse.nrRegions < 1) {
      printerror(application.verbose_state.debug, "ERROR pspec: region is not defined");
      return 0;
    }
    if(application.onpulse.bins_defined[0] == 0) {
      printerror(application.verbose_state.debug, "ERROR pspec: region not defined in bins");
      return 0;
    }

    twodfs   = (float *)malloc((1+fft_size/2)*(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1)*sizeof(float));
    s2dfs_p3 = (float *)malloc((1+fft_size/2)*(fin[0].NrSubints-fft_size+1)*sizeof(float));
    s2dfs_p2 = (float *)malloc((application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1)*(fin[0].NrSubints-fft_size+1)*sizeof(float));
    if(twodfs == NULL || s2dfs_p3 == NULL || s2dfs_p2 == NULL) {
      printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
      return 0;
    }
    verbose_definition noverbose;
    copyVerboseState(application.verbose_state, &noverbose);
    noverbose.verbose = 0;
    noverbose.nocounters = 1;
    /*    fin[0].NrSubints = fft_size + 300;  */
    for(i = 0; i < fin[0].NrSubints-fft_size+1; i++) {
      if(calc2DFS(&fin[0].data[i*fin[0].NrBins], fft_size, fin[0].NrBins, fft_size, twodfs, &application.onpulse, 0, noverbose) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate 2DFS");
	return 0;
      }
      for(l = 1; l < argc-1; l++) {
	if(strcmp(argv[l], "-p3zap") == 0) {
	  if(parse_command_string(application.verbose_state, argc, argv, l+1, 0, -1, "%f %f", &zapmin, &zapmax, NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	  //	  j = sscanf(argv[l+1], "%f %f", &zapmin, &zapmax);
	  //	  if(j != 2) {
	  //	    printerror(application.verbose_state.debug, "ERROR pspec: Cannot parse -p3zap option.");
	  //	    return 0;
	  //	  }
	  for(j = 0; j < (1+fft_size/2); j++) {
	    p3 = j/(float)fft_size;   // not sure if exactly right, but anyway 
	    if(p3 >= zapmin && p3 <= zapmax) {
	      for(k = 0; k < (application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1); k++) {
		twodfs[j*(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1)+k] = 0;
	      }
	    }
	  }
	}
      }
      if(s2dfs_p3_flag) {
	/* Collapse spectrum over P2 */
	for(j = 0; j < (1+fft_size/2); j++) {
	  s2dfs_p3[j*(fin[0].NrSubints-fft_size+1)+i] = 0;
	  for(k = 0; k < (application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1); k++) {
	    s2dfs_p3[j*(fin[0].NrSubints-fft_size+1)+i] += twodfs[j*(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1)+k];
	  }
	}
      }
      if(s2dfs_p2_flag) {
	/* Collapse spectrum over P3 */
	for(k = 0; k < (application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1); k++) {
	  s2dfs_p2[k*(fin[0].NrSubints-fft_size+1)+i] = 0;
	  for(j = 0; j < (1+fft_size/2); j++) {
	    s2dfs_p2[k*(fin[0].NrSubints-fft_size+1)+i] += twodfs[j*(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1)+k];
	  }
	}
      }
      if(!application.verbose_state.nocounters)
	printf("  Block %ld of the %ld  (%.2f%%)      \r", i, fin[0].NrSubints-fft_size+1, 100*(i+1)/(float)(fin[0].NrSubints-fft_size+1));
    }
    printf("  done                                             \n");
    if(s2dfs_p3_flag) {
      if(strcmp(s2dfs_p3_device, "?") == 0)
	printf("Specify plotting device to show the S2DFS P3 map: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, s2dfs_p3_device);
      strcpy(pgplot_options.box.xlabel, "Block number");
      strcpy(pgplot_options.box.ylabel, "P3 [cpp]");
      strcpy(pgplot_options.box.title, "S2DFS");
      /*      pgplotMap(&pgplot_options, s2dfs_p3, fin[0].NrSubints-fft_size+1, fft_size/2+1, 0, fin[0].NrSubints-fft_size+1, 0, fin[0].NrSubints-fft_size+1, 0, 0.5, 0, 0.5, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, "Block number", "P3 [cpp]", "S2DFS", 1, 1, 1, -0.5, -0.5, 1.0, 1.0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state); */
      pgplotMap(&pgplot_options, s2dfs_p3, fin[0].NrSubints-fft_size+1, fft_size/2+1, 0, fin[0].NrSubints-fft_size+1, 0, fin[0].NrSubints-fft_size+1, 0, 0.5, 0, 0.5, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state); 
      if(write_flag) {
	fout.NrSubints = fft_size/2+1;
	fout.NrBins = fin[0].NrSubints-fft_size+1;
	fout.gentype = GENTYPE_S2DFSP3;
	fout.tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
	if(fout.tsub_list != NULL)
	  free(fout.tsub_list);
	fout.tsub_list = (double *)malloc(sizeof(double));
	if(fout.tsub_list == NULL) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
	  return 0;
	}
	fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
	fout.yrangeset = 1;
	fout.yrange[0] = 0;
	fout.yrange[1] = 0.5;
	if(change_filename_extension(argv[argc-1], outputname, "s2dfs_p3", 1000, application.verbose_state) == 0) {
	  return 0;
	}
	if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	  return 0;
	if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
	  return 0;
	}
	//	appendHistoryLine(fout, argc, argv, application.verbose_state);
	if(writePSRData(&fout, s2dfs_p3, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
	  return 0;
	}
	closePSRData(&fout, 1, application.verbose_state); // Keep history, as header can be written out again
      }
    }
    if(s2dfs_p2_flag) {
      if(strcmp(s2dfs_p2_device, "?") == 0)
	printf("Specify plotting device to show the S2DFS P2 map: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, s2dfs_p2_device);
      strcpy(pgplot_options.box.xlabel, "Block number");
      strcpy(pgplot_options.box.ylabel, "P2 [cpp]");
      strcpy(pgplot_options.box.title, "S2DFS");
      /*	pgplotMap(&pgplot_options, s2dfs_p2, fin[0].NrSubints-fft_size+1, (application.onpulse.right[0]-application.onpulse.left[0]+1), 0, fin[0].NrSubints-fft_size+1, 0, fin[0].NrSubints-fft_size+1, -fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right[0]-application.onpulse.left[0]+1), fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right[0]-application.onpulse.left[0]+1), p2min, p2max, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, "Block number", "P2 [cpp]", "S2DFS", 1, 1, 1, -0.5, -0.5, 1.0, 1.0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state);  */
      if(p2range_set == 0) {
	p2min = -fin[0].NrBins/2-0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1);
	p2max = +fin[0].NrBins/2-0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1);
      }
      pgplotMap(&pgplot_options, s2dfs_p2, fin[0].NrSubints-fft_size+1, (application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1), 0, fin[0].NrSubints-fft_size+1, 0, fin[0].NrSubints-fft_size+1, -fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1), fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1), p2min, p2max, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state); 
      if(write_flag) {
	fout.NrSubints = (application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1);
	fout.NrBins = fin[0].NrSubints-fft_size+1;
	fout.gentype = GENTYPE_S2DFSP2;
	fout.tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
	if(fout.tsub_list != NULL)
	  free(fout.tsub_list);
	fout.tsub_list = (double *)malloc(sizeof(double));
	if(fout.tsub_list == NULL) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
	  return 0;
	}
	fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
	fout.yrangeset = 1;
	fout.yrange[0] = -fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1);
	fout.yrange[1] = fin[0].NrBins/2.0 -0.5*fin[0].NrBins/(float)(application.onpulse.right_bin[0]-application.onpulse.left_bin[0]+1);
	if(change_filename_extension(argv[argc-1], outputname, "s2dfs_p2", 1000, application.verbose_state) == 0) {
	  return 0;
	}
	if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	  return 0;
	if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
	  return 0;
	}
	//	appendHistoryLine(fout, argc, argv, application.verbose_state);
	if(writePSRData(&fout, s2dfs_p2, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
	  return 0;
	}
	closePSRData(&fout, 1, application.verbose_state); // Keep history, as header can be written out again
      }
    }
    free(twodfs);
    free(s2dfs_p3);
    free(s2dfs_p2);
  }  // End of if(s2dfs_p3_flag || s2dfs_p2_flag) 

//START REGION DEVELOP

  if(s_track_flag || s_modmap_flag) {
    printf("Calculating SLIDING SUBPULSE PHASE TRACK\n");
    lrfs = (float *)malloc((fft_size/2+1)*fin[0].NrBins*sizeof(float));
    phase_track = (float *)malloc(fin[0].NrBins*sizeof(float));
    s_phase_track = (float *)malloc(fin[0].NrBins*(fin[0].NrSubints-fft_size+1)*sizeof(float));
    s_mod_map = (float *)malloc(fin[0].NrBins*(fin[0].NrSubints-fft_size+1)*sizeof(float));
    stddev = (float *)malloc(fin[0].NrBins*sizeof(float));
    modindex = (float *)malloc(fin[0].NrBins*sizeof(float));
    rms_modindex = (float *)malloc(fin[0].NrBins*sizeof(float));
    if(lrfs == NULL || phase_track == NULL || stddev == NULL || modindex == NULL || rms_modindex == NULL || s_phase_track == NULL || s_mod_map == NULL) {
      printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
      return 0;
    }
    for(i = 0; i < fin[0].NrSubints-fft_size+1; i++) {
      verbose_definition noverbose;
      copyVerboseState(application.verbose_state, &noverbose);
      noverbose.verbose = 0;
      noverbose.nocounters = 1;
      if(calcLRFS(&fin[0].data[i*fin[0].NrBins], fft_size, fin[0].NrBins, fft_size, lrfs, subtractDC, NULL, phase_track, NULL, s_track_flag, freq_min, freq_max, track_only_first_region, NULL, 0, ftrack_mask, 0, &application.onpulse, &var_rms, argc, argv, noverbose) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate LRFS");
	return 0;
      }
      if(s_modmap_flag) {
	verbose_definition noverbose;
	cleanVerboseState(&noverbose);
	noverbose.verbose = 0;
	noverbose.nocounters = 1;
	read_partprofilePSRData(fin[0], profileI, NULL, 0, i, fft_size, noverbose);
	if(modSimple_flag) {
	  printerror(application.verbose_state.debug, "ERROR pspec: The -modSimple flag is not implemented for the modulation index map!");
	  return 0;
	}
	calcModindex(lrfs, profileI, fin[0].NrBins, fft_size, fft_size, stddev, rms_sigma, modindex, rms_modindex, &application.onpulse, var_rms, NULL, application.verbose_state);
	for(j = 0; j < fin[0].NrBins; j++) {
	  if(modindex[j] < 0 || modindex[j] < 3*rms_modindex[j])
	    modindex[j] = 0;
	  s_mod_map[j*(fin[0].NrSubints-fft_size+1)+i] = modindex[j];
	}
      }
      if(s_track_flag) {
	for(j = 0; j < fin[0].NrBins; j++) {
	  float x;
	  x = (xmax-xmin)*j/(float)(fin[0].NrBins-1);
	  s_phase_track[j*(fin[0].NrSubints-fft_size+1)+i] = phase_track[j] - slope*x + strack_dphase;
	  if(p3_subtract_flag)
	    s_phase_track[j*(fin[0].NrSubints-fft_size+1)+i] += 360.0*i/(float)p3_subtract;
	  s_phase_track[j*(fin[0].NrSubints-fft_size+1)+i] = derotate_deg(s_phase_track[j*(fin[0].NrSubints-fft_size+1)+i]);
	  if(strack_mask && s_modmap_flag) {
	    if(s_mod_map[j*(fin[0].NrSubints-fft_size+1)+i] <= 0) {
	      if(strack_mask == 2)
		s_phase_track[j*(fin[0].NrSubints-fft_size+1)+i] =  -1;
	      else
		s_phase_track[j*(fin[0].NrSubints-fft_size+1)+i] =  0;
	    }
	  }
	}
      }
      printf("  Block %ld of the %ld  (%.2f%%)      \r", i, fin[0].NrSubints-fft_size+1, 100*(i+1)/(float)(fin[0].NrSubints-fft_size+1));
    }
    if(s_modmap_flag) {
      if(strcmp(s_modmap_device, "?") == 0)
	printf("Specify plotting device to show the sliding modulation index map: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, s_modmap_device);
      strcpy(pgplot_options.box.xlabel, "Block number");
      strcpy(pgplot_options.box.ylabel, "Pulse longitude [deg]");
      strcpy(pgplot_options.box.title, "SLIDING MODULATION INDEX MAP");
      /*      pgplotMap(&pgplot_options, s_mod_map, fin[0].NrSubints-fft_size+1, fin[0].NrBins, 0, fin[0].NrSubints-fft_size+1, 0, fin[0].NrSubints-fft_size+1, xmin, xmax, xmin_zoom, xmax_zoom, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, "Block number", "Pulse longitude [deg]", "SLIDING MODULATION INDEX MAP", 1, 1, 1, -0.5, -0.5, 1.0, 1.0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state); */
      pgplotMap(&pgplot_options, s_mod_map, fin[0].NrSubints-fft_size+1, fin[0].NrBins, 0, fin[0].NrSubints-fft_size+1, 0, fin[0].NrSubints-fft_size+1, xmin, xmax, xmin_zoom, xmax_zoom, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state); 
      if(write_flag) {
	fout.NrSubints = fin[0].NrBins;
	fout.NrBins = fin[0].NrSubints-fft_size+1;
	if(change_filename_extension(argv[argc-1], outputname, "smod", 1000, application.verbose_state) == 0) {
	  return 0;
	}
	if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	  return 0;
	if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
	  return 0;
	}
	//	appendHistoryLine(fout, argc, argv, application.verbose_state);
	if(writePSRData(&fout, s_mod_map, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
	  return 0;
	}
	closePSRData(&fout, 1, application.verbose_state); // Keep history, as header can be written out again
      }
    }
    if(s_track_flag) {
      if(strcmp(s_track_device, "?") == 0)
	printf("Specify plotting device to show the sliding subpulse phase track: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, s_track_device);
      strcpy(pgplot_options.box.xlabel, "Block number");
      strcpy(pgplot_options.box.ylabel, "Pulse longitude [deg]");
      strcpy(pgplot_options.box.title, "SLIDING SUBPULSE PHASE TRACK");
      /*      pgplotMap(&pgplot_options, s_phase_track, fin[0].NrSubints-fft_size+1, fin[0].NrBins, 0, fin[0].NrSubints-fft_size+1, 0, fin[0].NrSubints-fft_size+1, xmin, xmax, xmin_zoom, xmax_zoom, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, "Block number", "Pulse longitude [deg]", "SLIDING SUBPULSE PHASE TRACK", 1, 1, 1, -0.5, -0.5, 1.0, 1.0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state); */
      pgplotMap(&pgplot_options, s_phase_track, fin[0].NrSubints-fft_size+1, fin[0].NrBins, 0, fin[0].NrSubints-fft_size+1, 0, fin[0].NrSubints-fft_size+1, xmin, xmax, xmin_zoom, xmax_zoom, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state); 
      if(write_flag) {
	fout.NrSubints = fin[0].NrBins;
	fout.NrBins = fin[0].NrSubints-fft_size+1;
	if(change_filename_extension(argv[argc-1], outputname, "strack", 1000, application.verbose_state) == 0) {
	  return 0;
	}
	if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	  return 0;
	if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
	  return 0;
	}
	//	appendHistoryLine(fout, argc, argv, application.verbose_state);
	if(writePSRData(&fout, s_phase_track, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
	  return 0;
	}
	closePSRData(&fout, 1, application.verbose_state); // Keep history, as header can be written out again
	if(strack_2ndfile_flag) {
	  fout.NrSubints = fin[0].NrBins;
	  fout.NrBins = fin[0].NrSubints-fft_size+1;
	  for(i = 0; i < fin[0].NrSubints-fft_size+1; i++) {
	    for(j = 0; j < fin[0].NrBins; j++) {
	      s_phase_track[j*(fin[0].NrSubints-fft_size+1)+i] = derotate_deg(s_phase_track[j*(fin[0].NrSubints-fft_size+1)+i] + strack_dphase2);
	      if(strack_mask && s_modmap_flag) {
		if(s_mod_map[j*(fin[0].NrSubints-fft_size+1)+i] <= 0) {
		  s_phase_track[j*(fin[0].NrSubints-fft_size+1)+i] =  0;
		}
	      }
	    }
	  }
	  if(change_filename_extension(argv[argc-1], outputname, "strack2", 1000, application.verbose_state) == 0) {
	    return 0;
	  }
	  if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	    return 0;
	  if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
	    return 0;
	  }
	  //	  appendHistoryLine(fout, argc, argv, application.verbose_state);
	  if(writePSRData(&fout, s_phase_track, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
	    return 0;
	  }
	  closePSRData(&fout, 1, application.verbose_state); // Keep history, as header can be written out again
	}
      }
    }
    free(lrfs);
    free(phase_track);
    free(s_phase_track);
    free(s_mod_map);
    free(stddev);
    free(modindex);
  } // End of if(s_track_flag || s_modmap_flag)
  
//START REGION DEVELOP
  if(lrac_flag) {
    float *lrac, *lrac2;
    long cclength, cclength2;
    int remove_wf;
    if(lrac_flag == 2)
      remove_wf = 0;
    else
      remove_wf = 1;
    // Note that the data is zero-padded with at least NrSubints zero's to avoid wrap problems.
    cclength = crosscorrelation_fft_padding_cclength(fin[i].NrSubints, fin[i].NrSubints);
    if(lrac_keepbothhalves == 0)
      cclength /= 2;
    lrac = (float *)malloc(originalNrPols*cclength*fin[0].NrBins*sizeof(float));
    if(lrac == NULL) {
      printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
      return 0;
    }

    for(i = 0; i < originalNrPols; i++) {
      if(calcLRAC(fin[i].data, fin[i].NrSubints, fin[i].NrBins, &lrac[i*cclength*fin[0].NrBins], remove_wf, lrac_zap0, lrac_keepbothhalves, &cclength2, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate LRAC");
	return 0;
      }
      if(cclength != cclength2) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR pspec: There is an unexpected mismatch between the expected and actual length of the auto-correlation output.\n");
	return 0;
      }
    }

    if(strcmp(lracdevice, "?") == 0)
      printf("Specify plotting device to show the LRAC: \n  ");
    strcpy(pgplot_options.viewport.plotDevice, lracdevice);
    strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
    strcpy(pgplot_options.box.ylabel, "Lag (subint number)");
    strcpy(pgplot_options.box.title, "LRAC");
    ret = get_period(fin[0], 0, &period, application.verbose_state);
    if(ret == 2) {
      printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
      return 0;
    }
    float y1, y2;
    if(lrac_keepbothhalves) {
      y1 = -cclength/2;
      y2 = cclength/2-1;
    }else {
      y1 = 0;
      y2 = cclength-1;
    }
    pgplotMap(&pgplot_options, lrac, fin[0].NrBins, cclength, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period, xmin_zoom, xmax_zoom, y1, y2, y1, y2, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, application.verbose_state); 
    if(write_flag) {
      fout.NrSubints = cclength;
      fout.NrBins = fin[0].NrBins;
      fout.NrPols = originalNrPols;
      fout.gentype = GENTYPE_LRAC;
      fout.yrangeset = 1;
      fout.yrange[0] = y1;
      fout.yrange[1] = y2;
      //      printf("XXXX set range to : %f %f\n", fout.yrange[0], fout.yrange[1]);
      fout.tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
      if(fout.tsub_list != NULL)
	free(fout.tsub_list);
      fout.tsub_list = (double *)malloc(sizeof(double));
      if(fout.tsub_list == NULL) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
	return 0;
      }
      fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
      //	fprintf(stderr, "XXXXX %lf\n", fout.fixedtsub);

      if(change_filename_extension(argv[argc-1], outputname, "lrac", 1000, application.verbose_state) == 0) {
	return 0;
      }
      if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	return 0;
      if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.\n");
	return 0;
      }
      //	appendHistoryLine(fout, argc, argv, application.verbose_state);

      /* Transpose the data for writing */
      lrac2 = (float *)malloc(originalNrPols*cclength*fin[0].NrBins*sizeof(float));
      if(lrac2 == NULL) {
	printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
	return 0;
      }

      for(i = 0; i < cclength; i++) {
	for(p = 0; p < originalNrPols; p++) {
	  memcpy(&lrac2[(originalNrPols*i+p)*fin[0].NrBins], &lrac[(p*cclength+i)*fin[0].NrBins], fin[0].NrBins*sizeof(float));
	}
      }
      if(writePSRData(&fout, lrac2, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.\n");
	return 0;
      }
      free(lrac2);
      closePSRData(&fout, 1, application.verbose_state); // Keep history, as header can be written out again
      fout.gentype = GENTYPE_UNDEFINED;
      fout.yrangeset = 0;
      fout.xrangeset = 0;
      fout.NrPols = 1;
    }
    free(lrac);
  } // End of if(lrac_flag)

//START REGION DEVELOP
  if(p3class) {

    float *p3class_data = malloc(fin[0].NrSubints*sizeof(float));
    if(p3class_data == NULL) {
      printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
      return 0;
    }

    // Only analyse first polarization channel
    if(p3classify(fin[0].data, fin[0].NrSubints, fin[0].NrBins, p3class_npad, p3class_nblock, p3class_nmin, p3class_cpp1, p3class_cpp2, p3class_data, &application.onpulse, argc, argv, application.verbose_state) == 0) {
      /* Plot only plot half the spectrum */
      printerror(application.verbose_state.debug, "ERROR pspec: P3 classification failed");
      return 0;
    }
  }



//START REGION DEVELOP

//  printf("XXXXXX After P3 fold: %p\n", fin[4].tsub_list);

  if(hrfs_flag) {
    float *hrfs;
    long hrfs_length;
    if(hrfs_flag) {
      printwarning(application.verbose_state.debug, "WARNING pspec: The frequency labelling might not be entirely correct and the code is untested for uneven nr of samples per period (actually it is not tested for non-power of two number of samples) and only with fftw3.");
    }
    if(calcHRFS(fin[0].data, fin[0].NrSubints, fin[0].NrBins, fft_size, &hrfs, &hrfs_length, &application.onpulse, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate HRFS");
      return 0;
    }
    if(strcmp(hrfsdevice, "?") == 0)
      printf("Specify plotting device to show the HRFS: \n  ");
    strcpy(pgplot_options.viewport.plotDevice, hrfsdevice);
    strcpy(pgplot_options.box.xlabel, "cpp");
    strcpy(pgplot_options.box.title, "HRFS");
    if(hrfs_flag == 1) {
      strcpy(pgplot_options.box.ylabel, "Spectral power");
      pgplotGraph1(&pgplot_options, hrfs, NULL, NULL, hrfs_length, 
		   0.0, fin[0].NrBins/2.0, 
		   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
    }else {
      strcpy(pgplot_options.box.ylabel, "Harmonic number");
      printwarning(application.verbose_state.debug, "WARNING pspec: The power corresponding to the shape of the average emission is masked.");
      for(i = 0; i < hrfs_length/fft_size; i++)
	hrfs[i*fft_size] = 0;
      pgplotMap(&pgplot_options, hrfs, fft_size, hrfs_length/fft_size,
		0.0, 1.0, 0.0, 1.0,
		0.0, fft_size/2, 0.0, fft_size/2,
		PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1.0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, application.verbose_state);
    }

    if(write_flag) {
      if(hrfs_flag == 1) {
	fout.NrSubints = 1;
	fout.NrBins = hrfs_length;
	fout.gentype = GENTYPE_HRFS_UNFOLDED;
	fout.tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
	if(fout.tsub_list != NULL)
	  free(fout.tsub_list);
	fout.tsub_list = (double *)malloc(sizeof(double));
	if(fout.tsub_list == NULL) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
	  return 0;
	}
	fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
	fout.xrangeset = 1;
	fout.yrangeset = 0;
	fout.xrange[0] = 0;
	fout.xrange[1] = fin[0].NrBins/2.0;
      }else {
	fout.NrSubints = hrfs_length/fft_size;
	fout.NrBins = fft_size;
	fout.gentype = GENTYPE_HRFS;
	fout.tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
	if(fout.tsub_list != NULL)
	  free(fout.tsub_list);
	fout.tsub_list = (double *)malloc(sizeof(double));
	if(fout.tsub_list == NULL) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
	  return 0;
	}
	fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
	fout.xrangeset = 1;
	fout.yrangeset = 0;
	fout.xrange[0] = 0;
	fout.xrange[1] = 1.0;
      }
      if(change_filename_extension(argv[argc-1], outputname, "hrfs", 1000, application.verbose_state) == 0) {
	return 0;
      }
      if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	return 0;
      if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.");
	return 0;
	}
      //	appendHistoryLine(fout, argc, argv, application.verbose_state);
      if(writePSRData(&fout, hrfs, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.");
	return 0;
      }
      closePSRData(&fout, 1, application.verbose_state); // Keep history, as header can be written out again
      fout.gentype = GENTYPE_UNDEFINED;
      fout.yrangeset = 0;
      fout.xrangeset = 0;
    }
    free(hrfs);
  }  // End of if(hrfs_flag)

//START REGION DEVELOP

  if(lrcc_flag) {
    lrcc = malloc(fin[0].NrBins*fin[0].NrBins*sizeof(float));
    if(lrcc == NULL) {
      printerror(application.verbose_state.debug, "ERROR pspec: Cannot allocate memory");
      return 0;
    }
    for(j = lrcc_lag; j <= lrcc_lagstop; j++) {
      if(calcLRCC(fin[0].data, fin[0].NrSubints, fin[0].NrBins, j, lrcc, lrcc_nosubmean, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR pspec: Cannot calculate LRCC");
	return 0;
      }
      if(strcmp(lrccdevice, "?") == 0)
	printf("Specify plotting device to show the LRCC: \n  ");
      strcpy(pgplot_options.viewport.plotDevice, lrccdevice);
      strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.ylabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.title, "LRCC");
      /*      pgplotMap(&pgplot_options, lrcc, fin[0].NrBins, fin[0].NrBins, 0, 360*(fin[0].NrBins-1)*fin[0].get_tsamp(fin[0], 0, application.verbose_state)/get_period(fin[0], 0, application.verbose_state), xmin_zoom, xmax_zoom, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/get_period(fin[0], 0, application.verbose_state), xmin_zoom, xmax_zoom, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, "Pulse phase [degrees]", "Pulse phase [degrees]", "LRCC", 1, 1, 1, -0.5, -0.5, 1.0, 1.0, 1, 1, 0, 1, 0, 1, 0, 0, application.verbose_state); */
      ret = get_period(fin[0], 0, &period, application.verbose_state);
      if(ret == 2) {
	printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
	return 0;
      }
      pgplotMap(&pgplot_options, lrcc, fin[0].NrBins, fin[0].NrBins, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period, xmin_zoom, xmax_zoom, 0, 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period, xmin_zoom, xmax_zoom, PPGPLOT_INVERTED_HEAT, application.itf, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, application.verbose_state); 
      if(write_flag) {
	fout.NrSubints = fin[0].NrBins;
	fout.NrBins = fin[0].NrBins;
	fout.gentype = GENTYPE_LRCC;
	fout.tsubMode = TSUBMODE_FIXEDTSUB;              // Only one "subint" which has length of total observation
	if(fout.tsub_list != NULL)
	  free(fout.tsub_list);
	fout.tsub_list = (double *)malloc(sizeof(double));
	if(fout.tsub_list == NULL) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR pspec: Memory allocation error");
	  return 0;
	}
	fout.tsub_list[0] = get_tobs(fin[0], application.verbose_state);
	fout.xrangeset = 0;
	fout.yrangeset = 1;
	fout.yrange[0] = 0;
	ret = get_period(fin[0], 0, &period, application.verbose_state);
	if(ret == 2) {
	  printerror(application.verbose_state.debug, "ERROR pspec (%s): Cannot obtain period", fin[0].filename);
	  return 0;
	}
	fout.yrange[1] = 360*(fin[0].NrBins-1)*get_tsamp(fin[0], 0, application.verbose_state)/period;
	
	if(j >= 0)
	  txt[0] = 'p';
	else
	  txt[0] = 'm';
	sprintf(txt+1, "%04ld.lrcc", labs(j));
	if(change_filename_extension(argv[argc-1], outputname, txt, 1000, application.verbose_state) == 0) {
	  return 0;
	}
	if(!openPSRData(&fout, outputname, application.oformat, 1, 0, 0, application.verbose_state))
	  return 0;
	if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write header.");
	  return 0;
	}
	//	appendHistoryLine(fout, argc, argv, application.verbose_state);
	if(writePSRData(&fout, lrcc, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspec: Unable to write data.");
	  return 0;
	}
	closePSRData(&fout, 1, application.verbose_state);  // Keep history, as header can be written out again
	fout.gentype = GENTYPE_UNDEFINED;
	fout.yrangeset = 0;
	fout.xrangeset = 0;
      }
    }
    free(lrcc);

  }  // End of if(lrcc_flag)

//START REGION RELEASE


  closePSRData(&fout, 0, application.verbose_state);  // Release header 
  for(i = 0; i < MaxNrPolarizations; i++)
    closePSRData(&fin[i], 0, application.verbose_state);
  free(profileI);

  ppgend();
//START REGION DEVELOP
  if(application.ppgplot_name != NULL) {
    pgcloseoutputfile();
  }

//START REGION RELEASE
  if(bootstrap > 0) {
    free(stddev_av);
    free(modindex_av);
    free(stddev_square);
    free(modindex_square);
    free(clone_profileI);
  }

  terminateApplication(&application);
  return 0;
}

int pgplotProfile(char *plotDevice, int windowwidth, int windowheight, float *profile, float *stddev, float *rms_stddev, float *modindex, float *rms_modindex, int nrx, float xmin, float xmax, char *xlabel, char *ylabel, char *title, int stddev_flag, int mod_flag, int zoom_flag, float xmin_zoom, float xmax_zoom, verbose_definition verbose)
{
  float min, max;
  float x, y, dy;
  int i, deviceID;
  deviceID = ppgopen(plotDevice);
  if(deviceID <= 0) {
    printerror(verbose.debug, "ERROR pspec: Cannot open plot device");
    return 0;
  }
  /*  fprintf(stderr, "nrx=%d   xmax=%f\n", nrx, xmax); */
  if(windowwidth > 0 && windowheight > 0) {
    x = windowwidth*0.01175548589341692789994673739445429916373;
    y = (windowheight-1)/(float)windowwidth;
    ppgpap(x,y);
  }
  ppgask(0);
  ppgslw(1);
  ppgpage();
  ppgsvp(0.1, 0.9, 0.1, 0.9);

  min = max = profile[0];
  for(i = 1; i < nrx; i++) {
    if(profile[i] > max)
      max = profile[i];
    if(profile[i] < min)
      min = profile[i];
    if(stddev_flag) {
      if(stddev[i] > max)
	max = stddev[i];
    }
    if(mod_flag) {
      if(modindex[i] > 3*rms_modindex[i]) {
	if(modindex[i]+rms_modindex[i] > max)
	  max = modindex[i]+rms_modindex[i];
      }
    }
  }

  ppgsci(1);
  ppgswin(xmin_zoom, xmax_zoom, min, max*1.03);
  ppgbox("bcnsti",0.0,0,"bcntsi",0.0,0);
  ppglab(xlabel, ylabel, title);

  ppgsci(1);

  ppgslw(5);
  for(i = 0; i < nrx; i++) {
    x = i*(xmax-xmin)/(float)nrx + xmin;
    y = profile[i];
    if(i == 0)
      ppgmove(x, y);
    else
      ppgdraw(x, y);
  }

  ppgslw(1);
  if(stddev_flag) {
    for(i = 0; i < nrx; i++) {
      x = i*(xmax-xmin)/(float)nrx + xmin;
      y = stddev[i];
      if(y > 3*rms_stddev[i]) {
	ppgpt1(x, y, 4);
      }
    }
  }
  if(mod_flag) {
    for(i = 0; i < nrx; i++) {
      x = i*(xmax-xmin)/(float)nrx + xmin;
      y = modindex[i];
      dy = rms_modindex[i];
      if(y > 3*dy) {
	ppgslw(5);
	ppgpt1(x, y, -1);
	ppgslw(1);
	ppgerr1(6, x, y, dy, 3);
      }
    }
  }

  ppgclos();
  return 1;
}
