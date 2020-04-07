//START REGION RELEASE
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "psrsalsa.h"

double f2_min, f2_max, f3_min, f3_max, fl_min, fl_max, f2n_min, f2n_max, f3n_min, f3n_max;
double twodsf_Imin;
double I_noise_max, I_noise_min;
float centroid_x, centroid_y;
float centroid_err_x, centroid_err_y;
int centroid_calculated;
int SelectedPlot;               /* 0 = 2dfs, 1 = lrfs, 2 = noise, 3 = line plots */
int PlotAvrgProfile;
int SelectedComponent;
int Centered;
int KeyCode;
datafile_definition twodfs, lrfs, AverageProfile, noise;

//START REGION DEVELOP
int AutoCalculateCentroid;
int NoiseCalculation;
int NrSelectedPulseRegions;
int PulseRegions[10][2];
//START REGION RELEASE

#define Max_nr_noise_patches 100
#define MaxNrNoiseBins 100

int f2npatch_min[Max_nr_noise_patches],f2npatch_max[Max_nr_noise_patches],f3npatch_min[Max_nr_noise_patches],f3npatch_max[Max_nr_noise_patches];
int nr_noise_patches;


void PlotWindow(verbose_definition verbose);
int MakeSelection(double sigma_noise);
void SelectRegion();
double calculate_noise_sigma(void);
void calculate_2dfs_Centroid(double sigma_noise, verbose_definition verbose);
//START REGION DEVELOP
int WhichRegion(int bin);
//START REGION RELEASE

int main(int argc, char **argv)
{
  int i, xi, yi, index;
  double I;
  char filename[1000], txt[1000], *input_filename_ptr;
//START REGION DEVELOP
  int SelectP3Region, SelectP2Region;
  double P3RegionLow, P3RegionHigh, P2RegionLow, P2RegionHigh;
//START REGION RELEASE
  double sigma_noise;
  psrsalsaApplication application;

  initApplication(&application, "pspecDetect", "[options] pulse_stack,\nwhere pulse_stack is the file name of the pulse stack that has been processed by\npspec to produce the 2DFS and LRFS.");

  //  application.switch_nocounters = 1;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_device = 1;

  cleanPSRData(&twodfs, application.verbose_state);
  cleanPSRData(&lrfs, application.verbose_state);
  cleanPSRData(&AverageProfile, application.verbose_state);
  cleanPSRData(&noise, application.verbose_state);
  closePSRData(&twodfs, 0, application.verbose_state);
  closePSRData(&lrfs, 0, application.verbose_state);
  closePSRData(&AverageProfile, 0, application.verbose_state);
  closePSRData(&noise, 0, application.verbose_state);
 
  SelectedPlot = 0;
  sigma_noise = 0;
  centroid_calculated = 0;
  nr_noise_patches = 0;
  SelectedComponent = 1;
  Centered = 1;
//START REGION DEVELOP
  NoiseCalculation = 1;
  SelectP3Region = 0;
  SelectP2Region = 0;
  NrSelectedPulseRegions = 0;
//START REGION RELEASE
  PlotAvrgProfile = 1;

//START REGION DEVELOP
  AutoCalculateCentroid = 0;
//START REGION RELEASE
  if(argc < 2) {
    printf("Interactive program designed to analyse features in the 2DFS to obtain\ncentroid P2 and P3 values and corresponding error-bars.\n\n");
    printApplicationHelp(&application);
//START REGION DEVELOP
    printf("Where optional options are:\n\n");
    printf("-p3 \"low high\"    Select vertical region 2dfs\n");
    printf("-p2 \"low high\"    Select horizontal region 2dfs\n");
    printf("-a                Automaticly calculate centroid\n");
    printf("-c                Select pulse longitude range\n");
    printf("-dn               Disable noise error\n");
    //    printf("-d  \"...\"         Plot device\n");
    //    printf("-v                Verbose mode\n");
//START REGION RELEASE
    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about the how to use the centroid information can be found in:\n");
    printf(" - Weltevrede et al. 2006, A&A, 445, 243\n");
    printf(" - Weltevrede et al. 2007, A&A, 469, 607.\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "-A") == 0) {
	AutoCalculateCentroid = 1;
      }else if(strcmp(argv[i], "-p3") == 0 || strcmp(argv[i], "-P3") == 0) {
	SelectP3Region = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &P3RegionLow, &P3RegionHigh, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspecDetect: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-p2") == 0 || strcmp(argv[i], "-P2") == 0) {
	SelectP2Region = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &P2RegionLow, &P2RegionHigh, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspecDetect: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "-C") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &PulseRegions[NrSelectedPulseRegions][0], &PulseRegions[NrSelectedPulseRegions][1], NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pspecDetect: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	NrSelectedPulseRegions++;
	i++;
      }else if(strcmp(argv[i], "-dn") == 0 || strcmp(argv[i], "-DN") == 0) {
	NoiseCalculation = 0;
//START REGION RELEASE
      }else {
	/* If the option is not recognized, assume it is a filename */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "pspecDetect: Unknown option: %s\n\nRun pspecDetect without command line arguments to show help", argv[i]);
	  terminateApplication(&application);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }


  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pspecDetect: No files specified");
    return 0;
  }


  input_filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state);

  sprintf(txt, "%d.2dfs", 1);
  if(change_filename_extension(input_filename_ptr, filename, txt, 1000, application.verbose_state) == 0) {
    return 0;
  }
  if(application.verbose_state.verbose)
    printf("Reading %s\n", filename);
  if(!openPSRData(&twodfs, filename, 0, 0, 1, 0, application.verbose_state))
    return 0;
  if(twodfs.NrPols > 1) {
    datafile_definition clone;
    if(preprocess_polselect(twodfs, &clone, 0, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "Cannot select first polarization channel\n");
      return 0;    
    }
    swap_orig_clone(&twodfs, &clone, application.verbose_state);
  }
  twodsf_Imin = 0;
  for(xi = 0; xi < twodfs.NrBins; xi++) {
    for(yi = 0; yi < twodfs.NrSubints; yi++) {
      I = twodfs.data[yi*twodfs.NrBins+xi];
      if(I < twodsf_Imin)
	twodsf_Imin = I;
    }    
  }
  if(application.verbose_state.verbose)
    printf("%ldx%ld points read from 2dfs\n", twodfs.NrBins, twodfs.NrSubints);
  /*
    Patrick: needs possibly fixing:
	p2min = -fin.NrBins/2-0.5*fin.NrBins/(float)(onpulseRegions.right[regionnr]-onpulseRegions.left[regionnr]+1);
	p2max = +fin.NrBins/2-0.5*fin.NrBins/(float)(onpulseRegions.right[regionnr]-onpulseRegions.left[regionnr]+1);
   */


  if(preprocess_polselect(twodfs, &noise, 0, application.verbose_state) != 1)
    exit(0);


  sprintf(txt, "lrfs");
  if(change_filename_extension(input_filename_ptr, filename, txt, 1000, application.verbose_state) == 0) {
    return 0;
  }
  if(application.verbose_state.verbose)
    printf("Reading %s\n", filename);
  if(!openPSRData(&lrfs, filename, 0, 0, 1, 0, application.verbose_state))
    return 0;
  if(application.verbose_state.verbose)
    printf("%ldx%ld points read from lrfs\n", lrfs.NrBins, lrfs.NrSubints);
  if(lrfs.NrPols > 1) {
    datafile_definition clone;
    if(preprocess_polselect(lrfs, &clone, 0, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "Cannot select first polarization channel\n");
      return 0;    
    }
    swap_orig_clone(&lrfs, &clone, application.verbose_state);
  }
  fl_min = 0;
  fl_max = lrfs.NrBins-1;


  if(application.verbose_state.verbose)
    printf("Reading %s\n", filename);
  if(!openPSRData(&AverageProfile, input_filename_ptr, 0, 0, 0, 0, application.verbose_state))
    return 0;
  if(!readHeaderPSRData(&AverageProfile, 0, 0, application.verbose_state))
    return 0;
  AverageProfile.data = malloc(AverageProfile.NrBins*sizeof(float));
  if(AverageProfile.data == NULL) {
    printerror(application.verbose_state.debug, "Memory allocation error\n");
    return 0;
  }
  if(!read_profilePSRData(AverageProfile, AverageProfile.data, NULL, 0, application.verbose_state))
    return 0;
  if(AverageProfile.NrBins != lrfs.NrBins) {
    printwarning(application.verbose_state.debug, "WARNING: It looks like data is rebinned? Check the units.");
    //    AverageProfile.fixedtsamp *= AverageProfile.NrBins/(double)lrfs.NrBins;
    //    AverageProfile.tsampMode = TSAMPMODE_FIXEDTSAMP;
    //    AverageProfile.NrBins = lrfs.NrBins;
    datafile_definition clone;
    AverageProfile.format = MEMORY_format;
    AverageProfile.NrPols = 1;
    AverageProfile.NrFreqChan = 1;
    AverageProfile.NrSubints = 1;
    if(preprocess_rebin(AverageProfile, &clone, lrfs.NrBins, application.verbose_state) == 0) {
      printwarning(application.verbose_state.debug, "WARNING: Rebinning of profile failed.");    
      return 0;
    }
    swap_orig_clone(&AverageProfile, &clone, application.verbose_state);
    printwarning(application.verbose_state.debug, "WARNING: Assuming the number of bins = %ld and the sampling time = %lf s.", AverageProfile.NrBins, AverageProfile.fixedtsamp);    
  }

  if(application.verbose_state.verbose)
    printf("%ld points read from avgprof\n\n", AverageProfile.NrBins);

  f2_min = f2n_min = -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins;
  f2_max = f2n_max = +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins;
  f3_min = 0;
  f3n_min = 0;
  f3_max = f3n_max = 0.5;

//START REGION DEVELOP
  if(SelectP3Region) {
    f3_min = P3RegionLow;
    f3_max = P3RegionHigh;
  }
  if(SelectP2Region) {
    f2_min = P2RegionLow;
    f2_max = P2RegionHigh;
  }
//START REGION RELEASE

  ppgopen(application.pgplotdevice);
  ppgask(0);
  ppgslw(1);
//START REGION DEVELOP
  if(!AutoCalculateCentroid) {
//START REGION RELEASE
    printf("Press h for help\n");
//START REGION DEVELOP
  }
//START REGION RELEASE
  SelectRegion();
  PlotWindow(application.verbose_state);
  do {
//START REGION DEVELOP
    if(AutoCalculateCentroid)
      KeyCode = 13;
    else
//START REGION RELEASE
      MakeSelection(sigma_noise);
    switch(KeyCode) {
    case 0: 
      PlotWindow(application.verbose_state);
      break;
    case 113:
    case 81:
    case 27: KeyCode = 27; break; 
    case 13: 
//START REGION DEVELOP
      if(NoiseCalculation)
//START REGION RELEASE
	sigma_noise = calculate_noise_sigma();
//START REGION DEVELOP
      else
        sigma_noise = 0;
//START REGION RELEASE
      calculate_2dfs_Centroid(sigma_noise, application.verbose_state);
      PlotWindow(application.verbose_state);
      break;
    case 67:
    case 99:
      if(Centered == 0) {
	printf("P2 Centring on\n");
	Centered = 1;
      }else {
	Centered = 0;
	printf("P2 Centring off\n");
      }
      printf("\nBe aware: P2 centring makes it less likely to detect a significant P2 offset just because the region included in the centroid calculation is offset. However, for clearly significant features you might want to turn P2 centring off to ensure the region included in the centroid calculation is centred on the feature, avoiding the corresponding bias.\n");
      break;
    case 82:
    case 114:
      sigma_noise = 0;
      centroid_calculated = 0;
      /*
      f2_min = f2n_min = twodfs.axes[2].bin_center(0).get(Unit());
      f2_max = f2n_max = twodfs.axes[2].bin_center(twodfs.axes[2].nelem-1).get(Unit());
      f3_min = twodfs.axes[1].bin_center(0).get(Unit());
      f3n_min = 0;
      f3_max = f3n_max = twodfs.axes[1].bin_center(twodfs.axes[1].nelem-1).get(Unit());
      fl_min = lrfs.axes[2].bin_center(0).get(Unit("deg"));
      fl_max = lrfs.axes[2].bin_center(lrfs.axes[2].nelem-1).get(Unit("deg"));
      */
      f2_min = f2n_min = -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins;
      f2_max = f2n_max = +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins;
      f3_min = 0;
      f3n_min = 0;
      f3_max = f3n_max = 0.5;
      fl_min = 0;
      fl_max = lrfs.NrBins-1;
//START REGION DEVELOP
      if(SelectP3Region) {
	f3_min = P3RegionLow;
	f3_max = P3RegionHigh;
      }
//START REGION RELEASE
      /* Reload noise to reset zapped spectral bins */
      if(noise.opened_flag)
        closePSRData(&noise, 0, application.verbose_state);
      cleanPSRData(&noise, application.verbose_state);
      if(preprocess_polselect(twodfs, &noise, 0, application.verbose_state) != 1)
        exit(0);

      nr_noise_patches = 0;
//START REGION DEVELOP
      NrSelectedPulseRegions = 0;
//START REGION RELEASE
      SelectRegion();
      PlotWindow(application.verbose_state);
      break;
    case 70:
    case 102:
      /*      centroid_calculated = 0; */   /* Can now unzoom while retaining the centroid fit */
      /*
      f2_min = twodfs.axes[2].bin_center(0).get(Unit());
      f2_max = twodfs.axes[2].bin_center(twodfs.axes[2].nelem-1).get(Unit());
      f3_min = twodfs.axes[1].bin_center(0).get(Unit());
      f3_max = twodfs.axes[1].bin_center(twodfs.axes[1].nelem-1).get(Unit());
      fl_min = lrfs.axes[2].bin_center(0).get(Unit("deg"));
      fl_max = lrfs.axes[2].bin_center(lrfs.axes[2].nelem-1).get(Unit("deg"));
      */
      f2_min = -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins;
      f2_max = +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins;
      f3_min = 0;
      f3_max = 0.5;
      fl_min = 0;
      fl_max = lrfs.NrBins-1;

//START REGION DEVELOP
      if(SelectP3Region) {
	f3_min = P3RegionLow;
	f3_max = P3RegionHigh;
      }
//START REGION RELEASE
      SelectRegion();
      PlotWindow(application.verbose_state);
      break;
    case 49:
    case 50:
    case 51:
    case 52:
    case 53:
    case 54:
    case 55:
    case 56:
    case 57: 
      SelectedComponent = KeyCode - 48;
      strcpy(filename, argv[argc - 1]);
      sprintf(txt, "%d.2dfs", SelectedComponent);
      if(change_filename_extension(input_filename_ptr, filename, txt, 1000, application.verbose_state) == 0) {
	return 0;
      }
      if(application.verbose_state.verbose)
	printf("Reading %s\n", filename);

      if(closePSRData(&twodfs, 0, application.verbose_state)) {
	printerror(0, "Closing file failed\n");
	return 0;
      }
      if(!openPSRData(&twodfs, filename, 0, 0, 1, 0, application.verbose_state))
	return 0;
      if(twodfs.NrPols > 1) {
	datafile_definition clone;
	if(preprocess_polselect(twodfs, &clone, 0, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "Cannot select first polarization channel\n");
	  return 0;    
	}
	swap_orig_clone(&twodfs, &clone, application.verbose_state);
      }

      if(application.verbose_state.verbose)
        printf("%ldx%ld points read from 2dfs\n", twodfs.NrBins, twodfs.NrSubints);
      twodsf_Imin = 0;
      for(xi = 0; xi < twodfs.NrBins; xi++) {
	for(yi = 0; yi < twodfs.NrSubints; yi++) {
	  I = twodfs.data[yi*twodfs.NrBins+xi];
	  if(I < twodsf_Imin)
	    twodsf_Imin = I;
	}    
      }
      sigma_noise = 0;
      centroid_calculated = 0; 
      nr_noise_patches = 0;
      f2_min = f2n_min = -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins;
      f2_max = f2n_max = +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins;
      f3_min = 0;
      f3n_min = 0;
      f3_max = f3n_max = 0.5;

      if(noise.opened_flag)
        closePSRData(&noise, 0, application.verbose_state);
      cleanPSRData(&noise, application.verbose_state);
      if(preprocess_polselect(twodfs, &noise, 0, application.verbose_state) != 1)
        exit(0);


      /*
     twodfs.read(filename);
      if(application.verbose_state.verbose)
        printf("%dx%d points read from 2dfs\n", twodfs.axes[2].nelem, twodfs.axes[1].nelem);
      for(xi = 0; xi < twodfs.axes[2].nelem; xi++) {
	for(yi = 0; yi < twodfs.axes[1].nelem; yi++) {
	  I = twodfs.data[yi*twodfs.axes[2].nelem+xi];
	  if(I < twodsf_Imin)
	    twodsf_Imin = I;
	}    
      }
      sigma_noise = 0;
      centroid_calculated = 0; 
      nr_noise_patches = 0;
      f2_min = f2n_min = twodfs.axes[2].bin_center(0).get(Unit());
      f2_max = f2n_max = twodfs.axes[2].bin_center(twodfs.axes[2].nelem-1).get(Unit());
      f3_min = twodfs.axes[1].bin_center(0).get(Unit());
      f3_max = twodfs.axes[1].bin_center(twodfs.axes[1].nelem-1).get(Unit());
      */



//START REGION DEVELOP
      if(SelectP3Region) {
	f3_min = P3RegionLow;
	f3_max = P3RegionHigh;
      }
//START REGION RELEASE
      SelectRegion();
      PlotWindow(application.verbose_state);
      break;
    case 32: 
      SelectedPlot++;
      if(SelectedPlot > 2
//START REGION DEVELOP
	 +2
//START REGION RELEASE
	 )
	SelectedPlot = 0;
//START REGION DEVELOP
      if(SelectedPlot == 3)   /* Skip plot 3, which is unnecessary */
	SelectedPlot = 4;
//START REGION RELEASE
      SelectRegion();
      PlotWindow(application.verbose_state);
      break;
    case 97: 
      PlotAvrgProfile++;
      if(PlotAvrgProfile > 1)
	PlotAvrgProfile = 0;
      SelectRegion();
      PlotWindow(application.verbose_state);
      break;
//START REGION DEVELOP
    case 115:
    case 83:   // S-option
      printf("Left and right bin:\n");
      scanf("%d %d", &PulseRegions[NrSelectedPulseRegions][0], &PulseRegions[NrSelectedPulseRegions][1]);
      NrSelectedPulseRegions++;
      SelectRegion();
      PlotWindow(application.verbose_state);
      break;
//START REGION RELEASE
    case 72: 
    case 104: 
    case 63:  // ?
      printf("\nGeneral strategy: Use space to show the 2dfs with the title \"For noise calculation ...\". Flag all signal only noise is visible. This sets the rms used in the error calculation. Use space to select the 2dfs and zoom in on the feature of interest. Depending on the situation, toggling 'c' would be beneficial (see help when using this option). Press return to calculate centroid of the zoomed in area (indicated by cross) and 1-sigma error box. The error is purely statistical, the systematic error resulting from the decision of what area to include in the centroid calculation can be assessed by selecting slightly different regions and repeat the calculation.\n\n");
      printf("h/?   = Help\n");
      printf("1..9  = Switch to component 1..9 (load new 2dfs)\n");
      printf("R     = Reset (zoom + noise-flagging settings)\n");
      printf("F     = Select new feature (but keep noise-flagging settings)\n");
//START REGION DEVELOP
      printf("S     = Select longitude range\n");
//START REGION RELEASE
      printf("C     = Toggle P2 centering (allows selection to be non-symmetric)\n");
      printf("A     = Toggle showing average profile superimposed over the LRFS\n");
      printf("SPACE = Switch between the 2dfs, lrfs, and the 2dfs samples included in the noise calculation\n");
      printf("ENTER = Calculate centroid\n");
      printf("Q/ESC = Quit\n");
      break;
    default: printf("Unknown key: %d\n", KeyCode); break;
    }
  }while(KeyCode != 27
//START REGION DEVELOP
	 && AutoCalculateCentroid == 0
//START REGION RELEASE
	 );

  ppgend();

  closePSRData(&twodfs, 0, application.verbose_state);
  closePSRData(&lrfs, 0, application.verbose_state);
  closePSRData(&AverageProfile, 0, application.verbose_state);
  if(noise.opened_flag)
    closePSRData(&noise, 0, application.verbose_state);

  terminateApplication(&application);
  return 0;
}


void SelectRegion()
{
  int i, xi, yi, offset;
  /*
  std::vector<IndexType> p0(2), p1(2);
  PhysDataArray<double, float> plot_subset2;

      if(SelectedPlot == 1)
        plot_subset2 = lrfs.get_hyperelem_ref(2, 0);
      else
        plot_subset2 = twodfs.get_hyperelem_ref(2, 0);
      if(SelectedPlot == 0) {
        p0[0] = plot_subset2.axes[0].ibin(PhysDataPoint<double>(f3_min, Unit()));
        p1[0] = plot_subset2.axes[0].ibin(PhysDataPoint<double>(f3_max, Unit()));
        p0[1] = plot_subset2.axes[1].ibin(PhysDataPoint<double>(f2_min, Unit()));
        p1[1] = plot_subset2.axes[1].ibin(PhysDataPoint<double>(f2_max, Unit()));
      }else if(SelectedPlot == 1) {
        p0[0] = plot_subset2.axes[0].ibin(PhysDataPoint<double>(f3_min, Unit()));
        p1[0] = plot_subset2.axes[0].ibin(PhysDataPoint<double>(f3_max, Unit()));
        p0[1] = plot_subset2.axes[1].ibin(PhysDataPoint<double>(fl_min, Unit("deg")));
        p1[1] = plot_subset2.axes[1].ibin(PhysDataPoint<double>(fl_max, Unit("deg")));
      }else if(SelectedPlot == 2 || SelectedPlot == 4) {
        p0[0] = plot_subset2.axes[0].ibin(PhysDataPoint<double>(f3n_min, Unit()));
        p1[0] = plot_subset2.axes[0].ibin(PhysDataPoint<double>(f3n_max, Unit()));
        p0[1] = plot_subset2.axes[1].ibin(PhysDataPoint<double>(f2n_min, Unit()));
        p1[1] = plot_subset2.axes[1].ibin(PhysDataPoint<double>(f2n_max, Unit()));
      }
      plot_subset2.get_subset(plot_subset, p0, p1);
      if(SelectedPlot == 2 || SelectedPlot == 4) {
	if(nr_noise_patches > 0) {
	  for(i = 0; i < nr_noise_patches; i++) {
	    for(xi = f2npatch_min[i]; xi <= f2npatch_max[i]; xi++) {
	      for(yi = f3npatch_min[i]; yi <= f3npatch_max[i]; yi++) {
		plot_subset.data[yi*plot_subset.axes[1].nelem+xi] = -10*fabs(twodsf_Imin);
	      }
	    }
	  }
	}
      }
  */

  if(SelectedPlot == 2 || SelectedPlot == 4) {
    if(nr_noise_patches > 0) {
      for(i = 0; i < nr_noise_patches; i++) {
	for(xi = f2npatch_min[i]; xi <= f2npatch_max[i]; xi++) {
	  for(yi = f3npatch_min[i]; yi <= f3npatch_max[i]; yi++) {
	    offset = yi*noise.NrBins+xi;
	    if(offset < 0 || offset >= noise.NrBins*noise.NrSubints) {
	      printerror(0, "Bug!!\n");
	    }else {
	      noise.data[offset] = -10*fabs(twodsf_Imin);
	    }
	  }
	}
      }
    }
  }
}

void PlotWindow(verbose_definition verbose)
{
  char txt[100];
  double Imax;
  double I;
  int i, xi, yi;
//START REGION DEVELOP
  double Imin;
  int NrPoints, nrpoints_flagged;
  int NrBins;
  long Bins[MaxNrNoiseBins], MaxBinValue;
//START REGION RELEASE

  /*  double x, y; */
//START REGION DEVELOP
//  if(AutoCalculateCentroid == 0) {
//START REGION RELEASE
    ppgpage();
    if(SelectedPlot == 0) {
      ppgsvp(0.1, 0.9, 0.1, 0.9);
      ppgswin(f2_min,f2_max,f3_min,f3_max);
      sprintf(txt, "2dfs feature component %d", SelectedComponent); 
      /*      ppglab("Fluctuation frequency (cycles/period)", "Fluctuation frequency (cycles/period)", txt); */
      /*      pgplot_colormap(plot_subset, cmap, 0.0, 0.0, Unit(), Unit()); */
      pgplot_options_definition pgplot_options;
      pgplot_clear_options(&pgplot_options);
      pgplot_options.viewport.dontopen = 1;
      pgplot_options.viewport.dontclose = 1;
      strcpy(pgplot_options.box.xlabel, "Fluctuation frequency (cycles/period)");
      strcpy(pgplot_options.box.ylabel, "Fluctuation frequency (cycles/period)");
      strcpy(pgplot_options.box.title, txt);
      /*      pgplotMap(&pgplot_options, twodfs.data, twodfs.NrBins, twodfs.NrSubints, -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins, +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins, f2_min, f2_max, 0, 0.5, f3_min, f3_max, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, "Fluctuation frequency (cycles/period)", "Fluctuation frequency (cycles/period)", txt, 1, 1, 1, -0.5, -0.5, 1.0, 1.0, 1, 1, 1, 1, 1, 0, 0, 0, application.verbose_state); */
      pgplotMap(&pgplot_options, twodfs.data, twodfs.NrBins, twodfs.NrSubints, -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins, +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins, f2_min, f2_max, 0, 0.5, f3_min, f3_max, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, verbose); 
      //      printf("XXXXX %f %f\n", -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins, +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)twodfs.NrBins);
      //	printf("XXXXX prof.nrbins=%ld\n", AverageProfile.NrBins);
      //	printf("XXXXX 2dfs.nrbins=%ld\n", twodfs.NrBins);
      if(centroid_calculated) {
	ppgsci(2);
	ppgslw(1);
	ppgmove(centroid_x,0.5);
	ppgdraw(centroid_x, -0.5);
        ppgmove(f2_min,centroid_y);
        ppgdraw(f2_max,centroid_y);  
        ppgmove(centroid_x-centroid_err_x,centroid_y-centroid_err_y);
        ppgdraw(centroid_x+centroid_err_x,centroid_y-centroid_err_y);
        ppgdraw(centroid_x+centroid_err_x,centroid_y+centroid_err_y);
        ppgdraw(centroid_x-centroid_err_x,centroid_y+centroid_err_y);
        ppgdraw(centroid_x-centroid_err_x,centroid_y-centroid_err_y);
        ppgsci(1);
      }
      /*      ppgbox("bcnsti",0.0,0,"bcnmsti",0.0,0); */
      /* Reset the viewport because of the side panels, so clicking works */
      /*      ppgswin(f2_min,f2_max,f3_min,f3_max); */

    }else if(SelectedPlot == 1) {
      ppgsvp(0.1, 0.9, 0.1, 0.9);
      ppgswin(fl_min,fl_max,f3_min,f3_max);
      /*      ppglab("Pulse longitude (degrees)", "Fluctuation frequency (cycles/period)", "lrfs feature"); */
      /*      pgplot_colormap(plot_subset, cmap, 0.0, 0.0, Unit("deg"), Unit()); */
      pgplot_options_definition pgplot_options;
      pgplot_clear_options(&pgplot_options);
      pgplot_options.viewport.dontopen = 1;
      pgplot_options.viewport.dontclose = 1;
      strcpy(pgplot_options.box.xlabel, "bins");
      strcpy(pgplot_options.box.ylabel, "Fluctuation frequency (cycles/period)");
      strcpy(pgplot_options.box.title, "lrfs feature");
      /*      pgplotMap(&pgplot_options, lrfs.data, lrfs.NrBins, lrfs.NrSubints, 0, lrfs.NrBins-1, fl_min, fl_max, 0, 0.5, f3_min, f3_max, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, "bins", "Fluctuation frequency (cycles/period)", "lrfs feature", 1, 1, 1, -0.5, -0.5, 1.0, 1.0, 1, 1, 1, 1, 1, 0, 0, 0, verbose); */
      pgplotMap(&pgplot_options, lrfs.data, lrfs.NrBins, lrfs.NrSubints, 0, lrfs.NrBins-1, fl_min, fl_max, 0, 0.5, f3_min, f3_max, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, verbose); 
      if(PlotAvrgProfile != 0) {
        i = 0;
	Imax = AverageProfile.data[0];
	for(xi=0; xi < AverageProfile.NrBins; xi++) {
	  if(AverageProfile.data[xi] > Imax)
	    Imax = AverageProfile.data[xi];
	}
	for(xi=0; xi < AverageProfile.NrBins; xi++) {
	  ppgsci(2
//START REGION DEVELOP
+WhichRegion(xi)
//START REGION RELEASE
);
	  if(xi >= fl_min && xi <= fl_max) {
	    I = AverageProfile.data[xi]*(f3_max-f3_min)/Imax + f3_min;
	    if(I < f3_min)
	      I = f3_min;
	    if(i == 0) {
	      ppgmove(xi,I);
	      i = 1;
	    }else {
	      ppgdraw(xi,I);
	    }
	  }
	}
	ppgsci(1);
      }
      if(centroid_calculated) {  
        ppgsci(2);
        ppgmove(0,centroid_y);
        ppgdraw(AverageProfile.NrBins,centroid_y);  
        ppgsci(1);
      }
      /*      ppgbox("bcnsti",0.0,0,"bcnmsti",0.0,0); */
    }else if(SelectedPlot == 2) {
      /*
      ppgsvp(0.1, 0.9, 0.1, 0.9);
      ppgswin(f2n_min,f2n_max,f3n_min,f3n_max);
      ppglab("Fluctuation frequency (cycles/period)", "Fluctuation frequency (cycles/period)", "2dfs noise");
      */

      for(xi = 0; xi < noise.NrBins; xi++) {
        for(yi = 0; yi < noise.NrSubints; yi++) {
          I = noise.data[yi*noise.NrBins+xi];
          if(I < twodsf_Imin)
            noise.data[yi*noise.NrBins+xi] = 0;
        }    
      }
      pgplot_options_definition pgplot_options;
      pgplot_clear_options(&pgplot_options);
      pgplot_options.viewport.dontopen = 1;
      pgplot_options.viewport.dontclose = 1;
      strcpy(pgplot_options.box.xlabel, "Fluctuation frequency (cycles/period)");
      strcpy(pgplot_options.box.ylabel, "Fluctuation frequency (cycles/period)");
      strcpy(pgplot_options.box.title, "For noise calulation: All signal should be flagged in this 2dfs plot");
      /*      pgplotMap(&pgplot_options, noise.data, noise.NrBins, noise.NrSubints, 
	      -AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)noise.NrBins, +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)noise.NrBins, f2_min, f2_max, 0, 0.5, f3_min, f3_max, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, "Fluctuation frequency (cycles/period)", "Fluctuation frequency (cycles/period)", "2dfs noise", 1, 1, 1, -0.5, -0.5, 1.0, 1.0, 1, 1, 1, 1, 1, 0, 0, 0, verbose); */
      pgplotMap(&pgplot_options, noise.data, noise.NrBins, noise.NrSubints, 
		-AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)noise.NrBins, +AverageProfile.NrBins/2.0-0.5*AverageProfile.NrBins/(float)noise.NrBins, f2_min, f2_max, 0, 0.5, f3_min, f3_max, PPGPLOT_INVERTED_HEAT, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, verbose); 
      SelectRegion();
      /*      ppgbox("bcnsti",0.0,0,"bcnsti",0.0,0); */
//START REGION DEVELOP
    }else if(SelectedPlot == 3) {
      /*
      ppgsvp(0.1, 0.9, 0.6, 0.9); 
      SelectedPlot = 0;
      SelectRegion();

      for(xi = 0; xi < plot_subset.axes[1].nelem; xi++) {
        I = 0;
        for(yi = 0; yi < plot_subset.axes[0].nelem; yi++) {
         I += plot_subset.data[yi*plot_subset.axes[1].nelem+xi];
        }
        if(xi == 0)
          Imin = Imax = I;
      if(I > Imax)
	Imax = I;
      if(I < Imin)
	Imin = I;
      }
      ppgswin(f2_min,f2_max,Imin,Imax);
      for(xi = 0; xi < plot_subset.axes[1].nelem; xi++) {
        I = 0;
        for(yi = 0; yi < plot_subset.axes[0].nelem; yi++) {
	  I += plot_subset.data[yi*plot_subset.axes[1].nelem+xi];
        }
        x = plot_subset.axes[1].bin_center(xi).get(Unit());
        if(xi == 0)
          ppgmove(x, I);
        else
          ppgdraw(x,I);  
      }
      if(centroid_calculated) {
        ppgsci(2);
        ppgmove(centroid_x, Imax);
        ppgdraw(centroid_x, Imin);
        ppgmove(centroid_x+centroid_err_x, Imax);
        ppgdraw(centroid_x+centroid_err_x, Imin);
        ppgmove(centroid_x-centroid_err_x, Imax);
        ppgdraw(centroid_x-centroid_err_x, Imin);
        ppgsci(1);
      }
      ppglab("1/P\\d2\\u (cycles/period)", "Strength [A.U.]", "P\\d2\\u distribution");
      ppgbox("bcnsti",0.0,0,"bcti",0.0,0);
  
      ppgsvp(0.1, 0.9, 0.1, 0.4); 
      for(yi = 0; yi < plot_subset.axes[0].nelem; yi++) {
        I = 0;
        for(xi = 0; xi < plot_subset.axes[1].nelem; xi++) {
  	  I += plot_subset.data[yi*plot_subset.axes[1].nelem+xi];
        }
        if(yi == 0)
          Imin = Imax = I;
        if(I > Imax)
    	  Imax = I;
        if(I < Imin)
	  Imin = I;
      }
      ppgswin(f3_min,f3_max,Imin,Imax);
      for(yi = 0; yi < plot_subset.axes[0].nelem; yi++) {
        I = 0;
        for(xi = 0; xi < plot_subset.axes[1].nelem; xi++) {
  	  I += plot_subset.data[yi*plot_subset.axes[1].nelem+xi];
        }
        x = plot_subset.axes[0].bin_center(yi).get(Unit());
        if(yi == 0)
          ppgmove(x, I);
        else
          ppgdraw(x, I);  
      }
      if(centroid_calculated) {
        ppgsci(2);
        ppgmove(centroid_y, Imax);
        ppgdraw(centroid_y, Imin);
        ppgmove(centroid_y+centroid_err_y, Imax);
        ppgdraw(centroid_y+centroid_err_y, Imin);
        ppgmove(centroid_y-centroid_err_y, Imax);
        ppgdraw(centroid_y-centroid_err_y, Imin);
        ppgsci(1);
      }
      ppglab("1/P\\d3\\u (cycles/period)", "Strength [A.U.]", "P\\d3\\u distribution");
      ppgbox("bcnsti",0.0,0,"bcti",0.0,0);
  
      SelectedPlot = 3;
      */
      fprintf(stderr, "PLOT 3 DISABLED FOR NOW\n");
    }else if(SelectedPlot == 4) {
      /*
      ppgsvp(0.1, 0.9, 0.6, 0.9);
      NrPoints = plot_subset.axes[1].nelem*plot_subset.axes[0].nelem;
      nrpoints_flagged = 0;
      Imin = Imax = 0;
      for(xi = 0; xi < plot_subset.axes[1].nelem; xi++) {
        for(yi = 0; yi < plot_subset.axes[0].nelem; yi++) {
  	  I = plot_subset.data[yi*plot_subset.axes[1].nelem+xi];
    	  if(I < twodsf_Imin) {
	    nrpoints_flagged++;
 	  }else if(I > Imax) {
	    Imax = I;
	  }else if(I < Imin) {
	    Imin = I;
  	  }
        }    
      }
      if(Imax < fabs(Imin))
        Imax = fabs(Imin);
      I = 0.01*(NrPoints - nrpoints_flagged);
      NrBins = (int)I;
      NrBins /= 2;
      NrBins *= 2;
      NrBins += 1;
      if(NrBins > MaxNrNoiseBins)
        NrBins = MaxNrNoiseBins;
      for(i = 0; i < NrBins; i++)
        Bins[i] = 0;
      MaxBinValue = 0;
      for(xi = 0; xi < plot_subset.axes[1].nelem; xi++) {
        for(yi = 0; yi < plot_subset.axes[0].nelem; yi++) {
  	  I = plot_subset.data[yi*plot_subset.axes[1].nelem+xi];
  	  if(I >= twodsf_Imin) {
	    I = NrBins*(I + Imax)/(2.0*Imax);
	    i = (int)I;
	    if(i == NrBins)
	      i = NrBins - 1;
	    Bins[i]++;
	    if(Bins[i] > MaxBinValue)
	      MaxBinValue = Bins[i];
	  }
        }
      } 
      ppgswin(-0.5,NrBins-0.5,0,log((double)MaxBinValue));
      ppgsci(2);
      ppgmove(-0.5, log((double)Bins[0]));
      ppgdraw(0.5, log((double)Bins[0]));
      for(i = 1; i < NrBins; i++) {
        ppgdraw(i-0.5, log((double)Bins[i]));
        ppgdraw(i+0.5, log((double)Bins[i]));
      }
      ppgsci(3);
      ppgmove(NrBins - NrBins/2 - 1, 0);
      ppgdraw(NrBins - NrBins/2 - 1, log((double)MaxBinValue));
      ppgsci(1);
      ppglab("bin number", "ln(Nr noise samples)", "Noise distribution");
      ppgbox("bcnsti",0.0,0,"bcti",0.0,0);
*/
      ppgsvp(0.1, 0.9, 0.6, 0.9);
      NrPoints = noise.NrBins * noise.NrSubints;
      nrpoints_flagged = 0;
      Imin = Imax = 0;
      for(xi = 0; xi < noise.NrBins; xi++) {
        for(yi = 0; yi < noise.NrSubints; yi++) {
  	  I = noise.data[yi*noise.NrBins+xi];
    	  if(I < twodsf_Imin) {
	    nrpoints_flagged++;
 	  }else if(I > Imax) {
	    Imax = I;
	  }else if(I < Imin) {
	    Imin = I;
  	  }
        }    
      }
      if(Imax < fabs(Imin))
        Imax = fabs(Imin);
      I = 0.01*(NrPoints - nrpoints_flagged);
      NrBins = (int)I;
      NrBins /= 2;
      NrBins *= 2;
      NrBins += 1;
      if(NrBins > MaxNrNoiseBins)
        NrBins = MaxNrNoiseBins;
      for(i = 0; i < NrBins; i++)
        Bins[i] = 0;
      MaxBinValue = 0;
      for(xi = 0; xi < noise.NrBins; xi++) {
        for(yi = 0; yi < noise.NrSubints; yi++) {
  	  I = noise.data[yi*noise.NrBins+xi];
  	  if(I >= twodsf_Imin) {
	    I = NrBins*(I + Imax)/(2.0*Imax);
	    i = (int)I;
	    if(i == NrBins)
	      i = NrBins - 1;
	    Bins[i]++;
	    if(Bins[i] > MaxBinValue)
	      MaxBinValue = Bins[i];
	  }
        }
      } 
      ppgswin(-0.5,NrBins-0.5,0,log((double)MaxBinValue));
      ppgsci(2);
      ppgmove(-0.5, log((double)Bins[0]));
      ppgdraw(0.5, log((double)Bins[0]));
      for(i = 1; i < NrBins; i++) {
        ppgdraw(i-0.5, log((double)Bins[i]));
        ppgdraw(i+0.5, log((double)Bins[i]));
      }
      ppgsci(3);
      ppgmove(NrBins - NrBins/2 - 1, 0);
      ppgdraw(NrBins - NrBins/2 - 1, log((double)MaxBinValue));
      ppgsci(1);
      ppglab("bin number", "ln(Nr noise samples)", "Noise distribution");
      ppgbox("bcnsti",0.0,0,"bcti",0.0,0);
//START REGION RELEASE
    }
//START REGION DEVELOP
  //}
//START REGION RELEASE
}

int MakeSelection(double sigma_noise)
{
  float x0, x1, y0, y1, dummy;
  char c;

  x0 = y0 = x1 = y1 = c = 0;

  KeyCode = 0;
  ppgsci(2);
  ppgband(0, 0, 0.0, 0.0, &x0, &y0, &c);
  if(c != 65) {
    ppgsci(1);
    KeyCode = c;
    return 0;
  }
  if(SelectedPlot == 3) {
    y0 = x0;
    x0 = (y0 - f3_min)/(f3_max-f3_min)*(f2_max-f2_min)+f2_min;
    printf("P3[P0]  = %lf\n", 1/y0);
    printf("P2[deg] = %lf\n", 360/x0);
    KeyCode = 0;
    ppgsci(1);
    return 0;
  }
  ppgband(2, 0, x0, y0, &x1, &y1, &c); 
  ppgsci(1);
  if(c != 65)
    return 0;

  if(y0 > y1) {
    dummy = y0;
    y0 = y1;
    y1 = dummy;
  }
  if(x0 > x1) {
    dummy = x0;
    x0 = x1;
    x1 = dummy;
  }

  if(SelectedPlot == 0) {
    if(Centered) {
      if(fabs(x0) > fabs(x1)) {
	x1 = fabs(x0);
	x0 = -fabs(x0);
      }else {
	x0 = -fabs(x1);
	x1 = fabs(x0);
      }
    }
    f2_min = x0;
    f2_max = x1;
    f3_min = y0;
    f3_max = y1;
    int nx, ny;
    float dx, dy;
    pgplotMapCoordinate_dbl(f2_min, f3_min, &nx, &ny);
    pgplotMapCoordinateInverse_dbl(&f2_min, &f3_min, nx, ny);
    pgplotMapCoordinate_dbl(f2_max, f3_max, &nx, &ny);
    pgplotMapCoordinateInverse_dbl(&f2_max, &f3_max, nx, ny);
    pgplotMapCoordinateBinSize(&dx, &dy);
    f2_min -= 0.49*dx;
    f2_max += 0.49*dx;
    f3_min -= 0.49*dy;
    f3_max += 0.49*dy;
  }else if(SelectedPlot == 1) {
    fl_min = x0;
    fl_max = x1;
    f3_min = y0;
    f3_max = y1;
  }else if(SelectedPlot == 2) {
    if(nr_noise_patches == Max_nr_noise_patches) {
      printf("Too many patches\n");
      nr_noise_patches--;
    }
    pgplotMapCoordinate(x0, y0, &(f2npatch_min[nr_noise_patches]), &(f3npatch_min[nr_noise_patches]));
    pgplotMapCoordinate(x1, y1, &(f2npatch_max[nr_noise_patches]), &(f3npatch_max[nr_noise_patches]));
    nr_noise_patches++;
    SelectRegion();
//START REGION DEVELOP
    if(NoiseCalculation)
//START REGION RELEASE
      sigma_noise = calculate_noise_sigma();
//START REGION DEVELOP
    else
      sigma_noise = 0;
//START REGION RELEASE
  }
  return 1;
}


double calculate_noise_sigma(void)
{
  double I, sigma_noise;
  int xi, yi;
  int nrpoints_flagged;
  sigma_noise = 0;
  nrpoints_flagged = 0;
  for(yi = 0; yi < noise.NrSubints; yi++) {
    for(xi = 0; xi < noise.NrBins; xi++) {
      I = noise.data[yi*noise.NrBins+xi];
      if(I < twodsf_Imin)
        nrpoints_flagged++;
      else
        sigma_noise += I*I;
    }
  }
  sigma_noise = sqrt(sigma_noise/(double)(noise.NrBins*noise.NrSubints-nrpoints_flagged));
  printf("Sigma noise = %e (%d points flagged)\n", sigma_noise, nrpoints_flagged);
  if(nrpoints_flagged == 0) {
    printf("\nNo signal is flagged, so the current error-bar is based on the rms of the noise + pulsar signal, overestimating the actual error. Press 'h' for a general help.\n\n");
  }
  return sigma_noise;
}

void calculate_2dfs_Centroid(double sigma_noise, verbose_definition verbose)
{
  //  fprintf(stderr, "calculate_2dfs_Centroid disabled for now\n");
  /*
  double I, Itot;
  int xi, yi;
  double x, y, xcent, ycent, xerr, yerr, xerr2, yerr2;

  std::vector<IndexType> p0(2), p1(2);
  PhysDataArray<double, float> feature;
  PhysDataArray<double, float> plot_subset2;

  plot_subset2 = twodfs.get_hyperelem_ref(2, 0);
  p0[0] = plot_subset2.axes[0].ibin(PhysDataPoint<double>(f3_min, Unit()));
  p1[0] = plot_subset2.axes[0].ibin(PhysDataPoint<double>(f3_max, Unit()));
  p0[1] = plot_subset2.axes[1].ibin(PhysDataPoint<double>(f2_min, Unit()));
  p1[1] = plot_subset2.axes[1].ibin(PhysDataPoint<double>(f2_max, Unit()));
  plot_subset2.get_subset(feature, p0, p1);

  plot_subset2 = twodfs.get_hyperelem_ref(2, 0);
  p0[0] = plot_subset2.axes[0].ibin(PhysDataPoint<double>(f3n_min, Unit()));
  p1[0] = plot_subset2.axes[0].ibin(PhysDataPoint<double>(f3n_max, Unit()));
  p0[1] = plot_subset2.axes[1].ibin(PhysDataPoint<double>(f2n_min, Unit()));
  p1[1] = plot_subset2.axes[1].ibin(PhysDataPoint<double>(f2n_max, Unit()));

  if(application.verbose_state || AutoCalculateCentroid)
    printf("Centroid:      "); 
  xcent = ycent = Itot = 0;
  for(yi = feature.axes[0].nelem-1; yi >= 0; yi--) {
    for(xi = 0; xi < feature.axes[1].nelem; xi++) {
      x = feature.axes[1].bin_center(xi).get(Unit());
      y = feature.axes[0].bin_center(yi).get(Unit());
      I = feature.data[yi*feature.axes[1].nelem+xi];
      Itot += I;
      xcent += I*x;
      ycent += I*y;
    }
  }
  xcent /= Itot;
  ycent /= Itot;
  if(application.verbose_state || AutoCalculateCentroid) {
    printf("(%lf, %lf)\n", xcent, ycent);
    printf("Noise error:   "); 
  }
  xerr = yerr = xerr2 = yerr2 = 0;
  for(yi = feature.axes[0].nelem-1; yi >= 0; yi--) {
    for(xi = 0; xi < feature.axes[1].nelem; xi++) {
      x = feature.axes[1].bin_center(xi).get(Unit());
      y = feature.axes[0].bin_center(yi).get(Unit());
      I = feature.data[yi*feature.axes[1].nelem+xi];
      xerr += (x - xcent)*(x - xcent);
      yerr += (y - ycent)*(y - ycent);
      xerr2 += I*I;
    }
  }
  xerr *= sigma_noise*sigma_noise/(Itot*Itot);
  yerr *= sigma_noise*sigma_noise/(Itot*Itot);
  xerr2 /= (Itot*Itot);
  yerr2 = xerr2;
  xerr2 *= 0.25*feature.axes[1].bin_size.get(Unit())*feature.axes[1].bin_size.get(Unit());
  yerr2 *= 0.25*feature.axes[0].bin_size.get(Unit())*feature.axes[0].bin_size.get(Unit());
  xerr = sqrt(xerr);
  yerr = sqrt(yerr);
  xerr2 = sqrt(xerr2);
  yerr2 = sqrt(yerr2);
  if(application.verbose_state || AutoCalculateCentroid) {
    printf("(%e, %e)\n", xerr, yerr);
    printf("Bin error:     (%e, %e)\n", xerr2, yerr2);
    printf("Half bin size: (%e, %e)\n", 0.5*feature.axes[1].bin_size.get(Unit()), 0.5*feature.axes[0].bin_size.get(Unit()));
  }
  centroid_x = xcent;
  centroid_y = ycent;
  centroid_calculated = 1;
  centroid_err_x = xerr; 
  centroid_err_y = yerr; 
  if(!AutoCalculateCentroid || application.verbose_state) {
    printf("P3[P0]  = %lf +- %lf\n", 1/ycent, (centroid_err_y)/(ycent*ycent));
    printf("P2[deg] = %lf +- %lf\n", 360/xcent, 360*(centroid_err_x)/(xcent*xcent));
  }
*/


  double I, Itot;
  int xi, yi, xstart, xend, ystart, yend, nrbins;
  double x, y, xcent, ycent, xerr, yerr;
//START REGION DEVELOP
  double xerr2, yerr2;
//START REGION RELEASE
  float xf, yf, binsizex, binsizey;


//START REGION DEVELOP
  if(AutoCalculateCentroid)
    verbose.verbose = 1;
//START REGION RELEASE

  if(SelectedPlot == 0) {    /* Only works when 2dfs is currently plotted, otherwise pgplotMapCoordinate isn't working */
    pgplotMapCoordinate(f2_min, f3_min, &xstart, &ystart);
    pgplotMapCoordinate(f2_max, f3_max, &xend, &yend);
    
    
    //    printf("XXXXX considering bins %d .. %d\n", xstart, xend);
    //    printf("XXXXX considering bins %d .. %d\n", ystart, yend);

    if(verbose.verbose)
      printf("Centroid:                            "); 
    xcent = ycent = Itot = 0;
    nrbins = 0;
    for(yi = ystart; yi <= yend; yi++) {
      for(xi = xstart; xi <= xend; xi++) {
	I = fabs(twodfs.data[yi*twodfs.NrBins+xi]);
	pgplotMapCoordinateInverse(&xf, &yf, xi, yi);
	//	if(yi == 0)
	//	  printf("XXXXXX xf=%f I=%f\n", xf, I);
	x = xf;
	y = yf;
	Itot += I;
	xcent += I*x;
	ycent += I*y;
	nrbins++;
      }
    }
    xcent /= Itot;
    ycent /= Itot;
    if(verbose.verbose) {
      printf("(%lf, %lf) cpp", xcent, ycent);
      printf(" based on %d selected bins\n", nrbins);
      printf("Statistical error caused by noise:   "); 
    }
    xerr = yerr = 0;
//START REGION DEVELOP
    xerr2 = yerr2 = 0;
//START REGION RELEASE
    for(yi = ystart; yi <= yend; yi++) {
      for(xi = xstart; xi <= xend; xi++) {
	I = twodfs.data[yi*twodfs.NrBins+xi];
	pgplotMapCoordinateInverse(&xf, &yf, xi, yi);
	x = xf;
	y = yf;
	xerr += (x - xcent)*(x - xcent);
	yerr += (y - ycent)*(y - ycent);
//START REGION DEVELOP
	xerr2 += I*I;
//START REGION RELEASE
      }
    }
    xerr *= sigma_noise*sigma_noise/(Itot*Itot);
    yerr *= sigma_noise*sigma_noise/(Itot*Itot);
//START REGION DEVELOP
    xerr2 /= (Itot*Itot);
    yerr2 = xerr2;
//START REGION RELEASE
    pgplotMapCoordinateBinSize(&binsizex, &binsizey);

//START REGION DEVELOP
    xerr2 *= 0.25*binsizex*binsizex;
    yerr2 *= 0.25*binsizey*binsizey;
//START REGION RELEASE
    xerr = sqrt(xerr);
    yerr = sqrt(yerr);
//START REGION DEVELOP
    xerr2 = sqrt(xerr2);
    yerr2 = sqrt(yerr2);
//START REGION RELEASE
    if(verbose.verbose) {
      printf("(%e, %e) cpp\n", xerr, yerr);
//START REGION DEVELOP
      printf("Bin error:                           (%e, %e)\n", xerr2, yerr2);
//START REGION RELEASE
      printf("Half bin size:                       (%e, %e) cpp\n", 0.5*binsizex, 0.5*binsizey);
    }
    centroid_x = xcent;
    centroid_y = ycent;
    centroid_calculated = 1;
    centroid_err_x = xerr; 
    centroid_err_y = yerr; 
//START REGION DEVELOP
    if(!AutoCalculateCentroid || verbose.verbose) {
//START REGION RELEASE
      printf("P3[P0]  = %lf +- %lf\n", 1/ycent, (centroid_err_y)/(ycent*ycent));
      printf("P2[deg] = %lf +- %lf\n", 360/xcent, 360*(centroid_err_x)/(xcent*xcent));
//START REGION DEVELOP
    }
//START REGION RELEASE
  }
}



//START REGION DEVELOP


int WhichRegion(int bin)
{
  int i;
  for(i = 0; i < NrSelectedPulseRegions; i++) {
    if(bin >= PulseRegions[i][0] && bin <= PulseRegions[i][1])
      return i+1;
  }
  return 0;
}

//START REGION DEVELOP

