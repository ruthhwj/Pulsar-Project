//START REGION RELEASE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "psrsalsa.h"

#define MaxNrJumps     100

//START REGION DEVELOP
void SHOWREVISIONINFO_prog() {
#include "ppol.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);

//START REGION RELEASE
int main(int argc, char **argv)
{
  int index, firstfiletoopen, nokeypresses, compute_PA, extendedPol, writeout, writeoutFilename;
  int manualOnpulseSelection, deviceID, normalize, correctLbias, nxsub, boxlinewidth;
  int titleset, titlelw, dashed, datalinewidth, noylabel, overlayPA, overlayFine, nrJumps;
  int outline, outline_color, extended, onlysignificantPA, twoprofiles, sigmaI;
  int doprojection;
  long i, j, f;
  float *profileI, loffset, correctQV, correctV, sigma_limit;
  float xsize, ysize, ysizepa, xtick, plotI1, plotI2, titlech;
  float plotl1, plotl2, plotp1, plotp2, PAoffset_data;
  float overlayalpha, overlaybeta, overlaypa0, overlayl0;
  float jump_longitude[MaxNrJumps], jump_offset[MaxNrJumps];
  char *filename_ptr, title[1000], txt[100], PlotDevice[100], PlotDevice2[100], *extension;
  char ofilename[1000];
  datafile_definition datain, dataout;
  psrsalsaApplication application;
  pgplot_options_definition pgplot_options;
//START REGION DEVELOP
  int subtract, subtractFile, elldist, sigmaset;
  int pa_dist_normalise, pa_dist_stat;
  int decompose_method, doselectmode;
  int show_pa_dist, projection_binnr, projection_weighting, proj_type, projection_nolabels;
  int fit_pa_dist_maxNrComp, fit_pa_dist_precision;
  float fit_pa_dist_sigmalimit, pamask_value, offpulsedata_rebinfactor, selectmode_cpa, selectmode_dpa, selectmode_cell, selectmode_dell;
  float *mapProjection, proj_rot_long, proj_rot_lat;
  int projection_nrx, projection_nry;
  datafile_definition padist_data, padist_fit, padist_resid, subtract_fin, pamask, data_offpulse;
  datafile_definition mode1, mode2;
  char pamask_name[MaxFilenameLength], offpulsedata_name[MaxFilenameLength];
//START REGION DEVELOP




//START REGION RELEASE
  initApplication(&application, "ppol", "[options] inputfile(s)");
//START REGION DEVELOP
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
//START REGION RELEASE
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_iformat = 1;
  application.switch_oformat = 1;
  application.switch_formatlist = 1;
  application.switch_header = 1;
  application.switch_headerlist = 1;
  application.switch_onpulse = 1;
  application.switch_onpulsef = 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_rebin = 1;
  application.switch_TSCR = 1;
  application.switch_FSCR = 1;
  application.switch_debase = 1;
  application.switch_rotateStokes = 1;
  application.switch_nocounters = 1;
  application.switch_stokes = 1;
  application.switch_deparang = 1;
  application.switch_deFaraday = 1;
  application.switch_filelist = 1;
  application.switch_dedisperse = 1;
  application.switch_deFaraday = 1;
  application.switch_changeRefFreq = 1;
  application.switch_fchan = 1;
  application.switch_history_cmd_only = 1;
  application.switch_nskip = 1;
  application.switch_nread = 1;
  application.switch_forceUniformFreqLabelling = 1;
//START REGION DEVELOP
  application.switch_rotateQU = 1;
  application.switch_ppgplot = 1;
  application.switch_showrevision = 1;
  application.switch_autoonpulse = 1;
  application.switch_subtractAvPA = 1;

//START REGION RELEASE
  application.oformat = PPOL_format;  // Default format for output
  nokeypresses = 0;
  extendedPol = 0;
  writeout = 0;
  writeoutFilename = -1;
  strcpy(PlotDevice, "?");
  strcpy(PlotDevice2, "?");
  manualOnpulseSelection = 0; // 0 = only if none specified, 1 = even if specified, -1 even if not specified
  extension = NULL;
  normalize = 1;
  correctLbias = 1;   /* Default: do the correct L debias method */
  correctQV = 1;
  correctV = 1;
  loffset = 0;
  xsize = 0.8;
  ysize = 1.0;
  ysizepa = 1.0;
  xtick = 0.0;
  nxsub = 0;
  sigma_limit = 3.0;
  boxlinewidth = 1;
  titleset = 0;
  title[0] = 0;
  titlech = 1;
  titlelw = 1;
  plotl1 = 0;
  plotl2 = 360;
  plotp1 = -180;
  plotp2 = 180;
  plotI1 = plotI2 = 0;
  PAoffset_data = 0;
  datalinewidth = 1;
  dashed = 0;
  noylabel = 0;
  overlayPA = 0;
  overlayFine = 1;
  outline = -1;
  nrJumps = 0;
  twoprofiles = 0;
  outline_color = 0;
  overlayalpha = 0;
  overlaybeta = 0;
  overlaypa0 = 0;
  overlayl0 = 0;
  sigmaI = 0;
  doprojection = 0;
//START REGION DEVELOP
  show_pa_dist = 0;
  subtract = 0;
  subtractFile = 0;
  pa_dist_normalise = 0;
  pa_dist_stat = 0;
  fit_pa_dist_maxNrComp = 2;
  fit_pa_dist_precision = 2;
  fit_pa_dist_sigmalimit = 5.0;
  pamask_name[0] = 0;
  elldist = 0;
  projection_nolabels = 0;
  offpulsedata_name[0] = 0;
  decompose_method = 0;
  sigmaset = 0;
  doselectmode = 0;

//START REGION RELEASE
  if(argc < 2) {
    printf("Program to convert Stokes parameters in other polarization products such as\n");
    printf("postition angle and linear intensity. The results can be written out and plotted\n");
    printf("using ppolFig or fitted using ppolFit.\n\n");
    printApplicationHelp(&application);
    //    fprintf(stdout, "-oformat           %d = PPOL_format, %d = PPOL_SHORT_format.\n", PPOL_format, PPOL_SHORT_format);
    fprintf(stdout, "General input/output options:\n");
    fprintf(stdout, "-ext               Write out polarised profile to file with this extension.\n");
    fprintf(stdout, "-stdout            Like -ofile, but now write to stdout.\n");
    fprintf(stdout, "-ofile             Like -stdout, but now write to the specified file.\n");

    fprintf(stdout, "\nOptions related to generating polarimetric data:\n");
//START REGION DEVELOP
    fprintf(stdout, "-addfile           Add the PA-swing from this datafile to the data.\n");
    fprintf(stdout, "                   This is the opposite from -subtractfile.\n");
    fprintf(stdout, "-decompose method  Decompose the polarimetric data into modes according to the\n");
    fprintf(stdout, "                   specified method, where method could be:\n");
    fprintf(stdout, "                   1 = split in two orthogonal modes\n");
    fprintf(stdout, "                   2 = split in two equal amplitude modes\n");
    fprintf(stdout, "                   3 = force one mode to follow the -paswing option\n");
    fprintf(stdout, "                   The data will be written to a filename equal to the input\n");
    fprintf(stdout, "                   filename, with a .mode1 and .mode2 extension appended.\n");
    //  -n          Nr of steps to use to resolve the different two modes before finding unique solution.
    //  -subV       Subtract Stokes V from Stokes I.
    //  -paswing    Overlay PA-swing with this alpha, beta, PA0 and L0.
//START REGION RELEASE
    fprintf(stdout, "-loffset           Shift longitudes by this amount.\n");
    fprintf(stdout, "-medianLdebias     Naively subtract the median L of offpulse region rather than\n");
    fprintf(stdout, "                   applying the better de-bias method of Wardle & Kronberg.\n");
    fprintf(stdout, "-noLdebias         Simply use L^2=Q^2+U^2, which is biased.\n");
    fprintf(stdout, "-nonorm            Do not normalize the output profile.\n");
//START REGION DEVELOP
    fprintf(stdout, "-offpulsedata      \"file rebinfactor\"\n");
    fprintf(stdout, "                   The specified file will be used to determine the noise\n");
    fprintf(stdout, "                   characteristics. This can be useful if the actual data to be\n");
    fprintf(stdout, "                   used to generate the polarization information is re-binned,\n");
    fprintf(stdout, "                   thereby reducing the number of offpulse bins, hence reducing\n");
    fprintf(stdout, "                   the precision of the rms. Set rebinfactor to the ratio of the\n");
    fprintf(stdout, "                   sampling times of the actual data and the offpulse data file\n");
    fprintf(stdout, "                   (>1 if actual data is rebinned). Any selection of the onpulse\n");
    fprintf(stdout, "                   region should be related to the file to be used for the noise\n");
    fprintf(stdout, "                   calculation. It is assumed that the rebinning has been done in\n");
    fprintf(stdout, "                   such a way that the mean remains the same, as is the case with\n");
    fprintf(stdout, "                   pmod -rebin.\n");
//START REGION RELEASE
    fprintf(stdout, "-paoffset          Add this angle to PA (in degrees).\n");
//START REGION DEVELOP
    fprintf(stdout, "-selectmode        \"cpa dpa cell dell\"\n");
    fprintf(stdout, "                   Filter data such that only polarizations with a PA within the\n");
    fprintf(stdout, "                   range defined by the centre PA cpa and offset dpa are accepted.\n");
    fprintf(stdout, "                   Likewise the ellipticity should be in the range defined by the\n");
    fprintf(stdout, "                   central value cell with offset dell. All values are in degree.\n");
//START REGION RELEASE
    fprintf(stdout, "-selectonpulse     Enables manual graphical selection of more on-pulse regions\n");
    fprintf(stdout, "                   in addition to any provided on the command line.\n");
//START REGION DEVELOP
    fprintf(stdout, "-noselectonpulse   Disables manual graphical selection of on-pulse regions, even\n");
    fprintf(stdout, "                   if not provided on the cmd line (possibly useful with -ao).\n");
//START REGION RELEASE
    fprintf(stdout, "-sigma             Set sigma limit on L required for PA calculation [def=%.1f].\n", sigma_limit);
    fprintf(stdout, "                   Here sigma is defined as the RMS of the unbiased off-pulse L.\n");    
//START REGION DEVELOP
    fprintf(stdout, "                   For the ellipticity angle the RMS in total polarization is used.\n");
    fprintf(stdout, "-sigmaI            Use in combination with -sigma to make sigma the RMS of I instead.\n");
    fprintf(stdout, "                   In addition, SJ sets the offpulse region by selecting that region\n");
    fprintf(stdout, "                   that gives the lowest offpulse RMS I, thereby biasing the RMS down.\n");
    fprintf(stdout, "-subtract          Subtract the RVM model defined by the -paswing option\n");
    fprintf(stdout, "                   from the data.\n");
    fprintf(stdout, "-subtractfile      Subtract the PA-swing from this datafile from the data.\n");
    fprintf(stdout, "                   This is the opposite from -addfile.\n");
//START REGION RELEASE
    fprintf(stdout, "-2                 Write out two pulse periods (so there is duplicated data).\n");
    fprintf(stdout, "-extendedpol       Also generate the total amount of polarization\n");
//START REGION DEVELOP
    fprintf(stdout, "                   and ellipticity angles.\n");
    fprintf(stdout, "-correctQV         Multiply stokes Q and V by -1.\n");
    fprintf(stdout, "-correctV          Multiply stokes V by -1.\n");

    fprintf(stdout, "\nOptions related to generating distributions:\n");
    fprintf(stdout, "-padist        Calculates the PA distribution with this number of PA bins.\n");
    fprintf(stdout, "-elldist       When specified combined with -padist, those options\n");
    fprintf(stdout, "               now use the ellipticity angle rather than the PA.\n");
    fprintf(stdout, "-padist_norm   If specified, the PA distribution is normalised such that the\n");
    fprintf(stdout, "               pulse longitude bin column with most counts will have a sum of 1.\n");
    fprintf(stdout, "-padist_stat   \"maxc limit prec\"\n");
    fprintf(stdout, "               Use in combination with -padist or -padist_norm. The computed\n");
    fprintf(stdout, "               distribution is fitted with von Mises functions. The maximum nr\n");
    fprintf(stdout, "               of functions used in the fit is maxc (plus a baseline), although\n");
    fprintf(stdout, "               a smaller number of functions may be fitted. How many components\n");
    fprintf(stdout, "               are used is controlled by the value of limit. Increasing this\n");
    fprintf(stdout, "               value requires a more significant effect to accept another\n");
    fprintf(stdout, "               function, making it less likely noise is being fitted, but more\n");
    fprintf(stdout, "               likely components are missed. Increasing prec (an integer)\n");
    fprintf(stdout, "               increases the number of trials performed to try to find the best\n");
    fprintf(stdout, "               possible fit. If -ext of -ofile is used the fit to the PA\n");
    fprintf(stdout, "               distribution is written out to a file with additional\n");
    fprintf(stdout, "               extension .fit, and the residual (the difference between the data\n");
    fprintf(stdout, "               and the fit) to a file with extension .resid.\n");
    fprintf(stdout, "               The default options are \"%d %.1f %d\".\n", fit_pa_dist_maxNrComp, fit_pa_dist_sigmalimit, fit_pa_dist_precision);    
    fprintf(stdout, "-padist_mask   \"maskfile maskvalue\" (see below for 1 inputs version)\n");
    fprintf(stdout, "               Use in combination with -padist. In addition to generating the\n");
    fprintf(stdout, "               pa distribution (which is not written out to a file), the input\n");
    fprintf(stdout, "               data is written out (converted into I, L, V, PA etc.) after\n");
    fprintf(stdout, "               applying a mask. The provided maskfile should be very similar to\n");
    fprintf(stdout, "               a pa-distribution calculated with the -padist option, but has in\n");
    fprintf(stdout, "               the first \"polarization channel\" the counts replaced with a\n");
    fprintf(stdout, "               value that acts to define a mask. The result is that if the\n");
    fprintf(stdout, "               provided maskvalue = the value in the mask for a particular\n");
    fprintf(stdout, "               sample in the input data, then the input data is unchanged.\n");
    fprintf(stdout, "               Otherwise the sample is set to zero. The result is that samples\n");
    fprintf(stdout, "               in the input data can be zapped depending on their PA values.\n");
    fprintf(stdout, "               Alternative:\n");
    fprintf(stdout, "-padist_mask   \"maskfile\"\n");
    fprintf(stdout, "               Same as above, but the data in the input file is replaced with\n");
    fprintf(stdout, "               the values in maskfile for the given PA values in the data, or\n");
    fprintf(stdout, "               zero when the PA is not significant.\n");
    fprintf(stdout, "-projection    Shows a projection of the Poincare sphere. Provide\n");
    fprintf(stdout, "               \"binnr hsize weighting type dlong dlat\". If binnr is set, only\n");
    fprintf(stdout, "               that pulse longitude bin is considered. If -1 it only uses the\n");
    fprintf(stdout, "               onpulse region and if -2 it will use all bins. Set the width in\n");
    fprintf(stdout, "               pixels of the map with hsize. If weighting is set to 1, instead\n");
    fprintf(stdout, "               of counts the points are weighted by their polarized power.\n");
    fprintf(stdout, "               Type 1 is the (equal area) Hammer-Aitoff, type 2 is a normal\n");
    fprintf(stdout, "               spherical projection, 3 is a long/lat map. The orientation of\n");
    fprintf(stdout, "               the projection can be set by dlong and dlat. No S/N limit is\n");
    fprintf(stdout, "               imposed when the input is Stokes data. ppol -extendedpol can\n");
    fprintf(stdout, "               be used first to calculate PA's and ellipticity, which allows\n");
    fprintf(stdout, "               ppol -sigma to be used.\n");


//START REGION RELEASE
    fprintf(stdout, "\nOptions affecting the plotting:\n");
    fprintf(stdout, "-1             Only plot PA-swing once (equivalent to -yrange \"0 180\").\n");
//START REGION DEVELOP
    fprintf(stdout, "-boxlw         Set linewidth of the box [def=%d].\n", boxlinewidth);
    fprintf(stdout, "-dash          Make L and V profiles dashed/dotted.\n");
//START REGION RELEASE
    fprintf(stdout, "-device        Specify plotting device for onpulse region selection.\n");
    fprintf(stdout, "-device2       Specify plotting device final plot.\n");
//START REGION DEVELOP
    fprintf(stdout, "-herrorbar     \"xleft xcentre xright y SizeOfMarkers lineWidth colourIndex\"\n");
    fprintf(stdout, "               Draw a horizontal errorbar. The position y is normalised\n");
    fprintf(stdout, "               between 0 and 1.\n");
    fprintf(stdout, "-herrorbar2    Identical to to -herrorbar option, but now for the PA panel,\n");
    fprintf(stdout, "               except that y is in degrees rather than being normalized.\n");
    fprintf(stdout, "-verrorbar     \"x ybottom ycentre ytop SizeOfMarkers lineWidth colourIndex\"\n");
    fprintf(stdout, "               Draw a vertical errorbar. The position y is normalised\n");
    fprintf(stdout, "               between 0 and 1.\n");
    fprintf(stdout, "-verrorbar2    Identical to to -verrorbar option, but now for the PA panel,\n");
    fprintf(stdout, "               except that y is in degrees rather than being normalized.\n");
    fprintf(stdout, "-lw            Set linewidth used to plot data [def=%d].\n", datalinewidth);
    fprintf(stdout, "-nokeypress    Do not wait for key presses to go to next page in output plot\n");
    fprintf(stdout, "               (by default this shouldn't happen if plotting to a file)\n");
    fprintf(stdout, "-noylabel      Disable labeling and numbers along y-axis.\n");
    fprintf(stdout, "-opm           Put OPM at this longitude and with this amount of degrees when\n");
    fprintf(stdout, "               using -paswing. You can use this option multiple times.\n"); 
    fprintf(stdout, "-outline       \"lw color\" Makes text of -text option appear in outline\n");
    fprintf(stdout, "-paswing       \"alpha beta pa0 l0\" Overlay PA-swing with these parameters.\n");
    fprintf(stdout, "-projection_nobox  Do not show a rectangular box with numerical labels around the\n");    
    fprintf(stdout, "                   plot produced by the -projection option.\n");    
    fprintf(stdout, "-text          \"x y ch lw font color\" \"text\" plots text at this position in\n");
    fprintf(stdout, "               the top panel. See -textkeywords for possible keywords to use.\n");
    fprintf(stdout, "-textkeywords  Lists to keywords you can use in the -text and -title options.\n");
    fprintf(stdout, "-title \"...\"   Set title (default is file name). The title supports keywords\n");
    fprintf(stdout, "               listed by the -textkeywords option.\n");
    fprintf(stdout, "-title_fmt     \"characterheight line_width\" (default being \"%.1f %d\").\n", titlech, titlelw);
//START REGION RELEASE
    fprintf(stdout, "-xrange        Specify longitude range covered in the plot in degrees.\n");
    fprintf(stdout, "-xrange_phase  Specify longitude range covered in the plot in phase.\n");
    fprintf(stdout, "-yrange        Specify PA-range covered in the plot.\n");
//START REGION DEVELOP
    fprintf(stdout, "-yrange2       Specify intensity range covered in the plot.\n");
    fprintf(stdout, "-xsize         Set xsize of the plot [def=%.1f].\n", xsize);
    fprintf(stdout, "-ysize         Set ysize of the plot [def=%.1f].\n", ysize);
    fprintf(stdout, "-ysizepa       Set relative ysize of PA plot [def=%.1f].\n", ysizepa);
    fprintf(stdout, "-xticks        \"XTICK NXSUB\"  Adjust PGPLOT tickmarks [def=\"%.0f %d\"].\n", xtick, nxsub);
//START REGION RELEASE

    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about fitting position-angle swings and using the beam-width information can be found in:\n");
    printf(" - 	Rookyard et al. 2015, MNRAS, 446, 3367\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else if(argc >= 2) {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-device") == 0) {
	strcpy(PlotDevice,argv[i+1]);
        i++;
      }else if(strcmp(argv[i], "-device2") == 0) {
	strcpy(PlotDevice2,argv[i+1]);
        i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-correctQV") == 0) {
	correctQV = -1;
      }else if(strcmp(argv[i], "-correctV") == 0) {
	correctV = -1;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-stdout") == 0) {
	writeout = 1;
      }else if(strcmp(argv[i], "-nonorm") == 0) {
	normalize = 0;
      }else if(strcmp(argv[i], "-ofile") == 0) {
	writeout = 1;
	writeoutFilename = i+1;
	i++;
      }else if(strcmp(argv[i], "-ext") == 0) {
	writeout = 1;
	extension = argv[i+1];
        i++;
      }else if(strcmp(argv[i], "-paoffset") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &PAoffset_data, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-loffset") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &loffset, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-nokeypress") == 0) {
	nokeypresses = 1;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-extendedpol") == 0) {
	extendedPol = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-dash") == 0) {
	dashed = 1;
      }else if(strcmp(argv[i], "-noylabel") == 0) {
	noylabel = 1;
      }else if(strcmp(argv[i], "-decompose") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &decompose_method, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-projection") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d %d %d %f %f", &projection_binnr, &projection_nrx, &projection_weighting, &proj_type, &proj_rot_long, &proj_rot_lat, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	projection_nry = 0.5*projection_nrx;
	doprojection = 1;
        i++;
      }else if(strcmp(argv[i], "-selectmode") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f %f", &selectmode_cpa, &selectmode_dpa, &selectmode_cell, &selectmode_dell, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	doselectmode = 1;
        i++;
      }else if(strcmp(argv[i], "-projection_nobox") == 0) {
	projection_nolabels = 1;
      }else if(strcmp(argv[i], "-xsize") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &xsize, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-ysize") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &ysize, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-sigma") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &sigma_limit, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
//START REGION DEVELOP
	sigmaset = 1;
//START REGION RELEASE
        i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-sigmaI") == 0) {
	sigmaI = 1;
      }else if(strcmp(argv[i], "-ysizepa") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &ysizepa, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-title") == 0) {
	strcpy(title, argv[i+1]);
	titleset = 1;
        i++;
      }else if(strcmp(argv[i], "-title_fmt") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %d", &titlech, &titlelw, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
	/*
      }else if(strcmp(argv[i], "-oformat") == 0) {
	j = sscanf(argv[i+1], "%d", &oformat);
	if(j != 1) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -oformat option");
	  return 0;
	}
        i++; */
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-boxlw") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &boxlinewidth, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-padist") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &show_pa_dist, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-elldist") == 0) {
	elldist = 1;
	extendedPol = 1;   // Need the ellipticities to get the ellipticity distribution
      }else if(strcmp(argv[i], "-padist_norm") == 0) {
	pa_dist_normalise = 1;
      }else if(strcmp(argv[i], "-padist_stat") == 0) {
	pa_dist_stat = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %f %d", &fit_pa_dist_maxNrComp, &fit_pa_dist_sigmalimit, &fit_pa_dist_precision, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-padist_mask") == 0 || strcmp(argv[i], "-padistmask") == 0 || strcmp(argv[i], "-pamask") == 0 || strcmp(argv[i], "-pa_mask") == 0) {
	int ret;
	ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%1000s %f", pamask_name, &pamask_value, NULL);
	if(ret == 1) {
	  pamask_value = sqrt(-1); // Set pamask_value to NaN to indicate that this value was not specified on the command line	
	}else if(ret != 2) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Error parsing '%s' option. Expected 2 or 1 parameters to be provided.", argv[i]);
	  return 0;
	}
	
	/*
	j = sscanf(argv[i+1], "%s %f", pamask_name, &pamask_value);
	if(j != 2) {
	  j = sscanf(argv[i+1], "%s", pamask_name);
	  pamask_value = sqrt(-1); // Set pamask_value to NaN to indicate that this value was not specified on the command line
	  if(j != 1) {
	    fflush(stdout);
	    printerror(application.verbose_state.debug, "ERROR ppol: Error parsing -padist_mask option. Expected 2 or 1 parameters to be provided.");
	    return 0;
	  }
	  }*/
        i++;
      }else if(strcmp(argv[i], "-offpulsedata") == 0) {
	int ret;
	ret = parse_command_string(application.verbose_state, argc, argv, i+1, 0, 1, "%1000s %f", offpulsedata_name, &offpulsedata_rebinfactor, NULL);
	if(ret != 2) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Error parsing '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-lw") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &datalinewidth, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-xrange") == 0 || strcmp(argv[i], "-xrange_phase") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &plotl1, &plotl2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	if(strcmp(argv[i], "-xrange_phase") == 0) {
	  plotl1 *= 360.0;
	  plotl2 *= 360.0;
	}
        i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-xticks") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %d", &xtick, &nxsub, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-yrange") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &plotp1, &plotp2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-yrange2") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &plotI1, &plotI2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-paswing") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f %f", &overlayalpha, &overlaybeta, &overlaypa0, &overlayl0, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	overlayPA = 1;
        i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-subtract") == 0) {
	subtract = 1;
      }else if(strcmp(argv[i], "-subtractfile") == 0) {
	subtract = 1;
	subtractFile = i+1;
	i++;
      }else if(strcmp(argv[i], "-addfile") == 0) {
	subtract = 2;
	subtractFile = i+1;
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-1") == 0) {
	plotp1 = 0.0;
	plotp2 = 180.0;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-text") == 0) {   /* Will be handled by pgplotPAplot */
	i += 2;
      }else if(strcmp(argv[i], "-herrorbar") == 0) { /* Will be handled by pgplotPAplot */
	i += 1;
      }else if(strcmp(argv[i], "-herrorbar2") == 0) { /* Will be handled by pgplotPAplot */
	i += 1;
      }else if(strcmp(argv[i], "-verrorbar") == 0) { /* Will be handled by pgplotPAplot */
	i += 1;
      }else if(strcmp(argv[i], "-verrorbar2") == 0) { /* Will be handled by pgplotPAplot */
	i += 1;
      }else if(strcmp(argv[i], "-outline") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &outline, &outline_color, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-selectonpulse") == 0) {
	manualOnpulseSelection = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-noselectonpulse") == 0) {
	manualOnpulseSelection = -1;
      }else if(strcmp(argv[i], "-opm") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &jump_longitude[nrJumps], &jump_offset[nrJumps], NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	nrJumps++;
	if(nrJumps == MaxNrJumps) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR ppol: Maximum number of %s options exceeded.", argv[i]);
	  return 0;
	}
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-textkeywords") == 0) {
	str_list_replace_keys(0);
	terminateApplication(&application);
	return 0;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-2") == 0) {
	twoprofiles = 1;
      }else if(strcmp(argv[i], "-medianLdebias") == 0) {
	correctLbias = 0;
      }else if(strcmp(argv[i], "-noLdebias") == 0) {
	correctLbias = -1;
      }else {
	/* If the option is not recognized, assume it is a filename */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "ppol: Unknown option: %s", argv[i]);
	  printerror(application.verbose_state.debug, "\nRun ppol without command line arguments to show help");
	  terminateApplication(&application);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }
 

//START REGION DEVELOP
  if(show_pa_dist == 0 && pamask_name[0] != 0) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR ppol: When using the -padist_mask option, the -padist option should be specified as well.");
    return 0;
  }

  if(subtract && overlayPA == 0 && subtractFile == 0) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR ppol: Cannot use the -subtract option without the -paswing option");
    return 0;
  }
//START REGION RELEASE

  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR ppol: No files specified");
    return 0;
  }

//START REGION DEVELOP
  if(application.ppgplot_name != NULL) {
    if(pgopenoutputfile(application.ppgplot_name) == 0) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR ppol: Cannot open %s", application.ppgplot_name);
      return 0;
    }
  }
//START REGION RELEASE

  firstfiletoopen = 1;
  while((filename_ptr = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {

    if(firstfiletoopen == 0) {
      if(nokeypresses == 0) {
	i = pgplot_device_type(PlotDevice2, application.verbose_state);
	if(i < 3 || i > 10) {                      // Do not ask for keypresses if creating final plots in a file (such as postscript)
	  printf("Press a key in the terminal to continue\n");
	  fflush(stdout);
	  pgetch();
	}
      }
    }
    
    //    cleanPSRData(&datain, application.verbose_state);
    cleanPSRData(&dataout, application.verbose_state);
    if(!openPSRData(&datain, filename_ptr, application.iformat, 0, 1, 0, application.verbose_state))
      return 0;
    if(application.verbose_state.verbose) {
      fflush(stdout);
      printwarning(application.verbose_state.debug, "Input data contains %ld bins, %ld pulses, %ld polarizations and %ld frequencies.", (datain.NrBins), datain.NrSubints, (datain.NrPols), datain.NrFreqChan);
    }
    
    /* Search commandline for header parameters to overwrite header
       parameters. */
    if(PSRDataHeader_parse_commandline(&datain, argc, argv, application.verbose_state) == 0)
      return 0;

    if(datain.poltype == POLTYPE_ILVPAdPA || datain.poltype == POLTYPE_PAdPA || datain.poltype == POLTYPE_ILVPAdPATEldEl) {
      if(doprojection == 0) {
	printerror(application.verbose_state.debug, "ERROR ppol: File already contains PA data. You can use ppolFig to plot this data.");
	return 0;
//START REGION DEVELOP
      }else {
	if(datain.poltype != POLTYPE_ILVPAdPATEldEl) {
	  printerror(application.verbose_state.debug, "ERROR ppol: The -projection option expects either IQUV data, or polarization information including a PA and ellipticity (i.e. generated with ppol -extendedpol).");
	  return 0;
	}
//START REGION RELEASE
      }
    }

    compute_PA = 1;
    if(datain.isFolded && datain.foldMode == FOLDMODE_FIXEDPERIOD) {
      if(datain.fixedPeriod <= 0.0) {
	printwarning(application.verbose_state.debug, "WARNING ppol: Period does not appear to be set, assuming it is 1 sec.");
	datain.fixedPeriod = 1.0;
      }
    }
    if(datain.isFolded && datain.tsampMode == TSAMPMODE_FIXEDTSAMP) {
      if(get_tsamp(datain, 0, application.verbose_state) <= 0.0) {
	printwarning(application.verbose_state.debug, "WARNING ppol: Assuming full period is stored.");
	double period;
	int ret;
	ret = get_period(datain, 0, &period, application.verbose_state);
	if(ret == 2) {
	  printerror(application.verbose_state.debug, "ERROR ppol (%s): Cannot obtain period", datain.filename);
	  return 0;
	}
	datain.fixedtsamp = period/(double)datain.NrBins;
      }
    }
    if(get_tsamp(datain, 0, application.verbose_state) < 0.0000001 || get_tsamp(datain, 0, application.verbose_state) >= 10) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR ppol: The sampling time does not appear to be set correctly in the header. Consider using the -header option.");
      return 0;
    }
    
        
    for(i = 1; i < argc; i++) {
      if(strcmp(argv[i], "-header") == 0) {
	fflush(stdout);
	printwarning(application.verbose_state.debug, "WARNING ppol: If using the -header option, be aware it applied BEFORE the preprocessing.");
	break;
      }
    }
    if(preprocessApplication(&application, &datain) == 0) {
      fflush(stdout);
      printwarning(application.verbose_state.debug, "WARNING ppol: Applying preprocess options failed.");
      return 0;
    }

    /* Check needs to be after preprocessing, as polarization parameters could been converted. */
    if(!(doprojection && datain.poltype == POLTYPE_ILVPAdPATEldEl)) {
      if(datain.poltype == POLTYPE_COHERENCY) {
	if(application.verbose_state.verbose) {
	  printf("Data is written as coherency parameters. Going to convert them into Stokes parameters first.\n");
	}
	if(preprocess_stokes(&datain, application.verbose_state) == 0) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR ppol: Conversion into Stokes parameters failed.");
	  return 0;
	}
      }else if(datain.poltype == POLTYPE_UNKNOWN) {
	fflush(stdout);
	printwarning(application.verbose_state.debug, "WARNING ppol: It is assumed that the data are Stokes parameters.");
      }else if(datain.poltype != POLTYPE_STOKES) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR ppol: The polarization type of the data cannot be handled by ppol.");
	return 0;
      }
      if(datain.NrPols != 4) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR ppol: Need Stokes IQUV data, but the number of polarization channels != 4!");
	return 0;
      }
    }
    
    if(datain.NrSubints > 1) {
      if(application.oformat == PPOL_format) {
	printf("Data contains multiple subints. Data will be written out in FITS format.\n");
	application.oformat = FITS_format;
      }
    }
    if(datain.NrFreqChan > 1) {
      if(application.oformat == PPOL_format) {
	printf("Data contains multiple frequency channels. Data will be written out in FITS format.\n");
	application.oformat = FITS_format;
      }
    }

//START REGION DEVELOP
    if(offpulsedata_name[0] != 0) {
      if(!openPSRData(&data_offpulse, offpulsedata_name, application.iformat, 0, 1, 0, application.verbose_state))
	return 0;
      if(application.verbose_state.verbose) {
	fflush(stdout);
	printwarning(application.verbose_state.debug, "Offpulse input data contains %ld bins, %ld pulses, %ld polarizations and %ld frequencies.", (data_offpulse.NrBins), data_offpulse.NrSubints, (data_offpulse.NrPols), data_offpulse.NrFreqChan);
      }


      if(data_offpulse.poltype == POLTYPE_COHERENCY) {
	if(application.verbose_state.verbose) {
	  printf("Offpulse data is written as coherency parameters. Going to convert them into Stokes parameters first.\n");
	}
	if(preprocess_stokes(&data_offpulse, application.verbose_state) == 0) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR ppol: Conversion into Stokes parameters failed.");
	  return 0;
	}
      }else if(data_offpulse.poltype == POLTYPE_UNKNOWN) {
	fflush(stdout);
	printwarning(application.verbose_state.debug, "WARNING ppol: It is assumed that the offpulse data are Stokes parameters.");
      }else if(data_offpulse.poltype != POLTYPE_STOKES) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR ppol: The polarization type of the offpulse data cannot be handled by ppol, it should be Stokes or coherency parameters.");
	return 0;
      }
      if(data_offpulse.NrPols != 4) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR ppol: Need Stokes IQUV data, but the number of polarization channels != 4 for the offpulse data!");
	return 0;
      }

      /* Convert the specified onpulse regions which were defined as
	 fractions on the commandline. This is now possible because we
	 know the number of bins. */
      region_frac_to_int(&(application.onpulse), data_offpulse.NrBins, 0);
    }else {
//START REGION RELEASE
      region_frac_to_int(&(application.onpulse), datain.NrBins, 0);
//START REGION DEVELOP
    }
//START REGION RELEASE


//START REGION DEVELOP
    if(doprojection) {
      extendedPol = 0;
      writeout = 0;
    }


    // Get onpulse region (or just show it)
//START REGION RELEASE
    long total_nr_offpulse_file_bins;
//START REGION DEVELOP
    if(doprojection == 0 || (doprojection == 1 && projection_binnr == -1)) {
//START REGION RELEASE
      strcpy(txt, "Select on-pulse region ");
      strcat(txt, datain.psrname);
//START REGION DEVELOP
      if(offpulsedata_name[0] == 0) {
//START REGION RELEASE
	total_nr_offpulse_file_bins = datain.NrBins;
//START REGION DEVELOP
      }else {
	total_nr_offpulse_file_bins = data_offpulse.NrBins;
      }
//START REGION RELEASE
      profileI = (float *)malloc(total_nr_offpulse_file_bins*sizeof(float));
      if(profileI == NULL) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR ppol: Memory allocation error.");
	return 0;
      }
      int ret;
//START REGION DEVELOP
      if(offpulsedata_name[0] == 0) {
//START REGION RELEASE
	ret = read_profilePSRData(datain, profileI, NULL, 0, application.verbose_state);
//START REGION DEVELOP
      }else {
	ret = read_profilePSRData(data_offpulse, profileI, NULL, 0, application.verbose_state);
      }
//START REGION RELEASE
      if(ret != 1) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR ppol: Cannot form pulse profile");
	return 0;
      }

      if(manualOnpulseSelection == 1 || (application.onpulse.nrRegions == 0 && manualOnpulseSelection != -1)) {
	pgplot_clear_options(&pgplot_options);
	strcpy(pgplot_options.viewport.plotDevice, PlotDevice);
	strcpy(pgplot_options.box.xlabel, "Bin");
	strcpy(pgplot_options.box.ylabel, "Intensity");
	strcpy(pgplot_options.box.title, txt);
	selectRegions(profileI, total_nr_offpulse_file_bins, &pgplot_options, 0, 0, 0, &(application.onpulse), application.verbose_state);
	if(firstfiletoopen == 0) {
	  ppgslct(deviceID);
	}
      }else {
	pgplot_clear_options(&pgplot_options);
	strcpy(pgplot_options.box.xlabel, "Bin");
	strcpy(pgplot_options.box.ylabel, "Intensity");
	strcpy(pgplot_options.box.title, txt);
	strcpy(pgplot_options.viewport.plotDevice, PlotDevice);
	pgplotGraph1(&pgplot_options, profileI, NULL, NULL, total_nr_offpulse_file_bins, 0, total_nr_offpulse_file_bins, 0, 0, total_nr_offpulse_file_bins, 0, 0, 0, 1, 0, 1, 0, 1, 1, &(application.onpulse), -1, application.verbose_state);
	if(firstfiletoopen == 0) {
	  ppgslct(deviceID);
	}
      }
      region_int_to_frac(&(application.onpulse), 1.0/(float)total_nr_offpulse_file_bins, 0);
      regionShowNextTimeUse(application.onpulse, "-onpulse", "-onpulsef", stdout);
      
//START REGION DEVELOP
    }
//START REGION RELEASE

//START REGION DEVELOP
    if(decompose_method == 0) {
      if(subtractFile) {
	//	  cleanPSRData(&subtract_fin, application.verbose_state);
	if(!openPSRData(&subtract_fin, argv[subtractFile], 0, 0, 1, 0, application.verbose_state))
	  return 0;
	if(subtract_fin.NrBins != datain.NrBins) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR ppol: The reference PA-swing has a different number of pulse phase bins compared to the data file.");
	  return 0;
	}
	if(subtract_fin.NrFreqChan != 1) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR ppol: The reference PA-swing should have a single frequency channel.");
	  return 0;
	}
      }
    }
//START REGION RELEASE


// Generate the PA data, or else do the projection stuff
//START REGION DEVELOP
    if(doprojection == 0) {
//START REGION RELEASE
      verbose_definition verbose2;
      copyVerboseState(application.verbose_state, &verbose2);
      if(application.verbose_state.verbose && datain.NrSubints == 1 && datain.NrFreqChan == 1)
	verbose2.verbose = 1;
      else
	verbose2.verbose = 0;
      
      int nolongitudes;
      if(datain.NrFreqChan > 1 || datain.NrSubints > 1) {
	nolongitudes = 1; // Since output will not be ascii, the pulse longitudes will be ignored and for some formats writing out will fail because tsamp is no longer defined.
      }else {
	nolongitudes = 0;
      }
//START REGION DEVELOP
      if(offpulsedata_name[0] == 0) {
//START REGION RELEASE
	if(make_paswing_fromIQUV(&datain, extendedPol, application.onpulse, normalize, correctLbias, correctQV, correctV, nolongitudes, loffset, PAoffset_data, NULL, 1.0, verbose2) == 0) {
	  return 0;
	}
//START REGION DEVELOP
      }else {
	if(make_paswing_fromIQUV(&datain, extendedPol, application.onpulse, normalize, correctLbias, correctQV, correctV, nolongitudes, loffset, PAoffset_data, &data_offpulse, offpulsedata_rebinfactor, verbose2) == 0) {
	  return 0;
	}
      }
      if(decompose_method) {
	decompose_polarization_modes(datain, &mode1, &mode2, decompose_method, verbose2);
	strcpy(ofilename, filename_ptr);
	strcat(ofilename, ".mode1");
	if(!openPSRData(&mode1, ofilename, FITS_format, 1, 0, 0, application.verbose_state))
	  return 0;
	if(!writeHeaderPSRData(&mode1, argc, argv, application.history_cmd_only, application.verbose_state))
	  return 0;
	if(writePSRData(&mode1, mode1.data, application.verbose_state) != 1)
	  return 0;
	closePSRData(&mode1, 0, application.verbose_state);
	strcpy(ofilename, filename_ptr);
	strcat(ofilename, ".mode2");
	if(!openPSRData(&mode2, ofilename, FITS_format, 1, 0, 0, application.verbose_state))
	  return 0;
	if(!writeHeaderPSRData(&mode2, argc, argv, application.history_cmd_only, application.verbose_state))
	  return 0;
	if(writePSRData(&mode2, mode2.data, application.verbose_state) != 1)
	  return 0;
	closePSRData(&mode2, 0, application.verbose_state);
      }
//START REGION RELEASE
//START REGION DEVELOP
    }else {    // if doprojection != 0
      mapProjection = malloc(projection_nrx*projection_nry*sizeof(float));
      if(mapProjection == NULL) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR ppol: Cannot allocate memory");
	return 0;
      }
      //      void makeHammerAitoffMap(float *map, int nrx, int nry, float background);
      //      makeHammerAitoffMap(mapHammerAitoff, 640, 320, -90);
      // Note that the sigma_limit is not used
      if(sigmaset) {
	if(datain.poltype == POLTYPE_ILVPAdPATEldEl) {
	  printwarning(application.verbose_state.debug, "WARNING ppol: A sigma limit is set by the user, but since the provided data has already a sigma limit applied when calculating the PA and ellipticity, it has no effect.");
	}else {
	  printwarning(application.verbose_state.debug, "WARNING ppol: A sigma limit is set by the user, but isn't used when producing a Poincare sphere from Stokes parameters. Consider generating PA's and ellipticities first with ppol -extendedpol after assigning a sigma limit. The generated output file can then be read in with a separate ppol -projection command to generated a Poincare sphere.");
	}
      }
      if(subtractFile) {
	if(subtract == 2) {
	  printerror(application.verbose_state.debug, "ERROR ppol: The -addfile option cannot be used together with the -projection option.");
	  return 0;
	}
	if(make_projection_map_fromIQUV(datain, mapProjection, projection_nrx, projection_nry, 0, projection_binnr, application.onpulse, projection_weighting, proj_type, proj_rot_long+2.0*(PAoffset_data), proj_rot_lat, &subtract_fin, application.verbose_state) == 0)
	  return 0;
      }else {
	if(make_projection_map_fromIQUV(datain, mapProjection, projection_nrx, projection_nry, 0, projection_binnr, application.onpulse, projection_weighting, proj_type, proj_rot_long+2.0*(PAoffset_data), proj_rot_lat, NULL, application.verbose_state) == 0)
	  return 0;
      }
    }

//START REGION RELEASE


//START REGION DEVELOP
    if(doprojection == 0 && decompose_method == 0) {
//START REGION RELEASE
      if(sigma_limit > 0
//START REGION DEVELOP
 || subtract
//START REGION RELEASE
	 ) {                                  // If not making a Poincare sphere projection, and a sigma limit and/or a PA-swing needs to be subtracted
//START REGION DEVELOP
//START REGION RELEASE
	// XXXXXXXXXXXXXXXXXXXXXXX   datain.NrPols-2 works NOT FOR ELLIPTICITY FORMAT, AS ASSUMES IT IS IN THE ONE BUT LAST POLARIZATION. THIS MADE SENSE, BECAUSE OT WOULD ALSO WORK FOR THE PAdPA FORMAT.
	int pachannel;
	if(datain.poltype == POLTYPE_ILVPAdPATEldEl) {
	  pachannel = 3;
	}else {
	  pachannel = datain.NrPols-2;
	}
	for(i = 0; i < datain.NrSubints; i++) {
	  for(f = 0; f < datain.NrFreqChan; f++) {
	    for(j = 0; j < datain.NrBins; j++) {
	      // Subtract something from PA's is required
//START REGION DEVELOP
	      if(subtract) {
		if(subtractFile == 0) {   // Subtract RVM model
		  datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] -= paswing(overlayalpha, overlaybeta, datain.tsamp_list[j], overlaypa0, overlayl0, nrJumps, jump_longitude, jump_offset, 0, 0);
		  datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] = derotate_180(datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))]);
		}else {
		  int pachannel_subtract_fin;
		  if(subtract_fin.poltype == POLTYPE_ILVPAdPATEldEl) {
		    pachannel_subtract_fin = 3;
		  }else {
		    pachannel_subtract_fin = subtract_fin.NrPols-2;
		  }
		  if(datain.data[j+datain.NrBins*(pachannel+1 + datain.NrPols*(f+datain.NrFreqChan*i))] > 0 && subtract_fin.data[j+subtract_fin.NrBins*(pachannel_subtract_fin+1)] > 0) {
		    if(subtract == 2) {  // Adding rather than subtracting
		      datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] += subtract_fin.data[j+subtract_fin.NrBins*(pachannel_subtract_fin)];
		    }else {
		      datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] -= subtract_fin.data[j+subtract_fin.NrBins*(pachannel_subtract_fin)];
		    }
		    datain.data[j+datain.NrBins*(pachannel+1 + datain.NrPols*(f+datain.NrFreqChan*i))] = sqrt(datain.data[j+datain.NrBins*(pachannel+1 + datain.NrPols*(f+datain.NrFreqChan*i))]*datain.data[j+datain.NrBins*(pachannel+1 + datain.NrPols*(f+datain.NrFreqChan*i))]+subtract_fin.data[j+subtract_fin.NrBins*(pachannel_subtract_fin+1)]*subtract_fin.data[j+subtract_fin.NrBins*(pachannel_subtract_fin+1)]);
		  }else {
		    datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] = 0;
		    datain.data[j+datain.NrBins*(pachannel+1 + datain.NrPols*(f+datain.NrFreqChan*i))] = -1;
		  }
		  // Make PA's from -90 to 90 when subtracting a model, so that residuals fluctuate around 0.
		  datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] = derotate_180(90+datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))]);
		  datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] -= 90;
		}
	      }
//START REGION RELEASE
	      // Apply L sigma limit to PA points
	      // Note: when showing file, offpulse_rms not defined for all formats
	      if(sigma_limit > 0) {
		if(datain.offpulse_rms == NULL) {
		  if(i == 0 && f == 0 && j == 0) {
		    fflush(stdout);
		    printwarning(application.verbose_state.debug, "WARNING ppol: Cannot apply sigma limit on PA-points for this type of input file. The used sigma limit will be the same as when the input file was generated.");
		  }
		}else {
		  if(sigmaI == 0) {
		    if(datain.data[j+datain.NrBins*(1+datain.NrPols*(f+i*datain.NrFreqChan))] < sigma_limit*datain.offpulse_rms[1+datain.NrPols*(f + datain.NrFreqChan*i)]) {
		      datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] = 0;
		      datain.data[j+datain.NrBins*(pachannel+1 + datain.NrPols*(f+datain.NrFreqChan*i))] = -1;
		    }
		  }else {
		    if(datain.data[j+datain.NrBins*(1+datain.NrPols*(f+i*datain.NrFreqChan))] < sigma_limit*datain.offpulse_rms[0+datain.NrPols*(f + datain.NrFreqChan*i)]) {
		      datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] = 0;
		      datain.data[j+datain.NrBins*(pachannel+1 + datain.NrPols*(f+datain.NrFreqChan*i))] = -1;
		    }
		  }
//START REGION DEVELOP
		  if(extendedPol) {
		    // The decision of ellipticity being significant depends on the rms in the total amount of polarization
		    if(sigmaI == 0) {
		      if(datain.data[j+datain.NrBins*(5+datain.NrPols*(f+i*datain.NrFreqChan))] < sigma_limit*datain.offpulse_rms[5+datain.NrPols*(f + datain.NrFreqChan*i)]) {
			datain.data[j+datain.NrBins*(6 + datain.NrPols*(f+datain.NrFreqChan*i))] = 0;
			datain.data[j+datain.NrBins*(7 + datain.NrPols*(f+datain.NrFreqChan*i))] = -1;
		      }
		    }else {
		      if(datain.data[j+datain.NrBins*(5+datain.NrPols*(f+i*datain.NrFreqChan))] < sigma_limit*datain.offpulse_rms[0+datain.NrPols*(f + datain.NrFreqChan*i)]) {
			datain.data[j+datain.NrBins*(6 + datain.NrPols*(f+datain.NrFreqChan*i))] = 0;
			datain.data[j+datain.NrBins*(7 + datain.NrPols*(f+datain.NrFreqChan*i))] = -1;
		      }
		    }
		  }
//START REGION RELEASE
		}
	      }
//START REGION DEVELOP
	      if(doselectmode) {
		if(extendedPol == 0) {
		  fflush(stdout);
		  printerror(application.verbose_state.debug, "ERROR ppol: -extendedpol should be used when the -selectmode option is used.");
		  return 0;
		}
		// Angle between -90 and 90 deg
		double paoffset = derotate_180_small_double(datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] - selectmode_cpa);
		// Derotates an angle between -45 and 45 degrees
		double elloffset = derotate_90_small_double(datain.data[j+datain.NrBins*(6 + datain.NrPols*(f+datain.NrFreqChan*i))] - selectmode_cell);	
		if(fabs(paoffset) > selectmode_dpa || fabs(elloffset) > selectmode_dell) {
		  datain.data[j+datain.NrBins*(pachannel + datain.NrPols*(f+datain.NrFreqChan*i))] = 0;
		  datain.data[j+datain.NrBins*(pachannel+1 + datain.NrPols*(f+datain.NrFreqChan*i))] = -1;
		  datain.data[j+datain.NrBins*(6 + datain.NrPols*(f+datain.NrFreqChan*i))] = 0;
		  datain.data[j+datain.NrBins*(7 + datain.NrPols*(f+datain.NrFreqChan*i))] = -1;
		}
	      }
//START REGION RELEASE
	    }
	  }
	}
      } // End of modifying the PA-data because of a sigma limit of PA-swing subtraction

      /*
	for(i = 0; i < datain.NrSubints; i++) {
	for(j = 0; j < datain.NrBins; j++) {
	printf("XXXXX %f %f\n", datain.data_pa[j+i*datain.NrBins], datain.data_dpa[j+i*datain.NrBins]);
	}
	}
      */
      
      if(datain.NrSubints > 1 || datain.NrFreqChan > 1) {
	fflush(stdout);
	printwarning(application.verbose_state.debug, "WARNING ppol: Only showing polarization of the first subint/frequency channel");
      }
      pgplot_clear_options(&pgplot_options);
      pgplot_options.viewport.dxplot = 0;
      pgplot_options.viewport.dyplot = 0;
      pgplot_options.viewport.xsize = xsize/0.8;  // Divide by 0.8 to make ppol backwards compatable.
      pgplot_options.viewport.ysize = ysize;

      ppgpage();
      strcpy(pgplot_options.viewport.plotDevice, PlotDevice2);
      pgplot_options.viewport.noclear = 1;  // Done somewhere else in the code
      pgplot_options.viewport.dontclose = 1;
      pgplot_options.viewport.dontopen = !firstfiletoopen;
      pgplot_options.box.box_lw = boxlinewidth;
      pgplot_options.box.label_lw = boxlinewidth;
      pgplot_options.box.box_xtick = xtick;
      pgplot_options.box.box_nxsub = nxsub;
      if(titleset)
	strcpy(pgplot_options.box.title, title);
      else
	strcpy(pgplot_options.box.title, datain.filename);
      pgplot_options.box.title_ch = titlech;
      pgplot_options.box.title_lw = titlelw;

      pgplotPAplot(datain, extendedPol, 0, 0
//START REGION DEVELOP
		   + extendedPol
//START REGION RELEASE
		   , &pgplot_options, "Pulse longitude [deg]", "I,Linear,V", "PA [deg]", "\\gx [deg]", plotl1, plotl2, 0, plotI1, plotI2, plotp1, plotp2, 0.0, sigma_limit, datalinewidth, ysizepa, dashed, noylabel, "-text", "-herrorbar", "-herrorbar2", "-verrorbar", "-verrorbar2", argc, argv, outline, outline, outline_color, overlayPA, overlayalpha, overlaybeta, overlaypa0, overlayl0, overlayFine, nrJumps, jump_longitude, jump_offset, NULL, -90, 90, 1.0, NULL, 1.0, 0, application.verbose_state);
      if(firstfiletoopen) {
	ppgqid(&deviceID);
      }
//START REGION DEVELOP
    }else if(decompose_method) { // No plot
    }else {  //End of if(doprojection == 0)                 // Make the actual Poincare plot
      /*
      strcpy(pgplot_options.box.xlabel, "Bin");
      strcpy(pgplot_options.box.ylabel, "Intensity");
      strcpy(pgplot_options.box.title, txt);
      strcpy(pgplot_options.viewport.plotDevice, PlotDevice);
      */
      pgplot_clear_options(&pgplot_options);
      strcpy(pgplot_options.viewport.plotDevice, PlotDevice2);
      pgplot_options.viewport.dontclose = 1;
      pgplot_options.viewport.dontopen = 1;
      ppgopen(pgplot_options.viewport.plotDevice);
      int cmaptype;
      cmaptype = pgplot_device_type(NULL, application.verbose_state);
      //      printf("cmaptype=%d\n", cmaptype);
      if(cmaptype <= 2)
	cmaptype = PPGPLOT_HEAT;
      else
	cmaptype = PPGPLOT_INVERTED_HEAT4;
      if(projection_nolabels) {
	pgplot_options.box.drawbox = 0;
	//	sprintf(pgplot_options.box.box_xopt, "bcsti");
	//	sprintf(pgplot_options.box.box_yopt, "bcsti");
      }else {
	pgplot_options.box.drawbox = 1;
	//	sprintf(pgplot_options.box.box_xopt, "bcnsti");
	//	sprintf(pgplot_options.box.box_yopt, "bcnsti");
      }
      if(titleset) {
	strcpy(pgplot_options.box.title, title);
      }else {
	strcpy(pgplot_options.box.title, datain.filename);
      }
      pgplot_options.box.title_ch = titlech;
      pgplot_options.box.title_lw = titlelw;
      pgplot_options.box.box_lw = boxlinewidth;
      pgplot_options.box.label_lw = boxlinewidth;
      pgplot_options.box.box_xtick = xtick;
      pgplot_options.box.box_nxsub = nxsub;
      if(proj_type == 3) {
	if(pgplotMap(&pgplot_options, mapProjection, projection_nrx, projection_nry, -180, 180, -180, 180, -90, 90, -90, 90, cmaptype, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, application.verbose_state) != 1) {
	  return 0;
	}
      }else {
	if(pgplotMap(&pgplot_options, mapProjection, projection_nrx, projection_nry, -2.25, 2.25, -2.25, 2.25, -1.125, 1.125, -1.125, 1.125, cmaptype, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, application.verbose_state) != 1) {
	  return 0;
	}
      }
      drawSphericalGrid(30, 30,  proj_rot_long, proj_rot_lat, 3, proj_type);
    }
//START REGION RELEASE
   
//START REGION DEVELOP
    if(show_pa_dist && doprojection == 0) {
      if(nokeypresses == 0) {
	printf("Press a in the terminal key to continue\n");
	pgetch();
      }
      pgplot_clear_options(&pgplot_options);
      pgplot_options.viewport.dontopen = 1;
      pgplot_options.viewport.noclear = 0;
      strcpy(pgplot_options.box.xlabel, "Pulse longitude bin");
      if(elldist == 0)
	strcpy(pgplot_options.box.ylabel, "PA [deg]");
      else
	strcpy(pgplot_options.box.ylabel, "\\gx [deg]");
      if(titleset)
	strcpy(pgplot_options.box.title, title);
      else
	strcpy(pgplot_options.box.title, datain.filename);

      datafile_definition *pamask_ptr;
      pamask_ptr = NULL;
      if(pamask_name[0] != 0) { // Check if mask is provided
	//	cleanPSRData(&pamask, application.verbose_state);
	if(!openPSRData(&pamask, pamask_name, 0, 0, 1, 0, application.verbose_state))
	  return 0;
	pamask_ptr = &pamask;
      }

      if(make_pa_distribution(datain, &padist_data, show_pa_dist, pa_dist_normalise, pamask_ptr, pamask_value, elldist, application.verbose_state) == 0)
	return 0;
      if(pamask_name[0] != 0) { // Check if mask is provided
	closePSRData(&pamask, 0, application.verbose_state);
      }
      int ret;
      if(elldist) {
	ret = pgplotMap(&pgplot_options, padist_data.data, padist_data.NrBins, padist_data.NrSubints, 0.5, padist_data.NrBins-0.5, padist_data.NrBins*plotl1/360.0, padist_data.NrBins*plotl2/360.0, -45+0.5*(90.0/(float)show_pa_dist), 45-0.5*(90.0/(float)show_pa_dist), -45, 45, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, application.verbose_state);
      }else {
	ret = pgplotMap(&pgplot_options, padist_data.data, padist_data.NrBins, padist_data.NrSubints, 0.5, padist_data.NrBins-0.5, padist_data.NrBins*plotl1/360.0, padist_data.NrBins*plotl2/360.0, -90+0.5*(180.0/(float)show_pa_dist), 90-0.5*(180.0/(float)show_pa_dist), -90, 90, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, application.verbose_state);
      }
      if(ret == 0)
	return 0;
      if(pa_dist_stat) {
	if(elldist) {
	  printerror(application.verbose_state.debug, "ERROR ppol: The -padist_stat option is not implemented for the ellipticity angle distribution.");
	  return 0;	  
	}
	if(fit_pa_distribution(padist_data, &padist_fit, &padist_resid, fit_pa_dist_maxNrComp, fit_pa_dist_sigmalimit, fit_pa_dist_precision, application.verbose_state) == 0) {
	  return 0;
	}
	if(nokeypresses == 0) {
	  printf("Press a in the terminal key to continue\n");
	  pgetch();
	}
	if(pgplotMap(&pgplot_options, padist_fit.data, padist_data.NrBins, padist_data.NrSubints, 0.5, padist_data.NrBins-0.5, 0, padist_data.NrBins, -90+0.5*(180.0/(float)show_pa_dist), 90-0.5*(180.0/(float)show_pa_dist), -90, 90, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, application.verbose_state) == 0)
	  return 0;
	if(nokeypresses == 0) {
	  printf("Press a in the terminal key to continue\n");
	  pgetch();
	}
	if(pgplotMap(&pgplot_options, padist_resid.data, padist_data.NrBins, padist_data.NrSubints, 0.5, padist_data.NrBins-0.5, 0, padist_data.NrBins, -90+0.5*(180.0/(float)show_pa_dist), 90-0.5*(180.0/(float)show_pa_dist), -90, 90, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, application.verbose_state) == 0)
	  return 0;
      }
      if(writeout && pamask_name[0] == 0) { // Do not write out the pa distribution when masking is done, as otherwise the masked data will be written to the same file, as multiple extensions cannot be provided on the command line
	// First write out the derived pa distribution from the data
	if(extension != NULL) {
	  if(change_filename_extension(filename_ptr, ofilename, extension, MaxFilenameLength, application.verbose_state) == 0) {
	    fflush(stdout);
	    printerror(application.verbose_state.debug, "ERROR ppol: Changing filename failed");
	    return 0;
	  }
	}
	if(writeoutFilename != -1) {
	  strcpy(ofilename, argv[writeoutFilename]);
	}
	if(writeoutFilename != -1 || extension != NULL) {
	  int oformat;
	  oformat = application.oformat;
	  if(oformat == PPOL_format) {
	    oformat = FITS_format;
	  }
	  if(!openPSRData(&padist_data, ofilename, oformat, 1, 0, 0, application.verbose_state))
	    return 0;
	  if(!writeHeaderPSRData(&padist_data, argc, argv, application.history_cmd_only, application.verbose_state))
	    return 0;
	  if(writePSRData(&padist_data, padist_data.data, application.verbose_state) != 1)
	    return 0;
	}

	// Write out the fit to the pa distribution
	if(pa_dist_stat) {
	  if(extension != NULL) {
	    if(change_filename_extension(filename_ptr, ofilename, extension, MaxFilenameLength, application.verbose_state) == 0) {
	      fflush(stdout);
	      printerror(application.verbose_state.debug, "ERROR ppol: Changing filename failed");
	      return 0;
	    }
	  }
	  if(writeoutFilename != -1) {
	    strcpy(ofilename, argv[writeoutFilename]);
	  }
	  strcat(ofilename, ".fit");
	  if(writeoutFilename != -1 || extension != NULL) {
	    if(!openPSRData(&padist_fit, ofilename, FITS_format, 1, 0, 0, application.verbose_state))
	      return 0;
	    if(!writeHeaderPSRData(&padist_fit, argc, argv, application.history_cmd_only, application.verbose_state))
	      return 0;
	    if(writePSRData(&padist_fit, padist_fit.data, application.verbose_state) != 1)
	      return 0;
	  }

	  if(extension != NULL) {
	    if(change_filename_extension(filename_ptr, ofilename, extension, MaxFilenameLength, application.verbose_state) == 0) {
	      fflush(stdout);
	      printerror(application.verbose_state.debug, "ERROR ppol: Changing filename failed");
	      return 0;
	    }
	  }
	  if(writeoutFilename != -1) {
	    strcpy(ofilename, argv[writeoutFilename]);
	  }
	  strcat(ofilename, ".resid");
	  if(writeoutFilename != -1 || extension != NULL) {
	    if(!openPSRData(&padist_resid, ofilename, FITS_format, 1, 0, 0, application.verbose_state))
	      return 0;
	    if(!writeHeaderPSRData(&padist_resid, argc, argv, application.history_cmd_only, application.verbose_state))
	      return 0;
	    if(writePSRData(&padist_resid, padist_resid.data, application.verbose_state) != 1)
	      return 0;
	  }
	}


      }
      closePSRData(&padist_data, 0, application.verbose_state);
      if(pa_dist_stat) {
	closePSRData(&padist_fit, 0, application.verbose_state);
	closePSRData(&padist_resid, 0, application.verbose_state);
      }
    }
//START REGION RELEASE

//START REGION DEVELOP
	if(subtractFile) {
	  //	  printf("XXXX close\n");
	  closePSRData(&subtract_fin, 0, application.verbose_state);
	}
//START REGION RELEASE

    /* Copy header parameters to output header */ 
    if(writeout
//START REGION DEVELOP
       && (show_pa_dist == 0 || pamask_name[0] != 0) && doprojection == 0    // Pa-distribution is written out above. This writes out the standard ppol output
//START REGION RELEASE
) {
      copy_params_PSRData(datain, &dataout, application.verbose_state);
      if(extension != NULL) {
	if(change_filename_extension(filename_ptr, ofilename, extension, MaxFilenameLength, application.verbose_state) == 0) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR ppol: Changing filename failed");
	  return 0;
	}
      }
      if(writeoutFilename != -1) {
	strcpy(ofilename, argv[writeoutFilename]);
      }
      if(writeoutFilename != -1 || extension != NULL) {
	if(!openPSRData(&dataout, ofilename, application.oformat, 1, 0, 0, application.verbose_state))
	  return 0;
      }else {
	if(datain.NrSubints > 1 || datain.NrFreqChan > 1) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot write data to stdout when there are multiple subints and/or frequency channels");
	  return 0;
	}
	dataout.fptr = stdout;
	dataout.fptr_hdr = stdout;
	dataout.format = application.oformat;
      }
//START REGION DEVELOP
      if(pamask_name[0] != 0 && !isnormal(pamask_value))
	dataout.poltype = POLTYPE_UNKNOWN;  // When applying the pa-mask and writing out the classification (i.e. -1, 0 and +1 for example), then set the poltype to be unknown, so that some of the pre-process options will work on this type of data.
//START REGION RELEASE
      if(!writeHeaderPSRData(&dataout, argc, argv, application.history_cmd_only, application.verbose_state))
	return 0;
      if(application.oformat == PPOL_format) {
	extended = 1;
	onlysignificantPA = 0;
      }else {
	extended = 0;
	onlysignificantPA = 1;
      }
      if(dataout.format == PPOL_format || dataout.format == PPOL_SHORT_format) {
	dataout.offpulse_rms = datain.offpulse_rms;
	if(writePPOLfile(dataout, datain.data, extended, onlysignificantPA, twoprofiles, 0.0, application.verbose_state) == 0) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR ppol: Cannot write data");
	  return 0;
	}
      }else {
	for(j = 0; j < datain.NrPols; j++) {
	  for(i = 0; i < datain.NrSubints; i++) {
	    for(f = 0; f < datain.NrFreqChan; f++) {
	      if(writePulsePSRData(&dataout, i, j, f, 0, datain.NrBins, &datain.data[datain.NrBins*(j+datain.NrPols*(f+i*datain.NrFreqChan))], application.verbose_state) != 1) {
		fflush(stdout);
		printerror(application.verbose_state.debug, "ERROR ppol: Cannot write data");
		return 0;
	      }
	    }
	  }
	}
      }
    }
    
    closePSRData(&datain, 0, application.verbose_state);
    if(dataout.format == PPOL_format || dataout.format == PPOL_SHORT_format) {
      dataout.offpulse_rms = NULL;  // Was set to be the same as datain, so avoid double free
    }
    closePSRData(&dataout, 0, application.verbose_state);
    firstfiletoopen = 0;
    if(compute_PA)
      free(profileI);
//START REGION DEVELOP
    if(doprojection)
      free(mapProjection);
    if(offpulsedata_name[0] != 0) {
      closePSRData(&data_offpulse, 0, application.verbose_state);
    }
//START REGION RELEASE
  }  // End of loop over input files
  ppgend();
//START REGION DEVELOP
  if(application.ppgplot_name != NULL) {
    pgcloseoutputfile();
  }
//START REGION RELEASE
  terminateApplication(&application);
  return 0;
}


//START REGION DEVELOP
