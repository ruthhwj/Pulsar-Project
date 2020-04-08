//START REGION RELEASE
/* 
Things to do:
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
//#ifndef __APPLE__
//#include <malloc.h>
//#endif
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "psrsalsa.h"
#include "cpgplot.h"

#define MaxNrJumps     100


void PrintHelp();
void convertAlphaBeta(double *alpha, double *beta
//START REGION DEVELOP
		      , int useslope, int usezeta
//START REGION RELEASE
		      );
void calcBeamWidths(int nalpha, int nbeta, double alphastart, double alphaend, double betastart, double betaend, 
//START REGION DEVELOP
		    int usezeta, int useslope, double elliptical_beam, 
//START REGION RELEASE
		    int calculate_beam_widths, int calculate_interpulse_widths, double pulse_width2, float *rhogrid, float *rhogrid2, int nocounters);

int internal_fit_pa_or_l0;

void PlotGrid(float *chigrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, double level, double suppress_fac, int GridDeviceID, double alpha, double beta, double lwbox, double labelcharheight, double boxlabelcharheight, int drawCross, int draw_title, int drawcontours, int nogray, 
//START REGION DEVELOP
	      int usezeta, int useslope, int use_sigma, 
//START REGION RELEASE
	      double chimax, double chimin, int maptype, int showwedge, char *showwedge_label);

void DoFitting(double alpha0, double beta0, double pa0, double dpa0, double l0, double dl0, double dh0, double ddh0, double ftol, double *fit_pa0, double *fit_l0, double *fit_a, double *fit_b, double *fit_dh, double *chi, int *nfunk, int searchAll, int report, FILE *reportStream, int finderrors, double nrofsigmas, int *nfitparameters, int amoeba_algorithm, verbose_definition verbose);

void PlotPAswing(double alpha, double beta, double pa0, double l0, int PlotFit, double leftPulseLongitude, double rightPulseLongitude, double dh);

void PlotContours(float *rhogrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, int nrlevels, float *TR, int GridDeviceID, int contour_txt, int contourcolor, int fixedContours, int nruserContours, float *userContours, int lwbox, int dotted);

void calcIntersectionRhoAndBanana(float *rhogrid, float *chigrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, int nruserContours, float *userContours, double chimax, double chimin, verbose_definition verbose);

void print_steepness(double alpha, double beta, double l0, double pa0, int verbose, double dh, double *sina_b);

double funk(double x[]);

double internal_funk_gsl(double *x, void *params);

double dy_180(double y1, double y2);
double dy_90(double y1, double y2);


//START REGION DEVELOP
void SHOWREVISIONINFO_prog() {
#include "ppolFit.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);
//START REGION RELEASE

/* Not all used anymore, especially fitalpha etc. ??? */
struct {
  /* The data */
  double *data_l;
  float *data_pa, *data_dpa;
  int NrDataPoints;
  /* Possible OPM jumps to take into account */
  double jump_longitude[MaxNrJumps], jump_offset[MaxNrJumps];
  int nrJumps, autojump;
  /* Possible longitude of a height difference */
  double add_height_longitude;
  /* The start position of the steepest gradient of the PA-swing and
     the maximum allowed deviation. */
  double l0_start, max_l0_diff;
  /* Possible point you want to force the PA-swing through */
  int force_set;
  double force_l, force_pa, force_dpa;
//START REGION DEVELOP
  int math_setting;
  int do_rho_width;
//START REGION RELEASE
  double pulse_width, sigma_width, rho_bcw, sigma_rho;
}fitterinfo;


int main(int argc, char **argv)
{
  char dumpfile[1000], c, device1[100], device2[100], txt[MaxStringLength], *txtptr, *showwedge_label;
  char prefix[1000];
  int i, j, nfunk, loadresults, GridDeviceID, PADeviceID, PSDeviceID, macrofilename, ret;
  int calculate_beam_widths, nrcontourlevels, nrcontourlevels2, redraw;
  int fixedContours, printsmooth, nogray, drawcontours, boxlw, nrfitparams;
  int contour_plot, contour_txt, calculate_interpulse_widths, iformat;
  int alphaset, betaset, lset, paset, printnofit, printfit, GridSearch, nalpha, nbeta;
  int onlyshowbest, showgraphics, drawCross, title_txt;
  int enableerrors, nruserContours, doBruteForce; //, nocounters;
  int showwedge, version;
  int gridDevice_resx, gridDevice_resy, paDevice_resx, paDevice_resy;
  int invertGrayscale, amoeba_algorithm, doMC, devicenores, fixseed;
  int suppress_greyscale, beamwidth_params_only_w;
  float *chigrid, *l0grid, *pa0grid, *dhgrid, *rhogrid, *rhogrid2, TR[6], userContours[500];
  double ftol;
  double l0, dl0, pa0, dpa0, dh0, ddh0, alpha0, beta0, chi, chi_old, chimax, chimin;
  double fit_pa0, fit_l0, fit_dh0, fit_alpha, fit_beta, bestalpha, bestbeta;
  double labelcharheight, boxlabelcharheight, suppress_fac;
  double alphastart, alphaend, betastart, betaend, add_longitude_shift;
  double leftPulseLongitude, rightPulseLongitude, nrofsigmas;
  double level, pulse_width2, bestl0, bestpa0; //, beshdh;
  double paErrorFac;
  FILE *fin, *macrofile, *beamwidth_params_fin;
  double sina_b;
//START REGION DEVELOP
  int calc_prob = 0;
  int use_sigma_level;
  double *prob_a, *prob_b;
//START REGION RELEASE
  int contourcolor, contourcolorIP;
  int doContourRange, contcol_nr_edge_steps, contcol_nr_fid_plane_steps;
  double contcol_edge1, contcol_edge1plus, contcol_edge1min, contcol_edge2, contcol_edge2plus, contcol_edge2min, contcol_phi0, contcol_phi0plus, contcol_phi0min, contcol_l0, contcol_l0plus, contcol_l0min, contcol_p;
  datafile_definition datain;
  int usegsl_error_grid_search;
  psrsalsaApplication application;
  verbose_definition beamwidth_params_vebose;
//START REGION DEVELOP
  int usezeta, useslope, doQuarterRange, dumpgrids;
  double elliptical_beam;
//START REGION RELEASE

  initApplication(&application, "ppolFit", "[options] paswing_file");
//START REGION DEVELOP
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
//START REGION RELEASE

  application.switch_verbose = 1;
  application.switch_nocounters = 1;
  application.switch_debug = 1;
  application.switch_history_cmd_only = 1;


//START REGION DEVELOP
  use_sigma_level = 0;
  fitterinfo.math_setting = 0;
  elliptical_beam = sqrt(-1);
//START REGION RELEASE
  suppress_greyscale = 6;
  fitterinfo.pulse_width = 0;
  nogray = 0;
  ftol = 1e-3;
  alphaset = betaset = lset = paset = printnofit = printfit = GridSearch = loadresults = 0;
  alphastart = 0;
  alphaend = 180;
  betastart = -90;
  betaend = 90;
  fitterinfo.max_l0_diff = 45;
  calculate_beam_widths = 0;
  calculate_interpulse_widths = 0;
  strcpy(dumpfile, "dump.dat");
  contour_txt = 1;
  nrcontourlevels = 19;
  nrcontourlevels2 = 10;
  fixedContours = 0;
  printsmooth = 0;
  onlyshowbest = 0;
  showgraphics = 1;
  iformat = -1;
  fitterinfo.nrJumps = 0;
  fitterinfo.autojump = 0;
  add_longitude_shift = 0;
  fitterinfo.add_height_longitude = 1e9;
  drawCross = 1;
  fitterinfo.force_set = 0;
  title_txt = 1;
  boxlw = 4;
  drawcontours = 3;
  labelcharheight = 1.7;
  boxlabelcharheight = 1;
  suppress_fac = 1;
  dh0 = -1e9;
  ddh0 = -1e9;
  enableerrors = 0;
  contour_plot = 0;
  doBruteForce = 0;
  sprintf(device1, "?");
  sprintf(device2, "?");
  nruserContours = 0;
  showwedge_label = "Reduced-\\gx\\u2";
  macrofilename = 0;
  beamwidth_params_fin = NULL;
  showwedge = 0;
  nrofsigmas = 3;
  contourcolor = 2;
  contourcolorIP = 3;
  gridDevice_resx = 680;
  gridDevice_resy = 680;
  paDevice_resx = 680;
  paDevice_resy = 340;
  amoeba_algorithm = 0;
  //  nocounters = 0;
  usegsl_error_grid_search = 1;
  devicenores = 0;
  fixseed = 0;
  cleanVerboseState(&beamwidth_params_vebose);
  /* To get rid of warnings, so assume I know what I'm doing. */
  fin = NULL;
  rhogrid =  rhogrid2 = NULL;
  bestalpha = bestbeta = -1e10;
  chimin = chimax = -1e10;
  PADeviceID = GridDeviceID = 0;
  doMC = 0;
  paErrorFac = 1;
  sprintf(prefix, "chi2");
  doContourRange = 0;
//START REGION DEVELOP
  usezeta = 0;
  useslope = 0;
  doQuarterRange = 0;
  dumpgrids = 0;
//START REGION RELEASE

  if(argc < 2) {
    printf("Program to fit the rotating vector model to a position angle swing as generated by ppol. Usage:\n\n");
    printApplicationHelp(&application);
    fprintf(stdout, "Define a grid search:\n");
    fprintf(stdout, "  -g \"a b\"      Do gridsearch of a alpha and b beta beta points\n"); 
    fprintf(stdout, "  -A \"min max\"  Specify search range in alpha [def=\"0 180\" degrees]\n");
    fprintf(stdout, "  -B \"min max\"  Specify search range in beta  [def=\"-90 90\" degrees]\n");
    fprintf(stdout, "  -l \"l0 dl\"    Set start value and stepsize for the longitude of magnetic axis\n");
    fprintf(stdout, "                [def=\"180 90\" degrees], same as -l0.\n");
    fprintf(stdout, "  -pa \"pa0 dpa\" Set start value and stepsize for the PA at the position of the\n");
    fprintf(stdout, "                magnetic axis [def=\"0 90\" degrees], same as -pa0.\n\n");
    fprintf(stdout, "Alternative: fixed alpha and beta fitting (only fit for pa_0 and l_0):\n");
    fprintf(stdout, "  -a            Fixed alpha value in deg (use together with -b)\n");
    fprintf(stdout, "  -b            Fixed beta value in deg (use together with -a)\n");
    fprintf(stdout, "  -printfit     Print out fit of pa-swing to the terminal\n"); 
    fprintf(stdout, "  -printnofit   Like -printfit, but pa_0 and l_0 are taken from the -l0 and -pa0\n");
    fprintf(stdout, "                options, so no fit will be done\n"); 
    fprintf(stdout, "  -printsmooth  Phase-wraps are tried to be avoided with -printfit\n"); 
    fprintf(stdout, "\nBeam-width (rho) contour options, which can be used in grid search mode.\n              (Consider using the -cont or 'r' options):\n");
    fprintf(stdout, "  -wmp        Define pulse width of main-pulse in degrees\n");
    fprintf(stdout, "  -wmp2       Alternative contours are derived using this width in degrees\n");
    fprintf(stdout, "  -wip        Same as -wmp, but now for the interpulse\n");
//START REGION DEVELOP
    fprintf(stdout, "  -elliptical For the beam widths contours, take into account this ellipticity.\n");
    fprintf(stdout, "              A negative value means elongation in the rotation direction, while\n");
    fprintf(stdout, "              positive means elongation in the direction of rotational meridian.\n");
    fprintf(stdout, "              The opening angle contours show the semi minor axis of the beam.\n");
//START REGION RELEASE
//START REGION DEVELOP
    fprintf(stdout, "\nSpecify a different type of search-grid (the meaning of -B changes):\n");
    fprintf(stdout, "  -zeta       Use alpha+beta values (in deg) instead of beta values in grid\n"); 
    fprintf(stdout, "  -slope      Use sin(alpha)/sin(beta) values instead of beta values in grid\n");
    fprintf(stdout, "  -islope     Use sin(beta)/sin(alpha) values instead of beta values in grid\n"); 
//START REGION RELEASE
    fprintf(stdout, "\nOptions affecting grid-search operation:\n");
//START REGION DEVELOP
    fprintf(stdout, "  -4          Do normal search + another with 1/4 of the l0 and pa0 range\n");
    fprintf(stdout, "              (better chance to find solutions). This will take ~twice as long.\n"); 
//START REGION RELEASE
    fprintf(stdout, "  -best       Show best solution and quit program\n"); 
    fprintf(stdout, "  -brute      Do normal search, then create a lattice of 16 points in l0, pa0\n");
    fprintf(stdout, "              space with point separations given by dl0 and dpa0 centered at l0\n");
    fprintf(stdout, "              and pa0.  This has a better chance to find solutions, but it will\n");
    fprintf(stdout, "              take ~17 times as long to run.\n"); 
    fprintf(stdout, "  -macro      Execute commands in specified file rather than via manual input\n");
    fprintf(stdout, "\nGeneral fit options:\n");   
//START REGION DEVELOP
    fprintf(stdout, "  -algorithm          Specify type of downhill-simplex algorithm to use\n");
    fprintf(stdout, "                      (0=def, 1=NR)\n"); 
//START REGION RELEASE
    fprintf(stdout, "  -autoopm            Shift each pa point by 90 degrees to see if it fits better\n"); 
    fprintf(stdout, "  -fitdh \"l h dh\"     Fit for emission height difference (as fraction of light\n");
    fprintf(stdout, "                      cylinder) at this pulse longitude l (deg), initial guess h\n"); 
    fprintf(stdout, "                      and initial step size dh.\n"); 
    fprintf(stdout, "  -maxdl              Set maximum allowed deviation of the position of the\n");
    fprintf(stdout, "                      magnetic axis from start value [def=%.1f deg].\n", fitterinfo.max_l0_diff); 
    fprintf(stdout, "  -forcepa \"l PA dPA\" Forces fit to go through PA with error bar dPA\n");
    fprintf(stdout, "                      at longitude l.\n");
    fprintf(stdout, "  -opm \"l o\"          Put an OPM at longitude l with offset o in degrees in\n");
    fprintf(stdout, "                      the model. You can use this option multiple times.\n"); 
    fprintf(stdout, "  -tol                Set tolerance for fitting (default is %f, and 1000x\n", ftol); 
    fprintf(stdout, "                      better when not doing a gridsearch).\n"); 
    fprintf(stdout, "\nOptions operating on input file:\n");
    fprintf(stdout, "  -dh \"l h\"   Put at pulse longitude l (in degrees) a height difference h\n");
    fprintf(stdout, "              (as fraction of the light cylinder) in the model.\n"); 
    fprintf(stdout, "  -dl         Apply this longitude shift in degrees to data\n"); 
    fprintf(stdout, "  -mc         When reading in the PA-values, each PA is taken to be a value\n");
    fprintf(stdout, "              from a Gaussian distribution defined by its error-bar. The\n");
    fprintf(stdout, "              errorbar is set to a fixed value. This allows Monte-Carlo\n");
    fprintf(stdout, "              type of analysis.\n");
//START REGION DEVELOP
    fprintf(stdout, "  -fixseed    Makes the random numbers reproducable from run to run\n");
//START REGION RELEASE
    fprintf(stdout, "  -paerrfac   Multiply all PA errorbars with this factor\n");
//START REGION DEVELOP
    fprintf(stdout, "\nOptions included by Math:\n");   
    //    fprintf(stdout, "  -nocounters Do not show progress counters\n"); 
    fprintf(stdout, "  -math       Math's secret package option to include :\n"); 
    fprintf(stdout, "                 :zoom in pulse phase in PPA fitting window\n");
    fprintf(stdout, "                 : -sigma \n");
    fprintf(stdout, "                 :fit OPM automatically\n");
    fprintf(stdout, "                 : -iformat 10\n");
    fprintf(stdout, "                 :set device 1&2 = /xs\n");
    fprintf(stdout, "  -mathnoautojump        Same as -math but no auto OPM jump fit\n");
    fprintf(stdout, "  -wmp_math   Use input W, W_err, rho, rho_err in the least chi2 fit (all in\n");
    fprintf(stdout, "              degrees).\n");
    fprintf(stdout, "  -makeps     Set device 1&2 = chi2.ps and pa.ps.\n\n");
    fprintf(stdout, "  -calcprob   Convert chi2 into (unnormalised) probability, using sum of\n");
    fprintf(stdout, "              exp(-0.5*relative_chi2). Output file 'prob.dat'\n");
    fprintf(stdout, "  -calcprob2  Convert chi2 into (unnormalised) probability, using best\n");
    fprintf(stdout, "              value of exp(-0.5*relative_chi2). Output file 'prob.dat'\n");
    fprintf(stdout, "  -sigma      Use 1,2,3 sigma level for contour level***requires Math's\n");
    fprintf(stdout, "              calcchi program. Default is 2,3,4 times the best chi2\n");
//START REGION RELEASE
    fprintf(stdout, "\nOther functionality:\n");   
    fprintf(stdout, "  -contcol    \"edg1 edg1+ edg1- edg2 edg2+ edg2- fi fi+ fi- l0 l0+ l0- P N1 N2\"\n");
    fprintf(stdout, "              Rather than doing a fit/showing a chi^2 grid, construct a\n"); 
    fprintf(stdout, "              collection of contours which can be used with the R option in\n");
    fprintf(stdout, "              interactive mode. The observed pulse starts at pulse longitude\n");
    fprintf(stdout, "              edg1 with pos/neg error edg1+/edg1- in deg. Similarly the end of\n");
    fprintf(stdout, "              the pulse is defined. The fiducial plane is at fi (in deg) with \n");
    fprintf(stdout, "              errorbars and the inflection point at l0 (in deg) with errorbars. \n");
    fprintf(stdout, "              The period of the pulsar is P seconds. All numbers are expected to\n");
    fprintf(stdout, "              be positive. N1 and N2 are the number of W and rho values generated.\n");
    fprintf(stdout, "              the overall range of pulse longitudes covered by the open field line\n");
    fprintf(stdout, "              region is based on the furthest edge from the fiducial plane position.\n");
//START REGION DEVELOP
    fprintf(stdout, "  -contcol2   Requires the same input parameters as -contcol and it behaves\n");
    fprintf(stdout, "              similarly. The output will be less efficient in the sense that\n");
    fprintf(stdout, "              duplicated or nearly duplicated entries are not avoided. There are N1\n");
    fprintf(stdout, "              steps made in each position of the edge of the pulse and N2 steps in\n");
    fprintf(stdout, "              the fiducial plane and inflection point position.\n");
//START REGION RELEASE
    fprintf(stdout, "\nPlot options:\n");   
    fprintf(stdout, "  -device1         Specify the pgplot device for the chi^2 grid\n"); 
    fprintf(stdout, "  -device2         Specify the pgplot device for the PA-swing\n"); 
    fprintf(stdout, "  -device1res      Change the resolution of device 1 [def=\"%d %d\"]\n", gridDevice_resx, gridDevice_resy); 
    fprintf(stdout, "  -device2res      Change the resolution of device 2 [def=\"%d %d\"]\n", paDevice_resx, paDevice_resy); 
    fprintf(stdout, "  -devicenores     Use the default resolution decided by pgplot\n"); 
    fprintf(stdout, "  -cont            Enables the plotting of rho contours, without waiting for\n"); 
    fprintf(stdout, "                   user input in interactive mode\n"); 
//START REGION DEVELOP
    fprintf(stdout, "  -nographics      Do not show any graphics\n"); 
    fprintf(stdout, "  -pgplot          Run additional pgplot commands\n"); 
    fprintf(stdout, "  -pgplotlist      Show supported pgplot commands for the -pgplot option\n"); 
    fprintf(stdout, "  -prefix          Default is chi2 (so output is chi2.ps and chi2.pgplot)\n"); 
//START REGION RELEASE
    fprintf(stdout, "  -showwedge       Show a colour wedge indicating the reduced chi^2 scale\n");
    fprintf(stdout, "  -showwedge_label Specify label to be shown next to the wedge\n");
    fprintf(stdout, "\nFile options:\n");   
//START REGION DEVELOP
    fprintf(stdout, "  -iformat    %d = PPOL_format, %d = PPOL_SHORT_format\n", PPOL_format, PPOL_SHORT_format);
//START REGION RELEASE
    fprintf(stdout, "  -load       Load the specified dump file containing the chi^2 grid and\n"); 
    fprintf(stdout, "              PA-swing. A dumpfile is automatically generated by ppolFit\n"); 
    fprintf(stdout, "              after a grid-search. The default output name is %s.\n", dumpfile); 
    fprintf(stdout, "  -save       Change the default name of the output dump file\n"); 
//START REGION DEVELOP
    fprintf(stdout, "  -savegrids  Write out the various grids (chi2 etc) to the specified file\n"); 
//START REGION RELEASE
    printf("\n");
    printf("Please use the appropriate citation when using results of this software in your publications:\n\n");
    printf("More information about fitting position-angle swings and using the beam-width information can be found in:\n");
    printf(" - 	Rookyard et al. 2015, MNRAS, 446, 3367\n\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {
    int index;
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-a") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &alpha0, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	alphaset = 1;
	i++;
      }else if(strcmp(argv[i], "-b") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &beta0, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	betaset = 1;
	i++;
      }else if(strcmp(argv[i], "-device1") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%100s", device1, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-device2") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%100s", device2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;      
      }else if(strcmp(argv[i], "-device1res") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &gridDevice_resx, &gridDevice_resy, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-device2res") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &paDevice_resx, &paDevice_resy, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-devicenores") == 0) {
	devicenores = 1;
      }else if(strcmp(argv[i], "-wmp") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &fitterinfo.pulse_width, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	calculate_beam_widths = 1;
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-wmp_math") == 0) {
	fitterinfo.do_rho_width = 1;
	contour_plot = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %lf %lf",&fitterinfo.pulse_width,&fitterinfo.sigma_width,&fitterinfo.rho_bcw,&fitterinfo.sigma_rho, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	calculate_beam_widths = 1;
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-macro") == 0) {
	macrofilename = i+1;
	i++;
      }else if(strcmp(argv[i], "-dl") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &add_longitude_shift, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-dh") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &fitterinfo.add_height_longitude, &dh0, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-fitdh") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %lf", &fitterinfo.add_height_longitude, &dh0, &ddh0, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-iformat") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &iformat, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-algorithm") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &amoeba_algorithm, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-wip") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &pulse_width2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	calculate_interpulse_widths = 1;
	i++;
      }else if(strcmp(argv[i], "-wmp2") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &pulse_width2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	calculate_interpulse_widths = 2;
	i++;
      }else if(strcasecmp(argv[i], "-paerrfac") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &paErrorFac, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-elliptical") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &elliptical_beam, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-contcol") == 0
//START REGION DEVELOP
	       || strcmp(argv[i], "-contcol2") == 0
//START REGION RELEASE
	       ) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", &contcol_edge1, &contcol_edge1plus, &contcol_edge1min, &contcol_edge2, &contcol_edge2plus, &contcol_edge2min, &contcol_phi0, &contcol_phi0plus, &contcol_phi0min, &contcol_l0, &contcol_l0plus, &contcol_l0min, &contcol_p, &contcol_nr_edge_steps, &contcol_nr_fid_plane_steps, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	//	j = sscanf(argv[i+1], "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", &contcol_edge1, &contcol_edge1plus, &contcol_edge1min, &contcol_edge2, &contcol_edge2plus, &contcol_edge2min, &contcol_phi0, &contcol_phi0plus, &contcol_phi0min, &contcol_l0, &contcol_l0plus, &contcol_l0min, &contcol_p, &contcol_nr_edge_steps, &contcol_nr_fid_plane_steps);
	//	if(j != 15) {
	//	  printerror(application.verbose_state.debug, "Cannot parse %s option, need 15 values.", argv[i]);
	//	  return 0;
	//	}
	doContourRange = 1;
//START REGION DEVELOP
	if(strcmp(argv[i], "-contcol2") == 0) {
	  doContourRange = 2;
	}
//START REGION RELEASE
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-l") == 0 || strcmp(argv[i], "-l0") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &l0, &dl0, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	lset = 1;
	i++;
      }else if(strcmp(argv[i], "-maxdl") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &fitterinfo.max_l0_diff, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-4") == 0) {
	doQuarterRange = 1;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-brute") == 0) {
	doBruteForce = 1;
      }else if(strcmp(argv[i], "-opm") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &fitterinfo.jump_longitude[fitterinfo.nrJumps], &fitterinfo.jump_offset[fitterinfo.nrJumps], NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	fitterinfo.nrJumps++;
	i++;
      }else if(strcmp(argv[i], "-forcepa") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf %lf", &fitterinfo.force_l, &fitterinfo.force_pa, &fitterinfo.force_dpa, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	fitterinfo.force_set = 1;
	i++;
      }else if(strcmp(argv[i], "-pa") == 0 || strcmp(argv[i], "-pa0") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &pa0, &dpa0, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	paset = 1;
	i++;
      }else if(strcmp(argv[i], "-A") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &alphastart, &alphaend, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-B") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf %lf", &betastart, &betaend, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-tol") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%lf", &ftol, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "showwedge_label") == 0) {
	showwedge_label = argv[i+1];
	i++;
      }else if(strcasecmp(argv[i], "-mc") == 0) {
	doMC = 1;
      }else if(strcasecmp(argv[i], "-fixseed") == 0) {
	fixseed = 1;
      }else if(strcmp(argv[i], "-cont") == 0) {
	contour_plot = 1;
	//      }else if(strcmp(argv[i], "-nocounters") == 0) {
	//	nocounters = 1;
      }else if(strcmp(argv[i], "-best") == 0) {
	onlyshowbest = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-nographics") == 0) {
	showgraphics = 0;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-autoopm") == 0) {
	fitterinfo.autojump = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-sigma") == 0) {
	use_sigma_level = 1;
	suppress_greyscale = 2;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-math") == 0) {
	fitterinfo.math_setting = 1;
	nruserContours = -1;
	use_sigma_level = 1;
	suppress_greyscale = 2;
	fitterinfo.autojump = 1;
	iformat = 10;
	sprintf(device1, "/xs");
	sprintf(device2, "/xs");
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-mathnoautojump") == 0) {
	fitterinfo.math_setting = 1;
	nruserContours = -1;
	use_sigma_level = 1;
	suppress_greyscale = 2;
	fitterinfo.autojump = 0;
	iformat = 10;
	sprintf(device1, "/xs");
	sprintf(device2, "/xs");
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-makeps") == 0) {
	sprintf(device1, "chi2.ps/ps");
	sprintf(device2, "pa.ps/ps");
	fitterinfo.math_setting = 2;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-zeta") == 0) {
	usezeta = 1;
      }else if(strcmp(argv[i], "-slope") == 0) {
	useslope = 1;
	usezeta = 0;
      }else if(strcmp(argv[i], "-islope") == 0) {
	useslope = 2;
	usezeta = 0;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-showwedge") == 0) {
	showwedge = 1;
      }else if(strcmp(argv[i], "-showwedge_label") == 0) {
	showwedge_label = argv[i+1];
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i],"-calcprob")==0){
	calc_prob = 1;
      }else if(strcmp(argv[i],"-calcprob2")==0){
	calc_prob = 2;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-load") == 0) {
	loadresults = 1;
	strcpy(dumpfile, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-save") == 0) {
	strcpy(dumpfile, argv[i+1]);
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-savegrids") == 0 || strcmp(argv[i], "-savegrid") == 0) {
	dumpgrids = i+1;  // So this variable points to the index of the specified file name
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-printnofit") == 0) {
	printnofit = 1;
      }else if(strcmp(argv[i], "-printfit") == 0) {
	printfit = 1;
      }else if(strcmp(argv[i], "-printsmooth") == 0) {
	printsmooth = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-pgplot") == 0) {   /* Ignore pgplot commants */
	i++;
      }else if(strcmp(argv[i], "-pgplotlist") == 0) {   
	printPgplotInterpratorHelp();
	terminateApplication(&application);
	return 0;
      }else if(strcmp(argv[i], "-prefix") == 0) {
	strcpy(prefix, argv[i+1]);
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-g") == 0) {
	GridSearch = 1;
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &nalpha, &nbeta, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else {
	/* If the option is not recognized, assume it is a filename */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "ppolFit: Unknown option: %s\n\nRun ppolFit without command line arguments to show help", argv[i]);	
	  terminateApplication(&application);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
      /*
      }else {
	if(i == argc - 1 && loadresults == 0) {  // If not loading the data, the last parameter is the PA file, so that is ok
	  break;
	}else if(i == argc - 1 && loadresults) {
	  printerror(application.verbose_state.debug, "Unknown option \"%s\".", argv[i]);
	  printerror(application.verbose_state.debug, "Note: when using the -load you're not allowed to specify the PA-swing file, as it is already stored in dump file.");
	  return 0;
	}else {
	  printerror(application.verbose_state.debug, "Unknown option \"%s\".", argv[i]);
	  return 0;
	}
      }
      */
    }
  } 

  // Make a list of contours to be used with the 'R' option.
  if(doContourRange) {
    double wopen, hem, theta_pc, rho_pc, shiftBCW;
    fprintf(stderr, "For most the most likely input values:\n");
    fprintf(stderr, "  Pulse starts/ends at:   %lf (+%f -%f) and %lf (+%f -%f) deg\n", contcol_edge1, contcol_edge1plus, contcol_edge1min, contcol_edge2, contcol_edge2plus, contcol_edge2min);
    fprintf(stderr, "  Fiducial plane at:      %lf (+%f -%f) deg\n", contcol_phi0, contcol_phi0plus, contcol_phi0min);
    wopen = fabs(contcol_edge2 - contcol_phi0);
    if(fabs(contcol_edge1 - contcol_phi0) > wopen)
      wopen = fabs(contcol_edge1 - contcol_phi0);
    wopen *= 2.0;
    fprintf(stderr, "  This make Wopen:        %lf deg\n", wopen);
    fprintf(stderr, "  Inflection point at:    %lf (+%f -%f) deg\n", contcol_l0, contcol_l0plus, contcol_l0min);
    shiftBCW = contcol_l0 - contcol_phi0;
    fprintf(stderr, "  BCW shift:              %lf deg\n", contcol_l0 - contcol_phi0);
    shiftBCW *= M_PI/180.0;
    if(contcol_edge1plus < 0 || contcol_edge1min < 0 || contcol_edge2plus < 0 || contcol_edge2min < 0 || contcol_l0min < 0 || contcol_l0plus < 0 || contcol_phi0min < 0 || contcol_phi0plus < 0) {
      printerror(application.verbose_state.debug, "ERROR ppolFit: All errors defined for the -contcol option should be positive.\n");
      return 0;
    }
    if(shiftBCW > 0) {
      hem = contcol_p*299792458.0*shiftBCW/(8.0*M_PI);
      fprintf(stderr, "  Emission height:        %lf km\n", 0.001*hem);
      theta_pc = asin(sqrt(2.0*M_PI*hem/(contcol_p*299792458.0)));
      rho_pc = theta_pc + atan(0.5*tan(theta_pc));
      fprintf(stderr, "  Opening angle beam:     %lf deg (half-opening angle)\n", rho_pc*180.0/M_PI);
    }else {
      fprintf(stderr, "  Emission height:        Inconsistent (BCW shift is negative)\n");
    }

    double rho_min, rho_max, shiftBCW_min, shiftBCW_max;
    // Max BCW shift possible
    shiftBCW = contcol_l0 + contcol_l0plus - (contcol_phi0 - contcol_phi0min);
    if(shiftBCW < 0) {
      printerror(application.verbose_state.debug, "ERROR ppolFit: No positive values for the BCW shift allowed, so rho is undefined\n");
      return 0;
    }
    shiftBCW_max = shiftBCW;
    shiftBCW *= M_PI/180.0;
    hem = contcol_p*299792458.0*shiftBCW/(8.0*M_PI);
    theta_pc = asin(sqrt(2.0*M_PI*hem/(contcol_p*299792458.0)));
    rho_pc = theta_pc + atan(0.5*tan(theta_pc));
    rho_pc *= 180.0/M_PI;
    rho_max = rho_pc;
      
    // Min BCW shift possible
    shiftBCW = contcol_l0 - contcol_l0min - (contcol_phi0 + contcol_phi0plus);
    if(shiftBCW < 0) {
      shiftBCW = 0;
    }
    shiftBCW_min = shiftBCW;
    shiftBCW *= M_PI/180.0;
    hem = contcol_p*299792458.0*shiftBCW/(8.0*M_PI);
    theta_pc = asin(sqrt(2.0*M_PI*hem/(contcol_p*299792458.0)));
    rho_pc = theta_pc + atan(0.5*tan(theta_pc));
    rho_pc *= 180.0/M_PI;
    rho_min = rho_pc;
 
    fprintf(stderr, "Predicted range in rho:    %lf ... %lf deg\n", rho_min, rho_max);

//START REGION DEVELOP
    // Simple solution: just use every posibility, which results in many duplicate, or near duplicate output lines
    if(doContourRange == 2) {
      double edge1, edge2, phi0, wopen2, l0, wopen_min, wopen_max;
      long i, j, k, l, failed_combi, rho_minmax_set;
      failed_combi = 0;
      rho_minmax_set = 0;
      for(i = 0; i < contcol_nr_edge_steps; i++) {
	edge1 = contcol_edge1 - contcol_edge1min + (contcol_edge1plus+contcol_edge1min)*(double)i/(double)(contcol_nr_edge_steps-1);
	for(j = 0; j < contcol_nr_edge_steps; j++) {
	  edge2 = contcol_edge2 - contcol_edge2min + (contcol_edge2plus+contcol_edge2min)*(double)j/(double)(contcol_nr_edge_steps-1);
	  for(k = 0; k < contcol_nr_fid_plane_steps; k++) {
	    phi0 = contcol_phi0 - contcol_phi0min + (contcol_phi0plus+contcol_phi0min)*(double)k/(double)(contcol_nr_fid_plane_steps-1);
	    for(l = 0; l < contcol_nr_fid_plane_steps; l++) {
	      l0 = contcol_l0 - contcol_l0min + (contcol_l0plus+contcol_l0min)*(double)l/(double)(contcol_nr_fid_plane_steps-1);
	      
	      wopen = edge2 - phi0;  // Right edge
	      if(wopen < 0) // Right edge left from fid plane
		wopen *= -1.0; // Right edge profile is effectively an left inner edge of beam.
	      
	      wopen2 = phi0 - edge1;  // Left edge
	      if(wopen2 < 0) // Left edge right from fid plane
		wopen2 *= -1.0; // Left edge profile is effectively an right inner edge of beam.
	      
	      if(wopen2 > wopen)
		wopen = wopen2;   // Effective profile width set by left edge
	      
	      wopen *= 2.0; // Full open field line width is twice the max offset
	      
	      if((i == 0 && j == 0 && k == 0) || wopen < wopen_min) {
		wopen_min = wopen;
	      }
	      if((i == 0 && j == 0 && k == 0) || wopen > wopen_max) {
		wopen_max = wopen;
	      }
	      
	      shiftBCW = l0 - phi0;
	      shiftBCW *= M_PI/180.0;
	      if(shiftBCW > 0) {
		hem = contcol_p*299792458.0*shiftBCW/(8.0*M_PI);
		theta_pc = asin(sqrt(2.0*M_PI*hem/(contcol_p*299792458.0)));
		rho_pc = theta_pc + atan(0.5*tan(theta_pc));
		rho_pc *= 180.0/M_PI;
		
		if(rho_minmax_set == 0 || rho_pc < rho_min) {
		  rho_min = rho_pc;
		  rho_minmax_set = 1;
		}
		if(rho_minmax_set == 0 || rho_pc > rho_max) {
		  rho_max = rho_pc;
		  rho_minmax_set = 1;
		}
		
		printf("%e %e\n", wopen, rho_pc);
	      }else {
		failed_combi++;
	      }
	    }
	  }
	}
      }
      fprintf(stderr, "\nRange in Wopen covered:     %lf ... %lf deg\n", wopen_min, wopen_max);
      fprintf(stderr, "Actual range in rho covered:    %lf ... %lf deg\n", rho_min, rho_max);
      if(failed_combi) {
	fprintf(stderr, "Some of the combinations in the pulse edges/fiducial plane/inflection point positions did not result in positive BCW effect.\n");
      }
    }else {
//START REGION RELEASE
      // Smarter solution: define regular rho grid (use actually regular BCW steps, since finding BCW shift from given rho is not particularly easy) and explore each possibility in turn
      long k, l;
      double phi0, edge1, edge2, wopen2, wopen_min, wopen_max;
      for(l = 0; l < contcol_nr_fid_plane_steps; l++) {
	shiftBCW = shiftBCW_min + (shiftBCW_max-shiftBCW_min)*(double)l/(double)(contcol_nr_fid_plane_steps-1);
	//This defines rho
	shiftBCW *= M_PI/180.0;
	hem = contcol_p*299792458.0*shiftBCW/(8.0*M_PI);
	theta_pc = asin(sqrt(2.0*M_PI*hem/(contcol_p*299792458.0)));
	rho_pc = theta_pc + atan(0.5*tan(theta_pc));
	rho_pc *= 180.0/M_PI;

	//	fprintf(stderr, "Try out rho=%e\n", rho_pc);

	// Now rho is fixed, explore each possible value of the fiducial plane position
	// Main thing to work out now the possible range in Wopen, which depends on the fiducial plane position and the pulse edges
	wopen_min = wopen_max = -1.0;
	for(k = 0; k < contcol_nr_fid_plane_steps*10; k++) {
	  phi0 = contcol_phi0 - contcol_phi0min + (contcol_phi0plus+contcol_phi0min)*(double)k/(double)(contcol_nr_fid_plane_steps*10-1);
	  // For a given fiducial plane position and rho, the inflection point is defined
	  l0 = phi0 + shiftBCW*180.0/M_PI;
	  //	  fprintf(stderr, "Try out phi0=%e, which implies l0=%e\n", phi0, l0);
	  // Check if for this fiducial plane position the rho value is allowed, or if required BCW shift is outside allowed range
	  if(l0 >= contcol_l0 - contcol_l0min && l0 <= contcol_l0 + contcol_l0plus) {
	    // With a fixed rho and fiducial plane position, work out the possible range in Wopen
	    for(i = 0; i < contcol_nr_edge_steps; i++) {
	      edge1 = contcol_edge1 - contcol_edge1min + (contcol_edge1plus+contcol_edge1min)*(double)i/(double)(contcol_nr_edge_steps-1);
	      for(j = 0; j < contcol_nr_edge_steps; j++) {
		edge2 = contcol_edge2 - contcol_edge2min + (contcol_edge2plus+contcol_edge2min)*(double)j/(double)(contcol_nr_edge_steps-1);

		wopen = edge2 - phi0;  // Right edge
		if(wopen < 0) // Right edge left from fid plane
		  wopen *= -1.0; // Right edge profile is effectively an left inner edge of beam.
		
		wopen2 = phi0 - edge1;  // Left edge
		if(wopen2 < 0) // Left edge right from fid plane
		  wopen2 *= -1.0; // Left edge profile is effectively an right inner edge of beam.
		
		if(wopen2 > wopen)
		  wopen = wopen2;   // Effective profile width set by left edge
		
		wopen *= 2.0; // Full open field line width is twice the max offset

		if(wopen_min < 0 || wopen < wopen_min) {
		  wopen_min = wopen;
		}
		if(wopen_max < 0 || wopen > wopen_max) {
		  wopen_max = wopen;
		}
	      }
	    }
	  }
	}
	if(wopen_min < 0 || wopen_max < 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Determining range in wopen failed. Possibly N2 is too small, which can especially be an issue if the error on l0 is small.\n");
	  return 0;
	}
	// So for the current value of rho we worked out the range in wopen. Now simply write out equally sampled W values
	for(i = 0; i < contcol_nr_edge_steps; i++) {
	  wopen = wopen_min + (wopen_max-wopen_min)*(double)i/(double)(contcol_nr_edge_steps-1);
	  printf("%e %e\n", wopen, rho_pc);
	}
      }
//START REGION DEVELOP
    }
//START REGION RELEASE
    fprintf(stderr, "\nThe output can be re-directed to a file and used with the R option of ppolFit in interactive mode.\n");

    terminateApplication(&application);
    return 0;
  }   // End of if doContourRange
//START REGION RELEASE

  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) > 0 && loadresults) {
    printerror(application.verbose_state.debug, "ERROR ppolFit: Unknown option \"%s\".", getNextFilenameFromList(&application, argv, application.verbose_state));
    printerror(application.verbose_state.debug, "Note: when using the -load you're not allowed to specify the PA-swing file, as it is already stored in dump file.");
    return 0;
  }
  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) > 1) {
    printerror(application.verbose_state.debug, "ERROR ppolFit: Only one input file can be selected on the command line.");
    return 0;
  }


  if(lset == 0) {
    if(loadresults == 0) {
      printwarning(application.verbose_state.debug, "  Assumed -l '180 90'");
      l0 = 180;
      dl0 = 90;
    }
  }
  if(paset == 0) {
    if(loadresults == 0) {
      printwarning(application.verbose_state.debug, "  Assumed -pa '0 90'");
      pa0 = 0;
      dpa0 = 90;
    }
  }

  gsl_rng *rand_num_gen;
  const gsl_rng_type *rand_num_gen_type;
  long idnum;
  gsl_rng_env_setup();    /* Set the default generators, can be influenced by environment variables */
  rand_num_gen_type = gsl_rng_default;
  rand_num_gen = gsl_rng_alloc(rand_num_gen_type);
  if(fixseed)
    idnum = 1;
  else
    randomize_idnum(&idnum);
  gsl_rng_set(rand_num_gen, idnum);


  fit_dh0 = dh0;
  if(printnofit == 1) {
    for(i = 0; i < 360; i++)
      printf("%d %lf\n", i, paswing_double(alpha0, beta0, i, pa0, l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, dh0));
    terminateApplication(&application);
    gsl_rng_free(rand_num_gen);
    return 0;
  }

  if(loadresults == 0) {
    //    cleanPSRData(&datain, application.verbose_state);
    iformat = guessPSRData_format(argv[argc-1], 0, application.verbose_state);
    if(iformat == -2 || iformat == -3)
      return 0;
    if(isValidPSRDATA_format(iformat) == 0) {
      printerror(application.verbose_state.debug, "ERROR ppolFit: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", argv[argc-1]);
      return 0;
    }
    if(iformat != PPOL_format && iformat != PPOL_SHORT_format) {
      printerror(application.verbose_state.debug, "ppolFit: Input file does not appear to be in a ppol format.");
      return 0;
    }
    if(!openPSRData(&datain, argv[argc-1], iformat, 0, 1, 0, application.verbose_state))
      return 0;
    fitterinfo.NrDataPoints = filterPApoints(&datain, application.verbose_state);
    if(fitterinfo.NrDataPoints == 0) {
      printerror(application.verbose_state.debug, "ppolFit: No significant pa points in data.");
      return 0;
    }

    if(paErrorFac != 1.0) {
      if(application.verbose_state.verbose)
	printf("ppolFit: Multiplying input PA errors with %lf.\n", paErrorFac);
      for(i = 0; i < datain.NrBins; i++) {
	//	datain.data_dpa[i] *= paErrorFac;
	if(datain.poltype == POLTYPE_ILVPAdPA || datain.poltype == POLTYPE_ILVPAdPATEldEl) {
	  datain.data[i+4*datain.NrBins] *= paErrorFac;
	}else {
	  datain.data[i+1*datain.NrBins] *= paErrorFac;
	}
      }
    }

    if(doMC) {
      if(application.verbose_state.verbose)
	printf("ppolFit: Randomizing input PAs and set errors to one.\n");
      for(i = 0; i < datain.NrBins; i++) {
	if(datain.poltype == POLTYPE_ILVPAdPA || datain.poltype == POLTYPE_ILVPAdPATEldEl) {
	  datain.data[i+3*datain.NrBins] += datain.data[i+4*datain.NrBins]*gsl_ran_gaussian(rand_num_gen, 1.0);
	  datain.data[i+4*datain.NrBins] = 1;
	}else {
	  datain.data[i                ] += datain.data[i+1*datain.NrBins]*gsl_ran_gaussian(rand_num_gen, 1.0);
	  datain.data[i+1*datain.NrBins] = 1;
	}
      }
    }


    if(datain.tsampMode != TSAMPMODE_LONGITUDELIST || datain.tsamp_list == NULL) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "ERROR ppolFit: Pulse longitudes were expected to be defined.");
      return 0;
    }
    fitterinfo.data_l = datain.tsamp_list;
    if(datain.poltype == POLTYPE_ILVPAdPA || datain.poltype == POLTYPE_ILVPAdPATEldEl) {
      fitterinfo.data_pa  = &(datain.data[3*datain.NrBins]);
      fitterinfo.data_dpa = &(datain.data[4*datain.NrBins]);
    }else {
      fitterinfo.data_pa  = datain.data;
      fitterinfo.data_dpa = &(datain.data[datain.NrBins]);
    }

    if(add_longitude_shift != 0.0) {
      if(application.verbose_state.verbose)
	printf("ppolFit: Add longitude offset of %lf to PA points.\n", add_longitude_shift);
      for(i = 0; i < fitterinfo.NrDataPoints; i++)
	fitterinfo.data_l[i] += add_longitude_shift;
    }

  }else {
    if(paErrorFac != 1.0) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "The -paerrfac option cannot be used together with the -load option");
      return 0;
    }
    if(add_longitude_shift != 0.0) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "The -dl option cannot be used together with the -load option");
      return 0;
    }
    if(doMC) {
      fflush(stdout);
      printerror(application.verbose_state.debug, "The -mc option cannot be used together with the -load option");
      return 0;
    }
  }



  if(loadresults) {
    fin = fopen(dumpfile, "rb");
    if(fin == NULL) {
      printerror(application.verbose_state.debug, "Cannot open %s", dumpfile); 
      return 0;
    }
    if(fread(txt, sizeof(char), 12, fin) != 12) {printerror(application.verbose_state.debug, "Read error from %s.", dumpfile); return 0; }
    txt[13] = 0;
    if(strcmp(txt, "fit_paswing") == 0) {
      if(fread(&version, sizeof(int), 1, fin) != 1) { printerror(application.verbose_state.debug, "Read error from %s.", dumpfile); return 0; }
    }else {
      printwarning(application.verbose_state.debug, "WARNING Unrecognized file version in %s. I assume it is written by an obsolete version of ppolFit.", dumpfile);
      rewind(fin);
      version = 0;
    }
    if(application.verbose_state.verbose)
      printf("File version of %s: %d\n", dumpfile, version);
    if(version < 0 || version > 5) {
      printerror(application.verbose_state.debug, "ERROR: ppolFit does not recognize file version %d", version);
    }
    if(version >= 2) {
      fread(&alphastart, sizeof(double), 1, fin);
      fread(&alphaend, sizeof(double), 1, fin);
      fread(&betastart, sizeof(double), 1, fin);
      fread(&betaend, sizeof(double), 1, fin);
      fread(&nalpha, sizeof(int), 1, fin);
      fread(&nbeta, sizeof(int), 1, fin);

      fread(&l0, sizeof(double), 1, fin);
      fread(&dl0, sizeof(double), 1, fin);
      fread(&pa0, sizeof(double), 1, fin);
      fread(&dpa0, sizeof(double), 1, fin);
    }else {
      float dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      alphastart = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      alphaend = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      betastart = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      betaend = dummyf;
      fread(&nalpha, sizeof(int), 1, fin);
      fread(&nbeta, sizeof(int), 1, fin);

      fread(&dummyf, sizeof(float), 1, fin);
      l0 = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      dl0 = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      pa0 = dummyf;
      fread(&dummyf, sizeof(float), 1, fin);
      dpa0 = dummyf;
    }

    if(version >= 2) {
//START REGION DEVELOP
      if(1 == 2) {
//START REGION RELEASE
	int dummy;
	fread(&dummy, sizeof(int), 1, fin);
	fread(&dummy, sizeof(int), 1, fin);
//START REGION DEVELOP
      }else {
	fread(&usezeta, sizeof(int), 1, fin);
	fread(&useslope, sizeof(int), 1, fin);
      }
//START REGION RELEASE
    }
    if(version >= 3) {
      if(fitterinfo.nrJumps != 0) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "You cannot use -opm together with the -load option as information is stored in the dump file");
	return 0;
      }
      fread(&fitterinfo.nrJumps, sizeof(int), 1, fin);
      for(i = 0; i < fitterinfo.nrJumps; i++) {
	fread(&fitterinfo.jump_longitude[i], sizeof(double), 1, fin);
	fread(&fitterinfo.jump_offset[i], sizeof(double), 1, fin);
      }
      if(fitterinfo.autojump != 0) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "You cannot use -autoopm together with the -load option as information is stored in the dump file");
	return 0;
      }
      fread(&fitterinfo.autojump, sizeof(int), 1, fin);
    }

    fread(&fitterinfo.NrDataPoints, sizeof(int), 1, fin);
    if(fitterinfo.NrDataPoints < 0) {
      printerror(application.verbose_state.debug, "ppolFit: Cannot load %d data points, something appears to be wrong with dumpfile\n", fitterinfo.NrDataPoints);
      return 0;
    }
    if(fitterinfo.NrDataPoints > 10000) {
      printwarning(application.verbose_state.debug, "WARNING ppolFit: Going to load %d PA data points, which seems a large number. Maybe something is wrong with dumpfile?\n", fitterinfo.NrDataPoints);
    }
    if(application.verbose_state.debug) {
      fprintf(stdout, "  Going to load %d data points\n", fitterinfo.NrDataPoints);
      fflush(stdout);
    }
    fitterinfo.data_l = (double *)malloc(fitterinfo.NrDataPoints*sizeof(double));
    fitterinfo.data_pa = (float *)malloc(fitterinfo.NrDataPoints*sizeof(float));
    fitterinfo.data_dpa = (float *)malloc(fitterinfo.NrDataPoints*sizeof(float));  
    if(fitterinfo.data_l == NULL || fitterinfo.data_pa == NULL || fitterinfo.data_dpa == NULL) {
      printerror(application.verbose_state.debug, "Cannot allocate memory"); 
      return 0;    
    } 
    for(i = 0; i < fitterinfo.NrDataPoints; i++) {
      float dummy_float;
      if(version >= 4) {
	fread(&fitterinfo.data_l[i], sizeof(double), 1, fin);
      }else {
	fread(&dummy_float, sizeof(float), 1, fin);
	fitterinfo.data_l[i] = dummy_float;
      }
      fread(&fitterinfo.data_pa[i], sizeof(float), 1, fin);
      fread(&fitterinfo.data_dpa[i], sizeof(float), 1, fin);
    }

    if(version >= 5) {
      int dummyint;
      if(fread(&dummyint, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Read error."); return 0; }
      if(dummyint < 0) {
	printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot load a string of %d bytes, something appears to be wrong with dumpfile\n", dummyint);
	return 0;
      }
      if(dummyint > MaxStringLength) {
	printerror(application.verbose_state.debug, "ERROR ppolFit: Command line stored in dump file appears to be very long (%d bytes). Maybe something is wrong with dumpfile?\n", dummyint);
      }
      if(fread(txt, 1, dummyint, fin) != dummyint) {printerror(application.verbose_state.debug, "Read error."); return 0; }
      txt[dummyint] = 0;
      if(application.verbose_state.verbose) {
	printf("Dump file created with command: %s\n", txt);
      }
    }

    GridSearch = 1; 
  }

  if(alphaset == 0 || betaset == 0) {
    if(GridSearch == 0) {
      printerror(application.verbose_state.debug, "Need to specify the -g option to do a grid search over alpha and beta, or use the -a and -b option to fix their values.");
      return 0;
    }
  }
  /*
  if(lset == 0 || paset == 0) {
    if(loadresults == 0) {
      printerror(application.verbose_state.debug, "Need to specify search range with -l and -pa.");
      return 0;
    }
    }*/

//START REGION DEVELOP
  if( (calc_prob == 1) ||  (calc_prob == 2) )  {
    prob_a = (double *)malloc((nalpha)*sizeof(double));
    prob_b = (double *)malloc((nbeta)*sizeof(double));
    if( (prob_a == NULL)||(prob_b == NULL) ) {
      printerror(application.verbose_state.debug, "Cannot allocate memory"); 
      return 0;    
    }
  }
//START REGION RELEASE


  if(GridSearch == 1) {
    chigrid = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    l0grid = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    pa0grid = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    dhgrid = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    rhogrid = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    rhogrid2 = (float *)malloc((nalpha)*(nbeta)*sizeof(float));
    if(chigrid == NULL || l0grid == NULL || pa0grid == NULL || dhgrid == NULL || rhogrid == NULL || rhogrid2 == NULL) {
      printerror(application.verbose_state.debug, "Cannot allocate memory"); 
      return 0;    
    }

    
    if(loadresults == 0) {
      for(i = 0; i < nalpha; i++) {
	for(j = 0; j < nbeta; j++) {
	  chigrid[nalpha*j+i] = 1e10;
	  l0grid[nalpha*j+i] = 1e10;
	  pa0grid[nalpha*j+i] = 1e10;
	  dhgrid[nalpha*j+i] = 1e10;
	}
      }
      for(i = 0; i < nalpha; i++) {
	for(j = 0; j < nbeta; j++) {
	  alpha0 = i*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
	  beta0 = j*(betaend-betastart)/(double)(nbeta-1)+betastart;
	  convertAlphaBeta(&alpha0, &beta0
//START REGION DEVELOP
			   , useslope, usezeta
//START REGION RELEASE
			   );
	  /* Patrick: Here you can do mrq stuff of
             ppolFit_mrq_stuff_doesnt_work (except it doesn't at
             this point).  There is probably no point in doing it
             ever..... */
	  DoFitting(alpha0, beta0, pa0, dpa0, l0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 0, stdout, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
	  chi_old = chi;
	  bestl0 = fit_l0;
	  bestpa0 = fit_pa0;
	  //	  beshdh = fit_dh0;
//START REGION DEVELOP
	  if(doQuarterRange) {
	    DoFitting(fit_alpha, fit_beta, pa0, 0.25*dpa0, l0, 0.25*dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 0, stdout, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
	    /* This should be a > sign. If the new solution is worse, take old one. */
	    if(chi > chi_old) {
	      chi = chi_old;
	    }
	  }
//START REGION RELEASE
	  if(doBruteForce) {
	    int bruteForceSignPa0, bruteForceSignL0, pa0step, l0step;
	    double newpa0, newl0;
	    for(bruteForceSignPa0 = -1; bruteForceSignPa0 <= 1; bruteForceSignPa0 += 2) {
	      for(bruteForceSignL0 = -1; bruteForceSignL0 <= 1; bruteForceSignL0 += 2) {
		for(pa0step = 1; pa0step <= 3; pa0step += 2) {
		  for(l0step = 1; l0step <= 3; l0step += 2) {
		    newpa0 = pa0+0.5*(double)bruteForceSignPa0*((double)pa0step*dpa0);
		    newl0 = l0+0.5*(double)bruteForceSignL0*((double)l0step*dl0);
		    /*		    printf("\nXXXX %f %f\n", newl0, newpa0); */
		    DoFitting(fit_alpha, fit_beta, newpa0, dpa0, newl0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 0, stdout, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
		    /* This should be a > sign. If the new solution is worse, take old one. */
		    if(chi < chi_old) {
		      chi_old = chi;
		      bestl0 = fit_l0;
		      bestpa0 = fit_pa0;
		      //		      beshdh = fit_dh0;
		    }
		  }
		}
	      }
	    }
	  }
	  chi = chi_old;
	  chigrid[nalpha*j+i] = chi/(double)(fitterinfo.NrDataPoints-nrfitparams);
	  l0grid[nalpha*j+i] = bestl0;
	  /*
	  double alpha_tst, beta_tst;
	  alpha_tst = i*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
	  beta_tst = j*(betaend-betastart)/(double)(nbeta-1)+betastart;
	  if(alpha_tst > 95 && alpha_tst < 105 && beta_tst < -7 && beta_tst > -8) {
	    fprintf(stderr, "XXXXXX %f\n", l0grid[nalpha*j+i]);
	    }*/
	  pa0grid[nalpha*j+i] = bestpa0;
	  dhgrid[nalpha*j+i] = fit_dh0;
	}
	if(application.verbose_state.nocounters == 0)
	  fprintf(stderr, "%.1f%%     \r",(100.0*(i+1))/(double)(nalpha));
      }
      fin = fopen(dumpfile, "wb");
      if(fin == NULL) {
	printerror(application.verbose_state.debug, "Cannot open %s", dumpfile); 
	return 0;
      }
      sprintf(txt, "fit_paswing");
      txt[12] = 0;
      if(fwrite(txt, sizeof(char), 12, fin) != 12) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      version = 5;
      if(fwrite(&version, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&alphastart, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&alphaend, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&betastart, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&betaend, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&nalpha, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&nbeta, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }

      if(fwrite(&l0, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&dl0, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&pa0, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(&dpa0, sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }

//START REGION DEVELOP
      if(1 == 2) {
//START REGION RELEASE
	int dummy;
	dummy = 0;
	if(fwrite(&dummy, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
	if(fwrite(&dummy, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
//START REGION DEVELOP
      }else {
	if(fwrite(&usezeta, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
	if(fwrite(&useslope, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      }
//START REGION RELEASE
      if(fwrite(&fitterinfo.nrJumps, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      for(i = 0; i < fitterinfo.nrJumps; i++) {
	if(fwrite(&fitterinfo.jump_longitude[i], sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
	if(fwrite(&fitterinfo.jump_offset[i], sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      }
      if(fwrite(&fitterinfo.autojump, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }

      if(fwrite(&fitterinfo.NrDataPoints, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      for(i = 0; i < fitterinfo.NrDataPoints; i++) {
	if(fwrite(&fitterinfo.data_l[i], sizeof(double), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
	if(fwrite(&fitterinfo.data_pa[i], sizeof(float), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
	if(fwrite(&fitterinfo.data_dpa[i], sizeof(float), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      }

      // Write command line
      int dummyint;
      constructCommandLineString(txt, MaxStringLength, argc, argv, application.verbose_state);
      dummyint = strlen(txt);
      if(fwrite(&dummyint, sizeof(int), 1, fin) != 1) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(txt, 1, dummyint, fin) != dummyint) {printerror(application.verbose_state.debug, "Write error."); return 0; }


      if(fwrite(chigrid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(l0grid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(pa0grid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      if(fwrite(dhgrid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Write error."); return 0; }
      fprintf(stderr, "Dumped data to %s\n", dumpfile); 
      if(fclose(fin) != 0) {printerror(application.verbose_state.debug, "File close error."); return 0; }

    }else {   /* Load grid instead of calculating it */
      if(fin == NULL) {
	printerror(application.verbose_state.debug, "ppolFit: BUG!!!!!!!!!!!!");
	return 0;
      }
      if(fread(chigrid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Read error."); return 0; }
      if(version >= 1) {
	if(fread(l0grid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Read error."); return 0; }
	if(fread(pa0grid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Read error."); return 0; }
	if(fread(dhgrid, sizeof(float), nalpha*nbeta, fin) != nalpha*nbeta) {printerror(application.verbose_state.debug, "Read error."); return 0; }
      }
      fclose(fin);
      if(application.verbose_state.verbose) fprintf(stderr, "Loaded data from %s\n", dumpfile); 
    }

    if(calculate_beam_widths != 0 || calculate_interpulse_widths != 0) {
      calcBeamWidths(nalpha, nbeta, alphastart, alphaend, betastart, betaend, 
//START REGION DEVELOP
		     usezeta, useslope, elliptical_beam, 
//START REGION RELEASE
		     calculate_beam_widths, calculate_interpulse_widths, pulse_width2, rhogrid, rhogrid2, application.verbose_state.nocounters);
    }
    
//START REGION DEVELOP
      if(dumpgrids) {  // Dump the found grids to a pulsar format so it can be plotted etc.
	int nrgrids = 6, p;
	float *data;
	data = malloc(nrgrids*sizeof(float)*nalpha*nbeta);
	if(data == NULL) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Cannot allocate memory");
	  return 0;
	}
	datafile_definition fout;
	cleanPSRData(&fout, application.verbose_state);
	if(loadresults == 0) {
	  copy_params_PSRData(datain, &fout, application.verbose_state);
	  fout.tsubMode = TSUBMODE_FIXEDTSUB;
	  if(fout.tsub_list != NULL)
	    free(fout.tsub_list);
	  fout.tsub_list = malloc(sizeof(double));
	  fout.tsub_list[0] = get_tobs(datain, application.verbose_state)/(double)nbeta;
	  if(fout.freqMode == FREQMODE_UNKNOWN) {
	    fout.freqMode = FREQMODE_UNIFORM;
	    if(fout.freqlabel_list != NULL) {
	      free(fout.freqlabel_list);
	      fout.freqlabel_list = NULL;
	    }
	    //	    fout.freq_list = malloc(2*sizeof(double));
	    //	    if(fout.freq_list == NULL) {
	    //	      fflush(stdout);
	    //	      printerror(application.verbose_state.debug, "ERROR ppolFit: Memory allocation error.");
	    //	      return 0;
	    //	    }
	    set_centre_frequency(&fout, M_PI, application.verbose_state);
	    if(set_bandwidth(&fout, -1.0, application.verbose_state) == 0) {
	      printerror(application.verbose_state.debug, "ERROR ppolFit: Bandwidth changing failed.");
	      return 0;
	    }
	  }
	}else {
	  fout.isFolded = 1;
	  fout.foldMode = FOLDMODE_FIXEDPERIOD;
	  fout.fixedPeriod = 1.0;
	  fout.freqMode = FREQMODE_UNIFORM;
	  if(fout.freqlabel_list != NULL) {
	    free(fout.freqlabel_list);
	    fout.freqlabel_list = NULL;
	  }
	  //	  fout.freq_list = malloc(2*sizeof(double));
	  //	  if(fout.freq_list == NULL) {
	  //	    fflush(stdout);
	  //	    printerror(application.verbose_state.debug, "ERROR ppolFit: Memory allocation error.");
	  //	    return 0;
	  //	  }
	  set_centre_frequency(&fout, M_PI, application.verbose_state);
	  if(set_bandwidth(&fout, -1.0, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR ppolFit: Bandwidth changing failed.");
	    return 0;
	  }
	  fout.feedtype = FEEDTYPE_UNKNOWN;
	  fout.tsubMode = TSUBMODE_FIXEDTSUB;
	  fout.tsub_list = malloc(sizeof(double));
	  fout.tsub_list[0] = M_PI;
	}
	fout.poltype = POLTYPE_UNKNOWN;
	fout.NrFreqChan = 1;
	fout.NrSubints = nbeta;
	fout.NrBins = nalpha;
	fout.NrPols = nrgrids;
	fout.xrangeset = 1;
	fout.xrange[0] = alphastart;
	fout.xrange[1] = alphaend;
	fout.yrangeset = 1;
	fout.yrange[0] = betastart;
	fout.yrange[1] = betaend;
	if(!openPSRData(&fout, argv[dumpgrids], FITS_format, 1, 0, 0, application.verbose_state))
	  return 0;
	if(writeHeaderPSRData(&fout, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Unable to write header.\n");
	  return 0;
	}
	printf("Writing grids to: %s\n", argv[dumpgrids]);
	printf("  The \"polarization channels\" are: reduced chi2, l0 (deg), PA0 (deg), dh (frac. RLC), rho (deg), 2nd calculated rho (deg)\n");
	for(i = 0; i < nbeta; i++) {
	  for(p = 0; p < nrgrids; p++) {
	    switch(p) {
	    case 0: memcpy(&data[(nrgrids*i+p)*nalpha], &chigrid[i*nalpha], nalpha*sizeof(float)); break;
	    case 1: memcpy(&data[(nrgrids*i+p)*nalpha], &l0grid[i*nalpha], nalpha*sizeof(float)); break;
	    case 2: memcpy(&data[(nrgrids*i+p)*nalpha], &pa0grid[i*nalpha], nalpha*sizeof(float)); break;
	    case 3: memcpy(&data[(nrgrids*i+p)*nalpha], &dhgrid[i*nalpha], nalpha*sizeof(float)); break;
	    case 4: memcpy(&data[(nrgrids*i+p)*nalpha], &rhogrid[i*nalpha], nalpha*sizeof(float)); break;
	    case 5: memcpy(&data[(nrgrids*i+p)*nalpha], &rhogrid2[i*nalpha], nalpha*sizeof(float)); break;
	    default: 
	      printerror(application.verbose_state.debug, "ERROR ppolFit: Bug! Undefined grid number.\n");
	      return 0;
	    }
	  }
	}
	if(writePSRData(&fout, data, application.verbose_state) == 0) {
	  printerror(application.verbose_state.debug, "ERROR ppolFit: Unable to write data.\n");
	  return 0;
	}
	closePSRData(&fout, 0, application.verbose_state);

	free(data);
      }
//START REGION RELEASE

    
    /* Find point with lowest chi^2 */
    for(i = 0; i < nalpha; i++) {
      for(j = 0; j < nbeta; j++) {
	if(i == 0 && j == 0) {
	  chimax = chigrid[nalpha*j+i];
	  chimin = chigrid[nalpha*j+i];
	}
	if(chigrid[nalpha*j+i] > chimax)
	  chimax = chigrid[nalpha*j+i];
	if(chigrid[nalpha*j+i] < chimin) {
	  chimin = chigrid[nalpha*j+i];
	  bestalpha = i*(alphaend-alphastart)/(double)(nalpha-1)+alphastart; 
	  bestbeta = j*(betaend-betastart)/(double)(nbeta-1)+betastart; 
	}
      }
    }

//START REGION DEVELOP
/* chi2-prob conversion */
    double prob, sum_prob, best_prob;
    if(calc_prob ==1) {
      for(i = 0; i < nalpha; i++) {
	sum_prob = 0;
	for(j = 0; j < nbeta; j++) {
	  prob = exp(-0.5*(chigrid[nalpha*j+i] -chimin));
	  sum_prob += prob;
	}
	prob_a[i] = sum_prob;
      }
      for(j = 0; j < nbeta; j++) {
	sum_prob = 0;
	for(i = 0; i < nalpha; i++) {
	  prob = exp(-0.5*(chigrid[nbeta*j+i] - chimin));
	  sum_prob += prob;
	}
	prob_b[j] = sum_prob;
      }
    }else     if(calc_prob == 2) {
      for(i = 0; i < nalpha; i++) {
	best_prob = -99.99;
	for(j = 0; j < nbeta; j++) {
	  prob = exp(-0.5*(chigrid[nalpha*j+i] -chimin));
	  if(best_prob < prob) best_prob = prob;
	}
	prob_a[i] = best_prob;
      }
      for(j = 0; j < nbeta; j++) {
	best_prob = -99.99;
	for(i = 0; i < nalpha; i++) {
	  prob = exp(-0.5*(chigrid[nbeta*j+i] - chimin));
	  if(best_prob < prob) best_prob = prob;
	}
	prob_b[j] = best_prob;
      }
    }
//START REGION RELEASE

/*********************/

    for(i = 0; i < nalpha; i++) {
      for(j = 0; j < nbeta; j++) {
	/* Normalize */
	chigrid[nalpha*j+i] = (chigrid[nalpha*j+i]-chimin)/(chimax-chimin);
	/* Invert */
	chigrid[nalpha*j+i] = 1-chigrid[nalpha*j+i];
      }
    }
    fprintf(stderr, "Found chi2 values between %e and %e\n", chimin, chimax);
    fprintf(stderr, "Best fit in grid: alpha = %f and beta = %f\n", bestalpha, bestbeta);
    alpha0 = bestalpha;
    beta0 = bestbeta;

    /* Set TR matrix for contour plotting */
    TR[0] = alphastart-0.5*(alphaend-alphastart)/(double)nalpha;
    TR[1] = (alphaend-alphastart)/(double)nalpha;
    TR[2] = 0;
    TR[3] = betastart-0.5*(betaend-betastart)/(double)nbeta;
    TR[4] = 0;
    TR[5] = (betaend-betastart)/(double)nbeta;

//START REGION DEVELOP
    if(	use_sigma_level == 1){
      leftPulseLongitude = fitterinfo.data_l[0]-10;
      rightPulseLongitude = fitterinfo.data_l[fitterinfo.NrDataPoints-1]+10;
    }else {
//START REGION RELEASE
      leftPulseLongitude = 0;
      rightPulseLongitude = 360;
//START REGION DEVELOP
    }
//START REGION RELEASE

    if(showgraphics) {
      if(onlyshowbest == 0) {
	PADeviceID = ppgopen(device2);
	//	ppgpap(8,0.5);
	if(devicenores == 0)
	  pgplot_setWindowsize(paDevice_resx, paDevice_resy, -1);
	PlotPAswing(alpha0, beta0, 0, 0, 0, leftPulseLongitude, rightPulseLongitude, fit_dh0);
	GridDeviceID = ppgopen(device1);
	if(devicenores == 0)
	  pgplot_setWindowsize(gridDevice_resx, gridDevice_resy, -1);
	//ppgpap(8,1);
      }
      level = suppress_greyscale*chimin;
      level = (chimax - level)/(chimax - chimin);
      PrintHelp();
    }

    if(macrofilename) {
      macrofile = fopen(argv[macrofilename], "r");
      if(macrofile == NULL) {
	printerror(application.verbose_state.debug, "Cannot open %s", argv[macrofilename]); 
	return 0;
      }
    }


    redraw = 1;
    do {

      if(showgraphics) {
	if(redraw) {
	  if(onlyshowbest == 2) {
	    GridDeviceID = ppgopen(device1);
	    if(devicenores == 0)
	      pgplot_setWindowsize(gridDevice_resx, gridDevice_resy, -1);
	    //	    ppgpap(8,1);
	  }

	  if(redraw == 2) {
	    if(onlyshowbest != 1) {
	      cpgslct(GridDeviceID);
	      PlotGrid(chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, level, suppress_fac, GridDeviceID, bestalpha, bestbeta, boxlw, labelcharheight, boxlabelcharheight, drawCross, title_txt, drawcontours, nogray, 
//START REGION DEVELOP
		       usezeta, useslope, use_sigma_level, 
//START REGION RELEASE
		       chimax, chimin, PPGPLOT_GRAYSCALE, showwedge, showwedge_label);
	    }
	  }else {
	    if(onlyshowbest != 1) {
	      cpgslct(GridDeviceID);
	      PlotGrid(chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, level, suppress_fac, GridDeviceID, -1, -1, boxlw, labelcharheight, boxlabelcharheight, drawCross, title_txt, drawcontours, nogray, 
//START REGION DEVELOP
		       usezeta, useslope, use_sigma_level, 
//START REGION RELEASE
		       chimax, chimin, PPGPLOT_GRAYSCALE, showwedge, showwedge_label);
	    }
	  }
//START REGION DEVELOP
	  if(commandline_PgplotInterprator(argc, argv, application.verbose_state.debug) == 0)
	    return 0;
//START REGION RELEASE
	  if(calculate_beam_widths != 0 && contour_plot) {
	    if(onlyshowbest != 1) {
	      cpgslct(GridDeviceID);
	      int dotted;
	      if(beamwidth_params_fin != NULL && beamwidth_params_only_w == 0) {  // If beamwidth_params_only_w is set, don't draw contour collection to current device
		rewind(beamwidth_params_fin);
		dotted = 0;
		int ret;
		do {
		  if(calculate_interpulse_widths == 0)
		    ret = fscanf(beamwidth_params_fin, "%lf %f", &fitterinfo.pulse_width, &userContours[0]);
		  else
		    ret = fscanf(beamwidth_params_fin, "%lf %f", &pulse_width2, &userContours[0]);
		  if(ret == 2) {
		    nruserContours = 1;
		    calculate_beam_widths = 1;
		    if(calculate_interpulse_widths == 0)
		      printf("Calculating contours with W=%f deg and rho=%f deg\n", fitterinfo.pulse_width, userContours[0]);
		    else
		      printf("Calculating contours for interpulse with W=%f deg and rho=%f deg\n", pulse_width2, userContours[0]);
		    calcBeamWidths(nalpha, nbeta, alphastart, alphaend, betastart, betaend, 
//START REGION DEVELOP
				   usezeta, useslope, elliptical_beam, 
//START REGION RELEASE
				   calculate_beam_widths, calculate_interpulse_widths, pulse_width2, rhogrid, rhogrid2, application.verbose_state.nocounters);
		    PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, GridDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, dotted);
		    calcIntersectionRhoAndBanana(rhogrid, chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nruserContours, userContours, chimax, chimin, beamwidth_params_vebose);
		  }
		}while(ret == 2);
	      }else {
		if(onlyshowbest != 1) {
		  dotted = 1;
		  PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, GridDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, dotted);
		}
	      }
	    }
	  }
	  if(calculate_interpulse_widths != 0 && contour_plot) {
	    if(onlyshowbest != 1) {
	      cpgslct(GridDeviceID);
	      int dotted;
	      dotted = 1;
	      if(beamwidth_params_fin != NULL && beamwidth_params_only_w == 0) {
		rewind(beamwidth_params_fin);
		dotted = 0;
		int ret;
		do {
		  if(calculate_interpulse_widths == 0)
		    ret = fscanf(beamwidth_params_fin, "%lf %f", &fitterinfo.pulse_width, &userContours[0]);
		  else
		    ret = fscanf(beamwidth_params_fin, "%lf %f", &pulse_width2, &userContours[0]);
		  if(ret == 2) {
		    nruserContours = 1;
		    calculate_interpulse_widths = 1;
		    if(calculate_interpulse_widths == 0)
		      printf("Calculating contours with W=%f deg and rho=%f deg\n", fitterinfo.pulse_width, userContours[0]);
		    else
		      printf("Calculating contours for interpulse with W=%f deg and rho=%f deg\n", pulse_width2, userContours[0]);
		    calcBeamWidths(nalpha, nbeta, alphastart, alphaend, betastart, betaend, 
//START REGION DEVELOP
				   usezeta, useslope, elliptical_beam, 
//START REGION RELEASE
				   calculate_beam_widths, calculate_interpulse_widths, pulse_width2, rhogrid, rhogrid2, application.verbose_state.nocounters);
		    PlotContours(rhogrid2, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels2, TR, GridDeviceID, contour_txt, contourcolorIP, fixedContours, nruserContours, userContours, boxlw, dotted);
		    calcIntersectionRhoAndBanana(rhogrid2, chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nruserContours, userContours, chimax, chimin, beamwidth_params_vebose);
		  }
		}while(ret == 2);
	      }else {
		dotted = 1;
		PlotContours(rhogrid2, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels2, TR, GridDeviceID, contour_txt, contourcolorIP, fixedContours, nruserContours, userContours, boxlw, dotted);
	      }
	    }
	  }
	  redraw = 0;
	  if(onlyshowbest == 2) {
	    ppgclos();
	  }
	}
	if(onlyshowbest == 0)
	   cpgslct(GridDeviceID);
      }
//START REGION DEVELOP
      if(  (calc_prob == 1) || (calc_prob == 2) ){
	FILE *fcalc;
	fcalc = fopen("prob.dat","w");
	for(i = 0; i < nalpha; i++) {
	  fprintf(fcalc,"%f %f %f %f\n", (alphaend-alphastart)*i/((double)(nalpha-1)) + alphastart, prob_a[i] ,(betaend-betastart)*i/((double)(nbeta-1)) + betastart, prob_b[i]  );
	}
	fclose(fcalc);
	calc_prob = 0;
      }else
//START REGION RELEASE
	if(onlyshowbest == 1) {
	  c = 98;
	  onlyshowbest = 2;
      }else if(onlyshowbest == 2) {
	c = 27;
      }else {
	if(showgraphics) {
	  if(macrofilename == 0) {
	    float dummyf1, dummyf2;
	    cpgband(7, 0, 0, 0, &dummyf1, &dummyf2, &c);
	    alpha0 = dummyf1;
	    beta0 = dummyf2;
	  }else {
	    int ret;
	    do {
	      ret = fscanf(macrofile, "%c", &c);
	      if(ret != 1) {
		printf("Reached EOF of macro\n");
		c = 27;
	      }
	    }while(c == 10);   /* Ignore returns */
	  }
	}else {
	  c = 27;
	}
      }
      
      switch(c) {
      case 27: 
      case 'q': break;             /* ESC */
      case 108:                   /* l */
	printf("Set new grayscale level to this number of sigma's (it is now %f): ", (chimax - level*(chimax - chimin))/chimin - 1);
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &level);
	}else {
	  fscanf(macrofile, "%lf", &level);
	  printf("%lf\n", level);
	}
	level = (level+1)*chimin;
	level = (chimax - level)/(chimax - chimin);
	printf("chi2 value corresponding to black = %f and white = %f\n", chimax - level*(chimax - chimin), chimin);
	redraw = 1;
	break;
//START REGION DEVELOP
      case 76:                   /* L */
	printf("Set new level to this fraction (it is now %f): ", level);
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &level);
	}else {
	  fscanf(macrofile, "%lf", &level);
	  printf("%lf\n", level);
	}
	printf("chi2 value corresponding to black = %f and white = %f\n", chimax - level*(chimax - chimin), chimin);
	redraw = 1;
	break;
//START REGION RELEASE
      case 65:                    /* Left mouse button */
	{
	  convertAlphaBeta(&alpha0, &beta0
//START REGION DEVELOP
			   , useslope, usezeta
//START REGION RELEASE
			   );
	  i = (nalpha-1)*(alpha0 - alphastart)/(double)(alphaend-alphastart);
	  j = (nbeta-1)*(beta0 - betastart)/(double)(betaend-betastart);
	  double optimum_l0, optimum_pa0;
	  optimum_l0 = l0grid[nalpha*j+i];
	  optimum_pa0 = pa0grid[nalpha*j+i];
	  if(application.verbose_state.verbose) {
	    printf("Clicked on grid position: %d %d\n", i, j);
	    printf("  alpha: %lf deg\n", alpha0);
	    printf("  beta:  %lf deg\n", beta0);
	    printf("  l0:    %lf deg\n", optimum_l0);
	    printf("  pa0:   %lf deg\n", optimum_pa0);
	  }

	  cpgslct(PADeviceID);
	  DoFitting(alpha0, beta0, pa0, dpa0, l0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 1, stdout, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
	  if(doBruteForce || loadresults) {
	    PlotPAswing(alpha0, beta0, optimum_pa0, optimum_l0, 1, leftPulseLongitude, rightPulseLongitude, fit_dh0);
	    if(doBruteForce) {
	      printwarning(application.verbose_state.debug, "WARNING: Showing solution as found during the grid-search. The fit done after clicking might not be accurate as it did not considered the range of initial parameters.  Run ppolFit with -v to see the parameters used in the plotted RVM curve.");
	    }else {
	      printwarning(application.verbose_state.debug, "WARNING: Showing solution as found during the grid-search when using the -load option. The fit done after clicking might not be accurate as it did not considered the range of initial parameters if the -brute option was used to generate the grid. Run ppolFit with -v to see the parameters used in the plotted RVM curve.");
	    }
	  }else {
	    PlotPAswing(alpha0, beta0, fit_pa0, fit_l0, 1, leftPulseLongitude, rightPulseLongitude, fit_dh0);
	  }
	}
	break;
      case 122:                    /* z */
	printf("Set left edge pulse longitude window: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &leftPulseLongitude);
	}else {
	  fscanf(macrofile, "%lf", &leftPulseLongitude);
	  printf("%lf\n", leftPulseLongitude);
	}
	printf("Set right edge pulse longitude window: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &rightPulseLongitude);
	}else {
	  fscanf(macrofile, "%lf", &rightPulseLongitude);
	  printf("%lf\n", rightPulseLongitude);
	}
	cpgslct(PADeviceID);
	PlotPAswing(alpha0, beta0, 0, 0, 0, leftPulseLongitude, rightPulseLongitude, fit_dh0);
	break;
      case 98:                     /* b */
	printf("Plotting best solution\n");
	if(onlyshowbest == 0)
	  cpgslct(PADeviceID);
	alpha0 = bestalpha;
	beta0 = bestbeta;
	convertAlphaBeta(&alpha0, &beta0
//START REGION DEVELOP
			 , useslope, usezeta
//START REGION RELEASE
			 );
	DoFitting(alpha0, beta0, pa0, dpa0, l0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 1, stdout, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
	DoFitting(alpha0, beta0, fit_pa0, 10, fit_l0, 10, dh0, ddh0, ftol*0.01, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 1, 1, stdout, enableerrors, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state); 
	alpha0 = fit_alpha;
	beta0 = fit_beta;
	printf("     Could do -paswing '%lf %lf %lf %lf' ", alpha0, beta0, fit_pa0, fit_l0);
	int jumpnr;
	if(fitterinfo.nrJumps > 0) {
	  for(jumpnr = 0; jumpnr < fitterinfo.nrJumps; jumpnr++) {
	    printf("-opm '%f %f' ", fitterinfo.jump_longitude[jumpnr], fitterinfo.jump_offset[jumpnr]);
	  }
	}
	printf("in ppolFig\n");
	print_steepness(alpha0, beta0, fit_l0, fit_pa0, application.verbose_state.verbose, fit_dh0, &sina_b);
//START REGION DEVELOP
	/*  output fitting file "fps.out" */
	//	FILE *fpsout;
	//	fpsout = fopen("fps.out","w");
	//	fprintf(fpsout,"%s %f %f %f %f  %f %f\n",argv[argc-1], alpha0, beta0, fit_pa0, fit_l0, sina_b,  chi/(double)(fitterinfo.NrDataPoints-(nrfitparams)));
	//	fclose(fpsout);
//START REGION RELEASE

	if(showgraphics) {
	  if(onlyshowbest != 0) {
	    PADeviceID = ppgopen(device2);
	    if(devicenores == 0)
	      pgplot_setWindowsize(paDevice_resx, paDevice_resy, -1);
	    //	    ppgpap(8,0.5);
	  }
	  PlotPAswing(alpha0, beta0, fit_pa0, fit_l0, 1, leftPulseLongitude, rightPulseLongitude, fit_dh0);
	  if(onlyshowbest != 0) {
	    ppgclos();
	  }
	}
	redraw = 2;
	break;	
//START REGION DEVELOP
      case 112:                    /* p */
	printf("Plot a RVM model with given parameters:\n");
	printf("alpha: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &alpha0);
	}else {
	  fscanf(macrofile, "%lf", &alpha0);
	  printf("%f\n", alpha0);
	}
	printf("beta: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &beta0);
	}else {
	  fscanf(macrofile, "%lf", &beta0);
	  printf("%f\n", beta0);
	}
	printf("l0: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &fit_l0);
	}else {
	  fscanf(macrofile, "%lf", &fit_l0);
	  printf("%f\n", fit_l0);
	}
	printf("pa0: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &fit_pa0);
	}else {
	  fscanf(macrofile, "%lf", &fit_pa0);
	  printf("%f\n", fit_pa0);
	}

	cpgslct(PADeviceID);
	PlotPAswing(alpha0, beta0, fit_pa0, fit_l0, 1, leftPulseLongitude, rightPulseLongitude, fit_dh0);
	break;
//START REGION RELEASE
      case 116:                    /* t */
	if(contour_txt == 0) {
	  printf("Drawing of labels for countours is now switched on\n");
	  contour_txt = 1;
	}else {
	  printf("Drawing of labels for countours is now switched off\n");
	  contour_txt = 0;
	}
	redraw = 1;
	break;	
      case 101:                    /* e */
	if(enableerrors != 1) {
	  printf("\nThe error estimation is done by starting at the best solution. Then one by one each fit parameter is changed until the chi2 becomes (nr of sigma+1) times bigger. In this process a downhill-simplex search is done for all parameters to take into account the covariance between the parameters.\n\n");
	  printf("Error estimation on best fit is now switched on\n");
	  enableerrors = 1;
	}else {
	  printf("Error estimation on best fit is now switched off\n");
	  enableerrors = 0;
	}
	redraw = 1;
	break;	
//START REGION DEVELOP
      case 69:                    /* E */
	if(enableerrors != 2) {
	  printf("\nThe error estimation is done by starting at the best solution. Then one by one each fit parameter is changed until the chi2 becomes (nr of sigma+1) times bigger. In this process a downhill-simplex search is done for all parameters (except PA_0) to take into account the covariance between the parameters.\n\n");
	  printf("Error estimation on best fit is now switched on (while keeping PA_0 fixed)\n");
	  enableerrors = 2;
	}else {
	  printf("Error estimation on best fit is now switched off (while keeping PA_0 fixed)\n");
	  enableerrors = 0;
	}
	redraw = 1;
	break;	
//START REGION RELEASE
      case 84:                    /* T */
	if(title_txt == 0) {
	  printf("Drawing title is switched on\n");
	  title_txt = 1;
	}else {
	  printf("Drawing title is switched off\n");
	  title_txt = 0;
	}
	redraw = 1;
	break;	
      case 103:                    /* g */
	if(nogray == 0) {
	  printf("Drawing grayscale plot switched off\n");
	  nogray = 1;
	}else {
	  printf("Drawing grayscale plot switched on\n");
	  nogray = 0;
	}
	redraw = 1;
	break;	
      case 115:                    /* s */
	printf("Suppress factor for grayscale (>= 1, currently is set to %lf): ", suppress_fac);
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &suppress_fac);
	}else {
	  fscanf(macrofile, "%lf", &suppress_fac);
	  printf("%f\n", suppress_fac);
	}
	redraw = 1;
	break;		
      case 83:                    /* S */
	printf("Set the number of sigmas in error calculation to (it is now %.1f): ", nrofsigmas);
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &nrofsigmas);
	}else {
	  fscanf(macrofile, "%lf", &nrofsigmas);
	  printf("%f\n", nrofsigmas);
	}
	redraw = 1;
	break;
      case 71:                    /* G */
	{
	double minchi2, maxchi2, chi2, minbeta, maxbeta, beta, minalpha, maxalpha, alpha, minl0, maxl0, minpa0, maxpa0, mindh, maxdh, minslope, maxslope;
	double minl0_2, maxl0_2, minpa0_2, maxpa0_2, mindh_2, maxdh_2, optimum_l0, optimum_pa0, optimum_chi2;
	int chi2unset, betaunset, fullsearch, reducel0, reducepa0, method;
	chi2unset = 1;
	betaunset = 1;
	fullsearch = 0;
	optimum_l0 = 0;
	do {
	  printf("What algorithm do you want to use: \n");
	  printf("  1. Look for range parameters by exploring best solutions for each alpha/beta pair (quick).\n");
	  printf("  2. For each alpha/beta pair, explore how far l0 and pa0 can be pushed while solving for the other parameter (slow, but more precise).\n");
	  printf("  3. Similar to option 2, but using downhill-simplex rather than Brent method.\n\n");
	  printf("Type number in terminal: ");
	  fflush(stdout);

	  if(macrofilename == 0) {
	    scanf("%d", &method);
	  }else {
	    fscanf(macrofile, "%d", &method);
	  }
	}while(method < 1 || method > 3);
	if(method == 1) {
	  fullsearch = 0;
	  printf("\nDoing a quick search\n");
	}else if(method == 2) {
	  fullsearch = 1;
	  usegsl_error_grid_search = 1;
	  printf("\nDoing a slow search with Brent method\n");
	}else if(method == 3) {
	  fullsearch = 1;
	  usegsl_error_grid_search = 0;
	  printf("\nDoing a slow search with downhill-simplex method\n");
	}

	/*
	printf("Do you want to make range in l0 smaller than 180 deg  (y/n): \n");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%s", txt);
	}else {
	  fscanf(macrofile, "%s", txt);
	}
	if(txt[0] == 'y' || txt[0] == 'Y')
	  reducel0 = 1;
	if(reducel0)
	  printf("Reducing l0 angle\n");
	else
	  printf("Not reducing l0 angle\n");
	fprintf(stderr, "\n");

	printf("\nWARNING: PA0 range is not adjusted for 180 deg jumps of l0, so PA0 range is probably overestimated\n\n");
	if(reducel0) {
	  printf("Prevent l0 to deviate more than 90 deg from this angle: \n");
	  fflush(stdout);
	  if(macrofilename == 0) {
	    scanf("%lf", &optimum_l0);
	  }else {
	    fscanf(macrofile, "%lf", &optimum_l0);
	  }
	  printf("Keeping l0 close to %lf degrees\n", optimum_l0);
	  fprintf(stderr, "\n");
	}
	*/

	reducel0 = 1;
	reducepa0 = 1;
	for(i = 0; i < nalpha; i++) {
	  for(j = 0; j < nbeta; j++) {
	    chi2 = 1.0 - chigrid[nalpha*j+i];
	    chi2 = chi2*(chimax - chimin) + chimin;
	    if((i == 0 && j == 0) || chi2 < optimum_chi2) {
	      optimum_chi2 = chi2;
	      optimum_l0 = l0grid[nalpha*j+i];
	      optimum_pa0 = pa0grid[nalpha*j+i];
	      //	      alpha = i*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
	      //	      beta = j*(betaend-betastart)/(double)(nbeta-1)+betastart;
	      //	      convertAlphaBeta(&alpha, &beta, useslope, usezeta);
	      //	      printf("a=%f b=%f chi=%f\n", alpha, beta, chi2);
	    }
	  }
	}
	printf("Found angles for l0 will be de-wrapped to an angle close to %lf degrees\n", optimum_l0);
	printf("Found angles for pa0 will be de-wrapped to an angle close to %lf degrees for chi2=%f\n", optimum_pa0, optimum_chi2);
	
	
	/*    int ret;
    double xfit[5], pa0, errorbar;
    xfit[0] = -166.844391;
    xfit[1] = 174.098953;
    fitterinfo.l0_start = xfit[1];
    xfit[2] = 3.053257;
    xfit[2] = 90.0;
    xfit[3] = 0.187013;
    xfit[3] = 3.11;
    xfit[4] = 0;
    //double internal_funk_gsl(double *x, void *params)
    //    ret = minimize_1D_double(1, internal_PA0_funk, xfit, 0, 180, 10, &pa0, 100, 0.0001, 0.0001, 1);
    ret = find_1D_error(internal_funk_gsl, xfit, 2, 5, 10, NULL, 3.0, 100, 0.0001, 0.0001, &errorbar, 1);
    if(ret != 0) {
      if(ret == 1) {
	printwarning(application.verbose_state.debug, "WARNING: Maximum nr of itterations exceeded\n");
      }else if(ret == 2) {
	printwarning(application.verbose_state.debug, "WARNING: Did not found root\n");
      }else if(ret == 3) {
	printwarning(application.verbose_state.debug, "WARNING: Lower and upper limit do not bracket a root\n");
      }else {
	printerror(application.verbose_state.debug, "ERROR: Unknown error in minimize_1D_double\n");
	exit(0);
      }
    }
    exit(0);
	*/
	for(i = 0; i < nalpha; i++) {
	  for(j = 0; j < nbeta; j++) {
	    chi2 = 1.0 - chigrid[nalpha*j+i];
	    chi2 = chi2*(chimax - chimin) + chimin;
	    if(chi2 <= (nrofsigmas+1)*chimin) {
	      if(chi2unset) {
		maxchi2 = minchi2 = chi2;
		chi2unset = 0;
	      }
	      if(chi2 > maxchi2)
		maxchi2 = chi2;
	      if(chi2 < minchi2)
		minchi2 = chi2;
	      alpha = i*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
	      beta = j*(betaend-betastart)/(double)(nbeta-1)+betastart;
	      convertAlphaBeta(&alpha, &beta
//START REGION DEVELOP
			       , useslope, usezeta
//START REGION RELEASE
			       );
	      double tmp_l0, tmp_pa0, tmp_slope;
	      tmp_l0 = l0grid[nalpha*j+i];
	      if(reducel0) {
		tmp_l0 = derotate_180_small_double(tmp_l0 - optimum_l0);
		//		if(tmp_l0 > 80) 
		//		  printf("XXXXXX: alpha=%f beta=%f\n", alpha, beta);
	      }
	      tmp_pa0 = pa0grid[nalpha*j+i];
	      if(reducepa0)
		tmp_pa0 = derotate_180_small_double(tmp_pa0 - optimum_pa0);

	      tmp_slope = sin(alpha*M_PI/180.0)/sin(beta*M_PI/180.0);
	      if(betaunset) {
		minbeta = maxbeta = beta;
		minalpha = maxalpha = alpha;
		minl0 = maxl0 = tmp_l0;
		minpa0 = maxpa0 = tmp_pa0;
		mindh = maxdh = dhgrid[nalpha*j+i];
		minl0_2 = maxl0_2 = tmp_l0;
		minpa0_2 = maxpa0_2 = tmp_pa0;
		mindh_2 = maxdh_2 = dhgrid[nalpha*j+i];
		minslope = maxslope = tmp_slope;
		betaunset = 0;
	      }
	      /*
	      if(alpha > 95 && alpha < 105 && beta < -7 && beta > -8) {
		fprintf(stderr, "XXXXXX %f\n", l0grid[nalpha*j+i]);
		}*/
	      if(beta < minbeta)
		minbeta = beta;
	      if(beta > maxbeta)
		maxbeta = beta;
	      if(alpha < minalpha)
		minalpha = alpha;
	      if(alpha > maxalpha)
		maxalpha = alpha;
	      if(tmp_l0 < minl0)
		minl0 = tmp_l0;    
	      if(tmp_l0 > maxl0)
		maxl0 = tmp_l0;    
	      if(tmp_pa0 < minpa0)
		minpa0 = tmp_pa0;    
	      if(tmp_pa0 > maxpa0)
		maxpa0 = tmp_pa0;    
	      if(dhgrid[nalpha*j+i] < mindh)
		mindh = dhgrid[nalpha*j+i];    
	      if(dhgrid[nalpha*j+i] > maxdh)
		maxdh = dhgrid[nalpha*j+i];
	      if(tmp_slope > maxslope)
		maxslope = tmp_slope;
	      if(tmp_slope < minslope)
		minslope = tmp_slope;

	      if(fullsearch) {
		double dx[5], xfit[5], dplus, dmin;
		int fixed[5], param_i, param_j;
		for(param_i = 0; param_i < 3; param_i++) {
		  for(param_j = 0; param_j < 5; param_j++) {
		    dx[param_j] = 1;
		    fixed[param_j] = 1;
		  }
		  xfit[0] = pa0grid[nalpha*j+i];
		  xfit[1] = l0grid[nalpha*j+i];
		  //		  fflush(stdout); fprintf(stderr, "Full search for l0=%f\n", fitterinfo.l0_start);
		  fitterinfo.l0_start = xfit[1];
		  xfit[2] = alpha;
		  xfit[3] = beta;
		  xfit[4] = dhgrid[nalpha*j+i];
		  fixed[0] = 0;   /* Allow other parameter than alpha and beta (which are fixed by grid search) to vary */
		  fixed[1] = 0;
		  fixed[4] = 0;
		  nrfitparams = 3;   /* In grid search, the number of fit parameters is 3, unless dh is not fitted for */
		  if(dh0 < -2 ) {
		    /* If not fitting for emission height difference,
		       fit for alpha instead, but with zero step size
		       so it is not really a search. There have to be
		       2 free fit parameters (and the parameter of
		       which we will determine the error bar is not a
		       free parameter), otherwise the amoeba algorithm
		       will not work. Note, all this will only be done
		       when using amoeba algorithm.*/
		    fixed[4] = 1;
		    fixed[2] = 0;
		    dx[2] = 0;
		    nrfitparams -= 1;
		    /* Don't try to find dh errors if not fitted for */
		    if(param_i == 2)
		      break;
		  }else {
		    if(usegsl_error_grid_search) {
		      printerror(application.verbose_state.debug, "Searching for errorbars when allowing emission height difference will not work without using downhill-simplex, as it is no longer a 1D fitting process.");
		      exit(0);
		    }
		  }

		  int debug_verbose;
		  debug_verbose = 0;
		  /*
		  if(debug_verbose) {
		    printf("XXXXX Grid claimed reduced chi2=%f for alpha=%f, beta=%f, pa0=%f, l0=%f and dh=%f\n", chi2, xfit[2], xfit[3], xfit[0], xfit[1], xfit[4]);
		    printf("XXXXX Grid claimed chi2 = %f (i.e. %f)\n", chi2, chi2*(fitterinfo.NrDataPoints-nrfitparams));
		    printf("XXXXX funk claims: chi2 = %f\n", funk(xfit));
		  }
		  */
		  
		  /* Calculate total chi2 limit: */
		  double chi2_notreduced;
		  chi2_notreduced = chimin*(fitterinfo.NrDataPoints-nrfitparams);

		  if(param_i == 0) {
		    if(usegsl_error_grid_search == 0) {
		      if(find_errors_amoeba_d(amoeba_algorithm, dx, fixed, xfit, chi2_notreduced, 5, funk, ftol, 0, &dplus, &dmin, nrofsigmas) != 0)
			return 0;
		    }else {
		      //		      fflush(stdout);
		      //		      fprintf(stderr, "DO A TRIAL: PA0 +\n");
		      internal_fit_pa_or_l0 = 1;
		      ret = find_1D_error(internal_funk_gsl, xfit, 0, 5, 1, -1, NULL,  nrofsigmas, chi2_notreduced, 2000, 0.01, 0.0, &dplus, debug_verbose);
		      if(ret != 0) {
			fflush(stdout);
			printwarning(application.verbose_state.debug, "WARNING: Finding error failed for this trial");
			if(ret == 1) {
			  printwarning(application.verbose_state.debug, "WARNING: Maximum nr of itterations exceeded");
			}else if(ret == 2) {
			  printwarning(application.verbose_state.debug, "WARNING: Did not found root");
			}else if(ret == 3) {
			  printwarning(application.verbose_state.debug, "WARNING: Lower and upper limit do not bracket a root");
			}else {
			  printwarning(application.verbose_state.debug, "WARNING: Unknown error in minimize_1D_double");
			}			
		      }
		      //		      fflush(stdout);
		      //		      fprintf(stderr, "DO A TRIAL: PA0 -\n");
		      ret = find_1D_error(internal_funk_gsl, xfit, 0, 5, 1, -1, NULL, -nrofsigmas, chi2_notreduced, 2000, 0.01, 0.0, &dmin, debug_verbose);
		      dmin *= -1.0;
		      if(ret != 0) {
			fflush(stdout);	
			printwarning(application.verbose_state.debug, "WARNING: Finding error failed for this trial");
			if(ret == 1) {
			  printwarning(application.verbose_state.debug, "WARNING: Maximum nr of itterations exceeded");
			}else if(ret == 2) {
			  printwarning(application.verbose_state.debug, "WARNING: Did not found root");
			}else if(ret == 3) {
			  printwarning(application.verbose_state.debug, "WARNING: Lower and upper limit do not bracket a root");
			}else {
			  printwarning(application.verbose_state.debug, "WARNING: Unknown error in minimize_1D_double");
			}			
		      }
		    }
		    tmp_pa0 = xfit[0]+dplus;
		    if(reducepa0)
		      tmp_pa0 = derotate_180_small_double(tmp_pa0 - optimum_pa0);
		    if(tmp_pa0 > maxpa0_2)
		      maxpa0_2 = tmp_pa0;
		    tmp_pa0 = xfit[0]+dmin;
		    if(reducepa0)
		      tmp_pa0 = derotate_180_small_double(tmp_pa0 - optimum_pa0);
		    if(tmp_pa0 < minpa0_2)
		      minpa0_2 = tmp_pa0;
		  }else if(param_i == 1) {
		    if(usegsl_error_grid_search == 0) {
		      if(find_errors_amoeba_d(amoeba_algorithm, dx, fixed, xfit, chi2_notreduced, 5, funk, ftol, 1, &dplus, &dmin, nrofsigmas) != 0)
			return 0;
		    }else {
		      //		      fflush(stdout);
		      //		      fprintf(stderr, "DO A TRIAL: l0 +\n");
		      internal_fit_pa_or_l0 = 2;
		      ret = find_1D_error(internal_funk_gsl, xfit, 1, 5, 1, fitterinfo.max_l0_diff, NULL, nrofsigmas, chi2_notreduced, 2000, 0.01, 0.0, &dplus, debug_verbose);
		      if(ret != 0) {
			fflush(stdout);
			printwarning(application.verbose_state.debug, "WARNING: Finding error failed for this trial");
			if(ret == 1) {
			  printwarning(application.verbose_state.debug, "WARNING: Maximum nr of itterations exceeded");
			}else if(ret == 2) {
			  printwarning(application.verbose_state.debug, "WARNING: Did not found root");
			}else if(ret == 3) {
			  printwarning(application.verbose_state.debug, "WARNING: Lower and upper limit do not bracket a root");
			}else {
			  printwarning(application.verbose_state.debug, "WARNING: Unknown error in minimize_1D_double");
			}			
		      }
		      //		      fflush(stdout);
		      //		      fprintf(stderr, "DO A TRIAL: l0 -\n");
		      ret = find_1D_error(internal_funk_gsl, xfit, 1, 5, 1, fitterinfo.max_l0_diff, NULL, -nrofsigmas, chi2_notreduced, 2000, 0.01, 0.0, &dmin, debug_verbose);
		      dmin *= -1.0;
		      if(ret != 0) {
			fflush(stdout);
			printwarning(application.verbose_state.debug, "WARNING: Finding error failed for this trial");
			if(ret == 1) {
			  printwarning(application.verbose_state.debug, "WARNING: Maximum nr of itterations exceeded");
			}else if(ret == 2) {
			  printwarning(application.verbose_state.debug, "WARNING: Did not found root");
			}else if(ret == 3) {
			  printwarning(application.verbose_state.debug, "WARNING: Lower and upper limit do not bracket a root");
			}else {
			  printwarning(application.verbose_state.debug, "WARNING: Unknown error in minimize_1D_double");
			}			
		      }
		    }
		    tmp_l0 = xfit[1]+dplus;
		    if(reducel0)
		      tmp_l0 = derotate_180_small_double(tmp_l0 - optimum_l0);
		    if(tmp_l0 > maxl0_2)
		      maxl0_2 = tmp_l0;
		    tmp_l0 = xfit[1]+dmin;
		    if(reducel0)
		      tmp_l0 = derotate_180_small_double(tmp_l0 - optimum_l0);
		    if(tmp_l0 < minl0_2)
		      minl0_2 = tmp_l0;
		  }else if(param_i == 2) {
		    if(find_errors_amoeba_d(amoeba_algorithm, dx, fixed, xfit, chi2_notreduced, 5, funk, ftol, 4, &dplus, &dmin, nrofsigmas) != 0)
		      return 0;
		    if(xfit[4]+dplus > maxdh_2)
		      maxdh_2 = xfit[4]+dplus;
		    if(xfit[4]+dmin < mindh_2)
		      mindh_2 = xfit[4]+dmin;
		
		  }
		}
	      }  // End of if(fullsearch)
	    }
	    if(application.verbose_state.nocounters == 0)
	      fprintf(stderr, "%.1f%%     \r",(100.0*(i+1))/(double)(nalpha));
	  }
	}

	fprintf(stderr, "\n");
	printf("After considering all the best solutions in the alpha-beta grid the following %.1f sigma ranges were found:\n", nrofsigmas);
	printf("The range in chi2                 = %f to %f", minchi2, maxchi2);
//START REGION DEVELOP
	printf(" (minimum should be %f, maximum ~ %f)", chimin, chimin*(1+nrofsigmas));
//START REGION RELEASE
	printf("\n");
	printf("The range in alpha                = %f to %f deg\n", minalpha, maxalpha);
	printf("The range in beta                 = %f to %f deg\n", minbeta, maxbeta);
	printf("The range in l0                   = %f to %f deg", minl0, maxl0);
	if(reducel0) {
	  printf(" offset with respect to %lf\n", optimum_l0);
	}else {
	  printf("\n");
	}
	printf("The range in pa0                  = %f to %f", minpa0, maxpa0);
	if(reducepa0) {
	  printf(" offset with respect to %lf\n", optimum_pa0);
	}else {
	  printf("\n");
	}
	printf("The range in sin(alpha)/sin(beta) = %f to %f\n", minslope, maxslope);
	printf("The range in dh                   = %f to %f (fraction of light cylinder radius)\n", mindh, maxdh);
	if(fullsearch) {
	  printf("\nEach of those solutions were taken as a start point of a calculation similar to the 'e' option. For each solution the allowed ranges for l0 (and pa0 and dh if used) were explored while fitting for pa0 and dh if used (or the other permutations in case the range in pa0 or dh are solved for) the following ranges were found:\n");
	  printf("The range in l0    = %f to %f deg", minl0_2, maxl0_2);
	  if(reducel0) {
	    printf(" offset with respect to %lf\n", optimum_l0);
	  }else {
	    printf("\n");
	  }
	  printf("The range in pa0   = %f to %f deg", minpa0_2, maxpa0_2);
	  if(reducepa0) {
	    printf(" offset with respect to %lf\n", optimum_pa0);
	  }else {
	    printf("\n");
	  }
	  printf("The range in dh    = %f to %f deg\n", mindh_2, maxdh_2);
	  redraw = 1;
	}
	}
	break;		
      case 99:                    /* c */
	printf("Nr contour levels: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%d", &drawcontours);
	}else {
	  fscanf(macrofile, "%d", &drawcontours);
	  printf("%d\n", drawcontours);
	}
	redraw = 1;
	break;		
      case 111:                    /* o */
	printf("Change output plot parameters\n");
	printf("line width of box (now set to %d): ", boxlw);
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%d", &boxlw);
	}else {
	  fscanf(macrofile, "%d", &boxlw);
	  printf("%d\n", boxlw);
	}
	printf("character height labels (now set to %f): ", labelcharheight);
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &labelcharheight);
	}else {
	  fscanf(macrofile, "%lf", &labelcharheight);
	  printf("%f\n", labelcharheight);
	}
	printf("label numbers character height (now set to %f): ", boxlabelcharheight);
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &boxlabelcharheight);
	}else {
	  fscanf(macrofile, "%lf", &boxlabelcharheight);
	  printf("%f\n", boxlabelcharheight);
	}
	printf("Colour index for contours (now set to %d): ", contourcolor);
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%d", &contourcolor);
	}else {
	  fscanf(macrofile, "%d", &contourcolor);
	  printf("%d\n", contourcolor);
	}
	printf("Colour index for contours of interpulse (now set to %d): ", contourcolorIP);
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%d", &contourcolorIP);
	}else {
	  fscanf(macrofile, "%d", &contourcolorIP);
	  printf("%d\n", contourcolorIP);
	}
	redraw = 1;
	break;		
      case 114:                    /* r */
	if(contour_plot == 0) {
	  printf("Countour plot for beam radius is switched on\n");
	  contour_plot = 1;
	}else {
	  printf("Countour plot for beam radius is switched off\n");
	  contour_plot = 0;
	}
	redraw = 1;
	break;		
      case 110:                  /* N */
	fixedContours = 0;
	nruserContours = 0;
	printf("Nr contour levels: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%d", &nrcontourlevels);
	}else {
	  fscanf(macrofile, "%d", &nrcontourlevels);
	  printf("%d\n", nrcontourlevels);
	}
	printf("nr contour levels interpulse: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%d", &nrcontourlevels2);
	}else {
	  fscanf(macrofile, "%d", &nrcontourlevels2);
	  printf("%d\n", nrcontourlevels2);
	}
	redraw = 1;
	break;			
      case 102:                  /* f */
	if(fixedContours) {
	  fixedContours = 0;
	  printf("Don't use fixed contours.\n");
	}else {
	  fixedContours = 1;
	  printf("Use fixed contours.\n");
	}
	redraw = 1;
	break;			
      case 70:                  /* F */
	printf("Specify contours (e.g. 5,10,20, hit return without specifying anything to disable drawing contours):\n");
	if(macrofilename == 0) {
	  fgets(txt, 990, stdin);
	}else {
	  int index, ret;
	  char character;
	  index = 0;
	  while((ret = fscanf(macrofile, "%c", &character))==1) {
	    txt[index++] = character;
	    if(index == 990 || character == '\n' || character == '\r') {
	      txt[index] = 0;
	      break;
	    }
	  }
	}
	if(strlen(txt) == 0) {
	  printwarning(application.verbose_state.debug, "WARNING: empty string. If not correct, try F again.");
	}
	//	strcat(txt, ",");
	txtptr = strtok(txt, ",");
	nruserContours = 0;
	do {
	  if(txtptr != NULL) {
	    //	    printf("XXXXX %s\n", txtptr);
	    sscanf(txtptr, "%f", &userContours[nruserContours++]);
	    printf("contour %d = %f\n", nruserContours, userContours[nruserContours-1]);
	  }
	  txtptr = strtok(NULL, ",");
	}while(txtptr != NULL);
	redraw = 1;
	break;			
      case 82:                  /* R */
	{
	  char beamwidth_params_file[1000];
	  printf("Specify ascii file with pulse a width and rho contour (could use the -contcol option): ");
	  fflush(stdout);
	  if(macrofilename == 0) {
	    scanf("%s", beamwidth_params_file);
	  }else {
	    int ret;
	    ret = fscanf(macrofile, "%s", beamwidth_params_file);
	    if(ret != 1) {
	      printerror(application.verbose_state.debug, "ERROR ppolFit: failed to read in string from macrofile");
	      return 0;
	    }
	  }
	  printf("Will read parameters from %s\n", beamwidth_params_file);
	  beamwidth_params_fin = fopen(beamwidth_params_file, "r");
	  if(beamwidth_params_fin == NULL) {
	    printerror(application.verbose_state.debug, "ERROR ppolFit: cannot open '%s'", beamwidth_params_file);
	    return 0;
	  }
	}
	printf("Specify MP or IP: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%s", txt);
	}else {
	  ret = fscanf(macrofile, "%s", txt);
	  if(ret != 1) {
	    printerror(application.verbose_state.debug, "ERROR ppolFit: failed to read in string from macrofile");
	    return 0;
	  }
	}
	if(strcasecmp(txt, "IP") == 0) {
	  calculate_interpulse_widths = 1;
	  calculate_beam_widths = 0;
	  printf("Take interpulse\n");
	}else {
	  calculate_interpulse_widths = 0;
	  calculate_beam_widths = 1;
	  printf("Take main pulse\n");
	}
	
	printf("Plot intersections between chi2 surface and contours (y/n): ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%s", txt);
	}else {
	  ret = fscanf(macrofile, "%s", txt);
	  if(ret != 1) {
	    printerror(application.verbose_state.debug, "ERROR ppolFit: cannot read in string from macrofile");
	    return 0;
	  }
	}
	if(strcasecmp(txt, "y") == 0) {
	  beamwidth_params_vebose.verbose = 1;
	  printf("Do show intersections between chi2 surface and contours\n");
	}else {
	  beamwidth_params_vebose.verbose = 0;
	  printf("Do not show intersections between chi2 surface and contours\n");
	}

	printf("Draw contours ONLY with w option, i.e. not to current device (y/n): ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%s", txt);
	}else {
	  ret = fscanf(macrofile, "%s", txt);
	  if(ret != 1) {
	    printerror(application.verbose_state.debug, "ERROR ppolFit: cannot read in string from macrofile");
	    return 0;
	  }
	}
	if(strcasecmp(txt, "y") == 0) {
	  beamwidth_params_only_w = 1;
	  printf("Drawing contours ONLY with w option\n");
	}else {
	  beamwidth_params_only_w = 0;
	  printf("Always drawing contours\n");
	}

	redraw = 1;
	break;
      case 67:                  /* C */
	if(drawCross == 2) {
	  drawCross = 0;
	  printf("Disabled drawing of the best-fit cross\n");
	}else if(drawCross == 1) {
	  drawCross = 2;
	  printf("Enabled drawing of the best-fit cross (Black)\n");
	}else if(drawCross == 0) {
	  drawCross = 1;
	  printf("Enabled drawing of the best-fit cross (Red)\n");
	}
	redraw = 1;
	break;			
      case 119:                /* w */
	invertGrayscale = 1;
	printf("Do you want to specify pgplot device for chi2 grid yourself (yes/no)? ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%s", txt);
	}else {
	  ret = fscanf(macrofile, "%s", txt);
	  if(ret != 1) {
	    printerror(application.verbose_state.debug, "ERROR ppolFit: cannot read in string from macrofile");
	    return 0;
	  }
	}
	if(strcasecmp(txt, "yes") == 0) {
	  printf("Specify pgplot device: ");
	  fflush(stdout);
	  if(macrofilename == 0) {
	    scanf("%s", txt);
	  }else {
	    ret = fscanf(macrofile, "%s", txt);
	    if(ret != 1) {
	      printerror(application.verbose_state.debug, "ERROR ppolFit: cannot read in string from macrofile");
	      return 0;
	    }
	  }	  
	  printf("Writing output to device %s\n", txt);
	  PSDeviceID = ppgopen(txt);
	  if(devicenores == 0)
	    pgplot_setWindowsize(gridDevice_resx, gridDevice_resy, -1);

	  printf("Do you want to invert grayscale (yes/no)? ");
	  fflush(stdout);
	  if(macrofilename == 0) {
	    scanf("%s", txt);
	  }else {
	    ret = fscanf(macrofile, "%s", txt);
	    if(ret != 1) {
	      printerror(application.verbose_state.debug, "ERROR ppolFit: cannot read in string from macrofile");
	      return 0;
	    }
	  }
	  if(strcasecmp(txt, "yes") == 0) {
	    invertGrayscale = 1;
	    printf("Inverting grayscale\n");
	  }else {
	    invertGrayscale = 0;
	    printf("Not inverting grayscale\n");
	  }

	}else {
	  printf("Writing output to %s.ps\n", prefix);
	  sprintf(txt, "%s.ps/cps", prefix);
	  PSDeviceID = ppgopen(txt);
	}
	ppgask(0);
	cpgslct(PSDeviceID);
	if(invertGrayscale) {
	  PlotGrid(chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, level, suppress_fac, PSDeviceID, bestalpha, bestbeta, boxlw, labelcharheight, boxlabelcharheight, drawCross, title_txt, drawcontours, nogray, 
//START REGION DEVELOP
		   usezeta, useslope, use_sigma_level, 
//START REGION RELEASE
		   chimax, chimin, PPGPLOT_INVERTED_GRAYSCALE, showwedge, showwedge_label);
	}else {
	  PlotGrid(chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, level, suppress_fac, PSDeviceID, bestalpha, bestbeta, boxlw, labelcharheight, boxlabelcharheight, drawCross, title_txt, drawcontours, nogray, 
//START REGION DEVELOP
		   usezeta, useslope, use_sigma_level, 
//START REGION RELEASE
		   chimax, chimin, PPGPLOT_GRAYSCALE, showwedge, showwedge_label);
	}

//START REGION DEVELOP
	if(commandline_PgplotInterprator(argc, argv, application.verbose_state.debug) == 0)
	  return 0;
//START REGION RELEASE

	if(calculate_beam_widths != 0 && contour_plot) {
	  int dotted;
	  if(beamwidth_params_fin != NULL) {
	    rewind(beamwidth_params_fin);
	    dotted = 0;
	    int ret;
	    do {
	      if(calculate_interpulse_widths == 0)
		ret = fscanf(beamwidth_params_fin, "%lf %f", &fitterinfo.pulse_width, &userContours[0]);
	      else
		ret = fscanf(beamwidth_params_fin, "%lf %f", &pulse_width2, &userContours[0]);
	      if(ret == 2) {
		nruserContours = 1;
		calculate_beam_widths = 1;
		if(calculate_interpulse_widths == 0)
		  printf("Calculating contours with W=%lf deg and rho=%f deg\n", fitterinfo.pulse_width, userContours[0]);
		else
		  printf("Calculating contours for ip with W=%f deg and rho=%f deg\n", pulse_width2, userContours[0]);
		calcBeamWidths(nalpha, nbeta, alphastart, alphaend, betastart, betaend, 
//START REGION DEVELOP
			       usezeta, useslope, elliptical_beam, 
//START REGION RELEASE
			       calculate_beam_widths, calculate_interpulse_widths, pulse_width2, rhogrid, rhogrid2, application.verbose_state.nocounters);
		PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, PSDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, dotted);
	      }
	    }while(ret == 2);
	  }else {
	    dotted = 1;
	    PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, PSDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, dotted);
	    PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, PSDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, dotted);
	  }
	}
	  /*	if(calculate_beam_widths != 0 && contour_plot) {
	  PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, PSDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, 1);
	  }*/

	/*	if(calculate_interpulse_widths != 0 && contour_plot) {
	  PlotContours(rhogrid2, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels2, TR, PSDeviceID, contour_txt, contourcolorIP, fixedContours, nruserContours, userContours, boxlw, 1);
	  } */

	if(calculate_interpulse_widths != 0 && contour_plot) {
	  int dotted;
	  if(beamwidth_params_fin != NULL) {
	    rewind(beamwidth_params_fin);
	    dotted = 0;
	    int ret;
	    do {
	      if(calculate_interpulse_widths == 0)
		ret = fscanf(beamwidth_params_fin, "%lf %f", &fitterinfo.pulse_width, &userContours[0]);
	      else
		ret = fscanf(beamwidth_params_fin, "%lf %f", &pulse_width2, &userContours[0]);
	      if(ret == 2) {
		nruserContours = 1;
		calculate_interpulse_widths = 1;
		if(calculate_interpulse_widths == 0)
		  printf("Calculating contours with W=%lf deg and rho=%f deg\n", fitterinfo.pulse_width, userContours[0]);
		else
		  printf("Calculating contours for ip with W=%f deg and rho=%f deg\n", pulse_width2, userContours[0]);
		calcBeamWidths(nalpha, nbeta, alphastart, alphaend, betastart, betaend, 
//START REGION DEVELOP
			       usezeta, useslope, elliptical_beam, 
//START REGION RELEASE
			       calculate_beam_widths, calculate_interpulse_widths, pulse_width2, rhogrid, rhogrid2, application.verbose_state.nocounters);
		PlotContours(rhogrid2, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels2, TR, GridDeviceID, contour_txt, contourcolorIP, fixedContours, nruserContours, userContours, boxlw, dotted);
	      }
	    }while(ret == 2);
	  }else {
	    dotted = 1;
	    PlotContours(rhogrid2, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels2, TR, GridDeviceID, contour_txt, contourcolorIP, fixedContours, nruserContours, userContours, boxlw, dotted);
	  }
	}


	ppgclos();
	printf("Writing chi2 grid done. Making a PA-swing plot is next.\n");
	printf("alpha: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &alpha0);
	}else {
	  fscanf(macrofile, "%lf", &alpha0);
	  printf("%f\n", alpha0);
	}
	printf("beta: ");
	fflush(stdout);
	if(macrofilename == 0) {
	  scanf("%lf", &beta0);
	}else {
	  fscanf(macrofile, "%lf", &beta0);
	  printf("%f\n", beta0);
	}
	DoFitting(alpha0, beta0, pa0, dpa0, l0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 1, stdout, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
	printf("Writing output to pa.ps\n");
	ppgopen("pa.ps/cps");
	PlotPAswing(alpha0, beta0, fit_pa0, fit_l0, 1, leftPulseLongitude, rightPulseLongitude, fit_dh0);
	ppgclos();
	cpgslct(GridDeviceID);
	printf("Making plots done.\n");
	break;				
//START REGION DEVELOP
      case 87:                /* W */
	printf("Writing output to %s.pgplot\n", prefix);
	sprintf(txt, "%s.pgplot", prefix);
	if(pgopenoutputfile(txt)) {
	  sprintf(txt, "%s.ps/cps", prefix);
	  PSDeviceID = ppgopen(txt);
	  cpgslct(PSDeviceID);
	  ppgask(0);
	  PlotGrid(chigrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, level, suppress_fac, PSDeviceID, bestalpha, bestbeta, boxlw, labelcharheight, boxlabelcharheight, drawCross, title_txt, drawcontours, nogray, 
		   usezeta, useslope, use_sigma_level, 
		   chimax, chimin, PPGPLOT_INVERTED_GRAYSCALE, showwedge, showwedge_label);
	  if(commandline_PgplotInterprator(argc, argv, application.verbose_state.debug) == 0)
	    return 0;
	  if(calculate_beam_widths != 0 && contour_plot) 
	    PlotContours(rhogrid, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels, TR, PSDeviceID, contour_txt, contourcolor, fixedContours, nruserContours, userContours, boxlw, 1);
	  if(calculate_interpulse_widths != 0 && contour_plot) 
	    PlotContours(rhogrid2, alphastart, alphaend, betastart, betaend, nalpha, nbeta, nrcontourlevels2, TR, PSDeviceID, contour_txt, contourcolorIP, fixedContours, nruserContours, userContours, boxlw, 1);
	  ppgclos();
	  pgcloseoutputfile();
	}
	cpgslct(GridDeviceID);
	break;				
//START REGION RELEASE
      case 104:
	PrintHelp();
	break;
      default:      printf("Unknown key: %d\n", c);
      }
    }while(c != 27 && c != 'q');
    if(showgraphics)
      ppgend();
    free(chigrid);
    free(l0grid);
    free(pa0grid);
    free(dhgrid);
    free(rhogrid);
    free(rhogrid2);
  }else {   /* No gridsearch */
    DoFitting(alpha0, beta0, pa0, dpa0, l0, dl0, dh0, ddh0, ftol, &fit_pa0, &fit_l0, &fit_alpha, &fit_beta, &fit_dh0, &chi, &nfunk, 0, 1, stderr, 0, nrofsigmas, &nrfitparams, amoeba_algorithm, application.verbose_state);
    alpha0 = fit_alpha;
    beta0 = fit_beta;
    if(application.verbose_state.verbose) printf("\n\n#Fitting procedure: ");
    for(i = 0; i < argc; i++)
      if(application.verbose_state.verbose) printf("%s ", argv[i]);
    if(application.verbose_state.verbose) printf("\n#\n#After %d steps the downhill-simplex found:\n#pa0 = %f and l0 = %f with chi2 = %e\n", nfunk, fit_pa0, fit_l0, chi);
    
    if(application.verbose_state.verbose) printf("#\n");
      
    if(printfit == 1) {
      dpa0 = paswing_double(alpha0, beta0, i, fit_pa0, fit_l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, fit_dh0);
      for(i = 0; i < 360; i++) {
	pa0 = paswing_double(alpha0, beta0, i, fit_pa0, fit_l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, fit_dh0);
	if(printsmooth) {
	  if(pa0 - dpa0 > 90)
	    pa0 -= 180;
	  else if(dpa0 - pa0 > 90)
	    pa0 += 180;
	  if(pa0 - dpa0 > 90)
	    pa0 -= 180;
	  else if(dpa0 - pa0 > 90)
	    pa0 += 180;
	  dpa0 = pa0;
	}
	if(i == 0)
	  printf("# BinNr PA[deg] alpha[deg] beta[deg] longitude_0[deg] PA_0[deg]\n");
	printf("%d %f %f %f %f %f\n", i, pa0, alpha0, beta0, fit_l0, fit_pa0);
      }
    }
  }

  if(loadresults == 0) {
    closePSRData(&datain, 0, application.verbose_state);
  }
  if(macrofilename)
    fclose(macrofile);
  terminateApplication(&application);
  gsl_rng_free(rand_num_gen);
  return 0;
}  // End of main

/* parameters are: PA_0 (deg), l_0 (deg), a(deg), b(deg), dh */
double funk(double x[])
{
  double pa, y1, dy, chi2, alpha, beta, heightshift;
  int i;
//START REGION DEVELOP
  double rho, upp, downn, dchi2;
  double r2d = 180.0/M_PI;
  double d2r = M_PI/180.0;
//START REGION RELEASE

  //  fflush(stdout); fprintf(stderr, "funk: pa0=%lf l0=%lf a=%lf b=%lf dh=%lf\n", x[0], x[1], x[2], x[3], x[4]);
  
  alpha = x[2];
  beta = x[3];
  heightshift = x[4];
  chi2 = 0;

  if(alpha == 0.0000)
    alpha = 0.00001;
  if(beta == 0.000)
    beta = 0.00001;

  /* If magnetic axis is not in allowed range, set chi2 to large number */
  if(fabs(fitterinfo.l0_start-x[1]) > fitterinfo.max_l0_diff) {
    //    fflush(stdout); fprintf(stderr, "fitterinfo.l0_start=%f\n", fitterinfo.l0_start);
    return 1e10;
  }

  /* Check if PA-swing goes through defined point */
  if(fitterinfo.force_set) {
    pa = paswing_double(alpha, beta, fitterinfo.force_l, x[0], x[1], fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, heightshift);
    if(fitterinfo.autojump)
      dy = dy_90(pa, fitterinfo.force_pa);
    else
      dy = dy_180(pa, fitterinfo.force_pa);
    if(dy > fitterinfo.force_dpa)
      return 1e10;
  }

  for(i = 0; i < fitterinfo.NrDataPoints; i++) {
    //    printf("XXXX %lf\n", heightshift);
    //    printf("XXXX %lf\n", fitterinfo.add_height_longitude);
    pa = paswing_double(alpha, beta, fitterinfo.data_l[i], x[0], x[1], fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, heightshift);
    y1 = fitterinfo.data_pa[i];
    if(fitterinfo.autojump)
      dy = dy_90(pa, y1);
    else
      dy = dy_180(pa, y1);
    chi2 += dy*dy/((double)(fitterinfo.data_dpa[i])*(double)(fitterinfo.data_dpa[i]));
    
//START REGION DEVELOP
    if(fitterinfo.do_rho_width){
      /*     upp, downn = error propagation from W to rho*/
      upp = 0.25*sin(d2r*alpha)*sin(d2r*(alpha+beta))*sin(d2r*(fitterinfo.pulse_width*0.5))*sin(d2r*alpha)*sin(d2r*(alpha+beta))*sin(d2r*(fitterinfo.pulse_width*0.5))*fitterinfo.sigma_width*fitterinfo.sigma_width;
      downn = 1.0 - ( (cos(d2r*alpha)*cos(d2r*(alpha+beta)) + sin(d2r*alpha)*sin(d2r*(alpha+beta))*cos(d2r*(fitterinfo.pulse_width*0.5))) * (cos(d2r*alpha)*cos(d2r*(alpha+beta)) + sin(d2r*alpha)*sin(d2r*(alpha+beta))*cos(d2r*(fitterinfo.pulse_width*0.5))));

      /*       chi2 from w-rho */
      rho = r2d*acos( cos(alpha*d2r)*cos((alpha+beta)*d2r) + sin(alpha*d2r)*sin((alpha +beta)*d2r)*cos(fitterinfo.pulse_width*d2r*0.5));
      // Patrick: abs -> fabs? It is an option written by Math, so I didn't change it. It is probably ok as it is compared with 1.
      // Patrick: eventually changed it to fabs to avoid compiler warnings
      if (fabs( cos(alpha*d2r)*cos((alpha+beta)*d2r) + sin(alpha*d2r)*sin((alpha +beta)*d2r)*cos(fitterinfo.pulse_width*d2r*0.5)) > 1)
	printwarning(0, "something wrong in chi2_rho");
 
      dchi2 = (fitterinfo.rho_bcw - rho)*(fitterinfo.rho_bcw - rho)/( (fitterinfo.sigma_rho*fitterinfo.sigma_rho) + (upp/downn)*(upp/downn) ) ;
     
/*       printf("Value a,b = %f, %f, ccc=  %f\n", alpha, beta,  (cos(alpha*d2r)*cos((alpha+beta)*d2r) + sin(alpha*d2r)*sin((alpha +beta)*d2r)*cos(fitterinfo.pulse_width*d2r*0.5)) ); */
      chi2 += dchi2;
      
    }
//START REGION RELEASE

  }


  //  fflush(stdout); fprintf(stderr, "funk:      chi2 = %f\n", chi2);
  return chi2;
}

//START REGION DEVELOP
//START REGION RELEASE


// parameters are: PA_0 (deg), l_0 (deg), a(deg), b(deg), dh 
double internal_PA0_funk(double pa0, void *params)
{
  double chi2;
  double *x;

  x = params;
  x[0] = pa0;
  
  chi2 = funk(x);
  return chi2;
}

//START REGION DEVELOP
//START REGION RELEASE

// parameters are: PA_0 (deg), l_0 (deg), a(deg), b(deg), dh 
double internal_L0_funk(double l0, void *params)
{
  double chi2;
  double *x;

  x = params;
  x[1] = l0;
  
  chi2 = funk(x);

  //  printf("YYYYY trial %f: chi2=%f\n", x[1], chi2);

  return chi2;
}

//START REGION DEVELOP
//START REGION RELEASE

// parameters are: PA_0 (deg), l_0 (deg), a(deg), b(deg), dh 
double internal_funk_gsl(double *x, void *params)
{
  int ret, nrpoints, debug_verbose, debug_verbose2, trial, higher_res_step;
  double chi2, chi2_before, pa0, l0, xnew[5], ftol;

  nrpoints = 200;  // 25 works often works though, 100 not always enough
  ftol =  0.00001;


  debug_verbose = 0;
  debug_verbose2 = 0;

  for(higher_res_step = 0; higher_res_step < 2; higher_res_step++) {
    memcpy(xnew, x, 5*sizeof(double));
    if(higher_res_step != 0) {
      fflush(stdout);
      printwarning(0, "Failed to find global minimum, trying to use a higher resolution grid");
      nrpoints *= 100;
    }
    //    if(x[2] >41.8 && x[2] < 41.9 && x[3] > 0.10 && x[3] < 0.11) {
    //      debug_verbose = 1;
    //      debug_verbose2 = 1;
    //      nrpoints *= 100;
    //    }

    chi2_before = funk(xnew);
    //    if(debug_verbose) {
    //      printf("XXXX internal_funk_gsl called with: pa0=%f l0=%f a=%f b=%f dh=%f\n", x[0], x[1], x[2], x[3], x[4]);
    //      printf("XXXX chi2 before optimization: %f\n", chi2_before);
    //    }


    for(trial = 0; trial < 2; trial++) {
      if(internal_fit_pa_or_l0 == 1)
	ret = minimize_1D_double(0, internal_L0_funk, x, fitterinfo.l0_start-fitterinfo.max_l0_diff, fitterinfo.l0_start+fitterinfo.max_l0_diff, nrpoints, 1, 2, &l0, 2000, ftol, 0.0, debug_verbose, debug_verbose2);
      else
	ret = minimize_1D_double(0, internal_PA0_funk, x, 0, 180, nrpoints, 0, 2, &pa0, 2000, ftol, 0.0, debug_verbose, debug_verbose2);
      if(ret == 1) {
	fflush(stdout);
	ftol *= 100;
	printwarning(debug_verbose, "WARNING: Maximum nr of itterations exceeded in internal_funk_gsl (alpha=%f, beta=%f)\nTrying to lower tolerance to %f", x[2], x[3], ftol);
      }else {
	break;
      }
    }
    if(ret != 0) {
      if(ret == 1) {
	fflush(stdout);
	printwarning(debug_verbose, "WARNING: Maximum nr of itterations exceeded in internal_funk_gsl (alpha=%f, beta=%f)\nTrying to lower tolerance", x[2], x[3]);
      }else if(ret == 2) {
	fflush(stdout);
	printwarning(debug_verbose, "WARNING: Did not found root");
      }else if(ret == 3) {
	fflush(stdout);
	printwarning(debug_verbose, "WARNING: Lower and upper limit do not bracket a root in internal_funk_gsl");
      }else {
	fflush(stdout);
	printwarning(debug_verbose, "WARNING: Unknown error in minimize_1D_double in internal_funk_gsl");
	exit(0);
      }
      //  }else {
      //      fprintf(stderr, "Success\n");
    }
    if(internal_fit_pa_or_l0 == 1) {
      xnew[1] = l0;
    }else {
      xnew[0] = pa0;
    }
    chi2 = funk(xnew);
    //    if(debug_verbose) {
    //      printf("XXXX chi2 after optimization: %f\n", chi2);
    //    }
    // Put in a bit of a margin, you don't want it to fail if exactly the same
    if(1.05*chi2_before < chi2) {
      fflush(stdout);
      printwarning(0, "Optimization failed: %f > %f for alpha=%f beta=%f", chi2, chi2_before, x[2], x[3]);
      if(higher_res_step > 0) // Allow the routine to try again with a higher resolution grid first
	exit(0);
    }else {
      break;
    }
  }
  return chi2;
}

//START REGION DEVELOP
//START REGION RELEASE

void DoFitting(double alpha0, double beta0, double pa0, double dpa0, double l0, double dl0, double dh0, double ddh0, double ftol, double *fit_pa0, double *fit_l0, double *fit_a, double *fit_b, double *fit_dh, double *chi, int *nfunk, int searchAll, int report, FILE *reportStream, int finderrors, double nrofsigmas, int *nfitparameters, int amoeba_algorithm, verbose_definition verbose)
{
  double xstart[5], dx[5], xfit[5], dplus[5], dmin[5], chi_d;
  int fixed[5];

  fitterinfo.l0_start = l0;

  xstart[0] = pa0;
  xstart[1] = l0;
  xstart[2] = alpha0;
  xstart[3] = beta0;
  xstart[4] = dh0;

  dx[0] = dpa0;
  dx[1] = dl0;
  dx[2] = 10;
  dx[3] = 10;
  dx[4] = ddh0;

  fixed[0] = 0;
  fixed[1] = 0;
  fixed[2] = 0;
  fixed[3] = 0;
  fixed[4] = 0;

  if(searchAll) {
    *nfitparameters = 4;
  }else {
    fixed[2] = 1;
    fixed[3] = 1;
    *nfitparameters = 2;
  }

  //  printf("XXXXX init: %lf\n", dh0);
  if(ddh0 >= -2) {
    (*nfitparameters) += 1;  // We actually want to get the emission height
  }else {         // Just a fixed parameter
    fixed[4] = 1;
    xstart[4] = dh0;
  }
//START REGION DEVELOP
  /* Keep PA fixed when pressing E */
  if(finderrors == 2) {
    fixed[0] = 1;
  }
//START REGION RELEASE

  if(finderrors && *nfitparameters <= 2) {
    printwarning(verbose.debug, "WARNING: Need at least three fit parameters to do error estimation.");
    finderrors = 0;
  }
  /* Try amoeba search, if fails to converge in reasonable time it
     tries again with a smaller tollerance. */
  do {
    if(finderrors) {
      printf("param0 = pa0, param1 = l0, param2 = alpha, param3 = beta, param4=dh\n");
    }
    if(doAmoeba_d(amoeba_algorithm, xstart, dx, fixed, xfit, &chi_d, 5, funk, ftol, nfunk, 0, finderrors, nrofsigmas, dplus, dmin) == 1) {
      printwarning(verbose.debug, "WARNING: Adjusting downhill-simplex tollerance to try to converge.");
      ftol *= 10;
    }else {
      *chi = chi_d;
      break;
    }
  }while(ftol < 0.01);

  *fit_pa0 = xfit[0];
  *fit_l0  = xfit[1];
  *fit_a   = xfit[2];
  *fit_b   = xfit[3];
  *fit_dh  = xfit[4];
  if(fixed[4]) { 
    if(fitterinfo.add_height_longitude <= 360)
      *fit_dh = dh0;
    else 
      *fit_dh = 0;
  }

  if(report) {
    if(searchAll == 0) {
      fprintf(reportStream, "Fit using specified alpha and beta value:");
    }else {
      fprintf(reportStream, "Fit keeping alpha and beta as a free parameter:");
    }
    fprintf(reportStream, "\n     alpha = %15f deg", *fit_a);
    if(finderrors) fprintf(reportStream, " (%.1f sigma error: %+15f %15f -> range = %f to %f)", nrofsigmas, dplus[2], dmin[2], *fit_a+dmin[2], *fit_a+dplus[2]);
    fprintf(reportStream, "\n     beta  = %15f deg", *fit_b);
    if(finderrors) fprintf(reportStream, " (%.1f sigma error: %+15f %15f -> range = %f to %f)", nrofsigmas, dplus[3], dmin[3], *fit_b+dmin[3], *fit_b+dplus[3]);
    fprintf(reportStream, "\n           = %15f for the other pole", *fit_b+2.0*(*fit_a)-180.0);
    fprintf(reportStream, "\n     l0    = %15f deg", *fit_l0);
    if(finderrors) fprintf(reportStream, " (%.1f sigma error: %+15f %15f -> range = %f to %f)", nrofsigmas, dplus[1], dmin[1], *fit_l0+dmin[1], *fit_l0+dplus[1]);
    fprintf(reportStream, "\n     pa0   = %15f deg", *fit_pa0);
    if(finderrors) fprintf(reportStream, " (%.1f sigma error: %+15f %15f -> range = %f to %f)", nrofsigmas, dplus[0], dmin[0], *fit_pa0+dmin[0], *fit_pa0+dplus[0]);
    fprintf(reportStream, "\n     dh    = %15f", *fit_dh);
    if(finderrors) fprintf(reportStream, " (%.1f sigma error: %+15f %15f -> range = %f to %f)", nrofsigmas, dplus[4], dmin[4], *fit_dh+dmin[4], *fit_dh+dplus[4]);
    fprintf(reportStream, "\n     reduced chi^2=%f (tot=%f) %d params and %d points\n", *chi/(double)(fitterinfo.NrDataPoints-(*nfitparameters)), *chi, *nfitparameters, fitterinfo.NrDataPoints);
  }
}

//START REGION DEVELOP
//START REGION RELEASE

void PlotGrid(float *chigrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, double level, double suppress_fac, int GridDeviceID, double alpha, double beta, double lwbox, double labelcharheight, double boxlabelcharheight, int drawCross, int draw_title, int drawcontours, int nogray, 
//START REGION DEVELOP
	      int usezeta, int useslope, int use_sigma, 
//START REGION RELEASE
	      double chimax, double chimin, int maptype, int showwedge, char *showwedge_label)
{
  int i;
//START REGION DEVELOP
  int j;
//START REGION RELEASE
  char txt1[100], txt2[100], txt3[100];
  float contours[200];

  verbose_definition noverbose;
  cleanVerboseState(&noverbose);
  noverbose.nocounters = 1;

/*   for system call calcchi(.f) for 1-2-3 sigma level */
//START REGION DEVELOP
  if(use_sigma){
    double sig_level;
    char call_calcchi[999];
    char rm_temp[999];
    sprintf(rm_temp,"rm temp_calcchi");
    system(rm_temp);
    //  fprintf(stderr,"calling calcchi\n");
    for( j=1; j<4; j++){
      sprintf(call_calcchi,"calcchi %d %d >> temp_calcchi",(fitterinfo.NrDataPoints-4), j);
      system(call_calcchi);
    }
    FILE *input_call;
    input_call = fopen("temp_calcchi", "r");
    if (input_call == NULL) {
      printerror(noverbose.debug, "Can't open \"temp_calcchi\" !!!");
      exit(1);
    }
    for( i=0 ; i<3 ; i++) {
      fscanf(input_call,"%lf", &sig_level);
      contours[i] = (chimax - (chimin+(sig_level/fitterinfo.NrDataPoints)) )/(chimax - chimin);
    }
    fclose(input_call);
  }else{
//START REGION RELEASE
    for(i = 0; i < drawcontours; i++) {
      contours[i] = (chimax - (i+2)*chimin)/(chimax - chimin);
    }
//START REGION DEVELOP
  }    
//START REGION RELEASE
  

  if(drawcontours > 200) {
    printwarning(noverbose.debug, "Maximum allowed number contours is 200");
    drawcontours = 200;
  }
  if(draw_title) {
    sprintf(txt1, "\\(2148)\\u2\\d grid");
  }else {
    txt1[0] = 0;
  }
  sprintf(txt2, "\\(2128) [deg]");
//START REGION DEVELOP
  if(usezeta)
    sprintf(txt2, "\\(2132) [deg]");
  if(useslope == 1)
    sprintf(txt2, "sin(\\(2127))/sin(\\(2128))");
  if(useslope == 2)
    sprintf(txt2, "sin(\\(2128))/sin(\\(2127))");
//START REGION RELEASE
  sprintf(txt3, "\\(2127) [deg]");
 
//START REGION DEVELOP
  if( fitterinfo.math_setting == 2)
    maptype = PPGPLOT_INVERTED_GRAYSCALE;
//START REGION RELEASE

  ppgask(0);
  ppgpage();
  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);

  pgplot_options.viewport.dontopen  = 1;
  pgplot_options.viewport.dontclose = 1;
  strcpy(pgplot_options.box.xlabel, txt3);
  strcpy(pgplot_options.box.ylabel, txt2);
  strcpy(pgplot_options.box.title, txt1);
  pgplot_options.box.label_ch = labelcharheight;
  pgplot_options.box.box_labelsize = boxlabelcharheight;
  pgplot_options.box.box_lw = lwbox;
  pgplot_options.box.label_lw = lwbox;

  /*  pgplotMap(&pgplot_options, chigrid, nalpha, nbeta, alphastart, alphaend, alphastart, alphaend, betastart, betaend, betastart, betaend, maptype, 0, nogray, drawcontours, contours, lwbox, 0, 1, 1, level, suppress_fac, 0, 0, txt3, txt2, txt1, 1, 1, 1, 0, 0, labelcharheight, boxlabelcharheight, lwbox, 0, 0, 0, 0, 0, 0, 0, 0); */
  pgplotMap(&pgplot_options, chigrid, nalpha, nbeta, alphastart, alphaend, alphastart, alphaend, betastart, betaend, betastart, betaend, maptype, 0, nogray, drawcontours, contours, lwbox, 0, 1, 1, level, suppress_fac, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, noverbose);


  if(showwedge) {
    /*
    printf("XXXXX %f %f\n", level, suppress_fac);
    printf("XXXXX %f %f\n", (1.0-level)*(chimax-chimin)+chimin, (1.0-suppress_fac)*(chimax-chimin)+chimin);
    */
    ppgsch(boxlabelcharheight*0.5);
    ppgslw(lwbox);
    ppgwedg("RI", 0, 6, (1.0-suppress_fac)*(chimax-chimin)+chimin, (1.0-level)*(chimax-chimin)+chimin, showwedge_label);
    ppgslw(1);
    ppgsch(1);
  }

  if(alpha >= 0 && drawCross) {
    if(drawCross == 1)
      ppgsci(2);
    else 
      ppgsci(0);
    ppgmove(alpha-2, beta);
    ppgdraw(alpha+2, beta);
    ppgmove(alpha, beta-2);
    ppgdraw(alpha, beta+2);
    ppgsci(1);
  }
  ppgslw(1);
}

//START REGION DEVELOP
//START REGION RELEASE

/* Sorts the alpha/beta points into contours. A contour is defined to
   be a collection of points seperated by less than maxsep from each
   other. The contour number of each point is stored in contnr and the
   number of found contours is returned. */
int define_contours_from_list_points(float *alpha, float *beta, int *contnr, int npoints, double maxsep)
{
  int i, j, curcont, found, findnextcontour;
  double sep;
  if(npoints < 1) {
    return 0;
  }
  for(i = 0; i < npoints; i++) {
    contnr[i] = 0;
  }


  /* Define first point with undefined contour nr to be part of a new contour */
  curcont = 0;
  do {
    curcont++;
    found = 0;
    findnextcontour = 0;
    for(i = 0; i < npoints; i++) {
      if(contnr[i] == 0) {
	contnr[i] = curcont;
	found = 1;
	break;
      }
    }
    if(found == 0)
      break;              /* All points are already asigned to a contour */
    else
      findnextcontour = 1;
    while(found) {
      found = 0;
      for(i = 0; i < npoints; i++) {
	if(contnr[i] == 0) {
	  for(j = 0; j < npoints; j++) {
	    if(contnr[j] == curcont) {
	      sep = sqrt((alpha[j]-alpha[i])*(alpha[j]-alpha[i])+(beta[j]-beta[i])*(beta[j]-beta[i]));
	      if(sep <= maxsep) {
		contnr[i] = curcont;
		found = 1;
		break;
	      }
	    }
	  }
	}
      }
    }
    found = 1;
  }while(findnextcontour);   

  return curcont-1;
}

//START REGION DEVELOP
//START REGION RELEASE

void calcIntersectionRhoAndBanana(float *rhogrid, float *chigrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, int nruserContours, float *userContours, double chimax, double chimin, verbose_definition verbose)
{
  int i, j, k, ialpha, ibeta, ialpha2, ibeta2, nfound, nfoundbest, *contnr, nrcontoursfound;
  float *alpha, *beta, *chi, chibest, maxsep;

  alpha = malloc(2*nalpha*nbeta*sizeof(float));
  beta = malloc(2*nalpha*nbeta*sizeof(float));
  chi = malloc(2*nalpha*nbeta*sizeof(float));
  contnr = malloc(2*nalpha*nbeta*sizeof(int));
  if(alpha == NULL || beta == NULL || chi == NULL || contnr == NULL) {
    printerror(verbose.debug, "ERROR ppolFit: Cannot allocate memory");
    return;
  }


  if(nruserContours > 0) {
    for(i = 0; i < nruserContours; i++) {
      nfound = 0;
      for(ialpha = 0; ialpha < nalpha; ialpha++) {
	for(ibeta = 0; ibeta < nbeta; ibeta++) {
	  /* See when making step to right we cross the contour */
	  ialpha2 = ialpha + 1;
	  ibeta2 = ibeta;
	  if(ialpha2 >= nalpha)
	    ialpha2 = ialpha;
	  if((rhogrid[nalpha*ibeta+ialpha] < userContours[i] && rhogrid[nalpha*ibeta2+ialpha2] > userContours[i]) || (rhogrid[nalpha*ibeta+ialpha] > userContours[i] && rhogrid[nalpha*ibeta2+ialpha2] < userContours[i])) {
	    alpha[nfound] = 0.5*(ialpha+ialpha2)*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
	    beta[nfound] = 0.5*(ibeta+ibeta2)*(betaend-betastart)/(double)(nbeta-1)+betastart;
	    chi[nfound] = 0.5*(chigrid[nalpha*ibeta+ialpha] + chigrid[nalpha*ibeta2+ialpha2]);
	    if(chi[nfound] > (chimax - 4*chimin)/(chimax - chimin))
	       nfound++;
	  }

	  /* See when making step to bottom we cross the contour */
	  ialpha2 = ialpha;
	  ibeta2 = ibeta+1;
	  if(ibeta2 >= nbeta)
	    ibeta2 = ibeta;
	  if((rhogrid[nalpha*ibeta+ialpha] < userContours[i] && rhogrid[nalpha*ibeta2+ialpha2] > userContours[i]) || (rhogrid[nalpha*ibeta+ialpha] > userContours[i] && rhogrid[nalpha*ibeta2+ialpha2] < userContours[i])) {
	    alpha[nfound] = 0.5*(ialpha+ialpha2)*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
	    beta[nfound] = 0.5*(ibeta+ibeta2)*(betaend-betastart)/(double)(nbeta-1)+betastart;
	    chi[nfound] = 0.5*(chigrid[nalpha*ibeta+ialpha] + chigrid[nalpha*ibeta2+ialpha2]);
	    if(chi[nfound] > (chimax - 4*chimin)/(chimax - chimin))
	       nfound++;
	  }

	}
      }
      if(nfound <= 0) {
	printf("No intersections found between 3-sigma banana and contour\n");
      }else {
	maxsep = sqrt(2.0*(alphaend-alphastart)*(alphaend-alphastart)/(double)((nalpha)*(nalpha)) + 2.0*(betaend-betastart)*(betaend-betastart)/(double)((nbeta)*(nbeta)));
	nrcontoursfound = define_contours_from_list_points(alpha, beta, contnr, nfound, maxsep);
	printf("Found %d cross points between contour %d and 3-sigma chi2 contour\n", nrcontoursfound, i+1);
	if(verbose.verbose) {
	  for(j = 0; j < nfound; j++) {
	    ppgsci(contnr[j]+4);
	    ppgpt1(alpha[j], beta[j], -1);	
	  }
	}
	for(k = 1; k <= nrcontoursfound; k++) {
	  nfoundbest = 0;
	  chibest = NAN;
	  for(j = 0; j < nfound; j++) {
	    if(contnr[j] == k) {
	      if(chi[j] > chibest || isnan(chibest)) {
		chibest = chi[j];
		nfoundbest = j;
	      }
	    }
	  }
	  if(verbose.verbose) {
	    ppgsci(8);
	    ppgslw(5);
	    ppgpt1(alpha[nfoundbest], beta[nfoundbest], 2);
	  }
	  chibest = 1.0-chibest;
	  chibest = chibest*(chimax-chimin) + chimin;
	  printf("Cross point %d: alpha=%f beta=%f chi=%f\n", k, alpha[nfoundbest], beta[nfoundbest], chibest);
	}
      }
    }
  }
  free(alpha);
  free(beta);
  free(chi);
  free(contnr);
  ppgsci(1);
  ppgslw(1);
}

//START REGION DEVELOP
//START REGION RELEASE

void PlotContours(float *rhogrid, double alphastart, double alphaend, double betastart, double betaend, int nalpha, int nbeta, int nrlevels, float *TR, int GridDeviceID, int contour_txt, int contourcolor, int fixedContours, int nruserContours, float *userContours, int lwbox, int dotted)
{
  float C[500];
  char txt[100];
  int i;

  if(fixedContours == 0) {
    for(i = 0; i < nrlevels; i++) 
      C[i] = 180.0*i/(double)(nrlevels-1.0);
  }else {
    C[0] = 1;
    C[1] = 2;
    C[2] = 3;
    C[3] = 4;
    C[4] = 5;
    C[5] = 10;
    C[6] = 20;
    C[7] = 30;
    C[8] = 40;
    C[9] = 50;
    C[10] = 60;
    C[11] = 70;
    C[12] = 80;
    C[13] = 90;
    C[14] = 100;
    C[15] = 110;
    C[16] = 120;
    C[17] = 130;
    C[18] = 140;
    C[19] = 150;
    C[20] = 160;
    C[21] = 170;
    C[22] = 180;
    nrlevels = 23;
  }
  if(nruserContours > 0) {
    for(i = 0; i < nruserContours; i++) 
      C[i] = userContours[i];
    nrlevels = nruserContours;
  }else if( nruserContours == -1){      /* Go into Math mode */
    C[0] = 5;
    C[1] = 10;
    C[2] = 15;
    C[3] = 20;
    C[4] = 25;
    C[5] = 30;
    nrlevels = 6;
  }

  ppgsci(contourcolor);
  ppgslw(lwbox);
  if(dotted)
    ppgsls(4);
  ppgcont(rhogrid, nalpha, nbeta, 1, nalpha, 1, nbeta, C, -nrlevels, TR);
  ppgsls(1);
  ppgslw(1);
  if(contour_txt) {
    for(i = 0; i < nrlevels; i++) {
      sprintf(txt, "%.0f", C[i]);    
      ppgconl(rhogrid, nalpha, nbeta, 1, nalpha, 1, nbeta, C[i], TR, txt, nalpha, 0.01*nalpha);
    }
  }
  ppgsci(1); 
}

//START REGION DEVELOP
//START REGION RELEASE

double dy_180(double y1, double y2)
{
  double dy;
  y1 = derotate_180_double(y1);
  y2 = derotate_180_double(y2);

  dy = fabs(y1-y2);
  if(fabs(y1-y2-180) < dy) {
    y2 += 180.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-180) < dy) {
    y2 += 180.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-180) < dy) {
    y2 += 180.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+180) < dy) {
    y2 -= 180.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+180) < dy) {
    y2 -= 180.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+180) < dy) {
    y2 -= 180.0;
    dy = fabs(y1-y2);
  }
  return dy;
}

double dy_90(double y1, double y2)
{
  double dy;
  y1 = derotate_180_double(y1);
  y2 = derotate_180_double(y2);

  dy = fabs(y1-y2);
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2-90) < dy) {
    y2 += 90.0;
    dy = fabs(y1-y2);
  }

  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  if(fabs(y1-y2+90) < dy) {
    y2 -= 90.0;
    dy = fabs(y1-y2);
  }
  return dy;
}

//START REGION DEVELOP
//START REGION RELEASE

void PlotPAswing(double alpha, double beta, double pa0, double l0, int PlotFit, double leftPulseLongitude, double rightPulseLongitude, double dh)
{
  int i;
  double oldpa, newpa;

  ppgpage();
  ppgask(0);
  ppgslw(1);
  ppgsvp(0.1, 0.9, 0.1, 0.9);
  ppgswin(leftPulseLongitude, rightPulseLongitude, 0, 180);
  ppglab("Pulse longitude", "PA", "PA-fit");
  ppgbox("bcnsti",0.0,0,"bcnmsti",0.0,0);
  
  for(i = 0; i < fitterinfo.NrDataPoints; i++) {
    ppgerr1(6, fitterinfo.data_l[i], derotate_180(fitterinfo.data_pa[i]), fitterinfo.data_dpa[i], 3); 
//START REGION DEVELOP
    if( fitterinfo.math_setting == 1){
      ppgerr1(6, fitterinfo.data_l[i], derotate_180(fitterinfo.data_pa[i])+90, fitterinfo.data_dpa[i], 3); 
      ppgerr1(6, fitterinfo.data_l[i], derotate_180(fitterinfo.data_pa[i])-90, fitterinfo.data_dpa[i], 3); 
    }
//START REGION RELEASE
  }

  if(PlotFit) { 
    ppgsci(2);
    /*    oldpa = paswing(0,pa0,l0,0); */
    oldpa = paswing_double(alpha, beta, 0, pa0, l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, dh);
    ppgmove(0, oldpa);
    for(i = 1; i < 3600; i++) {
      /*      newpa = paswing(0.1*i,pa0,l0,0); */
      newpa = paswing_double(alpha, beta, 0.1*i, pa0, l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, dh);
      if(fabs(newpa-oldpa) < 100)
	ppgdraw(0.1*i, newpa);
      else
	ppgmove(0.1*i, newpa);
      oldpa = newpa;
    }
    if(l0 > 360)
      l0 -= 360;
    if(l0 < 0)
      l0 += 360;
    ppgsci(3);
    ppgmove(l0, 0);
    ppgdraw(l0, 180);
    l0 += 180;
    ppgsci(4);
    if(l0 > 360)
      l0 -= 360;
    ppgmove(l0, 0);
    ppgdraw(l0, 180);
    ppgsci(1);
  }
}

//START REGION DEVELOP
//START REGION RELEASE

void PrintHelp()
{
  printf("General options:\n");
  printf("  ESC (or q): Exit\n");
  printf("  w:   Write out postscript files or other file\n");
//START REGION DEVELOP
  printf("  W:   Write out pgplotInterpretor files\n");
//START REGION RELEASE
  printf("Fit options:\n");
  printf("  b:   Find and plot best solution\n");
  printf("  e:   Toggle error estimation on best fit\n");
//START REGION DEVELOP
  printf("  E:   Toggle error estimation on best fit (while keeping PA_0 fixed)\n");
//START REGION RELEASE
  printf("  S:   Set the number of sigma's used for error calculation\n");
  printf("  G:   Find allowed range of fit parameters by considering the alpha-beta grid search\n");
  printf("Plot refinement options:\n");
  printf("  C:   Toggle showing of the cross of the best fit\n");
  printf("  g:   Toggle drawing chi2 surface in grayscale\n");
  printf("  l:   Set upper chi^2 level (for grayscale)\n");
//START REGION DEVELOP
  printf("  L:   Set upper chi^2 level as fraction of max (for grayscale)\n");
//START REGION RELEASE
  printf("  o:   Change output plot parameters (such as line widths)\n");
//START REGION DEVELOP
  printf("  p:   Plot a RVM model with given parameters\n");
//START REGION RELEASE
  printf("  s:   Set suppress factor of grayscale\n");
  printf("  T:   Toggle showing of the title in the chi^2 plot\n");
  printf("  z:   Zoom in on pulse longitude\n");
  printf("Contour options:\n");
  printf("  c:   Plot chi^2 contours\n");
  printf("  f:   Use build-in non-equally spaced specification of contour levels for beam radius\n");
  printf("  F:   Manually specify shown contour levels for beam radius\n");
  printf("  n:   Set nr of equally spaced contour levels for beam radius\n");
  printf("  r:   Toggle plotting contours of the beam radius\n");
  printf("  R:   Input pulse widths and rho values from ascii file which might help to get an idea of the allowed spread in the contours\n");
  printf("  t:   Switch drawing of labels for countours\n");
  printwarning(0, "WARNING: -maxdl is set to %f. Changing this parameter might change the results considerably.\n\n", fitterinfo.max_l0_diff);
}

//START REGION DEVELOP
//START REGION RELEASE

void print_steepness(double alpha, double beta, double l0, double pa0, int verbose, double dh, double *sina_b)
{
  int i;
  double l, l1, l2, psi, psi_old, dpsi, dpsidphi, max_dpsidphi, resolution;
  *sina_b = sin(alpha*M_PI/180.0)/sin(beta*M_PI/180.0);
  printf("\nsin(a)/sin(b)   = %f\n", sin(alpha*M_PI/180.0)/sin(beta*M_PI/180.0));
  max_dpsidphi = 0;

  l1 = l2 = -1;

  resolution = fitterinfo.data_l[1] - fitterinfo.data_l[0];
  for(i = 1; i < fitterinfo.NrDataPoints; i++) {
    if(fitterinfo.data_l[i] - fitterinfo.data_l[i-1] < resolution)
      resolution = fitterinfo.data_l[i] - fitterinfo.data_l[i-1];
  }
  if(verbose) printf("Resolution = %f degrees\n", resolution);
  for(i = 0; i < fitterinfo.NrDataPoints; i++) {
    l = fitterinfo.data_l[i];

/*    psi = paswing(l, pa0, l0, 0); */
    psi = paswing_double(alpha, beta, l, pa0, l0, fitterinfo.nrJumps, fitterinfo.jump_longitude, fitterinfo.jump_offset, fitterinfo.add_height_longitude, dh);
    if(i > 0) {
      dpsi = psi-psi_old;
      dpsidphi = dpsi/(fitterinfo.data_l[i]-fitterinfo.data_l[i-1]);
      if(fabs(dpsidphi) > fabs(max_dpsidphi)) {
	if(fabs(dpsi) > fabs(dpsi-180)) {
	  dpsi -= 180;
	}
	if(fabs(dpsi) > fabs(dpsi+180)) {
	  dpsi += 180;
	}
	dpsidphi = dpsi/(fitterinfo.data_l[i]-fitterinfo.data_l[i-1]);
	if(fabs(dpsidphi) > fabs(max_dpsidphi) && (fitterinfo.data_l[i] - fitterinfo.data_l[i-1] < 2.5*resolution)) {
	  max_dpsidphi = dpsidphi;
	  l1 = fitterinfo.data_l[i-1];
	  l2 = fitterinfo.data_l[i];
	  /*	  printf("l=%f %f\n", l, dpsi); */
	}
      }
    }
    psi_old = psi;
  }
  printf("measured steepest part = %f (between longitude %f - %f)\n", max_dpsidphi, l1, l2);
}

//START REGION DEVELOP

/*
  Calculates the correction factor to go from beam width rho (magnetic
  latitude of end points of pulse width along line of sight) to the
  beam width (in the direction of the minor axis) of an elliptical
  beam with ellipticy ell. A positive ellipticity means elongation
  along the rotational meridian. A negative value an elongation in the
  rotational direction.

  Angles are in deg.
 */
double beamcorrection(double zeta, double width, double rho, double ell, double alpha)
{
  double phi, corr;

  zeta *= M_PI/180.0;
  width *= M_PI/180.0;
  rho *= M_PI/180.0;

  /* phi is the magnetic longitude of the intersection points. For the
     cos^2 and sin^2 used below the sign and the quadrant make no
     difference. */
  phi = asin(sin(zeta)*sin(0.5*width)/sin(rho));
  if(ell >= 0)
    corr = sqrt((1.0-ell*ell)*cos(phi)*cos(phi)+sin(phi)*sin(phi));
  else
    corr = sqrt((1.0-ell*ell)*sin(phi)*sin(phi)+cos(phi)*cos(phi));
  //  fprintf(stderr, "XXXXX corr=%f phi=%f W=%f rho=%f at alpha=%f and beta=%f\n", corr, phi*180/M_PI, width*180/M_PI, rho*180/M_PI, alpha*180/M_PI, (zeta-alpha)*180/M_PI);
  return corr;
}

//START REGION DEVELOP
//START REGION RELEASE

// Converts from given alpha + "beta which can be different than beta"
// to actual alpha and beta.
void convertAlphaBeta(double *alpha, double *beta
//START REGION DEVELOP
		      , int useslope, int usezeta
//START REGION RELEASE
		      )
{
  double alpha0, beta0;
  alpha0 = *alpha;
  beta0 = *beta;
//START REGION DEVELOP
  /* beta = zeta - alpha */
  if(usezeta) {                                                    
    beta0 -= alpha0;
  }
  if(useslope == 1) {
    /* sin(beta) =  sin(alpha)/slope */
    if(beta0 == 0)
      beta0 = 1e-5;
    beta0 = sin(M_PI*alpha0/180.0)/beta0;
    if(beta0 < -1)
      beta0 = -1;
    if(beta0 > 1)
      beta0 = 1;
    beta0 = asin(beta0)*180.0/M_PI;
  }else if(useslope == 2) {
    /* sin(beta) =  sin(alpha)*y */
    beta0 = sin(M_PI*alpha0/180.0)*beta0;
    if(beta0 < -1)
      beta0 = -1;
    if(beta0 > 1)
      beta0 = 1;
    beta0 = asin(beta0)*180.0/M_PI;
  }
//START REGION RELEASE
  *alpha = alpha0;
  *beta = beta0;
}

//START REGION DEVELOP
//START REGION RELEASE

void calcBeamWidths(int nalpha, int nbeta, double alphastart, double alphaend, double betastart, double betaend, 
//START REGION DEVELOP
		    int usezeta, int useslope, double elliptical_beam, 
//START REGION RELEASE
		    int calculate_beam_widths, int calculate_interpulse_widths, double pulse_width2, float *rhogrid, float *rhogrid2, int nocounters)
{
  int i, j;
  double dummyf, chi;
  double alpha0, beta0;
  for(i = 0; i < nalpha; i++) {
    for(j = 0; j < nbeta; j++) {
      alpha0 = i*(alphaend-alphastart)/(double)(nalpha-1)+alphastart;
      beta0 = j*(betaend-betastart)/(double)(nbeta-1)+betastart;
      convertAlphaBeta(&alpha0, &beta0
//START REGION DEVELOP
		       , useslope, usezeta
//START REGION RELEASE
		       );
      if(calculate_beam_widths) {
	dummyf =  cos(alpha0*M_PI/180.0)*cos((alpha0+beta0)*M_PI/180.0) + sin(alpha0*M_PI/180.0)*sin((alpha0+beta0)*M_PI/180.0)*cos(0.5*fitterinfo.pulse_width*M_PI/180.0);
	chi = acos(dummyf)*180.0/M_PI;
//START REGION DEVELOP
	if(!isnan(elliptical_beam)) {
	  chi *= beamcorrection(alpha0+beta0, fitterinfo.pulse_width, chi, elliptical_beam, alpha0*M_PI/180.0);
	}
//START REGION RELEASE
	rhogrid[nalpha*j+i] = chi;
      }
      if(calculate_interpulse_widths == 1) {
	dummyf = -cos(alpha0*M_PI/180.0)*cos((alpha0+beta0)*M_PI/180.0) + sin(alpha0*M_PI/180.0)*sin((alpha0+beta0)*M_PI/180.0)*cos(0.5*pulse_width2*M_PI/180.0);
	chi = acos(dummyf)*180.0/M_PI;
//START REGION DEVELOP
	if(!isnan(elliptical_beam)) {
	  chi *= beamcorrection(alpha0+beta0, pulse_width2, chi, elliptical_beam, alpha0*M_PI/180.0);
	}
//START REGION RELEASE
	rhogrid2[nalpha*j+i] = chi;
      }
      if(calculate_interpulse_widths == 2) {
	beta0 = -2.0*alpha0-beta0;
	dummyf =  cos(alpha0*M_PI/180.0)*cos((alpha0+beta0)*M_PI/180.0) + sin(alpha0*M_PI/180.0)*sin((alpha0+beta0)*M_PI/180.0)*cos(0.5*pulse_width2*M_PI/180.0);
//START REGION DEVELOP
	if(!isnan(elliptical_beam)) {
	  chi *= beamcorrection(alpha0+beta0, pulse_width2, chi, elliptical_beam, alpha0*M_PI/180.0);
	}
//START REGION RELEASE
	chi = acos(dummyf)*180.0/M_PI;
	rhogrid2[nalpha*j+i] = chi;
      }
    }
    if(nocounters == 0)
      fprintf(stderr, "%.1f%%     \r",(100.0*(i+1))/(double)(nalpha));
  }
}

//START REGION DEVELOP
