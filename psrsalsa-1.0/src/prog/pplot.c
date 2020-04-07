//START REGION RELEASE
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "psrsalsa.h"
//START REGION DEVELOP
#include <stdlib.h>

//START REGION RELEASE

#define max_nr_zap         10250 // Max nr of zapped subints/frequency channels that are supported
#define max_nr_stack       1000 // Max nr of states stored in the stack used in interactive mode
#define MaxNrPointsInPolygon 20 // Max nr of points put in a single polygon. Used to be 500, but this wasn't working for png's
//START REGION DEVELOP

#define max_nrlines          20 // Max nr of lines that can be specified on the command line
#define max_nrtext           20 // Max nr of strings that can be specified on the command line

//START REGION RELEASE

enum xUnitsSwitch_values {
  XUNIT_BINS,       // 0
  XUNIT_DEG,        // 1
  XUNIT_PHASE,      // 2
  XUNIT_TIME,       // 3
  XUNIT_ENDOFLIST   // 4
};

//START REGION DEVELOP
void SHOWREVISIONINFO_prog() {
#include "pplot.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);

//START REGION RELEASE
float ypos(float *stackI, int bin, long PulseNr, float scale, int NrBins, long pulse_bottom, long pulse_top, long subint_start, float yUnitCmdLine, float dyshift);
int setBaselineParams(datafile_definition fin, float *baseline, float *dxshift, int *xUnitsSwitch, verbose_definition verbose);
int zapVectors(int nrZappedVectors, int *zappedVectors, long subint_start, datafile_definition fin, float *dataSubset, int onlysinglechannel, verbose_definition verbose);
int zapSubints(int nrZappedSubints, int *zappedSubints, long subint_start, int didtranspose_orig_nrbin, datafile_definition fin, float *dataSubset, int onlysinglesubint, verbose_definition verbose);
// Type = 0 = Vector centre
//        1 = Bottom of vector
//        2 = Top of vector
// If inverse is set, the vector indicated in proper units is converted into a vector nr.
void yvec2unit(datafile_definition fin, int yUnitsSwitch, int vectorMinMaxSpecified, float vectorMin, float vectorMax, int didtranspose_orig_nrbin, int type, float valuein, float *valueout, int mapmode, int inverse, verbose_definition verbose);

typedef struct {
  long subint_start, subint_end;
  long x1, x2;
  int grayscalemode;     // -1 = undefined, 0=no map, 1=colour map, 2=contourmap, 3=both a map/contour plot
  int nrcontours;        // Number of contours in contour mode
  int nrZappedVectors, nrZappedSubints;
}plot_state_def;

void copystackstate(plot_state_def state1, plot_state_def *state2)
{
  state2->subint_start = state1.subint_start;
  state2->subint_end = state1.subint_end;
  state2->x1 = state1.x1;
  state2->x2 = state1.x2;
  state2->grayscalemode = state1.grayscalemode;
  state2->nrcontours = state1.nrcontours;
  state2->nrZappedVectors = state1.nrZappedVectors;
  state2->nrZappedSubints = state1.nrZappedSubints;
}

typedef struct {
  float characterheight;
  int linewidth;
  int font;
}pgplot_text_state_def;

int main(int argc, char **argv)
{
  long subint_start, subint_end, subint_range_defined;   // The subint range being read in (controlled with -prange command line option)
  int bin1_start, bin2_start;      // The bin-range to be plotted specified on command line
  float dxshift_start;             // The horizontal pulse longitude shift specified with the -dl option
  float dyshift;                   // Vertical subint nr shift
  int change_pulse_numbering_flag, change_pulse_numbering_value;  // Flag (and value) set if subint numbering is changed on command line
  int longitudeRangeSet_flag;            // Flag set if a pulse longitude range is set (either via the -long, -time or -phase options)
  float longitude_start, longitude_end;  // The corresponding longitude range
  int xUnitsSwitch;   /* Set to enum values defined in xUnitsSwitch_values */
  int yUnitsSwitch;   /* 0=whatever units are define on command line (-yunit) or bins, 1=predefined units if specified in file */
  float yUnitCmdLine;               // Scale the verical axis (including subint numbering)
  int vectorMinMaxSpecified;        // Set if value corresponding to first and last vector are set
  float vectorMin, vectorMax;       // Corresponding values
  int yUnitsTobs;                   // Set if tobs is requested rather than subint nr.
  int yUnitsMHz;                   // Set if frequency in MHz is requested rather than channel nr.
  int fixverticalscale_flag;   // Flag indicating that the plotted yrange (in line plot mode) shouldn't be adjusted to fit all pulses in plotted range
  int current_polnr;                     // The current polarization channel which is selected
  int interactive_flag;                  // Flag indicating if in interactive mode
  //  int nrZappedVectors, nrZappedSubints;  // The number of vectors/subints currently being zapped
  int *zappedVectors, *zappedSubints;    // The list of channel/subint numbers currently being zapped
  char *inputfilename;                   // The current input file name
  int viewportOptionsSet;                // Set if the viewport options (-dx, -ys, -x and -y) are defined
  float viewport_startx, viewport_starty, viewport_endx, viewport_endy;  // The corresponding specified viewport options
  float ymargin_both, ymargin_top;       // Set vertical margins at top+bottom and top only
  int nrpanelsx, nrpanelsy; // Nr of panels specified with the -N option
  int appendframes_flag;    // Flag is set if the panels next to each other should touch each other
  int nokeypresses_flag;        // Set if no key presses are desired after finishing a page of panels
  float scale_start, scale2_start;  // The scaling of the intensities to map to colour or subint nr
  plot_state_def *stack_state;  // The stack containing the current and past plot settings
  int current_stack_pos;        // The nr of states stored in stack
  int grayscalemode_start;      // Indicates if a color map, contour map, or a line plot should be generated
  int nrcontours_start;         // Indicates the number of contours
  int showTwice_flag;           // Flag indicating if plot is shown twice next to each other
  int nonumside_flag;           // Disable numbering of side panels flag
  int disable_x_numbers, disable_y_numbers; // Set if numbering x/y axis is disabled
  int showtop, showright;       // Flags indicating if side panels of colour maps are requested
  int showwedge;                // Flag indicating if a wedge at side of colour map should appear indicating color scale
  int histswitch;               // Flag indicating if a histogram is requested
  int polymode_flag;            // Flag set if want to use polygons rather than line drawings
  int noboxx, noboxy;           // Flag set if disable drawing x/y axis
  int heading_string_index;              // Points to the command line parameter which contains the string to be used as the heading
  pgplot_text_state_def heading_font;    // Sets the pgplot commands to set the font etc. for the heading text
  pgplot_frame_def_internal pgplot_frame;
  pgplot_text_state_def title_font;    // Sets the pgplot commands to set the font etc. for the title text
  pgplot_text_state_def label_font;    // Sets the pgplot commands to set the font etc. for the label text
  pgplot_text_state_def box_font;      // Sets the pgplot commands to set the font etc. for the box text
//START REGION DEVELOP
  pgplot_text_state_def copy_font;      // Sets the pgplot commands to set the font etc. for the box text
//START REGION RELEASE
  int plotlw;                          // Linewidth of the line drawings (not in map mode)
  int notitleset;                      // Set if no title is set on command line
  char title[1000];
  int xtitle_set, ytitleset;           // Set if title of x/y axis is set on command line
  char xtitle[1000], ytitle[1000];     // The titles themselves
  int wedgelabel_set;                  // Set if a label next to the wedge in a colour map is set
  char wedgelabel[1000];               // The label itself
  long i, j, k;
  int ignorebins, ignorebins2;
  int ignorelastpulses, scalerange;
  float xmax_oscil, xmin_oscil, scalerange_min, scalerange_max;
  psrsalsaApplication application;
  datafile_definition fin;

//START REGION DEVELOP
  int plotdot, plotdotonlyonce;
  int markerwidth, markercolor, markerwidth2, markercolor2;
  char copyright[100];
  float markerPosition, markerPosition2;
  float copyrightalignment, dycopyright;
  float errorbar_x, errorbar_y, errorbar_s;
  int errorbar_set;
  float line_x1[max_nrlines], line_y1[max_nrlines], line_x2[max_nrlines], line_y2[max_nrlines];
  int nrlines, line_lw[max_nrlines], line_ci[max_nrlines];
  float text_x[max_nrtext], text_y[max_nrtext], text_ch[max_nrtext];
  int nrtext, text_f[max_nrtext], text_lw[max_nrtext], text_ci[max_nrtext], text_txt[max_nrtext], outline, outline_color;

//START REGION RELEASE
  initApplication(&application, "pplot", "[options] inputfile(s)");
//START REGION DEVELOP
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
//START REGION RELEASE
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_rebin = 1;
  application.switch_device = 1;
  application.switch_tscr = 1;
  application.switch_tscr_complete = 1;
  application.switch_TSCR = 1;
  application.switch_fscr = 1;
  application.switch_FSCR = 1;
  application.switch_rot = 1;
  application.switch_rotdeg = 1;
  application.switch_debase = 1;
  application.switch_headerlist = 1;
  application.switch_header = 1;
  application.switch_itf = 1;
  application.switch_formatlist = 1;
  application.switch_iformat = 1;
  application.switch_noweights = 1;
  application.switch_useweights = 1;
  application.switch_uniformweights = 1;
  application.switch_polselect = 1;
  application.switch_nocounters = 1;
  application.switch_conshift= 1;
  application.switch_circshift= 1;
  application.switch_dedisperse= 1;
  application.switch_deFaraday= 1;
  application.switch_stokes = 1;
  application.switch_coherence = 1;
  application.switch_filelist = 1;
  application.switch_size = 1;
  application.switch_macro = 1;
  application.switch_cmaplist = 1;
  application.switch_cmap = 1;
  application.switch_onpulse = 1;
  application.switch_onpulsef = 1;
  application.switch_onpulsegr = 1;
  application.switch_deparang = 1;
  application.switch_changeRefFreq = 1;
  application.switch_norm = 1;
  application.switch_normglobal = 1;
  application.switch_clip = 1;
  application.switch_align = 1;
  application.switch_templatedata = 1;
  application.switch_template = 1;
  application.switch_libversions = 1;
//START REGION DEVELOP
  application.switch_dedeFaraday= 1;
  application.switch_dededisperse= 1;
  application.switch_insertparang = 1;
  /* Overwrite default to inverted color map, because the levels are
     inverted when calling pgplotMap. This is required to make the
     alternative itf's make the signal weaker instead of stronger. */
//START REGION RELEASE
  application.cmap = PPGPLOT_INVERTED_HEAT;
  /* Overide default plotdevice */
#ifdef HAVEPGPLOT2 
  strcpy(application.pgplotdevice, "/pw");
#else
  strcpy(application.pgplotdevice, "/xs");
#endif
//START REGION DEVELOP
  application.switch_noplotsubset = 1;   // Seems to do strange things in postscript mode
  application.switch_invertxy = 1;
  application.switch_gate = 1;
  application.switch_ppgplot = 1;
  application.switch_invertf = 1;
  application.switch_invertfp = 1;
  application.switch_invertxy = 1;
  application.switch_invertfx = 1;
  application.switch_normRMS = 1;
  application.switch_normRMSi = 1;
  application.switch_fft = 1;
  application.switch_testsignal = 1;
  application.switch_autot = 1;
  application.switch_centert = 1;
  application.switch_fftsmooth = 1;
  application.switch_template_onpulse = 1;
  application.switch_absweights = 1;
  application.switch_checknan = 1;
  application.switch_removenan = 1;
  application.switch_showrevision = 1;
  application.switch_autoonpulse = 1;
  application.switch_shuffle = 1;
  application.switch_rotslope = 1;
  application.switch_zerodming = 1;
  application.switch_zerodming_adv = 1;
//START REGION RELEASE
  application.switch_forceUniformFreqLabelling = 1;

  interactive_flag = 0;
  dxshift_start = 0;
  dyshift = 0;
  xUnitsSwitch = XUNIT_BINS;    /* Default is showing x-axis in bins */
  yUnitsSwitch = 0;             /* Default is showing y-axis in vectors */
  yUnitCmdLine = 1;
  longitudeRangeSet_flag = 0;
  current_polnr = 0;
  subint_range_defined = 0;
  heading_font.characterheight = 1;
  heading_font.linewidth = 1;
  heading_font.font = 1;
  heading_string_index = 0;
  title_font.characterheight = 1;
  title_font.linewidth = 1;
  title_font.font = 1;
  box_font.characterheight = 1;
  box_font.linewidth = 1;
  box_font.font = 1;
  label_font.characterheight = 1;
  label_font.linewidth = 1;
  label_font.font = 1;
//START REGION DEVELOP
  copy_font.characterheight = 1.3;
  copy_font.linewidth = 1;
  copy_font.font = 3;
//START REGION RELEASE
  notitleset = 1;
  title[0] = 0;
  xtitle_set = 0;
  ytitleset = 0;
  wedgelabel_set = 0;
  viewportOptionsSet = 0;
  viewport_startx = 0.15;
  viewport_starty = 0.15;
  viewport_endx = 0.9;
  viewport_endy = 0.9;
  ymargin_both = 0;
  ymargin_top = 0;
  appendframes_flag = 0;
  showTwice_flag = 0;
  fixverticalscale_flag = 0;
  nokeypresses_flag = 0;
  grayscalemode_start = -1;   // -1 = undefined on command line
  nrcontours_start = 0;
  scale_start = 1;
  scale2_start = 0;
  nrpanelsx = 1;
  nrpanelsy = 1;
  vectorMinMaxSpecified = 0;
  bin1_start = -1;
  bin2_start = -1;
  polymode_flag = 1;
  ignorebins = 0;
  ignorebins2 = 0;
  showtop = 0;
  showright = 0;
  showwedge = 0;
  change_pulse_numbering_flag = 0;
  change_pulse_numbering_value = 0;
  plotlw = 1;
  ignorelastpulses = 0;
  histswitch = 0;
  xmax_oscil = 1e10;
  xmin_oscil = -1e10;
  noboxx = noboxy = 0;
  disable_x_numbers = 0;
  disable_y_numbers = 0;
  nonumside_flag = 0;
  yUnitsTobs = 0;
  yUnitsMHz = 0;
  clear_pgplot_frame(&pgplot_frame);
  scalerange = 0;
//START REGION DEVELOP
  nrlines = 0;
  nrtext = 0;
  outline = -1;
  markerPosition = -1;
  markerPosition2 = -1;
  copyright[0] = 0;
  copyrightalignment = 1.0;
  dycopyright = 0;
  errorbar_set = 0;
  plotdot = 0;
  plotdotonlyonce = 0;
//START REGION DEVELOP
//START REGION RELEASE

  if(argc < 2) {
    printf("Program to plot pulsar data in various ways.\n");
    printApplicationHelp(&application);
    printf("Other options:\n");
    printf("-interactive    Turn interactive mode on.\n");
    printf("-ia             short for -interactive.\n");
    
    printf("\nData selection options:\n");
    printf("-prange     Specify subint range to be plotted, default is all.\n");
    printf("-b          Specify bin range to be plotted, default is all.\n");
    //START REGION DEVELOP
    /*     printf("-B          Force baseline to be this number of degrees.\n"); Superseded by -header option */  
    //START REGION RELEASE
    printf("-long       Specify longitude range in degrees.\n");
    printf("-time       Specify longitude range in seconds.\n");
    printf("-phase      Specify longitude range in phase.\n");
    //START REGION DEVELOP
    printf("-xmax       Specify max time to be plotted in seconds (allows oscillescope \n");
    printf("              effect, doesn't work with map, affects last subint only).\n");
    printf("-xmin       Specify min time to be plotted in seconds (allows oscillescope \n");
    printf("              effect, doesn't work with map, affects last subint only).\n");
    //START REGION RELEASE
    printf("\nPanel options:\n");
    printf("-appendframes Remove space between the frames of multiple adjacent panels\n");
    printf("-N            \"nrx nry\" Create nrx by nry panels, rather than a single plot\n");
    printf("              per page\n");
    
    printf("\nPlot options:\n");
    //START REGION DEVELOP
    printf("-cont n       Show graphics as a contour map with n contours (can use -lw and -map).\n");
    printf("-dot          Plot a dot at value specified by -xmax with this linewidth.\n");
    printf("-error        \"x y sigma\" Add vertical errorbar.\n");
    //START REGION RELEASE
    printf("-hist         Turn on histogram mode. The bin centre is plotted at an integer\n");
    printf("              bin nr. If plotting versus pulse longitude/time/phase, you can\n");
    printf("              use the -dl option.\n");
    //START REGION DEVELOP
    printf("-marker       \"pulse linewidth color\" plot marker at this position.\n");
    printf("-line         \"x1 y1 x2 y2 lw color\" plot line at this position.\n");
    //START REGION RELEASE
    printf("-lw           set line width of line plot (not used in map mode).\n");
    printf("-lp           Show graphics as a \"Joy Division\" line plot rather than a colour\n");
    printf("              map (-map). For a single row of data -lp is the default.\n");
    printf("-map          Show graphics as a colour map (default if multiple rows in data).\n");
    printf("              The bin centre is plotted at an integer bin nr. If plotting versus\n");
    printf("              pulse longitude/time/phase, you can use the -dl option.\n");
    printf("-nokeypress   Do not wait for key presses to go to next page in output plot\n");
    printf("-nopoly       Do not use filled polygons but only line drawing, which sometimes\n");
    printf("              could give better results (or worse).\n");
    //START REGION DEVELOP
    printf("-outline      \"lw color\", which affects the -line and -text options\n");
    //START REGION RELEASE
    printf("-showwedge    Plot an annotated wedge to show color scale (use with -map).\n");
    printf("-showtop      Show a panel at top of the map (use with -map).\n");
    printf("-showright    Show a panel at right of the map (use with -map).\n");
    printf("-showtwice    Plot the map twice above each other (use with -map).\n");
    //START REGION DEVELOP
    printf("-text         \"x y ch lw font color\" \"text\" plot text at this position.\n");
    printf("              (keywords listed with the -textkeywords are supported).\n");
    //START REGION RELEASE
    printf("-textkeywords Lists of keywords you can use with");
    //START REGION DEVELOP
    printf(" -text and");
    //START REGION RELEASE
    printf(" -title.\n");
    //START REGION DEVELOP
    printf("-vmarker      \"x linewidth color\" plot vertical marker at this position.\n");
    //START REGION RELEASE
    printf("-dx           Set viewport x-range (start) to this value (default is %.2f).\n", viewport_startx);
    printf("-x            Set viewport x-range (end) to this value (default is %.2f).\n", viewport_endx);
    printf("-ys           Set viewport y-range (start) to this value (default is %.2f).\n", viewport_starty);
    printf("-y            Set viewport y-range (end) to this value (default is %.2f).\n", viewport_endy);
  
    //START REGION RELEASE
    printf("\nScaling and offset of data:\n");
    printf("-scale      Specify scale, default is 1. This option multiplies the intensity by\n");
    printf("            this factor. In -map mode, a value > 1 results in clipping, making\n");
    printf("            weak features clearer.\n");
    printf("-scale2     Specify scale2 (a value between 0 and 1), default is 0, only used in\n");
    printf("            map mode. This clips the low intensity values to emphasize the\n");
    printf("            bright features.\n");
    printf("-scalerange \"min max\" In map mode: The used color range is between min and max.\n");
    printf("-dl         Shift plot in pulse longitude.\n");
    //START REGION DEVELOP
    printf("-fix        Fix range of vertical axis in line plot (i.e. assume pulses are\n");
    printf("            within range 0-1, which can be modified with the -ym and\n");
    printf("            -ym2 options).\n");
    //START REGION RELEASE
    printf("-half       Shift pulses by half a pulse.\n");
    printf("-0          Force first pulse to be plotted as pulse zero (ignored with -map).\n");
    printf("-1          Force first pulse to be plotted as pulse one (ignored with -map).\n");
    //START REGION DEVELOP
    printf("-ip         Ignore this number of pulses to be plotted at end in line plot mode.\n");
    printf("-ib         Ignore this number of longitude bins to be plotted.\n");
    printf("-ib2        Ignore this extra number of bins at the right to be plotted.\n");
    printf("-ym         Set ymargins (add this amount to both sides of the vertical axis).\n");
    printf("-ym2        Set extra ymargin on top of plot (add this amount to only top side\n");
    printf("            of the vertical axis).\n");
    //START REGION RELEASE
    printf("-yunit      Multiply the y-axis with this number\n");
    printf("-vecRange   Specify value corresponding to first and last vector (vertical axis).\n");
    printf("-tobs       For subint/pulse phase plots, use observing time as y-axis.\n");
    printf("-MHz        For frequency/pulse phase plots, use MHz as y-axis.\n");
  
    printf("\nLabel options:\n");
    printf("-title      \"...\" Specify title (supports keywords listed with -textkeywords)\n");
    //START REGION RELEASE
    printf("-heading    \"...\" Specify heading.\n");
    //START REGION DEVELOP
    printf("-copy       \"...\" Specify copyright text.\n");
    printf("-copya            Specify copyright text alignment (default is %.1f).\n", copyrightalignment);
    printf("-copydy           Specify copyright dy (default is %.1f).\n", dycopyright);
    //START REGION RELEASE
    printf("-xtitle     \"...\" Set title of xaxis.\n");
    printf("-ytitle     \"...\" Set title of yaxis.\n");
    printf("-wedgetitle \"...\" Set title of wedge generated with -showwedge.\n");
    printf("-nonumx           Disable numbering x-as.\n");
    printf("-nonumy           Disable numbering y-as.\n");
    printf("-nonumside        Disable numbering on the side panels.\n");
    printf("-noboxx           Disable drawing the x-axis.\n");
    printf("-noboxy           Disable drawing the y-axis.\n");
    printf("-labels           \"heading_ch heading_lw heading_f title_ch title_lw title_f label_ch label_lw label_f box_ch box_lw box_f");
//START REGION DEVELOP
    printf(" copy_ch copy_lw copy_f");
    //START REGION RELEASE
    printf("\".\n");
    printf("                  default is \"%.1f %d %d %.1f %d %d %.1f %d %d %.1f %d %d", heading_font.characterheight, heading_font.linewidth, heading_font.font, title_font.characterheight, title_font.linewidth, title_font.font, label_font.characterheight, label_font.linewidth, label_font.font, box_font.characterheight, box_font.linewidth, box_font.font);
//START REGION DEVELOP
    printf(" %.1f %d %d", copy_font.characterheight, copy_font.linewidth, copy_font.font);
    //START REGION RELEASE
    printf("\".\n");
    printf("\n");
    printCitationInfo();
    terminateApplication(&application);
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      int index;
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-nokeypress") == 0) {
	nokeypresses_flag = 1;
      }else if(strcmp(argv[i], "-prange") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld %ld", &subint_start, &subint_end, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	if(subint_end < subint_start) {
	  printerror(application.verbose_state.debug, "ERROR pplot: The selected subint range seems to be such that the first subint exceeds the final.");
	  return 0;
	}
	subint_range_defined = 1;
        i++;
      }else if(strcmp(argv[i], "-appendframes") == 0) {
	appendframes_flag = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-ip") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &ignorelastpulses, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-ib") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &ignorebins, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-ib2") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &ignorebins2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-lw") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &plotlw, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-N") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &nrpanelsx, &nrpanelsy, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-dl") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &dxshift_start, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-dot") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &plotdot, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-xmax") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &xmax_oscil, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	polymode_flag = 0;
	printwarning(application.verbose_state.debug, "WARNING pplot: enabled the -nopoly flag to allow the -xmax option to work");
	i++;
      }else if(strcmp(argv[i], "-xmin") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &xmin_oscil, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-vecRange") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &vectorMin, &vectorMax, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	vectorMinMaxSpecified = 1;
	yUnitsSwitch = 1;
        i++;
      }else if(strcmp(argv[i], "-tobs") == 0) {
	yUnitsTobs = 1;	
      }else if(strcmp(argv[i], "-MHz") == 0) {
	yUnitsMHz = 1;	
      }else if(strcmp(argv[i], "-b") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &bin1_start, &bin2_start, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-marker") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %d %d", &markerPosition, &markerwidth, &markercolor, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-vmarker") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %d %d", &markerPosition2, &markerwidth2, &markercolor2, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-line") == 0) {
	if(nrlines < max_nrlines) {
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f %f %d %d", &line_x1[nrlines], &line_y1[nrlines], &line_x2[nrlines], &line_y2[nrlines], &line_lw[nrlines], &line_ci[nrlines], NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	  nrlines++;
	}
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-textkeywords") == 0) {
	str_list_replace_keys(0);
	terminateApplication(&application);
	return 0;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-text") == 0) {
	if(nrtext < max_nrtext) {
	  if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f %d %d %d", &text_x[nrtext], &text_y[nrtext], &text_ch[nrtext], &text_lw[nrtext], &text_f[nrtext], &text_ci[nrtext], NULL) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	    return 0;
	  }
	  text_txt[nrtext] = i+2;
	  nrtext++;
	}else {
	  printerror(application.verbose_state.debug, "pplot: Too much text specified");
	  return 0;
	}
	i += 2;
      }else if(strcmp(argv[i], "-outline") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d", &outline, &outline_color, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-yunit") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &yUnitCmdLine, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-ytitle") == 0) {
	strcpy(ytitle, argv[i+1]);
	ytitleset = 1;
	i++;
      }else if(strcmp(argv[i], "-xtitle") == 0) {
	strcpy(xtitle, argv[i+1]);
	xtitle_set = 1;
	i++;
      }else if(strcmp(argv[i], "-wedgetitle") == 0) {
	strcpy(wedgelabel, argv[i+1]);
	wedgelabel_set = 1;
	i++;
      }else if(strcmp(argv[i], "-nonumside") == 0) {
	nonumside_flag = 1;
      }else if(strcmp(argv[i], "-nonumx") == 0) {
	disable_x_numbers = 1;
      }else if(strcmp(argv[i], "-nonumy") == 0) {
	disable_y_numbers = 1;
      }else if(strcmp(argv[i], "-showtwice") == 0) {
	showTwice_flag = 1;
      }else if(strcmp(argv[i], "-labels") == 0) {
//START REGION DEVELOP
#if 0
//START REGION RELEASE
// Not executed in develop mode: ask for 3 parameters less (last copyright-related variables)
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %d %d %f %d %d %f %d %d %f %d %d", &heading_font.characterheight, &heading_font.linewidth, &heading_font.font, &title_font.characterheight, &title_font.linewidth, &title_font.font, &label_font.characterheight, &label_font.linewidth, &label_font.font, &box_font.characterheight, &box_font.linewidth, &box_font.font, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
//START REGION DEVELOP
#else
	// Executed when in develop mode, but not in release mode. Has the three extra parameters at the end
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %d %d %f %d %d %f %d %d %f %d %d %f %d %d", &heading_font.characterheight, &heading_font.linewidth, &heading_font.font, &title_font.characterheight, &title_font.linewidth, &title_font.font, &label_font.characterheight, &label_font.linewidth, &label_font.font, &box_font.characterheight, &box_font.linewidth, &box_font.font, &copy_font.characterheight, &copy_font.linewidth, &copy_font.font, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
#endif
//START REGION RELEASE
	/*
	j = sscanf(argv[i+1], "%f %d %d %f %d %d %f %d %d %f %d %d %f %d %d", &heading_font.characterheight, &heading_font.linewidth, &heading_font.font, &title_font.characterheight, &title_font.linewidth, &title_font.font, &label_font.characterheight, &label_font.linewidth, &label_font.font, &box_font.characterheight, &box_font.linewidth, &box_font.font, &copy_font.characterheight, &copy_font.linewidth, &copy_font.font);
	if(j != 12 
//START REGION DEVELOP
	   + 3
//START REGION RELEASE
	   ) {
	  printerror(application.verbose_state.debug, "pplot: Error parsing option %s", argv[i]);
	  return 0;
	}
	*/
	i++;
      }else if(strcmp(argv[i], "-scale") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &scale_start, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-scale2") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &scale2_start, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-scalerange") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &scalerange_min, &scalerange_max, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	scalerange = 1;
	i++;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-copy") == 0) {
	sprintf(copyright, "\\(0274) ");
	strcat(copyright, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-copya") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &copyrightalignment, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-copydy") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &dycopyright, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-ym") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &ymargin_both, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-ym2") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &ymargin_top, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-title") == 0) {
	strcpy(title, argv[i+1]);
	notitleset = 0;
	i++;
      }else if(strcmp(argv[i], "-heading") == 0) {
	heading_string_index = i+1;
	i++;
      }else if(strcmp(argv[i], "-x") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &viewport_endx, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	viewportOptionsSet = 1;
	i++;
      }else if(strcmp(argv[i], "-dx") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &viewport_startx, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	viewportOptionsSet = 1;
	i++;
      }else if(strcmp(argv[i], "-ys") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &viewport_starty, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	viewportOptionsSet = 1;
	i++;
      }else if(strcmp(argv[i], "-y") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &viewport_endy, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	viewportOptionsSet = 1;
	i++;
      }else if(strcmp(argv[i], "-long") == 0 || strcmp(argv[i], "-time") == 0 || strcmp(argv[i], "-phase") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &longitude_start, &longitude_end, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	xUnitsSwitch = XUNIT_DEG;
	if(strcmp(argv[i], "-phase") == 0) {
	  xUnitsSwitch = XUNIT_PHASE;
	}else if(strcmp(argv[i], "-time") == 0) {
	  xUnitsSwitch = XUNIT_TIME;
	}
	longitudeRangeSet_flag = 1;
	i++;
      }else if(strcmp(argv[i], "-0") == 0) {
	change_pulse_numbering_flag = 1;
	change_pulse_numbering_value = 0;
      }else if(strcmp(argv[i], "-1") == 0) {
	change_pulse_numbering_flag = 1;
	change_pulse_numbering_value = 1;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-fix") == 0) {
	fixverticalscale_flag = 1;
//START REGION RELEASE
      }else if(strcmp(argv[i], "-half") == 0) {
	dyshift = 0.5;
      }else if(strcmp(argv[i], "-hist") == 0) {
	histswitch = 1;
	grayscalemode_start = 0;
	//	polymode_flag = 0;  // Histogram mode only supported when not using polygons
      }else if(strcmp(argv[i], "-noboxx") == 0) {
	noboxx = 1;
      }else if(strcmp(argv[i], "-noboxy") == 0) {
	noboxy = 1;
      }else if(strcmp(argv[i], "-lp") == 0) {
	grayscalemode_start = 0;
      }else if(strcmp(argv[i], "-map") == 0) {
	if(grayscalemode_start == 2) {   // If -cont is specified, do both
	  grayscalemode_start = 3;
	}else {
	  grayscalemode_start = 1;
	}
      }else if(strcmp(argv[i], "-cont") == 0) {
	if(grayscalemode_start == 1) {   // If -map is specified, do both
	  grayscalemode_start = 3;
	}else {
	  grayscalemode_start = 2;
	}
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &nrcontours_start, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-showtop") == 0) {
	showtop = 1;
      }else if(strcmp(argv[i], "-showright") == 0) {
	showright = 1;
      }else if(strcmp(argv[i], "-showwedge") == 0) {
	showwedge = 1;
      }else if(strcmp(argv[i], "-interactive") == 0 || strcmp(argv[i], "-ia") == 0) {
	interactive_flag = 1;
      }else if(strcmp(argv[i], "-nopoly") == 0) {
	polymode_flag = 0;
//START REGION DEVELOP
      }else if(strcmp(argv[i], "-error") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f", &errorbar_x, &errorbar_y, &errorbar_s, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	errorbar_set = 1;
        i++;
//START REGION RELEASE
      }else {
	/* If the option is not recognized, assume it is a filename */
	if(argv[i][0] == '-') {
	  printerror(application.verbose_state.debug, "pplot: Unknown option: %s\n\nRun pplot without command line arguments to show help", argv[i]);
	  terminateApplication(&application);
	  return 0;
	}else {
	  if(applicationAddFilename(i, application.verbose_state) == 0)
	    return 0;
	}
      }
    }
  }

  if(histswitch && grayscalemode_start) {
    printerror(application.verbose_state.debug, "ERROR pplot: Cannot use the -hist and -map/-cont options simultaneously.");
    return 0;
  }

  if(applicationFilenameList_checkConsecutive(argv, application.verbose_state) == 0) {
    return 0;
  }

  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR pplot: No files specified");
    return 0;
  }

  if(numberInApplicationFilenameList(&application, argv, application.verbose_state) > 1 && interactive_flag) {
    printerror(application.verbose_state.debug, "ERROR pplot: Cannot plot multiple files in interactive mode.");
    return 0;
  }

  if(application.macro_ptr != NULL) {
    interactive_flag = 1;  // Interactive mode is implied when using macro    
  }

  //START REGION DEVELOP
  // Overwrite default viewport placement options to different values if multiple panels are selected
  //START REGION RELEASE
  if(viewportOptionsSet == 0) {
    if(nrpanelsx > 1 || nrpanelsy > 1) {
      viewport_endx = 0.95;
      viewport_startx = 0.05;
      viewport_starty = 0.1;
      viewport_endy = 0.95;
    }
  }

//START REGION DEVELOP
  if(application.ppgplot_name != NULL) {
    if(pgopenoutputfile(application.ppgplot_name) == 0) {
      printerror(application.verbose_state.debug, "ERROR pplot: Cannot open %s", application.ppgplot_name);
      return 0;
    }
  }

//START REGION RELEASE
  ppgopen(application.pgplotdevice);
  pgplot_setWindowsize(application.windowwidth, application.windowheight, -1);
  
  ppgask(0);
  ppgslw(1);

  // Allocate memory to store zapped vectors/subints
  zappedVectors = malloc(max_nr_zap*sizeof(int));
  zappedSubints = malloc(max_nr_zap*sizeof(int));
  // Allocate memory to store stack state
  stack_state = malloc(max_nr_stack*sizeof(plot_state_def));
  if(zappedVectors == NULL || zappedSubints == NULL || stack_state == NULL) {
    printerror(application.verbose_state.debug, "ERROR pplot: Memory allocation error.");
    return 0;
  }

  int curpanelnrx, curpanelnry;   // The current panel nr being used
  curpanelnrx = 0;
  curpanelnry = 0;
  int goingToPlotFirstPanel;                  // Is set if nothing is plotted yet
  goingToPlotFirstPanel = 1;
  int dataSubset_allocated;                   // Keeps track if memory is already allocated
  dataSubset_allocated = 0;
  int stack_poly_allocated;                   // Keeps track if memory is already allocated
  stack_poly_allocated = 0;
  int maxy_allocated;                        // Keeps track if memory is already allocated
  maxy_allocated = 0;

  //  cleanPSRData(&clone, application.verbose_state);  // Initialise clone, but release all associated memory
  //  closePSRData(&clone, 0, application.verbose_state);

  /* If there are files at the end of the command line, then process them one by one. */
  while((inputfilename = getNextFilenameFromList(&application, argv, application.verbose_state)) != NULL) {
    int data_read; // Flag indicated if data is read in (0=no, 1=yes, 2=read different polarization, but leave most variables as they are)
    data_read = 0;
    int nrpolarizations; // This nr will be set to the nr of polarizations in the current file BEFORE only the current polarization channel is selected
    nrpolarizations = 0;
    int didtranspose_orig_nrbin;    // This variable will be set if the input data needs a transpose allowing to make a freq/subint plot. It contains the nr of bins in the data (before transpose, so not nrsubints * nrbins
    didtranspose_orig_nrbin = 0;
    float baseline;      // Will be set to the range covered by the bins in a subint.
    float viewport_startxCurrentPanel; // Horizontal start position the current pgplot viewport
    float viewport_endxCurrentPanel;   // Horizontal end position the current pgplot viewport
    float scale, scale2;               // The used scaling of the data
    long *maxy;
    float dxshift;    // The horizontal shift taken from command-line or from file
    int bin1, bin2;   // Bin range to be plotted
    int zapmode;      // 0=zapping of vectors, 1=zapping of subints in subint/freq mode
    int nx_last_clicked_zap, ny_last_clicked_zap;  // Used when selecting a range of subints with the right mouse button
    nx_last_clicked_zap = 0;
    ny_last_clicked_zap = 0;

    float *dataSubset;    // This array contains the subint range being plotted. When plotting colour maps, the data which is plotted comes directly from variable fin.
    float *xstack_poly, *ystack_poly;
    float dummyf1, dummyf2, dummyf3, dummyf4;
    int dummyi;
    char txt1[1000], txt2[1000], txt3[1000];

    maxy = NULL;        // To avoid compiler warning
    dataSubset = NULL;  // To avoid compiler warning
    xstack_poly = ystack_poly = NULL; // To avoid compiler warning
    dxshift = dxshift_start;
    zapmode = 0;

    // Reset selected subint range if going to the next file if range was not defined on command line
    if(subint_range_defined == 0) {
      subint_start = -1;
    }
    if(goingToPlotFirstPanel == 0) {
      if(curpanelnrx == 0 && curpanelnry == 0) {
	if(nokeypresses_flag == 0) {
	  printf("Press a key to continue\n");
	  fflush(stdout);
	  pgetch_macro(&application, application.verbose_state);
	}
      }
    }
    goingToPlotFirstPanel = 0;

    /* Empty stack by setting current position back to zero */
    current_stack_pos = 0;

    printf("Plotting: %s\n", inputfilename);

    //    nrZappedVectors = 0;
    //    nrZappedSubints = 0;
    stack_state[0].nrZappedVectors = 0;
    stack_state[0].nrZappedSubints = 0;


    cleanPSRData(&fin, application.verbose_state);   // Not necessary, since already done when reading in data, but useful so it can be closed before reading in the file, which is useful in interactive mode where the data can be read in multiple times when switching polarization.

    bin1 = bin1_start;
    bin2 = bin2_start;
    scale = scale_start;    // Set the intensity scaling to the command-line values
    scale2 = scale2_start;
    stack_state[0].grayscalemode = grayscalemode_start;
    stack_state[1].grayscalemode = grayscalemode_start;  // Not entirely sure why required, but gets rid of a valgrind error.
    stack_state[0].nrcontours = nrcontours_start;
    stack_state[1].nrcontours = nrcontours_start;

    /* If no title is specified, use filename instead */
    if(notitleset) {
      strcpy(title, inputfilename);
    }

//START REGION DEVELOP
    // After removing -i option, this should have become obsolete, as without -i option used, NrFiles was set to 1
    //    viewport_endx /= (float)NrFiles;   
//START REGION RELEASE
    int firstPulsedataSubset;
    int ytitle_set_to_intensity, ytitle_set_to_pulsenumber, ytitle_set_to_subint, ytitle_set_to_freqchannel, ytitle_set_to_fluctuationFreq, ytitle_set_to_padist, ytitle_set_to_elldist, ytitle_set_to_pulselongitude, ytitle_set_to_p3fold, ytitle_set_to_spectral_power, ytitle_set_to_harmonic_number, ytitle_set_to_lag_number;
    ytitle_set_to_intensity = 0;
    ytitle_set_to_pulsenumber = 0;
    ytitle_set_to_subint = 0;
    ytitle_set_to_freqchannel = 0;
    ytitle_set_to_fluctuationFreq = 0;
    ytitle_set_to_padist = 0;
    ytitle_set_to_elldist = 0;
    ytitle_set_to_pulselongitude = 0;
    ytitle_set_to_p3fold = 0;
    ytitle_set_to_spectral_power = 0;
    ytitle_set_to_harmonic_number = 0;
    ytitle_set_to_lag_number = 0;
    firstPulsedataSubset = 0; // To avoid compile warnings
    do { // Corresponds with while(interactive_flag);
      float currentPanelScaling, currentPanelScalingx, currentPanelScalingy;
      float viewport_startyCurrentPanel, viewport_endyCurrentPanel;
      float xmin, xmax, ymin, ymax, xminshow, xmaxshow;
      float x, x2, y, y2;
      int redraw;
      char *newtext, ch;
      x = x2 = y = y2 = ch = 0; 
//START REGION DEVELOP
      // WORK OUT THE VIEWPORT RANGE TO USE FOR PLOT
      // viewport available to cover with all plots: viewport_endx-viewport_startx
      // This is split between nrpanelsx panels and (nrpanelsx-1) separations. The separations are a left plus right margin.
      // The separation (relative to 1) is viewport_startx and the panel size is (viewport_endx-viewport_startx) and the other margin is (1-viewport_endx)
      // The total plot viewport range necessary (if not scaling things) would therefore be:
      // (nrpanelsx-1)*viewport_startx + nrpanelsx*(viewport_endx-viewport_startx) + (nrpanelsx-1)*(1-viewport_endx) =
      //              -viewport_startx + nrpanelsx*viewport_endx                  + (nrpanelsx-1)*(1-viewport_endx) =
      //              -viewport_startx + nrpanelsx*viewport_endx                  + nrpanelsx*(1-viewport_endx) - (1-viewport_endx) =
      //              -viewport_startx                                            + nrpanelsx                   - (1-viewport_endx) =
      //              -viewport_startx                                            + nrpanelsx                   -1+viewport_endx =
      // viewport_endx-viewport_startx + nrpanelsx -1
      // This means that all separations/sizes need to be scaled with the available space (viewport_endx-viewport_startx) divided by this number
//START REGION RELEASE

      // Default, if there are not multiple panels
      currentPanelScaling = 1;
      // Start of left margin
      viewport_startxCurrentPanel = viewport_startx;
      if(nrpanelsx > 1) {
	if(appendframes_flag == 0)
	  currentPanelScalingx = (viewport_endx-viewport_startx)/(viewport_endx-viewport_startx + nrpanelsx -1);
	else
	  currentPanelScalingx = 1.0/nrpanelsx;
	currentPanelScaling = currentPanelScalingx;
	//	  printf("XXXX scaling = %f\n", currentPanelScalingx);
	// Add left margins
	if(appendframes_flag == 0 || curpanelnrx == 0)
	  viewport_startxCurrentPanel += curpanelnrx*currentPanelScalingx*(viewport_startx);
	// Add panels
	viewport_startxCurrentPanel += curpanelnrx*currentPanelScalingx*(viewport_endx-viewport_startx);
	// Add right margins
	if(appendframes_flag == 0)
	  viewport_startxCurrentPanel += curpanelnrx*currentPanelScalingx*(1-viewport_endx);
      }
      // End of plot
      viewport_endxCurrentPanel = viewport_endx;
      if(nrpanelsx > 1) {
	// Start of plot
	viewport_endxCurrentPanel = viewport_startxCurrentPanel;
	// Add size of plot
	viewport_endxCurrentPanel += currentPanelScalingx*(viewport_endx-viewport_startx);
      }
      //	printf("XXXX vp_start=%f vp_end=%f\n", viewport_startxCurrentPanel, viewport_endxCurrentPanel);
      //	printf("XXXX vp_start=%f vp_end=%f (settings)\n", viewport_startx, viewport_endx);
      
      // Start of bottom margin
      viewport_startyCurrentPanel = viewport_starty;
      if(nrpanelsy > 1) {
	currentPanelScalingy = (viewport_endy-viewport_starty)/(viewport_endy-viewport_starty + nrpanelsy -1);
	if(currentPanelScalingy < currentPanelScaling)
	  currentPanelScaling = currentPanelScalingy;
	//	  printf("XXXX scaling = %f\n", currentPanelScalingy);
	// Add bottom margins
	viewport_startyCurrentPanel += (nrpanelsy-curpanelnry-1)*currentPanelScalingy*(viewport_starty);
	// Add panels
	viewport_startyCurrentPanel += (nrpanelsy-curpanelnry-1)*currentPanelScalingy*(viewport_endy-viewport_starty);
	// Add top margins
	viewport_startyCurrentPanel += (nrpanelsy-curpanelnry-1)*currentPanelScalingy*(1-viewport_endy);
      }
      // End of plot
      viewport_endyCurrentPanel = viewport_endy;
      if(nrpanelsy > 1) {
	// Start of plot
	viewport_endyCurrentPanel = viewport_startyCurrentPanel;
	// Add size of plot
	viewport_endyCurrentPanel += currentPanelScalingy*(viewport_endy-viewport_starty);
      }
      //	printf("XXXX vp_start=%f vp_end=%f\n", viewport_startyCurrentPanel, viewport_endyCurrentPanel);
      //	printf("XXXX vp_start=%f vp_end=%f (settings)\n", viewport_starty, viewport_endy);
    
      if(curpanelnrx == 0 && curpanelnry == 0)
	ppgpage();
      
      if(data_read == 0 || data_read == 2) {  // If data is not yet read in (or reading in other polarization channel)
	int iformat;
	iformat = application.iformat;
	if(access(inputfilename, F_OK) != 0) {
	  fflush(stdout);
	  printerror(application.verbose_state.debug, "ERROR pplot: Cannot open file %s for reading.", inputfilename);
	  free(zappedVectors);
	  free(zappedSubints);
	  free(stack_state);
	  closePSRData(&fin, 0, application.verbose_state);
	  terminateApplication(&application);
	  return 0;	
	}
	if(iformat <= 0) {
	  iformat = guessPSRData_format(inputfilename, 0, application.verbose_state);
	  if(iformat == -2 || iformat == -3)
	    return 0;
	}
	if(isValidPSRDATA_format(iformat) == 0) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Input file cannot be opened. Please check if file %s exists and otherwise specify the correct input format with the -iformat option if the format is supported, but not automatically recognized.\n\n", inputfilename);
	  free(zappedVectors);
	  free(zappedSubints);
	  free(stack_state);
	  closePSRData(&fin, 0, application.verbose_state);
	  terminateApplication(&application);
	  return 0;
	}
	closePSRData(&fin, 0, application.verbose_state);  // Close/release memory if already loaded before in interactive mode (when switching polarizations).
	if(!openPSRData(&fin, inputfilename, iformat, 0, 1, 0, application.verbose_state)) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Error opening file.\n");
	  return 0;
	}

	/* Search commandline for header parameters to overwrite header
	   parameters. */
	if(PSRDataHeader_parse_commandline(&fin, argc, argv, application.verbose_state) == 0)
	  return 0;
	for(i = 1; i < argc; i++) {
	  if(strcmp(argv[i], "-header") == 0) {
	    printwarning(application.verbose_state.debug, "WARNING pplot: If using the -header option, be aware it applied BEFORE the preprocessing.");
	    break;
	  }
	}
	if(preprocessApplication(&application, &fin) == 0) {
	  return 0;
	}

	if(fin.gentype == GENTYPE_PENERGY) {
	  if(application.verbose_state.verbose) {
	    printwarning(application.verbose_state.debug, "The file appears to be penergy output. The \"polarization channels\" correspond to peak intensity (on & off-pulse), integrated energy (on & off-pulse), rms (on & off-pulse) and the S/N of the on-pulse region. The \"subintegrations\" correspond to the different polarizations in the original file.\n");
	  }
	}
	if((fin.poltype == POLTYPE_ILVPAdPA || fin.poltype == POLTYPE_ILVPAdPATEldEl) && fin.NrSubints == 1 && fin.NrFreqChan == 1) {
	  if(application.verbose_state.verbose) {
	    printwarning(application.verbose_state.debug, "The file appears to be ppol output. Plotting this data with ppol allows the profile and PA-swing to be shown in a single figure.\n");
	  }
	}

	if(fin.isFolded && fin.foldMode == FOLDMODE_FIXEDPERIOD) {
	  if(fin.fixedPeriod <= 0.0) {
	    printwarning(application.verbose_state.debug, "WARNING pplot: Period does not appear to be set, assuming it is 1 sec.");
	    fin.fixedPeriod = 1.0;
	  }
	}
	if(fin.tsampMode == TSAMPMODE_LONGITUDELIST) {
	  if(convert_to_fixed_tsamp(&fin, application.verbose_state) != 1) {
	    printerror(application.verbose_state.debug, "ERROR pplot: Cannot convert sampling to a regular grid.");
	    return 0;
	  }
	}
	if(fin.isFolded && fin.tsampMode == TSAMPMODE_FIXEDTSAMP) {
	  if(get_tsamp(fin, 0, application.verbose_state) <= 0.0) {
	    printwarning(application.verbose_state.debug, "WARNING pplot: Assuming full period is stored.");
	    double period;
	    int ret;
	    ret = get_period(fin, 0, &period, application.verbose_state);
	    if(ret == 2) {
	      printerror(application.verbose_state.debug, "ERROR pplot (%s): Cannot obtain period", fin.filename);
	      return 0;
	    }
	    fin.fixedtsamp = period/(double)fin.NrBins;
	  }
	}
	if(fin.NrSubints > 1) {
	  if(fin.freqMode != FREQMODE_UNIFORM) {
	    if(fin.freqMode == FREQMODE_FREQTABLE) {
	      printwarning(application.verbose_state.debug, "WARNING pplot: Different subints appear to have different effective frequencies. The subints will not be aligned properly if they are dedispersed with different reference frequencies.");
	    }else {
	      printwarning(application.verbose_state.debug, "WARNING pplot: The subints might not be aligned properly if they are dedispersed with different reference frequencies.");
	    }
	    if(showtop) {
	      printwarning(application.verbose_state.debug, "WARNING pplot: As a result, the -shoptop option might not result in a proper pulse profile.");
	    }
	  }
	}

	/* Convert the specified onpulse regions which were defined as
	   fractions on the commandline. This is now possible because we
	   know the number of bins. */
	region_frac_to_int(&(application.onpulse), fin.NrBins, 0);
	if(data_read == 0) {  /* Don't do this when switching polarization, only the first time the data is read */
	  if(setBaselineParams(fin, &baseline, &dxshift, &xUnitsSwitch, application.verbose_state) == 0)
	    return 0;
	  if(fin.yrangeset) {
	    yUnitsSwitch = 1;
	  }
	  if(longitudeRangeSet_flag == 0) {
	    longitude_start = 0 + dxshift;
	    longitude_end = baseline + dxshift;
	  }
	}   /* End of if(data_read == 0) */
	/* Get rid of polarization if it is there, so need to store the original nr of polarizations for future use, unless it is a make_paswing file */
	nrpolarizations = fin.NrPols;
	//	if(current_polnr != 0 && (fin.poltype == POLTYPE_ILVPAdPA || fin.poltype == POLTYPE_PAdPA)) {
	//	  fflush(stdout);
	//	  printwarning(application.verbose_state.debug, "WARNING pplot: Polarization select is not implemented for PA-swing files");
	//	}
	//	if(fin.NrPols > 1 && (fin.poltype != POLTYPE_ILVPAdPA && fin.poltype != POLTYPE_PAdPA)) { // Select polarization channel if there are multiple channels present
	if(fin.NrPols > 1) { // Select polarization channel if there are multiple channels present
	  datafile_definition clone;
	  if(current_polnr >= fin.NrPols)
	    current_polnr = 0;
	  //	  closePSRData(&clone, 0, application.verbose_state);  // Release header related memory if previously used in interactive mode when pressing 'P'
	  if(preprocess_polselect(fin, &clone, current_polnr, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pplot: Error selecting polarization channel %d.", current_polnr);
	    return 0;
	  }
	  swap_orig_clone(&fin, &clone, application.verbose_state); 
	}
	// Did zapping here     
	didtranspose_orig_nrbin = 0;

	ytitle_set_to_intensity = 0;
	ytitle_set_to_pulsenumber = 0;
	ytitle_set_to_subint = 0;
	ytitle_set_to_freqchannel = 0;
	ytitle_set_to_fluctuationFreq = 0;
	ytitle_set_to_padist = 0;
	ytitle_set_to_elldist = 0;
	ytitle_set_to_pulselongitude = 0;
	ytitle_set_to_p3fold = 0;
	ytitle_set_to_spectral_power = 0;
	ytitle_set_to_harmonic_number = 0;
	ytitle_set_to_lag_number = 0;
	/* If there is frequency and subint data, then transpose the data to make a subint/frequency plot */
	if(fin.NrFreqChan > 1) {
	  datafile_definition clone;
	  didtranspose_orig_nrbin = fin.NrBins;
	  //	  closePSRData(&clone, 0, application.verbose_state);  // Release header related memory if previously used
	  if(preprocess_transposeRawFBdata(fin, &clone, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pplot: Error transposing data.");
	    return 0;
	  }
	  /* Horizontal scale will be in subint number */
	  if(longitudeRangeSet_flag == 0) {
	    longitude_start = 0;
	    longitude_end = clone.NrSubints;
	    if(fin.NrSubints > 1)
	      xUnitsSwitch = XUNIT_PHASE;
	    else
	      xUnitsSwitch = XUNIT_BINS;
	  }
	  /* After the transpose, interpret the plot like a pulsestack */
	  clone.NrSubints = fin.NrFreqChan;
	  clone.NrBins = fin.NrSubints * fin.NrBins;
	  clone.NrFreqChan = 1;
	  if(ytitleset == 0)
	    ytitle_set_to_freqchannel = 1;
	  if(yUnitsMHz)
	    yUnitsSwitch = 1;

	  /*	  if(xtitle_set == 0) {
	    //	    xtitle_set = 1;
	    if(fin.NrSubints > 1) {
	      sprintf(xtitle, "Subint");
	    }else {
	      sprintf(xtitle, "Pulse longitude");
	    }
	    if(xUnitsSwitch == XUNIT_BINS) {
	      strcat(xtitle, " (bins)");
	    }else if(xUnitsSwitch == XUNIT_DEG) {
	      strcat(xtitle, " (deg)");
	    }else if(xUnitsSwitch == XUNIT_TIME) {
	      strcat(xtitle, " (sec)");
	    }
	    }*/
	  if(xUnitsSwitch == XUNIT_BINS) {
	    baseline = fin.NrSubints * fin.NrBins;
	  }else if(xUnitsSwitch == XUNIT_DEG) {
	    baseline = fin.NrSubints * 360;
	  }else if(xUnitsSwitch == XUNIT_PHASE) {
	    baseline = fin.NrSubints;
	  }else if(xUnitsSwitch == XUNIT_TIME) {
	    baseline = fin.NrSubints * fin.NrBins * get_tsamp(fin, 0, application.verbose_state);
	  }
	  //	  fprintf(stderr, "XXXXX Transpose: baseline = %f\n", baseline);
	  swap_orig_clone(&fin, &clone, application.verbose_state); 
	}else {                     // If there are no frequency channels
	  if(ytitleset == 0) {
	    if((fin.gentype == GENTYPE_2DFS || fin.gentype == GENTYPE_LRFS || fin.gentype == GENTYPE_S2DFSP3 || fin.gentype == GENTYPE_S2DFSP2) && (vectorMinMaxSpecified == 0 && fin.yrangeset == 1)) {
	      ytitle_set_to_fluctuationFreq = 1;
	    }else if(fin.gentype == GENTYPE_P3FOLD) {
	      ytitle_set_to_p3fold = 1;
	    }else if(fin.gentype == GENTYPE_HRFS_UNFOLDED) {
	      ytitle_set_to_spectral_power = 1;
	    }else if(fin.gentype == GENTYPE_HRFS) {
	      ytitle_set_to_harmonic_number = 1;
	    }else if(fin.gentype == GENTYPE_LRAC) {
	      ytitle_set_to_lag_number = 1;
	    }else if((fin.gentype == GENTYPE_LRCC) && (vectorMinMaxSpecified == 0 && fin.yrangeset == 1)) {
	      ytitle_set_to_pulselongitude = 1;
	    }else if(fin.gentype == GENTYPE_PADIST&& (vectorMinMaxSpecified == 0 && fin.yrangeset == 1)) {
	      ytitle_set_to_padist = 1;
	    }else if(fin.gentype == GENTYPE_ELLDIST&& (vectorMinMaxSpecified == 0 && fin.yrangeset == 1)) {
	      ytitle_set_to_elldist = 1;
	    }else if((fin.gentype == GENTYPE_RMMAP) && (vectorMinMaxSpecified == 0 && fin.yrangeset == 1)) {
	      strcpy(ytitle, "RM (rad/m\\u2\\d)");
	    }else if((fin.gentype == GENTYPE_PULSESTACK || fin.gentype == GENTYPE_PROFILE) && fin.NrSubints == 1) {
	      ytitle_set_to_intensity = 1;   // Handled lower down in the code
	    }else if((fin.gentype == GENTYPE_SEARCHMODE) && fin.NrSubints == 1 && fin.NrFreqChan == 1) { // A time-series
	      ytitle_set_to_intensity = 1;   // Handled lower down in the code
	    }else if(fin.gentype == GENTYPE_PULSESTACK) {
	      ytitle_set_to_pulsenumber = 1;
	      if(yUnitsTobs && fin.NrSubints > 1)
		yUnitsSwitch = 1;
	    }else {
	      ytitle_set_to_subint = 1;
	      if(yUnitsTobs && fin.NrSubints > 1)
		yUnitsSwitch = 1;
	    }
	  }
	}

	if(data_read == 0) {
	  if(subint_start >= fin.NrSubints) {
	    subint_start = fin.NrSubints -1;
	    subint_end = fin.NrSubints -1;
	  }
	  if(subint_start < 0) {
	    subint_start = 0;
	    subint_end = fin.NrSubints -1;
	  }
	  if(subint_end >= fin.NrSubints) {
	    subint_end = fin.NrSubints -1;
	  }
	  if(bin1 < 0) {
	    bin1 = 0;
	    bin2 = fin.NrBins -1;
	  }
	  if(xUnitsSwitch != XUNIT_BINS) {
	    bin1 = fin.NrBins*(longitude_start-dxshift)/baseline;
	    bin2 = fin.NrBins*(longitude_end-dxshift)/baseline - 1;
	    if(bin2 >= fin.NrBins)
	      bin2 = fin.NrBins-1;
	    if(bin1 < 0)
	      bin1 = 0;
	  }
	  if(application.verbose_state.verbose) printf("Specified binrange: %d - %d\n", bin1, bin2);
	}

	// Allocate temporary memory to produce plot and for the subset of the data to be plotted
	if(dataSubset_allocated) {
	  free(dataSubset);
	}
	if(subint_end < subint_start) {
	  printerror(application.verbose_state.debug, "ERROR pplot: The selected subint range seems to be such that the first subint exceeds the final.");
	  return 0;
	}
	dataSubset = (float *)malloc((subint_end-subint_start+1)*fin.NrBins*sizeof(float));
	//	printf("XXXXX dataSubset allocated for %ld X %ld floats\n", subint_end-subint_start+1, fin.NrBins);
	dataSubset_allocated = 1;
	if(stack_poly_allocated) {
	  free(xstack_poly);
	  free(ystack_poly);
	}
	if(maxy_allocated)
	  free(maxy);
	if(polymode_flag) {
	  if(histswitch) {  // In a histogram there are twice as many points drawn
	    xstack_poly = (float *)malloc(2*(fin.NrBins+2)*sizeof(float));
	    ystack_poly = (float *)malloc(2*(fin.NrBins+2)*sizeof(float));
	  }else {
	    xstack_poly = (float *)malloc((fin.NrBins+2)*sizeof(float));
	    ystack_poly = (float *)malloc((fin.NrBins+2)*sizeof(float));
	  }
	  stack_poly_allocated = 1;
	}else {
	  maxy = (long *)malloc(fin.NrBins*sizeof(long));
	  maxy_allocated = 1;
	}
	if(dataSubset == NULL 
	   || (maxy_allocated && maxy == NULL)
	   || (stack_poly_allocated && (xstack_poly == NULL || ystack_poly == NULL))) {
	  printerror(application.verbose_state.debug, "ERROR pplot: Memory allocation error.");
	  if(dataSubset == NULL && application.verbose_state.debug) {
	    printerror(application.verbose_state.debug, "Cannot allocate %ld bytes", (subint_end-subint_start+1)*fin.NrBins*sizeof(float));
	    printerror(application.verbose_state.debug, "subint_end=%ld subint_start=%ld fin.NrBins=%ld", subint_end, subint_start, fin.NrBins);
	  }
	  return 0;
	}

	/* Copy the relevant subint range of data in dataSubset */
	for(i = subint_start; i <= subint_end; i++) {
	  if(!readPulsePSRData(&fin, i, 0, 0, 0, fin.NrBins, &dataSubset[(i-subint_start)*fin.NrBins], application.verbose_state)) {
	    printerror(application.verbose_state.debug, "ERROR pplot: Cannot read data.");
	    return 0;
	  }
	}

	// If no -map or -lp is specified on command line, do a line plot if there is one subintegration, or a map if there are multiple
	if(stack_state[current_stack_pos].grayscalemode == -1) {
	  if(fin.NrSubints > 1) {
	    stack_state[current_stack_pos].grayscalemode = 1;
	  }else {
	    stack_state[current_stack_pos].grayscalemode = 0;
	  }
	}


	if(data_read == 0) {
	  if(change_pulse_numbering_flag) {
	    if(interactive_flag != 0 || stack_state[current_stack_pos].grayscalemode != 0) {
	      printwarning(application.verbose_state.debug, "WARNING: -map/-cont mode and interactive mode are incompatable with the -0 and -1 options. Commands are ignored.");
	    }else {
	      subint_end -= subint_start - change_pulse_numbering_value;
	      subint_start = change_pulse_numbering_value;
	    }
	  }
	  
	  firstPulsedataSubset = subint_start;

	  stack_state[current_stack_pos].subint_start = subint_start;
	  stack_state[current_stack_pos].subint_end = subint_end;
	  stack_state[current_stack_pos].x1 = bin1;
	  stack_state[current_stack_pos].x2 = bin2;
	  stack_state[current_stack_pos].nrZappedVectors = 0;
	  stack_state[current_stack_pos].nrZappedSubints = 0;
	  current_stack_pos++;
	}

	// See if zapping needs to be applied.
	// When running for the first time, current_stack_pos is zero at this point.
	if(current_stack_pos > 0) {
	  if(stack_state[current_stack_pos-1].nrZappedSubints > 0 || stack_state[current_stack_pos-1].nrZappedVectors) {
	    //	    printf("XXXXX %d %d (stack=%d)\n", stack_state[current_stack_pos-1].nrZappedSubints, stack_state[current_stack_pos-1].nrZappedVectors, current_stack_pos);
	    /* If the file is read again, make sure to apply zapping again. */
	    if(zapVectors(stack_state[current_stack_pos-1].nrZappedVectors, zappedVectors, subint_start, fin, dataSubset, -1, application.verbose_state) == 0)
	      return 0;
	    
	    // If there are freqs and subints, the zapped subints are separate from the zapped vectors (which are the freq. channels).
	    if(zapSubints(stack_state[current_stack_pos-1].nrZappedSubints, zappedSubints, subint_start, didtranspose_orig_nrbin, fin, dataSubset, -1, application.verbose_state) == 0)
	      return 0;
	  }
	}

      }   /* End of if(data_read == 0 || data_read == 2) */
      
      if(setBaselineParams(fin, &baseline, &dxshift, &xUnitsSwitch, application.verbose_state) == 0)
	return 0;

      if(heading_string_index && curpanelnrx == 0 && curpanelnry == 0) {
	ppgsch(heading_font.characterheight);
	ppgslw(heading_font.linewidth);
	ppgscf(heading_font.font);
	ppgsvp(0.1, 0.9, 0.1, 0.9);    // Set a default viewport first to generage the heading, important if multiple plots
	//	fprintf(stderr, "XXXXXX pgswin1 0, 1, 0, 1\n");
	ppgswin(0, 1, 0, 1);    // Set a default viewport first to generage the heading, important if multiple plots
	//	ppgmtxt("t", 0.15, 0.5, 0.5, argv[heading_string_index]);
	ppgptxt(0.5, 1.07, 0, 0.5, argv[heading_string_index]);
      }

      //START REGION DEVELOP
      if(copyright[0] != 0 && curpanelnrx == 0 && curpanelnry == 0) {
	ppgsvp(0.1, 0.9, 0.1, 0.9);    // Set a default viewport first to generage the heading, important if multiple plots
	//	fprintf(stderr, "XXXXXX pgswin2 0, 1, 0, 1\n");
	ppgswin(0, 1, 0, 1);    // Set a default viewport first to generage the heading, important if multiple plots
	ppgsch(copy_font.characterheight);
	ppgslw(copy_font.linewidth);
	ppgscf(copy_font.font);
	//	printf("XXXXX %f %f (%f %f)\n", y, y2, pgplot_frame.swin_y1, pgplot_frame.swin_y2);
	ppgptxt(1, -0.07+dycopyright, 0, copyrightalignment, copyright);
	/*
	if(appendframes_flag == 0)
	  ppgptxt(x2, pgplot_frame.swin_y1-0.11*(pgplot_frame.swin_y2-pgplot_frame.swin_y1)+dycopyright, 0, copyrightalignment, copyright);
	else if()
	  ppgptxt(x2, pgplot_frame.swin_y1-0.11*(pgplot_frame.swin_y2-pgplot_frame.swin_y1)+dycopyright, 0, copyrightalignment, copyright);
	*/
	ppgsch(1);
      }
    
      /*
      ppgslw(1);
      ppgsch(0.75);
      ppgscf(3);
      */
      //START REGION RELEASE

      ppgsvp(viewport_startxCurrentPanel, viewport_endxCurrentPanel, viewport_startyCurrentPanel, viewport_endyCurrentPanel);
    
      xmin = stack_state[current_stack_pos-1].x1;
      xmax = stack_state[current_stack_pos-1].x2;
      if(histswitch || stack_state[current_stack_pos-1].grayscalemode) {  // A histogram/colour map covers one bin more of baseline. Note, colour maps are handled separately below.
	xmin -= 0.5;
	xmax += 0.5;
      }
      if(xUnitsSwitch != XUNIT_BINS) {
	xmin *= baseline/(float)fin.NrBins;
	xmin += dxshift;
	xmax *= baseline/(float)fin.NrBins;
	xmax += dxshift;
      }
      y = stack_state[current_stack_pos-1].subint_start;
      y2 = stack_state[current_stack_pos-1].subint_end;
      //      printf("XXXXXX y, y2 = %f %f\n", y, y2);
      /* If only showing one pulse, the pulse number is not added to the y value */
      if(stack_state[current_stack_pos-1].subint_start == stack_state[current_stack_pos-1].subint_end) {
	if(stack_state[current_stack_pos-1].grayscalemode == 0) {
	  /* If vertical scale is fixed and just showing one pulse, assume it is normalised. */
	  //	  if(fixverticalscale_flag) {
	  //	    y = 0;
	  //	    y2 = 1;
	  //	  }else {
	  y = y2 = ypos(dataSubset, stack_state[current_stack_pos-1].x1, stack_state[current_stack_pos-1].subint_start, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
	  //	  }
	}
      }
      // Moved it to below for stargazing live. I assume this is what you want anyway. It now can also be used for a pulse-stack
      if(fixverticalscale_flag) {
	y = 0;
	y2 = 1;
      }
      //      printf("XXXXXX y, y2 = %f %f\n", y, y2);
      if(fixverticalscale_flag == 0 && stack_state[current_stack_pos-1].grayscalemode == 0) {
	for(i = stack_state[current_stack_pos-1].subint_start; i <= stack_state[current_stack_pos-1].subint_end; i++) {
	  for(j = stack_state[current_stack_pos-1].x1; j <= stack_state[current_stack_pos-1].x2; j++) {
	    if(ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) < y)
	      y  = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
	    if(ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) > y2)
	      y2 = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
	  }
	}
      }else if(y == y2 && stack_state[current_stack_pos-1].grayscalemode == 0){
	y2 += 1;
      }
      //      fprintf(stderr, "XXXXXX y, y2 = %f %f\n", y, y2);
      ymin = y-ymargin_both;
      ymax = y2+ymargin_both+ymargin_top;

      xminshow = xmin;
      xmaxshow = xmax;

      pgplot_frame.swin_showtwice = showTwice_flag;
      pgplot_frame.swin_x1 = xminshow;
      pgplot_frame.swin_x2 = xmaxshow;
      pgplot_frame.swin_y1 = ymin;
      pgplot_frame.swin_y2 = ymax;
      //      fprintf(stderr, "XXXXXX pgswin3 %f %f %f %f\n", xminshow,xmaxshow,ymin,ymax);
      if(stack_state[current_stack_pos-1].grayscalemode == 0) {   
	ppgswin(xminshow,xmaxshow,ymin,ymax);
      }

      if(stack_state[current_stack_pos-1].grayscalemode) {   
	//	printf("XXXXXX x1=%ld x2=%ld\n", stack_state[current_stack_pos-1].x1, stack_state[current_stack_pos-1].x2);
	float max, min, oldmin, oldmax;
	float xleft, xright, xleft2, xright2;
	min = max = dataSubset[stack_state[current_stack_pos-1].x1+ignorebins];
	//	printf("XXXX dataSubset_allocated=%d\n", dataSubset_allocated);
	for(i = stack_state[current_stack_pos-1].x1; i <= stack_state[current_stack_pos-1].x2; i++) {
	  for(j = 0; j <= stack_state[current_stack_pos-1].subint_end-stack_state[current_stack_pos-1].subint_start; j++) {
	    //	    printf("XXXXX investigating subint %ld (%ld %ld %d)\n", j+stack_state[current_stack_pos-1].subint_start-firstPulsedataSubset, j, stack_state[current_stack_pos-1].subint_start, firstPulsedataSubset);
	    if(dataSubset[(j+stack_state[current_stack_pos-1].subint_start-firstPulsedataSubset)*fin.NrBins+i] > max)
	      max = dataSubset[(j+stack_state[current_stack_pos-1].subint_start-firstPulsedataSubset)*fin.NrBins+i];
	    if(dataSubset[(j+stack_state[current_stack_pos-1].subint_start-firstPulsedataSubset)*fin.NrBins+i] < min)
	      min = dataSubset[(j+stack_state[current_stack_pos-1].subint_start-firstPulsedataSubset)*fin.NrBins+i];
	  }
	}
	//	fprintf(stderr, "XXXXXX min/max = %f %f\n", min, max);
	oldmin = min;
	oldmax = max;
	max = min + (max-min)/scale;      
	min = oldmin + (oldmax-oldmin)*scale2;
	xleft = stack_state[current_stack_pos-1].x1+ignorebins;
	xright = stack_state[current_stack_pos-1].x2-ignorebins-ignorebins2;
	xleft2 = 0;
	xright2 = fin.NrBins-1;
	// Note that baseline covered is one bin more than range in bin-centres, so extend the viewed range. Note that the labels are not drawn here.
	xleft -= 0.5;
	xright += 0.5;
	if(xUnitsSwitch != XUNIT_BINS) {
	  xleft *= baseline/(float)(fin.NrBins);
	  xright *= baseline/(float)(fin.NrBins);
	  xleft += dxshift;
	  xright += dxshift;
	  xleft2 = 0+dxshift;
	  xright2 = baseline*(fin.NrBins-1)/(float)(fin.NrBins)+dxshift;
	}
	/*      printf("L=%f R=%f (L=%f R=%f) %d\n", xleft, xright, xleft2, xright2, xUnitsSwitch); */
	printf("Plotting vectors (vertical axis): %ld - %ld\n", stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end);
	if(nonumside_flag) {
	  if(showright)
	    showright = 2;
	  if(showtop)
	    showtop = 2;
	}
	ppgsvp(viewport_startxCurrentPanel, viewport_endxCurrentPanel, viewport_startyCurrentPanel, viewport_endyCurrentPanel);
	
	dummyf1 = (0+change_pulse_numbering_value)*yUnitCmdLine+dyshift;
	dummyf2 = (fin.NrSubints-1+change_pulse_numbering_value)*yUnitCmdLine+dyshift;
	dummyf3 = (stack_state[current_stack_pos-1].subint_start-0.5)*yUnitCmdLine+dyshift;
	dummyf4 = (stack_state[current_stack_pos-1].subint_end+0.5)*yUnitCmdLine+dyshift;
	dummyi = showright;
	if((fin.yrangeset && yUnitsSwitch) || vectorMinMaxSpecified) {
	  if(showright)
	    dummyi += 2;
	}
      
	pgplot_options_definition pgplot_options;
	pgplot_clear_options(&pgplot_options);
	pgplot_options.viewport.windowwidth = application.windowwidth;
	pgplot_options.viewport.windowheight = application.windowheight;

	//START REGION DEVELOP
	// If viewport_startx=0.15 (default), dxplot should be set to zero, which corresponds to a default shift of 0.15
	//START REGION RELEASE
	pgplot_options.viewport.dxplot = viewport_startxCurrentPanel-0.15;
	//START REGION DEVELOP
	// Default values correspond to xsize = 1, which makes plot go from 0.15 to 0.9
	//START REGION RELEASE
	pgplot_options.viewport.xsize = (viewport_endxCurrentPanel-viewport_startxCurrentPanel)/(0.9-0.15);
	
	//	printf("XXXX %d dx=%f xsize=%f\n", curpanelnrx, pgplot_options.viewport.dxplot, pgplot_options.viewport.xsize);
	
	pgplot_options.viewport.dyplot = viewport_startyCurrentPanel-0.15;
	pgplot_options.viewport.ysize = (viewport_endyCurrentPanel-viewport_startyCurrentPanel)/(0.9-0.15);
	pgplot_options.viewport.noclear = 1;  // Clearing is already done above where viewport ranges are set
	strcpy(pgplot_options.viewport.plotDevice, application.pgplotdevice); 
	pgplot_options.viewport.dontopen = 1;
	pgplot_options.viewport.dontclose = 1;
	newtext = str_replace_header_params(fin, title, application.verbose_state);
	if(newtext == NULL) {
	  fflush(stdout);
	  printwarning(application.verbose_state.debug, "WARNING pplot: Cannot substitute keyword in title");
	  pgplot_options.box.title[0] = 0;
	}else {
	  strcpy(pgplot_options.box.title, newtext);
	  free(newtext);
	}
	pgplot_options.box.title_ch = title_font.characterheight*currentPanelScaling;
	pgplot_options.box.title_lw = title_font.linewidth;
	pgplot_options.box.title_f = title_font.font;
	pgplot_options.box.label_f = title_font.font;
	pgplot_options.box.label_ch = label_font.characterheight*currentPanelScaling;
	pgplot_options.box.box_labelsize = box_font.characterheight*currentPanelScaling;
	pgplot_options.box.box_lw = box_font.linewidth;
	pgplot_options.box.label_lw = box_font.linewidth;
	if(wedgelabel_set)
	  strcpy(pgplot_options.box.wedgelabel, wedgelabel);
	int levelset = 1;
	if(didtranspose_orig_nrbin)   // If there are frequency channels, let pgplotmap determine levels, as something goes wrong with datasubset, as probably transpose is not applied when copying the subset.
	  levelset = 0;
	if(scalerange) {  // If -scalerange is specified, overwrite the found range
	  levelset = 1;
	  min = scalerange_min;
	  max = scalerange_max;
	}
	if(stack_state[current_stack_pos-1].grayscalemode == 1) {   // Colour map
	  if(pgplotMap(&pgplot_options, fin.data, fin.NrBins, fin.NrSubints, xleft2, xright2, xleft, xright, dummyf1, dummyf2, dummyf3, dummyf4, application.cmap, application.itf, 0, 0, NULL, 1, 0, 1, levelset, min, max, 1, 2, dummyi, 0, showtop, 0, plotlw, showwedge, !application.do_noplotsubset, showTwice_flag, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR pplot: Cannot plot data.");
	    return 0;
	  }
	}else if(stack_state[current_stack_pos-1].grayscalemode == 2) {    // Contour plot
	  if(pgplotMap(&pgplot_options, fin.data, fin.NrBins, fin.NrSubints, xleft2, xright2, xleft, xright, dummyf1, dummyf2, dummyf3, dummyf4, application.cmap, application.itf, 1, stack_state[current_stack_pos-1].nrcontours, NULL, plotlw, 0, 1, levelset, min, max, 1, 2, dummyi, 0, showtop, 0, plotlw, showwedge, !application.do_noplotsubset, showTwice_flag, application.verbose_state) == 0) {    
	    printerror(application.verbose_state.debug, "ERROR pplot: Cannot plot data.");
	    return 0;
	  }
	}else {     // Both
	  if(pgplotMap(&pgplot_options, fin.data, fin.NrBins, fin.NrSubints, xleft2, xright2, xleft, xright, dummyf1, dummyf2, dummyf3, dummyf4, application.cmap, application.itf, 0, stack_state[current_stack_pos-1].nrcontours, NULL, plotlw, 0, 1, levelset, min, max, 1, 2, dummyi, 0, showtop, 0, plotlw, showwedge, !application.do_noplotsubset, showTwice_flag, application.verbose_state) == 0) {    
	    printerror(application.verbose_state.debug, "ERROR pplot: Cannot plot data.");
	    return 0;
	  }
	}
      }else {   /* if(stack_state[current_stack_pos-1].grayscalemode): i.e. if not making a grayscale plot, make a line plot */
	printf("Plotting vectors (vertical axis): %ld - %ld\n", stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end);
	ppgbbuf();  // Start buffering, which greatly increases speed of drawing.
	if(polymode_flag) {
	  /* Set linewidth to 1 for fill */
	  ppgslw(1);
	  for(i = stack_state[current_stack_pos-1].subint_end - ignorelastpulses; i >= stack_state[current_stack_pos-1].subint_start; i--) {
	    k = 0;
	    x = 0;
	    // First draw a polygon to erase background, then draw line later blow
	    for(j = stack_state[current_stack_pos-1].x1+ignorebins; j <= stack_state[current_stack_pos-1].x2-ignorebins; j++) {
	      if(!histswitch) {
		x = j;
	      }else {
		x = j - 0.5;
	      }
	      if(xUnitsSwitch != XUNIT_BINS) {
		x *= baseline/(float)fin.NrBins;
		x += dxshift;
	      }
	      y = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
	      if(k == 0) {
		xstack_poly[k] = x-100*0;
		ystack_poly[k] = ymin-100;
		k++;
	      }
	      int loopnr;
	      for(loopnr = 0; loopnr < 2; loopnr++) {  // Histogram requires two points to be drawn
		if(loopnr == 0) {
		  xstack_poly[k] = x;
		  ystack_poly[k] = y;
		  k++;
		}else {
		  float x2;
		  x2 = j + 0.5;
		  if(xUnitsSwitch != XUNIT_BINS) {
		    x2 *= baseline/(float)fin.NrBins;
		    x2 += dxshift;
		  }
		  xstack_poly[k] = x2;
		  ystack_poly[k] = y;
		  k++;
		}
		//	      fprintf(stderr, "(x,y) = (%f,%f))\n", xstack_poly[k], ystack_poly[k]);
		if(k == MaxNrPointsInPolygon) {
		  xstack_poly[k] = x+100*0;
		  ystack_poly[k] = ymin-100;
		  k++;
		  ppgsfs(1);
		  ppgsci(0);
		  ppgpoly(k, xstack_poly, ystack_poly);
		  /* Remember last few points (plus one out outside window) to make sure everything is filled. */
		  xstack_poly[0] = xstack_poly[k-3];
		  ystack_poly[0] = ymin-100;
		  xstack_poly[1] = xstack_poly[k-3];
		  ystack_poly[1] = ystack_poly[k-3];
		  xstack_poly[2] = xstack_poly[k-2];
		  ystack_poly[2] = ystack_poly[k-2];
		  k = 3;
		}
		if(!histswitch) {
		  break;
		}
	      }
	    }
	    xstack_poly[k] = x+100*0;
	    ystack_poly[k] = ymin-100;
	    k++;
	    ppgsfs(1);
	    ppgsci(0);
	    ppgpoly(k, xstack_poly, ystack_poly);
	    /* Now draw the actual line on top of the polygon that is used to erase any background lines */
	    ppgslw(plotlw);
	    ppgsci(1);
	    for(j = stack_state[current_stack_pos-1].x1+ignorebins; j <= stack_state[current_stack_pos-1].x2-ignorebins; j++) {
	      if(!histswitch) {
		x = j;
	      }else {
		x = j - 0.5;
	      }
	      if(xUnitsSwitch != XUNIT_BINS) {
		x *= baseline/(float)fin.NrBins;
		x += dxshift;
	      }
	      y = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
	      if(j == stack_state[current_stack_pos-1].x1+ignorebins) {
		ppgmove(x, y);
	      }else {
		ppgdraw(x, y);
	      }
	      if(histswitch) {
		x = j + 0.5;
		if(xUnitsSwitch != XUNIT_BINS) {
		  x *= baseline/(float)fin.NrBins;
		  x += dxshift;
		}
		ppgdraw(x, y);
	      }

	    }
	  }
	}else {   // End of if(polymode_flag)
	  int skip, didfirstmove;
	  double xdouble, ydouble;
	  skip = 0;
	  ppgslw(plotlw);
	  for(j = 0; j < fin.NrBins; j++)
	    maxy[j] = stack_state[current_stack_pos-1].subint_start;
	  for(i = stack_state[current_stack_pos-1].subint_start; i <= stack_state[current_stack_pos-1].subint_end - ignorelastpulses; i++) {
	    x = 0;
	    if(xUnitsSwitch != XUNIT_BINS) {
	      x *= baseline/(float)fin.NrBins;
	      x += dxshift;
	    }
	    ppgmove(x, ypos(dataSubset, 0,i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift));
	    didfirstmove = 0;
	    for(j = stack_state[current_stack_pos-1].x1+ignorebins; j <= stack_state[current_stack_pos-1].x2-ignorebins; j++) {
	      y =    ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
	      if(y > ypos(dataSubset, j, maxy[j], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) || i == maxy[j]) {
		if(skip == 0) {
		  if(!histswitch) {
		    x = j;
		    if(xUnitsSwitch != XUNIT_BINS) {
		      x *= baseline/(float)fin.NrBins;
		      x += dxshift;
		    }
		    if((i != stack_state[current_stack_pos-1].subint_end - ignorelastpulses || x <= xmax_oscil) && (i != stack_state[current_stack_pos-1].subint_start || x >= xmin_oscil)) {
		      if(didfirstmove == 0) {
			ppgmove(x, y);
			didfirstmove = 1;
		      }else {
			/* This is the most common draw command. Happens when line is above lines of previous pulse */
			ppgdraw(x, y);
		      }
		    }
		    //START REGION DEVELOP
		    else if(plotdot && plotdotonlyonce == 0 && !(i != stack_state[current_stack_pos-1].subint_end - ignorelastpulses || x <= xmax_oscil)) {
		      ppgslw(plotdot);
		      ppgpt1(xmax_oscil, y, -1);
		      plotdotonlyonce = 1;
		      ppgslw(plotlw);
		    }
		    //START REGION RELEASE
		  }else { // If histswitch is set
		    //		    x = j-1;  // Modified this line, since first bin was missing
		    x = j-0.5;
		    y = ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
		    if(xUnitsSwitch != XUNIT_BINS) {
		      x *= baseline/(float)fin.NrBins;
		      x += dxshift;
		    }
		    ppgmove(x, y);
		    y = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
		    if((i != stack_state[current_stack_pos-1].subint_end - ignorelastpulses || x <= xmax_oscil) && (i != stack_state[current_stack_pos-1].subint_start || x >= xmin_oscil))
		      ppgdraw(x, y);
		    //		    x = j;
		    x = j+0.5; // Modified this line, since first bin was missing
		    if(xUnitsSwitch != XUNIT_BINS) {
		      x *= baseline/(float)fin.NrBins;
		      x += dxshift;
		    }
		    if((i != stack_state[current_stack_pos-1].subint_end - ignorelastpulses || x <= xmax_oscil)  && (i != stack_state[current_stack_pos-1].subint_start || x >= xmin_oscil))
		      ppgdraw(x, y);
		  }
		}else {
		  if(j != stack_state[current_stack_pos-1].x1) {
		    xdouble = ((double)ypos(dataSubset, j-1, maxy[j-1], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) - (double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift))/(double)((double)ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) - (double)((double)ypos(dataSubset, j, maxy[j-1], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, maxy[j-1], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)));
		    ydouble = (double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) + xdouble*((double)ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift));
		    x = xdouble;
		    y = ydouble;
//START REGION DEVELOP
		      /*  There seems to be a bug in first bin. x should be between 0 and 1. 
			  
I think now there are exceptions that are not handled. Do a:
fakePulsar -cone '1 4.5 1 9 600 0 -1 1000 0' -micro "1 10 3" -np 15 -gg test.gg
pplot -prange "0 1" -b "130 180" test.gg 
and mark the line segments that are drawn in with the next ppgdraw and you see that things sometimes go wrong (around bin 143 or so). This check makes sure no completely wrong lines are drawn. Happens when line goes behind and next point is visible again?
*/
//START REGION RELEASE
		    if(x >= 0 && x <= 1) {
		      x += j-1;
		      if(xUnitsSwitch != XUNIT_BINS) {
			x *= baseline/(float)fin.NrBins;
			x += dxshift;
		      }
		      ppgmove(x, y);
		    }
//START REGION DEVELOP
		      /*
			else {
		      ppgmove(j-1+x, y);
		      y2 = ypos(j-1, maxy[j]) + x*(ypos(j, maxy[j])-ypos(j-1, maxy[j]));
		      printf("%f -> %f and %f -> %f\n", ypos(j-1, i), ypos(j, i), ypos(j-1, maxy[j]), ypos(j, maxy[j]));  
		      printf("pulse=%ld bin=%ld, %ld x=%f y=%f y2=%f\n\n", i, j, maxy[j], x, y, y2);  
		      } */
//START REGION RELEASE
		  }
		  y = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
		  x = j;
		  if(xUnitsSwitch != XUNIT_BINS) {
		    x *= baseline/(float)fin.NrBins;
		    x += dxshift;
		  }
//START REGION DEVELOP
		  /* These are the line segments for which the first part in not visible? */
		  /*		  ppgsci(2); */
//START REGION RELEASE
		  if((i != stack_state[current_stack_pos-1].subint_end - ignorelastpulses || x <= xmax_oscil)  && (i != stack_state[current_stack_pos-1].subint_start || x >= xmin_oscil))
		    ppgdraw(x, y);
		  /*		  ppgsci(1); */
		  skip = 0;
		}
	      }else {
		if(skip == 1) {
		  ppgmove(j, y);
		}else {
		  skip = 1;
		  
		  xdouble = ((double)ypos(dataSubset, j-1, maxy[j], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) - (double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift))/(double)((double)ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) - ((double)ypos(dataSubset, j, maxy[j], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, maxy[j], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)));
		  ydouble = (double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) + xdouble*((double)ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift)-(double)ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift));
		  x = xdouble;
		  y = ydouble;
		  x += j-1;
		  if(xUnitsSwitch != XUNIT_BINS) {
		    x *= baseline/(float)fin.NrBins;
		    x += dxshift;
		  }
		  
//START REGION DEVELOP
		  /* This shouldn't be necesary, but there is somewhere a bug in my program that makes such a test necesary. */
		    /*		  ppgsci(2); */
//START REGION RELEASE
		  if(x > 0 && x < 10000 && y > -stack_state[current_stack_pos-1].subint_end && y < 2*stack_state[current_stack_pos-1].subint_end) {
		    if((i != stack_state[current_stack_pos-1].subint_end - ignorelastpulses || x <= xmax_oscil)  && (i != stack_state[current_stack_pos-1].subint_start || x >= xmin_oscil)) {
//START REGION DEVELOP
		      /* Here it goes wrong */
		      /* For some reason calculation goes wrong sometimes. Only draw a line when calculated y position is below at least one of both start/end points, as it should. This avoids weird lines that shoot up. */
//START REGION RELEASE
		      if(y < ypos(dataSubset, j-1, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift) || y < ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift))
			ppgdraw(x, y); 
		    }
		  }
		  /*		  ppgsci(1);*/
		}
	      }
	      if(y >= ypos(dataSubset, j, maxy[j], scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift))
		maxy[j] = i;
	      
	    }
	  }  /* This is the end of the pulse number loop */
	} /* This is the end of the else corresponding to if(polymode_flag) */
	ppgebuf();	    //Stop buffering and draw to device
      } // This is end of the else corresponding to not in grayscale mode
      ppgsci(1);
      ppgsch(box_font.characterheight);
      ppgslw(box_font.linewidth);
      //      ppgscf(box_font.font);

      sprintf(txt1, "bc");
      sprintf(txt2, "bc");
      if(noboxx) {
	txt1[0] = 0;
	sprintf(txt2, "b");
      }
      if(noboxy) {
	sprintf(txt1, "b");
	txt2[0] = 0;
      }
    
      pgplot_options_definition pgplot_options;
      pgplot_clear_options(&pgplot_options);
      pgplot_options.box.drawbox = 0;
      pgplot_options.box.drawtitle = 0;
      pgplot_options.box.drawlabels = 0;
      pgplot_options.box.title_ch = title_font.characterheight*currentPanelScaling;
      pgplot_options.box.title_lw = title_font.linewidth;
      pgplot_options.box.title_f = title_font.font;
      pgplot_options.box.label_f = title_font.font;
      pgplot_options.box.label_ch = label_font.characterheight*currentPanelScaling;
      pgplot_options.box.box_labelsize = box_font.characterheight*currentPanelScaling;
      pgplot_options.box.box_lw = box_font.linewidth;
      pgplot_options.box.box_f = box_font.font;
      pgplot_options.box.label_lw = box_font.linewidth;
      //      printf("XXXXXX %f %d %d\n", pgplot_options.box.title_ch, pgplot_options.box.title_lw, pgplot_options.box.title_f);

      if(noboxx == 0 || noboxy == 0) {
	if(curpanelnrx == 0 || appendframes_flag == 0) {  // This used to be if(filenr == 0)
	  if(disable_x_numbers && disable_y_numbers ) {
	    strcat(txt1, "st");
	    strcat(txt2, "t");
	  }else if(disable_x_numbers && !disable_y_numbers) {
	    strcat(txt1, "st");
	    strcat(txt2, "nt");
	  }else if(!disable_x_numbers && !disable_y_numbers) {
	    strcat(txt1, "nst");
	    strcat(txt2, "nt");
	  }else if(!disable_x_numbers && disable_y_numbers) {
	    strcat(txt1, "nst");
	    strcat(txt2, "t");
	  }
	}else {
	  if(disable_x_numbers) {
	    strcat(txt1, "st");
	    strcat(txt2, "t");
	  }else {
	    strcat(txt1, "nst");
	    strcat(txt2, "t");
	  }
	}
	
	pgplot_options.box.drawbox = 1;   // First: Plot only the box + numbers
	
	if((fin.yrangeset && yUnitsSwitch)  // Special format with defined yrange
	   || (vectorMinMaxSpecified  && yUnitsSwitch) // -vecRange is defined
	   || (((fin.gentype == GENTYPE_SUBINTEGRATIONS || fin.gentype == GENTYPE_PULSESTACK) && fin.NrFreqChan == 1 && fin.NrSubints > 1) && yUnitsSwitch) // Observing time should be specified
	   || (yUnitsSwitch && didtranspose_orig_nrbin)  // Frequency in MHz should be specified
	   ) {
	  pgplot_frame.swin_showtwice = showTwice_flag;
	  pgplot_frame.swin_x1 = xminshow;
	  pgplot_frame.swin_x2 = xmaxshow;
	  pgplot_frame.swin_y1 = ymin;
	  pgplot_frame.swin_y2 = ymax;
	  strcpy(pgplot_options.box.box_xopt, txt1);
	  strcpy(pgplot_options.box.box_yopt, txt2);
	  //	  printf("XXXXX pgplotbox = %s %s\n", pgplot_options.box.box_xopt, pgplot_options.box.box_yopt);
	  float vectorMin2, vectorMax2;
	  //	  fprintf(stderr, "XXXXXX ymin,ymax = %f %f\n", ymin, ymax);
	  yvec2unit(fin, yUnitsSwitch, vectorMinMaxSpecified, vectorMin, vectorMax, didtranspose_orig_nrbin, 1, ymin, &vectorMin2, stack_state[current_stack_pos-1].grayscalemode, 0, application.verbose_state);
	  yvec2unit(fin, yUnitsSwitch, vectorMinMaxSpecified, vectorMin, vectorMax, didtranspose_orig_nrbin, 2, ymax, &vectorMax2, stack_state[current_stack_pos-1].grayscalemode, 0, application.verbose_state);
	  //	  fprintf(stderr, "XXXXXX ymin2,ymax2 = %f %f\n", vectorMin2, vectorMax2);
	  if(showTwice_flag == 0)
	    ppgswin(xminshow, xmaxshow, yUnitCmdLine*vectorMin2, yUnitCmdLine*vectorMax2);   // Set window in other units
	  else
	    ppgswin(xminshow, xmaxshow, yUnitCmdLine*vectorMin2, yUnitCmdLine*(vectorMax2+(vectorMax2-vectorMin2)));   // Set window in other units
	  pgplot_drawbox(&pgplot_options.box);   // Plot only the box + numbers
	  if(stack_state[current_stack_pos-1].grayscalemode == 0)
	    ppgswin(xminshow, xmaxshow, ymin, ymax);  // Go back to integer vector numbers
	  else {
	    ppgswin(xminshow, xmaxshow, ymin-0.5, ymax+0.5);  // Go back to integer vector numbers
	  }
	}else {
	  pgplot_options.box.drawbox = 1;
	  strcpy(pgplot_options.box.box_xopt, txt1);
	  strcpy(pgplot_options.box.box_yopt, txt2);
	  //	  printf("XXXXX pgplotbox3 = %s %s (box=%d labels=%d)\n", pgplot_options.box.box_xopt, pgplot_options.box.box_yopt, pgplot_options.box.drawbox, pgplot_options.box.drawlabels);
	  pgplot_drawbox(&pgplot_options.box);  // Plot only the box + numbers
	}

	//	if(fin.yrangeset && yUnitsSwitch) {
	  /* Use to be x, x2, instead of xminshow and xmaxshow */
	  /*
	  float vectorMin2, vectorMax2, vectorMin3, vectorMax3;
	  //	  fprintf(stderr, "XXXXXX ymin,ymax = %f %f\n", ymin, ymax);
	  vectorMin2 = fin.yrange[0];   // These are the centres of the vectors 
	  vectorMax2 = fin.yrange[1];
	  //	    fprintf(stderr, "XXXXX %f %f\n", vectorMin2, vectorMax2);
	  vectorMin3 = vectorMin2 - 0.5*(vectorMax2-vectorMin2)/(float)(fin.NrSubints-1);   // There are nrpulses-1 full bins spanning this range 
	  vectorMax3 = vectorMax2 + 0.5*(vectorMax2-vectorMin2)/(float)(fin.NrSubints-1);
	  //	    fprintf(stderr, "XXXXX %f %f\n", vectorMin3, vectorMax3);
	  if(showTwice_flag)
	    vectorMax3 += (fin.yrange[1] - fin.yrange[0])*fin.NrSubints/(float)(fin.NrSubints-1);
	  //START REGION DEVELOP
	  //	    vectorMin3 = (ymin/(float)(fin.NrSubints-1))*(vectorMax2 - vectorMin2) + vectorMin2;
	  //	    vectorMax3 = (ymax/(float)(fin.NrSubints-1))*(vectorMax2 - vectorMin2) + vectorMin2;
	  //	    vectorMin3 = vectorMin2;
	  //	    vectorMax3 = vectorMax2;
//START REGION RELEASE
//	  fprintf(stderr, "XXXXXX pgswin4 %f %f %f %f\n", xminshow, xmaxshow, vectorMin3, vectorMax3);
	  ppgswin(xminshow, xmaxshow, vectorMin3, vectorMax3);   // Set window in other units
	  //	  pgplot_frame.swin_y1 = vectorMin2;
	  //	  pgplot_frame.swin_y2 = vectorMax2;
	  //	  printf("XXXXX pgplotbox1 = %s %s (box=%d labels=%d)\n", pgplot_options.box.box_xopt, pgplot_options.box.box_yopt, pgplot_options.box.drawbox, pgplot_options.box.drawlabels);
	  pgplot_drawbox(pgplotbox);  // Plot only the box + numbers
	  //	  ppgswin(x, x2, vectorMin3, vectorMax3);
	  //	  fprintf(stderr, "XXXXXX pgswin5 %f %f %f %f\n", x, x2, ymin, ymax);
	  ppgswin(xminshow, xmaxshow, ymin, ymax);  // Go back to integer vector numbers
	  //	  ppgswin(x, x2, ymin, ymax);    // Go back to integer vector numbers
	  */
	  // Copied essentially from below, but now range is defines at centre of bin
	/*	}else if(vectorMinMaxSpecified) {  // Note: specified range corresponds to bottom and top of range
	  float vectorMin2, vectorMax2;
	  // Use to be x, x2, instead of xminshow and xmaxshow 
	  vectorMin2 = (ymin/(float)(fin.NrSubints-1))*(vectorMax - vectorMin) + vectorMin;
	  vectorMax2 = (ymax/(float)(fin.NrSubints-1))*(vectorMax - vectorMin) + vectorMin;
	  //	  printf("XXXXX ymin/max=%f %f and other variables are %f %f\n",ymin, ymax, vectorMin2, vectorMax2);
	  //	  fprintf(stderr, "XXXXXX pgswin6 %f %f %f %f\n", xminshow, xmaxshow, vectorMin2, vectorMax2);
	  ppgswin(xminshow, xmaxshow, vectorMin2, vectorMax2);   // Set window in other units
	  //	  pgplot_frame.swin_y1 = vectorMin2;
	  //	  pgplot_frame.swin_y2 = vectorMax2;
	  //	  printf("XXXXX pgplotbox2 = %s %s (box=%d labels=%d)\n", pgplot_options.box.box_xopt, pgplot_options.box.box_yopt, pgplot_options.box.drawbox, pgplot_options.box.drawlabels);
	  pgplot_drawbox(pgplotbox);   // Plot only the box + numbers
	  //	  fprintf(stderr, "XXXXXX pgswin7 %f %f %f %f\n", xminshow, xmaxshow, ymin, ymax);
	  ppgswin(xminshow, xmaxshow, ymin, ymax);  // Go back to integer vector numbers
	  //	  ppgswin(x, x2, ymin, ymax);    // Go back to integer vector numbers
	  */
      }

      pgplot_options.box.drawbox = 0;

      ppgsch(label_font.characterheight);
      //      ppgslw(label_font.linewidth);
      //      ppgscf(label_font.font);
      pgplot_options.box.label_f = label_font.font;
      if(curpanelnrx == 0 || appendframes_flag == 0) {
	if(xtitle_set == 0) {
	  if(xUnitsSwitch != XUNIT_BINS) {
	    if(xUnitsSwitch == XUNIT_TIME) {
	      sprintf(xtitle, "Time (sec)");
	    }else if(xUnitsSwitch == XUNIT_PHASE) {
	      if(didtranspose_orig_nrbin == 0)   // Only one subint on horizontal axis
		sprintf(xtitle, "Pulse phase");
	      else
		sprintf(xtitle, "Subint");
	    }else if(xUnitsSwitch == XUNIT_DEG) {
	      sprintf(xtitle, "Pulse longitude (deg)");
	    }
	    // Following made no sense, if yrange is set should not make a difference for the x-label.
	    // Maybe want to include longitudeRangeSet_flag and it depends in xrangeset rather than yrangeset.
	    //	    if((fin.gentype == GENTYPE_2DFS) && (vectorMinMaxSpecified == 0 && fin.yrangeset == 1)) {
	    if((fin.gentype == GENTYPE_2DFS) && (longitudeRangeSet_flag == 0 && fin.xrangeset == 1)) {
	      strcpy(xtitle, "fluctuation frequency (cycles/period)");
	    }else if((fin.gentype == GENTYPE_HRFS_UNFOLDED || fin.gentype == GENTYPE_HRFS) && (longitudeRangeSet_flag == 0 && fin.xrangeset == 1)) {
	      strcpy(xtitle, "fluctuation frequency (cycles/period)");
	    }
	  }else {
	    if(fin.gentype == GENTYPE_PROFILE || fin.gentype == GENTYPE_PULSESTACK || fin.gentype == GENTYPE_SUBINTEGRATIONS || fin.gentype == GENTYPE_DYNAMICSPECTRUM || fin.gentype == GENTYPE_LRFS || fin.gentype == GENTYPE_P3FOLD || fin.gentype == GENTYPE_LRCC || fin.gentype == GENTYPE_PADIST || fin.gentype == GENTYPE_ELLDIST) {
	      sprintf(xtitle, "Pulse longitude (bins)");
	    }else if(fin.gentype == GENTYPE_S2DFSP3 || fin.gentype == GENTYPE_S2DFSP2) {
	      sprintf(xtitle, "Block number (pulses)");
	    }else {
	      sprintf(xtitle, "Bin");
	    }
	  }
	}
	
	//	printf("XXXXX %d %d\n", ytitle_set_to_intensity, stack_state[current_stack_pos-1].grayscalemode);
	if(ytitleset == 0) {
	  if(ytitle_set_to_intensity) {  // This bit of code is to ensure that if in interactive mode, title get changed correctly when switching to/from map mode
	    if(stack_state[current_stack_pos-1].grayscalemode == 0) {
	      strcpy(ytitle, "Intensity");
	    }else {
	      strcpy(ytitle, "Pulse number");
	    }
	  }else if(ytitle_set_to_freqchannel && yUnitsSwitch == 0) {
	    strcpy(ytitle, "Frequency channel");
	  }else if(ytitle_set_to_pulsenumber && yUnitsSwitch == 0) {
	    strcpy(ytitle, "Pulse number");
	    // Unless only one subint is shown, an in line mode
	    if(stack_state[current_stack_pos-1].subint_start == stack_state[current_stack_pos-1].subint_end && stack_state[current_stack_pos-1].grayscalemode == 0) {
	      strcpy(ytitle, "Intensity");
	    }
	  }else if(ytitle_set_to_subint && yUnitsSwitch == 0) {
	    strcpy(ytitle, "Subint number");
	    // Unless only one subint is shown, an in line mode
	    if(stack_state[current_stack_pos-1].subint_start == stack_state[current_stack_pos-1].subint_end && stack_state[current_stack_pos-1].grayscalemode == 0) {
	      strcpy(ytitle, "Intensity");
	    }
	  }else if(ytitle_set_to_p3fold && yUnitsSwitch) {
	    strcpy(ytitle, "Pulse number");
	  }else if(ytitle_set_to_spectral_power) {
	    strcpy(ytitle, "Spectral Power");
	  }else if(ytitle_set_to_harmonic_number) {
	    strcpy(ytitle, "Harmonic Number");
	  }else if(ytitle_set_to_lag_number) {
	    strcpy(ytitle, "Lag number");
	  }else if((ytitle_set_to_pulsenumber || ytitle_set_to_subint) && yUnitsSwitch) {
	    strcpy(ytitle, "Time (sec)");
	  }else if(ytitle_set_to_fluctuationFreq && yUnitsSwitch) {
	    strcpy(ytitle, "fluctuation frequency (cycles/period)");
	  }else if(ytitle_set_to_fluctuationFreq && yUnitsSwitch == 0) {
	    strcpy(ytitle, "fluctuation frequency (bin)");
	  }else if(ytitle_set_to_padist && yUnitsSwitch) {
	    strcpy(ytitle, "PA (deg)");
	  }else if(ytitle_set_to_padist && yUnitsSwitch == 0) {
	    strcpy(ytitle, "PA (bin)");
	  }else if(ytitle_set_to_elldist && yUnitsSwitch) {
	    strcpy(ytitle, "\\gx (deg)");
	  }else if(ytitle_set_to_elldist && yUnitsSwitch == 0) {
	    strcpy(ytitle, "\\gx (bin)");
	  }else if(ytitle_set_to_pulselongitude && yUnitsSwitch) {
	    strcpy(ytitle, "Pulse longitude (deg)");
	  }else if(ytitle_set_to_pulselongitude && yUnitsSwitch == 0) {
	    strcpy(ytitle, "Pulse longitude (bins)");
	  }else if(yUnitsSwitch && didtranspose_orig_nrbin) {
	    strcpy(ytitle, "Frequency (MHz)");
	  }
	}

	pgplot_options.box.drawlabels = 1;
	strcpy(pgplot_options.box.xlabel, xtitle);
	strcpy(pgplot_options.box.ylabel, ytitle);
	//	printf("XXXXX pgplotbox4 = %s %s (box=%d labels=%d)\n", pgplot_options.box.box_xopt, pgplot_options.box.box_yopt, pgplot_options.box.drawbox, pgplot_options.box.drawlabels);
	pgplot_drawbox(&pgplot_options.box);  // Plot only the axis labels
      }else { 
	pgplot_options.box.drawlabels = 1;
	strcpy(pgplot_options.box.xlabel, xtitle);
	strcpy(pgplot_options.box.ylabel, "");
	//	printf("XXXXX pgplotbox5 = %s %s (box=%d labels=%d)\n", pgplot_options.box.box_xopt, pgplot_options.box.box_yopt, pgplot_options.box.drawbox, pgplot_options.box.drawlabels);
	pgplot_drawbox(&pgplot_options.box);   // Plot only the axis labels
      }
    
      if(!stack_state[current_stack_pos-1].grayscalemode) {   
	pgplot_options.box.drawlabels = 0;
	pgplot_options.box.drawtitle = 1;
	newtext = str_replace_header_params(fin, title, application.verbose_state);
	if(newtext == NULL) {
	  printwarning(application.verbose_state.debug, "WARNING pplot: Cannot substitute keyword in title");
	  pgplot_options.box.title[0] = 0;
	}else {
	  strcpy(pgplot_options.box.title, newtext);
	  free(newtext);
	}
	//	printf("XXXXX pgplotbox6 = %s %s (box=%d labels=%d)\n", pgplot_options.box.box_xopt, pgplot_options.box.box_yopt, pgplot_options.box.drawbox, pgplot_options.box.drawlabels);
	//	printf("XXXXXX %f %d %d\n", pgplot_options.box.title_ch, pgplot_options.box.title_lw, pgplot_options.box.title_f);
	pgplot_drawbox(&pgplot_options.box);   // Plot only the title
	pgplot_options.box.drawtitle = 0;
      }


      //START REGION DEVELOP
      if(nrlines > 0) {
	if(outline > 0) {
	  for(i = 0; i < nrlines; i++) {
	    ppgslw(line_lw[i]+outline);
	    ppgsci(outline_color);
	    ppgmove(line_x1[i], line_y1[i]);
	    ppgdraw(line_x2[i], line_y2[i]);      
	  }
	}
	for(i = 0; i < nrlines; i++) {
	  ppgsci(line_ci[i]);
	  ppgslw(line_lw[i]);
	  ppgmove(line_x1[i], line_y1[i]);
	  ppgdraw(line_x2[i], line_y2[i]);      
	}
	ppgslw(1);
	ppgsci(1);
      }
    
      //START REGION DEVELOP	
      if(nrtext > 0) {
	for(i = 0; i < nrtext; i++) {
	  int dofree;
	  ppgsch(text_ch[i]);
	  ppgscf(text_f[i]);
	  newtext = str_replace_header_params(fin, argv[text_txt[i]], application.verbose_state);
	  dofree = 1;
	  if(newtext == NULL) {
	    printwarning(application.verbose_state.debug, "WARNING pplot: Cannot substitute keyword in -text option");
	    newtext = fin.filename;
	    dofree = 0;
	  }
	  if(outline > 0) {
	    ppgslw(text_lw[i]+outline);
	    ppgsci(outline_color);
	    ppgptxt(text_x[i], text_y[i], 0, 0, newtext);
	  }
	  ppgsci(text_ci[i]);
	  ppgslw(text_lw[i]);
	  ppgptxt(text_x[i], text_y[i], 0, 0, newtext); 
	  ppgscf(1);
	  ppgsch(1);
	  ppgslw(1);
	  ppgsci(1);
	  if(dofree)
	    free(newtext);
	}
      }
	
      //START REGION DEVELOP	
      if(markerPosition >= 0) {
	ppgsci(markercolor);
	ppgslw(markerwidth);
	x = stack_state[current_stack_pos-1].x1;
	y = stack_state[current_stack_pos-1].subint_start + markerPosition;
	if(xUnitsSwitch != XUNIT_BINS) {
	  x *= baseline/(float)fin.NrBins;
	  x += dxshift;
	}
	ppgmove(x, y);
	x = stack_state[current_stack_pos-1].x2;
	if(xUnitsSwitch != XUNIT_BINS) {
	  x *= baseline/(float)fin.NrBins;
	  x += dxshift;
	}
	ppgdraw(x, y);
	ppgsci(1);
	ppgslw(1);
      }
      
      //START REGION DEVELOP	
      if(markerPosition2 >= 0) {
	ppgsci(markercolor2);
	ppgslw(markerwidth2);
	x = markerPosition2;
	if(fabs(x-xmax)/(xmax-xmin) < 0.005)
	  x = (xmax-xmin)*0.995+xmin;
	if(fabs(x-xmin)/(xmax-xmin) < 0.005)
	  x = (xmax-xmin)*0.005+xmin;
	ppgmove(x, ymin);
	ppgdraw(x, ymax);
	ppgsci(1);
	ppgslw(1);
      }
      
      //START REGION DEVELOP	
      if(errorbar_set) {
	/*      printf("%f %f %f\n",  errorbar_x, errorbar_y, errorbar_s); */
	ppgslw(plotlw);
	ppgerr1(6, errorbar_x, errorbar_y, errorbar_s*scale, 2.0);
	ppgslw(1);
      }

//START REGION RELEASE

      redraw = 0;
      if(interactive_flag) {
	do { // Corresponds to the while(redraw == 0)
	  int key, key2;
	  FILE *fout;
	  /* If first time entered loop, print help */
	  if(data_read == 2)
	    data_read = 1;
	  if(data_read == 0) {
	    key = '?';
	    data_read = 1;
	  }else {
	    printf("Option: ");
	    fflush(stdout);
	    do {
	      key = pgetch_macro(&application, application.verbose_state);
	    }while(key == '\n' || key == '\r');
	    printf("%c\n", key);
	  }
	  switch(key) {
	  case '?':
	    printf("a       Auto-scale, only works with line plot\n");
	    printf("b       Toggle units of the x-axis between bins and other units\n");
	    printf("I       Info (statistical), limited by current x-axis selection\n");
	    printf("l       Specify left/right bounds, i.e. the x-range\n");
	    printf("M       Change between colour map/line plot mode\n");
	    printf("n       Step in vector range (vertical range)\n");
	    printf("p       Pop stack\n");
	    printf("P       Switch polarization channel\n");
	    printf("q       Quit\n");
	    printf("r       Redraw\n");
	    printf("R       Read cursor position\n");
	    printf("s       Set scaling of data (similar to -scale/-scale2 options).\n");
	    printf("v       Set vector range (vertical range)\n");
	    printf("y       Toggle units of the y-axis\n");
	    if(zapmode == 0)
	      printf("z       Toggle from zap vectors (vertical) to zap subints (horizontal) mode\n");
	    else
	      printf("z       Toggle from zap subints (horizontal) to zap vectors (vertical) mode\n");
	    if(zapmode == 0)
	      printf("Z       Zap vectors (vertical axis, change to subints with 'z'). This only affects the plot, not the input-data.\n");
	    else
	      printf("Z       Zap subints (horizontal axis, change to vectors with 'z'). This only affects the plot, not the input-data.\n");
	    if(change_filename_extension(inputfilename, txt3, "zap", 999, application.verbose_state) == 1) {
	      if(zapmode == 0)
		printf("W       Write out list of zapped vectors to %s (change to subints with 'z')\n", txt3);
	      else
		printf("W       Write out list of zapped subints to %s (change to vectors with 'z')\n", txt3);
	    }
	    printf(",/. Move to the left/right (if zoomed in)\n");
	    printf("?   This help\n");
	    break;
	  case 'a':
	    if(stack_state[current_stack_pos-1].grayscalemode) {
	      printf("Resetting scale to 1.\n");
	      scale = 1;
	      redraw = 1;
	    }else {
	      //	      printf("XXXXXX scale was: %e\n", scale);
	      //	      printf("XXXXXX subint_start=%ld\n", stack_state[current_stack_pos-1].subint_start);
	      //	      printf("XXXXXX subint_end  =%ld\n", stack_state[current_stack_pos-1].subint_end);
	      scale = 1;
	      y = y2 = -1;
	      float float_tmp;
	      for(i = stack_state[current_stack_pos-1].subint_start; i <= stack_state[current_stack_pos-1].subint_end; i++) {
		for(j = stack_state[current_stack_pos-1].x1; j <= stack_state[current_stack_pos-1].x2; j++) {
		  float_tmp = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
		  // As ypos() works differently when only one pulse is shown
		  if(stack_state[current_stack_pos-1].subint_start != stack_state[current_stack_pos-1].subint_end)
		    float_tmp -= i;
		  if(-float_tmp > y || (i == stack_state[current_stack_pos-1].subint_start && j == stack_state[current_stack_pos-1].x1)) {
		    y = float_tmp;
		  }
		  if(float_tmp > y2 || (i == stack_state[current_stack_pos-1].subint_start && j == stack_state[current_stack_pos-1].x1)) {
		    y2 = float_tmp;
		    //		  printf("XXXXX y2=%e\n", y2);
		    //		}else {
		    //		  printf("XXXXX i=%ld y2=%e\n", i, y2);
		  }
		}
	      }
	      //	      printf("XXXXX y =%e\n", y);
	      //	      printf("XXXXX y2=%e\n", y2);
	      if(y2 > 0) {
		scale = fabs(1.0/y2);
		//		printf("XXXXXX scale is: %e\n", scale);
	      }
	      if(fabs(1.0/y) < scale && y > 0) {
		scale = fabs(1.0/y);
		//		printf("XXXXXX scale is updated to: %e\n", scale);
	      }
	      redraw = 1;
	    }
	    break;
	  case 'M':
	    printf("  c = contour map\n");
	    printf("  m = colour map\n");
	    printf("  C = contour + colour map\n");
	    printf("  l = line drawing\n");
	    fflush(stdout);
	    if(current_stack_pos >= max_nr_stack) {
	      printwarning(application.verbose_state.debug, "WARNING: Stack is full");
	      current_stack_pos--;
	    }else {
	      copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
	    }
	    do {
	      key2 = pgetch_macro(&application, application.verbose_state);
	    }while(key2 == '\n' || key2 == '\r');
	    if(key2 == 'm') {
	      printf("set colour map mode\n");
	      stack_state[current_stack_pos].grayscalemode = 1;
	      scale = 1;   // Reset scale, as has a different meaning in grayscale mode
	    }else if(key2 == 'c') {
	      printf("set contour map mode\n");
	      stack_state[current_stack_pos].grayscalemode = 2;
	      printf("\nType the desired nr of contours in the terminal: ");
	      fflush(stdout);
	      scanf("%d", &(stack_state[current_stack_pos].nrcontours));
	      scale = 1;   // Reset scale, as has a different meaning in grayscale mode
	    }else if(key2 == 'C') {
	      printf("set contour + colour map mode\n");
	      stack_state[current_stack_pos].grayscalemode = 3;
	      printf("\nType the desired nr of contours in the terminal: ");
	      fflush(stdout);
	      scanf("%d", &(stack_state[current_stack_pos].nrcontours));
	      scale = 1;   // Reset scale, as has a different meaning in grayscale mode
	    }else {
	      printf("set line plot mode\n");
	      stack_state[current_stack_pos].grayscalemode = 0;
	    }
	    current_stack_pos++;
	    redraw = 1;
	    break;
	  case 'p':
	    current_stack_pos--;
	    if(current_stack_pos <= 0) {
	      printwarning(application.verbose_state.debug, "WARNING: Stack is empty");
	      current_stack_pos = 1;
	    }else {
	      if(stack_state[current_stack_pos].nrZappedVectors != stack_state[current_stack_pos-1].nrZappedVectors || stack_state[current_stack_pos].nrZappedSubints != stack_state[current_stack_pos-1].nrZappedSubints) {
		data_read = 2;   /* This means, do not initialise most variables. Need to reload, since got different zapping applied. */
	      }
	    }
	    redraw = 1;
	    break;
	  case 'P': 
	    current_polnr++;
	    if(current_polnr >= nrpolarizations)
	      current_polnr = 0;
	    printf("Selected polarization channel %d\n", current_polnr);
	    redraw = 1;
	    data_read = 2;   /* This means, do not initialise most variables */
	    break;
	  case 'q':
	  case 3:                     /* Control C */
	    interactive_flag = 0;
	    redraw = 1;
	    break;
	  case 'b':
	    xUnitsSwitch += 1;
	    if(xUnitsSwitch == XUNIT_ENDOFLIST || fin.gentype == GENTYPE_S2DFSP3 || fin.gentype == GENTYPE_S2DFSP2)
	      xUnitsSwitch = 0;
	    //	    printf("XXXX xUnitsSwitch = %d\n", xUnitsSwitch);
	    redraw = 1;
	    break;
	  case 'y':
	    yUnitsSwitch += 1;
	    if(yUnitsSwitch == 2)
	      yUnitsSwitch = 0;
	    redraw = 1;
	    break;
	  case 'r':
	    redraw = 1;
	    break;
	  case 's':
	    if(stack_state[current_stack_pos-1].grayscalemode) {
	      printf("Set scale (default=1, current value %f). A value > 1 results in clipping, making weak features clearer.\n", scale);
	    }else {
	      printf("Set scale (default=1, current value %f). All intensities are multiplied with this value.\n", scale);
	    }
	    printf("Scale = ");
	    fflush(stdout);
	    if(application.macro_ptr == NULL) {
	      scanf("%f", &dummyf1);
	    }else {
	      int ret;
	      ret = fscanf(application.macro_ptr, "%f", &dummyf1);
	      if(ret != 1) {
		fclose(application.macro_ptr);
		application.macro_ptr = NULL;
		printf("\nReached end of macro, switching to keyboard input: type in two numbers\n");
		scanf("%f", &dummyf1);
	      }else {
		printf("%f\n", dummyf1);
	      }
	    }
	    scale = dummyf1;

	    if(stack_state[current_stack_pos-1].grayscalemode) {
	      printf("Set second scale (default=0, current value %f). A value closer to 1 results in clipping of low intesity values to emphasize the bright features.\n", scale2);
	      printf("Scale = ");
	      fflush(stdout);
	      if(application.macro_ptr == NULL) {
		scanf("%f", &dummyf1);
	      }else {
		int ret;
		ret = fscanf(application.macro_ptr, "%f", &dummyf1);
		if(ret != 1) {
		  fclose(application.macro_ptr);
		  application.macro_ptr = NULL;
		  printf("\nReached end of macro, switching to keyboard input: type in two numbers\n");
		  scanf("%f", &dummyf1);
		}else {
		  printf("%f\n", dummyf1);
		}
	      }
	      scale2 = dummyf1;
	    }

	    redraw = 1;
	    break;
	  case 'R':
	    printf("Left click to read positions, press other key to quit.\n");
	    do {
	      int nx, ny;
	      //	      ppgcurs(&x, &y, &ch);
	      ppgband(7, 0, 0, 0, &x, &y, &ch);
	      if(ch == 65) {
		pgplotMapCoordinate(x, y, &nx, &ny);
		//		printf("XXXXXX clicked on vector coordinates %f %f\n", x, y);
		/*
		if(fin.yrangeset && yUnitsSwitch) {
		  if(showTwice_flag)
		    dummyf1 = 2*y*(fin.yrange[1]-fin.yrange[0])/(float)(fin.NrSubints-1)+fin.yrange[0];
		  else
		    dummyf1 = y*(fin.yrange[1]-fin.yrange[0])/(float)(fin.NrSubints-1)+fin.yrange[0];
		}else if(vectorMinMaxSpecified) {
		  dummyf1 = y*(vectorMax-vectorMin)/(float)(fin.NrSubints-1)+vectorMin;
		}else {
		  dummyf1 = y*yUnitCmdLine;
		  }*/
		yvec2unit(fin, yUnitsSwitch, vectorMinMaxSpecified, vectorMin, vectorMax, didtranspose_orig_nrbin, 0, y, &dummyf1, stack_state[current_stack_pos-1].grayscalemode, 0, application.verbose_state);
		printf("%f %f (%d %d)\n", x, dummyf1, nx, ny);
	      }
	    }while(ch == 65);
	    break;
	  case 'z':
	    if(zapmode == 0) {   
	      if(didtranspose_orig_nrbin) {    // Will only work if there are multiple subints
		printf("Swapping to zapping of subints\n");
		zapmode = 1;
	      }else {
		printf("Swapping to zapping of subints only works if there is more than one subint and more than one frequency channel\n");
	      }
	    }else {
	      printf("Swapping to zapping of vectors (vertical axis)\n");
	      zapmode = 0;
	    }
	    break;
	  case 'W':    // W : write out zap file
	    if(zapmode == 0) {
	      if(stack_state[current_stack_pos-1].nrZappedVectors == 0) {
		printwarning(application.verbose_state.debug, "WARNING pplot: No vectors are zapped, so writing out of a zapfile is ignored.");
	      }else {
		int ret;
		if(fin.NrFreqChan == 1 && didtranspose_orig_nrbin == 0) {
		  ret = change_filename_extension(inputfilename, txt3, "subint.zap", 999, application.verbose_state);
		}else {
		  ret = change_filename_extension(inputfilename, txt3, "freq.zap", 999, application.verbose_state);
		}
		if(ret == 1) {
		  fout = fopen(txt3, "w");
		  if(fout != NULL) {  // Write out sorted and uniq list of numbers
		    long largest, count;
		    count = 0;
		    largest = zappedVectors[0];
		    for(i = 0; i < stack_state[current_stack_pos-1].nrZappedVectors; i++) {
		      if(zappedVectors[i] > largest) {
			largest = zappedVectors[i];
		      }
		    }
		    for(j = 0; j <= largest; j++) {
		      for(i = 0; i < stack_state[current_stack_pos-1].nrZappedVectors; i++) {
			if(zappedVectors[i] == j) {
			  fprintf(fout, "%d\n", zappedVectors[i]);
			  count++;
			  break;
			}
		      }
		    }
		    printf("pplot: Written %ld vector numbers to %s\n", count, txt3);
		    fclose(fout);
		  }else {
		    printerror(application.verbose_state.debug, "ERROR pplot: Cannot open zapfile %s.", txt3);
		    return 0;
		  }
		}else {
		  printerror(application.verbose_state.debug, "ERROR pplot: Filename is too long.");
		  return 0;
		}
	      }
	    }else {
	      if(stack_state[current_stack_pos-1].nrZappedSubints == 0) {
		printwarning(application.verbose_state.debug, "WARNING pplot: No subints are zapped, so writing out of a zapfile is ignored.");
	      }else {
		if(change_filename_extension(inputfilename, txt3, "subint.zap", 999, application.verbose_state) == 1) {
		  fout = fopen(txt3, "w");
		  if(fout != NULL) {  // Write out sorted and uniq list of numbers
		    long largest, count;
		    count = 0;
		    largest = zappedSubints[0];
		    for(i = 0; i < stack_state[current_stack_pos-1].nrZappedSubints; i++) {
		      if(zappedSubints[i] > largest) {
			largest = zappedSubints[i];
		      }
		    }
		    for(j = 0; j <= largest; j++) {
		      for(i = 0; i < stack_state[current_stack_pos-1].nrZappedSubints; i++) {
			if(zappedSubints[i] == j) {
			  fprintf(fout, "%d\n", zappedSubints[i]);
			  count++;
			  break;
			}
		      }
		    }
		    printf("pplot: Written %ld subint numbers to %s\n", count, txt3);
		    fclose(fout);
		  }else {
		    printerror(application.verbose_state.debug, "ERROR pplot: Cannot open zapfile %s.", txt3);
		    return 0;
		  }
		}else {
		  printerror(application.verbose_state.debug, "ERROR pplot: Filename is too long.");
		  return 0;
		}
	      }
	    }
	    break;
	  case 'Z':
	    if(current_stack_pos >= max_nr_stack) {
	      printwarning(application.verbose_state.debug, "WARNING: Stack is full");
	      current_stack_pos--;
	    }else {
	      copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
	    }
	    if(zapmode == 0)
	      printf("Left click to zap vector (vertical axis), left click followed by a right click removes a range. Press other key to quit. The result only visible after you're done with clicking.\n");
	    else
	      printf("Left click to zap subint (horizontal axis), left click followed by a right click removes a range. Press other key to quit. The result only visible after you're done with clicking.\n");
	    do {
	      int nx, ny, nx_final, ny_final, ret, ch_value;
	      if(application.macro_ptr == NULL) { // If no macro defined, use the keyboard instead
		//		ppgcurs(&x, &y, &ch);
		if(zapmode == 0)
		  ppgband(5, 0, 0, 0, &x, &y, &ch);
		else
		  ppgband(6, 0, 0, 0, &x, &y, &ch);
	      }else {
		if(application.verbose_state.verbose) {
		  printf("In macro mode, add lines with \"keyvalue xvalue yvalue\"\n");
		  printf("keyvalue = 65 = Left click\n");
		  printf("keyvalue = 88 = Right click\n");
		}
		ret = fscanf(application.macro_ptr, "%d %f %f", &ch_value, &x, &y);
		if(ret != 3) {
		  printf("\nReached end of macro, switching to keyboard input\n");
		  fclose(application.macro_ptr);
		  application.macro_ptr = NULL;
		  ch_value = 27;
		}
		ch = ch_value;
		//		printf("XXXXXX %d %f %f", ch, x, y);
	      }
	      if(ch == 65 || ch == 88) {
		pgplotMapCoordinate(x, y, &nx, &ny);
		if(zapmode == 0) {
		  if(ch == 65) {
		    printf("Zapping vector %d\n", ny);
		    ny_last_clicked_zap = ny;
		    ny_final = ny;
		  }else {
		    if(ny < ny_last_clicked_zap) {
		      int oldint;
		      oldint = ny;
		      ny = ny_last_clicked_zap;
		      ny_last_clicked_zap = oldint;
		    }
		    printf("Zapping vectors %d - %d\n", ny_last_clicked_zap, ny);
		    ny_final = ny;
		  }
		  for(ny = ny_last_clicked_zap; ny <= ny_final; ny++) {
		    if(stack_state[current_stack_pos].nrZappedVectors < max_nr_zap-1) {
		      zappedVectors[stack_state[current_stack_pos].nrZappedVectors++] = ny;
		    }else {
		      printwarning(application.verbose_state.debug, "WARNING: Maximum number of interactively zapped vectors exceeded.");
		    }
		    if(zapVectors(stack_state[current_stack_pos].nrZappedVectors, zappedVectors, subint_start, fin, dataSubset, ny, application.verbose_state) == 0)
		      return 0;
		  }
		}else {  // If zapping subints in a subint/freqchan plot
		  nx /= didtranspose_orig_nrbin;  // Divide by the nr of bins per subint as stored in this variable set at the time of the transpose
		  if(ch == 65) {
		    printf("Zapping subint %d\n", nx);
		    nx_last_clicked_zap = nx;
		    nx_final = nx;
		  }else {
		    if(nx < nx_last_clicked_zap) {
		      int oldint;
		      oldint = nx;
		      nx = nx_last_clicked_zap;
		      nx_last_clicked_zap = oldint;
		    }
		    printf("Zapping subints %d - %d\n", nx_last_clicked_zap, nx);
		    nx_final = nx;
		  }
		  for(nx = nx_last_clicked_zap; nx <= nx_final; nx++) {
		    if(stack_state[current_stack_pos].nrZappedSubints < max_nr_zap-1) {
		      zappedSubints[stack_state[current_stack_pos].nrZappedSubints++] = nx;
		    }else {
		      printwarning(application.verbose_state.debug, "WARNING: Maximum number of interactively zapped subints exceeded.");
		    }
		    if(zapSubints(stack_state[current_stack_pos].nrZappedSubints, zappedSubints, subint_start, didtranspose_orig_nrbin, fin, dataSubset, nx, application.verbose_state) == 0)
		      return 0;
		  }
		}
		nx_final = nx;
		ny_final = ny;
	      }else {
		printf("Got key: %d\n", ch);
	      }
	    }while(ch == 65);
	    redraw = 1;
	    current_stack_pos++;
	    break;
	  case 'v':
	    printf("Specify vector (vertical) range (two numbers): ");
	    fflush(stdout);
	    if(current_stack_pos >= max_nr_stack) {
	      printwarning(application.verbose_state.debug, "WARNING: Stack is full");
	      current_stack_pos--;
	    }else {
	      copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
	    }
	    if(application.macro_ptr == NULL) {
	      scanf("%f %f", &dummyf1, &dummyf2);
	    }else {
	      int ret;
	      ret = fscanf(application.macro_ptr, "%f %f", &dummyf1, &dummyf2);
	      if(ret != 2) {
		printf("\nReached end of macro, switching to keyboard input: type in two numbers\n");
		scanf("%f %f", &dummyf1, &dummyf2);
		fclose(application.macro_ptr);
		application.macro_ptr = NULL;
	      }else {
		printf("%f %f\n", dummyf1, dummyf2);
	      }
	    }
	    if((fin.yrangeset && yUnitsSwitch)  // Special format with defined yrange
	       || (vectorMinMaxSpecified && yUnitsSwitch) // -vecRange is defined
	       || (((fin.gentype == GENTYPE_SUBINTEGRATIONS || fin.gentype == GENTYPE_PULSESTACK) && fin.NrFreqChan == 1 && fin.NrSubints > 1) && yUnitsSwitch) // Observing time should be specified
	       || (yUnitsSwitch && didtranspose_orig_nrbin)  // Frequency in MHz should be specified
	       ) {
	      float rangemin, rangemax;
	      yvec2unit(fin, yUnitsSwitch, vectorMinMaxSpecified, vectorMin, vectorMax, didtranspose_orig_nrbin, 0, dummyf1, &rangemin, stack_state[current_stack_pos-1].grayscalemode, 1, application.verbose_state);
	      yvec2unit(fin, yUnitsSwitch, vectorMinMaxSpecified, vectorMin, vectorMax, didtranspose_orig_nrbin, 0, dummyf2, &rangemax, stack_state[current_stack_pos-1].grayscalemode, 1, application.verbose_state);
	      stack_state[current_stack_pos].subint_start = round(rangemin);
	      stack_state[current_stack_pos].subint_end = round(rangemax);
	    }else {
	      stack_state[current_stack_pos].subint_start = dummyf1/yUnitCmdLine;
	      stack_state[current_stack_pos].subint_end = dummyf2/yUnitCmdLine;
	    }
	    //	    if(fin.yrangeset && yUnitsSwitch) {
	      /*
		y = a*v+b;
		fin.yrange[1] = a*(fin.NrSubints-1) + b;
		fin.yrange[0] = b;
		->;
		fin.yrange[1] = a*(fin.NrSubints-1) + fin.yrange[0];
		(fin.yrange[1]-fin.yrange[0])/(float)(fin.NrSubints-1) = a;
		(y-b)/a = v;
	      */
	    /*	      stack_state[current_stack_pos].subint_start = (dummyf1-fin.yrange[0])*(fin.NrSubints-1)/(fin.yrange[1]-fin.yrange[0]);
	      stack_state[current_stack_pos].subint_end = (dummyf2-fin.yrange[0])*(fin.NrSubints-1)/(fin.yrange[1]-fin.yrange[0]);
	    }else {
	      if(vectorMinMaxSpecified) {
		stack_state[current_stack_pos].subint_start = (dummyf1-vectorMin)*(fin.NrSubints-1)/(vectorMax-vectorMin);
		stack_state[current_stack_pos].subint_end = (dummyf2-vectorMin)*(fin.NrSubints-1)/(vectorMax-vectorMin);
	      }else {
		stack_state[current_stack_pos].subint_start = dummyf1/yUnitCmdLine;
		stack_state[current_stack_pos].subint_end = dummyf2/yUnitCmdLine;
	      }
	      }*/

	    
	    //	    printf("XXXXXXXX %ld %ld\n", stack_state[current_stack_pos].subint_start, subint_start);
	    //	    printf("XXXXXXXX %ld %ld\n", stack_state[current_stack_pos].subint_end, subint_end);
	    if(stack_state[current_stack_pos].subint_start < subint_start) {
	      printwarning(application.verbose_state.debug, "WARNING: set start vector to %ld", subint_start);
	      stack_state[current_stack_pos].subint_start = subint_start;
	    }
	    if(stack_state[current_stack_pos].subint_start > subint_end) {
	      printwarning(application.verbose_state.debug, "WARNING: set end vector to %ld", subint_end);
	      stack_state[current_stack_pos].subint_start = subint_end;
	    }
	    if(stack_state[current_stack_pos].subint_end < subint_start) {
	      printwarning(application.verbose_state.debug, "WARNING: set end vector to %ld", subint_start);
	      stack_state[current_stack_pos].subint_end = subint_start;
	    }
	    if(stack_state[current_stack_pos].subint_end > subint_end) {
	      printwarning(application.verbose_state.debug, "WARNING: set end vector to %ld", subint_end);
	      stack_state[current_stack_pos].subint_end = subint_end;
	    }
	    current_stack_pos++;
	    redraw = 1;
	    break;
	  case 'l':
	    printf("Specify x-range (two numbers, or R to use mouse): ");
	    fflush(stdout);
	    if(current_stack_pos >= max_nr_stack) {
	      printwarning(application.verbose_state.debug, "WARNING: Stack is full");
	      current_stack_pos--;
	    }else {
	      copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
	    }
	    if(application.macro_ptr == NULL) {
	      /*	    scanf("%f %f", &x, &x2); */
	      /*	    scanf("%s %s", txt1, txt2); */
	      fgets(txt3, 990, stdin);
	      i = sscanf(txt3, "%s %s", txt1, txt2);
	      if(i == 1 && strcmp(txt1, "R") != 0) {
		fgets(txt3, 990, stdin);
		sscanf(txt3, "%s", txt2);
	      }
	      if(strcmp(txt1, "R") == 0) {
		printf("Click on left edge\n");
		//		ppgcurs(&x, &y, &ch);
		ppgband(6, 0, 0, 0, &x, &y, &ch);
	      }else {
		sscanf(txt1, "%f", &x);
	      }
	      if(strcmp(txt1, "R") == 0) {
		printf("Click on right edge\n");
		x2 = x;
		//		ppgcurs(&x2, &y, &ch);
		ppgband(4, 0, x, y, &x2, &y, &ch);
	      }else {
		sscanf(txt2, "%f", &x2);
	      }
	    }else {
	      int ret;
	      ret = fscanf(application.macro_ptr, "%f %f", &x, &x2);
	      if(ret != 2) {
		printf("\nReached end of macro, switching to keyboard input: type in two numbers\n");
		scanf("%f %f", &dummyf1, &dummyf2);
		fclose(application.macro_ptr);
		application.macro_ptr = NULL;
	      }
	    }
	    printf("Select %f-%f\n", x, x2);
	    if(xUnitsSwitch != XUNIT_BINS) {
	      x -= dxshift;
	      x /= baseline/(float)fin.NrBins;
	      x2 -= dxshift;
	      x2 /= baseline/(float)fin.NrBins;
	    }
	    stack_state[current_stack_pos].x1 = x;
	    stack_state[current_stack_pos].x2 = x2;
	    if(stack_state[current_stack_pos].x1 < 0) {
	      printwarning(application.verbose_state.debug, "WARNING: set x1 to 0");
	      stack_state[current_stack_pos].x1 = 0;
	    }
	    if(stack_state[current_stack_pos].x1 >= fin.NrBins) {
	      printwarning(application.verbose_state.debug, "WARNING: set x1 to %ld", fin.NrBins-1);
	      stack_state[current_stack_pos].x1 = fin.NrBins-1;
	    }
	    if(stack_state[current_stack_pos].x2 < 0) {
	      printwarning(application.verbose_state.debug, "WARNING: set x2 to 0");
	      stack_state[current_stack_pos].x2 = 0;
	    }
	    if(stack_state[current_stack_pos].x2 >= fin.NrBins) {
	      printwarning(application.verbose_state.debug, "WARNING: set x2 to %ld", fin.NrBins-1);
	      stack_state[current_stack_pos].x2 = fin.NrBins-1;
	    }
	    current_stack_pos++;
	    redraw = 1;
	    break;
	  case '.':                    /* Move to the right */
	    if(current_stack_pos >= max_nr_stack) {
	      printwarning(application.verbose_state.debug, "WARNING: Stack is full");
	      current_stack_pos--;
	    }else {
	      copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
	    }
	    y = stack_state[current_stack_pos].x2 - stack_state[current_stack_pos].x1;
	    x = stack_state[current_stack_pos].x1 + y;
	    x2 = stack_state[current_stack_pos].x2 + y;
	    if(x < 0) {
	      printwarning(application.verbose_state.debug, "WARNING: set x1 to 0");
	      x = 0;
	      x2 = y;
	    }
	    if(x2 >= fin.NrBins) {
	      printwarning(application.verbose_state.debug, "WARNING: set x2 to %ld", fin.NrBins-1);
	      x2 = fin.NrBins-1;
	      x = x2 - y;
	      if(x < 0)
		x = 0;
	    }
	    stack_state[current_stack_pos].x1 = x;
	    stack_state[current_stack_pos].x2 = x2;
	    current_stack_pos++;
	    redraw = 1;
	    break;
	  case ',':                    /* Move to the left */
	    if(current_stack_pos >= max_nr_stack) {
	      printwarning(application.verbose_state.debug, "WARNING: Stack is full");
	      current_stack_pos--;
	    }else {
	      copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
	    }
	    y = stack_state[current_stack_pos].x2 - stack_state[current_stack_pos].x1;
	    x = stack_state[current_stack_pos].x1 - y;
	    x2 = stack_state[current_stack_pos].x2 - y;
	    if(x < 0) {
	      printwarning(application.verbose_state.debug, "WARNING: set x1 to 0");
	      x = 0;
	      x2 = y;
	    }
	    if(x2 >= fin.NrBins) {
	      printwarning(application.verbose_state.debug, "WARNING: set x2 to %ld", fin.NrBins-1);
	      x2 = fin.NrBins-1;
	      x = x2 - y;
	      if(x < 0)
		x = 0;
	    }
	    stack_state[current_stack_pos].x1 = x;
	    stack_state[current_stack_pos].x2 = x2;
	    current_stack_pos++;
	    redraw = 1;
	    break;
	  case 'n':
	    i = stack_state[current_stack_pos-1].subint_end-stack_state[current_stack_pos-1].subint_start+1;
	    j = stack_state[current_stack_pos-1].subint_start;
	    k = stack_state[current_stack_pos-1].subint_end;
	    if(current_stack_pos >= max_nr_stack) {
	      printwarning(application.verbose_state.debug, "WARNING: Stack is full");
	      current_stack_pos--;
	    }else {
	      copystackstate(stack_state[current_stack_pos-1], &stack_state[current_stack_pos]);
	    }
	    stack_state[current_stack_pos].subint_start = j+i;
	    stack_state[current_stack_pos].subint_end = k+i;
	    if(stack_state[current_stack_pos].subint_start < subint_start) {
	      printwarning(application.verbose_state.debug, "WARNING: set vector to %ld", subint_start);
	      stack_state[current_stack_pos].subint_start = subint_start;
	    }
	    if(stack_state[current_stack_pos].subint_start > subint_end) {
	      printwarning(application.verbose_state.debug, "WARNING: set vector to %ld", subint_end);
	      stack_state[current_stack_pos].subint_start = subint_end;
	    }
	    if(stack_state[current_stack_pos].subint_end > subint_end) {
	      printwarning(application.verbose_state.debug, "WARNING: set vector to %ld", subint_end);
	      stack_state[current_stack_pos].subint_end = subint_end;
	    }
	    current_stack_pos++;
	    redraw = 1;
	    break;
	  case 'I':
	    printf("Vector -   min            at binnr  max            at binnr  mean           rms            centroid bin\n");
	    for(i = stack_state[current_stack_pos-1].subint_start; i <= stack_state[current_stack_pos-1].subint_end; i++) {
	      double Imin, Imax, mean, rms, x_centroid;
	      long bin_min, bin_max;
	      mean = 0;
	      rms = 0;
	      Imin = Imax = 0;
	      bin_min = bin_max = 0;
	      x_centroid = 0;
	      for(j = stack_state[current_stack_pos-1].x1; j <= stack_state[current_stack_pos-1].x2; j++) {
		double float_tmp;
		float_tmp = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
		// As ypos() works differently when only one pulse is shown
		if(stack_state[current_stack_pos-1].subint_start != stack_state[current_stack_pos-1].subint_end)
		  float_tmp -= i;
		mean += float_tmp;
		x_centroid += j*float_tmp;
		if(j == stack_state[current_stack_pos-1].x1 || float_tmp > Imax) {
		  Imax = float_tmp;
		  bin_max = j;
		}
		if(j == stack_state[current_stack_pos-1].x1 || float_tmp < Imin) {
		  Imin = float_tmp;
		  bin_min = j;
		}
	      }
	      x_centroid /= (double)mean;
	      mean /= (double)(stack_state[current_stack_pos-1].x2-stack_state[current_stack_pos-1].x1+1);
	      for(j = stack_state[current_stack_pos-1].x1; j <= stack_state[current_stack_pos-1].x2; j++) {
		double float_tmp;
		float_tmp = ypos(dataSubset, j, i, scale, fin.NrBins, stack_state[current_stack_pos-1].subint_start, stack_state[current_stack_pos-1].subint_end, subint_start, yUnitCmdLine, dyshift);
		// As ypos() works differently when only one pulse is shown
		if(stack_state[current_stack_pos-1].subint_start != stack_state[current_stack_pos-1].subint_end)
		  float_tmp -= i;
		rms += (float_tmp-mean)*(float_tmp-mean);
	      }
	      rms /= (float)(stack_state[current_stack_pos-1].x2-stack_state[current_stack_pos-1].x1+1);
	      rms = sqrt(rms);
	      printf("%6ld - %14e %9ld %14e %9ld %14e %14e %14e\n", i, Imin, bin_min, Imax, bin_max, mean, rms, x_centroid);
	    }
	    break;
	  default:
	    printf("Unknown option: '%c' (%d)\nPress '?' for help.\n", key, key);
	  }
	}while(redraw == 0);
      } // End of if(interactive_flag)
      curpanelnrx++;
      if(curpanelnrx == nrpanelsx) {
	curpanelnrx = 0;
	curpanelnry++;
      }
      if(curpanelnry == nrpanelsy)
	curpanelnry = 0;
    }while(interactive_flag);
    if(dataSubset_allocated) {
      free(dataSubset);
      dataSubset_allocated = 0;
    }
    if(stack_poly_allocated) {
      free(xstack_poly);
      free(ystack_poly);
      stack_poly_allocated = 0;
    }
    if(maxy_allocated) {
      free(maxy);
      maxy_allocated = 0;
    }
    closePSRData(&fin, 0, application.verbose_state);

    printf("Finished plotting of: %s\n", inputfilename);
  }   // End of while loop over all input files

  ppgend();
//START REGION DEVELOP
  if(application.ppgplot_name != NULL) {
    pgcloseoutputfile();
  }
//START REGION RELEASE
  free(zappedVectors);
  free(zappedSubints);
  free(stack_state);
  terminateApplication(&application);
  return 0;
} // End of main 
 
float ypos(float *stackI, int bin, long PulseNr, float scale, int NrBins, long pulse_bottom, long pulse_top, long subint_start, float yUnitCmdLine, float dyshift)
{
  long i;
  //  fprintf(stderr, "DEBUG: PulseNr=%ld\n", PulseNr);
  //  fprintf(stderr, "DEBUG: subint_start=%ld\n", subint_start);
  // fprintf(stderr, "DEBUG: NrBins=%d\n", NrBins);
  //fprintf(stderr, "DEBUG: bin=%d\n", bin);
  i = (PulseNr-subint_start)*NrBins+bin;
  //  fprintf(stderr, "DEBUG: index=%ld\n", i);
  /* If only plotting one pulse, don't add pulsenumber to vertical direction */
  if(pulse_bottom != pulse_top)
    return stackI[i]*scale+PulseNr*yUnitCmdLine + dyshift;
  else
    return stackI[i]*scale + dyshift;
}

// Returns 0 on error
int setBaselineParams(datafile_definition fin, float *baseline, float *dxshift, int *xUnitsSwitch, verbose_definition verbose)
{
  double period;
  int ret;
  if(fin.isFolded) {
    ret = get_period(fin, 0, &period, verbose);
    if(ret == 2) {
      printerror(verbose.debug, "ERROR pplot (%s): Cannot obtain period", fin.filename);
      return 0;
    }
  }else {
    ret = 1;
    period = -1;
  }
  if(fin.xrangeset) {   /* This code is repeated later in the code (see lfjgkldfjgklrejkltg) */
    *baseline = (fin.xrange[1]-fin.xrange[0])*(fin.NrBins)/(float)(fin.NrBins-0.999);
    *dxshift = fin.xrange[0];
    *xUnitsSwitch = XUNIT_DEG;
  }else {
    // Unless pulse longitudes are defined, tsamp and period must be known
    if(fin.isFolded != 0 || fin.gentype != GENTYPE_SEARCHMODE) {
      if(period < 0.001) {
	fflush(stdout);
	printwarning(verbose.debug, "pplot: The period does not appear to be set in the header. Consider using the -header option.");
	if(xUnitsSwitch != XUNIT_BINS) {
	  printerror(verbose.debug, "       Terminating.");
	  return 0;
	}else {
	  printwarning(verbose.debug, "       (warning only)");
	}
      }
    }
    if((get_tsamp(fin, 0, verbose) < 0.0000001 || get_tsamp(fin, 0, verbose) >= 100)) {
      fflush(stdout);
      printwarning(verbose.debug, "pplot: The sampling time does not appear to be set correctly in the header. Consider using the -header option");
      if(xUnitsSwitch != XUNIT_BINS) {
	printerror(verbose.debug, "       Terminating.");
	return 0;
      }else {
	printwarning(verbose.debug, "       (warning only)");
      }
    }
    if(fin.isFolded != 0 || fin.gentype != GENTYPE_SEARCHMODE)
      *baseline = 360.0*get_tsamp(fin, 0, verbose)*fin.NrBins/period;
    else
      *baseline = 360.0;
    if(*xUnitsSwitch == XUNIT_PHASE) {
      if(fin.isFolded != 0 || fin.gentype != GENTYPE_SEARCHMODE)
	*baseline = get_tsamp(fin, 0, verbose)*fin.NrBins/period;
      else
	*baseline = 1.0;
    }else if(*xUnitsSwitch == XUNIT_TIME) {
      *baseline = get_tsamp(fin, 0, verbose)*fin.NrBins;
    }
    if(verbose.verbose && xUnitsSwitch != XUNIT_BINS) {
      printf("Based on sampling time and pulse period the baseline appears to be %f.\n", *baseline);
      //      fprintf(stderr, "XXXXX setBaselineParams: baseline = %f (nrbins=%ld, nrsub=%ld, nrfreq=%ld)\n", *baseline, fin.NrBins, fin.NrSubints, fin.NrFreqChan);
      //		printf("DEBUG: xUnitsSwitch=%d\n", xUnitsSwitch);
    }
  }
  return 1;
}

// If onlysinglechannel is non-negative, zap that particular channel rather than look in the zappedVectors list
// Return 0 on error
int zapVectors(int nrZappedVectors, int *zappedVectors, long subint_start, datafile_definition fin, float *dataSubset, int onlysinglechannel, verbose_definition verbose)
{
  int i, nx, ny;
  float I;
  if(onlysinglechannel >= 0)
    nrZappedVectors = 1;
  for(i = 0; i < nrZappedVectors; i++) {
    I = 0;
    if(onlysinglechannel < 0)
      ny = zappedVectors[i];
    else
      ny = onlysinglechannel;
    for(nx = 0; nx < fin.NrBins; nx++) {
//START REGION DEVELOP
      /* dataSubset is used for determining the levels,
	 while fin.data is used to plot the data. So the
	 vector has to be zapped in both! Not sure why
	 there is both a dataSubset and a fin.data. It might
	 be necessary, or it might have been in the
	 past. */
//START REGION RELEASE
      dataSubset[(ny-subint_start)*fin.NrBins + nx] = 0;
      if(writePulsePSRData(&fin, ny, 0, 0, nx, 1, &I, verbose) != 1) {
	printerror(verbose.debug, "ERROR pplot: Error writing data while zapping");
	return 0;
      }
    }
  }
  return 1;
}

// If there are freqs and subints, the subints are zapped separately, as they are not stored as vectors, which are now frequency channls.
// If onlysinglesubint is non-negative, zap that particular subint rather than look in the nrZappedSubints list
// Return 0 on error
int zapSubints(int nrZappedSubints, int *zappedSubints, long subint_start, int didtranspose_orig_nrbin, datafile_definition fin, float *dataSubset, int onlysinglesubint, verbose_definition verbose)
{
  float I; 
  long binnr;
  int i, nx, ny;
  if(didtranspose_orig_nrbin) {
    if(onlysinglesubint >= 0)
      nrZappedSubints = 1;
    for(i = 0; i < nrZappedSubints; i++) {
      I = 0;
      if(onlysinglesubint < 0)
	nx = zappedSubints[i];
      else
	nx = onlysinglesubint;
      for(binnr = nx*didtranspose_orig_nrbin; binnr < (nx+1)*didtranspose_orig_nrbin; binnr++) {
	for(ny = 0; ny < fin.NrSubints; ny++) {
	  dataSubset[(ny-subint_start)*fin.NrBins + binnr] = 0;
	  if(writePulsePSRData(&fin, ny, 0, 0, binnr, 1, &I, verbose) != 1) {
	    printerror(verbose.debug, "ERROR pplot: Error writing data while zapping");
	    return 0;
	  }
	}
      }
    }
  }
  return 1;
}

// Type = 0 = Vector centre
//        1 = Bottom of vector
//        2 = Top of vector
// If inverse is set, the vector indicated in proper units is converted into a vector nr.
void yvec2unit(datafile_definition fin, int yUnitsSwitch, int vectorMinMaxSpecified, float vectorMin, float vectorMax, int didtranspose_orig_nrbin, int type, float valuein, float *valueout, int mapmode, int inverse, verbose_definition verbose)
{
  //  printf("XXXXX %d\n", yUnitsSwitch);
  if(vectorMinMaxSpecified == 0) {
    if(fin.yrangeset && yUnitsSwitch) {
      vectorMin = fin.yrange[0];   // These are the centres of the vectors. Simulate -vecRange being used
      vectorMax = fin.yrange[1];
      vectorMinMaxSpecified = 1;
    }else if(((fin.gentype == GENTYPE_SUBINTEGRATIONS || fin.gentype == GENTYPE_PULSESTACK) && fin.NrFreqChan == 1 && fin.NrSubints > 1) && yUnitsSwitch && didtranspose_orig_nrbin == 0) {
      vectorMin = 0;   // These are the centres of the vectors. Simulate -vecRange being used
      vectorMax = get_tobs(fin, verbose);
      vectorMinMaxSpecified = 1;
    }else if(yUnitsSwitch && didtranspose_orig_nrbin) {  // Frequency on axis
      //    printf("XXXXX %ld %ld\n", fin.NrFreqChan, fin.NrSubints);
      //    vectorMin = fin.freq_cent - 0.5*fin.bw + (0+0.5)*fin.bw/(double)fin.NrSubints;
      //    vectorMax = fin.freq_cent - 0.5*fin.bw + (0+0.5)*fin.bw/(double)fin.NrSubints;
      if(fin.freqMode != FREQMODE_UNIFORM) {
	fflush(stdout);
	printwarning(verbose.debug, "WARNING pplot: Frequency range is for first subint.");
      }
      vectorMin = get_nonweighted_channel_freq(fin, 0, verbose);
      vectorMax = get_nonweighted_channel_freq(fin, fin.NrSubints-1, verbose);
      vectorMinMaxSpecified = 1;
      //    printf("XXXXX %f %f\n", vectorMin, vectorMax);
    }
  }
  if(vectorMinMaxSpecified) {  // Note: specified range corresponds to bin centres
    // Note that the bottom/top of vector-range is divided over N-1 vectors
    float binsize;
    binsize = (vectorMax - vectorMin)/(float)(fin.NrSubints-1);
    if(inverse == 0) {
      *valueout = valuein*binsize + vectorMin;
      if(mapmode) {
	if(type == 1)
	  *valueout -= 0.5*binsize;   // There are nrpulses-1 full bins spanning this range
	else if(type == 2)
	  *valueout += 0.5*binsize;   // There are nrpulses-1 full bins spanning this range
      }
    }else {
      if(mapmode) {
	if(type == 1)
	  valuein += 0.5*binsize;
	else if(type == 2)
	  valuein -= 0.5*binsize;
      }
      *valueout = (valuein - vectorMin)/binsize;
    }
  }else {
    *valueout = valuein;
  }
  //  printf("XXXXX valueout=%f (valuein=%f type=%d inverse=%d mapmode=%d)\n", *valueout, valuein, type, inverse, mapmode);
}
//START REGION DEVELOP
