/* Note that Gaussian blobs are Gaussian in the map projection, not in
   reality. In the map projection the opposite magnetic pole is not a
   point, but a circle with a radius of pi. Therefore in the most
   extreme case, a gaussian blob at the other magnetic pole should
   look like a ring in the used map projection, but it now looks like
   a gaussian blob.*/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_sf_legendre.h"
#include "gsl/gsl_sf_ellint.h"
#include "psrsalsa.h"

#define MaxNrCones     100            /* The maximum nr of emission cones */
/* Nr of int. steps to simulate variable carousel rotation */
#define NrCourouselIntegrationSteps   100       
#define MaxNrJumps     100
#define NrEllIntegralInterpol  1000   // The elliptical integral is evaluated this number of times for interpolation table

void SHOWREVISIONINFO_prog() {
#include "fakePulsar.c.svninfo"
}
extern void (*SHOWREVISIONINFO)(void);

psrsalsaApplication application;

/* Used for nummerical integration of carousel rotation phase */
double carouselRotation_last_t;    
double carouselRotation_last_cone_phi[MaxNrCones];
double carouselRotation_last_wave_phi[MaxNrCones];
double carouselRotation_last_oscillation_phi[MaxNrCones];
int reset_num_integration_p4;

typedef struct {
  long NrPulses, NrBins, NrChannels;
  float tstart, frequency, bw;
  float l1, l2;             /* Pulse longitude range to store [deg] */
}OutputOptionsStruct;

typedef struct {
  float rho, halfwidth;     /* Half openings angle carousel and beamlet [rad] */
  float I;                  /* The intensity */
  float P4;                 /* Rotation period carousel [sec] */
  float phi0;               /* Start (t=0) rotation angle [rad] */
  float phi;                /* Current rotation angle is this value plus phi0*/
  float fP4, t4, phi4;      /* The fraction which which P4 oscillates,
			       the timescale [sec] and offset [rad] */
  float e;                  // Ellipticity of the beam
  float or;                 // Orientation of the ellipse, in radians
  int spacing;              // 0=sparks move with indentical angular velocity, 1=equidistant
  int argc;                 /* The argument of the list with intensities of 
			       the individual spark. (0=all 1) */
  char *map_filename;       // Set to NULL if a carousel is simulated, rather than read in from a polar map
  char *p3fold_filename;    // Set to NULL is a carousel is simulated, rather than read in from a p3 fold
  datafile_definition map;  // The actual polar map data, if used.
  char **argv;              /* The arguments as used if argc != 0 */
  double *lookupE;          // Elliptical integral lookup table
  int nrsparks;
}ConeDefinition;

typedef struct {
  float lambda;              /* The wavelength [rad] */
  float I;                  /* The intensity */
  float P4;                 /* Rotation period carousel [sec] */
  float phi0;               /* Start (t=0) rotation angle [rad] */
  float phi;                /* Current rotation angle is this value plus phi0*/
  float fP4, t4, phi4;      /* The fraction which which P4 oscillates, the timescale [sec] and offset [rad] */
}WaveDefinition;

typedef struct {
  int l, m;                 /* The l and m value of the spherical harmonic */
  float I;                  /* The intensity */
  float P4;                 /* Rotation period carousel [sec] */
  float phi0;               /* Start (t=0) rotation angle [rad] */
  float phi;                /* Current rotation angle is this value plus phi0*/
  float fP4, t4, phi4;      /* The fraction which which P4 oscillates, the timescale [sec] and offset [rad] */
  int alwayspos;            /* Set to 1 rather than 0 to make signal always positive */
}OscillationDefinition;

typedef struct {
  float I;                  /* The intensity */
  float x, y;               /* Position [rad] */
  float fwhm;               /* The full width half max [rad] */
}PatchDefinition;

typedef struct {
  float I;                  /* The intensity */
  float fwhm;               /* The full width half max [rad] */
}TopHatDefinition;

typedef struct {
  float I;                  /* The intensity */
  float separation;         /* separation [rad] */
  float fwhm;               /* The full width half max [deg] */
  float phi0;               /* Start phase of the microstructure */
}MicroStructureDefinition;

typedef struct {
  float p3;                 /* P3 value in pulse periods */
  float p2;                 /* P2 value in phase */
}SinusoidDefinition;

typedef struct {
  float alpha, zeta, beta, rm;
  float P1, topHatWidth;
  int nrcones, nrwaves, nrpatches, nroscillations, nrmicrostructure, nrsinusoids, nrtophats;
  ConeDefinition cone[MaxNrCones];
  WaveDefinition wave[MaxNrCones];
  PatchDefinition patch[MaxNrCones];
  OscillationDefinition oscillation[MaxNrCones];
  MicroStructureDefinition microstructure[MaxNrCones];
  SinusoidDefinition sinusoid[MaxNrCones];
  TopHatDefinition tophat[MaxNrCones];
}PulsarDefinitionStruct;

int showmap(PulsarDefinitionStruct *pulsar, int nrx, int nry, int nrfieldlines, float markphase, char *mapdevice, int writefile, int projection, float projection_long_rot, float projection_lat_rot);
/* int openPulseStack(char *filename, FILE **fstack, OutputOptionsStruct outoptions, PulsarDefinitionStruct pulsar); */
float CalculateEmission(float t, PulsarDefinitionStruct *pulsar, OutputOptionsStruct outputOpt, verbose_definition verbose);
float CalculateP3foldEmission(PulsarDefinitionStruct *pulsar, long nb, long np, verbose_definition verbose);
void calculateLOS(int reqOutput, int IP, float alpha_input, float beta_input, float deltaPulseLongitude1, float deltaPulseLongitude2, float trackheight);

int main(int argc, char **argv)
{
  long i, np, nb, nf, npol, idx, nidx, idnum, birdie_number;
  int showmap_flag, write_flag, showpulsstack_flag, lrfs_flag, maciej_flag, fft_size, nrfieldlines, nrJumps;
  int clip_lrfs_flag, showprofile_flag, slope_flag, IP_flag, los_flag, index, projection, dorvm, overwritePA, nronpulsebins, nonormalise;
  int dobirdie, *birdie_list;
  datafile_definition fstack;
  char mapdevice[1000], stackdevice[1000], lrfsdevice[1000], trackdevice[1000], profiledevice[1000], *ppgplot_name;
  float *Istack, *Iprofile, t, phase0, I, I2, s2n, s2n_pp, *lrfs, *phase_track, markphase, los_x1, los_x2, los_height, projection_long_rot, projection_lat_rot, pa, L, phi, phi0, scatter_timescale, scatter_reffreq, scatter_pwrlaw, overwritealpha, overwritebeta, overwritepa0, overwritel0, emissionheightBCW, emissionheightPWRLAW, birdie_I, mark_noise_pulses_amplitude;
  float jump_longitude[MaxNrJumps], jump_offset[MaxNrJumps];
  double Itot, normfac;
  pgplot_options_definition pgplot_options;

  OutputOptionsStruct outputOpt;
  PulsarDefinitionStruct pulsar;

  initApplication(&application, "fakePulsar", "[options]");
  SHOWREVISIONINFO = SHOWREVISIONINFO_prog;
  application.switch_verbose = 1;
  application.switch_debug = 1;
  application.switch_showrevision = 1;
  application.switch_nocounters = 1;
  application.switch_oformat = 1;
  application.switch_formatlist = 1;
  application.switch_history_cmd_only = 1;
  application.switch_fixseed = 1;
  /* Default output format is PSRFITS */
  application.oformat = FITS_format;


  /* Set default values */
  los_flag = 0;
  outputOpt.NrPulses = 1024;
  outputOpt.NrBins = 300;
  outputOpt.NrChannels = 1;
  outputOpt.tstart = 0;
  outputOpt.frequency = 9999;
  outputOpt.bw = 1;
  outputOpt.l1 = 0;
  outputOpt.l2 = 360;
  pulsar.alpha = 9;
  pulsar.beta = 4.5;
  pulsar.P1 = 1.0;
  pulsar.rm = 0;
  pulsar.nrcones = 0;
  pulsar.cone[0].I = 1.0;
  pulsar.cone[0].rho = 4.5;
  pulsar.cone[0].halfwidth = 1.0;
  pulsar.cone[0].nrsparks = 9;
  pulsar.cone[0].phi = 0;
  pulsar.cone[0].phi0 = 0;
  pulsar.cone[0].P4 = 1000.0;
  pulsar.cone[0].fP4 = -1.0;
  pulsar.cone[0].t4 = 1000.0;
  pulsar.cone[0].phi4 = 0.0;
  pulsar.cone[0].argc = 0;
  pulsar.cone[0].e = 0;
  pulsar.cone[0].or = 0;
  pulsar.cone[0].spacing = 0;
  pulsar.cone[0].map_filename = NULL;
  pulsar.cone[0].p3fold_filename = NULL;
  pulsar.nrwaves = 0;
  pulsar.wave[0].I = 1.0;
  pulsar.wave[0].lambda = 3.0;
  pulsar.wave[0].phi = 0;
  pulsar.wave[0].phi0 = 0;
  pulsar.wave[0].P4 = 200.0;
  pulsar.wave[0].fP4 = -1.0;
  pulsar.wave[0].t4 = 1000.0;
  pulsar.wave[0].phi4 = 0.0;
  pulsar.nrpatches = 0;
  pulsar.patch[0].I = 1.0;
  pulsar.patch[0].x = 0.0;
  pulsar.patch[0].y = 0.0;
  pulsar.patch[0].fwhm = 20.0;
  pulsar.nroscillations = 0;
  pulsar.oscillation[0].l = 10;
  pulsar.oscillation[0].m = 0;
  pulsar.oscillation[0].I = 1.0;
  pulsar.oscillation[0].P4 = 0.1;
  pulsar.oscillation[0].phi0 = 0.0;
  pulsar.oscillation[0].phi = 0.0;
  pulsar.oscillation[0].fP4 = -1.0;
  pulsar.oscillation[0].t4 = 100.0;
  pulsar.oscillation[0].phi4 = 0.0;
  pulsar.microstructure[0].I = 1.0;
  pulsar.microstructure[0].separation = 10.0;
  pulsar.microstructure[0].fwhm = 3.0;
  pulsar.microstructure[0].phi0 = 0.0;
  pulsar.nrsinusoids = 0;
  pulsar.sinusoid[0].p3 = 10;
  pulsar.sinusoid[0].p2 = 0.1;
  pulsar.nrtophats = 0;
  pulsar.nrmicrostructure = 0;
  pulsar.tophat[0].I = 1.0;
  pulsar.tophat[0].fwhm = 180.0;

  fft_size = 512;
  showprofile_flag = 0;
  showmap_flag = 0;
  showpulsstack_flag = 0;
  write_flag = 0;
  maciej_flag = 0;
  clip_lrfs_flag = 0;
  sprintf(mapdevice, "?");
  sprintf(stackdevice, "?");
  sprintf(lrfsdevice, "?");
  sprintf(trackdevice, "?");
  sprintf(profiledevice, "?");
  phase0 = 180.0;
  s2n = -1;
  s2n_pp = -1;
  lrfs_flag = 0;
  slope_flag = 0;
  mark_noise_pulses_amplitude = -1;
  ppgplot_name = 0;
  nrfieldlines = 0;
  markphase = -999999;
  dorvm = 0;
  scatter_timescale = -1;
  overwritePA = 0;
  nrJumps = 0;
  emissionheightBCW = 0;
  emissionheightPWRLAW = 0;
  dobirdie = 0;
  nonormalise = 0;
  pgplot_clear_options(&pgplot_options);
  /*

Spherical harmonics: m=0 : It is not a sine wave! It has the largest
amplitude at the origin.

  float theta;
  for(i = 0; i < 1000; i++) {
    theta = i*0.001*2*M_PI;
    I = plgndr(10,1,cos(theta));
    printf("%f %f\n", theta, I);
  }
  return 0;
  */

  if(argc < 2) {
    printf("Program to generate artificial folded single pulse data.\n");
    printApplicationHelp(&application);
    printf("Data options:\n");
    printf("  -gg             Write data to this file.\n");
    printf("  -np             Number of pulses [def=%ld].\n", outputOpt.NrPulses);
    printf("  -nb             Number of bins per pulse. Default is [%ld].\n", outputOpt.NrBins);
    printf("  -nf             Number of frequency channels. Default is [%ld].\n", outputOpt.NrChannels);
    printf("  -freq           \"freq bw\" Set the frequency and bandwidth [%.1f,%.1f].\n", outputOpt.frequency, outputOpt.bw);
    printf("  -l              Set pulse longitude range covered in output [%.1f - %.1f deg].\n", outputOpt.l1, outputOpt.l2);
    printf("\nGraphics options:\n");
    printf("  -prof           Show pulse profile.\n");
    printf("  -profd          Set pgplot device for the pulse profile.\n");
    printf("  -stack          Show pulse-stack.\n");
    printf("  -stackd         Set pgplot device for the pulse-stack.\n");
    printf("  -map            Show emission map with \"type long lat\".\n");
    printf("                     type=0: Zoomed in map\n");
    printf("                     type=1: Hammer-Aitoff\n");
    printf("                     type=2: Spherical.\n");
    printf("                     long and lat allows you to rotate the unzoomed maps.\n");
    printf("  -mapd           Set pgplot device for the emission map.\n");
    printf("  -mapnrB         Show this number of fieldlines.\n");
    printf("  -mapmarkphase   Mark this pulse longitude in degrees in the map.\n");
    printf("  -lrfs           Show LRFS.\n");
    printf("  -clip0          Cut out the first spectral channel of the LRFS.\n");
    printf("  -lrfsd          Set pgplot device for the LRFS.\n");
    printf("  -trackd         Set pgplot device for the subpulse phase.\n");
    printf("  -ppgplot        Output pgplot commands to this file.\n");
    printf("\nGeometry options:\n");
    printf("  -a              Angle alpha between P and B axis [%.2f deg].\n", pulsar.alpha);
    printf("  -b              Angle beta between B axis and line of sight [%.2f deg].\n", pulsar.beta);
    printf("  -P              Rotation period of the pulsar [%.5f sec].\n", pulsar.P1);
    printf("  -phase0         Phase of the peak of the profile  [%.2f deg].\n", phase0);
    printf("\nSimulation options:\n");
    printf("  -rvm            Generate polarization data with 100%% linear polarization\n");
    printf("                  described by the rotating vector model.\n");
    printf("  -paswing        Like -rvm, but use this alpha, beta, pa0 and l0 value instead.\n");
    printf("                  The angles are in degrees, and l0 is w.r.t. -phase0 value.\n");
    printf("  -opm            Put OPM at this longitude and with this amount of degrees when\n");
    printf("                  using -paswing or -rvm. You can use this option multiple times\n"); 
    printf("  -bcw            Add a BCW shift to the PA swing. Give two numbers: the\n");
    printf("                  emission height (as fraction of RLC) at the centre freq and\n"); 
    printf("                  the powerlaw index of the frequency dependence.\n"); 
    printf("  -rm             Add Faraday rotation with this RM, used together with -rvm.\n");
    printf("  -scatter        \"time reffreq index\" Convolve the data with a scatter tail\n");
    printf("                  with e-fold timescale in sec at the reference frequency in\n");
    printf("                  MHz, with a powerlaw dependence on frequency with index\n");
    printf("                  index (i.e. -4). Each pulse is treated seperately, so the\n");
    printf("                  end of scatter tail ends up at the start of the same\n"); 
    printf("                  pulse/subint.\n"); 
    printf("  -los1           \"IP_flag long\" Output the footprint parameters as function of\n");
    printf("                  emission height. IP_flag=0 (MP) or 1 (IP). long is with\n");
    printf("                  respect to the location steepest gradient PA-swing.\n");
    printf("  -los2           \"IP_flag long1 long2 height\" Output the footprint parameters\n");
    printf("                  of the line of sight between two pulse longitudes.\n");
    printf("                  IP_flag=0 (MP) or 1 (IP). long is with respect to location\n");
    printf("                  steepest gradient PA-swing and the emission height is a\n");
    printf("                  fraction with the light cylinder radius.\n");
    printf("\nAdd emission regions:\n");
    printf("\nDefine emission cone (can use this option more than once):\n");
    printf("  -cone 'I rho w n p4 phi0 fP4 t4 phi4'\n");
    printf("  (example: -cone '%.2f %.2f %.2f %d %.2f %.2f %.2f %.2f %.2f')\n", pulsar.cone[0].I, pulsar.cone[0].rho, pulsar.cone[0].halfwidth, pulsar.cone[0].nrsparks, pulsar.cone[0].P4, pulsar.cone[0].phi0, pulsar.cone[0].fP4, pulsar.cone[0].t4, pulsar.cone[0].phi4);
    printf("    I:    The intensity of the beam [%.2f].\n", pulsar.cone[0].I);
    printf("    rho:  half openings angle of the beam [%.2f deg].\n", pulsar.cone[0].rho);
    printf("    w:    half openings angle of the beamlets (where power is half) [%.2f deg].\n", pulsar.cone[0].halfwidth);
    printf("    n:    number of sparks [%d]. Zero means no drifting subpulses.\n", pulsar.cone[0].nrsparks);
    printf("    p4:   rotation period carousel [%.2f sec].\n", pulsar.cone[0].P4);
    printf("    phi0: start azimuthal angle of carousel [%.2f deg].\n", pulsar.cone[0].phi0);
    printf("    fP4:  fraction with which P4 varies [%.2f] (negative is disable).\n", pulsar.cone[0].fP4);
    printf("    t4:   timescale on which P4 varies [%.2f sec].\n", pulsar.cone[0].t4);
    printf("    phi4: start offset of P4 variation cycle [%.2f deg].\n", pulsar.cone[0].phi4);
    printf("\nChange cone into an ellipse (specify directly AFTER each of the -cone commands):\n");
    printf("  -ellipse 'e or spacing'\n");
    printf("    e:  The eccentricity of the ellipse [%.2f].\n", pulsar.cone[0].e);
    printf("    or: The orientation of the semi-major axis (indicated by rho) [%.2f deg].\n", pulsar.cone[0].or);
    printf("    spacing: 0=sparks move with indentical angular velocity, 1=equidistant\n");
    printf("\nChange the relative intensity of the beamlets in the carousel (specify directly AFTER each of the -cone commands):\n");
    printf("  -coneI 'I1 I2 .... IN'\n");
    printf("    Ii:  The relative intensity of beamlet i.\n");
    printf("\nDefine a polar map (generated by pfold for example) read in from file with name filename:\n");
    printf("  -readmap 'filename I p4 phi0 fP4 t4 phi4'\n");
    printf("    I:    The recorded intensities in the map are multiplied with this value [%.2f].\n", pulsar.cone[0].I);
    printf("    p4:   rotation period carousel [%.2f sec].\n", pulsar.cone[0].P4);
    printf("    phi0: start azimuthal angle of carousel [%.2f deg].\n", pulsar.cone[0].phi0);
    printf("    fP4:  fraction with which P4 varies [%.2f] (negative is disable).\n", pulsar.cone[0].fP4);
    printf("    t4:   timescale on which P4 varies [%.2f sec].\n", pulsar.cone[0].t4);
    printf("    phi4: start offset of P4 variation cycle [%.2f deg].\n", pulsar.cone[0].phi4);
    printf("\nDefine a p3-fold (generated by pfold for example) read in from file with name filename:\n");
    printf("  -readp3fold 'filename I phi0'\n");
    printf("    I:    The recorded intensities in the map are multiplied with this value [%.2f].\n", pulsar.cone[0].I);
    printf("    phi0: start phase in the modulation cycle [%.2f deg].\n", pulsar.cone[0].phi0);
    printf("    No interpolation is done: each pulse will be a row in the P3 fold.\n");
    printf("\nDefine emission waves:\n");
    printf("  -wave 'I L p4 phi0 fP4 t4 phi4'\n");
    printf("    I:    The intensity of the beam [%.2f].\n", pulsar.wave[0].I);
    printf("    L:    The wavelength [%.2f deg].\n", pulsar.wave[0].lambda);
    printf("    p4:   wave period [%.2f sec].\n", pulsar.wave[0].P4);
    printf("    phi0: start phase of wave [%.2f deg].\n", pulsar.wave[0].phi0);
    printf("    fP4:  fraction with which P4 varies [%.2f] (negative is disable).\n", pulsar.wave[0].fP4);
    printf("    t4:   timescale on which P4 varies [%.2f sec].\n", pulsar.wave[0].t4);
    printf("    phi4: start offset of P4 variation cycle [%.2f deg].\n", pulsar.wave[0].phi4);
    printf("\nDefine non-radial pulsation:\n");
    printf("  -oscillation  'l m I p4 phi0 fP4 t4 phi4'\n");
    printf("  -oscillation2 'l m I p4 phi0 fP4 t4 phi4'\n");
    printf("  (example: -oscillation '%d %d %.2f %.2f %.2f %.2f %.2f %.2f')\n", pulsar.oscillation[0].l, pulsar.oscillation[0].m, pulsar.oscillation[0].I, pulsar.oscillation[0].P4, pulsar.oscillation[0].phi0, pulsar.oscillation[0].fP4, pulsar.oscillation[0].t4, pulsar.oscillation[0].phi4);
    printf("    Here -oscillation2 forces signal to be always possitive by adding an offset.\n");    
    printf("    l:    The l of the spherical harmonic [%d].\n", pulsar.oscillation[0].l);
    printf("    m:    The m of the spherical harmonic [%d].\n", pulsar.oscillation[0].m);
    printf("    I:    The intensity of the beam [%.2f].\n", pulsar.oscillation[0].I);
    printf("    p4:   oscillation period [%.2f sec].\n", pulsar.oscillation[0].P4);
    printf("    phi0: start phase of oscillation [%.2f deg].\n", pulsar.oscillation[0].phi0);
    printf("    fP4:  fraction with which P4 varies [%.2f] (negative is disable).\n", pulsar.oscillation[0].fP4);
    printf("    t4:   timescale on which P4 varies [%.2f sec].\n", pulsar.oscillation[0].t4);
    printf("    phi4: start offset of P4 variation cycle [%.2f deg].\n", pulsar.oscillation[0].phi4);
    printf("\nDefine simple 2D sinusoid wave pattern (superimposed on other emission):\n");
    printf("  -2dsinusoid 'p3 p2'   (not shown on map)\n");
    printf("  (example: -2dsinusoid '%.2f %.2f'\n", pulsar.sinusoid[0].p3, pulsar.sinusoid[0].p2);
    printf("    p3:   The P3 value in pulse number [%.2f].\n", pulsar.sinusoid[0].p3);
    printf("    p2:   The P2 value in phase [%.2f].\n", pulsar.sinusoid[0].p2);
    printf("\nDefine emission patch (modulation pattern is convolved with the patches):\n");
    printf("  -patch 'I x y fwhm'\n");
    printf("    I:    The intensity of the patch [%.2f].\n", pulsar.patch[0].I);
    printf("    x,y:  The position of the patch [(%.2f,%.2f) deg].\n", pulsar.patch[0].x, pulsar.patch[0].y);
    printf("    fwhm: The full width half max [%.2f deg].\n", pulsar.patch[0].fwhm);
    printf("\nDefine a top-hat profile (modulation pattern is convolved with this top-hat\n");
    printf("function, not shown on map):\n");
    printf("  -tophat 'I fwhm'\n");
    printf("    I:    The intensity of the patch [%.2f].\n", pulsar.tophat[0].I);
    printf("    fwhm: The full width half max [%.2f deg].\n", pulsar.tophat[0].fwhm);
    printf("\nDefine microstructure:\n");
    printf("  -micro 'I sep fwhm'\n");
    printf("    I:    The intensity of the microstructure [%.2f].\n", pulsar.microstructure[0].I);
    printf("    sep:  The separation of the microstructure [%.2f deg].\n", pulsar.microstructure[0].separation);
    printf("    fwhm: The full width half max [%.2f deg].\n", pulsar.microstructure[0].fwhm);
    printf("\nOther:\n");
    printf("  -s2n        Set average S/N ratio of the single pulses per freq. channel.\n");
    printf("  -s3n        Set S/N ratio of the pulse profile.\n");
    printf("  -nfft       Set size of fft's [%d].\n", fft_size);
    printf("  -maciej     Enable test mode for Maciej.\n");
    printf("  -slope      A simple linearly increasing signal is superimposed on Stokes I.\n"); 
    printf("  -marknoise I  Inject a recognizable signal of amplitude I which is different\n"); 
    printf("              for different subintegrations. I should be positive.\n");
    printf("  -birdie     \"number intensity\"  Set the number of frequency channels (chosen\n");
    printf("              at random) affected by RFI. RFI is simulated as extra white-noise \n");
    printf("              with intensity I.\n");
    return 0;
  }else {
    for(i = 1; i < argc; i++) {
      index = i;
      if(processCommandLine(&application, argc, argv, &index)) {
	i = index;
      }else if(strcmp(argv[i], "-np") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &(outputOpt.NrPulses), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-nf") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &(outputOpt.NrChannels), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-nb") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%ld", &(outputOpt.NrBins), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-a") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &(pulsar.alpha), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-b") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &(pulsar.beta), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-phase0") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &phase0, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-s2n") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &s2n_pp, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-s3n") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &s2n, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-freq") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &outputOpt.frequency, &outputOpt.bw, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-nfft") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &fft_size, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-mapnrB") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d", &nrfieldlines, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-mapmarkphase") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &markphase, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-P") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &(pulsar.P1), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-rm") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &(pulsar.rm), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
	dorvm = 1;
      }else if(strcmp(argv[i], "-l") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &(outputOpt.l1), &(outputOpt.l2), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-los1") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %f", &(IP_flag), &(los_x1), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
	los_flag = 1;
      }else if(strcmp(argv[i], "-los2") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %f %f %f", &(IP_flag), &(los_x1), &(los_x2), &(los_height), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
	los_flag = 2;
      }else if(strcmp(argv[i], "-scatter") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f", &scatter_timescale, &scatter_reffreq, &scatter_pwrlaw, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-cone") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f %d %f %f %f %f %f", &(pulsar.cone[pulsar.nrcones].I), &(pulsar.cone[pulsar.nrcones].rho), &(pulsar.cone[pulsar.nrcones].halfwidth), &(pulsar.cone[pulsar.nrcones].nrsparks), &(pulsar.cone[pulsar.nrcones].P4), &(pulsar.cone[pulsar.nrcones].phi0), &(pulsar.cone[pulsar.nrcones].fP4), &(pulsar.cone[pulsar.nrcones].t4), &(pulsar.cone[pulsar.nrcones].phi4), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	pulsar.cone[pulsar.nrcones].argc = 0;
	pulsar.cone[pulsar.nrcones].map_filename = NULL;
	pulsar.cone[pulsar.nrcones].p3fold_filename = NULL;
	pulsar.cone[pulsar.nrcones].e = pulsar.cone[pulsar.nrcones].or = pulsar.cone[pulsar.nrcones].spacing = 0;
	(pulsar.nrcones)++;
        i++;
      }else if(strcmp(argv[i], "-readmap") == 0) {
	pulsar.cone[pulsar.nrcones].p3fold_filename = NULL;
	pulsar.cone[pulsar.nrcones].map_filename = malloc(1001);
	if(pulsar.cone[pulsar.nrcones].map_filename == NULL) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Memory allocation error.");
	  return 0;
	}
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%1000s %f %f %f %f %f %f", pulsar.cone[pulsar.nrcones].map_filename, &(pulsar.cone[0].I), &(pulsar.cone[pulsar.nrcones].P4), &(pulsar.cone[pulsar.nrcones].phi0), &(pulsar.cone[pulsar.nrcones].fP4), &(pulsar.cone[pulsar.nrcones].t4), &(pulsar.cone[pulsar.nrcones].phi4), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	pulsar.cone[pulsar.nrcones].rho = pulsar.cone[pulsar.nrcones].halfwidth = 0;
	pulsar.cone[pulsar.nrcones].nrsparks = 2;  // Maybe useful later on so the code knows there is a rotating carousel, rather than a non-varying cone.
	pulsar.cone[pulsar.nrcones].e = pulsar.cone[pulsar.nrcones].or = pulsar.cone[pulsar.nrcones].spacing = 0;
	pulsar.cone[pulsar.nrcones].phi = 0;
	pulsar.cone[pulsar.nrcones].argc = 0;
	(pulsar.nrcones)++;
	nonormalise = 1;  // Do not normalise the pulse stack, which is important because separate polarization channels are (at the moment) dealt with in separate calls to fakePulsar
        i++;
      }else if(strcmp(argv[i], "-readp3fold") == 0) {
	pulsar.cone[pulsar.nrcones].map_filename = NULL;
	pulsar.cone[pulsar.nrcones].p3fold_filename = malloc(1001);
	if(pulsar.cone[pulsar.nrcones].p3fold_filename == NULL) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Memory allocation error.");
	  return 0;
	}
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%1000s %f %f", pulsar.cone[pulsar.nrcones].p3fold_filename, &(pulsar.cone[0].I), &(pulsar.cone[pulsar.nrcones].phi0), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	pulsar.cone[pulsar.nrcones].rho = pulsar.cone[pulsar.nrcones].halfwidth = 0;
	pulsar.cone[pulsar.nrcones].nrsparks = 2;  // Maybe useful later on so the code knows there is a rotating carousel, rather than a non-varying cone.
	pulsar.cone[pulsar.nrcones].fP4 = 0;
	pulsar.cone[pulsar.nrcones].t4 = 0;
	pulsar.cone[pulsar.nrcones].phi4 = 0;
	pulsar.cone[pulsar.nrcones].e = pulsar.cone[pulsar.nrcones].or = pulsar.cone[pulsar.nrcones].spacing = 0;
	pulsar.cone[pulsar.nrcones].phi = 0;
	pulsar.cone[pulsar.nrcones].argc = 0;
	(pulsar.nrcones)++;
	nonormalise = 1;  // Do not normalise the pulse stack, which is important because separate polarization channels are (at the moment) dealt with in separate calls to fakePulsar
        i++;
      }else if(strcmp(argv[i], "-ellipse") == 0) {
	if(pulsar.nrcones == 0) {
	  printerror(application.verbose_state.debug, "Use the '%s' AFTER the -cone option.", argv[i]);
	  return 0;
	}
	if(pulsar.cone[pulsar.nrcones-1].map_filename != NULL || pulsar.cone[pulsar.nrcones-1].p3fold_filename != NULL) {
	  printerror(application.verbose_state.debug, "The '%s' option cannot be used after the -readmap of -readp3fold option.", argv[i]);
	  return 0;
	}
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %d", &(pulsar.cone[pulsar.nrcones-1].e), &(pulsar.cone[pulsar.nrcones-1].or), &(pulsar.cone[pulsar.nrcones-1].spacing), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-coneI") == 0) {
	if(pulsar.nrcones < 1) {
	  printerror(application.verbose_state.debug, "fakePulsar: The -coneI option must be specified AFTER the -cone option.\n");
	  return 0;
	}
	pulsar.cone[pulsar.nrcones-1].argc = i+1;
	pulsar.cone[pulsar.nrcones-1].argv = argv;
	i++;
      }else if(strcmp(argv[i], "-wave") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f %f %f %f %f", &(pulsar.wave[pulsar.nrwaves].I), &(pulsar.wave[pulsar.nrwaves].lambda), &(pulsar.wave[pulsar.nrwaves].P4), &(pulsar.wave[pulsar.nrwaves].phi0), &(pulsar.wave[pulsar.nrwaves].fP4), &(pulsar.wave[pulsar.nrwaves].t4), &(pulsar.wave[pulsar.nrwaves].phi4), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	pulsar.wave[pulsar.nrwaves].phi = 0;
	(pulsar.nrwaves)++;
        i++;
      }else if(strcmp(argv[i], "-oscillation") == 0 || strcmp(argv[i], "-oscillation2") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %d %f %f %f %f %f %f", &(pulsar.oscillation[pulsar.nroscillations].l), &(pulsar.oscillation[pulsar.nroscillations].m), &(pulsar.oscillation[pulsar.nroscillations].I), &(pulsar.oscillation[pulsar.nroscillations].P4), &(pulsar.oscillation[pulsar.nroscillations].phi0), &(pulsar.oscillation[pulsar.nroscillations].fP4), &(pulsar.oscillation[pulsar.nroscillations].t4), &(pulsar.oscillation[pulsar.nroscillations].phi4), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	if(strcmp(argv[i], "-oscillation") == 0)
	  pulsar.oscillation[pulsar.nroscillations].alwayspos = 0;
	else
	  pulsar.oscillation[pulsar.nroscillations].alwayspos = 1;
	(pulsar.nroscillations)++;
        i++;
      }else if(strcmp(argv[i], "-patch") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f %f", &(pulsar.patch[pulsar.nrpatches].I), &(pulsar.patch[pulsar.nrpatches].x), &(pulsar.patch[pulsar.nrpatches].y), &(pulsar.patch[pulsar.nrpatches].fwhm), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	(pulsar.nrpatches)++;
        i++;
      }else if(strcmp(argv[i], "-tophat") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &(pulsar.tophat[pulsar.nrtophats].I), &(pulsar.tophat[pulsar.nrtophats].fwhm), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	(pulsar.nrtophats)++;
        i++;
      }else if(strcmp(argv[i], "-micro") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f", &(pulsar.microstructure[pulsar.nrmicrostructure].I), &(pulsar.microstructure[pulsar.nrmicrostructure].separation), &(pulsar.microstructure[pulsar.nrmicrostructure].fwhm), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	(pulsar.nrmicrostructure)++;
        i++;
      }else if(strcmp(argv[i], "-2dsinusoid") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &(pulsar.sinusoid[pulsar.nrsinusoids].p3), &(pulsar.sinusoid[pulsar.nrsinusoids].p2), NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	(pulsar.nrsinusoids)++;
        i++;
      }else if(strcmp(argv[i], "-ppgplot") == 0) {
	ppgplot_name = argv[i+1];
        i++;
      }else if(strcmp(argv[i], "-mapd") == 0) {
	strcpy(mapdevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-profd") == 0) {
	strcpy(profiledevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-stackd") == 0) {
	strcpy(stackdevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-lrfsd") == 0) {
	strcpy(lrfsdevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-trackd") == 0) {
	strcpy(trackdevice, argv[i+1]);
	i++;
      }else if(strcmp(argv[i], "-gg") == 0) {
	write_flag = i+1;
	i++;
      }else if(strcmp(argv[i], "-map") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %f %f", &projection, &projection_long_rot, &projection_lat_rot, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	projection_long_rot *= M_PI/180.0;
	projection_lat_rot *= M_PI/180.0;
	showmap_flag = 1;
        i++;
      }else if(strcmp(argv[i], "-prof") == 0) {
	showprofile_flag = 1;
      }else if(strcmp(argv[i], "-stack") == 0) {
	showpulsstack_flag = 1;
      }else if(strcmp(argv[i], "-lrfs") == 0) {
	lrfs_flag = 1;
      }else if(strcmp(argv[i], "-clip0") == 0) {
	clip_lrfs_flag = 1;
      }else if(strcmp(argv[i], "-rvm") == 0) {
	dorvm = 1;
      }else if(strcmp(argv[i], "-paswing") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f %f %f", &overwritealpha, &overwritebeta, &overwritepa0, &overwritel0, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	overwritePA = 1;
	dorvm = 1;
        i++;
      }else if(strcasecmp(argv[i], "-bcw") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &emissionheightBCW, &emissionheightPWRLAW, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
        i++;
      }else if(strcmp(argv[i], "-opm") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f %f", &jump_longitude[nrJumps], &jump_offset[nrJumps], NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	// Pulse is centred at 180 deg, which is calculated later as longitude 0
	jump_longitude[nrJumps] -= 180;
	nrJumps++;
	i++;
	dorvm = 1;
      }else if(strcmp(argv[i], "-birdie") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%d %f", &birdie_number, &birdie_I, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	dobirdie = 1;
	i++;
      }else if(strcmp(argv[i], "-marknoise") == 0) {
	if(parse_command_string(application.verbose_state, argc, argv, i+1, 0, -1, "%f", &mark_noise_pulses_amplitude, NULL) == 0) {
	  printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot parse '%s' option.", argv[i]);
	  return 0;
	}
	i++;
      }else if(strcmp(argv[i], "-maciej") == 0) {
	maciej_flag = 1;
      }else if(strcmp(argv[i], "-slope") == 0) {
	slope_flag = 1;
      }else {
	printerror(application.verbose_state.debug, "fakePulsar: Unknown option: %s\n\nRun fakePulsar without command line arguments to show help", argv[i]);
	return 0;
      }
    }
  }

  if(los_flag) {
    calculateLOS(los_flag, IP_flag, pulsar.alpha, pulsar.beta, los_x1, los_x2, los_height);
    return 0;
  }

  if(pulsar.nrpatches == 0 && pulsar.nrtophats == 0 && pulsar.nrcones == 0 && pulsar.nrwaves == 0 && pulsar.nroscillations == 0 && pulsar.nrmicrostructure == 0 && pulsar.nrsinusoids == 0 && slope_flag == 0 && mark_noise_pulses_amplitude == 0) {
    printerror(application.verbose_state.debug, "Nothing to do.");
    return 0;
  }


  for(i = 0; i < pulsar.nrcones; i++) {
    if(pulsar.cone[i].map_filename != NULL) {
      if(openPSRData(&(pulsar.cone[i].map), pulsar.cone[i].map_filename, 0, 0, 1, 0, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR fakePulsar: Opening of %s failed.", pulsar.cone[i].map_filename);
	return 0;
      }
      if(pulsar.cone[i].map.gentype != GENTYPE_POLARMAP) {
	printerror(application.verbose_state.debug, "ERROR fakePulsar: File %s does not appear to be of the correct type. This data is expected to be in the format written out by pfold with the -carousel option.", pulsar.cone[i].map_filename);
	return 0;
      }
      if(pulsar.cone[i].map.xrangeset == 0 || pulsar.cone[i].map.yrangeset == 0) {
	printerror(application.verbose_state.debug, "ERROR fakePulsar: File %s does not appear to have the dimensions of the beam map stored. This data is expected to be in the format written out by pfold with the -carousel option.", pulsar.cone[i].map_filename);
	return 0;
      }
      if(pulsar.cone[i].map.NrPols != 1) {
	printwarning(application.verbose_state.debug, "WARNING fakePulsar: File %s has multiple polarization channels defined. Only the first channel will be used.", pulsar.cone[i].map_filename);
      }
      if(pulsar.cone[i].map.NrFreqChan != 1) {
	printwarning(application.verbose_state.debug, "WARNING fakePulsar: File %s has multiple frequency channels defined. Only the first channel will be used.", pulsar.cone[i].map_filename);
      }
    }else if(pulsar.cone[i].p3fold_filename != NULL) {
      if(openPSRData(&(pulsar.cone[i].map), pulsar.cone[i].p3fold_filename, 0, 0, 1, 0, application.verbose_state) == 0) {
	printerror(application.verbose_state.debug, "ERROR fakePulsar: Opening of %s failed.", pulsar.cone[i].p3fold_filename);
	return 0;
      }
      if(pulsar.cone[i].map.gentype != GENTYPE_P3FOLD) {
	printerror(application.verbose_state.debug, "ERROR fakePulsar: File %s does not appear to be of the correct type. This data is expected to be in the format written out by pfold with the -p3fold option.", pulsar.cone[i].p3fold_filename);
	return 0;
      }
      if(pulsar.cone[i].map.yrangeset == 0) {
	printerror(application.verbose_state.debug, "ERROR fakePulsar: File %s does not appear to have the dimensions of the p3-fold (i.e. P3) stored. This data is expected to be in the format written out by pfold with the -p3fold option.", pulsar.cone[i].p3fold_filename);
	return 0;
      }
      if(pulsar.cone[i].map.NrPols != 1) {
	printwarning(application.verbose_state.debug, "WARNING fakePulsar: File %s has multiple polarization channels defined. Only the first channel will be used.", pulsar.cone[i].p3fold_filename);
      }
      if(pulsar.cone[i].map.NrFreqChan != 1) {
	printwarning(application.verbose_state.debug, "WARNING fakePulsar: File %s has multiple frequency channels defined. Only the first channel will be used.", pulsar.cone[i].p3fold_filename);
      }
      if(pulsar.cone[i].map.NrBins != outputOpt.NrBins) {
	printwarning(application.verbose_state.debug, "WARNING fakePulsar: File %s has %ld pulse longitude bins. Can only simulate a pulse-stack with the same number of bins (specify using the -nb option).", pulsar.cone[i].p3fold_filename);
      }
    }
  }

  if(application.verbose_state.verbose) {
      printf("Geometry:\n");
      printf("  P1        = %f sec\n", pulsar.P1);
      printf("  alpha     = %f deg\n", pulsar.alpha);
      printf("  beta      = %f deg\n", pulsar.beta);
    for(i = 0; i < pulsar.nrpatches; i++) {
      printf("Patch %ld:\n", i+1);
      printf("  I         = %f\n", pulsar.patch[i].I);
      printf("  x         = %f deg\n", pulsar.patch[i].x);
      printf("  y         = %f deg\n", pulsar.patch[i].y);
      printf("  fwhm      = %f deg\n", pulsar.patch[i].fwhm);
    }
    for(i = 0; i < pulsar.nrtophats; i++) {
      printf("Top-hat %ld:\n", i+1);
      printf("  I         = %f\n", pulsar.tophat[i].I);
      printf("  fwhm      = %f deg\n", pulsar.tophat[i].fwhm);
    }
    for(i = 0; i < pulsar.nrcones; i++) {
      printf("Cone %ld:\n", i+1);
      printf("  I         = %f\n", pulsar.cone[i].I);
      if(pulsar.cone[i].map_filename != NULL) {
	printf("  map file  = %s\n", pulsar.cone[i].map_filename);
      }else if(pulsar.cone[i].p3fold_filename != NULL) {
	printf("  p3foldfile= %s\n", pulsar.cone[i].p3fold_filename);
      }else {
	printf("  rho       = %f deg\n", pulsar.cone[i].rho);
	printf("  halfwidth = %f deg\n", pulsar.cone[i].halfwidth);
	printf("  nrsparks  = %d\n", pulsar.cone[i].nrsparks);
      }
      if(pulsar.cone[i].p3fold_filename == NULL) {
	printf("  P4        = %f sec\n", pulsar.cone[i].P4);
      } 
      printf("  phi0      = %f deg\n", pulsar.cone[i].phi0);
     if(pulsar.cone[i].p3fold_filename == NULL) {
	printf("  fP4       = %f\n", pulsar.cone[i].fP4);
	printf("  t4        = %f sec\n", pulsar.cone[i].t4);
	printf("  phi4      = %f deg\n", pulsar.cone[i].phi4);
      }
      if(pulsar.cone[i].map_filename == NULL || pulsar.cone[i].p3fold_filename == NULL) {
	printf("  e         = %f\n", pulsar.cone[i].e);
	printf("  or        = %f deg\n", pulsar.cone[i].or);
	printf("  spacing   = %d\n", pulsar.cone[i].spacing);
      }
    }
    for(i = 0; i < pulsar.nrwaves; i++) {
      printf("Wave %ld:\n", i+1);
      printf("  I         = %f\n", pulsar.wave[i].I);
      printf("  lambda    = %f deg\n", pulsar.wave[i].lambda);
      printf("  P4        = %f sec\n", pulsar.wave[i].P4);
      printf("  phi0      = %f deg\n", pulsar.wave[i].phi0);
      printf("  fP4       = %f\n", pulsar.wave[i].fP4);
      printf("  t4        = %f sec\n", pulsar.wave[i].t4);
      printf("  phi4      = %f deg\n", pulsar.wave[i].phi4);
    }
    for(i = 0; i < pulsar.nroscillations; i++) {
      printf("Non radial oscillation %ld:\n", i+1);
      printf("  l         = %d\n", pulsar.oscillation[i].l);
      printf("  m         = %d\n", pulsar.oscillation[i].m);
      printf("  I         = %f\n", pulsar.oscillation[i].I);
      printf("  P4        = %f sec\n", pulsar.oscillation[i].P4);
      printf("  phi0      = %f deg\n", pulsar.oscillation[i].phi0);
      printf("  fP4       = %f\n", pulsar.oscillation[i].fP4);
      printf("  t4        = %f sec\n", pulsar.oscillation[i].t4);
      printf("  phi4      = %f deg\n", pulsar.oscillation[i].phi4);
      if(pulsar.oscillation[i].alwayspos)
	printf("  Allow Neg = NO\n");
      else
	printf("  Allow Neg = YES\n");	
    }
    for(i = 0; i < pulsar.nrmicrostructure; i++) {
      printf("Microstructure %ld:\n", i+1);
      printf("  I         = %f\n", pulsar.microstructure[i].I);
      printf("  sep       = %f deg\n", pulsar.microstructure[i].separation);
      printf("  fwhm      = %f deg\n", pulsar.microstructure[i].fwhm);
    }
    for(i = 0; i < pulsar.nrsinusoids; i++) {
      printf("2D sinusoid %ld:\n", i+1);
      printf("  P3        = %f pulses\n", pulsar.sinusoid[i].p3);
      printf("  P2        = %f phase\n", pulsar.sinusoid[i].p2);
    }
    printf("\n");
  }

  /* Convert angles to radians */
  pulsar.alpha *= M_PI/180.0;
  pulsar.beta *= M_PI/180.0;
  for(i = 0; i < pulsar.nrcones; i++) {
    pulsar.cone[i].rho *= M_PI/180.0;
    pulsar.cone[i].halfwidth *= M_PI/180.0;
    pulsar.cone[i].phi0 *= M_PI/180.0;
    pulsar.cone[i].phi4 *= M_PI/180.0;
    pulsar.cone[i].or *= M_PI/180.0;
    pulsar.cone[i].lookupE = NULL;
  }
  for(i = 0; i < pulsar.nrwaves; i++) {
    pulsar.wave[i].lambda *= M_PI/180.0;
    pulsar.wave[i].phi0 *= M_PI/180.0;
    pulsar.wave[i].phi4 *= M_PI/180.0;
  }
  for(i = 0; i < pulsar.nroscillations; i++) {
    pulsar.oscillation[i].phi0 *= M_PI/180.0;
    pulsar.oscillation[i].phi4 *= M_PI/180.0;
  }
  for(i = 0; i < pulsar.nrpatches; i++) {
    pulsar.patch[i].x *= M_PI/180.0;
    pulsar.patch[i].y *= M_PI/180.0;
    pulsar.patch[i].fwhm *= M_PI/180.0;
  }
  for(i = 0; i < pulsar.nrtophats; i++) {
    pulsar.tophat[i].fwhm *= M_PI/180.0;
  }
  for(i = 0; i < pulsar.nrmicrostructure; i++) {
    pulsar.microstructure[i].separation *= M_PI/180.0;
    pulsar.microstructure[i].fwhm *= M_PI/180.0;
  }
  pulsar.zeta = pulsar.alpha + pulsar.beta;
  if(pulsar.zeta < 0) {
    printerror(application.verbose_state.debug, "Please ensure alpha+beta is positive. Maybe you want to use beta=%f?", (-pulsar.zeta - pulsar.alpha)*180.0/M_PI);
    return 0;
  }
  

  /* Calculate start time to put peak at right spot */
  outputOpt.tstart = -phase0*pulsar.P1/360.0;

  if(s2n_pp > 0) {
    // If S/N per channel and pulse is defined, scale it up to get expected s2n of profile
    s2n = s2n_pp*sqrt(outputOpt.NrPulses*outputOpt.NrChannels);
  }

  if(s2n > 0 && nonormalise) { // Actually, even if not normalised this should work I think, but it certainly isn't tested
    printerror(application.verbose_state.debug, "fakePulsar: Noise addition is not allowed when normalisation of the output is disabled.");
    return 0;
  }
  if(dobirdie && nonormalise) { // Actually, even if not normalised this should work I think, but it certainly isn't tested
    printerror(application.verbose_state.debug, "fakePulsar: The -birdie option is not allowed when normalisation of the output is disabled.");
    return 0;
  }

  if(ppgplot_name != NULL) {
    if(pgopenoutputfile(ppgplot_name) == 0) {
      printerror(application.verbose_state.debug, "fakePulsar: Cannot open %s", ppgplot_name);
      return 0;
    }
  }

  if(showmap_flag)
    showmap(&pulsar, 600, 600, nrfieldlines, markphase, mapdevice, 0, projection, projection_long_rot, projection_lat_rot);

  cleanPSRData(&fstack, application.verbose_state);
  fstack.NrPols = 1;
  if(dorvm)
    fstack.NrPols = 4;
  if(application.verbose_state.verbose) {
    printf("Output:\n");
    printf("  Nr pulses                = %ld\n", outputOpt.NrPulses);
    printf("  Nr bins                  = %ld\n", outputOpt.NrBins);
    printf("  Nr frequency channels    = %ld\n", outputOpt.NrChannels);
    printf("  Nr polarization channels = %ld\n", fstack.NrPols);
    if(write_flag) {
      printf("  Output file = %s\n", argv[write_flag]);
    }
  }
  if(write_flag) {
    if(openPSRData(&fstack, argv[write_flag], application.oformat, 1, 0, 0, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "ERROR fakePulsar: Error opening file.");
      return 0;
    }
  }else {
    printwarning(application.verbose_state.debug, "WARNING fakePulsar: Use the -gg option if you want to write out the data to a file.");
  }
  fstack.NrSubints = outputOpt.NrPulses;
  fstack.NrBins = outputOpt.NrBins;
  fstack.NrFreqChan = outputOpt.NrChannels;
  fstack.isFolded = 1;
  fstack.foldMode = FOLDMODE_FIXEDPERIOD;
  fstack.fixedPeriod = pulsar.P1;
  fstack.fixedtsamp = pulsar.P1*(outputOpt.l2-outputOpt.l1)/(360.0*fstack.NrBins);
  fstack.tsampMode = TSAMPMODE_FIXEDTSAMP;
  fstack.freqMode = FREQMODE_UNIFORM;
  if(fstack.freqlabel_list != NULL) {
    free(fstack.freqlabel_list);
    fstack.freqlabel_list = NULL;
  }
  //  fstack.freq_list = malloc(2*sizeof(double));
  //  if(fstack.freq_list == NULL) {
  //    fflush(stdout);
  //    printerror(application.verbose_state.debug, "ERROR fakePulsar: Memory allocation error.");
  //    return 0;
  //  }
  set_centre_frequency(&fstack, outputOpt.frequency, application.verbose_state);
  fstack.freq_ref = 1e10; // Infinite frequency
  if(set_bandwidth(&fstack, outputOpt.bw, application.verbose_state) == 0) {
    printerror(application.verbose_state.debug, "ERROR fakePulsar: Bandwidth changing failed.");
    return 0;
  }
  //    fstack.channelbw = outputOpt.bw/(float)fstack.NrFreqChan;
  fstack.isDeDisp = 1;
  fstack.isDeFarad = 0;
  fstack.isDePar = 1;
  fstack.isDebase = 1;
  if(slope_flag)
    fstack.isDebase = 0;
  fstack.poltype = POLTYPE_STOKES;
  fstack.feedtype = FEEDTYPE_LINEAR;
  fstack.mjd_start = 43843.5;
  fstack.gentype = GENTYPE_PULSESTACK;
  fstack.tsubMode = TSUBMODE_FIXEDTSUB;
  free(fstack.tsub_list); // Allocated something in cleanPSRData
  fstack.tsub_list = (double *)malloc(sizeof(double));
  if(fstack.tsub_list == NULL) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR fakePulsar: Memory allocation error");
    return 0;
  }
  fstack.tsub_list[0] = fstack.fixedPeriod;
  //    printf("XXXX %lf\n", fstack.tobs);
  if(set_psrname_PSRData(&fstack, "PSR ATNF-4651", application.verbose_state) == 0) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR fakePulsar: Setting pulsar name failed.");
    return 0;
  }
  if(set_observatory_PSRData(&fstack, "North Pole", application.verbose_state) == 0) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR fakePulsar: Setting observatory name failed.");
    return 0;
  }
  if(set_scanID_PSRData(&fstack, "fakePulsar", application.verbose_state) == 0) {
    fflush(stdout);
    printerror(application.verbose_state.debug, "ERROR fakePulsar: Setting scan ID failed.");
    return 0;
  }
  fstack.telescope_X = 0;
  fstack.telescope_Y = 0;
  fstack.telescope_Z = 6371000;
  if(write_flag) {
    if(writeHeaderPSRData(&fstack, argc, argv, application.history_cmd_only, application.verbose_state) == 0) {
      printerror(application.verbose_state.debug, "Error writing header.");
      return 0;
    }
    //    appendHistoryLine(fstack, argc, argv, application.verbose_state);
  }

  //  fprintf(stderr, "Allocating %ld bytes\n", outputOpt.NrPulses*outputOpt.NrBins*outputOpt.NrChannels*fstack.NrPols*sizeof(float));
  Istack = (float *)malloc(outputOpt.NrPulses*outputOpt.NrBins*outputOpt.NrChannels*fstack.NrPols*sizeof(float));
  Iprofile = (float *)calloc(outputOpt.NrBins*fstack.NrPols, sizeof(float));
  if(Istack == NULL || Iprofile == NULL) {
    printerror(application.verbose_state.debug, "Cannot allocate memory");
    return 0;
  }

  /* Old NR code:
     long idnum_microstructure;
     idnum_microstructure = 0;
  */
  gsl_rng *microstructure_num_gen;
  const gsl_rng_type *rand_num_gen_type;
  gsl_rng_env_setup();    /* Set the default generators, can be influenced by environment variables */
  rand_num_gen_type = gsl_rng_default;
  microstructure_num_gen = gsl_rng_alloc(rand_num_gen_type);
  if(application.fixseed)
    idnum = 1;
  else
    randomize_idnum(&idnum);
  gsl_rng_set(microstructure_num_gen, idnum);


  gsl_rng *rand_num_gen;
  rand_num_gen = gsl_rng_alloc(rand_num_gen_type);
  if(application.fixseed)
    idnum = 2;
  else
    randomize_idnum(&idnum);
  gsl_rng_set(rand_num_gen, idnum);

  if(dobirdie) {  // If required, make a selection of channels to be affected by RFI.
    if(birdie_number > outputOpt.NrChannels) {
      printerror(application.verbose_state.debug, "ERROR fakePulsar: The number of birdies cannot exceed the number of frequency channels generated");
      return 0;
    }
    birdie_list = malloc(outputOpt.NrChannels*sizeof(int));
    if(birdie_list == NULL) {
      printerror(application.verbose_state.debug, "ERROR fakePulsar: Cannot allocate memory");
      return 0;
    }
    for(nf = 0; nf < outputOpt.NrChannels; nf++) {
      birdie_list[nf] = 0;  // Default: not a birdie
    }
    for(nf = 0; nf < birdie_number; nf++) {
      //      int ok;
      long freqnr;
      do {   // Select a random frequency channel which is not yet selected before to ensure there are really birdie_number selected channels
	//	ok = 0;
	freqnr = outputOpt.NrChannels*gsl_rng_uniform(rand_num_gen);
	if(freqnr >= outputOpt.NrChannels)  // Shouldn't happen
	  freqnr = 0;
      }while(birdie_list[freqnr] != 0);
      birdie_list[freqnr] = 1;
    }
  }

  /* Used for nummerical integration of carousel rotation phase */
  if(application.verbose_state.verbose) printf("Calculating pulse stack\n");
  reset_num_integration_p4 = 1;
  Itot = 0;
  for(np = 0; np < outputOpt.NrPulses; np++) {
    for(i = 0; i < pulsar.nrmicrostructure; i++) {
      /* old NR code 
      pulsar.microstructure[i].phi0 = 2*M_PI*ran0(&idnum_microstructure)-M_PI;
      */
      pulsar.microstructure[i].phi0 = 2*M_PI*gsl_rng_uniform(microstructure_num_gen)-M_PI;
    }
    /* old NR code
       I2 = ran0(&idnum_microstructure)-0.5;
    */
    I2 = gsl_rng_uniform(microstructure_num_gen)-0.5;
    for(nb = 0; nb < outputOpt.NrBins; nb++) {
      /*      t = np*pulsar.P1 + nb*pulsar.P1/(float)outputOpt.NrBins + outputOpt.tstart; */
      t = np*pulsar.P1 + ((outputOpt.l1 + nb*(outputOpt.l2-outputOpt.l1)/(float)outputOpt.NrBins)/360.0)*pulsar.P1 + outputOpt.tstart;
      /* Calculates the emission of the polar cap, or returns 1 if no polar cap is defined */
      I = CalculateEmission(t, &pulsar, outputOpt, application.verbose_state);
      I += CalculateP3foldEmission(&pulsar, nb, np, application.verbose_state);
      if(maciej_flag == 1) {
	/*      Random spikes throughout whole pulse stacks 
	for(i = 0; i < 80; i++) {
      	  if(nb == (int)((outputOpt.NrBins)*ran0(&idnum2)))
	    I += 300.0*fabs(gasdev(&idnum));
	} 
	*/
	/*  Spikes at given longitude with random amplitude 
	I2 = 30.0*(0.5+0.5*fabs(gasdev(&idnum)));
	if(nb == 150 || nb == 151 || nb == 100 || nb == 101)
	I += I2; */
	/* Baseline fluctuations */
	I += 1000*I2;
      }
      if(pulsar.nrsinusoids != 0) {
	for(i = 0; i < pulsar.nrsinusoids; i++) {
	  I *= 0.5*(1+sin(2*M_PI*np/pulsar.sinusoid[i].p3 + 2*M_PI*((outputOpt.l1 + nb*(outputOpt.l2-outputOpt.l1)/(float)outputOpt.NrBins)/360.0)/pulsar.sinusoid[i].p2));
	}
      }
      if(pulsar.nrtophats != 0) {
	double f;
	f = 2.0*M_PI*((float)nb+0.5)/(float)outputOpt.NrBins;
	for(i = 0; i < pulsar.nrtophats; i++) {
	  //	  printf("XXXX %f %f\n", f, pulsar.tophat[i].fwhm );
	  if(f >= M_PI-0.5*pulsar.tophat[i].fwhm && f <= M_PI+0.5*pulsar.tophat[i].fwhm)
	    I *= pulsar.tophat[i].I;
	  else
	    I = 0;
	}
      }
      Istack[np*outputOpt.NrBins + nb] = I;

      // Generate individual frequency channels
      nidx = fstack.NrPols*outputOpt.NrPulses*outputOpt.NrBins;
      phi0 = 2.0*calcRMAngle(get_weighted_channel_freq(fstack, np, 0, application.verbose_state), -1, 1, pulsar.rm);
      for(nf = 0; nf < outputOpt.NrChannels; nf++) {
	double bcw_shift, bcw_height;
	bcw_height = emissionheightBCW*pow(get_nonweighted_channel_freq(fstack, nf, application.verbose_state)/get_centre_frequency(fstack, application.verbose_state), emissionheightPWRLAW);
	if(dorvm) {
	  bcw_shift = 4*bcw_height*180/M_PI;
	  if(overwritePA == 0)
	    pa = paswing(pulsar.alpha*180.0/M_PI, pulsar.beta*180.0/M_PI, t*360.0/pulsar.P1 - bcw_shift, 0.0, 0.0, nrJumps, jump_longitude, jump_offset, 0.0, 0.0);
	  else
	    pa = paswing(overwritealpha, overwritebeta, t*360.0/pulsar.P1 - bcw_shift, overwritepa0, overwritel0, nrJumps, jump_longitude, jump_offset, 0.0, 0.0);
	  pa *= M_PI/180.0;
	}

	// Make a copy of Stokes I
	npol = 0;
	idx = npol*outputOpt.NrPulses*outputOpt.NrBins + np*outputOpt.NrBins + nb;
	Istack[idx + nf*nidx] = Istack[idx];
	Itot += I;

	if(fstack.NrPols == 4) {
	  phi = 2.0*pa;
	  phi += 2.0*calcRMAngle(get_weighted_channel_freq(fstack, np, nf, application.verbose_state), -1, 1, pulsar.rm)-phi0;
	  L = I;    // 100% linear polarization
	  npol = 1;
	  idx = npol*outputOpt.NrPulses*outputOpt.NrBins + np*outputOpt.NrBins + nb;
	  Istack[idx + nf*nidx] = L*cos(phi);
	  npol = 2;
	  idx = npol*outputOpt.NrPulses*outputOpt.NrBins + np*outputOpt.NrBins + nb;
	  Istack[idx + nf*nidx] = L*sin(phi);
	  npol = 3;
	  idx = npol*outputOpt.NrPulses*outputOpt.NrBins + np*outputOpt.NrBins + nb;
	  Istack[idx + nf*nidx] = 0;  // No circular polarization
	}
      }
    }
    // Superimposed baseline signal later so it doesn't affect the S2N calculation and allowing it to be unpolarized.
    for(nb = 0; nb < outputOpt.NrBins; nb++) {
      for(nf = 0; nf < outputOpt.NrChannels; nf++) {
	if(slope_flag) {
	  Istack[np*outputOpt.NrBins + nb + nf*nidx] += (nb + np*outputOpt.NrBins)/(float)outputOpt.NrBins;
	}
	if(mark_noise_pulses_amplitude > 0) {
	  long type = np % 6;   // Number of different types are defined
	  if(type == 0) {   // Slope over pulse
	    Istack[np*outputOpt.NrBins + nb + nf*nidx] += mark_noise_pulses_amplitude*(nb - 0.5*(outputOpt.NrBins-1))/(float)(0.5*outputOpt.NrBins);
	  }else if(type == 1) {  // Positive impulses
	    if(nb % 10 == 0) {
	      Istack[np*outputOpt.NrBins + nb + nf*nidx] += mark_noise_pulses_amplitude;
	    }
	  }else if(type == 2) {  // top-hat deviations
	    if(nb % 20 <= 10) {
	      Istack[np*outputOpt.NrBins + nb + nf*nidx] += mark_noise_pulses_amplitude;
	    }
	  }else if(type == 3) {   // sinusoid
	    Istack[np*outputOpt.NrBins + nb + nf*nidx] += mark_noise_pulses_amplitude*sin(2.0*M_PI*nb/(float)(outputOpt.NrBins-1));
	  }else if(type == 4) {  // sawtooth (shallow rise, sharp drop)
	    if(nb % 20 <= 10) {
	      Istack[np*outputOpt.NrBins + nb + nf*nidx] += mark_noise_pulses_amplitude*(nb % 20)/9.0;
	    }
	  }else if(type == 5) {   // Negative impulses
	    if(nb % 10 == 0) {
	      Istack[np*outputOpt.NrBins + nb + nf*nidx] += -mark_noise_pulses_amplitude;
	    }
	  }
	}
      }
    }
    if(application.verbose_state.verbose && application.verbose_state.nocounters == 0)
      printf("Pulse %ld of the %ld     \r", np+1, outputOpt.NrPulses);
  }
  if(application.verbose_state.verbose) printf("  Done                       \n");
  
  if(scatter_timescale > 0) {
    if(application.verbose_state.verbose) printf("Adding scatter-tail to data\n");
    nidx = fstack.NrPols*outputOpt.NrPulses*outputOpt.NrBins;
    double efoldtime;
    for(nf = 0; nf < outputOpt.NrChannels; nf++) {
      for(np = 0; np < outputOpt.NrPulses; np++) {
	efoldtime = scatter_efoldtime_bins(fstack, np, nf, scatter_timescale, scatter_reffreq, scatter_pwrlaw, application.verbose_state);
	for(npol = 0; npol < fstack.NrPols; npol++) {
	  if(convolveScatterTail_singlepulse(&(Istack[nidx*nf+npol*outputOpt.NrPulses*outputOpt.NrBins + np*outputOpt.NrBins]), fstack.NrBins, efoldtime, application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "ERROR fakePulsar: convolve with scatter tail failed");
	    return 0;
	  }
	}
	if(application.verbose_state.verbose && application.verbose_state.nocounters == 0)
	  printf("  %.1f%%     \r", 100.0*(np+nf*outputOpt.NrPulses)/(outputOpt.NrPulses*outputOpt.NrChannels));
      }
    }
    if(application.verbose_state.verbose) printf("  Done                       \n");
  }

  if(nonormalise == 0) {
    if(application.verbose_state.verbose) printf("Normalising signal\n");
    normfac = 1.0/Itot;
    nidx = fstack.NrPols*outputOpt.NrPulses*outputOpt.NrBins;
    for(nf = 0; nf < outputOpt.NrChannels; nf++) {
      for(np = 0; np < outputOpt.NrPulses; np++) {
	for(nb = 0; nb < outputOpt.NrBins; nb++) {
	  for(npol = 0; npol < fstack.NrPols; npol++) {
	    Istack[nidx*nf+npol*outputOpt.NrPulses*outputOpt.NrBins + np*outputOpt.NrBins + nb] *= normfac;
	  }
	}
	if(application.verbose_state.verbose && application.verbose_state.nocounters == 0)
	  printf("  %.1f%%     \r", 100.0*(np+nf*outputOpt.NrPulses)/(outputOpt.NrPulses*outputOpt.NrChannels));
      }
    }
    if(application.verbose_state.verbose) printf("  Done                       \n");
  }

  if(showprofile_flag || s2n > 0) {
    if(application.verbose_state.verbose) printf("Generating profile\n");
    nidx = fstack.NrPols*outputOpt.NrPulses*outputOpt.NrBins;
    for(nf = 0; nf < outputOpt.NrChannels; nf++) {
      for(np = 0; np < outputOpt.NrPulses; np++) {
	for(nb = 0; nb < outputOpt.NrBins; nb++) {
	  for(npol = 0; npol < fstack.NrPols; npol++) {
	    Iprofile[nb+npol*outputOpt.NrBins] += Istack[nidx*nf+npol*outputOpt.NrPulses*outputOpt.NrBins + np*outputOpt.NrBins + nb];
	  }
	}
 	if(application.verbose_state.verbose && application.verbose_state.nocounters == 0)
	  printf("  %.1f%%     \r", 100.0*(np+nf*outputOpt.NrPulses)/(outputOpt.NrPulses*outputOpt.NrChannels));
     }
    }
    if(application.verbose_state.verbose) printf("  Done                       \n");

    if(s2n > 0) { 
      double f;
      //Expected noise-level needed to get required s2n of profile if integrating over full longitude range
      f = 1.0/(s2n*sqrt(outputOpt.NrBins));
      nronpulsebins = 0;
      // Count nr of bins in profile which has 3 sigma detection
      for(nb = 0; nb < outputOpt.NrBins; nb++) {
	if(Iprofile[nb] > 3*f)
	  nronpulsebins++;
      }
      if(application.verbose_state.verbose)
	printf("Estimated nr of on-pulse bins to determine signal to noise: %d\n", nronpulsebins);
      if(nronpulsebins == 0) {
	printwarning(application.verbose_state.debug, "WARNING fakePulsar: There is no signal sticking out of the specified noise level. The onpulse region is set to be 1 bin wide.");
	nronpulsebins = 1;
      }
      f = 1.0/(s2n*sqrt(nronpulsebins));
      nronpulsebins = 0;
      if(s2n > 0) {
	for(nb = 0; nb < outputOpt.NrBins; nb++) {
	  if(Iprofile[nb] > 3*f)
	    nronpulsebins++;
	}
      }
      if(nronpulsebins == 0) {
	nronpulsebins = 1;
      }
      if(application.verbose_state.verbose)
	printf("Estimated nr of on-pulse bins to determine signal to noise: %d (refined)\n", nronpulsebins);

      if(nronpulsebins == 1) { // This could happen if a very low signal-to-noise is requested, especially if there are a lot of pulse longitude bins
	// Base the rms on W10 rather than that part of the profile sticking out of the noise
	double max_intensity;
	max_intensity = Iprofile[0];
	for(nb = 0; nb < outputOpt.NrBins; nb++) {
	  if(Iprofile[nb] > max_intensity)
	    max_intensity = Iprofile[nb];
	}
	nronpulsebins = 0;
	for(nb = 0; nb < outputOpt.NrBins; nb++) {
	  if(Iprofile[nb] >= 0.1*max_intensity)
	    nronpulsebins++;
	}
	if(application.verbose_state.verbose)
	  printf("Estimated nr of on-pulse bins to determine signal to noise (based on W10): %d\n", nronpulsebins);
      }

    }
  }

  if(s2n > 0 || dobirdie) {
    double f;
    if(s2n > 0) {
      f = 1.0/(s2n*sqrt(nronpulsebins*outputOpt.NrPulses*outputOpt.NrChannels));
      if(application.verbose_state.verbose) printf("Adding noise to pulse stack with rms=%lf based on %d on-pulse bins\n", f, nronpulsebins);
    }else {
      f = 0;
    }
    //    f = s2n*sqrt(outputOpt.NrBins*outputOpt.NrPulses)/Itot;
    // Make total signal in the simulated data = 1
    nidx = fstack.NrPols*outputOpt.NrPulses*outputOpt.NrBins;
    for(nf = 0; nf < outputOpt.NrChannels; nf++) {
      for(np = 0; np < outputOpt.NrPulses; np++) {
	for(nb = 0; nb < outputOpt.NrBins; nb++) {
	  for(npol = 0; npol < fstack.NrPols; npol++) {
	    if(s2n > 0) {
	      Istack[nidx*nf+npol*outputOpt.NrPulses*outputOpt.NrBins + np*outputOpt.NrBins + nb] += f*gsl_ran_gaussian(rand_num_gen, 1.0);
	    }
	    if(dobirdie) {
	      if(birdie_list[nf]) {
		Istack[nidx*nf+npol*outputOpt.NrPulses*outputOpt.NrBins + np*outputOpt.NrBins + nb] += birdie_I*gsl_ran_gaussian(rand_num_gen, 1.0);
	      }
	    }
	  }
	}
	if(application.verbose_state.verbose && application.verbose_state.nocounters == 0)
	  printf("  %.1f%%     \r", 100.0*(np+nf*outputOpt.NrPulses)/(outputOpt.NrPulses*outputOpt.NrChannels));
      }
    }
    if(application.verbose_state.verbose) printf("  Done                       \n");
  }

  if(showprofile_flag && s2n > 0) {
    if(application.verbose_state.verbose) printf("Generating profile after noise is added to data noise\n");
    for(nb = 0; nb < outputOpt.NrBins*fstack.NrPols; nb++) {
      Iprofile[nb] = 0;
    }
    nidx = fstack.NrPols*outputOpt.NrPulses*outputOpt.NrBins;
    for(nf = 0; nf < outputOpt.NrChannels; nf++) {
      for(np = 0; np < outputOpt.NrPulses; np++) {
	for(nb = 0; nb < outputOpt.NrBins; nb++) {
	  for(npol = 0; npol < fstack.NrPols; npol++) {
	    Iprofile[nb+npol*outputOpt.NrBins] += Istack[nidx*nf+npol*outputOpt.NrPulses*outputOpt.NrBins + np*outputOpt.NrBins + nb];
	  }
	}
 	if(application.verbose_state.verbose && application.verbose_state.nocounters == 0)
	  printf("  %.1f%%     \r", 100.0*(np+nf*outputOpt.NrPulses)/(outputOpt.NrPulses*outputOpt.NrChannels));
      }
    }
    if(application.verbose_state.verbose) printf("  Done                       \n");
  }

  if(showpulsstack_flag) {
    strcpy(pgplot_options.viewport.plotDevice, stackdevice);
    strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
    strcpy(pgplot_options.box.ylabel, "Pulse number");
    strcpy(pgplot_options.box.title, "Pulse-stack");
    /*    pgplotMap(&viewport, Istack, outputOpt.NrBins, outputOpt.NrPulses, outputOpt.l1, outputOpt.l2, outputOpt.l1, outputOpt.l2, 0, outputOpt.NrPulses, 0-0.5, outputOpt.NrPulses-0.5, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, "Pulse phase [degrees]", "Pulse number", "Pulse-stack", 1, 1, 1, 0, 0, 1.0, 1.0, 1, 0, 0, 0, 0, 0, 0, 0, application.debug);*/
    pgplotMap(&pgplot_options, Istack, outputOpt.NrBins, outputOpt.NrPulses, outputOpt.l1, outputOpt.l2, outputOpt.l1, outputOpt.l2, 0, outputOpt.NrPulses, 0-0.5, outputOpt.NrPulses-0.5, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, application.verbose_state);
  }
  if(showprofile_flag) {
    strcpy(pgplot_options.viewport.plotDevice, profiledevice);
    strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
    strcpy(pgplot_options.box.ylabel, "Intensity");
    strcpy(pgplot_options.box.title, "Pulse profile");
    pgplotGraph1(&pgplot_options, Iprofile, NULL, NULL, outputOpt.NrBins, outputOpt.l1, outputOpt.l2, 0, outputOpt.l1, outputOpt.l2, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
  }
  
  if(lrfs_flag) {
    lrfs = (float *)malloc((fft_size/2+1)*outputOpt.NrBins*sizeof(float));
    phase_track = (float *)malloc(outputOpt.NrBins*sizeof(float));
    if(lrfs == NULL || phase_track == NULL) {
      printerror(application.verbose_state.debug, "Cannot allocate memory");
      return 0;
    }
    if(calcLRFS(Istack, outputOpt.NrPulses, outputOpt.NrBins, fft_size, lrfs, clip_lrfs_flag, NULL, phase_track, NULL, 1, 0, 0.5, 0, NULL, 0, 0, 0, NULL, NULL, 0, NULL, application.verbose_state)) {
      strcpy(pgplot_options.viewport.plotDevice, lrfsdevice);
      strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.ylabel, "P3 [cpp]");
      strcpy(pgplot_options.box.title, "LRFS");
      /*      pgplotMap(&viewport, lrfs, outputOpt.NrBins, fft_size/2+1, outputOpt.l1, outputOpt.l2, outputOpt.l1, outputOpt.l2, 0, 0.5, 0, 0.5, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, "Pulse phase [degrees]", "P3 [cpp]", "LRFS", 1, 1, 1, 0, 0, 1.0, 1.0, 1, 1, 1, 0, 0, 0, 0, 0, application.debug); */
      pgplotMap(&pgplot_options, lrfs, outputOpt.NrBins, fft_size/2+1, outputOpt.l1, outputOpt.l2, outputOpt.l1, outputOpt.l2, 0, 0.5, 0, 0.5, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, application.verbose_state);
      strcpy(pgplot_options.viewport.plotDevice, trackdevice);
      strcpy(pgplot_options.box.xlabel, "Pulse phase [degrees]");
      strcpy(pgplot_options.box.ylabel, "Subpulse phase");
      strcpy(pgplot_options.box.title, "Subpulse phase track");
      pgplotGraph1(&pgplot_options, phase_track, NULL, NULL, outputOpt.NrBins, outputOpt.l1, outputOpt.l2, 0, outputOpt.l1, outputOpt.l2, 0, 0, 0, 0, 0, 1, 0, 1, 1, NULL, -1, application.verbose_state);
    }
    free(lrfs);
    free(phase_track);
  }

  if(write_flag) {
    if(application.verbose_state.verbose) printf("Writing out pulse stack\n");
    nidx = fstack.NrPols*outputOpt.NrPulses*outputOpt.NrBins;
    for(np = 0; np < outputOpt.NrPulses; np++) {
      for(nf = 0; nf < outputOpt.NrChannels; nf++) {
	for(npol = 0; npol < fstack.NrPols; npol++) {
	  if(writePulsePSRData(&fstack, np, npol, nf, 0, outputOpt.NrBins, &Istack[nf*nidx+npol*outputOpt.NrPulses*outputOpt.NrBins + np*outputOpt.NrBins], application.verbose_state) == 0) {
	    printerror(application.verbose_state.debug, "Error writing data.");
	    return 0;
	  }
	}
	if(application.verbose_state.verbose && application.verbose_state.nocounters == 0)
	  printf("  %.1f%%     \r", 100.0*(nf+np*outputOpt.NrChannels)/(outputOpt.NrPulses*outputOpt.NrChannels));
      }
    }
    closePSRData(&fstack, 0, application.verbose_state);
    if(application.verbose_state.verbose) printf("  Done                       \n");
  }

  gsl_rng_free(microstructure_num_gen);
  gsl_rng_free(rand_num_gen); 
  free(Istack);
  free(Iprofile);
  for(i = 0; i < pulsar.nrcones; i++) {
    if(pulsar.cone[i].map_filename != NULL || pulsar.cone[i].p3fold_filename != NULL) {
      if(closePSRData(&(pulsar.cone[i].map), 0, application.verbose_state) != 0) {
	fflush(stdout);
	printerror(application.verbose_state.debug, "ERROR fakePulsar: Closing file failed");
	return 0;
      }
    }
    if(pulsar.cone[i].map_filename != NULL) {
      free(pulsar.cone[i].map_filename);
    }
    if(pulsar.cone[i].p3fold_filename != NULL) {
      free(pulsar.cone[i].p3fold_filename);
    }
  }
  terminateApplication(&application);
  return 0;
}


/*
void Transformation(float t, PulsarDefinitionStruct pulsar, float *phi, float *theta)
{
  float y = sin(2*M_PI*t/pulsar.P1);
  float x = cos(2*M_PI*t/pulsar.P1)*sin(0.5*M_PI - pulsar.alpha) - 
            tan(0.5*M_PI - pulsar.zeta)*cos(0.5*M_PI - pulsar.alpha);
  *phi = polar_angle_rad(x, y);
  *theta = 0.5*M_PI - (asin(sin(0.5*M_PI - pulsar.alpha)*sin(0.5*M_PI - pulsar.zeta) +
                                 cos(0.5*M_PI - pulsar.alpha)*cos(0.5*M_PI -
				 pulsar.zeta)*cos(2*M_PI*t/pulsar.P1)));
}
*/


void Transformation(float t, PulsarDefinitionStruct pulsar, float *phi, float *theta)
{
  pulsar_polar_coordinates(t, pulsar.P1, pulsar.alpha, pulsar.zeta, phi, theta);
}

/* Return value = 1 (MP) or -1 (IP) */
int Transformation2(float t, PulsarDefinitionStruct pulsar, float *x, float *y)
{
  float phi, theta;
  Transformation(t, pulsar, &phi, &theta);
  *x = theta*sin(phi);
  *y = -theta*cos(phi);
  if(cos(theta) >= 0)
    return 1;
  else
    return -1;
}

/* Radially propagating waves */
float calcWaveEmission(PulsarDefinitionStruct pulsar, float x, float y)
{
  float I, r; //, phi;
  int nw;
  I = 0;
  r = sqrt((x)*(x)+(y)*(y));
  //  phi = polar_angle_rad(x, y);
  for(nw = 0; nw < pulsar.nrwaves; nw++) {  
    I += pulsar.wave[nw].I*0.5*(1+sin(2*M_PI*r/pulsar.wave[nw].lambda + pulsar.wave[nw].phi0 + pulsar.wave[nw].phi)); 
    /*    I += pulsar.wave[nw].I*cos(pulsar.wave[nw].m*phi)*plgndr(pulsar.wave[nw].l,pulsar.wave[nw].m,cos(r+pulsar.wave[nw].phi0 + pulsar.wave[nw].phi)); */
  }
  return I;
}


/* Radially propagating waves */
float calcOscillationEmission(PulsarDefinitionStruct pulsar, float x, float y)
{
  float I, r, phi;
  int no;
  I = 0;
  r = sqrt((x)*(x)+(y)*(y));
  phi = polar_angle_rad(x, y);
  for(no = 0; no < pulsar.nroscillations; no++) {  
    /* Old nr code 
       I += pulsar.oscillation[no].I*cos(pulsar.oscillation[no].m*phi)*plgndr(pulsar.oscillation[no].l,pulsar.oscillation[no].m,cos(r))*cos(pulsar.oscillation[no].phi0 + pulsar.oscillation[no].phi);
    */
    I += pulsar.oscillation[no].I*cos(pulsar.oscillation[no].m*phi)*gsl_sf_legendre_Plm(pulsar.oscillation[no].l, pulsar.oscillation[no].m, cos(r))*cos(pulsar.oscillation[no].phi0 + pulsar.oscillation[no].phi);
    if(pulsar.oscillation[no].alwayspos)
      I += pulsar.oscillation[no].I;
  }
  return I;
}


double ellint_exp_internal;
double *ellint_lookup_table_x = NULL;

double funk_ellint_root_function(double x, void *params)
{
  float *e = (float *)params;
  //  printf("ell = %f (%p)\n", *e, params);
  return ellint_exp_internal - gsl_sf_ellint_E(x, *e, GSL_PREC_DOUBLE);
}

// Given the mean Anomaly (polar coordinate that for a fixed carousel circulation time increases linearly with time), calculate the polar coordinate corresponding to an elliptical path such that the velocity along the path is constant. This means that if the sparks are equally spaced in meanAnomaly (i.e. separated by 2pi/Nsparks), they come out to be equidistant on the ellipse when using the corresponding trueAnomaly. These angles are all measured from the centre of the ellipse and in radians, not the focus as would be done for Keplarian mechanics.
double trueAnomaly(double meanAnomaly, ConeDefinition *cone, int verbose)
{
  int quadrant, ret, i;
  double reducedMA, reducedMA2, quarter, phi, inputangle;
  verbose_definition verbosestate;
  cleanVerboseState(&verbosestate);
  //  verbosestate.verbose = verbose;

  // Make sure that meanAnomaly is a value between 0 and 2pi
  // reducedMA: angle between 0 and pi/2, quadrant counts how the angle was reduced
  reducedMA = derotate_rad_double(meanAnomaly);
  if(reducedMA < 0.5*M_PI) {
    quadrant = 0;
    reducedMA2 = reducedMA;
  }else if(reducedMA < M_PI) {
    quadrant = 1;
    reducedMA2 = M_PI - reducedMA;
  }else if(reducedMA < 1.5*M_PI) {
    quadrant = 2;
    reducedMA2 = reducedMA - M_PI;
  }else {
    quadrant = 3;
    reducedMA2 = 2.0*M_PI - reducedMA;
  }

  if(cone->lookupE == NULL) {
    if(verbose) {
      fprintf(stdout, "Creating elliptical integral lookup table: ");
      fflush(stdout);
    }
    if(ellint_lookup_table_x == NULL) {
      ellint_lookup_table_x = malloc(NrEllIntegralInterpol*sizeof(double));
    }
    cone->lookupE = malloc(NrEllIntegralInterpol*sizeof(double));
    if(cone->lookupE == NULL || ellint_lookup_table_x == NULL) {
      fflush(stdout);
      printerror(verbose, "ERROR trueAnomaly: Cannot allocate memory");
      exit(0);
    }

    // Calculate the complete elliptic integral of the second kind, i.e. circumference of ellipse would be 4*quarter*semimajor axis
    //  printf("Normalised circumference of ellipse is %lf\n", 4.0*quarter);
    //  printf("MA %lf reduced to %lf in quadrant %d\n", reducedMA, reducedMA2, quadrant);
    quarter = gsl_sf_ellint_Ecomp(cone->e, GSL_PREC_DOUBLE);

    for(i = 0; i < NrEllIntegralInterpol; i++) {
      inputangle = 0.5*M_PI*i/(double)(NrEllIntegralInterpol-1);
      ellint_lookup_table_x[i] = inputangle;
      // Condition for equi-spaced points: reducedMA/(pi/2) = E(phi,e)/quarter, where E is the incomplete elliptic integral of the second kind and the complete elliptic integral of the second kind is E(pi/2).
      // Hence we can find the expected solution if the elliptical integral for the correct choice of phi.
      ellint_exp_internal = quarter*inputangle/(0.5*M_PI);
      //  printf("Expected value E = %lf\n", ellint_exp_internal);
      
      //  printf("Actual ell = %lf (%p)\n", cone->e, &(cone->e));
      // Precision is 1e-3
      ret = minimize_1D_double(1, funk_ellint_root_function, &(cone->e), 0, 0.5*M_PI, 0, 0, 0, &phi, 100, 1e-10, 0.0, 0, 0);
      if(ret != 0) {
	if(ret == 1) {
	  fflush(stdout);
	  printerror(verbose, "ERROR: Maximum nr of itterations exceeded in minimize_1D_double");
	  exit(0);
	}else if(ret == 2) {
	  fflush(stdout);
	  printerror(verbose, "ERROR: Did not found root");
	  exit(0);
	}else if(ret == 3) {
	  fflush(stdout);
	  printerror(verbose, "ERROR: Lower and upper limit do not bracket a root in internal_funk_gsl");
	  exit(0);
	}else {
	  fflush(stdout);
	  printerror(verbose, "ERROR: Unknown error in minimize_1D_double in internal_funk_gsl");
	  exit(0);
	}
      }
      //  printf("Found value E = %lf\n", gsl_sf_ellint_E(phi, cone->e, GSL_PREC_DOUBLE));
      //  printf("Difference = %lf\n", funk_ellint_root_function(phi, &(cone->e)));
      cone->lookupE[i] = phi;
    }
    if(verbose) {
      fprintf(stdout, "done\n");
      fflush(stdout);
    }
  }

  if(interpolate_double(ellint_lookup_table_x, cone->lookupE, NrEllIntegralInterpol, reducedMA2, &phi, 2, verbosestate) != 1) {
    fflush(stdout);
    printerror(verbose, "ERROR trueAnomaly: Interpolation failed.");
    exit(0);
  }


  if(quadrant == 1) {
    phi = M_PI - phi;
  }else if(quadrant == 2) {
    phi = phi + M_PI;
  }else if(quadrant == 3) {
    phi = 2.0*M_PI - phi;
  }

  //  printf("Found phi=%lf\n", phi);
  return phi;
}

float CalculateP3foldEmission(PulsarDefinitionStruct *pulsar, long nb, long np, verbose_definition verbose)
{
  long nc;
  float I, Itot;
  Itot = 0;
  for(nc = 0; nc < pulsar->nrcones; nc++) {
    if(pulsar->cone[nc].p3fold_filename != NULL) {
      long p3foldbin_long;
      double p3foldbin;
      double p3;
      p3 = pulsar->cone[nc].map.yrange[1] + pulsar->cone[nc].map.yrange[0]; // This should be P3 in pulses
      p3foldbin = (2.0*M_PI)*np/p3+pulsar->cone[nc].phi0;
      p3foldbin = derotate_rad_double(p3foldbin);
      p3foldbin_long = pulsar->cone[nc].map.NrSubints*p3foldbin/(2*M_PI);
      if(p3foldbin_long < 0) {  // Not sure if possible, but prevent problems with rounding
	p3foldbin_long += pulsar->cone[nc].map.NrSubints;
      }
      if(p3foldbin_long >= pulsar->cone[nc].map.NrSubints) {  // Not sure if possible, but prevent problems with rounding
	p3foldbin_long -= pulsar->cone[nc].map.NrSubints;
      }
      // Could do some interpolation, but didn't do this yet.
      if(readPulsePSRData(&(pulsar->cone[nc].map), p3foldbin_long, 0, 0, nb, 1, &I, verbose) != 1) {
	printerror(application.verbose_state.debug, "fakePulsar: Cannot read data from p3fold map");
	exit(-1);
      }
      Itot += I;
    }
  }
  return Itot;
}

float calcCarouselEmission(PulsarDefinitionStruct *pulsar, float x, float y, verbose_definition verbose)
{
  float I, phi_spark, xspark, yspark, rr, scale, cosor, sinor, tanphi;
  float log2 = log(2);
  int ns, nc, nrwords;
  char *ptr;
  static int firsttime = 1;



  I = 0;
  for(nc = 0; nc < pulsar->nrcones; nc++) {
    if(pulsar->cone[nc].e != 0) {
      sinor = sin(pulsar->cone[nc].or);
      cosor = cos(pulsar->cone[nc].or);
    }
    if(pulsar->cone[nc].p3fold_filename != NULL) {  // Dealth with in a different place where we know bin number etc.
    }else if(pulsar->cone[nc].map_filename != NULL) {
      double phi;
      // Convert Cartesian coordinates of the line of sight to polar coordinates
      rr = sqrt(x*x+y*y);
      phi = atan2(y, x);
      // Rather than rotating the polar map with the carousel period, make the line of sight rotate
      phi += pulsar->cone[nc].phi0 + pulsar->cone[nc].phi;
      // Get the updated Cartesian coordinates
      x = rr*cos(phi)*180.0/M_PI;
      y = rr*sin(phi)*180.0/M_PI;
      int bin_x, bin_y;
      bin_x = (x-pulsar->cone[nc].map.xrange[0])*(pulsar->cone[nc].map.NrBins-1)   /(pulsar->cone[nc].map.xrange[1]-pulsar->cone[nc].map.xrange[0]);
      bin_y = (y-pulsar->cone[nc].map.yrange[0])*(pulsar->cone[nc].map.NrSubints-1)/(pulsar->cone[nc].map.yrange[1]-pulsar->cone[nc].map.yrange[0]);
      // Check if we're inside the defined map
      if(bin_x < 0 || bin_x >= pulsar->cone[nc].map.NrBins || bin_y < 0 || bin_y >= pulsar->cone[nc].map.NrSubints) {
	I = 0;
      }else {
	// Could do some interpolation, but didn't do this yet.
	if(readPulsePSRData(&(pulsar->cone[nc].map), bin_y, 0, 0, bin_x, 1, &I, verbose) != 1) {
	  printerror(application.verbose_state.debug, "fakePulsar: Cannot read data from polar map");
	  exit(-1);
	}
      }
    }else if(pulsar->cone[nc].nrsparks > 0) {
      //      trueAnomaly(M_PI/2.0, pulsar->cone[nc], 0);
      //      exit(0);
      for(ns = 0; ns < pulsar->cone[nc].nrsparks; ns++) {
	phi_spark = ns*2.0*M_PI/(float)pulsar->cone[nc].nrsparks + pulsar->cone[nc].phi0 + pulsar->cone[nc].phi;
	if(pulsar->cone[nc].e == 0) {
	  xspark = pulsar->cone[nc].rho*cos(phi_spark);
	  yspark = pulsar->cone[nc].rho*sin(phi_spark);
	}else {
	  float bb, dummyf;

	  if(pulsar->cone[nc].spacing == 1) {
	    //	    phi_spark -= pulsar->cone[nc].or;
	    phi_spark = trueAnomaly(phi_spark, &(pulsar->cone[nc]), verbose.debug);
	    phi_spark += pulsar->cone[nc].or;
	  }
	  phi_spark += 1e-5;  // Avoid phi_spark to be exactly pi/2, which makes the code extremely slow.
	  tanphi = tan(phi_spark);
	  bb = pulsar->cone[nc].rho*pulsar->cone[nc].rho;
	  dummyf = (cosor+tanphi*sinor);
	  rr = dummyf * dummyf / bb;
	  bb *= (1.0-pulsar->cone[nc].e*pulsar->cone[nc].e);
	  dummyf = (sinor-tanphi*cosor);
	  rr += dummyf * dummyf / bb;
	  //	  if(rr == 0)
	  //	    fprintf(stderr, "Oops, rr is zero!\n");
	  //	  if(rr < 0)
	  //	    fprintf(stderr, "Oops, rr is negative!\n");	    
	  xspark = sqrt(1.0/rr);
	  if(derotate_deg(phi_spark*180.0/M_PI) > 90 && derotate_deg(phi_spark*180.0/M_PI) < 270)
	    xspark *= -1.0;

	  yspark = xspark * tanphi;
	  //	  fprintf(stderr, "XXXXXX %d %d %f %f %f\n", nc, ns, xspark, yspark, rr);
	}
	rr = (xspark-x)*(xspark-x)+(yspark-y)*(yspark-y);
	scale = 1;

	if(pulsar->cone[nc].argc != 0) {
	  ptr = pickWordFromString(pulsar->cone[nc].argv[pulsar->cone[nc].argc], ns+1, &nrwords, 0, ' ', verbose);
	  if(ptr != NULL) {
	    nrwords = sscanf(ptr, "%f", &scale);
	    if(nrwords != 1) {
	      printerror(application.verbose_state.debug, "fakePulsar: Cannot read '%s' as a floating point number", ptr);
	      exit(-1);
	    }
	  }
	  if(application.verbose_state.verbose && firsttime) 
	    printf("Intensity of spark %d: %f\n", ns+1, scale);
	}
	//	I += scale*pulsar->cone[nc].I*2.0*sqrt(log(2))*exp(-4*rr*log(2)/(pulsar->cone[nc].halfwidth*pulsar->cone[nc].halfwidth))/(pulsar->cone[nc].halfwidth*sqrt(M_PI));
	// The spark width was given to be the half-width. Probably normalising by 2D integral makes more sense than the current 1D integral
	I += scale*pulsar->cone[nc].I*0.5*log2*exp(-rr*log2/(pulsar->cone[nc].halfwidth*pulsar->cone[nc].halfwidth))/(pulsar->cone[nc].halfwidth*pulsar->cone[nc].halfwidth*M_PI);
      }
      firsttime = 0;
    }else {
      if(pulsar->cone[nc].e != 0.0) {
	 printerror(application.verbose_state.debug, "The use of zero sparks for an elliptical cone is not implemented.");
	 exit(0);
      }
      phi_spark = polar_angle_rad(x, y);
      xspark = pulsar->cone[nc].rho*cos(phi_spark);
      yspark = pulsar->cone[nc].rho*sin(phi_spark);
      rr = (xspark-x)*(xspark-x)+(yspark-y)*(yspark-y);
      //      I += pulsar->cone[nc].I*2.0*sqrt(log(2))*exp(-4*rr*log(2)/(pulsar->cone[nc].halfwidth*pulsar->cone[nc].halfwidth))/(pulsar->cone[nc].halfwidth*sqrt(M_PI));
      //Keep normalisation the same as for a 2D Gaussian. Maybe should have done this for a ring instead.
      I += pulsar->cone[nc].I*0.5*log2*exp(-rr*log2/(pulsar->cone[nc].halfwidth*pulsar->cone[nc].halfwidth))/(pulsar->cone[nc].halfwidth*pulsar->cone[nc].halfwidth*M_PI);
    }
  }
  return I;
}

float calcMicroStructureEmission(PulsarDefinitionStruct pulsar, float x, float y)
{
  float I, phi_spark, rr, phi;
  int ns, nm, nrsparks;
  I = 0;
  for(nm = 0; nm < pulsar.nrmicrostructure; nm++) {  
    nrsparks = 2*M_PI/pulsar.microstructure[nm].separation;
    for(ns = 0; ns < nrsparks; ns++) {
      phi_spark = ns*pulsar.microstructure[nm].separation + pulsar.microstructure[nm].phi0;
      if(phi_spark > 2.0*M_PI)
	phi_spark -= 2.0*M_PI;
      if(phi_spark < 0)
	phi_spark += 2.0*M_PI;
      phi = polar_angle_rad(x, y);
      rr = (phi-phi_spark)*(phi-phi_spark);
      I += pulsar.microstructure[nm].I*2.0*sqrt(log(2))*exp(-4*rr*log(2)/(pulsar.microstructure[nm].fwhm*pulsar.microstructure[nm].fwhm))/(pulsar.microstructure[nm].fwhm*sqrt(M_PI));
    }
  }
  return I;
}

float calcPatchEmission(PulsarDefinitionStruct pulsar, float x, float y)
{
  float I, rr;
  int np;
  if(pulsar.nrpatches == 0) {
    I = 1;
  }else {
    I = 0;
    for(np = 0; np < pulsar.nrpatches; np++) {  
      rr = (pulsar.patch[np].x-x)*(pulsar.patch[np].x-x)+(pulsar.patch[np].y-y)*(pulsar.patch[np].y-y);
      I += pulsar.patch[np].I*0.5*log(2)*exp(-rr*log(2)/(0.25*pulsar.patch[np].fwhm*pulsar.patch[np].fwhm))/(pulsar.patch[np].fwhm*pulsar.patch[np].fwhm*M_PI); 
      /* Should do a transformation from magnetic coordinates to blob coordinates in order to make sure that the shape of the blobs are really gaussian and not only in projection space. It somehow doesn't work, probably because the transformation does something else than I think it should do. Plotting a phi or theta map produces something different than I would expect.*/
      /*      rr = x*x+y*y;      
	      rr_patch = (pulsar.patch[np].x)*(pulsar.patch[np].x)+(pulsar.patch[np].y)*(pulsar.patch[np].y);
	      phi_patch = polar_angle_rad(pulsar.patch[np].x, pulsar.patch[np].y);
	      pulsar_polar_coordinates(phi_patch, 2.0*M_PI, sqrt(rr_patch), sqrt(rr), &phi, &theta);
	      rr = theta*theta;
	      if(np == 0)
	      printf("x=%f y=%f rr=%f\n", pulsar.patch[np].x, pulsar.patch[np].y, rr_patch); 
	      I += pulsar.patch[np].I*2.0*sqrt(log(2))*exp(-4*rr*log(2)/(0.25*pulsar.patch[np].fwhm*pulsar.patch[np].fwhm))/(0.5*pulsar.patch[np].fwhm*sqrt(M_PI)); 
	      if(np == 0)
	      I += phi;
      */
    }
  }
  return I;
}



/*
int openPulseStack(char *filename, FILE **fstack, OutputOptionsStruct outoptions, PulsarDefinitionStruct pulsar)
{
  Header_type   hdr;

  if(application.verbose_state.verbose) fprintf(stderr, "Going to open %s\n", filename);

  *fstack = fopen(filename, "wb");
  if(*fstack == NULL) {
    printf("Cannot open %s\n\n", filename);
    return 0;
  }
    
  hdr.redn.NBins = outoptions.NrBins;
  hdr.redn.NTimeInts = outoptions.NrPulses;
    
  strcpy(hdr.gen.HdrVer, "DPC_1.3");
  strcpy(hdr.gen.Platform, "PuMa Crate 0");
  strcpy(hdr.gen.ThisFileName, "200899999.00000.0.puma");
  strcpy(hdr.gen.ScanNum, "200899999");
  strcpy(hdr.gen.Comment, "fakePulsar");
  hdr.gen.NFiles=1;
  hdr.gen.FileNum=0;
  strcpy(hdr.gen.TapeID, "");
  hdr.gen.NTapes = 0;
  hdr.gen.TapeNum = 0;
  hdr.gen.ParBlkSize = 0;
  hdr.gen.DataBlkSize = outoptions.NrBins*outoptions.NrPulses*sizeof(float);
  hdr.gen.Cluster[0] = TRUE;
  hdr.gen.Cluster[1] = FALSE;
  hdr.gen.Cluster[2] = FALSE;
  hdr.gen.Cluster[3] = FALSE;
  hdr.gen.Cluster[4] = FALSE;
  hdr.gen.Cluster[5] = FALSE;
  hdr.gen.Cluster[6] = FALSE;
  hdr.gen.Cluster[7] = FALSE;
  hdr.gen.DataMJD = 0;
  hdr.gen.DataTime = 0;  
    
  strcpy(hdr.src.Pulsar, "PSR API-7486");
    
  hdr.redn.Raw = 0;
  hdr.redn.MJDint = 53000;
  hdr.redn.MJDfrac = 0.5;
  hdr.redn.DM = 99;
  hdr.redn.NFreqs = 1;
  hdr.redn.DeltaTime = pulsar.P1/outoptions.NrBins;
  hdr.redn.FreqCent = 1380.000000;
  hdr.redn.DeltaFreq = 640.000000;
  hdr.redn.IsDedisp = TRUE;
  hdr.redn.IsCohDedisp = 0;
  hdr.redn.CohFFTSize = 0;
  hdr.redn.IsPwr = FALSE;
  strcpy(hdr.redn.Zapfile, "");
  hdr.redn.Folded = TRUE;
  hdr.redn.FoldPeriod = pulsar.P1;
  hdr.redn.Polyco = 0;
  hdr.redn.NCoef = 0;
  hdr.redn.PolycoSpan = 0;
  hdr.redn.IsAdjusted = 0;
  hdr.redn.PolycoStored = 0;
  hdr.redn.Bary = 0;
  hdr.redn.OI = TRUE;
  hdr.redn.OQ = FALSE;
  hdr.redn.OU = FALSE;
  hdr.redn.OV = FALSE;
  hdr.redn.OX = FALSE;
  hdr.redn.OY = FALSE;
  hdr.redn.OP = FALSE;
  hdr.redn.OTheta = FALSE;
  hdr.redn.Opoldeg = FALSE;
  hdr.redn.Op = FALSE;
  hdr.redn.Ov = FALSE;
  hdr.redn.TRedn = 0;
  strcpy(hdr.redn.Command, "fakePulsar");
  hdr.redn.RednSoftwareVer = 1;
  hdr.redn.GenType = 4;
    
  hdr.mode.Nr = 1;
  hdr.mode.XPolScaleFac[0] = 0.959410;
  hdr.mode.XPolScaleFac[1] = 1.092279;
  hdr.mode.XPolScaleFac[2] = 1.030063;
  hdr.mode.XPolScaleFac[3] = 0.953323;
  hdr.mode.XPolScaleFac[4] = 0.000000;
  hdr.mode.XPolScaleFac[5] = 0.000000;
  hdr.mode.XPolScaleFac[6] = 0.000000;
  hdr.mode.XPolScaleFac[7] = 0.000000;
  hdr.mode.ActCluster[0] = 1;
  hdr.mode.ActCluster[1] = 1;
  hdr.mode.ActCluster[2] = 1;
  hdr.mode.ActCluster[3] = 1;
  hdr.mode.ActCluster[4] = 0;
  hdr.mode.ActCluster[5] = 0;
  hdr.mode.ActCluster[6] = 0;
  hdr.mode.ActCluster[7] = 0;
  hdr.mode.FIRFactor = 0;
  hdr.mode.NSHARCsAdded = 1;
  hdr.mode.NSampsAdded = 64; 
  hdr.mode.NFreqInFile = 64;
  hdr.mode.Tsamp = 409600;
  hdr.mode.Iout = TRUE;
  hdr.mode.Qout = TRUE;
  hdr.mode.Vout = TRUE;
  hdr.mode.Uout = TRUE;
  hdr.mode.Xout = FALSE;
  hdr.mode.Yout = FALSE;
  hdr.mode.NDMs = 0;        
  hdr.mode.DM[0] = 0.000000;
  hdr.mode.DM[1] = 0.000000;
  hdr.mode.DM[2] = 0.000000;
  hdr.mode.DM[3] = 0.000000;
  hdr.mode.DM[4] = 0.000000;
  hdr.mode.DM[5] = 0.000000;
  hdr.mode.DM[6] = 40.000000;
  hdr.mode.DM[7] = 0.000000;
  hdr.mode.RM=0;     
  hdr.mode.DC_Dynamic = TRUE;
  hdr.mode.ScaleDynamic = TRUE;
  hdr.mode.AdjustInterval=40; 
  hdr.mode.FloatsOut=FALSE;    
  hdr.mode.BitsPerSamp = 4;
  hdr.mode.SigmaRange = 3;
    
  pwheader(&hdr,*fstack);

  if(application.verbose_state.verbose) fprintf(stderr, "Header written to %s\n", filename);
  return 1;
}
*/

void interpolateP4_core(float t, float P4, float fP4, float t4, float phi4, float *phi, double *last_phi)
{
  int i;
  double t2;
  if(fP4 <= 0)
    /* Get the rotation position of the carousel for requested time sample */
    *phi = 2.0*M_PI*t/P4;
  else {      /* If variable, it should be nummerically integrated */
    /*
      phi = 2 Pi t / P4 = 2 Pi Int_0^t (dt/P4)
      P4 = P4 + deltaP4 sin (2Pi t/t4 + phi4)
    */
    t2 = carouselRotation_last_t;
    for(i = 0; i < NrCourouselIntegrationSteps; i++) {
      t2 += (t-carouselRotation_last_t)/(double)(NrCourouselIntegrationSteps);
      *last_phi += (2*M_PI*(t-carouselRotation_last_t)/(P4*(1.0+fP4*sin(2.0*M_PI*t2/t4 + phi4))))/(double)(NrCourouselIntegrationSteps);
      }
      *phi = *last_phi; 
    }
}

void interpolateP4(float t, PulsarDefinitionStruct *pulsar, OutputOptionsStruct outputOpt)
{
  int nc;

  if(reset_num_integration_p4 == 1) {
    carouselRotation_last_t = outputOpt.tstart;
    for(nc = 0; nc < pulsar->nrcones; nc++)
      carouselRotation_last_cone_phi[nc] = 0;
    for(nc = 0; nc < pulsar->nrwaves; nc++)
      carouselRotation_last_wave_phi[nc] = 0;
    for(nc = 0; nc < pulsar->nroscillations; nc++)
      carouselRotation_last_oscillation_phi[nc] = 0;
    reset_num_integration_p4 = 0;
  }
  for(nc = 0; nc < pulsar->nrcones; nc++)
    interpolateP4_core(t, pulsar->cone[nc].P4, pulsar->cone[nc].fP4, pulsar->cone[nc].t4, pulsar->cone[nc].phi4, &(pulsar->cone[nc].phi), &(carouselRotation_last_cone_phi[nc]));
  for(nc = 0; nc < pulsar->nrwaves; nc++)
    interpolateP4_core(t, pulsar->wave[nc].P4, pulsar->wave[nc].fP4, pulsar->wave[nc].t4, pulsar->wave[nc].phi4, &(pulsar->wave[nc].phi), &(carouselRotation_last_wave_phi[nc]));
  for(nc = 0; nc < pulsar->nroscillations; nc++)
    interpolateP4_core(t, pulsar->oscillation[nc].P4, pulsar->oscillation[nc].fP4, pulsar->oscillation[nc].t4, pulsar->oscillation[nc].phi4, &(pulsar->oscillation[nc].phi), &(carouselRotation_last_oscillation_phi[nc]));
  
  carouselRotation_last_t = t;
}

float CalculateEmission_xy(float x, float y, PulsarDefinitionStruct *pulsar, verbose_definition verbose)
{
  float I, I2, x2, y2, theta, phi;
  /* Calculate contribution to the emission of every spark. The sparks
     have a gaussian distribution of emission. If no polar cap is
     defined, set the intensity to 1. */
  if(pulsar->nrcones == 0 && pulsar->nrwaves == 0 && pulsar->nroscillations == 0 && pulsar->nrmicrostructure == 0) {
    I = 1;
  }else {
    I = 1;
    I = calcCarouselEmission(pulsar, x, y, verbose);
    I += calcWaveEmission(*pulsar, x, y);
    I += calcOscillationEmission(*pulsar, x, y);
    I += calcMicroStructureEmission(*pulsar, x, y);

    /* Add emission of other pole */
    theta = sqrt(x*x+y*y);
    phi = atan2(x,y);   /* Didn't think about this, not sure about phase relation between two carousels */
    theta = M_PI-theta;
    x2 = theta*sin(phi);
    y2 = theta*cos(phi);
    I += calcCarouselEmission(pulsar, x2, y2, verbose);
  }
  I2 = calcPatchEmission(*pulsar, x, y);
  return I*I2;
}      

float CalculateEmission(float t, PulsarDefinitionStruct *pulsar, OutputOptionsStruct outputOpt, verbose_definition verbose)
{
  float x, y, I;

  /* Get the position of the vector from the pulsar to Earth with
     respect to the magnetic axis (in polar coordinates). This is then
     converted in (x,y), which are the cartesian coordinates of the
     l.o.s. with respect to the magnetic axis at unit height. So x and
     y corresponds to the cartesian coordinates in a Desh&Ranking like
     emission map. */
  Transformation2(t, *pulsar, &x, &y);

  interpolateP4(t, pulsar, outputOpt);

  I = CalculateEmission_xy(x, y, pulsar, verbose);


  /*
  if(subpulseModel == 0) {                // Rotating carousel
    for(n = 0; n < NrSparks; n++) {
      r1 = n*2.0*M_PI/(float)NrSparks + dfi;
      x1 = sin(rho)*cos(r1);
      y1 = sin(rho)*sin(r1);
      rr = (x1-x)*(x1-x)+(y1-y)*(y1-y);
      I += 2.0*sqrt(log(2))*exp(-4*rr*log(2)/(SparkWidth*SparkWidth))/(SparkWidth*sqrt(M_PI));
    }
  }else {                                // Radially propagating waves 
    rr = (x)*(x)+(y)*(y);
    //	I = sin(sqrt(rr)*2*M_PI*(1.1)/rho-5.1*t*2*M_PI/P1); 
    I = sin(sqrt(rr)*2*M_PI*(0.9)/rho-5.03*t*2*M_PI/P1); //  1702 like 
    // I = sin(sqrt(rr)*2*M_PI*(0.001)/rho-10.05*t*2*M_PI/P1);  // broad pulses 
    I = sin(sqrt(rr)*2*M_PI*(0.001)/rho-20.05*t*2*M_PI/P1);  //narrow pulses
    //			I = sin(sqrt(rr)*2*M_PI*(0.9)/rho-0.99*t*2*M_PI/P1);    // 0815+0939
    //I = sin(sqrt(rr)*2*M_PI*(0.001)/rho-20.05*t*2*M_PI/P1);
    I *= I*I*I*I*I*I*I;
    x1 = sin(rho)*cos(phi);
    y1 = sin(rho)*sin(phi);
    rr = (x1-x)*(x1-x)+(y1-y)*(y1-y);
    I *= 2.0*sqrt(log(2))*exp(-4*rr*log(2)/(SparkWidth*SparkWidth))/(SparkWidth*sqrt(M_PI));
  }
  */
  return I;
}


void longlat2xy(float longitude, float latitude, float *x, float *y, float alpha)
{
  double x2, y2, z2, theta, phi, r;

  spherical2Cartesian(1, 0.5*M_PI-latitude, longitude, &x2, &y2, &z2);
  rotateY_d(1, alpha, &x2, &y2, &z2);
  cartesian2spherical(&r, &theta, &phi, x2, y2, z2);

  *x = theta*sin(phi);
  *y = -theta*cos(phi);
}

void xy2longlat(float *longitude, float *latitude, float x, float y, float alpha)
{
  double x2, y2, z2, theta, phi, r;

  //  int ok;

  alpha *= 180.0/M_PI;

  theta = sqrt(x*x+y*y);
  phi = atan2(x, -y);
  r = 1;
  //  if(theta < 0.01)
  //    ok = 1;
  //  else
  //    ok = 0;
  //  ok = 0;
  //  if(ok) printf("XXXXX x=%f y=%f -> theta=%f phi=%f\n", x, y, theta*180/M_PI, phi*180/M_PI);
  spherical2Cartesian(1, theta, phi, &x2, &y2, &z2);
  //  if(ok) printf("      x=%lf y=%lf z=%lf\n", x2, y2, z2);
  rotateY_d(1, alpha, &x2, &y2, &z2);
  //  if(ok) printf("      x=%lf y=%lf z=%lf\n", x2, y2, z2);
  cartesian2spherical(&r, &theta, &phi, x2, y2, z2);
  //  if(ok) printf("      theta=%f phi=%f\n", theta*180/M_PI, phi*180/M_PI);

  *latitude = 0.5*M_PI-theta;
  *longitude = phi;
  //  if(ok) printf("      long=%f lat=%f\n", (*longitude)*180/M_PI, (*latitude)*180/M_PI);
}

PulsarDefinitionStruct *mapfunction_pulsar;


// Not used?????
float mapfunction(float longitude, float latitude)
{
  float x, y;
  verbose_definition noverbose;
  cleanVerboseState(&noverbose);

  longlat2xy(longitude, latitude, &x, &y, -(mapfunction_pulsar->alpha)*180.0/M_PI);

  //    return latitude*180.0/M_PI;
  //  return longitude*180/M_PI;
  return CalculateEmission_xy(x, y, mapfunction_pulsar, noverbose);
}

int showmap(PulsarDefinitionStruct *pulsar, int nrx, int nry, int nrfieldlines, float markphase, char *mapdevice, int writefile, int projection, float projection_long_rot, float projection_lat_rot)
{
  float *cmap, dx, x, y, I, t, lon, lat, phi, weight;
  long i, j;
  int n, side, side2, oldside, oldside2, first;
  unsigned char c;
  FILE *fout;

  if(application.verbose_state.verbose) fprintf(stderr, "Going to make emission map.\n");
  cmap = (float *)malloc(nrx*nry*sizeof(float));
  if(cmap == NULL) {
    printerror(application.verbose_state.debug, "Cannot allocate memory");
    return 0;
  }


  if(writefile) {
    fout = fopen("dump.dat", "w+");
    if(fout == NULL) {
      printerror(application.verbose_state.debug, "Cannot open dump.dat");
    }
  }

  /* Find out how much we can zoom in on polar region plot so we have
     all features in view. */
  dx = 0;
  if(pulsar->nrpatches != 0) {
    for(n = 0; n < pulsar->nrpatches; n++) {
      x = pulsar->patch[n].x+0.5*pulsar->patch[n].fwhm;
      if(fabs(x) > dx)
	dx = fabs(x);
      x = pulsar->patch[n].x-0.5*pulsar->patch[n].fwhm;
      if(fabs(x) > dx)
	dx = fabs(x);
      y = pulsar->patch[n].y+0.5*pulsar->patch[n].fwhm;
      if(fabs(y) > dx)
	dx = fabs(y);
      y = pulsar->patch[n].y-0.5*pulsar->patch[n].fwhm;
      if(fabs(y) > dx)
	dx = fabs(y);
    }
  }else {
    for(n = 0; n < pulsar->nrcones; n++) {
      x = pulsar->cone[n].rho+5*pulsar->cone[n].halfwidth;
      if(x > dx)
	dx = x;
    }
    for(n = 0; n < pulsar->nrwaves; n++) {
      x = 5*pulsar->wave[n].lambda;
      if(x > dx)
	dx = x;
    }
    for(n = 0; n < pulsar->nroscillations; n++) {
      x = 5*2*M_PI/pulsar->oscillation[n].l;
      if(x > dx)
	dx = x;
    }
  }


  if(projection) {
    mapfunction_pulsar = pulsar;
    make_projection_map(cmap, nrx, nry, 0, projection_long_rot, projection_lat_rot, &mapfunction, projection);
  }else {
    for(j = 0; j < nry; j++) {
      for(i = 0; i < nrx; i++) {
	x = dx*(2*i/(float)(nrx)-1.0);
	y = dx*(2*j/(float)(nry)-1.0);
	//	fprintf(stderr, "%f %f (i=%ld/%d) j=%ld/%d)\n", x, y, i, nrx, j, nry);
	I = CalculateEmission_xy(x, y, pulsar, application.verbose_state);
	cmap[nrx*j+i] = I;
	c = I*1;
	if(writefile)
	  fwrite(&c, 1, 1, fout); 
      }
    }
    if(writefile) {
      fclose(fout);
      printf("Emission map written as dump.dat\n");
    }
  }

  //  fprintf(stderr, "OK so far\n");

  pgplot_options_definition pgplot_options;
  pgplot_clear_options(&pgplot_options);
  if(projection)
    pgplot_options.viewport.windowwidth = 2*nrx;
  else
    pgplot_options.viewport.windowwidth = nrx;
  pgplot_options.viewport.windowheight = nry;
  strcpy(pgplot_options.viewport.plotDevice, mapdevice);
  pgplot_options.viewport.dontclose = 1; 
  strcpy(pgplot_options.box.title, "Emission map");
  /*  pgplotMap(&viewport, cmap, nrx, nry, -dx, dx, -dx, dx, -dx, dx, -dx, dx, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, "X [rad]", "Y [rad]", "Emission map", 1, 1, 1, 0, 0, 1.0, 1.0, 1, 0, 0, 0, 0, 0, 0, 0, application.debug); */
  if(projection == 0) {
    strcpy(pgplot_options.box.xlabel, "X [rad]");
    strcpy(pgplot_options.box.ylabel, "Y [rad]");
    pgplotMap(&pgplot_options, cmap, nrx, nry, -dx, dx, -dx, dx, -dx, dx, -dx, dx, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, application.verbose_state);
  }else {
    if(projection == 3) {
      strcpy(pgplot_options.box.xlabel, "longitude");
      strcpy(pgplot_options.box.ylabel, "latitude");
      pgplotMap(&pgplot_options, cmap, nrx, nry, -180, 180, -180, 180, -90, 90, -90, 90, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, application.verbose_state);
    }else {
      strcpy(pgplot_options.box.xlabel, "X");
      strcpy(pgplot_options.box.ylabel, "Y");
      pgplotMap(&pgplot_options, cmap, nrx, nry, -2.25, 2.25, -2.25, 2.25, -1.125, 1.125, -1.125, 1.125, PPGPLOT_GRAYSCALE, 0, 0, 0, NULL, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, application.verbose_state);
    }
    printf("Press a key to plot rotational grid\n");
    pgetch();
    /* Draw magnetic coordinates first */
    // Doesn't work with rotation in longitude direction. Order of rotations is wrong way around.
    //    ppgsci(4);
    //    drawSphericalGrid(30, 30, projection_long_rot*180.0/M_PI, (projection_lat_rot+pulsar->alpha)*180.0/M_PI, 1, projection);

    /* Draw rotational coordinates */
    ppgsci(1);
    drawSphericalGrid(30, 30, projection_long_rot*180.0/M_PI, projection_lat_rot*180.0/M_PI, 1, projection);

    printf("Press a key to plot magnetic grid\n");
    pgetch();
  }

  if(nrfieldlines > 0) {
    ppgsls(2);
    for(i = 0; i < nrfieldlines; i++) {
      ppgmove(0, 0);
      ppgdraw(10*cos(2.0*M_PI*i/(float)nrfieldlines), 10*sin(2.0*M_PI*i/(float)nrfieldlines));
    }
    ppgsls(1);
  }

  /* Draw line of sight */
  ppgsci(2);
  first = 1;
  oldside = -1;
  for(i = 0; i <= 1000; i++) {
    t = i*pulsar->P1/1000.0;
    if(Transformation2(t, *pulsar, &x, &y) == 1)
      ppgsci(2);
    else
      ppgsci(3);
    if(projection) {
      xy2longlat(&lon, &lat, x, y, pulsar->alpha);
      if(projection == 1) {
	projectionHammerAitoff_xy(lon, lat, projection_long_rot, projection_lat_rot, &x, &y);
	if(x < 0)
	  side = -1;
	else
	  side = 1;
	if(side != oldside)
	  first = 1;
	oldside = side;
      }else if(projection == 2) {
	side = projection_sphere_xy(lon, lat, projection_long_rot, projection_lat_rot, &x, &y, &weight);
	if(side != oldside)
	  first = 1;
	oldside = side;
      }else if(projection == 3) {
	projection_longlat_xy(lon, lat, projection_long_rot, projection_lat_rot, &x, &y);
	if(x < 0)
	  side = -1;
	else
	  side = 1;
	if(side != oldside)
	  first = 1;
	oldside = side;
      }
    }
    if(first) {
      ppgmove(x, y);
      first = 0;
    }else {
      ppgdraw(x, y);
    }
  }
  ppgsci(1);

  /* Draw magnetic meridians */
  ppgsci(4);
  first = 1;
  oldside = -1;
  for(phi = 0; phi <= 2.0*M_PI+0.001; phi += 30.0*M_PI/180.0) {
    for(i = 0; i <= 1000; i++) {
      x = cos(phi)*2*M_PI*i/1000.0;
      y = sin(phi)*2*M_PI*i/1000.0;
      if(projection) {
	xy2longlat(&lon, &lat, x, y, pulsar->alpha);
	if(projection == 1) {
	  projectionHammerAitoff_xy(lon, lat, projection_long_rot, projection_lat_rot, &x, &y);
	  if(x < 0)
	    side = -1;
	  else
	    side = 1;
	  if(side != oldside)
	    first = 1;
	  oldside = side;
	}else if(projection == 2) {
	  side = projection_sphere_xy(lon, lat, projection_long_rot, projection_lat_rot, &x, &y, &weight);
	  if(side != oldside)
	    first = 1;
	  oldside = side;
	}else if(projection == 3) {
	  projection_longlat_xy(lon, lat, projection_long_rot, projection_lat_rot, &x, &y);
	  if(x < 0)
	    side = -1;
	  else
	    side = 1;
	  if(y < 0)
	    side2 = 0;
	  if(y > 0)
	    side2 = -1;
	  if(side != oldside || side2 != oldside2)
	    first = 1;
	  oldside = side;
	  oldside2 = side2;
	}
      }
      if(first) {
	ppgmove(x, y);
	first = 0;
      }else {
	ppgdraw(x, y);
      }
    }
  }
  ppgsci(1);

  /* Draw magnetic equal latitude lines */
  ppgsci(4);
  first = 1;
  oldside = -1;
  oldside2 = -1;
  for(phi = 0; phi <= M_PI+0.001; phi += 30.0*M_PI/180.0) {
    for(i = 0; i <= 1000; i++) {
      x = phi*cos(2*M_PI*i/1000.0);
      y = phi*sin(2*M_PI*i/1000.0);
      if(projection) {
	xy2longlat(&lon, &lat, x, y, pulsar->alpha);
	if(projection == 1) {
	  projectionHammerAitoff_xy(lon, lat, projection_long_rot, projection_lat_rot, &x, &y);
	  if(x < 0)
	    side = -1;
	  else
	    side = 1;
	  if(side != oldside)
	    first = 1;
	  oldside = side;
	}else if(projection == 2) {
	  side = projection_sphere_xy(lon, lat, projection_long_rot, projection_lat_rot, &x, &y, &weight);
	  if(side != oldside)
	    first = 1;
	  oldside = side;
	}else if(projection == 3) {
	  projection_longlat_xy(lon, lat, projection_long_rot, projection_lat_rot, &x, &y);
	  if(x < 0)
	    side = -1;
	  else
	    side = 1;
	  if(y < 0)
	    side2 = 0;
	  if(y > 0)
	    side2 = -1;
	  if(side != oldside || side2 != oldside2)
	    first = 1;
	  oldside = side;
	  oldside2 = side2;
	}
      }
      if(first) {
	ppgmove(x, y);
	first = 0;
      }else {
	ppgdraw(x, y);
      }
    }
  }
  ppgsci(1);

  /* Place a marker to mark a certain pulse longitude */
  if(markphase > -1000) {
    t = markphase*pulsar->P1/360.0;
    Transformation2(t, *pulsar, &x, &y);
    if(projection) {
      xy2longlat(&lon, &lat, x, y, pulsar->alpha);
      if(projection == 1) {
	projectionHammerAitoff_xy(lon, lat, projection_long_rot, projection_lat_rot, &x, &y);
      }else if(projection == 2) {
	projection_sphere_xy(lon, lat, projection_long_rot, projection_lat_rot, &x, &y, &weight);
      }
    }

    ppgsci(6);
    ppgsch(2);
    ppgslw(10);
    ppgmove(x, y);
    ppgpt1(x, y, 17);
    ppgsci(1);
    ppgsch(1);
    ppgslw(1);
  }


  ppgclos();
  free(cmap);
  if(application.verbose_state.verbose) fprintf(stderr, "Making emission map finished.\n");
  return 0;
}


/* 
If reqOutput = 1:

  For a given alpha_input and beta_input, calculate the s and phi
  (footprint parameters for a pure dipole field) corresponding to the
  pulse longitude that is deltaPulseLongitude1 away from the steepest
  gradient of the PA swing. If IP is set, the alpha/beta values are
  first transformed to the other pole and the offsets are the offsets
  w.r.t. the steepest gradient in the IP pole. (deltaPulseLongitude2
  and trackheight) are not used.

If reqOutput = 2: 

  For a given alpha_input and beta_input, calculate the s and phi
  (footprint parameters for a pure dipole field) as a function of the
  pulse longitude offset from the steepest gradient of the PA swing
  (between deltaPulseLongitude1 and deltaPulseLongitude2) for a given
  emission height trackheight as fraction of the light cylinder
  radius. If IP is set, the alpha/beta values are first transformed to
  the other pole and the offsets are the offsets w.r.t. the steepest
  gradient in the IP pole.

Output is sent to the stdout.

*/
void calculateLOS(int reqOutput, int IP, float alpha_input, float beta_input, float deltaPulseLongitude1, float deltaPulseLongitude2, float trackheight)
{
  long nrpoints = 100000, i;
  float alpha, beta, zeta, deltaPulseLongitude, r, fd, w, rho, s, phi, x, y; //, sign

  for(i = 0; i < nrpoints; i++) {    /* The number of points in the line of sight */
    if(reqOutput == 1) {
      deltaPulseLongitude = deltaPulseLongitude1;
      /* Emission height (fraction RLC) */
      r = i/(float)nrpoints;
    }else {
      deltaPulseLongitude = deltaPulseLongitude1 + (deltaPulseLongitude2 - deltaPulseLongitude1)*i/(float)nrpoints;
      r = trackheight;
    }
    if(IP == 0) {
      beta = beta_input;
      alpha = alpha_input;
    }else {
      beta = 2.0*alpha_input+beta_input-180.0;
      alpha = 180.0-alpha_input;
    }
    
    deltaPulseLongitude *= M_PI/180.0;
    alpha *= M_PI/180.0;
    beta *= M_PI/180.0;
    zeta = alpha+beta;


    /* Shift of PA and emission w.r.t. the fiducial point */
    fd = 2.0*r;
    /* fiducial point is -fd of steepest gradient PA and centre
       profile is -2*fd from steepest gradient. Total pulse width is
       therefore: */
    w = 2*fd+deltaPulseLongitude;
    //    if(w >= 0)
    //      sign = 1;
    //    else
    //      sign = -1;
    /*      rho = sqrt(w*w+beta*beta); */
    /* Multiply by 2 to get the full width */
    w *= 2.0;
    /* from thesis */
    rho = acos(cos(alpha)*cos(zeta)+sin(alpha)*sin(zeta)*cos(w/2.0));
    /* rho = 1.5*s*sqrt(r/rlc) (e.g. page 3 dyks 2008) 
       s = rho/(1.5*sqrt(r/rlc))
    */
    float theta;
    /* Page 68 handbook */
    theta = atan(-3/(2.0*tan(rho))+sqrt(2+(3/(2.0*tan(rho)))*(3/(2.0*tan(rho)))));
    s = sin(theta)/sqrt(r);
    /*            s = rho/(1.5*sqrt(r)); */
    /* This is the WANG type of solution
       p = 0.5*(alpha + zeta + rho);
       phi = M_PI-2.0*atan(sqrt(sin(p-alpha)*sin(p-rho)/(sin(p-zeta)*sin(p))));
       phi *= sign;
    */
    phi = asin(sin(zeta)*sin(w/2.0)/sin(rho));
    /* phi in the formula is the whole phi range between leading and
       trailing edge. We're here interested in the half angle. */
    /*      phi *= 0.5; */
    if(reqOutput == 1) {
      x = 180.0*rho*sin(phi)/M_PI;
      y = -180.0*rho*cos(phi)/M_PI;
    }else {
      x = s*sin(phi);
      y = -s*cos(phi);
    }
    if(reqOutput == 1) {
      if(i == 0)
	printf("#r (fraction of RLC), footprint parameter s (fraction last open field line), phi (degrees) and X and Y (magnetic coordinates emission point in degrees)\n");
      printf("%f %f %f %f %f\n", r, s, phi*180/M_PI, x, y);
    }else {
      if(i == 0)
	printf("#Pulse longitude (w.r.t. steepest part PA swing) footprint parameter s (fraction last open field line), phi (degrees) and X and Y (magnetic coordinates footprint in fraction RLC)\n");
      printf("%f %f %f %f %f\n", deltaPulseLongitude, s, phi*180/M_PI, x, y);
    }
      /*      if(i == 8000 && (j == 0 || j == 2)) {
	fprintf(stderr, "r=%f\n", r);
	fprintf(stderr, "w/2=%f\n", w*90/M_PI);
	fprintf(stderr, "phi=%f\n", phi*180/M_PI);
	fprintf(stderr, "rho=%f\n", rho*180/M_PI);
	fprintf(stderr, "s=%f\n", s);
	fprintf(stderr, "fakePulsar -np 10 -phase0 0 -a %f -b %f -P 0.197 -patch \"1 %f %f 1\" ", alpha*180/M_PI, beta*180/M_PI, x, y);
      }
      if(i == 8000 && (j == 1 || j == 3)) {
	fprintf(stderr, "-patch \"1 %f %f 1\" -map -mapd 1/xs -prof -profd 2/xs\n", x, y);
	fprintf(stderr, "w/2=%f\n", w*90/M_PI);
	fprintf(stderr, "phi=%f\n", phi*180/M_PI);
	fprintf(stderr, "rho=%f\n", rho*180/M_PI);
	fprintf(stderr, "s=%f\n\n", s);
	}*/
  }
}

