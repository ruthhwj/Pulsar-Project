//START REGION RELEASE
//START REGION DEVELOP

// This is your GSL version number you're using, multiplied by 100. Should be an integer value
#ifndef GSL_VERSION_NUMBER
  //#define GSL_VERSION_NUMBER   115
#endif

// Comment this line out if you don't want to use Numerical Recipes
//#define NRAVAIL     1  

//START REGION RELEASE
#define DM_CONST                4.148808e3      /* Eq: dt(sec) = DM_CONST*DM/F(MHz)**2, used to be 4.15e3 */

// Maximum number of regions to be definable as for instance on-pulse regions.
#define MAX_pulselongitude_regions 200
// Maximum number of components in a von Mises decomposition.
#define maxNrVonMisesComponents   100
// The maximum length of a word in a string when handled by pickWordFromString()
#define MaxPickWordFromString_WordLength 1000
// Maximum length of a filename supported in various places
#define MaxFilenameLength 10000
// Maximum length of a pgplot device
#define MaxPgplotDeviceLength 2000
// Maximum string length of labels etc.
#define MaxStringLength 10000
// Maximum coefficients of the polynomial fits for receiver models
#define MaxNrfitReceiverModelFitParameters 11

/* The file formats which are presently (partially) working in this library. */
//START REGION RELEASE
#define PUMA_format                  1
//START REGION DEVELOP
#define AO_ASCII_1_format            2
#define PRESTO_format                3
#define PARKESFB_format              4
//START REGION RELEASE
#define PSRCHIVE_ASCII_format        5
#define EPN_format                   6
#define FITS_format                  7
#define SIGPROC_format               8
//START REGION RELEASE
#define PPOL_format                  9
#define PPOL_SHORT_format            10
#define SIGPROC_ASCII_format         11
//START REGION DEVELOP
#define GMRT_ASCII_format            12
#define AO_ASCII_2_format            13
//#define LAST_VALID_DATA_FORMAT       12
//START REGION RELEASE
#define PSRSALSA_BINARY_format       20
#define MEMORY_format                99
//START REGION DEVELOP

//START REGION RELEASE

#define PPGPLOT_GRAYSCALE 1
#define PPGPLOT_INVERTED_GRAYSCALE 2
#define PPGPLOT_RED 3
#define PPGPLOT_INVERTED_RED 4
#define PPGPLOT_GREEN 5
#define PPGPLOT_INVERTED_GREEN 6
#define PPGPLOT_BLUE 7
#define PPGPLOT_INVERTED_BLUE 8
#define PPGPLOT_CYAN 9
#define PPGPLOT_INVERTED_CYAN 10
#define PPGPLOT_HEAT 20
#define PPGPLOT_INVERTED_HEAT 21
#define PPGPLOT_COLD 22
#define PPGPLOT_INVERTED_COLD 23
//START REGION DEVELOP
#define PPGPLOT_PLASMA 24
#define PPGPLOT_INVERTED_PLASMA 25
#define PPGPLOT_FOREST 26
#define PPGPLOT_INVERTED_FOREST 27
#define PPGPLOT_ALIEN_GLOW 28
#define PPGPLOT_INVERTED_ALIEN_GLOW 29
//START REGION RELEASE
#define PPGPLOT_HEAT2 30
#define PPGPLOT_INVERTED_HEAT2 31
#define PPGPLOT_HEAT3 32
#define PPGPLOT_INVERTED_HEAT3 33
#define PPGPLOT_HEAT4 34
#define PPGPLOT_INVERTED_HEAT4 35
#define PPGPLOT_DIVREDBLUE 36
#define PPGPLOT_INVERTED_DIVREDBLUE 37
#define PPGPLOT_INFERNO 38
//START REGION DEVELOP

//START REGION RELEASE

/* If adding a type, update returnGenType_str() and printHeaderGentypeOptions() in psrio.c */
#define GENTYPE_UNDEFINED        0
#define GENTYPE_PROFILE          1  /* Intensity vs pulse phase, otherwise the same as GENTYPE_SUBINTEGRATIONS. Should have one subint, but can have frequency resolution */
#define GENTYPE_PULSESTACK       2  /* Each subint is a single pulse */
#define GENTYPE_SUBINTEGRATIONS  3  /* Each subint are multiple pulses, tobs must be set in header */
//#define GENTYPE_FREQPHASE        4  /* Intensity as function of pulse phase and frequency, otherwise the same as profile */
#define GENTYPE_SEARCHMODE       4  /* Data is not folded at period, i.e. a time series. */
#define GENTYPE_BANDPASS         5  /* Intensity as function of frequency. Each subint can only have one bin. */
#define GENTYPE_DYNAMICSPECTRUM  6  /* Pulse intensity or s2n as function of time (bins) and frequncy. Tsamp should be the pulse period and there should be only one subintegration. */
#define GENTYPE_PENERGY          7  /* There should be 7 polarizations: Peak intensity (on & off-pulse) Integrated energy (on & off-pulse) RMS (on & off-pulse) and S/N. Tsamp should be the pulse period and the "subintegrations" correspond to polarizations in the original file. */
#define GENTYPE_POLNCAL         10  /* Polarization calibration cal signal */
#define GENTYPE_LRFS            20  /* yrange should be in cpp */
#define GENTYPE_2DFS            21  /* x and yrange should be in cpp */
#define GENTYPE_S2DFSP3         22  /* Sliding S2DFS P3 map. The yrange should be in cpp */
#define GENTYPE_S2DFSP2         23  /* Sliding S2DFS P2 map. The yrange should be in cpp */
#define GENTYPE_P3FOLD          25  /* yrange should be in pulse numbers */
#define GENTYPE_POLARMAP        26  /* xrange & yrange should be in degrees */
#define GENTYPE_HRFS_UNFOLDED   30  /* xrange should be in cpp */
#define GENTYPE_HRFS            31  /* xrange should be in cpp */
#define GENTYPE_LRCC            33  /* yrange should be in degrees */
#define GENTYPE_LRAC            34  /* yrange should be in lags */
#define GENTYPE_RMMAP           50  /* Faraday depth as function of pulse longitude (bins), yrange should be inrad/m^2 */
//#define GENTYPE_ILVPA          100  /* Data contains: Stokes I, L, V, PA and PAerror, otherwise the same as GENTYPE_SUBINTEGRATIONS */
#define GENTYPE_PADIST         101  /* Data contains: Pulse number is PA bin */
#define GENTYPE_ELLDIST        102  /* Data contains: Pulse number is ellipticity bin */
#define GENTYPE_RECEIVERMODEL  200  /* Data contains npol parameters. If nbin=2 then there is also an errorbar */
#define GENTYPE_RECEIVERMODEL2 201  /* Same as GENTYPE_RECEIVERMODEL, except last 2 polarizations contain the chisq and nfree */

#define FOLDMODE_UNKNOWN        -1  // Unknown mode of folding
#define FOLDMODE_FIXEDPERIOD     1  // Fixed period folding throughout data, currently the only model supported. Could implement ephemeris mode, period per subint, etc.

#define TSAMPMODE_UNKNOWN        -1  // Unknown mode of sampling
#define TSAMPMODE_FIXEDTSAMP      1  // Fixed sampling time throughout data, currently the only model supported. Could implement ephemeris mode, period per subint, etc.
#define TSAMPMODE_LONGITUDELIST   2  // Each bin has a specified pulse longitude, in degrees. Stored in the array tsamp_list, which has NrBins elements.

#define TSUBMODE_UNKNOWN        -1  // Unknown mode of making subintegrations
#define TSUBMODE_FIXEDTSUB       1  // A fixed subintegration time is defined, so each subintegration has the same length.
#define TSUBMODE_TSUBLIST        2  // Each subint has a defined duration, stored in the array 

/* 0 = undefined 1=STOKES, 2=AABBCRCI */
#define POLTYPE_UNKNOWN          -1  // Unknown meaning of polarization channels
#define POLTYPE_STOKES            1  // Stokes parameters
#define POLTYPE_COHERENCY         2  // Coherency parameters
#define POLTYPE_ILVPAdPA          3  // I, L, V, Pa and error on PA
#define POLTYPE_PAdPA             4  // Pa and error on PA
#define POLTYPE_ILVPAdPATEldEl    5  // I, L, V, Pa, error on PA, Totpol, Ellipticity, error on ellipticity

#define FREQMODE_UNKNOWN        -1  // Unknown mode of specifying frequency labelling
#define FREQMODE_UNIFORM         1  // A uniform sampling of frequency channels, i.e. a fixed channel bw + centre freq.
#define FREQMODE_FREQTABLE       2  // Each subint and freq channel has a defined frequency, stored in the array 

/* +1 is LIN:A=X,B=Y, +2 is CIRC:A=L,B=R (I)  (or -1,-2 to change handiness) 0=undefined*/
#define FEEDTYPE_UNKNOWN         0
#define FEEDTYPE_LINEAR          1  // linear:   A=X,B=Y
#define FEEDTYPE_CIRCULAR        2  // circular: A=L,B=R
#define FEEDTYPE_INV_LINEAR     -1  // linear:   A=Y,B=X
#define FEEDTYPE_INV_CIRCULAR   -2  // circular: A=R,B=L

//START REGION DEVELOP

#define SHAPEPAR_PEAKPHASE           1  // Phase of the peak
#define SHAPEPAR_PEAKAMP             2  // Amplitude of the peak
#define SHAPEPAR_W10                 3  // Width at 10% of maximum in phase
#define SHAPEPAR_W25                 4  // Width at 25% of maximum in phase
#define SHAPEPAR_W50                 5  // Width at 50% of maximum in phase
#define SHAPEPAR_W75                 6  // Width at 75% of maximum in phase
#define SHAPEPAR_W90                 7  // Width at 90% of maximum in phase
#define SHAPEPAR_PEAKSEPPHASE        10  // Find the offset from the main peak to another peak separated in the range indicated by the first and second aux parameter in calcVonMisesProfile_shape_parameter().
#define SHAPEPAR_PEAKAMPRATIO        11  // As SHAPEPAR_PEAKSEPPHASE, but calculate ratio of amplitudes
#define SHAPEPAR_PEAKAMPRATIO_RECI   12  // The reciprocal (^-1) of SHAPEPAR_PEAKAMPRATIO
#define SHAPEPAR_W10_MAXAMP          23  // Width at 10% of maximum in phase for that vonMises component which has the maximum amplitude
#define SHAPEPAR_W25_MAXAMP          24  // Width at 25% of maximum in phase for that vonMises component which has the maximum amplitude
#define SHAPEPAR_W50_MAXAMP          25  // Width at 50% of maximum in phase for that vonMises component which has the maximum amplitude
#define SHAPEPAR_W75_MAXAMP          26  // Width at 75% of maximum in phase for that vonMises component which has the maximum amplitude
#define SHAPEPAR_W90_MAXAMP          27  // Width at 90% of maximum in phase for that vonMises component which has the maximum amplitude

#define COMBINERECEIVERMODELS_CHI2     1  // For each channel, store that solution with the smallest chi2
#define COMBINERECEIVERMODELS_APPEND   2  // All receiver solutions are combined in a single file with probably multiply defined frequency channels

//START REGION RELEASE

/* Maximum number of filenames application library can handle */
#define MaxNrApplicationFilenames 1025

#define maxNrRotateStokes 10  // Maximum number of times the -rotateStokes option can be used
//START REGION DEVELOP
//START REGION RELEASE

#define MaxNrFitFunctions  100
#define MaxNrFitParameters  10
//START REGION DEVELOP

/*
   Each function is defined with a number of parameters (fixed values) p1..pn and a
   number of values a1..an (which can be refined).  

   when a function is added the functions, updates should be made in the following functions in fitting.c:
     evaluate_fitfunc_collection,
     mrq_func_internal,
     mrq_func_deriv_internal,
     show_fitfunctions_commandline_options, 
     parse_commandline_fitfunctions,
     print_fitfunctions, 
     countnrparameters_internal_mrq,
     initialize_mrq, 

*/
//START REGION RELEASE
#define FUNC_POLYNOMAL   1      /* a1*x**p1 */
//START REGION DEVELOP
#define FUNC_ATAN        2      /* a1*atan(a2*(x+a3)) */
#define FUNC_SIN         3      /* a1*sin(a2*(x+a3)) */
#define FUNC_COS         4      /* a1*cos(a2*(x+a3)) */
#define FUNC_NORMAL      10     /* a1*exp(-a2*(x-a3)^2) */

//START REGION DEVELOP


#define MaxNrConesRefraction      2

//START REGION DEVELOP
//START REGION RELEASE

/* Somehow this is not always defined on all systems */

#ifndef NAN
  #define NAN (0.0/0.0)
#endif

#ifndef M_PI
  #define M_PI  3.14159265358979323846264338327950288
#endif

// Call without the trailing \n
// If debug_flag is nonzero, the source file and line number are added to the message
#define printerror(debug_flag, ...)	 \
  { fflush(stdout); fprintf_color(stderr, 2, __VA_ARGS__);		\
  if(debug_flag) fprintf_color(stderr, 2, " (message generated in %s line %d)", __FILE__, __LINE__);	\
  fprintf_color(stderr, 0, "\n"); }

// Same as the above, but in different colour
#define printwarning(debug_flag, ...)	 \
  fflush(stdout); fprintf_color(stderr, 7, __VA_ARGS__);			\
  if(debug_flag) fprintf_color(stderr, 7, " (message generated in %s line %d)", __FILE__, __LINE__);	\
  fprintf_color(stderr, 0, "\n")

//START REGION DEVELOP
