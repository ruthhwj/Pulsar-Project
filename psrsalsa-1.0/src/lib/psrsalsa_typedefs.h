//START REGION RELEASE

#ifndef PSRSALSA_TYPEDEFS_LOADED

typedef struct {
  int nrRegions;
  /* Specifies if bins are defined */
  int *bins_defined;
  int *left_bin, *right_bin;
  /* Specifies if frac of baseline are defined */
  int *frac_defined;
  float *left_frac, *right_frac;
}pulselongitude_regions_definition;

typedef struct {
  int verbose;      // Toggle verbose mode
  int debug;        // Toggle debug mode. Among other things, it enables line nr printing in error messages
  int nocounters;   // If set, counters are disabled, which is good for log files
  int indent;       // Request verbose messages to have this nr of spaces in front of output.
}verbose_definition;

//START REGION DEVELOP
//START REGION RELEASE

// Requires in principle unnecesary non-dynamic memory allocations, but all functions were re-written to use pointers only to avoid too much non-dynamic memory use.
typedef struct {
  /* centre is in phase */
  double centre[maxNrVonMisesComponents], concentration[maxNrVonMisesComponents], height[maxNrVonMisesComponents];
  int nrcomponents;
}vonMises_collection_definition;

//START REGION DEVELOP
//START REGION RELEASE

// Requires in principle unnecesary non-dynamic memory allocations, but all functions were re-written to use pointers only to avoid too much non-dynamic memory use.
typedef struct { /* Please change pgplot_clear_viewport_def as well when changing this definition. */
  char plotDevice[MaxPgplotDeviceLength];  /* default = "?" */
  int windowwidth, windowheight;      /* The size of the device (if set) */
  float aspectratio;                  /* if size of the device (resolution) is unset, this aspect ratio is used if > 0 */
  float dxplot;                       /* default = 0, otherwise shifted to right */
  float xsize;                        /* default = 1, otherwise wider by this fraction */
  float dyplot;                       /* default = 0, otherwise shifter up */
  float ysize;                        /* default = 1, otherwise taller by this fraction */
  int noclear;                        /* default = 0, current plot is overwritten (no pgpage) */
  int dontopen;                       /* default = 0, otherwise device is assumed to be opened and selected */
  int dontclose;                      /* default = 0, otherwise request function to not close pgplot device */
}pgplot_viewport_definition;

typedef struct {
  int svp;                                   /* Set if want to call pgsvp */
  float svp_x1, svp_x2, svp_y1, svp_y2;      /* The coordinates for pgsvp, normally generated from numbers inside pgplot_viewport_def */
  int swin;                                  /* Set if want to call pgwin */
  int swin_showtwice;                        /* Set if plotting a map which is plotted twice on top of each other */
  float swin_x1, swin_x2, swin_y1, swin_y2;  /* The coordinates for pgswin */
  float TR[6];                               /* Will be calculated by pgplot_makeframe from internal_pgplot_XXXX values */
}pgplot_frame_def_internal;

// Requires in principle unnecesary non-dynamic memory allocations, but all functions were re-written to use pointers only to avoid too much non-dynamic memory use.
typedef struct {  /* If you change any of this, also change clear_pgplot_box */
  int drawbox;                       /* Draws the axes */
  int box_lw, box_f;                 /* Define line width and font*/
  float box_labelsize;               /* and the size of the numbers */
  float box_xtick, box_ytick;        /* Separation between major tick marks (default is auto) */
  int box_nxsub, box_nysub;          /* The number of subintervals between major tick marks (default is auto) */
  char box_xopt[10], box_yopt[10];   /* Defines the options of pgbox, default is bcnsti for both */

  int drawtitle;                     /* Draws a title on top of plot */
  float title_ch;                    /* Defines the character size of the title */
  int title_lw, title_f;             /* Defines the linewidth and the font */
  char title[MaxStringLength];

  int drawlabels;                    /* Draws the labels along the axes */
  float label_ch;                    /* Define the character size */
  int label_lw, label_f;             /* Defines the linewidth and the font */
  float dxlabel, dylabel;            /* For fine tuning of the placement of the labels */
  char xlabel[MaxStringLength];
  char ylabel[MaxStringLength];
  char wedgelabel[MaxStringLength];
}pgplot_box_definition;

// Requires in principle unnecesary non-dynamic memory allocations, but all functions were re-written to use pointers only to avoid too much non-dynamic memory use.
typedef struct {
  pgplot_viewport_definition viewport;
  pgplot_box_definition box;
}pgplot_options_definition;

//START REGION DEVELOP
//START REGION RELEASE

typedef struct {
  char *timestamp;   // These strings are NULL initially, or point to a string otherwise
  char *cmd;
  char *user;
  char *hostname;
  void *nextEntry;  // Set to NULL if it is the last entry. Is pointing to a new datafile_history_entry_definition struct
}datafile_history_entry_definition;

typedef struct 
{
  // ***************************
  // Some general IO properties
  // ***************************
  FILE *fptr, *fptr_hdr;  // File pointer to the file and possibly to the header file if separate. There are NULL if not defined.
  fitsfile *fits_fptr;    // FITS file pointer, as defined in fitsio.h, used to access fits files
  char *filename;         // The name of the data file. Should be an empty string if not in use (1 byte memory is assigned to hold byte value 0 when structure is initialised with cleanPRSData). When setting the filename, new memory should be allocated. This way of doing things avoids having to define an a-priory string length. You can use the function set_filename_PSRData()
  int format;             // File format, as defined above as PUMA_format etc.
  int version;            // File format version, as used for certain file types
  int opened_flag, enable_write_flag; // Flags to set if file is open and if writing to file is allowed.
  int dumpOnClose;        // In write mode, the individual write statements are buffered in memory and data is written when file is closed. Useful for ascii formats which do not have a predictable file position for each sample and therefore should be written from start to finish in one go.

  // *********************************************************
  // Some general information strings about the observation. 
  // The strings are not necessarily in any particular format.
  // String handeling similar as for filename.
  // *********************************************************
  char *psrname;     // Name of the pulsar, use set_psrname_PSRData().
  char *observatory; // Name of the observatory, use set_observatory_PSRData().
  char *instrument;  // Name of the instrument/backend used, use set_instrument_PSRData().
  char *scanID;      // Additional string identifying observation, besides psrname. Use set_scanid_PSRData().
  char *observer;    // Name of observer, use set_observer_PSRData()
  char *projectID;   // Name of observing project, use set_projectID_PSRData()
  char *institute;   // Name of the institute, use set_institute_PSRData().

  // Some other general information
  double telescope_X, telescope_Y, telescope_Z; // ITRF location of telescope in meters
  //  double telescope_long, telescope_lat;  /* Longitude and latitude of telescope in radians. */
  int NrBits;                 /* NrBits defines the sample size in bits*/
  char isDeDisp, isDeFarad, isDePar, isDebase;    /* -1 = unknown, 0 = dispersion and Faraday rotation and parallactic angle effect and baseline is still present in data, 1 = effects removed. If removed, it is with respect to freq_ref. */
  double dm, rm;              /* DM in pc/cm^3 and RM in radians/meter^2.*/
  double freq_ref;            // The reference frequency in MHz for dedispersion/deFaraday rotation. 1e10=infinite frequency (used to be -1), -2=Undefined
  int feedtype;               // FEEDTYPE_UNKNOWN, FEEDTYPE_LINEAR, FEEDTYPE_CIRCULAR
  int poltype;                // POLTYPE_UNKNOWN, POLTYPE_STOKES, POLTYPE_COHERENCY or POLTYPE_ILVPAdPA
  long double mjd_start;      // The MJD of the first sample (left edge) to be defined, or 0 if undefined
  char cableSwap;             // Cable swap (of two polarizations) occured when recording data? -1/0/1 = Unknown/No/Yes
  char cableSwapcor;          // 1 = Cable swap is corrected for by reordering the data, 0 = not corrected, -1 unknown

  // The dimensions of the data
  long NrSubints, NrBins, NrPols, NrFreqChan;  /* The dimensions of the data. NrSubints is the number of subintegrations, which could be the number of pulses if in single pulse mode. */

  // Information about folding/time resolution
  char isFolded;              // 0=Not folded, i.e. a time series. 1=Folded, -1 is Unknown. Note that if isFolded = 0, you probably want to set gentype=GENTYPE_SEARCHMODE
  char foldMode;              // FOLDMODE_UNKNOWN or FOLDMODE_FIXEDPERIOD
  double fixedPeriod;         // For foldMode=1: the fold period in seconds. The period is assumed to be constant and no PEPOCH is defined, so this is an approximate number. Use get_period() to obtain the period from a data set.
  char tsampMode;             // TSAMPMODE_UNKNOWN, TSAMPMODE_FIXEDTSAMP or TSAMPMODE_LONGITUDELIST
  double fixedtsamp;          // Sampling time in seconds, assumed to be constant throughout observation. Note this refers to actual durations of the samples stored in the file, not necessarily the sampling time before folding for example. In general, the sampling time is not necessarily the period divided by the number of bins, as not necessarily the whole period is stored. This parameter should be set for both folded and non-folded data. Use the get_tsamp() function to get acual sampling time rather than accessing this information directly.
  double *tsamp_list;         // TSAMPMODE_LONGITUDELIST: Each bin has a specified pulse longitude, in degrees. Stored in the array tsamp_list, which has NrBins elements. This can be used, for instance, for PA swings which are not necessarily regularly sampled. Should be set to NULL if not allocated. You can use teh get_pulse_longitude() function to get the pulse longitude, which also works in fixed sampling mode.

  char tsubMode;               // TSUBMODE_UNKNOWN, TSUBMODE_FIXEDTSUB or TSUBMODE_TSUBLIST
  double *tsub_list;           // Length of each subintegration in seconds (or zero if undefined). If TSUBMODE_FIXEDTSUB, there is only one value. If TSUBMODE_TSUBLIST, each subint has a defined value. Should be set to NULL if not allocated. Use get_tsub(), get_tobs() amd get_mjd_subint() to get duration of subint/observation etc.

  // Information about pointing
  double ra, dec;             /* Position in radians */

  // Information about the frequency coverage
  double bandwidth, centrefreq;   // The centre frequency and the bandwidth, all in MHz. These should always be defined as these are used for plotting etc. and define what (always uniform) frequency range was used to produce a channel. Use these (and the set equivalents) get_centre_freq(), get_bw(), get_channelbw(), do not modify them directly in the code. 
  // Could implement: (2) channel list for obs (3) channel list for each subint. 
  char freqMode;             // Frequency labeling. FREQMODE_UNKNOWN, FREQMODE_UNIFORM or FREQMODE_FREQTABLE
                             // Note this does affect things in the sense that a psrchive fits file can have non-uniform sampling which affects doeing the dedispersion and stuff. This effectively acts as some sort of reference frequency.
  double *freqlabel_list;         // FREQMODE_FREQTABLE: this is list of frequencies for each frequency channel and subint (subint is the outer loop). Otherwise this has no meaning and the frequency label is determined from the bandwidth and centre frequency. This variable should be set to NULL if not allocated. Use get/set_channel_nonweighted/weighted_freq() to get/set the frequency of a particular channel. Do not edit them manually.

  // Data reduction information.
  // The following flag, potentially becomes a number as there could be different options.
  int gentype;         /* See #define GENTYPE_ above */
  char isTransposed;            // This flag is set after calling preprocess_transposeRawFBdata(), i.e. the "subints" are really the frequency channels. Cannot remember exactly what other properties are. Used in pplot to draw frequency along the vertical axis.

  float xrange[2];     /* xrange and yrange defines bin centres of first and last bin. Defined for some gentypes. */
  float yrange[2];
  char xrangeset, yrangeset;
  datafile_history_entry_definition history;  // The history regarding how this data-set was processed

  // EVERYTHING BELOW THIS LINE IS NOT WRITTEN IN PSRSALSA HEADER/FILE YET, AS WANT TO THINK FIRST ABOUT SENSIBLE STRUCTURE

  // The actual data, these are things which are not copied with copy_params_PSRData()
  // DEFINE DATA PER SUBINT? WOULD SIMPLIFY CODE. CAN HAVE SUBINT STRUCTURE.
  float *data;                /* Data is organized as subint, freq, pol, bin */
  float *offpulse_rms;  /* One (off-pulse rms) value for each subint, frequency and polarization channel, ordered in the same way as data (but only one bin), or set to NULL if unused. For POLTYPE_ILVPAdPA the array is long enough to hold 5 polarizations, although only I, L and V are useful. For POLTYPE_ILVPAdPATEldEl type data there are 8 polarizations, and I, L, V and totpol contain useful information. */


  // Not really used? Removed them for now, since it is not clear if/when it will be implemented
  //  float fd_sang, fd_xyph;   /* in degrees */ 

  // Use the following in my code as well (weights)? I assume scales/offsets only need to be known during reading of fits files.
  // ALSO, WEIGHTS COULD BE DIFFERENT IN DIFFERENT SUBINTS. THINK IT SHOULD BE DEFINED FOR EACH SUBINT.
  float *scales, *offsets, *weights;  /* These are used in PSRFITS files. */
  int weight_stats_set;
  int weight_stats_zeroweightfound;
  int weight_stats_differentweights;
  int weight_stats_negativeweights;
  float weight_stats_weightvalue;

  long long datastart;        /* Undefined in PSRFITS, otherwise it points to the first data point in the file after the header */
  // Following is a bit silly. Should be another polarization channels
  //  long NrPApoints;            /* Set if PA-points were calculated */
  //  float *data_pa, *data_dpa;  /* The PA's and errorbars in degrees */

  // Could be useful I guess. Better way? But this is only useful for a pulse profile, so data is small, so can just define all bins, but set data to be undefined (error pa -1 or something).
  //  int bins_defined;           /* Set if the bin numbers of the original file where PA points came from, which are not necessarily all present */
  //  long int *data_bin;         /* The bin numbers */

  // Doesn't it simply follow from samptime if defined? Can I use a function instead? Or is it related to unequal sampling, like data_bin.
  //  int longitudes_defined;     /* Set if the longitudes of the bins are stored */
  //  float *data_long;           /* The longitudes in degrees */

  // Affected by above
  //  int rmspoints_defined; /* Set to the number of polarization channels stored, i.e. I, L, V with ppol for instance */
  //  float *rmsvalues;  /* One (off-pulse rms) value for each subint, frequency and polarization channel */


/*
 After adding a parameter, make sure to update:
     - cleanPSRData
     - copy_paramsPSRData
     - closePSRData         (if data needs to be released)
     - writePSRSALSAHeader  (increase version)
     - readPSRSALSAHeader   (make backward compatible)
 */
}datafile_definition;

//START REGION DEVELOP

typedef struct {
  int phasemodel;                  // 1=linear dep. with freq
  int phasemodel_defined;
  long phaseRefmodel;               // Additional fixed model: 0=none
  long diffgainRefmodel;               // Additional fixed model: 0=none
  int diffgainmodel;               // 1=const, 2 = linear, ....
  int diffgainmodel_defined;
  int gainmodel;                   // 1=const, 2 = linear, ....
  int gainmodel_defined;
  int or1model;                   // 1=const, 2 = linear, ....
  int or1model_defined;
  int or2model;                   // 1=const, 2 = linear, ....
  int or2model_defined;
  int ell1model;                   // 1=const, 2 = linear, ....
  int ell1model_defined;
  int ell2model;                   // 1=const, 2 = linear, ....
  int ell2model_defined;
  double reffreq;                 // freq below is 0.001*(freq[MHz]-reffreq[MHz])
  double phasemodel_values[MaxNrfitReceiverModelFitParameters];     // phase[deg] = [0] + [1]*freq[GHz]
  double phasemodel_redchi2;
  double diffgainmodel_values[MaxNrfitReceiverModelFitParameters];     // phase[deg] = [0] + [1]*freq[GHz]
  double diffgainmodel_redchi2;
  double gainmodel_values[MaxNrfitReceiverModelFitParameters];     // phase[deg] = [0] + [1]*freq[GHz]
  double gainmodel_redchi2;
  double ell1model_values[MaxNrfitReceiverModelFitParameters];     // ellipticity[deg] = [0] + [1]*freq[GHz]
  double ell1model_redchi2;
  double ell2model_values[MaxNrfitReceiverModelFitParameters];     // ellipticity[deg] = [0] + [1]*freq[GHz]
  double ell2model_redchi2;
  double or1model_values[MaxNrfitReceiverModelFitParameters];     // orientation[deg] = [0] + [1]*freq[GHz]
  double or1model_redchi2;
  double or2model_values[MaxNrfitReceiverModelFitParameters];     // orientation[deg] = [0] + [1]*freq[GHz]
  double or2model_redchi2;
  double mjd;
}analyticReceiverSolution_def;
// If changing these definitions, also change clear_analyticReceiverSolution() and analyticReceiverSolution_diffphase() and analyticReceiverSolution_diffgain()  and analyticReceiverSolution_gain() and analyticReceiverSolution_ellipticity() and analyticReceiverSolution_orientation() and printAnalyticReceiverSolution() and writeAnalyticReceiverSolution()

//START REGION DEVELOP

typedef struct {
  long order;
  long nbreak;
  double x_start, x_end;
  double *coefficients;   // Constants obtained from fitting (nrcoefficients = nbreak + order - 2)
  double *breakpoints;    // The position of the breakpoints
}splineFit_def;

//START REGION DEVELOP

/*
y = exp(-ln16*(x-center)^2/fwhm^2)
*/
typedef struct {
  double center;
  double min;
  double max;
  /*  double peak; */
  double fwhm;
  int NrPoints;
}refraction_PeakStruct;

typedef struct {
  double N;              /* relative density of cone, should be 1 for cone 1 */
                        /* so that r = 1 is defined at chi_c of cone 1 */
  double eps1;  
  double eps2;  
  double chi_c;          /* [rad], defined at r = 1 */
  double chi_cn;         /* chi_c/chi_c cone nr 1 */
}refraction_PlasmaConeStruct;

typedef struct {
  double chi_s;                   /* [chi_c of cone 1] at r=1 */
  double chi_e;                   /* [chi_c of cone 1] at r=1 */
  double alpha1;                  /* W = N^alpha1 */
  double f0;                      
  double gamma_f;                 
}refraction_PulsarEmissionStruct;

typedef struct {
  double r;                       /* remember chi_c defined at r = 1 */
  double chi;                     /* [rad] */
  double theta;                   /* [rad] */
  double phi;                     /* [rad] */
  double psi;                     /* [rad] */
  double N;                      /* = 1 at (r=1, chi=chi_c of cone 1) */
  double NormN;                  /* = N/N0 */
  double NormChi;                 /* Normalized chi = chi/chi0 */
  double NormTheta;               /* Normalized theta = theta/theta0 */
  double eta[3];                  /* The three roots of the dispersion law */
  double dlnNdchi;               /* N derivative */
  int NrSolutions;               /* The number of real roots that were found */
  double Position[3];             /* x, y and z  */
  double Direction[3];            /* n_x, n_y and n_z */
}refraction_PositionStruct;

typedef struct {
  int NrCones;
  refraction_PlasmaConeStruct Cone[MaxNrConesRefraction];
  double chi_s;                   /* [chi_c of cone 1] at r=1 */
  double chi_e;                   /* [chi_c of cone 1] at r=1 */
  double alpha1;                  /* W = N^alpha1 */
  refraction_PeakStruct f0;                      
  refraction_PeakStruct chi0;               /* center and min undefined */
  double gamma_f;
  double LongitudeRange;          /* phi = -range ... range */
}refraction_PulsarPhysicalParameters;

typedef struct {  
  int NrCalculatedPoints;         /* Nr calculated points in trajectory */  
  refraction_PositionStruct StartPoint;      /* The startpoint */
  refraction_PositionStruct *Points;         /* The points */
  double f0;
}refraction_RayParameters;

typedef struct {
  int NrPoints;                      /* Nr of requested rays */
  int DisableVariableEmissionHeight; /* 0 = No variable emission height */
  double EmissionHeight;              /* At 1 are chi_s and chi_c defined */
  double *chis;                       /* [chi_c of cone 1] */
  double **chi0;                      /* [rad] */
  double *r0;                         /* At 1 are chi_s and chi_c defined */
  int   *Succesful;                  /* 0 = Ray could not be calculated */
  double ***theta_f;                  /* The final ray directions [rad] */
                                     /* theta_f[nr_f0][nr_chi_0][nr_chi_s] */
  double *f0;
  int Nrf0;  
  int Nrchi0;
  double gamma_f;
}refraction_FinalRayDirection;

typedef struct {
  double alpha;     /* [rad] Angle between magn. axis and the rotation axis */
  double zeta;      /* [rad] Angle between the l.o.s. and the rotation axis */
  int NrPoints;
  double *pl;       /* [rad] Pulse longitude */
  double *phi;      /* [rad] Longitude of ray pointing to Earth */
  double *theta;    /* [rad] Final ray direction of ray pointing to Earth */
}refraction_LineOfSight;

typedef struct {
  int NrPulses;
  int NrPoints;
  double *pl;
  double **Intensity;
}refraction_PulseProfileStack;

typedef struct {
  char NrRays;
  char EquationSet;
  char EmissionHeight;
  char chi_c;
  char chi_s;
  char chi_e;
  char alpha1;
  char f0center;
  char f0min;
  char f0max;
  char f0fwhm;
  char f0NrPoints;
  char chi0max;
  char chi0fwhm;
  char chi0NrPoints;
  char gamma;
  int  NrCones;
  char eps1[MaxNrConesRefraction];
  char eps2[MaxNrConesRefraction];
  char N[MaxNrConesRefraction];
  char chi_cn[MaxNrConesRefraction];
  char alpha;
  char zeta;
  char LongitudeRange;
  char FileNr;
}refraction_InfoReadedStruct;

//START REGION DEVELOP
//START REGION RELEASE

// Never used, except in the definition of fitfunction_collection_type
typedef struct {
  int type;                           /* Type of function */
  double param[MaxNrFitParameters];   /* Parameters defining function */
  double start[MaxNrFitParameters];   /* Start values */
  int fit_flag[MaxNrFitParameters];   /* If set, then value is fitted */
  double value[MaxNrFitParameters];   /* Fitted values */
  double error[MaxNrFitParameters];   /* Fitted value errors */
}fitfunction_type;

// Requires in principle unnecesary non-dynamic memory allocations, but all functions were re-written to use pointers only to avoid too much non-dynamic memory use.
typedef struct {
  int nrfuncs;
  fitfunction_type func[MaxNrFitFunctions];
  double chi2, chi2_red;
}fitfunc_collection_type;


//START REGION DEVELOP




//START REGION DEVELOP
//START REGION RELEASE

typedef struct {
  char progname[MaxFilenameLength], *genusage;
  int switch_verbose, switch_debug, switch_nocounters;
  verbose_definition verbose_state;
  int switch_formatlist;
  int switch_iformat, iformat; 
  int switch_oformat, oformat;
  int switch_header;
  int switch_headerlist;
  int switch_onpulse, switch_onpulsef;
  int switch_polselect, polselectnr;
  int switch_itf, itf;
  int switch_rebin, dorebin; long rebin;
  int switch_nread; long nread;
  int switch_nskip; long nskip;
  int switch_conshift, doconshift;
  int switch_circshift, docircshift;
  int switch_rot, switch_rotdeg, doshiftphase; float shiftPhase_cmdline, shiftPhase;
  int switch_filelist, filelist;
  int switch_device; char pgplotdevice[MaxPgplotDeviceLength];
  int switch_tscr; long dotscr;
  int switch_tscr_complete; int tscr_complete;
  int switch_TSCR, doTSCR;
  int switch_fscr; long dofscr;
  int switch_FSCR, doFSCR;
  int switch_dedisperse, do_dedisperse;
  int switch_deFaraday, do_deFaraday;
//START REGION DEVELOP
  int switch_deFaradayTable, do_deFaradayTable; char *deFaradayTable_filename;
//START REGION RELEASE
  int switch_changeRefFreq; double newRefFreq;
  int switch_stokes, dostokes;
  int switch_coherence, docoherence;
  int switch_noweights, noweights;
  int switch_useweights, useweights;
  int switch_uniformweights, uniformweights;
  int switch_scale, doscale; float scale_scale, scale_offset;
  int switch_debase, dodebase;    // dodebase=1 = use -onpulse, dodebase=2 = use -onpulse2
  int switch_debase_slope, dodebase_slope;  // The dodebase flag is set as well to identify which onpulse option is used.
  int switch_onpulsegr, doonpulsegr;
  int switch_size, windowwidth, windowheight;
  int switch_macro; FILE *macro_ptr;
  int switch_noplotsubset, do_noplotsubset;
  int switch_cmap, cmap;
  int switch_cmaplist;
  int switch_insertparang, switch_deparang, do_parang_corr;
  int switch_history_cmd_only, history_cmd_only;
  int switch_norm, do_norm; float normvalue;
  int switch_normglobal, do_normglobal;
  int switch_clip, do_clip; float clipvalue;
  int switch_fchan, fchan_select;
  int switch_fixseed, fixseed;
  int switch_templatedata, template_data_index; datafile_definition template_file;
  int switch_template, template_specified;
  int switch_align, doalign;
  int switch_blocksize, blocksize;
  pulselongitude_regions_definition onpulse;
  vonMises_collection_definition vonMises_components;
  int switch_ext; char *extension;
  int switch_output; char outputname[MaxFilenameLength];
  int switch_shuffle, doshuffle;
  int switch_rotateStokes; int nr_rotateStokes, rotateStokes1[maxNrRotateStokes], rotateStokes2[maxNrRotateStokes]; float rotateStokesAngle[maxNrRotateStokes];
  int switch_libversions;
//START REGION DEVELOP
// The following is used by getNextFilenameFromList(), so needs to be defined
//START REGION RELEASE
  int doautot;
  int switch_forceUniformFreqLabelling;
  int switch_onpulse2, switch_onpulsef2;
  pulselongitude_regions_definition onpulse2;
//START REGION DEVELOP
  int switch_rotateQU, dorotateQU; float rotateQUangle;  // Replaced by rotateStokes option
  int switch_rotateUV, dorotateUV; float rotateUVangle;  // Replaced by rotateStokes option
  int switch_output2;
  int switch_invertxy, invertxy, switch_invertfp, invertfp, switch_invertf, invertf, switch_invertfx, invertfx;
  int switch_ppgplot;
  int switch_template_onpulse, switch_centert;
  int switch_gate, switch_fftsmooth, switch_smooth, switch_fftzap, switch_fft;
  int switch_autot;
  int switch_normRMS, donormRMS;
  int switch_normRMSi, donormRMSi;
  int switch_testsignal, do_testsignal;
  int switch_absweights;
  int switch_checknan, switch_removenan, dochecknan, doremovenan;
  int switch_checkinf, switch_removeinf, docheckinf, doremoveinf;
  int switch_useweightedfreq, useweightedfreq;  /* Have to add the   if(application.useweightedfreq) psrfits_set_noweights(1); yourself in your program */
  int switch_autoonpulse, do_autoonpulse;
  int switch_showrevision;
  int switch_dedeFaraday;
  int switch_dededisperse;  
  int switch_swapcables, doswapcables;
  int doGate, dofft;
  int template_onpulse, docentert;
  int do_fftzap;
  int absweights;
  int switch_rotslope, dorotslope; float rotslope;
  int switch_subtractRVM, dosubtractRVM; float subtractRVM_alpha, subtractRVM_beta, subtractRVM_l0, subtractRVM_pa0;
  int switch_subtractAvPA, dosubtractAvPA; 
  int switch_zerodming, dozerodming;
  int switch_zerodming_adv, dozerodming_adv, zerodming_adv_mode, zerodming_adv_maxwidth, zerodming_adv_extrabins, zerodming_adv_checkdedispersed, zerodming_adv_minwidening;
  float zerodming_adv_possigma, zerodming_adv_negsigma;
  char *ppgplot_name;
  float do_fftsmooth, do_smooth, fftzap_low, fftzap_high;
//START REGION RELEASE
  int *fzapMask;   /* If != NULL, this array indicates which frequency channels should be ignored during frequency scrunch options. */
}psrsalsaApplication;


//START REGION DEVELOP
typedef struct {
  int defined;       // 0=no, 1=yes
  long double t0;    // Epoch of periastron (MJD)
  long double tasc;  // Epoch of ascending node (MJD)
  long double pb;    // Orbital period (days)
  long double ecc;   // Orbital eccentricity
  int ecc_defined;   // set if defined
  long double pbdot; // The first time derivative of binary period, same unit as tempo2 (s/s?)
  long double a1;    // Projected semimajor axis of orbit (lt-s)
  long double om;    // Longitude of periastron (deg)
  long double omdot; // Periastron advance (deg/yr) 
  long double a1dot; // Rate of change of semimajor axis (lt-s / s) 
  long double gamma; // Post-Keplerian "gamma" term (s) 
  long double eps1, eps2; // ECC x sin(OM) and ECC x cos(OM) for ELL1 model
  int eps1_defined;  // set if defined
}binary_def;
//START REGION DEVELOP

//START REGION DEVELOP
typedef struct {
  // WITH command line flags
  int first_pulse, last_pulse; // The first and last pulses to be plotted (inclusive) [-fp, -lp]
  float alpha, beta; // Parameters used in the cartographic transform [-a, -b]
  float phi0; // Offset values for pulse number and longitude (0 by default) [-phi0]
  double P1, P3_hat; // The carousel rotational period [-p1, -p3h]
  float sampling_period; // Period of sampling, and hence total period of one subint [-sp]
  float maxbrightness; // Set the maximum brightness of the image (useful for the movie) [-mb]
  int SIZE;// Size in pixels of the image to be generated (SIZExSIZE) [-s]
  float R_max; // Maximum value of R to be plotted (corresponds to a corner value) [-rmax]
  int movie; // Creates a movie [-movie]
  int pulsepf; // Pulses per frame for the movie [-ppf]
  int frame_limiter;// Plot every nth frame [-fl]
  int stack; // Plots a pulse stack [-stack]
  int car; // Renders a carousel [-car]
  int p3car; // Renders a carousel from P3-folded data [-p3car]
  int N; // Specify the number of pulses in the carousel [-n]
  
  // WITHOUT command line flags
  float P3; // Tertiary drift period, P3
  int verbose; // These two values are set to 1 automatically when either -v or -debug is used
  int debug;
  psrsalsaApplication *app_pointer; // Pointer to the application. I can't remember why this exists...
  int losk; // Pulse number of the line of sight
  int k0; // k0
  float brightest_pixel; // Record the brightest pixel in the movie, so the user can set the brightness levels second time round.
  int NrBins, NrSubints, NrPols, NrFreqs; // Number of bins, subints, polarisation and frequency channels
  float subint_width; // width of each subint in seconds
  int whichsubint; // Which subint to read from the input (used to distinguish between P3-folded and non-P3-folded data)
  int num_pulses; // How many pulses to plot
  int A, B; // Parameters used in the cartographic transform
}carousel_parameter_list;					
//START REGION RELEASE

//START REGION RELEASE
//START REGION DEVELOP

//START REGION RELEASE
#define PSRSALSA_TYPEDEFS_LOADED 1
#endif
